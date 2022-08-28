import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
    """
    Bracket the minimum of the function.

    Given a function and distinct initial points, search in the
    downhill direction (as defined by the initial points) and return
    new points xa, xb, xc that bracket the minimum of the function
    f(xa) > f(xb) < f(xc). It doesn't always mean that obtained
    solution will satisfy xa<=x<=xb.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function to minimize.
    xa, xb : float, optional
        Bracketing interval. Defaults `xa` to 0.0, and `xb` to 1.0.
    args : tuple, optional
        Additional arguments (if present), passed to `func`.
    grow_limit : float, optional
        Maximum grow limit.  Defaults to 110.0
    maxiter : int, optional
        Maximum number of iterations to perform. Defaults to 1000.

    Returns
    -------
    xa, xb, xc : float
        Bracket.
    fa, fb, fc : float
        Objective function values in bracket.
    funcalls : int
        Number of function evaluations made.

    Examples
    --------
    This function can find a downward convex region of a function:

    >>> import matplotlib.pyplot as plt
    >>> from scipy.optimize import bracket
    >>> def f(x):
    ...     return 10*x**2 + 3*x + 5
    >>> x = np.linspace(-2, 2)
    >>> y = f(x)
    >>> init_xa, init_xb = 0, 1
    >>> xa, xb, xc, fa, fb, fc, funcalls = bracket(f, xa=init_xa, xb=init_xb)
    >>> plt.axvline(x=init_xa, color="k", linestyle="--")
    >>> plt.axvline(x=init_xb, color="k", linestyle="--")
    >>> plt.plot(x, y, "-k")
    >>> plt.plot(xa, fa, "bx")
    >>> plt.plot(xb, fb, "rx")
    >>> plt.plot(xc, fc, "bx")
    >>> plt.show()

    """
    _gold = 1.618034  # golden ratio: (1.0+sqrt(5.0))/2.0
    _verysmall_num = 1e-21
    fa = func(*(xa,) + args)
    fb = func(*(xb,) + args)
    if (fa < fb):  # Switch so fa > fb
        xa, xb = xb, xa
        fa, fb = fb, fa
    xc = xb + _gold * (xb - xa)
    fc = func(*((xc,) + args))
    funcalls = 3
    iter = 0

    while (fc < fb):
        tmp1 = (xb - xa) * (fb - fc)
        tmp2 = (xb - xc) * (fb - fa)
        val = tmp2 - tmp1
        if np.abs(val) < _verysmall_num:
            denom = 2.0 * _verysmall_num
        else:
            denom = 2.0 * val
        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom
        wlim = xb + grow_limit * (xc - xb)
        if iter > maxiter:
            raise RuntimeError("Too many iterations.")
        iter += 1
        if (w - xc) * (xb - w) > 0.0:
            fw = func(*((w,) + args))
            funcalls += 1
            if (fw < fc):
                xa = xb
                xb = w
                fa = fb
                fb = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            elif (fw > fb):
                xc = w
                fc = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            w = xc + _gold * (xc - xb)
            fw = func(*((w,) + args))
            funcalls += 1
        elif (w - wlim) * (wlim - xc) >= 0.0:
            w = wlim
            fw = func(*((w,) + args))
            funcalls += 1
        elif (w - wlim) * (xc - w) > 0.0:
            fw = func(*((w,) + args))
            funcalls += 1
            if (fw < fc):
                xb = xc
                xc = w
                w = xc + _gold * (xc - xb)
                fb = fc
                fc = fw
                fw = func(*((w,) + args))
                funcalls += 1
        else:
            w = xc + _gold * (xc - xb)
            fw = func(*((w,) + args))
            funcalls += 1
        xa = xb
        xb = xc
        xc = w
        fa = fb
        fb = fc
        fc = fw

        print(iter, xa, xb, xc, fa, fb, fc)
    return xa, xb, xc, fa, fb, fc, funcalls

class Brent:
    #need to rethink design of __init__
    def __init__(self, func, args=(), tol=1.48e-8, maxiter=500, full_output=0, disp=0):
        self.func = func
        self.args = args
        self.tol = tol
        self.maxiter = maxiter
        self._mintol = 1.0e-11
        self._cg = 0.3819660
        self.xmin = None
        self.fval = None
        self.iter = 0
        self.funcalls = 0
        self.disp = disp

    # need to rethink design of set_bracket (new options, etc.)
    def set_bracket(self, brack=None):
        self.brack = brack

    def get_bracket_info(self):
        #set up
        func = self.func
        args = self.args
        brack = self.brack
        ### BEGIN core bracket_info code ###
        ### carefully DOCUMENT any CHANGES in core ##
        if brack is None:
            xa, xb, xc, fa, fb, fc, funcalls = bracket(func, args=args)
        elif len(brack) == 2:
            xa, xb, xc, fa, fb, fc, funcalls = bracket(func,
                                                       xa=brack[0],
                                                       xb=brack[1],
                                                       args=args)
        elif len(brack) == 3:
            xa, xb, xc = brack
            if (xa > xc):  # swap so xa < xc can be assumed
                xc, xa = xa, xc
            if not ((xa < xb) and (xb < xc)):
                raise ValueError("Not a bracketing interval.")
            fa = func(*((xa,) + args))
            fb = func(*((xb,) + args))
            fc = func(*((xc,) + args))
            if not ((fb < fa) and (fb < fc)):
                raise ValueError("Not a bracketing interval.")
            funcalls = 3
        else:
            raise ValueError("Bracketing interval must be " "length 2 or 3 sequence.")
        ### END core bracket_info code ###

        return xa, xb, xc, fa, fb, fc, funcalls

    def optimize(self):
        # set up for optimization
        func = self.func
        xa, xb, xc, fa, fb, fc, funcalls = self.get_bracket_info()

        _mintol = self._mintol
        _cg = self._cg
        #################################
        #BEGIN CORE ALGORITHM
        #################################
        x = w = v = xb
        fw = fv = fx = func(*((x,) + self.args))
        if (xa < xc):
            a = xa
            b = xc
        else:
            a = xc
            b = xa
        deltax = 0.0
        funcalls += 1
        iter = 0

        if self.disp > 2:
            print(" ")
            print(f"{'Func-count':^12} {'x':^12} {'f(x)': ^12}")
            print(f"{funcalls:^12g} {x:^12.6g} {fx:^12.6g}")

        while (iter < self.maxiter):
            tol1 = self.tol * np.abs(x) + _mintol
            tol2 = 2.0 * tol1
            xmid = 0.5 * (a + b)
            # check for convergence
            if np.abs(x - xmid) < (tol2 - 0.5 * (b - a)):
                break
            # XXX In the first iteration, rat is only bound in the true case
            # of this conditional. This used to cause an UnboundLocalError
            # (gh-4140). It should be set before the if (but to what?).
            if (np.abs(deltax) <= tol1):
                if (x >= xmid):
                    deltax = a - x  # do a golden section step
                else:
                    deltax = b - x
                rat = _cg * deltax
            else:  # do a parabolic step
                tmp1 = (x - w) * (fx - fv)
                tmp2 = (x - v) * (fx - fw)
                p = (x - v) * tmp2 - (x - w) * tmp1
                tmp2 = 2.0 * (tmp2 - tmp1)
                if (tmp2 > 0.0):
                    p = -p
                tmp2 = np.abs(tmp2)
                dx_temp = deltax
                deltax = rat
                # check parabolic fit
                if ((p > tmp2 * (a - x)) and (p < tmp2 * (b - x)) and
                    (np.abs(p) < np.abs(0.5 * tmp2 * dx_temp))):
                    rat = p * 1.0 / tmp2  # if parabolic step is useful.
                    u = x + rat
                    if ((u - a) < tol2 or (b - u) < tol2):
                        if xmid - x >= 0:
                            rat = tol1
                        else:
                            rat = -tol1
                else:
                    if (x >= xmid):
                        deltax = a - x  # if it's not do a golden section step
                    else:
                        deltax = b - x
                    rat = _cg * deltax

            if (np.abs(rat) < tol1):  # update by at least tol1
                if rat >= 0:
                    u = x + tol1
                else:
                    u = x - tol1
            else:
                u = x + rat
            fu = func(*((u,) + self.args))  # calculate new output value
            funcalls += 1

            if (fu > fx):  # if it's bigger than current
                if (u < x):
                    a = u
                else:
                    b = u
                if (fu <= fw) or (w == x):
                    v = w
                    w = u
                    fv = fw
                    fw = fu
                elif (fu <= fv) or (v == x) or (v == w):
                    v = u
                    fv = fu
            else:
                if (u >= x):
                    a = x
                else:
                    b = x
                v = w
                w = x
                x = u
                fv = fw
                fw = fx
                fx = fu

            if self.disp > 2:
                print(f"{funcalls:^12g} {x:^12.6g} {fx:^12.6g}")

            iter += 1
        #################################
        #END CORE ALGORITHM
        #################################

        self.xmin = x
        self.fval = fx
        self.iter = iter
        self.funcalls = funcalls

    def get_result(self, full_output=False):
        if full_output:
            return self.xmin, self.fval, self.iter, self.funcalls
        else:
            return self.xmin

def doubleEmaEnvelopeD0(n, k_A, k_D):
    A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
    D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
    return (1 - A) * D

def doubleEmaEnvelopeD0Negative(n, k_A, k_D):
    A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
    D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
    return (A - 1) * D

def samplesToKp(timeInSamples):
    y = 1 - np.cos(2 * np.pi / timeInSamples)
    return -y + np.sqrt(y * (y + 2))

def processBrent(func, brack=None, args=(), tol=1.48e-8, maxiter=500, disp=0):
    brent = Brent(func=func,
                  args=args,
                  tol=tol,
                  full_output=True,
                  maxiter=maxiter,
                  disp=disp)
    brent.set_bracket(brack)
    brent.optimize()
    x, fval, nit, nfev = brent.get_result(full_output=True)
    return (x, fval, nit, nfev)

length = 10000
n = np.arange(length)
sampleA = 1100000
sampleD = 1000000
k_A = samplesToKp(sampleA)
k_D = samplesToKp(sampleD)

d0 = doubleEmaEnvelopeD0(n, k_A, k_D)

result = optimize.minimize_scalar(
    doubleEmaEnvelopeD0Negative,
    args=(k_A, k_D),
)
peakTime = result.x
peakValue = -result.fun

x, fval, nit, nfev = processBrent(
    doubleEmaEnvelopeD0Negative,
    args=(k_A, k_D),
)

print(f"Iteration            : {result.nit}")
print(f"Value                : {peakValue}")
print(f"Peak   Time [sample] : {peakTime}")
print(f"Attack Time [sample] : {sampleA}")
print(f"Decay  Time [sample] : {sampleD}")
print(f"Max                  : {np.max(d0)}")
print("--- debug")
print(f"Diff Iteration       : {result.nit - nit}")
print(f"Diff Peak Time       : {peakTime - x}")
print(f"Diff Value           : {peakTime + fval}")
print(f"nfev                 : {nfev}")

d0 /= peakValue

plt.plot(d0, lw=1, alpha=0.5, color="blue", label="d0")
plt.axvline(peakTime, lw=1, ls="--", color="black", alpha=0.5, label="peakTime")
plt.grid(color="#f0f0f0")
plt.legend()
plt.show()
