"""
>.\hiir-designer.exe 1 0.005 140
Number of x2 stages: 1
Stopband attenuation: 140 dB

Coefficients for phase 1 with transition bandwidth 0.0050000000000000001 (Total 19 coefficients)
{ 0.019911761024506557, 0.076569065603139905, 0.16170648261075027, 0.26428227031893498, 0.37320978687920564, 0.47939467893641913, 0.57665589850082322, 0.66168172238942402, 0.7334355636406803, 0.79240315662949701, 0.83992271287611509, 0.87769279111118159, 0.90746017802851253, 0.93085009866291657, 0.94929377019349725, 0.96401566368781932, 0.97605397317065279, 0.98629782872833549, 0.99553233211505254 }

real: static constexpr std::array<T, 10> co{T(0.019911761024506557),T(0.16170648261075027),T(0.37320978687920564),T(0.5766558985008232),T(0.7334355636406803),T(0.8399227128761151),T(0.9074601780285125),T(0.9492937701934973),T(0.9760539731706528),T(0.9955323321150525)};
imag: static constexpr std::array<T, 9> co{T(0.0765690656031399),T(0.264282270318935),T(0.4793946789364191),T(0.661681722389424),T(0.792403156629497),T(0.8776927911111816),T(0.9308500986629166),T(0.9640156636878193),T(0.9862978287283355)};
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

from numpy.polynomial import polynomial

def plot_tf(b, a):
    w, h = signal.freqz(b, a, worN=2048)
    fig, ax1 = plt.subplots()
    ax1.set_title('Digital filter frequency response')
    ax1.plot(w, 20 * np.log10(abs(h)), 'b')
    ax1.set_ylabel('Amplitude [dB]', color='b')
    ax1.set_xlabel('Frequency [rad/sample]')
    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(h))
    ax2.plot(w, angles, 'g')
    ax2.set_ylabel('Angle (radians)', color='g')
    ax2.grid()
    ax2.axis('tight')
    plt.show()

def add_delay(sos):
    return np.vstack((sos, [0, 1, 0, 1, 0, 0]))

def formatCoefficients(co):
    co_str = ",".join([f"T({str(c)})" for c in co])
    return f"static constexpr std::array<T, {len(co)}> co{{{co_str}}};"

def plot():
    co = [
        0.019911761024506557, 0.076569065603139905, 0.16170648261075027,
        0.26428227031893498, 0.37320978687920564, 0.47939467893641913,
        0.57665589850082322, 0.66168172238942402, 0.7334355636406803, 0.79240315662949701,
        0.83992271287611509, 0.87769279111118159, 0.90746017802851253,
        0.93085009866291657, 0.94929377019349725, 0.96401566368781932,
        0.97605397317065279, 0.98629782872833549, 0.99553233211505254
    ]

    def section(a):
        return [a, 0, 1, 1, 0, a]

    print("real: " + formatCoefficients(co[0::2]))
    print("imag: " + formatCoefficients(co[1::2]))

    sos_real = np.array([section(a) for a in co[0::2]])
    sos_imag = np.array([section(a) for a in co[1::2]])
    sos_imag = add_delay(sos_imag)

    def toPoly(sos):
        b, a = signal.sos2tf(sos)
        return (polynomial.Polynomial(b[::-1]), polynomial.Polynomial(a[::-1]))

    num_real, den_real = toPoly(sos_real)
    num_imag, den_imag = toPoly(sos_imag)

    num_poly = num_real * den_imag + num_imag * den_real
    den_poly = 2 * den_real * den_imag

    plot_tf(num_poly.coef[::-1], den_poly.coef[::-1])

if __name__ == "__main__":
    plot()
