/*
This is a port of some functions in cephes.

- https://netlib.org/cephes/

Below is the original readme which mentions the license.

```
   Some software in this archive may be from the book _Methods and
Programs for Mathematical Functions_ (Prentice-Hall or Simon & Schuster
International, 1989) or from the Cephes Mathematical Library, a
commercial product. In either event, it is copyrighted by the author.
What you see here may be used freely but it comes with no support or
guarantee.

   The two known misprints in the book are repaired here in the
source listings for the gamma function and the incomplete beta
integral.


   Stephen L. Moshier
   moshier@na-net.ornl.gov
```

Note that this port is using math functions from standard library. This means that:

- the output value may not match to the original cephes implementation.
- the accuracy may degrade.
*/

#pragma once
#include <cmath>
#include <limits>

namespace cephes {

constexpr double MAXLOG = 7.09782712893383996732E2;   /* log(MAXNUM) */
constexpr double PI = 3.14159265358979323846;         /* pi */
constexpr double PIO2 = 1.57079632679489661923;       /* pi/2 */
constexpr double MAXNUM = 1.79769313486231570815E308; /* 2**1024*(1-epsilon) */

inline double polevl(double x, const double coef[], int N)
{
  double ans;
  int i;
  const double *p;

  p = coef;
  ans = *p++;
  i = N;

  do ans = ans * x + *p++;
  while (--i);

  return (ans);
}

inline double p1evl(double x, const double coef[], int N)
{
  double ans;
  const double *p;
  int i;

  p = coef;
  ans = x + *p++;
  i = N - 1;

  do ans = ans * x + *p++;
  while (--i);

  return (ans);
}

double igam(double a, double x);

// Lower incomplete gamma function. Both arguments must be positive.
double igamc(double a, double x)
{
  constexpr double big = 4.503599627370496e15;
  constexpr double biginv = 2.22044604925031308085e-16;
  double ans, ax, c, yc, r, t, y, z;
  double pk, pkm1, pkm2, qk, qkm1, qkm2;

  if ((x <= 0) || (a <= 0)) return (1.0);

  if ((x < 1.0) || (x < a)) return (1.0 - igam(a, x));

  ax = a * std::log(x) - x - std::lgamma(a);
  if (ax < -MAXLOG) {
    return (0.0);
  }
  ax = std::exp(ax);

  y = 1.0 - a;
  z = x + y + 1.0;
  c = 0.0;
  pkm2 = 1.0;
  qkm2 = x;
  pkm1 = x + 1.0;
  qkm1 = z * x;
  ans = pkm1 / qkm1;

  do {
    c += 1.0;
    y += 1.0;
    z += 2.0;
    yc = y * c;
    pk = pkm1 * z - pkm2 * yc;
    qk = qkm1 * z - qkm2 * yc;
    if (qk != 0) {
      r = pk / qk;
      t = std::fabs((ans - r) / r);
      ans = r;
    } else
      t = 1.0;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if (std::fabs(pk) > big) {
      pkm2 *= biginv;
      pkm1 *= biginv;
      qkm2 *= biginv;
      qkm1 *= biginv;
    }
  } while (t > std::numeric_limits<double>::epsilon());

  return (ans * ax);
}

// Upper incomplete gamma function. Both arguments must be positive.
double igam(double a, double x)
{
  double ans, ax, c, r;

  if ((x <= 0) || (a <= 0)) return (0.0);
  if ((x > 1.0) && (x > a)) return (1.0 - igamc(a, x));

  ax = a * std::log(x) - x - std::lgamma(a);
  if (ax < -MAXLOG) return (0.0);
  ax = std::exp(ax);

  r = a;
  c = 1.0;
  ans = 1.0;

  do {
    r += 1.0;
    c *= x / r;
    ans += c;
  } while (c / ans > std::numeric_limits<double>::epsilon());

  return (ans * ax / a);
}

namespace SICI {
constexpr double EUL = 0.57721566490153286061;
constexpr double SN[] = {
  -8.39167827910303881427E-11, 4.62591714427012837309E-8,  -9.75759303843632795789E-6,
  9.76945438170435310816E-4,   -4.13470316229406538752E-2, 1.00000000000000000302E0,
};
constexpr double SD[] = {
  2.03269266195951942049E-12, 1.27997891179943299903E-9, 4.41827842801218905784E-7,
  9.96412122043875552487E-5,  1.42085239326149893930E-2, 9.99999999999999996984E-1,
};
constexpr double CN[] = {
  2.02524002389102268789E-11, -1.35249504915790756375E-8, 3.59325051419993077021E-6,
  -4.74007206873407909465E-4, 2.89159652607555242092E-2,  -1.00000000000000000080E0,
};
constexpr double CD[] = {
  4.07746040061880559506E-12, 3.06780997581887812692E-9, 1.23210355685883423679E-6,
  3.17442024775032769882E-4,  5.10028056236446052392E-2, 4.00000000000000000080E0,
};
constexpr double FN4[] = {
  4.23612862892216586994E0,  5.45937717161812843388E0,  1.62083287701538329132E0,
  1.67006611831323023771E-1, 6.81020132472518137426E-3, 1.08936580650328664411E-4,
  5.48900223421373614008E-7,
};
constexpr double FD4[] = {
  8.16496634205391016773E0,  7.30828822505564552187E0,  1.86792257950184183883E0,
  1.78792052963149907262E-1, 7.01710668322789753610E-3, 1.10034357153915731354E-4,
  5.48900252756255700982E-7,
};
constexpr double FN8[] = {
  4.55880873470465315206E-1, 7.13715274100146711374E-1,  1.60300158222319456320E-1,
  1.16064229408124407915E-2, 3.49556442447859055605E-4,  4.86215430826454749482E-6,
  3.20092790091004902806E-8, 9.41779576128512936592E-11, 9.70507110881952024631E-14,
};
constexpr double FD8[] = {
  9.17463611873684053703E-1,  1.78685545332074536321E-1,  1.22253594771971293032E-2,
  3.58696481881851580297E-4,  4.92435064317881464393E-6,  3.21956939101046018377E-8,
  9.43720590350276732376E-11, 9.70507110881952025725E-14,
};
constexpr double GN4[] = {
  8.71001698973114191777E-2, 6.11379109952219284151E-1, 3.97180296392337498885E-1,
  7.48527737628469092119E-2, 5.38868681462177273157E-3, 1.61999794598934024525E-4,
  1.97963874140963632189E-6, 7.82579040744090311069E-9,
};
constexpr double GD4[] = {
  1.64402202413355338886E0,  6.66296701268987968381E-1, 9.88771761277688796203E-2,
  6.22396345441768420760E-3, 1.73221081474177119497E-4, 2.02659182086343991969E-6,
  7.82579218933534490868E-9,
};
constexpr double GN8[] = {
  6.97359953443276214934E-1, 3.30410979305632063225E-1,  3.84878767649974295920E-2,
  1.71718239052347903558E-3, 3.48941165502279436777E-5,  3.47131167084116673800E-7,
  1.70404452782044526189E-9, 3.85945925430276600453E-12, 3.14040098946363334640E-15,
};
constexpr double GD8[] = {
  1.68548898811011640017E0,  4.87852258695304967486E-1,  4.67913194259625806320E-2,
  1.90284426674399523638E-3, 3.68475504442561108162E-5,  3.57043223443740838771E-7,
  1.72693748966316146736E-9, 3.87830166023954706752E-12, 3.14040098946363335242E-15,
};
} // namespace SICI

// Sine and cosine integral.
int sici(double x, double *si, double *ci)
{
  double z, c, s, f, g;
  int sign;

  if (x < 0.0) {
    sign = -1;
    x = -x;
  } else
    sign = 0;

  if (x == 0.0) {
    *si = 0.0;
    *ci = -std::numeric_limits<double>::infinity();
    return (0);
  }

  if (x > 1.0e9) {
    if (std::isfinite(x)) {
      *si = PIO2 - std::cos(x) / x;
      *ci = std::sin(x) / x;
    } else {
      if (sign == -1) {
        *si = -PIO2;
        *ci = NAN;
      } else {
        *si = PIO2;
        *ci = 0.0;
      }
    }
    return (0);
  }

  if (x <= 4.0) {
    z = x * x;
    s = x * polevl(z, SICI::SN, 5) / polevl(z, SICI::SD, 5);
    c = z * polevl(z, SICI::CN, 5) / polevl(z, SICI::CD, 5);

    if (sign) s = -s;
    *si = s;
    *ci = SICI::EUL + std::log(x) + c;
    return (0);
  }

  s = std::sin(x);
  c = std::cos(x);
  z = 1.0 / (x * x);
  if (x < 8.0) {
    f = polevl(z, SICI::FN4, 6) / (x * p1evl(z, SICI::FD4, 7));
    g = z * polevl(z, SICI::GN4, 7) / p1evl(z, SICI::GD4, 7);
  } else {
    f = polevl(z, SICI::FN8, 8) / (x * p1evl(z, SICI::FD8, 8));
    g = z * polevl(z, SICI::GN8, 8) / p1evl(z, SICI::GD8, 9);
  }
  *si = PIO2 - f * c - g * s;
  if (sign) *si = -(*si);
  *ci = f * s - g * c;

  return (0);
}

// Sine integral.
double Si(double x)
{
  int sign;
  if (x < 0.0) {
    sign = -1;
    x = -x;
  } else
    sign = 0;

  if (x == 0.0) return 0.0;

  if (x > 1.0e9) {
    if (std::isfinite(x)) return PIO2 - std::cos(x) / x;
    if (sign == -1) return -PIO2;
    return PIO2;
  }

  if (x <= 4.0) {
    double z = x * x;
    double s = x * polevl(z, SICI::SN, 5) / polevl(z, SICI::SD, 5);
    return sign ? -s : s;
  }

  double z = 1.0 / (x * x);
  double f, g;
  if (x < 8.0) {
    f = polevl(z, SICI::FN4, 6) / (x * p1evl(z, SICI::FD4, 7));
    g = z * polevl(z, SICI::GN4, 7) / p1evl(z, SICI::GD4, 7);
  } else {
    f = polevl(z, SICI::FN8, 8) / (x * p1evl(z, SICI::FD8, 8));
    g = z * polevl(z, SICI::GN8, 8) / p1evl(z, SICI::GD8, 9);
  }
  double si = PIO2 - f * std::cos(x) - g * std::sin(x);
  return sign ? -si : si;
}

// Cosine integral.
double Ci(double x)
{
  int sign;
  if (x < 0.0) {
    sign = -1;
    x = -x;
  } else
    sign = 0;

  if (x == 0.0) return -std::numeric_limits<double>::infinity();

  if (x > 1.0e9) {
    if (std::isfinite(x)) return std::sin(x) / x;
    if (sign == -1) return NAN;
    return 0.0;
  }

  if (x <= 4.0) {
    const double z = x * x;
    const double c = z * polevl(z, SICI::CN, 5) / polevl(z, SICI::CD, 5);
    return SICI::EUL + std::log(x) + c;
  }

  const double z = 1.0 / (x * x);
  double f, g;
  if (x < 8.0) {
    f = polevl(z, SICI::FN4, 6) / (x * p1evl(z, SICI::FD4, 7));
    g = z * polevl(z, SICI::GN4, 7) / p1evl(z, SICI::GD4, 7);
  } else {
    f = polevl(z, SICI::FN8, 8) / (x * p1evl(z, SICI::FD8, 8));
    g = z * polevl(z, SICI::GN8, 8) / p1evl(z, SICI::GD8, 9);
  }
  return f * std::sin(x) - g * std::cos(x);
}

double spence(double x)
{
  constexpr double A[8] = {
    4.65128586073990045278E-5, 7.31589045238094711071E-3, 1.33847639578309018650E-1,
    8.79691311754530315341E-1, 2.71149851196553469920E0,  4.25697156008121755724E0,
    3.29771340985225106936E0,  1.00000000000000000126E0,
  };
  constexpr double B[8] = {
    6.90990488912553276999E-4, 2.54043763932544379113E-2, 2.82974860602568089943E-1,
    1.41172597751831069617E0,  3.63800533345137075418E0,  5.03278880143316990390E0,
    3.54771340985225096217E0,  9.99999999999999998740E-1,
  };
  double w, y, z;
  int flag;

  if (x < 0.0) return (0.0);
  if (x == 1.0) return (0.0);
  if (x == 0.0) return (PI * PI / 6.0);

  flag = 0;

  if (x > 2.0) {
    x = 1.0 / x;
    flag |= 2;
  }

  if (x > 1.5) {
    w = (1.0 / x) - 1.0;
    flag |= 2;
  } else if (x < 0.5) {
    w = -x;
    flag |= 1;
  } else
    w = x - 1.0;

  y = -w * polevl(w, A, 7) / polevl(w, B, 7);

  if (flag & 1) y = (PI * PI) / 6.0 - std::log(x) * std::log(1.0 - x) - y;

  if (flag & 2) {
    z = std::log(x);
    y = -0.5 * z * z - y;
  }

  return (y);
}

} // namespace  cephes
