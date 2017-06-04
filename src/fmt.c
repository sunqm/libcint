/*
 * This code computes incomplete gamma function.  It is based on Xin
 * Wu's implementation.
 *
 *
 * List of Abbreviation(s)
 *
 * THO:
 * Gaussian-Expansion Methods for Molecular Integrals,
 * Hiroshi Taketa, Sigeru Huzinaga, and Kiyosi O-ohata,
 * Journal of the Physical Society of Japan,
 * Vol. 21, No. 11, 1966, 2313 - 2324.
 *
 */

#include <math.h>
#define SML_FLOAT64   1.0e-16
#define SML_FLOAT80   1.0e-20
#define SQRTPIE4      .8862269254527580136490837416705725913987747280611935641069038949264
#define SQRTPIE4l     .8862269254527580136490837416705725913987747280611935641069038949264l

#ifdef HAVE_QUADMATH_H
#include <quadmath.h>
#define SQRTPIE4q     .8862269254527580136490837416705725913987747280611935641069038949264q
#define SML_FLOAT128  1.0e-34
#endif

/*
 * Name
 *
 * fmtpse
 *
 * Synopsis
 *
 * double fmtpse(int m, double t)
 *
 * Description
 *
 * This function evaluates the auxiliary integral,
 *
 *             _ 1           2
 *            /     2 m  -t u
 * F (t)  =   |    u    e      du,
 *  m        _/  0
 *
 * by a power series expansion
 *
 *                    _                    2                     3                _
 *           exp(-t) |       t            t                     t                  |
 * F (t)  =  ------- | 1 + ----- + --------------- + ----------------------- + ... |
 *  m          2 b   |_    b + 1   (b + 1) (b + 2)   (b + 1) (b + 2) (b + 3)      _|,
 *
 * where b = m + 1 / 2. This power series expansion converges fast, when t is less than b + 1,
 * namely t < m + 3 / 2.
 *
 * Argument(s)
 *
 * int m:
 * F_m(t), see the Description section.
 *
 * double t:
 * F_m(t), see the Description section.
 *
 * Return Value
 * double:
 * F_m(t), see the Description section.
 *
 */

/*
 * Name
 *
 * fmt
 *
 * Synopsis
 *
 * double fmt(int m, double t)
 *
 * Description
 *
 * This function evaluates the auxiliary integral, see Eq. 2.11 in THO,
 *
 *             _ 1           2
 *            /     2 m  -t u
 * F (t)  =   |    u    e      du,
 *  m        _/  0
 *
 * where m replaces ν for more convenient typesetting.
 *
 * If t is less than SML16 or equals 0, then
 *
 *              1
 * F (t)  =  -------.
 *  m        2 m + 1
 *
 * If t is less than m + 3 / 2, the auxiliary integral is evaluated by
 * a power series expansion (see fmtpse.c for details).
 *
 * Otherwise F (t) is calculated first
 *            0
 *                    1
 *                    -
 *           1 /  π  \2       _
 * F (t)  =  - | --- |  erf( /t ).
 *  0        2 \  t  /
 *
 * Then the upward recurrence relation is used for F (t) of higher m
 *                                                  m
 *
 *            (2 m - 1) F     (t) - exp( -t )
 *                       m - 1
 *  F (t)  =  -------------------------------.
 *   m                      2 t
 *
 * Argument(s)
 *
 * int m:
 * F_m(t), see the Description section.
 *
 * double t:
 * F_m(t), see the Description section.
 *
 * Return Value
 *
 * double:
 * F_m(t), see the Description section.
 *
 */
void gamma_inc_like(double *f, double t, FINT m)
{
        FINT i;
        if (t < m + 1.5) {
                double b = m + 0.5;
                double x = 1;
                double s = 1;
                double e = .5 * exp(-t);
                if (t < SML_FLOAT64) {
                        f[m] = .5 / b;
                } else {
                        //f[m] = fmtpse(m, t);
                        for (i = 1; x > SML_FLOAT64; i++)
                        {
                                x *= t / (b + i);
                                s += x;
                        }
                        f[m] = e * s / b;
                }
                if (m > 0) {
                        for (i = m; i > 0; i--) {
                                b -= 1;
                                f[i-1] = (e + t * f[i]) / b;
                        }
                }
        } else {
                double pi2 = SQRTPIE4;
                double tt = sqrt(t);
                f[0] = pi2 / tt * erf(tt);
                if (m > 0) {
                        double e = exp(-t);
                        double b = .5 / t;
                        for (i = 1; i <= m; i++)
                                f[i] = b * ((2*i-1) * f[i-1] - e);
                }
        }
}

void lgamma_inc_like(long double *f, long double t, FINT m)
{
        FINT i;
        if (t < m + 1.5) {
                long double b = m + 0.5l;
                long double x = 1;
                long double s = 1;
                long double e = .5l * expl(-t);
                if (t < SML_FLOAT80) {
                        f[m] = .5l / b;
                } else {
                        //f[m] = fmtpse(m, t);
                        for (i = 1; x > SML_FLOAT80; i++)
                        {
                                x *= t / (b + i);
                                s += x;
                        }
                        f[m] = e * s / b;
                }
                if (m > 0) {
                        for (i = m; i > 0; i--) {
                                b -= 1;
                                f[i-1] = (e + t * f[i]) / b;
                        }
                }
        } else {
                long double pi2 = SQRTPIE4l;
                long double tt = sqrtl(t);
                f[0] = pi2 / tt * erfl(tt);
                if (m > 0) {
                        long double e = expl(-t);
                        long double b = .5l / t;
                        for (i = 1; i <= m; i++)
                                f[i] = b * ((2*i-1) * f[i-1] - e);
                }
        }
}

#ifdef HAVE_QUADMATH_H
void qgamma_inc_like(__float128 *f, __float128 t, FINT m)
{
        FINT i;
        if (t < m + 1.5) {
                __float128 b = m + .5q;
                __float128 x = 1;
                __float128 s = 1;
                __float128 e = .5q * expq(-t);
                if (t < SML_FLOAT128) {
                        f[m] = .5q / b;
                } else {
                        //f[m] = fmtpse(m, t);
                        for (i = 1; x > SML_FLOAT128; i++)
                        {
                                x *= t / (b + i);
                                s += x;
                        }
                        f[m] = e * s / b;
                }
                if (m > 0) {
                        for (i = m; i > 0; i--) {
                                b -= 1;
                                f[i-1] = (e + t * f[i]) / b;
                        }
                }
        } else {
                __float128 pi2 = SQRTPIE4q;
                __float128 tt = sqrtq(t);
                f[0] = pi2 / tt * erfq(tt);
                if (m > 0) {
                        __float128 e = expq(-t);
                        __float128 b = .5q / t;
                        for (i = 1; i <= m; i++)
                                f[i] = b * ((2*i-1) * f[i-1] - e);
                }
        }
}
#endif
