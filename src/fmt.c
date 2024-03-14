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

#include <float.h>
#include <math.h>
#include "cint_config.h"
#include "rys_roots.h"
#define SML_FLOAT64   (DBL_EPSILON * .5)
#define SML_FLOAT80   2.0e-20
#define SQRTPIE4      .8862269254527580136490837416705725913987747280611935641069038949264
#define SQRTPIE4l     .8862269254527580136490837416705725913987747280611935641069038949264l
#define ERFC_bound    200

#ifdef HAVE_QUADMATH_H
#include <quadmath.h>
#define SQRTPIE4q     .8862269254527580136490837416705725913987747280611935641069038949264q
#define SML_FLOAT128  1.0e-35
#endif

/*
 * Relative errors of fmt1_erfc_like are of
 *      (2*t)**(m-1) / (2m-3)!! * machine_precision * fmt_val
 * Errors of the other choice are
 *      (2m-1)!! / (2*t)**(m-1) * machine_precision * fmt_val
 * Given m, the turn-over point for t should satisfy
 *      (2m-1)!! / (2*t)**(m-1) > (2m-1)**.5
 * t0 = .5 * ((2m-1)!!/(2m-1)**.5)**(1/(m-1))
 */
static double TURNOVER_POINT[] = {
        0.,
        0.,
        0.866025403784,
        1.295010032056,
        1.705493613097,
        2.106432965305,
        2.501471934009,
        2.892473348218,
        3.280525047072,
        3.666320693281,
        4.05033123037 ,
        4.432891808508,
        4.814249856864,
        5.194593501454,
        5.574069276051,
        5.952793645111,
        6.330860773135,
        6.708347923415,
        7.08531930745 ,
        7.461828891625,
        7.837922483937,
        8.213639312398,
        8.589013237349,
        8.964073695432,
        9.338846443746,
        9.713354153046,
        10.08761688545,
        10.46165248270,
        10.83547688448,
        11.20910439128,
        11.58254788331,
        11.95581900374,
        12.32892831326,
        12.70188542111,
        13.07469909673,
        13.44737736550,
        13.81992759110,
        14.19235654675,
        14.56467047710,
        14.93687515212
};

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
static void fmt1_gamma_inc_like(double *f, double t, int m)
{
        int i;
        double b = m + 0.5;
        double bi;
        double e = .5 * exp(-t);
        double x = e;
        double s = e;
        double tol = SML_FLOAT64 * e;
        for (bi = b + 1.; x > tol; bi += 1.) {
                x *= t / bi;
                s += x;
        }
        f[m] = s / b;
        for (i = m; i > 0; i--) {
                b -= 1.;
                f[i-1] = (e + t * f[i]) / b;
        }
}

void gamma_inc_like(double *f, double t, int m)
{
        if (t < TURNOVER_POINT[m]) {
                fmt1_gamma_inc_like(f, t, m);
        } else {
                int i;
                double tt = sqrt(t);
                f[0] = SQRTPIE4 / tt * erf(tt);
                if (m > 0) {
                        double e = exp(-t);
                        double b = .5 / t;
                        for (i = 1; i <= m; i++)
                                f[i] = b * ((2*i-1) * f[i-1] - e);
                }
        }
}

static void fmt1_lgamma_inc_like(long double *f, long double t, int m)
{
        long double b = m + 0.5l;
        long double bi;
        long double e = .5l * expl(-t);
        long double x = e;
        long double s = e;
        long double tol = SML_FLOAT80 * e;
        int i;
        for (bi = b + 1.; x > tol; bi += 1.) {
                x *= t / bi;
                s += x;
        }
        f[m] = s / b;
        for (i = m; i > 0; i--) {
                b -= 1;
                f[i-1] = (e + t * f[i]) / b;
        }
}

void lgamma_inc_like(long double *f, long double t, int m)
{
        if (t < TURNOVER_POINT[m]) {
                fmt1_lgamma_inc_like(f, t, m);
        } else {
                int i;
                long double tt = sqrtl(t);
                f[0] = SQRTPIE4l / tt * erfl(tt);
                if (m > 0) {
                        long double e = expl(-t);
                        long double b = .5l / t;
                        for (i = 1; i <= m; i++)
                                f[i] = b * ((2*i-1) * f[i-1] - e);
                }
        }
}

static inline double _pow(double base, int exponent)
{
        int i;
        double result = 1;
        for (i = 1; i <= exponent; i <<= 1) {
                if (i & exponent) {
                        result *= base;
                }
                base *= base;
        }
        return result;
}

static inline long double _powl(long double base, int exponent)
{
        int i;
        long double result = 1.l;
        for (i = 1; i <= exponent; i <<= 1) {
                if (i & exponent) {
                        result *= base;
                }
                base *= base;
        }
        return result;
}

/* This function evaluates the auxiliary integral,
 *
 *     2  _ 1           2
 *  t s  /     2 m  -t u
 * e     |    u    e      du,
 *      _/  s
 *
 * by a power series expansion
 *
 * F[m] = e^{t s^2} int_l^1 u^{2m} e^{-t u^2} du
 *      = e^{t s^2} /(2m+1) int e^{-t u^2} d u^{2m+1}
 *      = e^{t s^2} /(2m+1) [e^{-t u^2} u^{2m+1}]_l^1 + (2t)/(2m+1) int u^{2m+2} e^{-t u^2} du
 *      = e^{t s^2} /(m+.5) (.5*e^{-t} - .5*e^{-t l^2} l^{2m+1}) + t F[m+1])
 */
void fmt1_erfc_like(double *f, double t, double lower, int m)
{
        int i;
        double lower2 = lower * lower;
        double b = m + 0.5;
        double bi;
        double e = .5 * exp(-t);
        double e1 = .5 * exp(-t * lower2) * lower;
        e1 *= _pow(lower2, m);
        double x = e;
        double x1 = e1;
        double s = e - e1;
        double div = 1.;
        double delta = s;
        double tol = SML_FLOAT64 * fabs(delta);
        for (bi = b + 1.; fabs(delta) > tol; bi += 1.) {
                div *= t / bi;
                x1 *= lower2;
                delta = (x - x1) * div;
                s += delta;
        }
        double val = s / b;
        f[m] = val;
        for (i = m; i > 0; i--) {
                b -= 1.;
                e1 /= lower2;
                val = (e - e1 + t * val) / b;
                f[i-1] = val;
        }
}
void fmt_erfc_like(double *f, double t, double lower, int m)
{
        if (lower == 0) {
                return gamma_inc_like(f, t, m);
        }

        int i;
        double lower2 = lower * lower;
        // F[m] < .5*sqrt(pi/t) * erfc(low*tt)
        if (t * lower2 > ERFC_bound) {
                for (i = 0; i <= m; i++) {
                        f[i] = 0;
                }
                return;
        }

        if (t < TURNOVER_POINT[m]) {
                fmt1_erfc_like(f, t, lower, m);
        } else {
                double tt = sqrt(t);
                // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
                double val = SQRTPIE4 / tt * (erfc(lower * tt) - erfc(tt));
                f[0] = val;
                if (m > 0) {
                        double e = exp(-t);
                        double e1 = exp(-t * lower2) * lower;
                        double b = .5 / t;
                        for (i = 0; i < m; i++) {
                                val = b * ((2*i+1) * val - e + e1);
                                e1 *= lower2;
                                f[i+1] = val;
                        }
                }
        }
}

void fmt_lerfc_like(long double *f, long double t, long double lower, int m)
{
        if (lower == 0) {
                return lgamma_inc_like(f, t, m);
        }

        int i;
        long double lower2 = lower * lower;
        // F[m] < .5*sqrt(pi/t) * erfc(low*tt)
        if (t * lower2 > ERFC_bound) {
                for (i = 0; i <= m; i++) {
                        f[i] = 0;
                }
                return;
        }

        if (t < TURNOVER_POINT[m]) {
                fmt1_lerfc_like(f, t, lower, m);
        } else {
                long double tt = sqrtl(t);
                // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
                long double val = SQRTPIE4l / tt * (erfcl(lower * tt) - erfcl(tt));
                f[0] = val;
                if (m > 0) {
                        long double e = expl(-t);
                        long double e1 = expl(-t * lower2) * lower;
                        long double b = .5l / t;
                        for (i = 0; i < m; i++) {
                                val = b * ((2*i+1) * val - e + e1);
                                e1 *= lower2;
                                f[i+1] = val;
                        }
                }
        }
}

void fmt1_lerfc_like(long double *f, long double t, long double lower, int m)
{
        int i;
        long double lower2 = lower * lower;
        long double b = m + 0.5l;
        long double bi;
        long double e = .5l * expl(-t);
        long double e1 = .5l * expl(-t * lower2) * lower;
        e1 *= _powl(lower2, m);
        long double x = e;
        long double x1 = e1;
        long double s = e - e1;
        long double div = 1.l;
        long double delta = s;
        long double tol = SML_FLOAT80 * fabsl(delta);
        for (bi = b + 1.l; fabsl(delta) > tol; bi += 1.l) {
                div *= t / bi;
                x1 *= lower2;
                delta = (x - x1) * div;
                s += delta;
        }
        long double val = s / b;
        f[m] = val;
        for (i = m; i > 0; i--) {
                b -= 1.l;
                e1 /= lower2;
                val = (e - e1 + t * val) / b;
                f[i-1] = val;
        }
}

#ifdef HAVE_QUADMATH_H
static void fmt1_qgamma_inc_like(__float128 *f, __float128 t, int m)
{
        __float128 b = m + .5q;
        __float128 bi;
        __float128 e = .5q * expq(-t);
        __float128 x = e;
        __float128 s = e;
        __float128 tol = SML_FLOAT128 * e;
        int i;
        for (bi = b + 1.; x > tol; bi += 1.) {
                x *= t / bi;
                s += x;
        }
        f[m] = s / b;
        for (i = m; i > 0; i--) {
                b -= 1;
                f[i-1] = (e + t * f[i]) / b;
        }
}

void qgamma_inc_like(__float128 *f, __float128 t, int m)
{
        if (t < TURNOVER_POINT[m]) {
                fmt1_qgamma_inc_like(f, t, m);
        } else {
                int i;
                __float128 tt = sqrtq(t);
                f[0] = SQRTPIE4q / tt * erfq(tt);
                if (m > 0) {
                        __float128 e = expq(-t);
                        __float128 b = .5q / t;
                        for (i = 1; i <= m; i++)
                                f[i] = b * ((2*i-1) * f[i-1] - e);
                }
        }
}

static inline __float128 _powq(__float128 base, int exponent)
{
        int i;
        __float128 result = 1.q;
        for (i = 1; i <= exponent; i <<= 1) {
                if (i & exponent) {
                        result *= base;
                }
                base *= base;
        }
        return result;
}

void fmt_qerfc_like(__float128 *f, __float128 t, __float128 lower, int m)
{
        if (lower == 0) {
                return qgamma_inc_like(f, t, m);
        }

        int i;
        __float128 lower2 = lower * lower;
        // F[m] < .5*sqrt(pi/t) * erfc(low*tt)
        if (t * lower2 > ERFC_bound) {
                for (i = 0; i <= m; i++) {
                        f[i] = 0;
                }
                return;
        }

        if (t < TURNOVER_POINT[m]) {
                fmt1_qerfc_like(f, t, lower, m);
        } else {
                __float128 tt = sqrtq(t);
                // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
                __float128 val = SQRTPIE4q / tt * (erfcq(lower * tt) - erfcq(tt));
                f[0] = val;
                if (m > 0) {
                        __float128 e = expq(-t);
                        __float128 e1 = expq(-t * lower2) * lower;
                        __float128 b = .5q / t;
                        for (i = 0; i < m; i++) {
                                val = b * ((2*i+1) * val - e + e1);
                                e1 *= lower2;
                                f[i+1] = val;
                        }
                }
        }
}

void fmt1_qerfc_like(__float128 *f, __float128 t, __float128 lower, int m)
{
        int i;
        __float128 lower2 = lower * lower;
        __float128 b = m + .5q;
        __float128 bi;
        __float128 e = .5q * expq(-t);
        __float128 e1 = .5q * expq(-t * lower2) * lower;
        e1 *= _powq(lower2, m);
        __float128 x = e;
        __float128 x1 = e1;
        __float128 s = e - e1;
        __float128 div = 1.q;
        __float128 delta = s;
        __float128 tol = SML_FLOAT128 * fabsq(delta);
        for (bi = b + 1.q; fabsq(delta) > tol; bi += 1.q) {
                div *= t / bi;
                x1 *= lower2;
                delta = (x - x1) * div;
                s += delta;
        }
        __float128 val = s / b;
        f[m] = val;
        for (i = m; i > 0; i--) {
                b -= 1.q;
                e1 /= lower2;
                val = (e - e1 + t * val) / b;
                f[i-1] = val;
        }
}
#endif
