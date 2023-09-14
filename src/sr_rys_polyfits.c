// This function is implemented based on libslater library
// https://github.com/nubakery/libslater

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rys_roots.h"
#include "sr_roots_part0_w.dat"
#include "sr_roots_part0_x.dat"
#include "sr_roots_part1_w.dat"
#include "sr_roots_part1_x.dat"
#include "sr_roots_part2_w.dat"
#include "sr_roots_part2_x.dat"
#include "sr_roots_part3_w.dat"
#include "sr_roots_part3_x.dat"

void _CINT_clenshaw_d1(double *rr, const double *x, double u, FINT nroot);
void _CINT_clenshaw_dc(double *rr, const double *x, double u, FINT nroot);
void _CINT_matmul_14_14(double *imc, double *im, FINT nroot);

static int searchsorted(const double *a, double v, int n)
{
        int i;
        int m = n / 2;
        if (a[m] <= v) {
                for (i = m + 1; i < n; i++) {
                        if (v < a[i]) {
                                break;
                        }
                }
        } else {
                for (i = 1; i < m; i++) {
                        if (v < a[i]) {
                                break;
                        }
                }
        }
        return i - 1;
}

int CINTsr_rys_polyfits(int nroots, double x, double lower, double *u, double *w)
{
        double xll = x * lower * lower;
        if (xll >= 38.44) {
                int i;
                for (i = 0; i < nroots; i++) {
                        u[i] = 0.;
                        w[i] = 0.;
                }
                return 0;
        }

        int il, ix;
        double t0, t1, u0, u1;
        double *dx, *dw;
        if (lower < 0.1) {
                il = searchsorted(SR_DATA0_UBASE, lower, 5);
                u0 = SR_DATA0_UBASE[il];
                u1 = SR_DATA0_UBASE[il+1];
                if (x < 12) {
                        ix = searchsorted(SR_DATA0_TBASE, x, 4);
                        dw = SR_DATA0_W + ((nroots-1)*nroots/2*4*3 + nroots*(il*3+ix)) * 196;
                        dx = SR_DATA0_X + ((nroots-1)*nroots/2*4*3 + nroots*(il*3+ix)) * 196;
                        t0 = SR_DATA0_TBASE[ix];
                        t1 = SR_DATA0_TBASE[ix+1];
                } else if (x <= SR_DATA1_TBASE[19]) {
                        ix = searchsorted(SR_DATA1_TBASE, x, 20);
                        dw = SR_DATA1_W + ((nroots-1)*nroots/2*4*19 + nroots*(il*19+ix)) * 196;
                        dx = SR_DATA1_X + ((nroots-1)*nroots/2*4*19 + nroots*(il*19+ix)) * 196;
                        t0 = SR_DATA1_TBASE[ix];
                        t1 = SR_DATA1_TBASE[ix+1];
                } else {
                        return 1;
                }
        } else if (lower < 0.6) {
                il = searchsorted(SR_DATA2_UBASE, lower, 6);
                u0 = SR_DATA2_UBASE[il];
                u1 = SR_DATA2_UBASE[il+1];
                if (x < 196) {
                        ix = searchsorted(SR_DATA2_TBASE, x, 11);
                        dw = SR_DATA2_W + ((nroots-1)*nroots/2*5*10 + nroots*(il*10+ix)) * 196;
                        dx = SR_DATA2_X + ((nroots-1)*nroots/2*5*10 + nroots*(il*10+ix)) * 196;
                        t0 = SR_DATA2_TBASE[ix];
                        t1 = SR_DATA2_TBASE[ix+1];
                } else if (x <= SR_DATAL_TBASE[6]) {
                        x = sqrt(xll);
                        ix = searchsorted(SR_DATAL_TBASE, x, 7);
                        dw = SR_DATAL_W + ((nroots-1)*nroots/2*5*6 + nroots*(il*6+ix)) * 196;
                        dx = SR_DATAL_X + ((nroots-1)*nroots/2*5*6 + nroots*(il*6+ix)) * 196;
                        t0 = SR_DATAL_TBASE[ix];
                        t1 = SR_DATAL_TBASE[ix+1];
                } else {
                        return 1;
                }
        } else {
                return 1;
        }
        double uu = (lower - u0) * 2 / (u1 - u0) - 1;
        double tt = (x - t0) * 2 / (t1 - t0) - 1;
        double im[14*6];
        _CINT_clenshaw_dc(im, dx, uu, nroots);
        _CINT_clenshaw_d1(u, im, tt, nroots);
        _CINT_clenshaw_dc(im, dw, uu, nroots);
        _CINT_clenshaw_d1(w, im, tt, nroots);
        return 0;
}
