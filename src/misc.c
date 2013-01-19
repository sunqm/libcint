/*
 * File: util.c
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 * basic functions
 */

void dcmplx_re(const int n, double *z, const double *re)
{
        int i;
        for (i = 0; i < n; i++) {
                z[0] = re[i];
                z[1] = 0;
                z += 2;
        }
}

void dcmplx_im(const int n, double *z, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[0] = 0;
                z[1] = im[i];
                z += 2;
        }
}

void dcmplx_pp(const int n, double *z, const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[0] = re[i];
                z[1] = im[i];
                z += 2;
        }
}
void dcmplx_pn(const int n, double *z, const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[0] =  re[i];
                z[1] = -im[i];
                z += 2;
        }
}
void dcmplx_np(const int n, double *z, const double *re, const double *im)
{
        int i;
        for (i = 0; i < n; i++) {
                z[0] = -re[i];
                z[1] =  im[i];
                z += 2;
        }
}


inline double square_dist(const double *r1, const double *r2)
{
        double r12[3];

        r12[0] = r1[0] - r2[0];
        r12[1] = r1[1] - r2[1];
        r12[2] = r1[2] - r2[2];

        return r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
}

