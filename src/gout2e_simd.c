#include <immintrin.h>
#include "cint.h"

void CINTgout2e(double *gout, double *g, FINT *idx,
                CINTEnvVars *envs, FINT gout_empty)
{
        int nf = envs->nf;
        int nrys_roots = envs->nrys_roots;
        int i, ix, iy, iz, n;
        int jx, jy, jz;
        __m128d r0, r1, r2, r3;
        //_mm256_zeroupper();

        if (gout_empty) {
                switch (nrys_roots) {
                case 1:
                        for (n = 0; n < nf; n++) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                gout[n] = g[ix] * g[iy] * g[iz];
                        }
                        break;
                case 2:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r3 = _mm_mul_pd  (r0, r1);
                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);

                                r3 = _mm_hadd_pd (r3, r0);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r0 = _mm_hadd_pd (r0, r0);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 3:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r3 = _mm_mul_pd  (r0, r1);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r3 = _mm_hadd_pd (r3, r0);

                                r0 = _mm_loadl_pd(r0, g+ix+2);
                                r0 = _mm_loadh_pd(r0, g+jx+2);
                                r1 = _mm_loadl_pd(r1, g+iy+2);
                                r1 = _mm_loadh_pd(r1, g+jy+2);
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadl_pd(r1, g+iz+2);
                                r1 = _mm_loadh_pd(r1, g+jz+2);
                                r0 = _mm_mul_pd  (r0, r1);
                                r3 = _mm_add_pd  (r0, r3);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r0 = _mm_hadd_pd (r0, r0);
                                _mm_store_sd(gout+n, r0);
                                gout[n] += g[ix+2] * g[iy+2] * g[iz+2];
                        }
                        break;
                case 4:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);

                                r3 = _mm_hadd_pd (r3, r0);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 5:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r3 = _mm_hadd_pd (r3, r0);

                                r0 = _mm_loadl_pd(r0, g+ix+4);
                                r0 = _mm_loadh_pd(r0, g+jx+4);
                                r1 = _mm_loadl_pd(r1, g+iy+4);
                                r1 = _mm_loadh_pd(r1, g+jy+4);
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadl_pd(r1, g+iz+4);
                                r1 = _mm_loadh_pd(r1, g+jz+4);
                                r0 = _mm_mul_pd  (r0, r1);
                                r3 = _mm_add_pd  (r0, r3);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                _mm_store_sd(gout+n, r0);
                                gout[n] += g[ix+4] * g[iy+4] * g[iz+4];
                        }
                        break;
                case 6:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+jx+4);
                                r1 = _mm_loadu_pd(g+jy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);

                                r3 = _mm_hadd_pd (r3, r0);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 7:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+jx+4);
                                r1 = _mm_loadu_pd(g+jy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r3 = _mm_hadd_pd (r3, r0);

                                r0 = _mm_loadl_pd(r0, g+ix+6);
                                r0 = _mm_loadh_pd(r0, g+jx+6);
                                r1 = _mm_loadl_pd(r1, g+iy+6);
                                r1 = _mm_loadh_pd(r1, g+jy+6);
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadl_pd(r1, g+iz+6);
                                r1 = _mm_loadh_pd(r1, g+jz+6);
                                r0 = _mm_mul_pd  (r0, r1);
                                r3 = _mm_add_pd  (r0, r3);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                _mm_store_sd(gout+n, r0);
                                gout[n] += g[ix+6] * g[iy+6] * g[iz+6];
                        }
                        break;
                default:
                        for (n = 0; n < nf; n++) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+6);
                                r1 = _mm_loadu_pd(g+iy+6);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+6);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);
                                for (i = 8; i < nrys_roots-1; i+=2) {
                                        r0 = _mm_loadu_pd(g+ix+i);
                                        r1 = _mm_loadu_pd(g+iy+i);
                                        r0 = _mm_mul_pd  (r0, r1);
                                        r1 = _mm_loadu_pd(g+iz+i);
                                        r0 = _mm_mul_pd  (r0, r1);
                                        r3 = _mm_add_pd  (r3, r0);
                                }
                                r3 = _mm_hadd_pd(r3, r3);
                                _mm_store_sd(gout+n, r3);
                                if (i < nrys_roots) {
                                        gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                        }
                        break;
                } // end switch nroots
        } else {
                switch (nrys_roots) {
                case 1:
                        for (n = 0; n < nf; n++) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                gout[n] += g[ix] * g[iy] * g[iz];
                        }
                        break;
                case 2:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r3 = _mm_mul_pd  (r0, r1);
                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_hadd_pd (r3, r0);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r0 = _mm_hadd_pd (r0, r0);
                                r2 = _mm_load_sd (gout+n);
                                r0 = _mm_add_sd  (r0, r2);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 3:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r3 = _mm_mul_pd  (r0, r1);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r3 = _mm_hadd_pd (r3, r0);

                                r0 = _mm_loadl_pd(r0, g+ix+2);
                                r0 = _mm_loadh_pd(r0, g+jx+2);
                                r1 = _mm_loadl_pd(r1, g+iy+2);
                                r1 = _mm_loadh_pd(r1, g+jy+2);
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadl_pd(r1, g+iz+2);
                                r1 = _mm_loadh_pd(r1, g+jz+2);
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_add_pd  (r0, r3);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r0 = _mm_hadd_pd (r0, r0);
                                r2 = _mm_load_sd (gout+n);
                                r0 = _mm_add_sd  (r0, r2);
                                r2 = _mm_set_sd  (g[ix+2] * g[iy+2] * g[iz+2]);
                                r0 = _mm_add_sd  (r0, r2);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 4:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_hadd_pd (r3, r0);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                r2 = _mm_load_sd (gout+n);
                                r0 = _mm_add_sd  (r0, r2);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 5:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r3 = _mm_hadd_pd (r3, r0);

                                r0 = _mm_loadl_pd(r0, g+ix+4);
                                r0 = _mm_loadh_pd(r0, g+jx+4);
                                r1 = _mm_loadl_pd(r1, g+iy+4);
                                r1 = _mm_loadh_pd(r1, g+jy+4);
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadl_pd(r1, g+iz+4);
                                r1 = _mm_loadh_pd(r1, g+jz+4);
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_add_pd  (r0, r3);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                r2 = _mm_load_sd (gout+n);
                                r0 = _mm_add_sd  (r0, r2);
                                r2 = _mm_set_sd  (g[ix+4] * g[iy+4] * g[iz+4]);
                                r0 = _mm_add_sd  (r0, r2);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 6:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+jx+4);
                                r1 = _mm_loadu_pd(g+jy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_hadd_pd (r3, r0);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                r2 = _mm_load_sd (gout+n);
                                r0 = _mm_add_sd  (r0, r2);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                case 7:
                        for (n = 0; n < nf-1; n+=2) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                jx = idx[3+n*3];
                                jy = idx[4+n*3];
                                jz = idx[5+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);

                                r0 = _mm_loadu_pd(g+jx  );
                                r1 = _mm_loadu_pd(g+jy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+jz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+jx+2);
                                r1 = _mm_loadu_pd(g+jy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+jx+4);
                                r1 = _mm_loadu_pd(g+jy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+jz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r3 = _mm_hadd_pd (r3, r0);

                                r0 = _mm_loadl_pd(r0, g+ix+6);
                                r0 = _mm_loadh_pd(r0, g+jx+6);
                                r1 = _mm_loadl_pd(r1, g+iy+6);
                                r1 = _mm_loadh_pd(r1, g+jy+6);
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadl_pd(r1, g+iz+6);
                                r1 = _mm_loadh_pd(r1, g+jz+6);
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_add_pd  (r0, r3);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_storeu_pd(gout+n, r3);
                        }
                        if (n < nf) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r0 = _mm_hadd_pd (r0, r0);
                                r2 = _mm_load_sd (gout+n);
                                r0 = _mm_add_sd  (r0, r2);
                                r2 = _mm_set_sd  (g[ix+6] * g[iy+6] * g[iz+6]);
                                r0 = _mm_add_sd  (r0, r2);
                                _mm_store_sd(gout+n, r0);
                        }
                        break;
                default:
                        for (n = 0; n < nf; n++) {
                                ix = idx[0+n*3];
                                iy = idx[1+n*3];
                                iz = idx[2+n*3];
                                r0 = _mm_loadu_pd(g+ix  );
                                r1 = _mm_loadu_pd(g+iy  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r1 = _mm_loadu_pd(g+iz  );
                                r0 = _mm_mul_pd  (r0, r1);
                                r2 = _mm_loadu_pd(g+ix+2);
                                r1 = _mm_loadu_pd(g+iy+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+2);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+4);
                                r1 = _mm_loadu_pd(g+iy+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+4);
                                r2 = _mm_mul_pd  (r2, r1);
                                r0 = _mm_add_pd  (r0, r2);
                                r2 = _mm_loadu_pd(g+ix+6);
                                r1 = _mm_loadu_pd(g+iy+6);
                                r2 = _mm_mul_pd  (r2, r1);
                                r1 = _mm_loadu_pd(g+iz+6);
                                r2 = _mm_mul_pd  (r2, r1);
                                r3 = _mm_add_pd  (r0, r2);
                                for (i = 8; i < nrys_roots-1; i+=2) {
                                        r0 = _mm_loadu_pd(g+ix+i);
                                        r1 = _mm_loadu_pd(g+iy+i);
                                        r0 = _mm_mul_pd  (r0, r1);
                                        r1 = _mm_loadu_pd(g+iz+i);
                                        r0 = _mm_mul_pd  (r0, r1);
                                        r3 = _mm_add_pd  (r3, r0);
                                }
                                r2 = _mm_loadu_pd(gout+n);
                                r3 = _mm_hadd_pd (r3, r3);
                                r3 = _mm_add_pd  (r3, r2);
                                _mm_store_sd(gout+n, r3);
                                if (i < nrys_roots) {
                                        gout[n] += g[ix+i] * g[iy+i] * g[iz+i];
                                }
                        }
                        break;
                } // end switch nroots
        }
}
