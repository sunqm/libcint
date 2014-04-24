/*
 * File: c2f.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#define __CVAR_CALL__   atm, *natm, bas, *nbas, env
#define __FVAR_FUNC__   const int *atm, const int *natm, \
                        const int *bas, const int *nbas, const double *env
#define C2F_(NAME)      int NAME##_(double *op, const unsigned int *shls, \
                                    __FVAR_FUNC__) \
                        { return NAME(op, shls, __CVAR_CALL__); }
// C2Fo for 2e integrals with optimizer
#define C2Fo_(NAME)     int NAME##_(double *op, const unsigned int *shls, \
                                    __FVAR_FUNC__, unsigned long *optptr) \
                        { CINTOpt *opt = (CINTOpt *)*optptr; \
                            return NAME(op, shls, __CVAR_CALL__, opt); }

typedef long CINTOptPtrAsInteger8;
#define OPTIMIZER2F_(NAME) \
    void NAME##_(CINTOptPtrAsInteger8 *optptr, __FVAR_FUNC__) { \
        CINTOpt **opt = (CINTOpt **)optptr; \
        NAME(opt, __CVAR_CALL__); }
