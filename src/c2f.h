/*
 * Copyright (C) 2013  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include "config.h"

#define __CVAR_CALL__   atm, *natm, bas, *nbas, env
#define __FVAR_FUNC__   const FINT *atm, const FINT *natm, \
                        const FINT *bas, const FINT *nbas, const double *env
#define C2F_(NAME)      FINT NAME##_(double *op, const FINT *shls, \
                                    __FVAR_FUNC__) \
                        { return NAME(op, shls, __CVAR_CALL__); }
// C2Fo for 2e integrals with optimizer
#define C2Fo_(NAME)     FINT NAME##_(double *op, const FINT *shls, \
                                    __FVAR_FUNC__, unsigned long *optptr) \
                        { CINTOpt *opt = (CINTOpt *)*optptr; \
                            return NAME(op, shls, __CVAR_CALL__, opt); }

typedef long CINTOptPtrAsInteger8;
#define OPTIMIZER2F_(NAME) \
    void NAME##_(CINTOptPtrAsInteger8 *optptr, __FVAR_FUNC__) { \
        CINTOpt **opt = (CINTOpt **)optptr; \
        NAME(opt, __CVAR_CALL__); }
