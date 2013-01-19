/*
 * File: c2f.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#define __CVAR_CALL__   shls, atm, *natm, bas, *nbas, env
#define __FVAR_FUNC__   const int *shls, const int *atm, const int *natm, \
                        const int *bas, const int *nbas, const double *env
#define C2F_(NAME)      void NAME##_(double *op, __FVAR_FUNC__) \
                        { NAME(op, __CVAR_CALL__); }
