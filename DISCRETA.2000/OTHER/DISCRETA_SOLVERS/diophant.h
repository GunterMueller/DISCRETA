/*4:*/
#line 100 "diophant.w"

#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#if defined(MPREC)
#include "freelip/lip.h"
#undef BLAS
extern long diophant(verylong**a_input,verylong*b_input,verylong*upperbounds_input,
int no_columns,int no_rows,
verylong factor_input,verylong norm_input,
int silent,int iterate,int iterate_no,
int bkz_beta_input,int bkz_p_input,
int stop_after_sol_input,int stop_after_loops_input,
int free_RHS_input,int*org_col_input,int no_org_col_input);
#else
extern long diophant(long**a_input,long*b_input,long*upperbounds_input,
int no_columns,int no_rows,
long factor_input,long norm_input,
int silent,int iterate,int iterate_no,
int bkz_beta_input,int bkz_p_input,
int stop_after_sol_input,int stop_after_loops_input,
int free_RHS_input,int*org_col_input,int no_org_col_input);
#endif
#endif

/*:4*/
