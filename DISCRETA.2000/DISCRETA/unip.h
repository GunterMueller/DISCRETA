/* unip.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef UNIP_INCLUDED
#define UNIP_INCLUDED

#ifndef GFQ_INCLUDED
#include <DISCRETA/gfq.h>
#endif

class unipoly_ob : public VECTOR_OB {
public:
	INT realloc_z(INT l);
	INT m_il(INT l);
	INT m_il_n(INT l);
	INT freeself();
	INT m_v(VECTOR_OP coeffs);
	INT degree();
	INT nullp();
	INT einsp();
	INT zero();
	INT one();
	INT m_one();
	INT add(UNIPOLY_OP b, UNIPOLY_OP c);
	INT add_apply(UNIPOLY_OP b);
	INT addinvers(UNIPOLY_OP b);
	INT addinvers_apply();
	INT mult(UNIPOLY_OP b, UNIPOLY_OP c);
	INT mult_apply(UNIPOLY_OP b);
	INT sprint(BYTE *str);
	INT latex(FILE *fp);
	INT sprint_latex(BYTE *str);
};

extern INT unip_f_print_sub;
extern INT unip_f_use_variable_name;
extern BYTE unip_variable_name[128];

#endif /* UNIP_INCLUDED */
