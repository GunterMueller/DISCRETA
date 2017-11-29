/* unip.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef UNIPOLYTRUE

#include <DISCRETA/unip.h>

#if TEXDOCU
INT UNIPOLY_OB::realloc_z(INT l)
#endif
{
	c_obj_k(VECTOR);
	((VECTOR_OP) this)->realloc_z(l);
	c_obj_k(UNIPOLY);
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::m_il(INT l)
#endif
{
	((VECTOR_OP) this)->m_il(l);
	c_obj_k(UNIPOLY);
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::m_il_n(INT l)
#endif
{
	((VECTOR_OP) this)->m_il_n(l);
	c_obj_k(UNIPOLY);
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::freeself()
#endif
{
	if (s_obj_k() != UNIPOLY) {
		return ((SYM_OP) this)->freeself();
		}
	c_obj_k(VECTOR);
	return ((VECTOR_OP) this)->freeself();
}

#if TEXDOCU
INT UNIPOLY_OB::m_v(VECTOR_OP coeffs)
#endif
{
	coeffs->swap((VECTOR_OP) this);
	c_obj_k(UNIPOLY);
	if (s_li() == 0) {
		printf("UNIPOLY::m_v() warning: vector of length zero\n");
		}
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::degree()
#endif
{
	INT l, i;

	l = s_li();
	for (i = l - 1; i >= 0; i--) {
		if (!s_i(i)->nullp())
			break;
		}
	if (i < l) {
		if (i < 0)
			i = 0;
		realloc_z(i + 1);
		}
	return i;
}

#if TEXDOCU
INT UNIPOLY_OB::nullp()
#endif
{
	INT d;

	d = degree();
	if (d > 0)
		return FALSE;
	if (s_i(0)->nullp())
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT UNIPOLY_OB::einsp()
#endif
{
	INT d;

	d = degree();
	if (d > 0)
		return FALSE;
	if (s_i(0)->einsp())
		return TRUE;
	else
		return FALSE;
}

#if TEXDOCU
INT UNIPOLY_OB::zero()
#endif
{
	realloc_z(1);
	s_i(0)->zero();
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::one()
#endif
{
	realloc_z(1);
	s_i(0)->one();
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::m_one()
#endif
{
	realloc_z(1);
	s_i(0)->m_one();
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::add(UNIPOLY_OP b, UNIPOLY_OP c)
#endif
{
	INT d1, d2, d3, i;
	
	d1 = degree();
	d2 = b->degree();
	d3 = MAXIMUM(d1, d2);
	c->m_il(d3 + 1);
	for (i = 0; i <= d1; i++) {
		if (i <= d2) {
			s_i(i)->add(b->s_i(i), c->s_i(i));
			}
		else {
			s_i(i)->copy(c->s_i(i));
			}
		}
	for (i = d1 + 1; i <= d2; i++) {
		b->s_i(i)->copy(c->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::add_apply(UNIPOLY_OP b)
#endif
{
	UNIPOLY_OB c;
	
	add(b, &c);
	c.swap(b);
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::addinvers(UNIPOLY_OP b)
#endif
{
	INT d1, i;
	
	d1 = degree();
	b->m_il(d1);
	for (i = 0; i <= d1; i++) {
		s_i(i)->addinvers(b->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::addinvers_apply()
#endif
{
	UNIPOLY_OB c;
	
	addinvers(&c);
	swap(&c);
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::mult(UNIPOLY_OP b, UNIPOLY_OP c)
#endif
{
	INT d1, d2, d3, i, j, k;
	SYM_OB x1, x2;
	
	d1 = degree();
	d2 = b->degree();
	d3 = d1 + d2;
	
	c->m_il(d3 + 1);
	s_i(0)->copy(&x1);
	x1.zero();
	for (i = 0; i <= d3; i++) {
		x1.copy(c->s_i(i));
		}
	for (i = 0; i <= d1; i++) {
		for (j = 0; j <= d2; j++) {
			k = i + j;
			s_i(i)->mult(b->s_i(j), &x1);
			x1.add(c->s_i(k), &x2);
			x2.swap(c->s_i(k));
			}
		}
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::mult_apply(UNIPOLY_OP b)
#endif
{
	UNIPOLY_OB c;
	
	mult(b, &c);
	c.swap(b);
	return OK;
}

INT unip_f_print_sub = FALSE;
INT unip_f_use_variable_name = FALSE;
BYTE unip_variable_name[128];

#if TEXDOCU
INT UNIPOLY_OB::sprint(BYTE *str)
#endif
/* appends to str. writes to maximal strlength of 200. */
{
	INT d, i, f_print_k, k, f_prev = FALSE;
	BYTE str1[256];
	SYM_OP coef;
	BYTE *x, *y;
	
	if (unip_f_use_variable_name)
		x = unip_variable_name;
	else
		x = "x";
	if (unip_f_print_sub)
		y = "_";
	else
		y = "^";
	d = degree();
	strcat(str, "(");
	for (i = d; i >= 0; i--) {
		coef = s_i(i);
		str1[0] = 0;
		if (coef->s_obj_k() == INTEGER) {
			k = ((INTEGER_OP) coef)->s_i();
			if (k == 0)
				continue;
			if (k < 0) {
				strcat(str1, " -");
				}
			else if (f_prev) {
				strcat(str1, " +");
				}
			if (k < 0)
				k = -k;
			if (k != 1 || (i == 0 && !unip_f_use_variable_name)) {
				sprintf(str1 + strlen(str1), "%ld", k);
				}
			}
		else {
			f_print_k = TRUE; 
			if (coef->einsp() && i > 0)
				f_print_k = FALSE;
			if (f_prev)
				strcat(str1, " +");
			if (f_print_k) {
				coef->sprint(str1);
				} 
			}
		if (i == 0) {
			if (unip_f_use_variable_name) {
				strcat(str1, x);
				strcat(str1, y);
				strcat(str1, "0");
				}
			}
		else if (i == 1) {
			strcat(str1, x);
			if (unip_f_print_sub) {
				strcat(str1, y);
				strcat(str1, "1");
				}
			}
		else if (i > 1) {
			strcat(str1, x);
			strcat(str1, y);
			sprintf(str1 + strlen(str1), "%ld", i);
			}
		if (strlen(str) + strlen(str1) < 200)
			strcat(str, str1);
		else
			return OK;
		f_prev = TRUE;
		}
	strcat(str, ")");
	return(OK);
}

#if TEXDOCU
INT UNIPOLY_OB::latex(FILE *fp)
#endif
{
	BYTE s[10000];

	s[0] = 0;
	sprint_latex(s);
	fprintf(fp, "%s", s);
	return OK;
}

#if TEXDOCU
INT UNIPOLY_OB::sprint_latex(BYTE *str)
#endif
/* appends to str. writes to maximal strlength of 200. */
{
	INT d, i, f_print_k, k, f_prev = FALSE;
	BYTE str1[256];
	SYM_OP coef;
	BYTE *x, *y;
	
	if (unip_f_use_variable_name)
		x = unip_variable_name;
	else
		x = "x";
	if (unip_f_print_sub)
		y = "_";
	else
		y = "^";
	d = degree();
	// strcat(str, "(");
	for (i = d; i >= 0; i--) {
		coef = s_i(i);
		str1[0] = 0;
		if (coef->s_obj_k() == INTEGER) {
			k = ((INTEGER_OP) coef)->s_i();
			if (k == 0)
				continue;
			if (k < 0) {
				strcat(str1, " -");
				}
			else if (f_prev) {
				strcat(str1, " +");
				}
			if (k < 0)
				k = -k;
			if (k != 1 || (i == 0 && !unip_f_use_variable_name)) {
				sprintf(str1 + strlen(str1), "%ld", k);
				}
			}
		else {
			f_print_k = TRUE; 
			if (coef->einsp() && i > 0)
				f_print_k = FALSE;
			if (f_prev)
				strcat(str1, " +");
			if (f_print_k) {
				coef->sprint_latex(str1);
				} 
			}
		if (i == 0) {
			if (unip_f_use_variable_name) {
				strcat(str1, x);
				strcat(str1, y);
				strcat(str1, "0");
				}
			}
		else if (i == 1) {
			strcat(str1, x);
			if (unip_f_print_sub) {
				strcat(str1, y);
				strcat(str1, "1");
				}
			}
		else if (i > 1) {
			strcat(str1, x);
			strcat(str1, y);
			if (i < 10)
				sprintf(str1 + strlen(str1), "%ld", i);
			else
				sprintf(str1 + strlen(str1), "{%ld}", i);
			}
		if (strlen(str) + strlen(str1) < 200)
			strcat(str, str1);
		else
			return OK;
		f_prev = TRUE;
		}
	// strcat(str, ")");
	return(OK);
}

#endif /* UNIPOLYTRUE */


