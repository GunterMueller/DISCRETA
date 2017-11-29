/* poly.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#ifndef POLY_INCLUDED
#define POLY_INCLUDED

#ifndef LIST_INCLUDED
#include <DISCRETA/list.h>
#endif

#ifdef MONOMTRUE
class monom_ob : public SYM_OB {
public:
	SYM_OP s_k() { 
		return (SYM_OP)ob_self.ob_monom->mo_koeff; };
	INTEGER_OP s_k_i() { 
		return (INTEGER_OP)s_k(); };
	INT s_k_ii() { return s_k_i()->s_i(); };
	void c_k(SYM_OP koeff) { 
		ob_self.ob_monom->mo_koeff = (OP)koeff; };
	void m_k_i(INT koeff) { s_k_i()->m_i(koeff); };
	void c_k_i(INT koeff) { s_k_i()->c_i(koeff); };

	VECTOR_OP s_s() { 
		return (VECTOR_OP)ob_self.ob_monom->mo_self; };
	void c_s(VECTOR_OP self) { 
		ob_self.ob_monom->mo_self = (OP)self; };
	INTEGER_OP s_s_l() { return s_s()->s_l(); };
	INT s_s_li() { return s_s()->s_li(); };
	SYM_OP s_s_i(INT i) { return s_s()->s_i(i); };
	INTEGER_OP s_s_ii(INT i) { 
		return (INTEGER_OP)s_s_i(i); };
	INT s_s_iii(INT i) { return s_s_ii(i)->s_i(); };
	void c_s_i(INT i, SYM_OP p) { *s_s_i(i) = *p; };
	void c_s_ii(INT i, INT j) { s_s_ii(i)->c_i(j); };
	void m_s_ii(INT i, INT j) { s_s_ii(i)->m_i(j); };

INT m_sk(VECTOR_OP self, SYM_OP koeff);
INT m_ski(VECTOR_OP self, INT koeff);
INT m_siki(INT self, INT koeff);
INT b_sk(VECTOR_OP self, SYM_OP koeff);
/* build_self koeff_monom */
/* not do be mixed up with SYM_OB::b_ks() ! */
INT b_ski(VECTOR_OP self, INT koeff);
INT fprint(FILE *f);
INT sprint(BYTE *str);
/* appends to str. writes to maximal strlength of 200. */
INT freeself();
INT compare(MONOM_OP b);
/* Bei gleichheit von self wird koeff verglichen */ 
INT copy(MONOM_OP b);
INT mult_scalar(SYM_OP a, MONOM_OP c);
/* a ist skalar */
/* before: mult_scalar_monom(OP a, OP b, OP c) */
INT mult_apply_scalar(SYM_OP a);
/* before: mult_apply_scalar_monom(OP a, OP b) */
INT add(SYM_OP b, MONOM_OP c);
INT addinvers(MONOM_OP b);
INT addinvers_apply();
};
INT add_koeff(SYM_OP a, SYM_OP b);
/* eqhandle bei insert bei polynomen, monopoly.
 * Addiert die beiden koeffizienten der monome a und b
 * nach dem monom b */
INT comp_monomvector_monomvector(SYM_OP a, SYM_OP b);
/* vergleicht zwei monome von vectorform */
/* wenn keine vectorform dann der normale vergleich */
#endif /* MONOMTRUE */

#ifdef POLYTRUE
class polynom_ob : public SYM_OB {
public:
	MONOM_OP s_mo() { 
		return (MONOM_OP)((LIST_OP)this)->s_s(); };
	void c_mo(MONOM_OP mo) { 
		((LIST_OP)this)->c_s(mo); };
	POLYNOM_OP s_n() { 
		return (POLYNOM_OP)((LIST_OP)this)->s_n(); };
	void c_n(POLYNOM_OP next) { 
		((LIST_OP)this)->c_n((LIST_OP)next); };
	
	SYM_OP s_k() { return s_mo()->s_k(); };
	INTEGER_OP s_k_i() { return (INTEGER_OP)s_k(); };
	INT s_k_ii() { return s_k_i()->s_i(); };
	void c_k(SYM_OP koeff) { s_mo()->c_k(koeff); };
	void m_k_i(INT koeff) { s_k_i()->m_i(koeff); };
	void c_k_i(INT koeff) { s_k_i()->c_i(koeff); };

	VECTOR_OP s_s() { return s_mo()->s_s(); };
	void c_s(VECTOR_OP self) { s_mo()->c_s(self); };
	INTEGER_OP s_s_l() { return s_s()->s_l(); };
	INT s_s_li() { return s_s()->s_li(); };
	SYM_OP s_s_i(INT i) { return s_s()->s_i(i); };
	INTEGER_OP s_s_ii(INT i) { return (INTEGER_OP)s_s_i(i); };
	INT s_s_iii(INT i) { return s_s_ii(i)->s_i(); };
	void c_s_i(INT i, SYM_OP p) { *s_s_i(i) = *p; };
	void c_s_ii(INT i, INT j) { s_s_ii(i)->c_i(j); };
	void m_s_ii(INT i, INT j) { s_s_ii(i)->m_i(j); };
	
	INT copy(POLYNOM_OP to) { 
		return ((LIST_OP)this)->copy((LIST_OP)to); };
	INT freeself() { 
		return ((LIST_OP)this)->freeself(); };
	
INT Print();
INT b_mn(MONOM_OP m, POLYNOM_OP next);
INT b_skn(VECTOR_OP self, SYM_OP koeff, POLYNOM_OP next);
INT b_skin(VECTOR_OP self, INT koeff, POLYNOM_OP next);
INT m_skn(VECTOR_OP self, SYM_OP koeff, POLYNOM_OP next);
INT m_skin(VECTOR_OP self, INT koeff, POLYNOM_OP next);
INT m_iindex(INT i);
INT m_iindex_iexponent(INT i, INT j);
INT nullp();
INT einsp();
INT m_scalar(SYM_OP a);
/* before: m_scalar_polynom(OP a, OP b) 
 * a ist scalar, b wird polynom. 
 * a [0]. */
INT mult_scalar(SYM_OP a, POLYNOM_OP erg);
/* before: mult_scalar_polynom(OP a, OP poly, OP erg) */
INT mult_polynom(POLYNOM_OP zwei, POLYNOM_OP c);
/* mult_polynom_polynom(OP eins, OP zwei, OP c) */
INT mult(SYM_OP b, POLYNOM_OP d);
INT mult_apply_scalar(SYM_OP a);
/* before: mult_apply_scalar_polynom(OP a, OP b) */
INT mult_apply_polynom(POLYNOM_OP b);
/* b = a * b */
/* before: mult_apply_polynom_polynom(OP a, OP b) */
INT mult_apply(SYM_OP b);
/* before: mult_apply_polynom(OP a, OP b) */
INT addinvers(POLYNOM_OP b);
/* before: addinvers_polynom(OP a, OP b) */
INT addinvers_apply();
/* before: addinvers_apply_polynom(OP a) */
INT add_polynom(POLYNOM_OP b, POLYNOM_OP c);
/* before: add_polynom_polynom(OP a, OP b, OP c) */
INT add_scalar(SYM_OP a, POLYNOM_OP c);
/* before: add_scalar_polynom(OP a, OP b, OP c) */
/* a scalar, b polynom, c result */
INT add(SYM_OP b, POLYNOM_OP c);
INT add_apply_polynom(POLYNOM_OP b);
INT add_apply_scalar(SYM_OP a);
/* before: add_apply_polynom_scalar(OP a, OP b) 
 * a polynom, b scalar. berechnet b := a + b. */
INT add_apply(SYM_OP b);
INT cycle_ind_Cn(INT n);
INT cycle_ind_Dn(INT n);
INT cycle_ind_Sn(INT n);
INT cycle_ind_An(INT n);
INT cycle_index_direct_product(POLYNOM_OP p, POLYNOM_OP q);
INT cycle_index_direct_sum(POLYNOM_OP p, POLYNOM_OP q);
INT cycle_index_Sk_Snmk(INT n, INT k);
INT cycle_index_arbitrary(VECTOR_OP gen);
INT cycle_index_labra(LABRA_OP L);
INT sum_of_coefficients(SYM_OP sum);
INT PolyaSubst(POLYNOM_OP erzFkt);
INT cycle_ind_Sm_times_Sn(INT m, INT n);
INT LaTeX(FILE *outLatex, char c, INT indiziere);
INT nabla(POLYNOM_OP res, LABRA_OP G, INT n, INT f_v);
};
INT number_of_double_cosets_A_Sn_B(INT n, POLYNOM_OP ca, POLYNOM_OP cb, SYM_OP res);

void mw_Latex_ZykInd_und_erzFkt(POLYNOM_OP ZykInd, POLYNOM_OP erzFkt, INT m, INT n);
void mw_Zykelindizes_und_Anzahlen(INT m, INT n);
#endif /* POLYTRUE */

#endif /* POLY_INCLUDED */


