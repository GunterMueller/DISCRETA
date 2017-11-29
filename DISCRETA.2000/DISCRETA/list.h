/* list.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#ifndef LIST_INCLUDED
#define LIST_INCLUDED

#ifdef LISTTRUE
class list_ob : public SYM_OB {
public:	
	SYM_OP s_s() { 
		return((SYM_OP)ob_self.ob_list->l_self); };
	MONOM_OP s_s_mo() { 
		return((MONOM_OP)ob_self.ob_list->l_self); };
	void c_s(SYM_OP self) { 
		ob_self.ob_list->l_self = (OP)self; };
	LIST_OP s_n() { 
		return((LIST_OP)ob_self.ob_list->l_next); };
	void c_n(LIST_OP next) { 
		ob_self.ob_list->l_next = (OP)next; };

INT m_sn(SYM_OP self, LIST_OP next);
/* before: m_sn_l(OP self, OP nx, OP a)*/
INT b_sn(SYM_OP self, LIST_OP next);
/* build_self next_list */
/* before: b_sn_l(OP self, OP nx, OP a) */
INT emptyp();
INT lastp();
INT length();
INT freeself();
INT fprint(FILE *f);
INT sprint(BYTE *str);
/* appends to str. writes to maximal strlength of 200. */
INT compare(LIST_OP b);
/* before: comp_list(OP a,OP b); */
/* vergleich zweier listen, 
 * z.b. 1,1,3  < 1,2,2 z.b. 2,2,3  > 2/3  */
INT transform_apply(INT (*tf)(SYM_OP self));
/* before: transform_apply_list(
 * OP von, INT (*tf)(OP zeiger)) */
INT transform(LIST_OP nach, 
	INT (*tf)(SYM_OP self_von, SYM_OP self_nach));
/* before: transformlist(OP von, OP nach, 
 *   INT (*tf)(OP zeiger, OP nachzeiger)) */
/* nach wird vor dem Ueberschreiben freigegeben. */
INT copy(LIST_OP nach);
INT transform_by(SYM_OP ve, LIST_OP to,
	INT (*tf)(SYM_OP ve, SYM_OP self_from, SYM_OP self_to));
/* before: trans2formlist(OP ve, OP vz, OP nach,
 *    INT (*tf)(OP ve, OP zeiger, OP nachzeiger)) */
/* ve ist konstante , vz ist liste */
INT insert(SYM_OP insert_this, 
	INT (*equality_handler)(
		SYM_OP self_from, SYM_OP self_to), 
	INT (*compare_func)(SYM_OP from, SYM_OP to));
/* before: insert_list(OP von, OP nach, 
 * INT (*eh)(OP von, OP nn), INT (*cf)(OP von, OP nn)) */
/* fuegt das object von in die liste nach ein */
/* von ist keine liste */
INT insert_list(LIST_OP insert_this_list, 
	INT (*equality_handler)(
		SYM_OP self_from, SYM_OP self_to), 
	INT (*compare_func)(SYM_OP from, SYM_OP to));
/* before: insert_list_list_2(OP von, OP nach, 
	INT (*eh)(OP von, OP nn), 
	INT (*cf)(OP von, OP nn)) */
/* Das object insert_this_list wird verbraucht */
/* compare_func NULL means compare_ab() */
};
INT test_list();
INT list_cmp_func(SYM_OP a, SYM_OP b);
#endif /* LISTTRUE */

#endif /* LIST_INCLUDED */
