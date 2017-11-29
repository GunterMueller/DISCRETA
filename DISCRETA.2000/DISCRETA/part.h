/* part.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#ifndef PART_INCLUDED
#define PART_INCLUDED

#ifdef PARTTRUE
/* VECTOR:
 * schwach monoton WACHSENDE (franz. Notation) 
 * Folge der (echt positiven) Teile 
 * EXPONENT: 
 * for (i = s_li() - 1, j = 0; i >= 0; i--, j++)
 *   s_ii(i) ist Anzahl der (j + 1) Teile
 */
class partition_ob : public SYM_OB {
public:
	OBJECTKIND s_k() { 
		return ob_self.ob_partition->pa_kind; };
	void c_k(OBJECTKIND kind) { 
		ob_self.ob_partition->pa_kind = kind; };

	VECTOR_OP s_s() { 
		return (VECTOR_OP)
			ob_self.ob_partition->pa_self; };
	void c_s(VECTOR_OP self) { 
		ob_self.ob_partition->pa_self = (OP)self; };

	INTEGER_OP s_l() { return s_s()->s_l(); };
	INT s_li() { return s_s()->s_li(); };
	INTEGER_OP s_i(INT i) { 
		return (INTEGER_OP) s_s()->s_i(i); };
	INT s_ii(INT i) { return s_i(i)->s_i(); };
	void c_i(INT i, INTEGER_OP p) { *s_i(i) = *p; };
	void c_ii(INT i, INT j) { s_i(i)->c_i(j); };
	void m_ii(INT i, INT j) { s_i(i)->m_i(j); };

INT freeself();
INT m_ks(OBJECTKIND kind, VECTOR_OP self);
INT b_ks(OBJECTKIND kind, VECTOR_OP self);
INT m_kli(OBJECTKIND kind, INT l);
INT m_kl(OBJECTKIND kind, INTEGER_OP l) { 
	return m_kli(kind, l->s_i()); };
INT b_i(INTEGER_OP l);
/* Bsp: 5 --> [5] */
INT partitionp();
/* VECTOR: Eintraege muessen groesser null sein 
 * und schwach monoton wachsen. */
INT weight_i();
INT weight_augmented_i();
INT length_i();
INT strictp();
/* true if no equal parts */
INT copy(PARTITION_OP b);
INT fprint(FILE *f);
INT sprint(BYTE *str);
INT sprint_latex(BYTE *s);
INT latex(FILE *fp);
INT compare(PARTITION_OP b);
INT conjugate(PARTITION_OP b);
/* before: fastconjugate_partition(OP part, OP b) */
INT equal_parts(INTEGER_OP b);
INT dimension(INTEGER_OP b);
/* before: dimension_partition(OP a, OP b) */
/* es wird die Dimension des durch die 
 * Partition bezeichneten Characters der Sn
 * berechnet. Dazu wird die Hakenformel verwendet */
INT dimension_augmented(INTEGER_OP b);
/* before: dimension_augpart() */
INT hook_length_i(INT i, INT j);
INT hook_length_augmented_i(INT i, INT j);
INT t_VECTOR_EXPONENT(PARTITION_OP nach);
INT t_EXPONENT_VECTOR(PARTITION_OP b);
INT induce2(PARTITION_OP nach);
INT augment();
INT de_augment();
INT first(INT n);
INT next_VECTOR_apply();
INT next_VECTOR(PARTITION_OP next);
INT first_into_k_parts_VECTOR(INT n, INT k);
INT next_into_k_parts_VECTOR(INT n, INT k);
INT class_length_Sn(SYM_OP res);
INT centralizer_order_Sn(SYM_OP res);
INT exp_nb_i_parts(INT i);
INT equal_parts(VECTOR_OP v);
};
INT part_ende();
INT charvalue(PARTITION_OP rep, 
	PARTITION_OP part, INTEGER_OP res);
INT charvalue_ij(VECTOR_OP V, 
	INT rep_i, INT part_i, INTEGER_OP res);
INT chartafel(INT n, MATRIX_OP M);
INT chartafel_vop(INT n, MATRIX_OP M, VECTOR_OP V);
INT vector_of_part(VECTOR_OP V, INT n);
INT kostka_number(PARTITION_OP inh, PARTITION_OP umriss, INTEGER_OP res);
#endif /* PARTTRUE */

#endif /* PART_INCLUDED */

