/* perm.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef PERM_INCLUDED
#define PERM_INCLUDED

#define NONCOMPARABLE 2

extern INT f_perm_print_start_with_zero;

class permutation_ob : public SYM_OB {
public:
	INT freeself();

	INT m_l(INTEGER_OP len);
	INT m_il(INT l);
	INT m_ks(OBJECTKIND kind, VECTOR_OP self);
		/* vector self wird kopiert. */
	INT b_ks(OBJECTKIND kind, VECTOR_OP self);
	VECTOR_OP s_s() { 
		return((VECTOR_OP)
			ob_self.ob_permutation->p_self); };
		/* vormals: s_p_s(p).
		 * Gibt das Bildvektor-Array zurueck. */
	void c_s(VECTOR_OP s) { 
		ob_self.ob_permutation->p_self = (OP)s; };
		/* vormals c_p_s(OP a, OP s). 
		 * Setzt p_self auf s. */
	OBJECTKIND s_k() { 
		return(ob_self.ob_permutation->p_kind); };
		/* vormals: s_p_k(p). 
		 * Gibt p_kind zurueck. */
	void c_k(OBJECTKIND kind) { 
		ob_self.ob_permutation->p_kind = 
			kind; };
	INTEGER_OP s_l() { 
		return(s_s()->s_l()); };
		/* vormals: s_p_l(). 
		 * Gibt die Laenge als Objekt zurueck. */
	INT s_li() { 
		return(s_l()->s_i()); };
		/* vormals: s_p_li(p)
		 * Gibt die Laenge des Bildvektor Arrays als INT 
		 * zurueck (d.h. den Grad der Permutation). */
	SYM_OP s_i(INT i) { 
		return(s_s()->s_i(i)); };
		/* vormals s_p_i(p, i)
		 * Gibt i-tes Objekt im 
		 * Bildvektor-Array zurueck. */
	INT s_ii(INT i) { 
		return(s_s()->s_ii(i)); };
		/* Vormals: s_p_ii(p, i)
		 * Gibt den i-ten Eintrag im 
		 * Bildvektor-Array 
		 * als INT zurueck, d.h. Bild(i+1). */
	void m_ii(INT i, INT j) { 
		s_s()->m_ii(i, j); };
		/* Setzt den i-ten Eintrag im 
		 * Bildvektor-Array auf j:
		 * d.h. Bild(i+1) := j. */
	void c_ii(INT i, INT j) { 
		s_s()->c_ii(i, j); };
		/* Setzt den i-ten Eintrag im 
		 * Bildvektor-Array auf j:
		 * d.h. Bild(i+1) := j. */

	// m_il_p, m_ks_p, b_ks_p, 

/* in iof.C: */
INT calc_size_on_file();
INT write_mem(MEM_OP mem, INT debug_depth);
INT read_mem(MEM_OP mem, INT debug_depth);
INT one();
INT fprint_list(FILE *fp);
INT print_list();
INT sprint_list(BYTE *str);
INT sprint(BYTE *str);
INT fprint_GAP(FILE *fp);
INT sprint_GAP(BYTE *str);
INT latex(FILE *fp);
INT sprint_latex(BYTE *s);
INT Add2Cycle(INT i0, INT i1);
INT Add3Cycle(INT i0, INT i1, INT i2);
INT Add4Cycle(INT i0, INT i1, INT i2, 
	INT i3);
INT Add5Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4);
INT Add6Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4, INT i5);
INT Add7Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4, INT i5, INT i6);
INT Add8Cycle(INT i0, INT i1, INT i2, 
	INT i3, INT i4, INT i5, 
	INT i6, INT i7);
INT AddNCycle(INT first, INT len);
INT AddNCycleOffset(INT first, INT len, INT offset);
INT is_full_cycle();
INT order(INT *order);
/* Rueckgabe: order == 1 falls p == id. 
 * Sonst order = Elementordnung(p). */
INT order_if_prime(INT *order);
/* Rueckgabe: order == 0, wenn 
 * Elementordnung (p) nicht prim und nicht eins;
 * sonst in order die entsprechende Primzahl 
 * bzw eins (wenn p == id). */
INT order_if_prime_power(INT *order, INT *prime, INT *k);
/* Rueckgabe: order == 0, 
 * wenn Elementordnung keine Primpotenz;
 * sonst in order die entsprechende 
 * Primpotenz (order == prime ^ k). 
 * order == 1 falls p == id 
 * (und prime == 1 und k == 0). */
INT IsEven() { return(signum() == 1); };

INT compare(PERMUTATION_OP b);
INT invers(PERMUTATION_OP b);
INT mult(PERMUTATION_OP b, PERMUTATION_OP c);
/* before: mult_permutation(OP a, OP b, OP c) */
/* c := a(b(i)) - (erst b, dann a). */
INT copy(PERMUTATION_OP b);
INT einsp();
INT first_permutation(INTEGER_OP l);
INT next_permutation_lex(
	PERMUTATION_OP next_perm);
INT perm_lehmercode(VECTOR_OP vec);
/* before: lehmercode_permutation(
 * OP perm, OP vec) */
/* diese prozedur berechnet zur 
 * permutation perm = [p1,....,pn]
 den zugehoerigen lehmercode vec [v1,...,vn] */
INT perm_lehmercode2(VECTOR_OP vec);
/* before: lehmercode2_permutation(
 * OP perm, OP vec) */
INT next_permutation(PERMUTATION_OP next);
/* before: next_permutation(OP start, OP n) */
/* Erzeuge naechste Permutation via Lehmercode nach n. 
 * Rueckgabe LASTPERMUTATION, wenn Ende. n ist dann 
 * unveraendert. */
INT signum();
/* before: signum_permutation(OP perm, OP b) */
INT numberof_inversions();
/* before: numberof_inversionen(OP a, OP b) */
/* Rueckgabe ist die anzahl der 
 * Inversionen in der permutation a */
INT rz(VECTOR_OP c);
/* reduzierte zerlegung */
INT sscan(BYTE **str, INT f_v);
INT cycle_type(VECTOR_OP v);
INT embed(INT n);
INT conjugate(PERMUTATION_OP kappa, PERMUTATION_OP p_kappa);

// perm2.C:
INT m_n_cycle(INT n);
INT reordering_for_cyclic_permutation(PERMUTATION_OP q);
INT reorder_to_cycles(PERMUTATION_OP q, VECTOR_OP v, INT f_v);
INT cyclic_multiplicator(INT n, INT a);
INT number_of_fixpoints();
#ifdef PARTTRUE
INT class_rep(PARTITION_OP type);
#endif
INT induce_action_on_columns(PERMUTATION_OP gg, MATRIX_OP I);
INT induce_action_on_blocks(PERMUTATION_OP gg, VECTOR_OP B);
INT induce3(PERMUTATION_OP b);
INT induce2(PERMUTATION_OP b, INT n);
INT induce_on_2tuples(PERMUTATION_OP p, INT f_injective);
INT induce2_01(PERMUTATION_OP b, INT n);
INT embed_at(PERMUTATION_OP b, INT n, INT at);
INT join(PERMUTATION_OP b, PERMUTATION_OP c);
INT add_fixpoint_in_front(PERMUTATION_OP b);
INT add_n_fixpoints_in_front(PERMUTATION_OP b, INT n);
INT add_n_fixpoints_at_end(PERMUTATION_OP b, INT n);
INT remove_fixpoint(PERMUTATION_OP b, INT i);
};
INT test_perm();

INT support_on_2_sets(VECTOR_OP R, INT n, VECTOR_OP support, INT *supp);
INT tuple2_rank(INT rank, INT *i, INT *j, INT n, INT f_injective);
INT tuple2_unrank(INT i, INT j, INT *rank, INT n, INT f_injective);
INT bruhat_comp_perm(PERMUTATION_OP a, PERMUTATION_OP b);
/* 1 if a>b 0 if a=b -1 if a<b NONCOMPARABLE else */
INT bru_comp(PERMUTATION_OP a, PERMUTATION_OP c);


// perm_grp.C:
INT vec_generators_is_trivial_group(VECTOR_OP gen);
INT is_abelian(VECTOR_OP G);
INT read_file_of_generators(VECTOR_OP G, char *fname);
INT write_file_of_generators(VECTOR_OP G, char *fname);
INT read_generators(VECTOR_OP G, FILE *fp);
INT write_generators(VECTOR_OP G, FILE *fp);
INT vec_induced_group_on_subset(VECTOR_OP V, VECTOR_OP subset, VECTOR_OP W);
INT vec_subgroup_of_hol_of_cyclic_group(VECTOR_OP V, INT n, INT i);
// n must be a prime !
INT vec_hol_of_cyclic_group(VECTOR_OP V, INT n);
INT vec_conjugate(VECTOR_OP gen, PERMUTATION_OP p);
INT vec_induce_action_on_blocks(VECTOR_OP gen, VECTOR_OP B);
INT vec_induce_action_on_columns(VECTOR_OP gen, MATRIX_OP M);
INT vec_induce3(VECTOR_OP gen);
INT vec_induce2(VECTOR_OP gen);
INT vec_induce_on_2tuples(VECTOR_OP gen, INT f_injective);
INT vec_add_fixpoint_in_front(VECTOR_OP gen);
INT vec_add_fixpoint_at_end(VECTOR_OP gen);
INT vec_generators_stabilize_point(VECTOR_OP a, VECTOR_OP b);
INT vec_generators_stabilize_a_point(VECTOR_OP a, INT pt, VECTOR_OP b);
// 0 <= pt < deg
INT vec_generators_degree(VECTOR_OP a);
INT vec_generators_group_order(VECTOR_OP a, SYM_OP go);
INT vec_generators_remove_fixpoint(VECTOR_OP a, VECTOR_OP b, INT i);
INT young2_generators(VECTOR_OP G, INT deg, INT t);
INT trivial_generators(VECTOR_OP G, INT deg);
INT cyclic_generators(VECTOR_OP G, INT deg);
INT symmetric_generators(VECTOR_OP G, INT deg);
INT symmetric_generators_pp(VECTOR_OP G, INT deg);
INT symmetric_generators_transpositions(VECTOR_OP G, INT deg);
INT alternating_generators(VECTOR_OP G, INT deg);
INT dihedral_generators(VECTOR_OP G, INT deg);
INT dihedral_generators_pp(VECTOR_OP G, INT deg);
INT Higman_Sims_176_generators(VECTOR_OP G);
INT number_of_solvable_groups(INT n);
INT solvable_group_permutation_representations(VECTOR_OP V, INT maxindex, INT ordering);
INT vector_select_ith(VECTOR_OP V, INT i);
INT solvable_group_from_catalog(VECTOR_OP V, INT n, INT m);
INT Mathieu_generators(VECTOR_OP G, INT n);
INT M11_generators(VECTOR_OP G);
/* 4 ply transitive of order 11 10 9 8, 
 * stabilizer of 12 in M12 */
INT M12_generators(VECTOR_OP G);
INT M23_generators(VECTOR_OP G);
/* 4 ply transitive, 
 * order 23 22 21 20 16 3, stabilizer of 24 in M24 */
INT M24_generators(VECTOR_OP G);
INT cube_generators(VECTOR_OP G);
INT Zn_mult_generators(INT n, INT mu, VECTOR_OP G);
INT Sn_wreath_Sm_generators(INT n, INT m, VECTOR_OP G);

// perm_grp2.C:
// binary operations on permutation groups given by generators !
INT vec_generators_diagonal_sum(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c);
INT vec_generators_comma(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c);
INT vec_generators_direct_sum(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c);
INT vec_generators_direct_product(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c);
INT perm_on_cartesian_product(PERMUTATION_OP a, PERMUTATION_OP b, PERMUTATION_OP c);
// a on the rows, b on the columns of the cartesian product
INT vec_generators_exponentiation(VECTOR_OP a, VECTOR_OP b, 
	VECTOR_OP c, INT f_simultaneously);
INT perm_on_exponentiation(PERMUTATION_OP a, INT a_in_component, 
	PERMUTATION_OP b, PERMUTATION_OP c, INT f_simultaneously);
// a on the rows, b on the columns of the matrix of the mapping
INT vec_generators_wreath_product(VECTOR_OP W, VECTOR_OP G, VECTOR_OP H, INT f_v);
INT wreath_embedding(PERMUTATION_OP g, INT n, INT m, PERMUTATION_OP q);
INT wreath_embedding_component(PERMUTATION_OP g, INT n, INT m, INT j, PERMUTATION_OP q);



#endif /* PERM_INCLUDED */

