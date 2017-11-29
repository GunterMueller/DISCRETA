/* lb.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef LABRA_INCLUDED
#define LABRA_INCLUDED

#ifndef PERM_INCLUDED
#include <DISCRETA/perm.h>
#endif

class labra_ob : public VECTOR_OB {
public:
	INTEGER_OP s_degree() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_degree_i() { 
		return(s_degree()->s_i()); };

	INTEGER_OP s_nb_gen() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_nb_gen_i() { 
		return(s_nb_gen()->s_i()); };
	VECTOR_OP s_G() { 
		return((VECTOR_OP)s_i(2)); };
	PERMUTATION_OP s_G_i(INT i) { 
		return((PERMUTATION_OP) s_G()->s_i(i)); };

	VECTOR_OP s_V() { 
		return((VECTOR_OP)s_i(3)); };
	INTEGER_OP s_V_i(INT i) { 
		return((INTEGER_OP)s_V()->s_i(i)); };
	INT s_V_ii(INT i) { 
		return(s_V_i(i)->s_i()); };

	VECTOR_OP s_KM() { 
		return((VECTOR_OP)s_i(4)); };
	PERMUTATION_OP s_KM_i(INT i) { 
		return((PERMUTATION_OP) s_KM()->s_i(i)); };

	VECTOR_OP s_EM() { 
		return((VECTOR_OP)s_i(5)); };
	PERMUTATION_OP s_EM_i(INT i) { 
		return((PERMUTATION_OP) s_EM()->s_i(i)); };

	PERMUTATION_OP s_base() { 
		return((PERMUTATION_OP)s_i(6)); };
	INTEGER_OP s_base_i(INT i) {
		return((INTEGER_OP) s_base()->s_i(i)); };
	INT s_base_ii(INT i) {
		return(s_base_i(i)->s_i()); };
	PERMUTATION_OP s_inv_base() { 
		return((PERMUTATION_OP)s_i(7)); };
	INTEGER_OP s_inv_base_i(INT i) {
		return((INTEGER_OP) s_inv_base()->s_i(i)); };
	INT s_inv_base_ii(INT i) {
		return(s_inv_base_i(i)->s_i()); };

	VECTOR_OP s_T() { 
		return((VECTOR_OP)s_i(8)); };
	VECTOR_OP s_T_i(INT i) { 
		return((VECTOR_OP)s_T()->s_i(i)); };
	INTEGER_OP s_T_ij(INT i, INT j) { 
		return((INTEGER_OP) s_T_i(i)->s_i(j)); };
	INT s_T_iji(INT i, INT j) { 
		return(s_T_ij(i, j)->s_i()); };

	VECTOR_OP s_tidx() { 
		return((VECTOR_OP)s_i(9)); };
	INTEGER_OP s_tidx_i(INT i) { 
		return((INTEGER_OP) s_tidx()->s_i(i)); };
	INT s_tidx_ii(INT i) { 
		return(s_tidx_i(i)->s_i()); };
	
	INTEGER_OP s_f_verbose() { 
		return((INTEGER_OP)s_i(10)); };
	INT s_f_verbose_i() { 
		return(s_f_verbose()->s_i()); };
	INTEGER_OP s_f_very_verbose() { 
		return((INTEGER_OP)s_i(11)); };
	INT s_f_very_verbose_i() { 
		return(s_f_very_verbose()->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init_quick(VECTOR_OP gen);
INT Init_no_generators(INT degree, INT f_verbose, INT f_very_verbose);
INT Init(INT degree, VECTOR_OP theGenerators, INT nb_gen, 
	INT f_verbose, INT f_very_verbose);
INT build(VECTOR_OP gen);
INT sprint(BYTE *s);
INT check_consistency();
INT Print_Sims();
INT Print();
INT print_T();
INT latex_transversal_indices(FILE *fp);
INT calc_T();
INT calc_transversal_matrix(MATRIX_OP M, INT f_v);
INT is_cyclic(VECTOR_OP path, PERMUTATION_OP per);
INT transversal_tree(BYTE *group_label, INT f_full, INT f_labels);
INT elements_first(VECTOR_OP path, PERMUTATION_OP p);
INT elements_next(VECTOR_OP path, PERMUTATION_OP p);
INT element_from_path(PERMUTATION_OP p, INT *path, INT path_len);
INT element_from_short_path(PERMUTATION_OP p, VECTOR_OP path);
INT transversal_as_set(INT i, VECTOR_OP chi);
INT stab_order(INT i, SYM_OP ord);
INT stab_generators(INT i, VECTOR_OP gen, INT f_v);
INT generators(VECTOR_OP p);
INT generators_bottom_up(VECTOR_OP p, INT nb_fix);
INT reduced_generating_set(VECTOR_OP p, INT f_bottom_up, INT f_v);
INT strong_generating_set(VECTOR_OP p, INT f_v);
INT group_order(SYM_OP go);
INT print_group_order();
INT depth(INT i);
INT rep_ij(INT i, INT j, PERMUTATION_OP rep);
INT path_exists(INT j, INT k);
/* TRUE, iff path j -> k exists */
INT path_exists_el(INT j_el, INT k_el);
/* TRUE, iff path j -> k exists 
 * j := inv_base[j_el], 
 * k := inv_base[k_el], 
 */
INT child_vectors(VECTOR_OP cv);
INT new_edge(INT i_el, INT j_el, PERMUTATION_OP m);
INT update_EM();
INT orbit(INT omega, VECTOR_OP G, INT G_len, INT f_verbose);
INT in_test(PERMUTATION_OP g);
INT sift(PERMUTATION_OP g);
INT schreier(VECTOR_OP G, INT G_len, INT omega, INT f_verbose);
INT schreier_and_sift(VECTOR_OP G, INT G_len, 
	INT omega, INT f_verbose, 
	VECTOR_OP G1, INTEGER_OP G1_len);
INT get_generators(VECTOR_OP V, INTEGER_OP V_len, INT j);
INT jerrum(INT f_verbose);
INT Sym_n(INT deg, INT f_verbose);
INT init_by_generators(VECTOR_OP G, INT deg, INT f_verbose);
INT all_sons(INT i, VECTOR_OP sons);
// sons of i are immediate descendants of i
// (there exists an edge \pi_i -> \pi_j
// sons becomes an integer vectors of indices of base points: 
// { KM(\pi(j)) | j = sons[i], i < sons->s_li() }
// is the set of edges leading to sons.
INT the_son(INT j, PERMUTATION_OP p);
INT all_stab_generators(INT i, VECTOR_OP gens);
// the group fixing \pi_1, ... , \pi_i is returned.
// gens contains a vector of indices of base points 
// such that
// { KM(\pi_j) | j = gens[i], i < gens->s_li() }
// is a generating set.
INT cycle(INT r, INT s, INT f_v, INT f_vv);
INT cycle_id();
INT orbits(VECTOR_OP tda_ordering, 
	VECTOR_OP orbit_first, VECTOR_OP orbit_length, 
	INT f_calc_labras, VECTOR_OP O_labra, INT f_v);
INT one_orbit(VECTOR_OP chi);
INT set_canon(VECTOR_OP X, VECTOR_OP X0, INT f_v, INT f_vv);
INT subset_orbits(VECTOR_OP X, VECTOR_OP Orbits, 
	VECTOR_OP Orbit_lengths, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2, INT t, INT *p_k_idx);
INT group_table(MATRIX_OP table, VECTOR_OP table_inv, 
	VECTOR_OP elements, VECTOR_OP generators, INT f_v);
};
INT test_jerrum_Sym(INT n, INT f_verbose);
void my_srand(INT seed);
INT my_rand(void);
INT lb_test_generators_from_file(BYTE *fname);
INT orbits_cheap(VECTOR_OP gen, VECTOR_OP tda_ordering, 
	VECTOR_OP orbit_first, VECTOR_OP orbit_length, INT f_v);

/* single_cosets.C */

#define SINGLE_COSET_MAX_DEG 500

typedef struct single_coset_work SINGLE_COSET_WORK;

struct single_coset_work {
	INT trans_rep[SINGLE_COSET_MAX_DEG];
	INT trans_rep_last[SINGLE_COSET_MAX_DEG];
	LABRA_OP G, U;
	// MATRIX_OP MG, MU;
	INT MG[SINGLE_COSET_MAX_DEG][SINGLE_COSET_MAX_DEG];
	INT MU[SINGLE_COSET_MAX_DEG][SINGLE_COSET_MAX_DEG];
	INT deg;
};

INT single_coset_labra_file();
SINGLE_COSET_WORK *single_coset_open(LABRA_OP G, LABRA_OP U, INT f_v);
INT single_coset_free(SINGLE_COSET_WORK *scw);
INT single_coset_first(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep);
INT single_coset_next(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep);

#endif /* LABRA_INCLUDED */

