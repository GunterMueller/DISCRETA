/* solvable.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef SOLVABLE_INCLUDED
#define SOLVABLE_INCLUDED

#ifndef DIVS_INCLUDED
#include <DISCRETA/divs.h>
#endif
#ifndef PERM_INCLUDED
#include <DISCRETA/perm.h>
#endif
#ifndef CP_INCLUDED
#include <DISCRETA/cp.h>
#endif
#ifndef MA_INCLUDED
#include <DISCRETA/ma.h>
#endif
#ifndef LO_INCLUDED
#include <DISCRETA/lo.h>
#endif

#undef FG_ISO_USE_CANONIC

typedef struct autlog_info AUTLOG_INFO;

typedef struct group_table GROUP_TABLE;

typedef struct gt_canon_info GT_CANON_INFO;

typedef struct fg_io_channel FG_IO_CHANNEL;
typedef struct fg_io FG_IO;

struct fg_io_channel {
	INT f_ascii;
	BYTE *prefix;
	BYTE fname[1024];
	FILE *fp;
	DATABASE_OP db;
};

struct fg_io {
	FG_IO_CHANNEL from;
	FG_IO_CHANNEL to;
};

// fg.C:
INT fg_io_open_channel_r_ascii(FG_IO_CHANNEL *ch, BYTE *prefix);
INT fg_io_open_channel_w_ascii(FG_IO_CHANNEL *ch, BYTE *prefix);
INT fg_io_open_channel_rw_ascii(FG_IO_CHANNEL *ch, BYTE *prefix, BYTE *mode);
INT fg_io_open_channel_r_db(FG_IO_CHANNEL *ch, BYTE *prefix);
INT fg_io_open_channel_w_db(FG_IO_CHANNEL *ch, BYTE *prefix);
INT fg_io_open_channel_rw_db(FG_IO_CHANNEL *ch, BYTE *prefix);
INT fg_io_close_channel(FG_IO_CHANNEL *ch);


class ze_ob : public VECTOR_OB {
public:
	INTEGER_OP s_step() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_step_i() { 
		return(s_step()->s_i()); };
	/* former: i */
	INTEGER_OP s_p() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_p_i() { 
		return(s_p()->s_i()); };
	INTEGER_OP s_n0() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_n0_i() { 
		return(s_n0()->s_i()); };
	INTEGER_OP s_n() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_n_i() { 
		return(s_n()->s_i()); };
	INTEGER_OP s_P() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_P_i() { 
		return(s_P()->s_i()); };
	VECTOR_OP s_A() { 
		return((VECTOR_OP)s_i(5)); };
	INTEGER_OP s_A_i(INT i) { 
		return((INTEGER_OP)s_A()->s_i(i)); };
	INT s_A_ii(INT i) { 
		return(s_A_i(i)->s_i()); };
	VECTOR_OP s_Av() { 
		return((VECTOR_OP)s_i(6)); };
	INTEGER_OP s_Av_i(INT i) { 
		return((INTEGER_OP)s_Av()->s_i(i)); };
	INT s_Av_ii(INT i) { 
		return(s_Av_i(i)->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init(INT step, INT p, INT n, INT P);
INT sprint(BYTE *s);
INT Print(FG_OP fg);
};


/* please update fg_util.C if this class is changed !
 * (FG_IDX_N, FG_IDX_M, FG_IDX_O, FG_IDX_HASH)
 */

class fg_ob : public VECTOR_OB {
public:
	/*
	 * Integers:
	 * version
	 * n
	 * m
	 * o
	 * GmZ_n
	 * Gd_n
	 */
	INTEGER_OP s_version() { 
		return((INTEGER_OP) s_i(0)); };
	INT s_version_i() { 
		return(s_version()->s_i()); };
	/* group order: */
	INTEGER_OP s_n() { 
		return((INTEGER_OP) s_i(1)); };
	INT s_n_i() { 
		return(s_n()->s_i()); };
	/* m-th group of order n (m >= 1): */
	INTEGER_OP s_m() { 
		return((INTEGER_OP) s_i(2)); };
	INT s_m_i() { 
		return(s_m()->s_i()); };
	
	/* the family (o >= 1): */
	INTEGER_OP s_o() { 
		return((INTEGER_OP) s_i(3)); };
	INT s_o_i() { 
		return(s_o()->s_i()); };
	/* the simplest family invariants: */
	INTEGER_OP s_GmZ_n() { 
		return((INTEGER_OP) s_i(4)); };
	INT s_GmZ_n_i() { 
		return(s_GmZ_n()->s_i()); };
	INTEGER_OP s_Gd_n() { 
		return((INTEGER_OP) s_i(5)); };
	INT s_Gd_n_i() { 
		return(s_Gd_n()->s_i()); };
	
	/* 
	 * a label for the group 
	 */
	STRING_OP s_label() {
		return((STRING_OP) s_i(6)); };
	BYTE *s_label_s() {
		return(s_label()->s_str()); };

  /* 
   * the normal series: 
   */
	INTEGER_OP s_nb_ze() { 
		return((INTEGER_OP) s_i(7)); };
	INT s_nb_ze_i() { 
		return(s_nb_ze()->s_i()); };
	VECTOR_OP s_ze() { 
		return((VECTOR_OP) s_i(8)); };
	ZE_OP s_ze_i(INT i) { 
		return((ZE_OP) s_ze()->s_i(i)); };

	
	/* 
	 * the group table
	 * is always available: 
	 */
 	MATRIX_OP s_theG() {
		return((MATRIX_OP) s_i(9)); };
	INTEGER_OP s_theG_ij(INT i, INT j) {
		return((INTEGER_OP) s_theG()->s_ij(i, j)); };
	INT s_theG_iji(INT i, INT j) {
		return( s_theG_ij(i, j)->s_i()); };
	INTEGER_OP s_theG_inv(INT i) {
		return( s_theG_ij(i, s_n_i())); };
	INT s_theG_inv_i(INT i) {
		return( s_theG_inv(i)->s_i()); };
	INTEGER_OP s_dim_n1() { 
		return((INTEGER_OP) s_i(10)); };
	INT s_dim_n1_i() { 
		return(s_dim_n1()->s_i()); };
	/* dimension (row length) of theG (= n + 5) */
	
	
	/*
	 * hash key 
	 * computed using the canonical form */
	INTEGER_OP s_f_has_hash() { 
		return((INTEGER_OP) s_i(11)); };
	INT s_f_has_hash_i() {
		return(s_f_has_hash()->s_i()); };
	VECTOR_OP s_hash() {
		return((VECTOR_OP) s_i(12)); };
	
	/* permutation for the canonical form: */
	INTEGER_OP s_f_has_p0() { 
		return((INTEGER_OP) s_i(13)); };
	INT s_f_has_p0_i() {
		return(s_f_has_p0()->s_i()); };
	PERMUTATION_OP s_p0() {
		return((PERMUTATION_OP) s_i(14)); };
	
	/* 
	 * the automorphism group:
	 * ago
	 * AutM
	 * T
	 */
	INTEGER_OP s_f_has_aut() { 
		return((INTEGER_OP) s_i(15)); };
	INT s_f_has_aut_i() {
		return(s_f_has_aut()->s_i()); };
	SYM_OP s_ago() { 
		return(s_i(16)); };
	MATRIX_OP s_AutM() {
		return((MATRIX_OP) s_i(17)); };
	PERMUTATION_OP s_AutM_ij(INT i, INT j) {
		return((PERMUTATION_OP) 
			s_AutM()->s_ij(i, j)); };
	VECTOR_OP s_T() { 
		return((VECTOR_OP) s_i(18)); };
	VECTOR_OP s_T_i(INT i) { 
		return((VECTOR_OP) s_T()->s_i(i)); };
	INTEGER_OP s_T_ij(INT i, INT j) { 
		return((INTEGER_OP) s_T_i(i)->s_i(j)); };
	INT s_T_iji(INT i, INT j) { 
		return(s_T_ij(i, j)->s_i()); };
	
	/* 
	 * classes of the automorphism group 
	 */
	INTEGER_OP s_f_has_aut_classes() { 
		return((INTEGER_OP) s_i(19)); };
	INT s_f_has_aut_classes_i() {
		return(s_f_has_aut_classes()->s_i()); };
	CLASS_REP_OP s_aut_classes() {
		return((CLASS_REP_OP) s_i(20)); };

	/* 
	 * subgroup lattice
	 * the sgl is either unfold (a SGL_OB - object)
	 * or packed (as a MEM_OB object)
	 */
	INTEGER_OP s_f_has_sgl() { 
		return((INTEGER_OP) s_i(21)); };
	INT s_f_has_sgl_i() {
		return(s_f_has_sgl()->s_i()); };
	SYM_OP s_sgl() {
		return(s_i(22)); };
	
	/* 
	 * the sylow type
	 */
	INTEGER_OP s_f_has_sylow_type() { 
		return((INTEGER_OP) s_i(23)); };
	INT s_f_has_sylow_type_i() {
		return(s_f_has_sylow_type()->s_i()); };
	VECTOR_OP s_sylow_type() { 
		return((VECTOR_OP) s_i(24)); };
	INTEGER_OP s_sylow_type_i(INT i) { 
		return((INTEGER_OP) 
		s_sylow_type()->s_i(i)); };
	INT s_sylow_type_ii(INT i) { 
		return(s_sylow_type_i(i)->s_i()); };
	
	/* 
	 * space for more information: 
	 */
	CONTI_OP s_C() {
		return((CONTI_OP) s_i(25)); };

INT field_name(INT i, INT j, BYTE *str);
INT export_GAP(FILE *fp);
INT save_ascii2(FILE *fp, INT f_long);
INT save_ascii(FILE *fp);
INT save_ascii_short(FILE *fp);
INT read_ascii(FILE *fp);
INT read_ascii_tg(FILE *fp, INT n);
INT init(INT nb_ze, INT *primes);
INT update();
INT init_extension(FG_OP G0, INT p, INT P, 
	PERMUTATION_OP aut);
INT init_Zp(INT p);
INT sprint(BYTE *s);
INT sprint1(BYTE *s);
INT sprint2(BYTE *s);
INT sprint3(BYTE *s);
INT fprint_gen(FILE *fp);
INT fprint_pres(FILE *fp);
INT Print_gen(void);
INT Print(void);
INT print_AutM_T();
INT sprint_aut_on_base(PERMUTATION_OP p, BYTE *str);
INT sprint_ze(INT i, BYTE *s);
INT sprint_tex_ze(INT i, BYTE *s);
INT fprint_int(FILE *fp, INT i);
INT print_int(INT i);
INT sprint_int(INT i, BYTE *str);
INT sprint_GAP_int(INT i, BYTE *str);
INT sprint_tex_int(INT i, BYTE *str);
INT sprint_int_vec(VECTOR_OP V, BYTE *str);
INT sprint_tex_int_vec(VECTOR_OP V, BYTE *str);
INT do_Q8();
INT do_Q8_2();
INT calc_ago_longint(LONGINT_OP ago);
INT sprint_ago(BYTE *str);
INT is_nilpotent();
INT nb_gen();
INT g_i(INT i);
INT solvable_group_from_catalogue(INT n, INT m);
INT solvable_group_from_file(BYTE *fname, INT m);
INT solvable_group_of_order_n_first(INT n, FG_IO **fgio);
INT solvable_group_of_order_n_next(FG_IO **fgio);


/* fg_color.C: */
INT order_structure(VECTOR_OP ord);
/* INT is_central(INT i); */
INT center(VECTOR_OP p);
INT Classes(BYTE *str, INT f_v, INT f_vv);
INT Coloring_seldoms_prefered(VECTOR_OP col);
INT Coloring(VECTOR_OP col, INT f_v, INT f_vv);
INT prepare_class_info(VECTOR_OP V, INT f_v, INT f_vv);

/* fg_table.C: */
INT int_shrink(FG_OP G, VECTOR_OP embedding, INT i);
INT i2gen_idx(INT i);
INT int2nw(INT ii, SHORT *nw);
INT nw2int(SHORT *nw, INT *ii);
INT mult_op(INTEGER_OP a, INTEGER_OP b, INTEGER_OP c);
INT inv_op(INTEGER_OP a, INTEGER_OP b);
INT conj_op(INTEGER_OP a, INTEGER_OP b, INTEGER_OP c);
INT commutator_op(INTEGER_OP a, INTEGER_OP b, INTEGER_OP c);
INT power_op(INTEGER_OP a, INT exp, INTEGER_OP c);
INT one_op(INTEGER_OP a);
INT onep_op(INTEGER_OP a);
INT gt_mult(INT a, INT b);
INT gt_inv(INT a);
INT gt_conj(INT a, INT b);
/* b^-1 a b */
INT gt_commutator(INT a, INT b);
/* a^-1 b^-1 a b */
INT gt_power(INT a, INT exp);
INT gt_order(INT i, INT *ord);
INT gt_order_if_prime(INT i, INT *ord);
INT gt_order_if_prime_power(INT i, INT *ord, INT *prime, INT *k);
INT calc_aut_i(INT i, PERMUTATION_OP aut, INT f_use_table);
INT apply_aut(INT i, SHORT *base_im, INT f_use_table);
INT bi2aut(SHORT *base_im, PERMUTATION_OP aut, INT f_use_table);
INT inv(INT i);
INT mult(INT i, INT j);
INT conjugate(INT i, INT j);
/* j^-1 * i * j */
INT commutator(INT i, INT j);
/* i^-1 * j^-1 * i * j */
INT power(INT i, INT exp);
INT nw_inv(INT i, SHORT *a, SHORT *b);
INT nw_mult(INT i, SHORT *a, SHORT *b, SHORT *c);
INT nw_vorbei_v(SHORT *a, SHORT *b, INT i, INT vorne);
INT nw_vorbei(SHORT *a, SHORT *b, INT i);
INT nw_reduce(INT i, SHORT *a);
INT nw_power_ip(INT i, SHORT *a, INT exp);
/* ip = in place */
INT nw_power(INT i, SHORT *a, SHORT *res, INT exp);
INT left_regular_representation(INT i, PERMUTATION_OP p);
/* Multiplikationsdarstellung des Elementes 
 * i auf allen 
 * Gruppenelementen: 0 <= elemente < n.
 * Die Permutation hat Eintraege 1 <= eintraege <= n. */
INT right_regular_representation(INT i, PERMUTATION_OP p);
INT inn_generator(INT i, PERMUTATION_OP p, INT f_use_table);
INT theG(INT f_verbose);
INT recalc_table(INT f_v);
INT theG_Zp_extension(FG_OP G0, PERMUTATION_OP A, INT f_verbose);
INT mult_Zp_extension(FG_OP G0, 
	INT i, INT ii, INT j, INT jj, 
	PERMUTATION_OP A, PERMUTATION_OP Av);
INT mult_Zp_extension2(INT n0, INT p, INT P, 
	INT i, INT ii, INT j, INT jj, 
	PERMUTATION_OP A, PERMUTATION_OP Av);
/* berechnet i * z^ii * j * z^jj. */
INT is_associative(INT f_verbose);
INT Inn_generators(VECTOR_OP gen, INT f_use_table);




/* fg_iso.C: */
INT int2rep(INT ii, INT *coset_reps);
INT rep2int(INT *coset_reps, INT *ii);
INT rep2aut(INT *coset_reps, PERMUTATION_OP aut);
INT aut2rep(PERMUTATION_OP aut, INT *coset_reps);
INT int2aut(INT ii, PERMUTATION_OP aut);
INT aut2int(PERMUTATION_OP aut, INT *ii);
INT aut2bi(PERMUTATION_OP aut, INT *base_im);
INT rep2bi(INT *coset_reps, INT *base_im);
INT int2bi(INT ii, INT *base_im);
INT bi2rep(INT *base_im, INT *coset_rep);
INT bi2int(INT *base_im, INT *ii);
INT test_bi();
INT autlog_info_fill_in_colors(FILE *fp, AUTLOG_INFO *info, INT f_B);
INT autlog_info_fill_in_theG(AUTLOG_INFO *info, INT f_B);
INT GmZ_Gd_orders(FILE *fp_txt, INT f_verbose);
INT Isomorphic(FG_OP fg2, INT f_verbose);
INT Isomorphic2(FG_OP fg2, INT f_verbose);
INT Isoclinic(FG_OP fg2, 
	INT f_verbose, INT f_very_verbose);
INT find_isomorphic_group(
	VECTOR_OP V, INT f_verbose);
INT find_isoclinic_group(
	VECTOR_OP V, INT f_verbose);
INT Aut(INT f_verbose, INT f_very_verbose);
INT get_aut_generators(VECTOR_OP V);
INT aut_on_generators(INT f_v, INT f_vv);
INT canonicize(PERMUTATION_OP p0, PERMUTATION_OP p0v, INT f_v, INT f_vv);
INT canonicize_with_generators(VECTOR_OP gen, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, INT f_v);
INT canonicize_qad(FILE *fp_txt, INT f_v, INT f_vv);
INT canonicize_qad_with_generators(FILE *fp_txt, VECTOR_OP gen, INT f_v);
INT Isomorphic_by_canonic_form(FG_OP G2, INT f_verbose);

/* fg_syl.C: */
/* INT fgg(); */
INT sgl(INT f_v, INT f_vv, INT f_vvv);
INT sylow_type(INT f_verbose);
INT max_p_Element(INT p, 
	INT *idx, INT *the_k, INT *f_found);
INT p_DiminoExtend(VECTOR_OP p, 
	VECTOR_OP gen, INT g, 
	INT f_p, INT prime, INT *f_not_added);
INT DiminoExtend(VECTOR_OP p, 
	VECTOR_OP gen, INT g);
INT Dimino(VECTOR_OP p);
INT p_Normalizer(INT p, 
	VECTOR_OP H, VECTOR_OP H_gen, 
	VECTOR_OP G, VECTOR_OP N, 
	VECTOR_OP N_gen);
INT p_Sylow(INT p, VECTOR_OP G, 
	VECTOR_OP P, VECTOR_OP P_gen, 
	INT f_verbose);

/* fg_ext.C: */
INT print_extension_matrix(CLASS_REP_OP R, BYTE *ex, INT p);
INT fprint_extension_matrix(FILE *fp_txt, CLASS_REP_OP R, BYTE *ex, INT p);
INT calc_Extension_matrix(FILE *fp_txt, MATRIX_OP M, INT p, INT f_reduced, INT f_verbose);
INT calc_extension_matrix(FILE *fp_txt, BYTE **ex, INT p, INT f_verbose);
INT rpe_by_stab(FILE *fp_txt, BYTE *ex, INT p, INT f_verbose);
INT calc_aut_classes_using_file(FILE *fp_txt, INT f_v, INT f_vv);
INT init_extension_matrix(FILE *fp_txt, INT p, BYTE **ex, INT f_verbose);
INT calc_aut_sgl_syl(FILE *fp_txt, 
	INT f_calc_aut, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, 
	INT f_v, INT f_vv);
INT do_extension_Zp(INT p, 
	INT f_calc_aut, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family, 
	INT f_reduce_by_classes, 
	VECTOR_OP theGroups, VECTOR_OP theFamilies);
INT determine_family(VECTOR_OP theFamilies, 
	INT f_verbose, INT f_very_verbose);

/* fg_direct.C: */
INT gen_is_maximal(INT i);
INT gen_is_central(INT i);
INT gen_is_direct(INT i);
INT adf_decomposition(FILE *fp_txt, INT f_v, INT f_vv);
INT complete_residual_group(FILE *fp_txt, FG_OP G, 
	VECTOR_OP Ord, INT f_v);
INT residual_group(FG_OP G, VECTOR_OP embedding, INT f_v);
INT complete_residue(VECTOR_OP embedding, VECTOR_OP primes, VECTOR_OP Gidx, VECTOR_OP Ord);
INT gen_power_ord_and_index(INT i, VECTOR_OP Ord, VECTOR_OP Gidx);
INT is_adf_irreducible(FILE *fp_txt, INT f_v, INT f_vv);
INT abelian_direct_factors(VECTOR_OP exp, VECTOR_OP Ord, VECTOR_OP Gidx, 
	VECTOR_OP embedding, VECTOR_OP primes, INT f_v);

};

/* fg_table.C: */
INT NW_one(SHORT *a, INT i);
INT NW_is_one(SHORT *a, INT i);
INT NW_copy(SHORT *a, SHORT *b, INT i);
INT NW_2_str(SHORT *a, INT i, BYTE *str);
/* haengt an str an */
INT NW_2_str_GAP(SHORT *a, INT i, INT nb_gen, BYTE *str);
INT NW_2_str_tex(SHORT *a, INT i, BYTE *str);
INT NW_shrink(SHORT *a, INT i, VECTOR_OP embedding);
INT NW_pad(SHORT *a, INT j, INT i);
INT NW_print(SHORT *a, INT i);
INT print_vec_of_fg(VECTOR_OP V);
INT do_fg_test(void);

/* fg_ext.C: */
INT fg_cmd_file_for_order_n(INT n_from, INT n_to, 
	INT f_calc_aut_group, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family);
INT fg_pvm_cmd_file_for_order_n(INT n_from, INT n_to, INT f_pvm, 
	INT f_calc_aut_group, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family);
INT extend_Zp_vec_of_groups(VECTOR_OP V, INT p, 
	INT f_calc_aut, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family, 
	INT f_reduce_by_classes, 
	VECTOR_OP theGroups, VECTOR_OP theFamilies);

/* fg_color.C: */
INT group_init_from_fg(GROUP_TABLE *G, FG_OP fg);

/* 
 * cl_rep.C
 */

class class_rep_ob : public VECTOR_OB {
public:
	/* the number of class representatives 
	 * (= s_R()->s_li() */
	INTEGER_OP s_nb_classes() { 
		return((INTEGER_OP) s_i(0)); };
	INT s_nb_classes_i() { 
		return(s_nb_classes()->s_i()); };
	/* the base len, 
	 * this is the length of the integer vectors 
	 * describing an (automorphism-) group element 
	 * by its image "`on the base"'.
	 * used to represent elements in R and Stab.
	 */
	INTEGER_OP s_base_len() { 
		return((INTEGER_OP) s_i(1)); };
	INT s_base_len_i() { 
		return(s_base_len()->s_i()); };
	
	VECTOR_OP s_base() {
		return((VECTOR_OP) s_i(2)); };
	INTEGER_OP s_base_i(INT i) {
		return((INTEGER_OP) s_base()->s_i(i)); };
	INT s_base_ii(INT i) {
		return(s_base_i(i)->s_i()); };

	MATRIX_OP s_R() {
		return((MATRIX_OP) s_i(3)); };
	INTEGER_OP s_R_ij(INT i, INT j) {
		return((INTEGER_OP) s_R()->s_ij(i, j)); };
	INT s_R_iji(INT i, INT j) {
		return(s_R_ij(i, j)->s_i()); };
	
	/* generators for centralizer, 
	 * elements as base images. 
	 * Stab[i]       is centralizer of i-th class representative,
	 *               it is a vector of generators. 
	 *               the vector elements are the rows of a matrix.
	 * Stab[i][j]    j-th element of i-th centralizer, 
	 *               given as a vector of base images. */
	VECTOR_OP s_stab() {
		return((VECTOR_OP) s_i(4)); };
	MATRIX_OP s_stab_i(INT i) {
		return((MATRIX_OP) s_stab()->s_i(i)); };
	INTEGER_OP s_stab_ijk(INT i, INT j, INT k) {
		return((INTEGER_OP) s_stab_i(i)->s_ij(j, k)); };
	INT s_stab_ijki(INT i, INT j, INT k) {
		return(s_stab_ijk(i, j, k)->s_i()); };

	INTEGER_OP s_version() { 
		return((INTEGER_OP) s_i(5)); };
	INT s_version_i() { 
		return(s_version()->s_i()); };
		
	/* cl = class length
	 * eo = element order
	 */
	VECTOR_OP s_cl() {
		return((VECTOR_OP) s_i(6)); };
	SYM_OP s_cl_i(INT i) {
		return(s_cl()->s_i(i)); };
	VECTOR_OP s_eo() {
		return((VECTOR_OP) s_i(7)); };
	SYM_OP s_eo_i(INT i) {
		return(s_eo()->s_i(i)); };

INT field_name(INT i, INT j, BYTE *str);
INT init(FG_OP G);
INT sprint(BYTE *s);
INT Print(INT f_verbose, INT f_very_verbose);
INT recalc_cl_eo(FG_OP G);
INT calc_aut_classes_using_file(FILE *fp_txt, FG_OP G, INT f_v, INT f_vv);
INT trivial_aut_classes(FILE *fp_txt, FG_OP G, INT f_v, INT f_vv);
INT get_aut_class_info(FG_OP G, VECTOR_OP cl, VECTOR_OP eo, INT f_v);
INT get_aci_ordered(FG_OP G, VECTOR_OP val, VECTOR_OP mult, INT f_v, INT f_vv);
INT sprint_aut_class_structure(FG_OP G, BYTE *str, INT f_v, INT f_vv);
INT prepare_aut_class_info(FG_OP G, VECTOR_OP V);
INT calc_aut_classes(FILE *fp_txt, FG_OP G, INT f_v, INT f_vv);
};

#define GROUP_CANONIC_FORM_MAX_AGO 9999999

class group_canonic_form_ob : public VECTOR_OB {
public:
	INTEGER_OP s_nb_gen() { 
		return((INTEGER_OP) s_i(0)); };
	INT s_nb_gen_i() { 
		return(s_nb_gen()->s_i()); };
	PERMUTATION_OP s_p0() {
		return((PERMUTATION_OP) s_i(1)); };
	PERMUTATION_OP s_p0v() {
		return((PERMUTATION_OP) s_i(2)); };
	VECTOR_OP s_g() { 
		return((VECTOR_OP) s_i(3)); };
	INTEGER_OP s_g_i(INT i) { 
		return((INTEGER_OP) s_g()->s_i(i)); };
	INT s_g_ii(INT i) { 
		return(s_g_i(i)->s_i()); };
	VECTOR_OP s_go() { 
		return((VECTOR_OP) s_i(4)); };
	INTEGER_OP s_go_i(INT i) { 
		return((INTEGER_OP) s_go()->s_i(i)); };
	INT s_go_ii(INT i) { 
		return(s_go_i(i)->s_i()); };
	VECTOR_OP s_ro() { 
		return((VECTOR_OP) s_i(5)); };
	INTEGER_OP s_ro_i(INT i) { 
		return((INTEGER_OP) s_ro()->s_i(i)); };
	INT s_ro_ii(INT i) { 
		return(s_ro_i(i)->s_i()); };
	VECTOR_OP s_re() { 
		return((VECTOR_OP) s_i(6)); };
	INTEGER_OP s_re_i(INT i) { 
		return((INTEGER_OP) s_re()->s_i(i)); };
	INT s_re_ii(INT i) { 
		return(s_re_i(i)->s_i()); };
	MATRIX_OP s_GG() { 
		return((MATRIX_OP) s_i(7)); };
	INTEGER_OP s_GG_ij(INT i, INT j) { 
		return((INTEGER_OP) s_GG()->s_ij(i, j)); };
	INT s_GG_iji(INT i, INT j) { 
		return(s_GG_ij(i, j)->s_i()); };
	INTEGER_OP s_ago_trunc() { 
		return((INTEGER_OP) s_i(8)); };
	INT s_ago_trunc_i() { 
		return(s_ago_trunc()->s_i()); };
	INTEGER_OP s_f_has_auts() { 
		return((INTEGER_OP) s_i(9)); };
	INT s_f_has_auts_i() { 
		return(s_f_has_auts()->s_i()); };
	INTEGER_OP s_f_auts_on_canonic_form() { 
		return((INTEGER_OP) s_i(10)); };
	INT s_f_auts_on_canonic_form_i() { 
		return(s_f_auts_on_canonic_form()->s_i()); };
	VECTOR_OP s_auts() { 
		return((VECTOR_OP) s_i(11)); };
	PERMUTATION_OP s_auts_i(INT i) { 
		return((PERMUTATION_OP) s_auts()->s_i(i)); };
	VECTOR_OP s_first_moved() { 
		return((VECTOR_OP) s_i(12)); };
	INTEGER_OP s_first_moved_i(INT i) { 
		return((INTEGER_OP) s_first_moved()->s_i(i)); };
	INT s_first_moved_ii(INT i) { 
		return(s_first_moved_i(i)->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT init();
INT sprint(BYTE *s);
INT calc_hash(FILE *fp_txt, VECTOR_OP hash, INT n, INT f_v, INT f_vv);

};
INT cf_read_and_print_hash(VECTOR_OP hash);

class sgo_info_ob : public VECTOR_OB {
public:
	/* o_len = orbit length under Aut(G) */
	INTEGER_OP s_o_len() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_o_len_i() { 
		return(s_o_len()->s_i()); };
	
	/* so_len = orbit length under Inn(G) */
	INTEGER_OP s_so_len() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_so_len_i() { 
		return(s_so_len()->s_i()); };

	/* the subgroup itself: 
	 * n = order, m = isomorphism type (in the database) */
	INTEGER_OP s_n() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_n_i() { 
		return(s_n()->s_i()); };
	INTEGER_OP s_m() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_m_i() { 
		return(s_m()->s_i()); };
	
	/* vector of generators */
	VECTOR_OP s_generators() { 
		return((VECTOR_OP)s_i(4)); };
	INTEGER_OP s_generators_i(INT i) { 
		return((INTEGER_OP)s_generators()->s_i(i)); };
	INT s_generators_ii(INT i) { 
		return(s_generators_i(i)->s_i()); };
	
	/* normalizer (layer/orbit) */
	INTEGER_OP s_normalizer() { 
		return((INTEGER_OP)s_i(5)); };
	INT s_normalizer_i() { 
		return(s_normalizer()->s_i()); };
	
	/* information about maximal subgroups: 
	 * (orbits / alpha_sup) */
	VECTOR_OP s_max_sub_orbits() { 
		return((VECTOR_OP)s_i(6)); };
	INTEGER_OP s_max_sub_orbits_i(INT i) { 
		return((INTEGER_OP)s_max_sub_orbits()->s_i(i)); };
	INT s_max_sub_orbits_ii(INT i) { 
		return(s_max_sub_orbits_i(i)->s_i()); };

	VECTOR_OP s_max_sub_alpha() { 
		return((VECTOR_OP)s_i(7)); };
	INTEGER_OP s_max_sub_alpha_i(INT i) { 
		return((INTEGER_OP)s_max_sub_alpha()->s_i(i)); };
	INT s_max_sub_alpha_ii(INT i) { 
		return(s_max_sub_alpha_i(i)->s_i()); };
	
	/* layer / orbit */
	INTEGER_OP s_l() { 
		return((INTEGER_OP)s_i(8)); };
	INT s_l_i() { 
		return(s_l()->s_i()); };
	INTEGER_OP s_o() { 
		return((INTEGER_OP)s_i(9)); };
	INT s_o_i() { 
		return(s_o()->s_i()); };
	
	/* placement: (x,y) coordinates, o_dx */
	INTEGER_OP s_place_x() { 
		return((INTEGER_OP)s_i(10)); };
	INT s_place_x_i() { 
		return(s_place_x()->s_i()); };
	INTEGER_OP s_place_y() { 
		return((INTEGER_OP)s_i(11)); };
	INT s_place_y_i() { 
		return(s_place_y()->s_i()); };
	INTEGER_OP s_o_dx() { 
		return((INTEGER_OP)s_i(12)); };
	INT s_o_dx_i() { 
		return(s_o_dx()->s_i()); };
};

// aut.h

// #ifndef AUT_INCLUDED
// #include <DISCRETA/SOLVABLE/aut.h>
// #endif

#define AUTLOG_MAX_G 64

#define AUT_MODE_BREAK_AFTER_FST 0
#define AUT_MODE_ONLY_COSET_REPS 1
#define AUT_MODE_FULL 2

struct autlog_info {
	INT An;
	INT Bn;
	INT *AtheG; /* n x n; dimension n x Adim_n */
	INT *BtheG;
	INT Adim_n;
	INT Bdim_n;
	INT f_A_has_colors;
	VECTOR_OP A_colors;
	INT f_B_has_colors;
	VECTOR_OP B_colors;
	INT nb_isomorphisms;
	INT mode;
	INT f_autologisms;
	INT f_verbose;
	INT f_very_verbose;
	INT A_nb_gen;
	INT B_nb_gen;
	INT A_g[AUTLOG_MAX_G];
	INT B_g[AUTLOG_MAX_G];
		/* 0 .. nb_gen-1: the generators. */
	INT (*add_aut)(AUTLOG_INFO *info, INT nb_gen, 
		INT *gen, SPERM *cp);
	INT (*add_aut_perm)(AUTLOG_INFO *info, INT nb_gen, 
		INT *gen, PERMUTATION_OP p);
	void *data;
	FILE *fp_txt;
};

INT canonicize_group_qad(FILE *fp_txt, FG_OP G, INT f_v, INT f_vv);
/* qad stands for "`quick and dirty"' ! */
INT canonicize_group_cf(FILE *fp_txt, FG_OP G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv);
INT canonicize_group(FILE *fp_txt, FG_OP G, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut, 
	INT f_v, INT f_vv);
INT autlog_test(AUTLOG_INFO *info);
INT autlog_GmZ_Ad_orders(AUTLOG_INFO *info, 
	INT *GmZ_order, INT *Gd_order, FILE *fp_txt);
INT al_test();


// group-table related functions:

// #ifndef GT_INCLUDED
// #include <DISCRETA/SOLVABLE/gt.h>
// #endif

struct group_table {
	MATRIX_OP table;
	VECTOR_OP table_inv;
	INT n;
	VECTOR_OP generators;
		/* needed for group_is_central() */
};


/* gt_color.C: */
INT group_init_simple(GROUP_TABLE *G, 
	MATRIX_OP M, VECTOR_OP generators);
INT group_init_generators(GROUP_TABLE *G, VECTOR_OP gen);
INT group_free(GROUP_TABLE *G);
INT group_mult(GROUP_TABLE *G, INT a, INT b);
INT group_inv(GROUP_TABLE *G, INT a);
INT group_conj(GROUP_TABLE *G, INT a, INT b);
INT group_commutator(GROUP_TABLE *G, INT a, INT b);
INT group_power(GROUP_TABLE *G, INT a, INT exp);
INT group_order(GROUP_TABLE *G, INT i, INT *ord);
INT group_order_if_prime(GROUP_TABLE *G, INT i, INT *ord);
INT group_order_if_prime_power(GROUP_TABLE *G, INT i, 
	INT *ord, INT *prime, INT *k);
INT group_order_structure(GROUP_TABLE *G, VECTOR_OP ord);
INT group_is_central(GROUP_TABLE *G, INT i);
INT group_center(GROUP_TABLE *G, VECTOR_OP p);
INT group_classes(GROUP_TABLE *G, VECTOR_OP class_no, 
	VECTOR_OP sv_last, VECTOR_OP sv_gen, 
	VECTOR_OP class_rep, VECTOR_OP class_len, VECTOR_OP class_el_ord);
INT group_print_classes(GROUP_TABLE *G, 
	VECTOR_OP class_no, VECTOR_OP class_rep, VECTOR_OP class_len);
INT group_first_class_coloring(GROUP_TABLE *G, 
	VECTOR_OP class_no, VECTOR_OP class_rep, 
	VECTOR_OP class_len, VECTOR_OP class_el_ord, 
	VECTOR_OP col, MATRIX_OP Col, INT f_v);
INT group_Classes(GROUP_TABLE *G, BYTE *str, INT f_v, INT f_vv);
INT group_Coloring_seldoms_prefered(GROUP_TABLE *G, VECTOR_OP col, INT f_v);
INT group_Coloring(GROUP_TABLE *G, VECTOR_OP col, INT f_v, INT f_vv);
INT group_refine_by_roots(GROUP_TABLE *G, 
	VECTOR_OP col, MATRIX_OP Col, INT f_v);
INT group_refine_by_p_roots(GROUP_TABLE *G, 
	INT p, VECTOR_OP col, MATRIX_OP Col, INT f_v);
INT group_p_roots(GROUP_TABLE *G, INT p, VECTOR_OP nb_p_roots);
INT group_p_roots_in_color(GROUP_TABLE *G, INT p, INT color, 
	VECTOR_OP col, VECTOR_OP nb_p_roots);
INT group_prepare_class_info(GROUP_TABLE *G, 
	VECTOR_OP V, INT f_v, INT f_vv);

/* gt_col_util.C: */
INT color_join(VECTOR_OP col1, VECTOR_OP col2, 
	VECTOR_OP col3, MATRIX_OP Col, INT f_v, INT f_vv);
INT color_v_join(VECTOR_OP col_0, VECTOR_OP col_v, 
	VECTOR_OP new_col, MATRIX_OP new_Col, 
	INT f_v, INT f_vv);
INT number_of_colors(MATRIX_OP Col, INT f_v);
INT prefer_seldom_colors(VECTOR_OP col, VECTOR_OP new_col, 
	PERMUTATION_OP col_perm_new_old, PERMUTATION_OP col_perm_old_new);
INT get_ci_ordered(VECTOR_OP cl, VECTOR_OP eo, 
	VECTOR_OP val, VECTOR_OP mult, INT f_v, INT f_vv);
INT sprint_class_structure(VECTOR_OP class_len, 
	VECTOR_OP class_el_ord, BYTE *str, INT f_v);

/* gt_canon.C: */

#define GT_CANON_MAX_GEN 64

struct gt_canon_info {
	INT n;
	INT *theG; /* n x n; dimension n x dim_n */
	INT dim_n;
	INT f_has_colors;
	VECTOR_OP colors;
	INT nb_isomorphisms;
	/* INT mode; */
	/* INT f_autologisms; */
	INT f_verbose;
	INT f_very_verbose;
	FILE *fp_txt;
};

INT group_calc_hash(FILE *fp_txt, GROUP_TABLE *G, 
	VECTOR_OP hash, PERMUTATION_OP p0, INT f_v, INT f_vv);
#if 0
INT group_calc_aut(FILE *fp_txt, GROUP_TABLE *G, VECTOR_OP base, 
	MATRIX_OP M, VECTOR_OP T, SYM_OP ago, INT f_v, INT f_vv);
#endif
INT group_calc_cf(FILE *fp_txt, GROUP_TABLE *G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv);
INT group_calc_cf_tg(FILE *fp_txt, GROUP_TABLE *G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv);
INT group_calc_cf_ab(FILE *fp_txt, GROUP_TABLE *G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv);
INT group_canonicize(FILE *fp_txt, GROUP_TABLE *G, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut, 
	INT f_v, INT f_vv);




#endif /* SOLVABLE_INCLUDED */

