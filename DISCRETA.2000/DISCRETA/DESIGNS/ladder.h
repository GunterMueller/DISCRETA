/* ladder.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef LADDER_INCLUDED
#define LADDER_INCLUDED

#ifndef VDI_INCLUDED
#include <DISCRETA/graphics.h>
#endif
#ifndef MA_INCLUDED
#include <DISCRETA/ma.h>
#endif
#ifndef FGA_INCLUDED
#include <DISCRETA/fga.h>
#endif
#ifndef LABRA_INCLUDED
#include <DISCRETA/lb.h>
#endif
#ifndef DIVS_INCLUDED
#include <DISCRETA/divs.h> /* for STRING ->s_str() */
#endif
#ifndef UNIP_INCLUDED
#include <DISCRETA/unip.h>
#endif
#ifndef DB_INCLUDED
#include <DISCRETA/db.h>
#endif

typedef struct dc_lattice DC_LATTICE;
typedef struct orderly_sets ORDERLY_SETS;


/* ladder_info.C: */

typedef struct ladder_info LADDER_INFO;

struct ladder_info {
	VECTOR_OP generators;
	SYM_OP go;
	INT deg;
	BYTE *g_label, *g_label_tex;
	
	INT t, k;
	INT lambda;
	INT f_verbose;

	INT i, per_kind;
	INT t0, t1, user_time;
	BYTE str[256];
	BYTE txt_out[256];
	BYTE bin_out[256];
	
	BYTE generators_fname[1024];
	FILE *fp;
	INT up_to_step; /* := 2 * k - 1 */

	INT type;
	void *data;

	VECTOR_OP dc;
	SYM_OP id;
	
};

/* ladder.C: */

class dcy_ob : public VECTOR_OB {
public:
	INTEGER_OP s_step() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_step_i() { 
		return(s_step()->s_i()); };
	/* type = 0: standard - Young Leiter  
	 *        1: S_2 $ S_{n/2} Leiter */
	INTEGER_OP s_type() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_type_i() { 
		return(s_type()->s_i()); };
/* type 0:
 * k = 0 if step = 0, 1
 * k = step / 2 if EVEN(step) && step >= 2
 * k = (step + 1) / 2 if ODD(step) && step >= 2 
 *
 * type 1:
 * k = step / 4
 */
	INTEGER_OP s_k() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_k_i() { 
		return(s_k()->s_i()); };
/* type 0:
 * fDown TRUE if step == 1 or (step >= 2 and EVEN(step))
 * fDown FALSE if step == 0 or (step >= 3 and ODD(step)) 
 * 
 * type 1:
 * step % 4 = 0: fDown FALSE
 * step % 4 = 1: fDown TRUE
 * step % 4 = 2: fDown TRUE
 * step % 4 = 3: fDown FALSE
 */
	INTEGER_OP s_fDown() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_fDown_i() { 
		return(s_fDown()->s_i()); };
	INTEGER_OP s_reduce_generators_mode() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_reduce_generators_mode_i() { 
		return(s_reduce_generators_mode()->s_i()); };
/* |T| = [Lim1 : Li] if fDown
 * |T| = [Li : Lim1] else */
	VECTOR_OP s_T() { 
		return((VECTOR_OP)s_i(5)); };
	SYM_OP s_T_i(INT i) { 
		return(s_T()->s_i(i)); };
	MATRIX_OP s_TDidx() { 
		return((MATRIX_OP)s_i(6)); };
	INTEGER_OP s_TDidx_ij(INT i, INT j) { 
		return((INTEGER_OP)s_TDidx()->s_ij(i, j)); };
	INT s_TDidx_iji(INT i, INT j) { 
		return(s_TDidx()->s_iji(i, j)); };
	MATRIX_OP s_TDfusel() { 
		return((MATRIX_OP)s_i(7)); };
	SYM_OP s_TDfusel_ij(INT i, INT j) { 
		return(s_TDfusel()->s_ij(i, j)); };
/* TDidx / TDfusel sind im 
 * Aufwaertsschritt einzeilig. */
	
	VECTOR_OP s_D() { 
		return((VECTOR_OP)s_i(8)); };
	SYM_OP s_D_i(INT i) { 
		return(s_D()->s_i(i)); };
	VECTOR_OP s_Ad() { 
		return((VECTOR_OP)s_i(9)); };
	VECTOR_OP s_Ad_i(INT i) { 
		return((VECTOR_OP)s_Ad()->s_i(i)); };
	SYM_OP s_Ad_ij(INT i, INT j) { 
		return(s_Ad_i(i)->s_i(j)); };
	INTEGER_OP s_omega() { 
		return((INTEGER_OP)s_i(10)); };
	INT s_omega_i() { 
		return(s_omega()->s_i()); };
	VECTOR_OP s_oiti() { 
		return((VECTOR_OP)s_i(11)); };
	INTEGER_OP s_oiti_i(INT i) { 
		return((INTEGER_OP) s_oiti()->s_i(i)); };
	INT s_oiti_ii(INT i) { 
		return(s_oiti_i(i)->s_i()); };

INT field_name(INT i, INT j, BYTE *str);
INT initialize_Young(VECTOR_OP T, INT n, INT i, INT type, void *data);
INT initialize_Wreath(VECTOR_OP T, INT n, INT i, INT type, void *data);
INT initialize_arbitrary(VECTOR_OP T, INT n, INT i, 
	INT omega, VECTOR_OP oiti, INT type, void *data);
INT do_dc_young_downstep(SYM_OP m, SYM_OP dv, INT n, INT *idx, 
	INT type, void *data);
INT do_downstep(DCY_OP dc_last, SYM_OP id, INT deg, INT step, INT f_verbose, 
	INT type, void *data);
INT do_upstep(DCY_OP dc0, SYM_OP id, INT deg, INT step, INT f_verbose, 
	INT type, void *data);
};

/* graphical.C: */
INT graphical_design(INT t, INT k, INT print_nrow, INT print_ncol);
INT print_graphical_KM_matrix(MATRIX_OP N, MATRIX_OP N1, MATRIX_OP N2, INT t, INT k, 
	VECTOR_OP supp_t, VECTOR_OP aut_t, 
	VECTOR_OP supp_k, VECTOR_OP aut_k);
INT print_graphical_KM_matrix_splitted(MATRIX_OP N, MATRIX_OP N1, MATRIX_OP N2, INT t, INT k, 
	VECTOR_OP supp_t, VECTOR_OP aut_t, 
	VECTOR_OP supp_k, VECTOR_OP aut_k, 
	INT nrow, INT ncol);
INT calc_graphical_KM_matrix(MATRIX_OP M, 
	MATRIX_OP N, MATRIX_OP N1, MATRIX_OP N2, INT n, 
	VECTOR_OP stab_order_t, VECTOR_OP supp_t, VECTOR_OP aut_t, 
	VECTOR_OP stab_order_k, VECTOR_OP supp_k, VECTOR_OP aut_k);
INT graphical_orbits(DCY_OB *dc, SYM_OP go, INT k, INT n, 
	VECTOR_OP stab_order, VECTOR_OP supp, VECTOR_OP aut);

/* kramer_mesner.C */
extern char *discreta_copyright_text;
void compute_KM(VECTOR_OP gen, BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, 
	INT f_strong_generators, 
	INT f_orderly_generation, 
	INT f_TDO, 
	INT f_extension_construction, 
	INT f_canonical_representatives, 
	INT f_k_k2_design, 
	INT k2);
INT calc_kramer_mesner_matrix(LABRA_OP labG, VECTOR_OP gen, SYM_OP go, INT deg, 
	BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, 
	INT f_k_k2_design, INT k2, 
	INT f_TDO, 
	INT f_extension_construction, 
	INT f_canonical_representatives);
INT design_orbits(MATRIX_OP X, VECTOR_OP orbits, INT f_complement);
INT design_orbits_vector(VECTOR_OP X, VECTOR_OP orbits, INT f_complement);
INT dc_calc_orbit_length(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP orbit_length);
INT dc_calc_stab_go_and_K(DCY_OP dc, INT k_max, VECTOR_OP stab_go, VECTOR_OP K);
INT dc_RR_to_K(VECTOR_OP RR, VECTOR_OP K);
INT dc_step_to_k(INT step);
INT dc_k_to_step(INT k);
INT dc_print_k_set(VECTOR_OP R, SYM_OP stab_go, SYM_OP go);
INT dc_get_k_set(PERMUTATION_OP d, VECTOR_OP R, INT k, INT f_v);
INT apply_perm_set(VECTOR_OP R, PERMUTATION_OP p, INT f_reorder);
INT design_calc_blocks(VECTOR_OP gen, VECTOR_OP RR, INT k, MATRIX_OP X, 
	VECTOR_OP O, VECTOR_OP O_first, VECTOR_OP O_len, INT f_v, INT f_vv);
INT dc_calc_set_orbit(VECTOR_OP gen, VECTOR_OP R, 
	VECTOR_OP O, INT f_v, INT f_vv);
INT dc_calc_go(DCY_OP dc0, SYM_OP go);
INT dc_calc_stab_go(DCY_OP dc, INT k_max, VECTOR_OP stab_go);
INT dc_orbit_size(DCY_OP dc, SYM_OP go, INT step, INT dc_no, INT f_v);
INT dc_print_dc(DCY_OP dc, SYM_OP go, INT step, INT dc_no);
INT dc_print_dc_(DCY_OP dc, SYM_OP go, INT step, INT dc_no, 
	INT f_d0, PERMUTATION_OP d0, PERMUTATION_OP d1, INT step2);
INT dc_calc_Ainf_t(DCY_OP dc0, INT up_to_step, 
	MATRIX_OP Ainf_t, /* VECTOR_OP layer, */ INT f_verbose);
INT dc_dc_no_to_dc_idx(DCY_OP dc0, INT dc_no, INT *step, INT *dc_idx);
INT dc_dc_idx_to_dc_no(DCY_OP dc0, INT step, INT dc_idx, INT *dc_no);
INT dc_calc_nb_d(DCY_OP dc0, INT up_to_step);
INT dc_calc_nb_d_via_MM(VECTOR_OP MM, INT k_max);
INT dc_Mtk(DCY_OP dc0, INT t, INT k, MATRIX_OP M);
INT dc_Mttp1(DCY_OP dc0, INT t, MATRIX_OP M);
INT km_normalizer_action_on_orbits(LABRA_OP labra_G, 
	VECTOR_OP Reps, PERMUTATION_OP p, PERMUTATION_OP q, INT f_v);
INT dc_calc_canonical_representatives(LABRA_OP labra_G, 
	VECTOR_OP RR, VECTOR_OP stab_go, INT f_v);
INT fuse_orbits(LABRA_OP labra_G, 
	VECTOR_OP Reps, VECTOR_OP new_rep_idx, VECTOR_OP new_reps, INT f_v);
INT build_KM_matrix_s_sp1(MATRIX_OP M, VECTOR_OP Orbits, 
	VECTOR_OP Orbits_above1, VECTOR_OP Orbits_above2, INT s);
INT orbits_below_to_orbits_above(VECTOR_OP Orbit_lengths, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2, 
	VECTOR_OP Orbits_above1, VECTOR_OP Orbits_above2, 
	SYM_OP go, INT s);
INT print_orbit_info(INT t, INT k, 
	VECTOR_OP Orbits, VECTOR_OP Orbit_lengths, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2);


// km_file.C:

#define LOC_MCKAY "discreta_mckay"
#define LOC_LLL "discreta_lll"
#define LOC_LARGESET1 "discreta_ls1"
#define LOC_LARGESET2 "discreta_ls2"
#define LOC_SPREAD "discreta_spread"
#define LOC_DANCE "discreta_dance"

INT km_read_ascii(BYTE *KM_fname, MATRIX_OP M, INT *v, INT *t, INT *k, 
	VECTOR_OP G_gen, VECTOR_OP RR, VECTOR_OP MM, VECTOR_OP stab_go, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2);
INT km_read_ascii_vtk(BYTE *KM_fname, INT *v, INT *t, INT *k);
INT km_write_ascii_G_gen(BYTE *km_fname, VECTOR_OP G_gen);
INT km_write_ascii_representatives(BYTE *km_fname, VECTOR_OP RR);
INT km_write_ascii_stab_go(BYTE *km_fname, VECTOR_OP stab_go);
INT km_write_stab_go_k_sets(BYTE *km_fname, VECTOR_OP stab_go, INT k);
INT km_write_ascii_KM_matrices(BYTE *km_fname, VECTOR_OP MM);
INT km_write_ascii_below(BYTE *km_fname, VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2);
INT km_print_M_asc(BYTE *g_label, BYTE *g_label_tex, BYTE *km_fname, 
	VECTOR_OP G_gen, SYM_OP go, INT deg, 
	MATRIX_OP M, INT gl_t, INT gl_k, INT f_k2, INT k2);
INT km_init_file_names(BYTE *g_label, BYTE *km_fname, INT t, INT k);
INT km_compute_TDO_decomposition(BYTE *g_label, INT t, INT k, MATRIX_OP Mtk);
void get_solutions(BYTE *KM_fname, INT lambda);
INT km_nb_of_solutions(BYTE *KM_fname, INT lambda);
INT km_read_until_lambdaend(FILE *fp);
INT km_get_solutions(BYTE *KM_fname, 
	INT lambda, INT from, INT len, VECTOR_OP S);
void check_solutions(BYTE *KM_fname, INT lambda);
void do_mckay(BYTE *KM_fname, INT lambda);
void do_spread(INT f_silent, BYTE *KM_fname, INT lambda);
void do_dance(INT f_silent, BYTE *KM_fname, INT lambda);
void do_LLL(INT f_silent, BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT lambda, INT f_iterate, INT nb_iterate);
void do_largeset1(BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT N);
void do_largeset2(BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT max_rounds, INT N, 
	INT max_solutions, INT max_loops);
INT show_km_matrix(BYTE *KM_fname);

// plesken_mtx.C:
INT dc_plesken_layer_idx(VECTOR_OP K_first, VECTOR_OP K_len, INT i);
INT dc_plesken_prepare(INT k_min, INT k_max, VECTOR_OP MM, 
	MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT f_v, INT f_vv);
INT dc_plesken_matrix_prepare(INT k_min, INT k_max, VECTOR_OP MM, 
	MATRIX_OP Ainf_t, MATRIX_OP coeff, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT f_v, INT f_vv);
INT dc_plesken_matrices_prepare(INT k_min, INT k_max, VECTOR_OP MM, 
	MATRIX_OP Ainf_block_wise, MATRIX_OP Ainf_t, MATRIX_OP Ainf_t_inv, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT f_v, INT f_vv);
INT dc_Ainf_inv_coeff_ma(INT k_max, MATRIX_OP coeff, INT f_v);
INT dc_Ainf_block_wise_via_MM(INT k_min, INT k_max, MATRIX_OP Ainf_block_wise, VECTOR_OP MM, INT f_v);
INT dc_Mtk_via_Mtr_Mrk(INT t, INT r, INT k, MATRIX_OP Mtr, MATRIX_OP Mrk, MATRIX_OP Mtk, INT f_v);
INT dc_calc_K_first_len(VECTOR_OP K, INT k_max, VECTOR_OP K_first, VECTOR_OP K_len);
INT dc_calc_K_first_len_from_MM(VECTOR_OP MM, INT k_max, VECTOR_OP K_first, VECTOR_OP K_len);
INT dc_find_first_len(VECTOR_OP K, INT t, INT *first, INT *len);



// plesken_ring.C:
INT plesken_product_block_wise(MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, VECTOR_OP args, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT k_min, 
	VECTOR_OP type, INT f_v, INT f_vv);
INT plesken_decompose_layer(MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT k_min, INT k_max, INT k, INT i, SYM_OP ni, 
	VECTOR_OP aiaj, INT f_v, INT f_vv);
INT plesken_multiply_columns(MATRIX_OP Ainf_block_wise, 
	VECTOR_OP orbit_length, VECTOR_OP args, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT k_min, INT k_max, INT k, 
	VECTOR_OP aiaj, INT f_v, INT f_vv);
INT plesken_product(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, VECTOR_OP args, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT k_min, 
	VECTOR_OP type, INT f_v, INT f_vv);
INT plesken_decompose(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	VECTOR_OP aiaj, VECTOR_OP aijk, INT first, INT len, INT d0, INT f_v);
INT plesken_decompose_with_inverse(MATRIX_OP Ainf_inv, 
	VECTOR_OP aiaj, VECTOR_OP aijk, INT first, INT len, INT d0, INT f_v);
INT plesken_decomposition_collect_layer_wise(VECTOR_OP orbit_length, 
	VECTOR_OP aijk, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, VECTOR_OP type, INT f_v, INT f_vv);
INT plesken_decomposition_collect_layer(VECTOR_OP orbit_length, 
	VECTOR_OP aijk, INT first, INT len, SYM_OP sum, INT f_v);


// plesken_Ik2.C:
INT calc_intersections_Iknm(MATRIX_OP Ainf_block_wise, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	MATRIX_OP I, INT k_min, INT k, 
	VECTOR_OP Sel1, VECTOR_OP Sel2, INT f_v, INT f_vv);
INT calc_intersections_Iknn(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	MATRIX_OP Iknn, INT k_min, INT k, 
	VECTOR_OP Sel, INT f_v, INT f_vv);
INT calc_intersections_Ik2(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT k_min, 
	MATRIX_OP Ik2, INT k, INT f_v, INT f_vv);
INT get_block_intersection_type(MATRIX_OP Ik2, INT k, 
	VECTOR_OP design_orbits, INT do_idx, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP type);
INT intersections_Iknm_blow_up(MATRIX_OP I, VECTOR_OP Sel1, VECTOR_OP Sel2, 
	INT k, VECTOR_OP I_vec, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v);
INT intersections_Iknn_blow_up(MATRIX_OP I, VECTOR_OP Sel, 
	INT k, VECTOR_OP I_vec, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v);
INT Iknm_get_intersections_of_2(MATRIX_OP I, VECTOR_OP Sel1, VECTOR_OP Sel2, 
	INT i1, INT i2, 
	INT num_layers, 
	VECTOR_OP type);
INT Iknn_get_intersections_of_2(MATRIX_OP I, VECTOR_OP Sel, 
	INT k, INT i1, INT i2, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP type);
INT intersections_Ik2_blow_up(MATRIX_OP Ik2, INT k, VECTOR_OP Ik2_vec, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v);
INT get_intersections_of_2(MATRIX_OP Ik2, INT k, INT i1, INT i2, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP type);
INT Ik2_print(MATRIX_OP Ik2, INT k, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, INT f_v);
INT calc_intersections_Ik3(MATRIX_OP Ainf, MATRIX_OP coeff, 
	VECTOR_OP orbit_length, 
	VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	MATRIX_OP Ik3, INT k, INT f_v, INT f_vv);



// plesken.C:
INT dc_calc_Ainf_t_via_MM(INT k_max, MATRIX_OP Ainf_t, VECTOR_OP MM, INT f_verbose);
	// outdated !
INT dc_Mtk_via_MM(INT t, INT k, MATRIX_OP M, VECTOR_OP MM, INT f_v);
	// outdated !

INT design_2intersections(INT k, MATRIX_OP Ik2, MATRIX_OP Ainf, MATRIX_OP Ainf_inv, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, VECTOR_OP orbits, 
	INT f_multiplicities, 
	INT f_use_complement, VECTOR_OP orbits_c, 
	VECTOR_OP type2, INT f_print_block_types);
INT print_block_intersection_types(VECTOR_OP block_types_sorted, 
	VECTOR_OP orbit_idx, INT highest_layer, INT num_layers);
INT print_gl_type2_intersections(VECTOR_OP type, INT f_multiplicities, INT f_complements);
INT print_type_with_multiplicities(VECTOR_OP gl_type, 
	VECTOR_OP type_vec, VECTOR_OP mult_vec);
INT print_invariants_and_classes(VECTOR_OP types_sorted, VECTOR_OP classes, 
	INT f_multiplicities, INT f_complements);
INT global_intersection_2_sets(INT k_max, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP block_types_sorted, VECTOR_OP orbit_idx, 
	INT f_multiplicities, VECTOR_OP gl_type_2_sets, INT f_v);
INT all_block_intersection_types(INT k_max, 
	MATRIX_OP Ik2, MATRIX_OP Ainf, MATRIX_OP Ainf_inv, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	VECTOR_OP design_orbits, 
	VECTOR_OP block_types_sorted, VECTOR_OP orbit_idx, INT f_v, INT f_vv);
INT block_intersection_type(INT k_max, MATRIX_OP Ainf, MATRIX_OP Ainf_inv, 
	VECTOR_OP orbit_length, VECTOR_OP K_first, VECTOR_OP K_len, 
	INT highest_layer, INT num_layers, 
	INT do_idx, VECTOR_OP design_orbits, 
	VECTOR_OP type, INT f_v, INT f_vv);





/* parameter.C: */
INT calc_lambda_max(INT t, INT v, INT k, INT lambda, SYM_OP l_max);
INT calc_lambda_ijs_matrix(INT t, INT v, INT k, INT lambda, INT s, MATRIX_OP M);
INT calc_lambda_ijs(INT t, INT v, INT k, INT lambda, INT s, INT i, INT j, SYM_OP lijs);
INT calc_lambda_ij(INT t, INT v, INT k, INT lambda, INT i, INT j, SYM_OP lij);
INT calc_mendelsohn_coefficient_matrix(INT t, INT m, MATRIX_OP M);
INT calc_mendelsohn_rhs(INT v, INT t, INT k, INT lambda, INT m, INT s, VECTOR_OP rhs);
INT calc_and_print_design_parameter(INT t, INT v, INT k, INT lambda);
INT design_calc_b(INT t, INT v, INT k, INT lambda, SYM_OP b);
INT design_calc_r(INT v, INT k, SYM_OP b, SYM_OP r);
INT design_print_Mendelsohn_and_Koehler(INT v, INT t, INT k, MATRIX_OP Mendelsohn_mtx, 
	VECTOR_OP Mendel_RHS, MATRIX_OP LAmbda, VECTOR_OP alpha_i, 
	VECTOR_OP Coeff, VECTOR_OP Constant_term, 
	VECTOR_OP Koehler_Constant_term_v, VECTOR_OP Koehler_Coeff_v);
INT solve_Mendelsohn(INT v, INT t, MATRIX_OP M, 
	VECTOR_OP RHS, VECTOR_OP alpha_i);
INT Mendelsohn_generalized_RHS(INT v, INT t, INT k, INT m, INT s, INT lambda, 
	VECTOR_OP RHS, INT f_v);
INT Mendelsohn(INT v, INT t, INT k, INT m, INT lambda, 
	MATRIX_OP M, VECTOR_OP RHS, MATRIX_OP LAmbda, INT f_v);
INT calc_LAmbda(INT v, INT t, INT k, INT lambda, MATRIX_OP LAmbda);
INT calc_Lambda(INT v, INT t, INT k, INT lambda, VECTOR_OP Lambda);
INT calc_delta_lambda(INT v, INT t, INT k);
INT Koehler_eqns_for_blocks(INT v, INT t, INT k, INT lambda, INT s_max, 
	INT max_intersection, 
	VECTOR_OP Coeff, VECTOR_OP Constant_term, BYTE *variable_name, INT f_v);
INT Koehler(INT v, INT t, INT k, INT lambda, INT m, INT s_max, 
        VECTOR_OP Coeff, VECTOR_OP Constant_term, INT f_v);
INT Koehler_j(INT v, INT t, INT m, INT s_max, INT j, VECTOR_OP Lambda, 
	VECTOR_OP constant_term, VECTOR_OP coeff);
INT Koehler_eqn_print(INT v, INT t, INT m, INT s_max, INT j, 
	VECTOR_OP constant_term, VECTOR_OP coeff, 
	BYTE *variable_name, BYTE *argument);
void design_print_solution_vector(BYTE *buf, INT no, INT width, BYTE *offset);
void design_print_solution_vector_numerical(BYTE *buf, INT no, BYTE *offset);


#define DP_RULE_COMPLEMENT 1
#define DP_RULE_REDUCED_T 2
#define DP_RULE_DERIVED 3
#define DP_RULE_RESIDUAL 4
#define DP_RULE_ALLTOP 5

#define DP_RULE_COMPL_REDUCED_T 6
#define DP_RULE_COMPL_DERIVED 7
#define DP_RULE_COMPL_RESIDUAL 8
#define DP_RULE_COMPL_ALLTOP 9

#define DP_RULE_TRUNG_SUPPLEMENTARY 10
#define DP_RULE_SUPPLEMENTARY 11
#define DP_RULE_TRUNG 12

class design_parameter_source_ob : public VECTOR_OB {
public:
	INTEGER_OP s_prev() { 
		return((INTEGER_OP) s_i(0)); };
	INT s_prev_i() { 
		return(s_prev()->s_i()); };
	INTEGER_OP s_rule() { 
		return((INTEGER_OP) s_i(1)); };
	INT s_rule_i() { 
		return(s_rule()->s_i()); };
	STRING_OP s_comment() { 
		return((STRING_OP) s_i(2)); };
	BYTE *s_comment_s() { 
		return(s_comment()->s_str()); };
	VECTOR_OP s_references() { 
		return((VECTOR_OP) s_i(3)); };
	STRING_OP s_references_i(INT i) { 
		return((STRING_OP) s_references()->s_i(i)); };
	BYTE *s_references_is(INT i) { 
		return(s_references_i(i)->s_str()); };

	INT field_name(INT i, INT j, BYTE *str);
	INT init();
	INT sprint(BYTE *s);
	INT sprint_text012(BYTE *s0, BYTE *s1, BYTE *s2);
};

class design_parameter_ob : public VECTOR_OB {
public:
	INTEGER_OP s_id() { 
		return((INTEGER_OP) s_i(0)); };
	INT s_id_i() { 
		return(s_id()->s_i()); };
	INTEGER_OP s_v() { 
		return((INTEGER_OP) s_i(1)); };
	INT s_v_i() { 
		return(s_v()->s_i()); };
	INTEGER_OP s_t() { 
		return((INTEGER_OP) s_i(2)); };
	INT s_t_i() { 
		return(s_t()->s_i()); };
	INTEGER_OP s_k() { 
		return((INTEGER_OP) s_i(3)); };
	INT s_k_i() { 
		return(s_k()->s_i()); };
	INTEGER_OP s_lambda() { 
		return((INTEGER_OP) s_i(4)); };
	INT s_lambda_i() { 
		return(s_lambda()->s_i()); };
	VECTOR_OP s_source() { 
		return((VECTOR_OP) s_i(5)); };
	DESIGN_PARAMETER_SOURCE_OP s_source_i(INT i) { 
		return((DESIGN_PARAMETER_SOURCE_OP) s_source()->s_i(i)); };

#if 0
	INTEGER_OP s_prev() { 
		return((INTEGER_OP) s_i(5)); };
	INT s_prev_i() { 
		return(s_prev()->s_i()); };
	INTEGER_OP s_rule() { 
		return((INTEGER_OP) s_i(6)); };
	INT s_rule_i() { 
		return(s_rule()->s_i()); };
	STRING_OP s_comment() { 
		return((STRING_OP) s_i(7)); };
	BYTE *s_comment_s() { 
		return(s_comment()->s_str()); };
#endif

	INT field_name(INT i, INT j, BYTE *str);
	INT init();
	INT sprint(BYTE *s);
	INT selfcomplementary();
	INT lambda_compl();
	INT complement(DESIGN_PARAMETER_OP p);
	INT reduce_t(DESIGN_PARAMETER_OP p);
	INT compl_reduce_t(DESIGN_PARAMETER_OP p);
	INT derive(DESIGN_PARAMETER_OP p);
	INT compl_derive(DESIGN_PARAMETER_OP p);
	INT residual(DESIGN_PARAMETER_OP p);
	INT compl_residual(DESIGN_PARAMETER_OP p);
	INT trung_supplementary(DESIGN_PARAMETER_OP p);
	INT supplementary(DESIGN_PARAMETER_OP p);
	INT alltop(DESIGN_PARAMETER_OP p);
};

INT design_parameter_produce(VECTOR_OP V, 
	INT t, INT v, INT k, INT lambda, BYTE *comment, INT f_v);
INT dp_search(VECTOR_OP V, INT len, DESIGN_PARAMETER_OP p);


INT i_bitvec_db(DATABASE_OP *db, BYTE *db_prefix);
void e_bitvec_db(DATABASE_OP db);
INT do_db_bitvec_create_and_close(BYTE *db_prefix);
INT read_01_file(BYTE *db_prefix, BYTE *fname);

INT i_db_dp(DATABASE_OP *p_db, BYTE *path);
void e_db_dp(DATABASE_OP db);
INT do_db_dp_create_and_close(BYTE *path);
INT do_db_dp_info(BYTE *path, INT btree_idx);
INT db_dp_highest_id(DATABASE_OP db, INT *highest_id);
INT db_dp_load_id_data(DATABASE_OP db, INT id, DESIGN_PARAMETER_OP p, DATATYPE *data);
INT db_dp_load_id(DATABASE_OP db, INT id, DESIGN_PARAMETER_OP p);
INT db_dp_search(DATABASE_OP db, 
	INT t, INT v, INT k, INT lambda, 
	INT *idx, INT *len, INT *f_found);
INT db_dp_parameter_add(DATABASE_OP db, 
	DESIGN_PARAMETER_OP p, INT *id, INT f_verbose, 
	INT *f_added, FILE *fp_tex);
INT db_dp_design_parameter_produce(DATABASE_OP db, 
	INT t, INT v, INT k, INT lambda, BYTE *comment, INT f_v, FILE *fp_tex, 
	MATRIX_OP data);
INT db_dp_HTML(BYTE *path);
INT db_dp_read_design_txt(BYTE *fname_txt, BYTE *fname_tex, 
	BYTE *db_path, INT f_clear_always, VECTOR_OP vec_data);



/* ladder2.C: */
INT dc_get_transversals(DCY_OP dc, SYM_OP go, INT n, INT step, INT dc_no, 
	VECTOR_OP dc_rep, VECTOR_OP dc_idx, VECTOR_OP dc_trans_idx, 
	VECTOR_OP dc_omega, VECTOR_OP oiti, VECTOR_OP trans, VECTOR_OP trans_idx, INT f_v);
INT dc_find_path(DCY_OP dc, SYM_OP go, INT step, INT dc_no, 
	VECTOR_OP dc_rep, VECTOR_OP dc_idx, VECTOR_OP dc_trans_idx, INT f_v);
INT dc2_print_k_set(VECTOR_OP R, SYM_OP stab_go, SYM_OP ol, INT alpha);
INT dc2_calc_intersection_indicator(DCY_OP dc, SYM_OP go, INT step, 
	PERMUTATION_OP d0, PERMUTATION_OP d1, INT k2, UNIPOLY_OP ii);
INT dc2_print_orbits(DCY_OP dc, SYM_OP go, INT step, 
	PERMUTATION_OP d0, PERMUTATION_OP d1, INT k2);
INT dc2_get_orbit_length_and_intersection(DCY_OP dc, SYM_OP go, INT step, INT dc_no, 
	PERMUTATION_OP d0, PERMUTATION_OP d1, INT k2, 
	VECTOR_OP R, SYM_OP ago, SYM_OP ol, INT *alpha);
INT dc2_intersection(PERMUTATION_OP d, 
	PERMUTATION_OP d1, PERMUTATION_OP d2, INT k1, INT k2);
INT dc2_calc_intersection_M(DCY_OP dc, SYM_OP go, INT k_max, MATRIX_OP M);
INT dc2_calc_intersection_M_tk(DCY_OP dc, SYM_OP go, INT t, INT k, MATRIX_OP M);
INT intersections_by_enumeration(DCY_OP dc, SYM_OP go, INT n, 
	INT s1, INT dcno1, INT s2, INT dcno2, UNIPOLY_OP ii);
INT double_cosets2(DCY_OP dc, DCY_OP dc2, SYM_OP go, INT n, 
	INT s1, INT dcno1, INT s2, INT dcno2, UNIPOLY_OP ii);

/* intersection.C: */
INT manage_block_alpha_i_data(VECTOR_OP block_alpha_i, 
	INTEGER_OP bai_len, 
	VECTOR_OP bai_ref, VECTOR_OP bai_refv, 
	VECTOR_OP bai_data, INT *ref, INT f_v);
INT total_intersection_3_tuples(MATRIX_OP X, INT k,
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, 
	VECTOR_OP intersection_3_tuples);
INT total_intersection_3_sets(MATRIX_OP X, INT k,
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, 
	VECTOR_OP intersection_3_sets);
INT intersection_of3_and_add(INT i, INT j, INT k, 
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, 
	VECTOR_OP intersection_3_sets, INT sign);
INT total_intersection_2_tuples_to_2_sets(VECTOR_OP intersection_2_tuples, 
	VECTOR_OP intersection_2_sets);
INT total_intersections_2_tuples(MATRIX_OP X, INT k, 
	SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	VECTOR_OP block_types, VECTOR_OP intersection_2_tuples, 
	INTEGER_OP bai_len, VECTOR_OP bai_ref, VECTOR_OP bai_refv,
	VECTOR_OP bai_data, INT f_v);
INT vec_sum_of_all_elements(VECTOR_OP V, SYM_OP S);
INT check_intersection_eqns(MATRIX_OP X, MATRIX_OP Mendelsohn, 
	VECTOR_OP Eqns, VECTOR_OP RHS, VECTOR_OP block_types, 
	INTEGER_OP bai_len, VECTOR_OP bai_ref, VECTOR_OP bai_refv, 
	VECTOR_OP bai_data, INT f_v);
INT check_intersection_eqn(INT block_orbit, MATRIX_OP X, MATRIX_OP Mendelsohn, 
	MATRIX_OP eqn, VECTOR_OP RHS, VECTOR_OP block_alpha_i, INT f_v, INT f_vv);
INT eqns_apply(VECTOR_OP Eqns, MATRIX_OP M);
INT intersection_eqns(VECTOR_OP Eqns, MATRIX_OP Ikk, INT k, INT f_v);
INT Intersection_Mtk_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	INT t, INT k, MATRIX_OP Mtk, MATRIX_OP Ainf, INT f_v, INT f_vv);
INT Intersection_M_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP M, MATRIX_OP Ainf, INT f_v, INT f_vv);
INT Intersection_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, MATRIX_OP Ainf, 
	INT dcno1, INT dcno2, UNIPOLY_OP ii, INT f_v, INT f_vv);
INT Intersection_of3_via_Ainf(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, MATRIX_OP Ainf, 
	INT i1, INT i2, INT i3, UNIPOLY_OP ii, INT f_v, INT f_vv);
INT Plesken_decompose(MATRIX_OP Ainf, VECTOR_OP aiaj, VECTOR_OP aijk, INT f_v);
INT unipoly_divide_out_oli(UNIPOLY_OP ii, SYM_OP go, VECTOR_OP stab_go, INT i);
INT collect_intersections_by_layer(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP Ainf, VECTOR_OP aijk, UNIPOLY_OP ii, INT f_v);

/* intersection_aijk.C: */
INT Intersection_M_via_Aijk(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, 
	MATRIX_OP M, MATRIX_OP Aijk);
INT Intersection_via_Aijk(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP K, MATRIX_OP Aijk, 
	INT dcno1, INT dcno2, UNIPOLY_OP ii);
INT intersection_M_via_Aijk(DCY_OP dc, SYM_OP go, INT k_max, MATRIX_OP M, 
	MATRIX_OP Aijk, INT nb_d);
INT intersection_via_Aijk(DCY_OP dc, SYM_OP go, MATRIX_OP Aijk, INT nb_dc, 
	INT dcno1, INT dcno2, UNIPOLY_OP ii);
INT intersection_M_tk_via_Aijk(DCY_OP dc, SYM_OP go, INT t, INT k, 
	MATRIX_OP M, MATRIX_OP Aijk, INT nb_d);
INT calc_Aijk(MATRIX_OP Ainf, MATRIX_OP Aijk);

/* report.C: */
void design_report(BYTE *KM_fname, INT f_html, INT f_select, 
	INT design_select_lambda, INT select_first, INT select_length);
void KM_fname_extract_group_label(BYTE *KM_fname, BYTE *group_label);
void design_report_plesken(BYTE *KM_fname);
INT do_report(BYTE *KM_fname, INT f_html, INT f_selected, 
	INT select_lambda, INT select_first, INT select_len);

/* ladder.C */
INT dc_young_downstep(
	SYM_OP m, SYM_OP dv, 
	INT n, INT k, INT *idx, 
	INT type, void *data);
INT dc_in_test_young(
	SYM_OP pi, 
	INT n, INT k, INT fDown, 
	INT type, void *data);
INT dc_do_step(DCY_OP dc0, 
	INT step, SYM_OP id, INT n, INT f_verbose, 
	INT type, void *data);
INT dc_trace_coset(
	DCY_OP dc0, SYM_OP g, 
	INT i, INT *dc_idx, SYM_OP fusel, 
	INT type, void *data);
INT vec_translate(
	VECTOR_OP V, SYM_OP x, VECTOR_OP Vx, 
	INT type, void *data);
INT vec_multiples(
	VECTOR_OP V, SYM_OP p);

/* mathieu.C: */
INT construct_W24(VECTOR_OP W24_blocks, 
	VECTOR_OP M24_gen, VECTOR_OP stab_gen, INT *block_idx, INT f_v);
INT construct_W23(VECTOR_OP W23_blocks, 
	VECTOR_OP M23_gen, VECTOR_OP stab_gen, INT *block_idx, INT f_v);
INT construct_W22(VECTOR_OP W22_blocks);
INT construct_graph_from_design(INT v, VECTOR_OP blocks, MATRIX_OP G, INT i1);
INT construct_two_graph(MATRIX_OP G, VECTOR_OP three_blocks);
INT get_graph_value(MATRIX_OP G, INT i, INT j);
INT blocks_disjoint(VECTOR_OP B1, VECTOR_OP B2);
INT block_intersection(VECTOR_OP B1, VECTOR_OP B2);
INT McKayWreath();
INT virasoro(INT f_v);

/* higman_sims.C: */
INT construct_Higman_Sims_176(VECTOR_OP HS_gen);
INT construct_Higman_design(VECTOR_OP W24_blocks, INT point_a, INT point_b, 
	VECTOR_OP H_point_idx, VECTOR_OP H_block_idx, MATRIX_OP I);
INT construct_Higman_Sims_graph(MATRIX_OP G);
INT construct_Higman_Sims_graph2(VECTOR_OP W22_blocks, MATRIX_OP G);

/* conway.C: */
INT construct_Conway3_276(VECTOR_OP Co3_gen);

/* ladder_info.C: */
LADDER_INFO *init_ladder_info(
	VECTOR_OP generators, SYM_OP go, INT deg, 
	BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, INT lambda, 
	INT f_verbose, 
	BYTE *generators_fname);
void free_ladder_info(LADDER_INFO *li);
INT li_init_file_names(LADDER_INFO *li);
void li_message(LADDER_INFO *li);
INT li_init_transversals(LADDER_INFO *li);
INT li_leiterspiel(LADDER_INFO *li);
INT li_Mtk(LADDER_INFO *li, 
	MATRIX_OP M, INT t, INT k);
INT li_stab_go(LADDER_INFO *li, 
	VECTOR_OP V, INT t);
INT li_perm_rep(LADDER_INFO *li, INT t);
INT li_set_reps(LADDER_INFO *li, 
	VECTOR_OP V, INT t);
INT li_print_dc_info(LADDER_INFO *li, 
	VECTOR_OP stab_go, VECTOR_OP reps, 
	INT t, FILE *fp, INT f_verbose);
INT li_print_M_asc(LADDER_INFO *li, 
	MATRIX_OP M, INT t, INT k, INT f_k2, INT k2);

/* dc_draw.h */

#define MAX_ORBITS 512

struct dc_lattice {
	MEM_OP plaz;
	MEM_OP o_dx;
	double extrema[6];
	INT n;
	INT n2;
	INT up_to_step;
	INT nb_d;
	INT d_layer[MAX_ORBITS];
	INT d_layer_idx[MAX_ORBITS];
		/* d_layer[i] und d_layer_idx[i]:
		 * i durchlaeft die Doppelnebenklassen 
		 * aller layers (stufen):
		 * layer der DNK und 
		 * laufende Nummer innerhalb dieses 
		 * layers;
		 * die layer werden von 
		 * step 1 beginnend bis up_to_step 
		 * durchlaufen (nur ODD(step)). 
		 * d_layer[0] = 0; d_layer_idx[0] = 0; */
	
	INT type;
	void *data;
	DCY_OP dc;
	MATRIX_OP Ainf_t;
	VECTOR_OP nl;
	VECTOR_OP orbit_size;
	INT f_with_perm;
	VDEVICE vdev;
};

#define DC_MAX_GRAPH_N 16
#define DC_GRAPH_RAD 15 
#define DC_GRAPH_OFFSET_X 0
// #define DC_GRAPH_OFFSET_X -30
#define DC_GRAPH_OFFSET_AIJ 20
#define DC_X_PIX 1000
#define DC_Y_PIX 1000

DC_LATTICE *dc_graphs(INT deg, BYTE *path);
DC_LATTICE *open_dcl(DCY_OP dc0, INT n, INT n2, 
	INT up_to_step, INT f_with_perm, 
	INT f_verbose, INT width_in_pages, 
	BYTE *file_name, INT type, void *data);
INT free_dcl(DC_LATTICE *p);
INT dcl_draw_func(void *dc_lattice, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev);
INT dcl_draw_orbit(DC_LATTICE *p, INT orbit_idx, INT at_x, INT at_y, INT dx, 
	INT width_in_pages, INT *gra_x, INT *gra_y, VDEVICE *vdev);
INT dcl_print_graph(DC_LATTICE *p, VECTOR_OP R, INT pix_x, INT pix_y, 
	INT width_in_pages, INT *gra_x, INT *gra_y, VDEVICE *vdev);
INT dcl_draw_legende(DC_LATTICE *p, 
	INT width_in_pages, INT *gra_x, INT *gra_y, VDEVICE *vdev);
INT dcl_calc_xy(INT n, INT *gra_x, INT *gra_y);
void dcl_3text(VDEVICE *vdev, INT x, INT y, BYTE *text1, BYTE *text2, BYTE *text3);
void dcl_halb(VDEVICE *vdev, INT x, INT y, BYTE *text);
INT dcl_line(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1);

// orderly:

struct orderly_sets {
	INT n, k;
	INT *set; // [k+1]
	INT l; // current length of the set;
	INT *nb_digits; // [k+1]
	INT *cur_digit; // [k+1]
	INT *digits; // [(k + 1) * n]
	LABRA_OP G;
	MATRIX_OP TG;
	SYM_OP go;
	INT f_stabilizer_action_verbose;
};

ORDERLY_SETS *orderly_sets_open(INT k, LABRA_OP G, INT f_stabilizer_action_verbose);
INT orderly_sets_close(ORDERLY_SETS * O);
INT orderly_sets_is_canonical(ORDERLY_SETS *O, LABRA_OP aut, SYM_OP ago);
INT orderly_sets_first(ORDERLY_SETS *O, LABRA_OP stab, SYM_OP stab_order);
INT orderly_sets_next(ORDERLY_SETS *O, LABRA_OP stab, SYM_OP stab_order);
INT orderly_sets_first_set(ORDERLY_SETS *O);
INT orderly_sets_next_set(ORDERLY_SETS *O);
INT orderly_sets_next_set_no_descend(ORDERLY_SETS *O);
INT orderly_sets_first_digit(ORDERLY_SETS *O, INT i);
INT orderly_sets_next_digit(ORDERLY_SETS *O, INT i);
INT orderly_sets_compute_good_digits(ORDERLY_SETS *O, LABRA_OP aut);

#endif /* LADDER_INCLUDED */

