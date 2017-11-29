/* ma.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef MA_INCLUDED
#define MA_INCLUDED

#ifndef GFQ_INCLUDED
#include <DISCRETA/gfq.h> // for PERM_REP
#endif

#define MATRIX_OVERSIZE 16

class matrix_ob : public SYM_OB {
public:
	INT freeself();
	
	// creation:
	INT m_lh(INTEGER_OP len, INTEGER_OP height);
	INT m_ilih(INT len, INT height);
	INT m_lh_n(INTEGER_OP len, INTEGER_OP height);
		/* Mit 0 vorbesetzen. */
	INT m_ilih_n(INT len, INT height);
		/* make_intlength_intheight_null_matrix */
		/* Mit 0 vorbesetzen. */
	
	INT b_lh(INTEGER_OP len, INTEGER_OP height);
	INT b_lhs(INTEGER_OP len, INTEGER_OP height, 
		SYM_OP self);
	INT b_lh_n(INTEGER_OP len, INTEGER_OP height);
		/* Mit 0 vorbesetzen. */
	
	// access functions: 
	INTEGER_OP s_h() { 
		return((INTEGER_OP)(
			ob_self.ob_matrix->m_height)); };
		/* height = Anzahl Zeilen. */
	INTEGER_OP s_l() { 
		return((INTEGER_OP)(
			ob_self.ob_matrix->m_length)); };
		/* length = Anzahl Spalten. */
	INT s_hi() { return(s_h()->s_i()); };
	INT s_li() { return(s_l()->s_i()); };
	void c_h(INTEGER_OP height) { 
		ob_self.ob_matrix->m_height = 
			(OP)height; };
	void c_l(INTEGER_OP length) { 
		ob_self.ob_matrix->m_length = 
			(OP)length; };

	SYM_OP s_s() { 
		return((SYM_OP)(
			ob_self.ob_matrix->m_self)); };
		/* Zeiger auf Speicher fuer l * h Objekte. */
	void c_s(SYM_OP self) { 
		ob_self.ob_matrix->m_self = (OP)self; };
	
	// access elements:
	SYM_OP s_ij(INT i, INT j) { 
		return(s_s() + i * s_li() + j); };
	INT s_iji(INT i, INT j) { 
		return(((INTEGER_OP)s_ij(i, j))->s_i()); };
	void m_iji(INT i, INT j, INT val) { 
		((INTEGER_OP)s_ij(i, j))->m_i(val); };
	void c_iji(INT i, INT j, INT val) { 
		((INTEGER_OP)s_ij(i, j))->c_i(val); };
	
	INT quadraticp() { return(s_li() == s_hi()); };
	
/* in iof.C: */
INT hip();
INT hip1();
INT calc_size_on_file();
INT write_mem(MEM_OP mem, INT debug_depth);
INT read_mem(MEM_OP mem, INT debug_depth);

INT nullp();
INT einsp();
INT compare(MATRIX_OP b);
INT add_apply_matrix(MATRIX_OP b);
INT add_apply(SYM_OP b);
INT add(MATRIX_OP b, MATRIX_OP ergeb);
INT copy(MATRIX_OP b);
INT mult_apply_scalar_matrix(SYM_OP a);
/* before: mult_apply_scalar_matrix(
 * OP a, OP b), a is scalar, b matrix */
INT mult_scalar(SYM_OP scalar, MATRIX_OP res);
/* before: mult_scalar_matrix(
 * OP scalar, OP mat, OP res) 
 * ACHTUNG: mat is now this ! */
/* res != mat notwendig ! */
/* komponentenweise multiplikation */
INT mult_matrix(MATRIX_OP b, MATRIX_OP c);
/* c = a * b; c must be different from a and b. */
/* before: mult_matrix_matrix(OP a, OP b, OP c) */
INT mult_vector(VECTOR_OP b, VECTOR_OP c);
INT mult_apply_matrix(MATRIX_OP b);
/* b = a * b */
/* before: mult_apply_matrix_matrix(OP a, OP b) */
INT mult(SYM_OP b, MATRIX_OP d);
/* before: mult_matrix(OP a, OP b, OP d) */
INT mult_apply(SYM_OP b);
/* before: mult_apply_matrix(OP a, OP b) */

INT swap_col(INT i, INT j);
/* before: change_column_ij(OP a, INT i, INT j) */
INT swap_row(INT i, INT j);
/* before: change_row_ij(OP a, INT i, INT j) */
INT zero();
INT one();
INT m_one();
INT homo_z(INT z);
INT norm_projective();
/* Normalform der projektiven Matrixgruppen:
 * Bestimme minimales j, so dass 
 * a[n-1][j] != 0 aber a[n-1][k] = 0
 * fuer k = 0, 1, .. j - 1.
 * Dann skalare Multiplikation mit a[n-1][j]^-1 */ 
INT invers(MATRIX_OP b);
INT sprint(BYTE *str);
/* appends to str. 
 * writes to maximal strlength of 200. */
INT fPrint(FILE *fp);
INT Print();
INT fprint_GAP(FILE *fp);
INT calc_column_width(VECTOR_OP col_width, INT f_latex);
INT fprint_raw(FILE *f);
INT sprint_latex(BYTE *str);
INT latex(FILE *fp);
INT sprint_latex_decomposed_vector(BYTE *str, 
	VECTOR_OP p, VECTOR_OP q);
INT sprint_latex_decomposed(BYTE *str, 
	PARTITION_OP p, PARTITION_OP q);
INT latex_lower_tri(FILE *fp);
INT latex_upper_tri(FILE *fp);
INT latex_upper_tri_block_wise(FILE *fp, INT block_size);
INT integer_print_dense();
INT transpose(MATRIX_OP At);
INT Asup2Ainf(MATRIX_OP Ainf);
INT Ainf2Asup(MATRIX_OP Asup);
INT Asup2Acover(MATRIX_OP Acover);
INT Asup2Acover_arbitrary(MATRIX_OP Acover);
INT Acover2nl(VECTOR_OP nl);

#ifdef PERMTRUE
INT perm_cols_rows(PERMUTATION_OP p, PERMUTATION_OP q, MATRIX_OP M1);
INT perm_cols(PERMUTATION_OP p, MATRIX_OP M1);
INT perm_rows(PERMUTATION_OP p, MATRIX_OP M1);
#endif

INT determinante(SYM_OP erg);
INT inc_row(void);
INT realloc_row(INT new_len);
INT inc_col(void);
INT realloc_column(INT new_len);
INT insert_row_at(INTEGER_OP len, INT i, 
	VECTOR_OP V, INT V_len);
INT search_row(INT len, INT nb_keys, 
	INT f_ascending, VECTOR_OP V, 
	INT *idx, INT *f_found);
/* V ist Vector mit mindestens nb_keys Eintraegen;
 * die matrix muss mindestens 
 * len Zeilen und nb_keys Spalten haben. */
INT split_column_wise(VECTOR_OP S, INT col_width);
INT split_row_wise(VECTOR_OP S, INT nrow);
INT cyclic_block_matrix(VECTOR_OP blocks, INT f_row_wise, INT f_forward);
INT perm_rep(PERM_REP *P, PERMUTATION_OP p);
INT matrix_of_perm_rep(PERM_REP *P, PERMUTATION_OP p);
INT complement();
INT stirling_first(INT n, INT f_extended, INT f_signless, INT f_v);
INT stirling_second(INT n, INT f_extended, INT f_ordered, INT f_v);
INT binomial(INT n, INT f_extended, INT f_inverse, INT f_v);
INT rank_p(INT p);
INT calc_hash_key(INT key_len, BYTE *hash_key, INT f_v);

// ma_geo.C:
INT automorphism_group_on_cols(PERMUTATION_OP p, 
	PERMUTATION_OP q, LABRA_OP aut_on_canonical_form);
INT canonical_form(INT f_row_perms, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, 
	VECTOR_OP DDp, VECTOR_OP DDb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	LABRA_OP aut, INT f_aut_on_canonical_form, VECTOR_OP aut_gen, 
	INT f_apply_perms_to_canonical_form, INT f_v, INT f_vv);
INT build_incidence_matrix(VECTOR_OP B, INT v, INT f_entries_start_with_1);
INT apply_perms(PERMUTATION_OP p, PERMUTATION_OP q);
INT calc_TDO(INT f_calc_second_tdo, INT f_ddp, INT f_ddb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	VECTOR_OP DDp, VECTOR_OP DDb, 
	MATRIX_OP tdo_scheme, 
	INT f_v, INT f_vv);
INT calc_TDA(FILE *TEX_fp, INT nb_geo, INT f_row_perms, 
	INT f_has_aut, LABRA_OP aut, 
	INT *nb_row_blocks_TDO, INT *nb_col_blocks_TDO, 
	INT *nb_row_blocks_TDA, INT *nb_col_blocks_TDA, 
	INT f_v, INT f_vv, INT f_incidences);
INT calc_blocks(VECTOR_OP B, PERMUTATION_OP q, PERMUTATION_OP qv, INT f_v);
INT action_on_blocks(VECTOR_OP gen_row_perms, 
	VECTOR_OP gen_col_perms, INT f_v, INT f_vv);
void print_incidences();
INT configuration_graph(MATRIX_OP G, MATRIX_OP I_and_G, MATRIX_OP G_and_I);
INT find_ij_line(INT i, INT j);
INT output_GEO(FILE *fp);
INT output_inc_tex(FILE *fp);
INT print_char_TD(FILE *fp, INT f_tex, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp);
INT print_char_and_inv_TD(FILE *fp, INT f_tex, 
	VECTOR_OP row_decomp_char, VECTOR_OP col_decomp_char, 
	VECTOR_OP row_decomp_inv, VECTOR_OP col_decomp_inv);
INT incma_latex(FILE *fp, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_labelled, 
	VECTOR_OP labelling_points, VECTOR_OP labelling_blocks, INT offset);
INT calc_blocks(VECTOR_OP B, 
	PERMUTATION_OP p, 
	PERMUTATION_OP q, PERMUTATION_OP qv, INT f_v);
INT TDA(VECTOR_OP gen, 
	VECTOR_OP row_decomp_char, VECTOR_OP col_decomp_char, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp_inv, VECTOR_OP col_decomp_inv, 
	INT f_v, INT f_vv);
INT calc_decomposition_matrix(
	VECTOR_OP row_decomp, VECTOR_OP col_decomp,
	MATRIX_OP S);
INT test_2_design(INT *k, INT *lambda);
INT nb_incidences();
INT levi(INT f_blocked, MATRIX_OP Levi);
INT self_dual(INT f_self_polar, 
	INT *f_is_sd, INT *f_is_sp, VECTOR_OP aut_gen, PERMUTATION_OP duality, 
	INT f_v, INT f_vv);
INT get_base_block_if_cyclic(PERMUTATION_OP per, VECTOR_OP bb);
};
INT mtx_test();
INT print_tdo_scheme(MATRIX_OP tdo_scheme, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, INT f_tex, FILE *fp);

#endif /* MA_INCLUDED */

