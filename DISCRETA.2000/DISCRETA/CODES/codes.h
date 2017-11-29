/* bch.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef BCH_INCLUDED
#define BCH_INCLUDED

#ifndef GFQ_INCLUDED
#include <DISCRETA/gfq.h>
#endif
#ifndef UNIP_INCLUDED
#include <DISCRETA/unip.h>
#endif
#ifndef MA_INCLUDED
#include <DISCRETA/ma.h>
#endif
#ifndef DIVS_INCLUDED
#include <DISCRETA/divs.h>
#endif

#define CODE_INDEX_ID 0
#define CODE_INDEX_N 1
#define CODE_INDEX_K 2
#define CODE_INDEX_D 4
#define CODE_INDEX_AGO 6

class code_ob : public VECTOR_OB {
public:
	INTEGER_OP s_id() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_id_i() { 
		return(s_id()->s_i()); };
	INTEGER_OP s_n() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_n_i() { 
		return(s_n()->s_i()); };
	INTEGER_OP s_k() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_k_i() { 
		return(s_k()->s_i()); };
	INTEGER_OP s_q() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_q_i() { 
		return(s_q()->s_i()); };
	INTEGER_OP s_d() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_d_i() { 
		return(s_d()->s_i()); };
	SYM_OP s_ago() { 
		return((INTEGER_OP)s_i(5)); };
	STRING_OP s_ago_text() { 
		return((STRING_OP) s_i(6)); };
	BYTE *s_ago_text_s() { 
		return(s_ago_text()->s_str()); };
	MATRIX_OP s_M() { 
		return((MATRIX_OP)s_i(7)); };
	INTEGER_OP s_M_ij(INT i, INT j) { 
		return((INTEGER_OP)s_M()->s_ij(i, j)); };
	INT s_M_iji(INT i, INT j) { 
		return(s_M_ij(i, j)->s_i()); };
	VECTOR_OP s_W() { 
		return((VECTOR_OP)s_i(8)); };
	INTEGER_OP s_W_i(INT i) { 
		return((INTEGER_OP)s_W()->s_i(i)); };
	INT s_W_ii(INT i) { 
		return(s_W_i(i)->s_i()); };
	VECTOR_OP s_T() { 
		return((VECTOR_OP)s_i(9)); };
	STRING_OP s_T_i(INT i) { 
		return((STRING_OP)s_T()->s_i(i)); };
	BYTE *s_T_is(INT i) { 
		return(s_T_i(i)->s_str()); };

INT field_name(INT i, INT j, BYTE *str);
INT Init(INT n, INT k, INT q);
INT Init_M(MATRIX_OP M);
INT sprint(BYTE *s);
INT Print();
INT fprint(FILE *fp);
INT set_ago(VECTOR_OP aut_gen);
INT init_map(INT n, INT k, INT q, INT *map, PERM_REP *P);
INT mindist();
};

// bch.C:
INT unip_mult_gfq(ZECH_DATA *Z, UNIPOLY_OP a, UNIPOLY_OP b, UNIPOLY_OP c);
INT unip_init_Xma(ZECH_DATA *Z, INT a, UNIPOLY_OP p);
INT unip_minpoly(ZECH_DATA *Z, INT a, UNIPOLY_OP m, INT f_v);
INT unip_QR_code_generator_polynomial(UNIPOLY_OP g, INT n, INT p, INT f_v, INT f_vv);
INT unip_BCH_generator_polynomial(UNIPOLY_OP g, INT n, INT designed_distance, 
	INT *bose_distance, INT p, INT f_v, INT f_vv);
INT generator_mat_cyclic_code(UNIPOLY_OP g, INT n, MATRIX_OP G);

INT i_db_codes(DATABASE_OP *p_db, BYTE *path);
void e_db_codes(DATABASE_OP db);
INT db_codes_create_and_close(BYTE *path);
INT db_codes_info(BYTE *path, INT btree_idx);
INT db_codes_dump(BYTE *path, INT btree_idx, BYTE *fname);
INT db_codes_highest_id(DATABASE_OP db, INT *highest_id);
INT db_codes_add(DATABASE_OP db, CODE_OP C);

// mindist.C
INT code_Mindist_of_matrix(MATRIX_OP G, INT q);
INT code_Mindist_of_map(INT n, INT k, INT q, VECTOR_OP gamma);
INT code_mindist_of_map(INT n, INT k, INT q, INT *gamma, PERM_REP *P);
int mindist(int n, int k, int q, int *d, int *G, 
	int f_v, int f_vv, int idx_zero, int idx_one, ZECH_DATA *Z);

#endif /* BCH_INCLUDED */

