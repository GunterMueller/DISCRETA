/* gfq.h 
 * 
 * Anton Betten 
 * Jul 30, 1995 
 *
 * header file for 
 * gfq_zech.C
 * gfq_psu.C
 * gfq_sz.C
 * perm_rep.C
 * gfq_nb.C
 * singer.C
 */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#ifndef GFQ_INCLUDED
#define GFQ_INCLUDED

// singer.C:

typedef char G_INT_TYPE;

typedef G_INT_TYPE *G_POLYNOM;

#define g_deg(p) (p[-1])
#define g_set_deg(p, i) (p[-1] = i)

#undef TEST_DEG 
#undef DEBUG_DIV_REM
#undef DEBUG_MULT_MOD
#undef DEBUG_GCD

#undef DEBUG_BERLEKAMP
#undef DEBUG_IS_IRREDUCIBLE
#undef DEBUG_IS_SQUAREFREE
#undef DEBUG_SINGER



G_POLYNOM g_alloc(INT deg);
void g_free(G_POLYNOM p);
void g_zero(G_POLYNOM p);
void g_one(G_POLYNOM p);
void g_zero_coeffs(G_POLYNOM p);
void g_sprint(G_POLYNOM p, BYTE *str);
void g_sprint_latex(G_POLYNOM p, BYTE *str);
void g_print(G_POLYNOM p);
void g_print_latex(G_POLYNOM p);
void g_derive(G_POLYNOM p, G_POLYNOM q);
void g_calc_deg(G_POLYNOM p, G_INT_TYPE chi);
INT g_is_const(G_POLYNOM p, G_INT_TYPE chi);
INT g_is_zero(G_POLYNOM p, INT chi);
INT g_is_one(G_POLYNOM p, INT chi);
G_INT_TYPE g_as_const(G_POLYNOM p, G_INT_TYPE chi);
INT g_cmp(G_POLYNOM p, G_POLYNOM q);
G_POLYNOM g_copy(G_POLYNOM p);
void g_move(G_POLYNOM p, G_POLYNOM q);
INT g_add(G_POLYNOM a, G_POLYNOM b, G_POLYNOM c, G_INT_TYPE chi);
INT g_rank(G_INT_TYPE *M, INT m, INT n, G_INT_TYPE chi);
/* m \times n matrix M over GF(chi) */
G_INT_TYPE *g_berlekamp_matrix(G_POLYNOM m, G_INT_TYPE chi);
INT g_power_mod_apply(G_POLYNOM p, INT l, G_POLYNOM m, G_INT_TYPE chi);
INT g_russ_power_mod(G_POLYNOM a, 
	G_POLYNOM b, INT l, 
	G_POLYNOM m, G_INT_TYPE chi);
INT g_mult_mod(G_POLYNOM a, G_POLYNOM b, G_POLYNOM c, G_POLYNOM m, G_INT_TYPE chi);
/* LUENEBURG: 
 * On The Rational Normal Form Of Endomorphisms p. 36 */
INT g_mult(G_POLYNOM a, G_POLYNOM b, G_POLYNOM c, G_INT_TYPE chi);
INT g_mult_scalar(G_POLYNOM a, G_POLYNOM b, G_INT_TYPE s, INT chi);
INT g_mult_apply_scalar(G_POLYNOM a, G_INT_TYPE s, INT chi);
INT g_reduce(G_POLYNOM a, G_POLYNOM m, G_INT_TYPE chi);
INT g_div_rem(G_POLYNOM m, G_POLYNOM n, G_POLYNOM q, G_POLYNOM r, G_INT_TYPE chi);
INT g_is_squarefree(G_POLYNOM p, INT chi);
INT g_is_irreducible(G_POLYNOM p, INT chi);
INT g_is_primitive(G_POLYNOM p, G_INT_TYPE chi, INT m, VECTOR_OP vp);
INT g_gcd(G_POLYNOM m, G_POLYNOM n, G_POLYNOM g, INT chi);
INT g_bezout(
	G_POLYNOM m, G_POLYNOM n, 
	G_POLYNOM u, G_POLYNOM v, 
	G_POLYNOM g, INT chi);
/* g := gcd(m, n) = u * m + v * n */
void g_numeric_pol(INT i, G_POLYNOM p, G_INT_TYPE chi);
INT g_pol_numeric(G_POLYNOM p, G_INT_TYPE chi);
G_POLYNOM g_singer(INT deg, G_INT_TYPE chi);
INT g_singer_test(G_POLYNOM p, G_INT_TYPE chi, INT m, VECTOR_OP vp);
void g_singer_kandidat(G_POLYNOM p, INT k, G_INT_TYPE b, G_INT_TYPE a);
/* generates the polynomial 
 * X^k + X^{k-1} + b X + a
 */
G_INT_TYPE g_inv_mod(G_INT_TYPE a, G_INT_TYPE chi);
G_INT_TYPE g_asr(G_INT_TYPE a, G_INT_TYPE chi);

/* already declared in sym.h: */
/* int singer(long chi, long deg, long *pp); */

INT g_pol_r(G_POLYNOM a, G_POLYNOM b, G_POLYNOM r, G_INT_TYPE chi);
INT mtx_frobenius(G_INT_TYPE **M, G_POLYNOM m, G_INT_TYPE chi, INT f_v);






/* gfq_zech.C: */
#undef DEBUG_Z_INVERSE_NUM

/* gfq_psu.C: */
#undef DEBUG_MTX3_PERM
#undef DEBUG_MTX3_PERM_VERBOSE
#undef DEBUG_Z_MTX3_GEN_HK_NUM
#undef DEBUG_Z_MTX3_GEN_QAB_NUM

/* gfq_sz.C: */
#undef DEBUG_Z_MTX4_GEN_MLAMBDA_NUM
#undef DEBUG_Z_MTX4_GEN_SAB_NUM
#undef DEBUG_Z_MTX4_GEN_T_NUM
#undef DEBUG_MTX4_PERM_VERBOSE
#undef DEBUG_MTX4_PERM

typedef struct zech_data ZECH_DATA;
typedef struct isotropic_rep ISOTROPIC_REP;
typedef struct tits_ovoid_rep TITS_OVOID_REP;
typedef struct perm_rep PERM_REP;
typedef struct galois GALOIS;

#if TEXDOCU
#else
this header file is for the definitions and declarations 
of the following files:
//CODE
gfq_zech.C
gfq_psu.C
gfq_sz.C
gfq_nb.C
perm_rep.C
///CODE
gfq\_zech provides ``small'' finite fields via ZECH logarithms, 
gfq\_psu and gfq\_sz give generators for PSU and Sz groups, 
gfq\_nb provides normal bases and irreducible polynomials 
for GF(q), perm\_rep gives permutation representations 
for projective spaces over GF(q). 
There is also the file singer.C which is for 
producing irreducible primitve polynomials.
#endif

#if TEXDOCU
struct zech_data {
	INT chi, deg, q, m_order;
	INT idx_zero, idx_m1;
	G_POLYNOM m;
	G_POLYNOM alpha;
	G_POLYNOM *V;
	INT *Z;
	INT *Num;
	INT *Num_inv;
	INT *Frob;
	INT *Frob_inv;
};
#else
This is the data structure for finite fields $GF(q)$ where $q=p$ and $f > 1$ 
(e.g. true extension fields). The polynomial $m$ is computed 
by the singer program (see singer.C). 
It is a primitive irreducible polynomial of degree $f$ over $GF(p)$, 
e.g. its root $X$ is a primitive element for $GF(q)^*$. 
This means that the powers $X^i$, $i = 0, 1, \ldots q - 2$ 
run through all non-zero field elements. 
chi stands for $p$, deg stands for $f$ in the previous notation. 
$m_order$ stands for $q-1$. idx\_zero is always $q-1$. 
idx\_m1 is the index of -1 as a ZECH-log value.

Basically, we have two ways of enumerating the $GF(q)$-elements.
Both use the numbers $\{0, 1, \ldots q-1\}$ to stand for 
finite field elements.

The first one is  ``numeric'', i.e. we call the element 
$a_{f-1}X^{f-1} + \ldots + a_1 X + a_0$ by its corresponding 
numeric value $a_{f-1} p^{f-1} + \ldots + a_1 p + a_0$. 
In this representation in particular 0 stands for $0 \in GF(q)$, 
$1$ stands for $1 \in GF(q)$, $p$ stands for the primitive element $X$. 

The other one is the zech logarithm. 
Here, the elements of $GF(q)^*$ and $0$ are treated differently. 
The $GF(q)^*$ element $X^i$ ($i=0,1,\ldots, q-2$) is represented 
by the number $i$. The element $0$ is represented by the number $q-1$. 

The arithmetical functions indicated in their name which 
numbering of elements is used. 
For example, the function z\_mult\_zlog(a,b) 
multiplies two Zech-log numbers (and returns the result as a 
ZECH log number), the function z\_mult\_num(a,b) 
multiplies numbers in their numeric representation (and 
returns the result in the same representation). 
#endif

#if TEXDOCU
struct isotropic_rep {
	ZECH_DATA *Z;
	INT q1;
	INT *IV;
	INT nb_iv;
		/* number of isotropic vectors */
};
#endif

#if TEXDOCU
struct tits_ovoid_rep {
	ZECH_DATA *Z;
	INT *O;
	INT nb_o;
	INT f_has_omega;
	INT *omega_idx;
	INT *omega;
		/* omega is a omega_h \times omega_l matrix */
		/* the rows are the orbits of x \mapsto x^2 
		 * on the elements of O 
		 * omega_len[i] = length of i-th orbit = 
		 * number of valid entries in i-th row 
		 * of omega. */
	INT *omega_len;
	INT omega_l;
	INT omega_h;
};
#endif

#if TEXDOCU
struct perm_rep {
	INT f_affine;
	INT deg, n;
	ZECH_DATA *Z;
	INT q, p, f, f_verbose;
	INT idx_zero;
		/* = 0;
		 * the zero element: 0 * X^deg - 1 + ... 0 X + 0 */
	INT idx_one;
		/* = Z->Num\_inv[0];
		 * the element 1 = X^0 */
		/* idx_zero and idx_one are also set 
		 * in the GFp case (to 0 and 1). */
	INT *vec1, *vec2;
};
#endif

#if TEXDOCU
struct galois {
	INT chi;
	INT deg;
	G_POLYNOM m;
	G_INT_TYPE *Frob;
};
#endif

/* gfq_zech.C: */
ZECH_DATA *zech_open(INT deg, INT chi, INT f_verbose);
void zech_free(ZECH_DATA *Z);
void z_print_elt_num_verbose(ZECH_DATA *Z, INT i);

INT z_mult_zlog(ZECH_DATA *Z, INT i1, INT i2);
INT z_mult_num(ZECH_DATA *Z, INT i1, INT i2);
INT z_inverse_zlog(ZECH_DATA *Z, INT i);
INT z_inverse_num(ZECH_DATA *Z, INT i);
INT z_negate_zlog(ZECH_DATA *Z, INT i);
INT z_negate_num(ZECH_DATA *Z, INT i);
INT z_is_zero_zlog(ZECH_DATA *Z, INT i);
INT z_is_zero_num(ZECH_DATA *Z, INT i);
INT z_zero_zlog(ZECH_DATA *Z);
INT z_zero_num(ZECH_DATA *Z);
INT z_is_one_zlog(ZECH_DATA *Z, INT i);
INT z_is_one_num(ZECH_DATA *Z, INT i);
INT z_one_zlog(ZECH_DATA *Z);
INT z_one_num(ZECH_DATA *Z);
INT z_mone_zlog(ZECH_DATA *Z);
INT z_mone_num(ZECH_DATA *Z);
INT z_prim_elt_zlog(ZECH_DATA *Z);
INT z_prim_elt_num(ZECH_DATA *Z);
INT z_add_zlog(ZECH_DATA *Z, INT i1, INT i2);
INT z_add_num(ZECH_DATA *Z, INT i1, INT i2);
INT z_apply_frob_zlog(ZECH_DATA *Z, INT i1);
INT z_apply_frob_num(ZECH_DATA *Z, INT i1);
INT z_apply_frob2_zlog(ZECH_DATA *Z, INT i1);
INT z_apply_frob2_num(ZECH_DATA *Z, INT i1);
INT z_apply_frob_sz_zlog(ZECH_DATA *Z, INT i1);
INT z_apply_frob_sz_num(ZECH_DATA *Z, INT i1);
INT z_apply_frob_sz0_zlog(ZECH_DATA *Z, INT i1);
INT z_apply_frob_sz0_num(ZECH_DATA *Z, INT i1);
INT z_trace_zlog(ZECH_DATA *Z, INT a);
INT z_trace_num(ZECH_DATA *Z, INT a);
INT z_trace2_zlog(ZECH_DATA *Z, INT a);
INT z_trace2_num(ZECH_DATA *Z, INT a);
INT z_norm_zlog(ZECH_DATA *Z, INT a);
INT z_norm_num(ZECH_DATA *Z, INT a);
INT z_norm2_zlog(ZECH_DATA *Z, INT a);
INT z_norm2_num(ZECH_DATA *Z, INT a);
void z_print_perm(INT *p, INT len);
INT z_sprint_perm(INT *p, INT length, BYTE *str);
INT z_search(ZECH_DATA *Z, G_POLYNOM p);

void z_mtx_i_print(INT *mtx, INT dim);
void z_mtx_print_nm(INT *mtx, INT n, INT m);
void z_vec_i_print(INT *vec, INT dim);
void z_mtx3_print(ZECH_DATA *Z, INT *mtx);
void z_mtx4_print(ZECH_DATA *Z, INT *mtx);
void z_mtx_i_gen_diag_num(ZECH_DATA *Z, INT *mtx, INT *x, INT dim);
void z_mtx3_gen_diag_num(ZECH_DATA *Z, INT *mtx, INT *x);
void z_mtx4_gen_diag_num(ZECH_DATA *Z, INT *mtx, INT *x);

INT z_mtx_vec_mult_GFq_num(ZECH_DATA *Z, INT *mtx, INT n, 
	INT *vec1, INT *vec2);
INT z_mtx_vec_mult_GFp(INT p, INT *mtx, INT n, 
	INT *vec1, INT *vec2);
INT z_vec_norm_projective_GFq_num(ZECH_DATA *Z, INT *vec, INT n);
INT z_vec_norm_projective_GFp(INT p, INT *vec, INT n);
INT z_mtx_affine_rep_Zq_num(ZECH_DATA *Z, INT *mtx, INT n, 
	PERMUTATION_OP perm);
INT z_mtx_projective_rep_Zq_num(ZECH_DATA *Z, INT *mtx, INT n, 
	PERMUTATION_OP perm);
INT z_mtx_rep_Zq_num(INT f_affine, ZECH_DATA *Z, INT *mtx, INT n, 
	PERMUTATION_OP perm);
INT z_mtx_affine_rep_Zp(INT p, INT *mtx, INT n, 
	PERMUTATION_OP perm);
INT z_mtx_projective_rep_Zp(INT p, INT *mtx, INT n, 
	PERMUTATION_OP perm);
INT z_mtx_rep_Zp(INT f_affine, INT p, INT *mtx, INT n, 
	PERMUTATION_OP perm);
INT z_Frobenius_rep_Zq_num(INT f_affine, ZECH_DATA *Z, INT n, 
	PERMUTATION_OP perm);
INT z_translation_rep_Zq_num(ZECH_DATA *Z, INT n, INT ei, INT betaj, 
	PERMUTATION_OP perm);
INT z_translation_rep_Zp(INT p, INT n, INT ei, 
	PERMUTATION_OP perm);
INT z_GL_n_q_gen_num(ZECH_DATA *Z, INT ***mtx, INT *len, INT n);
INT z_GL_n_p_gen(INT ***mtx, INT *len, INT n, INT p);
INT z_GL_gen(INT ***mtx, INT *len, INT n, 
	INT idx_zero, INT idx_one, INT idx_alpha);
INT z_test_GL_perm_rep(INT f_affine, INT f_semi, INT f, INT p, INT n);
INT z_GL_n_q_perm_rep_info(INT f_affine, 
	INT f_semi, INT f_special, INT f, INT p, INT n, 
	INT *deg, BYTE *label, BYTE *label_tex);
INT z_GL_n_q_perm_rep(INT f_affine, INT f_semi, INT f, INT p, INT n, 
	VECTOR_OP V, INT f_verbose);
INT z_mtx_determinante_GFp(INT p, INT *mtx, INT n);
INT z_mtx_determinante_GFq_num(ZECH_DATA *Z, INT *mtx, INT n);
INT z_mtx_inverse_GFp(INT p, INT *mtx, INT *mtx_inv, INT n);
INT z_mtx_inverse_GFq_num(ZECH_DATA *Z, INT *mtx, INT *mtx_inv, INT n);
INT z_gauss_n_m_GFp(INT p, INT f_special, INT *mtx, INT n, INT m);
INT z_gauss_n_m_GFq_num(ZECH_DATA *Z, INT f_special, INT *mtx, INT n, INT m);
INT z_gauss_n_m_rectangular_GFp(INT p, INT f_special, INT *mtx, INT n, INT m, INT dim_m, 
	INT *base_cols);
INT z_gauss_n_m_rectangular_GFq_num(ZECH_DATA *Z, INT f_special, 
	INT *mtx, INT n, INT m, INT dim_m, INT *base_cols);

/* gfq_psu.C: */
void z_mtx3_gen_Qab_num(ZECH_DATA *Z, INT *mtx, INT a, INT b);
void z_mtx3_gen_Hk_num(ZECH_DATA *Z, INT *mtx, INT k);
void z_mtx3_gen_T_num(ZECH_DATA *Z, INT *mtx);
ISOTROPIC_REP *ir_open(INT deg, INT chi);
void ir_close(ISOTROPIC_REP *ir);
INT ir_search(ISOTROPIC_REP *ir, INT *x);
INT z_collect_tr_plus_nr_zero(ZECH_DATA *Z, 
	INT *T, INT *nb_sol);
INT ir_vec3_check_if_normalized_num(ISOTROPIC_REP *ir, INT *v);
INT ir_vec3_is_normalized_num(ISOTROPIC_REP *ir, INT *v);
void ir_vec3_normalize_num(ISOTROPIC_REP *ir, INT *v);
void ir_mtx3_normalize_num(ISOTROPIC_REP *ir, INT *v);
INT ir_cmp_vec(INT *v1, INT *v2);
INT ir_mtx3_perm(ISOTROPIC_REP *ir, 
	INT *mtx, PERMUTATION_OP perm);
INT ir_generator_T(ISOTROPIC_REP *ir, 
	PERMUTATION_OP T_perm, INT f_verbose);
INT ir_generator_Hk(ISOTROPIC_REP *ir, 
	PERMUTATION_OP H_perm, INT f_verbose);
INT ir_generators_Q(ISOTROPIC_REP *ir, 
	VECTOR_OP V, INT f_verbose);
INT PSU_3_q2_generators(INT q, VECTOR_OP V);
INT gen_PSU(INT q);

/* gfq_sz.C: */
void z_mtx4_gen_Mlambda_num(ZECH_DATA *Z, INT *mtx, INT lambda);
void z_mtx4_gen_Sab_num(ZECH_DATA *Z, INT *mtx, INT a, INT b);
void z_mtx4_gen_T_num(ZECH_DATA *Z, INT *mtx);
TITS_OVOID_REP *tr_open(INT deg, INT chi);
void tr_close(TITS_OVOID_REP *tr);
INT tr_gen_tits_ovoid(TITS_OVOID_REP *tr);
INT tr_gen_omega(TITS_OVOID_REP *tr);
INT tr_cmp_vec(INT *v1, INT *v2);
INT tr_search(TITS_OVOID_REP *tr, INT *x);
void tr_vec4_normalize_num(TITS_OVOID_REP *tr, INT *v);
INT tr_mtx4_perm(TITS_OVOID_REP *tr, 
	INT *mtx, PERMUTATION_OP perm);
INT tr_generator_T(TITS_OVOID_REP *tr, 
	PERMUTATION_OP T_perm, INT f_verbose);
INT tr_generator_Mlambda(TITS_OVOID_REP *tr, 
	PERMUTATION_OP M_perm, INT f_verbose);
INT tr_generators_S(TITS_OVOID_REP *tr, 
	VECTOR_OP V, INT f_verbose);
INT Sz_q_generators(INT q, VECTOR_OP V);
INT gen_Sz(INT q);

/* perm_rep.C: */
INT perm_rep_GL_n_q_degree(INT f_affine, 
	INT f_semi, INT f_special, INT q, INT n, 
	INT f_verbose);
INT perm_rep_T_n_q_generators(INT q, INT n, VECTOR_OP V, INT f_verbose);
INT perm_rep_GL_n_q_generators(INT f_affine, 
	INT f_semi, INT f_special, INT f_translations, INT q, INT n, 
	VECTOR_OP V, INT f_verbose);
PERM_REP *open_perm_rep(INT q, 
	INT f_affine, INT n, INT f_verbose);
void perm_rep_free(PERM_REP *P);
INT perm_rep_i2vec(PERM_REP *P, INT i, INT *vec);
INT perm_rep_vec2i(PERM_REP *P, INT *vec, INT *i);
INT perm_rep_vec2i_with_normalization(PERM_REP *P, INT *vec, INT *i);
INT z_test_affine_rep(INT q, INT n);
INT z_affine_degree(INT q, INT n);
INT z_affine_i2vec(INT q, INT n, INT i, INT *vec);
INT z_affine_vec2i(INT q, INT n, INT *vec, INT *i);
INT z_test_projective_rep(INT f, INT p, INT n);
INT z_projective_degree(INT q, INT n);
INT z_projective_i2vec(INT q, INT idx_zero, INT idx_one, 
	INT n, INT i, INT *vec);
INT z_projective_vec2i(INT q, INT idx_zero, INT idx_one, 
	INT n, INT *vec, INT *i);
INT perm_rep_index_of_ei(PERM_REP *P, INT i);
INT perm_rep_mtx(PERM_REP *P, INT *mtx, PERMUTATION_OP perm);
INT perm_rep_to_mtx(PERM_REP *P, PERMUTATION_OP perm, INT *mtx);
INT perm_rep_determinante(PERM_REP *P, PERMUTATION_OP perm);
INT perm_rep_determinante_v(PERM_REP *P, PERMUTATION_OP perm, INT f_verbose);
INT perm_rep_translation_generators(PERM_REP *P, VECTOR_OP V, INT f_verbose);
INT perm_rep_GL_generators(PERM_REP *P, 
	INT f_semi, INT f_special, INT f_translations, 
	VECTOR_OP V, INT f_verbose);
INT perm_rep_GL2SL(PERM_REP *P, 
	VECTOR_OP GL_gen, VECTOR_OP SL_gen, INT f_verbose);
INT perm_rep_print_vectors(PERM_REP *P);
INT perm_rep_mult(PERM_REP *P, INT x, INT y);
INT perm_rep_inverse(PERM_REP *P, INT x);
INT perm_rep_negate(PERM_REP *P, INT x);
INT perm_rep_is_zero(PERM_REP *P, INT x);
INT perm_rep_zero(PERM_REP *P);
INT perm_rep_is_one(PERM_REP *P, INT x);
INT perm_rep_one(PERM_REP *P);
INT perm_rep_mone(PERM_REP *P);
INT perm_rep_add(PERM_REP *P, INT x, INT y);
INT perm_rep_power_int(PERM_REP *P, INT x, INT l);
INT perm_rep_gauss_n_m_rectangular(PERM_REP *P, 
	INT f_special, INT f_complete, MATRIX_OP M, VECTOR_OP base_cols, INT f_v);
INT perm_rep_get_kernel(PERM_REP *P, MATRIX_OP M, MATRIX_OP K, VECTOR_OP base_cols);

/* gfq_nb.C: */
INT gfq_GL_classes_bi(VECTOR_OP R, INT chi, INT k, INT f_v, INT f_vv);
INT gfq_bi2mtx(MATRIX_OP M, VECTOR_OP bi, INT q, INT k);
INT gfq_mtx2bi(MATRIX_OP M, VECTOR_OP bi, INT q);
INT gfq_gl_classes(VECTOR_OP V, INT q, INT k);
INT gfq_gl_class_centralizer_order_kung(INT q, 
	VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, SYM_OP co);
INT gfq_gl_class2matrix(INT q, VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, INT k, MATRIX_OP M);
INT gfq_gl_class_first(VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, INT k);
INT gfq_gl_class_next(VECTOR_OP V, VECTOR_OP V_mult, VECTOR_OP V_part, INT k);
INT gfq_pol_part_first(VECTOR_OP V_mult, VECTOR_OP V_part);
INT gfq_pol_part_next(VECTOR_OP V_mult, VECTOR_OP V_part);
INT gfq_choose_pol(VECTOR_OP V, INT k);
INT gfq_choose_pol_first(VECTOR_OP V, VECTOR_OP V_mult, INT k);
INT gfq_choose_pol_next(VECTOR_OP V, VECTOR_OP V_mult);
INT gfq_irred_polynomials(VECTOR_OP V, INT chi, INT deg_min, INT deg_max, INT f_v, INT f_vv);
INT gfq_irred_pol(VECTOR_OP V, INT chi, INT deg, INT f_v, INT f_vv);
void print_mtx(INT *M, INT dim_M, INT deg);
void test_gfq_nb(INT chi, INT deg, INT f_v);
INT gfq_minpol(GALOIS *GFq, G_POLYNOM g, G_POLYNOM mue);
INT gfq_dep(GALOIS *GFq, INT *v, INT *a, INT dim_a, INT m, INT *rho);
INT gfq_ord_ideal(GALOIS *GFq, INT i, G_POLYNOM mue, INT *nb, INT dim_nb);
INT gfq_generator(GALOIS *GFq, INT *nb, INT dim_nb, INT f_v, INT f_vv);
INT gfq_calc_nb(GALOIS *GFq, G_POLYNOM p, INT *nb, INT dim_nb, INT f_v);
INT gfq_frob_apply(GALOIS *GFq, G_POLYNOM p);
INT gfq_is_regular_word(INT *v, INT n, INT q);
INT gfq_first_irred_pol(INT q, INT p, INT deg, INT *v);
INT gfq_next_irred_pol(INT q, INT p, INT deg, INT *v);

#endif /* GFQ_INCLUDED */


