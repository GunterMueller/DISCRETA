/* gtP.h 
 * private include file for gt_canon.C
 */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#ifndef CP_INCLUDED
#include <cp.h>
#endif

typedef struct ordered_set ORDERED_SET;
typedef struct gt_canon_local GT_CANON_LOCAL;
typedef struct gt_canon_grid GT_CANON_GRID;

#define GT_CANON_MAX_N 512
#define GT_CANON_MAX_SET_SIZE GT_CANON_MAX_N
#define GT_CANON_MAX_TYPE 512
#define GT_CANON_MAX_GRID 512

struct ordered_set {
	INT a[GT_CANON_MAX_SET_SIZE];
	INT size;
};

struct gt_canon_local {
	GT_CANON_INFO *info;
	ORDERED_SET *E;
		/* allocated in autlog_int(), 
		 * freed in autlog_ext(). */
	
	INT first_moved;
		/* only used if mode 
		 * is ONLY_COSET_REPS */
	INT f_going_back;
	
	/* for canonicize: */
	INT is_first;
		/* first_moved used also */
	SPERM p0;
	SPERM p0v;
	INT nb_gen;
	INT g0[GT_CANON_MAX_GEN];
	INT go0[GT_CANON_MAX_GEN];
	INT nb_aut_gens;
	INT dim_aut_gens;
	SPERM *aut_gens;
	INT *auts_first_moved;

	INT n;
	INT *theG;
		/* n x (n + 1); (n + 1) - th 
		 * column are inverse elements;
		 * dimension n x dim_n */
	INT dim_n; /* >= n + 1 */
	/* INT nb_gen; */
	INT g[GT_CANON_MAX_GEN];
		/* 0 .. nb_gen-1: the generatores; */
		/* group A: a copy from AUTLOG_INFO 
		 * group B: searched via backtracking. */
	INT go[GT_CANON_MAX_GEN];
		/* group order:
		 * go[i] = |<g_0, ... , g_i>|, i.e.
		 * go[0] = |<g_0>|, 
		 * go[nb_gen - 1] = full group order. */
	SPERM p; /* degree n */
	SPERM pv; /* p^-1 */
	SPERM q; /* degree n */
	SPERM qv; /* q^-1 */
	SPERM r; /* degree n */
	SPERM rv; /* r^-1 */
	INT nb_grid;
	GT_CANON_GRID *Grid;
	INT f_has_element_colors;
	INT *element_colors;
};

struct gt_canon_grid {
	SPERM r; /* degree m = n of autlog_the_group */
	SPERM rv; /* r^-1 */
	INT m; /* = n of autlog_the_group */
	INT n; /* = type_len */
	INT type[GT_CANON_MAX_N][GT_CANON_MAX_TYPE];
	INT G_max;
	INT first[GT_CANON_MAX_GRID];
	INT len[GT_CANON_MAX_GRID];
	INT type_idx[GT_CANON_MAX_GRID];
	INT grid_entry[GT_CANON_MAX_GRID];
};

static INT gt_canonicize_main(GT_CANON_INFO *info, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut);
static INT gt_canonicize(GT_CANON_LOCAL *L, INT level);
static INT recalc_Ro_Re_GG(GT_CANON_LOCAL *L, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG);
static INT get_aut_gen(GT_CANON_LOCAL *L);
static void free_aut_gen(GT_CANON_LOCAL *L);
static INT get_canonical_labelling(GT_CANON_LOCAL *L, INT level);
static void print_labelling_g(GT_CANON_LOCAL *L, INT level);
static void print_labelling_go(GT_CANON_LOCAL *L, INT level);
static INT canonicize_compare_tables(GT_CANON_LOCAL *L, INT level);
static void print_E_colors(GT_CANON_LOCAL *L, INT level);
static void print_element_colors(GT_CANON_LOCAL *L, INT level);
static INT get_aut_gens(GT_CANON_LOCAL *L, 
	VECTOR_OP auts, VECTOR_OP first_moved);
static INT get_aut_group(GT_CANON_LOCAL *L, 
	LABRA_OP Aut, SYM_OP ago);
static INT info_fill_in_theG(GT_CANON_INFO *info, GROUP_TABLE *G);
static INT info_fill_in_colors(FILE *fp, GT_CANON_INFO *info, 
	GROUP_TABLE *G, INT f_v);
static void the_group_nil(GT_CANON_LOCAL *L);
static void the_group_int(GT_CANON_LOCAL *L);
static INT the_group_init_colors(GT_CANON_LOCAL *L, VECTOR_OP col);
static INT the_group_open_grid(GT_CANON_LOCAL *L, INT nb_grid, INT n);
static void the_group_ext(GT_CANON_LOCAL *L);
static INT the_group_print_generators(GT_CANON_LOCAL *L);
static INT the_group_order(GT_CANON_LOCAL *L, INT k);
/* k = 0: 1 fuer Einsgruppe, sonst go[k - 1] */
static void grid_nil(GT_CANON_GRID *p);
static void grid_int(GT_CANON_GRID *p, INT n);
static void grid_ext(GT_CANON_GRID *p);
static void gt_canon_nil(GT_CANON_LOCAL *a);
static void group_print(GT_CANON_LOCAL *G, SPERM *p, SPERM *pv);

/* gt_canon_grid.C: */
static INT calc_grid(GT_CANON_LOCAL *L, 
	GT_CANON_GRID *G, INT k, INT go, INT f_verbose);
static INT radix_sort(GT_CANON_GRID *G, INT radix, INT first, INT last);
static INT insert_idx(GT_CANON_GRID *G, 
	INT first, INT len, INT radix, 
	INT search_this, INT *idx);
static void print_grid(GT_CANON_GRID *G, SPERM *pv, 
	INT k, INT go_k, INT go);
static void print_type_vector(GT_CANON_GRID *G, INT i, INT k);
static INT calc_n_grids(GT_CANON_LOCAL *L, 
	INT go, INT n, INT f_verbose);
static INT choose_first_into_E(GT_CANON_LOCAL *L, 
	GT_CANON_GRID *Grid, INT k, ORDERED_SET *E, INT f_verbose);
static void print_E(ORDERED_SET *E, INT k);

/* gt_canon_dimino.C: */
static INT gt_add_generator(GT_CANON_LOCAL *B, 
	INT i, INT g, INT f_verbose);
static INT gt_dimino(GT_CANON_LOCAL *C, 
	INT k, INT g0, INT *go);
static INT gt_add_coset(GT_CANON_LOCAL *C, 
	INT coset_size, INT g2, INT *go);
static INT gt_mult(GT_CANON_LOCAL *A, INT i, INT j);


