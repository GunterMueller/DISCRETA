/* autP.h 
 * private include file for aut.C
 */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

typedef struct ordered_set ORDERED_SET;
typedef struct autlog_the_group AUTLOG_THE_GROUP;
typedef struct autlog_grid AUTLOG_GRID;
typedef struct autlog_factorgroup AUTLOG_FACTORGROUP;
typedef struct autlog_commutator_subgroup 
	AUTLOG_COMMUTATOR_SUBGROUP;
typedef struct autlog_local AUTLOG_LOCAL;

#ifdef SYSTEMMAC
#define AUTLOG_MAX_N 50
#define AUTLOG_MAX_SET_SIZE AUTLOG_MAX_N
#define AUTLOG_MAX_TYPE 50
#define AUTLOG_MAX_GRID 50
#else
#define AUTLOG_MAX_N 512
#define AUTLOG_MAX_SET_SIZE AUTLOG_MAX_N
#define AUTLOG_MAX_TYPE 512
#define AUTLOG_MAX_GRID 512
#endif

struct ordered_set {
	INT a[AUTLOG_MAX_SET_SIZE];
	INT size;
};

struct autlog_the_group {
	INT n;
	INT *theG;
		/* n x (n + 1); (n + 1) - th 
		 * column are inverse elements;
		 * dimension n x dim_n */
	INT dim_n; /* >= n + 1 */
	INT nb_gen;
	INT g[AUTLOG_MAX_G];
		/* 0 .. nb_gen-1: the generatores; */
		/* group A: a copy from AUTLOG_INFO 
		 * group B: searched via backtracking. */
	INT go[AUTLOG_MAX_G];
		/* group order:
		 * go[i] = |<g_0, ... , g_i>|, i.e.
		 * go[0] = |<g_0>|, 
		 * go[nb_gen - 1] = full group order. */
	SPERM p; /* degree n */
	SPERM pv; /* p^-1 */
	SPERM q;
		/* degree n, temporary, 
		 * for aut_add_coset() */
	SPERM qv; /* q^-1 */
		/* q is not a column permutation 
		 * as in the geometrical 
		 * isomorphism program ! */
	SPERM tmp1; /* degree n */
	SPERM tmp2; /* degree n */
	INT nb_grid;
	AUTLOG_GRID *Grid;
	INT f_has_element_colors;
	INT *element_colors;
	FILE *fp_txt;
};

struct autlog_grid {
	SPERM r; /* degree m = n of autlog_the_group */
	SPERM rv; /* r^-1 */
	INT m; /* = n of autlog_the_group */
	INT n; /* = type_len */
	INT type[AUTLOG_MAX_N][AUTLOG_MAX_TYPE];
	INT G_max;
	INT first[AUTLOG_MAX_GRID];
	INT len[AUTLOG_MAX_GRID];
	INT type_idx[AUTLOG_MAX_GRID];
	INT grid_entry[AUTLOG_MAX_GRID];
};

struct autlog_factorgroup {

	/* computed in autlog_center(): */
	INT *H;
	INT H_len;
		/* H, H_len are allocated or set in 
		 * autlog_center(); H contains list of 
		 * central elements of the group
		 * (the center Z). */

	/* computed in autlog_do_factorgroup(): */
	INT *H_rep_idx;
		/* for any element g \in G:
		 * the number of the 
		 * coset in the factorgroup G mod Z
		 * that it lies in */
	INT *H_rep;
		/* for any coset: the representing 
		 * element (as an element number in G); */
	
	INT nb_gen;
	INT *g;
		/* the generators for G as coset numbers 
		 * in the factor group (via H_rep_idx) */
	INT *g_idx;
		/* g_idx[i] is index into the factorgroup 
		 * AUTLOG_THE_GROUP generator table, 
		 * where this generator lies. 
		 * Otherwise -1, if it is redundant. 
		 * allocated in autlog_group_gen(). */
};

struct autlog_commutator_subgroup {
	INT dim_n; /* = AmZ->n */
	INT *C_idx; /* indices into C[] */
	INT nb_C;
		/* the number of different, 
		 * ascendingly sorted 
		 * commutators in C[] */
	INT *C;
		/* the commutators, sorted, 
		 * as element numbers in G */
	INT nb_gen;
		/* = AmZ->nb_gen */
	INT *g_len;
/* g_len:
 * i = 0, ... nb_gen - 1:
 * number of commutators 
 * formed by elements from
 * <g_0, ... g_i> of AmZ; 
 * g_len[nb_gen - 1] == nb_C. */
	INT *g;
		/* g[0] = 0 = id */
	INT *g_idx;
		/* -1 if not used in extending the derived group; 
		 * otherwise its index into 
		 * the g[] table of Ad or Bd. 
		 * especially: g_idx[0] = -1 
		 * (id is not considered to be a generator 
		 * by autlog_group_gen). */
};

struct autlog_local {
	AUTLOG_INFO *info;
	ORDERED_SET *E;
		/* allocated in autlog_int(), 
		 * freed in autlog_ext(). */
	
	AUTLOG_THE_GROUP A, B;
	
	/* the center */
	AUTLOG_FACTORGROUP A_Z, B_Z;

	/* the factor groups A mod Z and B mod Z: */
	AUTLOG_THE_GROUP AmZ, BmZ;
	
	/* the commutators: */
	AUTLOG_COMMUTATOR_SUBGROUP AcA, BcB;

	AUTLOG_THE_GROUP Ad, Bd;
		/* commutator subgroup (derived group) */
	
	AUTLOG_THE_GROUP *pA, *pB;
		/* these are the groups which we
		 * would like to find isomorphisms for
		 * set in autlog_int_isomorphism()
		 * to A and B, 
		 * in autlog_int_autologism()
		 * to AmZ, BmZ */
	
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
	INT g[AUTLOG_MAX_G];
	INT go[AUTLOG_MAX_G];
	INT nb_aut_gens;
	INT dim_aut_gens;
	SPERM *aut_gens;
	INT *auts_first_moved;
};

/* aut.C: */
static void print_E_colors(AUTLOG_LOCAL *L, INT level);
static void print_element_colors(AUTLOG_LOCAL *L, INT level);
static INT autlog_canonicize(AUTLOG_INFO *info, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut);
static INT recalc_Ro_Re_GG(AUTLOG_LOCAL *L, VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG);
static INT get_aut_gens(AUTLOG_LOCAL *L, VECTOR_OP auts, VECTOR_OP first_moved);
static INT get_aut_group(AUTLOG_LOCAL *L, LABRA_OP Aut, SYM_OP ago);
static INT canonicize(AUTLOG_LOCAL *L, INT level);
static INT get_aut_gen(AUTLOG_LOCAL *L);
static void free_aut_gen(AUTLOG_LOCAL *L);
static INT get_canonical_labelling(AUTLOG_LOCAL *L, INT level);
static void print_labelling_g(AUTLOG_LOCAL *L, INT level);
static void print_labelling_go(AUTLOG_LOCAL *L, INT level);
static INT canonicize_compare_tables(AUTLOG_LOCAL *L, INT level);
static INT autlog_do_level(
	AUTLOG_LOCAL *L, INT i);
static void autlog_print_mapping(
	AUTLOG_LOCAL *L, INT i);

/* aut_util.C: */
static void print_size();
static INT al_test2(FG_OP G, INT f_verbose);
static void autlog_the_group_nil(AUTLOG_THE_GROUP *G);
static void autlog_the_group_int(AUTLOG_THE_GROUP *G);
static INT autlog_the_group_init_colors(AUTLOG_THE_GROUP *G, VECTOR_OP col);
static INT autlog_the_group_open_grid(AUTLOG_THE_GROUP *G, INT nb_grid, INT n);
static void autlog_the_group_ext(AUTLOG_THE_GROUP *G);
static INT autlog_the_group_init_generators(
	AUTLOG_THE_GROUP *G, 
	INT nb_gen, INT *gen);
static INT autlog_the_group_print_generators(
	AUTLOG_THE_GROUP *G);
static INT autlog_the_group_order(
	AUTLOG_THE_GROUP *G, INT k);
/* k = 0: 1 fuer Einsgruppe, sonst go[k - 1] */
static void alfg_nil(AUTLOG_FACTORGROUP *p);
static void alfg_ext(AUTLOG_FACTORGROUP *p);
static INT alfg_print_generators(
	AUTLOG_FACTORGROUP *G_Z);
static void alcs_nil(
	AUTLOG_COMMUTATOR_SUBGROUP *p);
static void alcs_ext(
	AUTLOG_COMMUTATOR_SUBGROUP *p);
static void algrid_nil(AUTLOG_GRID *p);
static void algrid_int(AUTLOG_GRID *p, INT n);
static void algrid_ext(AUTLOG_GRID *p);
static void autlog_nil(AUTLOG_LOCAL *a);
static void autlog_group_print(AUTLOG_THE_GROUP *G, 
	SPERM *p, SPERM *pv);

/* aut_init.C: */
static INT autlog_int(AUTLOG_LOCAL *a, 
	AUTLOG_INFO *info, INT f_B);
static INT autlog_ext(AUTLOG_LOCAL *a);
static INT autlog_int_isomorphism(
	AUTLOG_LOCAL *L, INT f_verbose);
static INT autlog_int_autologism(
	AUTLOG_LOCAL *L, INT f_B, 
	INT f_verbose, INT f_very_verbose);
static INT autlog_do_group_gen_from_factorgroup(
	AUTLOG_THE_GROUP *GmZ, 
	AUTLOG_FACTORGROUP *G_Z, 
	INT f_verbose);
static INT autlog_do_center(
	AUTLOG_THE_GROUP *G, 
	AUTLOG_FACTORGROUP *G_Z, 
	INT f_verbose);
static INT autlog_do_do_factorgroup(
	AUTLOG_THE_GROUP *G, 
	AUTLOG_FACTORGROUP *G_Z, 
	AUTLOG_THE_GROUP *GmZ, 
	INT f_verbose);
static INT autlog_do_derivedgroup(
	AUTLOG_THE_GROUP *G, 
	AUTLOG_FACTORGROUP *G_Z, 
	AUTLOG_THE_GROUP *GmZ, 
	AUTLOG_COMMUTATOR_SUBGROUP *GcG, 
	AUTLOG_THE_GROUP *Gd, 
	INT f_verbose);
static INT autlog_do_derivedgroup_generators(
	AUTLOG_THE_GROUP *GmZ, 
	AUTLOG_COMMUTATOR_SUBGROUP *GcG, 
	INT f_verbose);
static INT autlog_center(
	AUTLOG_THE_GROUP *G, INT **Z, INT *nb_Z, 
	INT f_verbose);
static INT autlog_do_factorgroup(
	AUTLOG_THE_GROUP *G, INT *H, INT H_len, 
	INT **H_rep_idx, INT **H_rep, 
	INT **GmH, INT **GmH_gen, 
	INT *GmH_n, INT *GmH_dim, 
	INT f_verbose);
static INT autlog_derivedgroup(AUTLOG_THE_GROUP *G, 
	INT *H_rep_idx, INT *H_rep, 
	INT *GmH, INT GmH_n, INT GmH_dim,
	INT **C_idx, INT **C, INT *nb_C, 
	INT f_verbose);
static INT autlog_derivedgroup_generators(
	AUTLOG_THE_GROUP *GmZ, 
	INT *C_idx, INT *C, INT nb_C, 
	INT **g_len, INT **g, INT f_verbose);
static INT aldg_gen(AUTLOG_THE_GROUP *GmZ, 
	INT i0, INT i1, INT j0, INT j1, 
	INT *C_idx, INT *C, 
	INT *g, INT *nb_gen, INT f_verbose);



/* aut_dimino.C: */
static INT autlog_group_gen(
	AUTLOG_THE_GROUP *C, 
	INT nb_gen, INT *g, 
	INT **g_idx, INT f_verbose);
static INT autlog_add_generator(
	AUTLOG_THE_GROUP *B, 
	INT i, INT g, INT f_verbose);
static INT autlog_reduce(AUTLOG_LOCAL *L, 
	AUTLOG_THE_GROUP *A, 
	AUTLOG_THE_GROUP *B, 
	AUTLOG_GRID *G1, AUTLOG_GRID *G2, 
	ORDERED_SET *E, INT k, INT f_verbose);
static INT autlog_commutator_test(AUTLOG_LOCAL *L, 
	AUTLOG_THE_GROUP *A, 
	AUTLOG_THE_GROUP *B, 
	INT k, INT *offset, 
	INT ib, INT ie, 
	INT jb, INT je, INT f_verbose);
static INT aut_table_test(
	AUTLOG_THE_GROUP *A, 
	AUTLOG_THE_GROUP *B, 
	INT go_last, INT go, 
	INT ib, INT ie, 
	INT jb, INT je, INT f_verbose);
static INT autlog_dimino(AUTLOG_THE_GROUP *C, 
	INT k, INT g0, 
	INT f_test_it, 
	INT go_soll /* before: A_N2 */, 
	AUTLOG_THE_GROUP *A, 
	AUTLOG_GRID *C_grid, AUTLOG_GRID *A_grid, 
	INT *go);
static INT autlog_add_coset(AUTLOG_THE_GROUP *C, 
	INT f_test_it, AUTLOG_THE_GROUP *A, 
	AUTLOG_GRID *C_grid, AUTLOG_GRID *A_grid, 
	INT coset_size, INT g2, INT *go);
static INT autlog_mult(
	AUTLOG_THE_GROUP *A, INT i, INT j);
static INT autlog_mult_G(
	AUTLOG_THE_GROUP *G, INT i, INT j);
static INT autlog_invers_G(
	AUTLOG_THE_GROUP *G, INT i);
static INT autlog_conjugate_G(
	AUTLOG_THE_GROUP *G, INT i, INT j);
/* j^-1 * i * j */
static INT autlog_commutator_G(
	AUTLOG_THE_GROUP *G, INT i, INT j);
/* i^-1 *  j^-1 * i * j */
#if 0
static INT aut_table_compare(
	AUTLOG_THE_GROUP *A, 
	INT go_last, INT go, 
	INT ib, INT ie, INT jb, INT je, 
	SPERM *p, SPERM *pv, SPERM *p_, SPERM *pv_, 
	INT f_verbose);
#endif


/* aut_grid.C: */
static INT autlog_calc_grid(AUTLOG_THE_GROUP *C, 
	AUTLOG_GRID *G, INT k, INT go, INT f_verbose);
static INT autlog_radix_sort(AUTLOG_GRID *G, 
	INT radix, INT first, INT last);
static INT autlog_insert_idx(AUTLOG_GRID *G, 
	INT first, INT len, INT radix, 
	INT search_this, INT *idx);
static void autlog_print_grid(AUTLOG_GRID *G, SPERM *pv, 
	INT k, INT go_k, INT go);
static void autlog_print_type_vector(AUTLOG_GRID *G, INT i, INT k);
static INT autlog_calc_all_grids(AUTLOG_THE_GROUP *C, 
	INT go, INT f_verbose);
static INT autlog_calc_n_grids(AUTLOG_THE_GROUP *C, 
	INT go, INT n, INT f_verbose);
static INT autlog_E(
	AUTLOG_THE_GROUP *C1, AUTLOG_THE_GROUP *C2, 
	AUTLOG_GRID *G1, AUTLOG_GRID *G2, 
	INT k, ORDERED_SET *E, INT f_verbose);
static INT autlog_check_iso(
	AUTLOG_GRID *G1, AUTLOG_GRID *G2);
static INT autlog_choose_first_into_E(AUTLOG_THE_GROUP *G, 
	AUTLOG_GRID *Grid, INT k, ORDERED_SET *E, INT f_verbose);
static INT autlog_choose_E(
	AUTLOG_THE_GROUP *C1, AUTLOG_THE_GROUP *C2, 
	AUTLOG_GRID *G1, AUTLOG_GRID *G2, 
	INT k, ORDERED_SET *E, INT f_verbose);
static void autlog_print_E(
	ORDERED_SET *E, INT k);


