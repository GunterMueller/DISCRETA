/* geo.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef GEO_INCLUDED
#define GEO_INCLUDED

#ifndef CP_INCLUDED
#include <DISCRETA/cp.h>
#endif
#ifndef MA_INCLUDED
#include <DISCRETA/ma.h>
#endif

#define MAX_V 80
#define MAX_B 80
#define MAX_R 10 // used in iso.C
#define MAX_VB 80   /* MAX(MAX_V, MAX_B) */


typedef struct geo_canon_info GEO_CANON_INFO;

typedef struct geo_data GEO_DATA;

typedef struct iso_geo_data ISO_GEO_DATA;
typedef struct iso_info ISO_INFO;
typedef struct grid GRID;
typedef struct frame FRAME;

typedef struct tdo_scheme TDO_SCHEME;
typedef SHORT TDOSS;

#define MAX_GRID 80

#define TDOSS_OFF   ((2 * sizeof(INT *))/sizeof(SHORT))

struct tdo_scheme {
	INT m, n;
	INT *a;
};

/* isot.C */

struct iso_geo_data {
	INT dim_n;
	INT *theX; /* [v][dim_n] */
	INT f_pc;
	BYTE *pctheX;
		/* only in ISO_INFO allowed */
	INT f_R_allocated;
	INT *R; 
	/* R[i] is number of incidences in row i */

	INT v, b; /* the size of this geometry */
	INT V, B; /* the size of the larger geometry (B) */
	/* INT max_r; */

	/* second derivatives 
	 * (a map from the unordered pairs
	 * of points (p) or blocks (b)
	 * into integers, 
	 * ij2k() used for indexing) */
	INT f_use_ddp;
	INT f_use_ddb;
	SHORT *ddp, *ddb;

	/* tdo data: */
	INT f_tdo;
	INT tdo_m;
	INT tdo_V[MAX_V];
	INT tdo_n;
	INT tdo_B[MAX_B];

	/* additional coloring of 
	 * the points (p) / blocks (b)
	 * (multiple color fields allowed) */
	INT f_colors;
	INT nb_bcol, nb_pcol;
	INT *bcol, *pcol;

	/* current row-permutation (degree V), 
	 * used in ISO2;
	 * pv = p^{-1} */
	CPERM p, pv;
	/* current column-permutation (degree B), 
	 * qv = q^{-1} */
	CPERM q, qv;

	INT f_transpose_it;

	INT hbar[MAX_VB];
	INT hlen[MAX_VB];
	INT hlen01[MAX_VB];
	INT hlen1[MAX_VB];
	INT grid_entry[MAX_VB];
	INT G_max;
};

struct iso_info {
	INT *AtheX;
	/* v x max_r; 
	 * dimension v x max_r or 
	 * v x MAX_R */
	INT *BtheX;
	INT Af_full;
	INT Bf_full;
	BYTE *pcAtheX;
	BYTE *pcBtheX;
	INT f_pcA;
	INT f_pcB;
	
	INT v;
	INT b;
	INT max_r;
	
	INT *R; /* [MAX_V] */

#if 0
	/* unused: */
	INT nb_i_vbar;
	INT *i_vbar; /* [MAX_JJ] */
	INT nb_i_hbar;
	INT *i_hbar; /* [MAX_II] */
#endif
	
	INT tdo_m;
	INT tdo_V[MAX_V];
	INT tdo_n;
	INT tdo_B[MAX_B];
	
	INT nb_isomorphisms;
	INT f_break_after_fst;
	INT f_verbose;
	INT f_very_verbose;
	INT f_use_d;
	INT f_use_ddp;
	INT f_use_ddb;
	INT f_transpose_it;
	
	/* optionally: */
	INT *Ar; /* v entries */
	INT *Br;
	INT *Ad; /* v entries */
	INT *Bd;
	SHORT *Addp; /* (v \atop 2) entries */
	SHORT *Bddp;
	SHORT *Addb; /* (b \atop 2) entries */
	SHORT *Bddb;

	INT f_igd; /* use igd instead */
	ISO_GEO_DATA A, B;
};


/* ntdo.h */


#define THE_Y_OFFSET (2 * sizeof(void *) / sizeof(SHORT))

#define NTDO_MAX_N 300

typedef struct ntdo_info NTDO_INFO;
typedef struct ntdo NTDO;
typedef struct ntdo_grid NTDO_GRID;
typedef struct tdo_scheme TDO_SCHEME;
typedef struct ntdo_frame NTDO_FRAME;
typedef struct tdo_grad TDO_GRAD;


struct ntdo_info {
	INT v, b, nb_X;
	INT *theX;
	INT f_multivalued;
	INT *theVal;
	INT llambda, lmue;
	INT *lambda_0;
	INT *mue_0;
	INT f_transposed;
};

struct ntdo {
	NTDO_INFO *info;
	
	/* the incidence matrix, 
	 * eventually transposed: */
	INT v, b;
	INT nb_X;
	INT *theXi;
	INT *theXj;
	INT f_multivalued;
	INT *theVal;
	INT max_width;
	BYTE format_string[64];
	INT llambda, lmue;
	INT *lambda_0, *mue_0;

	INT max_size;
	NTDO_GRID *G_last, *G, *G_next;
	
	SPERM p, pv, q, qv;
	
	TDO_SCHEME tdos;
	TDO_SCHEME tdos2;
};

struct ntdo_grid {
	INT max_size;
	INT f_points;
	INT m;
	/* # Zeilen in type 
	 * (= v if f_point, = b if !f_points) */
	INT n;
	/* # Eintraege pro Zeile in type, (= G->G_max) */
	INT G_max;
	INT *first;
	INT *len;
	INT *type_idx;
	INT *grid_entry;
	/* the index into first / len, where to find 
	 * the row / column. */
	INT *type;
	INT f_perms_allocated;
	SPERM p, pv; /* degree m */
};

struct ntdo_frame {
	INT G_max;
	INT first[NTDO_MAX_N + 1];
	INT len[NTDO_MAX_N];
	INT grid_entry[NTDO_MAX_N];
};

struct tdo_grad {
	INT f_points;
	INT n;
	/* = tdo->v if f_points, 
	 * = tdo->inc->B otherwise */
	INT f_dd;
	INT N;
	/* = n * (n-1) / 2 if f_dd, 
	 * = n otherwise */
	INT *type; /* type[N] */

	INT nb_tdos;
	TDO_SCHEME *tdos;
	INT *mult;
};

// tdo_scheme.h

typedef struct ntdo NTDO;
typedef struct ntdo_grid NTDO_GRID;


/* iso_sub.C: */
void igd_nil(ISO_GEO_DATA *geo);
void igd_exit(ISO_GEO_DATA *geo);
void igd_init_theX(
	ISO_GEO_DATA *geo, INT *theX, INT *R, INT dim_n, 
	INT v, INT b);
void igd_init_pctheX(
	ISO_GEO_DATA *geo, BYTE *theX, INT *R, INT dim_n, 
	INT v, INT b);
INT igd_init_dd(
	ISO_GEO_DATA *geo, 
	INT f_ddp, SHORT *ddp, 
	INT f_ddb, SHORT *ddb);
void igd_init_color(ISO_GEO_DATA *p, 
	INT nb_pcol, INT nb_bcol, INT *pcol, INT *bcol);
void igd_init_tdo(
	ISO_GEO_DATA *p, TDOSS *tdoss);
void igd_init_tdo_V_B(ISO_GEO_DATA *p, 
	INT V, INT B, INT *Vi, INT *Bj);
INT igd_init_igd(
	ISO_GEO_DATA *geo, ISO_GEO_DATA *geo_info);
INT igd_init_hbar(ISO_GEO_DATA *p);
void igd_init_cperms(ISO_GEO_DATA *p);
INT igd_find_incidence(ISO_GEO_DATA *p, 
	INT i, INT j);
INT igd_print_theX(ISO_GEO_DATA *p);
INT igd_add_point(ISO_GEO_DATA *p, 
	INT k, INT j0, INT f_verbose);
INT igd_join(ISO_GEO_DATA *A, ISO_GEO_DATA *B, 
	INT k, INT f_verbose);
INT igd_del_point(ISO_GEO_DATA *p, 
	INT k, INT f_A, INT f_verbose);
void igd_set_ge(ISO_GEO_DATA *p);
INT subiso_test(ISO_INFO *info);

/* iso.C: */
INT iso_info_init_A_INT(
	ISO_INFO *info, INT *theA, INT f_full);
INT iso_info_init_A_BYTE(
	ISO_INFO *info, BYTE *theA, INT f_full);
INT iso_info_init_B_INT(
	ISO_INFO *info, INT *theB, INT f_full);
INT iso_info_init_B_BYTE(
	ISO_INFO *info, BYTE *theB, INT f_full);
INT iso_info_init_ddp(
	ISO_INFO *info, INT f_ddp, 
	SHORT *Addp, SHORT *Bddp);
INT iso_info_init_ddb(
	ISO_INFO *info, INT f_ddb, 
	SHORT *Addb, SHORT *Bddb);
INT iso_info_init_tdo(
	ISO_INFO *info, TDOSS *tdoss);
INT iso_info_init_tdo_V_B(ISO_INFO *info, 
	INT V, INT B, INT *Vi, INT *Bj);
INT init_ISO2(void);
INT iso_test(ISO_INFO *info);
INT print_theX_pq(INT *theX, INT dim_n, 
	INT v, INT b, INT *R, 
	CPERM *pv, CPERM *qv);
INT print_theX(INT *theX, 
	INT dim_n, INT v, INT b, INT *R);

struct grid {
	INT f_points;
	INT m;
	/* # Zeilen in type 
	 * (= v if f_point, = b if !f_points) */
	INT n;
	/* # Eintraege pro Zeile 
	 * in type, (= G->G_max) */
	INT G_max;
	INT first[MAX_GRID + 1];
	INT len[MAX_GRID];
	INT type_idx[MAX_GRID];
	INT grid_entry[MAX_GRID];
	/* the index into first / len, where to find 
	 * the row / column. */
	INT type[MAX_GRID][MAX_GRID];
};

struct frame {
	INT G_max;
	INT first[MAX_GRID + 1];
	INT len[MAX_GRID];
	INT grid_entry[MAX_GRID];
};



/*
 * TDO stuff
 */

INT tdos_init_decomposition(TDO_SCHEME *tdos, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp);
INT tdos_init_grid(TDO_SCHEME *tdos, 
	INT v, INT b, 
	INT nb_i_hbar, INT *i_hbar, 
	INT nb_i_vbar, INT *i_vbar);
INT tdos_free(TDO_SCHEME *tdos);
void tdos_nil(TDO_SCHEME *tdos);
INT tdos_m(TDO_SCHEME *tdos);
INT tdos_n(TDO_SCHEME *tdos);
/* m x n ist die Dimension der TDO Matrix; 
 * abgespeichert wird es in 
 * einer (m + 1) x (n + 1) Matrix, 
 * wobei die unterste Zeile und 
 * die rechteste Spalte die 
 * Bereichslaengen Vi bzw. Bj angeben. */
INT tdos_Vi(TDO_SCHEME *tdos, INT i);
INT tdos_Bj(TDO_SCHEME *tdos, INT j);
INT tdos_ij(TDO_SCHEME *tdos, INT i, INT j);
void tdos_copy(TDO_SCHEME *t1, TDO_SCHEME *t2);
void tdos_print(TDO_SCHEME *tdos);
INT tdos_cmp(TDO_SCHEME *t1, TDO_SCHEME *t2);

TDOSS *tdoss_init_decomposition(VECTOR_OP row_decomp, VECTOR_OP col_decomp);
void tdoss_free(TDOSS *t);
void tdoss_print(TDOSS *t);
INT tdoss_cmp(TDOSS *t1, TDOSS *t2);
INT tdos2tdoss(TDO_SCHEME *tdos, 
	SHORT *ddp_mult, SHORT *ddb_mult, TDOSS **t);
SHORT *tdoss_get_ddp_mult(TDOSS *t);
SHORT *tdoss_get_ddb_mult(TDOSS *t);
void tdoss_set_ddp_mult(TDOSS *t, SHORT *ddp_mult);
void tdoss_set_ddb_mult(TDOSS *t, SHORT *ddb_mult);
INT tdoss_m(TDOSS *t);
INT tdoss_n(TDOSS *t);
INT tdoss_Vi(TDOSS *t, INT i);
INT tdoss_Bj(TDOSS *t, INT j);
INT tdoss_V(TDOSS *t);
INT tdoss_B(TDOSS *t);
INT tdoss_ij(TDOSS *t, INT i, INT j);
INT tdoss_2_first(TDOSS *t, VECTOR_OP row_first, VECTOR_OP block_first);

INT tdos_init_ntdo(TDO_SCHEME *tdos, NTDO *tdo, 
	NTDO_GRID *G0, NTDO_GRID *G1, INT f_derived);
INT tdos_print_theX_short(TDO_SCHEME *tdos, INT nb_X, SHORT *theX);
INT tdos_print(TDO_SCHEME *tdos, NTDO *tdo);
INT TD_char_inv_print(FILE *fp, INT f_tex, 
	INT nrow, INT ncol, INT nb_X, INT *X, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_char_first, VECTOR_OP col_char_first, 
	VECTOR_OP row_inv_first, VECTOR_OP col_inv_first);
INT TD_print(FILE *fp, INT nrow, INT ncol, INT nb_X, INT *X, INT f_transposed, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP P_first, VECTOR_OP Q_first);
void grid_print(INT v, INT *hbar, BYTE **str);
void grid_fprint(FILE *fp, INT f_tex, INT v, INT *hbar, BYTE **str);
void grid_alloc(BYTE ***str, INT num, INT len);
void grid_free(BYTE **str, INT num);
void grid_prepare(BYTE **str, INT v, INT b, INT dim_M, INT *M, INT *vbar, INT *hbar);
void grid_prepare2(BYTE **str, INT v, INT b, INT dim_M, INT *M, INT *vbar, INT *hbar);

INT tdoss_print_theX_short(TDOSS *tdoss, INT nb_X, SHORT *theX);
INT incma_latex_picture(FILE *fp, 
	INT width, INT width_10, 
	INT f_outline_thin, BYTE *unit_length, 
	BYTE *thick_lines, BYTE *thin_lines, BYTE *geo_line_width, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X, 
	INT f_labelling_points, BYTE **point_labels, 
	INT f_labelling_blocks, BYTE **block_labels);
INT incma_latex(FILE *fp, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X);
INT incma_latex_integer_labels(FILE *fp, 
	INT v, INT b, 
	INT V, INT B, INT *Vi, INT *Bj, 
	INT *R, INT *X, INT dim_X, 
	VECTOR_OP point_labels, VECTOR_OP block_labels, INT offset);


/* ntdo.C: */
TDOSS *calc_ntdo_(INT v, INT b, 
	INT nb_X, INT *theX, INT f_multivalued, INT *theVal, 
	INT f_calc_second_tdo, 
	INT f_ddp, VECTOR_OP DDp, 
	INT f_ddb, VECTOR_OP DDb, 
	INT f_v, INT f_vv, INT f_dd_v, INT f_dd_vv, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp);
TDOSS *calc_ntdo(INT v, INT b, INT nb_X, INT *theX, 
	INT f_calc_second_tdo, 
	INT f_ddp, VECTOR_OP DDp, 
	INT f_ddb, VECTOR_OP DDb, 
	INT f_v, INT f_vv, INT f_dd_v, INT f_dd_vv);
void ntdo_the_Y_free(SHORT *theY);
INT ntdo_calc_theY(NTDO *tdo, SHORT *theY, SHORT *ddp, SHORT *ddb, INT f_v);
void the_Y_set_ddp(SHORT *theY, SHORT *ddp);
void the_Y_set_ddb(SHORT *theY, SHORT *ddb);
SHORT *the_Y_get_ddp(SHORT *theY);
SHORT *the_Y_get_ddb(SHORT *theY);
SHORT *the_Y_get_the_X(SHORT *theY);
void ntdo_print(NTDO *tdo);
void ntdo_print_row(NTDO *tdo, INT *vbar, INT *hbar, INT i, INT *M, INT dim_n);
void ntdo_top_border_row(NTDO *tdo, INT *vbar, INT *hbar);
void ntdo_top_middle_border_row(NTDO *tdo, INT *vbar, INT *hbar);
NTDO *ntdo_init(NTDO_INFO *info, INT f_verbose);
void ntdo_free(NTDO *N);
void ntdo_print_incidence_matrix(NTDO *N);
INT ntdo_calc(NTDO *tdo, INT f_v);
INT ntdo_next(NTDO *tdo, INT f_v);
INT ntdo_collect_types(NTDO *tdo, NTDO_GRID *G0, NTDO_GRID *G1, INT f_v);
INT ntdo_calc_grid(NTDO *tdo, NTDO_GRID *G0, NTDO_GRID *G1);
INT ntdo_refine_types(NTDO *tdo, NTDO_GRID *Gm1, NTDO_GRID *G1);
INT ntdo_radix_sort(NTDO *tdo, NTDO_GRID *G, INT radix, INT first, INT last);
INT ntdo_insert_idx(NTDO_GRID *G, INT first, INT len, INT radix, 
	INT search_this, INT *idx);

/* ntdo_grid.C:*/
void ntdo_grid_nil(NTDO_GRID *G);
NTDO_GRID *ntdo_grid_init(INT max_size);
void ntdo_grid_free(NTDO_GRID *G);
INT ntdo_init_perms(NTDO_GRID *G, INT degree);
INT ntdo_grid_init_partition(NTDO_GRID *G, 
	INT f_points, INT m, INT n, INT len, INT *parts);
INT ntdo_grid_init0(NTDO *tdo, NTDO_GRID *Gpoints, NTDO_GRID *Gblocks);
INT ntdo_grid_print(NTDO_GRID *G);
INT ntdo_grid_init_derived_i_first(NTDO *tdo, 
	NTDO_GRID *G, NTDO_GRID *G_old, INT derive_at_i, INT f_v);
INT ntdo_grid_init_derived_ij_first(NTDO *tdo, 
	NTDO_GRID *G, NTDO_GRID *G_old, INT I, INT J);
void ntdo_grid_copy_frame(NTDO_GRID *G1, NTDO_GRID *G2, INT f_v);

void frame2ntdo_grid(NTDO_FRAME *frame, NTDO_GRID *grid);

void tdog_exit(TDO_GRAD *p);
void tdog_nil(TDO_GRAD *p);
INT tdog_add_tdos(TDO_GRAD *p, TDO_SCHEME *tdos, INT i);

/* ntdo2.C: */
INT ntdo_calc_second_tdo(NTDO *tdo, INT f_v, INT f_vv);
INT ntdo_refine(NTDO *tdo, TDO_GRAD *tdog, 
	NTDO_GRID *G, NTDO_GRID *G_next, 
	INT f_points, NTDO_FRAME *frame, 
	SPERM *p, SPERM *pv, SPERM *q, SPERM *qv, INT f_v);

/* ntdo_dd.C: */
INT ntdo_calc_dd(NTDO *tdo, INT f_points, INT f_blocks, 
	SHORT **ddp, INT *Np, SHORT **ddp_mult, 
	SHORT **ddb, INT *Nb, SHORT **ddb_mult, 
	INT f_v, INT f_vv);
INT ntdo_dd(NTDO *tdo, INT f_points, SHORT **dd, INT *N, SHORT **dd_mult, INT f_v, INT f_vv);




/*
 *
 */



// tda.C:


void tda_geo_tex_head(FILE *fp);
void tda_geo_tex_foot(FILE *fp);
INT calc_TDO(FILE *fp, INT nrow, INT ncol, INT nb_X, INT *X, 
	INT f_calc_second_tdo, TDOSS **tdoss);
INT calc_TDA(FILE *TDA_fp, FILE *TEX_fp, INT nb_geo, 
	INT nrow, INT ncol, INT nb_X, INT *X, 
	INT f_transposed, LABRA_OP A, 
	INT *nb_row_blocks_TDO, INT *nb_col_blocks_TDO, 
	INT *nb_row_blocks_TDA, INT *nb_col_blocks_TDA, 
	INT f_v, INT f_vv, INT f_incidences);
INT geo_output_char_inv_decomposition(FILE *fp, INT f_tex, INT nb_geo, 
	INT nrow, INT ncol, INT nb_X, INT *X, 
	PERMUTATION_OP pp, PERMUTATION_OP qq, 
	VECTOR_OP row_char_first, 
	VECTOR_OP col_char_first, 
	VECTOR_OP row_inv_first, 
	VECTOR_OP col_inv_first, 
	SYM_OP ago, VECTOR_OP generators, INT f_incidences);
INT calc_blocks(INT nrow, INT ncol, INT nb_X, INT *X, INT f_transposed, 
	PERMUTATION_OP p, PERMUTATION_OP pv, 
	PERMUTATION_OP q, PERMUTATION_OP qv, 
	VECTOR_OP B, INT f_v);


/* geo_canon.h */

struct geo_canon_info {
	INT v, b;
	INT nb_X;
	INT *theX;
	INT f_transposed;
	TDOSS *tdoss;
	INT f_ddp, f_ddb;
	SHORT *ddp, *ddb;
	
	INT f_col_group;
	LABRA_OP col_group;
	MATRIX_OP col_transversals;
	
	INT f_print_backtrack_points;
	INT backtrack_points_mod;
	INT f_verbose;
	INT f_very_verbose;
	// INT f_aut_gens; /* auts and first_moved */
	INT f_get_aut_group; /* Aut and ago */
	INT f_aut_v, f_aut_vv;

	/* Output: */
	// VECTOR_OB auts;
	// VECTOR_OB first_moved;
	LABRA_OP aut;
	SYM_OP ago;
};

INT geo_maxtest_set(INT n, INT k, INT *theX, 
	LABRA_OP G, MATRIX_OP TG, 
	LABRA_OP aut, SYM_OP ago, INT f_v, INT f_vv);
INT geo_canonicize_set(INT n, INT k, INT *theX, 
	LABRA_OP G, MATRIX_OP TG, INT f_v, INT f_vv, 
	INT f_get_aut_group, LABRA_OP aut, SYM_OP ago, 
	PERMUTATION_OP transporter);
INT geo_canon_with_initial_decomposition(INT f_maxtest, INT *back_to, 
	INT nrow, INT ncol, INT nb_X, INT f_print_backtrack_points, INT *theX, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_transposed, INT f_get_aut_group, LABRA_OP aut, 
	INT f_v, INT f_vv);
INT geo_canon_with_initial_decomposition_and_ddp_ddb(
	INT f_maxtest, INT *back_to, 
	INT nrow, INT ncol, INT nb_X, 
	INT f_print_backtrack_points, INT *theX, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb,
	VECTOR_OP DDp, VECTOR_OP DDb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_transposed, INT f_get_aut_group, LABRA_OP aut, 
	INT f_v, INT f_vv);
INT geo_canon_simple(INT f_maxtest, INT *back_to,
	INT nrow, INT ncol, INT nb_X, INT f_print_backtrack_points, INT *theX, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_transposed, INT f_get_aut_group, LABRA_OP aut, 
	INT f_v, INT f_vv);
INT geo_Canonicize(INT f_maxtest, INT *back_to, 
	INT nrow, INT ncol, INT nb_X, INT f_print_backtrack_points, 
	INT *theY, TDOSS *tdoss, INT f_transposed, 
	INT f_ddp, SHORT *ddp, 
	INT f_ddb, SHORT *ddb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_get_aut_group, LABRA_OP aut, SYM_OP ago, INT f_v, INT f_vv);
INT geo_canonicize(GEO_CANON_INFO *info, PERMUTATION_OP p, PERMUTATION_OP q);



// geo_data.h

#define GEO_MAX_N 300

struct geo_data {
	INT f_second_kind;
	INT f_not_yet_tested;
	INT nb_X;
	CHUNK_MEMH *geodata;
	INT chunk_size; /* for geodata */
	INT nb_cand;
	INT nb_geo;
	
	/* flags for computation: */
	INT f_dont_go_further;
	INT f_dont_go_further_without;
	INT f_transposed;
	INT f_calc_second_tdo;
	INT f_calc_ddp;
	INT f_calc_ddb;
	
	/* verbose-flags for computation: */
	INT f_tdo_v;
	INT f_tdo_vv;
	INT f_dd_v;
	INT f_dd_vv;
	INT f_canon_v;
	INT f_canon_vv;
	
	/* flags for geo_add: */
	INT f_cand_verbose;
	INT f_add_verbose;
	INT f_print_candidates, print_candidates_mod;
	INT f_print_canonicized_candidates, print_canonicized_candidates_mod;
	INT f_print_geo, print_geo_mod;
	INT f_save_into_geo_file;
	BYTE geo_fname[1024];
	

	INT f_range, range_first, range_len;
	INT geo_hdl;
};

/* geo_data.C: */

INT geo_data_cmp(INT nb_X, SHORT *theX, SHORT *theY);
INT geo_data_add(GEO_DATA *gd, SHORT *theX);
INT geo_data_search_and_add(GEO_DATA *gd, SHORT *theX, INT *geo_hdl);

GEO_DATA *init_geo_data(INT nb_incidences);
void geo_data_print_range(GEO_DATA *gd);
void geo_data_print_header(GEO_DATA *gd, INT f_v);
void geo_data_print(GEO_DATA *gd, INT nrow, INT ncol);
void geo_tree_print(GEO_DATA *gd, INT hdl, INT nrow, INT ncol);
void geo_tree_node_print(GEO_DATA *gd, INT hdl, INT nrow, INT ncol);
INT multivalued_geo_tdo(MATRIX_OP M, MATRIX_OP N, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp);
void geo_print(GEO_DATA *gd, INT nb_X, SHORT *theX, INT nrow, INT ncol);

/* dcc_orderly.C: */
INT calc_kramer_mesner_matrix_by_orderly_generation(
	VECTOR_OP gen, SYM_OP go, INT deg, 
	LABRA_OP labG, MATRIX_OP TG, 
	BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, INT f_TDO);
INT dcc_orderly(LABRA_OP labG, MATRIX_OP TG, 
	INT n, INT k, VECTOR_OP Orbits, VECTOR_OP Stabs);
INT calc_KM_column(LABRA_OP G, MATRIX_OP TG, 
	VECTOR_OP Reps_t, VECTOR_OP Stabs_t, VECTOR_OP k_set, LABRA_OP k_stab, 
	MATRIX_OP M, INT n, INT t, INT k, INT k_idx);
INT calc_KM_tk_column_wise(LABRA_OP G, MATRIX_OP TG, 
	VECTOR_OP Reps_t, VECTOR_OP Stabs_t, 
	VECTOR_OP Reps_k, VECTOR_OP Stabs_k, 
	MATRIX_OP M, INT n, INT t, INT k);
INT calc_KM_tk(LABRA_OP G, MATRIX_OP TG, 
	VECTOR_OP Reps_t, VECTOR_OP Reps_k, MATRIX_OP M, INT n, INT t, INT k);
INT extend_to_k_sets(LABRA_OP G, MATRIX_OP TG, INT n, INT k, 
	VECTOR_OP km1_sets, VECTOR_OP km1_stabs, 
	VECTOR_OP k_sets, VECTOR_OP k_stabs, 
	INT f_v, INT f_vv);
INT calc_k_sets(LABRA_OP G, MATRIX_OP TG, INT n, INT k, 
	VECTOR_OP Reps, VECTOR_OP Stab, INT f_v, INT f_vv);
INT calc_k_stab_orders(VECTOR_OP k_stabs, VECTOR_OP k_stab_orders);
INT n_choose_k_first(INT *X, INT n, INT k);
INT n_choose_k_next(INT *X, INT n, INT k);
INT n_choose_k_next_at(INT *X, INT n, INT k, INT a);
void print_set(INT *X, INT n, INT k);
INT KM_tk_print_asc(LABRA_OP G, BYTE *g_label, BYTE *g_label_tex, 
	MATRIX_OP M, INT gl_t, INT gl_k);


/* geo_store.h */


class geo_by_base_blocks_ob : public VECTOR_OB {
public:
	STRING_OP s_hash() {
		return((STRING_OP) s_i(0)); };
	MATRIX_OP s_I() {
		return((MATRIX_OP) s_i(1)); };
	INTEGER_OP s_I_ij(INT i, INT j) {
		return((INTEGER_OP) s_I()->s_ij(i, j)); };
	INT s_I_iji(INT i, INT j) {
		return(s_I_ij(i, j)->s_i()); };
	VECTOR_OP s_BB() {
		return((VECTOR_OP) s_i(2)); };
	VECTOR_OP s_BB_i(INT i) {
		return((VECTOR_OP) s_BB()->s_i(i)); };
	INTEGER_OP s_BB_ij(INT i, INT j) {
		return((INTEGER_OP) s_BB_i(i)->s_i(j)); };
	INT s_BB_iji(INT i, INT j) {
		return(s_BB_ij(i, j)->s_i()); };
	INTEGER_OP s_id() {
		return((INTEGER_OP) s_i(3)); };
	INT s_id_i() {
		return s_id()->s_i(); };

INT field_name(INT i, INT j, BYTE *str);
INT init(VECTOR_OP BB, INT id);
INT print();
INT sprint(BYTE *s);
INT fprint_incidences(FILE *fp);
INT nb_incidences();
INT border(INT n);
INT span_design(VECTOR_OP gen, INT f_v, INT f_vv);
INT canonical_form(INT f_tdo, INT f_tdo2, 
	INT f_V, VECTOR_OP V, 
	INT f_B, VECTOR_OP B, 
	INT f_v, INT f_vv);
INT calc_hash(INT f_v, INT f_vv);
INT Equal(GEO_BY_BASE_BLOCKS_OP G1, INT f_v);
};

INT i_geo_by_base_blocks_db(DATABASE_OP *db, BYTE *db_prefix);
void e_geo_by_base_blocks_db(DATABASE_OP db);
INT db_geo_by_base_blocks_create_and_close(BYTE *db_prefix);
INT db_geo_by_base_blocks_find_isomorphic_via_hash(DATABASE_OP db, 
	FILE *fp_txt, GEO_BY_BASE_BLOCKS_OP G, INT f_verbose, 
	INT *idx_btree_hash, INT *f_found);
INT db_geo_by_base_blocks_export(BYTE *db_prefix, BYTE *geo_fname);
INT db_geo_by_base_blocks_export_id(BYTE *db_prefix, BYTE *id_fname);
INT db_geo_by_base_blocks_add(BYTE *db_prefix, INT f_0, INT f_first, 
	BYTE *group_fname, BYTE *base_block_fname, 
	INT f_inc_file, BYTE *inc_file_name, 
	INT f_tdo, INT f_tdo2, 
	INT f_tda, INT f_range, INT first, INT len, 
	INT f_test2design, 
	INT f_V, VECTOR_OP decomp_V, 
	INT f_B, VECTOR_OP decomp_B, 
	INT f_do_not_add, 
	INT f_v, INT f_vv);





#endif /* GEO_INCLUDED */

