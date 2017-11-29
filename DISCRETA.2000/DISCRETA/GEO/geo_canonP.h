/* geo_canonP.h */

#ifndef GEO_CANONP_INCLUDED
#define GEO_CANONP_INCLUDED

#ifndef GEO_INCLUDED
#include <DISCRETA/geo.h>
#endif

#undef DEBUG_HBAR

#define MAX_GEO_POINTS 1000
#define MAX_SET_SIZE MAX_GEO_POINTS

#define GEO_CANON_ALLOC_M

typedef struct ordered_set ORDERED_SET;
typedef struct geo_canon GEO_CANON;

struct ordered_set {
	INT a[MAX_SET_SIZE];
	INT size;
};

struct geo_canon {
	GEO_CANON_INFO *info;
	INT v, b, max_size;
	INT nb_X;
	INT *theXi;
	INT *theXj;
#ifdef GEO_CANON_ALLOC_M
	INT *M; /* [v][b] */
#endif
	INT *theX;
	INT *theY;
	INT *nb_x;
	INT *nb_y;
	
	INT tdo_m;
	INT tdo_n;
	INT *tdo_V;
	INT *tdo_B;
	ORDERED_SET *E; /* b */
	
	INT f_ddp;
	INT f_ddb;

	SHORT *ddp; /* (v \atop 2) entries */
	SHORT *ddb; /* (b \atop 2) entries */

	INT nb_P, nb_Q;
	NTDO_GRID **P, **Q;
#if 0
	NTDO_GRID *P[MAX_GEO_POINTS]; /* decomposition of the rows (points) */
	NTDO_GRID *Q[MAX_GEO_POINTS]; /* decomposition of the columns (blocks) */
#endif
	
	INT is_first;
	SPERM p0, p0v;
	SPERM q0, q0v;
	SPERM tmp_q1, tmp_q2;
	SPERM tmp_q3, tmp_q4;

	INT first_moved;
		/* only used if mode 
		 * is ONLY_COSET_REPS */

		// used for maxtest:
	INT f_going_back;
	INT back_to;
	
	INT nb_aut_gens;
	INT dim_aut_gens;
	SPERM *aut_gens;
	INT *auts_first_moved;
	INT nb_backtrack;
	INT nb_backtrack_points;

};

// geo_maxtest.C:
static INT geo_maxtest(GEO_CANON_INFO *info);
INT geo_maxtest_recursion(GEO_CANON *geo, INT level);

// geo_canon2.C:
INT geo_canonicize_recursion(GEO_CANON *geo, INT level);

// geo_util.C:
GEO_CANON *geo_canon_init(GEO_CANON_INFO *info);
void geo_canon_free(GEO_CANON *GC);
void geo_canon_free_PQ(GEO_CANON *GC);
INT geo_print(GEO_CANON *geo, INT k);
INT geo_get_aut_gens(GEO_CANON *geo, 
	VECTOR_OP auts, VECTOR_OP first_moved);
INT geo_get_aut_group(GEO_CANON *geo, 
	LABRA_OP Aut, SYM_OP ago, INT f_v, INT f_vv);
INT geo_get_aut_gen(GEO_CANON *geo);
void geo_free_aut_gen(GEO_CANON *geo);
INT geo_get_canonical_labelling(GEO_CANON *geo);
INT geo_canon_apply(GEO_CANON *geo);
INT geo_get_M_ij(GEO_CANON *GC, INT i, INT j);
INT geo_compare(GEO_CANON *geo, INT k);
INT geo_compare_col(GEO_CANON *geo, INT k, INT col);
INT print_grid_decomp(NTDO_GRID *Q);


INT geo_calc_grid(GEO_CANON *geo, INT k, INT f_v);
INT geo_radix_sort(GEO_CANON *geo, NTDO_GRID *Q, INT radix, 
	INT first, INT last);
INT geo_insert_idx(NTDO_GRID *G, INT first, INT len, 
	INT radix, INT search_this, INT *idx);
INT geo_choose_first_into_E(GEO_CANON *geo, 
	NTDO_GRID *Grid, NTDO_GRID *Grid_last, INT k, ORDERED_SET *E, INT f_v);
INT geo_add(GEO_CANON *geo, INT k, INT j, INT f_v);
INT geo_check_P(GEO_CANON *geo, INT k);

#endif /* GEO_CANONP_INCLUDED */



