/* geo_data.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>
#include <DISCRETA/ma.h> // for multivalued_geo

#include <DISCRETA/geo.h>
#include <DISCRETA/lb.h>

#undef DEBUG_GOING_BACK
#define BACKTEST_POINTS

INT geo_data_cmp(INT nb_X, SHORT *theX, SHORT *theY)
{
	INT i, ret;
	
	for (i = 0; i < nb_X; i++) {
		ret = theY[i] - theX[i];
		if (ret == 0)
			continue;
		if (ret > 0)
			return 1;
		if (ret < 0)
			return -1;
		}
	return 0;
}

INT geo_data_add(GEO_DATA *gd, SHORT *theX)
{
	INT hdl, i;
	void *p;
	INT *smaller, *larger;
	SHORT *mem_theX;

	hdl = gd->geodata->chunk_alloc_hdl();
	p = gd->geodata->hdl2ptr(hdl);
	smaller = (INT *)p;
	larger = smaller + 1;
	mem_theX = (SHORT *)(larger + 1);
	*smaller = -1;
	*larger = -1;
	for (i = 0; i < gd->nb_X; i++) {
		mem_theX[i] = theX[i];
		}
	if (gd->f_add_verbose) {
		printf("nb_X = %ld hdl = %ld geo added !\n", gd->nb_X, hdl);
		}
	return hdl;
}

INT geo_data_search_and_add(GEO_DATA *gd, SHORT *theX, INT *geo_hdl)
/* return FALSE: geometry already there, not added, 
 * return TRUE: new geometry, added. */
{
	INT last_hdl = 1, new_hdl, ret;
	INT hdl = 1; /* the root */
	void *p;
	INT *smaller, *larger;
	SHORT *mem_theX;
	
	while (hdl > 0) {
		p = gd->geodata->hdl2ptr(hdl);
		smaller = (INT *)p;
		larger = smaller + 1;
		mem_theX = (SHORT *)(larger + 1);
		ret = geo_data_cmp(gd->nb_X, theX, mem_theX);
		if (ret == 0) {
			*geo_hdl = hdl;
			return FALSE;
			}
		if (ret > 0) {
			last_hdl = hdl;
			hdl = *larger;
			}
		else {
			last_hdl = hdl;
			hdl = *smaller;
			}
		}
	
	new_hdl = geo_data_add(gd, theX);
	*geo_hdl = new_hdl;
	p = gd->geodata->hdl2ptr(last_hdl);
	smaller = (INT *)p;
	larger = smaller + 1;
	if (ret > 0)
		*larger = new_hdl;
	else
		*smaller = new_hdl;
	return TRUE;
}

GEO_DATA *init_geo_data(INT nb_incidences)
{
	GEO_DATA *gd;
	
	gd = (GEO_DATA *) my_malloc(sizeof(GEO_DATA), "init_geo_data");
	if (gd == NIL) {
		error("init_geo_data() no memory");
		return NIL;
		}
	gd->nb_X = nb_incidences;
	gd->geodata = NIL;
	gd->chunk_size = 1L << 8;
	gd->nb_cand = 0;
	gd->nb_geo = 0;
	gd->f_second_kind = TRUE;
	gd->f_dont_go_further = FALSE;
	gd->f_dont_go_further_without = FALSE;
	gd->f_transposed = FALSE;
	gd->f_calc_second_tdo = FALSE;
	gd->f_calc_ddp = FALSE;
	gd->f_calc_ddb = FALSE;
	gd->f_tdo_v = FALSE;
	gd->f_tdo_vv = FALSE;
	gd->f_dd_v = FALSE;
	gd->f_dd_vv = FALSE;
	gd->f_canon_v = FALSE;
	gd->f_canon_vv = FALSE;
	gd->f_cand_verbose = FALSE;
	gd->f_add_verbose = FALSE;
	gd->f_print_candidates = FALSE;
	gd->print_candidates_mod = 50;
	gd->f_print_canonicized_candidates = FALSE;
	gd->print_canonicized_candidates_mod = 50;
	gd->f_print_geo = FALSE;
	gd->print_geo_mod = 1;
	gd->f_save_into_geo_file = FALSE;
	gd->f_range = FALSE;
	gd->range_first = 0;
	gd->range_len = 0;
	gd->geo_hdl = -1;
	return gd;
}

void geo_data_print_range(GEO_DATA *gd)
{
	if (gd->f_range) {
		printf("range [%ld-%ld] ", 
			gd->range_first, gd->range_first + gd->range_len - 1);
		}
}

void geo_data_print_header(GEO_DATA *gd, INT f_v)
{
	if (gd == NIL)
		return;
	printf("GEO_DATA, nb_X = %ld:\n", gd->nb_X);
	printf("nb_cand = %ld nb_geo = %ld\n", gd->nb_cand, gd->nb_geo);
	printf("f_second_kind = %ld "
		"f_dont_go_further = %ld "
		"f_dont_go_further_without = %ld "
		"f_transposed = %ld f_calc_second_tdo = %ld "
		"f_calc_ddp = %ld f_calc_ddb = %ld\n", 
		gd->f_second_kind, 
		gd->f_dont_go_further, 
		gd->f_dont_go_further_without, 
		gd->f_transposed, gd->f_calc_second_tdo, 
		gd->f_calc_ddp, gd->f_calc_ddb);
	if (gd->f_range) {
		geo_data_print_range(gd);
		printf("\n");
		}
	if (gd->f_save_into_geo_file) {
		printf("saving into %s\n", gd->geo_fname);
		}
	if (f_v) {
		printf("f_tdo_v = %ld f_tdo_vv = %ld f_dd_v = %ld f_dd_vv = %ld "
			"f_canon_v = %ld f_canon_vv = %ld\n", 
			gd->f_tdo_v, gd->f_tdo_vv, 
			gd->f_dd_v, gd->f_dd_vv, 
			gd->f_canon_v, gd->f_canon_vv);
		printf("f_cand_verbose = %ld\n", gd->f_cand_verbose);
		printf("f_add_verbose = %ld\n", gd->f_add_verbose);
		printf("f_print_candidates = %ld print_candidates_mod = %ld\n", 
			gd->f_print_candidates, gd->print_candidates_mod);
		printf("f_print_canonicized_candidates = %ld "
			"print_canonicized_candidates_mod = %ld\n", 
			gd->f_print_canonicized_candidates, 
			gd->print_canonicized_candidates_mod);
		printf("f_print_geo = %ld print_geo_mod = %ld\n", 
			gd->f_print_geo, gd->print_geo_mod);
		}
}

void geo_data_print(GEO_DATA *gd, INT nrow, INT ncol)
{
	INT hdl_min, hdl_max, hdl;
	
	if (gd == NIL)
		return;
	geo_data_print_header(gd, TRUE);
	if (gd->geodata == NIL) {
		return;
		}
	hdl_min = 1;
	hdl_max = (gd->geodata->times - 1) * gd->geodata->chunk_size + 
				gd->geodata->entries_used - 1;
	
	
#if 0
	printf("in lexicographic order:\n");
	geo_tree_print(gd, 1);
#endif

#if 1
	printf("in the order of generation:\n");
	for (hdl = hdl_min; hdl <= hdl_max; hdl++) {
		geo_tree_node_print(gd, hdl, nrow, ncol);
		}
#endif

}

void geo_tree_print(GEO_DATA *gd, INT hdl, INT nrow, INT ncol)
{
	void *p;
	INT *smaller, *larger;

	if (hdl <= 0)
		return;
	p = gd->geodata->hdl2ptr(hdl);
	smaller = (INT *)p;
	larger = smaller + 1;
	geo_tree_print(gd, *smaller, nrow, ncol);
	geo_tree_node_print(gd, hdl, nrow, ncol);
	geo_tree_print(gd, *larger, nrow, ncol);
}

void geo_tree_node_print(GEO_DATA *gd, INT hdl, INT nrow, INT ncol)
{
	void *p;
	INT *smaller, *larger;
	SHORT *mem_theX;
	INT i;

	p = gd->geodata->hdl2ptr(hdl);
	smaller = (INT *)p;
	larger = smaller + 1;
	mem_theX = (SHORT *)(larger + 1);
	printf("hdl = %ld smaller = %ld larger = %ld theX = ", 
		hdl, *smaller, *larger);
	for (i = 0; i < gd->nb_X; i++) {
		printf("%ld ", (INT) mem_theX[i]);
		}
	printf("\n");
	geo_print(gd, gd->nb_X, mem_theX, nrow, ncol);
}

INT multivalued_geo_tdo(MATRIX_OP M, MATRIX_OP N, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp)
{
	INT *X;
	INT *Val;
	// SHORT Y[GEO_MAX_N << 1];
	// SHORT *theX;
	TDOSS *tdoss;
	INT nrow, ncol;
	INT nb_X, i, j, val;
	
	nrow = M->s_hi();
	ncol = M->s_li();
	X = (INT *) my_malloc(nrow * ncol * sizeof(INT), "multivalued_tdo() X");
	Val = (INT *) my_malloc(nrow * ncol * sizeof(INT), "multivalued_tdo() Val");
	nb_X = 0;
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			val = M->s_iji(i, j);
			if (val == 0)
				continue;
			X[nb_X] = i * ncol + j;
			Val[nb_X] = val;
			nb_X++;
			}
		}
	
	tdoss = calc_ntdo_(nrow, ncol, nb_X, X, 
		TRUE /* f_multivalued */, Val, 
		FALSE /* f_calc_second_tdo */, 
		FALSE /* f_calc_theY */, NIL, 
		FALSE /* f_calc_ddp */, FALSE /* f_calc_ddb */, 
		FALSE /* f_tdo_v */, FALSE /* f_tdo_vv */,
		FALSE /* f_dd_v */, FALSE /* f_dd_vv */, 
		p, q, row_decomp, col_decomp);
	
	M->perm_cols_rows(p, q, N);
#if 0
	theX = the_Y_get_the_X(Y);
	N->m_ilih(ncol, nrow);
	k = 0;
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			N->m_iji(i, j, theX[k++]);
			}
		}
#endif

	tdoss_free(tdoss);
	my_free(X);
	my_free(Val);
	
	return OK;
}

void geo_print(GEO_DATA *gd, INT nb_X, SHORT *theX, INT nrow, INT ncol)
{
#if 0
	INT X[GEO_MAX_N << 1];
	SHORT Y[GEO_MAX_N << 1];
	SHORT *ddp, *ddb;
	SHORT *theY;
	TDOSS *tdoss;
	INT i;
	LABRA_OB aut;
	SYM_OB ago;
	PERMUTATION_OB p, q;
	
	for (i = 0; i < nb_X; i++) {
		X[i] = theX[i];
		}


	tdoss = calc_ntdo(nrow, ncol, nb_X, X, gd->f_calc_second_tdo, 
		TRUE /* f_calc_theY */, Y, 
		gd->f_calc_ddp, gd->f_calc_ddb, 
		gd->f_tdo_v, gd->f_tdo_vv,
		gd->f_dd_v, gd->f_dd_vv);

	ddp = the_Y_get_ddp(Y);
	ddb = the_Y_get_ddb(Y);
	theY = the_Y_get_the_X(Y);
	for (i = 0; i < nb_X; i++) {
		X[i] = theY[i];
		}
	
	geo_Canonicize(FALSE /* f_maxtest */, NIL /* back_to */, 
		nrow, ncol, nb_X, 
		TRUE /* f_print_backtrack_points */, X, tdoss, FALSE /* f_transposed */, 
		gd->f_calc_ddp, ddp, 
		gd->f_calc_ddb, ddb, 
		&p, &q, 
		TRUE /* f_get_aut_group */, &aut, &ago, 
		gd->f_canon_v, gd->f_canon_vv);

	ntdo_the_Y_free(Y); /* eventually free ddp, ddb ! */

	{
	SHORT *theX;

	theX = the_Y_get_the_X(Y);
	printf("\n");
	tdoss_print_theX_short(tdoss, nb_X, theX);
	printf("automorphism group order = ");
	ago.println();
	}
	
	tdoss_free(tdoss);
#endif
}


