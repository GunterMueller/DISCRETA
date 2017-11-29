/* ntdo.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1996
 * (Apr. 23, 1996)
 */

#include <DISCRETA/discreta.h>

#include <DISCRETA/perm.h>
#include <DISCRETA/geo.h>

TDOSS *calc_ntdo_(INT v, INT b, 
	INT nb_X, INT *theX, INT f_multivalued, INT *theVal, 
	INT f_calc_second_tdo, 
	INT f_ddp, VECTOR_OP DDp, 
	INT f_ddb, VECTOR_OP DDb, 
	INT f_v, INT f_vv, INT f_dd_v, INT f_dd_vv, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp)
{
	NTDO_INFO info;
	NTDO *N;
	INT *lambda_0, *mue_0;
	SHORT *ddp = NIL, *ddp_mult = NIL;
	SHORT *ddb = NIL, *ddb_mult = NIL;
	INT ddb_N, ddp_N;
	TDOSS *tdoss = NIL;
	INT i, j, l1, l2;
	NTDO_GRID *Gpoints, *Gblocks;
	
	if (row_decomp && row_decomp->s_obj_k() == VECTOR && row_decomp->s_li()) {
		l1 = row_decomp->s_li();
		}
	else {
		l1 = 1;
		}
	lambda_0 = (INT *) my_malloc(l1 * sizeof(INT), "calc_ntdo_() lambda_0");
	if (l1 > 1) {
		for (i = 0; i < l1; i++) {
			lambda_0[i] = row_decomp->s_ii(i);
			}
		}
	else 
		lambda_0[0] = v;

	if (col_decomp && col_decomp->s_obj_k() == VECTOR && col_decomp->s_li()) {
		l2 = col_decomp->s_li();
		}
	else {
		l2 = 1;
		}
	mue_0 = (INT *) my_malloc(l2 * sizeof(INT), "calc_ntdo_() mue_0");
	if (l2 > 1) {
		for (i = 0; i < l2; i++) {
			mue_0[i] = col_decomp->s_ii(i);
			}
		}
	else 
		mue_0[0] = b;

	info.v = v;
	info.b = b;
	info.nb_X = nb_X;
	info.theX = theX;
	info.f_multivalued = f_multivalued;
	info.theVal = theVal;
	info.llambda = l1;
	info.lmue = l2;
	info.lambda_0 = lambda_0;
	info.mue_0 = mue_0;
	info.f_transposed = FALSE;
	N = ntdo_init(&info, f_vv);
	
	ntdo_calc(N, f_vv);
	
	if (f_v) {
		printf("first TDO:\n");
		ntdo_print(N);
		}

	
	tdos_init_ntdo(&N->tdos, N, N->G, N->G_next, FALSE);

	if (f_v) {
		printf("first TDO-SCHEME:\n");
		tdos_print(&N->tdos, N);
		}
	
	if (f_calc_second_tdo) {
		ntdo_calc_second_tdo(N, f_v, f_vv);
		tdos_free(&N->tdos);
		N->tdos = N->tdos2;
		tdos_nil(&N->tdos2);
		}
	if (N->G->f_points) {
		Gpoints = N->G;
		Gblocks = N->G_next;
		}
	else {
		Gpoints = N->G_next;
		Gblocks = N->G;
		}
	if (row_decomp) {
		row_decomp->m_il(Gpoints->G_max);
		for (i = 0; i < Gpoints->G_max; i++) {
			j = Gpoints->len[i];
			row_decomp->m_ii(i, j);
			}
		}
	if (col_decomp) {
		col_decomp->m_il(Gblocks->G_max);
		for (i = 0; i < Gblocks->G_max; i++) {
			j = Gblocks->len[i];
			col_decomp->m_ii(i, j);
			}
		}
	if (p) {
		if (v != N->p.l){
			error("calc_ntdo_() v != N->p.l");
			return NIL;
			}
		p->m_il(v);
		for (i = 0; i < v; i++) {
			j = N->p.a[i];
			p->m_ii(i, j + 1);
			}
		}
	if (q) {
		if (b != N->q.l) {
			error("calc_ntdo_() b != N->q.l");
			return NIL;
			}
		q->m_il(b);
		for (i = 0; i < b; i++) {
			j = N->q.a[i];
			q->m_ii(i, j + 1);
			}
		}

	printf("calling ntdo_calc_dd()\n"); fflush(stdout);
	if (ntdo_calc_dd(N, 
		f_ddp, f_ddb, 
		&ddp, &ddp_N, &ddp_mult, 
		&ddb, &ddb_N, &ddb_mult, f_dd_v, f_dd_vv) != OK) {
		error("calc_ntdo() error in ntdo_calc_dd()");
		return NIL;
		}
	printf("finished with ntdo_calc_dd()\n"); fflush(stdout);
	if (f_ddp) {
		INT len = (v * (v - 1)) >> 1, i;
		
		DDp->m_il(len);
		for (i = 0; i < len; i++) {
			DDp->m_ii(i, ddp[i]);
			}
		my_free(ddp);
		}
	if (f_ddb) {
		INT len = (b * (b - 1)) >> 1, i;
		
		DDb->m_il(len);
		for (i = 0; i < len; i++) {
			DDb->m_ii(i, ddb[i]);
			}
		my_free(ddb);
		}

	ddp = NIL;
	ddb = NIL;
	
	if (!tdos2tdoss(&N->tdos, ddp_mult, ddb_mult, &tdoss)) {
		error("calc_ntdo() error in tdos2tdoss()");
		return NIL;
		}
	ddp_mult = NIL; /* now stored in tdoss */
	ddb_mult = NIL;

	ntdo_free(N);
	
	my_free(lambda_0);
	my_free(mue_0);
	
	return tdoss;
}

TDOSS *calc_ntdo(INT v, INT b, INT nb_X, INT *theX, 
	INT f_calc_second_tdo, 
	INT f_ddp, VECTOR_OP DDp, 
	INT f_ddb, VECTOR_OP DDb, 
	INT f_v, INT f_vv, INT f_dd_v, INT f_dd_vv)
{
	return calc_ntdo_(v, b, 
		nb_X, theX, FALSE /* f_multivalued */, NIL /* theVal */, 
		f_calc_second_tdo, 
		f_ddp, DDp, 
		f_ddb, DDb, 
		f_v, f_vv, f_dd_v, f_dd_vv, NIL, NIL, NIL, NIL);
}


void ntdo_the_Y_free(SHORT *theY)
{
	SHORT *dd;
	
	dd = the_Y_get_ddp(theY);
	if (dd) {
		my_free(dd);
		the_Y_set_ddp(theY, NIL);
		}
	dd = the_Y_get_ddb(theY);
	if (dd) {
		my_free(dd);
		the_Y_set_ddb(theY, NIL);
		}
}

INT ntdo_calc_theY(NTDO *tdo, SHORT *theY, SHORT *ddp, SHORT *ddb, INT f_v)
{
	INT *M;
	INT i, j, k, ii, jj, v;
	
	M = (INT *) my_malloc(tdo->v * tdo->b * sizeof(INT), "ntdo_calc_theY");
	for (i = 0; i < tdo->v; i++) {
		for (j = 0; j < tdo->b; j++) {
			M[i * tdo->b + j] = 0;
			}
		}
	for (k = 0; k < tdo->nb_X; k++) {
		i = tdo->theXi[k];
		j = tdo->theXj[k];
		ii = tdo->p.a[i];
		jj = tdo->q.a[j];
		if (tdo->f_multivalued) {
			v = tdo->theVal[k];
			M[ii * tdo->b + jj] = v;
			}
		else {
			M[ii * tdo->b + jj] = 1;
			}
		}
	if (f_v) {
		for (i = 0; i < tdo->v; i++) {
			for (j = 0; j < tdo->b; j++) {
				printf(tdo->format_string, M[i * tdo->b + j]);
				}
			printf("\n");
			}
		}
	/* k = 0; */
	k = THE_Y_OFFSET;
	for (i = 0; i < tdo->v; i++) {
		for (j = 0; j < tdo->b; j++) {
			if (tdo->f_multivalued) {
				theY[k] = M[i *tdo->b + j];
				k++;
				}
			else {
				if (M[i * tdo->b + j]) {
					theY[k] = i * tdo->b + j;
					k++;
					}
				}
			}
		}
	the_Y_set_ddp(theY, ddp);
	the_Y_set_ddb(theY, ddb);
	my_free(M);
	return OK;
}

void the_Y_set_ddp(SHORT *theY, SHORT *ddp)
{
	SHORT **they = (SHORT **) theY;
	
	they[0] = ddp;
}

void the_Y_set_ddb(SHORT *theY, SHORT *ddb)
{
	SHORT **they = (SHORT **) theY;
	
	they[1] = ddb;
}

SHORT *the_Y_get_ddp(SHORT *theY)
{
	SHORT **they = (SHORT **) theY;
	
	return they[0];
}

SHORT *the_Y_get_ddb(SHORT *theY)
{
	SHORT **they = (SHORT **) theY;
	
	return they[1];
}

SHORT *the_Y_get_the_X(SHORT *theY)
{
	return theY + THE_Y_OFFSET;
}

void ntdo_print(NTDO *tdo)
{
	NTDO_GRID *Gblocks, *Gpoints;
	INT i, j, ii, jj, first, k;
	INT *vbar;
	INT *hbar;
	INT *M;
	
	
	vbar = (INT *) my_malloc((tdo->b + 1) * sizeof(INT), "ntdo_print() vbar");
	hbar = (INT *) my_malloc((tdo->v + 1) * sizeof(INT), "ntdo_print() hbar");
	M = (INT *) my_malloc(tdo->v * tdo->b * sizeof(INT), "ntdo_print() M");
	if (tdo->G_next->f_points) {
		Gpoints = tdo->G_next;
		Gblocks = tdo->G;
		}
	else {
		Gpoints = tdo->G;
		Gblocks = tdo->G_next;
		}
	for (i = 0; i <= tdo->b; i++)
		vbar[i] = FALSE;
	for (i = 0; i <= tdo->v; i++)
		hbar[i] = FALSE;
	for (i = 0; i < Gpoints->G_max; i++) {
		first = Gpoints->first[i];
		hbar[first] = TRUE;
		}
	hbar[tdo->v] = TRUE;
	for (i = 0; i < Gblocks->G_max; i++) {
		first = Gblocks->first[i];
		vbar[first] = TRUE;
		}
	vbar[tdo->b] = TRUE;
	for (i = 0; i < tdo->v; i++) {
		for (j = 0; j < tdo->b; j++) {
			M[i * tdo->b + j] = 0;
			}
		}
	for (k = 0; k < tdo->nb_X; k++) {
		i = tdo->theXi[k];
		j = tdo->theXj[k];
		ii = tdo->p.a[i];
		jj = tdo->q.a[j];
		M[ii * tdo->b + jj] = 1;
		}
	ntdo_top_border_row(tdo, vbar, hbar);
	for (i = 0; i < tdo->v; i++) {
		if (i > 0 && hbar[i]) {
			ntdo_top_middle_border_row(tdo, vbar, hbar);
			}
		ntdo_print_row(tdo, vbar, hbar, i, M, tdo->b);
		}
	ntdo_top_border_row(tdo, vbar, hbar);
	my_free(M);
	my_free(vbar);
	my_free(hbar);
}

void ntdo_print_row(NTDO *tdo, INT *vbar, INT *hbar, INT i, INT *M, INT dim_n)
{
	INT j;
	
	for (j = 0; j <= tdo->b; j++) {
		if (vbar[j]) {
			if (j == 0 || j == tdo->b)
				printf("#");
			else
				printf("|");
			}
		if (j < tdo->b) {
			if (M[i * dim_n + j]) {
				printf("X");
				}
			else {
				printf(".");
				}
			}
		}
	printf("\n");
}

void ntdo_top_border_row(NTDO *tdo, INT *vbar, INT *hbar)
{
	INT j;
	
	for (j = 0; j <= tdo->b; j++) {
		if (vbar[j])
			printf("#");
		if (j < tdo->b)
			printf("#");
		}
	printf("\n");
}

void ntdo_top_middle_border_row(NTDO *tdo, INT *vbar, INT *hbar)
{
	INT j;
	
	for (j = 0; j <= tdo->b; j++) {
		if (vbar[j]) {
			if (j == 0 || j == tdo->b)
				printf("#");
			else
				printf("+");
			}
		if (j < tdo->b)
			printf("-");
		}
	printf("\n");
}

NTDO *ntdo_init(NTDO_INFO *info, INT f_verbose)
{
	NTDO *N;
	INT x, i, j, k, v, b, nb_X;
	INT val, w;
	
#if 0
	if (info->v > NTDO_MAX_N) {
		error("ntdo_init() v > NTDO_MAX_N, please increase NTDO_MAX_N !");
		return NIL;
		}
	if (info->b > NTDO_MAX_N) {
		error("ntdo_init() b > NTDO_MAX_N, please increase NTDO_MAX_N !");
		return NIL;
		}
#endif
	N = (NTDO *) my_malloc(sizeof(NTDO), "ntdo_init NTDO");
	if (N == NIL) {
		error("ntdo_init(): no memory for N");
		return NIL;
		}
	tdos_nil(&N->tdos);
	tdos_nil(&N->tdos2);
	N->info = info;
	if (info->f_transposed) {
		v = N->v = info->b;
		b = N->b = info->v;
		}
	else {
		v = N->v = info->v;
		b = N->b = info->b;
		}
	nb_X = N->nb_X = info->nb_X;
	N->theXi = (INT *) my_malloc(N->nb_X * sizeof(INT), "ntdo_init theXi");
	if (N->theXi == NIL) {
		error("ntdo_init(): no memory for theXi");
		return NIL;
		}
	N->theXj = (INT *) my_malloc(N->nb_X * sizeof(INT), "ntdo_init theXj");
	if (N->theXj == NIL) {
		error("ntdo_init(): no memory for theXj");
		return NIL;
		}
	N->f_multivalued = info->f_multivalued;
	N->theVal = info->theVal;
	N->max_width = 1;
	if (N->f_multivalued) {
		for (k = 0; k < nb_X; k++) {
			val = N->theVal[k];
			if (val > 1000)
				w = 4;
			else if (val > 100)
				w = 3;
			else if (val > 10)
				w = 2;
			else
				w = 1;
			N->max_width = MAXIMUM(N->max_width, w);
			}
		sprintf(N->format_string, "%%.ld ", N->max_width);
		// sprintf(N->format_string, "%% .%ldld ", N->max_width);
		}
	for (k = 0; k < nb_X; k++) {
		x = info->theX[k];
		j = x % info->b;
		x -= j;
		i = x / info->b;
		if (info->f_transposed) {
			N->theXi[k] = j;
			N->theXj[k] = i;
			}
		else {
			N->theXi[k] = i;
			N->theXj[k] = j;
			}
		}
	if (info->f_transposed) {
		N->llambda = info->lmue;
		N->lmue = info->llambda;
		N->lambda_0 = info->mue_0;
		N->mue_0 = info->lambda_0;
		}
	else {
		N->llambda = info->llambda;
		N->lmue = info->lmue;
		N->lambda_0 = info->lambda_0;
		N->mue_0 = info->mue_0;
		}
	if (f_verbose) {
		ntdo_print_incidence_matrix(N);
		}
	N->max_size = (N->b > N->v) ? N->b : N->v;
	N->G_last = ntdo_grid_init(N->max_size);
	N->G = ntdo_grid_init(N->max_size);
	N->G_next = ntdo_grid_init(N->max_size);
	sp_nil(&N->p);
	sp_nil(&N->pv);
	sp_nil(&N->q);
	sp_nil(&N->qv);
	sp_int(&N->p, N->v, "sp_int ntdo_int");
	sp_int(&N->pv, N->v, "sp_int ntdo_int");
	sp_int(&N->q, N->b, "sp_int ntdo_int");
	sp_int(&N->qv, N->b, "sp_int ntdo_int");
	sp_id(&N->p);
	sp_id(&N->pv);
	sp_id(&N->q);
	sp_id(&N->qv);
	ntdo_grid_init0(N, N->G, N->G_last);
	return N;
}

void ntdo_free(NTDO *N)
{
	if (N) {
		tdos_free(&N->tdos);
		tdos_free(&N->tdos2);
		if (N->theXi) {
			my_free(N->theXi);
			N->theXi = NIL;
			}
		if (N->theXj) {
			my_free(N->theXj);
			N->theXj = NIL;
			}
		if (N->G_last) {
			ntdo_grid_free(N->G_last);
			N->G_last = NIL;
			}
		if (N->G) {
			ntdo_grid_free(N->G);
			N->G = NIL;
			}
		if (N->G_next) {
			ntdo_grid_free(N->G_next);
			N->G_next = NIL;
			}
		sp_free(&N->p);
		sp_free(&N->pv);
		sp_free(&N->q);
		sp_free(&N->qv);
		my_free(N);
		}
}

void ntdo_print_incidence_matrix(NTDO *N)
{
	INT i, j, k, ii, val, nb_X, next_i, next_j;
	
	if (N->f_multivalued) {
		for (k = 0; k < N->nb_X; k++) {
			printf("at %ld,%ld: %ld\n", N->theXi[k], N->theXj[k], N->theVal[k]);
			}
		}
	k = 0;
	nb_X = N->nb_X;
	if (k < nb_X) {
		next_i = N->theXi[k];
		next_j = N->theXj[k];
		}
	else {
		next_i = -1;
		next_j = -1;
		}
	for (i = 0; i < N->v; i++) {
		for (j = 0; j < N->b; j++) {
			if (i == next_i && j == next_j) {
				if (N->f_multivalued) {
					val = N->theVal[k];
					printf(N->format_string, val);
					}
				else
					printf("X");
				if (k < nb_X) {
					k++;
					next_i = N->theXi[k];
					next_j = N->theXj[k];
					}
				else {
					next_i = -1;
					next_j = -1;
					}
				}
			else {
				if (N->f_multivalued) {
					for (ii = 0; ii < N->max_width; ii++)
						printf(" ");
					printf(". ");
					}
				else
					printf(".");
				}
			}
		printf("\n");
		}
}

INT ntdo_calc(NTDO *tdo, INT f_v)
/* before: tdo_calc2 */
{
	NTDO_GRID *G;
	INT steps = 0;
	
	if (f_v) {
		printf("ntdo_calc() step = %ld\n", steps);
		printf("G_last:\n");
		ntdo_grid_print(tdo->G_last);
		printf("G:\n");
		ntdo_grid_print(tdo->G);
		}
	while (TRUE) {
		ntdo_next(tdo, f_v);
		steps++;
		if (f_v) {
			printf("ntdo_calc() step = %ld\n", steps);
			ntdo_grid_print(tdo->G_next);
			}
		if (tdo->G_next->G_max == tdo->G_last->G_max && steps >= 2)
			break; /* this is a TDO */
		G = tdo->G_last;
		tdo->G_last = tdo->G;
		tdo->G = tdo->G_next;
		tdo->G_next = G;
		}
	return TRUE;
}

INT ntdo_next(NTDO *tdo, INT f_v)
/* computes tdo->G_next out of tdo->G_last and tdo->G.
 * The G_last decomposition will eventually be refined. */
{
	NTDO_GRID *G0 = tdo->G;
	NTDO_GRID *G1 = tdo->G_next;
	
	G1->max_size = G0->max_size;
	G1->f_points = !G0->f_points;
	if (G1->f_points)
		G1->m = tdo->v;
	else
		G1->m = tdo->b;
	if (tdo->f_multivalued)
		if (G1->f_points)
			G1->n = tdo->b;
		else
			G1->n = tdo->v;
	else
		G1->n = G0->G_max;
	ntdo_collect_types(tdo, tdo->G, tdo->G_next, f_v);
	ntdo_refine_types(tdo, tdo->G_last, tdo->G_next);
	
	return TRUE;
}

INT ntdo_collect_types(NTDO *tdo, NTDO_GRID *G0, NTDO_GRID *G1, INT f_v)
/* before: tdo_collect_types */
{
	INT i;
	
	ntdo_calc_grid(tdo, G0, G1);
	
	for (i = 0; i < G1->m; i++) {
		G1->type_idx[i] = i;
		}
	G1->G_max = 0;
	/* G1->first[0] = 0; */
	
	if (f_v) {
		printf("in ntdo_collect_types:\n");
		ntdo_grid_print(G1);
		}
	
	return TRUE;
}

INT ntdo_calc_grid(NTDO *tdo, NTDO_GRID *G0, NTDO_GRID *G1)
/* before: tdo_recollect_types */
{
	INT i, j, k, ii, jj, ge;
	INT first, next, v, w, a, aa;
	INT offset;
	
	for (i = 0; i < G1->m; i++)
		for (j = 0; j < G1->n; j++)
			G1->type[i * G1->max_size + j] = 0;
	
	for (k = 0; k < tdo->nb_X; k++) {
		i = tdo->theXi[k];
		j = tdo->theXj[k];
		/* we have the incidence theX[k] = (i,j) */
		ii = tdo->p.a[i];
		jj = tdo->q.a[j];
		/* (i,j) lies at (ii, jj) in the incidence matrix */
		if (!G1->f_points) {
			ge = G0->grid_entry[ii];
			if (tdo->f_multivalued) {
				w = tdo->theVal[k];
				first = G0->first[ge];
				next = first + G0->len[ge];
				offset = jj * G1->max_size;
				for (a = first; a < next; a++) {
					v = G1->type[offset + a];
					if (v < w) {
						if (G1->type[offset + next - 1] != 0) {
							printf("first = %ld next = %ld a = %ld\n", first, next, a);
							printf("v = %ld w = %ld\n", v, w);
							printf("max_size = %ld n = %ld m = %ld\n", G1->max_size, G1->n, G1->m);
							fflush(stdout);
							return error("ntdo_calc_grid() overflow !");
							}
						for (aa = next - 1; aa > a; aa--) {
							G1->type[offset + aa] = 
								G1->type[offset + aa - 1];
							}
						G1->type[offset + a] = w;
						break;
						}
					}
				}
			else
				G1->type[jj * G1->max_size + ge]++;
			}
		else {
			ge = G0->grid_entry[jj];
			if (tdo->f_multivalued) {
				w = tdo->theVal[k];
				first = G0->first[ge];
				next = first + G0->len[ge];
				offset = ii * G1->max_size;
				for (a = first; a < next; a++) {
					v = G1->type[offset + a];
					if (v < w) {
						if (G1->type[offset + next - 1] != 0)
							return error("ntdo_calc_grid() overflow !");
						for (aa = next - 1; aa > a; aa--) {
							G1->type[offset + aa] = 
								G1->type[offset + aa - 1];
							}
						G1->type[offset + a] = w;
						break;
						}
					}
				}
			else
				G1->type[ii * G1->max_size + ge]++;
			}
		}
	return TRUE;
}

INT ntdo_refine_types(NTDO *tdo, NTDO_GRID *Gm1, NTDO_GRID *G1)
{
	INT old_k, old_first, old_len;
	
	/* each part of the previous decomposition gets refined: */
	for (old_k = 0; old_k < Gm1->G_max; old_k++) {
		old_first = Gm1->first[old_k];
		old_len = Gm1->len[old_k];
		ntdo_radix_sort(tdo, G1, 0, old_first, old_first + old_len - 1);
		}
	return TRUE;
}

INT ntdo_radix_sort(NTDO *tdo, NTDO_GRID *G, INT radix, INT first, INT last)
/* radix sort of 
 * G->type[first...last][radix] */
/* the type array will not be changed, only the indices type_idx. 
 * All permutations will be accumulated into p, pv, q, qv. */
{
	INT f_found, idx, k, l, t;
	INT first0, first1, res, k1;
	SPERM *perm, *perm_inv;
	
	if (first == last || radix == G->n) {
		/* [first,...,last] is filled into G. 
		 * grid_entry[first..last] will be set. */
		k = G->G_max;
		for (l = first; l <= last; l++)
			G->grid_entry[l] = k;
		G->first[k] = first;
		G->len[k] = last - first + 1;
		G->G_max++;
		/* printf("ntdo_radix_sort()|new entry: first = %ld len = %ld : ", 
			G->first[k], G->len[k]);
		i1 = G->type_idx[first];
		for (j = 0; j < G->n; j++) {
			printf("%ld ", G->type[i1 * G->max_size + j]);
			}
		printf("\n"); */
		return TRUE;
		}
	for (k = first; k <= last; k++) {
		f_found = ntdo_insert_idx(G, first, k - first, radix, k, &idx);
		if (idx != k) {
			/* apply s = (idx idx+1 ... k). 
			 * if (G->f_points) 
			 *   p := p * s, pv := s^-1 * pv
			 * else
			 *   q := q * s, qv := s^-1 * qv */
			if (G->f_points) {
				perm = &tdo->p;
				perm_inv = &tdo->pv;
				}
			else {
				perm = &tdo->q;
				perm_inv = &tdo->qv;
				}
			sp_mult_apply_forwc_r(perm, idx, k - idx + 1);
			sp_mult_apply_backwc_l(perm_inv, idx, k - idx + 1);
			t = G->type_idx[k];
			for (l = k; l > idx; l--)
				G->type_idx[l] = G->type_idx[l - 1];
			G->type_idx[idx] = t;
			/* grid_entry not yet set */
			}
		}
	first0 = first;
	first1 = G->type_idx[first0];
	for (k = first + 1; k <= last; k++) {
		k1 = G->type_idx[k];
		res = G->type[k1 * G->max_size + radix] - G->type[first1 * G->max_size + radix];
		if (res > 0) {
			error("ntdo_radix_sort(): not descending");
			return FALSE;
			}
		if (res < 0) {
			ntdo_radix_sort(tdo, G, radix + 1, first0, k - 1);
			first0 = k;
			first1 = G->type_idx[first0];
			}
		if (k == last) {
			ntdo_radix_sort(tdo, G, radix + 1, first0, k);
			}
		}
	return (TRUE);
}

INT ntdo_insert_idx(NTDO_GRID *G, INT first, INT len, INT radix, 
	INT search_this, INT *idx)
{
	INT i, st1, cur, cur1, res;
	INT f_found;
	
	st1 = G->type_idx[search_this];
	f_found = FALSE;
	for (i = 0; i < len; i++) {
		cur = first + i;
		cur1 = G->type_idx[cur];
		res = G->type[cur1 * G->max_size + radix] - G->type[st1 * G->max_size + radix];
		if (res == 0)
			f_found = TRUE;
		if (res < 0) {
			*idx = cur;
			return f_found;
			}
		}
	*idx = first + len;
	return f_found;
}



// #include "ntdo_grid.C"
/* ntdo_grid.C */

void ntdo_grid_nil(NTDO_GRID *G)
{
	G->first = NIL;
	G->len = NIL;
	G->type_idx = NIL;
	G->grid_entry = NIL;
	G->type = NIL;
	G->f_perms_allocated = FALSE;
	sp_nil(&G->p);
	sp_nil(&G->pv);
}

NTDO_GRID *ntdo_grid_init(INT max_size)
{
	NTDO_GRID *G;
	
	G = (NTDO_GRID *) my_malloc(sizeof(NTDO_GRID), "ntdo_grid_init NTDO_GRID");
	if (G == NIL) {
		error("ntdo_grid_init() no memory");
		return NIL;
		}
	ntdo_grid_nil(G);
	G->max_size = max_size;
	G->first = (INT *) my_malloc(max_size * sizeof(INT), "ntdo_grid_init ntdo->first");
	G->len = (INT *) my_malloc(max_size * sizeof(INT), "ntdo_grid_init ntdo->len");
	G->type_idx = (INT *) my_malloc(max_size * sizeof(INT), "ntdo_grid_init ntdo->type_idx");
	G->grid_entry = (INT *) my_malloc(max_size * sizeof(INT), "ntdo_grid_init ntdo->grid_entry");
	G->type = (INT *) my_malloc(max_size * max_size * sizeof(INT), "ntdo_grid_init ntdo->type");
	return G;
}

void ntdo_grid_free(NTDO_GRID *G)
{
	if (G) {
		if (G->first) {
			my_free(G->first);
			G->first = NIL;
			}
		if (G->len) {
			my_free(G->len);
			G->len = NIL;
			}
		if (G->type_idx) {
			my_free(G->type_idx);
			G->type_idx = NIL;
			}
		if (G->grid_entry) {
			my_free(G->grid_entry);
			G->grid_entry = NIL;
			}
		if (G->type) {
			my_free(G->type);
			G->type = NIL;
			}
		if (G->f_perms_allocated) {
			sp_free(&G->p);
			sp_free(&G->pv);
			G->f_perms_allocated = FALSE;
			}
		my_free(G);
		}
}

INT ntdo_init_perms(NTDO_GRID *G, INT degree)
{
	sp_int(&G->p, degree, "sp_int ntdo_init_perms");
	sp_int(&G->pv, degree, "sp_int ntdo_init_perms");
	sp_id(&G->p);
	sp_id(&G->pv);
	G->f_perms_allocated = TRUE;
	return OK;
}

INT ntdo_grid_init_partition(NTDO_GRID *G, 
	INT f_points, INT m, INT n, INT len, INT *parts)
{
	INT i, j, ii, f, l;

	G->f_points = f_points;
	G->m = m;
	G->n = n;
	if (m > G->max_size) {
		printf("m = %ld max_size = %ld\n", m, G->max_size);
		return error("ntdo_grid_init_partition(): m > max_size");
		}
	if (n > G->max_size) {
		printf("n = %ld max_size = %ld\n", n, G->max_size);
		return error("ntdo_grid_init_partition(): n > max_size");
		}
	G->G_max = 0;
	f = 0;
	for (i = 0; i < len; i++) {
		l = parts[i];
		G->first[G->G_max] = f;
		G->len[G->G_max] = l;
		for (j = 0; j < G->n; j++)
			G->type[f * G->max_size + j] = 0;
		for (ii = 0; ii < l; ii++) {
			G->type_idx[f + ii] = f;
			G->grid_entry[f + ii] = i;
			}
		G->G_max++;
		f += l;
		}
	return OK;
}

INT ntdo_grid_init0(NTDO *tdo, NTDO_GRID *Gpoints, NTDO_GRID *Gblocks)
{
	if (tdo->f_multivalued) {
		ntdo_grid_init_partition(Gpoints, TRUE /* f_points */, 
			tdo->v /* m */, tdo->b /* n */, 
			tdo->llambda, tdo->lambda_0);
	
		ntdo_grid_init_partition(Gblocks, FALSE /* f_points */, 
			tdo->b /* m */, tdo->v /* n */, 
			tdo->lmue, tdo->mue_0);
		}
	else {
		ntdo_grid_init_partition(Gpoints, TRUE /* f_points */, 
			tdo->v /* m */, tdo->lmue /* n */, 
			tdo->llambda, tdo->lambda_0);
	
		ntdo_grid_init_partition(Gblocks, FALSE /* f_points */, 
			tdo->b /* m */, tdo->llambda /* n */, 
			tdo->lmue, tdo->mue_0);
		}

	return TRUE;
#if 0
	INT i, j, ii, jj, first, len;
	
	/* init Gpoints: */
	Gpoints->f_points = TRUE;
	Gpoints->m = tdo->v;
	Gpoints->n = tdo->lmue;
	Gpoints->G_max = 0;
	first = 0;
	for (i = 0; i < tdo->llambda; i++) {
		len = tdo->lambda_0[i];
		Gpoints->first[Gpoints->G_max] = first;
		Gpoints->len[Gpoints->G_max] = len;
		for (j = 0; j < Gpoints->n; j++)
			Gpoints->type[first * Gpoints->max_size + j] = 0;
		for (ii = 0; ii < len; ii++) {
			Gpoints->type_idx[first + ii] = first;
			Gpoints->grid_entry[first + ii] = i;
			}
		Gpoints->G_max++;
		first += len;
		}

	/* init Gblocks: */
	Gblocks->f_points = FALSE;
	Gblocks->m = tdo->b;
	Gblocks->n = tdo->llambda;
	Gblocks->G_max = 0;
	first = 0;
	for (j = 0; j < tdo->lmue; j++) {
		len = tdo->mue_0[j];
		Gblocks->first[Gblocks->G_max] = first;
		Gblocks->len[Gblocks->G_max] = len;
		for (i = 0; i < Gblocks->n; i++)
			Gblocks->type[first * Gblocks->max_size + i] = 0;
		for (jj = 0; jj < len; jj++) {
			Gblocks->type_idx[first + jj] = first;
			Gblocks->grid_entry[first + jj] = j;
			}
		Gblocks->G_max++;
		first += len;
		}
	return TRUE;
#endif
}

INT ntdo_grid_print(NTDO_GRID *G)
{
	INT i, j, first, len, i1, ii, i0;
	
	if (G->f_points)
		printf("the points:\n");
	else
		printf("the blocks:\n");
	for (i = 0; i < G->G_max; i++) {
		first = G->first[i];
		len = G->len[i];
		printf("at %ld %ld x (", first, len);
		i1 = G->type_idx[first];
		for (j = 0; j < G->n; j++) {
			printf("%ld", G->type[i1 * G->max_size + j]);
			if (j < G->n)
				printf(" ");
			}
		printf(")");
		if (G->f_perms_allocated) {
			printf(" { ");
			for (ii = 0; ii < len; ii++) {
				i1 = first + ii;
				i0 = G->pv.a[i1];
				printf("%ld ", i0);
				}
			printf("}\n");
			}
		else
			printf("\n");
		}
	printf("total: %ld\n", G->m);
	return TRUE;
}

INT ntdo_grid_init_derived_i_first(NTDO *tdo, 
	NTDO_GRID *G, NTDO_GRID *G_old, INT derive_at_i, INT f_v)
/* Erstes Element des grid-entry i 
 * von G_old auszeichnen. 
 * G->G_max wird zu G_old->G_max + 1.
 * Setzt auch type_idx[] und type[][]. */
{
	INT i, i0, j, k, first, len, old_type_idx, ti;
	
	G->max_size = G_old->max_size;
	if (f_v) {
		printf("ntdo_grid_init_derived_i_first()\n");
		}
	G->f_points = G_old->f_points;
	G->m = G_old->m;
	G->n = G_old->n;
	for (i = 0, i0 = 0; i < G_old->G_max; i++, i0++) {
		first = G_old->first[i];
		len = G_old->len[i];
		old_type_idx = G_old->type_idx[first];
		if (old_type_idx < 0 || old_type_idx >= G_old->m) {
			printf("old_type_idx = %ld\n", old_type_idx);
			printf("m = %ld\n", G_old->m);
			error("ntdo_grid_init_derived_i_first() old_type_idx illegal !");
			return FALSE;
			}
		if (f_v) {
			printf("ntdo_grid_init_derived_i_first() i = %ld "
				"first = %ld len = %ld old_type_idx = %ld\n", 
				i, first, len, old_type_idx);
			fflush(stdout);
			}
		ti = first;
		for (k = 0; k < G_old->n; k++) {
			G->type[ti * G->max_size + k] = 
					G_old->type[old_type_idx * G_old->max_size + k];
			}
		if (f_v) {
			printf("ntdo_grid_init_derived_i_first() after copy type[]\n");
			fflush(stdout);
			}
		if (i != derive_at_i) {
			G->first[i0] = first;
			G->len[i0] = len;
			}
		else {
			G->first[i0] = first;
			G->len[i0] = 1;
			G->type_idx[first] = first;
			G->grid_entry[first] = i0;
			len--;
			if (len > 0) {
				i0++;
				first++;
				G->first[i0] = first;
				G->len[i0] = len;
				}
			}
		for (j = 0; j < len; j++) {
			G->type_idx[first + j] = ti;
			G->grid_entry[first + j] = i0;
			}
		if (f_v) {
			printf("ntdo_grid_init_derived_i_first() after entry-fill\n");
			fflush(stdout);
			}
		}
	/* G->first[i0] = G_old->first[i]; */
	G->G_max = i0;
	return TRUE;
}

INT ntdo_grid_init_derived_ij_first(NTDO *tdo, 
	NTDO_GRID *G, NTDO_GRID *G_old, INT I, INT J)
/* Erste Elemente der grid-entrys 
 * I und J von G_old auszeichnen. 
 * G->G_max wird zu G_old->G_max + 1 (bzw. + 2).
 * Setzt auch type_idx[] und type[][]. */
{
	INT i, i0, j, k, first, len;
	INT old_type_idx, ti;

	G->max_size = G_old->max_size;
	G->f_points = G_old->f_points;
	G->m = G_old->m;
	G->n = G_old->n;
	i0 = 0; /* aktueller Index nach G->first[] etc. */
	first = 0;
	for (i = 0; i < G_old->G_max; i++) {
		if (first != G_old->first[i]) {
			return error("ntdo_grid_init_derived_ij_first() first != G_old->first[i]");
			}
		/* first = G_old->first[i]; */
		len = G_old->len[i];
		old_type_idx = G_old->type_idx[first];
		
		if (i == I && i == J) {
			if (len < 2) {
				return error("grid_init_derived_ij_first() i == I && i == J && len < 2");
				}
			/* Ein Paar desselben 
			 * Bereichs auszeichnen: 
			 * gemeinsamer Zweierbereich, 
			 * NICHT einzeln ! 
			 * Grund dafuer: 
			 * Bei einzelner Unterteilung 
			 * waeren die Punkte (Bloecke) 
			 * unterschiedlich behandelt, 
			 * je nachdem wer zuerst liegt. 
			 * Dies muss vermieden werden, 
			 * da die TDO Invariante 
			 * fuer jede Ausgangslage der Inzidenz 
			 * gleich sein muss. */
			G->first[i0] = first;
			G->len[i0] = 2;
			ti = first;
			for (k = 0; k < G_old->n; k++) {
				G->type[ti * G->max_size + k] = 
					G_old->type[old_type_idx * G_old->max_size + k];
				}
			G->type_idx[first + 0] = ti;
			G->type_idx[first + 1] = ti;
			G->grid_entry[first + 0] = i0;
			G->grid_entry[first + 1] = i0;
			first += 2;
			len -= 2;
			i0++;
			}
		else if (i == I || i == J) {
			G->first[i0] = first;
			G->len[i0] = 1;
			ti = first;
			for (k = 0; k < G_old->n; k++) {
				G->type[ti * G->max_size + k] = 
					G_old->type[old_type_idx * G_old->max_size + k];
				}
			G->type_idx[first + 0] = ti;
			G->grid_entry[first + 0] = i0;
			i0++;
			first++;
			len--;
			}
		
		/* Eventuellen Rest verarbeiten 
		 * bzw. i disjunkt von I, J: */
		if (len) {
			G->first[i0] = first;
			G->len[i0] = len;
			ti = first;
			for (k = 0; k < G_old->n; k++) {
				G->type[ti * G->max_size + k] = 
					G_old->type[old_type_idx * G_old->max_size + k];
				}
			for (j = 0; j < len; j++) {
				G->type_idx[first + j] = ti;
				G->grid_entry[first + j] = i0;
				}
			i0++;
			first += len;
			}
		}
	/* G->first[i0] = G_old->first[i]; */
	G->G_max = i0;
	return OK;
}

void ntdo_grid_copy_frame(NTDO_GRID *G1, NTDO_GRID *G2, INT f_v)
{
	INT i, j, ti;

	if (f_v) {
		printf("ntdo_grid_copy_frame()\n");
		printf("G1 = %lx m = %ld n = %ld\n", G1, G1->m, G1->n);
		printf("G2 = %lx m = %ld n = %ld\n", G2, G2->m, G2->n);
		}
	G2->max_size = G1->max_size;
	G2->f_points = G1->f_points;
	G2->m = G1->m;
	G2->n = G1->n;
	for (i = 0; i < G1->G_max; i++) {
		G2->first[i] = G1->first[i];
		G2->len[i] = G1->len[i];
		ti = G1->type_idx[G1->first[i]];
		if (ti < 0 || ti >= G1->m) {
			printf("ti = %ld\n", ti);
			printf("m = %ld\n", G1->m);
			error("ntdo_grid_copy_frame() ti illegal !");
			return;
			}
		for (j = 0; j < G2->n; j++) {
			G2->type[ti * G2->max_size + j] = 
					G1->type[ti * G1->max_size + j];
			}
		}
	/* G2->first[i] = G1->first[i]; */
	G2->G_max = i;
	for (i = 0; i < G2->m; i++) {
		G2->type_idx[i] = G1->type_idx[i];
		G2->grid_entry[i] = G1->grid_entry[i];
		}
}

/*
 *
 */

void frame2ntdo_grid(NTDO_FRAME *frame, NTDO_GRID *grid)
/* kopiert alles ausser type_idx[], 
 * type[][], f_points, m, n. */
{
	INT i, j, first, len;
	
	grid->G_max = frame->G_max;
	for (i = 0; i < frame->G_max; i++) {
		first = frame->first[i];
		len = frame->len[i];
		grid->first[i] = first;
		grid->len[i] = len;
		for (j = 0; j < len; j++) {
			grid->grid_entry[first + j] = i;
			}
		}
	/* grid->first[i] = frame->first[i]; */
}

/*
 *
 */

void tdog_exit(TDO_GRAD *p)
{
	INT i;
	
	if (p->type) {
		my_free(p->type);
		p->type = NIL;
		}
	if (p->tdos) {
		for (i = 0; i < p->nb_tdos; i++) {
			tdos_free(&p->tdos[i]);
			}
		my_free(p->tdos);
		p->tdos = NIL;
		p->nb_tdos = 0;
		}
	if (p->mult) {
		my_free(p->mult);
		p->mult = NIL;
		}
}

void tdog_nil(TDO_GRAD *p)
{
	p->type = NIL;
	p->tdos = NIL;
	p->mult = NIL;
	p->nb_tdos = 0;
}

INT tdog_add_tdos(TDO_GRAD *p, TDO_SCHEME *tdos, INT i)
{
	INT l, j, res;
	
	for (l = 0; l < p->nb_tdos; l++) {
		res = tdos_cmp(tdos, &p->tdos[l]);
		if (res == -99) {
			error("tdog_add_tdos() error in tdos_cmp()");
			printf("l = %ld p->nb_tdos = %ld\n", l, p->nb_tdos);
			return FALSE;
			}
		if (res == 0) {
			p->mult[l]++;
			tdos_free(tdos);
			p->type[i] = l;
			return TRUE;
			}
		if (res < 0) {
			for (j = p->nb_tdos - 1; j >= l; j--) {
				tdos_copy(&p->tdos[j], &p->tdos[j + 1]);
				p->mult[j + 1] = p->mult[j];
				}
			/* alle bereits eingetragenen 
			 * Punkt-Ableitungstypen 
			 * updaten: */
			for (j = 0; j < i; j++) {
				if (p->type[j] >= l)
					p->type[j]++;
				}
			p->type[i] = l;
			p->mult[l] = 1;
			tdos_copy(tdos, &p->tdos[l]);
			tdos_nil(tdos);
			p->nb_tdos++;
			return TRUE;
			}
		}
	p->type[i] = p->nb_tdos;
	tdos_copy(tdos, &p->tdos[p->nb_tdos]);
	tdos_nil(tdos);
	p->mult[p->nb_tdos] = 1;
	p->nb_tdos++;
	return TRUE;
}


// end of ntdo_grid.C


