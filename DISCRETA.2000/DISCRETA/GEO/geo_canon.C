/* geo_canon.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */




#include <DISCRETA/discreta.h>


#include <DISCRETA/geo.h>
#include <DISCRETA/GEO/geo_canonP.h>
#include <DISCRETA/lb.h>

#undef ONLY_ONE


#define PRINT_BACKTRACK_POINTS
#define BACKTRACK_MOD 1000

#define GEO_CHECK_P
#undef DEBUG_GEO_CANON_INIT
#undef DEBUG_GEO_CANON_INIT_MEMORY
#undef DEBUG_GEO_CANON_INIT_DECOMPOSITION
#undef DEBUG_CALC_GRID


INT geo_calc_grid(GEO_CANON *geo, INT k, INT f_v)
/* computes Q[k + 1] according to the 
 * row-decomposition in P[k] 
 * and the column permutation in Q[k]. 
 * computes the types of the elements (columns) [k, \ldots, b-1] 
 * sorts them thereby copying Q[k]->p into Q[k+1]->p 
 * and modifying it there.
 * Remember that initially Q[0] reflects the tdo-decomposition 
 * and Q[0]->p is identity.
 *
 * the permutation which chooses [0, \ldots , k-1] 
 * is already in Q[k + 1]->p 
 * (it is the permutation Q[k] and the transposition which 
 * selects the k-1 element (if k > 0). */
{
	NTDO_GRID *P, *Q, *Q_last, *Q0;
	INT i, j, ti, a, b, ii, jj, ge;
	INT nby, *they, *qtype;
	INT off_type_P, off_dd, type_len;
	INT first, len;
	INT f_col_group;
	
	P = geo->P[k];
	Q0 = geo->Q[0];
	Q_last = geo->Q[k];
	Q = geo->Q[k + 1];
	
#ifdef DEBUG_CALC_GRID
	printf("calc_grid() k = %ld\n", k); fflush(stdout);
#endif
	f_col_group = geo->info->f_col_group;
	if (f_col_group)
		off_type_P = 1; /* 0-th entry: if it lies in the col_group orbit */
	else
		off_type_P = 0;
	off_dd = off_type_P + P->G_max;
	type_len = off_dd;
	if (geo->f_ddb) {
		type_len += k;
		}
	if (type_len > Q->max_size) {
		printf("type_len = %ld Q->max_size = %ld\n", type_len, Q->max_size);
		return error("geo_calc_grid() type_len > Q->max_size");
		}
	Q->m = geo->b;
	Q->n = type_len;

	for (i = 0; i < k; i++)
		Q->type_idx[i] = -1;
	for (i = k; i < Q->m; i++) {
		for (j = 0; j < Q->n; j++) {
			Q->type[i * Q->max_size + j] = 0;
			}
		Q->type_idx[i] = i;
		}

#if 0
	/* we begin each type vector with the grid_entry value of Q_last.
	 * This means that we get a {\em refinement} of Q_last into Q */
	/* this is also achieved by sorting each Q_last part 
	 * on its own ! */
	for (i = k; i < Q->m; i++) {
		ti = Q->type_idx[i];
		if (ti < 0)
			return error("geo_calc_grid() ti < 0 (a)");
		a = Q_last->grid_entry[i];
		Q->type[ti * Q->max_size + 0] = -a;
		}
#endif	

	if (f_col_group) {
		MATRIX_OP TG;
		INT f_in_orbit;

		TG = geo->info->col_transversals;
		for (j = k; j < Q->m; j++) {
			ti = Q->type_idx[j];
			f_in_orbit = TG->s_iji(k, j);
			Q->type[ti * Q->max_size + 0] = f_in_orbit;
			}
		}

#if 0
	for (a = 0; a < geo->nb_X; a++) {
		i = geo->theXi[a];
		j = geo->theXj[a];
		ii = P->p.a[i];
		jj = Q_last->p.a[j];
		if (jj < k)
			continue;
		ge = P->grid_entry[ii];
		ti = Q->type_idx[jj];
		if (ti < 0)
			return error("geo_calc_grid() ti < 0 (b)");
		Q->type[ti * Q->max_size + off_type_P + ge]++;
		}
#endif

	for (jj = k; jj < Q->m; jj++) {
		ti = Q->type_idx[jj];
		if (ti < 0)
			return error("geo_calc_grid() ti[jj] < 0");
		j = Q_last->pv.a[jj];
		nby = geo->nb_y[j];

		// we evaluate some constants here (outside the loop)
		// note that GNU C++ optimizer with -O2 
		// crashes otherwise !
		// gcc version 2.7.2.1
		they = geo->theY;
		they += j * geo->v;
		qtype = Q->type;
		qtype += ti * Q->max_size + off_type_P;
		for (a = 0; a < nby; a++) {
			i = they[a];
			ii = P->p.a[i];
			ge = P->grid_entry[ii];
			qtype[ge]++;
			}
		}
		
	if (geo->f_ddb) {
		// printf("types from ddb: ");
		for (ii = k; ii < Q->m; ii++) {
			ti = Q->type_idx[ii];
			if (ti < 0)
				return error("geo_calc_grid() ti < 0 (c)");
			i = Q_last->pv.a[ii];
			for (jj = 0; jj < k; jj++) {
				j = Q_last->pv.a[jj];
				a = ij2k(j, i, geo->b);
				b = geo->ddb[a];
				// printf("%ld ", b);
				Q->type[ti * Q->max_size + off_dd + jj] = b;
				}
			}
		// printf("\n");
		}
		
	if (geo->info->f_very_verbose) {
		printf("in geo_calc_grid() (k = %ld):\n", k);
		printf("unsorted type vectors (column-wise):\n");
		for (i = k; i < Q->m; i++)
			printf("%2ld ", i);
		printf("\n");
		for (i = k; i < Q->m; i++)
			printf("---");
		printf("\n");
		for (j = 0; j < type_len; j++) {
			for (i = k; i < Q->m; i++) {
				ti = Q->type_idx[i];
				a = Q->type[ti * Q->max_size + j];
				printf("%2ld ", a);
				}
			printf("\n");
			}
		}
		
	sp_mv(&Q_last->p, &Q->p);
	sp_mv(&Q_last->pv, &Q->pv);
	
	if (geo->info->f_very_verbose) {
		printf("Q0:\n");
		print_grid_decomp(Q0);
		printf("Q_last:\n");
		print_grid_decomp(Q_last);
		}
	Q->G_max = 0;
	/* Q->first[0] = k; */
	for (j = 0; j < Q0 /* Q_last */ ->G_max; j++) {
		first = Q0 /* Q_last */ ->first[j];
		len = Q0 /* Q_last */ ->len[j];
		if (k < first + len) {
#ifdef DEBUG_CALC_GRID
			printf("calc_grid() radix_sort(%ld, %ld)\n", k, first + len - 1); fflush(stdout);
#endif
			geo_radix_sort(geo, Q, 0 /* radix */, k, first + len - 1);
			break;
			}
		}
	for (j++; j < Q0 /* Q_last */ ->G_max; j++) {
		first = Q0 /* Q_last */ ->first[j];
		len = Q0 /* Q_last */ ->len[j];
#ifdef DEBUG_CALC_GRID
		printf("calc_grid() radix_sort(%ld, %ld)\n", first, first + len - 1); fflush(stdout);
#endif
		geo_radix_sort(geo, Q, 0 /* radix */, first, first + len - 1);
		}
	
	if (f_v) {
		printf("geo_calc_grid() k = %ld after sorting:\n", k);
		printf("off_type_P=%ld off_dd=%ld type_len=%ld\n", 
			off_type_P, off_dd, type_len);
		ntdo_grid_print(Q);
		}
	return OK;
}

#if 0
INT geo_calc_grid(GEO_CANON *geo, INT k, INT f_v)
/* computes Q[k + 1] according to the 
 * row-decomposition in P[k] 
 * and the column permutation in Q[k]. 
 * computes the types of the elements (columns) [k, \ldots, b-1] 
 * sorts them thereby copying Q[k]->p into Q[k+1]->p 
 * and modifying it there.
 * Remember that initially Q[0] reflects the tdo-decomposition 
 * and Q[0]->p is identity.
 *
 * the permutation which chooses [0, \ldots , k-1] 
 * is already in Q[k + 1]->p 
 * (it is the permutation Q[k] and the transposition which 
 * selects the k-1 element (if k > 0). */
{
	NTDO_GRID *P, *Q, *Q_last, *Q0;
	INT i, j, ti, a, b, ii, jj, ge;
	INT nby, *they, *qtype;
	INT off_type_P, off_dd, type_len;
	INT first, len, last;
	INT f_col_group;
	INT qge;
	
	P = geo->P[k];
	Q0 = geo->Q[0];
	Q_last = geo->Q[k];
	Q = geo->Q[k + 1];
	
#ifdef DEBUG_CALC_GRID
	printf("calc_grid() k = %ld\n", k); fflush(stdout);
#endif
	f_col_group = geo->info->f_col_group;
	if (f_col_group)
		off_type_P = 1; /* 0-th entry: if it lies in the col_group orbit */
	else
		off_type_P = 0;
	off_dd = off_type_P + P->G_max;
	type_len = off_dd;
	if (geo->f_ddb) {
		type_len += k;
		}
	if (type_len > Q->max_size) {
		printf("type_len = %ld Q->max_size = %ld\n", type_len, Q->max_size);
		return error("geo_calc_grid() type_len > Q->max_size");
		}
	Q->m = geo->b;
	Q->n = type_len;
	
	sp_mv(&Q_last->p, &Q->p);
	sp_mv(&Q_last->pv, &Q->pv);
	
	if (geo->info->f_very_verbose) {
		printf("Q0:\n");
		print_grid_decomp(Q0);
		printf("Q_last:\n");
		print_grid_decomp(Q_last);
		}
	Q->G_max = 0;
	qge = Q_last->grid_entry[k];
	first = Q_last->first[qge];
	len = Q_last->len[qge];
	if (geo->info->f_very_verbose) {
		printf("qge = %ld Q_last->first[qge]=%ld Q_last->len[qge]=%ld\n", qge, first, len);
		}
	if (k < first)
		return error("geo_calc_grid() k < first");
	if (k >= first + len)
		return error("geo_calc_grid() k >= first + len");
	last = first + len - 1;
	if (k == last) {
		Q->first[0] = k;
		Q->len[0] = 1;
		Q->grid_entry[k] = 0;
		Q->G_max++;
		for (ge = qge + 1; ge < Q_last->G_max; ge++) {
			first = Q_last->first[ge];
			len = Q_last->len[ge];
			Q->first[Q->G_max] = first;
			Q->len[Q->G_max] = len;
			for (i = first; i < first + len; i++)
				Q->grid_entry[i] = Q->G_max;
			Q->G_max++;
			}
		if (geo->info->f_very_verbose) {
			printf("class of size one, we just copy the remaining classes, Q=\n");
			print_grid_decomp(Q);
			}
		return OK;
		}

	for (i = 0; i < k; i++)
		Q->type_idx[i] = -1;
	for (i = k; i <= last; i++) {
		for (j = 0; j < Q->n; j++) {
			Q->type[i * Q->max_size + j] = 0;
			}
		Q->type_idx[i] = i;
		}
	for (i = last + 1; i < Q->m; i++)
		Q->type_idx[i] = -1;

#if 0
	/* we begin each type vector with the grid_entry value of Q_last.
	 * This means that we get a {\em refinement} of Q_last into Q */
	/* this is also achieved by sorting each Q_last part 
	 * on its own ! */
	for (i = k; i < Q->m; i++) {
		ti = Q->type_idx[i];
		if (ti < 0)
			return error("geo_calc_grid() ti < 0 (a)");
		a = Q_last->grid_entry[i];
		Q->type[ti * Q->max_size + 0] = -a;
		}
#endif	

	if (f_col_group) {
		MATRIX_OP TG;
		INT f_in_orbit;

		TG = geo->info->col_transversals;
		for (j = k; j <= last; j++) {
			ti = Q->type_idx[j];
			f_in_orbit = TG->s_iji(k, j);
			Q->type[ti * Q->max_size + 0] = f_in_orbit;
			}
		}

#if 0
	for (a = 0; a < geo->nb_X; a++) {
		i = geo->theXi[a];
		j = geo->theXj[a];
		ii = P->p.a[i];
		jj = Q_last->p.a[j];
		if (jj < k)
			continue;
		ge = P->grid_entry[ii];
		ti = Q->type_idx[jj];
		if (ti < 0)
			return error("geo_calc_grid() ti < 0 (b)");
		Q->type[ti * Q->max_size + off_type_P + ge]++;
		}
#endif

	for (jj = k; jj <= last; jj++) {
		ti = Q->type_idx[jj];
		if (ti < 0)
			return error("geo_calc_grid() ti[jj] < 0");
		j = Q_last->pv.a[jj];
		nby = geo->nb_y[j];

		// we evaluate some constants here (outside the loop)
		// note that GNU C++ optimizer with -O2 
		// crashes otherwise !
		// gcc version 2.7.2.1
		they = geo->theY;
		they += j * geo->v;
		qtype = Q->type;
		qtype += ti * Q->max_size + off_type_P;
		for (a = 0; a < nby; a++) {
			i = they[a];
			ii = P->p.a[i];
			ge = P->grid_entry[ii];
			qtype[ge]++;
			}
		}
		
	if (geo->f_ddb) {
		for (ii = k; ii <= last; ii++) {
			ti = Q->type_idx[ii];
			if (ti < 0)
				return error("geo_calc_grid() ti < 0 (c)");
			i = Q_last->pv.a[ii];
			for (jj = 0; jj < k; jj++) {
				j = Q_last->pv.a[jj];
				a = ij2k(j, i, geo->b);
				b = geo->ddb[a];
				Q->type[ti * Q->max_size + off_dd + jj] = b;
				}
			}
		}
		
	if (geo->info->f_very_verbose) {
		printf("in geo_calc_grid() (k = %ld last=%ld):\n", k);
		printf("unsorted type vectors (column-wise):\n");
		for (i = k; i <= last; i++)
			printf("%2ld ", i);
		printf("\n");
		for (i = k; i <= last; i++)
			printf("---");
		printf("\n");
		for (j = 0; j < type_len; j++) {
			for (i = k; i < Q->m; i++) {
				ti = Q->type_idx[i];
				a = Q->type[ti * Q->max_size + j];
				printf("%2ld ", a);
				}
			printf("\n");
			}
		}
		
#if 0
	for (j = 0; j < /* Q0 / Q_last */ Q_last ->G_max; j++) {
		first = /* Q0 */ Q_last ->first[j];
		len = /* Q0 */ Q_last  ->len[j];
		if (k < first + len) {
#ifdef DEBUG_CALC_GRID
			printf("calc_grid() radix_sort(%ld, %ld)\n", k, first + len - 1); fflush(stdout);
#endif
			geo_radix_sort(geo, Q, 0 /* radix */, k, first + len - 1);
			break;
			}
		}
	for (j++; j < /* Q0 / Q_last */ Q_last  ->G_max; j++) {
		first = /* Q0 */ Q_last  ->first[j];
		len = /* Q0 */ Q_last  ->len[j];
#ifdef DEBUG_CALC_GRID
		printf("calc_grid() radix_sort(%ld, %ld)\n", first, first + len - 1); fflush(stdout);
#endif
		geo_radix_sort(geo, Q, 0 /* radix */, first, first + len - 1);
		}
#endif
	geo_radix_sort(geo, Q, 0 /* radix */, k, last);
	for (ge = qge + 1; ge < Q_last->G_max; ge++) {
		first = Q_last->first[ge];
		len = Q_last->len[ge];
		Q->first[Q->G_max] = first;
		Q->len[Q->G_max] = len;
		for (i = first; i < first + len; i++)
			Q->grid_entry[i] = Q->G_max;
		Q->G_max++;
		}
	if (geo->info->f_very_verbose) {
		printf("type computed, sorted, refined classes, copied remaining classes, Q=\n");
		print_grid_decomp(Q);
		}
	
	if (f_v) {
		printf("geo_calc_grid() k = %ld after sorting:\n", k);
		printf("off_type_P=%ld off_dd=%ld type_len=%ld\n", 
			off_type_P, off_dd, type_len);
		ntdo_grid_print(Q);
		}
	return OK;
}
#endif

INT geo_radix_sort(GEO_CANON *geo, NTDO_GRID *Q, INT radix, 
	INT first, INT last)
{
	INT f_found, idx, k, l, t;
	INT first0, first1, res, k1;
	SPERM *perm, *perm_inv;

	if (first == last || radix == Q->n) {
		/* sort the interval [first .. last] 
		 * insert as new grid-entry into G
		 * grid_entry[first..last] will be set accordingly. */
		k = Q->G_max;
		/* if (Q->first[k] != first) {
			return error("geo_radix_sort() Q->first[k] != first");
			} */
		for (l = first; l <= last; l++)
			Q->grid_entry[l] = k;
		Q->first[k] = first;
		Q->len[k] = last - first + 1;
		/* Q->first[k + 1] = last + 1; */
		Q->G_max++;
#if 0
		if (iso->info->f_very_verbose) {
			printf("iso2_radix_sort()|new entry: first = %ld len = %ld : ", 
				(INT)first, (INT)G->len[k]);
			i1 = G->type_idx[first];
			for (j = 0; j < G->n; j++) {
				printf("%ld ", (INT)G->type[i1][j]);
				}
			printf("\n");
			}
#endif
		return TRUE;
		}
	for (k = first; k <= last; k++) {
		f_found = geo_insert_idx(Q, first, k - first, radix, k, &idx);
		if (idx != k) {
			/* Apply 
			 * s = (idx idx+1 ... k) 
			 * to the current incidence matrix.
			 *   q := q * s, qv := s^-1 * qv */
			perm = &Q->p;
			perm_inv = &Q->pv;
			sp_mult_apply_forwc_r(perm, idx, k - idx + 1);
			sp_mult_apply_backwc_l(perm_inv, idx, k - idx + 1);
			t = Q->type_idx[k];
			for (l = k; l > idx; l--)
				Q->type_idx[l] = Q->type_idx[l - 1];
			Q->type_idx[idx] = t;
			/* grid_entry not yet set, 
			 * need not be copied yet. */
			}
		}
	first0 = first;
	first1 = Q->type_idx[first0];
	for (k = first + 1; k <= last; k++) {
		k1 = Q->type_idx[k];
		res = Q->type[k1 * Q->max_size + radix] - 
			Q->type[first1 * Q->max_size + radix];
		if (res > 0) {
			return error("geo_radix_sort() not descending");
			}
		if (res < 0) {
			geo_radix_sort(geo, Q, radix + 1, first0, k - 1);
			first0 = k;
			first1 = Q->type_idx[first0];
			}
		if (k == last) {
			geo_radix_sort(geo, Q, radix + 1, first0, k);
			}
		}
	return OK;
}

INT geo_insert_idx(NTDO_GRID *G, INT first, INT len, 
	INT radix, INT search_this, INT *idx)
{
	INT i, st1, cur, cur1, res;
	INT f_found;
	
	st1 = G->type_idx[search_this];
	f_found = FALSE;
	for (i = 0; i < len; i++) {
		cur = first + i;
		cur1 = G->type_idx[cur];
		res = G->type[cur1 * G->max_size + radix] - 
				G->type[st1 * G->max_size + radix];
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

INT geo_choose_first_into_E(GEO_CANON *geo, 
	NTDO_GRID *Grid, NTDO_GRID *Grid_last, INT k, ORDERED_SET *E, INT f_v)
{
	INT i, j, l, first, len, i1, i0, i00;
	INT f_col_group, f_in_transversal;
	MATRIX_OP TG;
		
#if 0
	if (geo->info->f_verbose)
		f_v = TRUE;
#endif
	if (k != Grid->first[0])
		return error("geo_choose_first_into_E(): k != Grid->first[0]");
	f_col_group = geo->info->f_col_group;
	if (f_col_group)
		TG = geo->info->col_transversals;
	first = Grid->first[0];
	len = Grid->len[0];
	E->size = 0;
	for (i = 0; i < len; i++) {
		i1 = first + i;
		i0 = Grid->pv.a[i1];
			/* the original element number */
		i00 = Grid_last->p.a[i0];
			// the position in the Grid_last sorting 
		if (f_col_group) {
			f_in_transversal = TG->s_iji(k, i00);
			if (!f_in_transversal) {
				// if (k < 3)
				//	printf("not in transversal: k = %ld i00 = %ld\n", k, i00);
				continue;
				}
			}
		for (j = 0; j < E->size; j++)
			if (E->a[j] > i0)
				break;
		for (l = E->size - 1; l >= j; l--)
			E->a[l + 1] = E->a[l];
		E->a[j] = i0;
		E->size++;
		if (E->size >= MAX_SET_SIZE)
			return error("geo_choose_first_into_E() E->size >= MAX_SET_SIZE");
		}
#if 0
	if (f_col_group) {
		sp_mv(&Grid_last->p, &Grid->p);
		sp_mv(&Grid_last->pv, &Grid->pv);
		}
#endif
	if (f_v) {
		printf("E_%ld={ ", k);
		for (i = 0; i < E->size; i++) {
			printf("%ld ", E->a[i]);
			}
		printf("}\n");
		}
	return OK;
}

INT geo_add(GEO_CANON *geo, INT k, INT j, INT f_v)
// j is an element of E[level], 
// we apply a permutation such that 
// jj = q[j] lies at k-th position.
/* changes Q[k+1], 
 * computes P[k+1] */
{
	ORDERED_SET *E;
	NTDO_GRID *Q_last, *Q, *P, *P_last;
	INT jj, I, i, i0, i1, first, len, l;

#if 0
	if (geo->info->f_verbose)
		f_v = TRUE;
#endif
	Q = geo->Q[k+1];
	Q_last = geo->Q[k];
	P = geo->P[k+1];
	P_last = geo->P[k];
	jj = Q->p.a[j];
	if (geo->info->f_col_group) {
		PERMUTATION_OB per;
		LABRA_OP G;
		INT a, b;

		jj = Q_last->p.a[j];
		G = geo->info->col_group;
		G->rep_ij(k, jj, &per);
		if (f_v)
			printf("rep %ld->%ld\n", k, jj);
		for (a = 0; a < geo->b; a++) {
			b = per.s_ii(a) - 1;
			geo->tmp_q1.a[b] = a;
			geo->tmp_q2.a[a] = b;
			}
		// now: tmp_q1 maps jj to k, check it:
		if (geo->tmp_q1.a[jj] != k)
			return error("geo_add() geo->tmp_q1.a[jj] != k");
		sp_mult(&Q_last->p, &geo->tmp_q1, &geo->tmp_q3);
		if (geo->tmp_q3.a[j] != k)
			return error("geo_add() geo->tmp_q3.a[j] != k");
		sp_mult(&geo->tmp_q2, &Q_last->pv, &geo->tmp_q4);
		sp_mv(&geo->tmp_q3, &Q->p);
		sp_mv(&geo->tmp_q4, &Q->pv);
		}
	else {
		/* apply the transposition (k jj) to Q->p: */
		sp_mult_apply_tau_r(&Q->p, k, jj);
		sp_mult_apply_tau_l(&Q->pv, k, jj);
		}
	
	sp_mv(&P_last->p, &P->p);
	sp_mv(&P_last->pv, &P->pv);
	if (j != Q->pv.a[k])
		return error("geo_add() j != Q->pv.a[k]");
	P->G_max = 0;
	for (I = 0; I < P_last->G_max; I++) {
		first = P_last->first[I];
		len = P_last->len[I];
		i0 = first;
		for (i = 0; i < len; i++) {
			i1 = first + i;
			if (geo_get_M_ij(geo, P->pv.a[i1], j)) {
			// if (geo->M[P->pv.a[i1] * geo->b + j]) {
				if (i1 != i0) {
					sp_mult_apply_forwc_r(&P->p, i0, i1 - i0 + 1);
					sp_mult_apply_backwc_l(&P->pv, i0, i1 - i0 + 1);
					}
				i0++;
				}
			}
		l = i0 - first;
		if (l > 0) {
			P->first[P->G_max] = first;
			P->len[P->G_max] = l;
			P->G_max++;
			}
		if (l < len) {
			P->first[P->G_max] = first + l;
			P->len[P->G_max] = len - l;
			P->G_max++;
			}
		}
	for (I = 0; I < P->G_max; I++) {
		first = P->first[I];
		len = P->len[I];
		for (i = 0; i < len; i++) {
			P->grid_entry[first + i] = I;
			}
		}
#ifdef GEO_CHECK_P
	geo_check_P(geo, k);
#endif
	if (f_v) {
		printf("geo_add() level = %ld added %ld\n", k, j);
		for (i = 0; i <= k; i++) {
			E = geo->E + i;
			printf("%ld ", E->size);
			}
		printf(" (size) \n");
		for (i = 0; i <= k; i++) {
			printf("%ld ", (INT) Q->pv.a[i]);
			}
		printf("\n");
		geo_print(geo, k);
		fflush(stdout);
		}
	return OK;
}

INT geo_check_P(GEO_CANON *geo, INT k)
{
	NTDO_GRID *P, *Q;
	INT I, i, j, c, first, len;
	
	P = geo->P[k+1];
	Q = geo->Q[k+1];
	j = Q->pv.a[k];
	for (I = 0; I < P->G_max; I++) {
		first = P->first[I];
		len = P->len[I];
		c = geo_get_M_ij(geo, P->pv.a[first], j);
		// c = geo->M[P->pv.a[first] * geo->b + j];
		for (i = 1; i < len; i++) {
			if (geo_get_M_ij(geo, P->pv.a[first + i], j) != c)
			// if (geo->M[P->pv.a[first + i] * geo->b + j] != c)
				return error("geo_check_P() illegal decomposition");
			}
		}
	return OK;
}


// #include "geo_maxtest.C"
// geo_maxtest.C 


INT geo_maxtest_set(INT n, INT k, INT *theX, 
	LABRA_OP G, MATRIX_OP TG, 
	LABRA_OP aut, SYM_OP ago, INT f_v, INT f_vv)
// does not change theX !
{
	GEO_CANON_INFO info;
	INT f_print_backtrack_points = FALSE;
	INT backtrack_points_mod = 1;
	INT res;
	
	info.v = 1;
	info.b = n;
	info.nb_X = k;
	info.theX = theX;
	info.f_transposed = FALSE;
	info.tdoss = NIL;
	info.f_ddp = FALSE;
	info.f_ddb = FALSE;
	info.ddp = NIL;
	info.ddb = NIL;
	info.f_print_backtrack_points = f_print_backtrack_points;
	info.backtrack_points_mod = backtrack_points_mod;
	info.f_verbose = f_v;
	info.f_very_verbose = f_vv;
	info.f_get_aut_group = TRUE;
	info.f_aut_v = FALSE;
	info.f_aut_vv = FALSE;
	info.f_col_group = TRUE;
	info.col_group = G;
	info.col_transversals = TG;
	info.aut = aut;
	info.ago = ago;
	if (G->s_degree_i() != n) {
		error("geo_maxtest_set() wrong degree of G");
		return -2;
		}
	if (TG->s_li() != n) {
		error("geo_maxtest_set() wrong length of TG");
		return -2;
		}
	if (TG->s_hi() != n) {
		error("geo_maxtest_set() wrong height of TG");
		return -2;
		}
	res = geo_maxtest(&info);
	return res;
}

static INT geo_maxtest(GEO_CANON_INFO *info)
{
	GEO_CANON *geo;
	INT res;
	
	geo = geo_canon_init(info);
	geo->f_going_back = FALSE;
	geo->is_first = FALSE; // !!!
		// geo->p0, p0v, q0, q0v are initialized 
		// by geo_canon_init() with the identity permutation.
	geo->first_moved = MAX_GEO_POINTS;
	geo->nb_aut_gens = 0;
	geo->dim_aut_gens = 0;
	geo->aut_gens = NIL;
	geo->auts_first_moved = NIL;
	geo->nb_backtrack = 0;
	geo->nb_backtrack_points = 0;

	/* if (info->f_print_backtrack_points)
		printf("["); */
	geo_maxtest_recursion(geo, 0);
	if (geo->f_going_back) {
		res = geo->back_to;
		// printf("geo_maxtest() not maximal at %ld\n", res);
		}
	else {
		res = -1;
		// printf("geo_maxtest() is maximal !\n");
		if (info->f_get_aut_group) {
			geo_get_aut_group(geo, info->aut, info->ago, 
				info->f_aut_v, info->f_aut_vv);
			}
		}
	if (info->f_print_backtrack_points && 
		geo->nb_backtrack >= geo->info->backtrack_points_mod) {
		printf("nb_backtrack = %ld]", geo->nb_backtrack);
		fflush(stdout);
		}
	
	// geo_canon_apply(geo);

	geo_free_aut_gen(geo);

	geo_canon_free(geo);
	return res;
}

INT geo_maxtest_recursion(GEO_CANON *geo, INT level)
{
	ORDERED_SET *E;
	INT i, j, ret, len, ii, bad_elt;
	
	geo->nb_backtrack++;
	if (geo->info->f_print_backtrack_points) {
		if ((geo->nb_backtrack % geo->info->backtrack_points_mod) ==  0) {
			if (geo->nb_backtrack == geo->info->backtrack_points_mod)
				printf("[");
			printf(".");
			geo->nb_backtrack_points++;
			if (geo->nb_backtrack_points >= 50) {
				printf("\n");
				geo->nb_backtrack_points = 0;
				}
			fflush(stdout);
			}
		}

	if (geo->info->f_very_verbose) {
		printf("geo_maxtest_recursion() level = %ld\n", level);
		fflush(stdout);
		}
	geo_calc_grid(geo, level, geo->info->f_very_verbose);
	E = geo->E + level;
	geo_choose_first_into_E(geo, geo->Q[level + 1], geo->Q[level], level, E, 
		geo->info->f_very_verbose);
	len = E->size;
	for (i = 0; i < len; i++) {
		j = E->a[i];
		geo_add(geo, level, j, geo->info->f_very_verbose);
		if (i > 0 && level < geo->first_moved)
			geo->first_moved = level;
		if (geo->is_first) { // never the case, because we are doing a maxtest !!
			return error("geo_maxtest_recursion(): is_first ");
			if (level == geo->b - 1) {
				geo->is_first = FALSE;
				geo_get_canonical_labelling(geo);
				}
			else {
				geo_maxtest_recursion(geo, level + 1);
				}
			}
		else { /* (!L->is_first) */
			if (level == geo->b - 1) {
				ret = geo_compare(geo, level);
				if (ret == 0) {
					/* we have found an automorphism ! */
					if (geo->info->f_very_verbose) {
						printf("found an automorphism ! in transversal %ld\n", 
							geo->first_moved);
						fflush(stdout);
						}
					if (geo->first_moved != MAX_GEO_POINTS)
						geo_get_aut_gen(geo);
					if (geo->info->f_very_verbose) {
						printf("no = %ld\n", geo->nb_aut_gens);
						fflush(stdout);
						}
					if (level > geo->first_moved)
						return 1;
					}
#if 1
			   	if (ret == 1) {
					if (geo->info->f_verbose) {
						printf("found a better labelling !\n");
						geo_print(geo, level);
						}
					/* we found a better 'canonical' labelling ! */
					if (geo->info->f_very_verbose) {
						printf("found a better labelling !\n");
						}
					geo_get_canonical_labelling(geo);
					geo->first_moved = MAX_GEO_POINTS;
					geo_free_aut_gen(geo);
					geo->f_going_back = TRUE;
					geo->back_to = geo->b - 1;
					return 0;
					}
#endif
				}
			else { // !level == geo->b - 1
				ret = geo_compare(geo, level);
				if (ret > 0) {
					if (geo->info->f_verbose) {
						printf("going back\n");
						geo_print(geo, level);
						}
					geo->f_going_back = TRUE;
					bad_elt = -1;
					for (ii = 0; ii <= level; ii++) {
						bad_elt = MAXIMUM(bad_elt, geo->Q[level + 1]->pv.a[ii]);
						}
					geo->back_to = bad_elt;
					return 0;
					}
				if (ret >= 0) {
					if (geo_maxtest_recursion(geo, level + 1) == 1 && level > geo->first_moved)
						return 1;
					if (geo->f_going_back)
						return 0;
					}
				}
			}
		} /* next i */
	return 0;
}



// #include "geo_canon2.C"
// geo_canon2.C 

INT geo_canonicize_set(INT n, INT k, INT *theX, 
	LABRA_OP G, MATRIX_OP TG, INT f_v, INT f_vv, 
	INT f_get_aut_group, LABRA_OP aut, SYM_OP ago, 
	PERMUTATION_OP transporter)
// changes theX !
// afterwards, theX contains the canonical form.
{
	GEO_CANON_INFO info;
	INT f_print_backtrack_points = FALSE;
	INT backtrack_points_mod = 1;
	PERMUTATION_OB p, q;
	
	info.v = 1;
	info.b = n;
	info.nb_X = k;
	info.theX = theX;
	info.f_transposed = FALSE;
	info.tdoss = NIL;
	info.f_ddp = FALSE;
	info.f_ddb = FALSE;
	info.ddp = NIL;
	info.ddb = NIL;
	info.f_print_backtrack_points = f_print_backtrack_points;
	info.backtrack_points_mod = backtrack_points_mod;
	info.f_verbose = f_v;
	info.f_very_verbose = f_vv;
	info.f_get_aut_group = f_get_aut_group;
	info.aut = aut;
	info.ago = ago;
	info.f_aut_v = FALSE;
	info.f_aut_vv = FALSE;
	info.f_col_group = TRUE;
	info.col_group = G;
	info.col_transversals = TG;
	if (G->s_degree_i() != n)
		return error("geo_canonicize_set() wrong degree of G");
	if (TG->s_li() != n)
		return error("geo_canonicize_set() wrong length of TG");
	if (TG->s_hi() != n)
		return error("geo_canonicize_set() wrong height of TG");
	geo_canonicize(&info, &p, &q);
	q.copy(transporter);
	return OK;
}

INT geo_canon_with_initial_decomposition(INT f_maxtest, INT *back_to, 
	INT nrow, INT ncol, INT nb_X, INT f_print_backtrack_points, INT *theX, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_transposed, INT f_get_aut_group, LABRA_OP aut, 
	INT f_v, INT f_vv)
{
	SYM_OB ago;
	INT r;
	TDOSS *tdoss;
	
	tdoss = tdoss_init_decomposition(row_decomp, col_decomp);
	// printf("geo_canon_simple()\n");
	// fflush(stdout);
	r = geo_Canonicize(f_maxtest, back_to, 
		nrow, ncol, nb_X, f_print_backtrack_points, 
		theX, tdoss, f_transposed, 
		FALSE /* f_ddp */, NIL /* ddp */, 
		FALSE /* f_ddb */, NIL /* ddb */, 
		p, q, 
		f_get_aut_group, aut, &ago, f_v, f_vv);
	// printf("geo_canon_simple() finished with geo_Canonicize()\n");
	// fflush(stdout);
	tdoss_free(tdoss);
	return r;
}
	
INT geo_canon_with_initial_decomposition_and_ddp_ddb(
	INT f_maxtest, INT *back_to, 
	INT nrow, INT ncol, INT nb_X, 
	INT f_print_backtrack_points, INT *theX, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb,
	VECTOR_OP DDp, VECTOR_OP DDb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_transposed, INT f_get_aut_group, LABRA_OP aut, 
	INT f_v, INT f_vv)
{
	INT v = nrow;
	INT b = ncol;
	SYM_OB ago;
	INT r, len, i, j, k, a;
	TDOSS *tdoss;
	SHORT *ddp = NULL;
	SHORT *ddb = NULL;
	
	if (f_ddp) {
		len = (v * (v - 1)) >> 1;
		ddp = (SHORT *) my_malloc(len * sizeof(SHORT), "geo_canon_with_initial_decomposition_and_ddp_ddb() ddp");
		for (i = 0; i < len; i++) {
			ddp[i] = DDp->s_ii(i);
			}
		if (f_vv) {
			printf("DDp:\n");
			for (i = 0; i < v; i++) {
				printf("i=%ld ", i);
				for (j = i + 1; j < v; j++) {
					k = ij2k(i, j, v);
					a = ddp[k];
					printf("j=%ld:%ld ", j, a);
					}
				printf("\n");
				}
			}
		}
	if (f_ddb) {
		len = (b * (b - 1)) >> 1;
		ddb = (SHORT *) my_malloc(len * sizeof(SHORT), "geo_canon_with_initial_decomposition_and_ddp_ddb() ddb");
		for (i = 0; i < len; i++) {
			ddb[i] = DDb->s_ii(i);
			}
		if (f_vv) {
			printf("DDb:\n");
			for (i = 0; i < b; i++) {
				printf("i=%ld ", i);
				for (j = i + 1; j < b; j++) {
					k = ij2k(i, j, b);
					a = ddb[k];
					printf("j=%ld:%ld ", j, a);
					}
				printf("\n");
				}
			}
		}
	tdoss = tdoss_init_decomposition(row_decomp, col_decomp);
	// printf("geo_canon_simple()\n");
	// fflush(stdout);
	r = geo_Canonicize(f_maxtest, back_to, 
		nrow, ncol, nb_X, f_print_backtrack_points, 
		theX, tdoss, f_transposed, 
		f_ddp, ddp, 
		f_ddb, ddb, 
		p, q, 
		f_get_aut_group, aut, &ago, f_v, f_vv);
	// printf("geo_canon_simple() finished with geo_Canonicize()\n");
	// fflush(stdout);
	tdoss_free(tdoss);
	
	if (f_ddp) 
		my_free(ddp);
	if (f_ddb) 
		my_free(ddb);
	
	return r;
}
	
INT geo_canon_simple(INT f_maxtest, INT *back_to,
	INT nrow, INT ncol, INT nb_X, INT f_print_backtrack_points, INT *theX, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_transposed, INT f_get_aut_group, LABRA_OP aut, 
	INT f_v, INT f_vv)
{
	SYM_OB ago;
	INT r;
	
	// printf("geo_canon_simple()\n");
	// fflush(stdout);
	r = geo_Canonicize(f_maxtest, back_to, 
		nrow, ncol, nb_X, f_print_backtrack_points, 
		theX, NIL /* tdoss */, f_transposed, 
		FALSE /* f_ddp */, NIL /* ddp */, 
		FALSE /* f_ddb */, NIL /* ddb */, 
		p, q, 
		f_get_aut_group, aut, &ago, f_v, f_vv);
	// printf("geo_canon_simple() finished with geo_Canonicize()\n");
	// fflush(stdout);
	return r;
}

INT geo_Canonicize(INT f_maxtest, INT *back_to, 
	INT nrow, INT ncol, INT nb_X, INT f_print_backtrack_points, 
	INT *theX, TDOSS *tdoss, INT f_transposed, 
	INT f_ddp, SHORT *ddp, 
	INT f_ddb, SHORT *ddb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	INT f_get_aut_group, LABRA_OP aut, SYM_OP ago, INT f_v, INT f_vv)
// changes theX !
// afterwards, theX contains the canonical form.
// p and q only computed if f_maxtest is FALSE
{
	GEO_CANON_INFO *info;
	INT backtrack_points_mod = 100;
	
	info = (GEO_CANON_INFO *) my_malloc(sizeof(GEO_CANON_INFO), "geo_Canonicize");
	if (info == NIL)
		return error("geo_Canonicize() no memory for info");
	// ddp = the_Y_get_ddp(theY);
	// ddb = the_Y_get_ddb(theY);
	// theX = the_Y_get_the_X(theY);
	info->v = nrow;
	info->b = ncol;
	info->nb_X = nb_X;
	info->theX = theX;
	info->f_transposed = f_transposed;
	info->tdoss = tdoss;
	info->f_ddp = f_ddp;
	info->f_ddb = f_ddb;
	info->ddp = ddp;
	info->ddb = ddb;
	info->f_print_backtrack_points = f_print_backtrack_points;
	info->backtrack_points_mod = backtrack_points_mod;
	info->f_verbose = f_v;
	info->f_very_verbose = f_vv;
	info->f_get_aut_group = f_get_aut_group;
	info->aut = aut;
	info->ago = ago;
	info->f_aut_v = FALSE;
	info->f_aut_vv = FALSE;
	info->f_col_group = FALSE;
	info->col_group = NIL;
	info->col_transversals = NIL;
	if (f_maxtest) {
		INT bt;

		bt = geo_maxtest(info);
		*back_to = bt;
		}
	else {
		geo_canonicize(info, p, q);
		}
	my_free(info);
#if 0
	printf("geo_Canonicize() finished with geo_canonicize/geo_maxtest()\n");
	fflush(stdout);
	printf("p=");
	p->println();
	fflush(stdout);
	printf("q=");
	q->println();
	fflush(stdout);
#endif
	return OK;
}

INT geo_canonicize(GEO_CANON_INFO *info, PERMUTATION_OP p, PERMUTATION_OP q)
{
	GEO_CANON *geo;
	INT i, a;
	
	if (info->v > MAX_GEO_POINTS) 
		return error("geo_canonicize() v > MAX_GEO_POINTS");
	if (info->b > MAX_GEO_POINTS) 
		return error("geo_canonicize() b > MAX_GEO_POINTS");
	geo = geo_canon_init(info);
	geo->f_going_back = FALSE;
	geo->is_first = TRUE;
	geo->first_moved = MAX_GEO_POINTS;
	geo->nb_aut_gens = 0;
	geo->dim_aut_gens = 0;
	geo->aut_gens = NIL;
	geo->auts_first_moved = NIL;
	geo->nb_backtrack = 0;
	geo->nb_backtrack_points = 0;

	/* if (info->f_print_backtrack_points)
		printf("["); */
	geo_canonicize_recursion(geo, 0);
	if (info->f_print_backtrack_points && 
		geo->nb_backtrack >= geo->info->backtrack_points_mod) {
		printf("nb_backtrack = %ld]", geo->nb_backtrack);
		fflush(stdout);
		}
	
	if (info->f_get_aut_group) {
		geo_get_aut_group(geo, info->aut, info->ago, 
			info->f_aut_v, info->f_aut_vv);
		}
	// printf("geo_canonicize() info->v = %ld\n", geo->info->v);
	// printf("geo_canonicize() p0.l = %ld\n", geo->p0.l);
	// printf("geo_canonicize() info->b = %ld\n", geo->info->b);
	// printf("geo_canonicize() q0.l = %ld\n", geo->q0.l);
	p->m_il(geo->info->v);
	q->m_il(geo->info->b);
	if (geo->info->f_transposed) {
		if (geo->info->v != geo->q0.l) {
			printf("geo_canonicize() info->v = %ld\n", geo->info->v);
			printf("geo_canonicize() p0.l = %ld\n", geo->p0.l);
			printf("geo_canonicize() info->b = %ld\n", geo->info->b);
			printf("geo_canonicize() q0.l = %ld\n", geo->q0.l);
			return error("geo_canonicize() geo->info->v != geo->q0.l");
			}
		for (i = 0; i < geo->info->v; i++) {
			a = geo->q0.a[i];
			p->m_ii(i, a + 1);
			}
		if (geo->info->b != geo->p0.l) {
			printf("geo_canonicize() info->v = %ld\n", geo->info->v);
			printf("geo_canonicize() p0.l = %ld\n", geo->p0.l);
			printf("geo_canonicize() info->b = %ld\n", geo->info->b);
			printf("geo_canonicize() q0.l = %ld\n", geo->q0.l);
			return error("geo_canonicize() geo->info->b != geo->p0.l");
			}
		for (i = 0; i < geo->info->b; i++) {
			a = geo->p0.a[i];
			q->m_ii(i, a + 1);
			}
		}
	else {
		for (i = 0; i < geo->info->v; i++) {
			a = geo->p0.a[i];
			p->m_ii(i, a + 1);
			}
		for (i = 0; i < geo->info->b; i++) {
			a = geo->q0.a[i];
			q->m_ii(i, a + 1);
			}
		}
	// printf("geo_canonicize() calling geo_canon_apply()\n");
	// fflush(stdout);
	geo_canon_apply(geo);
	// printf("geo_canonicize() finished with geo_canon_apply()\n");
	// fflush(stdout);

	// printf("geo_canonicize() calling geo_free_aut_gen()\n");
	// fflush(stdout);
	geo_free_aut_gen(geo);
	// printf("geo_canonicize() finished with geo_free_aut_gen()\n");
	// fflush(stdout);

	// printf("geo_canonicize() calling geo_canon_free()\n");
	// fflush(stdout);
	geo_canon_free(geo);
	// printf("geo_canonicize() finished with geo_canon_free()\n");
	// fflush(stdout);
	
	return OK;
}

INT geo_canonicize_recursion(GEO_CANON *geo, INT level)
{
	ORDERED_SET *E;
	INT i, j, ret, len;
	
	geo->nb_backtrack++;
	if (geo->info->f_print_backtrack_points) {
		if ((geo->nb_backtrack % geo->info->backtrack_points_mod) ==  0) {
			if (geo->nb_backtrack == geo->info->backtrack_points_mod)
				printf("[");
			printf(".");
			geo->nb_backtrack_points++;
			if (geo->nb_backtrack_points >= 50) {
				printf("\n");
				geo->nb_backtrack_points = 0;
				}
			fflush(stdout);
			}
		}

	if (geo->info->f_very_verbose) {
		printf("geo_canonicize_recursion() level = %ld\n", level);
		fflush(stdout);
		}
	geo_calc_grid(geo, level, geo->info->f_very_verbose);
	E = geo->E + level;
	geo_choose_first_into_E(geo, geo->Q[level + 1], geo->Q[level], level, E, 
		geo->info->f_very_verbose);
	len = E->size;
	for (i = 0; i < len; i++) {
		j = E->a[i];
		geo_add(geo, level, j, geo->info->f_very_verbose);
		if (i > 0 && level < geo->first_moved)
			geo->first_moved = level;
		if (geo->is_first) {
			if (level == geo->b - 1) {
				geo->is_first = FALSE;
				geo_get_canonical_labelling(geo);
				}
			else {
				geo_canonicize_recursion(geo, level + 1);
				}
			}
		else { /* (!L->is_first) */
#ifdef ONLY_ONE
			return 0;
#endif
			if (level == geo->b - 1) {
				ret = geo_compare(geo, level);
				if (ret == 0) {
					/* we have found an automorphism ! */
					if (geo->info->f_very_verbose) {
						printf("found an automorphism ! in transversal %ld\n", 
							geo->first_moved);
						fflush(stdout);
						}
					geo_get_aut_gen(geo);
					if (geo->info->f_very_verbose) {
						printf("no = %ld\n", geo->nb_aut_gens);
						fflush(stdout);
						}
					if (level > geo->first_moved)
						return 1;
					}
				if (ret == 1) {
					/* we found a better 'canonical' labelling ! */
					if (geo->info->f_very_verbose) {
						printf("found a better labelling !\n");
						}
					geo_get_canonical_labelling(geo);
					geo->first_moved = MAX_GEO_POINTS;
					geo_free_aut_gen(geo);
					}
				}
			else {
				ret = geo_compare(geo, level);
				if (ret >= 0) {
					if (geo_canonicize_recursion(geo, level + 1) == 1 && level > geo->first_moved)
						return 1;
					}
				}
			}
		} /* next i */
	return 0;
}


// #include "geo_util.C"
// geo_util.C

GEO_CANON *geo_canon_init(GEO_CANON_INFO *info)
{
	GEO_CANON *GC;
	INT i, j, k, a;

	GC = (GEO_CANON *) my_malloc(sizeof(GEO_CANON), "geo_canon_init GEO_CANON *");
	if (GC == NIL) {
		error("geo_canon_init() no memory");
		return NIL;
		}
	GC->info = info;
	GC->nb_X = info->nb_X;
	if (info->f_transposed) {
		GC->v = info->b;
		GC->b = info->v;
		GC->f_ddp = info->f_ddb;
		GC->f_ddb = info->f_ddp;
		GC->ddp = info->ddb;
		GC->ddb = info->ddp;
		}
	else {
		GC->v = info->v;
		GC->b = info->b;
		GC->f_ddp = info->f_ddp;
		GC->f_ddb = info->f_ddb;
		GC->ddp = info->ddp;
		GC->ddb = info->ddb;
		}
	if (info->f_very_verbose) {
		printf("geo_canon_init(): f_transposed = %ld v = %ld b = %ld "
			"f_ddp = %ld f_ddb = %ld\n", 
			info->f_transposed, GC->v, GC->b, GC->f_ddp, GC->f_ddb);
		}
	GC->max_size = (GC->b > GC->v) ? GC->b : GC->v;
	if (GC->f_ddb)
		GC->max_size <<= 1;
#ifdef DEBUG_GEO_CANON_INIT
	printf("geo_canon_init() max_size = %ld\n", GC->max_size); fflush(stdout);
#endif
#ifdef DEBUG_GEO_CANON_INIT_MEMORY
	printf("geo_canon_init() allocating theXi theXj\n"); fflush(stdout);
#endif
	GC->theXi = (INT *) my_malloc(GC->nb_X * sizeof(INT), "geo_canon_init theXi");
	GC->theXj = (INT *) my_malloc(GC->nb_X * sizeof(INT), "geo_canon_init theXj");
#ifdef DEBUG_GEO_CANON_INIT_MEMORY
	printf("geo_canon_init() theXi theXj allocated\n"); fflush(stdout);
#endif
	if (GC->theXi == NIL || GC->theXj == NIL) {
		error("geo_canon_init() no memory for theXi/theXj");
		return NIL;
		}
#ifdef DEBUG_GEO_CANON_INIT
	for (k = 0; k < GC->nb_X; k++) {
		printf("%ld ", (INT) info->theX[k]);
		}
	printf("\n");
#endif
	for (k = 0; k < GC->nb_X; k++) {
		a = info->theX[k];
		if (a < 0) {
			error("geo_canon_init() a < 0");
			return NIL;
			}
		if (a >= GC->v * GC->b) {
			error("geo_canon_init() a >= GC->v * GC->b");
			return NIL;
			}
		j = a % info->b;
		a -= j;
		i = a / info->b;
		if (info->f_transposed) {
			GC->theXi[k] = j;
			GC->theXj[k] = i;
			}
		else {
			GC->theXi[k] = i;
			GC->theXj[k] = j;
			}
#ifdef DEBUG_GEO_CANON_INIT
		printf("geo_canon_init() (%ld,%ld)\n", GC->theXi[k], GC->theXj[k]); fflush(stdout);
#endif
		}

	GC->theX = (INT *) my_malloc(sizeof(INT) * GC->v * GC->b, "geo_canon_init theX");
	GC->theY = (INT *) my_malloc(sizeof(INT) * GC->b * GC->v, "geo_canon_init theX");
	GC->nb_x = (INT *) my_malloc(sizeof(INT) * GC->v, "geo_canon_init nb_x");
	GC->nb_y = (INT *) my_malloc(sizeof(INT) * GC->b, "geo_canon_init nb_y");
	if (GC->theX == NIL || GC->theY == NIL || GC->nb_x == NIL || GC->nb_y == NIL) {
		error("geo_canon_init() no memory for theX/theY/nb_x/nb_y");
		return NIL;
		}
	for (i = 0; i < GC->v; i++)
		GC->nb_x[i] = 0;
	for (j = 0; j < GC->b; j++)
		GC->nb_y[j] = 0;
	for (k = 0; k < GC->nb_X; k++) {
		i = GC->theXi[k];
		j = GC->theXj[k];
		GC->theX[i * GC->b + GC->nb_x[i]++] = j;
		GC->theY[j * GC->v + GC->nb_y[j]++] = i;
		}
	if (info->f_very_verbose) {
		printf("i / nb_x / theX[]:\n");
		for (i = 0; i < GC->v; i++) {
			printf("%ld %ld ", i, GC->nb_x[i]);
			for (j = 0; j < GC->nb_x[i]; j++) {
				printf("%ld ", GC->theX[i * GC->b + j]);
				}
			printf("\n");
			}
		printf("j / nb_y / theY[]:\n");
		for (j = 0; j < GC->b; j++) {
			printf("%ld %ld ", j, GC->nb_y[j]);
			for (i = 0; i < GC->nb_y[j]; i++) {
				printf("%ld ", GC->theY[j * GC->v + i]);
				}
			printf("\n");
			}
		}


#ifdef GEO_CANON_ALLOC_M
	GC->M = (INT *) my_malloc(GC->v * GC->b * sizeof(INT), "geo_canon_init M");
	if (GC->M == NIL) {
		error("geo_canon_init() no memory for M");
		return NIL;
		}
	for (i = 0; i < GC->v * GC->b; i++)
		GC->M[i] = 0;
	for (k = 0; k < GC->nb_X; k++) {
		i = GC->theXi[k];
		j = GC->theXj[k];
		GC->M[i * GC->b + j] = 1;
		}
#ifdef DEBUG_GEO_CANON_INIT
	printf("geo_canon_init() M initialized\n"); fflush(stdout);
#endif
	if (info->f_very_verbose) {
		printf("M:\n");
		for (i = 0; i < GC->v; i++) {
			for (j = 0; j < GC->b; j++) {
				if (GC->M[i * GC->b + j])
					printf("X");
				else
					printf(".");
				}
			printf("\n");
			}
		}
#endif

	if (info->tdoss) {
		if (info->f_transposed) {
			GC->tdo_m = tdoss_n(info->tdoss);
			GC->tdo_n = tdoss_m(info->tdoss);
			}
		else {
			GC->tdo_m = tdoss_m(info->tdoss);
			GC->tdo_n = tdoss_n(info->tdoss);
			}
		GC->tdo_V = (INT *) my_malloc(GC->tdo_m * sizeof(INT), "geo_canon_init tdo_V");
		GC->tdo_B = (INT *) my_malloc(GC->tdo_n * sizeof(INT), "geo_canon_init tdo_B");
		if (GC->tdo_V == NIL || GC->tdo_B == NIL) {
			error("geo_canon_init() no memory for tdo_V/tdo_B");
			return NIL;
			}
		if (info->f_transposed) {
			for (i = 0; i < GC->tdo_m; i++)
				GC->tdo_V[i] = tdoss_Bj(info->tdoss, i);
			for (j = 0; j < GC->tdo_n; j++)
				GC->tdo_B[j] = tdoss_Vi(info->tdoss, j);
			}
		else {
			for (i = 0; i < GC->tdo_m; i++)
				GC->tdo_V[i] = tdoss_Vi(info->tdoss, i);
			for (j = 0; j < GC->tdo_n; j++)
				GC->tdo_B[j] = tdoss_Bj(info->tdoss, j);
			}
		}
	else {
		GC->tdo_m = 1;
		GC->tdo_n = 1;
		GC->tdo_V = (INT *) my_malloc(GC->tdo_m * sizeof(INT), "geo_canon_init tdo_V");
		GC->tdo_B = (INT *) my_malloc(GC->tdo_n * sizeof(INT), "geo_canon_init tdo_B");
		if (GC->tdo_V == NIL || GC->tdo_B == NIL) {
			error("geo_canon_init() no memory for tdo_V/tdo_B");
			return NIL;
			}
		GC->tdo_V[0] = GC->v;
		GC->tdo_B[0] = GC->b;
		}
#ifdef DEBUG_GEO_CANON_INIT
	printf("geo_canon_init() tdo initialized\n"); fflush(stdout);
#endif
	if (info->f_very_verbose) {
		printf("geo_canon_init() (after transpose) tdo_m = %ld tdo_n = %ld\n", 
			GC->tdo_m, GC->tdo_n);
		printf("V: ");
		for (i = 0; i < GC->tdo_m; i++) {
			printf("%ld ", GC->tdo_V[i]);
			}
		printf("\n");
		printf("B: ");
		for (i = 0; i < GC->tdo_n; i++) {
			printf("%ld ", GC->tdo_B[i]);
			}
		printf("\n");
		}
#ifdef DEBUG_GEO_CANON_INIT_MEMORY
	printf("geo_canon_init() allocating E\n"); fflush(stdout);
#endif
	GC->E = (ORDERED_SET *) my_malloc((GC->b + 1) * sizeof(ORDERED_SET), 
		"geo_canon_init E");
#ifdef DEBUG_GEO_CANON_INIT_MEMORY
	printf("geo_canon_init() E allocated\n"); fflush(stdout);
#endif
	if (GC->E == NIL) {
		error("geo_canon_init() no memory for E");
		return NIL;
		}
	
	/*
	 * initialization of decompositions:
	 */
	GC->P = (NTDO_GRID **) my_malloc((GC->max_size + 1) * sizeof(NTDO_GRID *), 
		"geo_canon_init P");
	if (GC->P == NIL) {
		error("geo_canon_init() no memory for P");
		return NIL;
		}
	for (i = 0; i <= GC->max_size; i++) {
		GC->P[i] = ntdo_grid_init(GC->max_size);
		ntdo_init_perms(GC->P[i], GC->v);
		}
	GC->nb_P = GC->max_size + 1;
	
	GC->Q = (NTDO_GRID **) my_malloc((GC->max_size + 1) * sizeof(NTDO_GRID *), 
		"geo_canon_init Q");
	if (GC->Q == NIL) {
		error("geo_canon_init() no memory for Q");
		return NIL;
		}
	for (i = 0; i <= GC->max_size; i++) {
		GC->Q[i] = ntdo_grid_init(GC->max_size);
		ntdo_init_perms(GC->Q[i], GC->b);
		}
	GC->nb_Q = GC->max_size + 1;
	
	/* initial point/block decomposition 
	 * according to tdo: */
	ntdo_grid_init_partition(GC->P[0], TRUE /* f_points */, 
		GC->v /* m */, GC->max_size /* n */, 
		GC->tdo_m, GC->tdo_V);
	if (info->f_very_verbose) {
		printf("P[0]:\n");
		for (i = 0; i < GC->P[0]->G_max; i++) {
			printf("%ld/%ld ", GC->P[0]->first[i], GC->P[0]->len[i]);
			}
		printf("\n");
		for (i = 0; i < GC->v; i++) {
			printf("%ld ", GC->P[0]->grid_entry[i]);
			}
		printf("\n");
		}
	
	ntdo_grid_init_partition(GC->Q[0], FALSE /* f_points */, 
		GC->b /* m */, GC->max_size /* n */, 
		GC->tdo_n, GC->tdo_B);
	if (info->f_very_verbose) {
		printf("Q[0]:\n");
		for (i = 0; i < GC->Q[0]->G_max; i++) {
			printf("%ld/%ld ", GC->Q[0]->first[i], GC->Q[0]->len[i]);
			}
		printf("\n");
		for (i = 0; i < GC->b; i++) {
			printf("%ld ", GC->Q[0]->grid_entry[i]);
			}
		printf("\n");
		}
	
#ifdef DEBUG_GEO_CANON_INIT
	printf("geo_canon_init() decompositions initialized\n"); fflush(stdout);
#endif
	GC->is_first = TRUE;
	sp_nil(&GC->p0);
	sp_nil(&GC->p0v);
	sp_nil(&GC->q0);
	sp_nil(&GC->q0v);
	sp_nil(&GC->tmp_q1);
	sp_nil(&GC->tmp_q2);
	sp_nil(&GC->tmp_q3);
	sp_nil(&GC->tmp_q4);
	sp_int(&GC->p0, GC->v, "sp_int geo_canon_init");
	sp_int(&GC->p0v, GC->v, "sp_int geo_canon_init");
	sp_int(&GC->q0, GC->b, "sp_int geo_canon_init");
	sp_int(&GC->q0v, GC->b, "sp_int geo_canon_init");
	sp_int(&GC->tmp_q1, GC->b, "sp_int geo_canon_init");
	sp_int(&GC->tmp_q2, GC->b, "sp_int geo_canon_init");
	sp_int(&GC->tmp_q3, GC->b, "sp_int geo_canon_init");
	sp_int(&GC->tmp_q4, GC->b, "sp_int geo_canon_init");
	sp_id(&GC->p0);
	sp_id(&GC->p0v);
	sp_id(&GC->q0);
	sp_id(&GC->q0v);
	sp_id(&GC->tmp_q1);
	sp_id(&GC->tmp_q2);
	sp_id(&GC->tmp_q3);
	sp_id(&GC->tmp_q4);
	return GC;
}

void geo_canon_free(GEO_CANON *GC)
{
	if (GC == NIL)
		return;
	if (GC->theXi) {
		my_free(GC->theXi);
		GC->theXi = NIL;
		}
	if (GC->theXj) {
		my_free(GC->theXj);
		GC->theXj = NIL;
		}
#ifdef GEO_CANON_ALLOC_M
	if (GC->M) {
		my_free(GC->M);
		GC->M = NIL;
		}
#endif
	my_free(GC->theX);
	my_free(GC->theY);
	my_free(GC->nb_x);
	my_free(GC->nb_y);
	if (GC->tdo_V) {
		my_free(GC->tdo_V);
		GC->tdo_V = NIL;
		}
	if (GC->tdo_B) {
		my_free(GC->tdo_B);
		GC->tdo_B = NIL;
		}
	if (GC->E) {
		my_free(GC->E);
		GC->E = NIL;
		}
	geo_canon_free_PQ(GC);
	if (GC->P) {
		my_free(GC->P);
		GC->P = NIL;
		}
	if (GC->Q) {
		my_free(GC->Q);
		GC->Q = NIL;
		}
	sp_free(&GC->p0);
	sp_free(&GC->p0v);
	sp_free(&GC->q0);
	sp_free(&GC->q0v);
	sp_free(&GC->tmp_q1);
	sp_free(&GC->tmp_q2);
	sp_free(&GC->tmp_q3);
	sp_free(&GC->tmp_q4);
	my_free(GC);
}

void geo_canon_free_PQ(GEO_CANON *GC)
{
	INT i;

	for (i = 0; i < GC->nb_P; i++) {
		ntdo_grid_free(GC->P[i]);
		GC->P[i] = NIL;
		}
	GC->nb_P = 0;
	for (i = 0; i < GC->nb_Q; i++) {
		ntdo_grid_free(GC->Q[i]);
		GC->Q[i] = NIL;
		}
	GC->nb_Q = 0;
}

INT geo_print(GEO_CANON *geo, INT k)
{
	INT i, j, ii, jj, kk, first, v, b;
	INT vbar[MAX_GEO_POINTS + 1];
	INT hbar[MAX_GEO_POINTS + 1];
	INT *M;
	BYTE **str;
	NTDO_GRID *Q, *P;
	
	M = (INT *) my_malloc(geo->v * geo->b * sizeof(INT), "geo_print() M");
	if (M == NIL)
		return error("geo_print() no memory for M");
	Q = geo->Q[k+1];
	P = geo->P[k+1];
	for (i = 0; i <= MAX_GEO_POINTS; i++) {
		vbar[i] = FALSE;
		hbar[i] = FALSE;
		}
	first = 0;
	for (i = 0; i < P->G_max; i++) {
		hbar[first] = TRUE;
		first += P->len[i];
		}
	v = first;
	hbar[first] = TRUE;
	first = 0;
	for (i = 0; i < geo->tdo_n; i++) {
		vbar[first] = TRUE;
		first += geo->tdo_B[i];
		}
	b = k + 1;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			M[i * geo->b + j] = 0;
			}
		}
	for (kk = 0; kk < geo->nb_X; kk++) {
		i = geo->theXi[kk];
		j = geo->theXj[kk];
		ii = P->p.a[i];
		jj = Q->p.a[j];
		M[ii * geo->b + jj] = 1;
		}
	grid_alloc(&str, 3 * MAX_GEO_POINTS, 3 * MAX_GEO_POINTS);
	grid_prepare(str, v, b, geo->b, (INT *) M, vbar, hbar);
	grid_print(v, hbar, str);
	grid_free(str, 3 * MAX_GEO_POINTS);
	my_free(M);
	return OK;
}

INT geo_get_aut_gens(GEO_CANON *geo, 
	VECTOR_OP auts, VECTOR_OP first_moved)
{
	SPERM *p, *q;
	PERMUTATION_OP perm;
	INT i, n, nb_gen, ii, k;

	n = geo->b;
	nb_gen = geo->nb_aut_gens;
	auts->m_il(nb_gen);
	first_moved->m_il(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		p = geo->aut_gens + i;
		q = p;
		k = geo->auts_first_moved[i];
		first_moved->m_ii(i, k);
		perm = (PERMUTATION_OP) auts->s_i(i);
		perm->m_il(n);
		for (ii = 0; ii < n; ii++) {
			k = q->a[ii];
			perm->m_ii(ii, k + 1);
			}
		}
	return OK;
}

INT geo_get_aut_group(GEO_CANON *geo, 
	LABRA_OP Aut, SYM_OP ago, INT f_v, INT f_vv)
{
	SPERM *p, *q;
	SPERM tmp1, tmp2;
	PERMUTATION_OP KM_p;
	INT n, i, nb_gen, j, k, fst_moved;
	INT ii, jj;
	// INT f_v = TRUE, f_vv = FALSE;
	FILE *fp_txt;
	BYTE str[1024];

	// printf("in geo_get_aut_group()\n");
	// fflush(stdout);
	
	fp_txt = stdout;
	n = geo->b;
	Aut->Init_no_generators(n, f_v, f_vv);
	nb_gen = geo->nb_aut_gens;
	sp_nil(&tmp1);
	sp_nil(&tmp2);
	sp_int(&tmp1, n, "sp_int geo_get_aut_group");
	sp_int(&tmp2, n, "sp_int geo_get_aut_group");
	for (i = 0; i < nb_gen; i++) {
		p = geo->aut_gens + i;
		q = p;
		sp_mult(&geo->q0v, p, &tmp1);
		sp_mult(&tmp1, &geo->q0, &tmp2);
		q = &tmp2;
		fst_moved = geo->auts_first_moved[i];
		j = fst_moved /* L->go[fst_moved] */;
			/* the automorphism moves j_ */
		/* j = L->p0v.a[j]; */
		k = q->a[j];
		if (f_vv) {
			printf("i = %ld: fst_moved = %ld j = %ld k = %ld\n", 
				i, fst_moved, j, k);
			sp_print(p);
			fflush(stdout);
			}
		if (Aut->s_V_ii(k) == k) {
			Aut->s_V_i(k)->m_i(j);
			KM_p = Aut->s_KM_i(k);
			for (ii = 0; ii < n; ii++) {
				jj = q->a[ii];
				KM_p->m_ii(ii, jj + 1);
				}
			}
		else {
			if (f_vv) {
				printf("skipped!\n");
				fflush(stdout);
				}
			}
		} /* next i */
	Aut->update_EM();
	Aut->calc_T();
	Aut->group_order(ago);
	if (f_v) {
		fprintf(fp_txt, "ago = ");
		fflush(fp_txt);
		str[0] = 0;
		ago->sprint(str);
		fprintf(fp_txt, "%s\n", str);
		/* ago->println(); */
		/* Aut->print_group_order(); */
		fflush(fp_txt);
		}
	sp_free(&tmp1);
	sp_free(&tmp2);
	
	// printf("finished with geo_get_aut_group()\n");
	// fflush(stdout);
	return OK;
}

#define REALLOC_OVERSIZE 4

INT geo_get_aut_gen(GEO_CANON *geo)
{
	INT new_dim, i, n;
	SPERM *new_auts, *p;
	INT *new_first;

	if (geo->nb_aut_gens >= geo->dim_aut_gens) {
		/* realloc: */
		new_dim = geo->dim_aut_gens + REALLOC_OVERSIZE;
		new_auts = (SPERM *) my_malloc(new_dim * sizeof(SPERM), 
			"geo_get_aut_gen new_auts");
		if (new_auts == NIL)
			return error("geo_get_aut_gen(): no memory for new_auts");
		new_first = (INT *) my_malloc(new_dim * sizeof(INT), 
			"geo_get_aut_gen new_first");
		if (new_first == NIL)
			return error("geo_get_aut_gen(): no memory for new_first");
		for (i = 0; i < geo->nb_aut_gens; i++)
			new_first[i] = geo->auts_first_moved[i];
		if (geo->dim_aut_gens) {
			memcpy(new_auts, geo->aut_gens, geo->nb_aut_gens * sizeof(SPERM));
			my_free(geo->aut_gens);
			my_free(geo->auts_first_moved);
			}
		geo->aut_gens = new_auts;
		geo->auts_first_moved = new_first;
		geo->dim_aut_gens = new_dim;
		}
	p = geo->aut_gens + geo->nb_aut_gens;
	n = geo->b;
	sp_nil(p);
	sp_int(p, n, "sp_int geo_get_aut_gen");
	sp_mult(&geo->q0, &geo->Q[geo->b]->pv /* &L->pA->pv */, p);

#if 0
	if (f_v) {
		printf("pv=");
		sp_print(&L->pA->pv);
		printf("p0=");
		sp_print(&L->p0);
		sp_print(p);
		}
#endif
	geo->auts_first_moved[geo->nb_aut_gens] = geo->first_moved;

	geo->nb_aut_gens++;
	return OK;
}

void geo_free_aut_gen(GEO_CANON *geo)
{
	INT i;
	SPERM *p;

	if (geo->dim_aut_gens) {
		for (i = 0; i < geo->nb_aut_gens; i++) {
			p = geo->aut_gens + i;
			sp_free(p);
			}
		my_free(geo->aut_gens);
		my_free(geo->auts_first_moved);
		}
	geo->aut_gens = NIL;
	geo->auts_first_moved = 0;
	geo->nb_aut_gens = 0;
	geo->dim_aut_gens = 0;
}

INT geo_get_canonical_labelling(GEO_CANON *geo)
{
	NTDO_GRID *Q, *P;
	
	Q = geo->Q[geo->b];
	P = geo->P[geo->b];
	sp_mv(&Q->p, &geo->q0);
	sp_mv(&Q->pv, &geo->q0v);
	sp_mv(&P->p, &geo->p0);
	sp_mv(&P->pv, &geo->p0v);
	return OK;
}

INT geo_canon_apply(GEO_CANON *geo)
// computes the new theX 
// according to the canonical form.
{
	GEO_CANON_INFO *info;
	INT *M;
	INT i, j, k;

	info = geo->info;
	M = (INT *) my_malloc(info->v * info->b * sizeof(INT), "geo_canon_apply() M");
	if (M == NIL)
		return error("geo_canon_apply() no memory for M");
	for (i = 0; i < info->v; i++) {
		for (j = 0; j < info->b; j++) {
			M[i * info->b + j] = 0;
			}
		}
	for (k = 0; k < geo->nb_X; k++) {
		i = geo->theXi[k];
		j = geo->theXj[k];
		i = geo->p0.a[i];
		j = geo->q0.a[j];
		if (info->f_transposed) {
			M[j * info->b + i] = 1;
			}
		else {
			M[i * info->b + j] = 1;
			}
		}
	k = 0;
	for (i = 0; i < info->v; i++) {
		for (j = 0; j < info->b; j++) {
			if (M[i * info->b + j]) {
				info->theX[k] = i * info->b + j;
				k++;
				}
			}
		}
	my_free(M);
	return OK;
}

INT geo_get_M_ij(GEO_CANON *GC, INT i, INT j)
{
	INT val = 0;
	
#ifdef GEO_CANON_ALLOC_M
	val = GC->M[i * GC->b + j];
#else
	INT k, jj, ii;
	
	for (k = 0; k < GC->nb_X; k++) {
		ii = GC->theXi[k];
		if (ii != i)
			continue;
		jj = GC->theXj[k];
		if (jj != j)
			continue;
		val = 1;
		break;
		}
#endif
	return val;
}

INT geo_compare(GEO_CANON *geo, INT k)
{
	INT col, ret;
	
	for (col = 0; col <= k; col++) {
		ret = geo_compare_col(geo, k, col);
		if (ret) {
			// printf("geo_compare() level = %ld ret = %ld\n", k, ret);
			return ret;
			}
		}
	// printf("geo_compare() level = %ld ret = 0\n", k);
	return 0;
}

#if 0
INT geo_compare(GEO_CANON *geo, INT k)
{
	INT i, i0, j0, i1, j1, ret;
	
	j0 = geo->q0v.a[k];
	j1 = geo->Q[k + 1]->pv.a[k];
	for (i = 0; i < geo->v; i++) {
		i0 = geo->p0v.a[i];
		i1 = geo->P[k + 1]->pv.a[i];
		ret = geo->M[i1 * geo->b + j1] - geo->M[i0 * geo->b + j0];
		if (ret)
			return ret;
		}
	return 0;
}
#endif

INT geo_compare_col(GEO_CANON *geo, INT k, INT col)
{
	INT i, i0, j0, i1, j1, ret;
	
	j0 = geo->q0v.a[col];
	j1 = geo->Q[k + 1]->pv.a[col];
	for (i = 0; i < geo->v; i++) {
		i0 = geo->p0v.a[i];
		i1 = geo->P[k + 1]->pv.a[i];
		ret = geo_get_M_ij(geo, i1, j1) - geo_get_M_ij(geo, i0, j0);
		// ret = geo->M[i1 * geo->b + j1] - geo->M[i0 * geo->b + j0];
		if (ret)
			return ret;
		}
	return 0;
}

INT print_grid_decomp(NTDO_GRID *Q)
{
	INT j;
	
	for (j = 0; j < Q->G_max; j++)
		printf("[%2ld,%2ld] ", Q->first[j], Q->len[j]);
	printf("\n");
	return OK;
}


