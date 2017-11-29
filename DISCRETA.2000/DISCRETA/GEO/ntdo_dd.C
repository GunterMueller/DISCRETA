/* ntdo_dd.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */



#include <DISCRETA/discreta.h>


#include <DISCRETA/geo.h>

INT ntdo_calc_dd(NTDO *tdo, INT f_points, INT f_blocks, 
	SHORT **ddp, INT *Np, SHORT **ddp_mult, 
	SHORT **ddb, INT *Nb, SHORT **ddb_mult, 
	INT f_v, INT f_vv)
{
	INT f_print = FALSE;
	SHORT *mult;
	INT i, s, a;
	
#if 0
	if (f_points || f_blocks) {
		f_print = TRUE;
		}
#endif
	if (f_print)
		printf("(");
	if (f_points) {
		if (ntdo_dd(tdo, TRUE /* f_points */, ddp, Np, ddp_mult, f_v, f_vv) != OK) {
			return error("ntdo_calc_dd() error in ntdo_dd()");
			}
		if (f_print) {
			mult = *ddp_mult;
			s = 0;
			for (i = 1; i <= mult[0]; i++) {
				a = mult[i];
				s += a;
				printf("%ld ", a);
				}
			printf(": %ld", s);
			}
		}
	if (f_print)
		printf("|");
	if (f_blocks) {
		if (ntdo_dd(tdo, FALSE /* f_points */, ddb, Nb, ddb_mult, f_v, f_vv) != OK) {
			return error("ntdo_calc_dd() error in ntdo_dd()");
			}
		if (f_print) {
			mult = *ddb_mult;
			s = 0;
			for (i = 1; i <= mult[0]; i++) {
				a = mult[i];
				s += a;
				printf("%ld ", a);
				}
			printf(": %ld", s);
			}
		}
	if (f_print)
		printf(")\n");
	return OK;
}

INT ntdo_dd(NTDO *tdo, INT f_points, 
	SHORT **dd, INT *N, SHORT **dd_mult, INT f_v, INT f_vv)
/* Allociert dd auf einen N Vector von SHORTs. 
 * dd enthaelt fuer alle Paare 
 * (in der aktuellen Inzidenzlage) 
 * den Index des Ableitungs TDOs.
 * Allociert dd_mult auf einen 
 * Multiplizitaetenvektor; 
 * dd_mult[0] = Anzahl der auftretenden 
 *   Ableitungs TDOs.
 * dd_mult ist dann ein Vector 
 * der Laenge (dd_mult[0] + 1) sizeof(SHORT).
 * setzt N := n * (n - 1) >> 1, 
 * wobei n = tdo->v if f_points, 
 * = tdo->inc->B sonst. */
{
	TDO_GRAD tdog_, *tdog = &tdog_;
	TDO_SCHEME tdos_, *tdos = &tdos_;
	NTDO_GRID *G_last = NIL;
	NTDO_GRID *G = NIL;
	NTDO_GRID *G_next = NIL;
	NTDO_GRID *G0, *G1;
	SPERM p_, *p = &p_;
	SPERM pv_, *pv = &pv_;
	SPERM q_, *q = &q_;
	SPERM qv_, *qv = &qv_;
	
	SPERM *perm, *perm_inv;
	SPERM *tdo_perm, *tdo_perm_inv;
	INT I, I_len, i, i0;
	INT i1, J, J_len, j, j0, j1;
	INT k, max_size;
	INT ret = FALSE;
	SHORT *dd1 = NIL;
	SHORT *dd_mult1 = NIL;

	if (f_v) {
		printf("ntdo_dd(\n");
		fflush(stdout);
		}
	*dd = NIL;
	*dd_mult = NIL;
	tdog_nil(tdog);
	tdos_nil(tdos);
	sp_nil(p);
	sp_nil(q);
	sp_nil(pv);
	sp_nil(qv);
	sp_int(p, tdo->v, "sp_int ntdo_dd");
	sp_int(pv, tdo->v, "sp_int ntdo_dd");
	sp_int(q, tdo->b, "sp_int ntdo_dd");
	sp_int(qv, tdo->b, "sp_int ntdo_dd");
	/* save tdo->p: */
	sp_mv(&tdo->p, p);
	sp_mv(&tdo->pv, pv);
	sp_mv(&tdo->q, q);
	sp_mv(&tdo->qv, qv);

	max_size = tdo->G->max_size;
	G_last = tdo->G_last;
	G = tdo->G;
	G_next = tdo->G_next;
	tdo->G_last = ntdo_grid_init(max_size);
	tdo->G = ntdo_grid_init(max_size);
	tdo->G_next = ntdo_grid_init(max_size);
	if (tdo->G_last == NIL ||
		tdo->G == NIL ||
		tdo->G_next == NIL) {
		return error("ntdo_dd() no memory for GRID");
		}
	if (G->f_points == f_points) {
		G0 = G_next;
		G1 = G;
		}
	else {
		G0 = G;
		G1 = G_next;
		}
	/* G1 ist die TDO - Einteilung, 
	 * G0 die Gegenrichtungs 
	 * TDO - Einteilung. */
	
	if (f_points) {
		perm = p;
		perm_inv = pv;
		tdo_perm = &tdo->p;
		tdo_perm_inv = &tdo->pv;
		}
	else {
		perm = q;
		perm_inv = qv;
		tdo_perm = &tdo->q;
		tdo_perm_inv = &tdo->qv;
		}
	
	/* tdog initialisieren: */
	tdog->f_points = f_points;
	if (f_points)
		tdog->n = tdo->v;
	else
		tdog->n = tdo->b;
	tdog->N = tdog->n * (tdog->n - 1) >> 1;
	tdog->type = (INT *) my_malloc(tdog->N * sizeof(INT), "ntdo_dd tdog->type");
	if (tdog->type == NIL) {
		return error("ntdo_dd() no memory for tdog->type");
		}
	tdog->tdos = (TDO_SCHEME *) my_malloc(tdog->N * sizeof(TDO_SCHEME), "ntdo_dd tdog->tdos");
	if (tdog->tdos == NIL) {
		return error("ntdo_dd() no memory for tdog->tdos");
		}
	for (i = 0; i < tdog->N; i++) {
		tdos_nil(&tdog->tdos[i]);
		}
	tdog->mult = (INT *) my_malloc(tdog->N * sizeof(INT), "ntdo_dd tdog->mult");
	if (tdog->mult == NIL) {
		return error("ntdo_dd() no memory for tdog->mult");
		}
	
	/* ueber alle Paare (i1, j1) 
	 * mit i1 < j1 d.h. geeignete 
	 * (I, i), (J, j) Kombinationen: */
	for (I = 0; I < G1->G_max; I++) {
		i0 = G1->first[I];
		I_len = G1->len[I];
		if (f_vv) {
			printf("I=%ld i0 = %ld I_len=%ld\n", I, i0, I_len);
			fflush(stdout);
			}
		for (i = 0; i < I_len; i++) {
			i1 = i0 + i;
			
			if (f_vv) {
				printf("i1 = %ld\n", i1);
				fflush(stdout);
				}
			for (J = I; J < G1->G_max; J++) {
				j0 = G1->first[J];
				J_len = G1->len[J];
				if (I == J)
					j = i + 1;
				else
					j = 0;
				for (; j < J_len; j++) {
					j1 = j0 + j;
					if (f_vv) {
						printf("ntdo_grid_copy_frame(");
						fflush(stdout);
						}
					ntdo_grid_copy_frame(G0, tdo->G_last, FALSE);
					/* Gegenrichtung aus G0 
					 * nach tdo->G_last uebertragen. */
					if (f_vv) {
						printf(")\n");
						fflush(stdout);
						}
					if (f_vv) {
						printf("ntdo_grid_init_derived_ij_first(");
						fflush(stdout);
						}
					if (ntdo_grid_init_derived_ij_first(
						tdo, tdo->G, G1 /* G_old */, I, J) != OK) {
						return error("ntdo_dd() error in ntdo_grid_init_derived_ij_first");
						}
					if (f_vv) {
						printf(")\n");
						fflush(stdout);
						}
					if (f_vv) {
						printf("if II == J(");
						fflush(stdout);
						}
					if (I == J) {
						if (i != 0) {
							sp_mult_apply_tau_r(tdo_perm, i0, i1);
							sp_mult_apply_tau_l(tdo_perm_inv, i0, i1);
							/* dies hat j1 nicht 
							 * bewegt, da i1 < j1 */
							}
						if (j != 1) {
							sp_mult_apply_tau_r(tdo_perm, j0 + 1, j1);
							sp_mult_apply_tau_l(tdo_perm_inv, j0 + 1, j1);
							}
						}
					else {
						if (i != 0) {
							sp_mult_apply_tau_r(tdo_perm, i0, i1);
							sp_mult_apply_tau_l(tdo_perm_inv, i0, i1);
							}
						if (j != 0) {
							sp_mult_apply_tau_r(tdo_perm, j0, j1);
							sp_mult_apply_tau_l(tdo_perm_inv, j0, j1);
							}
						}
					if (f_vv) {
						printf(")\n");
						fflush(stdout);
						}
					if (f_vv) {
						printf("ntdo_calc(");
						fflush(stdout);
						}
					if (!ntdo_calc(tdo, FALSE)) {
						return error("ntdo_dd() error in ntdo_calc()");
						}
					if (f_vv) {
						printf(")\n");
						fflush(stdout);
						}
					if (tdos_init_ntdo(tdos, tdo, tdo->G, tdo->G_next, TRUE) != OK) {
						return error("ntdo_dd() error in tdos_int()");
						}
					k = ij2k(i1, j1, tdog->n);
					if (f_vv) {
						printf("ntdo_dd(): derivation at %ld,%ld (%ld/%ld,%ld/%ld)(\n", i1, j1, I, i, J, j);
						fflush(stdout);
						// tdos_print(tdos, tdo);
						// fflush(stdout);
						}
					if (!tdog_add_tdos(tdog, tdos, k)) {
						return error("ntdo_dd() error in tdog_add_tdos()");
						}
					if (f_vv) {
						printf(")\n");
						fflush(stdout);
						}
					/* urspruengliche Zeilen / 
					 * Spaltenlage wiederherstellen: */
					sp_mv(p, &tdo->p);
					sp_mv(pv, &tdo->pv);
					sp_mv(q, &tdo->q);
					sp_mv(qv, &tdo->qv);
					
					if (f_vv) {
						printf("mv)\n");
						fflush(stdout);
						}
					} /* next j */
				} /* next J */
			} /* next i */
		} /* next I */
	dd1 = (SHORT *) my_malloc(tdog->N * sizeof(SHORT), "ntdo_dd dd1");
	if (dd1 == NIL) {
		return error("ntdo_dd() no memory for dd1");
		}
	dd_mult1 = (SHORT *) my_malloc((tdog->nb_tdos + 1) * sizeof(SHORT), "ntdo_dd dd_mult1");
	if (dd_mult1 == NIL) {
		return error("ntdo_dd() no memory for dd_mult1");
		}
	for (i = 0; i < tdog->N; i++)
		dd1[i] = (SHORT) tdog->type[i];
	dd_mult1[0] = tdog->nb_tdos;
	for (i = 0; i < tdog->nb_tdos; i++)
		dd_mult1[1 + i] = (SHORT) tdog->mult[i];
	if (f_v) {
		printf("dd:\n");
		for (i = 0; i < tdog->N; i++)
			printf("%ld ", (INT)dd1[i]);
		printf("\n");
		printf("mult:\n");
		for (i = 0; i < (INT) dd_mult1[0]; i++)
			printf("%ld ", (INT)dd_mult1[1 + i]);
		printf("\n");
		}
	*dd = dd1;
	dd1 = NIL;
	*N = tdog->N;
	*dd_mult = dd_mult1;
		
	ret = TRUE;
	tdog_exit(tdog);
	if (tdo->G_last) {
		ntdo_grid_free(tdo->G_last);
		tdo->G_last = NIL;
		}
	if (tdo->G) {
		ntdo_grid_free(tdo->G);
		tdo->G = NIL;
		}
	if (tdo->G_next) {
		ntdo_grid_free(tdo->G_next);
		tdo->G_next = NIL;
		}
	tdo->G_last = G_last;
	tdo->G = G;
	tdo->G_next = G_next;
	sp_free(p);
	sp_free(pv);
	sp_free(q);
	sp_free(qv);
	if (f_v) {
		printf(")ntdo_dd\n");
		fflush(stdout);
		}
	return OK;
}

