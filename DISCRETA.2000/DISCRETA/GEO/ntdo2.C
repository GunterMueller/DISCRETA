/* ntdo2.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>


#include <DISCRETA/geo.h>


INT ntdo_calc_second_tdo(NTDO *tdo, INT f_v, INT f_vv)
{
	TDO_GRAD tdog;
	NTDO_FRAME *frame = NIL;
	NTDO_GRID *G_last, *G, *G_next, *Gtmp;
	SPERM p, q, pv, qv;
	INT f_points, m, n, max_size;
	
	// printf("in ntdo_calc_second_tdo()\n"); fflush(stdout);
	max_size = tdo->G->max_size;
	if (max_size != tdo->G_last->max_size) {
		return error("ntdo_calc_second_tdo(): max_size != tdo->G_last->max_size");
		}
	if (max_size != tdo->G_next->max_size) {
		return error("ntdo_calc_second_tdo(): max_size != tdo->G_next->max_size");
		}
	sp_nil(&p);
	sp_nil(&q);
	sp_nil(&pv);
	sp_nil(&qv);
	// printf("max_size =%ld\n", max_size); fflush(stdout);
	// printf("v =%ld\n", tdo->v); fflush(stdout);
	// printf("b =%ld\n", tdo->b); fflush(stdout);
	sp_int(&p, tdo->v, "sp_int ntdo_calc_second_tdo");
	sp_int(&pv, tdo->v, "sp_int ntdo_calc_second_tdo");
	sp_int(&q, tdo->b, "sp_int ntdo_calc_second_tdo");
	sp_int(&qv, tdo->b, "sp_int ntdo_calc_second_tdo");
	frame = (NTDO_FRAME *) my_malloc(sizeof(NTDO_FRAME), "ntdo_calc_second_tdo NTDO_FRAME");
	if (frame == NIL) {
		return error("ntdo_calc_second_tdo(): no memory for NTDO_FRAME");
		}
	/* save for being restored at the end */
	G_last = tdo->G_last;
	G = tdo->G;
	G_next = tdo->G_next;
	if (f_vv) {
		printf("G_last = tdo->G_last = %lx\n", G_last);
		printf("G = tdo->G = %lx\n", G);
		printf("G_next = tdo->G_next = %lx\n", G_next);
		}
	
	tdo->G_last = ntdo_grid_init(max_size);
	tdo->G = ntdo_grid_init(max_size);
	tdo->G_next = ntdo_grid_init(max_size);
	if (tdo->G_last == NIL ||
		tdo->G == NIL ||
		tdo->G_next == NIL) {
		return error("ntdo_calc_second_tdo() no memory for NTDO_GRID");
		}
	if (f_vv) {
		printf("tdo->G_last = %lx\n", tdo->G_last);
		printf("tdo->G = %lx\n", tdo->G);
		printf("tdo->G_next = %lx\n", tdo->G_next);
		}

	tdog_nil(&tdog);
	
	while (TRUE) {
		/* aktuelles TDO jetzt in G, G_next; 
		 * aktuelle Permutation 
		 * in p, pv, q, qv. */
		f_points = TRUE;
		m = tdo->v;
		n = tdo->b;
		
		/* save tdo->p: */
		sp_mv(&tdo->p, &p);
		sp_mv(&tdo->pv, &pv);
		sp_mv(&tdo->q, &q);
		sp_mv(&tdo->qv, &qv);
		if (ntdo_refine(tdo, &tdog, G, G_next, 
			f_points, frame, &p, &pv, &q, &qv, f_vv) != OK) {
			return error("ntdo_calc_second_tdo(): ntdo_refine");
			}
		/* Das (zeilenweise) verfeinerte 
		 * TDO wird durch 
		 * frame beschrieben (nur die Zeilen). */
		/* Die (Zeilen-) Permutationen 
		 * auf das verfeinerte TDO 
		 * befinden sich jetzt in p, pv */
		if (G->f_points == f_points) {
			ntdo_grid_copy_frame(G_next, tdo->G_last, FALSE);
				/* Alte Block-TDO-Einteilung von G_next 
				 * nach tdo->G_last kopieren. */
				/* setzt f_points, m, n */
			}
		else {
			ntdo_grid_copy_frame(G, tdo->G_last, FALSE);
				/* Alte Block-TDO-Einteilung von G 
				 * nach tdo->G_last kopieren. */
			}
		frame2ntdo_grid(frame, tdo->G);
			/* Neue Punkt-TDO-Einteilung 
			 * nach tdo->G. */
		
		tdo->G->f_points = f_points;
		tdo->G->m = m;
		tdo->G->n = n;
		ntdo_calc_grid(tdo, tdo->G_last, tdo->G);
		sp_mv(&p, &tdo->p);
		sp_mv(&pv, &tdo->pv);
		sp_mv(&q, &tdo->q);
		sp_mv(&qv, &tdo->qv);
		if (!ntdo_calc(tdo, FALSE)) {
			return error("ntdo_calc_second_tdo(): error in ntdo_calc");
			}
		/* Die Permutationen auf das (zeilenweise) 
		 * verfeinerte TDO befinden sich jetzt in 
		 * tdo->p, tdo->pv, tdo->q, tdo->qv. */
		if (tdos_init_ntdo(&tdo->tdos2, tdo, tdo->G, tdo->G_next, TRUE) != OK) {
			return error("ntdo_calc_second_tdo() tdos_int");
			}
		if (f_v) {
			printf("ntdo_calc_second_tdo after point-refinement:\n");
			tdos_print(&tdo->tdos2, tdo);
			}
		
		
		/* 
		 * now comes the block refinement: 
		 */
		
		Gtmp = G;
		G = tdo->G;
		tdo->G = Gtmp;
		Gtmp = G_next;
		G_next = tdo->G_next;
		tdo->G_next = Gtmp;
		
		/* save tdo->p: */
		sp_mv(&tdo->p, &p);
		sp_mv(&tdo->pv, &pv);
		sp_mv(&tdo->q, &q);
		sp_mv(&tdo->qv, &qv);
		f_points = FALSE;
		m = tdo->b;
		n = tdo->v;
		if (ntdo_refine(tdo, &tdog, G, G_next, 
			f_points, frame, &p, &pv, &q, &qv, f_vv) != OK) {
			return error("ntdo_calc_second_tdo() ntdo_refine");
			}
		/* Das (spaltenweise) verfeinerte 
		 * TDO wird durch 
		 * frame beschrieben. */
		/* Die (Spalten-) Permutationen 
		 * auf das verfeinerte TDO 
		 * befinden sich jetzt in q, qv */
		if (G->f_points == f_points) {
			ntdo_grid_copy_frame(G_next, tdo->G_last, FALSE);
				/* Alte Punkt-TDO-Einteilung von G_next 
				 * nach tdo->G_last kopieren. */
				/* setzt f_points, m, n */
			}
		else {
			ntdo_grid_copy_frame(G, tdo->G_last, FALSE);
				/* Alte Punkt-TDO-Einteilung von G 
				 * nach tdo->G_last kopieren. */
			}
		frame2ntdo_grid(frame, tdo->G);
			/* Neue Block-TDO-Einteilung 
			 * nach tdo->G. */
		tdo->G->max_size = tdo->G_last->max_size;
		tdo->G->f_points = f_points;
		tdo->G->m = m;
		tdo->G->n = n;
		ntdo_calc_grid(tdo, tdo->G_last, tdo->G);
		sp_mv(&p, &tdo->p);
		sp_mv(&pv, &tdo->pv);
		sp_mv(&q, &tdo->q);
		sp_mv(&qv, &tdo->qv);
		if (!ntdo_calc(tdo, FALSE)) {
			return error("ntdo_calc_second_tdo() ntdo_calc");
			}
		/* Die Permutationen auf das (spaltenweise) 
		 * verfeinerte TDO befinden sich jetzt in 
		 * tdo->p, tdo->pv, tdo->q, tdo->qv. */
		if (tdos_init_ntdo(&tdo->tdos2, tdo, tdo->G, tdo->G_next, TRUE) != OK) {
			return error("ntdo_calc_second_tdo() tdos_init_ntdo");
			}
		if (f_v) {
			printf("ntdo_calc_second_tdo after block-refinement:\n");
			tdos_print(&tdo->tdos2, tdo);
			}
		break;
		}
	
	tdog_exit(&tdog);
	if (frame) {
		my_free(frame);
		frame = NIL;
		}
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
	/* restore: */
	tdo->G_last = G_last;
	tdo->G = G;
	tdo->G_next = G_next;
	sp_free(&p);
	sp_free(&pv);
	sp_free(&q);
	sp_free(&qv);
	return OK;
}

INT ntdo_refine(NTDO *tdo, TDO_GRAD *tdog, 
	NTDO_GRID *G, NTDO_GRID *G_next, 
	INT f_points, NTDO_FRAME *frame, 
	SPERM *p, SPERM *pv, SPERM *q, SPERM *qv, INT f_v)
{
	TDO_SCHEME tdos;
	NTDO_GRID *G0, *G1;
	SPERM *perm, *perm_inv;
	SPERM *tdo_perm, *tdo_perm_inv;
	INT first, len, first1, len1;
	INT i, i1, j, j1, k, l, t;
	INT f_vv = FALSE;
	
	tdos_nil(&tdos);
	if (f_vv) {
		printf("ntdo_refine():\n");
		printf("tdo->G_last = %lx\n", tdo->G_last);
		printf("tdo->G = %lx\n", tdo->G);
		printf("tdo->G_next = %lx\n", tdo->G_next);
		printf("G = %lx\n", G);
		printf("G_next = %lx\n", G_next);
		}
	if (G->f_points == f_points) {
		G0 = G_next;
		G1 = G;
		}
	else {
		G0 = G;
		G1 = G_next;
		}
	/* G1 ist die unverfeinerte 
	 * TDO - Einteilung, 
	 * G0 die Gegenrichtungs 
	 * TDO - Einteilung (auch die alte). */
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
	frame->G_max = 0;
	if (f_v) {
		printf("ntdo_refine(): f_points = %ld = G1->f_points  = %ld, "
			"G0->f_points = %ld\n", 
			f_points, G1->f_points, G0->f_points);
		printf("G1->max_size = %ld G0->max_size = %ld\n", 
			G1->max_size, G0->max_size);
		fflush(stdout);
		printf("G0=\n");
		ntdo_grid_print(G0);
		printf("G1=\n");
		ntdo_grid_print(G1);
		}
	for (i = 0; i < G1->G_max; i++) {
		first = G1->first[i];
		len = G1->len[i];
		if (f_v) {
			printf("in ntdo_refine(): i = %ld first = %ld len = %ld\n", 
				i, first, len);
			fflush(stdout);
			}
		tdog->f_points = f_points;
		tdog->n = len;
		tdog->N = len;
		tdog->type = (INT *) my_malloc(tdog->N * sizeof(INT), "tdo_refine tdog->type");
		if (tdog->type == NIL) {
			return error("tdo_refine() no memory for tdog->type");
			}
		tdog->tdos = (TDO_SCHEME *) my_malloc(tdog->N * sizeof(TDO_SCHEME), 
			"tdo_refine TDO_SCHEME");
		if (tdog->tdos == NIL) {
			return error("tdo_refine() no memory for tdog->tdos");
			}
		for(i1 = 0; i1 < tdog->N; i1++) {
			tdos_nil(&tdog->tdos[i1]);
			}
		tdog->mult = (INT *) my_malloc(tdog->N * sizeof(INT), 
			"tdo_refine tdog->mult");
		if (tdog->mult == NIL) {
			return error("ntdo_refine() no memory for tdog->mult");
			}
		
		for (j = 0; j < len; j++) {
			if (f_v) {
				printf("in ntdo_refine(): j = %ld\n", j);
				fflush(stdout);
				}
			sp_mv(p, &tdo->p);
			sp_mv(pv, &tdo->pv);
			sp_mv(q, &tdo->q);
			sp_mv(qv, &tdo->qv);
			if (f_vv) {
				printf("in ntdo_refine() calling ntdo_grid_copy_frame()\n");
				fflush(stdout);
				printf("G1=\n");
				ntdo_grid_print(G1);
				fflush(stdout);
				}
			ntdo_grid_copy_frame(G0, tdo->G_last, f_vv);
			/* Gegenrichtung aus G0 nach 
			 * tdo->G_last uebertragen. */
			if (j != 0) {
				if (f_vv) {
					printf("in ntdo_refine() calling sp_mult_apply_tau_r()\n");
					fflush(stdout);
					}
				sp_mult_apply_tau_r(tdo_perm, first, first + j);
				sp_mult_apply_tau_l(tdo_perm_inv, first, first + j);
				}
			if (f_vv) {
				printf("in ntdo_refine() calling ntdo_grid_init_derived_i_first()\n");
				fflush(stdout);
				printf("G1=\n");
				ntdo_grid_print(G1);
				fflush(stdout);
				}
			if (!ntdo_grid_init_derived_i_first(tdo, tdo->G, G1, i, f_vv)) {
				return error("tdo_refine() ntdo_grid_init_derived_i_first");
				}
			/* In tdo->G ist jetzt die erste Zeile des 
			 * i-ten Bereichs ausgezeichnet. */
			if (f_vv) {
				printf("in ntdo_refine() calling ntdo_calc()\n");
				fflush(stdout);
				}
			if (!ntdo_calc(tdo, FALSE)) {
				return error("ntdo_refine() error in ntdo_calc()");
				}
			if (f_vv) {
				printf("in ntdo_refine() calling tdos_init_ntdo()\n");
				fflush(stdout);
				}
			if (tdos_init_ntdo(&tdos, tdo, tdo->G, tdo->G_next, TRUE) != OK) {
				return error("ntdo_refine() error in tdos_init_ntdo()");
				}
			if (f_v) {
				printf("ntdo_refine(): derivation at %ld (%ld/%ld)\n", first + j, i, j);
				fflush(stdout);
				tdos_print(&tdos, tdo);
				fflush(stdout);
				}
			if (f_vv) {
				printf("in ntdo_refine() calling tdog_add_tdos()\n");
				fflush(stdout);
				}
			if (!tdog_add_tdos(tdog, &tdos, j)) {
				return error("ntdo_refine() error in tdog_add_tdos()");
				}
			} /* next j */
		
		if (f_v) {
			printf("we found the following types:\n");
			for (k = 0; k < len; k++) {
				printf("%ld ", tdog->type[k]);
				}
			printf("\n");
			fflush(stdout);
			}
#ifdef DEBUG_NTDO_REFINE
		/* tdog_print1(tdog); */
		printf("mult: ");
		for (k = 0; k < tdog->nb_tdos; k++) {
			printf("%ld ", tdog->mult[k]);
			}
		printf("\n");
		printf("type: ");
		for (k = 0; k < len; k++) {
			printf("%ld ", tdog->type[k]);
			}
		printf("\n");
#endif

		/* Verfeinerung des Bereichs (first, len)
		 * nach Typen der Ableitungen
		 * (jetzt in tdog), sortieren 
		 * und nach frame schreiben.
		 * Die erfolgten Permutationen werden nach 
		 * p/pv bzw. q/qv geschrieben: */
		j = 0;
		for (k = 0; k < tdog->nb_tdos; k++) {
			first1 = j;
			for (j1 = j; j1 < len; j1++) {
				if (tdog->type[j1] == k) {
					if (j1 != j) {
						/* (j j+1 ... j1) anwenden: */
						t = tdog->type[j1];
						for (l = j1; l > j; l--)
							tdog->type[l] = tdog->type[l - 1];
						tdog->type[j] = t;
						if (f_vv) {
							printf("cycle %ld %ld\n", first + j, j1 - j + 1);
							}

						sp_mult_apply_forwc_r(perm, first + j, j1 - j + 1);
						sp_mult_apply_backwc_l(perm_inv, first + j, j1 - j + 1);
						}
					j++;
					}
				} /* next j1 */
			
			len1 = j - first1;
			/* Eintrag (first + first1, len1): */
			for (l = 0; l < len1; l++)
				frame->grid_entry[first + first1 + l] = frame->G_max;
			frame->first[frame->G_max] = first + first1;
			frame->len[frame->G_max] = len1;
			
			/* printf("at %ld len = %ld\n", 
				frame->first[frame->G_max], 
				frame->len[frame->G_max]); */
			frame->G_max++;
			/* frame->first[frame->G_max] = first + first1 + len1; */
			}
		tdog_exit(tdog);
		
		} /* next i */
	return OK;
}

