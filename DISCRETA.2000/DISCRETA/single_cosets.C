/* single_cosets.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LABRA_TRUE

#include <DISCRETA/lb.h>
#include <DISCRETA/ma.h>

static INT single_cosets(SINGLE_COSET_WORK *scw, INT i0, VECTOR_OP Reps);
static INT sc_maxtest(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep);
static INT single_coset_first_rep(SINGLE_COSET_WORK *scw);
static INT single_coset_next_rep(SINGLE_COSET_WORK *scw, INT i0);
static INT single_coset_next_back(SINGLE_COSET_WORK *scw, INT i0, INT r);
static void print_trans_rep(SINGLE_COSET_WORK *scw);
static INT calc_rep(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep);


INT single_coset_labra_file()
{
	SINGLE_COSET_WORK *scw;
	BYTE *gen_G_fname = "M24_labra.dsc";
	BYTE *gen_U_fname = "MOG_labra.dsc";

	SYM_OB go, go_;
	LABRA_OB labG, labU;
	VECTOR_OB B, R;
	PERMUTATION_OB rep;
	INT no;
	
#if 0
	printf("calling construct_W24()\n"); fflush(stdout);
	construct_W24(&B, TRUE /* f_v */ );
	printf("finished with construct_W24()\n"); fflush(stdout);
#endif
	labG.load(gen_G_fname);
	labU.load(gen_U_fname);

	
	scw = single_coset_open(&labG, &labU, TRUE);
	
	// single_cosets(scw, 0, &R);
	
	no = 0;
	single_coset_first(scw, &rep);
	do {
		no++;
		printf("%ld ", no);
		if (no % 10 == 0)
			printf("\n");
		fflush(stdout);
		;
	} while (single_coset_next(scw, &rep));
	printf("\n");
	
	single_coset_free(scw);
	
	return OK;
}

SINGLE_COSET_WORK *single_coset_open(LABRA_OP G, LABRA_OP U, INT f_v)
{
	SINGLE_COSET_WORK *scw;
	MATRIX_OB TG, TU;
	INT i, j, n;
	
	if (f_v) {
		G->Print_Sims();
		U->Print_Sims();
		}
	
	G->calc_transversal_matrix(&TG, FALSE);
	U->calc_transversal_matrix(&TU, FALSE);

	scw = (SINGLE_COSET_WORK *) my_malloc(
		sizeof(SINGLE_COSET_WORK), "SINGLE_COSET_WORK");
	scw->G = G;
	scw->U = U;
	// scw->MG = &TG;
	// scw->MU = &TU;
	scw->deg = n = G->s_degree_i();
	if (n >= SINGLE_COSET_MAX_DEG) {
		error("single_coset_open() n >= SINGLE_COSET_MAX_DEG");
		return NIL;
		}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			scw->MG[i][j] = TG.s_iji(i, j);
			scw->MU[i][j] = TU.s_iji(i, j);
			}
		}
	return scw;
}

INT single_coset_free(SINGLE_COSET_WORK *scw)
{
	my_free(scw);
	return OK;
}

INT single_coset_first(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep)
{
	INT r;

	single_coset_first_rep(scw);
	while (TRUE) {
		// printf("cand %ld: ", cand);
		// print_trans_rep();
		calc_rep(scw, rep);
		// rep->println();
		
		r = sc_maxtest(scw, rep);
		if (r == -1) {
			return OK;
			}
		else {
			// printf("not canonic at position %ld\n", r);
			if (!single_coset_next_back(scw, 0, r))
				break;
			}
		}
	error("single_coset_first() no coset representative found !");
	return OK;
}

INT single_coset_next(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep)
{
	INT r;

	if (!single_coset_next_rep(scw, 0))
		return FALSE;
	while (TRUE) {
		calc_rep(scw, rep);
		
		r = sc_maxtest(scw, rep);
		if (r == -1) {
			return TRUE;
			}
		else {
			if (!single_coset_next_back(scw, 0, r))
				break;
			}
		}
	return FALSE;
}

static INT single_cosets(SINGLE_COSET_WORK *scw, INT i0, VECTOR_OP Reps)
{
	INT ii, a, r, dots = 0;
	INT cand = 0;
	INT nb_R = 0;
	PERMUTATION_OB rep;
	VECTOR_OB rep1;

	Reps->m_il(0);
	single_coset_first_rep(scw);
	cand++;
	while (TRUE) {
		// printf("cand %ld: ", cand);
		// print_trans_rep();
		calc_rep(scw, &rep);
		// rep.println();
		
		r = sc_maxtest(scw, &rep);
		if (r == -1) {
			nb_R++;
			// printf("maximal, is a canonical coset representative no %ld\n", nb_R);
			printf("%ld ", nb_R);
			fflush(stdout);
			rep1.m_il(scw->deg);
			for (ii = 0; ii < scw->deg; ii++) {
				a = scw->trans_rep[ii];
				rep1.m_ii(ii, a);
				}
			Reps->inc();
			rep1.swap((VECTOR_OP) Reps->s_i(nb_R - 1));
			if (!single_coset_next_rep(scw, i0))
				break;
			}
		else {
			// printf("not canonic at position %ld\n", r);
			printf(".");
			dots++;
			if (dots % 50 == 0) {
				printf("\n");
				dots = 0;
				}
			fflush(stdout);
			if (!single_coset_next_back(scw, i0, r))
				break;
			}
		cand++;
		}
	printf("\n");
	printf("found %ld new single cosets\n", nb_R);
	return OK;
}

static INT sc_maxtest(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep)
{
	PERMUTATION_OB R, p1, p2, p3;
	INT i, j, ai, b, aj;
	
	rep->copy(&R);
	for (i = 0; i < scw->deg; i++) {
		ai = R.s_ii(i) - 1;
		for (j = 0; j < scw->deg; j++) {
			b = scw->MU[i][j];
			if (b) {
				aj = R.s_ii(j) - 1;
				if (aj < ai) {
					return j;
					}
				}
			} /* next j */
		scw->G->rep_ij(i, ai, &p1);
		p1.invers(&p2);
		R.mult(&p2, &p3);
		p3.swap(&R);
		} /* next i */
	return -1;
		// OK
}

static INT single_coset_first_rep(SINGLE_COSET_WORK *scw)
{
	INT i, j, a;
	
	for (i = 0; i < scw->deg; i++) {
		scw->trans_rep_last[i] = -1;
		}
	for (i = 0; i < scw->deg; i++) {
		scw->trans_rep[i] = -1;
		for (j = 0; j < scw->deg; j++) {
			a = scw->MG[i][j];
			if (a) {
				scw->trans_rep[i] = j;
				break;
				}
			}
		if (scw->trans_rep[i] == -1)
			return error("single_coset_first_rep() empty row");
		}
	return OK;
}

static INT single_coset_next_rep(SINGLE_COSET_WORK *scw, INT i0)
{
	INT i, j, j0, ii, a;
	
	for (i = 0; i < scw->deg; i++) {
		scw->trans_rep_last[i] = scw->trans_rep[i];
		}
	for (i = scw->deg - 1; i >= i0; i--) {
		j0 = scw->trans_rep[i];
		scw->trans_rep[i] = -1;
		for (j = j0 + 1; j < scw->deg; j++) {
			a = scw->MG[i][j];
			if (a) {
				scw->trans_rep[i] = j;
				break;
				}
			}
		if (scw->trans_rep[i] != -1) {
			for (ii = i + 1; ii < scw->deg; ii++) {
				scw->trans_rep[ii] = -1;
				for (j = 0; j < scw->deg; j++) {
					a = scw->MG[ii][j];
					if (a) {
						scw->trans_rep[ii] = j;
						break;
						}
					}
				if (scw->trans_rep[ii] == -1) {
					error("single_coset_next_rep() empty row");
					return FALSE;
					}
				} /* next ii */
			return TRUE;
			}
		}
	return FALSE;
}

static INT single_coset_next_back(SINGLE_COSET_WORK *scw, INT i0, INT r)
{
	INT i, j, j0, ii, a;
	
	for (i = 0; i < scw->deg; i++) {
		scw->trans_rep_last[i] = scw->trans_rep[i];
		}
	for (i = r; i >= i0; i--) {
		j0 = scw->trans_rep[i];
		scw->trans_rep[i] = -1;
		for (j = j0 + 1; j < scw->deg; j++) {
			a = scw->MG[i][j];
			if (a) {
				scw->trans_rep[i] = j;
				break;
				}
			}
		if (scw->trans_rep[i] != -1) {
			for (ii = i + 1; ii < scw->deg; ii++) {
				scw->trans_rep[ii] = -1;
				for (j = 0; j < scw->deg; j++) {
					a = scw->MG[ii][j];
					if (a) {
						scw->trans_rep[ii] = j;
						break;
						}
					}
				if (scw->trans_rep[ii] == -1) {
					error("single_coset_next_back() empty row");
					return FALSE;
					}
				} /* next ii */
			return TRUE;
			}
		}
	return FALSE;
}

static void print_trans_rep(SINGLE_COSET_WORK *scw)
{
	INT i;
	
	for (i = 0; i < scw->deg; i++) {
		printf("%ld ", scw->trans_rep[i]);
		}
	printf("\n");
}

static INT calc_rep(SINGLE_COSET_WORK *scw, PERMUTATION_OP rep)
{
	PERMUTATION_OB p1, p2;
	INT i, j;

	for (i = 0; i < scw->deg; i++) {
		j = scw->trans_rep[i];
		scw->G->rep_ij(i, j, &p1);
		if (i == 0)
			p1.copy(rep);
		else {
			// pre-multiply: 
			p1.mult(rep, &p2);
			p2.swap(rep);
			}
		}
	return OK;
}


#endif /* LABRA_TRUE */

