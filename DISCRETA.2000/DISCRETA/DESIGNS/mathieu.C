/* mathieu.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>

// for construct_Higman_Sims_176():

#ifdef SYM_GEO
#include <DISCRETA/geo.h>
#endif

#undef FREE_STAB
#define MAX_STEP 64

#if TEXDOCU
INT construct_W24(VECTOR_OP W24_blocks, 
	VECTOR_OP M24_gen, VECTOR_OP stab_gen, INT *block_idx, INT f_v)
#endif
{
	VECTOR_OB T, A;
	SYM_OB id;
	DCY_OB dc[MAX_STEP];
	DCY_OP Dc;
	INT erg = OK, i, l;
	INT up_to_k = 8, up_to_step;
	INT type;
	void *data;
	INT deg;
	VECTOR_OP D;
	SYM_OB go, stab_go;
	PERMUTATION_OP d;
	VECTOR_OB R;
	LABRA_OB labra_A, labra_A1;
	INT idx, f_found;
	
	up_to_step = 2 * up_to_k - 1;
	if (up_to_step >= MAX_STEP)
		return error("construct_W24() up_to_step too large");
	deg = 24;
	
	type = DO_TYPE_SYM;
	data = NIL;
	erg += M24_generators(&A);
	// A.save("M24_generators.dsc");
	reduce_generators_labra(&A, &go, FALSE /* f_verbose */, &labra_A);
	printf("group order: ");
	go.println();
	A.copy(M24_gen);
#if 0
	labra_A.reduced_generating_set(M24_gen, 
		FALSE /* f_bottom_up */, TRUE /* f_v */);
#endif

	do_copy(A.s_i(0), &id, type, data);
	do_one(&id, type, data);
	for (i = 0; i <= up_to_step; i++) {
		dc[i].initialize_Young(&T, deg, i /* step */, 
			0 /* type */, NIL /* data */);
		T.swap(dc[i].s_T());
		if (i <= up_to_step)
			dc[i].s_reduce_generators_mode()->m_i(2);
				/* 0: no reduction 
				 * 1: dimino 
				 * 2: labra */
		}
		
	dc[0].s_D()->m_il(1);
	do_copy(&id, dc[0].s_D_i(0), type, data);
	
	dc[0].s_Ad()->m_il(1);
	((SYM_OP) &A)->copy(dc[0].s_Ad()->s_i(0));
	dc[0].s_D()->println();
	dc[0].s_Ad()->println();

	for (i = 1; i <= up_to_step; i++)
		erg += dc_do_step(&dc[0], i, &id, deg /* deg */, 
			TRUE /* f_verbose */, type, data);

	Dc = &dc[up_to_step];
	D = Dc->s_D();
	l = D->s_li();
	printf("found %ld double cosets\n", l);
	for (i = 0; i < l; i++) {
		printf("%ld: ", i);
		d = (PERMUTATION_OP) D->s_i(i);
		/* d->println();
		d->invers(&dv); */
		dc_get_k_set(d, &R, 8, FALSE);
		reduce_generators_labra(Dc->s_Ad_i(i), &stab_go, 
			FALSE /* f_verbose */, &labra_A1);
		dc_print_k_set(&R, &stab_go, &go);
		// labra_A.save("M24_labra.dsc");
		// labra_A1.save("MOG_labra.dsc");
		if (i == 2) {
			printf("d=\n");
			d->print_list();
			dc_calc_set_orbit(&A, &R, W24_blocks, TRUE /* f_v */, TRUE /* f_vv */);
			R.quicksort(R.s_li(), TRUE);
			W24_blocks->search(W24_blocks->s_li(), TRUE, &R, &idx, &f_found);
			if (!f_found)
				return error("construct_W24() original block not found");
			*block_idx = idx - 1;
			labra_A1.reduced_generating_set(stab_gen, 
				FALSE /* f_bottom_up */, TRUE /* f_v */);
			// vec_conjugate(stab_gen, d);
			} /* if i */
		}
	if (f_v) {
		printf("the steiner system S(5,8,24) (a 5-(24,8,1) design):\n");
		for (i = 0; i < W24_blocks->s_li(); i++) {
			printf("block %ld: ", i);
			W24_blocks->s_i(i)->println();
			}
		}
	return OK;
}

#if TEXDOCU
INT construct_W23(VECTOR_OP W23_blocks, 
	VECTOR_OP M23_gen, VECTOR_OP stab_gen, INT *block_idx, INT f_v)
#endif
{
	VECTOR_OB T, A;
	SYM_OB id;
	DCY_OB dc[MAX_STEP];
	DCY_OP Dc;
	INT erg = OK, i, l;
	INT up_to_k = 7, up_to_step;
	INT type;
	void *data;
	INT deg;
	VECTOR_OP D;
	SYM_OB go, stab_go;
	PERMUTATION_OP d;
	VECTOR_OB R;
	LABRA_OB labra_A, labra_A1;
	INT idx, f_found;
	
	up_to_step = 2 * up_to_k - 1;
	if (up_to_step >= MAX_STEP)
		return error("construct_W23() up_to_step too large");
	deg = 23;
	
	type = DO_TYPE_SYM;
	data = NIL;
	erg += M23_generators(&A);
	// A.save("M23_generators.dsc");
	reduce_generators_labra(&A, &go, FALSE /* f_verbose */, &labra_A);
	printf("group order: ");
	go.println();
	A.copy(M23_gen);
#if 0
	labra_A.reduced_generating_set(M23_gen, 
		FALSE /* f_bottom_up */, TRUE /* f_v */);
#endif

	do_copy(A.s_i(0), &id, type, data);
	do_one(&id, type, data);
	for (i = 0; i <= up_to_step; i++) {
		dc[i].initialize_Young(&T, deg, i /* step */, 
			0 /* type */, NIL /* data */);
		T.swap(dc[i].s_T());
		if (i <= up_to_step)
			dc[i].s_reduce_generators_mode()->m_i(2);
				/* 0: no reduction 
				 * 1: dimino 
				 * 2: labra */
		}
		
	dc[0].s_D()->m_il(1);
	do_copy(&id, dc[0].s_D_i(0), type, data);
	
	dc[0].s_Ad()->m_il(1);
	((SYM_OP) &A)->copy(dc[0].s_Ad()->s_i(0));
	dc[0].s_D()->println();
	dc[0].s_Ad()->println();

	for (i = 1; i <= up_to_step; i++)
		erg += dc_do_step(&dc[0], i, &id, deg /* deg */, 
			TRUE /* f_verbose */, type, data);

	Dc = &dc[up_to_step];
	D = Dc->s_D();
	l = D->s_li();
	printf("found %ld double cosets\n", l);
	for (i = 0; i < l; i++) {
		printf("%ld: ", i);
		d = (PERMUTATION_OP) D->s_i(i);
		/* d->println();
		d->invers(&dv); */
		dc_get_k_set(d, &R, 7, FALSE);
		reduce_generators_labra(Dc->s_Ad_i(i), &stab_go, 
			FALSE /* f_verbose */, &labra_A1);
		dc_print_k_set(&R, &stab_go, &go);
		if (i == 3) {
			printf("d=\n");
			d->print_list();
			dc_calc_set_orbit(&A, &R, W23_blocks, TRUE /* f_v */, TRUE /* f_vv */);
			R.quicksort(R.s_li(), TRUE);
			W23_blocks->search(W23_blocks->s_li(), TRUE, &R, &idx, &f_found);
			if (!f_found)
				return error("construct_W23() original block not found");
			*block_idx = idx - 1;
			labra_A1.reduced_generating_set(stab_gen, 
				FALSE /* f_bottom_up */, TRUE /* f_v */);
			// vec_conjugate(stab_gen, d);
			} /* if i */
		}
	if (f_v) {
		printf("the steiner system S(4,7,23) (a 4-(23,7,1) design):\n");
		for (i = 0; i < W23_blocks->s_li(); i++) {
			printf("block %ld: ", i);
			W23_blocks->s_i(i)->println();
			}
		}
	return OK;
}

#if TEXDOCU
INT construct_W22(VECTOR_OP W22_blocks)
#endif
{
	VECTOR_OB T, A;
	SYM_OB id;
	DCY_OB dc[MAX_STEP];
	DCY_OP Dc;
	INT erg = OK, i, l;
	INT up_to_k = 6, up_to_step;
	INT type;
	void *data;
	INT deg;
	VECTOR_OP D;
	SYM_OB go, stab_go, go1;
	PERMUTATION_OP d;
	VECTOR_OB R;
	LABRA_OB labra_A, labra_A1;
	VECTOR_OB O1, O2;
	VECTOR_OB gsel;
	GROUP_SELECTION_OP gs;
	BYTE g_label[1024];
	BYTE g_label_tex[1024];
	
	up_to_step = 2 * up_to_k - 1;
	if (up_to_step >= MAX_STEP)
		return error("construct_W22() up_to_step too large");
	deg = 22;
	
	type = DO_TYPE_SYM;
	data = NIL;
	
	// A := PGGL(3,4)+ of degree 21 + 1 = 22
	gsel.m_il(10);
	gs = (GROUP_SELECTION_OP) gsel.s_i(0);
	gs->init(FGA_GROUP_PSSL, 3, 4, NIL);
	gs = (GROUP_SELECTION_OP) gsel.s_i(1);
	gs->init(FGA_ADD_FIXPOINT, 0, 0, NIL);
	km_get_group_from_selection(&gsel, 2 /* nb_gsel */, 
		&A, g_label, g_label_tex);
	printf("got group %s\n", g_label);
	fflush(stdout);
	reduce_generators_labra(&A, &go, FALSE /* f_verbose */, &labra_A);
	printf("group order: ");
	go.println();

	do_copy(A.s_i(0), &id, type, data);
	do_one(&id, type, data);
	for (i = 0; i <= up_to_step; i++) {
		dc[i].initialize_Young(&T, deg, i /* step */, 
			0 /* type */, NIL /* data */);
		T.swap(dc[i].s_T());
		if (i <= up_to_step)
			dc[i].s_reduce_generators_mode()->m_i(2);
				/* 0: no reduction 
				 * 1: dimino 
				 * 2: labra */
		}
		
	dc[0].s_D()->m_il(1);
	do_copy(&id, dc[0].s_D_i(0), type, data);
	
	dc[0].s_Ad()->m_il(1);
	((SYM_OP) &A)->copy(dc[0].s_Ad()->s_i(0));
	dc[0].s_D()->println();
	dc[0].s_Ad()->println();

	for (i = 1; i <= up_to_step; i++)
		erg += dc_do_step(&dc[0], i, &id, deg /* deg */, 
			TRUE /* f_verbose */, type, data);

	Dc = &dc[up_to_step];
	D = Dc->s_D();
	l = D->s_li();
	printf("found %ld double cosets\n", l);
	for (i = 0; i < l; i++) {
		printf("%ld: ", i);
		d = (PERMUTATION_OP) D->s_i(i);
		/* d->println(); */
		dc_get_k_set(d, &R, up_to_k, FALSE);
		reduce_generators_labra(Dc->s_Ad_i(i), &stab_go, 
			FALSE /* f_verbose */, &labra_A1);
		dc_print_k_set(&R, &stab_go, &go);
		// labra_A.save("M24_labra.dsc");
		// labra_A1.save("MOG_labra.dsc");
		if (i == 4) {
			// dc_get_transversals(dc, &go, deg, up_to_step, i);
			dc_calc_set_orbit(&A, &R, &O1, TRUE /* f_v */, TRUE /* f_vv */);
			}
		if (i == 16) {
			dc_calc_set_orbit(&A, &R, &O2, TRUE /* f_v */, TRUE /* f_vv */);
			}
		}
	O1.append(&O2, W22_blocks);
	return OK;
}

#if TEXDOCU
INT construct_graph_from_design(INT v, VECTOR_OP blocks, MATRIX_OP G, INT i1)
#endif
{
	INT vb;
	VECTOR_OP B, B1, B2;
	INT i, j, k, kk, l, ii, a, b, inter;
	
	b = blocks->s_li();
	vb = v + b;
	G->m_ilih_n(vb, vb);

	for (i = 0; i < v; i++) {
		for (j = i + 1; j < v; j++) {
			G->m_iji(i, j, 1);
			G->m_iji(j, i, 1);
			}
		}
	for (k = 0; k < b; k++) {
		B = (VECTOR_OP) blocks->s_i(k);
		l = B->s_li();
		for (ii = 0; ii < l; ii++) {
			a = B->s_ii(ii);
			// block k contains point a
			i = a;
			j = v + k;
			G->m_iji(i, j, 1);
			G->m_iji(j, i, 1);
			}
		
		}
	for (k = 0; k < b; k++) {
		for (kk = k + 1; kk < b; kk++) {
			B1 = (VECTOR_OP) blocks->s_i(k);
			B2 = (VECTOR_OP) blocks->s_i(kk);
			inter = block_intersection(B1, B2);
			if (inter == i1) {
				i = v + k;
				j = v + kk;
				G->m_iji(i, j, 1);
				G->m_iji(j, i, 1);
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT construct_two_graph(MATRIX_OP G, VECTOR_OP three_blocks)
#endif
{
	INT choice[64];
	INT n, k = 3, a, b, s, i, j, l = 0;
	VECTOR_OB B;

	n = G->s_hi();
	if (G->s_hi() != G->s_li())
		return error("construct_two_graph() not symmetric");
	three_blocks->m_il(0);
	
	n_Choose_k_first(choice, n, k);
	do {
		s = 0;
		for (i = 0; i < k; i++) {
			a = choice[i];
			for (j = i + 1; j < k; j++) {
				b = choice[j];
				if (G->s_iji(a, b))
					s++;
				}
			}
		if (ODD(s)) {
			B.m_il(k);
			for (i = 0; i < k; i++) {
				B.m_ii(i, choice[i]);
				}
			three_blocks->inc();
			B.swap(three_blocks->s_i(l));
			l++;
			if (l % 10 == 0)
				printf(",");
			else
				printf(".");
			if (l % 50 == 0)
				printf("%ld (%ld,%ld,%ld)\n", l, choice[0], choice[1], choice[2]);
			fflush(stdout);
			}
		} while (n_Choose_k_next(choice, n, k));
	return OK;
}

#if TEXDOCU
INT get_graph_value(MATRIX_OP G, INT i, INT j)
#endif
{
	// INT n = G->s_li();

	if (i < j)
		return G->s_iji(i, j);
	else if (i == j)
		return 0;
	else
		return G->s_iji(j, i);
}

#if TEXDOCU
INT blocks_disjoint(VECTOR_OP B1, VECTOR_OP B2)
#endif
{
	INT i, l, idx, f_found;
	
	l = B1->s_li();
	for (i = 0; i < l; i++) {
		B2->search(B2->s_li(), TRUE, B1->s_i(i), &idx, &f_found);
		if (f_found)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT block_intersection(VECTOR_OP B1, VECTOR_OP B2)
#endif
{
	INT n, i, l, idx, f_found;
	
	n = 0;
	l = B1->s_li();
	for (i = 0; i < l; i++) {
		B2->search(B2->s_li(), TRUE, B1->s_i(i), &idx, &f_found);
		if (f_found)
			n++;
		}
	return n;
}

#if TEXDOCU
INT McKayWreath()
#endif
/* S_2 \wreath S_12 \ S_{24} / Mathieu 24 */
/* WARNING !
 * this program takes 60 h on a 
 * HP with 400 MB RAM !
 * it uses at least 300 MB RAM.
 * the number of double cosets are:
step|index to step -1 | # double| computing
    |(<0 means        | cosets  | time
INT McKayWreath()
/* S_2 \wreath S_12 \ S_{24} / Mathieu 24 */
/* WARNING !
 * this program takes 60 h on a 
 * HP with 400 MB RAM !
 * it uses at least 300 MB RAM.
 * the number of double cosets are:
step|index to step -1 | # double| computing
    |(<0 means        | cosets  | time
    | downstep)       |         | (in h)
 0     1
 1 -24 1
 2 -23 1
 3   2 1
 4   1 1
 5 -22 1
 6 -21 1
 7   2 1
 8   2 1
 9 -20 1
10 -19 1
11   2 2
12   3 2
13 -18 3
14 -17 11
15   2 10
16   4 7
17 -16 18
18 -15 119
19   2 84
20   5 31
21 -14 153
22 -13 1528
23   2 867
24   6 207
25 -12 1520
26 -11 15615
27   2 8184
28   7 1332
29 -10 11192
30  -9 98518
31   2 50248
32   8 6723
33  -8 49396
34  -7 342780
35   2 173421
36   9 20079
37  -6 114515 21
38  -5 570206 38
39   2 287954 45:30
40  10 29922
41  -4 114360
42  -3 342229
43   2 173728 59
44  11 16764 59:30
45  -2 31310
46  -1 31310
47   2 16764
48  12 1858 60
*/ 
{
	VECTOR_OB T, A;
	SYM_OB id;
	DCY_OB dc[MAX_STEP];
	INT erg = OK, i;
	INT up_to_k = 12, up_to_step;
	INT type;
	void *data;
	INT deg;
	
	up_to_step = 4 * up_to_k;
	if (up_to_step >= MAX_STEP)
		return error("up_to_step too large");
	deg = 24;
	
	type = DO_TYPE_SYM;
	data = NIL;
	erg += M24_generators(&A);

	do_copy(A.s_i(0), &id, type, data);
	do_one(&id, type, data);
	for (i = 0; i <= up_to_step; i++) {
		dc[i].initialize_Wreath(&T, deg, i /* step */, 0 /* type */, NIL /* data */);
		T.swap(dc[i].s_T());
		if (i <= 4)
			dc[i].s_reduce_generators_mode()->m_i(2);
				/* 0: no reduction 
				 * 1: dimino 
				 * 2: labra */
		}
		
	dc[0].s_D()->m_il(1);
	do_copy(&id, dc[0].s_D_i(0), type, data);
	
	dc[0].s_Ad()->m_il(1);
	((SYM_OP) &A)->swap(dc[0].s_Ad()->s_i(0));
	dc[0].s_D()->println();
	dc[0].s_Ad()->println();

	for (i = 1; i <= up_to_step; i++)
		erg += dc_do_step(&dc[0], i, &id, deg /* deg */, TRUE /* f_verbose */, type, data);

	return OK;
}

#if TEXDOCU
INT virasoro(INT f_v)
#endif
{
	VECTOR_OB O, Nk, R2, M24_gen, stab_gen;
	VECTOR_OP R1;
	INT j, k, jj, j0, O_len, block_idx;
	
	construct_W24(&O, &M24_gen, &stab_gen, &block_idx, f_v);
	O_len = O.s_li();
	Nk.m_il_n(5);
	for (j = 0; j < O_len; j++) {
		R1 = (VECTOR_OP) O.s_i(j);
		R2.m_il_n(24);
		for (jj = 0; jj < 8; jj++) {
			k = R1->s_ii(jj);
			R2.m_ii(k, 1);
			}
		k = 0;
		for (jj = 0; jj < 12; jj++) {
			j0 = jj << 1;
			if (R2.s_ii(j0) && R2.s_ii(j0 + 1))
				k++;
			}
		Nk.s_i(k)->inc();
		}
	printf("n[k]: ");
	Nk.println();
	return OK;
}

#endif /* LADDER_TRUE */

