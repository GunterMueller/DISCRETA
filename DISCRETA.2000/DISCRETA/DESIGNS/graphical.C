/* graphical.C */

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
#include <DISCRETA/bruch.h>

#define MAX_STEP 64

// #define GRAPHICAL_MATRIX_NROW 50
// #define GRAPHICAL_MATRIX_NCOL 9

#if TEXDOCU
INT graphical_design(INT t, INT k, INT print_nrow, INT print_ncol)
#endif
{
	INT n, n2;
	VECTOR_OB T, A;
	SYM_OB id;
	DCY_OB dc[MAX_STEP];
	INT erg = OK, i;
	INT up_to_k = k, up_to_step;
	INT type;
	void *data;
	INT deg;
	SYM_OB go, stab_go, go1;
	VECTOR_OB R;
	LABRA_OB labra_A, labra_A1;
	VECTOR_OB O1, O2;
	VECTOR_OB gsel;
	GROUP_SELECTION_OP gs;
	BYTE g_label[1024];
	BYTE g_label_tex[1024];
	BYTE km_fname[1024];
	MATRIX_OB M, N, N1, N2;
	VECTOR_OB stab_order_t, supp_t, supp_k;
	VECTOR_OB stab_order_k, aut_t, aut_k;
	INT width_in_pages;
	
	n = 2 * k;
	n2 = (n * (n - 1)) / 2;
	
	up_to_step = 2 * up_to_k - 1;
	if (up_to_step >= MAX_STEP)
		return error("graphical_design() up_to_step too large");
	deg = n2;
	
	type = DO_TYPE_SYM;
	data = NIL;
	
	// A := S_k^[2]
	gsel.m_il(10);
	gs = (GROUP_SELECTION_OP) gsel.s_i(0);
	gs->init(FGA_GROUP_SYM, n, 0, NIL);
	gs = (GROUP_SELECTION_OP) gsel.s_i(1);
	gs->init(FGA_INDUCE_2_SETS, 0, 0, NIL);
	km_get_group_from_selection(&gsel, 2 /* nb_gsel */, 
		&A, g_label, g_label_tex);
	printf("got group %s\n", g_label);
	fflush(stdout);
	reduce_generators_labra(&A, &go, FALSE /* f_verbose */, &labra_A);
	printf("group order: ");
	go.println();

	km_init_file_names(g_label, km_fname, t, k);
	
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
	
	dc_Mtk(dc, t, k, &M);
	
	km_print_M_asc(g_label, g_label_tex, km_fname, 
		&A, &go, deg, 
		&M, t, k, FALSE /*  f_k2 */, 0 /* k2 */);
	
	graphical_orbits(dc, &go, t, n, &stab_order_t, &supp_t, &aut_t);
	printf("stab_order %ld orbits: ", t);
	stab_order_t.println();
	printf("support %ld orbits: ", t);
	supp_t.println();
	printf("ago %ld orbits: ", t);
	aut_t.println();
	graphical_orbits(dc, &go, k, n, &stab_order_k, &supp_k, &aut_k);
	printf("stab_order %ld orbits: ", k);
	stab_order_k.println();
	printf("support %ld orbits: ", k);
	supp_k.println();
	printf("ago %ld orbits: ", k);
	aut_k.println();

	width_in_pages = aut_k.s_li() / 5;
	if (width_in_pages <= 0)
		width_in_pages = 1;
	
	calc_graphical_KM_matrix(&M, &N, &N1, &N2, n, 
		&stab_order_t, &supp_t, &aut_t, 
		&stab_order_k, &supp_k, &aut_k);
	
#if 0
	print_graphical_KM_matrix_splitted(&N, &N1, &N2, t, k, 
		&supp_t, &aut_t, 
		&supp_k, &aut_k, 
		GRAPHICAL_MATRIX_NROW,
		GRAPHICAL_MATRIX_NCOL);
#endif

	print_graphical_KM_matrix_splitted(&N, &N1, &N2, t, k, 
		&supp_t, &aut_t, 
		&supp_k, &aut_k, 
		print_nrow,
		print_ncol);

#if 0
	sprintf(fname_ps, "g_t%ld_k%ld.ps", t, k);
	dcl = open_dcl(dc, n, n2, 
		up_to_step, FALSE /* f_with_perm */, 
	 	FALSE /* f_verbose */, width_in_pages, 
		fname_ps, type, data);

	free_dcl(dcl);
	dcl = NIL;
#endif

	
	return OK;
}

#if TEXDOCU
INT print_graphical_KM_matrix(MATRIX_OP N, MATRIX_OP N1, MATRIX_OP N2, INT t, INT k, 
	VECTOR_OP supp_t, VECTOR_OP aut_t, 
	VECTOR_OP supp_k, VECTOR_OP aut_k)
#endif
{
	INT l, h, i, j, a, b;
	SYM_OP q;
	BYTE str[256];

	l = N->s_li();
	h = N->s_hi();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			N->s_ij(i, j)->print();
			printf(" [n-%ld]_%ld\n", N1->s_iji(i, j), N2->s_iji(i, j));
			}
		}

	printf("\\begin{tabular}{|c|c|c||*{%ld}{c}", h);
	// for (j = 0; j < h; j++) 
		// printf("c");
	printf("|}\n");
	printf("\\hline\n");
	
	printf("$%ld \\backslash %ld$ & & & ", k, t);
	for (j = 0; j < h; j++) {
		printf(" %ld ", j + 1);
		if (j < h - 1)
			printf(" & ");
		}
	printf("\\\\\n");
	printf("\\hline\n");
	
	printf(" & $|S|$ & & ");
	for (j = 0; j < h; j++) {
		printf(" %ld ", supp_t->s_ii(j));
		if (j < h - 1)
			printf(" & ");
		}
	printf("\\\\\n");
	printf("\\hline\n");
	
	printf(" & & $|A_S|$ & ");
	for (j = 0; j < h; j++) {
		printf(" %ld ", aut_k->s_ii(j));
		if (j < h - 1)
			printf(" & ");
		}
	printf("\\\\\n");
	printf("\\hline\n");
	printf("\\hline\n");
	
	for (i = 0; i < l; i++) {
		printf("%ld & %ld & %ld & ", i + 1, supp_k->s_ii(i), aut_k->s_ii(i));
		for (j = 0; j < h; j++) {
			q = N->s_ij(j, i);
			if (q->nullp()) {
				printf("0 ");
				}
			else {
				str[0] = 0;
				q->sprint(str);
				printf(" $%s ", str);
				a = supp_t->s_ii(j);
				b = supp_k->s_ii(i) - a;
				if (b > 0) {
					if (a > 0) {
						printf("[n-%ld]", a);
						}
					else {
						printf("[n]");
						}
					if (b >= 10) {
						printf("_{%ld}", b);
						}
					else if (b > 1) {
						printf("_%ld", b);
						}
					else if (b > 0) {
						;
						}
					}
				else {
					}
				printf("$");
				}
			if (j < h - 1)
				printf(" & ");
			}
		printf("\\\\\n");
		}
	
	printf("\\hline\n");
	printf("\\end{tabular}\n\n\n");
	return OK;
}

#if TEXDOCU
INT print_graphical_KM_matrix_splitted(MATRIX_OP N, MATRIX_OP N1, MATRIX_OP N2, INT t, INT k, 
	VECTOR_OP supp_t, VECTOR_OP aut_t, 
	VECTOR_OP supp_k, VECTOR_OP aut_k, 
	INT nrow, INT ncol)
#endif
{
	INT l, h, i, j, a, b, i0, j0, height, width;
	SYM_OP q;
	BYTE str[256];

	l = N->s_li();
	h = N->s_hi();
#if 0
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			N->s_ij(i, j)->print();
			printf(" [n-%ld]_%ld\n", N1->s_iji(i, j), N2->s_iji(i, j));
			}
		}
#endif

#if 0
	if (l < ncol)
		return OK;
#endif
	for (i0 = 0; i0 < l; i0 += nrow) {
		if (i0 + nrow > l)
			height = l - i0;
		else
			height = nrow;

		for (j0 = 0; j0 < h; j0 += ncol) {
			if (j0 + ncol > h)
				width = h - j0;
			else
				width = ncol;

	printf("%% i0 = %ld j0 = %ld\n", i0, j0);
	printf("\\begin{tabular}{|c|c|c||*{%ld}{c}", width);
	printf("|}\n");
	printf("\\hline\n");
	
	printf("$%ld \\backslash %ld$ & & & ", k, t);
	for (j = 0; j < width; j++) {
		printf(" %ld ", j0 + j + 1);
		if (j < width - 1)
			printf(" & ");
		}
	printf("\\\\\n");
	printf("\\hline\n");
	
	printf(" & $|S|$ & & ");
	for (j = 0; j < width; j++) {
		printf(" %ld ", supp_t->s_ii(j0 + j));
		if (j < width - 1)
			printf(" & ");
		}
	printf("\\\\\n");
	printf("\\hline\n");
	
	printf(" & & $|A_S|$ & ");
	for (j = 0; j < width; j++) {
		printf(" %ld ", aut_k->s_ii(j0 + j));
		if (j < width - 1)
			printf(" & ");
		}
	printf("\\\\\n");
	printf("\\hline\n");
	printf("\\hline\n");
	
	for (i = 0; i < height; i++) {
		printf("%ld & %ld & %ld & ", i0 + i + 1, supp_k->s_ii(i0 + i), aut_k->s_ii(i0 + i));
		for (j = 0; j < width; j++) {
			q = N->s_ij(j0 + j, i0 + i);
			if (q->nullp()) {
				printf("0 ");
				}
			else {
				str[0] = 0;
				q->sprint(str);
				printf(" $%s ", str);
				a = supp_t->s_ii(j0 + j);
				b = supp_k->s_ii(i0 + i) - a;
				if (b > 0) {
					if (a > 0) {
						printf("[n-%ld]", a);
						}
					else {
						printf("[n]");
						}
					if (b >= 10) {
						printf("_{%ld}", b);
						}
					else if (b > 1) {
						printf("_%ld", b);
						}
					else if (b > 0) {
						;
						}
					}
				else {
					}
				printf("$");
				}
			if (j < width - 1)
				printf(" & ");
			}
		printf("\\\\\n");
		}
	printf("\\hline\n");
	printf("\\end{tabular}\n\n\n");
	
		} // next i0
		} // next j0
	
	return OK;
}

#if TEXDOCU
INT calc_graphical_KM_matrix(MATRIX_OP M, 
	MATRIX_OP N, MATRIX_OP N1, MATRIX_OP N2, INT n, 
	VECTOR_OP stab_order_t, VECTOR_OP supp_t, VECTOR_OP aut_t, 
	VECTOR_OP stab_order_k, VECTOR_OP supp_k, VECTOR_OP aut_k)
#endif
{
	INT i, j, l, h, o, u;
	SYM_OP m_op;
	SYM_OB m1, tmp;
	BRUCH_OB q;

	h = M->s_hi();
	l = M->s_li();
	N->m_ilih(l, h);
	N1->m_ilih(l, h);
	N2->m_ilih(l, h);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			m_op = M->s_ij(i, j);
			m_op->mult(stab_order_k->s_i(j), &tmp);
			tmp.ganzdiv(stab_order_t->s_i(i), &m1);
			m1.mult(aut_t->s_i(i), &tmp);
			if (tmp.s_obj_k() != INTEGER)
				return error("calc_graphical_KM_matrix() not an integer 1");
			o = ((INTEGER_OP) &tmp)->s_i();
			if (aut_k->s_i(j)->s_obj_k() != INTEGER)
				return error("calc_graphical_KM_matrix() not an integer 2");
			u = ((INTEGER_OP) aut_k->s_i(j))->s_i();
			q.m_ioiu(o, u);
			q.kuerzen();
			((SYM_OP) &q)->copy(N->s_ij(i, j));
			N1->m_iji(i, j, supp_t->s_ii(i));
			N2->m_iji(i, j, supp_k->s_ii(j) - supp_t->s_ii(i));
			}
		}
	return OK;
}

#if TEXDOCU
INT graphical_orbits(DCY_OB *dc, SYM_OP go, INT k, INT n, 
	VECTOR_OP stab_order, VECTOR_OP supp, VECTOR_OP aut)
#endif
{
	INT up_to_step = 2 * k - 1;
	DCY_OP Dc;
	VECTOR_OP D;
	PERMUTATION_OP d;
	SYM_OB stab_go, go1, ago, ago1;
	LABRA_OB labra_A1;
	VECTOR_OB R;
	INT i, l;
	INT type;
	void *data;
	VECTOR_OB support;
	INTEGER_OB int_ob;
	INT j, s;
	
	type = DO_TYPE_SYM;
	data = NIL;
	
	Dc = &dc[up_to_step];
	D = Dc->s_D();
	l = D->s_li();
	printf("found %ld double cosets on %ld sets\n", l, k);
	supp->m_il(l);
	aut->m_il(l);
	stab_order->m_il(l);
	for (i = 0; i < l; i++) {
		printf("%ld: ", i);
		d = (PERMUTATION_OP) D->s_i(i);
		/* d->println(); */
		dc_get_k_set(d, &R, k, FALSE);
		reduce_generators_labra(Dc->s_Ad_i(i), &stab_go, 
			FALSE /* f_verbose */, &labra_A1);
		dc_print_k_set(&R, &stab_go, go);
		stab_go.copy(stab_order->s_i(i));
		support_on_2_sets(&R, n, &support, &s);
		printf("support = %ld\n", s);
		supp->m_ii(i, s);
		stab_go.copy(&ago);
		for (j = 0; j < n - s; j++) {
			int_ob.m_i(j + 1);
			ago.ganzdiv(&int_ob, &ago1);
			ago1.swap(&ago);
			}
		printf("aut=");
		ago.println();
		ago.swap(aut->s_i(i));
		reduce_generators_labra(Dc->s_Ad_i(i), &stab_go, 
			FALSE /* f_verbose */, &labra_A1);
		printf("stab_go (2nd time) = ");
		stab_go.println();
		}
	return OK;
}

#endif /* LADDER_TRUE */

/* end of graphical.C */
