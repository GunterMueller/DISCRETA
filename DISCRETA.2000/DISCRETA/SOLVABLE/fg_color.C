/* fg_color.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>

#if 0
INT FG_OB::order_structure(VECTOR_OP ord)
{
	INT i, n, o;
	
	n = s_n_i();
	ord->m_il(n);
	for (i = 0; i < n; i++) {
		gt_order(i, &o);
		ord->m_ii(i, o);
		}
	return OK;
}
#endif

#if TEXDOCU
INT FG_OB::order_structure(VECTOR_OP ord)
#endif
{
	GROUP_TABLE table;

	group_init_from_fg(&table, this);
	group_order_structure(&table, ord);
	group_free(&table);
	return OK;
}

#if 0
INT FG_OB::is_central(INT i)
{
	ZE_OP ze;
	INT j, g_j, r;
	
	r = s_nb_ze_i();
	for (j = 0; j < r; j++) {
		ze = s_ze_i(j);
		g_j = ze->s_n0_i();
		if (gt_conj(i, g_j) != i)
			return FALSE;
		}
	return TRUE;
}
#endif

#if 0
INT FG_OB::center(VECTOR_OP p)
{
	INT i, n, j;
	
	n = s_n_i();
	p->m_il(n);
	for (i = 0; i < n; i++) {
		j = is_central(i);
		p->m_ii(i, j);
		}
	return OK;
}
#endif

#if TEXDOCU
INT FG_OB::center(VECTOR_OP p)
#endif
{
	GROUP_TABLE table;

	group_init_from_fg(&table, this);
	group_center(&table, p);
	group_free(&table);
	return OK;
}

#if 0
INT FG_OB::classes(VECTOR_OP class_no, 
	VECTOR_OP sv_last, VECTOR_OP sv_gen, 
	VECTOR_OP class_rep, VECTOR_OP class_len, VECTOR_OP class_el_ord)
{
	INT nb_generators, gen_i, i, n, j, nb_classes = 0, first, os, i0, i1, ord;
	INT *Q = NIL, nb_Q;
	
	n = s_n_i();
	nb_generators = nb_gen();
	Q = (INT *) my_malloc(n * sizeof(INT), "fg.C: classes");
	class_no->m_il(n);
	sv_last->m_il(n);
	sv_gen->m_il(n);
	class_rep->m_il(0);
	class_len->m_il(0);
	class_el_ord->m_il(0);
	for (i = 0; i < n; i++) {
		class_no->m_ii(i, -1);
		sv_last->m_ii(i, -1);
		sv_gen->m_ii(i, -1);
		}
	for (first = 0; first < n; first++) {
		if (class_no->s_ii(first) != -1)
			continue;
		/* now first is an unprocessed element, 
		 * we define a new class. */
		os = 1;
		Q[0] = first;
		nb_Q = 1;
		class_no->m_ii(first, nb_classes); /* the class no */
		gt_order(first, &ord);
		while (nb_Q) {
			i0 = Q[0];
			for (i = 1; i < nb_Q; i++)
				Q[i - 1] = Q[i];
			nb_Q--;
			/* apply all generators to i0: */
			for (i = 0; i < nb_generators; i++) {
				gen_i = g_i(i);
				i1 = gt_conj(i0, gen_i);
				if (class_no->s_ii(i1) == -1) { /* new orbit element */
					sv_last->m_ii(i1, i0);
					sv_gen->m_ii(i1, i);
					class_no->m_ii(i1, nb_classes);
					Q[nb_Q++] = i1;
					os++;
					}
				} /* next i */
			} /* while nb_Q */
		class_rep->inc();
		class_len->inc();
		class_el_ord->inc();
		class_rep->m_ii(nb_classes, first);
		class_len->m_ii(nb_classes, os);
		class_el_ord->m_ii(nb_classes, ord);
		nb_classes++;
		} /* next first (next class) */
	my_free(Q);
	return OK;
}

INT FG_OB::print_classes(VECTOR_OP class_no, VECTOR_OP class_rep, VECTOR_OP class_len)
{
	INT n, i, j, o, first, nb_classes;
	
	n = s_n_i();
	nb_classes = class_rep->s_li();
	for (i = 0; i < nb_classes; i++) {
		first = class_rep->s_ii(i);
		printf("class %ld: { %ld ", first);
		for (j = first + 1; j < n; j++) {
			if (class_no->s_ii(j) == i) {
				printf("%ld ", j);
				}
			}
		gt_order(first, &o);
		printf("} len %ld (el.-order %ld)\n", class_len->s_ii(i), o);
		}
	return OK;
}

INT FG_OB::first_class_coloring(
	VECTOR_OP class_no, VECTOR_OP class_rep, 
	VECTOR_OP class_len, VECTOR_OP class_el_ord, 
	VECTOR_OP col, MATRIX_OP Col, INT f_v)
{
	VECTOR_OB cl, ceo;
	INT n, i, idx, a, b;
	
	n = s_n_i();
	cl.m_il(n);
	ceo.m_il(n);
	for (i = 0; i < n; i++) {
		idx = class_no->s_ii(i);
		a = class_len->s_ii(idx);
		b = class_el_ord->s_ii(idx);
		cl.m_ii(i, a);
		ceo.m_ii(i, b);
		}
	color_join(&cl, &ceo, col, Col, FALSE /* f_v */, FALSE /* f_vv */);
	if (f_v) {
		VECTOR_OB val, mult;
		BYTE str[10000];
		INT c, nc, a, b;
		
		col->multiplicities(&val, &mult);
		str[0] = 0;
		val.sprint_multiplicities(&mult, str);
		printf("first class coloring: %s\n", str);
		nc = Col->s_hi();
		for (c = 0; c < nc; c++) {
			a = Col->s_iji(c, 0);
			b = Col->s_iji(c, 1);
			printf("color %ld means class-length %ld / "
				"element-order %ld\n", c, a, b);
			}
		}
	return OK;
}
#endif

#if 0
INT FG_OB::Classes(BYTE *str, INT f_v, INT f_vv)
{
	VECTOR_OB class_no;
	VECTOR_OB sv_last, sv_gen;
	VECTOR_OB class_rep, class_len, class_el_ord;
	INT n;

	n = s_n_i();
	classes(&class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		print_classes(&class_no, &class_rep, &class_len);
		}
	sprint_class_structure(&class_len, &class_el_ord, str, f_v);
	return OK;
}
#endif

#if TEXDOCU
INT FG_OB::Classes(BYTE *str, INT f_v, INT f_vv)
#endif
{
	GROUP_TABLE table;
	VECTOR_OB class_no;
	VECTOR_OB sv_last, sv_gen;
	VECTOR_OB class_rep, class_len, class_el_ord;

	if (f_v) {
		printf("classes of group %s:\n", s_label_s());
		fflush(stdout);
		}
	group_init_from_fg(&table, this);
	group_classes(&table, &class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		group_print_classes(&table, 
			&class_no, &class_rep, &class_len);
		}
	sprint_class_structure(&class_len, &class_el_ord, str, f_v);
	group_free(&table);
	return OK;
}

#if TEXDOCU
INT FG_OB::Coloring_seldoms_prefered(VECTOR_OP col)
#endif
{
	VECTOR_OB col0;
	PERMUTATION_OB p1, p2;
	
	Coloring(&col0, FALSE, FALSE);
	prefer_seldom_colors(&col0, col, &p1, &p2);
	return OK;
}

#if 0
INT FG_OB::Coloring(VECTOR_OP col, INT f_v, INT f_vv)
{
	VECTOR_OB class_no;
	VECTOR_OB sv_last, sv_gen;
	VECTOR_OB class_rep, class_len, class_el_ord;
	MATRIX_OB Col;
	BYTE str[10000];
	INT n;

	n = s_n_i();
	classes(&class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		print_classes(&class_no, &class_rep, &class_len);
		}
	str[0] = 0;
	sprint_class_structure(&class_len, &class_el_ord, str, f_v);
	
	/* first class coloring 
	 * (according to class length / element order): */
	first_class_coloring(&class_no, 
		&class_rep, &class_len, &class_el_ord, 
		col, &Col, f_v /* f_v */);
	
	refine_by_roots(col, &Col, f_v /* f_v */);
	
	return OK;
}
#endif

#if TEXDOCU
INT FG_OB::Coloring(VECTOR_OP col, INT f_v, INT f_vv)
#endif
{
	GROUP_TABLE table;

	if (f_v) {
		printf("Coloring for group %s:\n", s_label_s());
		fflush(stdout);
		}
	group_init_from_fg(&table, this);
	group_Coloring(&table, col, f_v, f_vv);
	group_free(&table);
	return OK;
}

#if 0
INT FG_OB::refine_by_roots(VECTOR_OP col, MATRIX_OP Col, INT f_v)
{
	VECTOR_OB vp, ve;
	BYTE str[10000];
	INT n, i, l, p, nc0, nc;
	
	if (f_v)
		printf("refine_by_roots():\n");
	n = s_n_i();
	factor_integer(n, &vp, &ve);
	if (f_v) {
		str[0] = 0;
		print_factorization(&vp, &ve, str);
		printf("n = %ld = %s\n", n, str);
		}
	l = vp.s_li();
	while (TRUE) {
		nc0 = number_of_colors(Col, f_v);
		for (i = 0; i < l; i++) {
			p = vp.s_ii(i);
			refine_by_p_roots(p, col, Col, f_v /* f_v */);
			}
		nc = number_of_colors(Col, f_v);
		if (nc == nc0)
			break;
		}
	return OK;
}

INT FG_OB::refine_by_p_roots(INT p, VECTOR_OP col, MATRIX_OP Col, INT f_v)
{
	VECTOR_OB Pr, pr, val, mult;
	BYTE str[10000];
	INT color, nb_colors;
	VECTOR_OB new_col;
	MATRIX_OB new_Col;
	
	if (f_v) {
		printf("refine_by_p_roots() p = %ld\n", p);
		}
	p_roots(p, &pr);
	pr.multiplicities(&val, &mult);
	if (f_v) {
		str[0] = 0;
		val.sprint_multiplicities(&mult, str);
		printf("number of %ld-roots: %s\n", p, str);
		fflush(stdout);
		}
	
	nb_colors = Col->s_hi();
	Pr.m_il(nb_colors);
	for (color = 0; color < nb_colors; color++) {
		p_roots_in_color(p, color, col, &pr);
		pr.multiplicities(&val, &mult);
		if (f_v) {
			str[0] = 0;
			val.sprint_multiplicities(&mult, str);
			printf("number of %ld-roots of elements of color %ld: %s\n", 
				p, color, str);
			fflush(stdout);
			}
		pr.swap((VECTOR_OP) Pr.s_i(color));
		}
	color_v_join(col, &Pr, &new_col, &new_Col, f_v /* f_v */, FALSE /* f_vv */);
	if (f_v) {
		VECTOR_OB val, mult;
		BYTE str[10000];
		INT c, nc, a, b, nt, j;
		
		new_col.multiplicities(&val, &mult);
		if (f_v) {
			str[0] = 0;
			val.sprint_multiplicities(&mult, str);
			printf("next class coloring: %s\n", str);
			}
		nc = new_Col.s_hi();
		nt = new_Col.s_li() - 1;
		for (c = 0; c < nc; c++) {
			a = new_Col.s_iji(c, 0);
			if (f_v) {
				printf("refined color %ld means old color %ld and new type ", c, a);
				for (j = 0; j < nt; j++) {
					b = new_Col.s_iji(c, 1 + j);
					printf("%ld ", b);
					}
				printf("\n");
				}
			}
		}
	new_col.swap(col);
	new_Col.swap(Col);
	return OK;
}

INT FG_OB::p_roots(INT p, VECTOR_OP nb_p_roots)
{
	INT i, j, n;
	
	n = s_n_i();
	nb_p_roots->m_il_n(n);
	for (i = 0; i < n; i++) {
		j = gt_power(i, p);
		nb_p_roots->s_i(j)->inc();
		}
	return OK;
}

INT FG_OB::p_roots_in_color(INT p, INT color, 
	VECTOR_OP col, VECTOR_OP nb_p_roots)
/* determines (for each element of G) 
 * the number of p-th roots of elements of given color. */
{
	INT i, j, n;
	
	n = s_n_i();
	nb_p_roots->m_il_n(n);
	for (i = 0; i < n; i++) {
		if (col->s_ii(i) != color)
			continue;
		j = gt_power(i, p);
		nb_p_roots->s_i(j)->inc();
		}
	return OK;
}
#endif

#if 0
INT FG_OB::prepare_class_info(VECTOR_OP V, INT f_v, INT f_vv)
{
	VECTOR_OB class_no;
	VECTOR_OB sv_last, sv_gen;
	VECTOR_OB class_rep, class_len, class_el_ord;
	INT n;
	VECTOR_OB val, mult;
	VECTOR_OP key;
	STRING_OP s;
	INT i, l, m, cl, eo;
	BYTE str[256];

	n = s_n_i();
	if (f_v) {
		printf("classes of group %s:\n", s_label_s());
		fflush(stdout);
		}
	classes(&class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		print_classes(&class_no, &class_rep, &class_len);
		}
	get_ci_ordered(&class_len, &class_el_ord, &val, &mult, FALSE, FALSE);
	
	l = val.s_li();
	V->m_il(l);
	for (i = 0; i < l; i++) {
		key = (VECTOR_OP) val.s_i(i);
		m = mult.s_ii(i);
		cl = key->s_ii(0);
		eo = key->s_ii(1);
		sprintf(str, "%ld x %ld/%ld", m, cl, eo);
		if (i < l - 1)
			strcat(str, ", ");
		s = (STRING_OP) V->s_i(i);
		s->init(str);
		}
	return OK;
}
#endif

#if TEXDOCU
INT FG_OB::prepare_class_info(VECTOR_OP V, INT f_v, INT f_vv)
#endif
{
	GROUP_TABLE table;

	if (f_v) {
		printf("classes of group %s:\n", s_label_s());
		fflush(stdout);
		}
	group_init_from_fg(&table, this);
	group_prepare_class_info(&table, V, f_v, f_vv);
	group_free(&table);
	return OK;
}

#if TEXDOCU
INT group_init_from_fg(GROUP_TABLE *G, FG_OP fg)
#endif
{
	ZE_OP ze;
	INT i, j, k, n, nb_gen;

	if (fg->s_theG()->s_obj_k() != MATRIX)
		return error("group_init_from_fg(): fg->s_theG()->s_obj_k() != MATRIX");
	G->table = (MATRIX_OP) callocobject("fg_color.C: group_init_from_fg");
	G->table_inv = (VECTOR_OP) callocobject("fg_color.C: group_init_from_fg");
	G->generators = (VECTOR_OP) callocobject("fg_color.C: group_init_from_fg");
	G->n = n = fg->s_n_i();
	// printf("n=%ld ", n); fflush(stdout);
	nb_gen = fg->s_nb_ze_i();
	// printf("nb_gen=%ld ", nb_gen); fflush(stdout);
	G->table->m_ilih(n, n);
	G->table_inv->m_il(n);
	G->generators->m_il(nb_gen);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			k = fg->s_theG_iji(i, j);
			G->table->m_iji(i, j, k);
			}
		}
	for (i = 0; i < n; i++) {
		k = fg->s_theG_inv_i(i);
		G->table_inv->m_ii(i, k);
		}
	for (i = 0; i < nb_gen; i++) {
		ze = fg->s_ze_i(i);
		j = ze->s_n0_i();
		G->generators->m_ii(i, j);
		}
	return OK;
}


#endif /* SOLVABLE_TRUE */


