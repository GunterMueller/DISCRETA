/* gt_color.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/divs.h> /* for STRING_OP */
#include <DISCRETA/solvable.h>

#if TEXDOCU
INT group_init_simple(GROUP_TABLE *G, 
	MATRIX_OP M, VECTOR_OP generators)
#endif
{
	INT i, j, k, n, nb_gen;

	G->table = (MATRIX_OP) callocobject("fg_color.C group_init_simple()");
	G->table_inv = (VECTOR_OP) callocobject("fg_color.C group_init_simple()");
	G->generators = (VECTOR_OP) callocobject("fg_color.C group_init_simple()");
	G->n = n = M->s_hi();
	G->table->m_ilih(n, n);
	G->table_inv->m_il(n);
	generators->copy(G->generators);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			k = M->s_iji(i, j);
			G->table->m_iji(i, j, k);
			if (k == 0) { /* inv(i) = j */
				G->table_inv->m_ii(i, j);
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT group_init_generators(GROUP_TABLE *G, VECTOR_OP gen)
#endif
{
	gen->copy(G->generators);
	return OK;
}

#if TEXDOCU
INT group_free(GROUP_TABLE *G)
#endif
{
	freeall(G->table);
	freeall(G->table_inv);
	freeall(G->generators);
	G->table = NULL;
	G->table_inv = NULL;
	G->generators = NULL;
	G->n = 0;
	return OK;
}

#if TEXDOCU
INT group_mult(GROUP_TABLE *G, INT a, INT b)
#endif
{
	INT c;

	c = G->table->s_iji(a, b);
	return c;
}

#if TEXDOCU
INT group_inv(GROUP_TABLE *G, INT a)
#endif
{
	INT a_inv;

	a_inv = G->table_inv->s_ii(a);
	return a_inv;
}

#if TEXDOCU
INT group_conj(GROUP_TABLE *G, INT a, INT b)
#else
/* $b^{-1} a b$ */
#endif
{
	INT b_inv, c;
	
	b_inv = G->table_inv->s_ii(b);
	c = group_mult(G, b_inv, a);
	c = group_mult(G, c, b);
	return c;
}

#if TEXDOCU
INT group_commutator(GROUP_TABLE *G, INT a, INT b)
#else
/* $a^{-1} b^{-1} a b$ */
#endif
{
	INT a_inv, b_inv, c;
	
	a_inv = G->table_inv->s_ii(a);
	b_inv = G->table_inv->s_ii(b);
	c = group_mult(G, a_inv, b_inv);
	c = group_mult(G, c, a);
	c = group_mult(G, c, b);
	return c;
}

#if TEXDOCU
INT group_power(GROUP_TABLE *G, INT a, INT exp)
#endif
{
	INT c = 0, b = a;
	
	while (exp > 0) {
		/* res = b^exp * c */
		if (ODD(exp)) {
			c = group_mult(G, b, c);
			exp--;
			continue; /* exp == 0 possible */
			}
		if (EVEN(exp)) {
			b = group_mult(G, b, b);
			exp >>= 1;
			}
		}
	return c;
}

#if TEXDOCU
INT group_order(GROUP_TABLE *G, INT i, INT *ord)
#endif
{
	INT j, ord1;

	
	j = i;
	ord1 = 1;
	while (j != 0) {
		j = group_mult(G, i, j);
		ord1++;
		}
	*ord = ord1;
	return(OK);
}

#if TEXDOCU
INT group_order_if_prime(GROUP_TABLE *G, INT i, INT *ord)
#endif
{
	INT ord1;
	
	group_order(G, i, &ord1);
	if (ord1 == 1) {
		*ord = 0;
		return(OK);
		}
	if (smallest_primedivisor(ord1) == ord1)
		*ord = ord1;
	else
		*ord = 0;
	return(OK);
}

#if TEXDOCU
INT group_order_if_prime_power(GROUP_TABLE *G, INT i, 
	INT *ord, INT *prime, INT *k)
#endif
{
	INT ord1, prime1, prime2, k1;

	group_order(G, i, &ord1);
	if (ord1 == 1) {
		*ord = 0;
		return(OK);
		}
	*ord = ord1;
	prime1 = 0;
	k1 = 0;
	while (ord1 != 1) {
		prime2 = smallest_primedivisor(ord1);
		if (prime1) {
			if (prime2 != prime1) {
				/* different primes */
				*ord = 0;
				return(OK);
				}
			}
		prime1 = prime2;
		k1++;
		ord1 /= prime1;
		}
	*prime = prime1;
	*k = k1;
	return(OK);
}

#if TEXDOCU
INT group_order_structure(GROUP_TABLE *G, VECTOR_OP ord)
#endif
{
	INT i, n, o;
	
	n = G->n;
	ord->m_il(n);
	for (i = 0; i < n; i++) {
		group_order(G, i, &o);
		ord->m_ii(i, o);
		}
	return OK;
}

#if TEXDOCU
INT group_is_central(GROUP_TABLE *G, INT i)
#endif
{
	INT j, g_j, r;
	
	r = G->generators->s_li();
	for (j = 0; j < r; j++) {
		g_j = G->generators->s_ii(j);
		if (group_conj(G, i, g_j) != i)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT group_center(GROUP_TABLE *G, VECTOR_OP p)
#endif
{
	INT i, n, j;
	
	n = G->n;
	p->m_il(n);
	for (i = 0; i < n; i++) {
		j = group_is_central(G, i);
		p->m_ii(i, j);
		}
	return OK;
}

#if TEXDOCU
INT group_classes(GROUP_TABLE *G, VECTOR_OP class_no, 
	VECTOR_OP sv_last, VECTOR_OP sv_gen, 
	VECTOR_OP class_rep, VECTOR_OP class_len, VECTOR_OP class_el_ord)
#endif
{
	INT nb_generators, gen_i, i, n, j, nb_classes = 0, first, os, i0, i1, ord;
	INT *Q = NIL, nb_Q;
	
	n = G->n;
	nb_generators = G->generators->s_li();
	Q = (INT *) my_malloc(n * sizeof(INT), "gt_color.C: group_classes");
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
		group_order(G, first, &ord);
		while (nb_Q) {
			i0 = Q[0];
			for (i = 1; i < nb_Q; i++)
				Q[i - 1] = Q[i];
			nb_Q--;
			/* apply all generators to i0: */
			for (i = 0; i < nb_generators; i++) {
				gen_i = G->generators->s_ii(i);
				i1 = group_conj(G, i0, gen_i);
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

#if TEXDOCU
INT group_print_classes(GROUP_TABLE *G, 
	VECTOR_OP class_no, VECTOR_OP class_rep, VECTOR_OP class_len)
#endif
{
	INT n, i, j, o, first, nb_classes;
	
	n = G->n;
	nb_classes = class_rep->s_li();
	for (i = 0; i < nb_classes; i++) {
		first = class_rep->s_ii(i);
		printf("class %ld: { %ld ", i, first);
		for (j = first + 1; j < n; j++) {
			if (class_no->s_ii(j) == i) {
				printf("%ld ", j);
				}
			}
		group_order(G, first, &o);
		printf("} len %ld (el.-order %ld)\n", class_len->s_ii(i), o);
		}
	return OK;
}

#if TEXDOCU
INT group_first_class_coloring(GROUP_TABLE *G, 
	VECTOR_OP class_no, VECTOR_OP class_rep, 
	VECTOR_OP class_len, VECTOR_OP class_el_ord, 
	VECTOR_OP col, MATRIX_OP Col, INT f_v)
#endif
{
	VECTOR_OB cl, ceo;
	INT n, i, idx, a, b;
	
	n = G->n;
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

#if TEXDOCU
INT group_Classes(GROUP_TABLE *G, BYTE *str, INT f_v, INT f_vv)
#endif
{
	VECTOR_OB class_no;
	VECTOR_OB sv_last, sv_gen;
	VECTOR_OB class_rep, class_len, class_el_ord;
	INT n;

	n = G->n;
	group_classes(G, &class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		group_print_classes(G, &class_no, &class_rep, &class_len);
		}
	sprint_class_structure(&class_len, &class_el_ord, str, f_v);
	return OK;
}

#if TEXDOCU
INT group_Coloring_seldoms_prefered(GROUP_TABLE *G, VECTOR_OP col, INT f_v)
#endif
{
	VECTOR_OB col0;
	PERMUTATION_OB p1, p2;
	
	group_Coloring(G, &col0, f_v, f_v);
	prefer_seldom_colors(&col0, col, &p1, &p2);
	return OK;
}

#if TEXDOCU
INT group_Coloring(GROUP_TABLE *G, VECTOR_OP col, INT f_v, INT f_vv)
#endif
{
	VECTOR_OB class_no;
	VECTOR_OB sv_last, sv_gen;
	VECTOR_OB class_rep, class_len, class_el_ord;
	MATRIX_OB Col;
	BYTE str[10000];
	INT n;

	n = G->n;
	group_classes(G, &class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		group_print_classes(G, &class_no, &class_rep, &class_len);
		}
	str[0] = 0;
	sprint_class_structure(&class_len, &class_el_ord, str, f_v);
	
	/* first class coloring 
	 * (according to class length / element order): */
	group_first_class_coloring(G, &class_no, 
		&class_rep, &class_len, &class_el_ord, 
		col, &Col, f_v /* f_v */);
	
	group_refine_by_roots(G, col, &Col, f_v /* f_v */);
	
	return OK;
}

#if TEXDOCU
INT group_refine_by_roots(GROUP_TABLE *G, 
	VECTOR_OP col, MATRIX_OP Col, INT f_v)
#endif
{
	VECTOR_OB vp, ve;
	BYTE str[10000];
	INT n, i, l, p, nc0, nc;
	
	if (f_v)
		printf("group_refine_by_roots():\n");
	n = G->n;
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
			group_refine_by_p_roots(G, p, col, Col, f_v /* f_v */);
			}
		nc = number_of_colors(Col, f_v);
		if (nc == nc0)
			break;
		}
	return OK;
}

#if TEXDOCU
INT group_refine_by_p_roots(GROUP_TABLE *G, 
	INT p, VECTOR_OP col, MATRIX_OP Col, INT f_v)
#endif
{
	VECTOR_OB Pr, pr, val, mult;
	BYTE str[10000];
	INT color, nb_colors;
	VECTOR_OB new_col;
	MATRIX_OB new_Col;
	
	if (f_v) {
		printf("group_refine_by_p_roots() p = %ld\n", p);
		}
	group_p_roots(G, p, &pr);
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
		group_p_roots_in_color(G, p, color, col, &pr);
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

#if TEXDOCU
INT group_p_roots(GROUP_TABLE *G, INT p, VECTOR_OP nb_p_roots)
#endif
{
	INT i, j, n;
	
	n = G->n;
	nb_p_roots->m_il_n(n);
	for (i = 0; i < n; i++) {
		j = group_power(G, i, p);
		nb_p_roots->s_i(j)->inc();
		}
	return OK;
}

#if TEXDOCU
INT group_p_roots_in_color(GROUP_TABLE *G, INT p, INT color, 
	VECTOR_OP col, VECTOR_OP nb_p_roots)
#else
determines (for each element of G) 
the number of p-th roots of elements of given color. 
#endif
{
	INT i, j, n;
	
	n = G->n;
	nb_p_roots->m_il_n(n);
	for (i = 0; i < n; i++) {
		if (col->s_ii(i) != color)
			continue;
		j = group_power(G, i, p);
		nb_p_roots->s_i(j)->inc();
		}
	return OK;
}

#if TEXDOCU
INT group_prepare_class_info(GROUP_TABLE *G, 
	VECTOR_OP V, INT f_v, INT f_vv)
#endif
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

	n = G->n;
	group_classes(G, &class_no, &sv_last, &sv_gen, 
		&class_rep, &class_len, &class_el_ord);
	if (f_vv) {
		group_print_classes(G, &class_no, &class_rep, &class_len);
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


#endif /* SOLVABLE_TRUE */


