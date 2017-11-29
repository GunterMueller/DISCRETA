/* lb.C: */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LABRA_TRUE

#include <DISCRETA/lb.h>
#include <DISCRETA/ma.h> // for calc_transversal_matrix()

#undef DEBUG_LABRA
#undef DEBUG_SIFT

#undef ANGST

static void traverse_tree(FILE *fp, LABRA_OP L, INT *path, INT depth, INT f_full, INT f_labels);
static INT traverse_tree_cyclic(LABRA_OP L, INT *path, INT depth);

/* 
 * LABRA
 */
#if TEXDOCU
#else
LABRA stands for \lq labelled branching\rq. 
A labelled branching is a data structure for permutation groups. 
It has been invented by Mark Jerrum in~\cite{Jerrum86}.
An algorithm for changing the base of a permutation group 
is described by Brown, Finkelstein and Purdom~\cite{Brownetal89a}.
#endif

#if TEXDOCU
INT LABRA_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "degree"; break;
	case 1: s = "nb_gen"; break;
	case 2: s = "G"; break;
	case 3: s = "V"; break;
	case 4: s = "KM"; break;
	case 5: s = "EM"; break;
	case 6: s = "base"; break;
	case 7: s = "inv_base"; break;
	case 8: s = "T"; break;
	case 9: s = "tidx"; break;
	case 10: s = "f_verbose"; break;
	case 11: s = "f_very_verbose"; break;
	default:
		return error("LABRA::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::init_quick(VECTOR_OP gen)
#endif
{
	PERMUTATION_OP perm;
	INT l, deg;
	
	l = gen->s_li();
	if (l == 0) {
		return error("LABRA::init_quick() no generators !");
		}
	perm = (PERMUTATION_OP) gen->s_i(0);
	deg = perm->s_li();
	Init(deg, 
		gen, l, 
		FALSE /* f_verbose */, 
		FALSE /* f_very_verbose */);
	jerrum(FALSE /* f_verbose */);
	// group_order(&go);
#if 0
	if (f_verbose) {
		Print();
		Print_T(TRUE /* f_verbose */);
		}
#endif
	return OK;
}

#if TEXDOCU
INT LABRA_OB::Init_no_generators(INT degree, INT f_verbose, INT f_very_verbose)
#endif
{
	VECTOR_OB G;

	G.m_il(0);
	Init(degree, &G, 0, f_verbose, f_very_verbose);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::Init(INT degree, VECTOR_OP theGenerators, INT nb_gen, 
	INT f_verbose, INT f_very_verbose)
#endif
{
	PERMUTATION_OP p;
	INT i;
	
#if 0
	if (s_obj_k() != EMPTY) {
		printf("LABRA_OB::Init() warning: non-EMPTY object, freeing\n");
		fflush(stdout);
		freeself();
		}
#endif
	m_il(12);
	c_obj_k(LABRA_KIND);

	if (f_very_verbose) {
		printf("opening LABRA of degree %ld, nb_gen = %ld\n", 
			degree, nb_gen);
		fflush(stdout);
		}
	s_degree()->m_i(degree);
	s_nb_gen()->m_i(nb_gen);
	theGenerators->copy(s_G());

	s_V()->m_il_n(degree);
	s_KM()->m_il(degree);
	s_EM()->m_il(degree);
	for (i = 0; i < degree; i++) {
		p = s_KM_i(i);
		p->m_il(degree);
		p->one();
		p = s_EM_i(i);
		p->m_il(degree);
		p->one();
		s_V_i(i)->m_i(i);
		}
	p = s_base();
	p->m_il(degree);
	p->one();
	p = s_inv_base();
	p->m_il(degree);
	p->one();
	s_T()->m_il(degree);
	for (i = 0; i < degree; i++) {
		s_T_i(i)->m_il(1);
		s_T_ij(i, 0)->m_i(i);
		}
	s_tidx()->m_il(0);
	s_f_verbose()->m_i(f_verbose);
	s_f_very_verbose()->m_i(f_very_verbose);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::build(VECTOR_OP gen)
#endif
{
	SYM_OB go;
	PERMUTATION_OP p;
	VECTOR_OB gen1;
	INT deg;
	
	gen->copy(&gen1);
	p = (PERMUTATION_OP) gen->s_i(0);
	deg = p->s_li();
	Init(deg, &gen1, gen1.s_li(), 
		FALSE /* f_verbose */, 
		FALSE /* f_very_verbose */);
	jerrum(FALSE /* f_verbose */);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::sprint(BYTE *s)
#endif
{
	BYTE str1[256];
	
	sprintf(str1, "LABRA: degree = %ld", s_degree_i());
	if (strlen(s) + strlen(str1) < 200)
		strcat(s, str1);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::check_consistency()
#endif
{
	INT i, j, i_el, j_el, n;
	PERMUTATION_OB per;

	n = s_degree_i();
	for (i = 0; i < n; i++) {
		i_el = s_base_ii(i) - 1;
		for (j = 0; j < n; j++) {
			j_el = s_base_ii(j) - 1;
			if (path_exists(i, j)) {
				rep_ij(i, j, &per);
				if (per.s_ii(i_el) - 1 != j_el) {
					printf("i = %ld i_el = %ld j = %ld j_el = %ld\n", i, i_el, j, j_el);
					return error("LABRA_OB::check_consistency() per.s_ii(i_el) - 1 != j_el");
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::Print_Sims()
#endif
{
	printf("LABRA: degree = %ld nb_gen = %ld\n", s_degree_i(), s_nb_gen_i());
	print_T();
#if 0
	INT i, j, i_el, j_el, len, t;
	SYM_OB stab_ord;
	SYM_OB group_order;
	INTEGER_OB len1;
	
	((INTEGER_OP) &group_order)->m_i(1);
	printf("               ");
	for (j = 0; j < s_degree_i(); j++) {
		if (j % 10 == 0)
			printf("|");
		else
			printf(" ");
		}
	printf(" length, group index, group order\n");
	for (i = 0; i < s_degree_i(); i++) {
		i_el = s_base()->s_ii(i) - 1;
		printf("T %3ld ", i);
		printf("base %3ld ", i_el);
		len = 0;
		for (j = 0; j < s_degree_i(); j++) {
			j_el = s_base()->s_ii(j) - 1;
			if (path_exists(i, j)) {
				len++;
				printf("*");
				}
			else
				printf(".");
			}
		len1.m_i(len);
		len1.mult_apply(&group_order);
		// stab_order(i, &stab_ord);
		printf(" %2ld ", len);
		// group_order.print();
		// printf(", ");
		// stab_ord.print();
		printf("\n");
		}
#endif
	return OK;
}

#if TEXDOCU
INT LABRA_OB::Print()
#endif
{
	INT i;
	
	printf("LABRA: degree = %ld nb_gen = %ld\n", s_degree_i(), s_nb_gen_i());
	printf("base: \n");
	for (i = 0; i < s_degree_i(); i++)
		printf("%ld %ld\n", i, s_base()->s_ii(i) - 1);
	printf("at: V (prev)    KM ;     EM:\n");
	for (i = 0; i < s_degree_i(); i++) {
		printf("%ld: %ld ", i + 1, s_V_ii(i) + 1);
#if 1
		s_KM_i(i)->print();
		printf("; ");
		s_EM_i(i)->print();
#else
		s_KM_i(i)->print();
#endif
		printf("\n");
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::print_T()
#endif
{
	INT i, t, len;
	
	len = s_tidx()->s_li();
	for (t = 0; t < len; t++) {
		i = s_tidx_ii(t);
		printf("transversal %ld: ", i);
		s_T_i(i)->println();
		}

#if 0
	INT i, j, i_el, j_el, len;
	SYM_OB group_order;
	INTEGER_OB len1;
	PERMUTATION_OB m, m1, t;
	
	((INTEGER_OP) &group_order)->m_i(1);
	for (i = 0; i < s_degree_i() - 1; i++) {
		i_el = s_base()->s_ii(i) - 1;
		s_EM_i(i_el)->copy(&m1);
		m1.invers(&m);
		printf("T %ld:\n", i);
		printf("base = %ld\n", i_el);
		len = 0;
		for (j = 0; j < s_degree_i(); j++) {
			j_el = s_base()->s_ii(j) - 1;
			if (path_exists(i, j)) {
				len++;
				m.mult(s_EM_i(j_el), &t);
				printf("%ld ", j);
				if (f_verbose)
					t.println();
				}
			}
		len1.m_i(len);
		len1.mult_apply(&group_order);
		printf("of length %ld; (sub) group index: ", len);
		group_order.println();
		}
#endif
	return OK;
}

#if TEXDOCU
INT LABRA_OB::latex_transversal_indices(FILE *fp)
#endif
{
	INT i, i_el, j, j_el, t, tt, len, len_t;
	PERMUTATION_OB m, m1, mm;
	
	len = s_tidx()->s_li();
	for (t = 0; t < len; t++) {
		i = s_tidx_ii(t);
		len_t = s_T_i(i)->s_li();
		
		i_el = s_base()->s_ii(i) - 1;
		s_EM_i(i_el)->copy(&m1);
		m1.invers(&m);
		fprintf(fp, "orbit of %ld is $\\{ ", i_el + 1);
		for (tt = 0; tt < len_t; tt++) {
			j = s_T_iji(i, tt);
			j_el = s_base()->s_ii(j) - 1;
			fprintf(fp, "%ld \\,", j_el + 1);
			}
		fprintf(fp, "\\}$\\\\\n");
		}
	
	return OK;


#if 0
	INT i, j, i_el, j_el, len, l;
	SYM_OB group_order;
	VECTOR_OB ol;
	INTEGER_OB len1;
	PERMUTATION_OB m, m1, t;
	
	((INTEGER_OP) &group_order)->m_i(1);
	ol.m_il(0);
	for (i = 0; i < s_degree_i() - 1; i++) {
		i_el = s_base()->s_ii(i) - 1;
		s_EM_i(i_el)->copy(&m1);
		m1.invers(&m);
		fprintf(fp, "orbit of %ld is $\\{ ", i_el + 1);
		for (j = 0; j < s_degree_i(); j++) {
			j_el = s_base()->s_ii(j) - 1;
			if (path_exists(i, j)) {
				fprintf(fp, "%ld \\,", j_el + 1);
				}
			}
		fprintf(fp, "\\}$\\\\\n");
		len = 0;
		for (j = 0; j < s_degree_i(); j++) {
			j_el = s_base()->s_ii(j) - 1;
			if (path_exists(i, j)) {
				len++;
				m.mult(s_EM_i(j_el), &t);
				fprintf(fp, "$");
				t.latex(stdout);
				fprintf(fp, "$\\\\\n");
				}
			}
		if (len > 1) {
			ol.inc();
			ol.m_ii(ol.s_li() - 1, len);
			}
		len1.m_i(len);
		len1.mult_apply(&group_order);
		// printf("of length %ld; (sub) group index: ", len);
		// group_order.println();
		}
	fprintf(fp, "the group order is $");
	l = ol.s_li();
	for (i = 0; i < l; i++) {
		fprintf(fp, "%ld ", ol.s_ii(i));
		if (i < l - 1)
			fprintf(fp, "\\cdot ");
		}
	fprintf(fp, "= ");
	group_order.fprint(fp);
	fprintf(fp, "$\\\\\n");

	return OK;
#endif
}

#if TEXDOCU
INT LABRA_OB::calc_T()
#endif
{
	INT i, v, l;
	
	s_T()->m_il(s_degree_i());
	for (i = 0; i < s_degree_i(); i++) {
		s_T_i(i)->m_il(0);
		}
	s_tidx()->m_il(0);
	
	for (i = 0; i < s_degree_i(); i++) {
		s_T_i(i)->append_i(i);
		v = i;
		while (s_V_ii(v) != v) {
			v = s_V_ii(v);
			s_T_i(v)->append_i(i);
			}
		}
	for (i = 0; i < s_degree_i(); i++) {
		l = s_T_i(i)->s_li();
		if (l > 1)
			s_tidx()->append_i(i);
		}	
	return OK;
}

#if TEXDOCU
INT LABRA_OB::calc_transversal_matrix(MATRIX_OP M, INT f_v)
#endif
{
	INT i, j, k, n, len_t;
	
	n = s_degree_i();
	M->m_ilih_n(n, n);
	for (i = 0; i < n; i++) {
		len_t = s_T_i(i)->s_li();
		for (k = 0; k < len_t; k++) {
			j = s_T_iji(i, k);
			M->m_iji(i, j, 1);
			}
		}
	return OK;
#if 0
	INT i, j, i_el, len, n;
	SYM_OB group_order, stab_ord;
	INTEGER_OB len1;
	
	n = s_degree_i();
	M->m_ilih_n(n, n);
	((INTEGER_OP) &group_order)->m_i(1);
	if (f_v) {
		for (j = 0; j < n; j++) {
			if (j % 10 == 0)
				printf("|");
			else
				printf(" ");
			}
		printf(" length, group index, group order\n");
		}
	for (i = 0; i < n - 1; i++) {
		i_el = s_base()->s_ii(i) - 1;
		if (f_v) {
			printf("T %3ld ", i);
			printf("base %3ld ", i_el);
			}
		len = 0;
		for (j = 0; j < n; j++) {
			if (path_exists(i, j)) {
				len++;
				if (f_v)
					printf("*");
				M->m_iji(i, j, 1);
				}
			else {
				if (f_v)
					printf(".");
				}
			}
		if (f_v) {
			len1.m_i(len);
			len1.mult_apply(&group_order);
			stab_order(i, &stab_ord);
			printf(" %2ld, ", len);
			group_order.print();
			printf(", ");
			stab_ord.print();
			printf("\n");
			}
		}
	M->m_iji(n - 1, n - 1, 1);
	return OK;
#endif
}

#if TEXDOCU
INT LABRA_OB::is_cyclic(VECTOR_OP path, PERMUTATION_OP per)
#endif
{
	INT i, a, l, ret;
	INT *p;
	
	l = s_degree_i();
	p = (INT *) my_malloc(sizeof(INT) * l, "p");
	
	ret = traverse_tree_cyclic(this, p, 0);
	if (ret) {
		path->m_il(l);
		for (i = 0; i < l; i++) {
			a = p[i];
			path->m_ii(i, a);
			}
		element_from_path(per, p, l);
		}
	
	my_free(p);
	return ret;
}

static INT traverse_tree_cyclic(LABRA_OP L, INT *path, INT depth)
{
	INT j, i_el, j_el;
	PERMUTATION_OB p;
	
	if (depth == L->s_degree_i()) {
		L->element_from_path(&p, path, depth);
		// printf("is_cyclic: ");
		// p.println();
		if (p.is_full_cycle())
			return TRUE;
		else
			return FALSE;
		}
	i_el = L->s_base()->s_ii(depth) - 1;
	for (j = 0; j < L->s_degree_i(); j++) {
		j_el = L->s_base()->s_ii(j) - 1;
		if (L->path_exists(depth, j)) {
			path[depth] = j;
			if (traverse_tree_cyclic(L, path, depth + 1))
				return TRUE;
			} // if path_exists
		} // next j
	return FALSE;
}


#if 0
#if TEXDOCU
INT LABRA_OB::elements_first(MATRIX_OP transversals, VECTOR_OP path, PERMUTATION_OP p)
#endif
{
	INT i, n;
	
	n = s_degree_i();
	calc_transversal_matrix(transversals, FALSE /* f_v */);
	path->m_il(n);
	for (i = 0; i < n; i++) {
		path->m_ii(i, i);
		}
	element_from_path_vec(p, path);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::elements_next(MATRIX_OP transversals, VECTOR_OP path, PERMUTATION_OP p)
#endif
{
	INT i, j, n;
	
	n = s_degree_i();
	for (i = n - 1; i >= 0; i--) {
		j = path->s_ii(i);
		for (j++; j < n; j++) {
			if (transversals->s_iji(i, j)) {
				path->m_ii(i, j);
				element_from_path_vec(p, path);
				return TRUE;
				}
			}
		path->m_ii(i, i);
		}
	return FALSE;
}
#endif

#if TEXDOCU
INT LABRA_OB::transversal_tree(BYTE *group_label, INT f_full, INT f_labels)
#endif
{
	BYTE fname[1024];
	FILE *fp;
	INT *path;

	if (f_full)
		sprintf(fname, "%s_tree_full.txt", group_label);
	else
		sprintf(fname, "%s_tree.txt", group_label);
	fp = fopen(fname, "w");

	path = (INT *) my_malloc(sizeof(INT) * s_degree_i(), "path");
	
	fprintf(fp, "COORDS_RANGE 1000000 700000\n");
	fprintf(fp, "TREE\n");

	traverse_tree(fp, this, path, 0, f_full, f_labels);
	
	my_free(path);
	fprintf(fp, "TREE_END\n");
	
	fclose(fp);
	return OK;
}

static void traverse_tree(FILE *fp, LABRA_OP L, INT *path, INT depth, INT f_full, INT f_labels)
{
	INT j, i_el, j_el, k;
	PERMUTATION_OB p;
	BYTE str[1000];
	
	if (depth == L->s_degree_i()) {
		fprintf(fp, "%ld ", depth);
		for (k = 0; k < depth; k++) {
			fprintf(fp, "%ld ", path[k] + 1);
			}
		if (f_labels) {
			L->element_from_path(&p, path, depth);
			str[0] = 0;
			p.sprint_latex(str);
			fprintf(fp, " $%s$\n", str);
			}
		fprintf(fp, "\n");
		return;
		}
	i_el = L->s_base()->s_ii(depth) - 1;
	for (j = 0; j < L->s_degree_i(); j++) {
		j_el = L->s_base()->s_ii(j) - 1;
		if (L->path_exists(depth, j)) {
			path[depth] = j;
			if (f_full) {
				traverse_tree(fp, L, path, depth + 1, f_full, f_labels);
				}
			else {
				if (depth == L->s_degree_i() - 1 || j == depth) {
					traverse_tree(fp, L, path, depth + 1, f_full, f_labels);
					}
				else {
					fprintf(fp, "%ld ", depth + 1);
					for (k = 0; k <= depth; k++) {
						fprintf(fp, "%ld ", path[k] + 1);
						}
					if (f_labels) {
						L->element_from_path(&p, path, depth + 1);
						str[0] = 0;
						p.sprint_latex(str);
						fprintf(fp, " $%s$", str);
						}
					fprintf(fp, "\n");
					}
				} // else
			} // if path_exists
		} // next j
}

#if TEXDOCU
INT LABRA_OB::elements_first(VECTOR_OP path, PERMUTATION_OP p)
#endif
{
	INT i, l;
	
	l = s_tidx()->s_li();
	path->m_il(l);
	for (i = 0; i < l; i++) {
		path->m_ii(i, 0);
		}
	element_from_short_path(p, path);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::elements_next(VECTOR_OP path, PERMUTATION_OP p)
#endif
{
	INT i, j, l, t, t_len;
	
	l = s_tidx()->s_li();
	for (i = l - 1; i >= 0; i--) {
		t = s_tidx_ii(i);
		t_len = s_T_i(t)->s_li();
		j = path->s_ii(i) + 1;
		if (j < t_len) {
			path->m_ii(i, j);
			element_from_short_path(p, path);
			return TRUE;
			}
		path->m_ii(i, 0);
		}
	return FALSE;
}
#if TEXDOCU
INT LABRA_OB::element_from_path(PERMUTATION_OP p, INT *path, INT path_len)
#endif
{
	INT i, j;
	PERMUTATION_OB a, b, c;
	
	b.m_il(s_degree_i());
	b.one();
	for (i = 0; i < path_len; i++) {
		j = path[i];
		rep_ij(i, j, &a);
		a.mult(&b, &c); // pre-multiply !
		c.swap(&b);
		}
	b.swap(p);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::element_from_short_path(PERMUTATION_OP p, VECTOR_OP path)
#endif
{
	INT i, j, jj, t, n, len;
	PERMUTATION_OB a, b, c;
	
	n = s_degree_i();
	len = s_tidx()->s_li();
	if (len != path->s_li())
		return error("element_from_short_path(): len != path->s_li()");

	b.m_il(n);
	b.one();
	for (i = 0; i < len; i++) {
		j = path->s_ii(i);
		t = s_tidx_ii(i);
		jj = s_T_iji(t, j);
		rep_ij(t, jj, &a);
		a.mult(&b, &c); // pre-multiply !
		c.swap(&b);
		}
	b.swap(p);
	return OK;
}

#if 0
#if TEXDOCU
INT LABRA_OB::element_from_short_path(PERMUTATION_OP p, VECTOR_OP path)
#endif
{
	INT i, j, jj, t, n, len;
	PERMUTATION_OB a, b, c;
	VECTOR_OB path2;
	
	n = s_degree_i();
	path2.m_il(n);
	for (i = 0; i < n; i++)
		path2.m_ii(i, i);
	len = s_tidx()->s_li();
	if (len != path->s_li())
		return error("element_from_short_path(): len != path->s_li()");
	for (i = 0; i < len; i++) {
		j = path->s_ii(i);
		t = s_tidx_ii(i);
		jj = s_T_iji(t, j);
		path2.m_ii(t, jj);
		}
	b.m_il(n);
	b.one();
	for (i = 0; i < n; i++) {
		j = path2.s_ii(i);
		rep_ij(i, j, &a);
		a.mult(&b, &c); // pre-multiply !
		c.swap(&b);
		}
	b.swap(p);
	return OK;
}
#endif

#if TEXDOCU
INT LABRA_OB::transversal_as_set(INT i, VECTOR_OP chi)
#endif
// warning: base indices are returned, not element numbers ! 
// (as indices into chi)
{
	INT j, n;
	
	n = s_degree_i();
	chi->m_il(n);
	// i_el = s_base()->s_ii(i) - 1;
	for (j = 0; j < n; j++) {
		if (path_exists(i, j)) {
			chi->m_ii(j, 1);
			}
		else {
			chi->m_ii(j, 0);
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::stab_order(INT i, SYM_OP ord)
#endif
{
	INT ii, j, len;
	INTEGER_OB len1;
	
	ord->freeself();
	((INTEGER_OP) ord)->m_i(1);
	for (ii = i; ii < s_degree_i() - 1; ii++) {
		// i_el = s_base()->s_ii(ii) - 1;
		len = 0;
		for (j = 0; j < s_degree_i(); j++) {
			if (path_exists(ii, j))
				len++;
			}
		len1.m_i(len);
		len1.mult_apply(ord);
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::stab_generators(INT i, VECTOR_OP gen, INT f_v)
#endif
{
	PERMUTATION_OP g0, g;
	PERMUTATION_OB g0v, h;
	INT ii, i_el, ii_el, iii, l;
	INT idx, f_found;

	gen->m_il(0);
	l = s_degree_i();
	i_el = s_base()->s_ii(i) - 1;
	g0 = s_EM_i(i_el);
	g0->invers(&g0v);
	if (f_v) {
		printf("stab(%ld), base_el = %ld\n", i, i_el);
		}
	for (ii = 0; ii < l; ii++) {
		ii_el = s_base()->s_ii(ii) - 1;
		if (path_exists(i, ii)) {
			g = s_EM_i(ii_el);
			g0v.mult(g, &h);
			if (f_v) {
				printf("in %ld: ", ii);
				h.println();
				}
			if (!h.einsp()) {
				/* empty stabilizer generators 
				 * causes problems in VS_search */
				gen->do_search(gen->s_li(), TRUE, &h, &idx, &f_found, 0, NIL);
				if (!f_found) {
					gen->inc();
					for (iii = gen->s_li() - 1; iii > idx; iii--)
						gen->s_i(iii)->swap(gen->s_i(iii - 1));
					h.copy((PERMUTATION_OP) gen->s_i(idx));
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::generators(VECTOR_OP p)
#endif
{
	INT i;
	INTEGER_OB len1;
	
	if (!s_base()->einsp())
		return error("LABRA::generators(): error: this routine "
			"has a bug for non-trivial base !");

	p->m_il(0);
	for (i = 0; i < s_degree_i() /* - 1 */; i++) {
		if (s_EM_i(i)->s_obj_k() != EMPTY && !s_EM_i(i)->einsp()) {
			// i_el also possible, order does not matter
			p->inc();
			s_EM_i(i)->copy((PERMUTATION_OP) p->s_i(p->s_li() - 1));
			}
		}
	if (p->s_li() == 0) {
		for (i = 0; i < s_degree_i() - 1; i++) {
			if (s_EM_i(i)->s_obj_k() != EMPTY) {
				p->inc();
				s_EM_i(i)->copy((PERMUTATION_OP) p->s_i(p->s_li() - 1));
				break;
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::generators_bottom_up(VECTOR_OP p, INT nb_fix)
#endif
{
	INT i, j;
	PERMUTATION_OB rep;
	
	p->m_il(0);
	for (i = nb_fix; i < s_degree_i(); i++) {
		for (j = 0; j < s_degree_i(); j++) {
			if (path_exists(i, j)) {
				rep_ij(i, j, &rep);
				p->inc();
				rep.copy((PERMUTATION_OP) p->s_i(p->s_li() - 1));
				}
			}
		}
	if (p->s_li() == 0) {
		for (i = 0; i < s_degree_i(); i++) {
			if (s_EM_i(i)->s_obj_k() != EMPTY) {
				p->inc();
				s_EM_i(i)->copy((PERMUTATION_OP) p->s_i(p->s_li() - 1));
				break;
				}
			}
		}
	return OK;
}

#define GENERATORS_FEAR

#if TEXDOCU
INT LABRA_OB::reduced_generating_set(VECTOR_OP p, INT f_bottom_up, INT f_v)
#endif
{
	VECTOR_OB gen0, gen2;
	SYM_OB go0, go1;
	LABRA_OB L;
	INT i, l, n, l1;

	n = s_degree_i();
	generators(&gen0);
	group_order(&go0);
	l = gen0.s_li();
	if (f_v) {
		printf("reduced_generating_set() group_order = ");
		go0.println();
		printf("number of generators (not reduced) = %ld\n", l);
		}
	if (l == 0) {
		gen0.swap(p);
		return OK;
		}
	l1 = 1;
	gen2.m_il(l1);
	if (f_bottom_up) {
		gen0.s_i(l - 1)->copy(gen2.s_i(0));
		}
	else {
		gen0.s_i(0)->copy(gen2.s_i(0));
		}
	
	if (f_bottom_up) {
		for (i = l - 2; i >= 0; i--) {
			L.init_by_generators(&gen2, n, FALSE);
			if (!L.in_test((PERMUTATION_OP) gen0.s_i(i))) {
				gen2.inc();
				gen0.s_i(i)->copy(gen2.s_i(l1));
				l1++;
				}
			}
		}
	else {
		i = 1;
		while (TRUE) {
			// L.init_by_generators(&gen2, n, FALSE);
			L.init_quick(&gen2);
			L.group_order(&go1);
			if (go1.sym_comp(&go0) == 0)
				break;
			for (; i < l; i++) {
				if (!L.in_test((PERMUTATION_OP) gen0.s_i(i))) {
					gen2.inc();
					gen0.s_i(i)->copy(gen2.s_i(l1));
					l1++;
					break;
					}
				}
			if (i == l)
				return error("LB::reduced_generating_set() i == l");
			}
		}
#ifdef GENERATORS_FEAR
	L.init_by_generators(&gen2, n, FALSE);
	L.group_order(&go1);
	if (go0.sym_comp(&go1) != 0) {
		go0.println();
		go1.println();
		return error("reduced_generating_set() go0 != go1");
		}
#endif
	if (f_v) {
#ifdef GENERATORS_FEAR
		printf("reduced_generating_set() group_order = ");
		go1.println();
#endif
		printf("number of reduced generators = %ld\n", l1);
		}
	gen2.swap(p);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::strong_generating_set(VECTOR_OP p, INT f_v)
#endif
{
	VECTOR_OB gen0, gen2;
	SYM_OB go0, go1;
	LABRA_OB L;
	INT i, l, n, l1;

	n = s_degree_i();
	generators_bottom_up(&gen0, 0);
	group_order(&go0);
	l = gen0.s_li();
	if (f_v) {
		printf("strong_generating_set() group_order = ");
		go0.println();
		printf("number of generators (not reduced) = %ld\n", l);
		}
	if (l == 0) {
		gen0.swap(p);
		return OK;
		}
	l1 = 1;
	gen2.m_il(l1);
	gen0.s_i(l - 1)->copy(gen2.s_i(0));
	
		for (i = l - 2; i >= 0; i--) {
			L.init_by_generators(&gen2, n, FALSE);
			if (!L.in_test((PERMUTATION_OP) gen0.s_i(i))) {
				gen2.inc();
				gen0.s_i(i)->copy(gen2.s_i(l1));
				l1++;
				}
			}
#ifdef GENERATORS_FEAR
	L.init_by_generators(&gen2, n, FALSE);
	L.group_order(&go1);
	if (go0.sym_comp(&go1) != 0) {
		go0.println();
		go1.println();
		return error("strong_generating_set() go0 != go1");
		}
#endif
	if (f_v) {
#ifdef GENERATORS_FEAR
		printf("strong_generating_set() group_order = ");
		go1.println();
#endif
		printf("number of strong generators = %ld\n", l1);
		}
	gen2.swap(p);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::group_order(SYM_OP go)
#endif
{
	INT i, j, len;
	INTEGER_OB len1;
	
	((INTEGER_OP) go)->m_i(1);
	for (i = 0; i < s_degree_i() - 1; i++) {
		len = 0;
		for (j = 0; j < s_degree_i(); j++) {
			if (path_exists(i, j)) {
				len++;
				}
			}
		len1.m_i(len);
		len1.mult_apply(go);
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::print_group_order()
#endif
{
	SYM_OB go;
	
	group_order(&go);
	printf("has group order ");
	go.println();
	return OK;
}

#if TEXDOCU
INT LABRA_OB::depth(INT i)
#endif
{
	INT j1, j2;
	
	j1 = i;
	while ((j2 = s_V_ii(j1)) != j1)
		j1 = j2;
	return j1;
}

#if TEXDOCU
INT LABRA_OB::rep_ij(INT i, INT j, PERMUTATION_OP rep)
#endif
{
	INT i_el, j_el;
	PERMUTATION_OB m1, m;
	
	i_el = s_base()->s_ii(i) - 1;
	j_el = s_base()->s_ii(j) - 1;
	s_EM_i(i_el)->copy(&m1);
	m1.invers(&m);
	if (path_exists(i, j)) {
		m.mult(s_EM_i(j_el), rep);
		}
	else {
		printf("i = %ld i_el = %ld j = %ld j_el = %ld\n", i, i_el, j, j_el);
		return error("LABRA::rep_ij() no path i -> j ");
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::path_exists(INT j, INT k)
#else
/* TRUE, iff path $j -> k$ exists. */
#endif
{
	INT j1, j2;
	
	if (j == k)
		return TRUE;
	j1 = k;
	while (((j2 = s_V_ii(j1)) != j1) && (j2 != j))
		j1 = j2;
	if (j2 == j)
		return TRUE;
	return FALSE;
}

#if TEXDOCU
INT LABRA_OB::path_exists_el(INT j_el, INT k_el)
#else
TRUE, iff path $j -> k$ exists 
$j := inv_base[j_el]$, 
$k := inv_base[k_el]$ 
#endif
{
	INT j, k;

	j = s_inv_base_ii(j_el) - 1;
	k = s_inv_base_ii(k_el) - 1;
	return path_exists(j, k);
}

#if TEXDOCU
INT LABRA_OB::child_vectors(VECTOR_OP cv)
#else
$cv[i][j] = 1$ iff $\pi_j$ is a child of $\pi_i$
#endif
{
	VECTOR_OB c;
	INT i, j, n;

	n = s_degree_i();
	cv->m_il(n);
	for (i = 0; i < n; i++) {
		c.m_il_n(n);
		for (j = 0; j < n; j++) {
			if (path_exists(i, j))
				c.m_ii(j, 1);
			}
		c.swap((VECTOR_OP) cv->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::new_edge(INT i_el, INT j_el, PERMUTATION_OP m)
#endif
{
	PERMUTATION_OB tmp;
	INT i, j;
	
#ifdef DEBUG_LABRA
	printf("new edge (%ld %ld) KM = ", i_el, j_el);
	m->println();
#endif
	i = s_inv_base_ii(i_el) - 1;
	j = s_inv_base_ii(j_el) - 1;
	s_V_i(j)->m_i(i);
	if (m->s_ii(i_el) - 1 != j_el) {
		printf("new edge (%ld %ld) KM = ", i_el, j_el);
		m->println();
		return error("new_edge() wrong m");
		}
	m->copy(s_KM_i(j_el));
	return OK;
}

#if TEXDOCU
INT LABRA_OB::update_EM()
#endif
{
	INT i, i_el, v, v_el;
	
	for (i = 0; i < s_degree_i(); i++) {
		i_el = s_base()->s_ii(i) - 1;
		v = s_V_ii(i);
		if (v == i)
			s_KM_i(i_el)->copy(s_EM_i(i_el));
		else {
			v_el = s_base()->s_ii(v) - 1;
			s_EM_i(v_el)->mult(s_KM_i(i_el), s_EM_i(i_el));
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::orbit(INT omega, VECTOR_OP G, INT G_len, INT f_verbose)
#else
computes the orbit of omega under $G$.
computes fusels into $KM[]$, 
i.e. $KM[j]$ is a permutation in $<G>$ which maps $omega$ to $j$.
the back pointers in $V[]$ are set (in new\_edge). 
The all point back to omega (for all orbit elements).
#endif
{
	PERMUTATION_OP g_i;
	PERMUTATION_OB tmp;
	INT *queue = NIL;
	INT *f_reached = NIL;
	INT nb_queue = 0, i, j, k;
	INT erg = OK;

	queue = (INT *) my_malloc(s_degree_i() * sizeof(INT), "lb.C: orbit()");
	if (queue == NIL)
		return error("LABRA::orbit() no memory for queue");
	f_reached = (INT *) my_malloc(s_degree_i() * sizeof(INT), "lb.C: orbit()");
	if (f_reached == NIL)
		return error("LABRA::orbit() no memory for f_reached");
	for (i = 0; i < s_degree_i(); i++)
		f_reached[i] = FALSE;
	f_reached[omega] = TRUE;
	queue[0] = omega;
	nb_queue = 1;
	if (f_verbose) {
		printf("orbit: { %ld, ", omega);
		fflush(stdout);
		}
	while (nb_queue) {
		j = queue[0];
		for (k = 1; k < nb_queue; k++)
			queue[k - 1] = queue[k];
		nb_queue--;
		for (i = 0; i < G_len; i++) {
			g_i = (PERMUTATION_OP) G->s_i(i);
			if (g_i->s_obj_k() != PERMUTATION)
				return error("LABRA::orbit(): g_i->s_obj_k() != PERMUTATION");
			k = g_i->s_ii(j) - 1;
			if (!f_reached[k]) {
				if (f_verbose) {
					printf("%ld, ", k);
					fflush(stdout);
					}
				f_reached[k] = TRUE;
				queue[nb_queue++] = k;
				if (j == omega) {
					new_edge(omega, k, g_i);
					}
				else {
					s_KM_i(j)->mult(g_i, &tmp);
					new_edge(omega, k, &tmp);
					}
				}
			} /* next i */
		} /* while (nb_queue) */
	if (f_verbose) {
		printf("} orbit finished \n");
		fflush(stdout);
		}
	
	if (queue) {
		my_free(queue);
		queue = NIL;
		}
	if (f_reached) {
		my_free(f_reached);
		f_reached = NIL;
		}
	return erg;
}

#if TEXDOCU
INT LABRA_OB::in_test(PERMUTATION_OP g)
#endif
{
	PERMUTATION_OB h, g1, g2;
	INT i, j, k_el, deg;
	
#ifdef DEBUG_SIFT
	printf("in_test: ");
	g->println();
#endif
	g->copy(&g2);
	deg = s_degree_i();
	i = 0;
	j = s_base()->s_ii(i) - 1; // i-th base point b[i] is j
	while (i < deg - 1 && (k_el = g2.s_ii(j) - 1) == j) {
		/* while g2 stabilizes j */
		i++;
		j = s_base()->s_ii(i) - 1;
		}
	/* now: g2 stabilizes b_0, ... b_{i-1} 
	 * and b_i^g2 = k_el; j = b_i */
	while ((i < deg - 1) && path_exists_el(j, k_el)) {
		s_EM_i(k_el)->invers(&h);
		g2.mult(&h, &g1); // we sift the coset representative off
		if (i != 0)
			g1.mult(s_EM_i(j), &g2);
				/* now: g2 = g2 * EM[k_el]^-1 * EM[j] */
		else
			g1.copy(&g2);
				/* now: g2 = g2 * EM[k_el]^-1 */
		/* now: j^g2 = j */
		if (g2.s_ii(j) - 1 != j) {
			printf("i = %ld j = %ld k_el = %ld\n", i, j, k_el);
			printf("g = "); g2.println();
			printf("EM[k_el] = "); s_EM_i(k_el)->println();
			error("sift: j^g2 != j");
			return FALSE;
			}
		i++; // we fall deeper into the stabilizer chain
			// lets determine how deep:
		j = s_base()->s_ii(i) - 1;
		while (i < deg - 1 && (k_el = g2.s_ii(j) - 1) == j) {
			/* while g2 stabilizes b[i] = j */
			i++;
			j = s_base()->s_ii(i) - 1;
			}
		}
	if (i == deg - 1) {
		return TRUE;
		}
	return FALSE;
}

#if TEXDOCU
INT LABRA_OB::sift(PERMUTATION_OP g)
#endif
{
	PERMUTATION_OB h, g1, g2;
	INT i, j, k, k_el, m, deg;
	
#ifdef DEBUG_SIFT
	printf("sift: ");
	g->println();
#endif
	g->copy(&g2);
	deg = s_degree_i();
	i = 0;
	j = s_base()->s_ii(i) - 1; // i-th base point b[i] is j
	while (i < deg - 1 && (k_el = g2.s_ii(j) - 1) == j) {
		/* while g2 stabilizes j */
		i++;
		j = s_base()->s_ii(i) - 1;
		}
	/* now: g2 stabilizes b_0, ... b_{i-1} 
	 * and b_i^g2 = k_el; j = b_i */
	while ((i < deg - 1) && path_exists_el(j, k_el)) {
		s_EM_i(k_el)->invers(&h);
		g2.mult(&h, &g1); // we sift the coset representative off
		if (i != 0)
			g1.mult(s_EM_i(j), &g2);
				/* now: g2 = g2 * EM[k_el]^-1 * EM[j] */
		else
			g1.copy(&g2);
				/* now: g2 = g2 * EM[k_el]^-1 */
		/* now: j^g2 = j */
		if (g2.s_ii(j) - 1 != j) {
			printf("i = %ld j = %ld k_el = %ld\n", i, j, k_el);
			printf("g = "); g2.println();
			printf("EM[k_el] = "); s_EM_i(k_el)->println();
			return error("sift: j^g2 != j");
			}
		i++; // we fall deeper into the stabilizer chain
			// lets determine how deep:
		j = s_base()->s_ii(i) - 1;
		while (i < deg - 1 && (k_el = g2.s_ii(j) - 1) == j) {
			/* while g2 stabilizes b[i] = j */
			i++;
			j = s_base()->s_ii(i) - 1;
			}
		}
	if (i == deg - 1) {
		return OK;
		}
	/* now: no path exists; 
	 * have to fuse orbit of k_el with that of j */

	// the case that k_el is in a deeper orbit:
	k = s_inv_base_ii(k_el) - 1;
	while ( ((m = s_V_ii(k)) != k) && 
		(s_inv_base()->s_ii(m) > s_inv_base()->s_ii(j)) ) {

		// we change g2 in such a way that it no longer maps omega to k_el 
		// but to m, the base point of the other orbit.
		// (multiplication by elements of the branching does not 
		// remove information from the generator g2).
		s_KM_i(k_el)->invers(&h);
		g2.mult(&h, &g1);
		g1.copy(&g2);
		k_el = m;
		k = s_inv_base_ii(k_el) - 1;
		/* now j^g2 = k_el again */
		if (g2.s_ii(j) - 1 != k_el) {
			return error("sift (fusing): j^g2 != k_el");
			}
		}
	
	// the case that we meet the root of another orbit.
	if (s_V_ii(k) == k) {
		// join the root by an edge 
		new_edge(j, k_el, &g2);
		update_EM();
		}
	else {
		// now the other case: k_el comes from an earlier orbit
		s_KM_i(k_el)->copy(&h);
		new_edge(j, k_el, &g2);
		update_EM();
		g2.invers(&g1);
		h.mult(&g1, &g2);
			/* g2 = (old) KM[k_el] * g2^-1 */

			// we do not know whether one can reach j from within the orbit of V[k], 
			// but KM[k_el] * g2^-1 does it: (V[k]^KM[k_el])^(g2^-1) = k_el^(g2^-1) = j !!!
			// therefore we sift (add) this element into the branching:
		sift(&g2);
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::schreier(VECTOR_OP G, INT G_len, INT omega, INT f_verbose)
#endif
{
	PERMUTATION_OP g_i;
	PERMUTATION_OB g, g1, h, h1;
	INT i, j, j_el, k_el;
	INT f_very_verbose = FALSE;

	if (f_verbose) {
		printf("schreier-stab of %ld\n", omega);
		fflush(stdout);
		}
	/* schreier generators for first orbit element: */
	if (f_verbose) {
		printf("schreier-stab at first orbit element omega = %ld\n", omega);
		fflush(stdout);
		}
	for (i = 0; i < G_len; i++) {
		G->s_i(i)->copy(&g);
		if (g.s_obj_k() != PERMUTATION)
			return error("LABRA::schreier() g.s_obj_k() != PERMUTATION");
		if (f_very_verbose) {
			printf("generator %ld = \n", i);
			g.println();
			fflush(stdout);
			}
		k_el = g.s_ii(omega) - 1;
		if (k_el == omega) {
			if (!g.einsp())
				sift(&g);
			}
		else {
			s_EM_i(k_el)->invers(&h);
			if (omega != s_base()->s_ii(0) - 1)
				h.mult(s_EM_i(omega), &h1);
			else
				h.copy(&h1);
			g.mult(&h1, &g1);
			if (!g1.einsp())
				sift(&g1);
			}
		} /* next i */

	/* schreier generators for other orbit elements: */
	for (j = s_inv_base()->s_ii(omega) - 1 + 1; 
		j < s_degree_i(); j++) {
		j_el = s_base()->s_ii(j) - 1;
		if (path_exists_el(omega, j_el)) {
			if (f_verbose) {
				printf("schreier-stab at orbit element "
					"base(%ld) = %ld\n", j, j_el);
				fflush(stdout);
				}
			for (i = 0; i < G_len; i++) {
				g_i = (PERMUTATION_OP) G->s_i(i);
				if (g.s_obj_k() != PERMUTATION)
					return error("LABRA::schreier() g.s_obj_k() != PERMUTATION");
				if (f_very_verbose) {
					printf("generator %ld = \n", i);
					g_i->println();
					fflush(stdout);
					}

				s_EM_i(omega)->invers(&h1);
				h1.mult(s_EM_i(j_el), &g);
				k_el = g_i->s_ii(j_el) - 1;
				if (k_el == omega) {
					g.mult(g_i, &g1);
					if (!g1.einsp())
						sift(&g1);
					}
				else {
					s_EM_i(k_el)->invers(&h);
					if (omega != s_base()->s_ii(0) - 1)
						h.mult(s_EM_i(omega), &h1);
					else
						h.copy(&h1);
					g.mult(g_i, &g1);
					g1.mult(&h1, &g);
					if (!g.einsp())
						sift(&g);
					}
				} /* next i */
			}
		} /* next j */
	if (f_verbose) {
		printf("\n");
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::get_generators(VECTOR_OP V, INTEGER_OP V_len, INT j)
#endif
{
	INT i, i_el, m;

	V_len->m_i(0);
	for (i = j + 1; i < s_degree_i(); i++) {
		i_el = s_base()->s_ii(i) - 1;
		m = s_V_ii(i);
		if (m != i && m >= j) {
			if (!s_KM_i(i_el)->einsp()) {
				V->append_element(V_len, s_KM_i(i_el));
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::schreier_and_sift(VECTOR_OP G, INT G_len, 
	INT omega, INT f_verbose, 
	VECTOR_OP G1, INTEGER_OP G1_len)
#endif
{
	PERMUTATION_OP g_i;
	PERMUTATION_OB g, g1, h, h1;
	INT i, j, j_el, k_el;
	LABRA_OB L; // temporary labra for sifting the stabilizers of omega
		// sifting the generators reduces their number 
		// without loosing the generating property (generating the stabilizer)
	INT f_very_verbose = FALSE;

	if (f_verbose) {
		printf("schreier-stab of %ld\n", omega);
		printf("calling LABRA::Init_no_generators()\n");
		fflush(stdout);
		}
	L.Init_no_generators(s_degree_i(), s_f_verbose_i(), s_f_very_verbose_i());
	/* schreier generators for first orbit element: */
	if (f_verbose) {
		printf("schreier-stab at first orbit element omega = %ld\n", omega);
		fflush(stdout);
		}
	for (i = 0; i < G_len; i++) {
		G->s_i(i)->copy(&g);
		if (g.s_obj_k() != PERMUTATION)
			return error("LABRA::schreier_and_sift() g.s_obj_k() != PERMUTATION");
		if (f_very_verbose) {
			printf("generator %ld = \n", i);
			g.println();
			fflush(stdout);
			}
		k_el = g.s_ii(omega) - 1;
		if (k_el == omega) {
			if (!g.einsp())
				L.sift(&g);
			}
		else {
			s_EM_i(k_el)->invers(&h);
			if (omega != s_base()->s_ii(0) - 1)
				h.mult(s_EM_i(omega), &h1);
			else
				h.copy(&h1);
			g.mult(&h1, &g1);
				// g1 is the schreier stabilizer:
				// it leaves omega fixed.
			if (!g1.einsp())
				L.sift(&g1);
			}
		} /* next i */

	/* 
	 * schreier generators for other orbit elements: 
	 */

	// for all deeper base points b[j] which lie in the orbit of omega:
	for (j = s_inv_base()->s_ii(omega) - 1 + 1; j < s_degree_i(); j++) {
		j_el = s_base()->s_ii(j) - 1;
		if (path_exists_el(omega, j_el)) {
			if (f_verbose) {
				printf("schreier-stab at orbit element base(%ld) = %ld\n", j, j_el);
				fflush(stdout);
				}
			for (i = 0; i < G_len; i++) {
				g_i = (PERMUTATION_OP) G->s_i(i);
				if (g_i->s_obj_k() != PERMUTATION)
					return error("LABRA::schreier_and_sift() g_i->s_obj_k() != PERMUTATION");
				if (f_very_verbose) {
					printf("generator %ld = \n", i);
					g_i->println();
					fflush(stdout);
					}
				s_EM_i(omega)->invers(&h1);
				h1.mult(s_EM_i(j_el), &g);
					// g maps omega to j_el
				
				k_el = g_i->s_ii(j_el) - 1;
				if (k_el == omega) {
					g.mult(g_i, &g1);
						// g stabilizes omega
					if (!g1.einsp())
						L.sift(&g1);
					}
				else {
					s_EM_i(k_el)->invers(&h);
					if (omega != s_base()->s_ii(0) - 1)
						h.mult(s_EM_i(omega), &h1);
					else
						h.copy(&h1);
						// now: h1 moves k_el to omega
					
					g.mult(g_i, &g1);
					g1.mult(&h1, &g);
					if (!g.einsp())
						L.sift(&g);
					}
				} /* next i */
			}
		} /* next j */
	if (f_verbose) {
		printf("calling get_generators()\n");
		fflush(stdout);
		}
	G1->m_il(0);
	G1_len->m_i(0);
	L.get_generators(G1, G1_len, -1);
	if (f_verbose) {
		printf("finished with schreier_and_sift()\n");
		fflush(stdout);
		}
	return OK;
}

static INT jerrum2(LABRA_OP L, VECTOR_OP G1, INTEGER_OP G1_len, INT j, INT f_verbose);

#if TEXDOCU
INT LABRA_OB::jerrum(INT f_verbose)
#endif
{
	VECTOR_OB G1;
	INTEGER_OB G1_len;
	INT i, j, ii;
	INT erg = OK;
	INT f_verb = FALSE;
	
	erg += s_G()->copy(&G1);
	erg += s_nb_gen()->copy(&G1_len);
	if (G1.s_li() > 40)
		f_verb = TRUE;
	if (f_verb) {
		printf("labra_%ld(", G1.s_li());
		fflush(stdout);
		}
	for (i = 0; i < s_degree_i(); i++) {
		if (f_verb) {
			printf(".");
			if (((i + 1) % 10) == 0)
				printf("%ld ", i + 1);
			if (((i + 1) % 40) == 0)
				printf("\n");
			fflush(stdout);
			}
		j = s_base()->s_ii(i) - 1;
		if (G1_len.s_i() == 0)
			break;
		if (f_verbose) {
			printf("jerrum: base(%ld) = %ld; nb_gen = %ld\n", 
				i, j, G1_len.s_i());
			}
		if (s_f_very_verbose_i()) {
			for (ii = 0; ii < G1_len.s_i(); ii++) {
				G1.s_i(ii)->println();
				}
			}
		if (orbit(j, &G1, G1_len.s_i(), 
			f_verbose) != OK)
			return error("LABRA::jerrum() error in orbit()");
		update_EM();
#if 0
		/* this does NOT work for PSSL_2(25) 
		 * group order 15600 correct at first time, 
		 * then incorrectly 7800 
		 * after applying jerrum to the (reduced set of) 
		 * generators again !!! 
		 * the #else part works correctly. 
		 */
		if (schreier(&G1, G1_len.s_i(), j, 
			f_verbose) != OK)
			return error("LABRA::jerrum() error in schreier()");
		if (get_generators(&G1, &G1_len, i + 1) != OK)
			return error("LABRA::jerrum() error in get_generators()");
#else
		jerrum2(this, &G1, &G1_len, j, f_verbose);
#if 0
		{
		VECTOR_OB G2;
		INTEGER_OB G2_len;

		if (schreier_and_sift(&G1, G1_len.s_i(), j, f_verbose, 
			&G2, &G2_len) != OK)
			return error("LABRA::jerrum() error in schreier_and_sift()");
		G2.swap(&G1);
		G2_len.swap(&G1_len);
		}
#endif
#endif
		}
	calc_T();
	if (f_verb) {
		printf(")labra\n");
		fflush(stdout);
		}
	return erg;
}

static INT jerrum2(LABRA_OP L, VECTOR_OP G1, INTEGER_OP G1_len, INT j, INT f_verbose)
{
	VECTOR_OB G2;
	INTEGER_OB G2_len;

	if (L->schreier_and_sift(G1, G1_len->s_i(), j, f_verbose, 
		&G2, &G2_len) != OK)
		return error("jerrum2() error in schreier_and_sift()");
	G2.swap(G1);
	G2_len.swap(G1_len);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::Sym_n(INT deg, INT f_verbose)
#endif
{
	PERMUTATION_OB p;
	VECTOR_OB G;
	INT i;
	
	if (deg <= 1)
		return error("LABRA_OB::Sym_n(): deg <= 1");
	G.m_il(deg - 1);
	for (i = 2; i <= deg; i++) {
		p.m_il(deg);
		p.one();
		p.Add2Cycle(i - 1, i);
		p.swap(G.s_i(i - 2));
		}
	if (Init(deg, &G, G.s_li(), 
		TRUE /* f_verbose */, 
		FALSE /* f_very_verbose */) != OK)
		return error("LABRA_OB::Sym_n(): error in Init");
	if (jerrum(f_verbose) != OK)
		return error("LABRA_OB::Sym_n(): error in jerrum");
	return OK;
}

#if TEXDOCU
INT LABRA_OB::init_by_generators(VECTOR_OP G, INT deg, INT f_verbose)
#endif
{
	if (Init(deg, G, G->s_li(), 
		TRUE /* f_verbose */, 
		FALSE /* f_very_verbose */) != OK)
		return error("LABRA_OB::init_by_generators(): error in Init");
	if (jerrum(f_verbose) != OK)
		return error("LABRA_OB::init_by_generators(): error in jerrum");
	return OK;
}

#if TEXDOCU
INT test_jerrum_Sym(INT n, INT f_verbose)
#endif
{
	INT i;
	LABRA_OB L;
	SYM_OB go;
	INT t0, t1, user_time;
	BYTE str[256];
	
	for (i = 2; i <= n; i++) {
		t0 = os_ticks();
		L.Sym_n(i, f_verbose);
		L.group_order(&go);
		printf("|S_%ld| = ", i);
		go.println();
		t1 = os_ticks();
		user_time = t1 - t0;
		sprintf(str, "Running time for Sym_n(%ld): ", i);
		print_delta_time(user_time, str);
		printf("%s\n", str);
		}
	return OK;
}


#if 0
INT test_jerrum(void)
{
	LABRA_OB L;
	PERMUTATION_OB p1, p2, p3, p4, p5, p6, p7, p8;
	VECTOR_OB G;
	INTEGER_OB G_len;
	INT per_kind, deg;
	INT erg = OK;

	deg = 8;
	per_kind = ik_permutation(deg);
	G_len.m_i(0);
	G.m_il(VECTOR_OVERSIZE);
	((SYM_OP) &p1)->init(per_kind);
	((SYM_OP) &p2)->init(per_kind);
	((SYM_OP) &p3)->init(per_kind);
	((SYM_OP) &p4)->init(per_kind);
	((SYM_OP) &p5)->init(per_kind);
	((SYM_OP) &p6)->init(per_kind);
	((SYM_OP) &p7)->init(per_kind);
	((SYM_OP) &p8)->init(per_kind);
	p1.one();
	p2.one();
	p3.one();
	p4.one();
	p5.one();
	p6.one();
	p7.one();
	p8.one();
	p1.Add2Cycle(1, 2);
	p2.Add2Cycle(2, 3);
	p3.Add2Cycle(3, 4);
	p4.Add2Cycle(4, 5);
	p5.Add2Cycle(5, 6);
	p6.Add2Cycle(6, 7);
	p7.Add2Cycle(7, 8);
	G.append_element(&G_len, &p1);
	G.append_element(&G_len, &p2);
	G.append_element(&G_len, &p3);
	G.append_element(&G_len, &p4);
	G.append_element(&G_len, &p5);
	G.append_element(&G_len, &p6);
	G.append_element(&G_len, &p7);
	erg += L.Init(deg, per_kind, &G, G_len.s_i(), 
		TRUE /* f_verbose */, 
		FALSE /* f_very_verbose */, 
		FALSE /* f_permh */);
	L.jerrum();
	L.Print();
	L.Print_T(FALSE);
	return erg;
}
#endif

#define VERDREHE_N_MAL 10000

#if 0
INT rubicks_cube(void)
{
	LABRA_OB L;
	PERMUTATION_OB p1, p2, p3, p4, p5, p6, q, q1;
	VECTOR_OB G, V, D, Dv;
	INTEGER_OB V_len, D_len;
	INT per_kind, deg;
	INT erg = OK;
	INT r, i, j;

	deg = 54;
	per_kind = ik_permutation(deg);
	G.m_il(6);
	((SYM_OP) &p1)->init(per_kind);
	((SYM_OP) &p2)->init(per_kind);
	((SYM_OP) &p3)->init(per_kind);
	((SYM_OP) &p4)->init(per_kind);
	((SYM_OP) &p5)->init(per_kind);
	((SYM_OP) &p6)->init(per_kind);
	((SYM_OP) &q)->init(per_kind);
	((SYM_OP) &q1)->init(per_kind);
	p1.one();
	p2.one();
	p3.one();
	p4.one();
	p5.one();
	p6.one();
	q.one();
	q1.one();
	
	/* oben */
	p1.Add4Cycle(1, 7, 9, 3);
	p1.Add4Cycle(2, 4, 8, 6);
	p1.Add4Cycle(19, 10, 37, 46);
	p1.Add4Cycle(47, 20, 11, 38);
	p1.Add4Cycle(48, 21, 12, 39);
	
	/* vorne */
	p2.Add4Cycle(10, 16, 18, 12);
	p2.Add4Cycle(11, 13, 17, 15);
	p2.Add4Cycle(7, 27, 30, 37);
	p2.Add4Cycle(8, 24, 29, 40);
	p2.Add4Cycle(21, 28, 43, 9);
	
	/* links */
	p3.Add4Cycle(19, 25, 27, 21);
	p3.Add4Cycle(20, 22, 26, 24);
	p3.Add4Cycle(1, 54, 28, 10);
	p3.Add4Cycle(4, 51, 31, 13);
	p3.Add4Cycle(48, 34, 16, 7);
	
	/* unten */
	p4.Add4Cycle(28, 34, 36, 30);
	p4.Add4Cycle(29, 31, 35, 33);
	p4.Add4Cycle(16, 25, 52, 43);
	p4.Add4Cycle(17, 26, 53, 44);
	p4.Add4Cycle(27, 54, 45, 18);
	
	/* rechts */
	p5.Add4Cycle(37, 43, 45, 39);
	p5.Add4Cycle(38, 40, 44, 42);
	p5.Add4Cycle(9, 18, 36, 46);
	p5.Add4Cycle(6, 15, 33, 49);
	p5.Add4Cycle(12, 30, 52, 3);
	
	/* hinten */
	p6.Add4Cycle(46, 52, 54, 48);
	p6.Add4Cycle(47, 49, 53, 51);
	p6.Add4Cycle(3, 45, 34, 19);
	p6.Add4Cycle(2, 42, 35, 22);
	p6.Add4Cycle(39, 36, 25, 1);
	
	p1.invers(&q); q.println();
	p2.invers(&q); q.println();
	p3.invers(&q); q.println();
	p4.invers(&q); q.println();
	p5.invers(&q); q.println();
	p6.invers(&q); q.println();
	
	p1.copy((PERMUTATION_OP) G.s_i(0));
	p2.copy((PERMUTATION_OP) G.s_i(1));
	p3.copy((PERMUTATION_OP) G.s_i(2));
	p4.copy((PERMUTATION_OP) G.s_i(3));
	p5.copy((PERMUTATION_OP) G.s_i(4));
	p6.copy((PERMUTATION_OP) G.s_i(5));
	erg += L.Init(deg, per_kind, &G, G.s_li(), 
		TRUE /* f_verbose */, 
		FALSE /* f_very_verbose */);
	L.jerrum(TRUE /* f_verbose */);
	L.Print();
	L.Print_T(FALSE);
	L.print_group_order();
	
	printf("verdrehe Wuerfel:\n");
	r = os_ticks();
	my_srand(r);
	for (i = 0; i < VERDREHE_N_MAL; i++) {
		j = my_rand() % 6;
		if (j == 0)
			q.mult(&p1, &q1);
		else if (j == 1)
			q.mult(&p2, &q1);
		else if (j == 2)
			q.mult(&p3, &q1);
		else if (j == 3)
			q.mult(&p4, &q1);
		else if (j == 4)
			q.mult(&p5, &q1);
		else if (j == 5)
			q.mult(&p6, &q1);
		else
			printf("j = %ld\n", j);
		q1.copy(&q);
		}
	q.println();
#if 0
	L.factor(&q, &V, &V_len);
	printf("factored !\n");
	L.write_out(&V, &V_len, &D, &Dv, &D_len);
#endif
	return erg;
}
#endif

/* es wird die folgende Numerierung verwendet:

                     oben

      links       vorne      rechts

                     unten

                     hinten


                     1   2   3
                     4   5   6
                     7   8   9

19 20 21  |  10 11 12  |  37 38 39
22 23 24  |  13 14 15  |  40 41 42
25 26 27  |  16 17 18  |  43 44 45

                   28 29 30
                   31 32 33
                   34 35 36

                   54 53 52
                   51 50 49
                   48 47 46
*/

#ifndef RAND_MAX
#define RAND_MAX 32767
#endif

static INT next = 1;

#if TEXDOCU
void my_srand(INT seed)
#endif
{
	next = seed;
}

#if TEXDOCU
INT my_rand(void)
#endif
/* 0 <= r <= RAND_MAX */
{
	next = next * 1103515245 + 12345;
	return ((next >> 16) & RAND_MAX);
}

#if TEXDOCU
INT lb_test_generators_from_file(BYTE *fname)
#endif
{
	VECTOR_OB G, G1;
	INTEGER_OB G1_len;
	PERMUTATION_OB p;
	PERMUTATION_OP p_perm;
	INT nb_gen, deg;
	INT i, j, a;
	FILE *fp;
	LABRA_OB L;
	BYTE str[1024];

	fp = fopen(fname, "r");
	fscanf(fp, "%ld", &nb_gen);
	G.m_il(nb_gen);
	printf("nb_gen = %ld\n", nb_gen);
	fscanf(fp, "%ld", &deg);
	printf("deg = %ld\n", deg);
	for (i = 0; i < nb_gen; i++) {
		p.m_il(deg);
		for (j = 0; j < deg; j++) {
			fscanf(fp, "%ld", &a);
			p.m_ii(j, a);
			}
		p.println();
		p.swap((PERMUTATION_OP) G.s_i(i));
		}
	fclose(fp);
	L.Init(deg, &G, nb_gen, TRUE, FALSE);
	printf("calling jerrum...\n");
	fflush(stdout);
	L.jerrum(TRUE);
	L.print_group_order();
	G1.m_il(0);
	L.get_generators(&G1, &G1_len, -1);

	strcpy(str, fname);
	strcat(str, ".out");
	fp = fopen(str, "w");
	fprintf(fp, "%ld %ld\n", G1_len.s_i(), deg);
	for (i = 0; i < G1_len.s_i(); i++) {
		p_perm = (PERMUTATION_OP) G1.s_i(i);
		p_perm->fprint_list(fp);
		}
	fclose(fp);

	return OK;
}

#if TEXDOCU
INT LABRA_OB::all_sons(INT i, VECTOR_OP sons)
#else
//PRE
// sons of i are immediate descendants of i
// (there exists an edge $\pi_i -> \pi_j$
// sons becomes an integer vectors of indices of base points: 
// $\{ KM(\pi(j)) | j = sons[i], i < sons->s_li() \}$
// is the set of edges leading to sons.
///PRE
#endif
{
	INT n = s_degree_i();
	INT j, nb_s = 0;
	
	sons->m_il(n);
	for (j = i + 1; j < n; j++) {
		if (s_V_ii(j) == i) {
			sons->m_ii(nb_s, j);
			nb_s++;
			}
		}
	sons->realloc_z(nb_s);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::the_son(INT j, PERMUTATION_OP p)
#endif
{
	INT j_el;

	j_el = s_base_ii(j) - 1;
	s_KM_i(j_el)->copy(p);
	return OK;
}

#if TEXDOCU
INT LABRA_OB::all_stab_generators(INT i, VECTOR_OP gens)
#else
//PRE
// the group fixing $\pi_1, ... , \pi_i$ is returned.
// gens contains a vector of indices of base points 
// such that
// $\{ KM(\pi_j) | j = gens[i], i < gens->s_li() \}$
// is a generating set.
///PRE
#endif
{
	INT n = s_degree_i();
	INT j, nb = 0;
	
	gens->m_il(n);
	for (j = i + 1; j < n; j++) {
		if (s_V_ii(j) != j && s_V_ii(j) > i) {
			gens->m_ii(nb, j);
			nb++;
			}
		}
	gens->realloc_z(nb);
	return OK;
}

#if TEXDOCU
INT orbits_cheap(VECTOR_OP gen, VECTOR_OP tda_ordering, 
	VECTOR_OP orbit_first, VECTOR_OP orbit_length, INT f_v)
#endif
{
	INT set_size, *Q = NIL, Q_len;
	INT i, j, k, cur, next, os;
	PERMUTATION_OP g;
	VECTOR_OB SVorbit;
	INT nb_orbits = 0, n = 0, n0;

	if (gen->s_li() == 0)
		return error("LABRA::orbits_cheap(): gen->s_li == 0");
	g = (PERMUTATION_OP) gen->s_i(0);
	set_size = g->s_li();
	Q = (INT *) my_malloc(set_size * sizeof(INT), "LABRA::orbits_cheap()");
	if (Q == NIL)
		return error("Q == NIL");
	SVorbit.m_il(set_size);
	tda_ordering->m_il(set_size);
	for (i = 0; i < set_size; i++)
		SVorbit.m_ii(i, -1);
	orbit_first->m_il(set_size);
	orbit_length->m_il(set_size);

	for (i = 0; i < set_size; i++) {
		if (SVorbit.s_ii(i) != -1)
			continue;
		SVorbit.m_ii(i, nb_orbits);
		tda_ordering->m_ii(n, i);
		n0 = n;
		n++;
		
		os = 1;
		Q[0] = i;
		Q_len = 1;
		while (Q_len) {
			cur = Q[0];
			for (k = 1; k < Q_len; k++)
				Q[k - 1] = Q[k];
			Q_len--;
			for (j = 0; j < gen->s_li(); j++) {
				g = (PERMUTATION_OP) gen->s_i(j);
				next = g->s_ii(cur) - 1;
				if (SVorbit.s_ii(next) != -1)
					continue;
				os++;
				SVorbit.m_ii(next, nb_orbits);
				tda_ordering->m_ii(n, next);
				n++;
				Q[Q_len++] = next;
				}
			}
		orbit_first->m_ii(nb_orbits, n0);
		orbit_length->m_ii(nb_orbits, os);
		nb_orbits++;
		}
	orbit_first->realloc_z(nb_orbits);
	orbit_length->realloc_z(nb_orbits);
	
	if (Q) {
		my_free(Q);
		Q = NIL;
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::orbits(VECTOR_OP tda_ordering, 
	VECTOR_OP orbit_first, VECTOR_OP orbit_length, 
	INT f_calc_labras, VECTOR_OP O_labra, INT f_v)
#endif
{
	INT i, i_el, j, n, n0, ol, nb_orbits;
	VECTOR_OB chi;
	VECTOR_OB f_visited; 

	n = s_degree_i();
	tda_ordering->m_il(n);
	f_visited.m_il_n(n);
	orbit_first->m_il(n);
	orbit_length->m_il(n);
	if (f_calc_labras)
		O_labra->m_il(n);
	nb_orbits = 0;
	n0 = 0;
	printf("orbits: %ld = ", n); fflush(stdout);
	if (!s_base()->einsp())
		return error("LABRA::orbits() base is not id !");
	while (n0 < n) {
		if (f_v) {
			printf("get the first unvisited element; f_visited = ");
			f_visited.println();
			}
		for (i = 0; i < n; i++) {
			if (f_visited.s_ii(i) == 0) {
				if (f_v) {
					printf("next, we compute the orbit of %ld\n", i);
					fflush(stdout);
					}
				cycle(0, s_inv_base_ii(i) - 1, FALSE, FALSE);
				break;
				}
			}
		orbit_first->m_ii(nb_orbits, n0);
		if (f_v) {
			printf("new orbit\n");
			fflush(stdout);
			}
		one_orbit(&chi);
			// computes the base indices of the one-orbit 
			// (orbit of \pi_0) into chi: 
			// chi(i) is 1 iff \pi_i is in the orbit of \pi_0
			// chi(0) is always 1 !
		if (f_v) {
			printf("chi = ");
			chi.println();
			fflush(stdout);
			}
		
		ol = 0;
		for (i = 0; i < n; i++) {
			if (chi.s_ii(i)) {
				i_el = s_base_ii(i) - 1;
				tda_ordering->m_ii(n0, i_el);
				f_visited.m_ii(i_el, 1);
				n0++;
				if (i != ol) {
					cycle(ol, i, FALSE, FALSE);
					}
				ol++;
				}
			}
		printf(" + %ld", ol); fflush(stdout);
		if (f_v) {
			printf("found an orbit of length %ld\n", ol);
			for (i = 0; i < ol; i++) {
				printf("%ld ", s_base_ii(i) - 1);
				}
			printf("\n");
			fflush(stdout);
			}
		orbit_length->m_ii(nb_orbits, ol);
		if (f_calc_labras)
			copy((LABRA_OP) O_labra->s_i(nb_orbits));
		nb_orbits++;
		}
	printf("\n");
	orbit_first->realloc_z(nb_orbits);
	orbit_length->realloc_z(nb_orbits);
	if (f_calc_labras)
		O_labra->realloc_z(nb_orbits);
	
	for (i = 0; i < n; i++) {
		j = tda_ordering->s_ii(i);
		if (s_base_ii(i) - 1 != j) {
			cycle(i, s_inv_base_ii(j) - 1, FALSE, FALSE);
			}
		}
	return OK;
}

#if TEXDOCU
INT LABRA_OB::one_orbit(VECTOR_OP chi)
#endif
{
	INT j, n;

	n = s_degree_i();
	chi->m_il_n(n);
	for (j = 0; j < n; j++) {
		if (path_exists(0, j))
			chi->m_ii(j, 1);
		}
	return OK;
}

#include "lb_set_canon.C"

#if TEXDOCU
INT LABRA_OB::set_canon(VECTOR_OP X, VECTOR_OP X0, INT f_v, INT f_vv)
#endif
{
	INT *XX, *YY;
	INT i, l, ago;
	SYM_OB go;
	
	l = X->s_li();
	if (l == 0) {
		X->copy(X0);
		group_order(&go);
		ago = go.s_i_i();
		return ago;
		}
#if 0
	if (f_v) {
		printf("LABRA_OB::set_canon() l=%ld X=", l);
		X->println();
		fflush(stdout);
		}
#endif
	XX = (INT *) my_malloc(l * sizeof(INT), "LABRA::set_canon(): XX");
	YY = (INT *) my_malloc(l * sizeof(INT), "LABRA::set_canon(): YY");
	for (i = 0; i < l; i++) 
		XX[i] = X->s_ii(i);
	ago = calc_canonical_form_of_set(this, l, XX, YY, f_v, f_vv);
	X0->m_il(l);
	for (i = 0; i < l; i++) 
		X0->m_ii(i, YY[i]);
	my_free(XX);
	my_free(YY);
	return ago;
}

#if TEXDOCU
INT LABRA_OB::subset_orbits(VECTOR_OP X, VECTOR_OP Orbits, 
	VECTOR_OP Orbit_lengths, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2, INT t, INT *p_k_idx)
#endif
{
	INT i, j, k, n_idx, n, len, len1, ago, ii, iii;
	INT *choice;
	VECTOR_OB N, N0, X0;
	VECTOR_OP pOk, pOlenk, pObelow1k, pObelow2k, ppObelow1k, ppObelow2k;
	VECTOR_OP pO, pOlen, pObelow1, pObelow2, ppObelow1;
	INTEGER_OB int_ob;
	INT f_v = TRUE;
	INT f_vv = FALSE;
	INT k_idx, f_found_k;
	INT f_found;
	INT idx1, f_found1;
	
	k = X->s_li();
	ago = set_canon(X, &X0, f_vv, FALSE);
	pOk = (VECTOR_OP) Orbits->s_i(k);
	pOlenk = (VECTOR_OP) Orbit_lengths->s_i(k);
	pObelow1k = (VECTOR_OP) Orbits_below1->s_i(k);
	pObelow2k = (VECTOR_OP) Orbits_below2->s_i(k);
	pOk->search(pOk->s_li(), TRUE /* f_ascending */, &X0, &k_idx, &f_found_k);
	if (f_found_k) {
		*p_k_idx = k_idx - 1;
		return OK;
		}
	else {
		if (f_v) {
			printf("adding %ld orbit: X=", k);
			X->print();
			printf(" X0=");
			X0.println();
			}
		len = pOk->s_li();
		pOk->inc();
		pOlenk->inc();
		pObelow1k->inc();
		pObelow2k->inc();
		for (i = len; i > k_idx; i--) {
			pOk->s_i(i)->swap(pOk->s_i(i - 1));
			pOlenk->s_i(i)->swap(pOlenk->s_i(i - 1));
			pObelow1k->s_i(i)->swap(pObelow1k->s_i(i - 1));
			pObelow2k->s_i(i)->swap(pObelow2k->s_i(i - 1));
			}
		X0.copy((VECTOR_OP) pOk->s_i(k_idx));
		pOlenk->m_ii(k_idx, ago);
		ppObelow1k = (VECTOR_OP) pObelow1k->s_i(k_idx);
		ppObelow2k = (VECTOR_OP) pObelow2k->s_i(k_idx);
		ppObelow1k->m_il(0);
		ppObelow2k->m_il(0);
		if (f_vv) {
			printf("%ld orbits:\n", k);
			pOk->Print();
			}
		}
	*p_k_idx = k_idx;
	if (k == t) {
		return OK;
		}
	
	choice = (INT *) my_malloc(k * sizeof(INT), "LABRA_OB::subset_orbits(): choice");
	n = k - 1;
	if (f_v) {
		printf("%ld-subsets of ", n);
		X->print();
		printf(" = ");
		X0.println();
		}
	N.m_il(n);
	n_Choose_k_first(choice, k, n);
	while (TRUE) {
		for (i = 0; i < n; i++) {
			j = choice[i];
			N.m_ii(i, X0.s_ii(j));
			}
		if (f_v) {
			printf("N=");
			N.print();
			}
		ago = set_canon(&N, &N0, f_vv, FALSE);
		if (f_v) {
			printf(" N0=");
			N0.print();
			printf(" ago=%ld\n", ago);
			fflush(stdout);
			}
		pO = (VECTOR_OP) Orbits->s_i(n);
		pOlen = (VECTOR_OP) Orbit_lengths->s_i(n);
		pObelow1 = (VECTOR_OP) Orbits_below1->s_i(n);
		pObelow2 = (VECTOR_OP) Orbits_below2->s_i(n);
		pO->search(pO->s_li(), TRUE /* f_ascending */, 
			&N0, &n_idx, &f_found);
		if (f_found) {
			n_idx--;
			if (pOlen->s_ii(n_idx) != ago)
				return error("pOlen->s_ii(n_idx) != ago");
			}
		else {
			// recursively add orbits below N0:
			subset_orbits(&N0, Orbits, Orbit_lengths, 
				Orbits_below1, Orbits_below2, t, &n_idx);
			if (f_vv) {
				printf("N0 added with n_idx = %ld\n", n_idx);
				printf("N0=");
				N0.println();
				printf("%ld-orbits:\n", n);
				((VECTOR_OP) Orbits->s_i(n))->Print();
				printf("updating %ld orbits below1:\nbefore:", k);
				((VECTOR_OP) Orbits_below1->s_i(k))->Print();
				}
			// update the orbit numbers in the k orbit layer:
			len = pOk->s_li();
			for (ii = 0; ii < len; ii++) {
				ppObelow1 = (VECTOR_OP) pObelow1k->s_i(ii);
				len1 = ppObelow1->s_li();
				for (iii = 0; iii < len1; iii++) {
					if (ppObelow1->s_ii(iii) >= n_idx)
						ppObelow1->s_i(iii)->inc();
					}
				}
			if (f_vv) {
				printf("after:\n");
				((VECTOR_OP) Orbits_below1->s_i(k))->Print();
				}

			}

		// now n_idx is the index of the n-orbit.
		// update the below info of the k orbit:
		int_ob.m_i(n_idx);
		ppObelow1k->search(ppObelow1k->s_li(), TRUE, 
			&int_ob, &idx1, &f_found1);
		if (f_found1) {
			idx1--;
			ppObelow2k->s_i(idx1)->inc();
			}
		else {
			len = ppObelow1k->s_li();
			ppObelow1k->inc();
			ppObelow2k->inc();
			for (ii = len; ii > idx1; ii--) {
				ppObelow1k->s_i(ii)->swap(ppObelow1k->s_i(ii - 1));
				ppObelow2k->s_i(ii)->swap(ppObelow2k->s_i(ii - 1));
				}
			ppObelow1k->m_ii(idx1, n_idx);
			ppObelow2k->m_ii(idx1, 1);
			}

		if (!n_Choose_k_next(choice, k, n))
			break;
		}


	my_free(choice);
	return OK;
}

INT LABRA_OB::group_table(MATRIX_OP table, VECTOR_OP table_inv, 
	VECTOR_OP elements, VECTOR_OP generators, INT f_v)
{
	SYM_OB go;
	INT n, deg;
	VECTOR_OB path, E, Eidx;
	PERMUTATION_OB p;
	PERMUTATION_OP pp, pq;
	INT i, j, idx, f_found;
	
	group_order(&go);
	if (go.s_obj_k() != INTEGER) {
		return error("LABRA_OB::group_table() group order too large (not even an INTEGER)\n");
		}
	n = go.s_i_i();
	deg = s_degree_i();
	printf("LABRA_OB::group_table() a group of order %ld of degree %ld\n", n, deg);
	elements->m_il(n);
	i = 0;
	elements_first(&path, &p);
	while (TRUE) {
		p.copy((PERMUTATION_OP) elements->s_i(i));
		i++;
		if (!elements_next(&path, &p))
			break;
		}
	if (i != n)
		return error("LABRA_OB::group_table()  i != n");
	elements->copy(&E);
	E.quicksort(E.s_li(), TRUE);
	Eidx.m_il(n);
	for (i = 0; i < n; i++) {
		j = elements->search_linear(E.s_i(i));
		if (j < 0) 
			return error("LABRA_OB::group_table() element not found");
		Eidx.m_ii(i, j);
		}
	if (f_v) {
		printf("elements:\n");
		elements->Print();
		printf("elements sorted:\n");
		E.Print();
		printf("element indices:\n");
		Eidx.println();
		}
	table->m_ilih(n, n);
	for (i = 0; i < n; i++) {
		pp = (PERMUTATION_OP) elements->s_i(i);
		for (j = 0; j < n; j++) {
			pq = (PERMUTATION_OP) elements->s_i(j);
			pp->mult(pq, &p);
			E.search(n, TRUE /* f_ascending */, &p, &idx, &f_found);
			if (!f_found)
				return error("LABRA_OB::group_table() group element not found");
			idx--;
			table->m_iji(i, j, Eidx.s_ii(idx));
			}
		}
	if (f_v) {
		printf("table:\n");
		table->Print();
		}
	table_inv->m_il(n);
	for (i = 0; i < n; i++) {
		pp = (PERMUTATION_OP) elements->s_i(i);
		pp->invers(&p);
		E.search(n, TRUE /* f_ascending */, &p, &idx, &f_found);
		if (!f_found)
			return error("LABRA_OB::group_table() inverse group element not found");
		idx--;
		table_inv->m_ii(i, Eidx.s_ii(idx));
		}
	if (f_v) {
		printf("table_inv:\n");
		table_inv->Print();
		}
	
	INT nb_gen = s_nb_gen_i();
	generators->m_il(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		pp = (PERMUTATION_OP) s_G_i(i);
		E.search(n, TRUE /* f_ascending */, pp, &idx, &f_found);
		if (!f_found)
			return error("LABRA_OB::group_table() generator element not found");
		idx--;
		generators->m_ii(i, Eidx.s_ii(idx));
		}
	if (f_v) {
		printf("generators:\n");
		table_inv->println();
		}
	return OK;
}
 



#endif /* LABRA_TRUE */


