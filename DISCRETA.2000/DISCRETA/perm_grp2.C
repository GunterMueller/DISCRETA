/* perm_grp2.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef PERMTRUE

#include <DISCRETA/ma.h> /* for PGL_2_p_as_perm() */
#include <DISCRETA/perm.h>
#include <DISCRETA/part.h>
#include <DISCRETA/fga.h> /* for perm_vec_get_degree() */


#define VEC_GENERATORS_EXPONENTIATION_DEBUG
#undef PERM_ON_EXPONENTIATION_DEBUG

#if TEXDOCU
INT vec_generators_diagonal_sum(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c)
#else
given generators a and b for permutation groups $G$ and $G$ on disjoint sets 
$X$ and $Y$, this routine computes generators for the disjoint action of $G$ on 
$X \cup Y$.
#endif
{
	PERMUTATION_OP pa, pb, pc;
	INT i, la, lb, da, db;

	la = a->s_li();
	lb = b->s_li();
	if (la == 0)
		return error("vec_generators_diagonal_sum() la == 0");
	if (la != lb)
		return error("vec_generators_diagonal_sum() la != lb");
	
	pa = (PERMUTATION_OP) a->s_i(0);
	da = pa->s_li();
	if (lb == 0)
		return error("vec_generators_diagonal_sum() lb == 0");
	pb = (PERMUTATION_OP) b->s_i(0);
	db = pb->s_li();
	
	c->m_il(la);
	for (i = 0; i < la; i++) {
		pa = (PERMUTATION_OP) a->s_i(i);
		pb = (PERMUTATION_OP) b->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(i);
		pa->join(pb, pc);
		}
	return OK;
}

#if TEXDOCU
INT vec_generators_comma(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c)
#else
Concatenates the lists (vectors) of generators in a and b. 
Result is c. Looks for the degree: if needed, one of the 
generating sets is embedded to get the largest degree. 
The fixepoints are added at the end.
#endif
{
	PERMUTATION_OP pa, pb, pc;
	INT i, k, la, lb, lc, da, db, dc;

	la = a->s_li();
	lb = b->s_li();
	lc = la + lb;
	if (la == 0)
		return error("vec_generators_comma() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	da = pa->s_li();
	if (lb == 0)
		return error("vec_generators_comma() lb == 0");
	pb = (PERMUTATION_OP) b->s_i(0);
	db = pb->s_li();
	dc = MAXIMUM(da, db);
	c->m_il(lc);
	k = 0;
	for (i = 0; i < la; i++) {
		pa = (PERMUTATION_OP) a->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(k);
		pa->add_n_fixpoints_at_end(pc, dc - da);
		k++;
		}
	for (i = 0; i < lb; i++) {
		pb = (PERMUTATION_OP) b->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(k);
		pb->add_n_fixpoints_at_end(pc, dc - db);
		k++;
		}
	return OK;
}

#if TEXDOCU
INT vec_generators_direct_sum(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c)
#else
given generators a and b for permutation groups $G$ and $H$ on disjoint sets 
$X$ and $Y$, this routine computes generators for the disjoint action of 
$X \cup Y$ (the action of the direct product of $G$ and $H$). 
#endif
{
	PERMUTATION_OP pa, pb, pc;
	INT i, la, lb, lc, da, db;

	la = a->s_li();
	lb = b->s_li();
	lc = la + lb;
	if (la == 0)
		return error("vec_generators_direct_sum() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	da = pa->s_li();
	if (lb == 0)
		return error("vec_generators_direct_sum() lb == 0");
	pb = (PERMUTATION_OP) b->s_i(0);
	db = pb->s_li();
	c->m_il(lc);
	for (i = 0; i < la; i++) {
		pa = (PERMUTATION_OP) a->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(i);
		pa->add_n_fixpoints_at_end(pc, db);
		}
	for (i = 0; i < lb; i++) {
		pb = (PERMUTATION_OP) b->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(la + i);
		pb->add_n_fixpoints_in_front(pc, da);
		}
	return OK;
}

#if TEXDOCU
INT vec_generators_direct_product(VECTOR_OP a, VECTOR_OP b, VECTOR_OP c)
#else
Given generators for Permutation groups $G$ acting on $X$ 
and $H$ acting of $Y$, this routine computes the action of 
the direct product on $\{(x,y) \in X \times Y\}$ via 
$(g \times h) \cdot (x,y) \mapsto (x^g, y^h)$.
Each generator for $G$ or $H$ results in a generator  
for $G \times H$ which moves only rows (or columns) of $X \times Y$.
#endif
{
	PERMUTATION_OP pa, pb, pc;
	PERMUTATION_OB id_n, id_m;
	INT i, la, lb, lc, n, m;

	la = a->s_li();
	lb = b->s_li();
	lc = la + lb;
	if (la == 0)
		return error("vec_generators_direct_product() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	n = pa->s_li();
	if (lb == 0)
		return error("vec_generators_direct_product() lb == 0");
	pb = (PERMUTATION_OP) b->s_i(0);
	m = pb->s_li();
	c->m_il(lc);
	id_n.m_il(n);
	id_n.one();
	id_m.m_il(m);
	id_m.one();
	for (i = 0; i < la; i++) {
		pa = (PERMUTATION_OP) a->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(i);
		perm_on_cartesian_product(pa, &id_m, pc);
		}
	for (i = 0; i < lb; i++) {
		pb = (PERMUTATION_OP) b->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(la + i);
		perm_on_cartesian_product(&id_n, pb, pc);
		}
	return OK;
}

#if TEXDOCU
INT perm_on_cartesian_product(PERMUTATION_OP a, PERMUTATION_OP b, PERMUTATION_OP c)
#else
Utility routine for vec\_generators\_direct\_product().
a acts on the rows, b on the columns of the cartesian product.
#endif
{
	INT n, m, nm, i, j, rank, i1, j1, rank1;

	n = a->s_li();
	m = b->s_li();
	nm = n * m;
	c->m_il(nm);
	for (rank = 0; rank < nm; rank++) {
		i = rank / m;
		j = rank % m;
		i1 = a->s_ii(i) - 1;
		j1 = b->s_ii(j) - 1;
		rank1 = i1 * m + j1;
		c->m_ii(rank, rank1 + 1);
		}
	return OK;
}

#if TEXDOCU
INT vec_generators_exponentiation(VECTOR_OP a, VECTOR_OP b, 
	VECTOR_OP c, INT f_simultaneously)
#endif
{
	PERMUTATION_OP pa, pb, pc;
	PERMUTATION_OB id_n, id_m;
	INT ii, i, k, la, lb, lc, n, m;

	la = a->s_li();
	lb = b->s_li();
	if (la == 0)
		return error("vec_generators_exponentiation() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	n = pa->s_li();
	if (lb == 0)
		return error("vec_generators_exponentiation() lb == 0");
	pb = (PERMUTATION_OP) b->s_i(0);
	m = pb->s_li();
	if (f_simultaneously)
		lc = la + lb;
	else
		lc = la * m +  lb;
	c->m_il(lc);
	id_n.m_il(n);
	id_n.one();
	id_m.m_il(m);
	id_m.one();
	k = 0;
	for (i = 0; i < la; i++) {
		pa = (PERMUTATION_OP) a->s_i(i);
		if (f_simultaneously) {
			pc = (PERMUTATION_OP) c->s_i(k);
			perm_on_exponentiation(pa, 0, &id_m, pc, TRUE);
			k++;
			}
		else {
			for (ii = 0; ii < m; ii++) {
				pc = (PERMUTATION_OP) c->s_i(k);
				perm_on_exponentiation(pa, ii, &id_m, pc, FALSE);
				k++;
				}
			}
		}
	for (i = 0; i < lb; i++) {
		pb = (PERMUTATION_OP) b->s_i(i);
		pc = (PERMUTATION_OP) c->s_i(k);
		perm_on_exponentiation(&id_n, 0, pb, pc, FALSE);
		k++;
		}
#ifdef VEC_GENERATORS_EXPONENTIATION_DEBUG
	printf("vec_generators_exponentiation() f_simultaneously = %ld\n", f_simultaneously);
	k = 0;
	for (i = 0; i < la; i++) {
		pa = (PERMUTATION_OP) a->s_i(i);
		printf("original generator of group a:\n");
		pa->println();
		printf("becomes:\n");
		if (f_simultaneously) {
			pc = (PERMUTATION_OP) c->s_i(k);
			pc->println();
			k++;
			}
		else {
			for (ii = 0; ii < m; ii++) {
				pc = (PERMUTATION_OP) c->s_i(k);
				pc->println();
				k++;
				}
			}
		} // next i
	for (i = 0; i < lb; i++) {
		printf("generators of group b:\n");
		pb = (PERMUTATION_OP) b->s_i(i);
		pb->println();
		printf("becomes:\n");
		pc = (PERMUTATION_OP) c->s_i(k);
		pc->println();
		k++;
		}
#endif
	return OK;
}

#if TEXDOCU
INT perm_on_exponentiation(PERMUTATION_OP a, INT a_in_component, 
	PERMUTATION_OP b, PERMUTATION_OP c, INT f_simultaneously)
#else
a acts on the rows, b on the columns of the matrix of the mapping of m to n
#endif
{
	VECTOR_OB R, R1;
	INT n, m, nm, i, j, k, rank, rank1;

	n = a->s_li();
	m = b->s_li();
	nm = 1;
	for (i = 0; i < m; i++)
		nm *= n;
	c->m_il(nm);
	R.m_il(m);
	R1.m_il(m);
#ifdef PERM_ON_EXPONENTIATION_DEBUG
	printf("perm_on_exponentiation()\n");
	printf("a = ");
	a->println();
	printf("in component %ld\n", a_in_component);
	printf("b = ");
	b->println();
	printf("f_simultaneously = %ld\n", f_simultaneously);
#endif
	for (rank = 0; rank < nm; rank++) {
		rank1 = rank;
		for (i = 0; i < m; i++) {
			k = rank1 % n;
			R.m_ii(i, k);
			rank1 = rank1 / n;
			}
		for (i = 0; i < m; i++) {
			k = R.s_ii(i);
			j = b->s_ii(i) - 1;
			if (f_simultaneously || i == a_in_component)
				k = a->s_ii(k) - 1;
			R1.m_ii(j, k);
			}
		rank1 = 0;
		for (i = 0; i < m; i++) {
			rank1 *= n;
			k = R1.s_ii(m - 1 - i);
			rank1 += k;
			}
#ifdef PERM_ON_EXPONENTIATION_DEBUG
		printf("rank = %ld R=", rank);
		R.print();
		printf(" R1=");
		R1.print();
		printf("rank1 = %ld\n", rank1);
#endif
		c->m_ii(rank, rank1 + 1);
		}
	return OK;
}

#if TEXDOCU
INT vec_generators_wreath_product(VECTOR_OP W, VECTOR_OP G, VECTOR_OP H, INT f_v)
#else
Computes generators for the wreath product $\la G \ra \wr \la H \ra$ into $W$.
#endif
{
	PERMUTATION_OP p;
	PERMUTATION_OB per;
	INT n, m, i, j, lG, lH, nb_gen;

	lG = G->s_li();
	lH = H->s_li();
	if (lG == 0)
		return error("vec_generators_wreath_product() lG == 0");
	if (lH == 0)
		return error("vec_generators_wreath_product() lH == 0");
	p = (PERMUTATION_OP) G->s_i(0);
	n = p->s_li();
	p = (PERMUTATION_OP) H->s_i(0);
	m = p->s_li();
	nb_gen = 0;
	W->m_il(0);
	if (f_v)
		printf("embedding generators of H:\n");
	for (i = 0; i < lH; i++) {
		p = (PERMUTATION_OP) H->s_i(i);
		wreath_embedding(p, n, m, &per);
		if (f_v)
			per.println();
		W->inc();
		per.swap(W->s_i(nb_gen));
		nb_gen++;
		}
	for (j = 0; j < m; j++) {
		if (f_v)
			printf("embedding generators of G into %ld-component:\n", j);
		for (i = 0; i < lG; i++) {
			p = (PERMUTATION_OP) G->s_i(i);
			wreath_embedding_component(p, n, m, j, &per);
			if (f_v)
				per.println();
			W->inc();
			per.swap(W->s_i(nb_gen));
			nb_gen++;
			}
		}
	return OK;
}

#if TEXDOCU
INT Sn_wreath_Sm_generators(INT n, INT m, VECTOR_OP G)
#else
Computes generators for $S_n \wr S_m$ into $G$.
#endif
{
	PERMUTATION_OP p;
	PERMUTATION_OB per;
	VECTOR_OB Sn, Sm;
	INT i, j, nb_gen_Sn, nb_gen_Sm, nb_gen;
	
	printf("Sn_wreath_Sm_generators(): n = %ld m = %ld\n", n, m);
	if (n == 1) {
		return symmetric_generators(G, m);
		}
	if (m == 1) {
		return symmetric_generators(G, n);
		}
	symmetric_generators(&Sn, n);
	symmetric_generators(&Sm, m);
	nb_gen_Sn = Sn.s_li();
	nb_gen_Sm = Sm.s_li();
	nb_gen = 0;
	G->m_il(0);
	printf("embedding generators of Sm:\n");
	for (i = 0; i < nb_gen_Sm; i++) {
		p = (PERMUTATION_OP) Sm.s_i(i);
		wreath_embedding(p, n, m, &per);
		per.println();
		G->inc();
		per.swap(G->s_i(nb_gen));
		nb_gen++;
		}
	for (j = 0; j < m; j++) {
		printf("embedding generators of Sn into %ld-component:\n", j);
		for (i = 0; i < nb_gen_Sn; i++) {
			p = (PERMUTATION_OP) Sn.s_i(i);
			wreath_embedding_component(p, n, m, j, &per);
			per.println();
			G->inc();
			per.swap(G->s_i(nb_gen));
			nb_gen++;
			}
		}
	return OK;
}

#if TEXDOCU
INT wreath_embedding(PERMUTATION_OP g, INT n, INT m, PERMUTATION_OP q)
#else
utility function for computing wreath product generators.
#endif
{
	INT i, j, ii, first, to;
	INT nm;
	
	if (g->s_li() != m)
		return error("wreath_embedding() g->s_li() != m");
	nm = n * m;
	q->m_il(nm);
	for (i = 0; i < m; i++) {
		j = g->s_ii(i) - 1;
		first = i * n;
		to = j * n;
		for (ii = 0; ii < n; ii++) {
			q->m_ii(first + ii, to + ii + 1);
			}
		}
	return OK;
}

#if TEXDOCU
INT wreath_embedding_component(PERMUTATION_OP g, INT n, INT m, INT j, PERMUTATION_OP q)
#else
utility function for computing wreath product generators.
#endif
{
	INT i, i_im, first;
	INT nm;
	
	if (g->s_li() != n)
		return error("wreath_embedding_component() g->s_li() != n");
	nm = n * m;
	q->m_il(nm);
	q->one();
	first = j * n;
	for (i = 0; i < n; i++) {
		i_im = g->s_ii(i) - 1;
		q->m_ii(first + i, first + i_im + 1);
		}
	return OK;
}

#endif /* PERMTRUE */



