/* perm_grp.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef PERMTRUE

#include <DISCRETA/ma.h> /* for PGL_2_p_as_perm() */
#include <DISCRETA/perm.h>
#include <DISCRETA/part.h>
#include <DISCRETA/fga.h>
#include <DISCRETA/lb.h> // stabilize point
#include <DISCRETA/ladder.h> // for construct_Higman_Sims_176
#ifdef SOLVABLE_TRUE
#include <DISCRETA/solvable.h>
#endif
#ifdef SYM_DB_FG
#include <DISCRETA/fg1.h>
#endif

#if TEXDOCU
INT vec_generators_is_trivial_group(VECTOR_OP gen)
#else
TRUE if the generators are all the identity, FALSE otherwise.
#endif
{
	PERMUTATION_OP p;
	INT i, l;
	
	l = gen->s_li();
	for (i = 0; i < l; i++) {
		p = (PERMUTATION_OP) gen->s_i(i);
		if (!p->einsp())
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT is_abelian(VECTOR_OP G)
#else
True if the elements in the vector G generate an abelian group.
#endif
{
	SYM_OP p, q;
	SYM_OB a, b, c, d, e, f;
	INT i, j, l;
	
	l = G->s_li();
	for (i = 0; i < l; i++) {
		p = G->s_i(i);
		for (j = i + 1; j < l; j++) {
			q = G->s_i(j);
			p->mult(q, &a);
			p->invers(&c);
			q->invers(&d);
			a.mult(&c, &e);
			e.mult(&d, &f);
#if 0
			printf("p=");
			p->println();
			printf("q=");
			q->println();
			printf("[p,q]=");
			f.println();
#endif
			if (!f.einsp())
				return FALSE;
			}
		}
	return TRUE;
}


#if TEXDOCU
INT read_file_of_generators(VECTOR_OP G, char *fname)
#else
opens the file fname for reading and read generators via read\_generators().
#endif
{
	FILE *fp;
	
	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("read_file_of_generators(): can't open file '%s' for reading !\n", fname);
		return ERROR;
		}
	read_generators(G, fp);
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT write_file_of_generators(VECTOR_OP G, char *fname)
#else
opens the file fname for writing and writes generators via write\_generators().
#endif
{
	FILE *fp;
	
	fp = fopen(fname, "w");
	if (fp == NULL)
		return ERROR;
	write_generators(G, fp);
	fclose(fp);
	return OK;
}

#ifndef BUFSIZE
#define BUFSIZE 50000
#endif

#if TEXDOCU
INT read_generators(VECTOR_OP G, FILE *fp)
#else
reads the generators of a permutation group from the file fp. 
The format is:
//PRE
nb-of-generators degree
List-of-generators
///PRE
where List-of-generators contains each generator in a row. 
Each generator is given by listing its image vector. 
The elements which are permuted are the numbers $1,\ldots,d$ 
where $d$ is the degree.
So, $S_5$ could be generated this way:
//PRE
2 5 startwith1
2 1 3 4 5
2 3 4 5 1
///PRE
#endif
{
	BYTE buf[BUFSIZE], *pbuf;
	BYTE buf1[BUFSIZE];
	PERMUTATION_OB p;
	INT i, j, erg = OK, nb, deg, a, a1;
	INT f_startwith0 = FALSE;
	INT f_v = FALSE;
	
	if (fp == NULL) {
		printf("read_generators(): fp == NIL !\n");
		return ERROR;
		}
	// fscanf(fp, "%ld", &nb);
	// fscanf(fp, "%ld", &deg);
	fgets(buf, BUFSIZE, fp);
	pbuf = buf;
	s_scan_int(&pbuf, &nb);
	s_scan_int(&pbuf, &deg);
	s_scan_token(&pbuf, buf1);
	if (strcmp(buf1, "startwith0") == 0)
		f_startwith0 = TRUE;
		// sscanf(buf, "%ld %ld", &nb, &deg);
	
	printf("read_generators() nb=%ld deg=%ld\n", nb, deg);
	if (f_startwith0)
		printf("read_generators() entries start with 0\n");
	fflush(stdout);
	erg += G->m_il(nb);
	for (i = 0; i < nb; i++) {
		erg += p.m_il(deg);
		erg += p.one();
		for (j = 0; j < deg; j++) {
			fscanf(fp, "%ld", &a);
			a1 = a;
			if (f_startwith0)
				a++;
			if (a < 1) 
				printf("entry %ld too small\n", a1);
			if (a > deg) 
				printf("entry %ld too big\n", a1);
			if (f_v) {
				printf("a=%ld\n", a);
				fflush(stdout);
				}
			p.m_ii(j, a);
			}
		if (f_v) {
			p.println();
			}
		p.copy((PERMUTATION_OP) G->s_i(i));
		}
	return erg;
}

#if TEXDOCU
INT write_generators(VECTOR_OP G, FILE *fp)
#else
wites an ASCII file of generators in a file ready for 
reading it with read\_generators().
#endif
{
	PERMUTATION_OP p;
	INT i, j, nb, deg, a;
	
	if (fp == NULL) {
		printf("write_generators(): fp == NIL !\n");
		return ERROR;
		}
	nb = G->s_li();
	deg = perm_vec_get_degree(G);
	printf("write_generators: nb=%ld deg=%ld\n", nb, deg);
	fflush(stdout);
	fprintf(fp, "%ld %ld\n", nb, deg);
	for (i = 0; i < nb; i++) {
		p = (PERMUTATION_OP) G->s_i(i);
		for (j = 0; j < deg; j++) {
			a = p->s_ii(j);
			fprintf(fp, "%ld ", a);
			}
		fprintf(fp, "\n");
		}
	return OK;
}

#if TEXDOCU
INT vec_induced_group_on_subset(VECTOR_OP V, VECTOR_OP subset, VECTOR_OP W)
#else
Assume the elements of $V$ stabilize the set \lq subset\rq.
Then $W$ becomes the restriction of these generators to this set.
The entries in subset must be in $\{0,1,\ldots,deg-1\}$.
#endif
{
	PERMUTATION_OP p;
	PERMUTATION_OB q;
	INT i, j, k, a, b, r, l;

	r = V->s_li();
	W->m_il(r);
	l = subset->s_li();
	for (i = 0; i < r; i++) {
		p = (PERMUTATION_OP) V->s_i(i);
		q.m_il(l);
		for (j = 0; j < l; j++) {
			a = subset->s_ii(j);
			b = p->s_ii(a) - 1;
			for (k = 0; k < l; k++) {
				if (subset->s_ii(k) == b)
					break;
				}
			if (k == l)
				return error("vec_induced_group_on_subset(): the set is not invariant under the group");
				q.m_ii(j, k + 1);
			}
		q.swap((PERMUTATION_OP) W->s_i(i));
		}
	return OK;
}

#if TEXDOCU
INT vec_subgroup_of_hol_of_cyclic_group(VECTOR_OP V, INT n, INT i)
#else
Computes generators for the subgroup of $Hol(C_n)$ of index $i$.
$n$ must be a prime.
#endif
{
	INT j, jj, m, m1, g, alpha, alpha1;
	PERMUTATION_OB p;
	
	m = n - 1;
	if (!is_prime(n)) {
		return error("vec_subgroup_of_hol_of_cyclic_group() n must be a prime !");
		}
	if (i == 0) {
		printf("WARNING: vec_subgroup_of_hol_of_cyclic_group(): "
			"index is 0, setting index to 1!\n");
		i = 1;
		}
	m1 = m / i;
	ggt_iipi(m, i, &g);
	if (m1 * i != m) {
		i = g;
		printf("WARNING: vec_subgroup_of_hol_of_cyclic_group(): "
			"index %ld does not divide n - 1 = %ld!\n", i, m);
		printf("setting index to %ld !\n", g);
		}
	m1 = m / g;
	printf("vec_subgroup_of_hol_of_cyclic_group() creating group "
		"$%ld \\unlhd %ld$\n", n, m1);
	alpha = primitive_root(n);
	alpha1 = alpha;
	for (j = 2; j <= g; j++) {
		alpha1 = (alpha1 * alpha) % n;
		}
	printf("primitive root alpha = %ld\n", alpha);
	printf("alpha^{%ld} = %ld mod %ld\n", g, alpha1, n);
	fflush(stdout);
	
	// the generator of the $n$-cycle:
	p.m_il(n);
	for (j = 0; j < n - 1; j++)
		p.m_ii(j, j + 2);
	p.m_ii(n - 1, 1);
	V->m_il(2);
	p.copy((PERMUTATION_OP) V->s_i(0));
	printf("generator of the cycle: ");
	p.println();
	fflush(stdout);
	
	for (j = 0; j < n; j++) {
		jj = j * alpha1 % n;
		p.m_ii(j, jj + 1);
		printf("$%ld \\mapsto %ld$\\\\\n", j, jj);
		fflush(stdout);
		}
	p.copy((PERMUTATION_OP) V->s_i(1));
	printf("2nd generator: ");
	p.println();
	fflush(stdout);
	
	return OK;
}

#if TEXDOCU
INT vec_hol_of_cyclic_group(VECTOR_OP V, INT n)
#else
Computes generators for the holomorph of $C_n$, 
the cyclic group of $n$ elements. The vector V 
becomes the vector of generators.
#endif
{
	VECTOR_OB perm, perminv, red;
	PERMUTATION_OB p;
	INT flag, i, j, ii, a, l = 0;

	perm.m_il_n(n + 1);
	perminv.m_il_n(n + 1);
	red.m_il_n(n + 1);

	p.m_il(n);
	for (i = 0; i < n - 1; i++)
		p.m_ii(i, i + 2);
	p.m_ii(n - 1, 1);
	V->m_il(1);
	l++;
	p.copy((PERMUTATION_OP) V->s_i(l - 1));
	for (a = 2; a < n; a++) {

		// we test the map: i \mapsto a * i mod n
		
		if (red.s_ii(a))
			continue;
		
		for (ii = 0; ii <= n; ii++)
			perminv.m_ii(ii, 0);
		flag = TRUE;
		
		for (i = 1; i < n; i++) {
			j = (i * a) % n;
			if (perminv.s_ii(j)) { // not an automorphism !
				flag = FALSE;
				break;
				}
			perminv.m_ii(j, i);
			perm.m_ii(i, j + 1);
	
			} // next i

		if (flag) { // we found an automorphism !
			i = 1;
			j = 0;
			do { // cycle of 1 under a:
				j++;
				i = (i * a) % n;
				red.m_ii(i, 1);
				} while (i > 1 && j < n);
			p.m_ii(0, 1);
			for (ii = 1; ii < n; ii++)
				p.m_ii(ii, perm.s_ii(ii));
			V->inc();
			l++;
			p.copy((PERMUTATION_OP) V->s_i(l - 1));
			}
		
		} // next a

	V->realloc_z(l);
	
	return OK;
}

#if TEXDOCU
INT vec_conjugate(VECTOR_OP gen, PERMUTATION_OP p)
#else
conjugates a vector of generators with the element p. 
The new vector elements are $\{ p^{-1} g p \mid g \in gen\}$.
#endif
{
	PERMUTATION_OB pv, tmp1;
	PERMUTATION_OP per;
	INT i, l;

	p->invers(&pv);
	l = gen->s_li();
	for (i = 0; i < l; i++) {
		per = (PERMUTATION_OP) gen->s_i(i);
		if (per->s_obj_k() != PERMUTATION)
			return error("vec_conjugate() not a permutation");
		pv.mult(per, &tmp1);
		tmp1.mult(p, per);
		}
	return OK;
}

#if TEXDOCU
INT vec_induce_action_on_blocks(VECTOR_OP gen, VECTOR_OP B)
#else
calls induce\_action\_on\_blocks() for all elements in the vector gen.
the result replaces the original vector gen.
#endif
{
	VECTOR_OB g;
	PERMUTATION_OP p1, p2;
	INT i, l;
	
	l = gen->s_li();
	g.m_il(l);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		p2 = (PERMUTATION_OP) g.s_i(i);
		p1->induce_action_on_blocks(p2, B);
		}
	g.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_induce_action_on_columns(VECTOR_OP gen, MATRIX_OP M)
#else
calls induce\_action\_on\_blocks() for all elements in the vector gen.
the result replaces the original vector gen.
#endif
{
	VECTOR_OB gen2;
	VECTOR_OB B;
	PERMUTATION_OP p1, p2;
	PERMUTATION_OB q, qv, g, g1;
	INT i, l;
	INT f_v = FALSE;
	
	l = gen->s_li();
	gen2.m_il(l);

	M->calc_blocks(&B, &q, &qv, f_v);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		p2 = (PERMUTATION_OP) gen2.s_i(i);
		p1->induce_action_on_blocks(&g, &B);
		q.mult(&g, &g1);
		g1.mult(&qv, p2);
		}
	gen2.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_induce3(VECTOR_OP gen)
#else
calls induce3 for all elements in the vector gen.
#endif
{
	VECTOR_OB g;
	PERMUTATION_OP p1, p2;
	INT i, l;
	
	l = gen->s_li();
	g.m_il(l);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		p2 = (PERMUTATION_OP) g.s_i(i);
		p1->induce3(p2);
		}
	g.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_induce2(VECTOR_OP gen)
#else
calls induce2 for all elements in the vector gen.
#endif
{
	VECTOR_OB g;
	PERMUTATION_OP p1, p2;
	INT i, l, n;
	
	l = gen->s_li();
	g.m_il(l);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		p2 = (PERMUTATION_OP) g.s_i(i);
		if (i == 0)
			n = p1->s_li();
		p1->induce2(p2, n);
		}
	g.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_induce_on_2tuples(VECTOR_OP gen, INT f_injective)
#else
calls induce\_on\_2tuples for all elements in the vector gen.
#endif
{
	VECTOR_OB g;
	PERMUTATION_OP p1, p2;
	INT i, l;
	
	l = gen->s_li();
	g.m_il(l);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		p2 = (PERMUTATION_OP) g.s_i(i);
		p1->induce_on_2tuples(p2, f_injective);
		}
	g.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_add_fixpoint_in_front(VECTOR_OP gen)
#else
calls add\_fixpoint\_in\_front for all elements in the vector gen.
#endif
{
	VECTOR_OB g;
	PERMUTATION_OP p1, p2;
	INT i, l;
	
	l = gen->s_li();
	g.m_il(l);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		p2 = (PERMUTATION_OP) g.s_i(i);
		p1->add_fixpoint_in_front(p2);
		}
	g.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_add_fixpoint_at_end(VECTOR_OP gen)
#else
calls add\_fixpoint\_at\_end for all elements in the vector gen.
#endif
{
	VECTOR_OB g;
	PERMUTATION_OP p1, p2;
	INT i, l, d;
	
	l = gen->s_li();
	g.m_il(l);
	for (i = 0; i < l; i++) {
		p1 = (PERMUTATION_OP) gen->s_i(i);
		d = p1->s_li();
		p2 = (PERMUTATION_OP) g.s_i(i);
		p1->embed_at(p2, d + 1, 0);
		}
	g.swap(gen);
	return OK;
}

#if TEXDOCU
INT vec_generators_stabilize_point(VECTOR_OP a, VECTOR_OP b)
#else
Given a Permutation group $G$ actin on $X = \{1,2,\ldots, n\}$, 
this routine computes generators for the stabilizer (in $G$) of 
the first point: $G_1$.
A Sims chain for $G$ is computed. Then, only generators 
for $G_1$ are read out of this chain 
(by a call to L.stab\_generators(1, b, FALSE)).
#endif
{
	LABRA_OB L;
	INT la, da;
	PERMUTATION_OP pa;
	
	la = a->s_li();
	if (la == 0)
		return error("vec_generators_stabilize_point() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	da = pa->s_li();
	L.Init(da, a, a->s_li(), FALSE, FALSE);
	L.jerrum(FALSE);
	L.stab_generators(1, b, FALSE);
	return OK;
}

#if TEXDOCU
INT vec_generators_stabilize_a_point(VECTOR_OP a, INT pt, VECTOR_OP b)
#else
Similar to vec\_generators\_stabilize\_point(), 
here, generators for $G_{pt}$ are returned.
The sims chain for $G$ is built like before (with respect to the 
base $1,2,\ldots,n$). Then, L.cycle(0, pt); is performed 
to change to the base $pt,1,2,\ldots,pt-1,pt+1,\ldots,n$.
Then, generators for $G_{pt}$ are determined.
We assume: $0 \le pt < deg$
#endif
{
	LABRA_OB L;
	INT la, da;
	PERMUTATION_OP pa;
	
	la = a->s_li();
	if (la == 0)
		return error("vec_generators_stabilize_point() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	da = pa->s_li();
	L.Init(da, a, a->s_li(), FALSE, FALSE);
	L.jerrum(FALSE);
	if (pt != 0)
		L.cycle(0, pt, FALSE, FALSE);
	L.stab_generators(1, b, FALSE);
	return OK;
}

#if TEXDOCU
INT vec_generators_degree(VECTOR_OP a)
#else
Rerurns the degree of the (permutation-)group generated by the elements of a.
#endif
{
	INT l, la;
	PERMUTATION_OP pa;
	
	la = a->s_li();
	if (la == 0)
		return error("vec_generators_degree() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	if (pa->s_obj_k() != PERMUTATION) {
		error("vec_generators_degree() not applicable, not a vector of permutations !\n");
		}
	l = pa->s_li();
	return l;
}

#if TEXDOCU
INT vec_generators_group_order(VECTOR_OP a, SYM_OP go)
#else
Determines the group order of the group generated by the 
elements in the vector a. A Sims-chain is built to 
compute the group order. 
If the vector does not contain PERMUTATION objects, 
nothing is done (and a group order -1 is returned).
#endif
{
	LABRA_OB L;
	INT la, da;
	PERMUTATION_OP pa;
	
	la = a->s_li();
	if (la == 0)
		return error("vec_generators_group_order() la == 0");
	pa = (PERMUTATION_OP) a->s_i(0);
	if (pa->s_obj_k() != PERMUTATION) {
		printf("vec_generators_group_order() not applicable, not a vector of permutations !\n");
		fflush(stdout);
		go->m_i_i(-1);
		return OK;
		}
	da = pa->s_li();
	printf("vec_generators_group_order(): nb_gen=%ld, gens=\n", la);
	a->Print();
	fflush(stdout);
	
	// printf("vec_generators_group_order(): calling Labra.Init()\n"); fflush(stdout);
	L.Init(da, a, a->s_li(), FALSE, FALSE);
	// printf("vec_generators_group_order(): calling Labra.jerrum()\n"); fflush(stdout);
	L.jerrum(FALSE);
	// printf("vec_generators_group_order(): calling Labra.group_order()\n"); fflush(stdout);
	L.group_order(go);
	// printf("vec_generators_group_order(): calling Labra.freeself()\n"); fflush(stdout);
	L.freeself_debug(0);
	// printf("vec_generators_group_order(): finished\n"); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT vec_generators_remove_fixpoint(VECTOR_OP a, VECTOR_OP b, INT i)
#else
calls remove\_fixpoint($g,i$) for each generator $g$ in $a$. 
$b$ contains the new elements.
#endif
{
	PERMUTATION_OP p1, p2;
	INT j, l;
	
	l = a->s_li();
	b->m_il(l);
	for (j = 0; j < l; j++) {
		p1 = (PERMUTATION_OP) a->s_i(j);
		p2 = (PERMUTATION_OP) b->s_i(j);
		p1->remove_fixpoint(p2, i);
		}
	return OK;
}

#if TEXDOCU
INT young2_generators(VECTOR_OP G, INT deg, INT t)
#else
Produces generators for a young group of type $S_{[t,n-t]}$.
Needs $deg \ge 2$. $t=0$ or $t=deg$ are handled correctly.
#endif
{
	PERMUTATION_OB p, q;
	INT nb_gen = 0;
	
	if (deg < 2)
		return error("young2_generators()|deg < 2");
	if (t < 0 || t > deg)
		return error("young2_generators()|t < 0 || t > deg");
	if (t == 0 || t == deg)
		return symmetric_generators(G, deg);
	
	G->m_il(0);
	if (t > 1) {
		p.m_il(deg);
		p.one();
		p.Add2Cycle(1, 2);
		G->inc();
		p.copy((PERMUTATION_OP) G->s_i(nb_gen));
		nb_gen++;

		p.m_il(deg);
		p.one();
		p.AddNCycle(1, t);
		G->inc();
		p.copy((PERMUTATION_OP) G->s_i(nb_gen));
		nb_gen++;
		}
	
	if (t < deg - 1) {
		p.m_il(deg);
		p.one();
		p.Add2Cycle(t + 1, t + 2);
		G->inc();
		p.copy((PERMUTATION_OP) G->s_i(nb_gen));
		nb_gen++;

		p.m_il(deg);
		p.one();
		p.AddNCycle(t + 1, deg - t);
		G->inc();
		p.copy((PERMUTATION_OP) G->s_i(nb_gen));
		nb_gen++;
		}
	
	return OK;
}

#if TEXDOCU
INT trivial_generators(VECTOR_OP G, INT deg)
#else
Generator for the trivial permutation group on \lq deg\rq elements.
#endif
{
	PERMUTATION_OB p;
	
	if (deg < 1)
		return error("trivial_generators()|deg < 1");

	p.m_il(deg);
	p.one();
	G->m_il(1);
	p.swap(G->s_i(0));
	
	return OK;
}

#if TEXDOCU
INT cyclic_generators(VECTOR_OP G, INT deg)
#else
Generator for the cyclic group on \lq deg\rq elements.
#endif
{
	PERMUTATION_OB p;
	INT i;
	
	if (deg < 1)
		return error("cyclic_generators()|deg < 1");

	// the cycle (1, 2, 3 ... deg):
	p.m_il(deg);
	p.one();
	for (i = 0; i < deg - 1; i++)
		p.m_ii(i, i + 2);
	p.m_ii(deg - 1, 1);
	G->m_il(1);
	p.swap(G->s_i(0));
	
	return OK;
}

#if TEXDOCU
INT symmetric_generators(VECTOR_OP G, INT deg)
#else
Generators for the symmetric group on \lq deg\rq elements.
The generating system consists of the cycle $(1,2,\ldots, deg)$ 
and of the element $(1,2)$ if $deg \ge 2$.
$deg < 1$ is forbidden.
#endif
{
	PERMUTATION_OB p;
	INT i;
	
	if (deg < 1)
		return error("symmetric_generators()|deg < 1");

	// the cycle (1, 2, 3 ... deg):
	p.m_il(deg);
	p.one();
	for (i = 0; i < deg - 1; i++)
		p.m_ii(i, i + 2);
	p.m_ii(deg - 1, 1);
	G->m_il(1);
	p.swap(G->s_i(0));
	
	if (deg >= 2) {
		G->realloc_z(2);
		p.m_il(deg);
		p.one();
		p.m_ii(0, 2);
		p.m_ii(1, 1);
		p.swap(G->s_i(1));
		}
	return OK;
}

#if TEXDOCU
INT symmetric_generators_pp(VECTOR_OP G, INT deg)
#else
Generators of prime-power order. Calls symmetric\_generators() 
and vec\_to\_vec\_pp().
#endif
{
	VECTOR_OB G1;
	INT erg = OK;
	
	erg += symmetric_generators(&G1, deg);
	erg += vec_to_vec_pp(&G1, G, 
		1 /* type */, NIL /* data */);
	return erg;
}

#if TEXDOCU
INT symmetric_generators_transpositions(VECTOR_OP G, INT deg)
#else
$deg < 2$ forbidden. 
Computes the generating system consisting of the elementary 
transpositions $(i,i+1)$ for $S_n$.
#endif
{
	PERMUTATION_OB p;
	INT i, erg = OK;
	
	if (deg < 2)
		return error("symmetric_generators_transpositions()|deg < 2");
	G->m_il(0);
	for (i = 0; i < deg - 1; i++) {
		p.m_il(deg);
		p.one();
		p.m_ii(i, i + 2);
		p.m_ii(i + 1, i + 1);
		G->inc();
		p.swap(G->s_i(G->s_li() - 1));
		}
	return erg;
}

#if TEXDOCU
INT alternating_generators(VECTOR_OP G, INT deg)
#else
$deg < 1$ forbidden, $deg = 1,2$ gives trivial generator (id).
Otherwise, a generating system for $A_n$ consisting of all 
3-cycles is returned.
#endif
{
	PERMUTATION_OB p;
	INT i, j, k, l, erg = OK;
	
	if (deg < 1)
		return error("alternating_generators()|deg < 1");
	if (deg <= 2)
		return trivial_generators(G, deg);
	l = 0;
	for (i = 0; i < deg; i++) {
		for (j = i + 1; j < deg; j++) {
			for (k = j + 1; k < deg; k++) {
				l++;
				}
			}
		}
	erg += G->m_il(l);
	l = 0;
	for (i = 0; i < deg; i++) {
		for (j = i + 1; j < deg; j++) {
			for (k = j + 1; k < deg; k++) {
				/* 3 cycle (i j k): */
				erg += p.m_il(deg);
				p.one();
				p.m_ii(i, j + 1); /* i -> j */
				p.m_ii(j, k + 1); /* j -> k */
				p.m_ii(k, i + 1); /* k -> i */
				p.swap(G->s_i(l));
				l++;
				}
			}
		}
	return erg;
}

#if TEXDOCU
INT dihedral_generators(VECTOR_OP G, INT deg)
#else
$deg < 1$ forbidden.
Returns generators for $D_n$, the dihedral group on $n$ elements 
($n = deg$).
#endif
{
	PERMUTATION_OB p, q;
	INT i, deg_2, erg = OK;
	
	if (deg < 1)
		return error("dihedral_generators()|deg < 1");
	if (deg <= 2)
		return trivial_generators(G, deg);
	erg += p.m_il(deg);
	erg += q.m_il(deg);
	erg += p.one();
	erg += q.one();
	for (i = 0; i < deg; i++)
		if (i < deg - 1)
			p.m_ii(i, i + 2);
		else
			p.m_ii(i, 1);
	deg_2 = deg >> 1;
	for (i = 0; i < deg_2; i++) {
		q.m_ii(i, deg - i);
		q.m_ii(deg - i - 1, i + 1);
		}
	erg += G->m_il(2);
	p.swap(G->s_i(0));
	q.swap(G->s_i(1));
	return erg;
}

#if TEXDOCU
INT dihedral_generators_pp(VECTOR_OP G, INT deg)
#else
Calls dihedral\_generators and vec\_to\_vec\_pp.
#endif
{
	VECTOR_OB G1;
	INT erg = OK;
	
	erg += dihedral_generators(&G1, deg);
	erg += vec_to_vec_pp(&G1, G, 
		1 /* type */, NIL /* data */);
	return erg;
}

#if TEXDOCU
INT Higman_Sims_176_generators(VECTOR_OP G)
#else
Returns generators for the sporadic simple group of Higman and Sims on 176 
points. The generators are first tried to read them in from a file. 
If this does not work (for example, because discreta\_home is not set), 
construct\_Higman\_Sims\_176() is called.
#endif
{
	BYTE str[1024];

	if (discreta_home == NIL)
		printf("Higman_Sims_176_generators() discreta_home is not set !\n");
	else {
#if 0
		sprintf(str, "%s/lib/aut_higman_design_generators.dsc", discreta_home);
		if (G->load(str) == OK)
			return OK;
#endif
		sprintf(str, "%s/lib/aut_higman_design_generators.txt", discreta_home);
		if (read_file_of_generators(G, str) == OK)
			return OK;
		}
	return construct_Higman_Sims_176(G);
}

#if TEXDOCU
INT solvable_group_permutation_representations(VECTOR_OP V, INT maxindex, INT ordering)
#else
Constructs permutation representations of V[0] (which must be 
an array of strings, i.e. an abstract presentation of a group) 
of degree $\le$ maxindex (if possible).
Therefore, it determines all subgroups of index $\le$ maxindex.
The trivial representation is always the first.
In general, there will be several permutation representations fround, 
so one needs to apply select afterwards. 
But, sometimes no representation is constructed because the group 
misses a non-trivial subgroup of index $\le$ maxindex.
The FG-library must exist, e.g. SOLVABLE\_TRUE must be defined.
#endif
{
#ifdef SOLVABLE_TRUE
	VECTOR_OB perm_reps;
	INT l;
#endif
	INT f_v = TRUE;
	
	if (f_v) {
		printf("solvable_group_permutation_representations() "
			"V->s_li() = %ld\n", V->s_li());
		fflush(stdout);
		}

#ifdef SOLVABLE_TRUE
#if 0
	lowindex(V, maxindex, &perm_reps, ordering, f_v);
	// where is lowindex ?
#endif
	l = perm_reps.s_li();
	if (f_v) {
		printf("solvable_group_permutation_representations() "
			"found %ld permutation representations\n", l);
		fflush(stdout);
		}
	perm_reps.swap(V);
#endif
	return OK;
}

#if TEXDOCU
INT vector_select_ith(VECTOR_OP V, INT i)
#else
selects the i-the entry of V.
the result will be a vector of length 1  
containing just this object. 
#endif
{
	// VECTOR_OB V1;
	SYM_OB tmp;
	INT f_v = TRUE;

	if (f_v) {
		printf("vector_select_ith() selcting %ld-th vector element\n", i);
		fflush(stdout);
		}
	if (i <= 0 || V->s_li() < i) {
		printf("vector_select_ith() i = %ld, V->s_li() = %ld\n", i, V->s_li());
		fflush(stdout);
		return ERROR;
		}
	if (f_v) {
		printf("the selected element:\n");
		V->s_i(i - 1)->println();
		fflush(stdout);
		}
	V->s_i(i - 1)->swap(&tmp);
	tmp.swap(V);
	
#if 0
	V1.m_il(1);
	V->s_i(i)->swap(V1.s_i(0));
	if (f_v) {
		printf("the selected element:\n");
		V1.s_i(0)->println();
		fflush(stdout);
		}
	V1.swap(V);
#endif
	return OK;
}

static INT nb_sol_groups[] = {
/*   0 */ 0, 1, 1, 1, 2, 1, 2, 1, 5, 2, 
/*  10 */ 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 
/*  20 */ 5, 2, 2, 1, 15, 2, 2, 5, 4, 1, 
/*  30 */ 4, 1, 51, 1, 2, 1, 14, 1, 2, 2, 
/*  40 */ 14, 1, 6, 1, 4, 2, 2, 1, 52, 2, 
/*  50 */ 5, 1, 5, 1, 15, 2, 13, 2, 2, 1, 
/*  60 */ 12, 1, 2, 4, 267, 1, 4, 1, 5, 1, 
/*  70 */ 4, 1, 50, 1, 2, 3, 4, 1, 6, 1, 
/*  80 */ 52, 15, 2, 1, 15, 1, 2, 1, 12, 1, 
/*  90 */ 10, 1, 4, 2, 2, 1, 231, 1, 5, 2, 
/* 100 */ 16, 1, 4, 1, 14, 2, 2, 1, 45, 1, 
/* 110 */ 6, 2, 43, 1, 6, 1, 5, 4, 2, 1, 
/* 120 */ 44, 2, 2, 1, 4, 5, 16, 1 
};

#if TEXDOCU
INT number_of_solvable_groups(INT n)
#else
Returns the number of solvable groups of order $n$. 
$n \le 127$ is required.
#endif
{
	if (n >= 128)
		return error("number_of_solvable_groups() n >= 128");
	return nb_sol_groups[n];
}

#if TEXDOCU
INT solvable_group_from_catalog(VECTOR_OP V, INT n, INT m)
#else
reads the $m$-th group of order $n$ from the solvable group library 
$1 \le m \le N(n)$ where $N(n)$ is the number of solvable groups of order $n$.
The groups is read from DISCRETA\_HOME/lib/solvable\_groups/ASCII
$m == -1$ means: read all solvable groups of order $n$.
the groups are put into the vector V (as FG\_OB objects).
This function needs the FG-library.
#endif
{
	BYTE str[1024];
	BYTE str1[1024];
#ifdef SOLVABLE_TRUE
	FG_OB G;
#endif
#ifdef SYM_DB_FG
	INT mm;
	FG_IO io;
#endif
	
	if (discreta_home == NIL) {
		printf("solvable_group_from_catalog() discreta_home is not set !\n");
		fflush(stdout);
		return ERROR;
		}
#ifdef SOLVABLE_TRUE
	printf("solvable_group_from_catalog() FG library not available !\n");
	fflush(stdout);
	return ERROR;
#endif
#ifndef SYM_DB_FG
	printf("solvable_group_from_catalog() DB_FG library not available !\n");
	fflush(stdout);
	return ERROR;
#endif
	sprintf(str, "%s/lib/solvable_groups/ASCII/g_%ld_", discreta_home, n);
		// fg.txt added within fg_io_open_channel_r_ascii() !
	sprintf(str1, "%s/lib/solvable_groups/ASCII/g_%ld_fg.txt", discreta_home, n);
	if (file_size(str1) <= 0) {
		printf("solvable_group_from_catalog() file %s does not exist !\n", str);
		fflush(stdout);
		return ERROR;
		}
#ifdef SYM_DB_FG
	fg_io_open_channel_r_ascii(&io.from, str);
	
	// V->m_il(0);
	mm = 0;
	while (G.read_ascii(io.from.fp)) {
		mm++;
		if (m < 0 || (mm == m)) {
			G.generators_and_relations(V, TRUE);
			// V->inc();
			// G.swap((FG_OP) V->s_i(V->s_li() - 1));
			}
		}
	fg_io_close_channel(&io.from);
#endif
}

#if TEXDOCU
INT Mathieu_generators(VECTOR_OP G, INT n)
#else
Returns generators for $M_n$ ($n \in \{11,12,23,24\}$).
The generators are taken from Hall~\cite{Hall59}.
#endif
{
	INT erg = OK;
	
	if (n == 11)
		erg += M11_generators(G);
	else if (n == 12)
		erg += M12_generators(G);
	else if (n == 23)
		erg += M23_generators(G);
	else if (n == 24)
		erg += M24_generators(G);
	else
		return error("wrong Mathieu group.");
	return erg;
}

static INT M11_gen(VECTOR_OP G, INT deg);

#if TEXDOCU
INT M11_generators(VECTOR_OP G)
#else
Returns generators for $M_{11}$.
$M_{11}$ is a
4-ply transitive group of order $11 \cdot 10 \cdot 9 \cdot 8 = 7920$, 
it is equal to the stabilizer of 12 in $M_{12}$.
#endif
{
	return M11_gen(G, 11);
}

static INT M11_gen(VECTOR_OP G, INT deg)
{
	PERMUTATION_OB u, a, b, x, y;
	INT erg = OK;
	
	if (deg < 11)
		return error("M11_gen(): deg < 11");
	erg += u.m_il(deg);
	erg += a.m_il(deg);
	erg += b.m_il(deg);
	erg += x.m_il(deg);
	erg += y.m_il(deg);
	erg += u.one();
	erg += a.one();
	erg += b.one();
	erg += x.one();
	erg += y.one();
	/* u := (1 2 3)(4 5 6)(7 8 9) */
	u.AddNCycle(1, 3);
	u.AddNCycle(4, 3);
	u.AddNCycle(7, 3);
	/* a := (2 4 3 7)(5 6 9 8) */
	a.Add4Cycle(2, 4, 3, 7);
	a.Add4Cycle(5, 6, 9, 8);
	/* b := (2 5 3 9)(4 8 7 6) */
	b.Add4Cycle(2, 5, 3, 9);
	b.Add4Cycle(4, 8, 7, 6);
	/* x := (1 10)(4 5)(6 8)(7 9) */
	x.Add2Cycle(1, 10);
	x.Add2Cycle(4, 5);
	x.Add2Cycle(6, 8);
	x.Add2Cycle(7, 9);
	/* y := (1 11)(4 6)(5 9)(7 8) */
	y.Add2Cycle(1, 11);
	y.Add2Cycle(4, 6);
	y.Add2Cycle(5, 9);
	y.Add2Cycle(7, 8);

	erg += G->m_il(5);
	u.swap(G->s_i(0));
	a.swap(G->s_i(1));
	b.swap(G->s_i(2));
	x.swap(G->s_i(3));
	y.swap(G->s_i(4));
	return erg;
}

#if TEXDOCU
INT M12_generators(VECTOR_OP G)
#else
Returns generators for $M_{12}$.
$M_{12}$ is of order 95040.
#endif
{
	PERMUTATION_OB z;
	INT erg = OK;
	
	erg += M11_gen(G, 12);
	
	erg += z.m_il(12);
	erg += z.one();
	/* z := (1 12)(4 7)(5 6)(8 9) */
	z.Add2Cycle(1, 12);
	z.Add2Cycle(4, 7);
	z.Add2Cycle(5, 6);
	z.Add2Cycle(8, 9);
	G->inc();
	z.swap(G->s_i(G->s_li() - 1));
	return erg;
}

static INT M23_gen(VECTOR_OP G, INT deg);

#if TEXDOCU
INT M23_generators(VECTOR_OP G)
#else
Returns generators for $M_{23}$.
$M_{23}$ is 4-ply transitive of order $23 \cdot 22 \cdot 21 \cdot 20 \cdot 16 \cdot 3$.
It is the stabilizer of 24 in $M_{24}$.
#endif
{
	return M23_gen(G, 23);
}

static INT M23_gen(VECTOR_OP G, INT deg)
{
	PERMUTATION_OB a, b;
	INT erg = OK;
	
	if (deg < 23)
		return error("M23_gen(): deg < 23");
	erg += a.m_il(deg);
	erg += b.m_il(deg);
	erg += a.one();
	erg += b.one();
	/* a = (1 2 3 4 5 6 7 8 9 10 .. 22 23) */
	a.AddNCycle(1, 23);
	/* b = (3 17 10 7 9)(4 13 14 19 5)
	 * (8 18 11 12 23)(15 20 22 21 16) */
	b.Add5Cycle(3, 17, 10, 7, 9);
	b.Add5Cycle(4, 13, 14, 19, 5);
	b.Add5Cycle(8, 18, 11, 12, 23);
	b.Add5Cycle(15, 20, 22, 21, 16);

	erg += G->m_il(2);
	a.swap(G->s_i(0));
	b.swap(G->s_i(1));
	return erg;
}

#if TEXDOCU
INT M24_generators(VECTOR_OP G)
#else
Returns generators for $M_{24}$.
#endif
{
	PERMUTATION_OB c;
	INT erg = OK;
	
	erg += M23_gen(G, 24);
	
	erg += c.m_il(24);
	erg += c.one();
	/* c := (1 24)(2 23)(3 12)(4 16)
	 * (5 18)(6 10)(7 20)(8 14)(9 21)
	 * (11 17)(13 22)(15 19) */
	c.Add2Cycle(1, 24);
	c.Add2Cycle(2, 23);
	c.Add2Cycle(3, 12);
	c.Add2Cycle(4, 16);
	c.Add2Cycle(5, 18);
	c.Add2Cycle(6, 10);
	c.Add2Cycle(7, 20);
	c.Add2Cycle(8, 14);
	c.Add2Cycle(9, 21);
	c.Add2Cycle(11, 17);
	c.Add2Cycle(13, 22);
	c.Add2Cycle(15, 19);
	
	G->inc();
	c.swap(G->s_i(G->s_li() - 1));
	return erg;
}

#if TEXDOCU
INT cube_generators(VECTOR_OP G)
#else
Returns generators for the automorphism group of the cube 
$\simeq S_2 \wr S_3 \simeq S_4$ on 8 points.
#endif
{
	PERMUTATION_OB a, b, c;
	INT erg = OK, deg = 8;
	
	erg += a.m_il(deg);
	erg += b.m_il(deg);
	erg += c.m_il(deg);
	erg += a.one();
	erg += b.one();
	erg += c.one();
	/* a := (1)(2 3 4)(5)(6 7 8) */
	a.Add3Cycle(2, 3, 4);
	a.Add3Cycle(6, 7, 8);
	/* b := (1)(2 3)(4)(5)(6 7)(8) */
	b.Add2Cycle(2, 3);
	b.Add2Cycle(6, 7);
	/* c := (1 6)(2 5)(3 8)(4 7) */
	c.Add2Cycle(1, 6);
	c.Add2Cycle(2, 5);
	c.Add2Cycle(3, 8);
	c.Add2Cycle(4, 7);

	erg += G->m_il(3);
	a.swap(G->s_i(0));
	b.swap(G->s_i(1));
	c.swap(G->s_i(2));
	return erg;
}

#if TEXDOCU
INT Zn_mult_generators(INT n, INT mu, VECTOR_OP G)
#endif
{
	PERMUTATION_OB p;
	
	if (n < 2)
		return error("Zn_mult_generators()|n < 2");
	p.cyclic_multiplicator(n, mu);
	G->m_il(1);
	p.swap(G->s_i(0));
	return OK;
}

#endif /* PERMTRUE */


