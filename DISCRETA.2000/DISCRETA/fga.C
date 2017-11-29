/* fga.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef FGA_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h> // for report_group_latex_stdout
#include <DISCRETA/lb.h> // for report_group_latex_stdout
#include <DISCRETA/fga.h>
#include <DISCRETA/gfq.h> /* for PSU_3(q^2) and Sz(q) */
#include <DISCRETA/DESIGNS/solid.h>

INT orbits(VECTOR_OP G, 
	VECTOR_OP SVorbit, VECTOR_OP SVlast, 
	VECTOR_OP SVgen, 
	VECTOR_OP Ofirst, VECTOR_OP Osize)
/* G generators of a permutation group
 * set_size = deg G
 * SVlast: vector of length set_size over integer
 * SVgen: vector of length set_size over integer
 * SVorbit: vector of length set_size over integer
 *          gives for every element 
 *          of the set its orbit number
 * Ofirst: vector over integer
 *         chosen orbit representatives
 * Osize: vector over integer
 *        the number of elements of each orbit
 * the vectors SVlast and Ofirst count the set 
 * from 1 to set_size. 
 */
{
	INT set_size, *Q = NIL, Q_len;
	INT i, j, k, cur, next, os;
	PERMUTATION_OP g;
	INT erg = OK;
	
	if (G->s_li() == 0)
		return error("orbits(): G->s_li == 0");
	g = (PERMUTATION_OP) G->s_i(0);
	set_size = g->s_li();
	Q = (INT *) my_malloc(set_size * sizeof(INT), "fga.C: orbits()");
	if (Q == NIL)
		return error("Q == NIL");
	erg += SVorbit->m_il(set_size);
	erg += SVlast->m_il(set_size);
	erg += SVgen->m_il(set_size);
	for (i = 1; i <= set_size; i++)
		SVorbit->m_ii(i - 1, -1);
	Ofirst->m_il(0);
	Osize->m_il(0);
	for (i = 1; i <= set_size; i++) {
		if (SVorbit->s_ii(i - 1) != -1)
			continue;
		SVorbit->m_ii(i - 1, Ofirst->s_li());
		SVlast->m_ii(i - 1, -1);
		SVgen->m_ii(i - 1, -1);
		os = 1;
		Q[0] = i;
		Q_len = 1;
		while (Q_len) {
			cur = Q[0];
			for (k = 1; k < Q_len; k++)
				Q[k - 1] = Q[k];
			Q_len--;
			for (j = 0; j < G->s_li(); j++) {
				g = (PERMUTATION_OP) G->s_i(j);
				next = g->s_ii(cur - 1);
				if (SVorbit->s_ii(next - 1) != -1)
					continue;
				os++;
				SVlast->m_ii(
					next - 1, cur);
				SVgen->m_ii(
					next - 1, j);
				SVorbit->m_ii(
					next - 1, Ofirst->s_li());
				Q[Q_len++] = next;
				}
			}
		Ofirst->inc();
		Osize->inc();
		Ofirst->m_ii(Ofirst->s_li() - 1, i);
		Osize->m_ii(Osize->s_li() - 1, os);
		}
	
	if (Q) {
		my_free(Q);
		Q = NIL;
		}
	return erg;
}

INT trace_schreier_vectors(INT i, 
	VECTOR_OP SVlast, VECTOR_OP SVgen, 
	VECTOR_OP G, SYM_OP p)
/* SVlast, SVgen: vectors over INT; 
 * G: vector over group elements 
 * (multiplicatively written)
 * p: will become the element 
 * leading form orbit_first to i 
 * (a product in the generators).*/
{
	SYM_OP g_k;
	SYM_OB q;
	INT erg = OK;
	INT last, j, k;
	
	G->s_i(0)->copy(p);
	p->one();
	j = i;
	while ((last = SVlast->s_ii(j - 1)) != -1) {
		/* there is a previous element, last;
		 * follow the path leading to orbit_first, 
		 * put the generators found in front of p: */
		k = SVgen->s_ii(j - 1);
			/* this is the generator 
			 * leading from last to j */
		g_k = G->s_i(k);
		erg += g_k->mult(p, &q);
		q.swap(p);
		j = last;
		}
	return erg;
}

/* schreier vectors: 
 * SVlast[i - 1] = j, SVgen[i - 1] = k iff g_k(j) = i 
 * SVlast[i - 1] = -1 iff i is a orbit 
 *       representative (orbit_first)
 * where G = { g_0, g_1, ... } are the generators 
 * of the operating group. */

INT schreier_stabilizer(VECTOR_OP G, 
	VECTOR_OP SVorbit, VECTOR_OP SVlast, 
	VECTOR_OP SVgen, 
	INT which_orbit, VECTOR_OP stab)
/* computes the stabilizer of the 
 * chosen orbit representative 
 * of orbit number which_orbit. */
{
	PERMUTATION_OP g;
	SYM_OB p0, p1, q0, q1, q2;
	INT erg = OK;
	INT set_size;
	INT nb_gen;
	INT i, j, k, idx, f_found;
	
	g = (PERMUTATION_OP) G->s_i(0);
	G->s_i(0)->copy(&p0);
	p0.one();
	erg += stab->m_il(1);
	p0.copy(stab->s_i(0));
	set_size = g->s_li();
	nb_gen = G->s_li();
	for (i = 1; i <= set_size; i++) {
		if (SVorbit->s_ii(i - 1) != which_orbit)
			continue;
		erg += trace_schreier_vectors(i, 
			SVlast, SVgen, G, &p0);
		for (j = 0; j < nb_gen; j++) {
			erg += G->s_i(j)->mult(&p0, &p1);
			g = (PERMUTATION_OP) G->s_i(j);
			k = g->s_ii(i - 1);
			erg += trace_schreier_vectors(k, 
				SVlast, SVgen, G, &q0);
			erg += q0.invers(&q1);
			erg += p1.mult(&q1, &q2);
			if (q2.einsp())
				continue;
			erg += stab->search(stab->s_li(), 
				TRUE, &q2, &idx, &f_found);
			if (f_found)
				continue;
			stab->inc();
			for (k = stab->s_li() - 1; k > idx; k--)
				stab->s_i(k)->swap(stab->s_i(k - 1));
			q2.copy(stab->s_i(idx));
			}
		}
	return erg;
}

INT transitivity(VECTOR_OP G, INT *t, 
	INT *f_sharp, INT *f_point_five, INT f_verbose)
{
	PERMUTATION_OP g;
	VECTOR_OB G0, G1, SVlast, 
		SVgen, SVorbit, Ofirst, Osize;
	INT erg = OK, t1, degree, i, os;
	
	g = (PERMUTATION_OP) G->s_i(0);
	degree = g->s_li();
	erg += G->copy(&G0);
	t1 = 0;
	/* number of points that 
	 * are currently stabilized */
	*f_sharp = FALSE;
	*f_point_five = FALSE;
	while (TRUE) {
		erg += orbits(&G0, &SVorbit, &SVlast, &SVgen, 
			&Ofirst, &Osize);
		if (f_verbose) {
			printf("SVorbit: "); SVorbit.println();
			printf("SVlast: "); SVlast.println();
			printf("SVgen: "); SVgen.println();
			printf("Ofirst: "); Ofirst.println();
			printf("Osize: "); Osize.println();
			fflush(stdout);
			}
		if (Ofirst.s_li() > t1 + 1) {
			os = Osize.s_ii(t1);
			for (i = t1 + 1; i < Ofirst.s_li(); i++) {
				if (Osize.s_ii(i) != os)
					break;
				}
			if (i == Ofirst.s_li())
				*f_point_five = TRUE;
			break;
			}
		if (Ofirst.s_li() == degree) {
			break;
			}
		erg += schreier_stabilizer(&G0, 
			&SVorbit, &SVlast, &SVgen, 
			t1 /* which_orbit */, &G1);
		if (f_verbose) {
			G1.println();
			fflush(stdout);
			}
		t1++;
		if (G1.s_li() == 1) { /* always one entry id */
			*f_sharp = TRUE;
			break;
			}
		erg += G1.copy(&G0);
		}
	if (t1 == degree - 1)
		t1++;
	*t = t1;
	return erg;
}

#if TEXDOCU
INT report_group_latex_stdout(VECTOR_OP G_gen)
#else
#endif
{
	BYTE str[1024];
	INT i, l;
	SYM_OB ol, goG;
	LABRA_OB labG;
	MATRIX_OB TG;
	
	l = G_gen->s_li();
	printf("\\begin{eqnarray*}\n");
	for (i = 0; i < l; i++) {
		if (i == 0)
			printf("G & = & \\langle ");
		else
			printf(" & & ");
		G_gen->s_i(i)->latex(stdout);
		if ( i < l - 1)
			printf(", ");
		else
			printf("\\rangle");
		printf("\\\\\n");
		}
	printf("\\end{eqnarray*}\n");	
	
	reduce_generators_labra(G_gen, &goG, FALSE /* f_verbose */, &labG);
	str[0] = 0;
	goG.sprint(str);
	printf("$|G| = %s$\\\\\n", str);
	
	printf("{Sims}-chain:\n");
	printf("\\[\n");
	
	labG.calc_transversal_matrix(&TG, FALSE /* f_v */);
	// TG.Print();
	TG.latex_upper_tri(stdout);
	
	printf("\\]\n");

	printf("the transversal elements:\\\\\n");
	labG.latex_transversal_indices(stdout);
	return OK;
}

#if TEXDOCU
INT perm_vec_get_degree(VECTOR_OP V)
#endif
{
	PERMUTATION_OP p;
	INT l, deg;

	l = V->s_li();
	if (l < 1)
		return error("perm_vec_get_degree(): l < 1");
	p = (PERMUTATION_OP) V->s_i(0);
	deg = p->s_li();
	return deg;
}

#if TEXDOCU
void go(VECTOR_OP p)
#endif
{
	VECTOR_OB q;
	SYM_OB go;
	
	printf("go(): determining group order "
		"(number of generators: %ld)...\n", p->s_li());
	p->copy(&q);
	gen_reduce(&q, &go, TRUE, stdout);
	printf("go = ");
	go.println();
}

#if TEXDOCU
INT gen_reduce(VECTOR_OP gen, SYM_OP group_order, INT f_verbose, FILE *fp)
#endif
{
	INT i;
	BYTE str[256], str1[1024];
	LABRA_OB L;
	
	reduce_generators_labra(
		gen, group_order, 
		FALSE /* f_verbose */, &L);
	if (f_verbose) {
		printf("number of generators: %ld\n", 
			gen->s_li());
		fflush(stdout);
		}
	if (fp) {
		fprintf(fp, "%% generators:\n");
		for (i = 0; i < gen->s_li(); i++) {
			str1[0] = 0;
			gen->s_i(i)->sprint(str1);
			fprintf(fp, "%% %s\n", str1);
			}
		}
	str[0] = 0;
	group_order->sprint(str);
	sprintf(str1, "%% group order: %s", str);
	if (fp)
		fprintf(fp, "%s\n", str1); fflush(fp);
	printf("%s\n", str1); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT reduce_generators(VECTOR_OP p, SYM_OP group_order, INT f_verbose)
#endif
{
	VECTOR_OB p1, G;
	INT erg = OK;
	
	p->swap(&p1);
	erg += group_elements(&p1, &G);
	if (group_order) {
		if (G.s_li() == 0) 
			((INTEGER_OP) group_order)->m_i(1);
		else
			((INTEGER_OP) group_order)->m_i(G.s_li());
		}
	p1.swap(p);
	return erg;
}

#undef TEST_REDUCE_GENERATORS_LABRA

static INT reduce_generators_labra_(VECTOR_OP p, SYM_OP group_order, 
	INT f_verbose, LABRA_OP L, INT f_test_labra);

#if TEXDOCU
INT reduce_generators_labra(VECTOR_OP p, SYM_OP group_order, 
	INT f_verbose, LABRA_OP L)
#endif
{
	INT f_test_labra = FALSE;
#ifdef TEST_REDUCE_GENERATORS_LABRA
	f_test_labra = TRUE;
#endif
	return reduce_generators_labra_(p, group_order, f_verbose, L, f_test_labra);
}

#if TEXDOCU
static INT reduce_generators_labra_(VECTOR_OP p, SYM_OP go, 
	INT f_verbose, LABRA_OP L, INT f_test_labra)
#endif
{
	// LABRA_OB L;
	PERMUTATION_OP perm;
	VECTOR_OB p1, G;
	INT erg = OK;
	INT deg;
	
	if (p->s_li() == 0) {
		((INTEGER_OP) go)->m_i(1);
		return OK;
		}
	if (vec_generators_is_trivial_group(p)) {
		((INTEGER_OP) go)->m_i(1);
		return OK;
		}
	p->copy(&p1);
	perm = (PERMUTATION_OP) p1.s_i(0);
	deg = perm->s_li();
	erg += L->Init(deg, 
		&p1, p1.s_li(), 
		FALSE /* f_verbose */, 
		FALSE /* f_very_verbose */);
	L->jerrum(FALSE /* f_verbose */);
	L->group_order(go);
	if (f_verbose) {
		L->Print();
		L->print_T();
		}
#if 0
	L->reduced_generating_set(&p1, TRUE /* f_bottom_up */, FALSE /* f_v */);
	L->group_order(go);
#endif
	L->generators(&p1);
	if (p1.s_li() > p->s_li())
		return erg;
	p1.swap(p);
	return erg;
}

#if TEXDOCU
INT group_elements(VECTOR_OP G_gen, VECTOR_OP G)
#else
Computes all elements of the group generated by G\_gen. 
The list of elements is returned in G.
This function uses the algorithm of Dimino (cf. Butler~\cite{Butler91b}).
#endif
{
	VECTOR_OB Ggen;
	INT i;
	
	Ggen.m_il(0);
	G->m_il(0);
	for (i = 0; i < G_gen->s_li(); i++) {
		dimino_ex(G, &Ggen, G_gen->s_i(i));
		}
	Ggen.swap(G_gen);
	return OK;
}

#if TEXDOCU
INT dimino_ex(VECTOR_OP G, VECTOR_OP G_gen, SYM_OP g)
#endif
{
	VECTOR_OB G_tmp, reps;
	SYM_OP g1;
	SYM_OB id, h, h1;
	INT erg = OK;
	INT j, k, i1, l1, idx, f_found;
	
	g->copy(&id);
	id.one();
	if (G->emptyp() || G->s_li() == 0) {
		G->m_il(1);
		id.copy(G->s_i(0));
		}
	erg += G->search(G->s_li(), 
		TRUE /* f_ascending */, 
		g, &idx, &f_found);
	if (f_found)
		return OK;
	if (G_gen->emptyp())
		G_gen->m_il(1);
	else
		G_gen->inc();
	g->copy(G_gen->s_i(G_gen->s_li() - 1));
	erg += G->copy(&G_tmp);
	reps.m_il(1);
	erg += id.copy(reps.s_i(0));
	/* the dimino induction step: 
	 * apply all generators 
	 * to all coset reps. */
	for (j = 0; j < reps.s_li(); j++) {
		for (i1 = 0; 
			i1 < G_gen->s_li(); i1++) {
			if (j == 0 && 
				i1 < G_gen->s_li() - 1)
				continue;
			g1 = G_gen->s_i(i1);
			g1->mult(reps.s_i(j), &h);
			l1 = G->s_li();
			erg += G->search(l1, 
				TRUE /* f_ascending */, 
				&h, &idx, &f_found);
			if (f_found)
				continue;
			reps.inc();
			erg += h.copy(
				reps.s_i(reps.s_li() - 1));
			erg += G->v_realloc(
				l1 + G_tmp.s_li());
			for (k = 0; k < G_tmp.s_li(); k++) {
				h.mult(G_tmp.s_i(k), &h1);
				/* G hat jetzt l1 + k 
				 * benuetzte Elemente */
				erg += G->search_and_insert(
					l1 + k, &h1);
				} /* next k */
			} /* next i1 */
		} /* next j */
	return OK;
}

/*
---------------------------------------------------------------------------------------------
*/


#define DEBUG_GET_GROUP

#if TEXDOCU
INT km_get_group_generators1(INT type, INT val1, INT val2, BYTE *s, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len)
#else
Using type/val1/val2 as a group description, this function 
generates solvable, sporadic and well-known groups if one of these is requested. 
The result (a vector of generators) is put on top of S, 
where S is a vector which is used as a stack.
Returns TRUE if something has been done, FALSE otherwise.
In the following table $n$ is the parameter val1, $m$ is val2.
\begin{center}
\begin{tabular}{ll}
type & group \\
\hline
FGA\_GROUP\_MATHIEU & $M_{n}$, $n \in \{11,12,23,24\}$ \\
FGA\_GROUP\_HIGMAN\_SIMS\_176 & $\HS_{176}$: Higman-Sims group on 176 points \\
FGA\_SOLVABLE\_GROUP & solvable group $m\#n$ from catalogue ($m \le 127$) \\
FGA\_GROUP\_TRIVIAL & $Id_{n}:$ trivial group on $n$ points \\
FGA\_GROUP\_SYM & $S_n$ generated by $(1,2,\ldots,n)$ and $(1,2)$ \\
FGA\_GROUP\_SYM\_TRANSPOSITIONS & $S_n$ generated by $(i,i+1)$, $1 \le i \le n - 1$ \\
FGA\_GROUP\_ALT & $A_n$ \\
FGA\_GROUP\_DIHEDRAL & $D_n$ (of order $2n$) \\
FGA\_GROUP\_CYCLIC & $C_n$ \\
FGA\_GROUP\_HOLOMORPH\_\\
$\;$OF\_CYCLIC\_GROUP & $\Hol(C_n)$ \\
FGA\_SUBGROUP\_OF\_HOLOMORPH\_\\
$\;$OF\_CYCLIC\_GROUP & subgroup of $\Hol(C_n)$ of index $m$\\
FGA\_GROUP\_SN\_WREATH\_SM & $S_n \wr S_m$ \\
FGA\_GROUP\_BY\_PERMUTATION\_GENERATOR & $\la p \ra$, $p$ is described by a string in s\\
\end{tabular}
\end{center}
#endif
{
	VECTOR_OB gen;
	
	if (type == FGA_GROUP_MATHIEU) {
		sprintf(g_label + strlen(g_label), "M%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf M}_{%ld}", val1);
		Mathieu_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, MATHIEU %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_HIGMAN_SIMS_176) {
		sprintf(g_label + strlen(g_label), "HS_176");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf HS}_{176}");
		Higman_Sims_176_generators(&gen);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, HIGMAN_SIMS_176: %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLVABLE_GROUP) {
		sprintf(g_label + strlen(g_label), "sg%ld.%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf SG}{%ld\\#%ld}", val1, val2);
		if (solvable_group_from_catalog(&gen, val1, val2) != OK)
			return ERROR;
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SOLVABLE GROUP %ld#%ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_TRIVIAL) {
		sprintf(g_label + strlen(g_label), "Id%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Id}_{%ld}", val1);
		trivial_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, trivial group on %ld points: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_SYM) {
		sprintf(g_label + strlen(g_label), "S%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf S}_{%ld}", val1);
		symmetric_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SYM %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_SYM_TRANSPOSITIONS) {
		sprintf(g_label + strlen(g_label), "S%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf S}_{%ld}", val1);
		symmetric_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SYM (generated by transpositions) %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_ALT) {
		sprintf(g_label + strlen(g_label), "A%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf A}_{%ld}", val1);
		alternating_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, ALT %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_DIHEDRAL) {
		sprintf(g_label + strlen(g_label), "D%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf D}_{%ld}", val1);
		dihedral_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, DIHEDRAL %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_CYCLIC) {
		sprintf(g_label + strlen(g_label), "C%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf C}_{%ld}", val1);
		cyclic_generators(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, CYCLIC %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_HOLOMORPH_OF_CYCLIC_GROUP) {
		sprintf(g_label + strlen(g_label), "C%ldHol", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm Hol}({\\bf C}_{%ld})", val1);
		vec_hol_of_cyclic_group(&gen, val1);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, HOLOMORPH_OF_CYCLIC %ld: %s %s \n", 
			*S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SUBGROUP_OF_HOLOMORPH_OF_CYCLIC_GROUP) {
		sprintf(g_label + strlen(g_label), "C%ldHolindex%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm Hol}({\\bf C}_{%ld})index%ld", val1, val2);
		vec_subgroup_of_hol_of_cyclic_group(&gen, val1, val2);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SUBGROUP_OF_HOLOMORPH_OF_CYCLIC %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_SN_WREATH_SM) {
		sprintf(g_label + strlen(g_label), "S%ldwrS%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf S}_{%ld}\\wr{\\bf S}_{%ld}", val1, val2);
		Sn_wreath_Sm_generators(val1, val2, &gen);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SN_WREATH_SM %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_BY_PERMUTATION_GENERATOR) {
		BYTE *s1 = s;
		PERMUTATION_OB perm;
		
		sprintf(g_label + strlen(g_label), "perm");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\langle %s \\rangle}", s);
		perm.sscan(&s1, TRUE /* f_v */);
		gen.m_il(1);
		perm.swap(gen.s_i(0));
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, PERMUTATION_GENERATOR %s: %s %s \n", 
			*S_len, s, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else
		return FALSE;
	return TRUE;
}

#if TEXDOCU
INT km_get_group_generators2(INT type, INT val1, INT val2, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len)
#else
Using type/val1/val2 as a group description, this function 
generates linear groups if one of these is requested. 
The result (a vector of generators) is put on top of S, 
where S is a vector which is used as a stack.
Returns TRUE if something has been done, FALSE otherwise.
The parameter val1 is used for the dimension ($=n$), 
val2 is the size of the field ($=q$).
\begin{center}
\begin{tabular}{ll}
type & group \\
\hline
FGA\_GROUP\_PSL & $\PSL_n(q)$ \\
FGA\_GROUP\_PGL & $\PGL_n(q)$ \\
FGA\_GROUP\_PSSL & $\PSSL_n(q)$ \\ 
FGA\_GROUP\_PGGL & $\PGGL_n(q)$ \\
FGA\_GROUP\_SL & $\SL_n(q)$ \\
FGA\_GROUP\_GL & $\GL_n(q)$ \\
FGA\_GROUP\_SSL & $\SSL_n(q)$ \\
FGA\_GROUP\_GGL & $\GGL_n(q)$ \\
FGA\_GROUP\_ASL & $\ASL_n(q)$ \\
FGA\_GROUP\_AGL & $\AGL_n(q)$ \\
FGA\_GROUP\_ASSL & $\ASSL_n(q)$ \\
FGA\_GROUP\_AGGL & $\AGGL_n(q)$ \\
FGA\_GROUP\_AFFINE\_TRANSLATIONS & \\
FGA\_GROUP\_PSU\_3\_Q2 & \\
\end{tabular}
\end{center}
#endif
{
	VECTOR_OB gen;
	
	if (type == FGA_GROUP_PSL) {
		sprintf(g_label + strlen(g_label), "PSL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm PSL}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			FALSE /* f_affine */ , 
			FALSE /* f_semi */, 
			TRUE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, TRUE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, PSL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_PGL) {
		sprintf(g_label + strlen(g_label), "PGL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm PGL}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			FALSE /* f_affine */ , 
			FALSE /* f_semi */, 
			FALSE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, TRUE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, PGL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_PSSL) {
		sprintf(g_label + strlen(g_label), "PSSL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm P}\\Sigma {\\rm L}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			FALSE /* f_affine */ , 
			TRUE /* f_semi */, 
			TRUE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, TRUE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, PSSL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_PGGL) {
		sprintf(g_label + strlen(g_label), "PGGL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm P}\\Gamma {\\rm L}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			FALSE /* f_affine */ , 
			TRUE /* f_semi */, 
			FALSE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, TRUE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, PGGL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_SL) {
		sprintf(g_label + strlen(g_label), "SL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm SL}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			FALSE /* f_semi */, 
			TRUE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_GL) {
		sprintf(g_label + strlen(g_label), "GL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm GL}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			FALSE /* f_semi */, 
			FALSE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, GL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_SSL) {
		sprintf(g_label + strlen(g_label), "SSL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "\\Sigma {\\rm L}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			TRUE /* f_semi */, 
			TRUE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SSL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_GGL) {
		sprintf(g_label + strlen(g_label), "GGL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "\\Gamma {\\rm L}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			TRUE /* f_semi */, 
			FALSE /* f_special */ , 
			FALSE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, GGL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_ASL) {
		sprintf(g_label + strlen(g_label), "ASL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm ASL}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			FALSE /* f_semi */, 
			TRUE /* f_special */ , 
			TRUE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, ASL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_AGL) {
		sprintf(g_label + strlen(g_label), "AGL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm AGL}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			FALSE /* f_semi */, 
			FALSE /* f_special */ , 
			TRUE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, AGL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_ASSL) {
		sprintf(g_label + strlen(g_label), "ASSL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm A}\\Sigma {\\rm L}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			TRUE /* f_semi */, 
			TRUE /* f_special */ , 
			TRUE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, ASSL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_AGGL) {
		sprintf(g_label + strlen(g_label), "AGGL_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm A}\\Gamma {\\rm L}_{%ld}(%ld)", val1, val2);
		perm_rep_GL_n_q_generators(
			TRUE /* f_affine */ , 
			TRUE /* f_semi */, 
			FALSE /* f_special */ , 
			TRUE /* f_translations */, 
			val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */ );
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, AGGL %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_AFFINE_TRANSLATIONS) {
		sprintf(g_label + strlen(g_label), "T_%ld_%ld", val1, val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm T}_{%ld}(%ld)", val1, val2);
		perm_rep_T_n_q_generators(val2 /* q */, val1 /* n */, &gen, FALSE /* f_verbose */);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, T %ld %ld: %s %s \n", 
			*S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_GROUP_PSU_3_Q2) {
		sprintf(g_label + strlen(g_label), "PSU_3_%ld_2", val2);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm PSU}_{3}(%ld^2)", val2);
		PSU_3_q2_generators(val2, &gen);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, PSU_3_%ld_2: %s %s \n", 
			*S_len, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else
		return FALSE;
	return TRUE;
}

#define SOLID_RADIUS 10000

#if TEXDOCU
INT km_get_group_solid(INT type, INT val1, INT val2, BYTE *s, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len)
#else
\begin{center}
\begin{tabular}{ll}
type & group \\
\hline
FGA\_TETRAHEDRON &  \\
FGA\_CUBE &  \\
FGA\_OCTAHEDRON &  \\
FGA\_DODECAHEDRON &  \\
FGA\_ICOSAHEDRON &  \\
\end{tabular}
\end{center}
#endif
{
	SOLID_OB A;
	SOLID_OP p;
	VECTOR_OP gen;
	PERMUTATION_OP per;
	double r;
	
	if (type == FGA_TETRAHEDRON) {
		sprintf(g_label + strlen(g_label), "Tetra");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Tetra}");
		A.tetrahedron(SOLID_RADIUS);
		A.copy((SOLID_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Tetrahedron %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_CUBE) {
		sprintf(g_label + strlen(g_label), "Cube");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Cube}");
		A.cube(SOLID_RADIUS);
		A.copy((SOLID_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Cube %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_OCTAHEDRON) {
		sprintf(g_label + strlen(g_label), "Octa");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Octa}");
		A.octahedron(SOLID_RADIUS);
		A.copy((SOLID_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Octahedron %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_DODECAHEDRON) {
		sprintf(g_label + strlen(g_label), "Dode");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Dode}");
		printf("calling dodecahedron()\n");
		A.dodecahedron(SOLID_RADIUS);
		A.copy((SOLID_OP) S->s_i(*S_len));
		printf("after copy()\n");
		
		S->s_i(*S_len)->println();
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Dodecahedron %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_ICOSAHEDRON) {
		sprintf(g_label + strlen(g_label), "Ico");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Ico}");
		A.icosahedron(SOLID_RADIUS);
		A.copy((SOLID_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Icosahedron %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_CUBE_4D) {
		sprintf(g_label + strlen(g_label), "Cube4D");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Cube4D}");
		A.cube4D(SOLID_RADIUS, SOLID_RADIUS * 2);
		A.copy((SOLID_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Cube4D %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_CUBUS_SIMUS) {
		sprintf(g_label + strlen(g_label), "Cubussimus");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Cubussimus}");
		printf("calling cubus_simus()\n");
		A.cubus_simus(1000.);
		A.copy((SOLID_OP) S->s_i(*S_len));
		printf("after copy()\n");
		
		S->s_i(*S_len)->println();
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Cubussimus %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_DODE_SIMUM) {
		sprintf(g_label + strlen(g_label), "Dodesimum");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf Dodesimum}");
		printf("calling cubus_simus()\n");
		A.dode_simum(1000.);
		// A.cubus_simus(1000.);
		// A.russian_snub_cube(1000.);
		A.copy((SOLID_OP) S->s_i(*S_len));
		printf("after copy()\n");
		
		S->s_i(*S_len)->println();
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, Dodesimum %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_CUBUS_EE) {
		sprintf(g_label + strlen(g_label), "CubusEE");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf CubusEE}");
		printf("calling cubus_simus()\n");
		A.snub_cube(1000.);
		// A.russian_snub_cube(1000.);
		A.copy((SOLID_OP) S->s_i(*S_len));
		printf("after copy()\n");
		
		S->s_i(*S_len)->println();
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, CubusEE %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_CUBUS_EE_RUSSIAN) {
		sprintf(g_label + strlen(g_label), "CubusEErussian");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf CubusEErussian}");
		printf("calling cubus_simus()\n");
		A.russian_snub_cube(1000.);
		A.copy((SOLID_OP) S->s_i(*S_len));
		printf("after copy()\n");
		
		S->s_i(*S_len)->println();
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, CubusEE %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}

	else if (type == FGA_SOLID_TRUNCATE) {
		sprintf(g_label + strlen(g_label), "T");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf trunc}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("SOLID_TRUNCATE not applicable");
		r = 1. / 3.;
		p->cut_vertices(r, &A);
		A.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, truncate %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_DUAL) {
		sprintf(g_label + strlen(g_label), "D");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf dual}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("SOLID_DUAL not applicable");
		p->dual(&A);
		A.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, dual %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_TRUNCATE_DODE) {
		sprintf(g_label + strlen(g_label), "T");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf trunc}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("SOLID_TRUNCATE_DODE not applicable");
		r = 1. / (2. + sqrt(2. * (1. - cos_grad(108))));
		p->cut_vertices(r, &A);
		A.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, truncate %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_TRUNCATE_CUBE) {
		sprintf(g_label + strlen(g_label), "T");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf trunc}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("SOLID_TRUNCATE_CUBE not applicable");
		r = 1. / (2. + sqrt(2.));
		p->cut_vertices(r, &A);
		A.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, truncate %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_RELABEL_POINTS) {
		sprintf(g_label + strlen(g_label), "relabel");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf relabel}");
		p = (SOLID_OP) S->s_i(*S_len - 2);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("FGA_SOLID_RELABEL_POINTS not applicable: arg1 must be a solid");
		gen = (SOLID_OP) S->s_i(*S_len - 1);
		if (gen->s_obj_k() != VECTOR) 
			return error("FGA_SOLID_RELABEL_POINTS not applicable: arg2 not a permutation");
		if (gen->s_li() != 1)
			return error("FGA_SOLID_RELABEL_POINTS not applicable: gen->s_li() != 1");
		per = (PERMUTATION_OP) gen->s_i(0);
		if (per->s_obj_k() != PERMUTATION) 
			return error("FGA_SOLID_RELABEL_POINTS not applicable: per is not a permutation");
		p->relabel_points(&A, per, FALSE /* f_relabel_vertex_labels */);
		A.copy((SOLID_OP) S->s_i(*S_len - 2));
		(*S_len)--;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, relabel %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_INDUCED_GROUP_ON_EDGES) {
		VECTOR_OB gen;
		
		sprintf(g_label + strlen(g_label), "onedges");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf onedges}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("FGA_SOLID_INDUCED_GROUP_ON_EDGES not applicable");
		p->induced_group_on_edges(p->s_group_generators(), &gen);
		gen.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, onedges %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_EDGE_MIDPOINTS) {
		sprintf(g_label + strlen(g_label), "E");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf edgemidpoints}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("FGA_SOLID_EDGE_MIDPOINTS not applicable");
		p->edge_midpoints(&A);
		A.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, edgemidpoints %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLID_ADD_CENTRAL_POINT) {
		sprintf(g_label + strlen(g_label), "addcenter");
		sprintf(g_label_tex + strlen(g_label_tex), "{\\bf addcenter}");
		p = (SOLID_OP) S->s_i(*S_len - 1);
		if (p->s_obj_k() != SOLID_KIND) 
			return error("FGA_SOLID_ADD_CENTRAL_POINT not applicable");
		p->add_central_point(&A);
		A.copy((SOLID_OP) S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, addcenter %s %s \n", 
			*S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else
		return FALSE;
	return TRUE;
}


#if TEXDOCU
INT km_get_group_generators_from_file(INT type, INT val1, INT val2, BYTE *fname, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len)
#else
if type is FGA\_GROUP\_FROM\_FILE:
reads generators from file and puts a vector on top of the stack S.
returns TRUE if something has been done, FALSE otherwise.
#endif
{
	VECTOR_OB gen;

	if (type == FGA_GROUP_FROM_FILE) {
		sprintf(g_label + strlen(g_label), "file_%s", fname);
		sprintf(g_label_tex + strlen(g_label_tex), "{\\rm file %s}", fname);
		read_file_of_generators(&gen, fname);
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, FILE %s: %s %s \n", 
			*S_len, fname, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else
		return FALSE;
	return TRUE;
}

#if TEXDOCU
INT km_binary_operations(INT type, INT val1, INT val2, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len)
#else
\begin{center}
\begin{tabular}{ll}
type & action \\
\hline
FGA\_DIRECT\_SUM & direct sum of topmost two groups \\
FGA\_DIRECT\_PRODUCT & direct product of topmost two groups \\
FGA\_WREATH\_PRODUCT & wreath product of topmost two groups \\
FGA\_EXPONENTIATION & exponentiation of topmost two groups \\
FGA\_ON\_MAPPINGS & action on mappings of topmost two groups \\
FGA\_COMMA & concatenation of topmost two generating sets \\
\end{tabular}
\end{center}
Returns TRUE if something has been done.
#endif
{
	VECTOR_OB gen;
	SYM_OP pa, pb;
	SYM_OB c;
	VECTOR_OP a, b;
	
	if (type == FGA_DIRECT_SUM) {
		sprintf(g_label + strlen(g_label), "x");
		sprintf(g_label_tex + strlen(g_label_tex), "x");
		if (*S_len < 2)
			return error("get_group(): DIRECT_SUM: 2 args must be given !");
		pa = S->s_i(*S_len - 2);
		pb = S->s_i(*S_len - 1);
		if (pa->s_obj_k() == SOLID_KIND && pb->s_obj_k() == SOLID_KIND) {
			SOLID_OP A, B;
			
			A = (SOLID_OP) pa;
			B = (SOLID_OP) pb;
			A->direct_sum(B, (SOLID_OP) &c);
			}
		else if (pa->s_obj_k() == VECTOR && pb->s_obj_k() == VECTOR) {
			vec_generators_direct_sum((VECTOR_OP)pa, (VECTOR_OP)pb, (VECTOR_OP) &c);
			}
		else {
			return error("FGA_DIRECT_SUM wrong object kinds");
			}
		// vec_generators_direct_sum(a, b, &gen);
		S->s_i(*S_len - 2)->freeself();
		S->s_i(*S_len - 1)->freeself();
		*S_len = *S_len - 2;
		c.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, DIRECT_SUM: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_DIRECT_PRODUCT) {
		printf("X\n"); fflush(stdout);
		sprintf(g_label + strlen(g_label), "X");
		sprintf(g_label_tex + strlen(g_label_tex), "\\times");
		if (*S_len < 2)
			return error("get_group(): DIRECT_PRODUCT: 2 args must be given !");
		pa = S->s_i(*S_len - 2);
		pb = S->s_i(*S_len - 1);
		if (pa->s_obj_k() == VECTOR && pb->s_obj_k() == SOLID_KIND) {
			SOLID_OP B;
			
			B = (SOLID_OP) pb;
			B->direct_product((VECTOR_OP) pa, (SOLID_OP) &c);
			}
		else if (pa->s_obj_k() == VECTOR && pb->s_obj_k() == VECTOR) {
			vec_generators_direct_product((VECTOR_OP)pa, (VECTOR_OP)pb, (VECTOR_OP) &c);
			}
		else {
			return error("FGA_DIRECT_PRODUCT wrong object kinds");
			}
		S->s_i(*S_len - 2)->freeself();
		S->s_i(*S_len - 1)->freeself();
		*S_len = *S_len - 2;
		c.copy(S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, DIRECT_PRODUCT: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		printf("X finished\n"); fflush(stdout);
		}
	else if (type == FGA_WREATH_PRODUCT) {
		sprintf(g_label + strlen(g_label), "wr");
		sprintf(g_label_tex + strlen(g_label_tex), "\\wr");
		if (*S_len < 2)
			return error("get_group(): WREATH_PRODUCT: 2 groups must be given !");
		a = (VECTOR_OP) S->s_i(*S_len - 2);
		b = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_generators_wreath_product(&gen, a, b, TRUE /* f_v */);
		S->s_i(*S_len - 2)->freeself();
		S->s_i(*S_len - 1)->freeself();
		*S_len = *S_len - 2;
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, WREATH_PRODUCT: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_EXPONENTIATION) {
		sprintf(g_label + strlen(g_label), "exp");
		sprintf(g_label_tex + strlen(g_label_tex), "\\uparrow");
		if (*S_len < 2)
			return error("get_group(): EXPONENTIATION: 2 groups must be given !");
		a = (VECTOR_OP) S->s_i(*S_len - 2);
		b = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_generators_exponentiation(a, b, &gen, FALSE);
		S->s_i(*S_len - 2)->freeself();
		S->s_i(*S_len - 1)->freeself();
		*S_len = *S_len - 2;
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, EXPONENTIATION: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_ON_MAPPINGS) {
		sprintf(g_label + strlen(g_label), "on_mappings");
		sprintf(g_label_tex + strlen(g_label_tex), "on n^m");
		if (*S_len < 2)
			return error("get_group(): ON_MAPPINGS: 2 groups must be given !");
		a = (VECTOR_OP) S->s_i(*S_len - 2);
		b = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_generators_exponentiation(a, b, &gen, TRUE);
		S->s_i(*S_len - 2)->freeself();
		S->s_i(*S_len - 1)->freeself();
		*S_len = *S_len - 2;
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, ON_MAPPINGS: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_COMMA) {
		sprintf(g_label + strlen(g_label), ",");
		sprintf(g_label_tex + strlen(g_label_tex), ",");
		if (*S_len < 2)
			return error("get_group(): COMMA: 2 groups must be given !");
		a = (VECTOR_OP) S->s_i(*S_len - 2);
		b = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_generators_comma(a, b, &gen);
		S->s_i(*S_len - 2)->freeself();
		S->s_i(*S_len - 1)->freeself();
		*S_len = *S_len - 2;
		gen.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, COMMA: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else
		return FALSE;
	return TRUE;
}


#if TEXDOCU
INT km_unary_operations(INT type, INT val1, INT val2, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len)
#else
\begin{center}
\begin{tabular}{ll}
type & action (on topmost group) \\
\hline
FGA\_INDUCE\_3\_SETS & induced action on 3-subsets \\
FGA\_INDUCE\_2\_SETS & induced action on 2-subsets \\
FGA\_INDUCE\_2\_TUPLES & induced action on 2-tuples \\
FGA\_INDUCE\_INJECTIVE\_2\_TUPLES & induced action on injective 2-tuples \\
FGA\_ADD\_FIXPOINT & add a fixed point\\
FGA\_STABILIZE\_POINT & point stabilizer of first point \\
FGA\_SELECT\_ITH & select $i$-th element of topmost vector \\
FGA\_SOLVABLE\_GROUP\_\\
$\;$PERMUTATION\_REPRESENTATION & \\
FGA\_ABSTRACT\_PRESENTATION & \\
FGA\_POWER & generators for the group generated by the $i$-th powers\\
\end{tabular}
\end{center}
Returns TRUE if something has been done.
#endif
{
	VECTOR_OB gen, gen2, tmp;
	VECTOR_OP a;
	
	if (type == FGA_INDUCE_3_SETS) {
		sprintf(g_label + strlen(g_label), "on3sets");
		sprintf(g_label_tex + strlen(g_label_tex), "^{[3]}");
		if (*S_len == 0)
			return error("get_group(): INDUCE_3_SETS: no group defined");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_induce3(a);
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, INDUCE_3_SETS: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_INDUCE_2_SETS) {
		sprintf(g_label + strlen(g_label), "on2sets");
		sprintf(g_label_tex + strlen(g_label_tex), "^{[2]}");
		if (*S_len == 0)
			return error("get_group(): INDUCE_2_SETS: no group defined");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_induce2(a);
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, INDUCE_2_SETS: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_INDUCE_2_TUPLES) {
		sprintf(g_label + strlen(g_label), "on2tuples");
		sprintf(g_label_tex + strlen(g_label_tex), "^{(2)}");
		if (*S_len == 0)
			return error("get_group(): INDUCE_2_TUPLES: no group defined");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_induce_on_2tuples(a, FALSE /* f_injective */);
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, INDUCE_2_TUPLES: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_INDUCE_INJECTIVE_2_TUPLES) {
		sprintf(g_label + strlen(g_label), "oninj2tuples");
		sprintf(g_label_tex + strlen(g_label_tex), "^{(2)_i}");
		if (*S_len == 0)
			return error("get_group(): INDUCE_INJECTIVE_2_TUPLES: no group defined");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_induce_on_2tuples(a, TRUE /* f_injective */);
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, INDUCE_INJECTIVE_2_TUPLES: %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_ADD_FIXPOINT) {
		sprintf(g_label + strlen(g_label), "p");
		sprintf(g_label_tex + strlen(g_label_tex), "+");
		if (*S_len == 0)
			return error("get_group(): ADD_FIXPOINT: no group defined");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_add_fixpoint_in_front(a);
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, ADD_FIXPOINT (in front): %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_STABILIZE_POINT) {
		sprintf(g_label + strlen(g_label), "-");
		sprintf(g_label_tex + strlen(g_label_tex), "-");
		if (*S_len == 0)
			return error("get_group(): STABILIZE_POINT: no group defined");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		vec_generators_stabilize_point(a, &gen);
		vec_generators_remove_fixpoint(&gen, &gen2, 0);
		S->s_i(*S_len - 1)->freeself();
		(*S_len)--;
		gen2.copy((VECTOR_OP) S->s_i(*S_len));
		(*S_len)++;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, STABILIZE_POINT %s %s \n", *S_len, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SELECT_ITH) {
		sprintf(g_label + strlen(g_label), "%ld-th", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "%ld-th", val1);
		if (*S_len == 0)
			return error("get_group(): SELECT_ITH stack is empty");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		if (vector_select_ith(a, val1) != OK)
			return ERROR;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, SELECT_ITH %ld, %s %s \n", *S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_SOLVABLE_GROUP_PERMUTATION_REPRESENTATION) {
		sprintf(g_label + strlen(g_label), "on%ldpt", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "perm reps on $\\le$ %ld points", val1);
		if (*S_len == 0)
			return error("get_group(): FGA_SOLVABLE_GROUP_PERMUTATION_REPRESENTATION stack is empty");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		if (solvable_group_permutation_representations(a, val1, val2) != OK)
			return ERROR;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, FGA_SOLVABLE_GROUP_PERMUTATION_REPRESENTATION %ld ordering %ld, %s %s \n", *S_len, val1, val2, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else if (type == FGA_ABSTRACT_PRESENTATION) {
#ifdef RUCKDESCHEL_TRUE
		sprintf(g_label + strlen(g_label), "abstract%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "abstract presentation %ld", val1);
		if (*S_len == 0)
			return error("get_group(): FGA_ABSTRACT_PRESENTATION stack is empty");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		if (perms_to_abstract_presentation(a, &tmp, val1) != OK)
			return ERROR;
		tmp.swap(S->s_i(*S_len - 1));
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, FGA_ABSTRACT_PRESENTATION strategy %ld, %s %s \n", *S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
#else
		return error("FGA_ABSTRACT_PRESENTATION not available");
#endif
		}
	else if (type == FGA_POWER) {
		sprintf(g_label + strlen(g_label), "power%ld", val1);
		sprintf(g_label_tex + strlen(g_label_tex), "^{%ld}", val1);
		if (*S_len == 0)
			return error("get_group(): FGA_POWER stack is empty");
		a = (VECTOR_OP) S->s_i(*S_len - 1);
		if (a->power_elementwise_apply(val1) != OK)
			return ERROR;
#ifdef DEBUG_GET_GROUP
		printf("get_group: S_len = %ld, FGA_POWER %ld, %s %s \n", *S_len, val1, g_label, g_label_tex);
		print_deg(S, *S_len);
		fflush(stdout);
#endif
		}
	else 
		return FALSE;
	return TRUE;
}

#if TEXDOCU
VECTOR_OP km_get_generators(VECTOR_OP S, INT S_len)
#endif
{
	SYM_OP q;
	SOLID_OP A;
	VECTOR_OP gen = NULL;

	q = S->s_i(S_len - 1);
	if (q->s_obj_k() == SOLID_KIND) {
		A = (SOLID_OP) q;
		gen = A->s_group_generators();
		}
	else if (q->s_obj_k() == VECTOR) {
		gen = (VECTOR_OP) q;
		}
	else {
		printf("km_get_generators() unknown object kind\n");
		q->printobjectkind();
		error("km_get_generators()");
		}
	return gen;
}

#if TEXDOCU
INT km_get_deg(VECTOR_OP S, INT S_len)
#endif
{
	VECTOR_OP gen;
	PERMUTATION_OP p;
	
	gen = km_get_generators(S, S_len);
	if (gen->s_obj_k() != VECTOR)
		return error("km_get_deg() gen->s_obj_k() != VECTOR");
	if (gen->s_li() < 1)
		return error("km_get_deg(): no generators");
	p = (PERMUTATION_OP) gen->s_i(0);
	return p->s_li();
}

#if TEXDOCU
INT print_deg(VECTOR_OP S, INT S_len)
#endif
{
	INT i;
	
	i = km_get_deg(S, S_len);
	printf("get_group: deg = %ld\n", i);
	fflush(stdout);
	return OK;
}

#define MAX_STACK_SIZE 100

#if TEXDOCU
INT km_get_group_from_selection(VECTOR_OP gsel, INT nb_g_sel, 
	VECTOR_OP generators, BYTE *g_label, BYTE *g_label_tex)
#endif
{
	GROUP_SELECTION_OP gs;
	VECTOR_OB S; // stack
	VECTOR_OB gen;
	VECTOR_OP a;
	SYM_OB go;
	INT S_len = 0;
	INT i, type, val1, val2;
	BYTE *s;
	BYTE str[1000];
	SOLID_OP p;

	S.m_il(MAX_STACK_SIZE);
	*g_label = 0;
	*g_label_tex = 0;
	for (i = 0; i < nb_g_sel; i++) {
		gs = (GROUP_SELECTION_OP) gsel->s_i(i);
		str[0] = 0;
		gs->sprint(str);
		printf("km_get_group_from_selection: i=%ld gsel = %s\n", i, str);
		fflush(stdout);
		type = gs->s_type_i();
		val1 = gs->s_val1_i();
		val2 = gs->s_val2_i();
		s = gs->s_s_s();

		if (km_get_group_generators1(type, val1, val2, s, g_label, g_label_tex, &S, &S_len))
			;
		else if (km_get_group_generators2(type, val1, val2, g_label, g_label_tex, &S, &S_len))
			;
		else if (km_get_group_solid(type, val1, val2, s, g_label, g_label_tex, &S, &S_len))
			;
		else if (type == FGA_GROUP_FROM_FILE) {
			km_get_group_generators_from_file(type, val1, val2, s, g_label, g_label_tex, &S, &S_len);
			}
		// binary operations:
		else if (km_binary_operations(type, val1, val2, g_label, g_label_tex, &S, &S_len))
			;
		// unary operations:
		else if (km_unary_operations(type, val1, val2, g_label, g_label_tex, &S, &S_len))
			;

		else {
			return error("unknown group selection type (or not yet implemented\n");
			}
#ifdef DEBUG_GET_GROUP
		if (S_len) {
			a = km_get_generators(&S, S_len);
			// a->Print();
			// vec_generators_group_order(a, &go);
			// printf("km_get_group_from_selection(): group_order = ");
			// go.println();
			}
#endif
	
		}

	printf("km_get_group_from_selection() S_len = %ld\n", S_len);
	fflush(stdout);
	if (S_len != 1)
		return error("km_get_group_from_selection(): S_len != 1");
	i = km_get_deg(&S, S_len);
	a = km_get_generators(&S, S_len);
	if (S.s_i(0)->s_obj_k() == SOLID_KIND) {
		BYTE str[10000];
		FILE *fp;
		
		sprintf(str, "%s.graph", g_label);
		p = (SOLID_OP) S.s_i(0);
		p->write_graphfile(str);
		fp = fopen("solid_name", "w");
		fprintf(fp, "%s\n", str);
		fclose(fp);
		sprintf(str, "%s.solid", g_label);
		p->save(str);
		}
	else {
		call_system("rm solid_name");
		}
	a->copy(generators);
	printf("finished with km_get_group_from_selection()\n");
	fflush(stdout);
	return OK;
}

static void add_gsel(VECTOR_OP V, INT type, INT val1, INT val2, BYTE *s);

#if TEXDOCU
INT compose_gsel_from_strings(VECTOR_OP V, INT num_args, BYTE **args)
#else
\begin{center}
\begin{tabular}{ll}
description & mapped onto \\
\hline
PSL & FGA\_GROUP\_PSL \\
PGL & FGA\_GROUP\_PGL \\
PSSL & FGA\_GROUP\_PSSL \\
PGGL & FGA\_GROUP\_PGGL \\
SL & FGA\_GROUP\_SL \\
GL & FGA\_GROUP\_GL \\
SSL & FGA\_GROUP\_SSL \\
GGL & FGA\_GROUP\_GGL \\
ASL & FGA\_GROUP\_ASL \\
AGL & FGA\_GROUP\_AGL \\
ASSL & FGA\_GROUP\_ASSL \\
AGGL & FGA\_GROUP\_AGGL \\
T & FGA\_GROUP\_AFFINE\_TRANSLATIONS \\
PSU\_3\_Q2 & FGA\_GROUP\_PSU\_3\_Q2 \\
SYM or S & FGA\_GROUP\_SYM \\
SYM\_TRANSPOSITIONS & FGA\_GROUP\_SYM\_TRANSPOSITIONS \\
ALT or A & FGA\_GROUP\_ALT \\
DIHEDRAL or D & FGA\_GROUP\_DIHEDRAL \\
CYCLIC or C & FGA\_GROUP\_CYCLIC \\
HOLOMORPH\_OF\_ & FGA\_GROUP\_HOLOMORPH\_\\
$\;$CYCLIC\_GROUP or HolC & $\;$OF\_CYCLIC\_GROUP \\
TRIVIAL or Id & FGA\_GROUP\_TRIVIAL \\
SUBGROUP\_OF\_HOLOMORPH\_ & FGA\_SUBGROUP\_OF\_HOLOMORPH\_ \\
$\;$OF\_CYCLIC\_GROUP or HolCsubgroup & $\;$OF\_CYCLIC\_GROUP \\
SN\_WREATH\_SM or SnwrSm & FGA\_GROUP\_SN\_WREATH\_SM \\
PERM & FGA\_GROUP\_BY\_PERMUTATION\_GENERATOR \\
M11 or M12 or M23 or M24 & FGA\_GROUP\_MATHIEU \\
HS176 & FGA\_GROUP\_HIGMAN\_SIMS\_176 \\
SOLVABLE\_GROUP or solvable & FGA\_SOLVABLE\_GROUP \\
file & FGA\_GROUP\_FROM\_FILE \\
INDUCE\_2\_SETS or [2] & FGA\_INDUCE\_2\_SETS \\
INDUCE\_2\_TUPLES or (2) & FGA\_INDUCE\_2\_TUPLES \\
INDUCE\_INJECTIVE\_ & FGA\_INDUCE\_INJECTIVE\_\\
$\;$2\_TUPLES or  (2)i & $\;$2\_TUPLES \\
INDUCE\_3\_SETS or [3] & FGA\_INDUCE\_3\_SETS \\
ADD\_FIXPOINT or + & FGA\_ADD\_FIXPOINT \\
STABILIZE\_POINT or - & FGA\_STABILIZE\_POINT \\
POWER & FGA\_POWER \\
COMMA or , & FGA\_COMMA \\
DIRECT\_SUM or x & FGA\_DIRECT\_SUM \\
DIRECT\_PRODUCT or X & FGA\_DIRECT\_PRODUCT \\
WREATH\_PRODUCT or wr & FGA\_WREATH\_PRODUCT \\
EXPONENTIATION or exp & FGA\_EXPONENTIATION \\
ON\_MAPPINGS & FGA\_ON\_MAPPINGS \\
\end{tabular}
\end{center}
#endif
{
	INT i = 0, val1, val2;
	BYTE *p, *q, *fname;
	BYTE *str;
	
	V->m_il(0);
	while (i < num_args) {
		p = args[i++];
		val1 = 0;
		val2 = 0;
		str = NIL;
		
		// the linear groups:
		if (strcmp(p, "PSL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_PSL, val1, val2, NIL);
			}
		else if (strcmp(p, "PGL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_PGL, val1, val2, NIL);
			}
		else if (strcmp(p, "PSSL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_PSSL, val1, val2, NIL);
			}
		else if (strcmp(p, "PGGL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_PGGL, val1, val2, NIL);
			}
		else if (strcmp(p, "SL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_SL, val1, val2, NIL);
			}
		else if (strcmp(p, "GL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_GL, val1, val2, NIL);
			}
		else if (strcmp(p, "SSL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_SSL, val1, val2, NIL);
			}
		else if (strcmp(p, "GGL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_GGL, val1, val2, NIL);
			}
		else if (strcmp(p, "ASL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_ASL, val1, val2, NIL);
			}
		else if (strcmp(p, "AGL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_AGL, val1, val2, NIL);
			}
		else if (strcmp(p, "ASSL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_ASSL, val1, val2, NIL);
			}
		else if (strcmp(p, "AGGL") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_AGGL, val1, val2, NIL);
			}
		else if (strcmp(p, "T") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_AFFINE_TRANSLATIONS, val1, val2, NIL);
			}
		else if (strcmp(p, "PSU_3_Q2") == 0) {
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_PSU_3_Q2, val1, val2, NIL);
			}
		// missing:
		// FGA_GROUP_SZ_Q
		
		// the well known groups:
		else if (strcmp(p, "SYM") == 0 || strcmp(p, "S") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_SYM, val1, val2, NIL);
			}
		else if (strcmp(p, "SYM_TRANSPOSITIONS") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_SYM_TRANSPOSITIONS, val1, val2, NIL);
			}
		else if (strcmp(p, "ALT") == 0 || strcmp(p, "A") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_ALT, val1, val2, NIL);
			}
		else if (strcmp(p, "DIHEDRAL") == 0 || strcmp(p, "D") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_DIHEDRAL, val1, val2, NIL);
			}
		else if (strcmp(p, "CYCLIC") == 0 || strcmp(p, "C") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_CYCLIC, val1, val2, NIL);
			}
		else if (strcmp(p, "HOLOMORPH_OF_CYCLIC_GROUP") == 0 || strcmp(p, "HolC") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_HOLOMORPH_OF_CYCLIC_GROUP, val1, val2, NIL);
			}
		else if (strcmp(p, "TRIVIAL") == 0 || strcmp(p, "Id") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_GROUP_TRIVIAL, val1, val2, NIL);
			}
		else if (strcmp(p, "SUBGROUP_OF_HOLOMORPH_OF_CYCLIC_GROUP") == 0 || strcmp(p, "HolCsubgroup") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_SUBGROUP_OF_HOLOMORPH_OF_CYCLIC_GROUP, val1, val2, NIL);
			}
		else if (strcmp(p, "SN_WREATH_SM") == 0 || strcmp(p, "SnwrSm") == 0) {
			sscanf(args[i++], "%ld", &val1);
			sscanf(args[i++], "%ld", &val2);
			add_gsel(V, FGA_GROUP_SN_WREATH_SM, val1, val2, NIL);
			}
		else if (strcmp(p, "PERM") == 0) {
			// sscanf(args[i++], "%ld", &val1);
			// sscanf(args[i++], "%ld", &val2);
			val1 = 0;
			val2 = 0;
			q = args[i++];
			add_gsel(V, FGA_GROUP_BY_PERMUTATION_GENERATOR, val1, val2, q);
			}
			
		// solid:
		else if (strcmp(p, "TETRA") == 0 || strcmp(p, "TETRAHEDRON") == 0) {
			add_gsel(V, FGA_TETRAHEDRON, val1, val2, NIL);
			}
		else if (strcmp(p, "CUBE") == 0) {
			add_gsel(V, FGA_CUBE, val1, val2, NIL);
			}
		else if (strcmp(p, "OCTA") == 0 || strcmp(p, "OCTAHEDRON") == 0) {
			add_gsel(V, FGA_OCTAHEDRON, val1, val2, NIL);
			}
		else if (strcmp(p, "DODE") == 0 || strcmp(p, "DODECAHEDRON") == 0) {
			add_gsel(V, FGA_DODECAHEDRON, val1, val2, NIL);
			}
		else if (strcmp(p, "ICO") == 0 || strcmp(p, "ICOSAHEDRON") == 0) {
			add_gsel(V, FGA_ICOSAHEDRON, val1, val2, NIL);
			}
		else if (strcmp(p, "TRUNCATE") == 0) {
			add_gsel(V, FGA_SOLID_TRUNCATE, val1, val2, NIL);
			}
		else if (strcmp(p, "DUAL") == 0) {
			add_gsel(V, FGA_SOLID_DUAL, val1, val2, NIL);
			}
		else if (strcmp(p, "TRUNCATE_DODE") == 0) {
			add_gsel(V, FGA_SOLID_TRUNCATE_DODE, val1, val2, NIL);
			}
		else if (strcmp(p, "TRUNCATE_CUBE") == 0) {
			add_gsel(V, FGA_SOLID_TRUNCATE_CUBE, val1, val2, NIL);
			}
		else if (strcmp(p, "CUBE4D") == 0) {
			add_gsel(V, FGA_CUBE_4D, val1, val2, NIL);
			}
		else if (strcmp(p, "RELABEL_POINTS") == 0) {
			add_gsel(V, FGA_SOLID_RELABEL_POINTS, val1, val2, NIL);
			}
		else if (strcmp(p, "ON_EDGES") == 0) {
			add_gsel(V, FGA_SOLID_INDUCED_GROUP_ON_EDGES, val1, val2, NIL);
			}
		else if (strcmp(p, "EDGE_MIDPOINTS") == 0) {
			add_gsel(V, FGA_SOLID_EDGE_MIDPOINTS, val1, val2, NIL);
			}
		else if (strcmp(p, "ADD_CENTRAL_POINT") == 0) {
			add_gsel(V, FGA_SOLID_ADD_CENTRAL_POINT, val1, val2, NIL);
			}
		else if (strcmp(p, "CUBUS_SIMUS") == 0) {
			add_gsel(V, FGA_SOLID_CUBUS_SIMUS, val1, val2, NIL);
			}
		else if (strcmp(p, "DODE_SIMUM") == 0) {
			add_gsel(V, FGA_SOLID_DODE_SIMUM, val1, val2, NIL);
			}
		else if (strcmp(p, "CUBUS_EE") == 0) {
			add_gsel(V, FGA_SOLID_CUBUS_EE, val1, val2, NIL);
			}
		else if (strcmp(p, "CUBUS_EE_RUSSIAN") == 0) {
			add_gsel(V, FGA_SOLID_CUBUS_EE_RUSSIAN, val1, val2, NIL);
			}
		
		// sporadic groups:
		else if (strcmp(p, "M11") == 0) {
			val1 = 11;
			add_gsel(V, FGA_GROUP_MATHIEU, val1, val2, NIL);
			}
		else if (strcmp(p, "M12") == 0) {
			val1 = 12;
			add_gsel(V, FGA_GROUP_MATHIEU, val1, val2, NIL);
			}
		else if (strcmp(p, "M23") == 0) {
			val1 = 23;
			add_gsel(V, FGA_GROUP_MATHIEU, val1, val2, NIL);
			}
		else if (strcmp(p, "M24") == 0) {
			val1 = 24;
			add_gsel(V, FGA_GROUP_MATHIEU, val1, val2, NIL);
			}
		else if (strcmp(p, "HS176") == 0) {
			add_gsel(V, FGA_GROUP_HIGMAN_SIMS_176, val1, val2, NIL);
			}
		
		// solvable group:
		else if (strcmp(p, "SOLVABLE_GROUP") == 0 || strcmp(p, "solvable") == 0) {
			sscanf(args[i++], "%ld", &val1); // order
			sscanf(args[i++], "%ld", &val2); // catalogue number
			add_gsel(V, FGA_SOLVABLE_GROUP, val1, val2, NIL);
			}
		
		// from file:
		else if (strcmp(p, "file") == 0) {
			fname = args[i++];
			val1 = 0;
			val2 = 0;
			add_gsel(V, FGA_GROUP_FROM_FILE, val1, val2, fname);
			}
		
		// unary operators:
		else if (strcmp(p, "INDUCE_2_SETS") == 0 || strcmp(p, "[2]") == 0) {
			add_gsel(V, FGA_INDUCE_2_SETS, val1, val2, NIL);
			}
		else if (strcmp(p, "INDUCE_2_TUPLES") == 0 || strcmp(p, "(2)") == 0) {
			add_gsel(V, FGA_INDUCE_2_TUPLES, val1, val2, NIL);
			}
		else if (strcmp(p, "INDUCE_INJECTIVE_2_TUPLES") == 0 || strcmp(p, "(2)i") == 0) {
			add_gsel(V, FGA_INDUCE_INJECTIVE_2_TUPLES, val1, val2, NIL);
			}
		else if (strcmp(p, "INDUCE_3_SETS") == 0 || strcmp(p, "[3]") == 0) {
			add_gsel(V, FGA_INDUCE_3_SETS, val1, val2, NIL);
			}
		else if (strcmp(p, "ADD_FIXPOINT") == 0 || strcmp(p, "+") == 0) {
			add_gsel(V, FGA_ADD_FIXPOINT, val1, val2, NIL);
			}
		else if (strcmp(p, "STABILIZE_POINT") == 0 || strcmp(p, "-") == 0) {
			add_gsel(V, FGA_STABILIZE_POINT, val1, val2, NIL);
			}
		else if (strcmp(p, "POWER") == 0) {
			sscanf(args[i++], "%ld", &val1);
			add_gsel(V, FGA_POWER, val1, val2, NIL);
			}
		
		
		
		// binary operators:
		else if (strcmp(p, "COMMA") == 0 || strcmp(p, ",") == 0) {
			add_gsel(V, FGA_COMMA, val1, val2, NIL);
			}
		else if (strcmp(p, "DIRECT_SUM") == 0 || strcmp(p, "x") == 0) {
			add_gsel(V, FGA_DIRECT_SUM, val1, val2, NIL);
			}
		else if (strcmp(p, "DIRECT_PRODUCT") == 0 || strcmp(p, "X") == 0) {
			add_gsel(V, FGA_DIRECT_PRODUCT, val1, val2, NIL);
			}
		else if (strcmp(p, "WREATH_PRODUCT") == 0 || strcmp(p, "wr") == 0) {
			add_gsel(V, FGA_WREATH_PRODUCT, val1, val2, NIL);
			}
		else if (strcmp(p, "EXPONENTIATION") == 0 || strcmp(p, "exp") == 0) {
			add_gsel(V, FGA_EXPONENTIATION, val1, val2, NIL);
			}
		else if (strcmp(p, "ON_MAPPINGS") == 0) {
			add_gsel(V, FGA_ON_MAPPINGS, val1, val2, NIL);
			}
		
		}
	return OK;
}

static void add_gsel(VECTOR_OP V, INT type, INT val1, INT val2, BYTE *s)
{
	GROUP_SELECTION_OP gs;
	BYTE str[10000];
	INT i;
	
	i = V->s_li();
	V->inc();
	gs = (GROUP_SELECTION_OP) V->s_i(i);
	gs->init(type, val1, val2, s);
	
	str[0] = 0;
	// group_selection_print(gs, str);
	gs->sprint(str);
	printf("added item %ld: %s\n", i, str);
	fflush(stdout);
}

#if TEXDOCU
INT GROUP_SELECTION_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "type"; break;
	case 1: s = "val1"; break;
	case 2: s = "val2"; break;
	case 3: s = "s"; break;
	default:
		return error("GROUP_SELECTION_OB::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT GROUP_SELECTION_OB::init(INT type, INT val1, INT val2, BYTE *s)
#endif
{
	m_il(4);
	c_obj_k(GROUP_SELECTION_KIND);
	
	s_type()->m_i(type);
	s_val1()->m_i(val1);
	s_val2()->m_i(val2);
	if (s)
		s_s()->init(s);
	else
		s_s()->init("");
	return OK;
}

#if TEXDOCU
INT GROUP_SELECTION_OB::sprint(BYTE *s)
#endif
{
	INT type, val1, val2;
	BYTE str[1024];
	BYTE *p;
	
	type =s_type_i();
	val1 = s_val1_i();
	val2 = s_val2_i();
	p = s_s_s();
	if (type == FGA_GROUP_MATHIEU) {
		sprintf(str, "M%ld", val1);
		}
	else if (type == FGA_GROUP_HIGMAN_SIMS_176) {
		sprintf(str, "HS_176");
		}
	else if (type == FGA_GROUP_PSL) {
		sprintf(str, "PSL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_PGL) {
		sprintf(str, "PGL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_PSSL) {
		sprintf(str, "PSSL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_PGGL) {
		sprintf(str, "PGGL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_SL) {
		sprintf(str, "SL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_GL) {
		sprintf(str, "GL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_SSL) {
		sprintf(str, "SSL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_GGL) {
		sprintf(str, "GGL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_ASL) {
		sprintf(str, "ASL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_AGL) {
		sprintf(str, "AGL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_ASSL) {
		sprintf(str, "ASSL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_AGGL) {
		sprintf(str, "AGGL_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_AFFINE_TRANSLATIONS) {
		sprintf(str, "T_%ld_%ld", val1, val2);
		}
	else if (type == FGA_GROUP_PSU_3_Q2) {
		sprintf(str, "PSU_3_%ld", val2);
		}
	else if (type == FGA_GROUP_SZ_Q) {
		sprintf(str, "SZ_%ld", val2);
		}
	else if (type == FGA_GROUP_FROM_FILE) {
		sprintf(str, "from file '%s'", p);
		}
	else if (type == FGA_GROUP_SYM) {
		sprintf(str, "S%ld", val1);
		}
	else if (type == FGA_GROUP_SYM_TRANSPOSITIONS) {
		sprintf(str, "S%ld", val1);
		}
	else if (type == FGA_GROUP_ALT) {
		sprintf(str, "A%ld", val1);
		}
	else if (type == FGA_GROUP_DIHEDRAL) {
		sprintf(str, "D%ld", val1);
		}
	else if (type == FGA_GROUP_CYCLIC) {
		sprintf(str, "C%ld", val1);
		}
	else if (type == FGA_GROUP_TRIVIAL) {
		sprintf(str, "Id_%ld", val1);
		}
#if 0
	else if (type == FGA_GROUP_Zn_MULTIPLICATOR) {
		sprintf(str, "Znx%ld", val1);
		}
#endif
	else if (type == FGA_GROUP_SN_WREATH_SM) {
		sprintf(str, "S%ld_wr_S%ld", val1, val2);
		}
	else if (type == FGA_GROUP_HOLOMORPH_OF_CYCLIC_GROUP) {
		sprintf(str, "Hol(C%ld)", val1);
		}
	else if (type == FGA_TETRAHEDRON) {
		sprintf(str, "Tetra");
		}
	else if (type == FGA_CUBE) {
		sprintf(str, "Cube");
		}
	else if (type == FGA_OCTAHEDRON) {
		sprintf(str, "Octa");
		}
	else if (type == FGA_DODECAHEDRON) {
		sprintf(str, "Dode");
		}
	else if (type == FGA_ICOSAHEDRON) {
		sprintf(str, "Ico");
		}
	else if (type == FGA_SOLID_TRUNCATE) {
		sprintf(str, "truncate");
		}
	else if (type == FGA_SOLID_DUAL) {
		sprintf(str, "dual");
		}
	else if (type == FGA_SOLID_TRUNCATE_DODE) {
		sprintf(str, "truncate_dode");
		}
	else if (type == FGA_SOLID_TRUNCATE_CUBE) {
		sprintf(str, "truncate_cube");
		}
	else if (type == FGA_CUBE_4D) {
		sprintf(str, "cube4D");
		}
	else if (type == FGA_SOLID_RELABEL_POINTS) {
		sprintf(str, "relabel_points");
		}
	else if (type == FGA_SOLID_INDUCED_GROUP_ON_EDGES) {
		sprintf(str, "induced_group_on_edges");
		}
	else if (type == FGA_SOLID_EDGE_MIDPOINTS) {
		sprintf(str, "edge_midpoints");
		}
	else if (type == FGA_SOLID_ADD_CENTRAL_POINT) {
		sprintf(str, "add_central_point");
		}
	else if (type == FGA_SOLID_CUBUS_SIMUS) {
		sprintf(str, "cubus_simus");
		}
	else if (type == FGA_SOLID_DODE_SIMUM) {
		sprintf(str, "dode_simum");
		}
	else if (type == FGA_SOLID_CUBUS_EE) {
		sprintf(str, "cubus_ee");
		}
	else if (type == FGA_SOLID_CUBUS_EE_RUSSIAN) {
		sprintf(str, "cubus_ee_russian");
		}
	else if (type == FGA_DIRECT_SUM) {
		sprintf(str, "x");
		}
	else if (type == FGA_DIRECT_PRODUCT) {
		sprintf(str, "X");
		}
	else if (type == FGA_INDUCE_2_SETS) {
		sprintf(str, "^[2]");
		}
	else if (type == FGA_INDUCE_2_TUPLES) {
		sprintf(str, "^(2)");
		}
	else if (type == FGA_INDUCE_INJECTIVE_2_TUPLES) {
		sprintf(str, "^(2)i");
		}
	else if (type == FGA_INDUCE_3_SETS) {
		sprintf(str, "^[3]");
		}
	else if (type == FGA_ADD_FIXPOINT) {
		sprintf(str, "+");
		}
	else if (type == FGA_STABILIZE_POINT) {
		sprintf(str, "-");
		}
	else if (type == FGA_WREATH_PRODUCT) {
		sprintf(str, "wr");
		}
	else if (type == FGA_EXPONENTIATION) {
		sprintf(str, "exp");
		}
	else if (type == FGA_ON_MAPPINGS) {
		sprintf(str, "on n^m");
		}
	else if (type == FGA_HOLOMORPH) {
		sprintf(str, "hol");
		}
	else if (type == FGA_GROUP_BY_PERMUTATION_GENERATOR) {
		sprintf(str, "<%s>", p);
		}
	else if (type == FGA_COMMA) {
		sprintf(str, ",");
		}
	else if (type == FGA_SELECT_ITH) {
		sprintf(str, "%ld-th", val1);
		}
	else if (type == FGA_SOLVABLE_GROUP_PERMUTATION_REPRESENTATION) {
		sprintf(str, "on <= %ld points (lowindex order %ld)", val1, val2);
		}
	else if (type == FGA_SOLVABLE_GROUP) {
		sprintf(str, "solvable group %ld#%ld", val1, val2);
		}
	else if (type == FGA_ABSTRACT_PRESENTATION) {
		sprintf(str, "abstract presentation strategy=%ld", val1);
		}
	else if (type == FGA_SUBGROUP_OF_HOLOMORPH_OF_CYCLIC_GROUP) {
		sprintf(str, "Hol(C%ld) subgroup of index %ld", val1, val2);
		}
	else if (type == FGA_POWER) {
		sprintf(str, "power %ld ", val1);
		}
	else {
		sprintf(str, "unknown: type=%ld val1=%ld val2=%ld", type, val1, val2);
		}
	strcat(s, str);
	return OK;
	
}

#endif /* FGA_TRUE */


