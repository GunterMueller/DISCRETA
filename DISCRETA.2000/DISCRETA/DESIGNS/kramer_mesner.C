/* kramer_mesner.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#include <stdlib.h> // for system

#ifdef LADDER_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/ladder.h>

#ifdef SYM_GEO
#include <DISCRETA/geo.h>
#endif

char *discreta_copyright_text = 
"\n"
" #####      #     ####    ####   #####   ######   #####    ##\n"
" #    #     #    #       #    #  #    #  #          #     #  #\n"
" #    #     #     ####   #       #    #  #####      #    #    #\n"
" #    #     #         #  #       #####   #          #    ######\n"
" #    #     #    #    #  #    #  #   #   #          #    #    #\n"
" #####      #     ####    ####   #    #  ######     #    #    #\n"
"\n"
"  -  the design construction tool by\n"
"     Anton Betten <Anton.Betten@uni-bayreuth.de>\n"
"     Evi Haberberger <Evi.Haberberger@uni-bayreuth.de>\n"
"     Reinhard Laue <Reinhard.Laue@uni-bayreuth.de>\n"
"     Alfred Wassermann <Alfred.Wassermann@uni-bayreuth.de>\n"
"\n"
"  -  with additional programs and algorithms by\n"
"     Donald E. Knuth http://www-cs-faculty.stanford.edu/~knuth/\n"
"     Brendan D. McKay <bdm@cs.anu.edu.au>\n"
"\n"
"  -  for more information, see\n"
"     http://www.mathe2.uni-bayreuth.de/~discreta/\n" 
"\n";


#if TEXDOCU
void compute_KM(VECTOR_OP gen, BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, 
	INT f_strong_generators, 
	INT f_orderly_generation, 
	INT f_TDO, 
	INT f_extension_construction, 
	INT f_canonical_representatives, 
	INT f_k_k2_design, 
	INT k2)
#else
If f\_orderly\_generation is true, the function 
calc\_kramer\_mesner\_matrix\_by\_orderly\_generation() is called 
(this function is defined in dcc\_orderly.C).
Otherwise, the function calc\_kramer\_mesner\_matrix() is called.
#endif
{
	INT deg;
	SYM_OB go;
	MATRIX_OB TG;
	LABRA_OB labG;
	
	/* printf("compute_KM()\n"); */
	// deg = ((PERMUTATION_OP) gen->s_i(0))->s_li();
	deg = vec_generators_degree(gen);
	
	// vec_generators_group_order(gen, &go);
	
	printf("group = %s deg = %ld\n", g_label, deg);
	fflush(stdout);
	// go.println();
	

	reduce_generators_labra(gen, &go, FALSE /* f_verbose */, &labG);
	printf("group order of G: ");
	go.println();
	if (!go.einsp()) {
		printf("calling print_T: ");
		labG.print_T();
		printf("calling calc_transversal_matrix: ");
		labG.calc_transversal_matrix(&TG, TRUE /* f_v */);

		if (deg < 50) {
			TG.Print();
			TG.latex_upper_tri(stdout);
			labG.latex_transversal_indices(stdout);
			}	

		if (deg < 10) {
			printf("reducing the number of generators:\n");
			fflush(stdout);
			labG.reduced_generating_set(gen, FALSE /* f_bottom_up */, TRUE /* f_v */);
			printf("reduced generating set:\n");
			gen->Latex(stdout);
			printf("\n");
			}
		}

	// write_generators(gen, stdout);
	{ 
		BYTE fname[1000];
		FILE *fp;
		
		sprintf(fname, "gen_%s.txt", g_label);
		fp = fopen(fname, "w");
		write_generators(gen, fp);
		fclose(fp);
		sprintf(fname, "gen_%s.g", g_label);
		fp = fopen(fname, "w");
		gen->fprint_GAP(fp);
		fclose(fp);
	}
	
	if (f_strong_generators) {
		labG.strong_generating_set(gen, TRUE /* f_v */);
		printf("strong generating set:\n");
		gen->Latex(stdout);
		printf("\n");
		write_generators(gen, stdout);
		}
	

#if 0
	sprintf(gen_G_fname, "labra_%s.dsc", g_label);
	labG.save(gen_G_fname);
#endif

	printf("t = %ld k = %ld\n", t, k);
	if (f_strong_generators)
		printf("f_strong_generators\n");
	if (f_orderly_generation)
		printf("f_orderly_generation\n");
	if (f_TDO)
		printf("f_TDO\n");
	if (f_extension_construction)
		printf("f_extension_construction\n");
	if (f_canonical_representatives)
		printf("f_canonical_representatives\n");
	if (f_k_k2_design)
		printf("$\\{k,k2\\} = \\{%ld,%ld\\}$ design\n", k, k2);
	fflush(stdout);


	if (f_orderly_generation) {
		if (f_k_k2_design) {
			error("error: use k2 not yet implemented with orderly generation !\n"
				"use non-orderly generation instead !");
			return;
			}
#ifdef SYM_GEO
		calc_kramer_mesner_matrix_by_orderly_generation(gen, &go, deg, &labG, &TG, 
			g_label, g_label_tex, t, k, f_TDO);
			// declaration in geo_data.h !
#else
		error("calc_kramer_mesner_matrix_by_orderly_generation() not available");
#endif
		}
	else {
		calc_kramer_mesner_matrix(&labG, gen, &go, deg, g_label, g_label_tex, 
			t, k, f_k_k2_design, k2, 
			f_TDO, 
			f_extension_construction, 
			f_canonical_representatives);
		}
	
}


#if TEXDOCU
INT calc_kramer_mesner_matrix(LABRA_OP labG, VECTOR_OP gen, SYM_OP go, INT deg, 
	BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, 
	INT f_k_k2_design, INT k2, 
	INT f_TDO, 
	INT f_extension_construction, 
	INT f_canonical_representatives)
#else
This is the internal function for computing Kramer-Mesner matrices 
via Leiterspiel. The parameters are stored in an LADDER\_INFO structure 
of the file li.C.
Finally, a file 
$KM\_\langle group \rangle\_t\langle t \rangle\_k \langle k \rangle.txt$ 
is written containing the result.
This file contains 
\begin{enumerate}
\item
the generators in LaTeX, 
\item
the Kramer Mesner matrix in ASCII, 
\item
the generators in coded form, 
\item
a list of set representatives in coded form, 
\item
the corresponding stabilizer orders (coded)
\item
and a list of neighboured KM-matrices, i.e. $M_{i,i+1}$ 
for $0 \le i <k$ (coded), 
\item
(optionally) the TDO version of the Kramer Mesner matrix in ASCII (see below).
\end{enumerate}
In later stages of the KM-program (not this routine), the KM-file 
can be extended by
\begin{enumerate}
\item
solution vectors 
\item
intersection matrices called Ik2 (coded), 
\item
intersection matrices called Ikk (coded), 
\end{enumerate}
The flag f\_canonical\_representatives allows to switch over 
to lexicographically-least representatives of the sets of each orbit 
(so-called canonical representative). This may take some time 
in the end of the computation but the numbers may look nicer 
afterwards. Note that this does not imply that the orbits 
are sorted lexicographically: the ordering of the orbits 
is imposed by the algorithm Leiterspiel which is used for their construction. 
But apart from ordering the orbits, the canonical representatives 
can be used to compare results of different sources: 
if the permutation groups are equal, there should be a 1-1 mapping 
between the orbit representatives of different programs.
The flag f\_TDO determines if a TDO-decomposition of the KM-matrix 
should be computed in the end (and appended to the KM-file). 
Thid TDO decomposition may help to find large entries in the matrix 
since the algorithms reorders rows and colums lexicographically. 
If the flag f\_k\_k2\_design is set, construction of 
block designs with different values of $\lambda$ is supported 
(here $k$ and $k2$). Namely, the combined matrix 
$(M_{t,k}|M_{t,k2})$ is computed. Clearly, $k2$ has to be specified 
in this case.
The flag f\_extension\_construction is no longer in use.
We give an example of a KM-file. Starting with the group 
PGGL(2,32) one gets for $t=5$, $k=7$ and $\lambda=42$ 317 designs. 
The file KM\_PGGL\_2\_32\_t5\_k7.txt is (shortened):\par
{\mytt
\input ../REFERENCE/KM/KM_PGGL_2_32_t5_k7.tex
}
#endif
{
	LADDER_INFO *li;
	INT max_k;
	INT f_calc_stab_go = TRUE; // dv->f_the_stab;
	INT f_calc_reps = TRUE; // dv->f_the_sets;
	INT f_save_set_orbits = TRUE;
	INT v;
	MATRIX_OB Mtk, Mtk2, Mtkm1, Mtm1km1, Mtk3, Mttp1;
	INT i, j, m1, n1, m2, n2, m3, n3, a, mm, nn, kk;
	VECTOR_OB G_gen, RR, R, MM, stab_go;

	gen->copy(&G_gen);
	if (f_k_k2_design)
		max_k = MAXIMUM(k, k2);
	else
		max_k = k;
	// printf("max_k = %ld\n", max_k); fflush(stdout);
	li = init_ladder_info(gen, go, deg, g_label, g_label_tex, 
		t, max_k, 0 /* dv->lambda */, 
		TRUE /* f_verbose */, 
		NIL /* dv->generators_fname */ );
	
	li_init_file_names(li);
		/* compute filename of kramer mesner matrix file into li->txt_out 
		 * and writes this file into "group_label" */

	li_message(li);
		/* hi I am DCC */
	
	li_init_transversals(li);

	li_leiterspiel(li);

	v = li->deg;

	printf("computing KM-matrix for t=%ld k=%ld:\n", li->t, li->k);
	fflush(stdout);
	
	li_Mtk(li, &Mtk, t, k);

	printf("the matrix is of size %ld x %ld\n", Mtk.s_hi(), Mtk.s_li());
	fflush(stdout);
	
	if (f_k_k2_design) {
		printf("computing KM-matrix for t=%ld k=%ld:\n", t, k2);
		fflush(stdout);
	
		li_Mtk(li, &Mtk2, t, k2);
	
		m1 = Mtk.s_hi();
		n1 = Mtk.s_li();
		m2 = Mtk2.s_hi();
		n2 = Mtk2.s_li();
		if (m1 != m2)
			return error("m1 != m2");
		m3 = m1;
		n3 = n1 + n2;
		Mtk3.m_ilih(n3, m3);
		for (i = 0; i < m1; i++) {
			for (j = 0; j < n1; j++) {
				a = Mtk.s_iji(i, j);
				Mtk3.m_iji(i, j, a);
				}
			}
		for (i = 0; i < m2; i++) {
			for (j = 0; j < n2; j++) {
				a = Mtk2.s_iji(i, j);
				Mtk3.m_iji(i, n1 + j, a);
				}
			}
		Mtk3.swap(&Mtk);
		li_print_M_asc(li, &Mtk, t, k, TRUE /* f_k2 */, k2);
		}
	else if (f_extension_construction) {
		printf("extension-construction:\n");
		
		printf("computing KM-matrix for t=%ld k-1=%ld:\n", li->t, li->k - 1);
		fflush(stdout);
		li_Mtk(li, &Mtkm1, li->t, li->k - 1);
		mm = Mtkm1.s_hi();
		nn = Mtkm1.s_li();
		printf("%ld x %ld matrix\n", mm, nn);
		// Mtkm1.Print();
	
		printf("computing KM-matrix for t-1=%ld k-1=%ld:\n", li->t - 1, li->k - 1);
		fflush(stdout);
		li_Mtk(li, &Mtm1km1, li->t - 1, li->k - 1);
		mm = Mtm1km1.s_hi();
		nn = Mtm1km1.s_li();
		printf("%ld x %ld matrix\n", mm, nn);
		// Mtm1km1.Print();
	
		m1 = Mtk.s_hi();
		n1 = Mtk.s_li();
		n2 = Mtkm1.s_li();
		m2 = Mtm1km1.s_hi();
		m3 = m1 + m2;
		n3 = n1 + n2;
		Mtk3.m_ilih_n(n3, m3);
		for (i = 0; i < m1; i++) {
			for (j = 0; j < n1; j++) {
				a = Mtk.s_iji(i, j);
				Mtk3.m_iji(i, j, a);
				}
			}
		for (i = 0; i < m1; i++) {
			for (j = 0; j < n2; j++) {
				a = Mtkm1.s_iji(i, j);
				Mtk3.m_iji(i, n1 + j, a);
				}
			}
		for (i = 0; i < m2; i++) {
			for (j = 0; j < n2; j++) {
				a = Mtm1km1.s_iji(i, j);
				Mtk3.m_iji(m1 + i, n1 + j, a);
				}
			}
		Mtk3.swap(&Mtk);
		
		li_print_M_asc(li, &Mtk, t, k, FALSE /* f_k2 */, 0);
		}
	else {
		li_print_M_asc(li, &Mtk, t, k, FALSE /* f_k2 */, 0);
		}
	
#ifdef SYM_GEO
	km_write_ascii_G_gen(li->txt_out, &G_gen);
#endif

	// the set representatives:
	RR.m_il(k + 1);
	for (kk = 0; kk <= k; kk++) {
		li_set_reps(li, &R, kk);
		R.swap(RR.s_i(kk));
		}

	if (f_canonical_representatives) {
		dc_calc_canonical_representatives(labG, &RR, &stab_go, TRUE /* f_v */);
		}
	else {
		dc_calc_stab_go((DCY_OP) li->dc->s_s(), k, &stab_go);
		}

#ifdef SYM_GEO
	km_write_ascii_representatives(li->txt_out, &RR);
	km_write_ascii_stab_go(li->txt_out, &stab_go);
	km_write_stab_go_k_sets(li->txt_out, &stab_go, k);
#endif

	MM.m_il(k);
	for (kk = 0; kk < k; kk++) {
		printf("computing KM_%ld,%ld ", kk, kk + 1);
		fflush(stdout);
		dc_Mttp1((DCY_OP) li->dc->s_s(), kk, &Mttp1);
		Mttp1.swap((MATRIX_OP) MM.s_i(kk));
		printf("finished\n"); fflush(stdout);
		}
#ifdef SYM_GEO
	km_write_ascii_KM_matrices(li->txt_out, &MM);
#endif
	
#if 0
	printf("writing KM-matrix:\n");
	fflush(stdout);
	sprintf(KM_fname, "KM_%s_t%ld_k%ld.dsc", g_label, li->t, li->k);
	Mtk.save(KM_fname);
#endif
	
#ifdef SYM_GEO
	if (f_TDO) {
		printf("computing TDO decomposition\n");
		fflush(stdout);
		km_compute_TDO_decomposition(g_label, li->t, li->k, &Mtk);
		}
#endif
	
	
	if (f_calc_stab_go && f_calc_reps) {
		VECTOR_OB S_go, R;

		printf("computing stabilizers/representatives "
			"for i=%ld to k=%ld:\n", (INT)1, li->k);
		fflush(stdout);
		for (i = 1; i <= k; i++) {
			li_stab_go(li, &S_go, i);
			li_set_reps(li, &R, i);
			li_print_dc_info(li, &S_go, &R, i, stdout, FALSE /* f_verbose */);
			if (f_save_set_orbits) {
				DCY_OP dcy = (DCY_OP) li->dc->s_s();
				DCY_OP dc2;
				INT t1;
			
				if (i == 0)
					t1 = 0;
				else if (i == 1)
					t1 = 1;
				else
					t1 = 2 * i - 1;
				dc2 = dcy + t1;
				}
			}
		}


	free_ladder_info(li);

	return OK;
}

#if TEXDOCU
INT design_orbits(MATRIX_OP X, VECTOR_OP orbits, INT f_complement)
#else
This function converts a solution column $X$ to a list of orbit indices 
stored into orbits. 
In fact, the solution $X$ is given as a matrix of just one 
column. If f\_complement id true, exactly the non-selected orbit 
indices are returned. They form the orbits of the complementary design.
#endif
{
	INT i, l, j, nb = 0;

	l = X->s_hi();
	for (i = 0; i < l; i++) {
		if (f_complement) {
			if (X->s_iji(i, 0) == 1)
				continue;
			}
		else {
			if (X->s_iji(i, 0) == 0)
				continue;
			}
		nb++;
		}
	orbits->m_il(nb);
	j = 0;
	for (i = 0; i < l; i++) {
		if (f_complement) {
			if (X->s_iji(i, 0) == 1)
				continue;
			}
		else {
			if (X->s_iji(i, 0) == 0)
				continue;
			}
		orbits->m_ii(j, i);
		j++;
		}
	if (j != nb)
		return error("design_orbits() j != nb");
	return OK;
}

#if TEXDOCU
INT design_orbits_vector(VECTOR_OP X, VECTOR_OP orbits, INT f_complement)
#else
This function converts a solution vector $X$ to a list of orbit indices 
stored into orbits. 
If f\_complement id true, exactly the non-selected orbit 
indices are returned. 
They form the orbits of the complementary design.
#endif
{
	INT i, l, j, nb = 0;

	l = X->s_li();
	for (i = 0; i < l; i++) {
		if (f_complement) {
			if (X->s_ii(i) == 1)
				continue;
			}
		else {
			if (X->s_ii(i) == 0)
				continue;
			}
		nb++;
		}
	orbits->m_il(nb);
	j = 0;
	for (i = 0; i < l; i++) {
		if (f_complement) {
			if (X->s_ii(i) == 1)
				continue;
			}
		else {
			if (X->s_ii(i) == 0)
				continue;
			}
		orbits->m_ii(j, i);
		j++;
		}
	if (j != nb)
		return error("design_orbits_vector() j != nb");
	return OK;
}

#if TEXDOCU
INT dc_calc_orbit_length(SYM_OP go, VECTOR_OP stab_go, VECTOR_OP orbit_length)
#else
Computes the vector orbit\_length out of stab\_go. 
Each entry is computed as the index of the corresponding entry in 
stab\_go in the group-order go.
#endif
{
	VECTOR_OP p1, p2;
	INT i, j, k, l;
	SYM_OP pa, pb;
	
	k = stab_go->s_li();
	orbit_length->m_il(k);
	for (i = 0; i < k; i++) {
		p1 = (VECTOR_OP) stab_go->s_i(i);
		p2 = (VECTOR_OP) orbit_length->s_i(i);
		l = p1->s_li();
		p2->m_il(l);
		for (j = 0; j < l; j++) {
			pa = p1->s_i(j);
			pb = p2->s_i(j);
			go->ganzdiv(pa, pb);
			}
		}
	return OK;
}

#if TEXDOCU
INT dc_calc_stab_go_and_K(DCY_OP dc, INT k_max, VECTOR_OP stab_go, VECTOR_OP K)
#else
This function is used to get the results of the Leiterspiel out of 
the DCY data-structure. dc is actually a pointer to a lot of 
DCY structures, namely the complete ladder. The representatives 
of all sets corresponding to $i$ orbits $i \le k$ are taken 
out of the data-structure. 
They correspond to partitions of exactly 2 blocks, namely $[i,n-i]$ 
when $n$ is the degree.
The intermediate steps (corresponding to partitions 
$[i,1,n-i-1]$) are unimportant.
This function computes stab\_go and K. 
stab\_go contains the stabilizer orders of all orbit representatives 
for $i$-orbits, $i \le k\_max$. 
$K$ contains the size of the representative, i.e. the value $i$ itself.
#endif
{
	VECTOR_OP Ad, p_stab_go;
	LABRA_OB labra_A;
	INT up_to_step, i, k, s, dc_no, nb_d, nb_D;
	
	up_to_step = dc_k_to_step(k_max);
	nb_D = dc_calc_nb_d(dc, up_to_step);
	stab_go->m_il(k_max + 1);
	K->m_il(nb_D);
	dc_no = 0;
	for (k = 0; k <= k_max; k++) {
		s = dc_k_to_step(k);
		nb_d = dc[s].s_D()->s_li();
		p_stab_go = (VECTOR_OP) stab_go->s_i(k);
		p_stab_go->m_il(nb_d);
		for (i = 0; i < nb_d; i++) {
			Ad = dc[s].s_Ad_i(i);
			reduce_generators_labra(Ad, p_stab_go->s_i(i), 
				FALSE /* f_verbose */, &labra_A);
			K->m_ii(dc_no, k);
			dc_no++;
			}
		}
	return OK;
}

#if TEXDOCU
INT dc_RR_to_K(VECTOR_OP RR, VECTOR_OP K)
#else
This function computes the K-vector out of RR.
K is a vector of integers containing the 
size of the set for each representative (so, a sequence 
like $(0,1,2,2,3,...,3,4,....,4)$ if the number of orbits 
on 0-,1-,2-..sets is 1,1,2,.. ). 
RR is a vector of vectors (of vectors) holding the
representatives. 
#endif
{
	INT nb_d, nb_d1, i, j, k, l;
	VECTOR_OP pR;
	
	k = RR->s_li();
	nb_d = 0;
	for (i = 0; i < k; i++) {
		pR = (VECTOR_OP) RR->s_i(i);
		l = pR->s_li();
		nb_d += l;
		}
	K->m_il(nb_d);
	nb_d1 = 0;
	for (i = 0; i < k; i++) {
		pR = (VECTOR_OP) RR->s_i(i);
		l = pR->s_li();
		for (j = 0; j < l; j++) {
			K->m_ii(nb_d1++, i);
			}
		}
	return OK;
}

#if TEXDOCU
INT dc_step_to_k(INT step)
#else
computes k out of step.
#endif
{
	if (step == 0)
		return 0;
	if (step == 1)
		return 1;
	return (step + 1) / 2;
}

#if TEXDOCU
INT dc_k_to_step(INT k)
#else
computes the step number corresponding to the partition $[k,n-k]$.
#endif
{
	if (k == 0)
		return 0;
	else
		return 2 * k - 1;
}

#if TEXDOCU
INT dc_print_k_set(VECTOR_OP R, SYM_OP stab_go, SYM_OP go)
#endif
{
	SYM_OB go1;
	INT l, j;
	BYTE s1[256];
	BYTE s2[256];
	
	l = R->s_li();
	printf(" { ");
	for (j = 0; j < l; j++) {
		printf("%ld", R->s_ii(j));
		if (j < l - 1)
			printf(", ");
		}
	printf(" }");
	go->ganzdiv(stab_go, &go1);
	s1[0] = 0;
	s2[0] = 0;
	stab_go->sprint(s1);
	go1.sprint(s2);
	printf("_%s_%s\n", s1, s2);
	return OK;
}

#if TEXDOCU
INT dc_get_k_set(PERMUTATION_OP d, VECTOR_OP R, INT k, INT f_v)
#else
This function reads out the representative from a given permutation d. 
R becomes a vector of integers holding all elements of the set.
#endif
{
	INT j, j0, j1, n, kk;

	n = d->s_li();
	j0 = n - k;
	R->m_il(k);
	for (j = 0; j < k; j++) {
		j1 = j0 + j;
		kk = d->s_ii(j1) - 1;
		R->m_ii(j, kk);
		}
	if (f_v) {
		printf("%ld-set: { ", k);
		for (j = 0; j < k; j++) {
			printf("%ld", R->s_ii(j));
			if (j < k - 1)
				printf(", ");
			}
		printf(" }\n");
		// R->println();
		}
	return OK;
}

#if TEXDOCU
INT apply_perm_set(VECTOR_OP R, PERMUTATION_OP p, INT f_reorder)
#endif
{
	INT i, a, b, l;
	
	l = R->s_li();
	for (i = 0; i < l; i++) {
		a = R->s_ii(i);
		b = p->s_ii(a) - 1;
		R->m_ii(i, b);
		}
	if (f_reorder) {
		R->quicksort(R->s_li(), TRUE);
		}
	return OK;
}

#if TEXDOCU
INT design_calc_blocks(VECTOR_OP gen, VECTOR_OP RR, INT k, MATRIX_OP X, 
	VECTOR_OP O, VECTOR_OP O_first, VECTOR_OP O_len, INT f_v, INT f_vv)
#else
RR is a vector of vectors holding the representing sets 
for all $i \le k$ (of length $k+1$). 
At the beginning, the vector at position $k$ is selected.
From the set representatives and the selection vector (column) X 
and the generators gen of the group, the design is built. 
O contains a list of blocks of the design.
The function  dc\_calc\_set\_orbit (see below) is used to compute all sets 
of one particular orbit. These orbits are then concatenated to 
form O. O\_first and O\_len are integer vectors 
holding the positions of the start of the blocks of one particular orbit 
and the length of that orbit. 
#endif
{
	VECTOR_OB O1, R0;
	VECTOR_OP R, R1;
	INT a, b, i, j, l, o_first, o_len;
	
	O->m_il(0);
	R = (VECTOR_OP) RR->s_i(k);
	l = R->s_li(); // number of k-orbits 
	if (l != X->s_hi())
		return error("design_calc_blocks() l != X->s_hi()");
	O_first->m_il(l);
	O_len->m_il(l);
	o_first = 0;
	o_len = 0;
	b = 0;
	for (i = 0; i < l; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		R1 = (VECTOR_OP) R->s_i(i);
		R0.m_il(k);
		for (j = 0; j < k; j++) {
			a = R1->s_ii(j);
			R0.m_ii(j, a - 1);
			}
		dc_calc_set_orbit(gen, &R0, &O1, f_v, f_vv);
		o_len = O1.s_li();
		O->realloc_z(b + o_len);
		for (j = 0; j < o_len; j++) {
			O1.s_i(j)->swap(O->s_i(b + j));
			}
		b += o_len;
		O_first->m_ii(i, o_first);
		O_len->m_ii(i, o_len);
		o_first += o_len;
		}
	return OK;
}

#if TEXDOCU
INT dc_calc_set_orbit(VECTOR_OP gen, VECTOR_OP R, 
	VECTOR_OP O, INT f_v, INT f_vv)
#else
Computes the orbit $O$ of the set $R$ under the group generated by $gen$. 
$O$ becomes a (sorted) vector of vectors over integer.
Useful for determination of the complete block set of a design 
which was found via the Kramer-Mesner method.
Clearly, this routine is applicable only for reasonably {\em small} orbits.
#endif
{
	INTEGER_OB O_len, Q_len, int_ob;
	VECTOR_OB Q, R1;
	INT i, j, k, l, n, ii, a, b, idx, f_found;
	PERMUTATION_OP p;
	VECTOR_OP R0;

	if (f_v) {
		printf("dc_calc_set_orbit() calculating orbit of ");
		R->println();
		fflush(stdout);
		}
	l = gen->s_li();
	n = R->s_li();
	O_len.m_i(1);
	O->m_il(VECTOR_OVERSIZE);
	R->copy(&R1);
	R1.quicksort(R1.s_li(), TRUE /* f_ascending */);
	R1.copy((VECTOR_OP) O->s_i(0));
	
	Q_len.m_i(1);
	Q.m_il(VECTOR_OVERSIZE);
	Q.m_ii(0, 0);
	while (Q_len.s_i() > 0) {

		// get the index of the next unprocessed orbit element:
		// shrink the queue
		i = Q.s_ii(0);
		for (ii = 1; ii < Q_len.s_i(); ii++) {
			Q.m_ii(ii - 1, Q.s_ii(ii));
			}
		Q_len.dec();

		for (j = 0; j < l; j++) {
			
			// get the pointer to R0 each time new 
			// because O might have changed in between !
			R0 = (VECTOR_OP) O->s_i(i);
			p = (PERMUTATION_OP) gen->s_i(j);
			for (k = 0; k < n; k++) {
				a = R0->s_ii(k);
				b = p->s_ii(a) - 1;
				R1.m_ii(k, b);
				}
			R1.quicksort(n, TRUE /* f_ascending */);
			/* R1.println();
			fflush(stdout); */
			O->search(O_len.s_i(), TRUE /* f_ascending */, &R1, &idx, &f_found);
			if (f_found)
				continue;
			if (f_vv) {
				printf("%ld ", O_len.s_i());
				/* R1.println(); */
				if (O_len.s_i() % 10 == 0)
					printf("\n");
				fflush(stdout);
				}
			O->insert_at(&O_len, idx, &R1);
			if (i >= idx)
				i++;
			for (k = 0; k < Q_len.s_i(); k++)
				if (Q.s_ii(k) >= idx)
					Q.s_i(k)->inc();
			int_ob.m_i(idx);
			Q.insert_at(&Q_len, Q_len.s_i(), &int_ob);

			
			} // next j
		
		} // while Q not empty
	O->realloc_z(O_len.s_i());
	if (f_vv) {
		printf("\n");
		}
	if (f_v) {
		printf("dc_calc_set_orbit() found an orbit of length %ld\n", O_len.s_i());
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT dc_calc_go(DCY_OP dc0, SYM_OP go)
#else
Computes the group order of the stabilizer of the first double coset at step 0. 
Normally, this is the group itself (the group of the ladder and of the KM-system).
#endif
{
	VECTOR_OP Ad;
	LABRA_OB labra_A1;
	INT type;
	void *data;

	type = 0;
	data = NIL;
	Ad = dc0->s_Ad_i(0);
	reduce_generators_labra(Ad, go, FALSE /* f_verbose */, &labra_A1);
	return OK;
}

#if TEXDOCU
INT dc_calc_stab_go(DCY_OP dc, INT k_max, VECTOR_OP stab_go)
#else
Computes a vector stab\_go of all stabilizer orders of double cosets.
Only the double cosets for $j$-sets are taken ($0 \le j \le k\_max$). 
Some steps in the ladder correspond to partitions which do not belong 
to sets. They are left out.
#endif
{
	INT i, j, s, nb_d, up_to_step;
	VECTOR_OP Ad, p_stab_go;
	LABRA_OB labA;
	SYM_OB ago;
	
	up_to_step = dc_k_to_step(k_max);
	nb_d = dc_calc_nb_d(dc, up_to_step);
	stab_go->m_il(k_max + 1);
	for (i = 0; i <= k_max; i++) {
		s = dc_k_to_step(i);
		nb_d = dc[s].s_D()->s_li();
		p_stab_go = (VECTOR_OP) stab_go->s_i(i);
		p_stab_go->m_il(nb_d);
		for (j = 0; j < nb_d; j++) {
			Ad = dc[s].s_Ad_i(j);
			reduce_generators_labra(Ad, &ago, FALSE /* f_verbose */, &labA);
			ago.swap(p_stab_go->s_i(j));
			}
		}
	return OK;
}

#if TEXDOCU
INT dc_orbit_size(DCY_OP dc, SYM_OP go, INT step, INT dc_no, INT f_v)
#else
returns the orbit length (as an integer) of orbit $dcno$ in $step$. 
An error is raised if the orbit length is not an INTEGER 
(but a long integer instead).
#endif
{
	DCY_OP Dc;
	VECTOR_OP Ad;
	SYM_OB stab_go, ol;
	LABRA_OB labra_A1;
	INT type;
	void *data;
	INT oli;
	INT i, l;

	type = 0;
	data = NIL;
	Dc = dc + step;
	Ad = Dc->s_Ad_i(dc_no);
	reduce_generators_labra(Ad, &stab_go, FALSE /* f_verbose */, &labra_A1);
	go->ganzdiv(&stab_go, &ol);
	if (ol.s_obj_k() != INTEGER)
		return error("dc_orbit_size() ol not an INTEGER!");
	oli = ((INTEGER_OP) &ol)->s_i();
	if (f_v) {
		printf("dc_orbit_size() step = %ld nb dcs = %ld\n", 
			Dc->s_step_i(), Dc->s_D()->s_li());
		printf("step %ld dc_no %ld stab_order = ", step, dc_no);
		stab_go.print();
		printf(" orbit_length = %ld\n", oli);
		printf("stab generators:");
		l = Ad->s_li();
		for (i = 0; i < l; i++) {
			Ad->s_i(i)->println();
			}
		}
	return oli;
}

#if TEXDOCU
INT dc_print_dc(DCY_OP dc, SYM_OP go, INT step, INT dc_no)
#endif
{
	return dc_print_dc_(dc, go, step, dc_no, FALSE, NIL, NIL, 0);
}

#if TEXDOCU
INT dc_print_dc_(DCY_OP dc, SYM_OP go, INT step, INT dc_no, 
	INT f_d0, PERMUTATION_OP d0, PERMUTATION_OP d1, INT step2)
#endif
{
	DCY_OP Dc;
	PERMUTATION_OP d;
	VECTOR_OB R;
	SYM_OB stab_go;
	LABRA_OB labra_A1;
	INT type;
	void *data;
	PERMUTATION_OB d1v, q, tmp;
	INT n, i, l, alpha, a, k, k2;

	type = 0;
	data = NIL;
	Dc = dc + step;
	printf("step %ld dc no %ld: ", step, dc_no);
	d = (PERMUTATION_OP) Dc->s_D()->s_i(dc_no);
	/* d->println(); */
	reduce_generators_labra(Dc->s_Ad_i(dc_no), &stab_go, 
		FALSE /* f_verbose */, &labra_A1);
	if (step == 0)
		k = 0;
	else if (step == 1)
		k = 1;
	else
		k = step / 2 + 1;
	dc_get_k_set(d, &R, k, FALSE);
	if (f_d0) {
		d0->mult(d, &q);
		dc_get_k_set(&q, &R, k, FALSE);
		}
	dc_print_k_set(&R, &stab_go, go);
	// printf("stab_go = ");
	// stab_go.println();
	// printf("orbit length = ");
	// go1.println();
	if (f_d0) {
		n = d1->s_li();
		d1->invers(&d1v);
		apply_perm_set(&R, &d1v, TRUE /* f_reorder */);
		alpha = 0;
		l = R.s_li();
		k2 = step2 / 2 + 1;
		for (i = 0; i < l; i++) {
			a = R.s_ii(i);
			if (a >= n - k2)
				alpha++;
			}
		printf("intersection %ld\n", alpha);
		}
	return OK;
}

#if TEXDOCU
INT dc_calc_Ainf_t(DCY_OP dc0, INT up_to_step, 
	MATRIX_OP Ainf_t, INT f_verbose)
#else
This is the OLD version to compute the plesken matrix.
It computes the $A^\wedge$ matrix on the whole and is 
ineffective for large examples. All $M_{i,j}$ 
are computed for all $i < j \le k\_max$ 
where $k\_max$ is dc\_step\_to\_k(up\_to\_step). 
The function dc\_Mtk() is called to compute $M_{i,j}$.
The size of the matrix is $A^\wedge$ determined via dc\_calc\_nb\_d().
#endif
{
	MATRIX_OB A, Ainf;
	INT nb_d0, nb_d;
	INT nb_D, nb_D1, t_sets, k_sets;
	INT i, j, k, t, ii, jj, T;
	DCY_OP L;
	
	nb_d0 = dc_calc_nb_d(dc0, up_to_step);
	if (f_verbose)
		printf("nb_d = %ld\n", nb_d0);
	Ainf_t->m_ilih_n(nb_d0, nb_d0);
	// layer->m_il_n(nb_d0);
	// layer->m_ii(0, 0);
	nb_D = 1;
	for (i = 1; i <= up_to_step; i += 2) {
		L = dc0 + i;
		k = L->s_D()->s_li();
		j = i / 2;
		j++;
		// for (ii = 0; ii < k; ii++) 
	// 		layer->m_ii(nb_D + ii, j);
		nb_D += k;
		}
	
	for (i = 0; i < nb_d0; i++)
		Ainf_t->m_iji(i, i, 1);
	/* Ainf_t->m_iji(0, 1, n2); */
	T = (up_to_step + 1) >> 1;
	nb_D = 0; 
		/* number of dc's before the t-sets */
	for (t = 0; t <= T; t++) {
		nb_d = nb_D;
			/* number of dc's before the k sets */
		for (k = t + 1; k <= T; k++) {
			dc_Mtk(dc0, t, k, &A);
			if (f_verbose) {
				printf("M_%ld,%ld = \n", t, k);
				A.Print();
				}
			if (k == t + 1) {
				t_sets = A.s_hi();
				nb_D1 = nb_D + A.s_hi();
				nb_d += t_sets; /* important ! */
				}
			else {
				if (A.s_hi() != t_sets)
					return error("dc_calc_Ainf_t() A.s_hi() != t_sets");
				}
			k_sets = A.s_li();
			for (ii = 0; ii < t_sets; ii++) {
				for (jj = 0; jj < k_sets; jj++) {
					Ainf_t->m_iji(nb_D + ii, nb_d + jj, A.s_iji(ii, jj));
					}
				}
			nb_d += k_sets;
			}
		nb_D = nb_D1;
		}
	L = dc0 + up_to_step;
	nb_d += L->s_D()->s_li();
	if (nb_d != nb_d0)
		printf("dc_calc_Ainf_t(): nb_d = %ld != nb_d0 = %ld", nb_d, nb_d0);
	Ainf_t->latex_upper_tri(stdout);
	Ainf_t->transpose(&Ainf);
	Ainf.latex_lower_tri(stdout);
#if 0
	nb_d = 1;
	for (i = 1; i < up_to_step; i += 2) {
		t = ((i - 1) >> 1) + 1;
		dc_Mttp1(dc0, t, &A);
		if (f_verbose) {
			printf("M_%ld,%ld = \n", t, t + 1);
			A.Print();
			}
		nb_d1 = nb_d + A.s_hi();
		for (j = 0; j < A.s_hi(); j++) {
			for (k = 0; k < A.s_li(); k++) {
				Ainf_t->m_iji(
				nb_d + j, 
				nb_d1 + k, 
				A.s_iji(j, k));
				}
			}
		nb_d = nb_d1;
		}
	L = dc0 + up_to_step;
	nb_d += L->s_D()->s_li();
	if (nb_d != nb_d0)
		printf("dc_calc_Ainf_t(): "
		"nb_d = %ld != nb_d0 = %ld", 
		nb_d, nb_d0);
#endif
	return OK;
}

#if TEXDOCU
INT dc_dc_no_to_dc_idx(DCY_OP dc0, INT dc_no, INT *step, INT *dc_idx)
#else
Transforms to dc\_no (a global dc-number) to a 
step and dc\_idx pair.
#endif
{
	INT i, nb_d0, nb_d = 1;
	DCY_OP L;
	
	if (dc_no == 0) {
		*step = 0;
		*dc_idx = 0;
		return OK;
		}
	for (i = 1; ; i += 2) {
		L = dc0 + i;
		nb_d0 = nb_d;
		nb_d += L->s_D()->s_li();
		if (dc_no < nb_d) {
			*step = i;
			*dc_idx = dc_no - nb_d0;
			return OK;
			}
		}
}

#if TEXDOCU
INT dc_dc_idx_to_dc_no(DCY_OP dc0, INT step, INT dc_idx, INT *dc_no)
#else
inverse function to dc\_dc\_no\_to\_dc\_idx(). 
Computes dc\_no out of step and dc\_idx.
#endif
{
	INT i, nb_d = 1;
	DCY_OP L;
	
	if (step == 0) {
		*dc_no = dc_idx;
		return OK;
		}
	for (i = 1; i < step; i += 2) {
		L = dc0 + i;
		nb_d += L->s_D()->s_li();
		}
	*dc_no = nb_d + dc_idx;
	return OK;
}

#if TEXDOCU
INT dc_calc_nb_d(DCY_OP dc0, INT up_to_step)
#else
Computes the number of all orbits of the group on $i$-sets, 
$i \le k\_max$, k\_max determined via up\_to\_step.
#endif
{
	INT i, nb_d = 1;
	DCY_OP L;
	
	if (!ODD(up_to_step)) {
		printf("dc_calc_nb_d(): !ODD(up_to_step)");
		return -1;
		}
	for (i = 1; i <= up_to_step; i += 2) {
		L = dc0 + i;
		nb_d += L->s_D()->s_li();
		}
	printf("dc_calc_nb_d() found %ld double cosets\n", nb_d);
	return nb_d;
}

#if TEXDOCU
INT dc_calc_nb_d_via_MM(VECTOR_OP MM, INT k_max)
#else
computes the total number of double cosets corresponding 
to $i$-orbits out of the vector $MM$ of neighboured KM-matrices. 
The size of the matrices is added.
#endif
{
	MATRIX_OP pM;
	INT i, nb_d = 0;
	
	for (i = 0; i < k_max; i++) {
		pM = (MATRIX_OP) MM->s_i(i);
		nb_d += pM->s_hi();
		}
	pM = (MATRIX_OP) MM->s_i(k_max - 1);
	nb_d += pM->s_li();
	printf("dc_calc_nb_d_via_MM() found %ld double cosets\n", nb_d);
	return nb_d;
}

#if TEXDOCU
INT dc_Mttp1(DCY_OP dc0, INT t, MATRIX_OP M)
#else
Computes the neighboured KM-matrix $M_{t,t+1}$ out of the results 
of the Leiterspiel.  
#endif
{
	INT i, j, k, i1, nb_Dim1, nb_Di, nb_Dip1;
	INT idx_Di;
	INT erg = OK;
	DCY_OP Lim1, Li, Lip1;
	
	if (t == 0) {
		DCY_OP L0, L1;
		INT nb_d0, nb_d1;

		L0 = dc0;
		L1 = dc0 + 1;
		nb_d0 = L0->s_D()->s_li();
		nb_d1 = L1->s_D()->s_li();
		erg += M->m_ilih_n(nb_d1, nb_d0);
		for (i = 0; i < nb_d0; i++) {
			for (k = 0; k < L1->s_T()->s_li(); k++) {
				idx_Di = L1->s_TDidx_iji(k, i);
				/* no upstep here ! */
				M->s_ij(i, idx_Di)->inc();
				}
			}
		return erg;
		}
	i1 = t * 2;
	Lim1 = dc0 + (i1 - 1);
	Li = dc0 + i1;
	Lip1 = dc0 + (i1 + 1);
	if (i1 != 1 && i1 != 2 && Lim1->s_fDown_i())
		return error("dc_Mttp1()|i1 != 1 && i1 != 2 && Lim1->s_fDown_i()");
	nb_Dim1 = Lim1->s_D()->s_li();
	nb_Di = Li->s_D()->s_li();
	nb_Dip1 = Lip1->s_D()->s_li();
	erg += M->m_ilih_n(nb_Dip1, nb_Dim1);
	for (i = 0; i < nb_Dim1; i++) {
		for (k = 0; k < Li->s_T()->s_li(); k++) {
			idx_Di = Li->s_TDidx_iji(k, i);
			j = Lip1->s_TDidx_iji(0, idx_Di);
			M->s_ij(i, j)->inc();
			}
		}
	return erg;
}

#if TEXDOCU
INT dc_Mtk(DCY_OP dc0, INT t, INT k, MATRIX_OP M)
#else
Computes $M_{t,k}$ out of all $M_{i,i+1}$ 
for $t \le i < k$ via the recursion formula 
$M_{t,k} = M_{t,l} \cdot M_{l,k}$
#endif
{
	MATRIX_OB M1, M2;
	INT erg = OK;
	INT t1, i, j, s;
	SYM_OP h;
	SYM_OB h1, h2, h3, s_ob;
	
	printf("dc_Mtk() t = %ld k = %ld\n", t, k);
#if 0
	if (t == 0)
		return error("dc_Mtk() t == 0");
#endif
	if (k < t) 
		return error("dc_Mtk(): k < t !");
	if (k == t) {
		DCY_OP L;
		INT nb_d;

		if (t == 1)
			L = dc0 + 1;
		else
			L = dc0 + 2 * t - 1;
		nb_d = L->s_D()->s_li();
		erg += M->m_ilih_n(nb_d, nb_d);
		M->one();
		return OK;
		}
	
	erg += dc_Mttp1(dc0, t, M);
		/* M := M_t,t+1 */
	for (t1 = t + 1; t1 < k; t1++) {
		/* now M = M_t,t1 */
		erg += dc_Mttp1(dc0, t1, &M1);
			/* M1 := M_t1,t1+1 */
		erg += M->mult(&M1, &M2);
		
#if 0
		printf("M = ");
		M->Print();
		printf("M1 = ");
		M1.Print();
		printf("M2 = ");
		M2.Print();
#endif
		
		M2.copy(M);
			/* M := (k - t) atop 
			        (k - t1) * 
			                M_t,k   (k = t1 + 1)

			 *    = (t1 + 1 - t) atop 
			        (t1 + 1 - t1) * 
			                M_t,t1+1 

			 *    = (t1 + 1 - t) * M_t,t1+1  */
		s = t1 + 1 - t;
		for (i = 0; i < M->s_hi(); i++) {
			for (j = 0; j < M->s_li(); j++) {
				h = M->s_ij(i, j);
				s_ob.m_i_i(s);
				h->ganzdiv(&s_ob, &h1);
				// h1 = h / s;
				
				h1.mult(&s_ob, &h2);
				h2.addinvers_apply();
				h2.add(h, &h3);
				if (!h3.nullp()) {
				// if (h1 * s != h) {
					printf("t = %ld\n", t);
					printf("k = %ld\n", k);
					// printf("h = %ld\n", h);
					// printf("s = %ld\n", s);
					// printf("h1 = %ld\n", h1);
					return error("dc_Mtk()|h1 * s != h");
					}
				h1.copy(M->s_ij(i, j));
				// M->m_iji(i, j, h1);
				}
			}
		}
	return erg;
}

#if TEXDOCU
INT km_normalizer_action_on_orbits(LABRA_OP labra_G, 
	VECTOR_OP Reps, PERMUTATION_OP p, PERMUTATION_OP q, INT f_v)
#else
Computes the action of the permutation $p$ (which normalizes $G$) 
on the representatives of orbits. The result is written into $q$.
Reps must contain canonical orbit representatives, 
otherwise this function does not work.
#endif
{
	INT *theX = NIL;
	MATRIX_OB TG;
	VECTOR_OP pR;
	VECTOR_OB R1;
	LABRA_OB labra_A;
	SYM_OB ago;
	INT v, i, ii, j, l, ll, a, b, idx;
	PERMUTATION_OB transporter;
	
	l = Reps->s_li();
	v = p->s_li();
	if (v != labra_G->s_degree_i())
		return error("km_normalizer_action_on_orbits() wrong degree of permutation");
	labra_G->calc_transversal_matrix(&TG, FALSE);
	theX = (INT *) my_malloc(v * sizeof(INT), "km_normalizer_action_on_orbits");
	q->m_il(l);
	for (i = 0; i < l; i++) {
		pR = (VECTOR_OP) Reps->s_i(i);
		ll = pR->s_li();
		R1.m_il(ll);
		for (j = 0; j < ll; j++) {
			a = pR->s_ii(j) - 1;
			b = p->s_ii(a);
			R1.m_ii(j, b);
			}
		R1.quicksort(R1.s_li(), TRUE);
		for (ii = 0; ii < ll; ii++) {
			a = R1.s_ii(ii);
			theX[ii] = a - 1;
			}
#ifdef SYM_GEO
		geo_canonicize_set(v, ll, theX, 
			labra_G, &TG, FALSE /* f_v */, FALSE /* f_vv */, 
			TRUE /* f_get_aut_group */, &labra_A, &ago, &transporter);
#else
		ago.m_i_i(-1);
#endif
		for (ii = 0; ii < ll; ii++) {
			a = theX[ii] + 1;
			R1.m_ii(ii, a);
			}
		idx = Reps->search_linear(&R1);
		if (idx < 0) {
			printf("km_normalizer_action_on_orbits() set not found\n R1 = ");
			R1.println();
			return error("km_normalizer_action_on_orbits() not normalizing");
			}
		q->m_ii(i, idx + 1);
		}
	my_free(theX);
	return OK;
}

#if TEXDOCU
INT dc_calc_canonical_representatives(LABRA_OP labra_G, 
	VECTOR_OP RR, VECTOR_OP stab_go, INT f_v)
#else
Computes the canonical set-representatives via the function 
geo\_canonicize\_set() from geo\_canon.C.
Additionally, the stabilizer orders are stored in the vector stab\_go.
#endif
{
	INT a, i, j, k, l, ii, v;
	INT *theX = NIL;
	MATRIX_OB TG;
	// SYM_OB go;
	LABRA_OB labra_A;
	SYM_OB ago;
	VECTOR_OP pR, ppR;
	VECTOR_OP p_stab_go;
	PERMUTATION_OB transporter;
	
	v = labra_G->s_degree_i();
	theX = (INT *) my_malloc(v * sizeof(INT), "dc_calc_canonical_representatives");
	labra_G->calc_transversal_matrix(&TG, FALSE);
	// labra_G->group_order(&go);
	
	k = RR->s_li() - 1;
	stab_go->m_il(k + 1);
	for (i = 0; i <= k; i++) {
		p_stab_go = (VECTOR_OP) stab_go->s_i(i);
		if (f_v) {
			printf("%ld-sets:\n", i); fflush(stdout);
			}
		pR = (VECTOR_OP) RR->s_i(i);
		// pR->Print();
		l = pR->s_li();
		p_stab_go->m_il(l);
		for (j = 0; j < l; j++) {
			if (f_v) {
				printf("no %ld (of %ld): ", j + 1, l);
				}
			ppR = (VECTOR_OP) pR->s_i(j);
			if (f_v) {
				printf(" { ");
				}
			for (ii = 0; ii < i; ii++) {
				a = ppR->s_ii(ii);
				theX[ii] = a - 1;
				if (f_v)
					printf("%ld ", a);
				}
			if (f_v) {
				printf(" } ");
				printf("canonical: ");
				}
#ifdef SYM_GEO
			geo_canonicize_set(v, i, theX, 
				labra_G, &TG, FALSE /* f_v */, FALSE /* f_vv */, 
				TRUE /* f_get_aut_group */, &labra_A, &ago, &transporter);
#else
			ago.m_i_i(-1);
#endif
			// labra_A.group_order(&ago);
			printf(" { ");
			for (ii = 0; ii < i; ii++) {
				a = theX[ii] + 1;
				if (f_v) {
					printf("%ld ", a);
					}
				ppR->m_ii(ii, a);
				}
			if (f_v) {
				printf(" } ");
				printf("stab order ");
				ago.println();
				}
			ago.swap(p_stab_go->s_i(j));
			fflush(stdout);
			
			} // next j
		} // next i
	my_free(theX);
	return OK;
}

#if TEXDOCU
INT fuse_orbits(LABRA_OP labra_G, 
	VECTOR_OP Reps, VECTOR_OP new_rep_idx, VECTOR_OP new_reps, INT f_v)
#else
Fuses the orbits of a group $U$ under the group $G$. $U$ must be a subgroup of $G$. 
Reps must contain canonical orbit representatives, 
otherwise this function does not work.
new\_reps will contain the new canonical representatives of orbits under $G$ 
obtained from cononization of orbit representatives in Reps. 
new\_rep\_idx will contain the index of the new representative 
for each (old) representative in Reps.
All integers of representatives start with 1.
integers in new\_rep\_idx start with 0.  
#endif
{
	INT *theX = NIL;
	MATRIX_OB TG;
	VECTOR_OP pR;
	VECTOR_OB R1;
	LABRA_OB labra_A;
	INTEGER_OB new_reps_len;
	SYM_OB ago;
	INT v, i, ii, j, l, ll, a, idx, f_found;
	PERMUTATION_OB transporter;
	
	l = Reps->s_li();
	v = labra_G->s_degree_i();
	labra_G->calc_transversal_matrix(&TG, FALSE);
	theX = (INT *) my_malloc(v * sizeof(INT), "fuse_orbits");
	new_rep_idx->m_il(l);
	new_reps->m_il(0);
	new_reps_len.m_i(0);
	for (i = 0; i < l; i++) {
		pR = (VECTOR_OP) Reps->s_i(i);
		ll = pR->s_li();
		R1.m_il(ll);
		for (j = 0; j < ll; j++) {
			a = pR->s_ii(j) - 1;
			R1.m_ii(j, a);
			}
		R1.quicksort(R1.s_li(), TRUE);
		for (ii = 0; ii < ll; ii++) {
			a = R1.s_ii(ii);
			theX[ii] = a;
			}
#ifdef SYM_GEO
		geo_canonicize_set(v, ll, theX, 
			labra_G, &TG, FALSE /* f_v */, FALSE /* f_vv */, 
			TRUE /* f_get_aut_group */, &labra_A, &ago, &transporter);
#else
		ago.m_i_i(-1);
#endif
		for (ii = 0; ii < ll; ii++) {
			a = theX[ii] + 1;
			R1.m_ii(ii, a);
			}
		new_reps->search(new_reps_len.s_i(), TRUE /* f_ascending */, &R1, 
			&idx, &f_found);
		if (!f_found) {
			for (ii = 0; ii < i; ii++) {
				a = new_rep_idx->s_ii(ii);
				if (a >= idx)
					new_rep_idx->s_i(ii)->inc();
				}
			new_reps->insert_at(&new_reps_len, idx, &R1);
			new_rep_idx->m_ii(i, idx);
			if (f_v) {
				printf("fuse_orbits: i=%ld, found new orbit representative: ", i);
				R1.println();
				}
			}
		else {
			idx--;
			new_rep_idx->m_ii(i, idx);
			}
		}
	new_reps->realloc_z(new_reps_len.s_i());
	if (f_v) {
		printf("fuse_orbits: new representatives:\n");
		new_reps->Print();
		printf("fuse_orbits: new_rep_idx= ");
		new_rep_idx->println();
		}
	my_free(theX);
	return OK;
}

#if TEXDOCU
INT build_KM_matrix_s_sp1(MATRIX_OP M, VECTOR_OP Orbits, 
	VECTOR_OP Orbits_above1, VECTOR_OP Orbits_above2, INT s)
#endif
{
	VECTOR_OP pOabove1, pOabove2;
	VECTOR_OP ppOabove1, ppOabove2;
	VECTOR_OP Os, Osp1;
	INT m, n, i, j, h, l;
	
	Os = (VECTOR_OP) Orbits->s_i(s);
	Osp1 = (VECTOR_OP) Orbits->s_i(s + 1);
	m = Os->s_li();
	n = Osp1->s_li();
	M->m_ilih_n(n, m);
	pOabove1 = (VECTOR_OP) Orbits_above1->s_i(s + 1);
	pOabove2 = (VECTOR_OP) Orbits_above2->s_i(s + 1);
	for (j = 0; j < n; j++) {
		ppOabove1 = (VECTOR_OP) pOabove1->s_i(j);
		ppOabove2 = (VECTOR_OP) pOabove2->s_i(j);
		l = ppOabove1->s_li();
		for (h = 0; h < l; h++) {
			i = ppOabove1->s_ii(h);
			ppOabove2->s_i(h)->copy(M->s_ij(i, j));
			}
		}
	return OK;
}

#if TEXDOCU
INT orbits_below_to_orbits_above(VECTOR_OP Orbit_lengths, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2, 
	VECTOR_OP Orbits_above1, VECTOR_OP Orbits_above2, 
	SYM_OP go, INT s)
#endif
{
	INT len_s, len_sp1;
	VECTOR_OP Olen_s, Olen_sp1;
	VECTOR_OP Obelow1, Obelow2;
	VECTOR_OP Oabove1, Oabove2;
	VECTOR_OP pObelow1, pObelow2;
	VECTOR_OP pOabove1, pOabove2;
	INT i, j, h, l;
	SYM_OB a, a1, b, b1, c, d, e;
	
	Olen_s = (VECTOR_OP) Orbit_lengths->s_i(s);
	Olen_sp1 = (VECTOR_OP) Orbit_lengths->s_i(s + 1);
	len_s = Olen_s->s_li();
	len_sp1 = Olen_sp1->s_li();
	Obelow1 = (VECTOR_OP) Orbits_below1->s_i(s + 1);
	Obelow2 = (VECTOR_OP) Orbits_below2->s_i(s + 1);
	Oabove1 = (VECTOR_OP) Orbits_above1->s_i(s + 1);
	Oabove2 = (VECTOR_OP) Orbits_above2->s_i(s + 1);
	Oabove1->m_il(len_sp1);
	Oabove2->m_il(len_sp1);
	for (j = 0; j < len_sp1; j++) {
		Olen_sp1->s_i(j)->copy(&b);
		go->ganzdiv_integral(&b, &b1);
		pObelow1 = (VECTOR_OP) Obelow1->s_i(j);
		pObelow2 = (VECTOR_OP) Obelow2->s_i(j);
		pOabove1 = (VECTOR_OP) Oabove1->s_i(j);
		pOabove2 = (VECTOR_OP) Oabove2->s_i(j);
		l = pObelow1->s_li();
		pOabove1->m_il(l);
		pOabove2->m_il(l);
		for (h = 0; h < l; h++) {
			i = pObelow1->s_ii(h);
			pOabove1->m_ii(h, i);
			Olen_s->s_i(i)->copy(&a);
			go->ganzdiv_integral(&a, &a1);
			
			pObelow2->s_i(h)->copy(&c);
			b1.mult(&c, &d);
			d.ganzdiv_integral(&a1, &e);
			e.copy(pOabove2->s_i(h));
			}
		}
	return OK;
}

#if TEXDOCU
INT print_orbit_info(INT t, INT k, 
	VECTOR_OP Orbits, VECTOR_OP Orbit_lengths, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2)
#endif
{
	INT i, j, len;
	VECTOR_OP pO, pOlen, pObelow1, pObelow2;
	
	for (i = t; i <= k; i++) {
		printf("%ld-orbit info:\n", i);
		pO = (VECTOR_OP) Orbits->s_i(i);
		pOlen = (VECTOR_OP) Orbit_lengths->s_i(i);
		pObelow1 = (VECTOR_OP) Orbits_below1->s_i(i);
		pObelow2 = (VECTOR_OP) Orbits_below2->s_i(i);
		len = pO->s_li();
		for (j = 0; j < len; j++) {
			printf("%ld: ", j);
			pO->s_i(j)->print();
			printf(" ago=");
			pOlen->s_i(j)->print();
			if (i > 0) {
				printf(" below:");
				pObelow1->s_i(j)->print();
				pObelow2->s_i(j)->print();
				// printf(" above:");
				// pOabove1->s_i(j)->print();
				// pOabove2->s_i(j)->print();
				}
			printf("\n");
			}
		}
	return OK;
}


#endif /* LADDER_TRUE */

/* end of kramer_mesner.C */
