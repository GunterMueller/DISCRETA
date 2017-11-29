/* ma_geo.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef MATRIXTRUE

#include <DISCRETA/ma.h>

#ifdef PERMTRUE
#include <DISCRETA/perm.h>
#endif

#include <DISCRETA/lb.h>

#ifdef SYM_GEO
#include <DISCRETA/geo.h>
#endif

static INT mtx_compare(MATRIX_OP I, MATRIX_OP J);

#if TEXDOCU
INT MATRIX_OB::automorphism_group_on_cols(PERMUTATION_OP p, 
	PERMUTATION_OP q, LABRA_OP aut_on_canonical_form)
#else
computes the automorphism group {\em on the canonical form} !
p maps the rows of the matrix to the canonical form, 
q maps the columns of the matrix to the canonical form.

you can get the generators on the original matrix with
//PRE
aut.reduced\_generating\_set(aut\_gen, FALSE /* f\_bottom\_up */, TRUE /* f\_v */);
q.invers(qv);
vec\_conjugate(aut\_gen, qv);
aut\_gen.Print(); 
///PRE
#endif
{
	INT *theX;
	INT l, h, i, j, a, nb_X, n, back_to;
	INT t0, t1, user_time;
	BYTE str[1024];
	SYM_OB ago;
	
	t0 = os_ticks();
	l = s_li();
	h = s_hi();
	nb_X = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (s_iji(i, j))
				nb_X++;
			}
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a == 0)
				continue;
			theX[n++] = i * l + j;
			}	
		}
#ifdef SYM_GEO
	geo_canon_simple(FALSE, &back_to, h, l, nb_X, TRUE /* f_print_dots */, 
		theX, p, q, FALSE /* f_transposed */, 
		TRUE /* f_get_aut_group */, aut_on_canonical_form, 
		FALSE /* f_canon_v */,  FALSE /* f_canon_vv */);
#else
	return error("MATRIX_OB::automorphism_group_on_cols(): GEOLIB not available !");
#endif
	// changes theX !

#if 0
	canonical_form->m_ilih_n(l, h);
	for (k = 0; k < nb_X; k++) {
		a = theX[k];
		j = a % l;
		i = a / l;
		canonical_form->m_iji(i, j, 1);
		}
#endif
	printf("\n");
	
	printf("automorphism group order ");
	aut_on_canonical_form->group_order(&ago);
	ago.println();
	
	
	my_free(theX);

	t1 = os_ticks();
	user_time = t1 - t0;
	str[0] = 0;
	print_delta_time(user_time, str);
	printf("total computing time: %s\n", str);
	fflush(stdout);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::canonical_form(INT f_row_perms, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, 
	VECTOR_OP DDp, VECTOR_OP DDb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	LABRA_OP aut, INT f_aut_on_canonical_form, VECTOR_OP aut_gen, 
	INT f_apply_perms_to_canonical_form, INT f_v, INT f_vv)
#else
computes the canonical form 
via row-permutations if f\_row\_perms is true,
via col-permutations if f\_row\_perms is false.
p maps the rows of the matrix to the canonical form, 
q maps the columns of the matrix to the canonical form.
aut\_on\_canonical\_form is the automorphism group on the 
canonical form if f\_aut\_on\_canonical\_form is TRUE, 
on the original matrix if FALSE. 
aut contains row-permutations if f\_row\_perms is TRUE, 
col-permutations otherwise.
The matrix will be changed to contain the canonical form if 
f\_apply\_perms\_to\_canonical\_form is TRUE.
#endif
{
	INT *theX;
	INT l, h, i, j, k, a, nb_X, n, back_to;
	INT f_transpose;
	INT t0, t1, user_time;
	BYTE str[1024];
	SYM_OB ago;
	LABRA_OB aut_on_canonical_form;
	VECTOR_OB aut_gen1;
	PERMUTATION_OB pqv;
	MATRIX_OB canonical_form;
	
	t0 = os_ticks();
	f_transpose = f_row_perms;
	l = s_li();
	h = s_hi();
	p->m_il(h);
	q->m_il(l);
	aut_on_canonical_form.Init_no_generators(l, FALSE, FALSE);
	nb_X = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (s_iji(i, j))
				nb_X++;
			}
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a == 0)
				continue;
			theX[n++] = i * l + j;
			}	
		}
#ifdef SYM_GEO
	// printf("MA::canonical_form() calling geo_canon_simple()\n");
	// fflush(stdout);
	geo_canon_with_initial_decomposition_and_ddp_ddb(
		FALSE /* f_maxtest */, &back_to, 
		h, l, nb_X, 
		f_v /* f_print_dots */, theX, 
		row_decomp, col_decomp, 
		f_ddp, f_ddb, 
		DDp, DDb, 
		p, q, f_transpose, 
		TRUE /* f_get_aut_group */, &aut_on_canonical_form, 
		f_v /* f_canon_v */,  f_vv /* f_canon_vv */);
	// printf("MA::canonical_form() finished with geo_canon_simple()\n");
	// fflush(stdout);
#else
	return error("MATRIX_OB::canonical_form(): GEOLIB not available !");
#endif
	// changes theX !

	if (f_apply_perms_to_canonical_form) {
		canonical_form.m_ilih_n(l, h);
		for (k = 0; k < nb_X; k++) {
			a = theX[k];
			j = a % l;
			i = a / l;
			canonical_form.m_iji(i, j, 1);
			}
		swap(&canonical_form);
		}
	if (f_v) {
		printf("\n");
		}
	
	aut_on_canonical_form.group_order(&ago);
	if (f_v) {
		printf("automorphism group order ");
		ago.println();
		}
	// printf("MA::canonical_form() calling reduced_generating_set()\n");
	// fflush(stdout);
	if (!f_aut_on_canonical_form) {
		// aut_on_canonical_form.generators(&aut_gen1);
		aut_on_canonical_form.reduced_generating_set(&aut_gen1, FALSE /* f_bottom_up */, TRUE /* f_v */);
		if (f_transpose) {
			p->invers(&pqv);
			}
		else {
			q->invers(&pqv);
			}
		vec_conjugate(&aut_gen1, &pqv);
		aut->init_quick(&aut_gen1);
		aut_gen1.swap(aut_gen);
		// aut_gen.Print(); 
		}
	else {
		// aut_on_canonical_form.generators(aut_gen);
		aut_on_canonical_form.reduced_generating_set(aut_gen, FALSE /* f_bottom_up */, TRUE /* f_v */);
		aut_on_canonical_form.swap(aut);
		}
	// printf("MA::canonical_form() finished with reduced_generating_set()\n");
	// fflush(stdout);
	
	my_free(theX);

	t1 = os_ticks();
	user_time = t1 - t0;
	str[0] = 0;
	print_delta_time(user_time, str);
	if (f_v) {
		printf("total computing time: %s\n", str);
		fflush(stdout);
		}
	return OK;
}


#if TEXDOCU
INT MATRIX_OB::build_incidence_matrix(VECTOR_OP B, INT v, INT f_entries_start_with_1)
#else
This function builds the incidence matrix out of the block system in B. 
v is the number of points and will be the number of rows of the incidence matrix. 
if f\_entries\_start\_with\_1 is TRUE, a one is subtracted from each value in the blocks. 
#endif
{
	VECTOR_OP Bj;
	INT b;
	INT i, j, k, l;

	b = B->s_li();
	m_ilih_n(b, v);
	for (j = 0; j < b; j++) {
		Bj = (VECTOR_OP) B->s_i(j);
		if (Bj->s_obj_k() != VECTOR)
			return error("MA::build_incidence_matrix() Bj is not a vector!");
		l = Bj->s_li();
		for (k = 0; k < l; k++) {
			i = Bj->s_ii(k);
			if (f_entries_start_with_1)
				i--;
			if (i >= v)
				return error("MA::build_incidence_matrix() i >= v!");
			m_iji(i, j, 1);
			}
		}
	return OK;
}


#if TEXDOCU
INT MATRIX_OB::apply_perms(PERMUTATION_OP p, PERMUTATION_OP q)
#else
applies p and q to rows and columns of the matrix object. 
This means that $I[i][j]$ is mapped to $I[p(i)][q(i)]$.
#endif
{
	MATRIX_OB I;
	INT h, l, i, j, i1, j1, a;
	
	h = s_hi();
	l = s_li();
	I.m_ilih_n(l, h);
	for (i = 0; i < h; i++) {
		i1 = p->s_ii(i) - 1;
		for (j = 0; j < l; j++) {
			j1 = q->s_ii(j) - 1;
			a = s_iji(i, j);
			I.m_iji(i1, j1, a);
			}
		}
	swap(&I);
	return TRUE;
}

#if TEXDOCU
INT MATRIX_OB::calc_TDO(INT f_calc_second_tdo, INT f_ddp, INT f_ddb, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	VECTOR_OP DDp, VECTOR_OP DDb, 
	MATRIX_OP tdo_scheme, 
	INT f_v, INT f_vv)
#else
computes TDO (first or second, depending on f\_calc\_second\_tdo) of the current matrix. 
$p$ and $q$ contain the TDO permutation which take the original matrix into its TDO form. 
This function interprets any matrix as a 0/1 matrix. 
Each non-zero entry is regarded as a 1.
row\_decomp and col\_decomp are initial decompositions which are respected by the TDO. 
Afterwords, they are refined and contain the TDO decomposition. 
tdo\_scheme contains the matrix of row sums of the TDO decomposition. 
#endif
{
	INT m, n, l, h, i, j, a, nb_X, n1;
	INT *theX;
	VECTOR_OB aut_gen;
	// INT f_TDA_v = f_v;
	// INT f_TDA_vv = f_vv;
	TDOSS *tdoss;
	MATRIX_OB I1;

	l = s_li();
	h = s_hi();
	nb_X = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (s_iji(i, j))
				nb_X++;
			}
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n1 = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a == 0)
				continue;
			theX[n1++] = i * l + j;
			}	
		}

	// f_v = TRUE;
	// f_vv = TRUE;
		
	printf("calling calc_ntdo_()\n"); fflush(stdout);
	tdoss = calc_ntdo_(h, l, 
		nb_X, theX, FALSE /* f_multivalued */, NIL /* *theVal */, 
		f_calc_second_tdo, 
		f_ddp, DDp, 
		f_ddb, DDb, 
		f_v /* f_v */, 
		f_vv /* f_vv */, 
		TRUE /* f_dd_v */, 
		FALSE /* f_dd_vv */, 
		p, q, 
		row_decomp, col_decomp);
	printf("finished with calc_ntdo_()\n"); fflush(stdout);
	m = tdoss_m(tdoss);
	n = tdoss_n(tdoss);
	tdo_scheme->m_ilih_n(n, m);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			a = tdoss_ij(tdoss, i, j);
			tdo_scheme->m_iji(i, j, a);
			}
		}
	printf("TDO-scheme:\n");
	tdo_scheme->Print();
	fflush(stdout);
	
#if 0
	tdoss_free(tdoss);
	printf("tdoss freed\n");fflush(stdout);
#endif

	copy(&I1);
	// printf("nach copy\n");fflush(stdout);
	I1.apply_perms(p, q);
	// printf("nach I1.apply_perms\n");fflush(stdout);
	
	// printf("I1=\n");
	// I1.Print();
	// fflush(stdout);
	
#if 0
	I1.calc_decomposition_matrix(row_decomp, col_decomp, tdo_scheme);
	printf("TDO-scheme:\n");
	tdo_scheme->Print();
	fflush(stdout);
	
	
	my_free(theX);
#endif
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::calc_TDA(FILE *TEX_fp, INT nb_geo, INT f_row_perms, 
	INT f_has_aut, LABRA_OP aut, 
	INT *nb_row_blocks_TDO, INT *nb_col_blocks_TDO, 
	INT *nb_row_blocks_TDA, INT *nb_col_blocks_TDA, 
	INT f_v, INT f_vv, INT f_incidences)
#else
This function calculates the TDA of the matrix. Either f\_has\_aut is TRUE and aut 
contains the automorphism group and the matrix is in its canonical form 
or the automorphism group is recalculated via MA::canonical\_form().
#endif
{
	MATRIX_OB I;
	PERMUTATION_OB p, q;
	LABRA_OB aut1;
	INT l, h, i, j, a, nb_X, n;
	INT *theX;
	VECTOR_OB aut_gen;
	INT f_TDA_v = f_v;
	INT f_TDA_vv = f_vv;
	VECTOR_OB row_decomp, col_decomp;

	copy(&I);
	l = I.s_li();
	h = I.s_hi();
	row_decomp.m_il(1);
	row_decomp.m_ii(0, l);
	col_decomp.m_il(1);
	col_decomp.m_ii(0, h);
	if (f_has_aut) {
		aut->copy(&aut1);
		}
	else {
		I.canonical_form(f_row_perms, &row_decomp, &col_decomp, 
			FALSE /* f_ddp */, FALSE /* f_ddb */,
			NULL /* DDp */, NULL /* DDb */, 
			&p, &q, 
			&aut1, TRUE /* f_aut_on_canonical_form */, &aut_gen, 
			TRUE /* f_apply_perms_to_canonical_form */, 
			TRUE /* f_v */, FALSE /* f_vv */);
		}
	nb_X = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			if (I.s_iji(i, j))
				nb_X++;
			}
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = I.s_iji(i, j);
			if (a == 0)
				continue;
			theX[n++] = i * l + j;
			}	
		}

	::calc_TDA(NIL, TEX_fp, nb_geo, 
		h, l, nb_X, theX, 
		f_row_perms, &aut1,
		nb_row_blocks_TDO, nb_col_blocks_TDO, 
		nb_row_blocks_TDA, nb_col_blocks_TDA, 
		f_TDA_v, f_TDA_vv, f_incidences);
	
	my_free(theX);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::calc_blocks(VECTOR_OP B, PERMUTATION_OP q, PERMUTATION_OP qv, INT f_v)
#endif
{
	VECTOR_OP block;
	VECTOR_OB B1;
	INT i, j, v, b;
	INT nb_empty, idx, f_found;
	
	b = s_li();
	v = s_hi();
	B->m_il(b);
	for (j = 0; j < b; j++) {
		block = (VECTOR_OP) B->s_i(j);
		block->m_il(0);
		}
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_iji(i, j) == 0)
				continue;
			block = (VECTOR_OP) B->s_i(j);
			block->inc();
			block->m_ii(block->s_li() - 1, i);
			}
		}
	for (j = 0; j < b; j++) {
		block = (VECTOR_OP) B->s_i(j);
		block->quicksort(block->s_li(), TRUE /* f_ascending */);
		}
	B->copy(&B1);
	B->quicksort(B->s_li(), TRUE /* f_ascending */);
	if (f_v) {
		printf("the original blocks:\n");
		for (i = 0; i < b; i++) {
			block = (VECTOR_OP) B1.s_i(i);
			printf("block %ld: ", i);
			block->println();
			}
		printf("the lexicographic blocks:\n");
		for (i = 0; i < b; i++) {
			block = (VECTOR_OP) B->s_i(i);
			printf("block %ld: ", i);
			block->println();
			}
		}
	q->m_il(b);
	nb_empty = 0;
	for (i = 0; i < b; i++) {
		block = (VECTOR_OP) B1.s_i(i);
		if (block->s_li() == 0) {
			idx = nb_empty;
			nb_empty++;
			}
		else {
			B->search(b, TRUE, block, &idx, &f_found);
			if (!f_found)
				return error("calc_blocks(): not found");
			idx--;
			}
		q->m_ii(i, idx + 1);
		if (f_v) {
			printf("original block %ld is lexicographic block %ld\n", i, idx);
			}
		}
	q->invers(qv);
	return OK;
}

INT MATRIX_OB::action_on_blocks(VECTOR_OP gen_row_perms, 
	VECTOR_OP gen_col_perms, INT f_v, INT f_vv)
{
	INT i, l;
	VECTOR_OB B;
	PERMUTATION_OB qlex, qlexv, tmp1, tmp2;
	VECTOR_OB gen1, gen2;
	
	calc_blocks(&B, &qlex, &qlexv, f_vv);
	
	gen_row_perms->copy(&gen1);
	l = gen1.s_li();
	vec_induce_action_on_blocks(&gen1, &B);
	gen2.m_il(l);
	for (i = 0; i < l; i++) {
		qlex.mult((PERMUTATION_OP) gen1.s_i(i), &tmp1);
		tmp1.mult(&qlexv, &tmp2);
		tmp2.copy((PERMUTATION_OP) gen2.s_i(i));
		}
	if (f_v) {
		l = gen1.s_li();
		for (i = 0; i < l; i++) {
			printf("generator %ld on points: ", i);
			gen1.s_i(i)->println();
			printf("on blocks: ");
			gen2.s_i(i)->println();
			}
		}
	gen2.swap(gen_col_perms);
	return OK;
}

void MATRIX_OB::print_incidences()
{
	INT i, j, v, b;
	
	v = s_hi();
	b = s_li();
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_iji(i, j))
				printf("X");
			else
				printf(".");
			}
		printf("\n");
		}
	printf("\n");
}

#if TEXDOCU
INT MATRIX_OB::configuration_graph(MATRIX_OP G, MATRIX_OP I_and_G, MATRIX_OP G_and_I)
#endif
{
	VECTOR_OB B, Bj;
	INT i, j, k, a, v, b, l;

	B.m_il(0);
	l = 0;
	v = s_hi();
	b = s_li();
	for (i = 0; i < v; i++) {
		for (j = i + 1; j < v; j++) {
			k = find_ij_line(i, j);
			if (k >= 0)
				continue;
			Bj.m_il(2);
			Bj.m_ii(0, i);
			Bj.m_ii(1, j);
			B.inc();
			Bj.swap((VECTOR_OP) B.s_i(l));
			l++;
			}
		}
	G->build_incidence_matrix(&B, v, FALSE /* f_entries_start_with_1 */);
	I_and_G->m_ilih(b + l, v);
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = s_iji(i, j);
			I_and_G->m_iji(i, j, a);
			}
		}
	for (i = 0; i < v; i++) {
		for (j = 0; j < l; j++) {
			a = G->s_iji(i, j);
			I_and_G->m_iji(i, b + j, a);
			}
		}
	G_and_I->m_ilih(b + l, v);
	for (i = 0; i < v; i++) {
		for (j = 0; j < l; j++) {
			a = G->s_iji(i, j);
			G_and_I->m_iji(i, j, a);
			}
		}
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = s_iji(i, j);
			G_and_I->m_iji(i, l + j, a);
			}
		}
	
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::find_ij_line(INT i, INT j)
#endif
{
	INT k, b;

	b = s_li();
	for (k = 0; k < b; k++) {
		if (s_iji(i, k) && s_iji(j, k))
			return k;
		}
	return -1;
}

#if TEXDOCU
INT MATRIX_OB::output_GEO(FILE *fp)
#endif
{
	INT i, j, a, l, h, nb_X = 0;
	
	l = s_li();
	h = s_hi();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a) {
				fprintf(fp, "%ld ", i * l + j);
				nb_X++;
				}
			}
		}
	fprintf(fp, "\n");
	printf("MA::output_GEO(): nb_X = %ld\n", nb_X);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::output_inc_tex(FILE *fp)
#endif
{
	INT i, j, a, l, h;
	BYTE *alphabet = NIL;
	
	l = s_li();
	h = s_hi();
	alphabet = (BYTE *) my_malloc((l + 11) * sizeof(BYTE), "alphabet");
	for (i = 0; i < 10; i++)
		alphabet[i] = '0' + i;
	for (i = 10; i < l; i++) {
		j = i - 10;
		if (j < 26)
			alphabet[i] = 'a' + j;
		else {
			j -= 26;
			alphabet[i] = 'A' + j;
			}
		}
	alphabet[l] = 0;
	// printf("alphabet: %s\n", alphabet);
	
	fprintf(fp, "$\\begin{array}{*{%ld}{c}}\n", h);
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a) {
				fprintf(fp, "%c", alphabet[j]);
				}
			}
		if (i < h - 1) {
			fprintf(fp, "&");
			}
		else {
			fprintf(fp, "\\\\\n");
			}
		}
	fprintf(fp, "\\end{array}$\n");
	fprintf(fp, "\n");
	my_free(alphabet);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::print_char_TD(FILE *fp, INT f_tex, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp)
#endif
{
	INT v, b, nb_X, m, n, i, j, a, *theX, f;
	PERMUTATION_OB p, q;
	VECTOR_OB row_first, col_first;
	
	m = row_decomp->s_li();
	row_first.m_il(m);
	f = 0;
	for (i = 0; i < m; i++) {
		row_first.m_ii(i, f);
		f += row_decomp->s_ii(i);
		}
	n = col_decomp->s_li();
	col_first.m_il(n);
	f = 0;
	for (i = 0; i < n; i++) {
		col_first.m_ii(i, f);
		f += col_decomp->s_ii(i);
		}
	
	v = s_hi();
	b = s_li();
	p.m_il(v);
	q.m_il(b);
	p.one();
	q.one();
	nb_X = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_iji(i, j))
				nb_X++;
			}
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = s_iji(i, j);
			if (a == 0)
				continue;
			theX[n++] = i * b + j;
			}	
		}

	
	TD_char_inv_print(fp, f_tex, 
		v, b, nb_X, theX, 
		&p, &q, 
		&row_first, &col_first, 
		&row_first, &col_first);
	
	my_free(theX);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::print_char_and_inv_TD(FILE *fp, INT f_tex, 
	VECTOR_OP row_decomp_char, VECTOR_OP col_decomp_char, 
	VECTOR_OP row_decomp_inv, VECTOR_OP col_decomp_inv)
#endif
{
	INT v, b, nb_X, m, n, i, j, a, *theX, f;
	PERMUTATION_OB p, q;
	VECTOR_OB row_first_char, col_first_char;
	VECTOR_OB row_first_inv, col_first_inv;
	
	m = row_decomp_char->s_li();
	row_first_char.m_il(m);
	f = 0;
	for (i = 0; i < m; i++) {
		row_first_char.m_ii(i, f);
		f += row_decomp_char->s_ii(i);
		}
	n = col_decomp_char->s_li();
	col_first_char.m_il(n);
	f = 0;
	for (i = 0; i < n; i++) {
		col_first_char.m_ii(i, f);
		f += col_decomp_char->s_ii(i);
		}
	
	m = row_decomp_inv->s_li();
	row_first_inv.m_il(m);
	f = 0;
	for (i = 0; i < m; i++) {
		row_first_inv.m_ii(i, f);
		f += row_decomp_inv->s_ii(i);
		}
	n = col_decomp_inv->s_li();
	col_first_inv.m_il(n);
	f = 0;
	for (i = 0; i < n; i++) {
		col_first_inv.m_ii(i, f);
		f += col_decomp_inv->s_ii(i);
		}
	
	v = s_hi();
	b = s_li();
	p.m_il(v);
	q.m_il(b);
	p.one();
	q.one();
	nb_X = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (s_iji(i, j))
				nb_X++;
			}
		}
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	
	n = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = s_iji(i, j);
			if (a == 0)
				continue;
			theX[n++] = i * b + j;
			}	
		}

	
	TD_char_inv_print(fp, f_tex, 
		v, b, nb_X, theX, 
		&p, &q, 
		&row_first_char, &col_first_char, 
		&row_first_inv, &col_first_inv);
	
	my_free(theX);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::incma_latex(FILE *fp, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_labelled, 
	VECTOR_OP labelling_points, VECTOR_OP labelling_blocks, INT offset)
#endif
{
	INT v, b, m, n, i, j, a, nb_X, r, max_r;
	INT *Vi, *Bj, *X, *R;

	v = s_hi();
	b = s_li();
	m = row_decomp->s_li();
	n = col_decomp->s_li();
	nb_X = 0;
	max_r = 0;
	for (i = 0; i < v; i++) {
		r = 0;
		for (j = 0; j < b; j++) {
			if (s_iji(i, j))
				r++;
			}
		nb_X += r;
		max_r = MAXIMUM(max_r, r);
		}
	Vi = (INT *) my_malloc(m * sizeof(INT), "MA::incma_latex Vi");
	Bj = (INT *) my_malloc(n * sizeof(INT), "MA::incma_latex Bj");
	X = (INT *) my_malloc(v * max_r * sizeof(INT), "MA::incma_latex X");
	R = (INT *) my_malloc(v * sizeof(INT), "MA::incma_latex R");
	for (i = 0; i < m; i++) {
		a = row_decomp->s_ii(i);
		Vi[i] = a;
		}
	for (i = 0; i < n; i++) {
		a = col_decomp->s_ii(i);
		Bj[i] = a;
		}
	for (i = 0; i < v; i++) {
		r = 0;
		for (j = 0; j < b; j++) {
			if (s_iji(i, j)) {
				X[i * max_r + r] = j;
				r++;
				}
			}
		R[i] = r;
		}
	if (f_labelled) {
		::incma_latex_integer_labels(fp, v, b, m, n, 
			Vi, Bj, R, X, max_r /* dim_X */, 
			labelling_points, labelling_blocks, offset);
		}
	else {
		::incma_latex(fp, v, b, m, n, Vi, Bj, R, X, max_r /* dim_X */);
		}
	my_free(Vi);
	my_free(Bj);
	my_free(X);
	my_free(R);
	return OK;
}

INT print_tdo_scheme(MATRIX_OP tdo_scheme, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, INT f_tex, FILE *fp)
{
	INT m, n, i, j, a;
	
	m = row_decomp->s_li();
	n = col_decomp->s_li();
	if (!f_tex) {
		fprintf(fp, "      ");
		for (j = 0; j < n; j++) {
			a = col_decomp->s_ii(j);
			fprintf(fp, "%3ld ", a);
			}
		fprintf(fp, "\n");
		fprintf(fp, "------");
		for (j = 0; j < n; j++) {
			fprintf(fp, "----");
			}
		fprintf(fp, "\n");
		for (i = 0; i < m; i++) {
			a = row_decomp->s_ii(i);
			fprintf(fp, "%3ld | ", a);
			for (j = 0; j < n; j++) {
				a = tdo_scheme->s_iji(i, j);
				fprintf(fp, "%3ld ", a);
				}
			fprintf(fp, "\n");
			}
		}
	else {
		fprintf(fp, "\\[\n");
		fprintf(fp, "\\begin{array}{r|*{%ld}{r}}\n", n);
		fprintf(fp, "& ");
		for (j = 0; j < n; j++) {
			a = col_decomp->s_ii(j);
			fprintf(fp, "%3ld ", a);
			if (j < n - 1)
				fprintf(fp, "& ");
			}
		fprintf(fp, "\\\\\n");
		fprintf(fp, "\\hline\n");
		for (i = 0; i < m; i++) {
			a = row_decomp->s_ii(i);
			fprintf(fp, "%3ld & ", a);
			for (j = 0; j < n; j++) {
				a = tdo_scheme->s_iji(i, j);
				fprintf(fp, "%3ld ", a);
				if (j < n - 1)
					fprintf(fp, "& ");
				}
			fprintf(fp, "\\\\\n");
			}
		fprintf(fp, "\\end{array}\n");
		fprintf(fp, "\\]\n");
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::calc_blocks(VECTOR_OP B, 
	PERMUTATION_OP p, 
	PERMUTATION_OP q, PERMUTATION_OP qv, INT f_v)
#else
Calculates all blocks of the incidence system.
The incidence system is permuted by p on the rows before calculating the blocks.
Then, B is sorted. The original block no i 
becomes q[i] in the sorted vector B.
qv is q invers.
#endif
{
	VECTOR_OP block;
	VECTOR_OB B1;
	INT v, b, i, i1, j, idx, f_found, nb_empty;
	INT ii, l;

	v = s_hi();
	b = s_li();
	B->m_il(b);
	for (j = 0; j < b; j++) {
		block = (VECTOR_OP) B->s_i(j);
		block->m_il(0);
		for (i = 0; i < v; i++) {
			if (s_iji(i, j) == 0)
				continue;
			i1 = p->s_ii(i) - 1;
			block->inc();
			block->m_ii(block->s_li() - 1, i1);
			}
		}
	for (j = 0; j < b; j++) {
		block = (VECTOR_OP) B->s_i(j);
		block->quicksort(block->s_li(), TRUE /* f_ascending */);
		}
	B->copy(&B1);
	B->quicksort(B->s_li(), TRUE /* f_ascending */);
	if (f_v) {
		printf("the original blocks:\n");
		for (i = 0; i < b; i++) {
			block = (VECTOR_OP) B1.s_i(i);
			printf("block %ld: ", i);
			printf("[");
			l = block->s_li();
			for (ii = 0; ii < l; ii++) {
				printf("%ld ", block->s_ii(ii) + 1);
				}
			printf("]\n");
			// block->println();
			}
		printf("the lexicographic blocks:\n");
		for (i = 0; i < b; i++) {
			block = (VECTOR_OP) B->s_i(i);
			printf("block %ld: ", i);
			printf("[");
			l = block->s_li();
			for (ii = 0; ii < l; ii++) {
				printf("%ld ", block->s_ii(ii) + 1);
				}
			printf("]\n");
			// block->println();
			}
		}
	q->m_il(b);
	nb_empty = 0;
	for (i = 0; i < b; i++) {
		block = (VECTOR_OP) B1.s_i(i);
		if (block->s_li() == 0) {
			idx = nb_empty;
			nb_empty++;
			}
		else {
			B->search(b, TRUE, block, &idx, &f_found);
			if (!f_found)
				return error("MA::calc_blocks(): not found");
			idx--;
			}
		q->m_ii(i, idx + 1);
		// if (f_v) {
		//	printf("original block %ld is lexicographic block %ld\n", i, idx);
		//	}
		}
	q->invers(qv);
	if (f_v) {
		printf("found %ld blocks:\n", b);
		for (j = 0; j < b; j++) {
			block = (VECTOR_OP) B->s_i(j);
			printf("%ld (original block no %ld)\n", j, qv->s_ii(j) - 1);
			// block->println();
			}
		}
	
	
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::TDA(VECTOR_OP gen, 
	VECTOR_OP row_decomp_char, VECTOR_OP col_decomp_char, 
	PERMUTATION_OP p, PERMUTATION_OP q, 
	VECTOR_OP row_decomp_inv, VECTOR_OP col_decomp_inv, 
	INT f_v, INT f_vv)
#else
This function computes the TDA of a given incidence matrix.
gen must be a vector of generators for the automorphism group of that matrix. 
row\_decomp\_char and col\_decomp\_char are NOT USED by the function!
The matrix is in the original form and $p$ and $q$ map it onto the canonical form. 
gen are the generators for Aut on the canonical form. 
The function conjugates each generator g of gen to $p^{-1}gp$ 
so that it gets the automorphism group on the original matrix.
Afterwards, $p$ and $q$ take the matrix into the TDA form:
the original point i lies at p[i] in the TDA.
the original block j lies at q[j] in the TDA.
row\_decomp\_inv and col\_decomp\_inv describe the decomposition of the TDA.
This function does not use LABRA::orbits but uses orbits\_cheap instead.
#endif
{
	VECTOR_OB P_tda_ordering;
	VECTOR_OB Q_tda_ordering;
	VECTOR_OB P_orbit_first, P_orbit_length;
	VECTOR_OB Q_orbit_first, Q_orbit_length;
	// VECTOR_OB P_O_labra;
	// VECTOR_OB Q_O_labra;
	PERMUTATION_OB pv, qv;
	VECTOR_OB B;
	PERMUTATION_OB qlex, qlexv;
	INT v, b, i, j, l, first, length;
	PERMUTATION_OB tmp1, tmp2;
	VECTOR_OB gen0, gen1, gen2;
	// LABRA_OB A1;
	
	v = s_hi();
	b = s_li();
	// printf("MA::TDA() f_v = %ld f_vv = %ld\n", f_v, f_vv);
	
#if 0
	A->orbits(&P_tda_ordering, 
		&P_orbit_first, &P_orbit_length, 
		FALSE /* f_calc_labras */, &P_O_labra, FALSE /* f_v */);
#else
	orbits_cheap(gen, &P_tda_ordering, &P_orbit_first, &P_orbit_length, FALSE /* f_v */);
#endif
	if (P_tda_ordering.s_li() != v)
		return error("MA::TDA() P_tda_ordering.s_li() != v");
	// A->strong_generating_set(&gen, FALSE /* f_v */);
	
	p->m_il(v);
	for (i = 0; i < v; i++) {
		j = P_tda_ordering.s_ii(i);
		p->m_ii(j, i + 1);
		}
	p->invers(&pv);
	if (f_vv) {
		printf("tda-ordering of the points:\n");
		P_tda_ordering.println();
		} 
	gen->copy(&gen0);
	vec_conjugate(&gen0, p);
		// computes p^{-1} g p
	if (f_vv) {
		l = gen->s_li();
		for (i = 0; i < l; i++) {
			printf("generator %ld on points: ", i);
			gen->s_i(i)->println();
			printf("on TDA ordering of points: ");
			gen0.s_i(i)->println();
			}
		}
	
	l = P_orbit_first.s_li();
	for (i = 0; i < l; i++) {
		first = P_orbit_first.s_ii(i);
		length = P_orbit_length.s_ii(i);
		if (f_vv) {
			printf("%ld/%ld\n", first, length);
			}
		}
	if (f_vv) {
		printf("the incidence matrix:\n");
		Print();
		printf("the TDA ordering of points is obtained by the following permutation\n");
		p->println();
		}
	calc_blocks(&B, p, &qlex, &qlexv, f_vv);
	
	Q_orbit_first.m_il(1);
	Q_orbit_length.m_il(1);
	Q_orbit_first.m_ii(0, 0);
	// Q_orbit_len.m_ii(0, b);

#if 0
	l = gen.s_li();
	gen0.m_il(l);
	for (i = 0; i < l; i++) {
		pv.mult((PERMUTATION_OP) gen.s_i(i), &tmp1);
		tmp1.mult(p, &tmp2);
		tmp2.copy((PERMUTATION_OP) gen0.s_i(i));
		}
#endif
	gen0.copy(&gen1);
	vec_induce_action_on_blocks(&gen1, &B);
	l = gen1.s_li();
	gen2.m_il(l);
	for (i = 0; i < l; i++) {
		qlex.mult((PERMUTATION_OP) gen1.s_i(i), &tmp1);
		tmp1.mult(&qlexv, &tmp2);
		tmp2.copy((PERMUTATION_OP) gen2.s_i(i));
		}
	if (f_vv) {
		l = gen->s_li();
		for (i = 0; i < l; i++) {
			printf("generator %ld on points: ", i);
			gen->s_i(i)->println();
			printf("on TDA ordering of points: ");
			gen0.s_i(i)->println();
			printf("on lexicographic blocks: ");
			gen1.s_i(i)->println();
			printf("on blocks: ");
			gen2.s_i(i)->println();
			}
		}

#if 0
	A1.init_by_generators(&gen2, b, FALSE);
	if (f_vv) {
		printf("the group induced on the blocks ");
		A1.print_group_order();
		}
	A1.orbits(&Q_tda_ordering, &Q_orbit_first, &Q_orbit_length, 
		FALSE /* f_calc_labras */, &Q_O_labra, FALSE /* f_v */);
#else
	orbits_cheap(&gen2, &Q_tda_ordering, &Q_orbit_first, &Q_orbit_length, FALSE /* f_v */);
#endif
	
	if (b != Q_tda_ordering.s_li())
		return error("MA::TDA() b != Q_tda_ordering.s_li()");
	q->m_il(b);
	for (i = 0; i < b; i++) {
		j = Q_tda_ordering.s_ii(i);
		q->m_ii(j, i + 1);
		}
	q->invers(&qv);
	if (f_vv) {
		printf("tda-ordering of the blocks:\n");
		Q_tda_ordering.println();
		}
	l = Q_orbit_first.s_li();
	for (i = 0; i < l; i++) {
		first = Q_orbit_first.s_ii(i);
		length = Q_orbit_length.s_ii(i);
		if (f_vv) {
			printf("%ld/%ld\n", first, length);
			}
		}
	P_orbit_length.copy(row_decomp_inv);
	Q_orbit_length.copy(col_decomp_inv);


	return OK;
}

#if TEXDOCU
INT MATRIX_OB::calc_decomposition_matrix(
	VECTOR_OP row_decomp, VECTOR_OP col_decomp,
	MATRIX_OP S)
#else
This function computes the row decomposition matrix of the decomposition 
given by row\_decomp and col\_decomp. 
By the way, this function also tests if the matrix is really row-tactical.
#endif
{
	INT m, n, i, j, s;
	INT i0, i1, I;
	INT j0, j1, J;
	MATRIX_OB It, S1;

	m = row_decomp->s_li();
	n = col_decomp->s_li();
	S->m_ilih_n(n, m);
	i0 = 0;
	for (i = 0; i < m; i++) {
		I = row_decomp->s_ii(i);
		for (i1 = 0; i1 < I; i1++) {
			j0 = 0;
			for (j = 0; j < n; j++) {
				J = col_decomp->s_ii(j);
				s = 0;
				for (j1 = 0; j1 < J; j1++) {
					if (s_iji(i0 + i1, j0 + j1))
						s++;
					}
				if (i1 == 0) {
					S->m_iji(i, j, s);
					// printf("s[%ld][%ld]=%ld\n", i, j, s);
					// fflush(stdout);
					}
				else {
					if (s != S->s_iji(i, j)) {
						printf("the matrix:\n");
						Print();
						printf("row_decomp:\n");
						row_decomp->println();
						printf("col_decomp:\n");
						col_decomp->println();
						fflush(stdout);
						return error("MA::calc_decomposition() not row tactical !");
						}
					}
				j0 += J;
				} // next j
			}
		i0 += I;
		} // next i
	// warning: infinite recursion !
	// transpose(&It);
	// It.calc_decomposition_matrix(col_decomp, row_decomp, &S1);
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::test_2_design(INT *k, INT *lambda)
#else
This function tests if the matrix is the incidence matrix of a 2-design. 
This is done by evaluating the product 
$II^t$ and checking if this matrix has $r$ on the diagonal and $\lambda$ elsewhere.
$r$ is returned in $k$ and $\lambda$ in lambda 
($r = k$ for a symmetric design, therefore the mixture in notation).
#endif
{
	MATRIX_OB It, IIt;
	INT v, b, i, j, k1, l1, a;

	transpose(&It);
	mult(&It, &IIt);
	v = s_hi();
	b = s_li();
	if (IIt.s_hi() != v)
		return error("MA::test_2_design() IIt.s_hi() != v");
	if (IIt.s_li() != v)
		return error("MA::test_2_design() IIt.s_li() != v");
	for (i = 0; i < v; i++) {
		for (j = 0; j < v; j++) {
			a = IIt.s_iji(i, j);
			if (i == j) {
				if (i == 0)
					k1 = a;
				else if (k1 != a) {
					printf("at %ld %ld:\n", i, j);
					printf("IIt=\n");
					IIt.Print();
					printf("I=\n");
					Print();
					return FALSE;
					}
				}
			else {
				if (i == 0 && j == 1)
					l1 = a;
				else if (l1 != a) {
					printf("at %ld %ld:\n", i, j);
					printf("IIt=\n");
					IIt.Print();
					printf("I=\n");
					Print();
					return FALSE;
					}
				}
			}
		}
	*k = k1;
	*lambda = l1;
	return TRUE;
}

#if TEXDOCU
INT MATRIX_OB::nb_incidences()
#else
This function counts the number of non-zero entries of a given matrix.
#endif
{
	INT i, j, h, l, a, n = 0;

	h = s_hi();
	l = s_li();
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a)
				n++;
			}
		}
	return n;
}

#if TEXDOCU
INT MATRIX_OB::levi(INT f_blocked, MATRIX_OP Levi)
#else
This function computes the incidence system for the levi graph for the given incidence system. 
The levi graph is a bipartite graph and can be described as follows:
The levi graph has the points and blocks of I as its points and has an edge between 
point $p$ and (original) block $B$ if and only if $p$ lied on $B$ in $I$. 
if f\_blocked is TRUE, two additional columns are added which contain the 
incidence vectors of all the points and all the blocks. 
#endif
{
	INT i, j, h, l, a, nb_X, n;

	h = s_hi();
	l = s_li();
	nb_X = nb_incidences();
	if (f_blocked) {
		Levi->m_ilih_n(nb_X + 2, h + l);
		for (i = 0; i < h; i++) {
			Levi->m_iji(i, 0, 1);
			}
		for (i = 0; i < l; i++) {
			Levi->m_iji(h + i, 1, 1);
			}
		n = 2;
		}
	else {
		Levi->m_ilih_n(nb_X, h + l);
		n = 0;
		}
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = s_iji(i, j);
			if (a) {
				Levi->m_iji(i, n, 1);
				Levi->m_iji(h + j, n, 1);
				n++;
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT MATRIX_OB::self_dual(INT f_self_polar, 
	INT *f_is_sd, INT *f_is_sp, VECTOR_OP aut_gen, PERMUTATION_OP duality, 
	INT f_v, INT f_vv)
#else
This function checks if the incidence system (with $v=b$) is self dual 
and (if f\_self\_polar is TRUE) checks if it is self polar in the case that it is self dual.
f\_is\_sd is TRUE iff $I$ is self-dual, 
f\_is\_sp is TRUE iff $I$ is self-polar. 
#endif
{
	INT v, b, ret, i, j, o;
	VECTOR_OB row_decomp, col_decomp;
	PERMUTATION_OB p1, q1, p2, q2, p2v, q2v, aa, bb, c, d, p;
	VECTOR_OB aut_P, aut_B, aut_P_t, path;
	LABRA_OB aut, aut_t, Aut;
	MATRIX_OB Ican, J, Jcan;
	// MATRIX_OB transversals;
	
	v = s_hi();
	b = s_li();
	if (v != b)
		return error("MA::self_dual(): v != b");
	
	aut_gen->m_il(0);
	duality->m_il(0);
	
	// the trivial decomposition:
	row_decomp.m_il(1);
	row_decomp.m_ii(0, v);
	col_decomp.m_il(1);
	col_decomp.m_ii(0, b);
		
	canonical_form(TRUE /* f_row_perms */, 
		&row_decomp, &col_decomp, 
		FALSE /* f_ddp */, FALSE /* f_ddb */,
		NULL /* DDp */, NULL /* DDb */, 
		&p1, &q1, 
		&aut, FALSE /* f_aut_on_canonical_form */, &aut_P, 
		FALSE /* f_apply_perms_to_canonical_form */, 
		FALSE, FALSE);
	copy(&Ican);
	Ican.apply_perms(&p1, &q1);
	if (f_v) {
		printf("self_dual: canonical form:\n");
		Ican.print_incidences();
		printf("p1:\n");
		p1.println();
		printf("q1:\n");
		q1.println();
		}
	transpose(&J);
	J.canonical_form(TRUE /* f_row_perms */, 
		&col_decomp, &row_decomp, 
		FALSE /* f_ddp */, FALSE /* f_ddb */,
		NULL /* DDp */, NULL /* DDb */, 
		&p2, &q2, 
		&aut_t, FALSE /* f_aut_on_canonical_form */, &aut_P_t, 
		FALSE /* f_apply_perms_to_canonical_form */, 
		FALSE, FALSE);
	J.copy(&Jcan);
	Jcan.apply_perms(&p2, &q2);
	if (f_v) {
		printf("self_dual: canonical form of transposed:\n");
		Jcan.print_incidences();
		printf("p2:\n");
		p2.println();
		printf("q2:\n");
		q2.println();
		}
	*f_is_sd = FALSE;
	*f_is_sp = FALSE;
	if (mtx_compare(&Ican, &Jcan) == 0) {
		// nb_self_dual++;
		// printf("self_dual ! (no %ld of %ld)\n", nb_self_dual, nb_geo);
		// out_inc(OUT_fp, X, nb_X);
		*f_is_sd = TRUE;
		}
	
	p2.invers(&p2v);
	q2.invers(&q2v);
	p1.mult(&p2v, &aa);
	q1.mult(&q2v, &bb);
	
	action_on_blocks(&aut_P, &aut_B, f_vv, FALSE);
	vec_generators_diagonal_sum(&aut_P, &aut_B, aut_gen);
	
	if (!*f_is_sd) {
		*f_is_sp = FALSE;
		return OK;
		}
	
	c.m_il(v + v);
	for (i = 0; i < v; i++) {
		j = aa.s_ii(i);
		c.m_ii(i, v + j);
		}
	for (i = 0; i < v; i++) {
		j = bb.s_ii(i);
		c.m_ii(v + i, j);
		}
	c.copy(duality);
	if (f_v) {
		printf("c=\n");
		c.println();
		}
	
	if (!f_self_polar)
		return OK;
	
	Aut.init_by_generators(aut_gen, v + b, FALSE /* f_v */);
	Aut.print_group_order();
	Aut.elements_first(&path, &p);
#if 0
	if (f_v) {
		printf("MA::self_dual(): transversals:\n");
		transversals.Print();
		}
#endif
	i = 1;
	while (TRUE) {
		if (f_v) {
			printf("MA::self_dual(): Aut.element %ld: path=", i);
			path.println();
			printf("p=");
			p.println();
			}
		p.mult(&c, &d);
		if (f_v) {
			printf("d=");
			d.println();
			}
		d.order(&o);
		if (f_v) {
			printf("has order %ld\n", o);
			}
		if (o == 2) {
			*f_is_sp = TRUE;
			}
		
		
		ret = Aut.elements_next(&path, &p);
		if (!ret)
			break;
		i++;
		}

	return OK;
}

static INT mtx_compare(MATRIX_OP I, MATRIX_OP J)
{
	INT i, j, a, b, l, h;

	l = I->s_li();
	if (l != J->s_li())
		return error("mtx_compare(): l differs !");
	h = I->s_hi();
	if (h != J->s_hi())
		return error("mtx_compare(): h differs !");
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			a = I->s_iji(i, j);
			b = J->s_iji(i, j);
			if (a < b)
				return -1;
			if (a > b)
				return 1;
			}
		}
	return 0;
}

#if TEXDOCU
INT MATRIX_OB::get_base_block_if_cyclic(PERMUTATION_OP per, VECTOR_OP bb)
#else
#endif
{
	INT i, ii, idx, h, l;
	
	h = s_hi();
	l = 0;
	for (i = 0; i < h; i++) {
		if (s_iji(i, 0))
			l++;
		}
	bb->m_il(l);
	if (h != per->s_li())
		return error("MA::get_base_block_if_cyclic() h != per->s_li()");
	ii = 0;
	idx = 0;
	for (i = 0; i < h; i++) {
		if (s_iji(ii, 0)) {
			bb->m_ii(idx, i);
			idx++;
			}
		ii = per->s_ii(ii) - 1;
		}
	if (idx != l)
		return error("MA::get_base_block_if_cyclic() idx != l");
	if (ii != 0)
		return error("MA::get_base_block_if_cyclic() ii != 0");
	return OK;
}


#endif /* MATRIXTRUE */

