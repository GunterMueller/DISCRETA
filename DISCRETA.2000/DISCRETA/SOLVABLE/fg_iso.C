/* fg_iso.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>
#ifndef LABRA_INCLUDED
#include <DISCRETA/lb.h>
#endif

#define MAX_NW 64

#define FG_ISO_USE_CANONIC_FORM
#define USE_COLORING

#undef FG_ISO_USE_CANONIC_TG

#undef FG_ISO_DEBUG_CANONICIZE

#if TEXDOCU
INT FG_OB::int2rep(INT ii, INT *coset_reps)
#endif
{
	INT i, j;

	for (i = s_nb_ze_i() - 1; i >= 0; i--) {
		j = ii % s_T_i(i)->s_li();
		coset_reps[i] = s_T_iji(i, j);
		ii -= j;
		ii /= s_T_i(i)->s_li();
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::rep2int(INT *coset_reps, INT *ii)
#endif
{
	INT i, j, k;
	
	*ii = 0;
	for (i = 0; i < s_nb_ze_i(); i++) {
		(*ii) *= s_T_i(i)->s_li();
		k = coset_reps[i];
		for (j = 0; j < s_T_i(i)->s_li(); j++) {
			if (s_T_iji(i, j) == k)
				break;
			}
		if (j == s_T_i(i)->s_li()) {
			return error("FG_OB::rep2int(): coset_reps[i] not found");
			}
		(*ii) += j;
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::rep2aut(INT *coset_reps, PERMUTATION_OP aut)
#endif
{
	PERMUTATION_OB p1;
	PERMUTATION_OP q;
	INT i, j;

	aut->m_il(s_n_i());
	aut->one();
	for (i = 0; i < s_nb_ze_i(); i++) {
		j = coset_reps[i];
		q = s_AutM_ij(i, j);
		q->mult(aut, &p1);
		/* aut->mult(q, &p1); */
			/* q wird vorgeschaltet ! */
		p1.swap(aut);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::aut2rep(PERMUTATION_OP aut, INT *coset_reps)
#endif
{
	ZE_OP ze;
	PERMUTATION_OB p, q, rv;
	PERMUTATION_OP r;
	INT i, j, g_i;

	aut->copy(&p);
	for (i = 0; i < s_nb_ze_i(); i++) {
		ze = s_ze_i(i);
		g_i = ze->s_n0_i();
		j = p.s_ii(g_i) - 1;
		coset_reps[i] = j;
		r = s_AutM_ij(i, j);
		r->invers(&rv);
		p.mult(&rv, &q);
		/* rv.mult(&p, &q); */
			/* rv nachgeschaltet ! */
		q.swap(&p);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::int2aut(INT ii, PERMUTATION_OP aut)
#endif
{
	INT cr[MAX_NW];
	
	int2rep(ii, cr);
	rep2aut(cr, aut);
	return OK;
}

#if TEXDOCU
INT FG_OB::aut2int(PERMUTATION_OP aut, INT *ii)
#endif
{
	INT cr[MAX_NW];
	
	aut2rep(aut, cr);
	rep2int(cr, ii);
	return OK;
}

#if TEXDOCU
INT FG_OB::aut2bi(PERMUTATION_OP aut, INT *base_im)
#endif
{
	ZE_OP ze;
	INT i, j, g_i;

	for (i = 0; i < s_nb_ze_i(); i++) {
		ze = s_ze_i(i);
		g_i = ze->s_n0_i();
		j = aut->s_ii(g_i) - 1;
		base_im[i] = j;
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::rep2bi(INT *coset_reps, INT *base_im)
#endif
{
	PERMUTATION_OB aut;

	rep2aut(coset_reps, &aut);
	aut2bi(&aut, base_im);
	return OK;
}

#if TEXDOCU
INT FG_OB::int2bi(INT ii, INT *base_im)
#endif
{
	PERMUTATION_OB aut;

	int2aut(ii, &aut);
	aut2bi(&aut, base_im);
	return OK;
}

#if TEXDOCU
INT FG_OB::bi2rep(INT *base_im, INT *coset_rep)
#endif
{
	INT bi[MAX_NW];
	PERMUTATION_OB rv;
	PERMUTATION_OP r;
	INT i, j, ii;
	
	for (i = 0; i < s_nb_ze_i(); i++)
		bi[i] = base_im[i];
	for (i = 0; i < s_nb_ze_i(); i++) {
		j = bi[i];
		coset_rep[i] = j;
		r = s_AutM_ij(i, j);
		r->invers(&rv);
		for (ii = i + 1; ii < s_nb_ze_i(); ii++)
			bi[ii] = rv.s_ii(bi[ii]) - 1;
		}
	return OK;
	
}

#if TEXDOCU
INT FG_OB::bi2int(INT *base_im, INT *ii)
#endif
{
	INT cr[MAX_NW];

	bi2rep(base_im, cr);
	rep2int(cr, ii);
	return OK;
}

#if TEXDOCU
INT FG_OB::test_bi()
#endif
{
	PERMUTATION_OB aut;
	INT bi[MAX_NW];
	INT ago, i, j;
	
	if (s_ago()->s_obj_k() != INTEGER)
		return error("FG_OB::test_bi() ago not an INTEGER");
	ago = ((INTEGER_OP) s_ago())->s_i();
	for (i = 0; i < ago; i++) {
		int2bi(i, bi);
		int2aut(i, &aut);
		if (i == 0 && !aut.einsp())
			printf("FG_OB::test_bi(): i == 0 && !aut.einsp()\n");
		bi2int(bi, &j);
		if (j != i)
			return error("FG_OB::test_bi(): i != j nach bi2int()");
		aut2int(&aut, &j);
		if (j != i)
			return error("FG_OB::test_bi(): i != j nach aut2int()");
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::autlog_info_fill_in_colors(FILE *fp, AUTLOG_INFO *info, INT f_B)
#endif
{
	VECTOR_OP col = (VECTOR_OP) callocobject("fg_iso.C autlog_info_fill_in_colors");
	VECTOR_OB val, mult;
	BYTE str[10000];
		
	Coloring_seldoms_prefered(col);
	col->multiplicities(&val, &mult);
	str[0] = 0;
	val.sprint_multiplicities(&mult, str);
	fprintf(fp, "coloring: %s\n", str);
	fflush(fp);
	
	if (f_B) {
		info->f_B_has_colors = TRUE;
		info->B_colors = col;
		}
	else {
		info->f_A_has_colors = TRUE;
		info->A_colors = col;
		}	
	return OK;
}

#if TEXDOCU
INT FG_OB::autlog_info_fill_in_theG(AUTLOG_INFO *info, INT f_B)
#endif
{
	MATRIX_OP G;
	INT *theG = NIL;
	ZE_OP ze;
	INT n, i, j, size, dim_n1, len, l, g;
	
	/* fill in the group table: */
	n = s_n_i();
	dim_n1 = s_dim_n1_i();
	size = n * dim_n1 * sizeof(INT);
	theG = (INT *) my_malloc(size, "autlog_info_fill_in_theG");
	if (theG == NIL)
		return error("autlog_info_fill_in_theG(): no memory");
	G = s_theG();
	for (i = 0; i < n; i++) {
		for (j = 0; j < dim_n1; j++) {
			theG[i * dim_n1 + j] = G->s_iji(i, j);
			}
		}
	if (f_B) {
		info->Bn = n;
		info->BtheG = theG;
		info->Bdim_n = dim_n1;
		info->f_B_has_colors = FALSE;
		info->B_colors = NIL;
		}
	else {
		info->An = n;
		info->AtheG = theG;
		info->Adim_n = dim_n1;
		info->f_A_has_colors = FALSE;
		info->A_colors = NIL;
		}	

	/* fill in the generator indices: */
	len = s_nb_ze_i();
	if (len >= AUTLOG_MAX_G) {
		Srfs("FG::autlog_info_fill_in_theG", "len >= AUTLOG_MAX_G");
		return ERROR;
		}
	for (l = 0; l < len; l++) {
		ze = s_ze_i(l);
		g = ze->s_n0_i();
		if (f_B)
			info->B_g[l] = g;
		else
			info->A_g[l] = g;
		}
	if (f_B)
		info->B_nb_gen = len;
	else
		info->A_nb_gen = len;
	return OK;
}

#if TEXDOCU
INT FG_OB::GmZ_Gd_orders(FILE *fp_txt, INT f_verbose)
#endif
{
	AUTLOG_INFO info;
	INT GmZ_n, Gd_n;

	info.fp_txt = stdout;
	info.AtheG = NIL;
	autlog_info_fill_in_theG(&info, FALSE /* f_B */);
	info.f_verbose = f_verbose;
		/* to print the orders */
	info.f_very_verbose = FALSE;
	if (!autlog_GmZ_Ad_orders(&info, &GmZ_n, &Gd_n, fp_txt)) {
		Srff("FG::GmZ_Gd_orders", "autlog_GmZ_Ad_orders");
		return ERROR;
		}
	s_GmZ_n()->m_i(GmZ_n);
	s_Gd_n()->m_i(Gd_n);

	my_ptr_free( (void **) &info.AtheG);
	return OK;
}

#if TEXDOCU
INT FG_OB::Isomorphic(FG_OP fg2, INT f_verbose)
#endif
{
#ifdef FG_ISO_USE_CANONIC_FORM
	return Isomorphic_by_canonic_form(fg2, f_verbose);
#else
	return Isomorphic2(fg2, f_verbose);
#endif
}

#if TEXDOCU
INT FG_OB::Isomorphic2(FG_OP fg2, INT f_verbose)
#endif
{
	AUTLOG_INFO info;
	INT erg = OK, ret, n;
	
	info.fp_txt = stdout;
	info.AtheG = NIL;
	info.BtheG = NIL;
	n = s_n_i();
	if (n != fg2->s_n_i())
		return error("FG::Isomorphic()|different group orders!\n");
	if (s_dim_n1_i() != fg2->s_dim_n1_i())
		return error("FG::Isomorphic()|different dim_n1!\n");
	
	autlog_info_fill_in_theG(&info, FALSE /* f_B */);
	fg2->autlog_info_fill_in_theG(&info, TRUE /* f_B */);

	info.nb_isomorphisms = 0;
	info.mode = AUT_MODE_BREAK_AFTER_FST;
	info.f_autologisms = FALSE;
	info.f_verbose = FALSE;
	info.f_very_verbose = FALSE;

	info.add_aut = NIL;
	info.add_aut_perm = NIL;
	info.data = NIL;
	if (!autlog_test(&info)) {
		Srff("FG::isomorphic", "autlog_test");
		return ERROR;
		}
	if (info.nb_isomorphisms) {
		if (f_verbose)
			printf("FG::Isomorphic()|groups are isomorphic\n");
		ret = TRUE;
		}
	else {
		if (f_verbose)
			printf("FG::Isomorphic()|"
			"non isomorphic groups\n");
		ret = FALSE;
		}

l_exit:
	my_ptr_free( (void **) &info.AtheG);
	my_ptr_free( (void **) &info.BtheG);
	return ret;
}

#if TEXDOCU
INT FG_OB::Isoclinic(FG_OP fg2, INT f_verbose, INT f_very_verbose)
#endif
{
	AUTLOG_INFO info;
	INT erg = OK, ret;
	
	info.fp_txt = stdout;
	info.AtheG = NIL;
	info.BtheG = NIL;
	if (s_GmZ_n_i() != fg2->s_GmZ_n_i())
		return error("FG::Isoclinic()|different GmZ group orders!\n");
	if (s_Gd_n_i() != fg2->s_Gd_n_i())
		return error("FG::Isoclinic()|different Gd group orders!\n");
	
	autlog_info_fill_in_theG(&info, FALSE /* f_B */);
	fg2->autlog_info_fill_in_theG(&info, TRUE /* f_B */);

	info.nb_isomorphisms = 0;
	info.mode = AUT_MODE_BREAK_AFTER_FST;
	info.f_autologisms = TRUE;
	info.f_verbose = f_verbose;
	info.f_very_verbose = f_very_verbose;

	info.add_aut = NIL;
	info.add_aut_perm = NIL;
	info.data = NIL;
	if (!autlog_test(&info)) {
		Srff("FG::Isoclinic", "autlog_test");
		return ERROR;
		}
	if (info.nb_isomorphisms) {
		if (f_verbose)
			printf("FG::Isoclinic()|groups are isoclinic\n");
		ret = TRUE;
		}
	else {
		if (f_verbose)
			printf("FG::Isoclinic()|non isoclinic groups\n");
		ret = FALSE;
		}

l_exit:
	my_ptr_free( (void **) &info.AtheG);
	my_ptr_free( (void **) &info.BtheG);
	return ret;
}

#if TEXDOCU
INT FG_OB::find_isomorphic_group(VECTOR_OP V, INT f_verbose)
#else
/* Die Gruppenordnungen der Gruppen 
 * in V muessen stimmen. */
#endif
{
	FG_OP G2;
	INT i, l;
	
	l = V->s_li();
	printf("find_isomorphic_group(%ld:", l);
	fflush(stdout);
	for (i = 0; i < V->s_li(); i++) {
		G2 = (FG_OP) V->s_i(i);
		if (Isomorphic(G2, f_verbose)) {
			printf("isomorphic !)\n");
			fflush(stdout);
			return i;
			}
		printf(".");
		fflush(stdout);
		}
	printf(")\n");
	fflush(stdout);
	return -1;
}

#if TEXDOCU
INT FG_OB::find_isoclinic_group(VECTOR_OP V, INT f_verbose)
#endif
{
	FG_OP G2;
	INT i, GmZ_n, Gd_n;
	
	GmZ_n = s_GmZ_n_i();
	Gd_n = s_Gd_n_i();
	if (f_verbose) {
		printf("find_isoclinic_group(): GmZ_n = %ld Gd_n = %ld\n", GmZ_n, Gd_n);
		fflush(stdout);
		}
	for (i = 0; i < V->s_li(); i++) {
		G2 = (FG_OP) V->s_i(i);
		if (G2->s_GmZ_n_i() != GmZ_n)
			continue;
		if (G2->s_Gd_n_i() != Gd_n)
			continue;
		if (f_verbose) {
			printf("find_isoclinic_group() i = %ld, vor Isoclinic()\n", i);
			fflush(stdout);
			}
		if (Isoclinic(G2, FALSE /* f_verbose */, 
			FALSE /* f_very_verbose */))
			return i;
		}
	return -1;
}

static INT add_aut_perm(AUTLOG_INFO *info, INT nb_gen, INT *gen, PERMUTATION_OP p);

#undef DEBUG_AUT

#if TEXDOCU
INT FG_OB::Aut(INT f_verbose, INT f_very_verbose)
#endif
{
	MATRIX_OB M;
	PERMUTATION_OP perm;
	AUTLOG_INFO info;
	SYM_OB ago, b;
	INT n, i, j, l, len;
	INT erg = OK;
	
	info.fp_txt = stdout;
	info.AtheG = NIL;
	info.BtheG = NIL;
	n = s_n_i();
	s_ago()->freeself();
	((INTEGER_OP) s_ago())->m_i(0);
	autlog_info_fill_in_theG(&info, FALSE /* f_B */);
	autlog_info_fill_in_theG(&info, TRUE /* f_B */);
	
	info.nb_isomorphisms = 0;
	info.mode = AUT_MODE_ONLY_COSET_REPS;
	info.f_autologisms = FALSE;
	info.f_verbose = FALSE;
	info.f_very_verbose = FALSE;
#ifdef DEBUG_AUT
	info.f_verbose = TRUE;
	info.f_very_verbose = TRUE;
#endif
	
	info.add_aut = NIL;
	info.add_aut_perm = add_aut_perm;
	len = s_nb_ze_i();
	M.m_ilih(n, len);
	info.data = (void *) &M;
	if (!autlog_test(&info)) {
		Srff("FG::Aut", "autlog_test");
		return ERROR;
		}
	s_T()->m_il(len);
	/* calc aut_group_order */
	/* ago = 1; */
	((INTEGER_OP) &ago)->m_i(1);
	for (i = 0; i < len; i++) {
		if (f_very_verbose)
			printf("orbit of %ld:\n", info.A_g[i]);
		s_T_i(i)->m_il(0);
		l = 0;
		for (j = 0; j < n; j++) {
			if (M.s_ij(i, j)->s_obj_k() != EMPTY) {
				perm = (PERMUTATION_OP) M.s_ij(i, j);
				if (f_very_verbose)
					printf("%ld ", j);
				/* perm->println(); */
				s_T_i(i)->inc();
				s_T_ij(i, l)->m_i(j);
				l++;
				}
			}
		/* ago *= l; */
		((INTEGER_OP) &b)->m_i(l);
		b.mult_apply(&ago);
		if (f_very_verbose) {
			printf("length %ld ago = ", l);
			ago.println();
			}
		}
	/* s_ago()->m_i(ago); */
	ago.swap(s_ago());
	s_f_has_aut()->m_i(TRUE);
	M.swap(s_AutM());
	if (f_verbose) {
		printf("aut_group_order: ");
		s_ago()->println();
		}
	if (info.AtheG) {
		my_free(info.AtheG);
		info.AtheG = NIL;
		}
	if (info.BtheG) {
		my_free(info.BtheG);
		info.BtheG = NIL;
		}
	return erg;
}

#if TEXDOCU
static INT add_aut_perm(AUTLOG_INFO *info, INT nb_gen, INT *gen, PERMUTATION_OP p)
#endif
{
	MATRIX_OP M = (MATRIX_OP) info->data;
	INT base, i, j;
	
#ifdef DEBUG_AUT
	printf("add_aut_perm():");
	p->println();
	fflush(stdout);
#endif
	for (i = 0; i < nb_gen; i++) {
		base = gen[i];
		j = p->s_ii(base) - 1;
		if (M->s_ij(i, j)->s_obj_k() == EMPTY) {
			((SYM_OP) p)->copy(
				M->s_ij(i, j));
			}
		if (j != base)
			break;
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::get_aut_generators(VECTOR_OP V)
#endif
{
	VECTOR_OP Ti;
	INT i, j, k, l, Tlen;
	
	if (!s_f_has_aut_i())
		return error("FG_OB::get_aut_generators: !s_f_has_aut_i()");
	l = 0;
	Tlen = s_T()->s_li();
	for (i = 0; i < Tlen; i++) {
		l += s_T_i(i)->s_li();
		}
	V->m_il(l);
	l = 0;
	for (i = 0; i < Tlen; i++) {
		Ti = s_T_i(i);
		for (j = 0; j < Ti->s_li(); j++) {
			k = Ti->s_ii(j);
			s_AutM_ij(i, k)->copy((PERMUTATION_OP) V->s_i(l));
			l++;
			}
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::aut_on_generators(INT f_v, INT f_vv)
#endif
{
	GROUP_TABLE table;
	FILE *fp_txt = stdout;
	VECTOR_OB base;
	INT i, gi, l;

	group_init_from_fg(&table, this);
	if (f_vv) {
		fprintf(fp_txt, "FG::aut_on_generators() computing canonical form...\n");
		fflush(fp_txt);
		}
	l = s_nb_ze_i();
	base.m_il(l);
	for (i = 0; i < l; i++) {
		gi = g_i(i);
		base.m_ii(i, gi);
		}
	return error("FG::aut_on_generators() group_calc_aut not available at the moment");
#if 0
	group_calc_aut(fp_txt, &table, &base, 
		s_AutM(), s_T(), s_ago(), f_v, f_vv);
#endif
	group_free(&table);
	s_f_has_aut()->m_i(TRUE);
	return OK;
}

#if TEXDOCU
INT FG_OB::canonicize(PERMUTATION_OP p0, PERMUTATION_OP p0v, INT f_v, INT f_vv)
#endif
{
	GROUP_TABLE table;
	GROUP_CANONIC_FORM_OB cf;
	FILE *fp_txt = stdout;

#if 0
	if (canonicize_group_cf(stdout, this, &cf, f_v, f_vv) != OK)
		return error("FG::canonicize() error in canonicize_group_cf()");
#endif

	group_init_from_fg(&table, this);
	if (f_vv) {
		fprintf(fp_txt, "FG::canonicize() computing canonical form...\n");
		fflush(fp_txt);
		}
	group_calc_cf(fp_txt, &table, &cf, f_v, f_vv);
	if (f_vv) {
		fprintf(fp_txt, "FG::canonicize() computing hash value...\n");
		fflush(fp_txt);
		}
	cf.calc_hash(fp_txt, s_hash(), s_n_i(), f_v, f_vv);
	s_f_has_hash()->m_i(TRUE);
	cf.s_p0()->copy(s_p0());
	s_f_has_p0()->m_i(TRUE);

	cf.s_p0()->copy(p0);
	cf.s_p0v()->copy(p0v);
	if (f_vv) {
		fprintf(fp_txt, "finished!\n");
		fflush(fp_txt);
		}
	group_free(&table);
	return OK;
}

#if TEXDOCU
INT FG_OB::canonicize_with_generators(VECTOR_OP gen, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, INT f_v)
#endif
{
	GROUP_TABLE table;
	GROUP_CANONIC_FORM_OB cf;
	INT f_vv = FALSE;
	FILE *fp_txt = stdout;

	group_init_from_fg(&table, this);
	group_init_generators(&table, gen);
	if (f_vv) {
		fprintf(fp_txt, "FG::canonicize() computing canonical form...\n");
		fflush(fp_txt);
		}
	group_calc_cf(fp_txt, &table, &cf, f_v, f_vv);
	if (f_vv) {
		fprintf(fp_txt, "FG::canonicize() computing hash value...\n");
		fflush(fp_txt);
		}
	cf.calc_hash(fp_txt, s_hash(), s_n_i(), f_v, f_vv);
	s_f_has_hash()->m_i(TRUE);
	cf.s_p0()->copy(s_p0());
	s_f_has_p0()->m_i(TRUE);

	cf.s_p0()->copy(p0);
	cf.s_p0v()->copy(p0v);
	if (f_vv) {
		fprintf(fp_txt, "finished!\n");
		fflush(fp_txt);
		}
	group_free(&table);
	return OK;
}

#if TEXDOCU
INT FG_OB::canonicize_qad(FILE *fp_txt, INT f_v, INT f_vv)
#endif
{
	PERMUTATION_OB p, pv;

#ifdef FG_ISO_DEBUG_CANONICIZE
	f_v = TRUE;
	fprintf(fp_txt, "FG::canonicize_qad() calling canonicize_group_qad()\n");
	fflush(fp_txt);
#endif

#if 0
	if (canonicize_group_qad(fp_txt, this, f_v, f_vv) != OK)
		return error("FG::canonicize_qad() error in canonicize_group_qad()");
#endif
	if (canonicize(&p, &pv, f_v, f_vv) != OK)
		return error("FG::canonicize_qad() error in canonicize()");

#ifdef FG_ISO_DEBUG_CANONICIZE
	printf("FG::canonicize_qad() the group:\n");
	Print();
	fflush(stdout);
#endif
	return OK;
}

#if TEXDOCU
INT FG_OB::canonicize_qad_with_generators(FILE *fp_txt, VECTOR_OP gen, INT f_v)
#endif
{
	PERMUTATION_OB p, pv;
	INT f_vv = FALSE;

	if (canonicize_with_generators(gen, &p, &pv, f_v) != OK)
		return error("FG::canonicize_qad_with_generators() "
			"error in canonicize_with_generators()");
	return OK;
}

#if TEXDOCU
INT FG_OB::Isomorphic_by_canonic_form(FG_OP G2, INT f_verbose)
#endif
{
	PERMUTATION_OB p1, p2;
	PERMUTATION_OB p1v, p2v;
	MATRIX_OP theG1, theG2;
	INT n, i, j, i1, i2, j1, j2, a, b, aa, bb;
	
	n = s_n_i();
	if (n != G2->s_n_i())
		return error("FG::Isomorphic_by_canonic_form()|different group orders!\n");
	if (!s_f_has_p0_i()) {
		canonicize(&p1, &p1v, f_verbose, FALSE);
		}
	else {
		s_p0()->copy(&p1);
		p1.invers(&p1v);
		}
	if (!G2->s_f_has_p0_i()) {
		G2->canonicize(&p2, &p2v, f_verbose, FALSE);
		}
	else {
		G2->s_p0()->copy(&p2);
		p2.invers(&p2v);
		}

	theG1 = s_theG();
	theG2 = G2->s_theG();
	for (i = 0; i < n; i++) {
		i1 = p1v.s_ii(i) - 1;
		i2 = p2v.s_ii(i) - 1;
		for (j = 0; j < n; j++) {
			j1 = p1v.s_ii(j) - 1;
			j2 = p2v.s_ii(j) - 1;
			a = theG1->s_iji(i1, j1);
			b = theG2->s_iji(i2, j2);
			aa = p1.s_ii(a) - 1;
			bb = p2.s_ii(b) - 1;
			if (aa != bb)
				return FALSE;
			}
		}
	return TRUE;
}

#ifdef FG_ISO_USE_CANONIC_TG

#if TEXDOCU
INT FG_OB::canonic(PERMUTATION_OP p)
#else
uses TG routines.
#endif
{
	GROUP_CANON GC;
	ARRAY < VEKTOR < short > > tafel;
	INT n, i, j, a;
	INT t0, t1, user_time;
	BYTE str[256];
	
	n = s_n_i();
	tafel.REALLOC(n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = s_theG_iji(i, j);
			tafel[i + 1][j + 1] = a + 1;
			}
		}
	printf("FG::canonic(): calling group_can()...\n");
	fflush(stdout);
	
	t0 = os_ticks();
	GC.Init(n, &tafel);
	GC.Run();
	VEKTOR < short > & p_ = GC.GetCanNum();
	p->m_il(n);
	for (i = 0; i < n; i++) {
		a = p_[i + 1];
		p->m_ii(i, a);
		}
	t1 = os_ticks();
	user_time = t1 - t0;
	str[0] = 0;
	print_delta_time(user_time, str);
	printf("found canonic form %s\n", str);
	fflush(stdout);
	return OK;
}

#endif /* FG_ISO_USE_CANONIC_TG */


#endif /* SOLVABLE_TRUE */

