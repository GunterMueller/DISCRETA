/* aut.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/perm.h>

#include <DISCRETA/solvable.h>

#ifndef LABRA_INCLUDED
#include <DISCRETA/lb.h>
#endif

#include <DISCRETA/SOLVABLE/autP.h>
/* private include file ! */


#define GROUP_CANONICIZE_CALC_HASH

#undef PRINT_OCCASIONALLY
#define MIN_LEVEL_FOR_PRINT 3
#define PRINT_INTERVALL_SEC 5

#define PRINT_BACKTRACK_POINTS
#define BACKTRACK_MOD 250

static INT last_print_ticks, nb_backtrack, nb_backtrack_points;

#if TEXDOCU
INT canonicize_group_qad(FILE *fp_txt, FG_OP G, INT f_v, INT f_vv)
#else
/* qad stands for "`quick and dirty"' ! */
#endif
{
	GROUP_CANONIC_FORM_OB cf;

	if (canonicize_group_cf(fp_txt, G, &cf, f_v, f_vv) != OK)
		return f_error(fp_txt, "canonicize_group_qad() error in canonicize_group_cf()");
	return OK;
}

#if TEXDOCU
INT canonicize_group_cf(FILE *fp_txt, FG_OP G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv)
#else
#endif
{
	LABRA_OB L;
	SYM_OB ago;
	INT ago_trunc;
	INT f_aut = TRUE;
	INT f_aut_gens = TRUE;
	INT n, nb_gen;
	/* VECTOR_OB hash; */

	n = G->s_n_i();
	cf->init();
	if (canonicize_group(fp_txt, G, cf->s_p0(), cf->s_p0v(), 
		cf->s_go(), cf->s_g(), 
		cf->s_ro(), cf->s_re(), cf->s_GG(), 
		cf->s_auts(), cf->s_first_moved(), f_aut_gens, 
		&L, &ago, f_aut, f_v, f_vv) != OK)
		return f_error(fp_txt, "canonicize_group_cf() error in canonicize_group()");
	nb_gen = cf->s_ro()->s_li();
	cf->s_nb_gen()->m_i(nb_gen);
	if (f_vv) {
		printf("canonicize_group_cf():\n");
		printf("ago = ");
		fflush(stdout);
		ago.println();
		printf("nb_gen = ");
		cf->s_nb_gen()->println();
		printf("g = ");
		cf->s_g()->println();
		printf("go = ");
		cf->s_go()->println();
		printf("ro = ");
		cf->s_ro()->println();
		printf("re = ");
		cf->s_re()->println();
		printf("GG = \n");
		cf->s_GG()->Print();
		fflush(stdout);
		}
	cf->s_ago_trunc()->m_i(0);
	if (f_aut) {
		if (ago.s_obj_k() != INTEGER || 
			((INTEGER_OP) &ago)->s_i() > GROUP_CANONIC_FORM_MAX_AGO) 
			ago_trunc = GROUP_CANONIC_FORM_MAX_AGO;
		else
			ago_trunc = ((INTEGER_OP) &ago)->s_i();
		cf->s_ago_trunc()->m_i(ago_trunc);
		if (f_v) {
			fprintf(fp_txt, "ago_trunc = %ld\n", ago_trunc);
			fflush(fp_txt);
			}
		}
	if (f_aut_gens) {
		cf->s_f_has_auts()->m_i(TRUE);
		cf->s_f_auts_on_canonic_form()->m_i(FALSE);
		if (f_vv) {
			printf("has %ld auts\n", cf->s_auts()->s_li());
			printf("first_moved = ");
			cf->s_first_moved()->println();
			fflush(stdout);
			}
		}
#ifdef GROUP_CANONICIZE_CALC_HASH
	if (f_vv) {
		fprintf(fp_txt, "computing hash value...\n");
		fflush(fp_txt);
		}
	cf->calc_hash(fp_txt, G->s_hash(), n, f_v, f_vv);
	G->s_f_has_hash()->m_i(TRUE);
	if (f_vv) {
		fprintf(fp_txt, "finished!\n");
		fflush(fp_txt);
		}

	cf->s_p0()->copy(G->s_p0());
	G->s_f_has_p0()->m_i(TRUE);
#if 0
	if (f_v) {
		printf("canonical labelling p0 = \n");
		G->s_p0()->println();
		fflush(stdout);
		}
#endif
#endif
	return OK;
}

#if TEXDOCU
INT canonicize_group(FILE *fp_txt, FG_OP G, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut, 
	INT f_v, INT f_vv)
#endif
{
	AUTLOG_INFO info;
	INT t0, t1, user_time;
	BYTE str[256];

	if (f_v) {
		/* print_size(); */
		fprintf(fp_txt, "(");
		fflush(fp_txt);
		}

	t0 = os_ticks();
	info.fp_txt = fp_txt;
	G->autlog_info_fill_in_theG(&info, FALSE /* f_B */ );
	/* computes:
		info->An
		info->AtheG
		info->Adim_n
		info->A_nb_gen = s_nb_ze_i();
			info->A_g[l]
		info->f_A_has_colors = FALSE, 
		info->A_colors = NIL
	*/
	G->autlog_info_fill_in_colors(fp_txt, &info, FALSE /* f_B */);
		/* computes the coloring and fills it into info */

	info.nb_isomorphisms = 0;
	info.mode = AUT_MODE_BREAK_AFTER_FST;
	info.f_autologisms = FALSE;
	info.f_verbose = f_v;
	info.f_very_verbose = f_vv;
	if (f_vv) {
		printf("f_verbose = %ld f_very_verbose = %ld\n", f_v, f_vv);
		}

	info.add_aut = NIL;
	info.add_aut_perm = NIL;
	info.data = NIL;
	
	if (!autlog_canonicize(&info, p0, p0v, go, g, Ro, Re, GG, 
		auts, first_moved, f_aut_gens, Aut, ago, f_aut)) {
		return error("canonicize_group() error in autlog_canonicize()");
		}

	t1 = os_ticks();
	user_time = t1 - t0;
	if (f_v) {
		fprintf(fp_txt, "nb_backtrack = %ld\n", nb_backtrack);
		str[0] = 0;
		print_delta_time(user_time, str);
		fprintf(fp_txt, "user time:  %s", str);
		fprintf(fp_txt, ")\n");
		fflush(fp_txt);
		}
#if 0
	if (info.nb_isomorphisms) {
		if (f_verbose)
			printf("FG::Isomorphic()|"
			"groups are isomorphic\n");
		ret = TRUE;
		}
	else {
		if (f_verbose)
			printf("FG::Isomorphic()|"
			"non isomorphic groups\n");
		ret = FALSE;
		}
#endif

l_exit:
	my_ptr_free( (void **) &info.AtheG);
	if (info.f_A_has_colors) {
		freeall(info.A_colors);
		info.A_colors = NIL;
		info.f_A_has_colors = FALSE;
		}
	return OK;
}

#if TEXDOCU
static INT autlog_canonicize(AUTLOG_INFO *info, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut)
#endif
{
	AUTLOG_LOCAL *L = NIL;
	INT n, i, j, nb_gen;
	
	L = (AUTLOG_LOCAL *) my_malloc(sizeof(AUTLOG_LOCAL), "autlog_canonicize L");
	if (L == NIL)
		return error("autlog_canonicize(): no memory for L");
	autlog_nil(L);
	L->info = info;
	n = L->info->An;
	L->A.n = n;
	autlog_the_group_int(&L->A);
		/* allocates p, pv, q, qv, tmp1, tmp2 */
	L->A.fp_txt = info->fp_txt;
	if (info->f_A_has_colors) {
		autlog_the_group_init_colors(&L->A, info->A_colors);
		}
	
	L->A.theG = info->AtheG;
	L->A.dim_n = info->Adim_n;
	L->f_going_back = FALSE;
	
	autlog_the_group_open_grid(&L->A, 1 /* nb_grid */, L->A.n /* n */);
	if (info->f_very_verbose) {
		printf("autlog_canonicize(): "
		"nach autlog_the_group_open_grid(A) n = %ld\n", L->A.n);
		fflush(stdout);
		}

	L->pA = &L->A;

	L->E = (ORDERED_SET *) my_malloc(L->pA->n * sizeof(ORDERED_SET), "autlog_canonicize E");
	if (L->E == NIL)
		return error("autlog_canonicize(): no memory for L->E");


	L->is_first = TRUE;
	L->first_moved = AUTLOG_MAX_G;
	sp_nil(&L->p0);
	sp_nil(&L->p0v);
	sp_int(&L->p0, L->pA->n, "autlog_canonicize");
	sp_int(&L->p0v, L->pA->n, "autlog_canonicize");
	sp_id(&L->p0);
	sp_id(&L->p0v);
	L->nb_aut_gens = 0;
	L->dim_aut_gens = 0;
	L->aut_gens = NIL;
	L->auts_first_moved = NIL;
	
	last_print_ticks = os_ticks();
	nb_backtrack = 0;
	nb_backtrack_points = 0;
	
	canonicize(L, 0 /* level */);

	/* get the canonical labelling: */
	p0->m_il(n);
	p0v->m_il(n);
	for (i = 0; i < n; i++) {
		j = L->p0.a[i];
		p0->m_ii(i, j + 1);
		j = L->p0v.a[i];
		p0v->m_ii(i, j + 1);
		}
	
	/* get the generators and orders 
	 * of the canonical subgroup chain */
	nb_gen = L->nb_gen;
	go->m_il(nb_gen + 1);
	g->m_il(nb_gen + 1);
	for (i = 0; i <= nb_gen; i++) {
		j = L->go[i];
		go->m_ii(i, j);
		j = L->g[i];
		g->m_ii(i, j);
		}

	if (f_aut_gens) {
		get_aut_gens(L, auts, first_moved);
		}
	if (f_aut) {
		get_aut_group(L, Aut, ago);
		}
	recalc_Ro_Re_GG(L, Ro, Re, GG);

	sp_free(&L->p0);
	sp_free(&L->p0v);
	free_aut_gen(L);

	/* from autlog_ext(): */
	/* Kopien auf info, nicht freigeben: */
	L->A.theG = NIL;
	L->Ad.theG = NIL;
	my_ptr_free( (void **) &L->E);
	autlog_the_group_ext(&L->A);
	alfg_ext(&L->A_Z);
	autlog_the_group_ext(&L->AmZ);
		/* frees p, pv, q, qv, tmp1, tmp2 
		 * frees Grid if allocated. */
	alcs_ext(&L->AcA);
	autlog_the_group_ext(&L->Ad);


	my_free(L);
	return TRUE;
}

#if TEXDOCU
static INT recalc_Ro_Re_GG(AUTLOG_LOCAL *L, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG)
#endif
{
	INT nb_gen, i, j, g_i, g_j, k, k_, n;
	INT go_last, go, g, ro, re, re_;
	INT *theG, dim_n;

	n = L->pA->n;
	nb_gen = L->nb_gen;
	theG = L->pA->theG;
	dim_n = L->pA->dim_n;
	Ro->m_il(nb_gen);
	Re->m_il(nb_gen);
	GG->m_ilih(nb_gen, nb_gen);
	for (i = 0; i < nb_gen; i++) {
		g = L->g[i + 1];
		go_last = L->go[i];
		if (L->info->f_very_verbose) {
			printf("recalc_Ro_Re_G(): i = %ld go_last = %ld g = %ld\n", i, go_last, g);
			fflush(stdout);
			}
		ro = 1;
		re = g;
		while (L->p0.a[re] >= go_last) {
			re = theG[re * dim_n + g];
			ro++;
			}
		re_ = L->p0.a[re];
		Ro->m_ii(i, ro);
		Re->m_ii(i, re_);
		if (L->info->f_very_verbose) {
			printf("ro = %ld re = %ld re_ = %ld\n", ro, re, re_); 
			fflush(stdout);
			}
		}
	for (i = 0; i < nb_gen; i++) {
		g_i = L->g[i + 1];
		for (j = 0; j < nb_gen; j++) {
			g_j = L->g[j + 1];
			k = theG[g_i * dim_n + g_j];
			k_ = L->p0.a[k];
			GG->m_iji(i, j, k_);
			}
		}
	if (L->info->f_very_verbose) {
		printf("GG:\n");
		GG->Print();
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
static INT get_aut_gens(AUTLOG_LOCAL *L, VECTOR_OP auts, VECTOR_OP first_moved)
#endif
{
	SPERM *p, *q;
	PERMUTATION_OP perm;
	INT i, n, nb_gen, ii, k;

	n = L->pA->n;
	nb_gen = L->nb_aut_gens;
	auts->m_il(nb_gen);
	first_moved->m_il(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		p = L->aut_gens + i;
		q = p;
		k = L->auts_first_moved[i];
		first_moved->m_ii(i, k);
		perm = (PERMUTATION_OP) auts->s_i(i);
		perm->m_il(n);
		for (ii = 0; ii < n; ii++) {
			k = q->a[ii];
			perm->m_ii(ii, k + 1);
			}
		}
	return OK;
}

#if TEXDOCU
static INT get_aut_group(AUTLOG_LOCAL *L, LABRA_OP Aut, SYM_OP ago)
#endif
{
	SPERM *p, *q;
	SPERM tmp1, tmp2;
	PERMUTATION_OP KM_p;
	INT n, i, nb_gen, j_, j, k, fst_moved;
	INT ii, jj;
	INT f_v = TRUE, f_vv = FALSE;
	FILE *fp_txt;
	BYTE str[1024];

	fp_txt = L->info->fp_txt;
	n = L->pA->n;
	Aut->Init_no_generators(n, f_v, f_vv);
	nb_gen = L->nb_aut_gens;
	sp_nil(&tmp1);
	sp_nil(&tmp2);
	sp_int(&tmp1, n, "aut.C: get_aut_group");
	sp_int(&tmp2, n, "aut.C: get_aut_group");
	for (i = 0; i < nb_gen; i++) {
		p = L->aut_gens + i;
		q = p;
		sp_mult(&L->p0v, p, &tmp1);
		sp_mult(&tmp1, &L->p0, &tmp2);
		q = &tmp2;
		fst_moved = L->auts_first_moved[i];
		j = L->go[fst_moved];
			/* the automorphism moves j_ */
		/* j = L->p0v.a[j]; */
		k = q->a[j];
		if (L->info->f_very_verbose) {
			printf("i = %ld: fst_moved = %ld j = %ld k = %ld\n", 
				i, fst_moved, j, k);
			sp_print(p);
			fflush(stdout);
			}
		if (Aut->s_V_ii(k) == k) {
			Aut->s_V_i(k)->m_i(j);
			KM_p = Aut->s_KM_i(k);
			for (ii = 0; ii < n; ii++) {
				jj = q->a[ii];
				KM_p->m_ii(ii, jj + 1);
				}
			}
		else {
			if (L->info->f_very_verbose) {
				printf("skipped!\n");
				fflush(stdout);
				}
			}
		} /* next i */
	Aut->update_EM();
	Aut->group_order(ago);
	if (L->info->f_verbose) {
		fprintf(fp_txt, "ago = ");
		fflush(fp_txt);
		str[0] = 0;
		ago->sprint(str);
		fprintf(fp_txt, "%s\n", str);
		/* ago->println(); */
		/* Aut->print_group_order(); */
		fflush(fp_txt);
		}
	sp_free(&tmp1);
	sp_free(&tmp2);
	return OK;
}

#if TEXDOCU
static void print_E_colors(AUTLOG_LOCAL *L, INT level)
#endif
{
	ORDERED_SET *E;
	INT i, l, len, a, b;
	FILE *fp_txt;
	
	fp_txt = L->info->fp_txt;
	/* fprintf(fp_txt, "print_E_colors() level = %ld\n", level); */
	if (!L->pA->f_has_element_colors) {
		fprintf(fp_txt, "no colors !\n");
		return;
		}
	for (l = 0; l < level; l++) {
		E = L->E + l;
		len = E->size;
		fprintf(fp_txt, "E[%ld]={ ", l);
		for (i = 0; i < len; i++) {
			a = E->a[i];
			b = L->pA->element_colors[a];
			fprintf(fp_txt, "%ld ", b);
			}
		fprintf(fp_txt, "} %ld\n", len);
		}
	fflush(fp_txt);
}

#if TEXDOCU
static void print_element_colors(AUTLOG_LOCAL *L, INT level)
#endif
{
	FILE *fp_txt;
	INT i, ii, go, go_last = 0, i0, b;
	
	fp_txt = L->info->fp_txt;
	if (!L->pA->f_has_element_colors)
		return;
	for (ii = 0; ii <= level; ii++) {
		go = autlog_the_group_order(L->pA, ii);
		fprintf(fp_txt, "level %ld go = %ld elements: ", ii, go);
		for (i = go_last; i < go; i++) {
			i0 = L->pA->pv.a[i];
			b = L->pA->element_colors[i0];
			fprintf(fp_txt, "%ld ", b);
			}
		fprintf(fp_txt, "\n");
		go_last = go;
		}
	fflush(fp_txt);
}

#if TEXDOCU
static INT canonicize(AUTLOG_LOCAL *L, INT level)
#endif
{
	ORDERED_SET *E;
	INT i, ii, g, len, go, go_, go_last, ret;
	
	nb_backtrack++;
#ifdef PRINT_BACKTRACK_POINTS
	if ((nb_backtrack % BACKTRACK_MOD) ==  0) {
		fprintf(L->info->fp_txt, ".");
		nb_backtrack_points++;
		if (nb_backtrack_points >= 50) {
			fprintf(L->info->fp_txt, "\n");
			nb_backtrack_points = 0;
			}
		fflush(L->info->fp_txt);
		}
#endif
#ifdef PRINT_OCCASIONALLY
	if (level >= MIN_LEVEL_FOR_PRINT) {
		INT t, dt, dt_sec;
		
		t = os_ticks();
		dt = t - last_print_ticks;
		dt_sec = dt / os_ticks_per_second();
		if (dt_sec > PRINT_INTERVALL_SEC) {
			print_E_colors(L, level);
			print_element_colors(L, level);
			last_print_ticks = t;
			}
		}
#endif

	if (L->info->f_very_verbose) {
		printf("canonicize() level = %ld\n", level);
		/* printf("f_verbose = %ld f_very_verbose = %ld\n", 
			L->info->f_verbose, L->info->f_very_verbose); */
		fflush(stdout);
		if (L->info->f_very_verbose >= 2) {
			autlog_group_print(L->pA, &L->pA->p, &L->pA->pv);
			}
		}
	autlog_calc_grid(
		L->pA, L->pA->Grid, level, 
		L->pA->n /* go */, 
		L->info->f_very_verbose);
	E = L->E + level;
	autlog_choose_first_into_E(L->pA, L->pA->Grid, level, 
		E, L->info->f_very_verbose);
	len = E->size;
	for (i = 0; i < len; i++) {
		g = E->a[i];
		autlog_add_generator(L->pA, level, g, L->info->f_very_verbose);
		if (L->info->f_very_verbose) {
			if (L->info->f_very_verbose >= 2) {
				printf("canonicize() level = %ld added generator %ld\n", level, g);
				fflush(stdout);
				printf("canonicize() group order(%ld) = %ld\n", level + 1, 
					autlog_the_group_order(L->pA, level + 1));
				fflush(stdout);
				}
			print_labelling_g(L, level + 1);
			if (L->info->f_very_verbose >= 2) {
				autlog_group_print(L->pA, &L->pA->p, &L->pA->pv);
				}
			}
		if (i > 0 && level < L->first_moved)
			L->first_moved = level;
		if (L->is_first) {
			if (autlog_the_group_order(L->pA, level + 1) == L->pA->n) {
				L->is_first = FALSE;
				get_canonical_labelling(L, level + 1);
				}
			else {
				canonicize(L, level + 1);
				}
			}
		else { /* (!L->is_first) */
			go_last = autlog_the_group_order(L->pA, level);
			go = autlog_the_group_order(L->pA, level + 1);
			if (go == L->pA->n) {
				ret = canonicize_compare_tables(L, level + 1);
				if (ret == 0) {
					/* we have found an automorphism ! */
					if (L->info->f_very_verbose) {
						printf("found an automorphism ! in transversal %ld\n", 
							L->first_moved);
						fflush(stdout);
						print_labelling_g(L, level + 1);
						if (L->info->f_very_verbose >= 1) {
							autlog_group_print(L->pA, &L->pA->p, &L->pA->pv);
							}
						}
					get_aut_gen(L);
					if (L->info->f_very_verbose) {
						printf("no = %ld\n", L->nb_aut_gens);
						fflush(stdout);
						}
					if (level > L->first_moved)
						return 1;
					}
				if (ret == 1) {
					/* we found a better 'canonical' labelling ! */
					if (L->info->f_very_verbose) {
						printf("found a better labelling !\n");
						}
					get_canonical_labelling(L, level + 1);
					L->first_moved = AUTLOG_MAX_G;
					free_aut_gen(L);
					}
				}
			else {
				ret = canonicize_compare_tables(L, level + 1);
				if (ret >= 0) {
					if (canonicize(L, level + 1) == 1 && level > L->first_moved)
						return 1;
					}
				}
			}
		} /* next i */
	return 0;
}

#define REALLOC_OVERSIZE 4

#if TEXDOCU
static INT get_aut_gen(AUTLOG_LOCAL *L)
#endif
{
	INT new_dim, i, n;
	SPERM *new_auts, *p;
	INT *new_first;

	if (L->nb_aut_gens >= L->dim_aut_gens) {
		/* realloc: */
		new_dim = L->dim_aut_gens + REALLOC_OVERSIZE;
		new_auts = (SPERM *) my_malloc(new_dim * sizeof(SPERM), "aut.C: get_aut_gen()");
		if (new_auts == NIL)
			return error("get_aut_gen(): no memory for new_auts");
		new_first = (INT *) my_malloc(new_dim * sizeof(INT), "aut.C: get_aut_gen()");
		if (new_first == NIL)
			return error("get_aut_gen(): no memory for new_first");
		for (i = 0; i < L->nb_aut_gens; i++)
			new_first[i] = L->auts_first_moved[i];
		if (L->dim_aut_gens) {
			memcpy(new_auts, L->aut_gens, L->nb_aut_gens * sizeof(SPERM));
			my_free(L->aut_gens);
			my_free(L->auts_first_moved);
			}
		L->aut_gens = new_auts;
		L->auts_first_moved = new_first;
		L->dim_aut_gens = new_dim;
		}
	p = L->aut_gens + L->nb_aut_gens;
	n = L->pA->n;
	sp_nil(p);
	sp_int(p, n, "aut.C: get_aut_gen()");
	sp_mult(&L->p0, &L->pA->pv, p);


	if (L->info->f_very_verbose) {
		printf("pv=");
		sp_print(&L->pA->pv);
		printf("p0=");
		sp_print(&L->p0);
		sp_print(p);
#if 0
		PERMUTATION_OB perm;
		INT i, j, n;

		n = L->pA->n;
		perm.m_il(n);
		for (i = 0; i < n; i++) {
			j = p->a[i];
			perm.m_ii(i, j + 1);
			}
		perm.println();
#endif
		}
	L->auts_first_moved[L->nb_aut_gens] = L->first_moved;

	L->nb_aut_gens++;
	return OK;
}

#if TEXDOCU
static void free_aut_gen(AUTLOG_LOCAL *L)
#endif
{
	INT i;
	SPERM *p;

	if (L->dim_aut_gens) {
		for (i = 0; i < L->nb_aut_gens; i++) {
			p = L->aut_gens + i;
			sp_free(p);
			}
		my_free(L->aut_gens);
		my_free(L->auts_first_moved);
		}
	L->aut_gens = NIL;
	L->auts_first_moved = 0;
	L->nb_aut_gens = 0;
	L->dim_aut_gens = 0;
}

#if TEXDOCU
static INT get_canonical_labelling(AUTLOG_LOCAL *L, INT level)
#endif
{
	INT go, go_, ii;
	
	sp_mv(&L->pA->p, &L->p0);
	sp_mv(&L->pA->pv, &L->p0v);
	L->nb_gen = level;
	go_ = 0;
	for (ii = 0; ii <= L->nb_gen; ii++) {
		go = autlog_the_group_order(L->pA, ii);
		L->go[ii] = go;
		L->g[ii] = L->pA->pv.a[go_];
		go_ = go;
		}
	if (L->info->f_very_verbose) {
		printf("canonicize() canonic labelling:\n");
		printf("generators: ");
		print_labelling_g(L, level);
		printf("group orders: ");
		print_labelling_go(L, level);
		if (L->info->f_very_verbose >= 1) {
			autlog_group_print(L->pA, &L->pA->p, &L->pA->pv);
			}
		}
	return OK;
}

#if TEXDOCU
static void print_labelling_g(AUTLOG_LOCAL *L, INT level)
#endif
{
	INT go, go_, ii, g;
	
	printf("{ ");
	go_ = 0;
	for (ii = 0; ii <= level; ii++) {
		go = autlog_the_group_order(L->pA, ii);
		g = L->pA->pv.a[go_];
		printf("%ld ", g);
		go_ = go;
		}
	printf("}\n");
	fflush(stdout);
}

#if TEXDOCU
static void print_labelling_go(AUTLOG_LOCAL *L, INT level)
#endif
{
	INT go, ii;
	
	printf("go: { ");
	for (ii = 0; ii <= level; ii++) {
		go = autlog_the_group_order(L->pA, ii);
		printf("%ld ", go);
		}
	printf("}\n");
	fflush(stdout);
}

#if TEXDOCU
static INT canonicize_compare_tables(AUTLOG_LOCAL *L, INT level)
#endif
{
	INT i1, i0, i0_, j1, j0, j0_, k0, k0_, k1, k1_, i1_save;
	INT *theG, dim_n;
	INT ret;
	SPERM *p, *pv, *p_, *pv_;
	INT ii, go, go_;
	
	ret = 0;
	for (ii = 0; ii <= level; ii++) {
		go = autlog_the_group_order(L->pA, ii);
		go_ = L->go[ii];
		if (go > go_) {
			ret = 1;
			goto l_exit;
			}
		if (go < go_) {
			ret = -1;
			goto l_exit;
			}
		}
	go = autlog_the_group_order(L->pA, level);
	theG = L->pA->theG;
	dim_n = L->pA->dim_n;
	p = &L->pA->p;
	pv = &L->pA->pv;
	p_ = &L->p0;
	pv_ = &L->p0v; 
	for (i1 = 0; i1 < go; i1++) {
		for (j1 = 0; j1 <= i1; j1++) {
			i0 = pv->a[i1];
			i0_ = pv_->a[i1];
			j0 = pv->a[j1];
			j0_ = pv_->a[j1];
			k0 = theG[i0 * dim_n + j0];
			k0_ = theG[i0_ * dim_n + j0_];
			k1 = p->a[k0];
			k1_ = p_->a[k0_];
			if (k1 >= go) {
				return error("k1 >= go");
				}
			if (k1_ >= go) {
				return error("k1_ >= go");
				}
			if (k1_ > k1) {
				ret = 1;
				goto l_exit;
				}
			if (k1_ < k1) {
				ret = -1;
				goto l_exit;
				}
			} /* next j1 */
		i1_save = i1;
		j1 = i1;
		for (; i1 >= 0; i1--) {
			i0 = pv->a[i1];
			i0_ = pv_->a[i1];
			j0 = pv->a[j1];
			j0_ = pv_->a[j1];
			k0 = theG[i0 * dim_n + j0];
			k0_ = theG[i0_ * dim_n + j0_];
			k1 = p->a[k0];
			k1_ = p_->a[k0_];
			if (k1 >= go) {
				return error("k1 >= go");
				}
			if (k1_ >= go) {
				return error("k1_ >= go");
				}
			if (k1_ > k1) {
				ret = 1;
				goto l_exit;
				}
			if (k1_ < k1) {
				ret = -1;
				goto l_exit;
				}
			} /* next i1 */
		i1 = i1_save;
		}
#if 0
	if (level > 1) {
		ret = canonicize_compare_tables(L, level - 1);
		if (ret != 0)
			goto l_exit;
		}
	go_last = autlog_the_group_order(L->pA, level - 1);
	go = autlog_the_group_order(L->pA, level);
	ret = aut_table_compare(L->pA, 
		go_last, go, 
		go_last /* ib */, go /* ie */, 
		0 /* jb */, go /* je */, 
		&L->pA->p, &L->pA->pv, 
		&L->p0, &L->p0v, 
		L->info->f_verbose);
	if (ret == 0) {
		ret = aut_table_compare(L->pA, 
			go_last, go, 
			0 /* ib */, go_last /* ie */, 
			go_last /* jb */, go /* je */, 
			&L->pA->p, &L->pA->pv, 
			&L->p0, &L->p0v, 
			L->info->f_verbose);
		}
#endif

l_exit:
	if (L->info->f_very_verbose) {
		printf("compare tables: level = %ld ret = %ld\n", level, ret);
		}
	return ret;
}

#if TEXDOCU
INT autlog_test(AUTLOG_INFO *info)
#endif
{
	AUTLOG_LOCAL *L = NIL;
	
	L = (AUTLOG_LOCAL *) my_malloc(sizeof(AUTLOG_LOCAL), "autlog_test L");
	if (L == NIL)
		return error("autlog_test(): no memory for L");
	autlog_nil(L);
	if (info->f_verbose) {
		/* print_size(); */
		printf("("); fflush(stdout);
		}
	if (autlog_int(L, info, TRUE /* f_B */)) {
		if (info->f_verbose) {
			printf("|"); fflush(stdout);
			}
		autlog_do_level(L, 0);
		}
	else {
		info->nb_isomorphisms = 0;
		}
	
	
	autlog_ext(L);
	my_free(L);
	if (info->f_verbose) {
		printf(")\n"); fflush(stdout);
		}
	return TRUE;
}

#if TEXDOCU
static INT autlog_do_level(AUTLOG_LOCAL *L, INT i)
#else
/* Rueckgabe FALSE, falls wegen f\_break\_after\_first
 * abgebrochen wurde. */
#endif
{
	INT j;
	
	if (i == 0) {
		if (L->E == NIL)
			return error("autlog_do_level(): L->E == NIL");
		}
	
	if (L->info->f_very_verbose) {
		printf("*** AUTLOG: level %ld ***\n", i);
		}

	if (i == L->pA->nb_gen) {
		if (autlog_the_group_order(L->pA, i) != L->pA->n)
			return error("autlog_do_level(): "
			"autlog_the_group_order(L->pA, i) "
			"!= L->pA->n");
		L->info->nb_isomorphisms++;
		if (L->info->f_verbose) {
			printf("*** found automorphism nr. %ld !\n", 
				L->info->nb_isomorphisms);
			}
		autlog_print_mapping(L, i);
		
		if (L->info->mode == AUT_MODE_BREAK_AFTER_FST)
			return FALSE; /* terminate program */
		else {
			if (L->info->mode == AUT_MODE_ONLY_COSET_REPS) {
				L->f_going_back = TRUE;
				if (L->info->nb_isomorphisms == 1)
					L->first_moved = i;
				}
			return TRUE; /* proceed further */
			}
		}
		
	/* aut_calc_grid(a, TRUE, i); */
	
	autlog_calc_grid(
		L->pB, L->pB->Grid, i, 
		L->pA->n /* go */, 
		L->info->f_very_verbose);
	/* aut_calc_grid(a, FALSE, i); */
	
	autlog_E(L->pA, L->pB, 
		L->pA->Grid + i, 
		L->pB->Grid, i, L->E + i, 
		L->info->f_very_verbose);
	/* aut_check_iso(a, i, &a->E[i]); */

	if (L->info->f_very_verbose && 
		L->info->f_verbose) {
		autlog_print_E(L->E + i, i);
		}
	
	autlog_reduce(L, L->pA, L->pB, 
		L->pA->Grid + i, L->pB->Grid, 
		L->E + i, i, 
		L->info->f_very_verbose);
	/* aut_reduce(a, i); */

	if (L->info->f_very_verbose) {
		autlog_print_E(L->E + i, i);
		}

	if (L->E[i].size == 0L) {
		return TRUE;
		}
	
	while (L->E[i].size > 0) {
		
		autlog_add_generator(L->pB, i, L->E[i].a[0], 
			L->info->f_very_verbose);
		
		if (!autlog_do_level(L, i + 1)) {
			return FALSE;
			}
		
		if (L->info->mode == AUT_MODE_ONLY_COSET_REPS) {
			if (L->f_going_back && L->first_moved < i)
				return TRUE;
			L->f_going_back = FALSE;
			if (L->first_moved > i)
				L->first_moved = i;
			}
			
		/* delete first of E[i]: */
		for (j = 1; j < L->E[i].size; j++) {
			L->E[i].a[j - 1] = L->E[i].a[j];
			}
		L->E[i].size--;
		if (L->info->f_very_verbose) {
			autlog_print_E(L->E + i, i);
			}
		
		}

	return TRUE;
}

#if TEXDOCU
static void autlog_print_mapping(
	AUTLOG_LOCAL *L, INT i)
#endif
{
	INT l, g, g1, g2;
	BYTE s[256];
	PERMUTATION_OB p;
	INT j;
	
	sp_mult(&L->pA->p, &L->pB->pv, &L->pA->tmp1);
	if (L->info->f_verbose) {
		s[0] = 0;
		sp_sprint(&L->pA->tmp1, s);
		printf("%s\n", s);
		for (l = 0; l < L->pA->nb_gen; l++) {
			g = L->pA->g[l];
			/* g2 = iso->Bg[l]; */
			g2 = L->pA->tmp1.a[g];
			
			/* g1 = iso->Ap.a[g];
			g2 = iso->Bpv.a[g1]; */
			
			printf("%ld -> %ld", g, g2);
			if (l < L->pA->nb_gen - 1)
				printf(", ");
			}
		printf("\n");
		}
	if (L->info->add_aut) {
		if ((*L->info->add_aut)(
			L->info, L->pA->nb_gen, 
			L->pA->g, &L->pA->tmp1) != OK) {
			Srff("autlog_print_mapping", "add_aut");
			}
		}
	if (L->info->add_aut_perm) {
		p.m_il(L->pA->n);
		for (l = 0; l < L->pA->n; l++) {
			j = L->pA->tmp1.a[l] + 1;
			p.m_ii(l, j);
			}
		if ((*L->info->add_aut_perm)(
			L->info, L->pA->nb_gen, 
			L->pA->g, &p) != OK) {
			Srff("autlog_print_mapping", "add_aut_perm");
			}
		}
}

//#include "aut_util.C"
/* aut_util.C 
 * included from aut.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#if TEXDOCU
static void print_size()
#endif
{
	printf("AUTLOG_MAX_N = %ld\n", (INT) AUTLOG_MAX_N);
	printf("AUTLOG_MAX_TYPE = %ld\n", (INT) AUTLOG_MAX_TYPE);
	printf("AUTLOG_MAX_GRID = %ld\n", (INT) AUTLOG_MAX_GRID);
	printf("sizeof ORDERED_SET = %ld\n", (INT) sizeof(ORDERED_SET));
	printf("sizeof AUTLOG_THE_GROUP = %ld\n", (INT) sizeof(AUTLOG_THE_GROUP));
	printf("sizeof AUTLOG_GRID = %ld\n", (INT) sizeof(AUTLOG_GRID));
	printf("sizeof AUTLOG_FACTORGROUP = %ld\n", (INT) sizeof(AUTLOG_FACTORGROUP));
	printf("sizeof AUTLOG_COMMUTATOR_SUBGROUP = %ld\n", (INT) sizeof(AUTLOG_COMMUTATOR_SUBGROUP));
	printf("sizeof AUTLOG_LOCAL = %ld\n", (INT) sizeof(AUTLOG_LOCAL));
	fflush(stdout);
}

#if TEXDOCU
INT autlog_GmZ_Ad_orders(AUTLOG_INFO *info, 
	INT *GmZ_order, INT *Gd_order, FILE *fp_txt)
#else
/* the simplest family invariants. */
#endif
{
	AUTLOG_LOCAL *L = NIL;
	INT n1, n2;
	
	info->f_autologisms = TRUE;
	L = (AUTLOG_LOCAL *) my_malloc(sizeof(AUTLOG_LOCAL), "autlog_GmZ_Ad_orders");
	if (L == NIL)
		return error("autlog_GmZ_Ad_orders(): no memory for L");
	autlog_nil(L);
	if (info->f_verbose) {
		fprintf(fp_txt, "|GmZ|");
		fflush(fp_txt);
		}
	autlog_int(L, info, FALSE /* f_B */);
	if (info->f_verbose) {
		fprintf(fp_txt, " = ");
		fflush(fp_txt);
		}
	
	n1 = L->AmZ.n;
	n2 = autlog_the_group_order(&L->Ad, L->Ad.nb_gen);
	*GmZ_order = n1;
	*Gd_order = n2;
	if (info->f_verbose) {
		fprintf(fp_txt, "%ld |Gd| = %ld\n", n1, n2);
		fflush(fp_txt);
		}
	
	autlog_ext(L);
	my_free(L);
	return TRUE;
}

#if TEXDOCU
INT al_test()
#endif
{
	FG_OB Q8;
	ZE_OP ze0, ze1, ze2;
	INT primes[] = {2, 2, 2};
	
	Q8.init(3, primes);
	ze0 = Q8.s_ze_i(0);
	ze1 = Q8.s_ze_i(1);
	ze2 = Q8.s_ze_i(2);

	ze1->s_P()->m_i(1);
	ze1->s_A_i(0)->m_i(1);
	ze1->s_Av_i(0)->m_i(1);
	
	ze2->s_P()->m_i(1);
	ze2->s_A_i(0)->m_i(1);
	ze2->s_Av_i(0)->m_i(1);
	ze2->s_A_i(1)->m_i(3);
	ze2->s_Av_i(1)->m_i(3);
	
	printf("vor Q8.theG()\n");
	fflush(stdout);
	Q8.theG(TRUE /* f_verbose */);
	printf("nach Q8.theG()\n");
	fflush(stdout);
	
	al_test2(&Q8, TRUE);
	return OK;
}

#if TEXDOCU
static INT al_test2(FG_OP G, INT f_verbose)
#endif
{
	ZE_OP ze;
	INT *theG = NIL, size;
	AUTLOG_INFO info;
	INT l, len, g, g_i, i, j, k, n, dim_n, ago;
	INT erg = OK;
	
	n = G->s_n_i();
	dim_n = n + 1;
	((INTEGER_OP) G->s_ago())->m_i(0);
	size = n * dim_n * sizeof(INT);
	theG = (INT *) my_malloc(size, "al_test2");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			k = G->s_theG_iji(i, j);
			theG[i * dim_n + j] = k;
			if (k == 0)
				theG[i * dim_n + n] = j;
			if (f_verbose)
				printf("%ld ", theG[i * dim_n + j]);
			}
		if (f_verbose)
			printf("\n");
		}
	info.An = n;
	info.Bn = n;
	info.AtheG = theG;
	info.BtheG = theG;
	info.Adim_n = dim_n;
	info.Bdim_n = dim_n;
	info.nb_isomorphisms = 0;
	info.mode = AUT_MODE_FULL
		/* or AUT_MODE_ONLY_COSET_REPS */;
	info.f_autologisms = TRUE;
	info.f_verbose = TRUE;
	info.f_very_verbose = TRUE;
	len = G->s_nb_ze_i();
	if (len >= AUTLOG_MAX_G) {
		Srfs("al_test2", "len >= MAX_G");
		return ERROR;
		}
	for (l = 0; l < len; l++) {
		ze = G->s_ze_i(l);
		g_i = ze->s_n0_i();
		if (f_verbose)
			printf("%ld ", g_i);
		info.A_g[l] = g_i;
		info.B_g[l] = g_i;
		}
	if (f_verbose)
		printf("\n");
	info.A_nb_gen = len;
	info.B_nb_gen = len;

	info.add_aut = NIL;
	info.add_aut_perm = NIL;
#if 0
	M.m_ilih(n, len);
	info.data = (void *) &M;
#endif
	printf("vor autlog_test()\n");
	fflush(stdout);
	if (!autlog_test(&info))
		return error("al_test2(): autlog_test");
	printf("nach autlog_test()\n");
	fflush(stdout);
	return OK;
}

#if TEXDOCU
static void autlog_the_group_nil(AUTLOG_THE_GROUP *G)
#endif
{
	G->theG = NIL;
	G->Grid = NIL;
	sp_nil(&G->p);
	sp_nil(&G->pv);
	sp_nil(&G->q);
	sp_nil(&G->qv);
	sp_nil(&G->tmp1);
	sp_nil(&G->tmp2);
	G->nb_grid = 0;
	G->f_has_element_colors = FALSE;
	G->element_colors = NIL;
	G->fp_txt = NIL;
}

#if TEXDOCU
static void autlog_the_group_int(AUTLOG_THE_GROUP *G)
#endif
{
	sp_int(&G->p, G->n, "autlog_the_group_int");
	sp_int(&G->pv, G->n, "autlog_the_group_int");
	sp_int(&G->q, G->n, "autlog_the_group_int");
	sp_int(&G->qv, G->n, "autlog_the_group_int");
	sp_int(&G->tmp1, G->n, "autlog_the_group_int");
	sp_int(&G->tmp2, G->n, "autlog_the_group_int");
	sp_id(&G->p);
	sp_id(&G->pv);
	sp_id(&G->q);
	sp_id(&G->qv);
	sp_id(&G->tmp1);
	sp_id(&G->tmp2);
	G->f_has_element_colors = FALSE;
}

#if TEXDOCU
static INT autlog_the_group_init_colors(AUTLOG_THE_GROUP *G, VECTOR_OP col)
#endif
{
	INT n, i, a;
	
	G->f_has_element_colors = TRUE;
	n = G->n;
	if (n != col->s_li())
		return error("autlog_the_group_init_colors(): n != col->s_li()");
	G->element_colors = (INT *) my_malloc(n * sizeof(INT), "autlog_the_group_init_colors");
	for (i = 0; i < n; i++) {
		a = col->s_ii(i);
		G->element_colors[i] = a;
		}
	return OK;
}

#if TEXDOCU
static INT autlog_the_group_open_grid(AUTLOG_THE_GROUP *G, INT nb_grid, INT n)
#endif
{
	INT i;
	
	if (G->nb_grid)
		return error("autlog_the_group_open_grid"
		"(): G->nb_grid");
	G->Grid = (AUTLOG_GRID *) my_malloc(nb_grid * sizeof(AUTLOG_GRID), 
		"autlog_the_group_open_grid");
	if (G->Grid == NIL)
		return error("autlog_the_group_open_grid(): no memory for G->Grid");
	G->nb_grid = nb_grid;
	for (i = 0; i < G->nb_grid; i++) {
		algrid_nil(G->Grid + i);
		algrid_int(G->Grid + i, n);
		}
	return OK;
}

#if TEXDOCU
static void autlog_the_group_ext(AUTLOG_THE_GROUP *G)
#endif
{
	INT i;
	
	my_ptr_free( (void **) &G->theG);
	if (G->nb_grid) {
		for (i = 0; i < G->nb_grid; i++) {
			algrid_ext(G->Grid + i);
			}
		my_ptr_free( (void **) &G->Grid);
		G->nb_grid = 0;
		}
	sp_free(&G->p);
	sp_free(&G->pv);
	sp_free(&G->q);
	sp_free(&G->qv);
	sp_free(&G->tmp1);
	sp_free(&G->tmp2);
	if (G->f_has_element_colors) {
		my_free(G->element_colors);
		G->element_colors = NIL;
		G->f_has_element_colors = FALSE;
		}
}

#if TEXDOCU
static INT autlog_the_group_init_generators(
	AUTLOG_THE_GROUP *G, 
	INT nb_gen, INT *gen)
#endif
{
	INT i;
	
	for (i = 0; i < nb_gen; i++)
		G->g[i] = gen[i];
	G->nb_gen = nb_gen;
	return OK;
}

#if TEXDOCU
static INT autlog_the_group_print_generators(
	AUTLOG_THE_GROUP *G)
#endif
{
	INT i;
	
	for (i = 0; i < G->nb_gen; i++)
		printf("%ld ", G->g[i]);
	printf("\n");
	return OK;
}

#if TEXDOCU
static INT autlog_the_group_order(
	AUTLOG_THE_GROUP *G, INT k)
#else
/* k = 0: 1 fuer Einsgruppe, sonst go[k - 1] */
#endif
{
	if (k < 0)
		return error("autlog_the_group_order k < 0");
	if (k == 0)
		return 1;
	return G->go[k - 1];
}

#if TEXDOCU
static void alfg_nil(AUTLOG_FACTORGROUP *p)
#endif
{
	p->H = NIL;
	p->H_len = 0;
	p->H_rep_idx = NIL;
	p->H_rep = NIL;
	p->nb_gen = 0;
	p->g = NIL;
	p->g_idx = NIL;
}

#if TEXDOCU
static void alfg_ext(AUTLOG_FACTORGROUP *p)
#endif
{
	my_ptr_free( (void **) &p->H);
	my_ptr_free( (void **) &p->H_rep_idx);
	my_ptr_free( (void **) &p->H_rep);
	my_ptr_free( (void **) &p->g);
	my_ptr_free( (void **) &p->g_idx);
	p->H_len = 0;
	p->nb_gen = 0;
}

#if TEXDOCU
static INT alfg_print_generators(
	AUTLOG_FACTORGROUP *G_Z)
#endif
{
	INT i;
	
	for (i = 0; i < G_Z->nb_gen; i++)
		printf("%ld ", G_Z->g[i]);
	printf("\n");
	return OK;
}

#if TEXDOCU
static void alcs_nil(AUTLOG_COMMUTATOR_SUBGROUP *p)
#endif
{
	p->C_idx = NIL;
	p->nb_C = 0;
	p->C = NIL;
	p->nb_gen = 0;
	p->g_len = NIL;
	p->g = NIL;
	p->g_idx = NIL;
}

#if TEXDOCU
static void alcs_ext(AUTLOG_COMMUTATOR_SUBGROUP *p)
#endif
{
	my_ptr_free( (void **) &p->C_idx);
	my_ptr_free( (void **) &p->C);
	my_ptr_free( (void **) &p->g_len);
	my_ptr_free( (void **) &p->g);
	my_ptr_free( (void **) &p->g_idx);
	p->nb_C = 0;
	p->nb_gen = 0;
}

#if TEXDOCU
static void algrid_nil(AUTLOG_GRID *p)
#endif
{
	sp_nil(&p->r);
	sp_nil(&p->rv);
}

#if TEXDOCU
static void algrid_int(AUTLOG_GRID *p, INT n)
#endif
{
	p->m = n;
	sp_int(&p->r, p->m, "algrid_int");
	sp_int(&p->rv, p->m, "algrid_int");
	sp_id(&p->r);
	sp_id(&p->rv);
}

#if TEXDOCU
static void algrid_ext(AUTLOG_GRID *p)
#endif
{
	sp_free(&p->r);
	sp_free(&p->rv);
}

#if TEXDOCU
static void autlog_nil(AUTLOG_LOCAL *a)
#endif
{
	a->info = NIL;
	a->E = NIL;
	autlog_the_group_nil(&a->A);
	autlog_the_group_nil(&a->B);
	autlog_the_group_nil(&a->AmZ);
	autlog_the_group_nil(&a->BmZ);
	alfg_nil(&a->A_Z);
	alfg_nil(&a->B_Z);
	alcs_nil(&a->AcA);
	alcs_nil(&a->BcB);
	autlog_the_group_nil(&a->Ad);
	autlog_the_group_nil(&a->Bd);
}

#if TEXDOCU
static void autlog_group_print(AUTLOG_THE_GROUP *G, 
	SPERM *p, SPERM *pv)
#endif
{
	INT i, j, i0, j0, k0, k;
	
	for (i = 0; i < G->n; i++) {
		i0 = pv->a[i];
		for (j = 0; j < G->n; j++) {
			j0 = pv->a[j];
			k0 = G->theG[i0 * G->dim_n + j0];
			/* k = p->a[k0]; */
			k = k0;
			printf("%ld ", k);
			}
		printf("\n");
		}
	printf("\n");
	printf("p  = [");
	for (i = 0; i < p->l; i++) {
		j = p->a[i];
		printf("%ld ", j);
		}
	printf("]\n");
	printf("pv = [");
	for (i = 0; i < pv->l; i++) {
		j = pv->a[i];
		printf("%ld ", j);
		}
	printf("]\n");
	fflush(stdout);
}

//#include "aut_init.C"
/* aut_init.C 
 * included from aut.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#if TEXDOCU
static INT autlog_int(AUTLOG_LOCAL *a, AUTLOG_INFO *info, INT f_B)
#endif
{
	a->info = info;
	a->A.n = a->info->An;
	a->B.n = a->info->Bn;
	autlog_the_group_int(&a->A);
	if (info->f_A_has_colors)
		autlog_the_group_init_colors(&a->A, info->A_colors);
	if (f_B) {
		autlog_the_group_int(&a->B);
		if (info->f_B_has_colors)
			autlog_the_group_init_colors(&a->B, info->B_colors);
		}
	
	a->A.theG = info->AtheG;
	a->B.theG = info->BtheG;
	a->A.dim_n = info->Adim_n;
	a->B.dim_n = info->Bdim_n;
	a->f_going_back = FALSE;
	
	if (!info->f_autologisms)
		autlog_int_isomorphism(a, 
		info->f_very_verbose /* f_verbose */);
	else {
		if (!autlog_int_autologism(a, f_B, 
			info->f_very_verbose /* f_verbose */, 
			info->f_very_verbose /* f_very_verbose */))
			/* kein Autologismus moeglich: */
			return FALSE;
		}
	/* now a->pA and a->pB set. */

	a->E = (ORDERED_SET *) my_malloc(a->pA->nb_gen * sizeof(ORDERED_SET), "autlog_int E");
	if (a->E == NIL)
		return error("autlog_int(): no memory for a->E");

	return TRUE;
}

#if TEXDOCU
static INT autlog_ext(AUTLOG_LOCAL *a)
#endif
{
	/* Kopien auf info, nicht freigeben: */
	a->A.theG = NIL;
	a->B.theG = NIL;
	a->Ad.theG = NIL;
	a->Bd.theG = NIL;
	
	my_ptr_free( (void **) &a->E);
	
	autlog_the_group_ext(&a->A);
	autlog_the_group_ext(&a->B);
	alfg_ext(&a->A_Z);
	alfg_ext(&a->B_Z);
	autlog_the_group_ext(&a->AmZ);
	autlog_the_group_ext(&a->BmZ);
	
	alcs_ext(&a->AcA);
	alcs_ext(&a->BcB);

	autlog_the_group_ext(&a->Ad);
	autlog_the_group_ext(&a->Bd);
	return TRUE;
}

#if TEXDOCU
static INT autlog_int_isomorphism(
	AUTLOG_LOCAL *L, INT f_verbose)
#endif
{
	if (L->A.n != L->B.n)
		return error("autlog_int_isomorphism(): L->A.n != L->B.n - no iso !");
	{
		INT *g_idx = NIL;
	
		autlog_group_gen(&L->A, 
			L->info->A_nb_gen, 
			L->info->A_g, &g_idx, f_verbose);
		if (f_verbose) {
			printf("autlog_int_isomorphism(): nach autlog_group_gen(A)\n");
			fflush(stdout);
			}
		my_ptr_free( (void **) &g_idx);
		if (f_verbose) {
			printf("autlog_int_isomorphism(): nach my_ptr_free()\n");
			fflush(stdout);
			}

	}
	
	autlog_the_group_open_grid(&L->A, L->A.nb_gen /* nb_grid */, L->A.n /* n */);
	if (f_verbose) {
		printf("autlog_int_isomorphism(): nach autlog_the_group_open_grid(A)\n");
		fflush(stdout);
		}
	
	autlog_calc_all_grids(&L->A, L->A.n /* go */, f_verbose);
	if (f_verbose) {
		printf("autlog_int_isomorphism(): nach autlog_calc_all_grids(A)\n");
		fflush(stdout);
		}
	
	autlog_the_group_open_grid(&L->B, 1 /* nb_grid */, L->B.n /* n */);
	if (f_verbose) {
		printf("autlog_int_isomorphism(): nach autlog_the_group_open_grid(B)\n");
		fflush(stdout);
		}

	L->pA = &L->A;
	L->pB = &L->B;
	
	return OK;
}

#if TEXDOCU
static INT autlog_int_autologism(
	AUTLOG_LOCAL *L, INT f_B, 
	INT f_verbose, INT f_very_verbose)
#else
/* Rueckgabe FALSE: kein Autologismus moeglich. */
#endif
{
	INT i;
	
	autlog_the_group_init_generators(&L->A, L->info->A_nb_gen, L->info->A_g);
	if (f_B) {
		autlog_the_group_init_generators(&L->B, L->info->B_nb_gen, L->info->B_g);
		}
	
	if (f_verbose) {
		printf("autlog_int(): vor autlog_center()\n");
		fflush(stdout);
		}
	
	autlog_do_center(&L->A, &L->A_Z, f_very_verbose);
	if (f_B) {
		autlog_do_center(&L->B, &L->B_Z, f_very_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_center()\n");
		fflush(stdout);
		}
	
	autlog_do_do_factorgroup(
		&L->A, &L->A_Z, &L->AmZ, 
		f_very_verbose);
	if (f_B) {
		autlog_do_do_factorgroup(
			&L->B, &L->B_Z, &L->BmZ, 
			f_very_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_do_factorgroup()\n");
		fflush(stdout);
		}
	if (f_B) {
		if (L->AmZ.n != L->BmZ.n) {
			if (f_verbose) {
				printf("autlog_int(): L->AmZ.n != L->BmZ.n - no iso !\n");
				fflush(stdout);
				}
			return FALSE;
			}
		}
	
	if (f_verbose) {
		printf("generators for A mod Z:\n");
		alfg_print_generators(&L->A_Z);
		if (f_B) {
			printf("generators for B mod Z:\n");
			alfg_print_generators(&L->B_Z);
			}
		}
	
	autlog_do_group_gen_from_factorgroup(
		&L->AmZ, &L->A_Z, f_very_verbose);
	if (f_B) {
		autlog_do_group_gen_from_factorgroup(
			&L->BmZ, &L->B_Z, f_very_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_do_group_gen_from_factorgroup()\n");
		fflush(stdout);
		}
	
	if (f_verbose) {
		printf("reduced generators for A mod Z:\n");
		autlog_the_group_print_generators(&L->AmZ);
		if (f_B) {
			printf("reduced generators for B mod Z:\n");
			autlog_the_group_print_generators(&L->BmZ);
			}
		}

	autlog_the_group_open_grid(&L->AmZ, 
		L->AmZ.nb_gen /* (kann Null sein - 
			Grupppe abelsch) nb_grid */, 
		L->AmZ.n /* n */);
	if (f_B) {
		autlog_the_group_open_grid(&L->BmZ, 
			1 /* L->BmZ.nb_gen */ /* nb_grid */, 
			L->BmZ.n /* n */);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_the_group_open_grid(GmZ)\n");
		fflush(stdout);
		}
	
	autlog_calc_all_grids(&L->AmZ, L->AmZ.n /* go */, f_verbose);
	if (f_B) {
		autlog_calc_n_grids(&L->BmZ, L->BmZ.n /* go */, 1, f_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_calc_all(n)_grids(GmZ)\n");
		fflush(stdout);
		}
	if (f_B) {
		if (L->AmZ.nb_gen) {
			if (!autlog_check_iso(
				L->AmZ.Grid + 0, L->BmZ.Grid + 0)) {
				if (f_verbose) {
					printf("autlog_int(): no iso !\n");
					fflush(stdout);
					}
				return FALSE;
				}
			}
		}
	
	autlog_do_derivedgroup(&L->A, &L->A_Z, &L->AmZ, &L->AcA, &L->Ad, f_very_verbose);
	if (f_B) {
		autlog_do_derivedgroup(&L->B, &L->B_Z, &L->BmZ, &L->BcB, &L->Bd, f_very_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_do_derivedgroup()\n");
		fflush(stdout);
		}
	
	autlog_do_derivedgroup_generators(&L->AmZ, &L->AcA, f_very_verbose);
	if (f_B) {
		autlog_do_derivedgroup_generators(&L->BmZ, &L->BcB, f_very_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_do_derivedgroup_generators()\n");
		fflush(stdout);
		}
	
	autlog_group_gen(&L->Ad, 
		L->AcA.g_len[MAXIMUM(L->AcA.nb_gen - 1, 0)], 
		L->AcA.g, &L->AcA.g_idx, f_verbose);
	if (f_B) {
		autlog_group_gen(&L->Bd, 
			L->BcB.g_len[MAXIMUM(L->BcB.nb_gen - 1, 0)], 
			L->BcB.g, &L->BcB.g_idx, f_verbose);
		}
	if (f_verbose) {
		printf("autlog_int(): nach autlog_group_gen(Gd)\n");
		fflush(stdout);
		}
	
	{
		INT Ad_n = 0, Bd_n = 0;
		
		Ad_n = autlog_the_group_order(&L->Ad, L->Ad.nb_gen);
		if (f_B)
			Bd_n = autlog_the_group_order(&L->Bd, L->Bd.nb_gen);
		if (f_verbose) {
			printf("Ad_n = %ld Bd_n = %ld\n", 
				Ad_n, Bd_n);
			fflush(stdout);
			}
		if (f_B) {
			if (Ad_n != Bd_n) {
				if (f_verbose) {
					printf("no iso !\n");
					fflush(stdout);
					}
				return FALSE;
				}
			}

		autlog_the_group_open_grid(&L->Ad, 
			L->Ad.nb_gen /* kann 0 sein; nb_grid */, 
			L->Ad.n /* n */);
		if (f_B) {
			autlog_the_group_open_grid(&L->Bd, 
				1 /* L->Bd.nb_gen */ /* nb_grid */, 
				L->Bd.n /* n */);
			}
		if (f_verbose) {
			printf("autlog_int(): nach autlog_the_group_open_grid(Gd)\n");
			fflush(stdout);
			}
		
		autlog_calc_all_grids(&L->Ad, 
			Ad_n /* nicht: L->Ad.n */ /* go */, f_verbose);
		if (f_B) {
			autlog_calc_n_grids(&L->Bd, 
				Bd_n /* nicht: L->Bd.n */ /* go */, 1, f_verbose);
			}
		if (f_verbose) {
			printf("autlog_int(): nach autlog_calc_all(n)_grids(Gd)\n");
			fflush(stdout);
			}
		if (f_B) {
			if (L->Ad.nb_gen) {
				if (!autlog_check_iso(
					L->Ad.Grid + 0, L->Bd.Grid + 0)) {
					if (f_verbose) {
						printf("autlog_int(): no iso !\n");
						fflush(stdout);
						}
					return FALSE;
					}
				}
			}
	}
	
	L->pA = &L->AmZ;
	L->pB = &L->BmZ;
	
	return TRUE;
}

#if TEXDOCU
static INT autlog_do_group_gen_from_factorgroup(
	AUTLOG_THE_GROUP *GmZ, 
	AUTLOG_FACTORGROUP *G_Z, 
	INT f_verbose)
#endif
{
	autlog_group_gen(GmZ, 
		G_Z->nb_gen, G_Z->g, &G_Z->g_idx, f_verbose);
	return OK;
}

#if TEXDOCU
static INT autlog_do_center(
	AUTLOG_THE_GROUP *G, 
	AUTLOG_FACTORGROUP *G_Z, 
	INT f_verbose)
#endif
{
	autlog_center(G, &G_Z->H, &G_Z->H_len, f_verbose);
	return OK;
}

#if TEXDOCU
static INT autlog_do_do_factorgroup(
	AUTLOG_THE_GROUP *G, 
	AUTLOG_FACTORGROUP *G_Z, 
	AUTLOG_THE_GROUP *GmZ, 
	INT f_verbose)
#endif
{
	autlog_do_factorgroup(G, G_Z->H, G_Z->H_len, 
		&G_Z->H_rep_idx, &G_Z->H_rep, 
		&GmZ->theG, 
		&G_Z->g /* vorlaeufige Generatoren fuer GmZ */, 
		&GmZ->n, &GmZ->dim_n, 
		f_verbose);
	G_Z->nb_gen = G->nb_gen;
	/* SPERMs initialisieren: */
	autlog_the_group_int(GmZ);
	return OK;
}

#if TEXDOCU
static INT autlog_do_derivedgroup(
	AUTLOG_THE_GROUP *G, 
	AUTLOG_FACTORGROUP *G_Z, 
	AUTLOG_THE_GROUP *GmZ, 
	AUTLOG_COMMUTATOR_SUBGROUP *GcG, 
	AUTLOG_THE_GROUP *Gd, 
	INT f_verbose)
#endif
{
	INT erg = OK;
	
	erg += autlog_derivedgroup(G, 
		G_Z->H_rep_idx, G_Z->H_rep, 
		GmZ->theG, GmZ->n, GmZ->dim_n,
		&GcG->C_idx, &GcG->C, &GcG->nb_C, 
		f_verbose);
	if (erg) {
		printf("error in autlog_derivedgroup()\n");
		fflush(stdout);
		}
	GcG->dim_n = GmZ->n;
	Gd->n = G->n;
	Gd->theG = G->theG;
	Gd->dim_n = G->dim_n;
	autlog_the_group_int(Gd);
		/* SPERMs initialisieren */
	return erg;
}

#if TEXDOCU
static INT autlog_do_derivedgroup_generators(
	AUTLOG_THE_GROUP *GmZ, 
	AUTLOG_COMMUTATOR_SUBGROUP *GcG, 
	INT f_verbose)
#endif
{
	INT erg = OK;
	
	erg += autlog_derivedgroup_generators(GmZ, 
		GcG->C_idx, GcG->C, GcG->nb_C, 
		&GcG->g_len /* Anzahl der Generatoren 
			in g aus der i - ten Untergruppe von GmZ */, 
		&GcG->g /* vorlaeufige Generatoren fuer Gd, 
			erster Eintrag ist 0 (= id), 
			Abschlusseintrag -1. */, 
		f_verbose);
	if (erg) {
		printf("error in autlog_derivedgroup()\n");
		fflush(stdout);
		}
	GcG->nb_gen = GmZ->nb_gen;
		/* = Laenge von g_len[] */
	return erg;
}

#if TEXDOCU
static INT autlog_center(AUTLOG_THE_GROUP *G, INT **Z, INT *nb_Z, INT f_verbose)
#endif
{
	INT *Z1 = NIL;
	INT i, j, k, nb_Z1 = 0;
	
	Z1 = (INT *) my_malloc(G->n * sizeof(INT), "autlog_center");
	if (Z1 == NIL)
		return error("autlog_center(): no memory for Z1");
	for (i = 0; i < G->n; i++) {
		for (j = 0; j < G->nb_gen; j++) {
			k = autlog_conjugate_G(G, i, G->g[j]);
			if (k != i)
				break;
			}
		if (j == G->nb_gen) {
			Z1[nb_Z1] = i;
			nb_Z1++;
			}
		}
	if (f_verbose) {
		printf("autlog_center(): center of order %ld.\n", nb_Z1);
		}
	*Z = Z1;
	*nb_Z = nb_Z1;
	return OK;
}

#if TEXDOCU
static INT autlog_do_factorgroup(
	AUTLOG_THE_GROUP *G, INT *H, INT H_len, 
	INT **H_rep_idx, INT **H_rep, 
	INT **GmH, INT **GmH_gen, 
	INT *GmH_n, INT *GmH_dim, 
	INT f_verbose)
#else
//PRE
/* Input:
 * G->n group order, 
 * G->g[], G->nb\_gen generators.
 * H: the elements (their indices into G) 
 *    of the subgroup
 * H\_len: size of subgroup.
 * Output:
 * H\_rep\_idx: the index of the coset 
 *            containing the element of G.
 * H\_rep: index (into G) of coset 
 *        representatives (GmH\_n many).
 * GmH: the factor group table, 
 *      last column are inverse elements.
 * GmH\_gen: G->g[] as generators of the 
 *          factor group (0, 1, .. G->nb\_gen);
 *    GmH\_gen[i] = H\_rep\_idx[G->g[i]].
 * GmH\_n: order of the factor group 
 *        (= G->n / H\_len).
 * GmH\_dim: length of factor group 
 *     table (width), (= GmH\_n + 1) */
///PRE
#endif
{
	INT *H_rep_idx1 = NIL;
	INT *H_rep1 = NIL;
	INT *GmH1 = NIL;
	INT *GmH_gen1 = NIL;
	INT GmH_dim1 = 0;
	INT GmH_n1 = 0;
	INT i, j, k, i1, j1, k1, idx_inv;
	
	H_rep_idx1 = (INT *) my_malloc(G->n * sizeof(INT), "autlog_do_factorgroup");
	H_rep1 = (INT *) my_malloc(G->n * sizeof(INT), "autlog_do_factorgroup");
	if (H_rep_idx1 == NIL || H_rep1 == NIL)
		return error("autlog_do_factorgroup(): "
		"no memory for H_rep_idx1 / H_rep1");
	for (i = 0; i < G->n; i++) {
		H_rep_idx1[i] = -1;
		}
	for (i = 0; i < G->n; i++) {
		if (H_rep_idx1[i] != -1)
			continue;
		for (j = 0; j < H_len; j++) {
			k = autlog_mult_G(G, i, H[j]);
			if (H_rep_idx1[k] != -1)
				return error("autlog_do_factorgroup(): "
				"H_rep_idx1[k] != -1");
			H_rep_idx1[k] = GmH_n1;
			}
		H_rep1[GmH_n1++] = i;
		}
	*H_rep_idx = H_rep_idx1;
	*H_rep = H_rep1;
	*GmH_n = GmH_n1;
	if (GmH_n1 != G->n / H_len)
		return error("GmH_n != G->n / H_len");
	if (f_verbose) {
		printf("autlog_do_factorgroup(): "
		"found %ld cosets\n", GmH_n1);
		}
	
	GmH_gen1 = (INT *) my_malloc(G->nb_gen * sizeof(INT), "autlog_do_factorgroup");
	if (GmH_gen1 == NIL)
		return error("autlog_do_factorgroup(): "
		"no memory for GmH_gen1");
	for (i = 0; i < G->nb_gen; i++) {
		j = H_rep_idx1[G->g[i]];
		GmH_gen1[i] = j;
		}
	*GmH_gen = GmH_gen1;
	
	GmH_dim1 = GmH_n1 + 1;
	GmH1 = (INT *) my_malloc(
		GmH_n1 * GmH_dim1 * sizeof(INT), "autlog_do_factorgroup");
	if (GmH1 == NIL)
		return error("autlog_do_factorgroup(): "
		"no memory for GmH1");
	
	/* i, j, k: elements of G mod H, 
	 * i1, j1, k1: elements of G. */
	for (i = 0; i < GmH_n1; i++) {
		i1 = H_rep1[i];
		idx_inv = -1;
		for (j = 0; j < GmH_n1; j++) {
			j1 = H_rep1[j];
			k1 = autlog_mult_G(G, i1, j1);
			k = H_rep_idx1[k1];
			if (k == 0) {
				if (idx_inv != -1)
					return error("autlog_do_factorgroup"
					"(): idx_inv != -1");
				idx_inv = j;
				}
			GmH1[i * GmH_dim1 + j] = k;
			if (f_verbose) {
				printf("%ld ", k);
				}
			}
		if (idx_inv == -1)
			return error("autlog_do_factorgroup"
			"(): idx_inv == -1");
		GmH1[i * GmH_dim1 + GmH_n1] = idx_inv;
		if (f_verbose)
			printf("  | %ld\n", idx_inv);
		}
	*GmH_dim = GmH_dim1;
	*GmH = GmH1;
	return OK;
}

#if TEXDOCU
static INT autlog_derivedgroup(
	AUTLOG_THE_GROUP *G, 
	INT *H_rep_idx, INT *H_rep, 
	INT *GmH, INT GmH_n, INT GmH_dim,
	INT **C_idx, INT **C, INT *nb_C, 
	INT f_verbose)
#endif
{
	INT *C_idx1 = NIL;
	INT *C1 = NIL;
	INT nb_C1 = 0;
	INT i, j, k, i1, j1, k1, l, l1;
	
	C_idx1 = (INT *) my_malloc(GmH_n * GmH_n * sizeof(INT), "autlog_derivedgroup");
	C1 = (INT *) my_malloc(GmH_n * GmH_n * sizeof(INT), "autlog_derivedgroup");
	if (C_idx1 == NIL || C1 == NIL)
		return error("autlog_derivedgroup(): no memory for C_idx1/C1");
	for (i = 0; i < GmH_n; i++) {
		i1 = H_rep[i];
		for (j = 0; j < GmH_n; j++) {
			j1 = H_rep[j];
			k1 = autlog_commutator_G(G, i1, j1);
			for (l = 0; l < nb_C1; l++) {
				if (C1[l] >= k1)
					break;
				}
			if (nb_C1 == 0) {
				C1[nb_C1++] = k1;
				}
			else if (l == nb_C1) {
				C1[l] = k1;
				nb_C1++;
				}
			else if (C1[l] > k1) {
				/* raufschieben: */
				for (l1 = nb_C1; l1 > l; l1--)
					C1[l1] = C1[l1 - 1];
				C1[l] = k1;
				/* bisherige Eintraege aktualisieren: */
				for (l1 = 0; l1 < i * GmH_n + j; l1++)
					if (C_idx1[l1] >= l)
						C_idx1[l1]++;
				nb_C1++;
				}
			C_idx1[i * GmH_n + j] = l;
			}
		}
	if (f_verbose) {
		printf("Die Kommutatoren:\n");
		for (i = 0; i < nb_C1; i++)
			printf("%ld ", C1[i]);
		printf("\n");
		for (i = 0; i < GmH_n; i++) {
			for (j = 0; j < GmH_n; j++) {
				k = C_idx1[i * GmH_n + j];
				printf("%ld ", k);
				}
			printf("\n");
			}
		}
	*C_idx = C_idx1;
	*C = C1;
	*nb_C = nb_C1;
	return OK;
}

#if TEXDOCU
static INT autlog_derivedgroup_generators(
	AUTLOG_THE_GROUP *GmZ, 
	INT *C_idx, INT *C, INT nb_C, 
	INT **g_len, INT **g, INT f_verbose)
#endif
{
	INT *g_len1 = NIL;
	INT *g1 = NIL;
	INT nb_gen;
	INT k, go, go1, i;

	g_len1 = (INT *) my_malloc((GmZ->nb_gen + 1) * sizeof(INT), "autlog_derivedgroup_generators()");
	g1 = (INT *) my_malloc((nb_C + 1) * sizeof(INT), "autlog_derivedgroup_generators()");
	if (g_len1 == NIL || g1 == NIL)
		return error("autlog_derivedgroup_generators"
		"(): no memory for g_len1/g1");
	g1[0] = 0;
	nb_gen = 1;
	g_len1[0] = 1;
	for (k = 0; k < GmZ->nb_gen; k++) {
		go = autlog_the_group_order(GmZ, k);
		go1 = autlog_the_group_order(GmZ, k + 1);
		aldg_gen(GmZ, go, go1, 0, 
			go1, C_idx, C, g1, &nb_gen, f_verbose);
		aldg_gen(GmZ, 0, go, go, 
			go1, C_idx, C, g1, &nb_gen, f_verbose);
		g_len1[k] = nb_gen;
		if (f_verbose) {
			printf("k = %ld: ", k);
			for (i = (k == 0 ? 
				0 : g_len1[k - 1]); i < g_len1[k]; i++)
				printf("%ld ", g1[i]);
			printf("\n");
			}
		}
	if (nb_gen != nb_C)
		return error("autlog_derivedgroup_generators"
		"(): nb_gen != nb_C");
	g1[nb_gen] = -1;
		/* possibly addressed in 
		 * autlog_commutator_test() */
		/* therefore the + 1 in my_malloc() */
	*g_len = g_len1;
	*g = g1;
	return OK;
}

#if TEXDOCU
static INT aldg_gen(
	AUTLOG_THE_GROUP *GmZ, 
	INT i0, INT i1, INT j0, INT j1, 
	INT *C_idx, INT *C, INT *g, 
	INT *nb_gen, INT f_verbose)
#endif
{
	INT i, j, ii, jj, k, l;
	INT c_idx, c;
	
	for (i = i0; i < i1; i++) {
		ii = GmZ->pv.a[i];
		for (j = j0; j < j1; j++) {
			jj = GmZ->pv.a[j];
			c_idx = C_idx[ii * GmZ->n + jj];
			c = C[c_idx];
#if 0
			if (f_verbose) {
				printf("i = %ld j = %ld, "
					"ii = %ld jj = %ld, c_idx = %ld c = %ld\n", 
					i, j, ii, jj, c_idx, c);
				}
#endif
			for (k = 0; k < *nb_gen; k++) {
				if (g[k] == c)
					break;
				}
			if (k == *nb_gen) {
				/* nicht gefunden bzw. neue Liste */
				g[*nb_gen] = c;
				(*nb_gen)++;
				}
			}
		}
	return OK;
}


// #include "aut_grid.C"
/* aut_grid.C 
 * included from aut.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


/*
 *
 */

#define DEBUG_CALC_GRID
#undef DEBUG_CALC_GRID_VERY_VERBOSE

#if TEXDOCU
static INT autlog_calc_grid(
	AUTLOG_THE_GROUP *C, 
	AUTLOG_GRID *G, INT k, INT go, INT f_verbose)
#else
/* Relativordnungen bezueglich 
 * $E = U_{0}, U_1, ... U_{k}$
 * nach G berechnen
 * (k == 0 heisst nur bezueglich Einsgruppe). 
 * $U_k = <U_0, ... U_{k - 1}>$.
 * Es wird die Tabelle fuer alle Elemente 
 * zwischen go\_k und go aufgestellt.
 * Diese muessen mittels C-$>$p 
 * an den Beginn der Gruppentafel 
 * permutiert sein. */
#endif
{
	INT len, k1, i, j, i1, j1;
	INT ro, re, o, o1, t_i, n0 = 0;
	INT dim_n, *theG;
	INT go_k, go_k1;
	SPERM *p, *pv;
	
	dim_n = C->dim_n;
	theG = C->theG;
	p = &C->p;
	pv = &C->pv;
	
	if (G->m != C->n)
		return error("autlog_calc_grid(): G->m != C->n");
	sp_id(&G->r);
	sp_id(&G->rv);
	if (C->f_has_element_colors)
		n0 = 1;
	G->n = n0 + 1 + k * 2; /* (n0 for element coloring) */
	go_k = autlog_the_group_order(C, k);
	/* go_k war vorher a->n1[k]. */
	len = go - go_k;
	for (i = 0; i < len; i++)
		for (j = 0; j < G->n; j++)
			G->type[i][j] = 0;
	
	for (i = 0; i < go_k; i++)
		G->type_idx[i] = - 1;
	for ( ; i < go; i++)
		G->type_idx[i] = i - go_k;
#if 0
	for (i = 0; i < go_k; i++) {
		if (i >= go_km1 /* a->n1[k] */)
			G->type_idx[i] = i - go_km1 /* a->n1[k] */;
		else
			G->type_idx[i] = - 1;
		}
#endif
	
	for (i = go_k /* a->n1[k] */; i < go; i++) {
		i1 = pv->a[i];
		/* (i1/j1) der Eintrag in die 
		 * urspruengliche (unpermutierte) Gruppentafel. */
		t_i = G->type_idx[i] 
			/* i - go_km1 */ /* a->n1[k] */;
#ifdef DEBUG_CALC_GRID_VERY_VERBOSE
		if (f_verbose) {
			printf("go_k = %ld i = %ld i1 = %ld t_i = %ld\n", 
				go_k, i, i1, t_i);
			fflush(stdout);
			}
#endif
		if (C->f_has_element_colors)
			G->type[t_i][0] = - C->element_colors[i1];
		
		/* Neue Zeile in AUTLOG_GRID */
		/* Die relativen Elementordnungen
		 * mit den zugehoerigen Potenzelementen in allen 
		 * bereits festgelegten Untergruppen 
		 * (die Elementnummern dort sind also 
		 * kompatibel links wie rechts):
		 *
		 * type[][0]: element color
		 *
		 * type[][0 + n0] relative order ro bzgl. der 
		 *  i^ro < go(0) = 1 (also i^ro = id)
		 *
		 * type[][1 + n0] re = i^ro
		 * type[][2 + n0] ro bzgl. der 
		 *         re = i^ro < go(1) = |<g_0>|
		 *
		 *       [3 + n0]
		 *       [4 + n0]        go(2) = |<g_0, g_1>|
		 *  ...
		 * type[][2k-1 + n0] re  go(k) = |<g_0, ... , g_{k-1}>|
		 * type[][2k + n0]   ro
		 * 
		 * das sind zusammen (pro Element i) n0 + 2k + 1 Eintraege
		 * (= G->n), nach denen spaeter sortiert wird.
		 */
		j = i;
		ro = 1;
		for (k1 = k; k1 >= 0; k1--) {
			go_k1 = autlog_the_group_order(C, k1);
#ifdef DEBUG_CALC_GRID_VERY_VERBOSE
			if (f_verbose) {
				printf("k1 = %ld go_k1 = %ld\n", 
					k1, go_k1);
				fflush(stdout);
				}
#endif
			while (TRUE) {
				j1 = pv->a[j];
				o = theG[i1 * dim_n + j1];
				o1 = p->a[o];
				ro++; /* Invariante: o1 = i^ro */
				if (o1 < go_k1 /* a->n1[k1] */) {
					/* (muss irgendwann eintreten) */
					break;
					}
				/* o1 >= go_k1, also weiterpotenzieren: */
				j = o1;
				} /* while */
			/* jetzt: ro ist Relativordnung */
			re = o1; /* kompatible numerierung. */
			
			/* re/ro eintragen: */
			G->type[t_i][n0 + 2 * k1] = ro;
			if (k1 > 0)
				G->type[t_i][n0 + 2 * k1 - 1] = re;
			
#ifdef DEBUG_CALC_GRID_VERY_VERBOSE
			if (f_verbose) {
				printf("ro = %ld re = %ld\n", ro, re);
				fflush(stdout);
				}
#endif
			/* wir sind evtl. bereits in eine tiefere 
			 * Untergruppe gefallen: */
			while (k1 > 0 && 
				re < autlog_the_group_order(C, k1 - 1)) {
				k1--;
				
				/* re/ro eintragen: */
				G->type[t_i][n0 + 2 * k1] = ro;
				if (k1 > 0)
					G->type[t_i][n0 + 2 * k1 - 1] = re;
				}
			/* jetzt ist k1 die minimale Untergruppe 
			 * mit o1 < iso->n1[k1]. 
			 * wir machen weiter. */
			j = o1;
			} /* next k1 */
		} /* next i */
	G->G_max = 0;
	G->first[0] = go_k;

/* if (a->info->f_very_verbose && a->info->f_verbose) {
		printf("in autlog_calc_grid() (k = %ld):\n", k);
		for (i = 0; i < len; i++) {
			for (j = 0; j < G->n; j++) {
				printf("%ld ", G->type[i][j]);
				}
			printf("\n");
			}
		} */
	
#ifdef DEBUG_CALC_GRID
#if 0
	if (f_verbose) {
		printf("vor autlog_radix_sort():\n");
		autlog_print_grid(G, pv, k, go_k, go);
		}
#endif
#endif

	autlog_radix_sort(G, 0 /* radix */, 
		go_k /* a->n1[k] */, 
		go - 1 /* a->n1[k] + len - 1 */);
	
#ifdef DEBUG_CALC_GRID
	if (f_verbose) {
		printf("nach autlog_radix_sort():\n");
		autlog_print_grid(G, pv, k, go_k, go);
		}
#endif

	return OK;
}

#if TEXDOCU
static INT autlog_radix_sort(AUTLOG_GRID *G, INT radix, INT first, INT last)
#else
/* Die Typvektoren werden aufsteigend sortiert.
 * Most significant ist dabei die 0-te Stelle von type (radix = 0). */
#endif
{
	INT f_found, idx, k, l, t, i1, j, x, y;
	INT first0, first1, res, k1;
	
	if (first == last || radix == G->n) {
		/* Bereich first .. last 
		 * als neuen grid entry nach G eintragen. 
		 * grid_entry[first..last] setzen. */
		k = G->G_max;
		if (G->first[k] != first)
			return error("autlog_radix_sort(): G->first[k] != first");
		for (l = first; l <= last; l++)
			G->grid_entry[l] = k;
		G->len[k] = last - first + 1;
		G->first[k + 1] = last + 1;
		G->G_max++;
		/* if (a->info->f_very_verbose && 
			a->info->f_verbose) {
			printf("autlog_radix_sort()|new entry: first = %ld len = %ld : ", 
				(INT) first, (INT) G->len[k]);
			i1 = G->type_idx[first];
			for (j = 0; j < G->n; j++) {
				printf("%ld ", (INT) G->type[i1][j]);
				}
			printf("\n");
			} */
		return OK;
		}
	
	for (k = first; k <= last; k++) {
		f_found = autlog_insert_idx(G, first, k - first, radix, k, &idx);
		if (idx != k) {
			/* s = (idx idx+1 ... k) auf den aktuellen Stand 
			 * der Matrix anwenden. 
			 *   p := p * s, pv := s^-1 * pv */
			sp_mult_apply_forwc_r(&G->r, idx, k - idx + 1);
			sp_mult_apply_backwc_l(&G->rv, idx, k - idx + 1);
			t = G->type_idx[k];
			for (l = k; l > idx; l--)
				G->type_idx[l] = G->type_idx[l - 1];
			G->type_idx[idx] = t;
			/* grid_entry ist noch nicht gesetzt, 
			 * braucht nicht mitgeschoben zu werden. */
			}
		}
	first0 = first;
	first1 = G->type_idx[first0];
	for (k = first + 1; k <= last; k++) {
		k1 = G->type_idx[k];
		x = G->type[k1][radix];
		y = G->type[first1][radix];
		res = x - y;
		if (res < 0) /* y > x */
			return error("autlog_radix_sort(): descending");
		if (res > 0) { /* x > y */
			autlog_radix_sort(G, radix + 1, first0, k - 1);
			first0 = k;
			first1 = G->type_idx[first0];
			}
		if (k == last) {
			autlog_radix_sort(G, radix + 1, first0, k);
			}
		}
	return OK;
}

#if TEXDOCU
static INT autlog_insert_idx(AUTLOG_GRID *G, 
	INT first, INT len, INT radix, 
	INT search_this, INT *idx)
#else
/* aufsteigende Sortierung, 
 * bestimme Einfuegeposition fuer "search\_this" Element 
 */
#endif
{
	INT i, st1, cur, cur1, res, a, b;
	INT f_found;
	
	st1 = G->type_idx[search_this];
	f_found = FALSE;
	for (i = 0; i < len; i++) {
		cur = first + i;
		cur1 = G->type_idx[cur];
		a = G->type[cur1][radix];
		b = G->type[st1][radix];
		res = a - b;
		if (res == 0) /* a == b */
			f_found = TRUE;
		if (res > 0) { /* a > b */ /* aufsteigende Sortierung ! */
			*idx = cur;
			return f_found;
			}
		}
	*idx = first + len;
	return f_found;
}

#if TEXDOCU
static void autlog_print_grid(
	AUTLOG_GRID *G, SPERM *pv, 
	INT k, INT go_k, INT go)
#endif
{
	INT i, j, i0, i1, k1, t_i, ro, re, first;
	
#if 0
	printf("the grid (ro/re):\n");
	for (i = go_k; i < go; i++) {
		i1 = G->rv.a[i];
		i0 = pv->a[i1];
		t_i = G->type_idx[i];
		printf("i = %ld i1 = %ld i0 = %ld: %ld", 
			i, i1, i0, G->type[t_i][0]);
		for (k1 = 1; k1 <= k; k1++) {
			ro = G->type[t_i][2 * k1];
			re = G->type[t_i][2 * k1 - 1];
			printf(", (%ld/%ld)", ro, re);
			}
		printf("\n");
		}
#endif
	printf("G_max = %ld\n", G->G_max);
	for (i = 0; i < G->G_max; i++) {
		first = G->first[i];
		printf("first/len = %ld/%ld: ", G->first[i], G->len[i]);
		autlog_print_type_vector(G, first, k);
		}
}

#if TEXDOCU
static void autlog_print_type_vector(AUTLOG_GRID *G, INT i, INT k)
#endif
{
	INT t_i, k1, ro, re;
	
	t_i = G->type_idx[i];
	printf("%ld", G->type[t_i][0]);
	for (k1 = 1; k1 <= k; k1++) {
		ro = G->type[t_i][2 * k1];
		re = G->type[t_i][2 * k1 - 1];
		printf(", (%ld/%ld)", ro, re);
		}
	printf("\n");
}

#if TEXDOCU
static INT autlog_calc_all_grids(
	AUTLOG_THE_GROUP *C, 
	INT go, INT f_verbose)
#endif
{
	return autlog_calc_n_grids(
		C, go, C->nb_gen, f_verbose);
}

#if TEXDOCU
static INT autlog_calc_n_grids(
	AUTLOG_THE_GROUP *C, 
	INT go, INT n, INT f_verbose)
#endif
{
	INT k;
	
	if (C->nb_grid < n)
		return error("autlog_calc_n_grids(): "
		"C->nb_grid < n");
	for (k = 0; k < n; k++)
		autlog_calc_grid(C, C->Grid + k, 
			k, go, f_verbose);
	return OK;
}

#if TEXDOCU
static INT autlog_E(
	AUTLOG_THE_GROUP *C1, AUTLOG_THE_GROUP *C2, 
	AUTLOG_GRID *G1, AUTLOG_GRID *G2, 
	INT k, ORDERED_SET *E, INT f_verbose)
#endif
{
	if (!autlog_check_iso(G1, G2)) {
		E->size = 0;
		return OK;
		}
	autlog_choose_E(C1, C2, G1, G2, 
		k, E, f_verbose);
	return OK;
}

#if TEXDOCU
static INT autlog_check_iso(
	AUTLOG_GRID *G1, AUTLOG_GRID *G2)
#else
/* TRUE, wenn beide GRIDs gleich. */
#endif
{
	INT i, j, first, i1, i2;
	
	if (G1->G_max != G2->G_max)
		return FALSE;
	for (i = 0; i < G1->G_max; i++) {
		if (G1->first[i] != G2->first[i])
			return FALSE;
		if (G1->len[i] != G2->len[i])
			return FALSE;
		}
	for (i = 0; i < G1->G_max; i++) {
		first = G1->first[i];
		i1 = G1->type_idx[first];
		i2 = G2->type_idx[first];
		for (j = 0; j < G1->n; j++) {
			if (G1->type[i1][j] != G2->type[i2][j])
				return FALSE;
			}
		}
	/* now: a bijection exists. */
	return TRUE;
}

#if TEXDOCU
static INT autlog_choose_first_into_E(AUTLOG_THE_GROUP *G, 
	AUTLOG_GRID *Grid, INT k, ORDERED_SET *E, INT f_verbose)
#endif
{
	INT i, j, l, first, len;
	INT g0, g1, g2, ge, go, i0, i1, i2;
		
	go = autlog_the_group_order(G, k); /* before: n1 */
	if (go != Grid->first[0])
		return error("autlog_choose_first_into_E(): go != Grid->first[0]");
	ge = Grid->G_max - 1;
	first = Grid->first[ge];
	len = Grid->len[ge];
	E->size = 0;
	for (i = 0; i < len; i++) {
		i2 = first + i;
		i1 = Grid->rv.a[i2];
		i0 = G->pv.a[i1];
			/* die urspruengliche Elementnummer. */
		for (j = 0; j < E->size; j++)
			if (E->a[j] > i0)
				break;
		for (l = E->size - 1; l >= j; l--)
			E->a[l + 1] = E->a[l];
		E->a[j] = i0;
		E->size++;
		/* E->a[E->size++] = i2; */
		}
	if (f_verbose) {
		autlog_print_E(E, k);
		}
	return OK;
}

#if TEXDOCU
static INT autlog_choose_E(
	AUTLOG_THE_GROUP *C1, AUTLOG_THE_GROUP *C2, 
	AUTLOG_GRID *G1, AUTLOG_GRID *G2, 
	INT k, ORDERED_SET *E, INT f_verbose)
#endif
{
	INT i, j, l, first;
	INT g0, g1, g2, ge, go, i0, i1, i2;
		
	go = autlog_the_group_order(C1, k); /* before: n1 */
	if (go != autlog_the_group_order(C2, k))
		return error("C1->go[] != C2->go[]");
	/* n1 = a->n1[k]; */
	if (go != G1->first[0])
		return error("autlog_choose_E(): go != G1->first[0]");
	if (go != G2->first[0])
		return error("autlog_choose_E(): go != G2->first[0]");
	if (C1->nb_gen <= k)
		return error("autlog_choose_E(): C1->nb_gen <= k");
#if 0
	if (C1->nb_gen <= k) {
		g1 = go;
			/* das naechste freie Element 
			 * (in der C1->p Lage) */
		g2 = G1->r.a[g1];
		g0 = C1->pv.a[g1];
		if (f_verbose) {
			printf("choosing generator "
			"%ld now in %ld\n", g0, g2);
			}
		C1->g[C1->nb_gen] = g0;
		C1->nb_gen++;
		}
#endif
	g0 = C1->g[k];
	g1 = C1->p.a[g0];
	g2 = G1->r.a[g1];
	/* das g0 - te Gruppenelement liegt momentan 
	 * in der g2 - ten Zeile / Spalte
	 * (bezogen auf die sortierte 
	 * Reihenfolge des grids). */
	ge = G1->grid_entry[g2];
	/* links (A): nimm den Block, in dem g1 liegt. 
	 * rechts (B): alle Bloecke des Grideintrags ge
	 * nach E[k]. */
	E->size = 0;
	first = G2->first[ge];
	for (i = 0; i < G2->len[ge]; i++) {
		i2 = first + i;
		i1 = G2->rv.a[i2];
		i0 = C2->pv.a[i1];
			/* die urspruengliche Elementnummer. */
		for (j = 0; j < E->size; j++)
			if (E->a[j] > i0)
				break;
		for (l = E->size - 1; l >= j; l--)
			E->a[l + 1] = E->a[l];
		E->a[j] = i0;
		E->size++;
		/* E->a[E->size++] = i2; */
		}
	if (f_verbose) {
		autlog_print_E(E, k);
		}
	return TRUE;
}

#if TEXDOCU
static void autlog_print_E(ORDERED_SET *E, INT k)
#endif
{
	INT l;
	
	printf("E[%ld] = {", k);
	for (l = 0; l < E->size; l++) {
		printf("%ld", E->a[l]);
		if (l < E->size - 1)
			printf(", ");
		}
	printf("}\n");
}


// #include "aut_dimino.C"
/* aut_dimino.C 
 * included from aut.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#if TEXDOCU
static INT autlog_group_gen(AUTLOG_THE_GROUP *C, 
	INT nb_gen, INT *g, INT **g_idx, INT f_verbose)
#endif
{
	INT k, last_go, go, g_k, g_k1;
	INT *g_idx1 = NIL;
	
	sp_id(&C->p);
	sp_id(&C->pv);
	sp_id(&C->q);
	sp_id(&C->qv);
	C->nb_gen = 0;
	g_idx1 = (INT *) my_malloc(nb_gen * sizeof(INT), "autlog_group_gen()");
	if (g_idx1 == NIL)
		return error("autlog_group_gen(): "
		"no memory for g_idx1");
	for (k = 0; k < nb_gen; k++) {
		g_k = g[k];
		g_k1 = C->p.a[g_k];
		last_go = autlog_the_group_order(
			C, C->nb_gen);
		if (g_k1 < last_go) {
			g_idx1[k] = -1;
				/* redundant generator */
			if (f_verbose) {
				printf("%ld - th generator "
					"%ld skipped\n", 
					k, g_k);
				fflush(stdout);
				}
			continue;
			}
		
		autlog_dimino(C, C->nb_gen, g_k, 
			FALSE /* f_test_it */, 
			0 /* go_soll */ /* before: A_N2 */, 
			NIL /* A */, NIL /* C_grid */, 
			NIL /* A_grid */, 
			&go);
		
		C->g[C->nb_gen] = g_k;
		C->go[C->nb_gen] = go;
		g_idx1[k] = C->nb_gen;
			/* new generator number */
		
		if (f_verbose) {
			printf("autlog_group_gen(): "
				"g[%ld] = %ld, go = %ld\n", 
				C->nb_gen, C->g[C->nb_gen], 
				C->go[C->nb_gen]);
			fflush(stdout);
			}
		C->nb_gen++;
		sp_mult(&C->p, &C->q, &C->tmp1);
		sp_mult(&C->qv, &C->pv, &C->tmp2);
		sp_mv(&C->tmp1, &C->p);
		sp_mv(&C->tmp2, &C->pv);
		sp_id(&C->q);
		sp_id(&C->qv);
		}
	*g_idx = g_idx1;
	return OK;
}

#if TEXDOCU
static INT autlog_add_generator(
	AUTLOG_THE_GROUP *B, 
	INT i, INT g, INT f_verbose)
#endif
{
	INT g1, go, last_go;
	
	if (f_verbose) {
		printf("autlog_add_generator(): "
		"i = %ld generator %ld\n", i, g);
		}
	g1 = B->p.a[g];
	last_go = autlog_the_group_order(B, i);
	if (g1 < last_go) {
		printf("i = %ld g = %ld g1 = %ld "
			"last_go = %ld B->nb_gen = %ld "
			"|B| = %ld\n", 
			i, g, g1, last_go, B->nb_gen, 
			autlog_the_group_order(B, B->nb_gen) );
		printf("error in autlog_add_generator(): "
		"g1 < last_go"); fflush(stdout);
		return ERROR;
		}
	sp_id(&B->q);
	sp_id(&B->qv);

	autlog_dimino(B, i, g, 
		FALSE /* f_test_it */, 
		0 /* go_soll */, 
		NIL /* A */, NIL /* C_grid */, 
		NIL /* A_grid */, 
		&go);
	B->g[i] = g;
	B->go[i] = go;

	sp_mult(&B->p, &B->q, &B->tmp1);
	sp_mult(&B->qv, &B->pv, &B->tmp2);
	sp_mv(&B->tmp1, &B->p);
	sp_mv(&B->tmp2, &B->pv);
	sp_id(&B->q);
	sp_id(&B->qv);
	
	return OK;
}

#if TEXDOCU
static INT autlog_reduce(AUTLOG_LOCAL *L, 
	AUTLOG_THE_GROUP *A, 
	AUTLOG_THE_GROUP *B, 
	AUTLOG_GRID *G1, AUTLOG_GRID *G2, 
	ORDERED_SET *E, INT k, INT f_verbose)
#endif
{
	INT go_A, go_B, go_A0, go_B0;
	INT g0, l, m, f_ok, offset, ret;

	if (E->size == 0)
		return TRUE;
	go_A = autlog_the_group_order(A, k + 1);
	if (f_verbose) {
		printf("autlog_reduce(): go_A = %ld\n", go_A);
		}
	go_A0 = autlog_the_group_order(A, k);
	go_B0 = autlog_the_group_order(B, k);
	if (go_A0 != go_B0)
		return error("autlog_reduce(): "
		"go_A0 != go_B0");
	m = 0;
	while (m < E->size) {
		g0 = E->a[m];
		sp_id(&B->q);
		sp_id(&B->qv);

		f_ok = FALSE;
		if (autlog_dimino(B, k, g0, TRUE /* f_test_it */, 
			go_A /* go_soll */ /* before: A_N2 */, 
			A, G2 /* C_grid */, G1 /* A_grid */, 
			&go_B)) {
			/* now test the group table: */
			if (aut_table_test(A, B, go_A0, go_A, 
				go_A0 /* ib */, go_A /* ie */, 
				0 /* jb */, go_A /* je */, f_verbose)) {
				if (aut_table_test(A, B, go_A0, go_A, 
					0 /* ib */, go_A0 /* ie */, 
					go_A0 /* jb */, go_A /* je */, f_verbose)) {
					
					
					if (L->info->f_autologisms) {
						offset = 0;
						ret = autlog_commutator_test(
							L, A, B, k, &offset, 
							go_A0 /* ib */, go_A /* ie */, 
							0 /* jb */, 
							go_A /* je */, f_verbose);
						if (ret == -1) {
	printf("autlog_reduce() error in "
	"autlog_commutator_test() (first time)\n");
	printf("k = %ld offset = %ld\n", k, offset);
	fflush(stdout);
	return error("1");
							}
						if (ret) {
							ret = autlog_commutator_test(
								L, A, B, k, &offset, 
								0 /* ib */, go_A0 /* ie */, 
								go_A0 /* jb */, 
								go_A /* je */, f_verbose);
							if (ret == -1) {
	printf("autlog_reduce() error in "
	"autlog_commutator_test() (second time)\n");
	printf("k = %ld offset = %ld\n", k, offset);
	fflush(stdout);
	return error("2");
								}
							if (ret) {
								f_ok = TRUE;
								}
							}
						}
					else
						f_ok = TRUE;
					
					}
				}
			}
		if (!f_ok) {
			for (l = m + 1; l < E->size; l++)
				E->a[l - 1] = E->a[l];
			E->size--;
			}
		else {
			m++;
			}
		}
	if (f_verbose) {
		autlog_print_E(E, k);
		}
	return TRUE;
}

#if TEXDOCU
static INT autlog_commutator_test(AUTLOG_LOCAL *L, 
	AUTLOG_THE_GROUP *A, 
	AUTLOG_THE_GROUP *B, 
	INT k, INT *offset, 
	INT ib, INT ie, INT jb, INT je, INT f_verbose)
#else
/* A = AmZ, B = BmZ. 
 * returns FALSE, if no autologism possible. */
#endif
{
	INT Ai0, Bi0, Aj0, Bj0;
	INT A_C_idx, B_C_idx, A_c, B_c, A_c1, B_c1;
	INT go_Ad, go_Bd;
	INT i, i1, j1;
	INT g_next, last_g_idx, last_go, go;
	/* g_next + offset is index of next (unprocessed) 
	 * generator in AcA->g[]. 
	 * This may be equal to AcA->g_len[AcA->nb_gen - 1], 
	 * in this case it is not a useful generator. */
	
	if (k == 0)
		g_next = 1;
		/* skip first generator id in AcA->g[]. */
	else
		g_next = L->AcA.g_len[k - 1];
	last_g_idx = -1;
	for (i = 0; i < g_next + *offset; i++) {
		if (L->AcA.g_idx[i] != -1) {
			if (L->AcA.g_idx[i] != last_g_idx + 1)
				return error("autlog_commutator_test(): "
				"L->AcA.g_idx[i] != last_g_idx + 1");
			last_g_idx = L->AcA.g_idx[i];
			if (L->Ad.g[last_g_idx] != L->AcA.g[i])
				return error("autlog_commutator_test(): "
				"L->Ad.g[last_g_idx] != L->AcA.g[i]");
			}
		}
	last_go = autlog_the_group_order(
		&L->Ad, last_g_idx + 1);
		/* works fine even if last_g_idx is -1 */
	for (i1 = ib; i1 < ie; i1++) {
		Ai0 = A->pv.a[ A->qv.a[i1] ];
		Bi0 = B->pv.a[ B->qv.a[i1] ];
		for (j1 = jb; j1 < je; j1++) {
			Aj0 = A->pv.a[ A->qv.a[j1] ];
			Bj0 = B->pv.a[ B->qv.a[j1] ];
			A_C_idx = L->AcA.C_idx[Ai0 * L->AcA.dim_n + Aj0];
			B_C_idx = L->BcB.C_idx[Bi0 * L->BcB.dim_n + Bj0];
			A_c = L->AcA.C[A_C_idx];
			B_c = L->BcB.C[B_C_idx];
			A_c1 = L->Ad.p.a[A_c];
			B_c1 = L->Bd.p.a[B_c];
			if (f_verbose) {
				printf("autlog_commutator_test() "
				"A_c = %ld A_c1 = %ld "
				"B_c = %ld B_c1 = %ld\n", 
				A_c, A_c1, B_c, B_c1);
				fflush(stdout);
				}
			if (A_c1 < last_go) {
				if (A_c1 != B_c1) {
					if (f_verbose) {
						printf("autlog_commutator_test() "
						"A_c1 != %ld B_c1 = %ld - no iso !\n", 
						A_c1, B_c1);
						fflush(stdout);
						}
					return FALSE;
					}
				/* so far OK, now have to check 
				 * if this generator 
				 * has to be skipped in AcA.g[]
				 * (it is a non generator, but it could 
				 * have been occured earlier in g[]): */
				if (L->AcA.g[g_next + *offset] == A_c) {
				
					/* it is a non generator: */
					if (L->AcA.g_idx[g_next + *offset] != -1)
						return error("autlog_commutator_test(): "
						"L->AcA.g_idx[g_next + *offset] != -1");

					/* skip it: */
					(*offset)++;
					}
				}
			else {
				/* a new generator for 
				 * Ad, Bd has been found. */
				if (L->AcA.g[g_next + *offset] != A_c)
					return error("autlog_commutator_test(): "
					"L->AcA.g[g_next + *offset] != A_c");
				if (L->Ad.g[last_g_idx + 1] != A_c) {
					return error("autlog_commutator_test(): "
					"L->Ad.g[last_g_idx + 1] != A_c");
					}
				if (f_verbose) {
					printf("autlog_commutator_test(): "
					"a new generator\n");
					fflush(stdout);
					}
				if (B_c1 < autlog_the_group_order(
					&L->Bd, last_g_idx + 1)) {
					if (f_verbose) {
						printf("autlog_commutator_test(): "
						"B_c1 < last_go - no - iso !");
						fflush(stdout);
						}
					return FALSE;
					}
				
				if (autlog_add_generator(&L->Bd, 
					last_g_idx + 1, 
					B_c, f_verbose) != OK) {
	printf("autlog_commutator_test(): "
		"error in autlog_add_generator()\n");
	printf("i1 = %ld j1 = %ld Ai0 = %ld Aj0 = %ld "
		"Bi0 = %ld Bj0 = %ld\n", 
		i1, j1, Ai0, Aj0, Bi0, Bj0);
	printf("A_C_idx = %ld B_C_idx = %ld "
		"A_c = %ld B_c = %ld A_c1 = %ld B_c1 = %ld\n", 
		A_C_idx, B_C_idx, A_c, B_c, A_c1, B_c1);
	printf("A->n = %ld B->n = %ld\n", A->n, B->n);
					fflush(stdout);
					return -1;
					}
				go_Ad = autlog_the_group_order(
					&L->Ad, last_g_idx + 2);
				go_Bd = autlog_the_group_order(
					&L->Bd, last_g_idx + 2);
				if (go_Ad != go_Bd) {
					if (f_verbose) {
						printf("go_Ad = %ld != "
							"go_Bd = %ld - no iso !\n", 
							go_Ad, go_Bd);
						fflush(stdout);
						}
					return FALSE;
					}
				if (!aut_table_test(
					&L->Ad, &L->Bd, last_go, go_Ad, 
					last_go /* ib */, go_Ad /* ie */, 
					0 /* jb */, go_Ad /* je */, f_verbose)) {
					if (f_verbose) {
						printf("autlog_commutator_test(): "
						"aut_table_test() - no iso !\n");
						fflush(stdout);
						}
					return FALSE;
					}
				if (!aut_table_test(
					&L->Ad, &L->Bd, last_go, go_Ad, 
					0 /* ib */, last_go /* ie */, 
					last_go /* jb */, 
					go_Ad /* je */, f_verbose)) {
					if (f_verbose) {
						printf("autlog_commutator_test(): "
						"aut_table_test() - no iso !\n");
						fflush(stdout);
						}
					return FALSE;
					}
				(*offset)++;
				last_g_idx++;
				last_go = autlog_the_group_order(
					&L->Ad, last_g_idx + 1);
				}
			}
		}
	return TRUE;
}

#if TEXDOCU
static INT aut_table_test(
	AUTLOG_THE_GROUP *A, 
	AUTLOG_THE_GROUP *B, 
	INT go_last, INT go, 
	INT ib, INT ie, INT jb, INT je, INT f_verbose)
#else
/* Rueckgabe FALSE, wenn Differenzen 
 * in den Gruppentafeln von A und B 
 * im Rechteck (ib, jb) bis (ie, je) auftreten. 
 * Es werden als Elemente im Rechteck die durch 
 * p und q permutierten Elemente angenommen. */
#endif
{
	INT Ai0, Aj0, Ak0, Bi0, Bj0, Bk0, i1, j1, Ak1, Bk1;
	INT *AtheG, *BtheG, Adim_n, Bdim_n;
	SPERM *Ap, *Apv, *Aq, *Aqv;
	SPERM *Bp, *Bpv, *Bq, *Bqv;
	
	AtheG = A->theG;
	BtheG = B->theG;
	Adim_n = A->dim_n;
	Bdim_n = B->dim_n;
	Ap = &A->p;
	Apv = &A->pv;
	Aq = &A->q;
	Aqv = &A->qv;
	Bp = &B->p;
	Bpv = &B->pv;
	Bq = &B->q;
	Bqv = &B->qv;
	for (i1 = ib; i1 < ie; i1++) {
		Ai0 = Apv->a[ Aqv->a[i1] ];
		Bi0 = Bpv->a[ Bqv->a[i1] ];
		for (j1 = jb; j1 < je; j1++) {
			Aj0 = Apv->a[ Aqv->a[j1] ];
			Bj0 = Bpv->a[ Bqv->a[j1] ];
			Ak0 = AtheG[Ai0 * Adim_n + Aj0];
			Bk0 = BtheG[Bi0 * Bdim_n + Bj0];
			Ak1 = Aq->a[ Ap->a[Ak0] ];
			Bk1 = Bq->a[ Bp->a[Bk0] ];
			if (Ak1 >= go) {
				if (f_verbose)
					printf("(%ld/%ld): "
					"Ak1 = %ld > go = %ld - no iso !\n", 
						i1, j1, Ak1, go);
				return FALSE;
				}
			if (Bk1 >= go) {
				if (f_verbose)
					printf("(%ld/%ld): "
					"Bk1 = %ld > go = %ld - no iso !\n", 
						i1, j1, Bk1, go);
				return FALSE;
				}
			if (Ak1 != Bk1) {
				if (f_verbose)
					printf("(%ld/%ld): "
					"Ak1 = %ld != Bk1 = %ld - no iso !\n", 
						i1, j1, Ak1, Bk1);
				return FALSE;
				}
			}
		}
	return TRUE;
}

#if TEXDOCU
static INT autlog_dimino(
	AUTLOG_THE_GROUP *C, 
	INT k, INT g0, 
	INT f_test_it, INT go_soll /* before: A_N2 */, 
	AUTLOG_THE_GROUP *A, 
	AUTLOG_GRID *C_grid, AUTLOG_GRID *A_grid, 
	INT *go)
#else
/* Rueckgabe FALSE nur wenn f\_test\_it TRUE. */
#endif
{
	INT g1, coset_size, rep_pos;
	INT go1, h, s0, s1, s2, l;
	
	coset_size = 
		autlog_the_group_order(C, k);
		/* old_coset_size */
	go1 = coset_size;
		/* the old group size, 
		 * will grow to the new group size. */
	rep_pos = 0;
	g1 = C->p.a[g0];
	if (C->q.a[g1] != g1) {
		printf("coset_size = %ld\n", coset_size);
		printf("g0 = %ld\n", g0);
		printf("g1 = %ld\n", g1);
		printf("C->q.a[g1] = %ld\n", (INT) C->q.a[g1]);
		return error("aut_dimino(): C->q.a[g1] != g1");
		}
	if (!autlog_add_coset(C, f_test_it, 
		A, C_grid, A_grid,  
		coset_size, g1, &go1)) {
		/* only if f_test_it TRUE */
		return FALSE;
		}
	rep_pos += coset_size;
	while (rep_pos < go1) {
		for (l = 0; l <= k; l++) {
			if (l == k) {
				s0 = g0;
				}
			else {
				s0 = C->g[l];
				}
			s1 = C->p.a[s0];
			s2 = C->q.a[s1];
			h = autlog_mult(C, rep_pos, s2);
			/* printf("%ld * %ld = %ld\n", 
				rep_pos, s2, h); */
			if (h >= go1) {
				if (f_test_it && go1 >= go_soll) {
					/* printf("dimino produces too "
					"many cosets - no iso !\n"); */
					return FALSE;
					}
				if (!autlog_add_coset(C, f_test_it, 
					A, C_grid, A_grid,  
					coset_size, h, &go1)) {
					/* only if f_test_it TRUE */
					return FALSE;
					}
				}
			} /* for l */
		rep_pos += coset_size;
		}
	if (f_test_it && go1 != go_soll) {
		/* printf("end of dimino go_soll = %ld, "
		"go = %ld - no iso !\n", 
				go_soll, go1); */
		return FALSE;
		}
	*go = go1;
	return TRUE;
}

#if TEXDOCU
static INT autlog_add_coset(
	AUTLOG_THE_GROUP *C, 
	INT f_test_it, AUTLOG_THE_GROUP *A, 
	AUTLOG_GRID *C_grid, AUTLOG_GRID *A_grid, 
	INT coset_size, INT g2, INT *go)
#else
/* Annahme, dass die grids auf die Elemente 
 * in der Reihenfolge vor Anwendung von 
 * (C-$>$q sowie A-$>$q)
 * berechnet worden sind 
 * (nur falls f\_test\_it TRUE). 
 * Die erfolgten Permutationen 
 * werden nach q akkumuliert. 
 * Rueckgabe FALSE kann nur bei 
 * f\_test\_it TRUE auftreten. */
#endif
{
	INT i, j, x, y, x1, y1;
	INT A_ge, B_ge, g0, g1;
	
	g1 = C->qv.a[g2];
	g0 = C->pv.a[g1];
	if (g2 < *go) {
		printf("g2 = %ld g1 = %ld "
		"g0 = %ld *go = %ld\n", 
			g2, g1, g0, *go);
		return error("autlog_add_coset(): g2 < *go");
		}
	for (i = 0; i < coset_size; i++) {
		j = autlog_mult(C, i, g2);
		if (j < *go) {
			printf("i = %ld g2 = %ld g1 = %ld "
				"g0 = %ld qv->a[g2] = %ld\n", 
				i, g2, g1, g0, (INT) C->qv.a[g2]);
			printf("j = %ld *go = %ld\n", j, *go);
			if (C->fp_txt) {
				fprintf(C->fp_txt, "autlog_add_coset(): j < *go"
					"i = %ld g2 = %ld g1 = %ld "
					"g0 = %ld qv->a[g2] = %ld\n", 
					i, g2, g1, g0, (INT) C->qv.a[g2]);
				}
			return error("autlog_add_coset(): j < *go");
			}
		if (f_test_it) {
			/* q - Permutation in A und C 
			 * rueckgaengig machen: */
			x = A->qv.a[*go];
			y = C->qv.a[j];
			/* r - Permutation anwenden, welche 
			 * zur Sortierung im grid gehoert 
			 * (und von der p - Lage ausgeht). */
			x1 = A_grid->r.a[x];
			y1 = C_grid->r.a[y];
			A_ge = A_grid->grid_entry[x1];
			B_ge = C_grid->grid_entry[y1];
			if (A_ge < 0)
				return error("autlog_add_coset(): "
					"A_ge < 0");
			if (B_ge < 0)
				return error("autlog_add_coset(): "
					"B_ge < 0");
			if (A_ge != B_ge) {
				/* printf("A_ge = %ld != B_ge = "
					"%ld - no iso\n", A_ge, B_ge); */
				return FALSE;
				}
			}
		if (j != *go) {
			sp_mult_apply_tau_r(&C->q /* p */, j, *go);
			sp_mult_apply_tau_l(&C->qv /* pv */, j, *go);
			/* have to update g2, too: */
			if (g2 == j)
				g2 = *go;
			else if (g2 == *go)
				g2 = j;
			}
		(*go)++;
		}
	return TRUE;
}

#if TEXDOCU
static INT autlog_mult(
	AUTLOG_THE_GROUP *A, INT i, INT j)
#endif
{
	INT i0, j0, k0, i1, j1, k1, k2;
	
	i1 = A->qv.a[i];
	j1 = A->qv.a[j];
	i0 = A->pv.a[i1];
	j0 = A->pv.a[j1];
	k0 = A->theG[i0 * A->dim_n + j0];
	k1 = A->p.a[k0];
	k2 = A->q.a[k1];
	return k2;
}

#if TEXDOCU
static INT autlog_mult_G(
	AUTLOG_THE_GROUP *G, INT i, INT j)
#endif
{
	return G->theG[i * G->dim_n + j];
}

#if TEXDOCU
static INT autlog_invers_G(
	AUTLOG_THE_GROUP *G, INT i)
#endif
{
	return G->theG[i * G->dim_n + G->n];
}

#if TEXDOCU
static INT autlog_conjugate_G(
	AUTLOG_THE_GROUP *G, INT i, INT j)
#else
/* $j^{-1} * i * j$ */
#endif
{
	INT jv, k, l;
	
	jv = autlog_invers_G(G, j);
	k = autlog_mult_G(G, jv, i);
	l = autlog_mult_G(G, k, j);
	return l;
}

#if TEXDOCU
static INT autlog_commutator_G(
	AUTLOG_THE_GROUP *G, INT i, INT j)
#else
/* $i^{-1} *  j^{-1} * i * j$ */
#endif
{
	INT iv, jv, k, l;
	
	iv = autlog_invers_G(G, i);
	jv = autlog_invers_G(G, j);
	k = autlog_mult_G(G, iv, jv);
	l = autlog_mult_G(G, k, i);
	k = autlog_mult_G(G, l, j);
	return k;
}

#if 0
static INT aut_table_compare(
	AUTLOG_THE_GROUP *A, 
	INT go_last, INT go, 
	INT ib, INT ie, INT jb, INT je, 
	SPERM *p, SPERM *pv, SPERM *p_, SPERM *pv_, 
	INT f_verbose)
/* Rueckgabe FALSE, wenn Differenzen 
 * in den Gruppentafeln von A und B 
 * im Rechteck (ib, jb) bis (ie, je) auftreten. 
 * Es werden als Elemente im Rechteck die durch 
 * p und q permutierten Elemente angenommen. */
{
	INT i1, i0, i0_, j1, j0, j0_, k0, k0_, k1, k1_;
	INT *theG, dim_n;
	
	theG = A->theG;
	dim_n = A->dim_n;
	for (i1 = ib; i1 < ie; i1++) {
		i0 = pv->a[i1];
		i0_ = pv_->a[i1];
		for (j1 = jb; j1 < je; j1++) {
			j0 = pv->a[j1];
			j0_ = pv_->a[j1];
			k0 = theG[i0 * dim_n + j0];
			k0_ = theG[i0_ * dim_n + j0_];
			k1 = p->a[k0];
			k1_ = p_->a[k0_];
			if (k1 >= go) {
				return error("k1 >= go");
				}
			if (k1_ >= go) {
				return error("k1_ >= go");
				}
			if (k1_ > k1) {
				return 1;
				}
			if (k1_ < k1) {
				return -1;
				}
			}
		}
	return TRUE;
}
#endif


#endif /* SOLVABLE_TRUE */




