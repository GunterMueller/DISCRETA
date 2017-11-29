/* gt_canon.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/solvable.h>
#include <DISCRETA/lb.h>

#include <DISCRETA/SOLVABLE/gtP.h>
/* private include file ! */

#undef GROUP_CANON_TG

#ifdef GROUP_CANON_TG
#include <TG/group_canon.h>
#include <TG/sims.h>
#endif


#undef PRINT_OCCASIONALLY
#define MIN_LEVEL_FOR_PRINT 3
#define PRINT_INTERVALL_SEC 5

#define PRINT_BACKTRACK_POINTS
#define BACKTRACK_MOD 250

static INT last_print_ticks, nb_backtrack, nb_backtrack_points;

#if TEXDOCU
INT group_calc_hash(FILE *fp_txt, GROUP_TABLE *G, 
	VECTOR_OP hash, PERMUTATION_OP p0, INT f_v, INT f_vv)
#endif
{
	GROUP_CANONIC_FORM_OB cf;
	INT n;

	n = G->n;
	group_calc_cf(fp_txt, G, &cf, f_v, f_vv);
	if (f_vv) {
		fprintf(fp_txt, "computing hash value...\n");
		fflush(fp_txt);
		}
	cf.calc_hash(fp_txt, hash, n, f_v, f_vv);
	if (f_vv) {
		fprintf(fp_txt, "finished!\n");
		fflush(fp_txt);
		}

	cf.s_p0()->copy(p0);
#if 0
	if (f_v) {
		printf("canonical labelling p0 = \n");
		p0->println();
		fflush(stdout);
		}
#endif
	return OK;
}

#ifdef GROUP_CANON_TG
#if TEXDOCU
static INT calc_aut_M(INT n, LABRA_TG& GCAUT, VECTOR_OP base, 
	MATRIX_OP M, VECTOR_OP T, SYM_OP ago, INT f_v, INT f_vv)
#endif
{
	INT i, j, k, l, a, len, g_i, g_i_image;
	SYM_OB ago1, b;
	VECTOR_OP Ti;
	PERMUTATION_OP perm;
	PERMUTATION_OB p;

	LABRA_TG TEST_LAB;
	TEST_LAB.Init(GCAUT.Dim());
	TEST_LAB = GCAUT;
	TEST_LAB.CycleId();
	SIMS SI;
	SI.Init(GCAUT.Dim());
	SI = TEST_LAB;
	if (f_vv)
		SI.Print();
	
	len = base->s_li();
	M->m_ilih(n, len);
	for (i = 0; i < len; i++) {
		g_i = base->s_ii(i);
		if (f_vv) {
			printf("\ntransversal %ld, g_i = %ld:\n", i, g_i);
			fflush(stdout);
			}
		ARRAY < PERMUT < short > > & trans = SI[g_i + 1];
		l = trans.Used();
		for (j = 0; j < l; j++) {
			PERMUT <short > & p_tg = trans[j + 1];
			p.m_il(n);
			for (k = 0; k < n; k++) {
				a = p_tg[k + 1];
				p.m_ii(k, a);
				}
			if (f_vv) {
				p.println();
				fflush(stdout);
				}
			g_i_image = p.s_ii(g_i) - 1;
			if (f_vv) {
				printf("%ld \\mapsto %ld\n", g_i, g_i_image);
				fflush(stdout);
				}
			p.swap((PERMUTATION_OP) M->s_ij(i, g_i_image));
			}
		}

	T->m_il(len);
	((INTEGER_OP) &ago1)->m_i(1);
	for (i = 0; i < len; i++) {
		if (f_vv)
			printf("orbit of %ld:\n", base->s_ii(i));
		Ti = (VECTOR_OP) T->s_i(i);
		Ti->m_il(0);
		l = 0;
		for (j = 0; j < n; j++) {
			if (M->s_ij(i, j)->s_obj_k() != EMPTY) {
				perm = (PERMUTATION_OP) M->s_ij(i, j);
				if (f_vv)
					printf("%ld ", j);
				/* perm->println(); */
				Ti->inc();
				Ti->m_ii(l, j);
				l++;
				}
			}
		/* ago *= l; */
		((INTEGER_OP) &b)->m_i(l);
		b.mult_apply(&ago1);
		if (f_vv) {
			printf("length %ld ago = ", l);
			ago1.println();
			}
		}
	ago1.swap(ago);
	
	return OK;
}

#if TEXDOCU
INT group_calc_aut(FILE *fp_txt, GROUP_TABLE *G, VECTOR_OP base, 
	MATRIX_OP M, VECTOR_OP T, SYM_OP ago, INT f_v, INT f_vv)
#endif
{
	GROUP_CANON GC;
	ARRAY < VEKTOR < short > > tafel;
	INT n, i, j, a;
	INT t0, t1, user_time;
	BYTE str[256];

	n = G->n;
	tafel.REALLOC(n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = G->table->s_iji(i, j);
			tafel[i + 1][j + 1] = a + 1;
			}
		}
	if (f_v) {
		printf("group_calc_aut(): calling GC.Init(), GC.Run()...\n");
		fflush(stdout);
		}
	
	t0 = os_ticks();
	GC.Init(n, &tafel);
	GC.Run();
	
	LABRA_TG& GCAUT = GC.Aut();
	
	calc_aut_M(n, GCAUT, base,  M, T, ago, f_v, f_vv);

	t1 = os_ticks();
	user_time = t1 - t0;
	if (f_v) {
		str[0] = 0;
		print_delta_time(user_time, str);
		printf("found aut group in %s\n", str);
		fflush(stdout);
		}

	return OK;
}
#endif

#if TEXDOCU
INT group_calc_cf(FILE *fp_txt, GROUP_TABLE *G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv)
#endif
{
#ifdef GROUP_CANON_TG
	/* group_calc_cf_tg(fp_txt, G, cf, f_v, f_vv); */
	group_calc_cf_ab(fp_txt, G, cf, f_v, f_vv);
#else
	group_calc_cf_ab(fp_txt, G, cf, f_v, f_vv);
#endif
	return OK;
}

static INT total_canon_time = 0;

#ifdef GROUP_CANON_TG
#if TEXDOCU
INT group_calc_cf_tg(FILE *fp_txt, GROUP_TABLE *G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv)
#endif
{
	GROUP_CANON GC;
	ARRAY < VEKTOR < short > > tafel;
	INT n, i, j, g_i, g_j, g_i0, g_j0, aa, a, nb_gen;
	INT t0, t1, user_time;
	BYTE str[256];
	PERMUTATION_OB p, pv;
	VECTOR_OB go, g, ro, re;
	MATRIX_OB GG;

	n = G->n;
	cf->init();
	tafel.REALLOC(n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a = G->table->s_iji(i, j);
			tafel[i + 1][j + 1] = a + 1;
			}
		}
	if (f_v) {
		printf("group_calc_cf_tg(): calling GC.Init(), GC.Run()...\n");
		fflush(stdout);
		}
	
	t0 = os_ticks();
	GC.Init(n, &tafel);
	GC.Run();
	VEKTOR < short > & p_ = GC.GetCanNum();
	p.m_il(n);
	for (i = 0; i < n; i++) {
		a = p_[i + 1];
		p.m_ii(i, a);
		}
	p.copy(cf->s_p0v());
	p.invers(&pv);
	pv.copy(cf->s_p0());

	VEKTOR < short > & go_ = GC.GetGO();
	nb_gen = go_.Used();
	cf->s_nb_gen()->m_i(nb_gen);
	go.m_il(nb_gen + 1);
	for (i = 0; i < nb_gen; i++) {
		a = go_[i + 1];
		go.m_ii(i, a);
		}
	go.m_ii(nb_gen, n);
	go.copy(cf->s_go());
	
	g.m_il_n(nb_gen + 1);
	g.copy(cf->s_g()); /* not used ! */
	
	VEKTOR < short > & ro_ = GC.GetRelOrders();
	ro.m_il(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		a = ro_[i + 1];
		ro.m_ii(i, a);
		}
	ro.copy(cf->s_ro());
	
	VEKTOR < short > & re_ = GC.GetImages();
	re.m_il(nb_gen);
	for (i = 0; i < nb_gen; i++) {
		a = re_[i + 1] - 1; /* element no ab = element no tg - 1 !!! */
		re.m_ii(i, a);
		}
	re.copy(cf->s_re());
	
	GG.m_ilih(nb_gen, nb_gen);
	for (i = 0; i < nb_gen; i++) {
		g_i = go.s_ii(i);
		g_i0 = p.s_ii(g_i) - 1;
		for (j = 0; j < nb_gen; j++) {
			g_j = go.s_ii(j);
			g_j0 = p.s_ii(g_j) - 1;
			aa = G->table->s_iji(g_i0, g_j0);
			a = pv.s_ii(aa) - 1;
			GG.m_iji(i, j, a);
			}
		}
	GG.swap(cf->s_GG());

	cf->s_ago_trunc()->m_i(0);
	if (f_vv) {
		printf("group_calc_cf_tg():\n");
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
	
	t1 = os_ticks();
	user_time = t1 - t0;
	total_canon_time += user_time;
	if (f_v) {
		str[0] = 0;
		print_delta_time(user_time, str);
		printf("found canonic form %s\n", str);
		str[0] = 0;
		print_delta_time(total_canon_time, str);
		printf("total for canonic form: %s\n", str);
		fflush(stdout);
		}

#if 0
	if (group_canonicize(fp_txt, G, cf->s_p0(), cf->s_p0v(), 
		cf->s_go(), cf->s_g(), 
		cf->s_ro(), cf->s_re(), cf->s_GG(), 
		cf->s_auts(), cf->s_first_moved(), f_aut_gens, 
		&L, &ago, f_aut, f_v, f_vv) != OK)
		return f_error(fp_txt, "group_calc_cf() error in group_canonicize()");
#endif
	return OK;
}
#endif

#if TEXDOCU
INT group_calc_cf_ab(FILE *fp_txt, GROUP_TABLE *G, 
	GROUP_CANONIC_FORM_OP cf, INT f_v, INT f_vv)
#endif
{
	LABRA_OB L;
	SYM_OB ago;
	INT ago_trunc;
	INT f_aut = TRUE;
	INT f_aut_gens = TRUE;
	INT n, nb_gen;

	n = G->n;
	cf->init();
	if (group_canonicize(fp_txt, G, cf->s_p0(), cf->s_p0v(), 
		cf->s_go(), cf->s_g(), 
		cf->s_ro(), cf->s_re(), cf->s_GG(), 
		cf->s_auts(), cf->s_first_moved(), f_aut_gens, 
		&L, &ago, f_aut, f_v, FALSE) != OK)
		return f_error(fp_txt, "group_calc_cf() error in group_canonicize()");
	nb_gen = cf->s_ro()->s_li();
	cf->s_nb_gen()->m_i(nb_gen);
	if (f_vv) {
		printf("group_calc_cf():\n");
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
	return OK;
}

#if TEXDOCU
INT group_canonicize(FILE *fp_txt, GROUP_TABLE *G, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut, 
	INT f_v, INT f_vv)
#endif
{
	GT_CANON_INFO info;
	INT t0, t1, user_time;
	BYTE str[256];

	if (f_v) {
		/* print_size(); */
		fprintf(fp_txt, "(");
		fflush(fp_txt);
		}

	t0 = os_ticks();
	info.fp_txt = fp_txt;
	if (f_v) {
		/* print_size(); */
		fprintf(fp_txt, "info_fill_in_theG ");
		fflush(fp_txt);
		}
	info_fill_in_theG(&info, G);
	/* computes:
		info->n
		info->theG
		info->dim_n
		info->f_has_colors = FALSE, 
		info->colors = NIL
	*/
	if (f_v) {
		/* print_size(); */
		fprintf(fp_txt, "fill_in_colors ");
		fflush(fp_txt);
		}
	info_fill_in_colors(fp_txt, &info, G, f_v);
		/* computes the coloring and fills it into info */

	info.nb_isomorphisms = 0;
	/* info.mode = GT_CANON_MODE_BREAK_AFTER_FST; */
	/* info.f_autologisms = FALSE; */
	info.f_verbose = f_v;
	info.f_very_verbose = f_vv;
	if (f_vv) {
		printf("f_verbose = %ld f_very_verbose = %ld\n", f_v, f_vv);
		fflush(stdout);
		}

	/* info.add_aut = NIL;
	info.add_aut_perm = NIL;
	info.data = NIL; */
	
	if (!gt_canonicize_main(&info, p0, p0v, go, g, Ro, Re, GG, 
		auts, first_moved, f_aut_gens, Aut, ago, f_aut)) {
		return error("group_canonicize() error in gt_canonicize_main()");
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

l_exit:
	my_ptr_free( (void **) &info.theG);
	if (info.f_has_colors) {
		freeall(info.colors);
		info.colors = NIL;
		info.f_has_colors = FALSE;
		}
	return OK;
}

#if TEXDOCU
static INT gt_canonicize_main(GT_CANON_INFO *info, 
	PERMUTATION_OP p0, PERMUTATION_OP p0v, 
	VECTOR_OP go, VECTOR_OP g, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG, 
	VECTOR_OP auts, VECTOR_OP first_moved, INT f_aut_gens, 
	LABRA_OP Aut, SYM_OP ago, INT f_aut)
#endif
{
	GT_CANON_LOCAL *L = NIL;
	INT n, i, j, nb_gen;
	
	L = (GT_CANON_LOCAL *) my_malloc(sizeof(GT_CANON_LOCAL), 
		"gt_canonicize_main L");
	if (L == NIL)
		return error("gt_canonicize_main(): no memory for L");
	gt_canon_nil(L);
	L->info = info;
	n = L->info->n;
	L->n = n;
	the_group_int(L);
		/* allocates p, pv, q, qv, tmp1, tmp2 */
	/* L->A.fp_txt = info->fp_txt; */
	if (info->f_has_colors) {
		the_group_init_colors(L, info->colors);
		}
	
	L->theG = info->theG;
	L->dim_n = info->dim_n;
	L->f_going_back = FALSE;
	
	the_group_open_grid(L, 1 /* nb_grid */, L->n /* n */);
	if (info->f_very_verbose) {
		printf("gt_canonicize_main(): "
		"nach the_group_open_grid n = %ld\n", L->n);
		fflush(stdout);
		}

	L->E = (ORDERED_SET *) my_malloc(L->n * sizeof(ORDERED_SET), 
		"gt_canonicize_main E");
	if (L->E == NIL)
		return error("gt_canonicize_main(): no memory for L->E");


	L->is_first = TRUE;
	L->first_moved = AUTLOG_MAX_G;
	sp_nil(&L->p0);
	sp_nil(&L->p0v);
	sp_int(&L->p0, L->n, "gt_canonicize_main L->p0");
	sp_int(&L->p0v, L->n, "gt_canonicize_main L->p0v");
	sp_id(&L->p0);
	sp_id(&L->p0v);
	L->nb_aut_gens = 0;
	L->dim_aut_gens = 0;
	L->aut_gens = NIL;
	L->auts_first_moved = NIL;
	
	last_print_ticks = os_ticks();
	nb_backtrack = 0;
	nb_backtrack_points = 0;
	
	gt_canonicize(L, 0 /* level */);

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
		j = L->go0[i];
		go->m_ii(i, j);
		j = L->g0[i];
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
	L->theG = NIL;
	my_ptr_free( (void **) &L->E);
	the_group_ext(L);

	my_free(L);
	return TRUE;
}

#if TEXDOCU
static INT gt_canonicize(GT_CANON_LOCAL *L, INT level)
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
		printf("gt_canonicize() level = %ld\n", level);
		/* printf("f_verbose = %ld f_very_verbose = %ld\n", 
			L->info->f_verbose, L->info->f_very_verbose); */
		fflush(stdout);
		if (L->info->f_very_verbose >= 2) {
			group_print(L, &L->p, &L->pv);
			}
		}
	calc_grid(L, L->Grid, level, L->n /* go */, 
		L->info->f_very_verbose);
	E = L->E + level;
	choose_first_into_E(L, L->Grid, level, 
		E, L->info->f_very_verbose);
	len = E->size;
	for (i = 0; i < len; i++) {
		g = E->a[i];
		gt_add_generator(L, level, g, L->info->f_very_verbose);
		if (L->info->f_very_verbose) {
			if (L->info->f_very_verbose >= 2) {
				printf("gt_canonicize() level = %ld added generator %ld\n", level, g);
				fflush(stdout);
				printf("gt_canonicize() group order(%ld) = %ld\n", level + 1, 
					the_group_order(L, level + 1));
				fflush(stdout);
				}
			print_labelling_g(L, level + 1);
			if (L->info->f_very_verbose >= 2) {
				group_print(L, &L->p, &L->pv);
				}
			}
		if (i > 0 && level < L->first_moved)
			L->first_moved = level;
		if (L->is_first) {
			if (the_group_order(L, level + 1) == L->n) {
				L->is_first = FALSE;
				get_canonical_labelling(L, level + 1);
				}
			else {
				gt_canonicize(L, level + 1);
				}
			}
		else { /* (!L->is_first) */
			go_last = the_group_order(L, level);
			go = the_group_order(L, level + 1);
			if (go == L->n) {
				ret = canonicize_compare_tables(L, level + 1);
				if (ret == 0) {
					/* we have found an automorphism ! */
					if (L->info->f_very_verbose) {
						printf("found an automorphism ! in transversal %ld\n", 
							L->first_moved);
						fflush(stdout);
						print_labelling_g(L, level + 1);
						if (L->info->f_very_verbose >= 1) {
							group_print(L, &L->p, &L->pv);
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
					L->first_moved = GT_CANON_MAX_GEN;
					free_aut_gen(L);
					}
				}
			else {
				ret = canonicize_compare_tables(L, level + 1);
				if (ret >= 0) {
					if (gt_canonicize(L, level + 1) == 1 && level > L->first_moved)
						return 1;
					}
				}
			}
		} /* next i */
	return 0;
}

#if TEXDOCU
static INT recalc_Ro_Re_GG(GT_CANON_LOCAL *L, 
	VECTOR_OP Ro, VECTOR_OP Re, MATRIX_OP GG)
#endif
{
	INT nb_gen, i, j, g_i, g_j, k, k_, n;
	INT go_last, go, g, ro, re, re_;
	INT *theG, dim_n;

	n = L->n;
	nb_gen = L->nb_gen;
	theG = L->theG;
	dim_n = L->dim_n;
	Ro->m_il(nb_gen);
	Re->m_il(nb_gen);
	GG->m_ilih(nb_gen, nb_gen);
	for (i = 0; i < nb_gen; i++) {
		g = L->g0[i + 1];
		go_last = L->go0[i];
		if (L->info->f_very_verbose) {
			printf("recalc_Ro_Re_G(): "
				"i = %ld go_last = %ld g = %ld\n", 
				i, go_last, g);
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
			printf("ro = %ld re = %ld re_ = %ld\n", 
				ro, re, re_); 
			fflush(stdout);
			}
		}
	for (i = 0; i < nb_gen; i++) {
		g_i = L->g0[i + 1];
		for (j = 0; j < nb_gen; j++) {
			g_j = L->g0[j + 1];
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

#define REALLOC_OVERSIZE 4

#if TEXDOCU
static INT get_aut_gen(GT_CANON_LOCAL *L)
#endif
{
	INT new_dim, i, n;
	SPERM *new_auts, *p;
	INT *new_first;

	if (L->nb_aut_gens >= L->dim_aut_gens) {
		/* realloc: */
		new_dim = L->dim_aut_gens + REALLOC_OVERSIZE;
		new_auts = (SPERM *) my_malloc(new_dim * sizeof(SPERM), 
			"gt_canon.C: get_aut_gen()");
		if (new_auts == NIL)
			return error("get_aut_gen(): no memory for new_auts");
		new_first = (INT *) my_malloc(new_dim * sizeof(INT), 
			"gt_canon.C: get_aut_gen()");
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
	n = L->n;
	sp_nil(p);
	sp_int(p, n, "gt_canon.C: get_aut_gen()");
	sp_mult(&L->p0, &L->pv, p);


	if (L->info->f_very_verbose) {
		printf("pv=");
		sp_print(&L->pv);
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
static void free_aut_gen(GT_CANON_LOCAL *L)
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
static INT get_canonical_labelling(GT_CANON_LOCAL *L, INT level)
#endif
{
	INT go, go_, ii;
	
	sp_mv(&L->p, &L->p0);
	sp_mv(&L->pv, &L->p0v);
	L->nb_gen = level;
	go_ = 0;
	for (ii = 0; ii <= L->nb_gen; ii++) {
		go = the_group_order(L, ii);
		L->go0[ii] = go;
		L->g0[ii] = L->pv.a[go_];
		go_ = go;
		}
	if (L->info->f_very_verbose) {
		printf("canonicize() canonical labelling:\n");
		printf("generators: ");
		print_labelling_g(L, level);
		printf("group orders: ");
		print_labelling_go(L, level);
		if (L->info->f_very_verbose >= 1) {
			group_print(L, &L->p, &L->pv);
			}
		}
	return OK;
}

#if TEXDOCU
static void print_labelling_g(GT_CANON_LOCAL *L, INT level)
#endif
{
	INT go, go_, ii, g;
	
	printf("{ ");
	go_ = 0;
	for (ii = 0; ii <= level; ii++) {
		go = the_group_order(L, ii);
		g = L->pv.a[go_];
		printf("%ld ", g);
		go_ = go;
		}
	printf("}\n");
	fflush(stdout);
}

#if TEXDOCU
static void print_labelling_go(GT_CANON_LOCAL *L, INT level)
#endif
{
	INT go, ii;
	
	printf("go: { ");
	for (ii = 0; ii <= level; ii++) {
		go = the_group_order(L, ii);
		printf("%ld ", go);
		}
	printf("}\n");
	fflush(stdout);
}

#if TEXDOCU
static INT canonicize_compare_tables(GT_CANON_LOCAL *L, INT level)
#endif
{
	INT i1, i0, i0_, j1, j0, j0_, k0, k0_, k1, k1_, i1_save;
	INT *theG, dim_n;
	INT ret;
	SPERM *p, *pv, *p_, *pv_;
	INT ii, go, go_;
	
	ret = 0;
	for (ii = 0; ii <= level; ii++) {
		go = the_group_order(L, ii);
		go_ = L->go0[ii];
		if (go > go_) {
			ret = 1;
			goto l_exit;
			}
		if (go < go_) {
			ret = -1;
			goto l_exit;
			}
		}
	go = the_group_order(L, level);
	theG = L->theG;
	dim_n = L->dim_n;
	p = &L->p;
	pv = &L->pv;
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
				return error("canonicize_compare_tables() "
					"k1 >= go");
				}
			if (k1_ >= go) {
				return error("canonicize_compare_tables() "
					"k1_ >= go");
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
				return error("canonicize_compare_tables() "
					"k1 >= go");
				}
			if (k1_ >= go) {
				return error("canonicize_compare_tables() "
					"k1_ >= go");
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

l_exit:
	if (L->info->f_very_verbose) {
		printf("compare tables: level = %ld ret = %ld\n", level, ret);
		}
	return ret;
}

#if TEXDOCU
static void print_E_colors(GT_CANON_LOCAL *L, INT level)
#endif
{
	ORDERED_SET *E;
	INT i, l, len, a, b;
	FILE *fp_txt;
	
	fp_txt = L->info->fp_txt;
	/* fprintf(fp_txt, "print_E_colors() level = %ld\n", level); */
	if (!L->f_has_element_colors) {
		fprintf(fp_txt, "no colors !\n");
		return;
		}
	for (l = 0; l < level; l++) {
		E = L->E + l;
		len = E->size;
		fprintf(fp_txt, "E[%ld]={ ", l);
		for (i = 0; i < len; i++) {
			a = E->a[i];
			b = L->element_colors[a];
			fprintf(fp_txt, "%ld ", b);
			}
		fprintf(fp_txt, "} %ld\n", len);
		}
	fflush(fp_txt);
}

#if TEXDOCU
static void print_element_colors(GT_CANON_LOCAL *L, INT level)
#endif
{
	FILE *fp_txt;
	INT i, ii, go, go_last = 0, i0, b;
	
	fp_txt = L->info->fp_txt;
	if (!L->f_has_element_colors)
		return;
	for (ii = 0; ii <= level; ii++) {
		go = the_group_order(L, ii);
		fprintf(fp_txt, "level %ld go = %ld elements: ", ii, go);
		for (i = go_last; i < go; i++) {
			i0 = L->pv.a[i];
			b = L->element_colors[i0];
			fprintf(fp_txt, "%ld ", b);
			}
		fprintf(fp_txt, "\n");
		go_last = go;
		}
	fflush(fp_txt);
}

#if TEXDOCU
static INT get_aut_gens(GT_CANON_LOCAL *L, 
	VECTOR_OP auts, VECTOR_OP first_moved)
#endif
{
	SPERM *p, *q;
	PERMUTATION_OP perm;
	INT i, n, nb_gen, ii, k;

	n = L->n;
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
static INT get_aut_group(GT_CANON_LOCAL *L, 
	LABRA_OP Aut, SYM_OP ago)
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
	n = L->n;
	Aut->Init_no_generators(n, f_v, f_vv);
	nb_gen = L->nb_aut_gens;
	sp_nil(&tmp1);
	sp_nil(&tmp2);
	sp_int(&tmp1, n, "gt_canon.C: get_aut_group");
	sp_int(&tmp2, n, "gt_canon.C: get_aut_group");
	for (i = 0; i < nb_gen; i++) {
		p = L->aut_gens + i;
		q = p;
		sp_mult(&L->p0v, p, &tmp1);
		sp_mult(&tmp1, &L->p0, &tmp2);
		q = &tmp2;
		fst_moved = L->auts_first_moved[i];
		j = L->go0[fst_moved];
			/* the automorphism moves j_ */
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
static INT info_fill_in_theG(GT_CANON_INFO *info, GROUP_TABLE *G)
#endif
{
	MATRIX_OP table;
	INT *theG = NIL;
	INT n, i, j, size, dim_n1;
	
	/* fill in the group table: */
	n = G->n;
	dim_n1 = n + 1;
	size = n * dim_n1 * sizeof(INT);
	theG = (INT *) my_malloc(size, "info_fill_in_theG");
	if (theG == NIL)
		return error("info_fill_in_theG(): no memory");
	table = G->table;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			theG[i * dim_n1 + j] = table->s_iji(i, j);
			}
		}
	for (i = 0; i < n; i++) {
		theG[i * dim_n1 + n] = G->table_inv->s_ii(i);
		}
	info->n = n;
	info->theG = theG;
	info->dim_n = dim_n1;
	info->f_has_colors = FALSE;
	info->colors = NIL;
	return OK;
}

#if TEXDOCU
static INT info_fill_in_colors(FILE *fp, GT_CANON_INFO *info, 
	GROUP_TABLE *G, INT f_v)
#endif
{
	VECTOR_OP col = (VECTOR_OP) callocobject("gt_canon.C: info_fill_in_colors");
	VECTOR_OB val, mult;
	BYTE str[10000];
		
	group_Coloring_seldoms_prefered(G, col, FALSE);
	col->multiplicities(&val, &mult);
	if (f_v) {
		str[0] = 0;
		val.sprint_multiplicities(&mult, str);
		fprintf(fp, "coloring: %s\n", str);
		fflush(fp);
		}
	
	info->f_has_colors = TRUE;
	info->colors = col;
	return OK;
}

#if TEXDOCU
static void the_group_nil(GT_CANON_LOCAL *L)
#endif
{
	L->theG = NIL;
	L->Grid = NIL;
	sp_nil(&L->p);
	sp_nil(&L->pv);
	sp_nil(&L->q);
	sp_nil(&L->qv);
	sp_nil(&L->r);
	sp_nil(&L->rv);
	L->nb_grid = 0;
	L->f_has_element_colors = FALSE;
	L->element_colors = NIL;
}

#if TEXDOCU
static void the_group_int(GT_CANON_LOCAL *L)
#endif
{
	sp_int(&L->p, L->n, "the_group_int p");
	sp_int(&L->pv, L->n, "the_group_int pv");
	sp_int(&L->q, L->n, "the_group_int q");
	sp_int(&L->qv, L->n, "the_group_int qv");
	sp_int(&L->r, L->n, "the_group_int r");
	sp_int(&L->rv, L->n, "the_group_int rv");
	sp_id(&L->p);
	sp_id(&L->pv);
	sp_id(&L->q);
	sp_id(&L->qv);
	sp_id(&L->r);
	sp_id(&L->rv);
	L->f_has_element_colors = FALSE;
}

#if TEXDOCU
static INT the_group_init_colors(GT_CANON_LOCAL *L, VECTOR_OP col)
#endif
{
	INT n, i, a;
	
	L->f_has_element_colors = TRUE;
	n = L->n;
	if (n != col->s_li())
		return error("the_group_init_colors(): n != col->s_li()");
	L->element_colors = (INT *) my_malloc(n * sizeof(INT), 
		"the_group_init_colors");
	for (i = 0; i < n; i++) {
		a = col->s_ii(i);
		L->element_colors[i] = a;
		}
	return OK;
}

#if TEXDOCU
static INT the_group_open_grid(GT_CANON_LOCAL *L, INT nb_grid, INT n)
#endif
{
	INT i;
	
	if (L->nb_grid)
		return error("the_group_open_grid(): L->nb_grid");
	L->Grid = (GT_CANON_GRID *) my_malloc(nb_grid * sizeof(GT_CANON_GRID), 
		"the_group_open_grid");
	if (L->Grid == NIL)
		return error("the_group_open_grid(): no memory for L->Grid");
	L->nb_grid = nb_grid;
	for (i = 0; i < L->nb_grid; i++) {
		grid_nil(L->Grid + i);
		grid_int(L->Grid + i, n);
		}
	return OK;
}

#if TEXDOCU
static void the_group_ext(GT_CANON_LOCAL *L)
#endif
{
	INT i;
	
	my_ptr_free( (void **) &L->theG);
	if (L->nb_grid) {
		for (i = 0; i < L->nb_grid; i++) {
			grid_ext(L->Grid + i);
			}
		my_ptr_free( (void **) &L->Grid);
		L->nb_grid = 0;
		}
	sp_free(&L->p);
	sp_free(&L->pv);
	sp_free(&L->q);
	sp_free(&L->qv);
	sp_free(&L->r);
	sp_free(&L->rv);
	if (L->f_has_element_colors) {
		my_free(L->element_colors);
		L->element_colors = NIL;
		L->f_has_element_colors = FALSE;
		}
}

#if TEXDOCU
static INT the_group_print_generators(GT_CANON_LOCAL *L)
#endif
{
	INT i;
	
	for (i = 0; i < L->nb_gen; i++)
		printf("%ld ", L->g[i]);
	printf("\n");
	return OK;
}

#if TEXDOCU
static INT the_group_order(GT_CANON_LOCAL *L, INT k)
#else
/* k = 0: 1 fuer Einsgruppe, sonst go[k - 1] */
#endif
{
	if (k < 0)
		return error("the_group_order k < 0");
	if (k == 0)
		return 1;
	return L->go[k - 1];
}

#if TEXDOCU
static void grid_nil(GT_CANON_GRID *p)
#endif
{
	sp_nil(&p->r);
	sp_nil(&p->rv);
}

#if TEXDOCU
static void grid_int(GT_CANON_GRID *p, INT n)
#endif
{
	p->m = n;
	sp_int(&p->r, p->m, "grid_int");
	sp_int(&p->rv, p->m, "grid_int");
	sp_id(&p->r);
	sp_id(&p->rv);
}

#if TEXDOCU
static void grid_ext(GT_CANON_GRID *p)
#endif
{
	sp_free(&p->r);
	sp_free(&p->rv);
}

#if TEXDOCU
static void gt_canon_nil(GT_CANON_LOCAL *L)
#endif
{
	L->info = NIL;
	L->E = NIL;
	the_group_nil(L);
}

#if TEXDOCU
static void group_print(GT_CANON_LOCAL *G, SPERM *p, SPERM *pv)
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

// #include "gt_canon_grid.C"
/* gt_canon_grid.C 
 * included from gt_canon.C */

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
static INT calc_grid(GT_CANON_LOCAL *L, 
	GT_CANON_GRID *G, INT k, INT go, INT f_verbose)
#else
Relativordnungen bezueglich 
$E = U_{0}, U_1, ... U_{k}$
nach G berechnen
(k == 0 heisst nur bezueglich Einsgruppe). 
$U_k = <U_0, ... U_{k - 1}>$.
Es wird die Tabelle fuer alle Elemente 
zwischen go\_k und go aufgestellt.
Diese muessen mittels C-$>$p 
an den Beginn der Gruppentafel 
permutiert sein.
#endif
{
	INT len, k1, i, j, i1, j1;
	INT ro, re, o, o1, t_i, n0 = 0;
	INT dim_n, *theG;
	INT go_k, go_k1;
	SPERM *p, *pv;
	
	dim_n = L->dim_n;
	theG = L->theG;
	p = &L->p;
	pv = &L->pv;
	
	if (G->m != L->n)
		return error("calc_grid(): G->m != L->n");
	sp_id(&G->r);
	sp_id(&G->rv);
	if (L->f_has_element_colors)
		n0 = 1;
	G->n = n0 + 1 + k * 2; /* (n0 for element coloring) */
	go_k = the_group_order(L, k);
	/* go_k war vorher a->n1[k]. */
	len = go - go_k;
	for (i = 0; i < len; i++)
		for (j = 0; j < G->n; j++)
			G->type[i][j] = 0;
	
	for (i = 0; i < go_k; i++)
		G->type_idx[i] = - 1;
	for ( ; i < go; i++)
		G->type_idx[i] = i - go_k;
	
	for (i = go_k; i < go; i++) {
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
		if (L->f_has_element_colors)
			G->type[t_i][0] = - L->element_colors[i1];
		
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
			go_k1 = the_group_order(L, k1);
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
				re < the_group_order(L, k1 - 1)) {
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
		printf("in calc_grid() (k = %ld):\n", k);
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
		printf("vor radix_sort():\n");
		print_grid(G, pv, k, go_k, go);
		}
#endif
#endif

	radix_sort(G, 0 /* radix */, go_k, go - 1);
	
#ifdef DEBUG_CALC_GRID
	if (f_verbose) {
		printf("nach radix_sort():\n");
		print_grid(G, pv, k, go_k, go);
		}
#endif

	return OK;
}

#if TEXDOCU
static INT radix_sort(GT_CANON_GRID *G, INT radix, INT first, INT last)
#else
Die Typvektoren werden aufsteigend sortiert.
Most significant ist dabei die 0-te Stelle von type (radix = 0).
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
			return error("radix_sort(): G->first[k] != first");
		for (l = first; l <= last; l++)
			G->grid_entry[l] = k;
		G->len[k] = last - first + 1;
		G->first[k + 1] = last + 1;
		G->G_max++;
		/* if (a->info->f_very_verbose && 
			a->info->f_verbose) {
			printf("radix_sort()|new entry: first = %ld len = %ld : ", 
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
		f_found = insert_idx(G, first, k - first, radix, k, &idx);
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
			return error("radix_sort(): descending");
		if (res > 0) { /* x > y */
			radix_sort(G, radix + 1, first0, k - 1);
			first0 = k;
			first1 = G->type_idx[first0];
			}
		if (k == last) {
			radix_sort(G, radix + 1, first0, k);
			}
		}
	return OK;
}

#if TEXDOCU
static INT insert_idx(GT_CANON_GRID *G, 
	INT first, INT len, INT radix, 
	INT search_this, INT *idx)
#else
aufsteigende Sortierung, 
bestimme Einfuegeposition fuer "search\_this" Element 
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
static void print_grid(GT_CANON_GRID *G, SPERM *pv, 
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
		print_type_vector(G, first, k);
		}
}

#if TEXDOCU
static void print_type_vector(GT_CANON_GRID *G, INT i, INT k)
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
static INT calc_n_grids(GT_CANON_LOCAL *L, 
	INT go, INT n, INT f_verbose)
#endif
{
	INT k;
	
	if (L->nb_grid < n)
		return error("calc_n_grids(): "
		"L->nb_grid < n");
	for (k = 0; k < n; k++)
		calc_grid(L, L->Grid + k, 
			k, go, f_verbose);
	return OK;
}

#if TEXDOCU
static INT choose_first_into_E(GT_CANON_LOCAL *L, 
	GT_CANON_GRID *Grid, INT k, ORDERED_SET *E, INT f_verbose)
#endif
{
	INT i, j, l, first, len;
	INT g0, g1, g2, ge, go, i0, i1, i2;
		
	go = the_group_order(L, k); /* before: n1 */
	if (go != Grid->first[0])
		return error("choose_first_into_E(): go != Grid->first[0]");
	ge = Grid->G_max - 1;
	first = Grid->first[ge];
	len = Grid->len[ge];
	E->size = 0;
	for (i = 0; i < len; i++) {
		i2 = first + i;
		i1 = Grid->rv.a[i2];
		i0 = L->pv.a[i1];
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
		print_E(E, k);
		}
	return OK;
}

#if TEXDOCU
static void print_E(ORDERED_SET *E, INT k)
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


// #include "gt_canon_dimino.C"
/* gt_canon_dimino.C 
 * included from gt_canon.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */

#if TEXDOCU
static INT gt_add_generator(GT_CANON_LOCAL *B, 
	INT i, INT g, INT f_verbose)
#endif
{
	INT g1, go, last_go;
	
	if (f_verbose) {
		printf("gt_add_generator(): "
		"i = %ld generator %ld\n", i, g);
		}
	g1 = B->p.a[g];
	last_go = the_group_order(B, i);
	if (g1 < last_go) {
		printf("i = %ld g = %ld g1 = %ld "
			"last_go = %ld B->nb_gen = %ld "
			"|B| = %ld\n", 
			i, g, g1, last_go, B->nb_gen, 
			the_group_order(B, B->nb_gen) );
		printf("error in gt_add_generator(): "
		"g1 < last_go"); fflush(stdout);
		return ERROR;
		}
	sp_id(&B->q);
	sp_id(&B->qv);

	gt_dimino(B, i, g, &go);
	B->g[i] = g;
	B->go[i] = go;

	sp_mult(&B->p, &B->q, &B->r);
	sp_mult(&B->qv, &B->pv, &B->rv);
	sp_mv(&B->r, &B->p);
	sp_mv(&B->rv, &B->pv);
	sp_id(&B->q);
	sp_id(&B->qv);
	
	return OK;
}

#if TEXDOCU
static INT gt_dimino(GT_CANON_LOCAL *C, 
	INT k, INT g0, INT *go)
#else
/* Rueckgabe FALSE nur wenn f\_test\_it TRUE. */
#endif
{
	INT g1, coset_size, rep_pos;
	INT go1, h, s0, s1, s2, l;
	
	coset_size = the_group_order(C, k);
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
		return error("gt_dimino(): C->q.a[g1] != g1");
		}
	if (!gt_add_coset(C, coset_size, g1, &go1)) {
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
			h = gt_mult(C, rep_pos, s2);
			/* printf("%ld * %ld = %ld\n", 
				rep_pos, s2, h); */
			if (h >= go1) {
				if (!gt_add_coset(C, 
					coset_size, h, &go1)) {
					return FALSE;
					}
				}
			} /* for l */
		rep_pos += coset_size;
		}
	*go = go1;
	return TRUE;
}

#if TEXDOCU
static INT gt_add_coset(GT_CANON_LOCAL *C, 
	INT coset_size, INT g2, INT *go)
#else
Annahme, dass die grids auf die Elemente 
in der Reihenfolge vor Anwendung von 
(C-$>$q sowie A-$>$q)
berechnet worden sind 
(nur falls f\_test\_it TRUE). 
Die erfolgten Permutationen 
werden nach q akkumuliert. 
Rueckgabe FALSE kann nur bei 
f\_test\_it TRUE auftreten. 
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
		return error("gt_add_coset(): g2 < *go");
		}
	for (i = 0; i < coset_size; i++) {
		j = gt_mult(C, i, g2);
		if (j < *go) {
			printf("i = %ld g2 = %ld g1 = %ld "
				"g0 = %ld qv->a[g2] = %ld\n", 
				i, g2, g1, g0, (INT) C->qv.a[g2]);
			printf("j = %ld *go = %ld\n", j, *go);
			if (C->info->fp_txt) {
				fprintf(C->info->fp_txt, "gt_add_coset(): j < *go"
					"i = %ld g2 = %ld g1 = %ld "
					"g0 = %ld qv->a[g2] = %ld\n", 
					i, g2, g1, g0, (INT) C->qv.a[g2]);
				}
			return error("gt_add_coset(): j < *go");
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
static INT gt_mult(GT_CANON_LOCAL *A, INT i, INT j)
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


#endif /* SOLVABLE_TRUE */


