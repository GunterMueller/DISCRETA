/* dcc_orderly.C */

#include <DISCRETA/discreta.h>


#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/cp.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/fga.h>
#include <DISCRETA/ladder.h>

#include <DISCRETA/geo.h>

INT calc_kramer_mesner_matrix_by_orderly_generation(
	VECTOR_OP gen, SYM_OP go, INT deg, 
	LABRA_OP labG, MATRIX_OP TG, 
	BYTE *g_label, BYTE *g_label_tex, 
	INT t, INT k, INT f_TDO)
{
	LADDER_INFO *li;
	INT n;
	MATRIX_OB Mtk, Mtk2;
	VECTOR_OB Orbits, Stabs;
	
	li = init_ladder_info(gen, go, deg, g_label, g_label_tex, 
		t, k, 0 /* lambda */, 
		TRUE /* f_verbose */, 
		NIL );
	
	li_init_file_names(li);
		/* compute filename of kramer mesner matrix file into li->txt_out 
		 * and writes this file into "group_label" */

	li_message(li);
		/* hi I am DCC */
	
	n = li->deg;
	
	printf("computing KM-matrix for t=%ld k=%ld:\n", t, k);
	fflush(stdout);
	
	dcc_orderly(labG, TG, n, k, &Orbits, &Stabs);
	
#if 0
	calc_KM_tk(labG, TG, 
		(VECTOR_OP) Orbits.s_i(t), 
		(VECTOR_OP) Orbits.s_i(k), 
		&Mtk, n, t, k);
	Mtk.Print();
	printf("writing KM-matrix:\n");
	fflush(stdout);
	KM_tk_print_asc(labG, g_label, g_label_tex, &Mtk, t, k);
#endif
	
#if 1
	calc_KM_tk_column_wise(labG, TG, 
		(VECTOR_OP) Orbits.s_i(t), 
		(VECTOR_OP) Stabs.s_i(t), 
		(VECTOR_OP) Orbits.s_i(k), 
		(VECTOR_OP) Stabs.s_i(k), &Mtk, n, t, k);
	Mtk.Print();
	printf("writing KM-matrix:\n");
	fflush(stdout);
	KM_tk_print_asc(labG, g_label, g_label_tex, &Mtk, t, k);
#endif

	if (f_TDO) {
		km_compute_TDO_decomposition(g_label, t, k, &Mtk);
		}
	
	return OK;
}

INT dcc_orderly(LABRA_OP labG, MATRIX_OP TG, INT n, INT k, 
	VECTOR_OP Orbits, VECTOR_OP Stabs)
{
	VECTOR_OP k_orbits, k_stabs;
	VECTOR_OP km1_orbits, km1_stabs;
	VECTOR_OB val, mult;
	BYTE str[10000];
	INT kk;
	
	Orbits->m_il(k + 1);
	Stabs->m_il(k + 1);
	for (kk = 0; kk <= k; kk++) {
		VECTOR_OB k_stab_order;
		
		k_orbits = (VECTOR_OP) Orbits->s_i(kk);
		k_stabs = (VECTOR_OP) Stabs->s_i(kk);
		if (kk > 0) {
			km1_orbits = (VECTOR_OP) Orbits->s_i(kk - 1);
			km1_stabs = (VECTOR_OP) Stabs->s_i(kk - 1);
			
			extend_to_k_sets(labG, TG, n, kk, 
				km1_orbits, km1_stabs, k_orbits, k_stabs, TRUE, FALSE);
			}
		else {
			calc_k_sets(labG, TG, n, kk, k_orbits, k_stabs, TRUE, FALSE);
			}
	
		calc_k_stab_orders(k_stabs, &k_stab_order);
		
		k_stab_order.multiplicities(&val, &mult);
		str[0] = 0;
		val.sprint_multiplicities(&mult, str);
		printf("stabilizer orders: %s\n", str);
		}
	return OK;
}

INT calc_KM_column(LABRA_OP G, MATRIX_OP TG, 
	VECTOR_OP Reps_t, VECTOR_OP Stabs_t, VECTOR_OP k_set, LABRA_OP k_stab, 
	MATRIX_OP M, INT n, INT t, INT k, INT k_idx)
{
	INT X[10000];
	INT Y[10000];
	INT ZZ[10000];
	INT i, ii, iii, j, kmt, nt, a, idx, f_found;
	VECTOR_OB t_orbit;
	SYM_OB go, aa, aa1, bb, cc, cc1;
	LABRA_OP t_stab;
	PERMUTATION_OB transporter;

	G->group_order(&go);
	nt = Reps_t->s_li();
	if (nt != M->s_hi())
		return error("calc_KM_column(): nt != M->s_hi()");
	for (i = 0; i < nt; i++) 
		M->m_iji(i, k_idx, 0);
	kmt = k - t;
	for (i = 0; i < k; i++) {
		j = k_set->s_ii(i);
		X[i] = j;
		}
	n_choose_k_first(Y, k, kmt);
	while (TRUE) {
		ii = 0;
		iii = 0;
		for (i = 0; i < k; i++) {
			if (iii < kmt && i == Y[iii]) {
				iii++;
				}
			else {
				ZZ[ii] = X[i];
				ii++;
				}
			}
		if (iii != kmt)
			return error("calc_KM_column() iii != kmt");
		if (ii != t)
			return error("calc_KM_column() ii != t");
		
		geo_canonicize_set(n, t, ZZ, G, TG, FALSE, FALSE, 
			FALSE /* f_get_aut_group */, NIL, NIL, &transporter);
		t_orbit.m_il(t);
		for (i = 0; i < t; i++) {
			a = ZZ[i];
			t_orbit.m_ii(i, a);
			}
		Reps_t->search(nt, TRUE, &t_orbit, &idx, &f_found);
		if (!f_found) {
			return error("calc_KM_tk() k_orbit not found");
			}
		else {
			idx--;
			M->s_ij(idx, k_idx)->inc();
			}
		
		if (!n_choose_k_next(Y, k, kmt))
			break;
		}
	
	k_stab->group_order(&aa);
  	go.ganzdiv(&aa, &aa1);
	// aa1.println();
	for (i = 0; i < nt; i++) {
		// printf("M_%ld,%ld=", i, k_idx);
		// M->s_ij(i, k_idx)->println();
		M->s_ij(i, k_idx)->mult(&aa1, &bb);
		// bb.println();
		t_stab = (LABRA_OP) Stabs_t->s_i(i);
		t_stab->group_order(&cc);
  		go.ganzdiv(&cc, &cc1);
		// cc1.println();
  	 	bb.ganzdiv(&cc1, M->s_ij(i, k_idx));
		// printf("M_%ld,%ld=", i, k_idx);
		// M->s_ij(i, k_idx)->println();
		// printf("\n");
		}
	return OK;
}

INT calc_KM_tk_column_wise(LABRA_OP G, MATRIX_OP TG, 
	VECTOR_OP Reps_t, VECTOR_OP Stabs_t, 
	VECTOR_OP Reps_k, VECTOR_OP Stabs_k, 
	MATRIX_OP M, INT n, INT t, INT k)
{
	INT nt, nk, k_idx;

	nt = Reps_t->s_li();
	nk = Reps_k->s_li();
	M->m_ilih_n(nk, nt);
	printf("("); fflush(stdout);
	for (k_idx = 0; k_idx < nk; k_idx++) {
		calc_KM_column(G, TG, 
			Reps_t, Stabs_t, 
			(VECTOR_OP) Reps_k->s_i(k_idx), 
			(LABRA_OP) Stabs_k->s_i(k_idx), 
			M, n, t, k, k_idx);
		if (k_idx + 1 % 10 == 0)
			if (k_idx + 1 % 50 == 0)
				if (k_idx + 1 % 100 == 0)
					printf(".\n");
				else
					printf(":");
			else
				printf(",");
		else
			printf(".");
		fflush(stdout);
		}
	printf(")\n"); fflush(stdout);
	return OK;
}

INT calc_KM_tk(LABRA_OP G, MATRIX_OP TG, 
	VECTOR_OP Reps_t, VECTOR_OP Reps_k, MATRIX_OP M, INT n, INT t, INT k)
{
	INT nt, nk, i, ii, iii, kmt, a, b, idx, f_found;
	VECTOR_OP t_orbit;
	VECTOR_OB k_orbit;
	INT X[10000];
	INT Y[10000];
	INT Z[10000];
	INT ZZ[10000];
	INT free_positions[10000];
	PERMUTATION_OB transporter;

	kmt = k - t;
	nt = Reps_t->s_li();
	nk = Reps_k->s_li();
	M->m_ilih_n(nk, nt);
	for (i = 0; i < nt; i++) {
		t_orbit = (VECTOR_OP) Reps_t->s_i(i);
		for (ii = 0; ii < n; ii++)
			X[ii] = 0;
		for (ii = 0; ii < t; ii++) {
			iii = t_orbit->s_ii(ii);
			X[iii] = 1;
			}
		iii = 0;
		for (ii = 0; ii < n; ii++) {
			if (X[ii] == 0) {
				free_positions[iii++] = ii;
				}
			}
		n_choose_k_first(Y, n - t, kmt);
		while (TRUE) {
			for (ii = 0; ii < n; ii++) {
				Z[ii] = X[ii];
				}
			for (ii = 0; ii < kmt; ii++) {
				a = Y[ii];
				b = free_positions[a];
				Z[b] = 1;
				}
			iii = 0;
			for (ii = 0; ii < n; ii++) {
				if (Z[ii]) {
					ZZ[iii++] = ii;
					}
				}
			if (iii != k)
				return error("calc_KM_tk() iii != k");
			
			geo_canonicize_set(n, k, ZZ, G, TG, FALSE, FALSE, 
				FALSE /* f_calc_aut_group */, NIL, NIL, &transporter);
			k_orbit.m_il(k);
			for (ii = 0; ii < k; ii++) {
				a = ZZ[ii];
				k_orbit.m_ii(ii, a);
				}
			Reps_k->search(nk, TRUE, &k_orbit, &idx, &f_found);
			if (!f_found) {
				return error("calc_KM_tk() k_orbit not found");
				}
			else {
				idx--;
				M->s_ij(i, idx)->inc();
				}
			
			if (!n_choose_k_next(Y, n - t, kmt))
				break;
			}
		
		}
	return OK;
}

INT extend_to_k_sets(LABRA_OP G, MATRIX_OP TG, INT n, INT k, 
	VECTOR_OP km1_sets, VECTOR_OP km1_stabs, 
	VECTOR_OP k_sets, VECTOR_OP k_stabs, 
	INT f_v, INT f_vv)
{
	INT i, j, l, ii, iii, a, idx, f_found, res;
	VECTOR_OP km1_set;
	VECTOR_OB k_orbit;
	INT X[10000];
	INT free_positions[10000];
	INT ZZ[10000];
	LABRA_OB aut;
	SYM_OB ago;
	INTEGER_OB k_len;
	PERMUTATION_OB transporter;

	k_len.m_i(0);
	k_sets->m_il(0);
	k_stabs->m_il(0);
	l = km1_sets->s_li();
	for (i = 0; i < l; i++) {
		km1_set = (VECTOR_OP) km1_sets->s_i(i);
		for (ii = 0; ii < n; ii++) {
			X[ii] = 0;
			}
		for (ii = 0; ii < k - 1; ii++) {
			j = km1_set->s_ii(ii);
			X[j] = 1;
			}
		iii = 0;
		for (ii = 0; ii < n; ii++) {
			if (X[ii] == 0) {
				free_positions[iii++] = ii;
				}
			}
		if (iii != n - k + 1)
			return error("extend_to_k_sets(): iii != n - k + 1");
		for (j = 0; j < n - k + 1; j++) {
			if (f_vv) {
				printf("."); fflush(stdout);
				}
			X[free_positions[j]] = 1;
			
			iii = 0;
			for (ii = 0; ii < n; ii++) {
				if (X[ii]) {
					ZZ[iii++] = ii;
					}
				}
			if (iii != k)
				return error("extend_to_k_sets() iii != k");
			
			// printf("canonicize("); fflush(stdout);
			geo_canonicize_set(n, k, ZZ, G, TG, FALSE, FALSE, 
				FALSE /* f_get_aut_group */, NIL, NIL, &transporter);
			// printf(")\ngeo_maxtest_set("); fflush(stdout);
			res = geo_maxtest_set(n, k, ZZ, G, TG, &aut, &ago, FALSE, FALSE);
			// printf(")\n"); fflush(stdout);
			if (res != -1)
				return error("extens_to_k_sets(): canonical set not maximal !");
			k_orbit.m_il(k);
			for (ii = 0; ii < k; ii++) {
				a = ZZ[ii];
				k_orbit.m_ii(ii, a);
				}
			k_sets->search(k_len.s_i(), TRUE, &k_orbit, &idx, &f_found);
			if (!f_found) {
				if (f_vv) {
					printf("new orbit (extension of %ld by %ld): ", i, j);
					// print_set(Y, n, k);
					k_orbit.println();
					}
				k_sets->insert_at(&k_len, idx, &k_orbit);
				k_len.dec();
				k_stabs->insert_at(&k_len, idx, &aut);
				}
			
			X[free_positions[j]] = 0;
			
			} // next j
		if (f_vv) {
			printf("\n"); fflush(stdout);
			}
		} // next i
	k_sets->realloc(&k_len);
	k_stabs->realloc(&k_len);
	if (f_v) {
		printf("found %ld orbits on %ld sets\n", k_sets->s_li(), k);
		}
	// Reps->println();
	return OK;
}

INT calc_k_sets(LABRA_OP G, MATRIX_OP TG, INT n, INT k, 
	VECTOR_OP Reps, VECTOR_OP Stab, INT f_v, INT f_vv)
{
	INT i, l = 1, res;
	INT X[10000];
	INTEGER_OB len_ob;
	VECTOR_OB the_set;
	INT idx, f_found;
	LABRA_OB aut;
	SYM_OB ago;
	// PERMUTATION_OB transporter;

	len_ob.m_i(0);
	Reps->m_il(0);
	Stab->m_il(0);
	n_choose_k_first(X, n, k);
#if 0
	for (i = 0; i < k; i++) 
		Y[i] = X[i];
#endif
	while (TRUE) {
		// printf("%3ld: ", l);
		// print_set(Y, n, k);
		// geo_canonicize_set(n, k, Y, G, TG, FALSE, FALSE, &transporter);
		res = geo_maxtest_set(n, k, X, G, TG, &aut, &ago, FALSE, FALSE);
		if (res == -1) {
			if (f_vv) {
				printf("new orbit !\n");
				printf("ago = ");
				ago.println();
				}
		
			// print_set(X, n, k);
			the_set.m_il(k);
			for (i = 0; i < k; i++) {
				the_set.m_ii(i, X[i]);
				}
			Reps->search(len_ob.s_i(), TRUE, &the_set, &idx, &f_found);
			if (!f_found) {
				if (f_vv) {
					printf("new orbit (cand %ld): ", l);
					print_set(X, n, k);
					the_set.println();
					}
				Reps->insert_at(&len_ob, idx, &the_set);
				len_ob.dec();
				Stab->insert_at(&len_ob, idx, &aut);
				}
			else {
				return error("calc_k_sets() new orbit already there !");
				}
			}
		if (res == -1)
			res = n;
		if (!n_choose_k_next_at(X, n, k, res))
			break;
#if 0
		for (i = 0; i < k; i++) 
			Y[i] = X[i];
#endif
		l++;
		}

	// note:
	// \sum_orbits go / ago(orbit) =  \left( n \atop k \right)
	// where go is the order of the group G
	Reps->realloc(&len_ob);
	Stab->realloc(&len_ob);
	if (f_v) {
		printf("found %ld orbits on %ld sets\n", Reps->s_li(), k);
		}
	// Reps->println();
	return OK;
}

INT calc_k_stab_orders(VECTOR_OP k_stabs, VECTOR_OP k_stab_orders)
{
	LABRA_OP stab;
	SYM_OP go;
	INT i, l;

	l = k_stabs->s_li();
	k_stab_orders->m_il(l);
	for (i = 0; i < l; i++) {
		go = k_stab_orders->s_i(i);
		stab = (LABRA_OP) k_stabs->s_i(i);
		stab->group_order(go);
		}
	return OK;
}

INT n_choose_k_first(INT *X, INT n, INT k)
{
	INT i;
	
	for (i = 0; i < k; i++) 
		X[i] = i;
	return TRUE;
}

INT n_choose_k_next(INT *X, INT n, INT k)
{
	INT i, j, r, ii;
	
	for (i = k - 1; i >= 0; i--) {
		j = X[i];
		r = k - i;
		if (j + 1 <= n - r) {
			for (ii = 0 ; ii < r; ii++)
				X[i + ii] = j + 1 + ii;
			return TRUE;
			}
		}
	return FALSE;
}

INT n_choose_k_next_at(INT *X, INT n, INT k, INT a)
{
	INT i, j, r, ii;
	
	for (i = k - 1; i >= 0; i--) {
		j = X[i];
		if (j > a)
			continue;
		r = k - i;
		if (j + 1 <= n - r) {
			for (ii = 0 ; ii < r; ii++)
				X[i + ii] = j + 1 + ii;
			return TRUE;
			}
		}
	return FALSE;
}

void print_set(INT *X, INT n, INT k)
{
	INT i, ii;
	
	ii = 0;
	for (i = 0; i < n; i++) {
		if (X[ii] == i) {
			printf("X");
			ii++;
			}
		else
			printf(".");
		}
	printf("\n");
}

INT KM_tk_print_asc(LABRA_OP G, BYTE *g_label, BYTE *g_label_tex, MATRIX_OP M, INT gl_t, INT gl_k)
{
	BYTE str1[256];
	FILE *fp;
	BYTE fname[256];
	INT i, m, n, l, deg;
	VECTOR_OP gen;
	SYM_OB go;

	deg = G->s_degree_i();
	gen = G->s_G();
	sprintf(fname, "KM_%s_t%ld_k%ld.txt", g_label, gl_t, gl_k);
	fp = fopen(fname, "w");
	G->group_order(&go);
	str1[0] = 0;
	go.sprint(str1);
	fprintf(fp, "%% this file:    %s\n", fname);
	fprintf(fp, "%% group:        %s\n", g_label);
	fprintf(fp, "%% group (tex):  %s\n", g_label_tex);
	fprintf(fp, "%% order:        %s\n", str1);
	fprintf(fp, "%% degree:       %ld\n", deg);
	l = gen->s_li();
	fprintf(fp, "%% # generators: %ld\n", l);
	for (i = 0; i < l; i++) {
		str1[0] = 0;
		gen->s_i(i)->sprint(str1);
		fprintf(fp, "%% %s\n", str1);
		}
	fflush(fp);
	
	fprintf(fp, "%% t, k:\n");
	fprintf(fp, "%% %ld %ld\n", gl_t, gl_k);

	m = M->s_hi();
	n = M->s_li();

	fprintf(fp, "%% m, n:\n");
	fprintf(fp, "%ld %ld\n", m, n);
	M->fprint_raw(fp);
	
	fclose(fp);

	/* printf("written file %s of size %ld\n", fname, file_size(fname));
	fflush(stdout); */
	return TRUE;
}

