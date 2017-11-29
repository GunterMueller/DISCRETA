/* fg_ext.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>
// #include <TG/tg.h> /* for rpe_by_stab() */

#define MAX_NW 64

#define USE_OBJ_LIB
#define DEBUG_BASE_IMAGE

#undef malloc
#undef free

static INT fg_reduce_possible_extensions(FG_OP G, 
	INT p, INT nb_cr, INT *cr_idx, 
	BYTE *ex, INT f_verbose);

#define BUF_SIZE 1024

#if TEXDOCU
INT fg_cmd_file_for_order_n(INT n_from, INT n_to, 
	INT f_calc_aut_group, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family)
#endif
{
	fg_pvm_cmd_file_for_order_n(n_from, n_to, FALSE /* f_pvm */, 
		f_calc_aut_group, f_calc_aut_classes, 
		f_calc_sgl, f_calc_sylow_type, f_calc_family);
	return OK;
}

#if TEXDOCU
INT fg_pvm_cmd_file_for_order_n(INT n_from, INT n_to, INT f_pvm, 
	INT f_calc_aut_group, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family)
#endif
{
	BYTE str[BUF_SIZE];
	FILE *fp, *fp1;
	INT i, j, l, p, e, i1;
	VECTOR_OB vp, ve;
	
	fp = fopen("fg.cmd", "w");
#if 0
	system("date >a");
	fp1 = fopen("a", "r");
	fgets(str, BUF_SIZE, fp1);
	fclose(fp1);
#endif
	date_as_string(str);
	if (n_from < 2)
		n_from = 2;
	fprintf(fp, "REM\n");
	fprintf(fp, "fg.cmd written at %s\n", str);
	fprintf(fp, "by function fg_pvm_cmd_file_for_order_n(), "
		"order: [%ld,\\ldots,%ld] f_pvm = %ld\n", n_from, n_to, f_pvm);
	fprintf(fp, "MER\n\n");
	fprintf(fp, "compute_aut_group %ld\n", f_calc_aut_group);
	fprintf(fp, "compute_aut_classes %ld\n", f_calc_aut_classes);
	fprintf(fp, "compute_sgl %ld\n", f_calc_sgl);
	fprintf(fp, "compute_sylow_type %ld\n", f_calc_sylow_type);
	fprintf(fp, "compute_family %ld\n", f_calc_family);
	fprintf(fp, "reduce_by_classes 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	for (i = n_from; i <= n_to; i++) {
		fprintf(fp, "\n");
		fprintf(fp, "REM\n");
		fprintf(fp, "********** ORDER %ld\n", i);
		fprintf(fp, "MER\n");
		factor_integer(i, &vp, &ve);
		l = vp.s_li();
		p = vp.s_ii(0);
		e = ve.s_ii(0);
		if (l == 1 && e == 1) {
			fprintf(fp, "Z %ld \"%ld#1\"\n", i, i);
			fprintf(fp, "\n");	
			}
		else {
			for (j = l - 1; j >= 0; j--) {
				p = vp.s_ii(j);
				e = ve.s_ii(j);
				i1 = i / p;
				if (f_pvm) {
					fprintf(fp, "extendall_pvm %ld %ld\n", i1, p);
					}
				else {
					fprintf(fp, "extendall %ld %ld\n", i1, p);
					}
				}
			}
		}
	fprintf(fp, "END\n");
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT FG_OB::print_extension_matrix(CLASS_REP_OP R, BYTE *ex, INT p)
#endif
{
	fprint_extension_matrix(stdout, R, ex, p);
	return OK;
}

#if TEXDOCU
INT FG_OB::fprint_extension_matrix(FILE *fp_txt, 
	CLASS_REP_OP R, BYTE *ex, INT p)
#endif
{
	INT base_im[MAX_NW];
	INT coset_rep[MAX_NW];
	PERMUTATION_OB aut;
	INT i1, P, n, o, k;
	INT nb_classes, base_len;

	fprintf(fp_txt, "p-extension matrix for group %s and prime %ld:\n", 
		s_label_s(), p);
	fflush(fp_txt);
	n = s_n_i();
	nb_classes = R->s_nb_classes_i();
	base_len = R->s_base_len_i();
	fprintf(fp_txt, "class  aut: [base-images] element-order extension-matrix\n");
	for (i1 = 0; i1 < nb_classes; i1++) {
		fprintf(fp_txt, "%3ld [", i1);
		for (k = 0; k < base_len; k++) {
			base_im[k] = R->s_R_iji(i1, k);
			fprintf(fp_txt, "%3ld", base_im[k]);
			if (k < base_len - 1)
				fprintf(fp_txt, ", ");
			}
		bi2rep(base_im, coset_rep);
		rep2aut(coset_rep, &aut);
		aut.order(&o);
		fprintf(fp_txt, "] %5ld ", o);
		for (P = 0; P < n; P++) {
			if (ex[i1 * n + P]) {
				fprintf(fp_txt, "X");
				}
			else {
				fprintf(fp_txt, ".");
				}
			}
		fprintf(fp_txt, "\n");
		}
	fflush(fp_txt);
	return OK;
}

#if TEXDOCU
INT FG_OB::calc_Extension_matrix(FILE *fp_txt, 
	MATRIX_OP M, INT p, INT f_reduced, INT f_verbose)
#endif
{
	CLASS_REP_OP R;
	BYTE *ex;
	INT n, nb_classes, i, j, a, nb_poss = 0;

	if (!s_f_has_aut_classes_i())
		return f_error(fp_txt, "FG_OB::calc_Extension_matrix(): !s_f_has_aut_classes_i()");
	R = s_aut_classes();
	nb_classes = R->s_nb_classes_i();
	n = s_n_i();
	calc_extension_matrix(fp_txt, &ex, p, f_verbose);
	if (f_reduced) {
		// rpe_by_stab(fp_txt, ex, p, f_verbose);
		// we dont have classes at the moment
		}
	M->m_ilih_n(n, nb_classes);
	for (i = 0; i < nb_classes; i++) {
		for (j = 0; j < n; j++) {
			if (ex[i * n + j]) {
				M->m_iji(i, j, 1);
				nb_poss++;
				}
			}
		}
	my_free(ex);
	return nb_poss;
}

#if TEXDOCU
INT FG_OB::calc_extension_matrix(FILE *fp_txt, BYTE **ex, INT p, INT f_verbose)
#endif
{
	CLASS_REP_OP R;
	INT base[MAX_NW];
	INT base_im[MAX_NW];
	INT coset_rep[MAX_NW];
	INT P_base_im[MAX_NW];
	ZE_OP ze;
	PERMUTATION_OB aut;
	BYTE *ex1 = NIL;
	INT n, size;
	INT i, j, k, g_i, P, i1, o;
	INT nb_classes, base_len;
	INT nb_poss = 0;
	INT f_very_verbose = FALSE;

	if (!s_f_has_aut_classes_i())
		return f_error(fp_txt, "FG_OB::calc_extension_matrix(): !s_f_has_aut_classes_i()");
	R = s_aut_classes();
	nb_classes = R->s_nb_classes_i();
	base_len = R->s_base_len_i();
	/* ago = s_ago_i(); */
	n = s_n_i();
	size = nb_classes * n * sizeof(BYTE);
	if (f_verbose) {
		if (TRUE)
			f_very_verbose = TRUE;
		}
	if (base_len != s_nb_ze_i())
		return f_error(fp_txt, "FG::calc_extension_matrix(): base_len != s_nb_ze_i()");
	if (f_very_verbose) {
		fprintf(fp_txt, "in fg_possible_extensions() ago = ");
		s_ago()->fprintln(fp_txt);
		fprintf(fp_txt, "in %ld classes\n", nb_classes);
		}
	ex1 = (BYTE *) my_malloc(size, "fg_ext.C: calc_extension_matrix");
	if (ex1 == NIL)
		return f_error(fp_txt, "FG::calc_extension_matrix(): "
		"no memory for ex1");
	*ex = ex1;
	for (i = 0; i < base_len; i++) {
		ze = s_ze_i(i);
		g_i = ze->s_n0_i();
		base[i] = g_i;
		}
	for (i = 0; i < nb_classes * n; i++)
		ex1[i] = 0;
	
	/* Bestimme alle moeglichen Kombinationen von 
	 * Automorphismen (ii) 
	 * und p-ten Potenzen (P): */

	if (f_very_verbose) {
		fprintf(fp_txt, "class  aut: [base-images] element-order\n");
		}
	for (i1 = 0; i1 < nb_classes; i1++) {
		/* ii = cr_idx[i1];
		int2aut(ii, &aut); */
		if (f_very_verbose) {
			fprintf(fp_txt, "%3ld [", i1);
			}
		for (k = 0; k < base_len; k++) {
			base_im[k] = R->s_R_iji(i1, k);
			if (f_very_verbose) {
				fprintf(fp_txt, "%3ld", base_im[k]);
				if (k < base_len - 1)
					fprintf(fp_txt, ", ");
				}
			}
		bi2rep(base_im, coset_rep);
		rep2aut(coset_rep, &aut);
		aut.order(&o);
		if (f_very_verbose) {
			fprintf(fp_txt, "] %5ld ", o);
			fflush(fp_txt);
			}
		
		for (P = 0; P < n; P++) {
		
			/* Bestimme den inneren Automorphismus 
			 * von P: */
			for (i = 0; i < base_len; i++) {
				g_i = base[i];
				P_base_im[i] = gt_conj(g_i, P);
				}
			/* bi2int(P_base_im, &P_aut); */
			
			/* P muss Fixpunkt sein: */
			if (aut.s_ii(P) - 1 != P)
				goto loop;

			/* Die p-te Potenz des Automorphismus 
			 * muss mit dem inneren Autom. von P 
			 * uebereinstimmen: */
			for (i = 0; i < base_len; i++) {
				g_i = base[i];
				for (j = 0; j < p; j++)
					g_i = aut.s_ii(g_i) - 1;
				if (g_i != P_base_im[i])
					goto loop;
				}
			/* Automorphismen stimmen genau 
			 * dann ueberein, wenn sie auf den 
			 * Generatoren uebereinstimmen. */
			ex1[i1 * n + P] = 1;
			nb_poss++;
			if (f_very_verbose) {
				fprintf(fp_txt, "X");
				}
			goto loop1;
loop:
			if (f_very_verbose) {
				fprintf(fp_txt, ".");
				}
loop1:
			;
			}
			if (f_very_verbose) {
				fprintf(fp_txt, "\n");
				}
		}
	if (f_verbose) {
		fprintf(fp_txt, "fg_possible_extensions(): "
		"%ld = %ld * %ld possible: %ld\n", 
		nb_classes * n, nb_classes, n, nb_poss);
		fflush(fp_txt);
		}
	
	return nb_poss;
}

#if TEXDOCU
INT FG_OB::rpe_by_stab(FILE *fp_txt, BYTE *ex, INT p, INT f_verbose)
#else
does not work at the moment as it uses LABRA\_TG
#endif
{
	return error("error: FG_OB::rpe_by_stab called");
#if 0
	CLASS_REP_OP R;
	INT base_im[MAX_NW];
	INT coset_rep[MAX_NW];
	PERMUTATION_OB aut;
	INT nb_poss1, nb_poss, i1, P, P1, n, j, k, a;
	INT nb_classes, base_len, deg;
	INT f_very_verbose = FALSE;
	/* ARRAY < LABRA_TG > Stab(R->s_nb_classes_i()); */
	
	if (!s_f_has_aut_classes_i())
		return f_error(fp_txt, "FG_OB::reduce(): !s_f_has_aut_classes_i()");
	R = s_aut_classes();
	if (f_verbose) {
		fprintf(fp_txt, "FG_OB::reduce possible extensions by stabilizers\n");
		fflush(fp_txt);
		}
	nb_classes = R->s_nb_classes_i();
	base_len = R->s_base_len_i();
	n = s_n_i();
	nb_poss = 0;
	for (i1 = 0; i1 < nb_classes; i1++) {
		for (P = 0; P < n; P++) {
			if (ex[i1 * n + P])
				nb_poss++;
			}
		}
	nb_poss1 = nb_poss;

	
	for (i1 = 0; i1 < nb_classes; i1++) {
		/* ii = cr_idx[i1]; */
		for (P = 0; P < n; P++) {
			if (ex[i1 * n + P] == 0)
				continue;
			{
			int nb_gen, *gen;

			/* Stab.Add(i1, n);
			LABRA_TG &L = Stab[i1]; */
			
			LABRA_TG L;
			
			nb_gen = R->s_stab_i(i1)->s_hi();
			gen = (int *) my_malloc(nb_gen * n * sizeof(int), "fg_ext.C: rpe_by_stab");

#if 0
			fprintf(fp_txt, "extension (class,P)=(%ld,%ld), "
				"found %ld generators for centralizer of\n", i1, P, nb_gen);
			for (k = 0; k < base_len; k++) {
				a = R->s_R_iji(i1, k);
				fprintf(fp_txt, "%ld ", a);
				}
			fprintf(fp_txt, "\n");
			fflush(fp_txt);
#endif
			for (j = 0; j < nb_gen; j++) {
				for (k = 0; k < base_len; k++) {
					a = R->s_stab_ijki(i1, j, k);
					base_im[k] = a;
					}
				bi2rep(base_im, coset_rep);
				rep2aut(coset_rep, &aut);
				/* aut.println(); */
				for (k = 0; k < n; k++) {
					a = aut.s_ii(k);
					gen[j * n + k] = a;
					}
#ifdef DEBUG_BASE_IMAGE
				{
				INT g_k, b;

				for (k = 0; k < base_len; k++) {
					g_k = R->s_base_ii(k);
					a = aut.s_ii(g_k) - 1;
					b = R->s_stab_ijki(i1, j, k);
					if (a != b) {
						printf("a = %ld != b = %ld", a, b);
						return f_error(fp_txt, "error in rpe_by_stab");
						}
					}
				}
#endif
				}

			/* printf("calling L.init()\n");
			fflush(stdout); */
			L.Init(n);
			/* printf("calling L.InitErz()\n");
			fflush(stdout); */
			L.InitErz(nb_gen, gen);
			/* printf("calling L.gens_to_labra()\n");
			fflush(stdout); */
			/* L.gens_to_labra(); */
			/* L.Ordnung().Print(0); */
			
			/* printf("Gruppe=\n");
			L.Print(); */
			 
			if (L.pi_i(P+1) != 1) {
			 /* printf("cycle %ld %ld\n", (INT) 1, (INT) L.pi_i(P+1)); */
			 L.Cycle(1,L.pi_i(P+1));
			 }
			L.Nachfolger_start(1);
			if (f_very_verbose) {
				fprintf(fp_txt, "row %ld orbit of %ld: ", i1, (INT) P);
				fflush(fp_txt);
				}
			while (!L.IsLastNachfolger()) {
				P1 = L.Nachfolger(1);
				if (f_verbose) {
					fprintf(fp_txt, "%ld ", P1 - 1);
					fflush(fp_txt);
					}
				if (P1 - 1 != P) {
					ex[i1 * n + P1 - 1] = 0;
					nb_poss--;
					}
				}
			if (f_very_verbose) {
				fprintf(fp_txt, "\n");
				}

			my_free(gen);
			}
			}
		}
	if (f_verbose) {
		fprintf(fp_txt, "reduced from %ld to %ld\n", nb_poss1, nb_poss);
		fflush(fp_txt);
		}
	return nb_poss;
#endif
}

#if TEXDOCU
static INT fg_reduce_possible_extensions(FG_OP G, 
	INT p, INT nb_cr, INT *cr_idx, 
	BYTE *ex, INT f_verbose)
#endif
{
	INT base[MAX_NW];
	INT base_im[MAX_NW];
	INT P_base_im[MAX_NW];
	ZE_OP ze;
	PERMUTATION_OB aut, aut1, beta, beta_v, tmp, aut_beta;
	INT ago, n, size;
	INT i, j, k, g_i, P, P1, P_beta, i1, ii, iii, iv, p1;
	INT aut1_idx, aut_beta_idx, nb_poss = 0;
	INT f_very_verbose = FALSE;

	if (G->s_ago()->s_obj_k() != INTEGER)
		return error("fg_reduce_possible_extensions() ago not an integer");
	ago = ((INTEGER_OP) G->s_ago())->s_i();
	n = G->s_n_i();

	nb_poss = 0;
	for (i1 = 0; i1 < nb_cr; i1++) {
		for (P = 0; P < n; P++) {
			if (ex[i1 * n + P])
				nb_poss++;
			}
		}
	for (i1 = 0; i1 < nb_cr; i1++) {
		ii = cr_idx[i1];
		G->int2aut(ii, &aut);
		if (f_very_verbose) {
			if (i1 % 100 == 0) {
				printf("%ld %ld: ", i1, ii);
				if ((i1 > 0) && (i1 % (100 * 10) == 0))
					printf("\n");
				fflush(stdout);
				}
			}

		for (P = 0; P < n; P++) {
			if (!ex[i1 * n + P])
				continue;
			for (p1 = 1; p1 <= 1 /* < p */; p1++) {
				/* p1 - te Potenz von aut nach aut1, 
				 * von P nach P1: */
				aut.copy(&aut1);
				k = 1;
				while (k < p1) {
					aut.mult(&aut1, &tmp);
					tmp.swap(&aut1);
					k++;
					}
				G->aut2int(&aut1, &aut1_idx);
				if (cr_idx[aut1_idx] != aut1_idx)
					return error("fg_reduce_possible_extensions(): "
						"reduce not allowed together "
						"with conjugacy classes !");

				P1 = G->gt_power(P, p1);
				if (ex[aut1_idx * n + P1]) {
			/* Nebenklasse streichen: */
			for (iii = 0; iii < ago; iii++) {
				G->int2aut(iii, &beta);
				G->aut2int(&beta, &iv);
				if (iv != iii) {
					printf("iii = %ld "
					"iv = %ld\nbeta = ", iii, iv);
					beta.println();
					return error("iv != iii");
					}
				beta.invers(&beta_v);
				beta_v.mult(&aut1, &tmp);
				tmp.mult(&beta, &aut_beta);
				P_beta = beta.s_ii(P1) - 1;
				G->aut2int(&aut_beta, &aut_beta_idx);
				if (aut_beta_idx != ii || P_beta != P) {
					if (ex[aut_beta_idx * n + P_beta]) {
						ex[aut_beta_idx * n + P_beta] = 0;
						nb_poss--;
						}
					}
				}
					} /* if (ex[]) */
				} /* for p1 */
			}
		}
	if (f_verbose) {
		printf("fg_reduce_possible_extensions()|"
			" reduced to %ld.\n", nb_poss);
		fflush(stdout);
		}
	return nb_poss;
}

#if TEXDOCU
INT FG_OB::calc_aut_classes_using_file(FILE *fp_txt, INT f_v, INT f_vv)
#endif
{
	CLASS_REP_OB R;
	
	fprintf(fp_txt, "FG::calc_aut_classes_using_file...");
	fflush(fp_txt);
	R.calc_aut_classes_using_file(fp_txt, this, f_v, f_vv);
	R.swap(s_aut_classes());
	s_f_has_aut_classes()->m_i(TRUE);
	fprintf(fp_txt, "FG::calc_aut_classes_using_file...finished!");
	fflush(fp_txt);
	return OK;
}

#if TEXDOCU
INT FG_OB::init_extension_matrix(FILE *fp_txt, INT p, BYTE **ex, INT f_verbose)
#endif
{
	INT nb_poss;
	BYTE *ex1;

	if (!s_f_has_aut_classes_i()) {
		calc_aut_classes_using_file(fp_txt, f_verbose, FALSE);
		}
	nb_poss = calc_extension_matrix(fp_txt, &ex1, p, f_verbose);
	// nb_poss = rpe_by_stab(fp_txt, ex1, p, f_verbose);
	// we dont have classes at the moment
	if (f_verbose) {
		CLASS_REP_OP R;
		
		R = s_aut_classes();
		fprint_extension_matrix(fp_txt, R, ex1, p);
		}
	*ex = ex1;
	return nb_poss;
}

#if TEXDOCU
INT FG_OB::calc_aut_sgl_syl(FILE *fp_txt, 
	INT f_calc_aut, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, 
	INT f_v, INT f_vv)
#endif
{
	if (f_calc_aut) {
		fprintf(fp_txt, "aut:\n");
		fflush(fp_txt);
		Aut(f_v, f_vv);
		}
	if (f_calc_aut_classes) {
		fprintf(fp_txt, "aut-classes:\n");
		fflush(fp_txt);
		calc_aut_classes_using_file(fp_txt, f_v, f_vv);
		}
	if (f_calc_sylow_type) {
		fprintf(fp_txt, "Sylow-type: ");
		fflush(fp_txt);
		sylow_type(f_vv);
		s_sylow_type()->println();
		fflush(stdout);
		}

	if (f_calc_sgl) {
		fprintf(fp_txt, "SGL:\n");
		fflush(fp_txt);
		sgl(FALSE /* f_v */, FALSE /* f_vv */, FALSE /* f_vvv */ );
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::do_extension_Zp(INT p, 
	INT f_calc_aut, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family, 
	INT f_reduce_by_classes, 
	VECTOR_OP theGroups, VECTOR_OP theFamilies)
#endif
{
	INT base_im[MAX_NW];
	INT coset_rep[MAX_NW];
	FG_OB G1;
	FG_OP G2;
	PERMUTATION_OB aut;
	INT n, m, /* ago, */ old_nb_groups;
	INT i, i1, P, cnt, nb_poss, k;
	INT nb_classes, base_len;
	INT iso_idx;
	/* INT nb_cr, *cr_idx; */
	BYTE *ex = NIL;
	BYTE str[256];
	INT f_verbose = TRUE;
	CLASS_REP_OP R;

	if (theGroups->emptyp())
		theGroups->m_il(0);
	if (theFamilies->emptyp())
		theFamilies->m_il(0);
	old_nb_groups = theGroups->s_li();
	/* ago = s_ago_i(); */
	n = s_n_i();
	printf("do_extension_Zp() p = %ld "
		"f_calc_aut = %ld f_calc_aut_classes = %ld "
		"f_calc_sgl = %ld f_calc_sylow_type = %ld f_calc_family = %ld "
		"f_reduce_by_classes = %ld\n", 
		p, f_calc_aut, f_calc_aut_classes, 
		f_calc_sgl, f_calc_sylow_type, f_calc_family, 
		f_reduce_by_classes);
	fflush(stdout);
	
	nb_poss = init_extension_matrix(stdout, p, &ex, f_verbose);
	R = s_aut_classes();
	nb_classes = R->s_nb_classes_i();
	base_len = R->s_base_len_i();

	cnt = 0;
	for (i1 = 0; i1 < nb_classes; i1++) {
		for (P = 0; P < n; P++) {
			if (!ex[i1 * n + P])
				continue;
			for (k = 0; k < base_len; k++) {
				base_im[k] = R->s_R_iji(i1, k);
				}
			bi2rep(base_im, coset_rep);
			rep2aut(coset_rep, &aut);
			
			if (f_verbose) {
				printf("group %s %ld-extension %ld/%ld: "
					"(aut, P) = ([", 
					s_label_s(), p, cnt, nb_poss);
				for (k = 0; k < base_len; k++) {
					printf("%ld", (INT) base_im[k]);
					if (k < base_len - 1)
						printf(", ");
					}
				printf("], %ld) = (", P);
				printf(", ");
				print_int(P);
				printf(")\n");
				fflush(stdout);
				}
			G1.init_extension(this, p, P, &aut);
			cnt++;

			G1.theG_Zp_extension(this, &aut, FALSE /* f_verbose */);
			if (!G1.is_associative(FALSE /* f_verbose */))
				goto loop;
			iso_idx = G1.find_isomorphic_group(theGroups, FALSE /* f_verbose */);
			if (iso_idx >= 0)
				goto loop;
			
			G1.calc_aut_sgl_syl(stdout, f_calc_aut, f_calc_aut_classes, 
				f_calc_sgl, f_calc_sylow_type, 
				TRUE /* f_v */, FALSE /* f_vv */);

#if 0
			if (f_calc_aut) {
				printf("aut:\n");
				fflush(stdout);
				G1.Aut(TRUE /* f_verbose */, 
					FALSE /* f_very_verbose */);
				}
			if (f_calc_aut_classes) {
				printf("aut-classes:\n");
				fflush(stdout);
				G1.calc_aut_classes_using_file(stdout, f_verbose, FALSE /* f_vv */);
				}
			if (f_calc_sylow_type) {
				printf("Sylow-type: ");
				fflush(stdout);
				G1.sylow_type(FALSE /* f_verbose */);
				G1.s_sylow_type()->println();
				fflush(stdout);
				}
			
			if (f_calc_sgl) {
				printf("SGL:\n");
				fflush(stdout);
				G1.sgl(FALSE /* f_v */, FALSE /* f_vv */, FALSE /* f_vvv */ );
				}
#endif

			m = theGroups->s_li() + 1;
			sprintf(str, "%ld#%ld", G1.s_n_i(), m);
			G1.s_label()->init(str);
			G1.s_m()->m_i(m);
			G1.GmZ_Gd_orders(stdout, TRUE /* f_verbose: print orders */);
			printf("FG::do_extension_Zp()|"
				"found new group %s\n", 
				G1.s_label_s());
			
			if (f_calc_family) {
				printf("Family:\n");
				fflush(stdout);
				G1.determine_family(theFamilies, 
					TRUE /* f_verbose */, 
					FALSE /* f_very_verbose */);
				}
			
			theGroups->inc();
			G1.swap(theGroups->s_i(theGroups->s_li() - 1));
			
loop:
			;
			}
		}
	printf("\nerzeugt: %ld; "
		"neu gefundene Gruppen: %ld total %ld\n", 
		cnt, theGroups->s_li() - old_nb_groups, 
		theGroups->s_li());
	/* my_free(cr_idx); */
	if (ex) {
		my_free(ex);
		ex = NIL;
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::determine_family(VECTOR_OP theFamilies, 
	INT f_verbose, INT f_very_verbose)
#endif
{
	INT isoc_idx, o;
	
	isoc_idx = find_isoclinic_group(theFamilies, f_very_verbose);
	if (isoc_idx >= 0) {
		o = ((FG_OP) theFamilies->s_i(isoc_idx))->s_o_i();
		if (f_verbose) {
			printf("determine_family() "
			"found isoclinism class o = %ld\n", o);
			fflush(stdout);
			}
		s_o()->m_i(o);
		}
	else {
		/* new family: */
		o = theFamilies->s_li() + 1;
			/* first family gets numero 1 */
		if (f_verbose) {
			printf("determine_family() "
			"new isoclinism class o = %ld\n", o);
			fflush(stdout);
			}
		s_o()->m_i(o);
		theFamilies->inc();
		((SYM_OP) this)->copy(
			theFamilies->s_i(theFamilies->s_li() - 1));
		}
	return OK;
}

#if TEXDOCU
INT extend_Zp_vec_of_groups(VECTOR_OP V, INT p, 
	INT f_calc_aut, INT f_calc_aut_classes, 
	INT f_calc_sgl, INT f_calc_sylow_type, INT f_calc_family, 
	INT f_reduce_by_classes, 
	VECTOR_OP theGroups, VECTOR_OP theFamilies)
#endif
{
	FG_OP G;
	INT i;
	
	for (i = 0; i < V->s_li(); i++) {
		G = (FG_OP) V->s_i(i);
		G->do_extension_Zp(p, 
			f_calc_aut, f_calc_aut_classes, f_calc_sgl, f_calc_sylow_type, f_calc_family, 
			f_reduce_by_classes, theGroups, theFamilies);
		}
	return OK;
}

#endif /* SOLVABLE_TRUE */



