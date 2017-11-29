/* report.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>
#include <stdlib.h> // for system()

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

#define BUFSIZE 16000
#define INTERSECTION_MATRIX_COL_WIDTH_FOR_SPLIT 2
#define INTERSECTION_TOO_LARGE_IF_N_COLS 50
#define INTERSECTION_EQNS_COL_WIDTH_FOR_SPLIT 25
#define KM_MATRIX_COL_WIDTH_FOR_SPLIT 50
#define SOLUTION_VECTOR_WIDTH_FOR_SPLIT 100

#define F_USE_MULTIPLICITIES FALSE
#define F_USE_COMPLEMENTS FALSE

static INT report_plesken_matrix(INT k, VECTOR_OP MM, VECTOR_OP RR, VECTOR_OP stab_go, SYM_OP go);
static INT read_until_lambdaend(FILE *fp);
static INT print_and_split_matrix(MATRIX_OP M, INT width_in_cols, INT f_clearpage);
static INT km_report_group(INT k, SYM_OP go, VECTOR_OP G_gen, INT f_sims_chain);
static INT km_report_orbits(INT k, SYM_OP go, VECTOR_OP stab_go, VECTOR_OP RR);
static INT km_report_orbits_html(FILE *fp, BYTE *group_label, 
	INT k, SYM_OP go, VECTOR_OP stab_go, VECTOR_OP RR);
static INT km_report_KM_matrices_html(FILE *fp_html, BYTE *group_label, INT k, VECTOR_OP MM);
static INT km_report_KM_matrix_html(FILE *fp, BYTE *group_label, INT t, INT k, MATRIX_OP Mtk);
static INT design_output_geo(VECTOR_OP gen, INT k, VECTOR_OP RR, MATRIX_OP X);
static INT design_print_base_blocks(VECTOR_OP gen, 
	VECTOR_OP RR, INT k, MATRIX_OP X, SYM_OP go, VECTOR_OP stab_go);
static INT base_block_fname(BYTE *KM_name, INT lambda, 
	BYTE *base_fname, BYTE *fname);
static INT design_clear_base_blocks_file(BYTE *KM_name, INT lambda);
static INT design_print_base_blocks_file(BYTE *KM_name, INT lambda, INT no,
	VECTOR_OP RR, INT k, MATRIX_OP X);
static INT print_nice_plesken_matrix(MATRIX_OP Ainf_block_wise, 
	INT k_min, INT k_max, 
	VECTOR_OP RR, VECTOR_OP stab_go, 
	SYM_OP go, FILE *fp);
static INT print_nice_matrix(MATRIX_OP M, INT t, INT k, 
	VECTOR_OP RR, VECTOR_OP stab_go, 
	SYM_OP go, FILE *fp);
static INT prepare_set_tex(VECTOR_OP R, SYM_OP go, SYM_OP stab_go, BYTE *s);
static INT write_r0_tex_file();

#if TEXDOCU
void design_report(BYTE *KM_fname, INT f_html, INT f_select, 
	INT design_select_lambda, INT select_first, INT select_length)
#endif
{
	BYTE base_fname[1000];
	BYTE html_fname[1000];
	BYTE group_label[1000];
	BYTE str[100000];
	FILE *fp_html;
	BYTE buf[BUFSIZE];
	INT v, t, k, lambda, m1;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	VECTOR_OB orbit_length, Eqns;
	MATRIX_OB Mtk, Mendelsohn_mtx, LAmbda;
	VECTOR_OB Mendelsohn_RHS, alpha_i;
	VECTOR_OB Mendel_RHS; // vector of mendelsohn RHS for s=1,2,3... (s=0 is empty)
	VECTOR_OB Lambda, Coeff, Constant_term;
	VECTOR_OB Koehler_Constant_term_v, Koehler_Coeff_v;
	VECTOR_OB S;
	INT s, s_max;
	SYM_OB go, ol;
	INT i, m, n, l, a, no;
	FILE *fp;
	// INT f_output_geo = TRUE;
	VECTOR_OB block_types, bai_ref, bai_refv, bai_data;
	INTEGER_OB bai_len;
	MATRIX_OB X, Y;
	BYTE *offset;
	SYM_OB Sum;
	VECTOR_OB orbits, orbits_c;
	INT highest_layer, num_layers;


	write_r0_tex_file();
	lambda = 0;
	printf("%% report() f_select = %ld design_select_lambda = %ld\n", 
		f_select, design_select_lambda);
	printf("%% select_first = %ld, select_length = %ld\n", 
		select_first, select_length);
	printf("KM-file: %s\n", KM_fname);
	fflush(stdout);
	
	strcpy(base_fname, KM_fname);
	l = strlen(base_fname);
	if (l > 4 && strcmp(base_fname + l - 4, ".txt") == 0)
		base_fname[l - 4] = 0;
	printf("%% base_fname = %s\n", base_fname);
	fflush(stdout);

	KM_fname_extract_group_label(KM_fname, group_label);
	printf("%% group_label = %s\n", group_label);
	fflush(stdout);
	
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m1 = k;
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	printf("v=%ld t=%ld k=%ld lambda=%ld m=%ld\n", v, t, k, lambda, m1);
	fflush(stdout);
	
	if (f_html) {
		sprintf(html_fname, "%s.html", base_fname);
		fp_html = fopen(html_fname, "w");
		fprintf(fp_html, "<html>\n");
		fprintf(fp_html, "<head><title>DISCRETA report: "
			"%ld-(%ld,%ld,\\lambda) designs</title></head>\n", t, v, k);
		fprintf(fp_html, "<body>\n\n<p>\n\n\n");
		fprintf(fp_html, "<center>\n<h2> <a href=\""
			"http://www.mathe2.uni-bayreuth.de/~discreta"
			"\"> DISCRETA </a> report on </h2>\n"
			"<h2>%ld-(%ld,%ld,\\lambda) designs</h2>\n", t, v, k);
		fprintf(fp_html, "<p>\n<h2>generated from the group %s</h2>\n<p>\n\n", group_label);
		fprintf(fp_html, "<a href=\"http://www.mathe2.uni-bayreuth.de/betten/anton.html\"> Anton Betten, </a> \n");
		fprintf(fp_html, "<a href=\"http://www.mathe2.uni-bayreuth.de/~evi\"> Evi Haberberger, </a> \n");
		fprintf(fp_html, "<a href=\"http://www.mathe2.uni-bayreuth.de/people/laue.html\"> Reinhard Laue, </a>\n");
		fprintf(fp_html, "<a href=\"http://did.mat.uni-bayreuth.de/wassermann/wassermann.html\"> Alfred Wassermann, </a>\n<p>\n");
		fprintf(fp_html, "Mathematical department, University of Bayreuth\n\n");
		fprintf(fp_html, "</center>\n\n<p>\n\n");
		}
	
	if (f_html) {
		BYTE str[BUFSIZE];
		FILE *fp1;
		
		l = file_size(KM_fname);
		fprintf(fp_html, "This report is generated from the file "
			"<a href=\"%s\"> %s </a>. This file is %ld bytes long.\n<p>\n", 
			KM_fname, KM_fname, l);
		system("rm a");
		system("date >a");
		fp1 = fopen("a", "r");
		fgets(str, BUFSIZE, fp1);
		fclose(fp1);
		fprintf(fp_html, "creation date: %s<p>\n", str);
		}
	

	highest_layer = k - 1;
	num_layers = k - 1 - t;
	printf("highest_layer = %ld num_layers = %ld\n", highest_layer, num_layers);
	
	printf("stab_go=\n");
	stab_go.Print();
	((VECTOR_OP) stab_go.s_i(0))->s_i(0)->copy(&go);
	str[0] = 0;
	go.sprint(str);
	// go.println();
	printf("the group order is %s\n", str);
	printf("computing orbit length vector:\n");
	fflush(stdout);
	dc_calc_orbit_length(&go, &stab_go, &orbit_length);
	if (f_html) {
		fprintf(fp_html, "the group %s of order %s is "
			"generated by the following permutations:<p>\n", group_label, str);
		l = G_gen.s_li();
		for (i = 0; i < l; i++) {
			str[0] = 0;
			G_gen.s_i(i)->sprint(str);
			fprintf(fp_html, "%s<br>\n", str);
			}
		fprintf(fp_html, "<p>\n");
		fprintf(fp_html, "jump to the Kramer Mesner matrix "
			"<a href=\"#matrix%ldsets%ldsets\"> M<sub> %ld,%ld</sub> </a>\n<p>\n", t, k, t, k);
		fprintf(fp_html, "jump to "
			"<a href=\"#solutions\"> solution vectors </a>\n<p>\n");
		}


	printf("\n\n\n#if TEXDOCU\n");
	printf("#else\n");
	printf("\\section{$%ld$-$(%ld,%ld,\\lambda)$ designs}\n", t, v, k);

	
	printf("\\begin{verbatim}\n");
	printf("%s\n", KM_fname);
	printf("\\end{verbatim}\n");


	km_report_group(k, &go, &G_gen, FALSE /* f_sims_chain */);
	km_report_orbits(k, &go, &stab_go, &RR);
	if (f_html) {
		km_report_orbits_html(fp_html, group_label, k, &go, &stab_go, &RR);
		}
	printf("%% after report orbits: group_label = %s\n", group_label);
	fflush(stdout);


	printf("The Kramer-Mesner matrix for %ld- and %ld-orbits "
		"(size $%ld \\times %ld$)\n", t, k, Mtk.s_hi(), Mtk.s_li());
#if 0
	printf("{\\arraycolsep0pt\n");
	printf("{\\tiny\n");
	printf("\\[\n");
	Mtk.latex(stdout);
	printf("\\]\n");
	printf("}\n");
	printf("}\n");
#endif
	if (Mtk.s_li() < 40) {
		printf("%% print_nice_matrix\n"); fflush(stdout);
		print_nice_matrix(&Mtk, t, k, 
			&RR, &stab_go, 
			&go, stdout);
		// report_plesken_matrix(k, &MM, &RR, &stab_go, &go);

		}
	print_and_split_matrix(&Mtk, KM_MATRIX_COL_WIDTH_FOR_SPLIT, TRUE /* f_clearpage */);

#if 0
	printf("$A^\\vee$ block-wise:\\\\\n");
	for (j = 1; j <= k; j++) {
		printf("%ld-th non diagonal elements:\\\\\n", j);
		for (i = 0; i <= k; i++) {
			i1 = i + j;
			if (i1 > k)
				continue;
			printf("$A_{%ld,%ld}$:\\\\\n", i, i1);
			Aij = (MATRIX_OP) Ainf_block_wise.s_ij(i, i1);
			printf("\\[\n");
			Aij->latex(stdout);
			printf("\\]\n");
			}
		}
#endif
	if (f_html) {
		km_report_KM_matrices_html(fp_html, group_label, k, &MM);
		km_report_KM_matrix_html(fp_html, group_label, t, k, &Mtk);
		}
	

	if (f_html) {
		fprintf(fp_html, "<p>\n\n<a name=\"solutions\"><h2>Solutions:</h2></a>\n<p>\n\n");
		fflush(stdout);
		}
	printf("%% group_label=%s\n", group_label);

	fp = fopen(KM_fname, "r");
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strncmp(buf, "LAMBDA", 6) != 0)
			continue;
		sscanf(buf, "LAMBDA %ld", &lambda);
		if (f_select) {
			if (lambda != design_select_lambda) {
				read_until_lambdaend(fp);
				continue;
				}
			printf("\n\n\n#if SELECTLAMBDA\n");
			printf("#else\n");
			}


		printf("\\clearpage\n");
		{
		BYTE cmd[1000];
		sprintf(cmd, "t122.out -tex KM_%s_t%ld_k%ld.txt %ld 3 >aa", 
			group_label, t, k, lambda);
		system(cmd);
		}
		printf("\\input report_KM_%s_t%ld_k%ld.txt_%ld.tex\\\\\n", 
			group_label, t, k, lambda);
		printf("\\underline{{\\large $\\lambda = %ld$}}\\\\\n", lambda);
		calc_and_print_design_parameter(t, v, k, lambda);
	
		if (f_html) {
			fprintf(fp_html, "<p>\n\n<a name=\"lambda%ld\"><h2>\\lambda=%ld:</h2></a>\n<p>\n\n", lambda, lambda);
			fflush(fp_html);
			}
		
		Mendelsohn(v, t, k, k, lambda, 
			&Mendelsohn_mtx, &Mendelsohn_RHS, 
			&LAmbda, FALSE /* f_v */);
		s_max = 3;
		Mendel_RHS.m_il(s_max + 1);
		Mendelsohn_RHS.copy((VECTOR_OP) Mendel_RHS.s_i(1));
		for (s = 2; s <= s_max; s++) {
			Mendelsohn_generalized_RHS(v, t, k, k, s, lambda, 
				(VECTOR_OP) Mendel_RHS.s_i(s), FALSE /* f_v */);
			}
		calc_Lambda(v, t, k, lambda, &Lambda);

		// Koehler equations for m1 = k (generalized up to s_max):
		Koehler(v, t, k, lambda, k, s_max, 
			&Coeff, 
			&Constant_term, FALSE /* f_v */);
		
		// Koehler equations for v (generalized up to s_max):
		Koehler(v, t, k, lambda, v, s_max, 
			&Koehler_Coeff_v, 
			&Koehler_Constant_term_v, FALSE /* f_v */);
		
		design_print_Mendelsohn_and_Koehler(v, t, k, &Mendelsohn_mtx, 
			&Mendel_RHS, &LAmbda, &alpha_i, 
			&Coeff, &Constant_term, 
			&Koehler_Constant_term_v, &Koehler_Coeff_v);
	


		printf("solution vectors for $\\lambda = %ld$: \\\\\n", lambda);
		// printf("found solutions, checking...\n");
		bai_len.m_i(0);
		
		design_clear_base_blocks_file(KM_fname, lambda);
		no = 0;
		while (TRUE) {
			if (fgets(buf, BUFSIZE, fp) == NULL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0;
			while (l >= 0) {
				if (buf[l - 1] == ' ') {
					buf[l - 1] = 0;
					l--;
					}
				else 
					break;
				}
			if (strncmp(buf, "LAMBDAEND", 9) == 0) {
				if (f_select) {
					printf("#endif\n");
					}
				break;
				}
			// printf("checking solution %ld: %d\n", no, buf);
			l = strlen(buf);
			if (l != n) {
				printf("l = %ld != n = %ld\n", l, n);
				continue;
				}
			no++;
			printf("%% ${\\cal D}_{%ld}$: \n", no);
			if (f_select) {
				if (no < select_first)
					continue;
				if (no >= select_first + select_length) {
					read_until_lambdaend(fp);
					break;
					}
				}

#if 0
			if (no > BREAK_AFTER_N_SOLUTIONS) {
				read_until_lambdaend(fp);
				break;
				}
#endif
			offset = "0pt";
			
			design_print_solution_vector(buf, no, SOLUTION_VECTOR_WIDTH_FOR_SPLIT, offset);
			design_print_solution_vector_numerical(buf, no, offset);

			if (f_html) {
				fprintf(fp_html, "%s<br>\n", buf);
				}
		
			
			X.m_ilih(1, n);
			for (i = 0; i < n; i++) {
				if (buf[i] == '1')
					X.m_iji(i, 0, 1);
				else
					X.m_iji(i, 0, 0);
				}
			design_orbits(&X, &orbits, FALSE /* f_complement */ );
			if (F_USE_COMPLEMENTS) {
				design_orbits(&X, &orbits_c, TRUE /* f_complement */ );
				}
			Mtk.mult(&X, &Y);
			for (i = 0; i < m; i++) {
				a = Y.s_iji(i, 0);
				if (a != lambda)
					printf("WARNING !!! Y[%ld] = %ld != lambda = %ld\n", i, a, lambda);
				}
			if (0 /* f_output_geo && 0 */) {
				printf("$");
				design_output_geo(&G_gen, k, &RR, &X);
				printf("$");
				}
			if (1) {
				design_print_base_blocks(&G_gen, &RR, k, &X, 
					&go, &stab_go);
				}
			design_print_base_blocks_file(KM_fname, lambda, no,
				&RR, k, &X);
			
			}
		printf("found %ld designs \n", no);
		} // end while
	fclose(fp);
	printf("#endif\n\n");
	if (f_html) {
		fprintf(fp_html, "</body>\n");
		fprintf(fp_html, "</html>\n");
		fclose(fp_html);
		}

}

void KM_fname_extract_group_label(BYTE *KM_fname, BYTE *group_label)
{
	BYTE str[1000];
	INT l, i;
	
	strcpy(str, KM_fname);
	l = strlen(str);
	if (l > 4 && strcmp(str + l - 4, ".txt") == 0)
		str[l - 4] = 0;
	group_label[0] = 0;
	l = strlen(str);
	if (l > 3) {
		strcpy(group_label, str + 3);
		}
	else {
		strcpy(group_label, str);
		}
	l = strlen(group_label);
	for (i = l - 1; i >= 0; i--) {
		if (group_label[i] == '_') {
			group_label[i] = 0;
			break;
			}
		}
	for (i--; i >= 0; i--) {
		if (group_label[i] == '_') {
			group_label[i] = 0;
			break;
			}
		}
}

INT report_plesken_matrix(INT k, VECTOR_OP MM, VECTOR_OP RR, VECTOR_OP stab_go, SYM_OP go)
{
	MATRIX_OB Ainf, Ainf_inv, Ainf_block_wise;
	VECTOR_OB K_first, K_len;
	INT f_v = TRUE;
	INT f_vv = FALSE;
	
	// dc_calc_Ainf_t_via_MM(k, &Ainf, &MM, FALSE /* f_verbose */);
	dc_plesken_matrices_prepare(0 /* k_min */, k /* k_max */, MM, 
		&Ainf_block_wise, &Ainf, &Ainf_inv, 
		&K_first, &K_len, f_v, f_vv);
	
	printf("%% print_nice_plesken_matrix\n"); fflush(stdout);
	print_nice_plesken_matrix(&Ainf_block_wise, 
		0 /* k_min */, k /* k_max */, 
		RR, stab_go, 
		go, stdout);
	printf("%% finished\n"); fflush(stdout);
	return OK;
}

#if TEXDOCU
void design_report_plesken(BYTE *KM_fname)
#endif
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n, i, j, i1, l, h;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	VECTOR_OB orbit_length, K, Eqns;
	MATRIX_OB Ainf, Ainf_inv, Ainf_block_wise;
	VECTOR_OB K_first, K_len;
	INT highest_layer, num_layers;
	SYM_OB go;
	MATRIX_OP Aij;
	INT f_v = TRUE;
	INT f_vv = FALSE;

	printf("KM-file: %s\n", KM_fname);
	fflush(stdout);
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	printf("v=%ld t=%ld k=%ld\n", v, t, k);
	fflush(stdout);
	
	// dc_calc_Ainf_t_via_MM(k, &Ainf, &MM, FALSE /* f_verbose */);
	dc_plesken_matrices_prepare(0 /* k_min */, k /* k_max */, &MM, 
		&Ainf_block_wise, &Ainf, &Ainf_inv, 
		&K_first, &K_len, f_v, f_vv);
	
	highest_layer = k - 1;
	num_layers = k - 1 - t;
	printf("highest_layer = %ld num_layers = %ld\n", highest_layer, num_layers);
	
	stab_go.s_i(0)->copy(&go);
	printf("the group order is ");
	go.println();
	printf("computing orbit length vector:\n");
	fflush(stdout);
	dc_calc_orbit_length(&go, &stab_go, &orbit_length);


	printf("\n\n\n#if TEXDOCU\n");
	printf("#else\n");
	printf("\\section{$%ld$-$(%ld,%ld,\\lambda)$ designs}\n", t, v, k);

	
	printf("\\begin{verbatim}\n");
	printf("%s\n", KM_fname);
	printf("\\end{verbatim}\n");


	stab_go.s_i(0)->copy(&go);
	km_report_group(k, &go, &G_gen, FALSE /* f_sims_chain */);
	km_report_orbits(k, &go, &stab_go, &RR);

	printf("$A^\\vee$ block-wise:\\\\\n");
	for (j = 1; j <= k; j++) {
		printf("%ld-th non diagonal elements:\\\\\n", j);
		for (i = 0; i <= k; i++) {
			i1 = i + j;
			if (i1 > k)
				continue;
			printf("$A_{%ld,%ld}$:\\\\\n", i, i1);
			Aij = (MATRIX_OP) Ainf_block_wise.s_ij(i, i1);
			l = Aij->s_li();
			h = Aij->s_hi();
			if (l > KM_MATRIX_COL_WIDTH_FOR_SPLIT) {
				print_and_split_matrix(Aij, 
					KM_MATRIX_COL_WIDTH_FOR_SPLIT, 
					TRUE /* f_clearpage */);
				}
			else {
				printf("{\\arraycolsep1pt\n");
				printf("{\\tiny\n");
				printf("\\[\n");
				Aij->latex(stdout);
				printf("\\]\n");
				printf("}%%\n");
				printf("}%%\n");
				}
			}
		}
	
	return ;
}

#if TEXDOCU
static INT read_until_lambdaend(FILE *fp)
#endif
{
	BYTE buf[BUFSIZE];
	INT l;
	
	// search for LAMBDAEND:
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strncmp(buf, "LAMBDAEND", 9) == 0)
			break;
		}
	return OK;
}


#if TEXDOCU
static INT print_and_split_matrix(MATRIX_OP M, INT width_in_cols, INT f_clearpage)
#endif
{
	MATRIX_OP eqn;
	VECTOR_OB S;
	INT j, l;
	
	M->split_column_wise(&S, width_in_cols);
	l = S.s_li();
	for (j = 0; j < l; j++) {
		eqn = (MATRIX_OP) S.s_i(j);
		printf("{\\arraycolsep1pt\n");
		printf("{\\tiny\n");
		printf("\\[\n");
		eqn->latex(stdout);
		printf("\\]\n");
		printf("}%%\n");
		printf("}%%\n");
		if (f_clearpage) {
			;
			}
		else {
			printf("%%");
			}
		printf("\\clearpage\n");
		}
	return OK;
}

#if TEXDOCU
static INT km_report_group(INT k, SYM_OP go, VECTOR_OP G_gen, INT f_sims_chain)
#else
#endif
{
	BYTE str[1024];
	INT i, l;
	SYM_OB goG;
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
	
	str[0] = 0;
	go->sprint(str);
	printf("$|G| = %s$\\\\\n", str);
	
	if (f_sims_chain) {
		printf("{Sims}-chain:\n");
		printf("\\[\n");
	
		reduce_generators_labra(G_gen, &goG, FALSE /* f_verbose */, &labG);
		labG.calc_transversal_matrix(&TG, FALSE /* f_v */);
		// TG.Print();
		TG.latex_upper_tri(stdout);
	
		printf("\\]\n");

		printf("the transversal elements:\\\\\n");
		labG.latex_transversal_indices(stdout);
		}
	
	
	return OK;
}

#if TEXDOCU
static INT km_report_orbits(INT k, SYM_OP go, 
	VECTOR_OP stab_go, VECTOR_OP RR)
#else
//PRE
input: k, go, stab\_go, RR
computes K from RR (RR are the set representatives for all $i \le k$)
///PRE
#endif
{
	BYTE str[1024];
	BYTE str1[1024];
	INT i, j, l, jj, ll, a;
	VECTOR_OP p_stab_go, pR, ppR;
	SYM_OB ol, goG;
	
	printf("now the orbit representatives, format is: \\\\\n");
	printf("local orbit number (global orbit number): "
		"$\\{$ set representative $\\}$ stabilizer order\\\\\n");
	
	printf("\\begin{multicols}{2}\n");
	
	for (i = 0; i <= k; i++) {
		p_stab_go = (VECTOR_OP) stab_go->s_i(i);
		pR = (VECTOR_OP) RR->s_i(i);
		l = pR->s_li();
		printf("\\underline{{\\bf %ld-orbits:}} \\\\\n", i);
		for (j = 0; j < l; j++) {
			ppR = (VECTOR_OP) pR->s_i(j);
			ll = ppR->s_li();
			printf("%ld: $\\{", j + 1);
			for (jj = 0; jj < ll; jj++) {
				a = ppR->s_ii(jj);
				printf("%ld", a);
				if (jj < ll - 1)
					printf(", ");
				}
			str[0] = 0;
			str1[0] = 0;
			go->ganzdiv(p_stab_go->s_i(j), &ol);
			p_stab_go->s_i(j)->sprint(str);
			// ol.sprint(str1);
			// printf("\\}_{%s,%s}$\\\\\n", str1, str);
			printf("\\}_{%s}$\\\\\n", str);
			}
		}
	printf("\\end{multicols}\n");
	return OK;
}

#if TEXDOCU
static INT km_report_orbits_html(FILE *fp, BYTE *group_label, 
	INT k, SYM_OP go, VECTOR_OP stab_go, VECTOR_OP RR)
#else
//PRE
input: k, go, stab\_go, RR
computes K from RR (RR are the set representatives for all $i \le k$)
///PRE
#endif
{
	BYTE str[1024];
	BYTE str1[1024];
	INT i, j, jj, l, ll, a, nb_d;
	VECTOR_OP p_stab_go, pR, ppR;
	SYM_OB ol, goG;
	
	fprintf(fp, "<table border=\"1\">\n");
	fprintf(fp, "<caption>\n");
	fprintf(fp, "<center>\n");
	fprintf(fp, "orbits of %s on i-sets, i less than or equal to %ld\n", group_label, k);
	fprintf(fp, "</center>\n");
	fprintf(fp, "</caption>\n");
	fprintf(fp, "<tr><th> i </th> <th> # of orbits </th>  <th> jump </th> \n");
	nb_d = 0;
	for (i = 0; i <= k; i++) {
		p_stab_go = (VECTOR_OP) stab_go->s_i(i);
		pR = (VECTOR_OP) RR->s_i(i);
		l = pR->s_li();
		fprintf(fp, "<tr><td> %ld </td> <td> %ld </td> <td> "
			"jump to <a href=\"#orbits%ldsets\"> orbits </a> ", i, l, i);
		if (i < k) {
			fprintf(fp, "/  <a href=\"#matrix%ldsets%ldsets\"> "
				"KM-matrix M<sub>%ld,%ld</sub> </a> </td> \n", i, i + 1, i, i + 1);
			}
		else {
			fprintf(fp, "</td>\n");
			}
		nb_d += l;
		}
	fprintf(fp, "</table>\n");
	fprintf(fp, "altogether %ld orbits\n\n<p><p>\n\n", nb_d);
	
	for (i = 0; i <= k; i++) {
		p_stab_go = (VECTOR_OP) stab_go->s_i(i);
		pR = (VECTOR_OP) RR->s_i(i);
		l = pR->s_li();
		fprintf(fp, "<a name=\"orbits%ldsets\"> </a>\n", i);
		fprintf(fp, "<table border=\"1\">\n");
		fprintf(fp, "<caption>\n");
		fprintf(fp, "<center>\n");
		fprintf(fp, "orbit representatives of %s on %ld-sets\n", group_label, i);
		fprintf(fp, "</center>\n");
		fprintf(fp, "</caption>\n");
		fprintf(fp, "<tr><th> # </th> <th> representative </th> <th> order of the set-stabilizer </th> \n");
		for (j = 0; j < l; j++) {
			ppR = (VECTOR_OP) pR->s_i(j);
			ll = ppR->s_li();
			fprintf(fp, "<tr> <td> %ld </td> <td> {", j + 1);
			for (jj = 0; jj < ll; jj++) {
				a = ppR->s_ii(jj);
				fprintf(fp, "%ld", a);
				if (jj < ll - 1)
					fprintf(fp, ", ");
				}
			str[0] = 0;
			str1[0] = 0;
			go->ganzdiv(p_stab_go->s_i(j), &ol);
			p_stab_go->s_i(j)->sprint(str);
			// ol.sprint(str1);
			// printf("\\}_{%s,%s}$\\\\\n", str1, str);
			fprintf(fp, "} </td> <td> %s </td>\n", str);
			}
		fprintf(fp, "</table>\n\n<p>\n\n");
		}
	return OK;
}

#if TEXDOCU
static INT km_report_KM_matrices_html(FILE *fp, BYTE *group_label, INT k, VECTOR_OP MM)
#endif
{
	INT i, l, h;
	MATRIX_OP pM;
	
	for (i = 0; i < k; i++) {
		pM = (MATRIX_OP) MM->s_i(i);
		l = pM->s_li();
		h = pM->s_hi();
		fprintf(fp, "<a name=\"matrix%ldsets%ldsets\"> </a>\n", i, i + 1);
		fprintf(fp, "<h2>\n");
		fprintf(fp, "Kramer Mesner matrix M <sub>%ld,%ld</sub> "
			"for %s of size %ld x %ld\n", i, i + 1, group_label, h, l);
		fprintf(fp, "</h2><p>\n");
		fprintf(fp, "<pre>\n");
		pM->fPrint(fp);
		fprintf(fp, "</pre>\n\n<p>\n\n");
		}
	return OK;
}

#if TEXDOCU
static INT km_report_KM_matrix_html(FILE *fp, BYTE *group_label, INT t, INT k, MATRIX_OP Mtk)
#endif
{
	INT l, h;
	
	if (k != t + 1) {
		l = Mtk->s_li();
		h = Mtk->s_hi();
		fprintf(fp, "<a name=\"matrix%ldsets%ldsets\"> </a>\n", t, k);
		fprintf(fp, "<h2>\n");
		fprintf(fp, "Kramer Mesner matrix M <sub>%ld,%ld</sub> "
			"for %s of size %ld x %ld\n", t, k, group_label, h, l);
		fprintf(fp, "</h2><p>\n");
		fprintf(fp, "<pre>\n");
		Mtk->fPrint(fp);
		fprintf(fp, "</pre>\n\n<p>\n\n");
		fflush(fp);
		}
	return OK;
}

static INT design_output_geo(VECTOR_OP gen, INT k, VECTOR_OP RR, MATRIX_OP X)
{
	VECTOR_OB O, O_first, O_len;
	VECTOR_OP pO;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	INT v, b, *theX, nb_X;
	PERMUTATION_OP p;
	PERMUTATION_OB pp, qq;
	MATRIX_OB I;
	INT i, j, ii;
	LABRA_OB aut;
	SYM_OB ago;
	VECTOR_OB aut_gen;
	
	design_calc_blocks(gen, RR, k, X, &O, &O_first, &O_len, f_v, f_vv);
	b = O.s_li();
	p = (PERMUTATION_OP) gen->s_i(0);
	v = p->s_li();
	nb_X = b * k;
	theX = (INT *) my_malloc(nb_X * sizeof(INT), "theX");
	I.m_ilih_n(b, v);
	for (j = 0; j < b; j++) {
		pO = (VECTOR_OP) O.s_i(j);
		if (pO->s_li() != k)
			return error("design_output_geo() pO->s_li() != k");
		for (ii = 0; ii < k; ii++) {
			i = pO->s_ii(ii);
			I.m_iji(i, j, 1);
			}
		}
	printf("I:\\\\\n");
	I.latex(stdout);
	ii = 0;
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			if (I.s_iji(i, j))
				theX[ii++] = i * b + j;
			}
		}
	if (ii != nb_X)
		return error("design_output_geo() ii != nb_X");
	
	printf("%% ");
	for (ii = 0; ii < nb_X; ii++) {
		printf("%ld ", theX[ii]);
		}
	printf(" %% GEO %ld %ld %ld\n", v, b, nb_X);
	printf("%% ");

#if 0
#ifdef SYM_GEO
	geo_canon_simple(FALSE, &back_to, v, b, nb_X, TRUE /* f_print_dots */, 
		theX, &pp, &qq, FALSE /* f_transposed */, 
		TRUE /* f_get_aut_group */, &aut, 
		FALSE /* f_canon_v */,  FALSE /* f_canon_vv */);
#else
	return error("design_output_geo() GEOLIB not available !");
#endif
	// changes theX !
	printf("\n");

	
	printf("%% ");
	for (ii = 0; ii < nb_X; ii++) {
		printf("%ld ", theX[ii]);
		}
	printf(" %% G_CANON %ld %ld %ld\n", v, b, nb_X);

	printf("\n");
	printf("automorphism group order ");
	aut.group_order(&ago);
        ago.println();
	// aut.save("aut_higman_design.dsc");
	// aut.generators(&aut_gen, &ago);
	// aut.reduced_generating_set(&aut_gen, FALSE /* f_bottom_up */, TRUE /* f_v */);
	// aut_gen.save("aut_higman_design_generators.dsc");
	// write_file_of_generators(&aut_gen, "aut_higman_design_generators.txt");
	// fflush(stdout);
	// aut_gen.swap(HS_gen);
	aut.strong_generating_set(&aut_gen, TRUE /* f_v */);
	printf("strong generators: \\\\\n");
	aut_gen.Latex(stdout);
#endif

	my_free(theX);
	return OK;
}

#if TEXDOCU
static INT design_print_base_blocks(VECTOR_OP gen, 
	VECTOR_OP RR, INT k, MATRIX_OP X, SYM_OP go, VECTOR_OP stab_go)
#else
#endif
{
	SYM_OB ol;
	VECTOR_OP R, R1, p_stab_go;
	INT i, jj, l, a;
	BYTE str[1024];
	BYTE str1[1024];
	
	R = (VECTOR_OP) RR->s_i(k);
	p_stab_go = (VECTOR_OP) stab_go->s_i(k);
	l = R->s_li(); // number of k-orbits 
	if (l != X->s_hi())
		return error("design_print_base_blocks() l != X->s_hi()");
	
	printf("\\begin{multicols}{2}\n");
	printf("\\noindent\n");
	
	for (i = 0; i < l; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		R1 = (VECTOR_OP) R->s_i(i);
		printf("$\\{");
		for (jj = 0; jj < k; jj++) {
			a = R1->s_ii(jj);
			printf("%ld", a);
			if (jj < k - 1)
				printf(", ");
			}
		str[0] = 0;
		str1[0] = 0;
		go->ganzdiv(p_stab_go->s_i(i), &ol);
		p_stab_go->s_i(i)->sprint(str);
		ol.sprint(str1);
		printf("\\}_{%s,%s}$\\\\\n", str1, str);
		}
	
	printf("\\end{multicols}\n");
	
	printf("base blocks:\n");
	for (i = 0; i < l; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		R1 = (VECTOR_OP) R->s_i(i);
		for (jj = 0; jj < k; jj++) {
			a = R1->s_ii(jj);
			printf("%ld ", a);
			}
		printf("\n");
		}
	return OK;
}

static INT base_block_fname(BYTE *KM_name, INT lambda, 
	BYTE *base_fname, BYTE *fname)
{
	BYTE *p;
	
	strcpy(fname, KM_name);
	if ((p = strrchr(fname, '.')) != NULL) {
		*p = 0;
		}
	strcpy(base_fname, fname);
	sprintf(fname + strlen(fname), "_%ld.base_blocks", lambda);
	return OK;
}

static INT design_clear_base_blocks_file(BYTE *KM_name, INT lambda)
{
	BYTE base_fname[1000];
	BYTE fname[1000];
	BYTE cmd[1000];
	
	base_block_fname(KM_name, lambda, base_fname, fname);
	sprintf(cmd, "rm %s", fname);
	printf("%%");
	call_system(cmd);
	return OK;
}

#if TEXDOCU
static INT design_print_base_blocks_file(BYTE *KM_name, INT lambda, INT no,
	VECTOR_OP RR, INT k, MATRIX_OP X)
#else
#endif
{
	SYM_OB ol;
	VECTOR_OP R, R1, p_stab_go;
	INT i, jj, l, a;
	BYTE str[1024];
	BYTE str1[1024];
	BYTE base_fname[1000];
	BYTE fname[1000];
	FILE *fp;
	
	
	base_block_fname(KM_name, lambda, base_fname, fname);
	fp = fopen(fname, "a");
	
	fprintf(fp, "GEOMETRY %ld %s_%ld\n", no, base_fname, lambda);
	fprintf(fp, "k=%ld\n", k);
	R = (VECTOR_OP) RR->s_i(k);
	for (i = 0; i < l; i++) {
		if (X->s_iji(i, 0) == 0)
			continue;
		R1 = (VECTOR_OP) R->s_i(i);
		for (jj = 0; jj < k; jj++) {
			a = R1->s_ii(jj);
			fprintf(fp, "%ld ", a);
			}
		fprintf(fp, "\n");
		}
	fprintf(fp, "END\n\n");
	fclose(fp);
	
	return OK;
}

static INT print_nice_plesken_matrix(MATRIX_OP Ainf_block_wise, 
	INT k_min, INT k_max, 
	VECTOR_OP RR, VECTOR_OP stab_go, 
	SYM_OP go, FILE *fp)
{
	MATRIX_OB M1;
	MATRIX_OP M;
	VECTOR_OP R, R1, R2, p_stab_go;
	INT m, n, i, ii, j, I, J, a;
	BYTE str[1024];
	
	fprintf(fp, "\\[\n");
	fprintf(fp, "\\begin{array}{r|");
	for (J = k_min; J <= k_max; J++) {
		R = (VECTOR_OP) RR->s_i(J);
		n = R->s_li();
		fprintf(fp, "*{%ld}{r}", n);
		if (J < k_max)
			fprintf(fp, "|");
		}
	fprintf(fp, "|}\n");
	
	fprintf(fp, " & \n");	
	for (J = k_min; J <= k_max; J++) {
		R = (VECTOR_OP) RR->s_i(J);
		p_stab_go = (VECTOR_OP) stab_go->s_i(J);
		n = R->s_li();
		for (j = 0; j < n; j++) {
			R1 = (VECTOR_OP) R->s_i(j);
			str[0] = 0;
			prepare_set_tex(R1, go, 
				p_stab_go->s_i(j), str);
			fprintf(fp, "\\begin{rotate}{90}$%s$\\end{rotate}\n", str);	

			if (j < n - 1)
				fprintf(fp, " & ");
			}
		if (J < k_max)
			fprintf(fp, " & ");
		else {
			fprintf(fp, "\\\\\n");
			fprintf(fp, "\\hline\n");
			}
		}
	
	
	for (I = k_min; I <= k_max; I++) {
		R = (VECTOR_OP) RR->s_i(I);
		p_stab_go = (VECTOR_OP) stab_go->s_i(I);
		m = R->s_li();
		for (i = 0; i < m; i++) {
			R1 = (VECTOR_OP) R->s_i(i);
			str[0] = 0;
			prepare_set_tex(R1, go, 
				p_stab_go->s_i(i), str);
			fprintf(fp, "%s & ", str);	

			for (J = k_min; J <= k_max; J++) {
				R2 = (VECTOR_OP) RR->s_i(J);
				n = R2->s_li();
				if (J > I)
					M = (MATRIX_OP) Ainf_block_wise->s_ij(I, J);
				else {
					M1.m_ilih_n(n, m);
					if (J == I) {
						for (ii = 0; ii < n; ii++)
							M1.m_iji(ii, ii, 1);
						}
					M = &M1;
					}
				for (j = 0; j < n; j++) {
					a = M->s_iji(i, j);
					fprintf(fp, "%ld ", a);	
					if (j < n - 1)
						fprintf(fp, " & ");
					}
				if (J < k_max)
					fprintf(fp, " & ");
				else {
					fprintf(fp, "\\\\\n");
					}
				
				} // next J
			} // next i
		fprintf(fp, "\\hline\n");
		} // next I

	fprintf(fp, "\\end{array}\n");
	fprintf(fp, "\\]\n");
	return OK;
}

static INT print_nice_matrix(MATRIX_OP M, INT t, INT k, 
	VECTOR_OP RR, VECTOR_OP stab_go, 
	SYM_OP go, FILE *fp)
{
	VECTOR_OP R, R1, p1, p2;
	INT m, n, i, j;
	BYTE str[1024];
	
	m = M->s_hi();
	n = M->s_li();
	p1 = (VECTOR_OP) stab_go->s_i(t);
	p2 = (VECTOR_OP) stab_go->s_i(k);
	fprintf(fp, "\\[\n");
	fprintf(fp, "\\begin{array}{*{%ld}{r}}\n", n + 1);
	fprintf(fp, " & ");
	R = (VECTOR_OP) RR->s_i(k);
	for (j = 0; j < n; j++) {
		R1 = (VECTOR_OP) R->s_i(j);
		str[0] = 0;
		prepare_set_tex(R1, go, p2->s_i(j), str);
		fprintf(fp, "\\begin{rotate}{90}$%s$\\end{rotate}\n", str);

		if (j < n - 1)
			fprintf(fp, " & ");
		else
			fprintf(fp, "\\\\\n");
		}
	R = (VECTOR_OP) RR->s_i(t);
	for (i = 0; i < m; i++) {
		R1 = (VECTOR_OP) R->s_i(i);
		str[0] = 0;
		prepare_set_tex(R1, go, p1->s_i(i), str);
		fprintf(fp, "%s & ", str);

		for (j = 0; j < n; j++) {
			fprintf(fp, "%ld ", M->s_iji(i, j));
			if (j < n - 1)
				fprintf(fp, " & ");
			else
				fprintf(fp, "\\\\\n");
			}
		}
	fprintf(fp, "\\end{array}\n");
	fprintf(fp, "\\]\n");
	return OK;
}

static INT prepare_set_tex(VECTOR_OP R, SYM_OP go, SYM_OP stab_go, BYTE *s)
{
	SYM_OB ol;
	BYTE str[1024];
	BYTE str1[1024];
	BYTE str2[1024];
	INT jj, a, l;
	
	str[0] = 0;
	sprintf(str + strlen(str), "\\{");
	l = R->s_li();
	for (jj = 0; jj < l; jj++) {
		a = R->s_ii(jj);
		sprintf(str + strlen(str), "%ld", a);
		if (jj < l - 1)
			sprintf(str + strlen(str), ", ");
		}
	str1[0] = 0;
	str2[0] = 0;
	go->ganzdiv(stab_go, &ol);
	stab_go->sprint(str1);
	ol.sprint(str2);
	sprintf(str + strlen(str), "\\}_{%s,%s}", str1, str2);
	strcat(s, str);
	return OK;
}

#if TEXDOCU
INT do_report(BYTE *KM_fname, INT f_html, INT f_selected, 
	INT select_lambda, INT select_first, INT select_len)
#else
This is the main entry point from the user interface. 
This routine starts the report and produce and translates 
and displays the latex report immediately on the creen using 
unix exec commands. 
The following extenal calls are made:
//PRE
discreta\_report [options] KM\_fname
grepdocu [options] -no\_filename -no\_underscore\_translation report\_fname
latex r
xdvi r.dvi
///PRE
Additionally, the file "all\_reports.tex" is extended by one line 
containing the latex input command.
#endif
{
	BYTE report_fname[1024];
	BYTE cmd[1024];
	BYTE options[1024];
	BYTE grep_options[1024];
	FILE *fp;
	
	sprintf(report_fname, "%s_report", KM_fname);
	options[0] = 0;
	grep_options[0] = 0;
	if (f_html) {
		sprintf(options + strlen(options), "-html ");
		}
	if (f_selected) {
		sprintf(options + strlen(options), "-select %ld %ld %ld ", select_lambda, select_first, select_len);
		sprintf(grep_options + strlen(grep_options), "-string SELECTLAMBDA ");
		}
	sprintf(cmd, "discreta_report %s %s | tee %s", options, KM_fname, report_fname);
	call_system(cmd);
	
	fp = fopen(report_fname, "a");
	fprintf(fp, "#endif\n\n");
	fclose(fp);
	printf("written file %s of size %ld\n", 
		report_fname, file_size(report_fname));
	fflush(stdout);

	
	sprintf(cmd, "grepdocu %s -no_filename -no_underscore_translation %s", 
		grep_options, report_fname);
	call_system(cmd);
	strcat(report_fname, ".tex");
	printf("written file %s of size %ld\n", 
		report_fname, file_size(report_fname));
	fflush(stdout);

	
	call_system("cp r0.tex r.tex");
	fp = fopen("r.tex", "a");
	fprintf(fp, "\\input %s\n\n", report_fname);
	fprintf(fp, "\\end{document}\n\n");
	fclose(fp);
	
	fp = fopen("all_reports.tex", "a");
	fprintf(fp, "\\input %s\n\n", report_fname);
	fclose(fp);
	
	sprintf(cmd, "latex r.tex");
	call_system(cmd);
	
	sprintf(cmd, "xdvi r.dvi");
	call_system(cmd);
	
	return OK;
}

static INT write_r0_tex_file()
{
	FILE *fp;
	
	fp = fopen("r0.tex", "w");
	fprintf(fp, "\n");
fprintf(fp, "\\documentclass[]{article}\n");
fprintf(fp, "%%12pt,titlepage\n");
fprintf(fp, "\\usepackage{fancyheadings}\n");
fprintf(fp, "%%\\usepackage{amstex}\n");
fprintf(fp, "\\usepackage{amsmath}\n");
fprintf(fp, "\\usepackage{amssymb}\n");
fprintf(fp, "\\usepackage{multicol}\n");
fprintf(fp, "\\usepackage{epsf}\n");
fprintf(fp, "\\usepackage{supertabular}\n");
fprintf(fp, "\\usepackage{wrapfig}\n");
fprintf(fp, "%%\\usepackage{blackbrd}\n");
fprintf(fp, "\\usepackage{epic,eepic}\n");
fprintf(fp, "\\usepackage{rotating}\n");
fprintf(fp, "%%\\usepackage{concmath}\n");
fprintf(fp, "%%\\usepackage{beton}\n");
fprintf(fp, "%%\\documentstyle[12pt,epsf]{article}\n");
fprintf(fp, "%%\\pagestyle{headings}\n");
fprintf(fp, "%%\\pagestyle{empty}\n");
fprintf(fp, "%%\\thispagestyle{empty}\n");
fprintf(fp, "%%\\oddsidemargin=-2cm\n");
fprintf(fp, "%%\\evensidemargin=0pt\n");
fprintf(fp, "\\newfont{\\mytt}{cmtt10 scaled 700}\n");
fprintf(fp, "%%\\input dina4.sty\n");
fprintf(fp, "%%This is documentsubstyle DINA4 for DIN A4 pagesize. GMD Z1.BN  12.06.85\n");
fprintf(fp, "\\evensidemargin 0in\n");
fprintf(fp, "\\oddsidemargin 0in\n");
fprintf(fp, "\\marginparwidth 0pt\n");
fprintf(fp, "\\marginparsep 0pt\n");
fprintf(fp, "\\topmargin -1in\n");
fprintf(fp, "\\headheight 0.7cm\n");
fprintf(fp, "\\headsep 1.8cm\n");
fprintf(fp, "%%\\footheight 0.7cm\n");
fprintf(fp, "\\footskip 2cm\n");
fprintf(fp, "\\textheight 22cm\n");
fprintf(fp, "\\textwidth 6.2in\n");
fprintf(fp, "\\marginparpush 0pt\n");
fprintf(fp, "%%\\title{}\n");
fprintf(fp, "%%\\author{{\\sc Anton Betten}}\n");
fprintf(fp, "%%\\date{}\n");
fprintf(fp, "\\begin{document}\n");
fprintf(fp, "%%\\maketitle\n");
	fclose(fp);
	return OK;
}

#endif /* LADDER_TRUE */

