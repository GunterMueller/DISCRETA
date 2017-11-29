//cb.C

#if TEXDOCU
void calc_lambda_i(INT t, INT v, INT k, INT lambda)
#else
clips the output directly into the text window via out\_str.
Apart from that, the function is nearly equivalent to 
calc\_delta\_lambda() in parameter.C.
#endif
{
	BYTE str[1024];
	INT i, a, b, c, g, rhs_a, rhs_b, delta_lambda, dl;

	if (v == 0) {
		printf("v not set !\n");
		return;
		}
	if (t == 0) {
		printf("t not set !\n");
		return;
		}
	if (k == 0) {
		printf("k not set !\n");
		return;
		}
	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	sprintf(str, "v=%ld t=%ld k=%ld lambda=%ld\n", v, t, k, lambda);
	out_str(str);
	fflush(stdout);
	for (i = t; i >= 0; i--) {
		if (i == t) {
			rhs_a = lambda;
			rhs_b = 1;
			delta_lambda = 1;
			}
		else {
			a = rhs_a * (v - i);
			b = rhs_b * (k - i);
			ggt_iipi(a, b, &g);
			a /= g;
			b /= g;
			kgv_iipi(delta_lambda, b, &dl);
			delta_lambda = dl;
			sprintf(str, "t'=%ld lambda'=%ld/%ld delta_lambda=%ld\n", 
				i, a, b, delta_lambda);
			out_str(str);
			fflush(stdout);
			rhs_a = a;
			rhs_b = b;
			}
		}
}

INT km_fname_available()
{
	if (file_size("km_fname") > 0)
		return TRUE;
	else
		return FALSE;
}

INT read_km_fname(BYTE *KM_fname)
{
	FILE *fp;
	BYTE buf[BUFSIZE];
	INT l;
	
	fp = fopen("km_fname", "r");
	if (fgets(buf, BUFSIZE, fp) == NIL) {
		fclose(fp);
		return error("error reading KM-fname");
		}
	fclose(fp);
#if 0
	l = strlen(buf);
	if (l > 0)
		buf[l - 1] = 0;
#endif
	strcpy(KM_fname, buf);
	text_field_set_string(Text_file, KM_fname);
	return OK;
}

INT solid_fname_available()
{
	if (file_size("solid_name") > 0)
		return TRUE;
	else
		return FALSE;
}

INT read_solid_fname(BYTE *fname)
{
	FILE *fp;
	BYTE buf[BUFSIZE];
	INT l;
	
	fp = fopen("solid_name", "r");
	if (fgets(buf, BUFSIZE, fp) == NIL) {
		fclose(fp);
		return error("error reading solid_name");
		}
	fclose(fp);
	l = strlen(buf);
	if (l > 0 && buf[l - 1] == '\n')
		buf[l - 1] = 0;
	strcpy(fname, buf);
	// text_field_set_string(Text_file, KM_fname);
	return OK;
}


void do_mendelsohn(INT t, INT v, INT k, INT lambda)
{
	BYTE str[1024];
	MATRIX_OB M;
	MATRIX_OB LAmbda;
	VECTOR_OB RHS;
	INT s_max = 3;

	if (v == 0) {
		printf("v not set !\n");
		return;
		}
	if (t == 0) {
		printf("t not set !\n");
		return;
		}
	if (k == 0) {
		printf("k not set !\n");
		return;
		}
	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	sprintf(str, "v=%ld t=%ld k=%ld lambda=%ld\n", v, t, k, lambda);
	out_str(str);
	fflush(stdout);
	Mendelsohn(v, t, k, s_max, lambda, &M, &RHS, &LAmbda, TRUE /* f_v */);
	
}

void do_koehler(INT t, INT v, INT k, INT lambda, INT m)
{
	BYTE str[1024];
	// VECTOR_OB RHS;
	VECTOR_OB Coeff, Constant_term;
	INT s_max = 3;

	if (v == 0) {
		printf("v not set !\n");
		return;
		}
	if (t == 0) {
		printf("t not set !\n");
		return;
		}
	if (k == 0) {
		printf("k not set !\n");
		return;
		}
	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	if (m == 0) {
		printf("m not set !\n");
		return;
		}
	sprintf(str, "v=%ld t=%ld k=%ld lambda=%ld m=%ld\n", v, t, k, lambda, m);
	out_str(str);
	fflush(stdout);
	
	// calc_Lambda(v, t, k, lambda, &RHS);
	Koehler(v, t, k, lambda, m, s_max, 
		&Coeff, &Constant_term, TRUE /* f_v */);
}

#if 0
void km_preproc_cb(Widget w, void *w1, void *call_data)
{
	KM_DESIGN_VALUES dv;
	INT v, t, k, lambda, m, n;
	MATRIX_OB Mtk, KM, Ainf, I;
	VECTOR_OB G_gen, RR, MM, go, stab_go, K;
	VECTOR_OB RHS;
	LABRA_OB labraG;
	MATRIX_OB TG;
	INT i, j, l, ii, nb_d, nb_d1, deg, a, i1, i2;
	VECTOR_OP pR;
	SYM_OB ago, g;
	VECTOR_OB OL, OL1, Lambda;
	INT k_first, k_len;
	BYTE *fname;

	km_dial_to_dv(&dv);
	lambda = dv.lambda;
	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	printf("KM-file: %s lambda = %ld\n", dv.KM_fname, lambda);
	fflush(stdout);
	
	km_read_ascii(dv.KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &I);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n); fflush(stdout);

	if (G_gen.s_obj_k() != VECTOR) {
		error("km_preproc_cb(): no generators");
		return;
		}
	if (RR.s_obj_k() != VECTOR) {
		error("km_preproc_cb(): no set representatives");
		return;
		}
	if (MM.s_obj_k() != VECTOR) {
		error("km_preproc_cb(): no KM-matrices (t,t+1)");
		return;
		}
	if (MM.s_li() < k) {
		error("km_preproc_cb(): too few KM-matrices (t,t+1)");
		return;
		}

	reduce_generators_labra(&G_gen, &go, FALSE /* f_verbose */, &labraG);
	printf("the group order is ");
	go.println(); fflush(stdout);
	labraG.calc_transversal_matrix(&TG, FALSE);
	deg = labraG.s_degree_i();
	if (deg != v) {
		error("km_preproc_cb() deg != v");
		return;
		}
		
	dc_calc_Ainf_t_via_MM(k, &Ainf, &MM, FALSE /* f_verbose */);
	
	nb_d = stab_go.s_li();
	K.m_il(nb_d);
	nb_d1 = 0;
	for (i = 0; i <= k; i++) {
		pR = (VECTOR_OP) RR.s_i(i);
		l = pR->s_li();
		for (j = 0; j < l; j++) {
			K.m_ii(nb_d1++, i);
			}
		}

	// printf("computing Aijk "); fflush(stdout);
	// calc_Aijk(&Ainf, &Aijk);
	// printf("finished !\n"); fflush(stdout);
		


	OL.m_il(n + 1);
	dc_find_first_len(&K, k, &k_first, &k_len);
	if (k_len != n) {
		error("km_preproc_cb(): k_len != n");
		return;
		}
	printf("computing vector of orbit length (we print the stabilizer orders):\n");
	for (j = 0; j < n; j++) {
		go.ganzdiv(stab_go.s_i(k_first + j), OL.s_i(j));
		stab_go.s_i(k_first + j)->print();
		// OL.s_i(j)->print();
		printf(" ");
		}
	printf("\n");
	calc_Lambda(v, t, k, lambda, &Lambda);
	printf("b = ");
	Lambda.s_i(0)->println();
	Lambda.s_i(0)->copy(OL.s_i(n));
	OL.gcd_all_elements_and_divide_out(&g, &OL1);
	printf("the gcd of all orbit length and b is ");
	g.println();
	for (j = 0; j <= n; j++) {
		OL1.s_i(j)->print();
		printf(" ");
		}
	printf("\n");
	KM.m_ilih(n, m + 1);
	for (j = 0; j < n; j++) {
		OL1.s_i(j)->copy(KM.s_ij(0, j));
		}
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			Mtk.s_ij(i, j)->copy(KM.s_ij(i + 1, j));
			}
		}
	
	
	printf("preprocessing with lambda = %ld\n", lambda); fflush(stdout);
	RHS.m_il(m + 1);
	for (i = 0; i < m; i++) {
		RHS.m_ii(i + 1, lambda);
		}
	OL1.s_i(n)->copy(RHS.s_i(0));
	
	km_preprocess(&KM, &RHS, lambda);

}
#endif

void show_KM_matrix(BYTE *KM_fname)
{
	char *fname;
	FILE *fp;
	char buf[BUFSIZE];
	XmTextPosition tp;

	fname = KM_fname;
	fp = fopen(fname, "r");
	if (fp == NIL) {
		printf("cannot open file %s!\n", fname);
		return;
		}
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NIL)
			break;
		if (strncmp(buf, "ASCII", 5) == 0) {
			while (TRUE) {
				if (fgets(buf, BUFSIZE, fp) == NIL)
					break;
				if (strncmp(buf, "ASCIIEND", 8) == 0) {
					break;
					}
				}
			}
		else {
			tp = XmTextGetLastPosition(Discreta_scrolled_text);
			XmTextInsert(Discreta_scrolled_text, tp, buf);
			}
		}
	fclose(fp);
}

void cb_do_LLL(INT f_silent, BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT lambda, INT f_iterate, INT nb_iterate)
{
	char *fname;
	FILE *fp;
	char s[1024];
	char buf[BUFSIZE];
	XmTextPosition tp;

	do_LLL(f_silent, KM_fname, c0_factor, beta, p, lambda, f_iterate, nb_iterate);
	fname = "lll.log.1";
	fp = fopen(fname, "r");
	if (fp == NIL) {
		printf("cannot open file %s !\n", fname);
		return;
		}
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NIL)
			break;
		tp = XmTextGetLastPosition(Discreta_scrolled_text);
		XmTextInsert(Discreta_scrolled_text, tp, buf);
		}
	fclose(fp);

}

void cb_do_iso_test(BYTE *KM_fname, INT lambda)
{
	BYTE base_name[1000];
	BYTE group_label[1000];
	BYTE cmd[1000];
	INT l;
	
	printf("isomorphism test for %s, lambda=%ld\n", KM_fname, lambda);
	
	strcpy(base_name, KM_fname);
	l = strlen(base_name);
	if (l > 4 && strcmp(base_name + l - 4, ".txt") == 0)
		base_name[l - 4] = 0;
	printf("base_name = %s\n", base_name);
	KM_fname_extract_group_label(KM_fname, group_label);
	printf("group_label = %s\n", group_label);
	
	sprintf(cmd, "t146.out gen_%s.txt %s_%ld.base_blocks", group_label, base_name, lambda);
	system(cmd);
	
	sprintf(cmd, "t149.out %s_%ld.geo", base_name, lambda);
	system(cmd);
	
	sprintf(cmd, "t153.out %s_%ld_c.geo", base_name, lambda);
	system(cmd);
	
	sprintf(cmd, "t147.out %s_%ld_c_merge.geo", base_name, lambda);
	system(cmd);
	
	sprintf(cmd, "latex %s_%ld_c_merge.tex", base_name, lambda);
	system(cmd);
	
	sprintf(cmd, "xdvi %s_%ld_c_merge.dvi", base_name, lambda);
	system(cmd);
	
	
}

void cb_do_iso_test_sylow(BYTE *KM_fname, INT p, INT lambda, INT fv)
{
	BYTE cmd[1000];
	FILE *fp;
	
	printf("isomorphism test by Sylow subgroup for %s, p=%ld, lambda=%ld\n", KM_fname, p, lambda);
	
	if (nb_primes(p) != 1) {
		printf("please enter a prime number for p!\n");
		return;
		}
	
	fp = fopen("iso.g", "w");
	
fprintf(fp, "#\n");
fprintf(fp, "# iso.g \n");
fprintf(fp, "#\n");
fprintf(fp, "# Evi Haberberger\n");
fprintf(fp, "# November 1999\n");
fprintf(fp, "#\n");
fprintf(fp, "# classifies the isomorphism types of designs with\n");
fprintf(fp, "# prescribed automorphism group and p-Sylowgroup for prescribed p\n");
fprintf(fp, "# main function is \"extract_iso\"\n");
fprintf(fp, "#\n");
fprintf(fp, "\n");
fprintf(fp, "Read(\"%s/lib/gap_programs/discreta.g\");\n", discreta_home);
fprintf(fp, "\n");
fprintf(fp, "Read(\"%s/lib/gap_programs/d1.g\");\n", discreta_home);
fprintf(fp, "Read(\"%s/lib/gap_programs/d2.g\");\n", discreta_home);
fprintf(fp, "Read(\"%s/lib/gap_programs/d3.g\");\n", discreta_home);
fprintf(fp, "Read(\"%s/lib/gap_programs/d4.g\");\n", discreta_home);
fprintf(fp, "Read(\"%s/lib/gap_programs/d5.g\");\n", discreta_home);
fprintf(fp, "Read(\"%s/lib/gap_programs/d6.g\");\n", discreta_home);
fprintf(fp, "\n");
fprintf(fp, "Read(\"%s/lib/gap_programs/d0.g\");\n", discreta_home);	
fprintf(fp, "\n");
fprintf(fp, "\n");
fprintf(fp, "km := \"%s\";\n", KM_fname);
fprintf(fp, "p := %ld;\n", p);
fprintf(fp, "lambda := %ld;\n", lambda);
fprintf(fp, "\n");
fprintf(fp, "extract_iso(km, p, lambda, %ld);\n", fv);
fprintf(fp, "\n");
fprintf(fp, "# the last number indicates, how extensive the report should be:\n");
fprintf(fp, "# 0 : no report\n");
fprintf(fp, "# 1 : short report\n");
fprintf(fp, "# 2 : long report\n");
fprintf(fp, "\n");
fprintf(fp, "\n");

	fclose(fp);
	
	
	sprintf(cmd, "gap iso.g");
	system(cmd);
	
}


