/* fg.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef SOLVABLE_TRUE

#include <DISCRETA/solvable.h>
#include <DISCRETA/sgl.h> /* for PrintOrbits() */
#ifdef SYM_DB_FG
#include <DISCRETA/fg1.h> // for FG_IO
#endif

#define MAX_NW 64

/*
 * ZE
 */

#if TEXDOCU
INT ZE_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "step"; break;
	case 1: s = "p"; break;
	case 2: s = "n0"; break;
	case 3: s = "n"; break;
	case 4: s = "P"; break;
	case 5: s = "A"; break;
	case 6: s = "Av"; break;
	default:
		return error("ZE::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}

#if TEXDOCU
INT ZE_OB::init(INT step, INT p, INT n, INT P)
#endif
{
	INT erg = OK;
	INT n0 = n / p;
	
	if (step <= 0)
		return error("ZE::init()|step <= 0");
	if (n0 * p != n)
		return error("ZE::init()|n0 * p != n");
	erg += m_il(7);
	c_obj_k(ZE);
	s_step()->m_i(step);
	s_p()->m_i(p);
	s_n0()->m_i(n0);
	s_n()->m_i(n);
	s_P()->m_i(P);
	erg += s_A()->m_il_n(step - 1);
	erg += s_Av()->m_il_n(step - 1);
	return erg;
}

#if TEXDOCU
INT ZE_OB::sprint(BYTE *s)
#endif
{
	BYTE str[512];
	
	sprintf(str, "p=%ld P=%ld A: ", s_p_i(), s_P_i());
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	s_A()->sprint(s);
	sprintf(str, " Av: ");
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	s_Av()->sprint(s);
	return OK;
}

#if TEXDOCU
INT ZE_OB::Print(FG_OP fg)
#endif
{
	BYTE str1[512], str2[512];
	
	str1[0] = 0; str2[0] = 0;
	fg->sprint_int_vec(s_A(), str1);
	fg->sprint_int_vec(s_Av(), str2);
	printf("ze %ld: p=%ld P=%ld n0 = %ld ", 
		s_step_i(), s_p_i(), s_P_i(), s_n0_i());
	printf("n = %ld A: (%s) Av: (%s)\n", 
		s_n_i(), str1, str2);
	return OK;
}

/*
 * FG
 */

#if TEXDOCU
INT FG_OB::field_name(INT i, INT j, BYTE *str)
#endif
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "version"; break;
	case 1: s = "n"; break;
	case 2: s = "m"; break;
	case 3: s = "o"; break;
	case 4: s = "GmZ_n"; break;
	case 5: s = "Gd_n"; break;
	case 6: s = "label"; break;
	case 7: s = "nb_ze"; break;
	case 8: s = "ze"; break;
	case 9: s = "theG"; break;
	case 10: s = "dim_n1"; break;
	case 11: s = "f_has_hash"; break;
	case 12: s = "hash"; break;
	case 13: s = "f_has_p0"; break;
	case 14: s = "p0"; break;
	case 15: s = "f_has_aut"; break;
	case 16: s = "ago"; break;
	case 17: s = "AutM"; break;
	case 18: s = "T"; break;
	case 19: s = "f_has_aut_classes"; break;
	case 20: s = "aut_classes"; break;
	case 21: s = "f_has_sgl"; break;
	case 22: s = "sgl"; break;
	case 23: s = "f_has_sylow_type"; break;
	case 24: s = "sylow_type"; break;
	case 25: s = "C"; break;
	default:
		return error("FG::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}
#if 0
/* version 1.00: */
INT FG_OB::field_name(INT i, INT j, BYTE *str)
{
	BYTE *s;
	
	switch (i) {
	case 0: s = "version"; break; /* |-> 0 */
	case 1: s = "nb_ze"; break; /* |-> 7 */
	case 2: s = "n"; break; /* |-> 1 */
	case 3: s = "m"; break; /* |-> 2 */
	case 4: s = "o"; break; /* |-> 3 */
	case 5: s = "GmZ_n"; break; /* |-> 4 */
	case 6: s = "Gd_n"; break; /* |-> 5 */
	case 7: s = "f_has_hash"; break; /* |-> 11 */
	case 8: s = "hash"; break; /* |-> 12 */
	case 9: s = "f_has_p0"; break; /* |-> 13 */
	case 10: s = "p0"; break; /* |-> 14 */
	case 11: s = "dim_n1"; break; /* |-> 10 */
	case 12: s = "ze"; break; /* |-> 8 */
	case 13: s = "theG"; break; /* |-> 9 */
	case 14: s = "ago"; break; /* |-> 16 */
	case 15: s = "AutM"; break; /* |-> 17 */
	case 16: s = "T"; break; /* |-> 18 */
	case 17: s = "label"; break; /* |-> 6 */
	case 18: s = "C"; break; /* |-> 25 */
	case 19: s = "sylow_type"; break; /* |-> 24 */
	default:
		return error("FG::field_name()|i out of range");
	}
	strcpy(str, s);
	return OK;
}
		INT embedding[] = { 0, 7, 1, 2, 3, 4, 5, 
			11, 12, 13, 14, 10, 8, 9, 
			16, 17, 18, 6, 25, 24  };
#endif

#if TEXDOCU
INT FG_OB::export_GAP(FILE *fp)
#endif
{
	BYTE s1[1024];
	INT i, j, gj, a, b, c, l, p, P;
	INT f_first = TRUE;
	ZE_OP ze;

	if (s_theG()->s_obj_k() == EMPTY)
		recalc_table(FALSE);
	l = s_nb_ze_i();
	for (i = 0; i < l; i++) {
		fprintf(fp, "g%ld := AbstractGenerator( \"g%ld\" );\n", i + 1, i + 1);
		}
	fprintf(fp, "GROUP := rec( generators := [");
	fflush(fp);
	for (i = 0; i < l; i++) {
		fprintf(fp, "g%ld", i + 1);
		if (i < l - 1)
			fprintf(fp, ", ");
		}
	fprintf(fp, "],\nrelators := [\n");
	fflush(fp);
	for (i = 0; i < l; i++) {
		if (!f_first)
			fprintf(fp, ", \n");
		ze = s_ze_i(i);
		p = ze->s_p_i();
		P = ze->s_P_i();
		s1[0] = 0;
		sprint_GAP_int(P, s1);
		fprintf(fp, "g%ld^%ld / ( %s )", l - i, p, s1);
		fflush(fp);
		f_first = FALSE;
		}
	for (i = 0; i < l; i++) {
		ze = s_ze_i(i);
		for (j = 0; j < i; j++) {
			if (!f_first)
				fprintf(fp, ", \n");
			gj = g_i(j);
			a = gt_inv(gj);
			b = ze->s_A_ii(j); // b := g_i^{-1} * g_j * g_i
			c = gt_mult(a, b); // c := g_j^{-1} * b = [g_j, g_i]
			s1[0] = 0;
			sprint_GAP_int(c, s1);
			fprintf(fp, "Comm( g%ld, g%ld ) / ( %s )", l - j, l - i, s1);
			fflush(fp);
			f_first = FALSE;
			}
		}
	fprintf(fp, "\n] );\n");
	fprintf(fp, "Add( Glist, AgGroupFpGroup( GROUP ) );\n\n");
	fflush(fp);
	return OK;
}

#if TEXDOCU
INT FG_OB::save_ascii2(FILE *fp, INT f_long)
#endif
{
	ZE_OP ze;
	INT i, j, k, l;
	BYTE *head;
	
	if (f_long)
		head = "GROUPSTART";
	else
		head = "GROUP";
	
	fprintf(fp, "%s %ld %ld\n", head, s_n_i(), s_m_i());
	l = s_nb_ze_i();
	for (i = 0; i < l; i++) {
		ze = s_ze_i(i);
		fprintf(fp, "%ld %ld", ze->s_p_i(), ze->s_P_i());
		for (j = 0; j < i; j++) {
			k = ze->s_A_ii(j);
			fprintf(fp, " %ld", k);
			}
		fprintf(fp, "\n");
		}
	
	if (!f_long)
		return OK;
	
	if (s_f_has_hash_i()) {
		fprintf(fp, "HASH\n");
		s_hash()->save_ascii(fp);
		}
	
	if (s_f_has_p0_i()) {
		fprintf(fp, "P0\n");
		s_p0()->save_ascii(fp);
		}
	
	if (s_f_has_aut_i()) {
		fprintf(fp, "AUT\n");
		s_ago()->save_ascii(fp);
		s_AutM()->save_ascii(fp);
		s_T()->save_ascii(fp);
		}
	
	if (s_f_has_aut_classes_i()) {
		fprintf(fp, "AUTCLASSES\n");
		s_aut_classes()->save_ascii(fp);
		}
	
	if (s_f_has_sgl_i()) {
		fprintf(fp, "SGL\n");
		s_sgl()->save_ascii(fp);
		}
	
	if (s_f_has_sylow_type_i()) {
		fprintf(fp, "SYLOW_TYPE\n");
		s_sylow_type()->save_ascii(fp);
		}
	
	fprintf(fp, "GROUPEND\n");
	return OK;
}

#if TEXDOCU
INT FG_OB::save_ascii(FILE *fp)
#endif
{
	return save_ascii2(fp, TRUE);
}

#if TEXDOCU
INT FG_OB::save_ascii_short(FILE *fp)
#endif
{
	return save_ascii2(fp, FALSE);
}

#define BUFSIZE 1024

#if TEXDOCU
INT FG_OB::read_ascii(FILE *fp)
#endif
{
	INT primes[MAX_NW];
	INT P[MAX_NW];
	INT A[MAX_NW][MAX_NW];
	BYTE buf[BUFSIZE];
	BYTE str[256];
	BYTE *p;
	INT n, m, l, i, j, k, n0, f_with_info;
	ZE_OP ze, ze1;
	PERMUTATION_OB aut, autv;
	
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			return FALSE;
			}
		p = buf;
		if (!s_scan_token(&p, str)) {
			printf("FG_OB::read_ascii() can't read label 'GROUP' !\n");
			return FALSE;
			}
		if (strcmp(str, "GROUP") == 0) {
			f_with_info = FALSE;
			break;
			}
		if (strcmp(str, "GROUPSTART") == 0) {
			f_with_info = TRUE;
			break;
			}
		}
	if (!s_scan_int(&p, &n)) {
		printf("FG_OB::read_ascii() can't read n !\n");
		return FALSE;
		}
	if (!s_scan_int(&p, &m)) {
		printf("FG_OB::read_ascii() can't read m !\n");
		return FALSE;
		}
	l = nb_primes(n);
	printf("reading group %ld#%ld l = %ld\n", n, m, l);
	fflush(stdout);
	for (i = 0; i < l; i++) {
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			printf("FG_OB::read_ascii() GROUP not complete !\n");
			return FALSE;
			}
		p = buf;
		if (!s_scan_int(&p, primes + i)) {
			printf("FG_OB::read_ascii() can't read p !\n");
			return FALSE;
			}
		if (!s_scan_int(&p, P + i)) {
			printf("FG_OB::read_ascii() can't read P !\n");
			return FALSE;
			}
		for (j = 0; j < i; j++) {
			if (!s_scan_int(&p, &k)) {
				printf("FG_OB::read_ascii() can't read A[i][j] !\n");
				return FALSE;
				}
			A[i][j] = k;
			}
		}
#if 0
	printf("primes:\n");
	for (i = 0; i < l; i++) {
		printf("%ld ", primes[i]);
		}
	printf("\n");
	printf("P:\n");
	for (i = 0; i < l; i++) {
		printf("%ld ", P[i]);
		}
	printf("\n");
	printf("A:\n");
	for (i = 0; i < l; i++) {
		for (j = 0; j < i; j++) {
			printf("%ld ", A[i][j]);
			}
		printf("\n");
		}
	printf("\n");
#endif
	printf("calling init()\n");
	fflush(stdout);
	init(l, primes);
	if (s_n_i() != n) {
		error("FG_OB::read_ascii() s_n_i() != n");
		return FALSE;
		}
	s_m()->m_i(m);
	sprintf(str, "%ld#%ld", n, m);
	s_label()->init(str);
	printf("read group %s\n", str);
	fflush(stdout);
	for (i = 0; i < l; i++) {
		ze = s_ze_i(i);
		ze->s_P()->m_i(P[i]);
		for (j = 0; j < i; j++) {
			ze->s_A_i(j)->m_i(A[i][j]);
			}
		printf("calling calc_aut_i() i = %ld\n", i);
		calc_aut_i(i, &aut, FALSE /* f_use_table */);
		aut.invers(&autv);
		for (j = 0; j < i; j++) {
			ze1 = s_ze_i(j);
			n0 = ze1->s_n0_i();
			k = autv.s_ii(n0) - 1;
			ze->s_Av_i(j)->m_i(k);
			}
		}
	/* recalc_table(FALSE); */
	/* s_theG()->freeself(); */

	if (!f_with_info) {
		printf("read group %s\n", str);
		fflush(stdout);
		// Print();
		// fflush(stdout);
		return TRUE;
		}
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			return FALSE;
			}
		p = buf;
		if (!s_scan_token(&p, str)) {
			printf("FG_OB::read_ascii() error reading footer !\n");
			return FALSE;
			}
		if (strcmp(str, "GROUPEND") == 0)
			return TRUE;
		
		if (strcmp(str, "HASH") == 0) {
			printf("FG_OB::read_ascii() reading HASH !\n");
			fflush(stdout);
			s_hash()->load_ascii(fp);
			s_f_has_hash()->m_i(TRUE);
			printf("FG_OB::read_ascii() finished reading HASH !\n");
			fflush(stdout);
			continue;
			}

		if (strcmp(str, "P0") == 0) {
			printf("FG_OB::read_ascii() reading P0 !\n");
			fflush(stdout);
			s_p0()->load_ascii(fp);
			s_f_has_p0()->m_i(TRUE);
			printf("FG_OB::read_ascii() finished reading P0 !\n");
			fflush(stdout);
			continue;
			}

		if (strcmp(str, "AUT") == 0) {
			printf("FG_OB::read_ascii() reading AUT !\n");
			fflush(stdout);
			s_ago()->load_ascii(fp);
			s_AutM()->load_ascii(fp);
			s_T()->load_ascii(fp);
			s_f_has_aut()->m_i(TRUE);
			printf("FG_OB::read_ascii() finished reading AUT !\n");
			fflush(stdout);
			continue;
			}

		if (strcmp(str, "AUTCLASSES") == 0) {
			printf("FG_OB::read_ascii() reading AUTCLASSES !\n");
			fflush(stdout);
			s_aut_classes()->load_ascii(fp);
			s_f_has_aut_classes()->m_i(TRUE);
			printf("FG_OB::read_ascii() finished reading AUTCLASSES !\n");
			fflush(stdout);
			continue;
			}
		
		if (strcmp(str, "SGL") == 0) {
			printf("FG_OB::read_ascii() reading SGL !\n");
			fflush(stdout);
			s_sgl()->load_ascii(fp);
			s_f_has_sgl()->m_i(TRUE);
			printf("FG_OB::read_ascii() finished reading SGL !\n");
			fflush(stdout);
			continue;
			}

		if (strcmp(str, "SYLOW_TYPE") == 0) {
			printf("FG_OB::read_ascii() reading SYLOW_TYPE !\n");
			fflush(stdout);
			s_sylow_type()->load_ascii(fp);
			s_f_has_sylow_type()->m_i(TRUE);
			printf("FG_OB::read_ascii() finished reading SYLOW_TYPE !\n");
			fflush(stdout);
			continue;
			}

		}

	return TRUE;
}

#if TEXDOCU
INT FG_OB::read_ascii_tg(FILE *fp, INT n)
#endif
{
	INT primes[MAX_NW];
	INT P[MAX_NW];
	INT A[MAX_NW][MAX_NW];
	INT A_[10000];
	BYTE buf[BUFSIZE];
	BYTE str[256];
	BYTE *p;
	INT m, l, i, j, jj, k, n0;
	ZE_OP ze, ze1;
	PERMUTATION_OB aut, autv;
	
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			return FALSE;
			}
		p = buf;
		if (!s_scan_token(&p, str)) {
			printf("FG_OB::read_ascii_tg() can't read label 'GROUP' !\n");
			return FALSE;
			}
		if (strcmp(str, "GROUP") == 0)
			break;
		}
	m = -1;
	l = nb_primes(n);
	printf("reading group %ld#%ld l = %ld\n", n, m, l);
	fflush(stdout);
	for (i = 0; i < l; i++) {
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			printf("FG_OB::read_ascii_tg() can't read p= z= line!\n");
			return FALSE;
			}
		p = buf;
		sscanf(p, "p=%ld z=%ld", primes + i, P + i);
		P[i]--;
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			printf("FG_OB::read_ascii_tg() can't read permutation !\n");
			return FALSE;
			}
		n0 = 1;
		for (j = 0; j < i; j++) {
			n0 *= primes[j];
			}
		
		p = buf;
		for (j = 0; j < n0; j++) {
			if (!s_scan_int(&p, &k)) {
				printf("FG_OB::read_ascii_tg() can't read A_[j] !\n");
				return FALSE;
				}
			A_[j] = k - 1;
			}

		for (j = 0; j < i; j++) {
			n0 = 1;
			for (jj = 0; jj < j; jj++) {
				n0 *= primes[jj];
				}
			A[i][j] = A_[n0];
			}
		}
	if (fgets(buf, BUFSIZE, fp) == NULL) {
		printf("FG_OB::read_ascii_tg() can't read canonical form !\n");
		return FALSE;
		}
	printf("calling init()\n");
	fflush(stdout);
	init(l, primes);
	if (s_n_i() != n) {
		error("FG_OB::read_ascii_tg() s_n_i() != n");
		return FALSE;
		}
	sprintf(str, "%ld#%ld", n, m);
	s_label()->init(str);
	printf("read group %s\n", str);
	fflush(stdout);
	for (i = 0; i < l; i++) {
		ze = s_ze_i(i);
		ze->s_P()->m_i(P[i]);
		for (j = 0; j < i; j++) {
			ze->s_A_i(j)->m_i(A[i][j]);
			}
		printf("calling calc_aut_i() i = %ld\n", i);
		calc_aut_i(i, &aut, FALSE /* f_use_table */);
		aut.invers(&autv);
		for (j = 0; j < i; j++) {
			ze1 = s_ze_i(j);
			n0 = ze1->s_n0_i();
			k = autv.s_ii(n0) - 1;
			ze->s_Av_i(j)->m_i(k);
			}
		}
	/* recalc_table(FALSE); */
	s_theG()->freeself();
	return TRUE;
}

#if TEXDOCU
INT FG_OB::init(INT nb_ze, INT *primes)
#endif
{
	INT erg = OK;
	INT n, i;
	BYTE str[256];
	CONTI_OP C;
	
	erg += m_il(26);
	c_obj_k(FG);
	s_version()->m_i(110); /* version 1.10 */
	s_nb_ze()->m_i(nb_ze);
	s_n()->m_i(0);
	s_m()->m_i(0);
	s_o()->m_i(0);
	s_GmZ_n()->m_i(0);
	s_Gd_n()->m_i(0);
	s_f_has_hash()->m_i(FALSE);
	s_f_has_p0()->m_i(FALSE);
	s_f_has_aut()->m_i(FALSE);
	s_f_has_aut_classes()->m_i(FALSE);
	s_f_has_sgl()->m_i(FALSE);
	s_f_has_sylow_type()->m_i(FALSE);
	s_dim_n1()->m_i(0);
	erg += s_ze()->m_il(nb_ze);
	n = 1;
	str[0] = 0;
	for (i = 0; i < nb_ze; i++) {
		n *= primes[i];
		erg += s_ze_i(i)->init(i + 1 /* step */, 
			primes[i] /* p */, 
			n /* n */, 0 /* P */);
		sprintf(str + strlen(str), "%ld", primes[i]);
		if (i < nb_ze - 1)
			sprintf(str + strlen(str), ".");
		}
	erg += s_C()->init();
	C = s_C();
	erg += C->add_label(str);
	sprintf(str, "%ld", n);
	erg += s_label()->init(str);
	s_n()->m_i(n);
	((INTEGER_OP) s_ago())->m_i(0);
	/* erg += nw_enum(); */
	erg += s_sylow_type()->m_il(0);
	return erg;
}

#if TEXDOCU
INT FG_OB::update()
#endif
{
	if (s_li() == 15) {
		FG_OB G;
		INT embedding[] = { 1, 2, 3, 4, 5, 6, 
			11, 12, 13, 14, 15, 16, 17, 18, 19 };
		INT i, l;

		printf("FG::update() updating from old version (with 15 entries) to 1.00\n");
		fflush(stdout);
		G.m_il(20);
		c_obj_k(FG);
		l = sizeof(embedding) / sizeof(INT);
		for (i = 0; i < l; i++) {
			s_i(i)->swap(G.s_i(embedding[i]));
			}
		G.s_version()->m_i(100);
		G.s_f_has_hash()->m_i(FALSE);
		G.s_f_has_p0()->m_i(FALSE);
		G.swap(this);
		printf("FG::update() updated to version 1.00\n");
		fflush(stdout);
		}
	if (s_version_i() == 100) {
		FG_OB G;
		INT embedding[] = { 0, 7, 1, 2, 3, 4, 5, 
			11, 12, 13, 14, 10, 8, 9, 
			16, 17, 18, 6, 25, 24  };
		INT i, l;

		printf("FG::update() updating from 1.00 version to 1.10\n");
		fflush(stdout);
		G.m_il(26);
		c_obj_k(FG);
		l = sizeof(embedding) / sizeof(INT);
		for (i = 0; i < l; i++) {
			s_i(i)->swap(G.s_i(embedding[i]));
			}
		if (G.s_AutM()->s_obj_k() != EMPTY)
			s_f_has_aut()->m_i(TRUE);
		else
			s_f_has_aut()->m_i(FALSE);
		s_f_has_aut_classes()->m_i(FALSE);
		s_f_has_sgl()->m_i(FALSE);
		if (G.s_sylow_type()->s_li() != 0)
			s_f_has_sylow_type()->m_i(TRUE);
		else
			s_f_has_sylow_type()->m_i(FALSE);
		G.s_version()->m_i(110);
		G.swap(this);
		printf("FG::update() updated to version 1.10\n");
		fflush(stdout);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::init_extension(FG_OP G0, INT p, INT P, PERMUTATION_OP aut)
#endif
{
	PERMUTATION_OB aut_inv;
	INT primes[MAX_NW];
	ZE_OB ze;
	ZE_OP p_ze;
	INT n0, n1, i, g_i;
	INT g_i_bild, g_i_bild_inv, old_nb_ze;
	
	aut->invers(&aut_inv);
	n0 = G0->s_n_i();
	n1 = n0 * p;
	ze.init(G0->s_nb_ze_i() + 1 /* step */, p, n1 /* n */, P);
	for (i = 0; i < G0->s_nb_ze_i(); i++) {
		p_ze = G0->s_ze_i(i);
		g_i = p_ze->s_n0_i();
		g_i_bild = aut->s_ii(g_i) - 1;
		g_i_bild_inv = aut_inv.s_ii(g_i) - 1;
		ze.s_A()->m_ii(i, g_i_bild);
		ze.s_Av()->m_ii(i, g_i_bild_inv);
		}

	old_nb_ze = G0->s_nb_ze_i();
	for (i = 0; i < old_nb_ze; i++) {
		p_ze = G0->s_ze_i(i);
		primes[i] = p_ze->s_p_i();
		}
	primes[old_nb_ze] = p;
	init(old_nb_ze + 1, primes);
	for (i = 1; i < old_nb_ze; i++)
		G0->s_ze_i(i)->copy(s_ze_i(i));

	ze.swap(s_ze_i(s_ze()->s_li() - 1));
	return OK;
}

#if TEXDOCU
INT FG_OB::init_Zp(INT p)
#endif
{
	init(1, &p);
	theG(FALSE /* f_verbose */);
	Aut(FALSE /* f_verbose */, FALSE /* f_very_verbose */);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint(BYTE *s)
#endif
{
	BYTE str[10000];
	
#if FALSE
	CONTI_OP C;

	C = s_C();
	if (C->s_nb_L_i()) {
		s1 = C->s_L_i_str(0);
		}
	else
		s1 = "";
#endif
	sprintf(str, "FG: %s ", s_label_s());
	sprint1(str);
	sprint2(str);
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	sprint3(str);
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint1(BYTE *s)
#endif
{
	BYTE str[512];
	BYTE str1[512];
	
	str1[0] = 0;
	s_ago()->sprint(str1);
	sprintf(str, "ago = %s o: %ld GmZ_n = %ld Gd_n = %ld ", 
		str1, s_o_i(), s_GmZ_n_i(), s_Gd_n_i());
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint2(BYTE *s)
#endif
{
	BYTE str[512];
	INT i, len;
	
	s_sylow_type()->sprint(s);
	strcat(s, " ");
	len = s_nb_ze_i();
	for (i = 0; i < len; i++) {
		str[0] = 0;
		sprint_ze(i, str);
		if (strlen(s) + strlen(str) < 200)
			strcat(s, str);
		if (i < len - 1)
			strcat(s, " ");
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint3(BYTE *s)
#endif
{
	// BYTE str[1024];
	
	sprintf(s, 
		"has_hash = %ld has_p0 = %ld "
		"has_aut = %ld has_aut_classes = %ld"
		"has_sgl = %ld has_sylow_type = %ld", 
		s_f_has_hash_i(), s_f_has_p0_i(), 
		s_f_has_aut_i(), s_f_has_aut_classes_i(), 
		s_f_has_sgl_i(), s_f_has_sylow_type_i());
	return OK;
}

#if TEXDOCU
INT FG_OB::fprint_gen(FILE *fp)
#endif
{
	INT j, l;
	BYTE str[1024];

	str[0] = 0;
	sprint1(str);
	fprintf(fp, "GROUP %s %s\n", s_label_s(), str);
	l = s_nb_ze_i();
	for (j = 0; j < l; j++) {
		str[0] = 0;
		sprint_ze(j, str);
		fprintf(fp, "%s\n", str);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::fprint_pres(FILE *fp)
#endif
{
	INT j, l;
	BYTE str[1024];

	l = s_nb_ze_i();
	for (j = 0; j < l; j++) {
		str[0] = 0;
		sprint_ze(j, str);
		fprintf(fp, "%s\n", str);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::Print_gen(void)
#endif
{
	return fprint_gen(stdout);
}

#if TEXDOCU
INT FG_OB::Print(void)
#endif
{
	BYTE str[512], str1[256], str3[1024];
	INT i, len;
	
	str1[0] = 0;
	s_ago()->sprint(str1);
	printf("FG: %s ago = %s\n", s_label_s(), str1);
	fflush(stdout);
	len = s_nb_ze_i();
	for (i = 0; i < len; i++) {
		str[0] = 0;
		sprint_ze(i, str);
		printf("%s\n", str);
		}
	printf("m: ");
	fflush(stdout);
	s_m()->print();
	fflush(stdout);
	printf(" n: ");
	fflush(stdout);
	s_n()->print();
	fflush(stdout);
	printf(" o: ");
	fflush(stdout);
	s_o()->print();
	fflush(stdout);
	printf(" GmZ_n: ");
	fflush(stdout);
	s_GmZ_n()->print();
	fflush(stdout);
	printf(" Gd_n: ");
	fflush(stdout);
	s_Gd_n()->print ();
	fflush(stdout);
	printf(" ago: ");
	fflush(stdout);
	s_ago()->println();
	fflush(stdout);

	str3[0] = 0;
	sprint3(str3);
	printf("%s.\n", str3);
	fflush(stdout);
#if 0
	f_has_hash = s_f_has_hash_i();
	f_has_p0 = s_f_has_p0_i();
	printf("f_has_hash = %ld f_has_p0 = %ld\n", 
		f_has_hash, f_has_p0);
	fflush(stdout);
	if (f_has_hash)
		s_hash()->println();
	fflush(stdout);
	if (f_has_p0)
		s_p0()->println();
	fflush(stdout);
#endif
	if (s_f_has_sylow_type_i()) {
		printf("sylow type: ");
		fflush(stdout);
		s_sylow_type()->println();
		}
	if (s_f_has_aut_i()) {
		printf("AutM: \n");
		fflush(stdout);
		/* s_AutM()->Print();
		printf("\n"); */
		print_AutM_T();
		fflush(stdout);
		}
	if (s_f_has_sgl_i()) {
		printf("SGL: \n");
		fflush(stdout);
		{
			SGL_OP L;
			
			L = (SGL_OP) s_sgl();
			L->PrintOrbits(FALSE, DO_TYPE_FG, this);
		}
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::print_AutM_T()
#endif
{
	MATRIX_OP M;
	INT i, j, k, h, l;
	PERMUTATION_OP p;
	VECTOR_OP Ti;
	BYTE str[2048];

	M = s_AutM();
	if (M->s_obj_k() != MATRIX)
		return OK;
	h = M->s_hi();
	l = M->s_li();
	k = 0;
	for (i = 0; i < h; i++) {
		for (j = 0; j < l; j++) {
			p = (PERMUTATION_OP) M->s_ij(i, j);
			if (p->s_obj_k() != PERMUTATION) 
				printf(" * ");
			else {
				printf("%2ld ", k);
				k++;
				}
			}
		printf("\n");
		}
	k = 0;
	for (i = 0; i < h; i++) {
		printf("transversal:\n");
		Ti = s_T_i(i);
		Ti->println();
		for (j = 0; j < l; j++) {
			p = (PERMUTATION_OP) M->s_ij(i, j);
			if (p->s_obj_k() == PERMUTATION)  {
				printf("%ld (in col %ld): ", k, j);
				/* p->println(); */
				str[0] = 0;
				sprint_aut_on_base(p, str);
				printf("%s\n", str);
				k++;
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_aut_on_base(PERMUTATION_OP p, BYTE *str)
#endif
{
	INT i, gi, j, l;
	BYTE s1[256], s2[256];
		
	l = s_nb_ze_i();
	for (i = 0; i < l; i++) {
		gi = g_i(i);
		s1[0] = 0;
		s2[0] = 0;
		sprint_int(gi, s1);
		j = p->s_ii(gi) - 1;
		sprint_int(j, s2);
		sprintf(str + strlen(str), "%s -> %s", s1, s2);
		if (i < l - 1)
			sprintf(str + strlen(str), ", ");
			
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_ze(INT i, BYTE *s)
#endif
{
	BYTE str[512];
	ZE_OP ze;
	
	str[0] = 0;
	ze = s_ze_i(i);
	sprint_int(ze->s_n0_i(), str);
	sprintf(str + strlen(str), "^%ld=", ze->s_p_i());
	sprint_int(ze->s_P_i(), str);
	if (i != 0) {
		sprintf(str + strlen(str), " (");
		if (strlen(s) + strlen(str) < 200)
			strcat(s, str);
		str[0] = 0;
		sprint_int_vec(ze->s_A(), str);
		strcat(str, ")");
		}
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_tex_ze(INT i, BYTE *s)
#endif
{
	BYTE str[512];
	ZE_OP ze;
	
	str[0] = 0;
	ze = s_ze_i(i);
	sprint_tex_int(ze->s_n0_i(), str);
	sprintf(str + strlen(str), "^{%ld}=", ze->s_p_i());
	sprint_tex_int(ze->s_P_i(), str);
	if (i != 0) {
		sprintf(str + strlen(str), " (");
		if (strlen(s) + strlen(str) < 200)
			strcat(s, str);
		str[0] = 0;
		sprint_tex_int_vec(ze->s_A(), str);
		strcat(str, ")");
		}
	if (strlen(s) + strlen(str) < 200)
		strcat(s, str);
	return OK;
}

#if TEXDOCU
INT FG_OB::fprint_int(FILE *fp, INT i)
#endif
{
	BYTE str[1024];
	
	str[0] = 0;
	sprint_int(i, str);
	fprintf(fp, "%s", str);
	return OK;
}

#if TEXDOCU
INT FG_OB::print_int(INT i)
#endif
{
	BYTE str[1024];
	
	str[0] = 0;
	sprint_int(i, str);
	printf("%s", str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_int(INT i, BYTE *str)
#endif
{
	SHORT nw[MAX_NW];
	
	int2nw(i, nw);
	NW_2_str(nw, s_nb_ze_i() - 1, str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_tex_int(INT i, BYTE *str)
#endif
{
	SHORT nw[MAX_NW];
	
	int2nw(i, nw);
	NW_2_str_tex(nw, s_nb_ze_i() - 1, str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_GAP_int(INT i, BYTE *str)
#endif
{
	SHORT nw[MAX_NW];
	
	int2nw(i, nw);
	NW_2_str_GAP(nw, s_nb_ze_i() - 1, s_nb_ze_i(), str);
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_int_vec(VECTOR_OP V, BYTE *str)
#endif
{
	INT i, l, ii;
	
	l = V->s_li();
	for (i = 0; i < l; i++) {
		ii = V->s_ii(i);
		sprint_int(ii, str);
		if (i < l - 1)
			strcat(str, " ");
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_tex_int_vec(VECTOR_OP V, BYTE *str)
#endif
{
	INT i, l, ii;
	
	l = V->s_li();
	for (i = 0; i < l; i++) {
		ii = V->s_ii(i);
		sprint_tex_int(ii, str);
		if (i < l - 1)
			strcat(str, ", \\, ");
		}
	return OK;
}

#if TEXDOCU
INT do_fg_test(void)
#endif
{
	FG_OB Q8;
	INT primes[] = { 2, 2, 2 };

	Q8.init(3, primes);
	Q8.do_Q8();
	Q8.sylow_type(TRUE /* f_verbose */);
	return OK;
}

#if TEXDOCU
INT FG_OB::do_Q8()
#endif
{
	INT erg = OK;
	// ZE_OP ze0 = s_ze_i(0);
	ZE_OP ze1 = s_ze_i(1);
	ZE_OP ze2 = s_ze_i(2);
	
	ze1->s_P()->m_i(1);
	ze1->s_A_i(0)->m_i(1);
	ze1->s_Av_i(0)->m_i(1);
	
	ze2->s_P()->m_i(1);
	ze2->s_A_i(0)->m_i(1);
	ze2->s_Av_i(0)->m_i(1);
	ze2->s_A_i(1)->m_i(3);
	ze2->s_Av_i(1)->m_i(3);
	
	printf("vor theG()\n");
	fflush(stdout);
	erg += theG(TRUE /* f_verbose */);
	printf("nach theG()\n");
	fflush(stdout);
	
	erg += Aut(TRUE /* f_verbose */, TRUE /* f_very_verbose */);
	
	return erg;
}

#if TEXDOCU
INT FG_OB::do_Q8_2()
#endif
{
	INT erg = OK;
	// ZE_OP ze0 = s_ze_i(0);
	ZE_OP ze1 = s_ze_i(1);
	ZE_OP ze2 = s_ze_i(2);
	ZE_OP ze3 = s_ze_i(3);
	
	ze1->s_P()->m_i(1);
	ze1->s_A_i(0)->m_i(1);
	ze1->s_Av_i(0)->m_i(1);
	
	ze2->s_P()->m_i(1);
	ze2->s_A_i(0)->m_i(1);
	ze2->s_Av_i(0)->m_i(1);
	ze2->s_A_i(1)->m_i(3);
	ze2->s_Av_i(1)->m_i(3);
	
	ze3->s_P()->m_i(0);
	ze3->s_A_i(0)->m_i(1);
	ze3->s_Av_i(0)->m_i(1);
	ze3->s_A_i(1)->m_i(2);
	ze3->s_Av_i(1)->m_i(2);
	ze3->s_A_i(2)->m_i(4);
	ze3->s_Av_i(2)->m_i(4);
	
	erg += theG(TRUE /* f_verbose */);
	
	erg += Aut(TRUE /* f_verbose */, TRUE /* f_very_verbose */);
	
	return erg;
}

#if TEXDOCU
INT FG_OB::calc_ago_longint(LONGINT_OP ago)
#endif
{
	LONGINT_OB b, c;
	INT i, j, l;

	ago->m_i(1);
	if (!s_f_has_aut_i()) {
		printf("FG::calc_ago_longint(): f_has_aut is FALSE\n");
		return OK;
		}
	if (s_T()->s_obj_k() != VECTOR) {
		printf("FG::calc_ago_longint(): T (and aut ?) not calculated\n");
		return OK;
		}
	l = s_T()->s_li();
	for (i = 0; i < l; i++) {
		j = s_T_i(i)->s_li();
		b.m_i(j);
		ago->mult_longint(&b, &c);
		c.swap(ago);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::sprint_ago(BYTE *str)
#endif
{
	LONGINT_OB ago;
	BYTE s[256], *s1;
	INT i, j, l;

	calc_ago_longint(&ago);
	s[0] = 0;
	ago.sprint(s);
	/* LONGINT->sprint() eventually begins with spaces */
	s1 = s;
	while (*s1 == ' ')
		s1++;
	sprintf(str + strlen(str), "%s =", s1);
	l = s_T()->s_li();
	for (i = 0; i < l; i++) {
		j = s_T_i(i)->s_li();
		sprintf(str + strlen(str), " %ld", j);
		}
	return OK;
}

#if TEXDOCU
INT FG_OB::is_nilpotent()
#endif
{
	VECTOR_OP V;
	INT i, len;

	if (!s_f_has_sylow_type_i()) {
		printf("FG::is_nilpotent(): no sylow type\n");
		fflush(stdout);
		return -1;
		}
	V = s_sylow_type();
	len = V->s_li();
	for (i = 0; i < len; i++) {
		if (V->s_ii(i) != 1)
			return FALSE;
		}
	return TRUE;
}

#if TEXDOCU
INT FG_OB::nb_gen()
#endif
{
	return s_nb_ze_i();
}

#if TEXDOCU
INT FG_OB::g_i(INT i)
#endif
{
	ZE_OP ze;
	
	ze = s_ze_i(i);
	return ze->s_n0_i();
}

#if TEXDOCU
INT FG_OB::solvable_group_from_catalogue(INT n, INT m)
#else
//PRE
// reads the m-th group of order n from the solvable group library 
// $1 \le m \le$ number of solvable groups of order n.
// in DISCRETA\_HOME/lib/solvable\_groups/ASCII
// m == -1 means: read all groups of order n.
///PRE
#endif
{
	BYTE str[1024];
	FILE *fp;
	INT mm;
	
	if (discreta_home == NIL) {
		printf("FG::solvable_group_from_catalogue() discreta_home is not set !\n");
		fflush(stdout);
		return error("");
		}
	sprintf(str, "%s/lib/solvable_groups/ASCII/g_%ld_fg.txt", discreta_home, n);
	if (file_size(str) <= 0) {
		printf("FG::solvable_group_from_catalogue() file %s does not exist !\n", str);
		fflush(stdout);
		return error("");
		}
	fp = fopen(str, "r");
	mm = 0;
	while (read_ascii(fp)) {
		mm++;
		if (m < 0 || (mm == m)) {
			break;
			}
		}
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT FG_OB::solvable_group_from_file(BYTE *fname, INT m)
#else
//PRE
// reads the m-th group of the given file.
///PRE
#endif
{
	FILE *fp;
	INT mm;
	
	if (file_size(fname) <= 0) {
		printf("FG::solvable_group_from_file() file %s does not exist !\n", fname);
		fflush(stdout);
		return error("");
		}
	fp = fopen(fname, "r");
	mm = 0;
	while (read_ascii(fp)) {
		mm++;
		if (m < 0 || (mm == m)) {
			break;
			}
		}
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT FG_OB::solvable_group_of_order_n_first(INT n, FG_IO **fgio)
#endif
{
	BYTE fname[1024];
	BYTE str1[1024];
	FG_IO *io;

	if (discreta_home == NIL) {
		printf("FG_OB::solvable_group_of_order_n_first() discreta_home is not set !\n");
		fflush(stdout);
		return FALSE;
		}
	sprintf(fname, "%s/lib/solvable_groups/ASCII/g_%ld_", discreta_home, n);
		// fg.txt added within fg_io_open_channel_r_ascii() !
	sprintf(str1, "%s/lib/solvable_groups/ASCII/g_%ld_fg.txt", discreta_home, n);
	if (file_size(str1) <= 0) {
		printf("FG_OB::solvable_group_of_order_n_first() file %s does not exist !\n", str1);
		fflush(stdout);
		return FALSE;
		}
	io = (FG_IO *) my_malloc(sizeof(FG_IO), "FG_IO solvable_group_of_order_n_first");
	*fgio = io;
	fg_io_open_channel_r_ascii(&io->from, fname);
	if (!read_ascii(io->from.fp))
		return FALSE;
	return TRUE;
}

#if TEXDOCU
INT FG_OB::solvable_group_of_order_n_next(FG_IO **fgio)
#endif
{
	FG_IO *io = *fgio;
	
	if (!read_ascii(io->from.fp)) {
		fg_io_close_channel(&io->from);
		my_free(io);
		*fgio = NIL;
		return FALSE;
		}
	return TRUE;
}

INT fg_io_open_channel_r_ascii(FG_IO_CHANNEL *ch, BYTE *prefix)
{
	return fg_io_open_channel_rw_ascii(ch, prefix, "r");
}

INT fg_io_open_channel_w_ascii(FG_IO_CHANNEL *ch, BYTE *prefix)
{
	return fg_io_open_channel_rw_ascii(ch, prefix, "w");
}

INT fg_io_open_channel_rw_ascii(FG_IO_CHANNEL *ch, BYTE *prefix, BYTE *mode)
{
	ch->f_ascii = TRUE;
	ch->prefix = prefix;
	sprintf(ch->fname, "%sfg.txt", prefix);
	ch->fp = fopen(ch->fname, mode);
	printf("fg_io_open_channel_rw_ascii(): %s %s\n", mode, ch->fname);
	fflush(stdout);
	return OK;
}

INT fg_io_open_channel_r_db(FG_IO_CHANNEL *ch, BYTE *prefix)
{
	return fg_io_open_channel_rw_db(ch, prefix);
}

INT fg_io_open_channel_w_db(FG_IO_CHANNEL *ch, BYTE *prefix)
{
#if 0
	BYTE str[1024];
	
	sprintf(str, "%sfg.db", prefix);
	if (file_size(str) <= 0)
		do_db_fg_create_and_close(prefix);
	return fg_io_open_channel_rw_db(ch, prefix);
#endif
	return error("fg_io_open_channel_w_db(): not yet implemented");
}

INT fg_io_open_channel_rw_db(FG_IO_CHANNEL *ch, BYTE *prefix)
{
#if 0
	ch->f_ascii = FALSE;
	ch->prefix = prefix;
	sprintf(ch->fname, "%sfg.db", prefix);
	if (i_fg(&ch->db, prefix) != OK) {
		return error("fg_io_open_channel_w_db() error in i_fg()");
		}
#endif
	return error("fg_io_open_channel_rw_db(): not yet implemented");
	return OK;
}

INT fg_io_close_channel(FG_IO_CHANNEL *ch)
{
	if (ch->f_ascii) {
		fclose(ch->fp);
		}
	else {
#if 0
		e_fg(ch->db);
#endif
		return error("fg_io_close_channel() e_fg() not yet implemented");
		}
	return OK;
}


#endif /* SOLVABLE_TRUE */


