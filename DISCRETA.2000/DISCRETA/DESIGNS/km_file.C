/* km_file.C */

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



#define BUFSIZE 10000

#if TEXDOCU
INT km_read_ascii(BYTE *KM_fname, MATRIX_OP M, INT *v, INT *t, INT *k, 
	VECTOR_OP G_gen, VECTOR_OP RR, VECTOR_OP MM, VECTOR_OP stab_go, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2)
#endif
{
	BYTE buf[BUFSIZE];
	FILE *fp;
	INT v1 = -1;
	INT t1 = -1, k1 = -1;
	INT m = -1, n = -1;
	INT i, j, a, b, l;
	INT f_compact = FALSE;
	VECTOR_OB B;

	printf("reading file %s of size %ld\n", KM_fname, file_size(KM_fname));
	fflush(stdout);
	fp = fopen(KM_fname, "r");
	if (fp == NIL)
		return error("km_read_ascii() cannot open file");

	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (buf[0] != '%')
			break;
		// printf("%s\n", buf);
		if (strcmp(buf, "% t, k:") == 0) {
			printf("reading t and k\n"); fflush(stdout);
			if (fgets(buf, BUFSIZE, fp) == NULL)
				return error("km_read_ascii() error reading t and k");
			sscanf(buf, "%% %ld %ld", &t1, &k1);
			}
		if (strcmp(buf, "% compact") == 0 || strcmp(buf, "%compact") == 0) {
			f_compact = TRUE;
			}
		if (strcmp(buf, "% m, n:") == 0) {
			printf("reading m and n\n"); fflush(stdout);
			if (fgets(buf, BUFSIZE, fp) == NULL)
				return error("km_read_ascii() error reading m and n");
			sscanf(buf, "%ld %ld", &m, &n);
			break;
			}
		if (strncmp(buf, "% degree:", 9) == 0) {
			printf("reading v\n"); fflush(stdout);
			sscanf(buf + 9, "%ld", &v1);
			}
		}
	if (t1 == -1)
		return error("km_read_ascii() t and k not found");
	if (m == -1)
		return error("km_read_ascii() m and n not found");
	if (v1 == -1)
		return error("km_read_ascii() v not found");
	if (f_compact) {
		printf("reading compact KM-matrix\n"); fflush(stdout);
		M->m_ilih_n(n, m);
		for (i = 0; i < m; i++) {
			fscanf(fp, "%ld", &l);
			for (j = 0; j < l; j++) {
				fscanf(fp, "%ld", &a);
				fscanf(fp, "%ld", &b);
				M->m_iji(i, a, b);
				}
			}
		}
	else {
		printf("reading KM-matrix\n"); fflush(stdout);
		M->m_ilih(n, m);
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				fscanf(fp, "%ld", &a);
				M->m_iji(i, j, a);
				}
			}
		}
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		// printf("%s\n", buf);
		if (strcmp(buf, "GENERATORS") == 0) {
			printf("reading generators\n"); fflush(stdout);
			G_gen->load_ascii(fp);
			}
		if (strcmp(buf, "SET-REPRESENTATIVES") == 0) {
			printf("reading set-representatives\n"); fflush(stdout);
			RR->load_ascii(fp);
			}
		if (strcmp(buf, "KM-MATRICES") == 0) {
			printf("reading KM-matrices\n"); fflush(stdout);
			MM->load_ascii(fp);
			}
		if (strcmp(buf, "STABILIZER-ORDERS") == 0) {
			printf("reading stabilizer orders\n"); fflush(stdout);
			stab_go->load_ascii(fp);
			}
		if (strcmp(buf, "ORBITS-BELOW") == 0) {
			printf("reading orbits below\n"); fflush(stdout);
			B.load_ascii(fp);
			B.s_i(0)->swap(Orbits_below1);
			B.s_i(1)->swap(Orbits_below2);
			}
		}
	*t = t1;
	*k = k1;
	*v = v1;
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT km_read_ascii_vtk(BYTE *KM_fname, INT *v, INT *t, INT *k)
#endif
{
	BYTE buf[BUFSIZE];
	FILE *fp;
	INT v1 = -1;
	INT t1 = -1, k1 = -1;
	INT m = -1, n = -1;
	INT l;

	printf("reading file %s of size %ld\n", KM_fname, file_size(KM_fname));
	fflush(stdout);
	fp = fopen(KM_fname, "r");
	if (fp == NIL)
		return error("km_read_ascii_vtk() cannot open file");

	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (buf[0] != '%')
			break;
		// printf("%s\n", buf);
		if (strcmp(buf, "% t, k:") == 0) {
			printf("reading t and k\n"); fflush(stdout);
			if (fgets(buf, BUFSIZE, fp) == NULL)
				return error("km_read_ascii_vtk() error reading t and k");
			sscanf(buf, "%% %ld %ld", &t1, &k1);
			}
		if (strcmp(buf, "% m, n:") == 0) {
			printf("reading m and n\n"); fflush(stdout);
			if (fgets(buf, BUFSIZE, fp) == NULL)
				return error("km_read_ascii_vtk() error reading m and n");
			sscanf(buf, "%ld %ld", &m, &n);
			break;
			}
		if (strncmp(buf, "% degree:", 9) == 0) {
			printf("reading v\n"); fflush(stdout);
			sscanf(buf + 9, "%ld", &v1);
			}
		}
	if (t1 == -1)
		return error("km_read_ascii_vtk() t and k not found");
	if (m == -1)
		return error("km_read_ascii_vtk() m and n not found");
	if (v1 == -1)
		return error("km_read_ascii_vtk() v not found");
	*t = t1;
	*k = k1;
	*v = v1;
	fclose(fp);
	return OK;
}

#if TEXDOCU
INT km_write_ascii_G_gen(BYTE *km_fname, VECTOR_OP G_gen)
#endif
{
	FILE *fp;
	
	printf("appending generators for G to KM-file %s\n", km_fname); fflush(stdout);
	fp = fopen(km_fname, "a");
	fprintf(fp, "GENERATORS\n");
	G_gen->save_ascii(fp);
	fclose(fp);
	printf("km file size %ld\n", file_size(km_fname)); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_write_ascii_representatives(BYTE *km_fname, VECTOR_OP RR)
#endif
{
	FILE *fp;

	printf("appending set-representatives to KM-file %s\n", km_fname); fflush(stdout);
	fp = fopen(km_fname, "a");
	fprintf(fp, "SET-REPRESENTATIVES\n");
	RR->save_ascii(fp);
	fclose(fp);
	printf("km file size %ld\n", file_size(km_fname)); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_write_ascii_stab_go(BYTE *km_fname, VECTOR_OP stab_go)
#endif
{
	FILE *fp;

	printf("appending stabilizer orders to KM-file %s\n", km_fname); fflush(stdout);
	fp = fopen(km_fname, "a");
	fprintf(fp, "STABILIZER-ORDERS\n");
	stab_go->save_ascii(fp);
	fclose(fp);
	printf("km file size %ld\n", file_size(km_fname)); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_write_stab_go_k_sets(BYTE *km_fname, VECTOR_OP stab_go, INT k)
#endif
{
	FILE *fp;
	INT i, len;
	SYM_OP p;
	BYTE str[1000];
	VECTOR_OP p_stab_go;

	printf("appending stabilizer orders k-sets to KM-file %s\n", km_fname); fflush(stdout);
	fp = fopen(km_fname, "a");
	fprintf(fp, "STABILIZER-ORDER-K-SETS\n");
	p_stab_go = (VECTOR_OP) stab_go->s_i(k);
	len = p_stab_go->s_li();
	for (i = 0; i < len; i++) {
		p = p_stab_go->s_i(i);
		str[0] = 0;
		p->sprint(str);
		fprintf(fp, "%s ", str);
		}
	fprintf(fp, "\n");
	
	fclose(fp);
	printf("km file size %ld\n", file_size(km_fname)); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_write_ascii_KM_matrices(BYTE *km_fname, VECTOR_OP MM)
#endif
{
	FILE *fp;

	printf("appending KM-file (t,t+1) to KM-file %s\n", km_fname); fflush(stdout);
	fp = fopen(km_fname, "a");
	fprintf(fp, "KM-MATRICES\n");
	MM->save_ascii(fp);
	fclose(fp);
	printf("km file size %ld\n", file_size(km_fname)); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_write_ascii_below(BYTE *km_fname, VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2)
#endif
{
	FILE *fp;
	VECTOR_OB B;
	
	B.m_il(2);
	Orbits_below1->swap((VECTOR_OP) B.s_i(0));
	Orbits_below2->swap((VECTOR_OP) B.s_i(1));

	printf("appending orbits_below to km_file %s\n", km_fname); fflush(stdout);
	fp = fopen(km_fname, "a");
	fprintf(fp, "ORBITS-BELOW\n");
	B.save_ascii(fp);
	Orbits_below1->swap((VECTOR_OP) B.s_i(0));
	Orbits_below2->swap((VECTOR_OP) B.s_i(1));
	fclose(fp);
	printf("km file size %ld\n", file_size(km_fname)); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_print_M_asc(BYTE *g_label, BYTE *g_label_tex, BYTE *km_fname, 
	VECTOR_OP G_gen, SYM_OP go, INT deg, 
	MATRIX_OP M, INT gl_t, INT gl_k, INT f_k2, INT k2)
#endif
{
	BYTE str1[10000];
	FILE *fp;
	BYTE *fname = km_fname;
	INT i, m, n, l;

	fp = fopen(fname, "w");
	str1[0] = 0;
	go->sprint(str1);
	fprintf(fp, "%% this file:    %s\n", km_fname);
	fprintf(fp, "%% group:        %s \n", g_label);
	fprintf(fp, "%% group:        %s \n", g_label_tex);
	fprintf(fp, "%% order:        %s \n", str1);
	fprintf(fp, "%% degree:       %ld \n", deg);
	l = G_gen->s_li();
	fprintf(fp, "%% # generators: %ld \n", l);
	for (i = 0; i < l; i++) {
		str1[0] = 0;
		G_gen->s_i(i)->sprint(str1);
		fprintf(fp, "%% %s\n", str1);
		}
	fflush(fp);
	
	if (f_k2) {
		fprintf(fp, "%% 2 values of k !\n");
		fprintf(fp, "%% t, k1, k2:\n");
		fprintf(fp, "%% %ld %ld %ld\n", gl_t, gl_k, k2);
		}
	else {
		fprintf(fp, "%% t, k:\n");
		fprintf(fp, "%% %ld %ld\n", gl_t, gl_k);
		}
	fflush(fp);

	if (M) {
		m = M->s_hi();
		n = M->s_li();

		fprintf(fp, "%% m, n:\n");
		fprintf(fp, "%ld %ld\n", m, n);
		fflush(fp);
		M->fprint_raw(fp);
		}
	
	fclose(fp);

	printf("written file %s of size %ld\n", km_fname, file_size(km_fname));
	fflush(stdout);
	return TRUE;
}

#if TEXDOCU
INT km_init_file_names(BYTE *g_label, BYTE *km_fname, INT t, INT k)
#endif
{
	FILE *fp;
	
	sprintf(km_fname, "KM_");
	strcat(km_fname, g_label);
	sprintf(km_fname + strlen(km_fname), "_t%ld_k%ld", t, k);
	strcat(km_fname, ".txt");
	
	fp = fopen("km_fname", "w");
	fputs(km_fname, fp);
	fclose(fp);

	printf("%s\n", km_fname); fflush(stdout);
	return OK;
}

#if TEXDOCU
INT km_compute_TDO_decomposition(BYTE *g_label, INT t, INT k, MATRIX_OP Mtk)
#else
This function computes the TDO decomposition of the Kramer-Mesner matrix.
The function multivalued\_geo\_tdo() from geo\_data.C is called.
#endif
{
#ifdef SYM_GEO
	MATRIX_OB Mtk2;
	BYTE str[1024];
	FILE *fp;
	INT i, m, n, a, l;
	PERMUTATION_OB p, q;
	VECTOR_OB row_decomp, col_decomp;
	
	sprintf(str, "KM_%s_t%ld_k%ld.txt", g_label, t, k);
	multivalued_geo_tdo(Mtk, &Mtk2, &p, &q, &row_decomp, &col_decomp);
	printf("multivalued TDO computed\n"); fflush(stdout);
	// Mtk2.Print();
	fp = fopen(str, "a");
	m = Mtk2.s_hi();
	n = Mtk2.s_li();
	fprintf(fp, "%%\n");
	fprintf(fp, "%% TDO decomposition:\n");
	fprintf(fp, "%ld %ld\n", m, n);
	Mtk2.fprint_raw(fp);
	fprintf(fp, "%%\n");
	fprintf(fp, "%% row decomposition:\n");
	l = row_decomp.s_li();
	fprintf(fp, "%% ");
	for (i = 0; i < l; i++) {
		a = row_decomp.s_ii(i);
		fprintf(fp, "%ld ", a);
		}
	fprintf(fp, "\n");

	fprintf(fp, "%%\n");
	fprintf(fp, "%% column decomposition:\n");
	l = col_decomp.s_li();
	fprintf(fp, "%% ");
	for (i = 0; i < l; i++) {
		a = col_decomp.s_ii(i);
		fprintf(fp, "%ld ", a);
		}
	fprintf(fp, "\n");

	fprintf(fp, "%%\n");
	fprintf(fp, "%% row permutation:\n");
	l = p.s_li();
	fprintf(fp, "%% ");
	for (i = 0; i < l; i++) {
		a = p.s_ii(i);
		fprintf(fp, "%ld ", a);
		}
	fprintf(fp, "\n");

	fprintf(fp, "%%\n");
	fprintf(fp, "%% col permutation:\n");
	l = q.s_li();
	fprintf(fp, "%% ");
	for (i = 0; i < l; i++) {
		a = q.s_ii(i);
		fprintf(fp, "%ld ", a);
		}
	fprintf(fp, "\n");

	fclose(fp);
	return OK;
#else
	return error("km_compute_TDO_decomposition(): SYM_GEO not available !");
#endif /* SYM_GEO */
}
	
#if TEXDOCU
void get_solutions(BYTE *KM_fname, INT lambda)
#else
This function tries to read the file "solutions" which contains 01-solutions 
of KM-systems produced for instance with the LLL-solver of 
A. Wassermann. If such a file is found, its contents is  
appended to the KM-file (whose filename is KM\_fname). 
It is surrounded by a line 
//PRE
LAMBDA 42
///PRE
at the top (here 42 stands for the chosen value of lambda) and 
//PRE
LAMBDAEND 42
///PRE
at the bottom.
#endif
{
	BYTE buf[BUFSIZE];
	INT v, t, k, m, n;
	MATRIX_OB Mtk;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	VECTOR_OB RHS;
	INT l;
	FILE *fp, *fp_sol;

	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	printf("KM-file: %s lambda = %ld\n", KM_fname, lambda);
	fflush(stdout);
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	
	fp_sol = fopen("solutions", "r");
	if (fp_sol == NIL) {
		error("get_solutions() cannot open file solutions");
		return;
		}

	fp = fopen(KM_fname, "a");
	fprintf(fp, "LAMBDA %ld\n", lambda);
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp_sol) == NULL)
			break;
		l = strlen(buf);
		if (l) 
			buf[l - 1] = 0;
		while (l > 0 && !(buf[l - 1] == '0' || buf[l - 1] == '1')) {
			buf[l - 1] = 0;
			l--;
			}
		if (l != n) {
			printf("WARNING: the solutions has wrong length "
				"(has: %ld should have %ld\n", l, n);
			fflush(stdout);
			}
		// printf("%s\n", buf);
		fprintf(fp, "%s\n", buf);
		}
	fprintf(fp, "LAMBDAEND %ld\n\n", lambda);
	

	fclose(fp);
	fclose(fp_sol);

}

#if TEXDOCU
INT km_nb_of_solutions(BYTE *KM_fname, INT lambda)
#endif
{
	FILE *fp_KM;
	BYTE buf[BUFSIZE];
	INT l, lamb, no;
	
	fp_KM = fopen(KM_fname, "r");
	no = -1;
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp_KM) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strncmp(buf, "LAMBDA", 6) != 0)
			continue;
		sscanf(buf, "LAMBDA %ld", &lamb);
		if (lamb != lambda) {
			km_read_until_lambdaend(fp_KM);
			continue;
			}
		no = 0;
		while (TRUE) {
			if (fgets(buf, BUFSIZE, fp_KM) == NULL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0;
			if (strncmp(buf, "LAMBDAEND", 9) == 0) {
				break;
				}
			no++;
			}
		break;
		} // while 
	fclose(fp_KM);
	if (no == -1) {
		printf("km_nb_of_solutions(): "
			"no soutions for lambda=%ld found in file %s\n", 
			lambda, KM_fname);
		return error("");
		}
	return no;
	
}

INT km_read_until_lambdaend(FILE *fp)
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
INT km_get_solutions(BYTE *KM_fname, 
	INT lambda, INT from, INT len, VECTOR_OP S)
#endif
{
	FILE *fp_KM;
	BYTE buf[BUFSIZE];
	INT l, lamb, no, i;
	VECTOR_OB s;
	
	S->m_il(len);
	fp_KM = fopen(KM_fname, "r");
	no = -1;
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp_KM) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strncmp(buf, "LAMBDA", 6) != 0)
			continue;
		sscanf(buf, "LAMBDA %ld", &lamb);
		if (lamb != lambda) {
			km_read_until_lambdaend(fp_KM);
			continue;
			}
		no = 0;
		while (TRUE) {
			if (fgets(buf, BUFSIZE, fp_KM) == NULL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0;
			if (strncmp(buf, "LAMBDAEND", 9) == 0) {
				break;
				}
			if (no >= from && no < from + len) {
				for (i = 0; i < l; i++) {
					if (!(buf[i] == '1' || buf[i] == '0'))
						break;
					}
				l = i;
				s.m_il_n(l);
				for (i = 0; i < l; i++) {
					if (buf[i] == '1')
						s.m_ii(i, 1);
					}
				s.swap(S->s_i(no - from));
				}
			no++;
			}
		break;
		} // while 
	if (no == -1) {
		printf("no soutions for lambda=%ld found in file %s\n", 
			lambda, KM_fname);
		return error("");
		}
		
	fclose(fp_KM);
	return OK;
	
}


#if TEXDOCU
void check_solutions(BYTE *KM_fname, INT lambda)
#else
This function searches for a LAMBDA section in the KM-file 
with appropriate lambda (specified as argument). 
If it finds such a section, the solutions are checked. 
There is the possibility to check intersection equations which 
is deactivated at the moment.
#endif
{
	BYTE buf[BUFSIZE];
	INT v, t, k, m, n;
	MATRIX_OB Mtk, X, Y;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	INT i, a, l, no;
	FILE *fp;
	BYTE str[1024];
	VECTOR_OB Eqns, RHS;
	MATRIX_OB Mendel, LAmbda;
	VECTOR_OB block_types, bai_ref, bai_refv, bai_data;
	INTEGER_OB bai_len;
	INT idx;

	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	printf("KM-file: %s lambda = %ld\n", KM_fname, lambda);
	fflush(stdout);
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	
	sprintf(str, "LAMBDA %ld", lambda);
	fp = fopen(KM_fname, "r");
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strcmp(buf, str) != 0)
			continue;
		printf("found solutions, checking...\n");
		no = 0;
		while (TRUE) {
			if (fgets(buf, BUFSIZE, fp) == NULL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0;
			if (strncmp(buf, "LAMBDAEND", 9) == 0)
				break;
			no++;
			// printf("checking solution %ld: %d\n", no, buf);
			l = strlen(buf);
			if (l != n) {
				printf("l = %ld != n = %ld\n", l, n);
				continue;
				}
			X.m_ilih(1, n);
			for (i = 0; i < n; i++) {
				if (buf[i] == '1')
					X.m_iji(i, 0, 1);
				else
					X.m_iji(i, 0, 0);
				}
			Mtk.mult(&X, &Y);
			for (i = 0; i < m; i++) {
				a = Y.s_iji(i, 0);
				if (a != lambda)
					printf("WARNING !!! Y[%ld] = %ld != lambda = %ld\n", i, a, lambda);
				}
			if (no % 10 == 0) {
				printf(",");
				}
			else
				printf(".");
			if (no % 50 == 0)
				printf("%ld\n", no);
			fflush(stdout);
			}
		printf("\nfinished with solution check for lambda = %ld\n", lambda);

		}

	fclose(fp);
	

}

#if TEXDOCU
void do_mckay(BYTE *KM_fname, INT lambda)
#else
This is the DISCRETA interface to B.D. McKay-s equation solver. 
This solver is a backtrack solver with some intelligent tests 
during the processing. The equations are limited in size 
which is due to a constant in the McKay program. 
Also the solver is not practical for too large equations. 
It proves to be useful to show the nonexistence of solutions, 
especially if there are only few equations, i.e. if the 
Kramer-Mesner matrix has only few rows.
#endif
{
	BYTE s[1000];
	char *logfile = "mckay.log";

	if (lambda == 0) {
		printf("lambda not set !\n");
		return;
		}
	sprintf(s, LOC_MCKAY);
	sprintf(s + strlen(s), " %s %ld <%s", 
		logfile, lambda, KM_fname);
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);
	printf("finished !\n");
	fflush(stdout);
	
}

#if TEXDOCU
void do_spread(INT f_silent, BYTE *KM_fname, INT lambda)
#endif
{
	BYTE s[1024];

	sprintf(s, LOC_SPREAD);

	sprintf(s + strlen(s), " %ld ", lambda);

	if (f_silent) {
		sprintf(s + strlen(s), " silent ");
		}
	
	sprintf(s + strlen(s), " %s", KM_fname);

	sprintf(s + strlen(s), " | tee spread.log");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	printf("finished !\n");
	fflush(stdout);
	
}

#if TEXDOCU
void do_dance(INT f_silent, BYTE *KM_fname, INT lambda)
#endif
{
	BYTE s[1024];

	sprintf(s, "stripcolumnsdance %ld %s | %s ", lambda, KM_fname, LOC_DANCE);

	if (f_silent) {
		sprintf(s + strlen(s), " silent ");
		}
	
	sprintf(s + strlen(s), " | tee dance.log");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	printf("finished !\n");
	fflush(stdout);
	
}

#if TEXDOCU
void do_LLL(INT f_silent, BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT lambda, INT f_iterate, INT nb_iterate)
#else
This is the DISCRETA interface to Alfred Wassermanns LLL-equation solver. 
The flag f\_with determines if the solutions should be stored in a 
file "solutions". The text output is captured in a file "lll.log".
Afterwards, the last 10 lines of this file are copied into a 
file "lll.log.1". This file is useful for a quick inspection 
of the result of the solver (i.e. number of solutions). 
It is also very important to check that the solver really has terminated 
correctly so the the set of solutions is complete. 
Especially during large computations, some users terminate by hand 
using ctrl-c. In this case, one gets only a lower bound for 
the correct number of solutions of the system.
#endif
{
	char s[1024];

	/* printf("km_lll_callback\n"); */
	sprintf(s, LOC_LLL);

	sprintf(s + strlen(s), " %ld %ld ", 
		(INT) (c0_factor * lambda), 
		lambda);
		
	if (f_silent) {
		sprintf(s + strlen(s), " silent ");
		}
		
	if (f_iterate) {
		sprintf(s + strlen(s), " iterate %ld ", nb_iterate);
		}
	else {
		sprintf(s + strlen(s), " bkz %ld %ld ", beta, p);
		}

	sprintf(s + strlen(s), " %s", KM_fname);

#if 0
	if (f_selection) {
		sprintf(s + strlen(s), " %ld", selection);
		}
#endif
	sprintf(s + strlen(s), " | tee lll.log");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	sprintf(s, "tail -10 lll.log >lll.log.1");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	printf("finished !\n");
	fflush(stdout);
	
}

#if TEXDOCU
void do_largeset1(BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT N)
#endif
{
	char s[1024];

	sprintf(s, LOC_LARGESET1);
	sprintf(s + strlen(s), " %s %ld %ld %ld %ld", 
		KM_fname, N, 
		c0_factor, beta, p);

	sprintf(s + strlen(s), " | tee ls1.log");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	sprintf(s, "tail -10 ls1.log >ls1.log.1");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	printf("finished !\n");
	fflush(stdout);
	
}

#if TEXDOCU
void do_largeset2(BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT max_rounds, INT N, 
	INT max_solutions, INT max_loops)
#endif
{
	char s[1024];

	sprintf(s, LOC_LARGESET2);
	sprintf(s + strlen(s), " %s %ld %ld %ld %ld %ld %ld %ld", 
		KM_fname, N, max_solutions, max_loops, max_rounds, 
		c0_factor, beta, p);

	sprintf(s + strlen(s), " | tee ls2.log");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	sprintf(s, "tail -10 ls2.log >ls2.log.1");
	printf("calling ...\n%s\n", s);
	fflush(stdout);
	system(s);

	printf("finished !\n");
	fflush(stdout);
	
}

INT show_km_matrix(BYTE *KM_fname)
{
	FILE *fp;
	char buf[BUFSIZE];

	fp = fopen(KM_fname, "r");
	if (fp == NIL) {
		printf("cannot open file %s!\n", KM_fname);
		fflush(stdout);
		return error("error");
		}
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NIL)
			break;
		printf("%s", buf);
		}
	fclose(fp);
	return OK;
}




#endif /* LADDER_TRUE */


