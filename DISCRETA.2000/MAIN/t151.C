// t151.C
//
// 12.12.1999
//
// computes TDA
// 
// 
//
// Anton Betten
// Bayreuth


#include <DISCRETA/discreta.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/geo.h>
#include <stdlib.h>

INT Compute_TDA(BYTE *geo_file1, BYTE *geo_file2, 
	INT f_range, INT range_f, INT range_l,
	INT f_reset_labels, INT f_v, INT f_vv);
INT compute_TDA(MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT *f_ddp, INT *f_ddb, VECTOR_OP DDp, VECTOR_OP DDb, 
	INT f_aut_gen, VECTOR_OP aut_gen, 
	INT f_v, INT f_vv);
INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, VECTOR_OP DDp, VECTOR_OP DDb, 
	INT f_aut_gen, VECTOR_OP aut_gen);

void usage()
{
	printf("usage: t151.out [options] geo_file\n");
	printf("computes the TDA\n");
	printf("requires generators for the automorphism group\n");
	printf("in the input file\n");
	printf("available options: \n");
	printf("-v                          : verbose\n");
	printf("-vv                         : very verbose\n");
	printf("-reset_labels               : verbose\n");
	printf("-range   first len          : range [first ... first - len - 1]\n");
}

int main(int argc, char **argv)
{
	if (argc < 2) {
		usage();
		exit(1);
		}
	
	discreta_init();
	{
	INT t0, t1, user_time;
	BYTE str[256], *p;
	INT i, l;
	BYTE *geo_file1;
	BYTE geo_file2[1000];
	INT f_range = FALSE, range_f, range_l;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	INT f_reset_labels = FALSE;
	INT f_input_geo = FALSE;

	for (i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			f_v = TRUE;
			}
		else if (strcmp(argv[i], "-vv") == 0) {
			f_vv = TRUE;
			}
		else if (strcmp(argv[i], "-reset_labels") == 0) {
			f_reset_labels = TRUE;
			}
		else if (strcmp(argv[i], "-range") == 0) {
			f_range = TRUE;
			i++;
			sscanf(argv[i], "%ld", &range_f);
			i++;
			sscanf(argv[i], "%ld", &range_l);
			}
		}
	geo_file1 = argv[argc - 1];
	
	t0 = os_ticks();

	l = strlen(geo_file1);
	if (l > 4 && strcmp(geo_file1 + l - 4, ".geo") == 0) {
		f_input_geo = TRUE;
		}
	else {
		return error("input file has unknown extension (should be .geo)");
		}

	
	strcpy(geo_file2, geo_file1);
	geo_file2[l - 4] = 0;

	if (f_input_geo) {
		strcat(geo_file2, "_tda.geo");
		Compute_TDA(geo_file1, geo_file2, 
			f_range, range_f, range_l, f_reset_labels, f_v, f_vv);
		}
	else  {
		return error("unknown input file format.");
		}
	printf("written file %s\n", geo_file2);
	




	
	
	t1 = os_ticks();
	
	user_time = t1 - t0;
	strcpy(str, "Running time for main(): ");
	print_delta_time(user_time, str);
	printf("%s\n", str);

	// my_malloc_dump();
	
	}
	discreta_exit();
	return 0;
}

#define BUFSIZE 50000

BYTE buf[BUFSIZE];


INT Compute_TDA(BYTE *geo_file1, BYTE *geo_file2, 
	INT f_range, INT range_f, INT range_l,
	INT f_reset_labels, INT f_v, INT f_vv)
{
	FILE *in_fp, *out_fp;
	INT v = -1, b = -1;
	INT i, j, k, a, a1, l, geo_nr, len;
	BYTE geo_label[1000];
	BYTE *p_str, str[1024];
	VECTOR_OB labelling_P, labelling_B;
	MATRIX_OB I;
	VECTOR_OB row_decomp, col_decomp;
	VECTOR_OB DDp, DDb;
	VECTOR_OB aut_gen;
	INT nb_geo = 0;
	INT f_ddp, f_ddb, f_aut_gen;
	
	
	printf("computing TDA for %s, result written to %s\n", geo_file1, geo_file2); fflush(stdout);
	in_fp = fopen(geo_file1, "r");
	out_fp = fopen(geo_file2, "w");
	
	
	
	while (TRUE) {
		if (fgets(buf, BUFSIZE, in_fp) == NIL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0; /* delete newline */
		if (buf[0] == '#')
			continue;
		
		if (strncmp(buf, "GEOMETRY", 8) != 0)
			continue;
		
		nb_geo++;
		p_str = &buf[9];
		s_scan_int(&p_str, &geo_nr);
		geo_label[0] = 0;
		s_scan_token(&p_str, geo_label);
		DDp.m_il(0);
		DDb.m_il(0);
		f_ddp = FALSE;
		f_ddb = FALSE;
		f_aut_gen = FALSE;
		
		if (f_range) {
			if (nb_geo < range_f || nb_geo >= range_f + range_l) {
				while (TRUE) {
					if (fgets(buf, BUFSIZE, in_fp) == NIL)
						break;
					l = strlen(buf);
					if (l)
						buf[l - 1] = 0; /* delete newline */
					if (buf[0] == '#')
						continue;
					if (strncmp(buf, "END", 3) == 0) {
						break;
						}
					}
				continue;
				}
			}


		while (TRUE) {
			if (fgets(buf, BUFSIZE, in_fp) == NIL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0; /* delete newline */
			if (buf[0] == '#')
				continue;
			
			if (strncmp(buf, "v=", 2) == 0) {
				sscanf(buf, "v=%ld b=%ld", &v, &b);
				labelling_P.m_il(0);
				labelling_B.m_il(0);
				}
			else if (strncmp(buf, "INCIDENCE_MATRIX", 16) == 0) {
				I.m_ilih_n(b, v);
				for (i = 0; i < v; i++) {
					if (fgets(buf, BUFSIZE, in_fp) == NIL)
						return error("error reading INCIDENCE_MATRIX\n");
					for (j = 0; j < b; j++) {
						if (buf[j] == 'X') {
							I.m_iji(i, j, 1);
							}
						}
					}
				row_decomp.m_il(1);
				row_decomp.m_ii(0, v);
				col_decomp.m_il(1);
				col_decomp.m_ii(0, b);
				printf("INCIDENCE_MATRIX of size %ld x %ld\n", 
					I.s_hi(), I.s_li());
				fflush(stdout);
				}
			else if (strncmp(buf, "LABELLING_OF_POINTS", 19) == 0) {
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading LABELLING_OF_POINTS");
				p_str = buf;
				labelling_P.m_il(v);
				for (i = 0; i < v; i++) {
					s_scan_int(&p_str, &a);
					labelling_P.m_ii(i, a);
					}
				printf("LABELLING_OF_POINTS: ");
				labelling_P.println();
				fflush(stdout);
				}
			else if (strncmp(buf, "LABELLING_OF_BLOCKS", 19) == 0) {
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading LABELLING_OF_BLOCKS");
				p_str = buf;
				labelling_B.m_il(b);
				for (i = 0; i < b; i++) {
					s_scan_int(&p_str, &a);
					labelling_B.m_ii(i, a);
					}
				printf("LABELLING_OF_BLOCKS: ");
				labelling_B.println();
				fflush(stdout);
				}
			else if (strncmp(buf, "DECOMPOSITION_OF_POINTS", 23) == 0) {
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading DECOMPOSITION_OF_POINTS");
				p_str = buf;
				s_scan_int(&p_str, &a);
				row_decomp.m_il(a);
				for (i = 0; i < a; i++) {
					s_scan_int(&p_str, &a1);
					row_decomp.m_ii(i, a1);
					}
				printf("DECOMPOSITION_OF_POINTS: ");
				row_decomp.println();
				fflush(stdout);
				}
			else if (strncmp(buf, "DECOMPOSITION_OF_BLOCKS", 23) == 0) {
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading DECOMPOSITION_OF_BLOCKS");
				p_str = buf;
				s_scan_int(&p_str, &a);
				col_decomp.m_il(a);
				for (i = 0; i < a; i++) {
					s_scan_int(&p_str, &a1);
					col_decomp.m_ii(i, a1);
					}
				printf("DECOMPOSITION_OF_BLOCKS: ");
				col_decomp.println();
				fflush(stdout);
				}
			else if (strncmp(buf, "DDP", 3) == 0) {
				f_ddp = TRUE;
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading DDP");
				p_str = buf;
				len = (v * (v - 1)) >> 1;
				DDp.m_il(len);
				for (i = 0; i < len; i++) {
					s_scan_int(&p_str, &a);
					DDp.m_ii(i, a);
					}
				printf("DDP: ");
				DDp.println();
				fflush(stdout);
				}
			else if (strncmp(buf, "DDB", 3) == 0) {
				f_ddb = TRUE;
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading DDB");
				p_str = buf;
				len = (b * (b - 1)) >> 1;
				DDb.m_il(len);
				for (i = 0; i < len; i++) {
					s_scan_int(&p_str, &a);
					DDb.m_ii(i, a);
					}
				printf("DDB: ");
				DDb.println();
				fflush(stdout);
				}
			else if (strncmp(buf, "AUT_GENS", 8) == 0) {
				PERMUTATION_OB pp;
				
				f_aut_gen = TRUE;
				if (fgets(buf, BUFSIZE, in_fp) == NIL)
					return error("error reading AUT_GENS");
				p_str = buf;
				s_scan_int(&p_str, &l);
				aut_gen.m_il(l);
				for (i = 0; i < l; i++) {
					if (fgets(buf, BUFSIZE, in_fp) == NIL)
						return error("error reading AUT_GENS");
					p_str = buf;
					pp.m_il(v);
					for (j = 0; j < v; j++) {
						s_scan_int(&p_str, &a1);
						pp.m_ii(j, a1);
						}
					pp.copy((PERMUTATION_OP) aut_gen.s_i(i));
					aut_gen.s_i(i)->println();
					}
				}
			else if (strncmp(buf, "END", 3) == 0) {
				compute_TDA(&I, 
					&labelling_P, &labelling_B, 
					&row_decomp, &col_decomp, 
					&f_ddp, &f_ddb, &DDp, &DDb, 
					f_aut_gen, &aut_gen, 
					f_v, f_vv);
				if (f_reset_labels) {
					INT v, b, i;
					v = I.s_hi();
					for (i = 0; i < v; i++) 
						labelling_P.m_ii(i, i + 1);
					b = I.s_li();
					for (i = 0; i < b; i++) 
						labelling_B.m_ii(i, i + 1);
					}
				write_geo(out_fp, geo_nr, geo_label, &I, 
					&labelling_P, &labelling_B, 
					&row_decomp, &col_decomp, 
					f_ddp, f_ddb, &DDp, &DDb, 
					f_aut_gen, &aut_gen);
				printf("geo %ld %s written\n", geo_nr, geo_label);
				fflush(out_fp);
				break;
				}
			}
		
		}
	
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}

INT compute_TDA(MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT *f_ddp, INT *f_ddb, VECTOR_OP DDp, VECTOR_OP DDb, 
	INT f_aut_gen, VECTOR_OP aut_gen, 
	INT f_v, INT f_vv)
{
	PERMUTATION_OB pp, qq;
	VECTOR_OB aut_gen2;
	INT v = I->s_hi();
	INT b = I->s_li();
	LABRA_OB A1, A2;
	VECTOR_OB P_tda_ordering, P_orbit_first, P_orbit_length, P_O_labra;
	VECTOR_OB B_tda_ordering, B_orbit_first, B_orbit_length, B_O_labra;
	INT i, a;
	
	printf("computing TDA:\n");
	fflush(stdout);
	
	if (!f_aut_gen) {
		printf("WARING: no generators for the automorphism group, skipping\n");
		return OK;
		}

	aut_gen->copy(&aut_gen2);
	vec_induce_action_on_columns(&aut_gen2, I);
	
	A1.init_quick(aut_gen);
	A2.init_quick(&aut_gen2);
	A1.orbits(&P_tda_ordering, &P_orbit_first, &P_orbit_length, 
		FALSE /* f_calc_labras */, &P_O_labra, FALSE /* f_v */);
	printf("P_orbits_first:\n");
	P_orbit_first.println();
	printf("P_orbit_length:\n");
	P_orbit_length.println();
	A2.orbits(&B_tda_ordering, &B_orbit_first, &B_orbit_length, 
		FALSE /* f_calc_labras */, &B_O_labra, FALSE /* f_v */);
	printf("B_orbits_first:\n");
	B_orbit_first.println();
	printf("B_orbit_length:\n");
	B_orbit_length.println();
	
	printf("tda ordering of points:\n");
	P_tda_ordering.println();
	pp.m_il(v);
	for (i = 0; i < v; i++) {
		a = P_tda_ordering.s_ii(i);
		pp.m_ii(a, i + 1);
		}
	printf("pp=\n");
	pp.println();

	printf("tda ordering of blocks:\n");
	B_tda_ordering.println();
	qq.m_il(b);
	for (i = 0; i < b; i++) {
		a = B_tda_ordering.s_ii(i);
		qq.m_ii(a, i + 1);
		}
	printf("qq=\n");
	qq.println();

	I->apply_perms(&pp, &qq);
	labelling_P->apply_perm(&pp);
	labelling_B->apply_perm(&qq);
	
	if (*f_ddp) {
		DDp->apply_perm_to_vector_of_pairs(&pp, v);
		}
	if (*f_ddb) {
		DDb->apply_perm_to_vector_of_pairs(&qq, b);
		}
	vec_conjugate(aut_gen, &pp);
	
	P_orbit_length.copy(row_decomp);
	B_orbit_length.copy(col_decomp);
		
	return OK;
}

INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, VECTOR_OP DDp, VECTOR_OP DDb, 
	INT f_aut_gen, VECTOR_OP aut_gen)
{
	INT v = I->s_hi();
	INT b = I->s_li();
	INT i, j, a, l;
	PERMUTATION_OP p;
	SYM_OB go;
	BYTE str[1000];
	LABRA_OB aut;
	
	fprintf(fp, "GEOMETRY %ld %s\n", nb_geo, label);
	fprintf(fp, "v=%ld b=%ld\n", v, b);
	fprintf(fp, "INCIDENCE_MATRIX\n");
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = I->s_iji(i, j);
			if (a == 0) {
				fprintf(fp, ".");
				}
			else
				fprintf(fp, "X");
			}
		fprintf(fp, "\n");
		}
	fprintf(fp, "LABELLING_OF_POINTS\n");
	for (i = 0; i < v; i++) {
		labelling_P->s_i(i)->fprint(fp);
		fprintf(fp, " ");
		}
	fprintf(fp, "\n");
	fprintf(fp, "LABELLING_OF_BLOCKS\n");
	for (j = 0; j < b; j++) {
		labelling_B->s_i(j)->fprint(fp);
		fprintf(fp, " ");
		}
	fprintf(fp, "\n");
	fprintf(fp, "DECOMPOSITION_OF_POINTS\n");
	l = row_decomp->s_li();
	fprintf(fp, "%ld ", l);
	for (i = 0; i < l; i++) {
		row_decomp->s_i(i)->fprint(fp);
		fprintf(fp, " ");
		}
	fprintf(fp, "\n");
	fprintf(fp, "DECOMPOSITION_OF_BLOCKS\n");
	l = col_decomp->s_li();
	fprintf(fp, "%ld ", l);
	for (j = 0; j < l; j++) {
		col_decomp->s_i(j)->fprint(fp);
		fprintf(fp, " ");
		}
	fprintf(fp, "\n");
	

	if (f_ddp) {
		fprintf(fp, "DDP\n");
		l = (v * (v - 1)) >> 1;
		for (j = 0; j < l; j++) {
			DDp->s_i(j)->fprint(fp);
			fprintf(fp, " ");
			}
		fprintf(fp, "\n");
		}
	if (f_ddb) {
		fprintf(fp, "DDB\n");
		l = (b * (b - 1)) >> 1;
		for (j = 0; j < l; j++) {
			DDb->s_i(j)->fprint(fp);
			fprintf(fp, " ");
			}
		fprintf(fp, "\n");
		}

	if (f_aut_gen) {
		aut.init_quick(aut_gen);
		aut.group_order(&go);
		str[0] = 0;
		go.sprint(str);
		fprintf(fp, "AUT_GENS (group order %s)\n", str);
		l = aut_gen->s_li();
		fprintf(fp, "%ld\n", l);
		for (i = 0; i < l; i++) {
			p = (PERMUTATION_OP) aut_gen->s_i(i);
			for (j = 0; j < v; j++) {
				a = p->s_ii(j);
				fprintf(fp, "%ld ", a);
				}
			fprintf(fp, "\n");
			}
		}
	//fprintf(fp, "\n");
	
	fprintf(fp, "END\n\n");
}

