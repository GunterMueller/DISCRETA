// t148.C
//
// 09.12.1999
//
// computes TDO
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

INT Compute_tdo(BYTE *geo_file1, BYTE *geo_file2, INT f_tdo2, INT f_ddp, INT f_ddb, 
	INT f_range, INT range_f, INT range_l);
INT compute_tdo(INT f_tdo2, INT f_ddp, INT f_ddb, 
	MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	VECTOR_OP DDp, VECTOR_OP DDb);
INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, VECTOR_OP DDp, VECTOR_OP DDb);

void usage()
{
	printf("usage: t148.out [options] geo_file\n");
	printf("computes TDO of geometries. \n");
	printf("available options: \n");
	printf("-tdo2                       : compute second TDO\n");
	printf("-ddp                        : compute second point derivatives\n");
	printf("-ddb                        : compute second block derivatives\n");
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
	INT f_tdo2 = FALSE;
	INT f_ddp = FALSE;
	INT f_ddb = FALSE;
	INT f_range = FALSE, range_f, range_l;
	INT f_input_geo = FALSE;

	for (i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-tdo2") == 0) {
			f_tdo2 = TRUE;
			}
		else if (strcmp(argv[i], "-ddp") == 0) {
			f_ddp = TRUE;
			}
		else if (strcmp(argv[i], "-ddb") == 0) {
			f_ddb = TRUE;
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
		strcat(geo_file2, "_tdo.geo");
		Compute_tdo(geo_file1, geo_file2, f_tdo2, f_ddp, f_ddb, f_range, range_f, range_l);
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


INT Compute_tdo(BYTE *geo_file1, BYTE *geo_file2, INT f_tdo2, INT f_ddp, INT f_ddb, 
	INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT v = -1, b = -1;
	INT i, j, k, a, a1, l, geo_nr;
	BYTE geo_label[1000];
	BYTE *p_str, str[1024];
	VECTOR_OB labelling_P, labelling_B;
	// INT f_permute_points, f_permute_blocks;
	// PERMUTATION_OB p0, q0;
	MATRIX_OB I;
	VECTOR_OB row_decomp, col_decomp;
	VECTOR_OB DDp, DDb;
	INT nb_geo = 0;
	
	
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
		p_str = &buf[8];
		s_scan_int(&p_str, &geo_nr);
		geo_label[0] = 0;
		s_scan_token(&p_str, geo_label);
		printf("GEOMETRY %ld %s\n", geo_nr, geo_label);
		fflush(stdout);


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
				printf("v=%ld b=%ld\n", v, b);
				fflush(stdout);
				labelling_P.m_il(0);
				labelling_B.m_il(0);
				printf("labelling_P/B allocated\n");
				fflush(stdout);
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
				// I.Print();
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
			else if (strncmp(buf, "END", 3) == 0) {
				
				compute_tdo(f_tdo2, f_ddp, f_ddb, 
					&I, &labelling_P, &labelling_B, 
					&row_decomp, &col_decomp, 
					&DDp, &DDb);
				
				write_geo(out_fp, geo_nr, geo_label, &I, 
					&labelling_P, &labelling_B, 
					&row_decomp, &col_decomp, 
					f_ddp, f_ddb, &DDp, &DDb);
				printf("geo %ld %s written\n", geo_nr, geo_label);
				fflush(stdout);
				break;
				}
			}
		
		}
	
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}

INT compute_tdo(INT f_tdo2, INT f_ddp, INT f_ddb, 
	MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	VECTOR_OP DDp, VECTOR_OP DDb)
{
	PERMUTATION_OB pp, qq;
	MATRIX_OB tdo_scheme;
	INT v = I->s_hi();
	INT b = I->s_li();
	
	printf("computing TDO:\n");
	fflush(stdout);
	printf("initial decomposition:\n");
	printf("row_decomp:\n");
	row_decomp->println();
	printf("col_decomp:\n");
	col_decomp->println();
	I->calc_TDO(f_tdo2, f_ddp, f_ddb, 
		&pp, &qq, 
		row_decomp, col_decomp, 
		DDp, DDb, 
		&tdo_scheme, 
		FALSE /* f_v */, FALSE /* f_vv */);
	printf("t148: TDO finished!\n"); fflush(stdout);
	fflush(stdout);
	printf("pp=");
	pp.println();
	printf("qq=");
	qq.println();
	if (f_ddp) {
		printf("DDp=");
		DDp->println();
		}
	if (f_ddb) {
		printf("DDb=");
		DDb->println();
		}
	
	I->apply_perms(&pp, &qq);
	labelling_P->apply_perm(&pp);
	labelling_B->apply_perm(&qq);
#if 0
	if (f_ddp) {
		DDp->apply_perm_to_vector_of_pairs(&pp, v);
		}
	if (f_ddb) {
		DDb->apply_perm_to_vector_of_pairs(&qq, b);
		}
#endif

	printf("row_decomp=");
	row_decomp->println();
	printf("col_decomp=");
	col_decomp->println();
	return OK;
}

INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, VECTOR_OP DDp, VECTOR_OP DDb)
{
	INT v = I->s_hi();
	INT b = I->s_li();
	INT i, j, a, l;
	
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
	fprintf(fp, "END\n\n");
}

