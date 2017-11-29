// t153.C
//
// 13.12.1999
//
// merger for geometries
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
#include <DISCRETA/divs.h>
#include <DISCRETA/db.h>
#include <stdlib.h>

INT Merge(BYTE *geo_file1, BYTE *geo_file2, BYTE *base_fname, 
	INT f_range, INT range_f, INT range_l,
	INT f_v, INT f_vv);
INT merge(MATRIX_OP I, INT geo_nr, BYTE *geo_label, BYTE *geo_file2, 
	VECTOR_OP hash_table, VECTOR_OP hash_geo_nr, VECTOR_OP hash_geo_label, 
	INT f_v, INT f_vv);
INT equal_to(MATRIX_OP I, BYTE *geo_file2, INT geo_nr, BYTE *geo_label);
INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, 
	VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, VECTOR_OP DDp, VECTOR_OP DDb, 
	INT f_aut_gen, VECTOR_OP aut_gen);

void usage()
{
	printf("usage: t153.out [options] geo_file\n");
	printf("merger for geometries\n");
	printf("available options: \n");
	printf("-v                          : verbose\n");
	printf("-vv                         : very verbose\n");
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
	BYTE base_fname[1000];
	INT f_range = FALSE, range_f, range_l;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	INT f_input_geo = FALSE;

	for (i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			f_v = TRUE;
			}
		else if (strcmp(argv[i], "-vv") == 0) {
			f_vv = TRUE;
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
	strcpy(base_fname, geo_file2);

	if (f_input_geo) {
		strcat(geo_file2, "_merge.geo");
		Merge(geo_file1, geo_file2, base_fname, f_range, range_f, range_l, f_v, f_vv);
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


INT Merge(BYTE *geo_file1, BYTE *geo_file2, BYTE *base_fname, 
	INT f_range, INT range_f, INT range_l,
	INT f_v, INT f_vv)
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
	BYTE fname_hash[1000];
	VECTOR_OB hash_table, hash_geo_nr, hash_geo_label;
	INT nb_written = 0, nb_skip = 0;
	
	strcpy(fname_hash, base_fname);
	strcat(fname_hash, ".hash");
	printf("merging %s, result written to %s\nusing hash file %s\n", 
		geo_file1, geo_file2, fname_hash); fflush(stdout);
	
	in_fp = fopen(geo_file1, "r");
	
	// delete the output file:
	out_fp = fopen(geo_file2, "w");
	fclose(out_fp);
	
	hash_table.m_il(0);
	hash_geo_nr.m_il(0);
	hash_geo_label.m_il(0);
	
	
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
		printf("GEOMETRY %ld %s\n", geo_nr, geo_label);
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
				if (merge(&I, geo_nr, geo_label, geo_file2, 
					&hash_table, &hash_geo_nr, &hash_geo_label, 
					f_v, f_vv)) {
					printf("geo %ld %s is new\n", geo_nr, geo_label);
					out_fp = fopen(geo_file2, "a");
					write_geo(out_fp, geo_nr, geo_label, &I, 
						&labelling_P, &labelling_B, 
						&row_decomp, &col_decomp, 
						f_ddp, f_ddb, &DDp, &DDb, 
						f_aut_gen, &aut_gen);
					fclose(out_fp);
					printf("geo %ld %s written\n", geo_nr, geo_label);
					fflush(out_fp);
					nb_written++;
					}
				else {
					printf("already there, skipping\n");
					nb_skip++;
					fflush(stdout);
					}
				break;
				}
			}
		
		}
	// write hash file:
	{
	FILE *fp = fopen(fname_hash, "w");
	INT l = hash_table.s_li();
	INT i;
	BYTE *h, *label;
	INT nr;
	
	fprintf(fp, "%ld\n", l);
	for (i = 0; i < l; i++) {
		h = ((STRING_OP) hash_table.s_i(i))->s_str();
		nr = hash_geo_nr.s_ii(i);
		label = ((STRING_OP) hash_geo_label.s_i(i))->s_str();
		fprintf(fp, "%s %ld %s\n", h, nr, label);
		}
	fclose(fp);
	}
	printf("merge finished; found %ld geometries, skipped %ld geometries\n", 
		nb_written, nb_skip);
	printf("geometries written to %s\nhash file is %s\n", geo_file2, fname_hash);
	fflush(stdout);
	
	fclose(in_fp);
	return OK;
}

INT merge(MATRIX_OP I, INT geo_nr, BYTE *geo_label, BYTE *geo_file2, 
	VECTOR_OP hash_table, VECTOR_OP hash_geo_nr, VECTOR_OP hash_geo_label, 
	INT f_v, INT f_vv)
{
	STRING_OB h, label;
	BYTE hash_key[BTREEMAXKEYLEN + 1];
	INT key_len = BTREEMAXKEYLEN;
	INT l, idx, f_found;
	INT nr_cand;
	BYTE *label_cand;
	
	printf("hashing:\n");
	fflush(stdout);

	I->calc_hash_key(key_len, hash_key, f_v);
	h.init(hash_key);
	
	l = hash_table->s_li();
	hash_table->search(l, TRUE, &h, &idx, &f_found);
	if (f_found) {
		idx--;
		while (TRUE) {
			nr_cand = hash_geo_nr->s_ii(idx);
			label_cand = ((STRING_OP) hash_geo_label->s_i(idx))->s_str();
			if (equal_to(I, geo_file2, nr_cand, label_cand)) {
				printf("equal to %ld %s, skipping\n", nr_cand, label_cand);
				return FALSE;
				}
			else {
				printf("not equal to %ld %s\n", nr_cand, label_cand);
				}
			idx++;
			if (idx == l)
				break;
			if (h.sym_comp(hash_table->s_i(idx)) != 0)
				break;
			}
		}
	printf("adding new hash key at position %ld (new len = %ld)\n", idx, l + 1);
	hash_table->inc();
	hash_geo_nr->inc();
	hash_geo_label->inc();
	for (INT i = l; i > idx; i--) {
		hash_table->s_i(i)->swap(hash_table->s_i(i - 1));
		hash_geo_nr->s_i(i)->swap(hash_geo_nr->s_i(i - 1));
		hash_geo_label->s_i(i)->swap(hash_geo_label->s_i(i - 1));
		}
	h.swap((STRING_OP) hash_table->s_i(idx));
	hash_geo_nr->m_ii(idx, geo_nr);
	label.init(geo_label);
	label.swap((STRING_OP) hash_geo_label->s_i(idx));
	
	return TRUE;
}

INT equal_to(MATRIX_OP I, BYTE *geo_file2, INT geo_nr, BYTE *geo_label)
{
	FILE *fp;
	INT l;
	BYTE *p_str, str[1000];
	INT geo_nr1;
	BYTE geo_label1[1000];
	MATRIX_OB I1;
	INT v = I->s_hi();
	INT b = I->s_li();
	INT v1, b1;
	
	fp = fopen(geo_file2, "r");
	
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NIL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0; /* delete newline */
		if (buf[0] == '#')
			continue;
		
		if (strncmp(buf, "GEOMETRY", 8) != 0)
			continue;
		
		p_str = &buf[9];
		s_scan_int(&p_str, &geo_nr1);
		geo_label1[0] = 0;
		s_scan_token(&p_str, geo_label1);
		if (geo_nr1 == geo_nr && strcmp(geo_label1, geo_label) == 0) {
			
			while (TRUE) {
				if (fgets(buf, BUFSIZE, fp) == NIL)
					break;
				l = strlen(buf);
				if (l)
					buf[l - 1] = 0; /* delete newline */
				if (buf[0] == '#')
					continue;
			
				if (strncmp(buf, "v=", 2) == 0) {
					sscanf(buf, "v=%ld b=%ld", &v1, &b1);
					if (v1 != v || b1 != b) {
						fclose(fp);
						return FALSE;
						}
					}
				else if (strncmp(buf, "INCIDENCE_MATRIX", 16) == 0) {
					I1.m_ilih_n(b, v);
					for (INT i = 0; i < v; i++) {
						if (fgets(buf, BUFSIZE, fp) == NIL)
							return error("error reading INCIDENCE_MATRIX\n");
						for (INT j = 0; j < b; j++) {
							if (buf[j] == 'X') {
								I1.m_iji(i, j, 1);
								}
							}
						}
					printf("INCIDENCE_MATRIX of size %ld x %ld\n", 
						I1.s_hi(), I1.s_li());
					if (I1.sym_comp(I) == 0) {
						fclose(fp);
						return TRUE;
						}
					else {
						fclose(fp);
						return FALSE;
						}
					fflush(stdout);
					}
				else if (strncmp(buf, "END", 3) == 0) {
					return error("error: INCIDENCE_MATRIX not found\n");
					}
				
				}
			
			
			
			}
		}
	printf("compare_with() could not find GEOMETRY %ld %s in file %s\n", 
		geo_nr, geo_label, geo_file2);
	return error("error");
	fclose(fp);
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

