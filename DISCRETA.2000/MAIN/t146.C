// t146.C
//
// 09.12.1999
//
// 
// converts from .shrink / .base_blocks / .blocks / .inc  (or .txt)  to .geo
// and from .geo to .inc
// Examples for the various file formats are given below
// 
//
// Anton Betten
// Bayreuth




#include <DISCRETA/discreta.h>
#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/geo.h>
#include <DISCRETA/ladder.h> // for dc_calc_set_orbit()
#include <stdlib.h>

INT convert_shrink_geo(BYTE *shrink_file, BYTE *geo_file, BYTE *base_fname, INT f_range, INT range_f, INT range_l);
INT get_col(BYTE *alphabet, BYTE c);
INT convert_blocks_geo(BYTE *base_file, BYTE *geo_file, INT f_range, INT range_f, INT range_l);
INT convert_base_blocks_geo(BYTE *base_file, BYTE *generators_file, 
	BYTE *geo_file, INT f_range, INT range_f, INT range_l);
INT convert_inc_geo(BYTE *inc_file, BYTE *geo_file, BYTE *base_fname, 
	INT f_range, INT range_f, INT range_l);
INT convert_geo_inc(BYTE *geo_file, BYTE *inc_file, INT f_range, INT range_f, INT range_l);
INT write_inc(FILE *fp, MATRIX_OP I);
INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, VECTOR_OP labelling_P, VECTOR_OP labelling_B);

void usage()
{
	printf("usage: t146.out [options] [generator_file] input_file\n");
	printf("converts .shrink / .base_blocks / .blocks / .inc to .geo\n");
	printf("or .geo to .inc\n");
	printf("available options: \n");
	printf("-range   first len          : range [first ... first - len - 1]\n");
	printf("the file of generators is needed only if the input is of type .base_blocks\n");
}

int main(int argc, char **argv)
{
	if (argc <= 1) {
		usage();
		exit(1);
		}
	
	discreta_init();
	{
	INT t0, t1, user_time;
	BYTE str[256], *p;
	INT i, l;
	BYTE *input_file;
	BYTE base_fname[1000];
	BYTE output_file[1000];
	INT f_range = FALSE, range_f, range_l;
	INT f_input_shrink = FALSE;
	INT f_input_base_blocks = FALSE;
	INT f_input_blocks = FALSE;
	INT f_input_inc = FALSE;
	INT f_input_geo = FALSE;
	
	for (i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-range") == 0) {
			f_range = TRUE;
			i++;
			sscanf(argv[i], "%ld", &range_f);
			i++;
			sscanf(argv[i], "%ld", &range_l);
			}
		}
	input_file = argv[argc - 1];
	
	t0 = os_ticks();

	l = strlen(input_file);
	
	strcpy(output_file, input_file);

	if (l > 4 && strcmp(input_file + l - 4, ".inc") == 0) {
		f_input_inc = TRUE;
		printf("input is a .inc file\n");
		output_file[l - 4] = 0;
		}
	else if (l > 4 && strcmp(input_file + l - 4, ".txt") == 0) {
		f_input_inc = TRUE;
		printf("input is a .txt file\n");
		output_file[l - 4] = 0;
		}
	else if (l > 4 && strcmp(input_file + l - 4, ".geo") == 0) {
		f_input_geo = TRUE;
		printf("input is a .geo file\n");
		output_file[l - 4] = 0;
		}
	else if (l > 12 && strcmp(input_file + l - 12, ".base_blocks") == 0) {
		f_input_base_blocks = TRUE;
		printf("input is a .base_blocks file\n");
		output_file[l - 12] = 0;
		}
	else if (l > 7 && strcmp(input_file + l - 7, ".blocks") == 0) {
		f_input_blocks = TRUE;
		printf("input is a .blocks file\n");
		output_file[l - 7] = 0;
		}
	else if (l > 7 && strcmp(input_file + l - 7, ".shrink") == 0) {
		f_input_shrink = TRUE;
		printf("input is a .shrink file\n");
		output_file[l - 7] = 0;
		}
	else {
		return error("input file has unknown extension (should be .shrink or .blocks or .base_blocks or .inc or .geo)");
		}
	fflush(stdout);
	strcpy(base_fname, output_file);

	if (f_input_inc) {
		strcat(output_file, ".geo");
		convert_inc_geo(input_file, output_file, base_fname, f_range, range_f, range_l);
		}
	else if (f_input_blocks) {
		strcat(output_file, ".geo");
		convert_blocks_geo(input_file, output_file, f_range, range_f, range_l);
		}
	else if (f_input_base_blocks) {
		BYTE *group_generators_file = argv[argc - 2];
		strcat(output_file, ".geo");
		convert_base_blocks_geo(input_file, group_generators_file, output_file, f_range, range_f, range_l);
		}
	else if (f_input_shrink) {
		strcat(output_file, ".geo");
		convert_shrink_geo(input_file, output_file, base_fname, f_range, range_f, range_l);
		}
	else  {
		strcat(output_file, ".inc");
		convert_geo_inc(input_file, output_file, f_range, range_f, range_l);
		}
	printf("written file %s\n", output_file);
	
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

INT convert_shrink_geo(BYTE *shrink_file, BYTE *geo_file, BYTE *base_fname, INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT nb_geo = 0;
	INT nrow, ncol, nb_X, nbX;
	INT a, i, j, ii, last_i;
	INT l, ll;
	BYTE *p;
	BYTE *alphabet = NIL;
	INT line, next_line_first_X, col;
	INT f_equal;
	BYTE s[64], s_line[1024], c, cc;
	BYTE **old_lines;
	INT f_v = FALSE;
	MATRIX_OB I;
	VECTOR_OB labelling_P, labelling_B;
	INT v, b;
	
	in_fp = fopen(shrink_file, "r");
	out_fp = fopen(geo_file, "w");



	if (fgets(buf, BUFSIZE, in_fp) == NIL)
		return error("empty file");
	
	sscanf(buf, "%ld %ld %ld", &nrow, &ncol, &nb_X);
	printf("reading geometries of size %ld %ld with nb_X=%ld\n", nrow, ncol, nb_X);
	
	alphabet = (BYTE *) my_malloc((ncol + 11) * sizeof(BYTE), "alphabet");
	
	for (i = 0; i < 10; i++)
		alphabet[i] = '0' + i;
	for (i = 10; i < ncol; i++) {
		j = i - 10;
		if (j < 26)
			alphabet[i] = 'a' + j;
		else {
			j -= 26;
			alphabet[i] = 'A' + j;
			}
		}
	alphabet[ncol] = 0;
	printf("alphabet: %s\n", alphabet);
	old_lines = (BYTE **) my_malloc(nrow * sizeof(BYTE *), "old_lines");
	for (i = 0; i < nrow; i++) {
		old_lines[i] = (BYTE *) my_malloc((ncol + 1) * sizeof(BYTE), "line");
		old_lines[i][0] = 0;
		}
	
	while (TRUE) {
		if (fgets(buf, BUFSIZE, in_fp) == NIL)
			break;
		if (buf[0] == '-')
			break;
		if (f_v) {
			printf("nb_geo = %ld: reading %s\n", nb_geo, buf);
			fflush(stdout);
			}
		I.m_ilih_n(ncol, nrow);
		nb_geo++;
		
		l = strlen(buf);
		last_i = 0;
		line = 0;
		nbX = 0;
		for (i = 0; i < l; i++) {
			if (buf[i] == 0)
				return error("primature end of line !");
			if (buf[i] == ',' || buf[i] == ';') {
				cc = buf[i];
				buf[i] = 0;
				strcpy(s_line, buf + last_i);
				if (strlen(s_line) == 0) {
					strcpy(s_line, old_lines[line]);
					}
				ll = strlen(s_line);
				if (f_v) {
					printf("line %ld: %s\n", line, s_line);
					fflush(stdout);
					}
				ll = strlen(s_line);
				for (ii = 0; ii < ll; ii++) {
					c = s_line[ii];
					col = get_col(alphabet, c);
					// printf("col = %ld\n", col);
					// X[nbX++] = line * ncol + col;
					I.m_iji(line, col, 1);
					nbX++;
					}
				strcpy(old_lines[line], s_line);
				last_i = i + 1;
				line++;
				if (line == nrow) {
					if (cc != ';')
						return error("last line not terminated by ';'!");
					break;
					}
				}
			}
		if (nbX != nb_X) {
			printf("nbX = %ld nb_X = %ld\n", nbX, nb_X);
			fflush(stdout);
			return error("nbX != nbX");
			}
		
		v = nrow;
		b = ncol;
		
		labelling_P.m_il(v);
		for (i = 0; i < v; i++) {
			labelling_P.m_ii(i, i + 1);
			}
		labelling_B.m_il(b);
		for (j = 0; j < b; j++) {
			labelling_B.m_ii(j, j + 1);
			}


		write_geo(out_fp, nb_geo, base_fname, &I, &labelling_P, &labelling_B);
		
		
		}

	for (i = 0; i < nrow; i++) {
		my_free(old_lines[i]);
		}
	my_free(old_lines);
	my_free(alphabet);



	
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}

INT get_col(BYTE *alphabet, BYTE c)
{
	INT i, l;

	l = strlen(alphabet);
	for (i = 0; i < l; i++) {
		if (alphabet[i] == c)
			return i;
		}
	return error("unknown character !");
}

INT convert_blocks_geo(BYTE *base_file, BYTE *geo_file, INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT v, b;
	INT i, j, k, a, l;
	BYTE *p_str, str[1024];
	VECTOR_OB labelling_P, labelling_B;
	MATRIX_OB I;
	INT nb_geo = 0;
	
	
	in_fp = fopen(base_file, "r");
	out_fp = fopen(geo_file, "w");
	
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
		printf("GEOMETRY %ld \n", nb_geo);
		
		
		j = -1;
		v = -1;
		b = -1;
		k = -1;
		while (TRUE) {
			if (fgets(buf, BUFSIZE, in_fp) == NIL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0; /* delete newline */
			if (buf[0] == '#')
				continue;
		
			if (strncmp(buf, "v=", 2) == 0) {
				p_str = &buf[2];
				sscanf(buf, "v=%ld b=%ld", &v, &b);
				I.m_ilih_n(b, v);
				labelling_P.m_il(v);
				for (i = 0; i < v; i++) {
					labelling_P.m_ii(i, i + 1);
					}
				labelling_B.m_il(b);
				for (j = 0; j < b; j++) {
					labelling_B.m_ii(j, j + 1);
					}
				j = 0;
				}
			else if (strncmp(buf, "k=", 2) == 0) {
				p_str = &buf[2];
				s_scan_int(&p_str, &k);
				}
			else if (strncmp(buf, "END", 3) == 0) {
				write_geo(out_fp, nb_geo, "", &I, &labelling_P, &labelling_B);
				break;
				}
			else {
				if (v < 0)
					return error("v not defined\n");
				if (k < 0)
					return error("k not defined\n");
				p_str = buf;
				for (i = 0; i < k; i++) {
					s_scan_int(&p_str, &a);
					I.m_iji(a - 1, j, 1);
					}
				j++;
				}
			
			}
		}
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}

INT convert_base_blocks_geo(BYTE *base_file, BYTE *generators_file, 
	BYTE *geo_file, INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT v, b;
	INT i, k, a, l;
	BYTE *p_str, str[1024], geo_label[1000];
	VECTOR_OB labelling_P, labelling_B;
	MATRIX_OB I;
	INT geo_nr = 0;
	VECTOR_OB gen;
	VECTOR_OB BB, BB1, B, O;
	
	
	in_fp = fopen(base_file, "r");
	out_fp = fopen(geo_file, "w");
	
	read_file_of_generators(&gen, generators_file);
	v = vec_generators_degree(&gen);
	
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
		
	
		p_str = &buf[8];
		s_scan_int(&p_str, &geo_nr);
		geo_label[0] = 0;
		s_scan_token(&p_str, geo_label);
		printf("GEOMETRY %ld %s (v=%ld)\n", geo_nr, geo_label, v);
		BB.m_il(0);
		
		k = -1;
		while (TRUE) {
			if (fgets(buf, BUFSIZE, in_fp) == NIL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0; /* delete newline */
			if (buf[0] == '#')
				continue;
		
			if (strncmp(buf, "k=", 2) == 0) {
				p_str = &buf[2];
				s_scan_int(&p_str, &k);
				}
			else if (strncmp(buf, "END", 3) == 0) {
				b = BB.s_li();
				INT j, ii;
				VECTOR_OP pO;
				
				I.m_ilih_n(b, v);
				for (j = 0; j < b; j++) {
					pO = (VECTOR_OP) BB.s_i(j);
					if (pO->s_li() != k)
						return error("design_output_geo() pO->s_li() != k");
					for (ii = 0; ii < pO->s_li(); ii++) {
						i = pO->s_ii(ii);
						I.m_iji(i, j, 1);
						}
					}


				labelling_P.m_il(v);
				for (i = 0; i < v; i++) {
					labelling_P.m_ii(i, i + 1);
					}
				labelling_B.m_il(b);
				for (j = 0; j < b; j++) {
					labelling_B.m_ii(j, j + 1);
					}

				write_geo(out_fp, geo_nr, geo_label, &I, &labelling_P, &labelling_B);
				break;
				}
			else {
				if (k < 0)
					return error("k not defined\n");
				p_str = buf;
				B.m_il(k);
				for (i = 0; i < k; i++) {
					s_scan_int(&p_str, &a);
					B.m_ii(i, a - 1);
					}
				dc_calc_set_orbit(&gen, &B, &O, FALSE /*f_v*/, FALSE /*f_vv*/);
				BB.append(&O, &BB1);
				BB1.swap(&BB);
				}
			
			}
		}
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}

INT convert_inc_geo(BYTE *inc_file, BYTE *geo_file, BYTE *base_fname, 
	INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT v, b, nb_X;
	INT i, j, k, a;
	INT *X;
	BYTE *p_str, str[1024];
	VECTOR_OB labelling_P, labelling_B;
	INT f_permute_points, f_permute_blocks;
	PERMUTATION_OB p0, q0;
	MATRIX_OB I;
	INT nb_geo = 0;
	
	
	in_fp = fopen(inc_file, "r");
	out_fp = fopen(geo_file, "w");
	
	fscanf(in_fp, "%ld %ld %ld", &v, &b, &nb_X);
	printf("reading geometries of size %ld %ld with nb_X=%ld\n", 
		v, b, nb_X);


	X = (INT *) my_malloc((nb_X + 1) * sizeof(INT), "do_it X");
	if (X == NIL)
		return error("pp() no memory for X");


	while (TRUE) {
		// I.freeself();
		fscanf(in_fp, "%ld", &a);
		if (a == -1)
			break;
		X[0] = a;
		for (i = 1; i < nb_X; i++) {
			fscanf(in_fp, "%ld", &a);
			X[i] = a;
			}
		fgets(buf, BUFSIZE, in_fp);
		p_str = buf;
		f_permute_points = FALSE;
		f_permute_blocks = FALSE;
		while (s_scan_token(&p_str, str)) {
			if (strcmp(str, "permute_points") == 0) {
				f_permute_points = TRUE;
				p0.m_il(v);
				for (i = 0; i < v; i++) {
					if (!s_scan_int(&p_str, &a)) 
						return error("error while scanning permutation of points");
					if (a < 1 || a > v)
						return error("error while scanning permutation of points: a < 1 || a > v");
					p0.m_ii(i, a);
					}
				printf("initial permutation of points:\n");
				p0.println();
				}
			else if (strcmp(str, "permute_points_to") == 0) {
				f_permute_points = TRUE;
				p0.m_il(v);
				for (i = 0; i < v; i++) {
					if (!s_scan_int(&p_str, &a)) 
						return error("error while scanning permutation of points");
					if (a < 1 || a > v)
						return error("error while scanning permutation of points: a < 1 || a > v");
					p0.m_ii(i, a);
					}
				p0.invers_apply();
				printf("initial permutation of points:\n");
				p0.println();
				}
			else if (strcmp(str, "permute_blocks") == 0) {
				f_permute_blocks = TRUE;
				q0.m_il(b);
				for (i = 0; i < b; i++) {
					if (!s_scan_int(&p_str, &a)) 
						return error("error while scanning permutation of blocks");
					if (a < 1 || a > b)
						return error("error while scanning permutation of blocks: a < 1 || a > b");
					q0.m_ii(i, a);
					}
				printf("initial permutation of blocks:\n");
				q0.println();
				}
			else if (strcmp(str, "permute_blocks_to") == 0) {
				f_permute_blocks = TRUE;
				q0.m_il(b);
				for (i = 0; i < b; i++) {
					if (!s_scan_int(&p_str, &a)) 
						return error("error while scanning permutation of blocks");
					if (a < 1 || a > b)
						return error("error while scanning permutation of blocks: a < 1 || a > b");
					q0.m_ii(i, a);
					}
				q0.invers_apply();
				printf("initial permutation of blocks:\n");
				q0.println();
				}
			else {
				printf("skipping unknown command: %s\n", str);
				}
			}
		I.m_ilih_n(b, v);
		for (k = 0; k < nb_X; k++) {
			a = X[k];
			i = a / b;
			j = a % b;
			I.m_iji(i, j, 1);
			}
		
		nb_geo++;
		if (f_range) {
			if (nb_geo < range_f || nb_geo >= range_f + range_l) {
				continue;
				}
			}
		printf("\ngeo no %5ld\n", nb_geo);
		
		labelling_P.m_il(v);
		for (i = 0; i < v; i++) {
			labelling_P.m_ii(i, i + 1);
			}
		labelling_B.m_il(b);
		for (j = 0; j < b; j++) {
			labelling_B.m_ii(j, j + 1);
			}
		if (!f_permute_points) {
			p0.m_il(v);
			p0.one();
			}
		if (!f_permute_blocks) {
			q0.m_il(b);
			q0.one();
			}
		I.apply_perms(&p0, &q0);
		labelling_P.apply_perm(&p0);
		labelling_B.apply_perm(&q0);
		
		write_geo(out_fp, nb_geo, base_fname, &I, &labelling_P, &labelling_B);
		
		} // while
	
	my_free(X);
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}


INT convert_geo_inc(BYTE *geo_file, BYTE *inc_file, INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT v = -1, b = -1, nb_X = -1, v1, b1, nb_X1;
	INT i, j, k, a, l;
	INT *X;
	BYTE *p_str, str[1024];
	// VECTOR_OB labelling_P, labelling_B;
	// INT f_permute_points, f_permute_blocks;
	// PERMUTATION_OB p0, q0;
	MATRIX_OB I;
	INT nb_geo = 0;
	
	printf("convert_geo_inc\n"); fflush(stdout);
	
	in_fp = fopen(geo_file, "r");
	out_fp = fopen(inc_file, "w");
	

	X = (INT *) my_malloc((nb_X + 1) * sizeof(INT), "do_it X");
	if (X == NIL)
		return error("pp() no memory for X");

	
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
		
		printf("."); fflush(stdout);
		nb_geo++;
		p_str = &buf[9];
		s_scan_int(&p_str, &l);
		s_scan_token(&p_str, str);
		while (TRUE) {
			if (fgets(buf, BUFSIZE, in_fp) == NIL)
				break;
			l = strlen(buf);
			if (l)
				buf[l - 1] = 0; /* delete newline */
			if (buf[0] == '#')
				continue;
			
			if (strncmp(buf, "v=", 2) == 0) {
				sscanf(buf, "v=%ld b=%ld", &v1, &b1);
				if (v == -1) {
					v = v1;
					b = b1;
					// nb_X = nb_X1;
					}
				else {
					if (v1 != v)
						return error("different v\n");
					if (b1 != b)
						return error("different b\n");
					}
				}
			if (strncmp(buf, "INCIDENCE_MATRIX", 16) == 0) {
				I.m_ilih_n(b, v);
				nb_X1 = 0;
				for (i = 0; i < v; i++) {
					if (fgets(buf, BUFSIZE, in_fp) == NIL)
						return error("error reading INCIDENCE_MATRIX\n");
					for (j = 0; j < b; j++) {
						if (buf[j] == 'X') {
							I.m_iji(i, j, 1);
							nb_X1++;
							}
						}
					}
				if (nb_X == -1) {
					nb_X = nb_X1;
					fprintf(out_fp, "%ld %ld %d\n", v, b, nb_X);
					}
				else {
					if (nb_X1 != nb_X)
						return error("different nb_X\n");
					}
				}
			if (strncmp(buf, "END", 3) == 0) {
				write_inc(out_fp, &I);
				break;
				}
			}
		
		}
	fprintf(out_fp, "-1 %ld geometries\n", nb_geo);
	return OK;
}


INT write_inc(FILE *fp, MATRIX_OP I)
{
	INT v = I->s_hi();
	INT b = I->s_li();
	INT i, j, a;
	
	for (i = 0; i < v; i++) {
		for (j = 0; j < b; j++) {
			a = I->s_iji(i, j);
			if (a) {
				fprintf(fp, "%ld ", i * b + j);
				}
			}
		}
	fprintf(fp, "\n");
}

INT write_geo(FILE *fp, INT nb_geo, BYTE *label, MATRIX_OP I, VECTOR_OP labelling_P, VECTOR_OP labelling_B)
{
	INT v = I->s_hi();
	INT b = I->s_li();
	INT i, j, a;
	
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
	fprintf(fp, "END\n\n");
}


