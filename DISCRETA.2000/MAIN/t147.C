// t147.C
//
// 09.12.1999
//
// write geometries in tex
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

INT latex_it(BYTE *geo_file, BYTE *tex_file, INT f_range, INT range_f, INT range_l);
INT latex_geo(FILE *fp, INT geo_nr, BYTE *geo_label, 
	MATRIX_OP I, VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, 
	VECTOR_OP DDp, VECTOR_OP DDb, 
	VECTOR_OP aut_gen);
INT print_aut_gen_TEX(FILE *TEX_fp, VECTOR_OP aut_gen);
INT latex_header(FILE *fp);
INT latex_footer(FILE *fp);

void usage()
{
	printf("usage: t147.out [options] geo_file\n");
	printf("write geometries in tex\n");
	printf("available options: \n");
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
	BYTE *geo_file;
	BYTE tex_file[1000];
	INT f_range = FALSE, range_f, range_l;
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
	geo_file = argv[argc - 1];
	
	t0 = os_ticks();
	
	l = strlen(geo_file);
	if (l > 4 && strcmp(geo_file + l - 4, ".geo") == 0) {
		f_input_geo = TRUE;
		}
	else {
		return error("input file has unknown extension (should be .geo)");
		}


	strcpy(tex_file, geo_file);
	tex_file[l - 4] = 0;

	if (f_input_geo) {
		strcat(tex_file, ".tex");
		latex_it(geo_file, tex_file, f_range, range_f, range_l);
		}
	else  {
		return error("unknown input file format.");
		}
	printf("written file %s\n", tex_file);
	
	
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


INT latex_it(BYTE *geo_file, BYTE *tex_file, INT f_range, INT range_f, INT range_l)
{
	FILE *in_fp, *out_fp;
	INT v = -1, b = -1;
	INT i, j, k, a, a1, l, geo_nr, len;
	BYTE geo_label[1000];
	BYTE *p_str, str[1024];
	VECTOR_OB labelling_P, labelling_B;
	// INT f_permute_points, f_permute_blocks;
	// PERMUTATION_OB p0, q0;
	MATRIX_OB I;
	VECTOR_OB row_decomp, col_decomp;
	INT f_ddp, f_ddb;
	VECTOR_OB DDp, DDb;
	INT nb_geo = 0;
	LABRA_OB aut;
	VECTOR_OB aut_gen;
	
	
	in_fp = fopen(geo_file, "r");
	out_fp = fopen(tex_file, "w");
	
	latex_header(out_fp);
	
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
		aut_gen.m_il(0);
		DDp.m_il(0);
		DDb.m_il(0);
		f_ddp = FALSE;
		f_ddb = FALSE;
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
				latex_geo(out_fp, geo_nr, geo_label, &I, 
					&labelling_P, &labelling_B, 
					&row_decomp, &col_decomp, 
					f_ddp, f_ddb, &DDp, &DDb, 
					&aut_gen);
				printf("geo %ld %s written\n", geo_nr, geo_label);
				fflush(out_fp);
				break;
				}
			}
		
		}
	latex_footer(out_fp);
	
	fprintf(out_fp, "\n");
	
	fclose(in_fp);
	fclose(out_fp);
	return OK;
}

INT latex_geo(FILE *fp, INT geo_nr, BYTE *geo_label, 
	MATRIX_OP I, VECTOR_OP labelling_P, VECTOR_OP labelling_B, 
	VECTOR_OP row_decomp, VECTOR_OP col_decomp, 
	INT f_ddp, INT f_ddb, 
	VECTOR_OP DDp, VECTOR_OP DDb, 
	VECTOR_OP aut_gen)
{
	BYTE geo_label_tex[1000];
	BYTE str[1000];
	INT v = I->s_hi();
	INT b = I->s_li();
	INT i, j, k, a, l;
	BYTE c;
	MATRIX_OB DDP, DDB;
	INT block_size = 10;
	
	l = strlen(geo_label);
	geo_label_tex[0] = 0;
	for (i = 0; i < l; i++) {
		c = geo_label[i];
		if (c == '_') {
			strcat(geo_label_tex, "\\_");
			}
		else {
			str[0] = c;
			str[1] = 0;
			strcat(geo_label_tex, str);
			}
		}
	
	fprintf(fp, "\\subsubsection{geo no %ld %s}\n", geo_nr, geo_label_tex);

	fprintf(fp, "{\\mytt\n");
	fprintf(fp, "\\noindent%%\n");
	I->print_char_TD(fp, TRUE /* f_tex */, 
		row_decomp, col_decomp);
	fprintf(fp, "}%% \n");
	fprintf(fp, "labelling of points: \n$");
	labelling_P->latex(fp);
	fprintf(fp, "$\\\\\n");
	fprintf(fp, "labelling of blocks: \n$");
	labelling_B->latex(fp);
	fprintf(fp, "$\\\\\n");
	fprintf(fp, "\\par\\noindent\n");
	
	if (f_ddp) {
		fprintf(fp, "DDp: \n");
		// printf(fp, "$");
		// DDp->latex(fp);
		// fprintf(fp, "$\\\\\n");
		DDP.m_ilih_n(v, v);
		for (i = 0; i < v; i++) {
			for (j = i + 1; j < v; j++) {
				k = ij2k(i, j, v);
				a = DDp->s_ii(k);
				DDP.m_iji(i, j, a);
				}
			}
		DDP.latex_upper_tri_block_wise(fp, block_size);
		// fprintf(fp, "\\[\n");
		// fprintf(fp, "\\mbox{DDp=}\n");
		// DDP.latex_upper_tri(fp);
		// fprintf(fp, "\\]\n");
		}
	if (f_ddb) {
		fprintf(fp, "DDb: \n");
		// fprintf(fp, "$");
		// DDb->latex(fp);
		// fprintf(fp, "$\\\\\n");
		DDB.m_ilih_n(b, b);
		for (i = 0; i < b; i++) {
			for (j = i + 1; j < b; j++) {
				k = ij2k(i, j, b);
				a = DDb->s_ii(k);
				DDB.m_iji(i, j, a);
				}
			}
		DDB.latex_upper_tri_block_wise(fp, block_size);
		// fprintf(fp, "\\[\n");
		// fprintf(fp, "\\mbox{DDb=}\n");
		// DDB.latex_upper_tri(fp);
		// fprintf(fp, "\\]\n");
		}
	fprintf(fp, "{\n");
	I->incma_latex(fp, row_decomp, col_decomp, 
		TRUE /* f_labelled */, 
		labelling_P, labelling_B, 0 /* offset */);
	fprintf(fp, "}%%\n");
	fflush(fp);
	if (aut_gen->s_li()) {
		print_aut_gen_TEX(fp, aut_gen);
		}
	fprintf(fp, "\\clearpage\n");
}

INT print_aut_gen_TEX(FILE *TEX_fp, VECTOR_OP aut_gen)
{
	INT i, l;
	LABRA_OB aut;
	SYM_OB go;
	BYTE str[1000];
	
	aut.init_quick(aut_gen);
	aut.group_order(&go);
	str[0] = 0;
	go.sprint(str);
	fprintf(TEX_fp, "\\par\n$|Aut| = %s$\\\\\n", str);
	fprintf(TEX_fp, "\\[\n");
	fprintf(TEX_fp, "\\begin{array}{rl}\n");
	l = aut_gen->s_li();
	for (i = 0; i < l; i++) {
		if (i == 0) {
			fprintf(TEX_fp, "Aut \\; = & \\langle ");
			}
		else {
			fprintf(TEX_fp, " & ");
			}
		aut_gen->s_i(i)->latex(TEX_fp);
		if (i == l - 1) {
			fprintf(TEX_fp, " \\rangle");
			}
		else {
			fprintf(TEX_fp, ", ");
			}
		fprintf(TEX_fp, "\\\\[3pt]\n");
		}
	fprintf(TEX_fp, "\\end{array}\n");
	fprintf(TEX_fp, "\\]\n");
	return OK;
}

INT latex_header(FILE *fp)
{
	fprintf(fp, "\\documentclass[11pt,a4paper]{report}\n");
	fprintf(fp, "\\usepackage{rotating}\n");
	fprintf(fp, "\\newfont{\\mytt}{cmtt12 scaled 1000}\n");
	fprintf(fp, "\\textheight=640pt\n");
	fprintf(fp, "\\textwidth=440pt\n");
	fprintf(fp, "\\topmargin=0pt\n");
	fprintf(fp, "\\headsep=0pt\n");
	fprintf(fp, "\\footskip=10pt\n");
	fprintf(fp, "\\mathsurround=1pt\n");
	fprintf(fp, "\\evensidemargin=15pt\n");
	fprintf(fp, "\\oddsidemargin=15pt\n");
	fprintf(fp, "\n");
	fprintf(fp, "\\pagestyle{empty}\n");
	fprintf(fp, "\n");
	fprintf(fp, "\\begin{document}\n");
	fprintf(fp, "\n");
}

INT latex_footer(FILE *fp)
{

	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "\\end{document}\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
}

