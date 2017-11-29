#include <stdio.h>

#define BUFSIZE 10000
#define FALSE 0
#define TRUE 1
#define WIDTH 80

int f_filename = TRUE;
int f_underscore_translation = TRUE;

int do_it(FILE *fin, FILE *fout, char *fname_in, char *search_string);
void output_verb_string(FILE *fout, char *s);
int convert(FILE *fin, FILE *fout, char *fname_in);
void convert_to_tex(char *p, char *q);

int main(int argc, char **argv)
{
	int i, j, l;
	char fname_in[1024];
	char fname[1024];
	FILE *fin, *fout;
	int first_file = 1;
	char *p;
	int f_out_prefix = FALSE;
	char *out_prefix = NULL;
	int f_string = FALSE;
	char *search_string = "TEXDOCU";
	int f_convert = FALSE;
		
	for (i = 1; i < argc; i++) {
		p = argv[i];
		if (p[0] != '-') {
			first_file = i;
			break;
			}
		if (strcmp(p, "-no_filename") == 0)
			f_filename = FALSE;
		else if (strcmp(p, "-convert") == 0)
			f_convert = TRUE;
		else if (strcmp(p, "-no_underscore_translation") == 0)
			f_underscore_translation = FALSE;
		else if (strcmp(p, "-outprefix") == 0) {
			f_out_prefix = TRUE;
			out_prefix = argv[i + 1];
			i++;
			}
		else if (strcmp(p, "-string") == 0) {
			f_string = TRUE;
			search_string = argv[i + 1];
			i++;
			}
		}
	for (i = first_file; i < argc; i++) {
		strcpy(fname_in, argv[i]);
		if (f_out_prefix) {
			strcpy(fname, out_prefix);
			strcat(fname, argv[i]);
			}
		else {
			strcpy(fname, argv[i]);
			}
		printf("in-file %s\n", fname_in);
		fflush(stdout);
		l = strlen(fname);
#if 0
		for (j = l - 1; j >= 0; j--) {
			if (fname[j] == '.') {
				strcpy(fname + j + 1, "tex");
				break;
				}
			}
		if (j < 0) {
			strcpy(fname + l, ".tex");
			}
#endif
		strcpy(fname + l, ".tex");
		printf("out-file %s\n", fname);
		fflush(stdout);
		fin = fopen(fname_in, "r");
		if (fin == NULL) {
			printf("can't open %s for reading, skipping\n", fname_in);
			continue;
			}
		fout = fopen(fname, "w");
		if (fout == NULL) {
			printf("can't open %s for writing, skipping\n", fname);
			fclose(fin);
			continue;
			}
		if (f_convert) {
			convert(fin, fout, fname_in);
			}
		else {
			do_it(fin, fout, fname_in, search_string);
			}

		fclose(fin);
		fclose(fout);
		}
	
	return 0;
}

int do_it(FILE *fin, FILE *fout, char *fname_in, char *search_string)
{
	char buf[BUFSIZE];
	int l, i, ii, j;
	int line, f_else_seen, f_code, f_preformatted;
	int save_f_no_underscore_translation;
	char search_text[1000];
	int search_text_len;
	char fname2[1024];
	char c;
	
	sprintf(search_text, "#if %s", search_string);
	search_text_len = strlen(search_text);
	if (f_filename) {
		for (i = 0, ii = 0; i < strlen(fname_in); i++) {
			c = fname_in[i];
			if (c == '_') {
				fname2[ii++] = '\\';
				}
			fname2[ii++] = c;
			}
		fname2[ii] = 0;
		fprintf(fout, "\\section{{\\tt %s}}\n", fname2);
		/* fprintf(fout, "\\subsection{file \\verb'%s':}\\\\\n", fname2);*/
		}
	while (1) {
		if (fgets(buf, BUFSIZE, fin) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strncmp(buf, search_text, search_text_len) == 0) {
			fprintf(fout, "\\begin{flushleft}\n");
			line = 0;
			f_else_seen = FALSE;
			f_code = FALSE;
			f_preformatted = FALSE;
			while (1) {
				if (fgets(buf, BUFSIZE, fin) == NULL)
					break;
				l = strlen(buf);
				if (l)
					buf[l - 1] = 0;
				if (strncmp(buf, "#else", 5) == 0) {
					f_else_seen = TRUE;
					fprintf(fout, "\\end{flushleft}\n");
					fprintf(fout, "\\begin{quote}\n");
					continue;
					}
				if (strncmp(buf, "//CODE", 6) == 0) {
					fprintf(fout, "\\par\n");
					f_code = TRUE;
					continue;
					}
				if (strncmp(buf, "///CODE", 7) == 0) {
					f_code = FALSE;
					continue;
					}
				if (strncmp(buf, "//PRE", 5) == 0) {
					f_preformatted = TRUE;
					continue;
					}
				if (strncmp(buf, "///PRE", 6) == 0) {
					f_preformatted = FALSE;
					continue;
					}
				if (strncmp(buf, "//TEX", 5) == 0) {
					f_underscore_translation = FALSE;
					continue;
					}
				if (strncmp(buf, "///TEX", 6) == 0) {
					f_underscore_translation = TRUE;
					continue;
					}
				if (strncmp(buf, "#endif", 6) == 0) {
					if (f_else_seen) {
						fprintf(fout, "\\end{quote}\n");
						fprintf(fout, "\\par\n");
						}
					else {
						fprintf(fout, "\\end{flushleft}\n");
						}
					break;
					}
#if 0
				if (f_underscore_translation) {
					l = strlen(buf);
					for (i = l - 1; i >= 0; i--) {
						if (buf[i] == '_') {
							for (j = l - 1; j >= i; j--)
								buf[j + 1] = buf[j];
							l++;
							buf[i] = '\\';
							}
						}
					}
#endif
				if (!f_else_seen)
					output_verb_string(fout, buf);
				else {
					if (f_code) {
						/* output_verb_string(fout, buf);*/
						fprintf(fout, "\\verb'%s'\\\\\n", buf);
						}
					else if (f_preformatted) {
						fprintf(fout, "%s\\\\\n", buf);
						}
					else {
						fprintf(fout, "%s\n", buf);
						}
					}
				line++;
				}
			}
		}
	fprintf(fout, "\n");
	return 1;
}

void output_verb_string(FILE *fout, char *s)
{
	int i, l, ll, first;
	char str[BUFSIZE];
	char str_tex[BUFSIZE];
	char c;
	
	l = strlen(s);
	ll = l / WIDTH;
	if (ll * WIDTH < l)
		ll++;
	for (i = 0; i < ll; i++) {
		first = i * WIDTH;
		strncpy(str, s + first, WIDTH);
		str[WIDTH] = 0;
		convert_to_tex(str, str_tex);
		fprintf(fout, "{\\tt %s}\\\\\n", str_tex);
		
#if 0
		if (first + WIDTH > l)
			fprintf(fout, "\\noindent\\verb'%s'\\\\\n", s + first);
		else {
			c = s[first + WIDTH];
			s[first + WIDTH] = 0;
			fprintf(fout, "\\noindent\\verb'%s'\\\\\n", s + first);
			s[first + WIDTH] = c;
			}
#endif
		}
}

int convert(FILE *fin, FILE *fout, char *fname_in)
{
	char buf[BUFSIZE];
	char buf2[BUFSIZE];
	int l, i, ii, j;
	int line;
	char fname2[1024];
	char c;
	
	if (f_filename) {
		for (i = 0, ii = 0; i < strlen(fname_in); i++) {
			c = fname_in[i];
			if (c == '_') {
				fname2[ii++] = '\\';
				}
			fname2[ii++] = c;
			}
		fname2[ii] = 0;
		fprintf(fout, "\\par file {\\tt %s}:\n\\par\n", fname2);
		/* fprintf(fout, "\\subsection{file \\verb'%s':}\\\\\n", fname2);*/
		}
	fprintf(fout, "\\begin{flushleft}\n");
	fprintf(fout, "{\\tt\n");
	line = 0;
	while (1) {
		if (fgets(buf, BUFSIZE, fin) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		convert_to_tex(buf, buf2);
#if 0
		ii = 0;
		for (i = 0, ii = 0; i < strlen(buf); i++) {
			c = buf[i];
			if (c == '_') {
				buf2[ii++] = '\\';
				}
			buf2[ii++] = c;
			}
		buf2[ii] = 0;
#endif
		/* fprintf(fout, "\\par\n"); */
		/* output_verb_string(fout, buf2);*/
		if (strlen(buf2) == 0)
			strcpy(buf2, "\\ ");
		fprintf(fout, "%s\\\\\n", buf2);
		line++;
		}
	fprintf(fout, "}% tt\n");
	fprintf(fout, "\\end{flushleft}\n");
	fprintf(fout, "\n");
	return 1;
}

void convert_to_tex(char *p, char *q)
{
	int i, j, l;
	char c;
	
	i = 0;
	q[0] = 0;
	while ((c = p[i++]) != 0) {
		if (c == '\t') {
			strcat(q, "\\ \\ ");
			}
		else if (c == ' ') {
			strcat(q, "\\ ");
			}
		else if (c == '\\') {
			strcat(q, "\\symbol{92}");
			}
		else if (c == '\'') { /* Akut */
			strcat(q, "\\symbol{19}");
			}
		else if (c == ',') {
			strcat(q, "\\symbol{44}");
			}
		else if (c == '!') {
			strcat(q, "\\symbol{33}");
			}
		else if (c == '\"') {
			strcat(q, "\\symbol{34}");
			}
		else if (c == '.') {
			strcat(q, "\\symbol{46}");
			}
		else if (c == '-') {
			strcat(q, "\\symbol{45}");
			}
		else if (c == '#') {
			strcat(q, "\\symbol{35}");
			}
		else if (c == '$') {
			strcat(q, "\\symbol{36}");
			}
		else if (c == '&') {
			strcat(q, "\\symbol{38}");
			}
		else if (c == '~') {
			strcat(q, "\\symbol{126}");
			}
		else if (c == '_') {
			strcat(q, "\\symbol{95}");
			}
		else if (c == '^') {
			strcat(q, "\\symbol{94}");
			}
		else if (c == '%') {
			strcat(q, "\\symbol{37}");
			}
		else if (c == '{') {
			strcat(q, "\\symbol{123}");
			}
		else if (c == '}') {
			strcat(q, "\\symbol{125}");
			}
		else if (c == '@') { /* 64 */
			strcat(q, "\\symbol{64}");
			}
		else {
			l = strlen(q);
			q[l] = c;
			q[l + 1] = 0;
			}
		}
	
}
