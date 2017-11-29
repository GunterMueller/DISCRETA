
#include <stdlib.h>

#define BUFSIZE 50000


DATA_3D *open_data_3D(int nb_V)
{
	DATA_3D *p;
	int i;
	
	p = (DATA_3D *)my_malloc(sizeof(DATA_3D), "DATA_3D");
	if (p == NULL)
		error("open_data_3D() no memory for DATA_3D");
	p->nb_V = nb_V;
	p->x = (int *) my_malloc(nb_V * sizeof(int), "DATA_3D x");
	p->y = (int *) my_malloc(nb_V * sizeof(int), "DATA_3D y");
	p->z = (int *) my_malloc(nb_V * sizeof(int), "DATA_3D z");
	p->labels = (char **) my_malloc(nb_V * sizeof(char *), "DATA_3D labels");
	if (p->x == NULL)
		error("open_data_3D() no memory for DATA_3D, x");
	if (p->y == NULL)
		error("open_data_3D() no memory for DATA_3D, y");
	if (p->z == NULL)
		error("open_data_3D() no memory for DATA_3D, z");
	for (i = 0; i < nb_V; i++) {
		p->labels[i] = NULL;
		}
	p->f_labels = FALSE;
	return p;
}

void free_data_3D(DATA_3D *p)
{
	int i, l;
	
	l = p->nb_V;
	for (i = 0; i < l; i++) {
		if (p->labels[i]) {
			my_free(p->labels[i]);
			p->labels[i] = NULL;
			}
		}
	my_free(p->labels);
	my_free(p->x);
	my_free(p->y);
	my_free(p->z);
	my_free(p);
}

GRAPH_MULTI *open_graph_multi(int nb_V, int nb_E, int f_multi)
{
	GRAPH_MULTI *p;
	
	p = (GRAPH_MULTI *) my_malloc(sizeof(GRAPH_MULTI), "GRAPH_MULTI");
	if (p == NULL)
		error("open_graph_multi() no memory for GRAPH_MULTI");
	p->nb_V = nb_V;
	p->nb_E = nb_E;
	p->f_multi = f_multi;
	p->e1 = (int *) my_malloc(nb_E * sizeof(int), "GRAPH_MULTI e1");
	p->e2 = (int *) my_malloc(nb_E * sizeof(int), "GRAPH_MULTI e2");
	p->mult = (int *) my_malloc(nb_E * sizeof(int), "GRAPH_MULTI m");
	if (p->e1 == NULL)
		error("open_graph_multi() no memory for GRAPH_MULTI, e1");
	if (p->e2 == NULL)
		error("open_graph_multi() no memory for GRAPH_MULTI, e2");
	if (p->mult == NULL)
		error("open_graph_multi() no memory for GRAPH_MULTI, mult");
	return p;
}

void free_graph_multi(GRAPH_MULTI *p)
{
	my_free(p->e1);
	my_free(p->e2);
	my_free(p->mult);
	my_free(p);
}

int save_graph_multi_to_file(char *fname, GRAPH_MULTI *gm, DATA_3D *d3)
{
	char name[1000], str[256], *p;
	FILE *fp, *fp1;
	int i, l;
	int nb_V = gm->nb_V;
	int nb_E = gm->nb_E;
	
	printf("SOLID::write_graphfile() fname=%s nb_V=%d nb_E=%d\n", 
		fname, nb_V, nb_E); fflush(stdout);
	strcpy(name, fname);
	if ((p = strrchr(name, '.')) != NULL) {
		*p = 0;
		}

	system("rm a");
	system("date >a");
	fp1 = fopen("a", "r");
	fgets(str, 256, fp1);
	fclose(fp1);
	if ((l = strlen(str)) > 0) {
		str[l - 1] = 0;
		}


	fp = fopen(fname, "w");
	fprintf(fp, "<GRAPH NAME=\"%s\" NUMBER_OF_VERTICES=%d NUMBER_OF_EDGES=%d>\n", 
		name, nb_V, nb_E);
	fprintf(fp, "<!-- created by DISCRETA, %s -->\n", str);
	fprintf(fp, "<COORDS3D_INT>\n");
	for (i = 0; i < nb_V; i++) 
	{
		fprintf(fp, "%d %d %d\n", 
			(int) d3->x[i], 
			(int) d3->y[i], 
			(int) d3->z[i]); 
	}
	fprintf(fp, "</COORDS3D_INT>\n");
	fprintf(fp, "<EDGELIST>\n");
	for (i = 0; i < nb_E; i++) 
	{
		fprintf(fp, "%d %d\n", gm->e1[i], gm->e2[i]);
	}
	fprintf(fp, "</EDGELIST>\n");
	if (d3->f_labels) {
		fprintf(fp, "<VERTEXLABELS>\n");
		for (i = 0; i < nb_V; i++) 
		{
			fprintf(fp, "\"%s\"\n", d3->labels[i]); 
		}
		fprintf(fp, "</VERTEXLABELS>\n");
	}

	fprintf(fp, "</GRAPH>\n\n");
	fclose(fp);
	return TRUE;
}


#if 0
int save_graph_multi_to_file(char *fname, GRAPH_MULTI *gm, DATA_3D *d3)
{
	FILE *fp;
	int i;
	
	fp = fopen(fname, "w");
	fprintf(fp, "NUMBER_OF_VERTICES %d\n", gm->nb_V);
	fprintf(fp, "NUMBER_OF_EDGES %d\n", gm->nb_E);
	if (gm->f_multi)
		fprintf(fp, "MULTIGRAPH\n");
	fprintf(fp, "VERTICES\n");
	for (i = 0; i < gm->nb_V; i++) {
		fprintf(fp, "%d %d %d\n", d3->x[i], d3->y[i], d3->z[i]);
		}
	fprintf(fp, "EDGES\n");
	for (i = 0; i < gm->nb_E; i++) {
		fprintf(fp, "%d %d", gm->e1[i], gm->e2[i]);
		if (gm->f_multi) {
			fprintf(fp, " %d\n", gm->mult[i]);
			}
		else {
			fprintf(fp, "\n");
			}
		}
	if (d3->f_labels) {
		fprintf(fp, "VERTEX-LABELS\n");
		for (i = 0; i < gm->nb_V; i++) {
			fprintf(fp, "%s\n", d3->labels[i]);
			}
		}

	fclose(fp);
	return TRUE;
}
#endif


static int read_graph2(FILE *fp, char *buf, char *pbuf, char *name, GRAPH_MULTI **gm, DATA_3D **d3);
static int read_until_closing_command(FILE *fp, char *command);

int read_graph_multi_from_file(char *fname, GRAPH_MULTI **gm, DATA_3D **d3)
{
	char name[BUFSIZE], buf[BUFSIZE], *pbuf, token[BUFSIZE];
	FILE *fp;
	int l;
	
	printf("read_graph_multi_from_file() opening file %s\n", fname); fflush(stdout);
	fp = fopen(fname, "r");
	while (fgets(buf, BUFSIZE, fp) != NULL) {
		l = strlen(buf);
		if (l > 0 && buf[l - 1] == '\n') {
			buf[l - 1] = 0;
			l--;
			}
		pbuf = buf;
		while (TRUE) {
			if (!s_scan_token(&pbuf, token)) // end of line
				break;
			if (strcmp(token, "<") == 0) {
				if (!s_scan_token(&pbuf, token)) // end of line
					error("end of line after <");
				if (strcmp(token, "GRAPH") == 0) {
					if (!read_graph2(fp, buf, pbuf, name, gm, d3))
						error("error while reading GRAPH");
					else {
						fclose(fp);
						return TRUE;
						}
					}
				}
			}
		}
	fclose(fp);
	error("error while reading GRAPH");
	return FALSE;
#if 0
	fscanf(fp, "%d", &nb_V);
	fscanf(fp, "%d", &nb_E);
	fscanf(fp, "%d", &f_multi);
	GM = open_graph_multi(nb_V, nb_E, f_multi);
	D3 = open_data_3D(nb_V);
	for (i = 0; i < nb_V; i++) {
		fscanf(fp, "%d", &x);
		fscanf(fp, "%d", &y);
		fscanf(fp, "%d", &z);
		D3->x[i] = x;
		D3->y[i] = y;
		D3->z[i] = z;
		}
	for (i = 0; i < nb_E; i++) {
		fscanf(fp, "%d", &e1);
		fscanf(fp, "%d", &e2);
		if (f_multi)
			fscanf(fp, "%d", &m);
		else
			m = 1;
		GM->e1[i] = e1;
		GM->e2[i] = e2;
		GM->mult[i] = m;
		}
#endif
}

#undef VERBOSE

static int read_graph2(FILE *fp, char *buf, char *pbuf, char *name, GRAPH_MULTI **gm, DATA_3D **d3)
{
	char token[BUFSIZE];
	GRAPH_MULTI *GM;
	DATA_3D *D3;
	int f_multi, nb_V, nb_E;
	int i, x, y, z, e1, e2, m, l;
	INT ii;

	GM = NULL;
	D3 = NULL;
	nb_V = -1;
	nb_E = -1;
	f_multi = FALSE;
	name[0] = 0;
	while (TRUE) {
		if (!s_scan_token(&pbuf, token)) // end of line
			error("end of line encountered");
		if (strcmp(token, "NUMBER_OF_VERTICES") == 0) {
			if (!s_scan_token(&pbuf, token)) // end of line
				error("end of line encountered, expecting '='");
			s_scan_int(&pbuf, &ii);
			nb_V = ii;
#ifdef VERBOSE
			printf("nb_V = %d\n", nb_V);
#endif
			}
		else if (strcmp(token, "NUMBER_OF_EDGES") == 0) {
			if (!s_scan_token(&pbuf, token)) // end of line
				error("end of line encountered, expecting '='");
			s_scan_int(&pbuf, &ii);
			nb_E = ii;
#ifdef VERBOSE
			printf("nb_E = %d\n", nb_E);
#endif
			}
		else if (strcmp(token, "NAME") == 0) {
			if (!s_scan_token(&pbuf, token)) // end of line
				error("end of line encountered, expecting '='");
			s_scan_str(&pbuf, name);
#ifdef VERBOSE
			printf("name = %s\n", name);
#endif
			}
		else if (strcmp(token, ">") == 0) {
			break;
			}
		} // while
	if (nb_V < 0)
		error("NUMBER_OF_VERTICES not specified"); 
	if (nb_E < 0)
		error("NUMBER_OF_VERTICES not specified");
	GM = open_graph_multi(nb_V, nb_E, f_multi);
	D3 = open_data_3D(nb_V);
#ifdef VERBOSE
	printf("reading graph, nb_V = %d, nb_E=%d\n", nb_V, nb_E);
#endif
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			error("end of line encountered");
		l = strlen(buf);
		if (l > 0 && buf[l - 1] == '\n') {
			buf[l - 1] = 0;
			l--;
			}
		pbuf = buf;
		if (!s_scan_token(&pbuf, token)) // end of line
			continue;
		if (strcmp(token, "<") == 0) {
			if (!s_scan_token(&pbuf, token)) // end of line
				error("end of line after <");
			if (strcmp(token, "/") == 0) {
				if (!s_scan_token(&pbuf, token)) // end of line
					error("end of line after </");
				if (strcmp(token, "GRAPH") == 0)
					break;
				}
			else if (strcmp(token, "COORDS3D_INT") == 0) {
#ifdef VERBOSE
				printf("reading coordinates\n");
#endif
				for (i = 0; i < nb_V; i++) {
					fscanf(fp, "%d", &x);
					fscanf(fp, "%d", &y);
					fscanf(fp, "%d", &z);
#ifdef VERBOSE
					printf("%d: %d %d %d\n", i, x, y, z);
#endif
					D3->x[i] = x;
					D3->y[i] = y;
					D3->z[i] = z;
					}
				read_until_closing_command(fp, "COORDS3D_INT");
				} // end else if COORDS3D_INT
			else if (strcmp(token, "EDGELIST") == 0) {
				for (i = 0; i < nb_E; i++) {
					fscanf(fp, "%d", &e1);
					fscanf(fp, "%d", &e2);
					if (f_multi)
						fscanf(fp, "%d", &m);
					else
						m = 1;
					GM->e1[i] = e1;
					GM->e2[i] = e2;
					GM->mult[i] = m;
					}
				read_until_closing_command(fp, "EDGELIST");
				} // end else if EDGELIST
			else if (strcmp(token, "VERTEXLABELS") == 0) {
				for (i = 0; i < nb_V; i++) {
					char buf1[10000], *pbuf1, str[1000];
					int l;
				
					fgets(buf1, 10000, fp);
					l = strlen(buf1);
					if (l && buf1[l - 1] == '\n') {
						buf1[l - 1] = 0;
						l--;
						}
					pbuf1 = buf1;
					s_scan_str(&pbuf1, str);
					l = strlen(str);
					D3->labels[i] = (char *)malloc((l + 1) * sizeof(char));
					strcpy(D3->labels[i], str);
#ifdef VERBOSE
					printf("label %d = %s\n", i, D3->labels[i]);
#endif
					}
				D3->f_labels = TRUE;
				read_until_closing_command(fp, "VERTEXLABELS");
				} // end else if EDGELIST
			}
		} // while
	// finished reading <GRAPH> ... </GRAPH>
	*gm = GM;
	*d3 = D3;
	printf("read graph %s with nb_V = %d nb_E = %d\n", name, GM->nb_V, GM->nb_E);
	fflush(stdout);
	return TRUE;
}

static int read_until_closing_command(FILE *fp, char *command)
{
	char buf[BUFSIZE], *pbuf, token[BUFSIZE];
	int l;
 
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL) {
			printf("expecting /%s\n", command);
			return error("read_until_closing_command(): end of line encountered");
			}
		l = strlen(buf);
		if (l > 0 && buf[l - 1] == '\n') {
			buf[l - 1] = 0;
			l--;
			}
		pbuf = buf;
		if (!s_scan_token(&pbuf, token)) // end of line
			continue;
		if (strcmp(token, "<") == 0) {
			if (!s_scan_token(&pbuf, token)) { // end of line
				printf("expecting /%s\n", command);
				return error("end of line after <");
				}
			if (strcmp(token, "/") == 0) {
				if (!s_scan_token(&pbuf, token)) { // end of line
					printf("expecting /%s\n", command);
					return error("end of line after </");
					}
				if (strcmp(token, command) == 0)
					return TRUE;
				}
			}
		}
	return FALSE;
}

#if 0
int read_graph_multi_from_file(char *fname, GRAPH_MULTI **gm, DATA_3D **d3)
{
	char buf[BUFSIZE], *pbuf, token[BUFSIZE];
	FILE *fp;
	GRAPH_MULTI *GM = NULL;
	DATA_3D *D3 = NULL;
	int f_multi = FALSE, nb_V = -1, nb_E = -1;
	int i, x, y, z, e1, e2, m, l;
	
	printf("read_graph_multi_from_file() opening file %s\n", fname); fflush(stdout);
	fp = fopen(fname, "r");
	while (fgets(buf, BUFSIZE, fp) != NULL) {
		l = strlen(buf);
		if (l > 0 && buf[l - 1] == '\n') {
			buf[l - 1] = 0;
			l--;
			}
		pbuf = buf;
		while (TRUE) {
			if (!my_s_scan_token(&pbuf, token)) // end of line
				break;
			if (strcmp(token, "#") == 0) { // comment line
				break;
				}
			if (strcmp(token, "NUMBER_OF_VERTICES") == 0) {
				nb_V = my_s_scan_int(&pbuf);
				printf("nb_V = %d\n", nb_V);
				}
			else if (strcmp(token, "NUMBER_OF_EDGES") == 0) {
				nb_E = my_s_scan_int(&pbuf);
				printf("nb_E = %d\n", nb_E);
				}
			else if (strcmp(token, "MULTIGRAPH") == 0) {
				f_multi = TRUE;
				}
			else if (strcmp(token, "VERTICES") == 0) {
				GM = open_graph_multi(nb_V, nb_E, f_multi);
				D3 = open_data_3D(nb_V);
				for (i = 0; i < nb_V; i++) {
					fscanf(fp, "%d", &x);
					fscanf(fp, "%d", &y);
					fscanf(fp, "%d", &z);
					D3->x[i] = x;
					D3->y[i] = y;
					D3->z[i] = z;
					}
				}
			else if (strcmp(token, "VERTEX-LABELS") == 0) {
				for (i = 0; i < nb_V; i++) {
					char buf1[10000];
					int l;
					
					fgets(buf1, 10000, fp);
					l = strlen(buf1);
					if (l && buf1[l - 1] == '\n') {
						buf1[l - 1] = 0;
						l--;
						}
					D3->labels[i] = (char *)my_malloc((l + 1) * sizeof(char), 
						"read_graph_multi_from_file label");
					strcpy(D3->labels[i], buf1);
					printf("label %ld = %s\n", i, D3->labels[i]);
					}
				D3->f_labels = TRUE;
				}
			else if (strcmp(token, "EDGES") == 0) {
				if (GM == NULL) {
					return error("please specify VERTICES before EDGES\n");
					}
				for (i = 0; i < nb_E; i++) {
					fscanf(fp, "%d", &e1);
					fscanf(fp, "%d", &e2);
					if (f_multi)
						fscanf(fp, "%d", &m);
					else
						m = 1;
					GM->e1[i] = e1;
					GM->e2[i] = e2;
					GM->mult[i] = m;
					}
				}
			}
		}
#if 0
	fscanf(fp, "%d", &nb_V);
	fscanf(fp, "%d", &nb_E);
	fscanf(fp, "%d", &f_multi);
	GM = open_graph_multi(nb_V, nb_E, f_multi);
	D3 = open_data_3D(nb_V);
	for (i = 0; i < nb_V; i++) {
		fscanf(fp, "%d", &x);
		fscanf(fp, "%d", &y);
		fscanf(fp, "%d", &z);
		D3->x[i] = x;
		D3->y[i] = y;
		D3->z[i] = z;
		}
	for (i = 0; i < nb_E; i++) {
		fscanf(fp, "%d", &e1);
		fscanf(fp, "%d", &e2);
		if (f_multi)
			fscanf(fp, "%d", &m);
		else
			m = 1;
		GM->e1[i] = e1;
		GM->e2[i] = e2;
		GM->mult[i] = m;
		}
#endif
	*gm = GM;
	*d3 = D3;
	fclose(fp);
	printf("read_graph_multi_from_file() read file %s "
		"with nb_V = %d nb_E = %d\n", fname, GM->nb_V, GM->nb_E);
	fflush(stdout);
	return TRUE;
}
#endif

void graph_standard_coordinates(GRAPH_MULTI *gm, DATA_3D *d3)
{
	DATA_3D_BOUNDS b;
	int x0, y0, z0;
	double sx, sy, sz;
	
	graph_center(gm, d3, &x0, &y0, &z0);
	data_3d_translate(d3, -x0, -y0, -z0);
	data_3d_compute_bounds(d3, &b);
	sx = (double) STANDARD_SCALE_X / (double) (b.xsize >> 1);
	sy = (double) STANDARD_SCALE_Y / (double) (b.ysize >> 1);
	sz = (double) STANDARD_SCALE_Z / (double) (b.zsize >> 1);
	data_3d_scale(d3, sx, sy, sz);
}

void graph_center(GRAPH_MULTI *gm, DATA_3D *d3, int *x0, int *y0, int *z0)
{
	int x = 0, y = 0, z = 0, i;
	double s;
	
	for (i = 0; i < gm->nb_V; i++) {
		x += d3->x[i];
		y += d3->y[i];
		z += d3->z[i];
		}
	s = 1. / (double)gm->nb_V;
	*x0 = (int)((double)x * s);
	*y0 = (int)((double)y * s);
	*z0 = (int)((double)z * s);
}

void data_3d_compute_bounds(DATA_3D *p, DATA_3D_BOUNDS *b)
{
	int xmin, xmax;
	int ymin, ymax;
	int zmin, zmax;
	int i, x, y, z;
	
	if (p->nb_V <= 0)
		error("data_3d_compute_bounds() nb_V <= 0");
	xmin = xmax = p->x[0];
	ymin = ymax = p->y[0];
	zmin = zmax = p->z[0];
	for (i = 1; i < p->nb_V; i++) {
		x = p->x[i];
		y = p->y[i];
		z = p->z[i];
		xmin = MINIMUM(xmin, x);
		ymin = MINIMUM(ymin, y);
		zmin = MINIMUM(zmin, z);
		xmax = MAXIMUM(xmax, x);
		ymax = MAXIMUM(ymax, y);
		zmax = MAXIMUM(zmax, z);
		}
	b->xmin = xmin;
	b->ymin = ymin;
	b->zmin = zmin;
	b->xsize = xmax - xmin;
	b->ysize = ymax - ymin;
	b->zsize = zmax - zmin;
}

void bounds_3D_print(DATA_3D_BOUNDS *b)
{
	printf("xmin=%d xsize=%d\n", b->xmin, b->xsize);
	printf("ymin=%d ysize=%d\n", b->ymin, b->ysize);
	printf("zmin=%d zsize=%d\n", b->zmin, b->zsize);
}

void bounds_3D_center(DATA_3D_BOUNDS *b, int *x0, int *y0, int *z0)
{
	*x0 = b->xmin + (b->xsize >> 1);
	*y0 = b->ymin + (b->ysize >> 1);
	*z0 = b->zmin + (b->zsize >> 1);
}

void data_3d_translate(DATA_3D *p, int dx, int dy, int dz)
{
	int i;
	
	for (i = 0; i < p->nb_V; i++) {
		p->x[i] += dx;
		p->y[i] += dy;
		p->z[i] += dz;
		}
}

void data_3d_scale(DATA_3D *p, double sx, double sy, double sz)
{
	int i;
	
	for (i = 0; i < p->nb_V; i++) {
		p->x[i] = (int)((double)p->x[i] * sx);
		p->y[i] = (int)((double)p->y[i] * sy);
		p->z[i] = (int)((double)p->z[i] * sz);
		}
}

