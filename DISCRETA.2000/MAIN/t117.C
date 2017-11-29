// t117.C
//
//
// draws a graph from a .graph file
// 
//
// Anton Betten
// Bayreuth, 17.03.1999

#include <DISCRETA/discreta.h>
#include <DISCRETA/graphics.h>

#include <stdio.h>		// printf
#include <string.h>		// memset
#include <time.h>		// clock

#include <stdlib.h>

typedef struct draw_local DRAW_LOCAL;

// t_graph.C:

#define STANDARD_SCALE_X 1000
#define STANDARD_SCALE_Y 1000
#define STANDARD_SCALE_Z 1000

typedef struct data_3D DATA_3D;
typedef struct data_3D_bounds DATA_3D_BOUNDS;
typedef struct graph_multi GRAPH_MULTI;

DATA_3D *open_data_3D(int nb_V);
void free_data_3D(DATA_3D *p);
GRAPH_MULTI *open_graph_multi(int nb_V, int nb_E, int f_multi);
void free_graph_multi(GRAPH_MULTI *p);
int save_graph_multi_to_file(char *fname, GRAPH_MULTI *gm, DATA_3D *d3);
int read_graph_multi_from_file(char *fname, GRAPH_MULTI **gm, DATA_3D **d3);
void graph_standard_coordinates(GRAPH_MULTI *gm, DATA_3D *d3);
void graph_center(GRAPH_MULTI *gm, DATA_3D *d3, int *x0, int *y0, int *z0);
void data_3d_compute_bounds(DATA_3D *p, DATA_3D_BOUNDS *b);
void bounds_3D_print(DATA_3D_BOUNDS *b);
void bounds_3D_center(DATA_3D_BOUNDS *b, int *x0, int *y0, int *z0);
void data_3d_translate(DATA_3D *p, int dx, int dy, int dz);
void data_3d_scale(DATA_3D *p, double sx, double sy, double sz);

struct data_3D {
	int nb_V;
        int *x, *y, *z;
	int f_labels;
	char **labels;
};

struct data_3D_bounds {
	int xmin, ymin, zmin;
	int xsize, ysize, zsize;
};

struct graph_multi {
	int nb_V;
	int nb_E;
	int f_multi;
	int *e1, *e2, *mult;
};


struct draw_local {
	double extrema[6];
	GRAPH_MULTI *GM;
	DATA_3D *D3;
	INT f_dots;
	double scale_x, scale_y;
	INT f_select, select_n;
	INT *selection;
};

typedef struct labelstruct {
	INT point_index;
	BYTE *align;
	BYTE *label;
} LABEL;

#define BUFSIZE 50000

#define SIZE_X 300
#define SIZE_Y 300

#define RADIUS (12)
#define RADIUS1 (180)

#define SCALE_FACTOR_X 1.
#define SCALE_FACTOR_Y 1.

#define MAX_POINTS 500

#define DASH 30
#define EPSILON 0.0001

#define SCALE 1.3
#define ROTATE 10

#define FNAME_MASK "%s.%s"


INT draw_pic(void *draw_data, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev);
INT is_selected(DRAW_LOCAL *p, INT a);

INT do_it(BYTE *base_fname, BYTE *output_fname, INT f_dots, INT f_append, 
	INT f_coordinates, INT xmin, INT ymin, INT xmax, INT ymax,
	double scale_x, double scale_y, 
	INT f_select, INT select_n, INT *selection)
{
	BYTE fname1[10000];
	BYTE fname2[10000];
	VDEVICE vdev;
	DRAW_LOCAL *data = NIL;

	data = (DRAW_LOCAL *) my_malloc(sizeof(DRAW_LOCAL), "do_it");
	data->extrema[0] =  -1000;
	data->extrema[1] = 1000;
	data->extrema[2] = -1000;
	data->extrema[3] = 1000;
	data->extrema[4] = -1800;
	data->extrema[5] = 1800;
	data->f_dots = f_dots;
	data->scale_x = scale_x;
	data->scale_y = scale_y;
	data->f_select = f_select;
	data->select_n = select_n;
	data->selection = selection;
	
	printf("scale_x = %lf scale_y = %lf\n", scale_x, scale_y);
	printf("base_fname=%s output_fname=%s\n", base_fname, output_fname);
	sprintf(fname1, "%s.graph", base_fname);
	if (!read_graph_multi_from_file(fname1, &data->GM, &data->D3))
		return FALSE;
	graph_standard_coordinates(data->GM, data->D3);

	sprintf(fname2, FNAME_MASK, output_fname, "mp");
	if (f_append) {
		FILE *fp;
		INT xmin_, ymin_, xmax_, ymax_;
		
		printf("append drawing to file %s:\n", fname2); fflush(stdout);
		fp = fopen(fname2, "a");
		mp_begin_drawing(fp);
		if (f_coordinates) {
			xmin_ = xmin;
			ymin_ = ymin;
			xmax_ = xmax;
			ymax_ = ymax;
			}
		else {
			mp_get_output_coordinates(&xmin_, &ymin_, &xmax_, &ymax_);
			}
		draw_mp_picture(fp, data, &vdev, 
			xmin_, ymin_, xmax_, ymax_, draw_pic);
		mp_end_drawing(fp);
		fclose(fp);
		// draw_mp_file(data, &vdev, fname, 1000 /* factor 1000 */, draw_pic);
		printf("finished.\n"); fflush(stdout);
		}
	else {
		printf("draw into file %s:\n", fname2); fflush(stdout);
		if (f_coordinates) {
			mp_set_output_coordinates(xmin, ymin, xmax, ymax);
			}
		draw_mp_file(data, &vdev, fname2, 1000 /* factor 1000 */, draw_pic);
		printf("finished.\n"); fflush(stdout);
		}

	return OK;
}

int main(int argc, char **argv)
{
	INT t0, t1, user_time;

	if( argc < 3 ) {
		printf("usage: t117.out [options] fname output_fname\n");
		printf("reads the file fname.graph and produces output_fname.mp\n");
		printf("options:\n");
		printf("-dots           : draw the subgroups as dots\n");
		printf("-scale fx fy    : scale the coordinates by the factor (fx,fy)\n");
		printf("                : where fx,fy are floating point numbers\n");
		printf("                : default output area is 1500 pixels which \n");
		printf("                : is OK for large graphs\n");
		printf("-append         : appends to an existing file fname.mp\n");
		printf("                : in this case, no header or beginfig\n");
		printf("                : is written into the file\n");
		printf("-coordinates xmin ymin xmax ymax    : use given coordinates for drawing\n");
		printf("-select n i_1 ... i_n               : draw the points (i_1,\\ldots,i_n) special\n");
		printf("                : the numbering of points starts with 0\n");
		return 1;
		}
	discreta_init();
	{
	// INT t0, t1, user_time;
	INT i, j;
	BYTE s[256], *base_fname, *output_fname;
	INT f_dots = FALSE;
	INT f_append = FALSE;
	INT f_coordinates = FALSE;
	INT xmin, ymin, xmax, ymax;
	INT f_select = FALSE, select_n;
	INT *selection = NULL;
	
	t0 = os_ticks();
	{
	double scale_x = SCALE_FACTOR_X;
	double scale_y = SCALE_FACTOR_Y;
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-dots") == 0) {
			f_dots = TRUE;
			}
		else if (strcmp(argv[i], "-append") == 0) {
			f_append = TRUE;
			}
		else if (strcmp(argv[i], "-coordinates") == 0) {
			f_coordinates = TRUE;
			if (argc < i + 4)
				return error("not enough arguments for coordinates");
			sscanf(argv[++i], "%ld", &xmin);
			sscanf(argv[++i], "%ld", &ymin);
			sscanf(argv[++i], "%ld", &xmax);
			sscanf(argv[++i], "%ld", &ymax);
			}
		else if (strcmp(argv[i], "-scale") == 0) {
			if (argc < i + 2)
				return error("not enough arguments for scaling factors");
			sscanf(argv[++i], "%lf", &scale_x);
			sscanf(argv[++i], "%lf", &scale_y);
			}
		else if (strcmp(argv[i], "-select") == 0) {
			f_select = TRUE;
			sscanf(argv[++i], "%ld", &select_n);
			if (argc < i + select_n)
				return error("not enough arguments for selection");
			selection = (INT *) my_malloc(select_n * sizeof(INT), "select");
			for (j = 0; j < select_n; j++) {
				sscanf(argv[++i], "%ld", &selection[j]);
				}
			}
		else
			break;
		}
	base_fname = argv[i++];
	output_fname = argv[i++];
	
	do_it(base_fname, output_fname, f_dots, f_append, 
		f_coordinates, xmin, ymin, xmax, ymax, scale_x, scale_y, 
		f_select, select_n, selection);
	}

	t1 = os_ticks();
	user_time = t1 - t0;
	s[0] = 0;
	print_delta_time_100(user_time, s);
	printf("total computing time: %s\n", s);
	fflush(stdout);
	
	}
	discreta_exit();
	return 0;
}

INT draw_pic(void *draw_data, 
	INT x_, INT y_, INT w_, INT h_, VDEVICE *vdev)
{
	DRAW_LOCAL *p = (DRAW_LOCAL *) draw_data;
	INT *Px, *Py, *Pz;
	INT x0, y0, x1, y1;
	INT *x, *y;
	INT *pts;
	INT rad, rad1;
	INT i;
	INT nb_nodes = p->GM->nb_V + 8;
	INT e1, e2;
	DATA_3D_BOUNDS b;
	
	// printf("vdev->co: %ld %ld %ld %ld\n", vdev->co[0], vdev->co[1], vdev->co[2], vdev->co[3]);

	user2vdev(vdev, 0, 0, 0, &x0, &y0, p->extrema);
	user2vdev(vdev, RADIUS, 0, 0, &x1, &y1, p->extrema);
	rad = ABS(x1 - x0);
	printf("rad = %ld\n", rad);
	
	user2vdev(vdev, 0, 0, 0, &x0, &y0, p->extrema);
	user2vdev(vdev, RADIUS1, 0, 0, &x1, &y1, p->extrema);
	rad1 = ABS(x1 - x0);
	printf("rad1 = %ld\n", rad1);
	
	Px = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic Px");
	Py = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic Py");
	Pz = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic Pz");
	x = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic x");
	y = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic y");
	pts = (INT *) my_malloc(nb_nodes * sizeof(INT), "draw_pic pts");
	

	for (i = 0; i < p->GM->nb_V; i++) {
		Px[i] = (INT) ((double) p->D3->x[i] * p->scale_x);
		Py[i] = (INT) ((double) p->D3->y[i] * p->scale_y);
		}
	Px[p->GM->nb_V + 0] = (INT) ((double) b.xmin * p->scale_x);
	Py[p->GM->nb_V + 0] = (INT) ((double) b.ymin * p->scale_y);
	Px[p->GM->nb_V + 1] = (INT) ((double) (b.xmin + b.xsize) * p->scale_x);
	Py[p->GM->nb_V + 1] = (INT) ((double) b.ymin * p->scale_y);
	Px[p->GM->nb_V + 2] = (INT) ((double) (b.xmin + b.xsize) * p->scale_x);
	Py[p->GM->nb_V + 2] = (INT) ((double) (b.ymin + b.ysize) * p->scale_y);
	Px[p->GM->nb_V + 3] = (INT) ((double) b.xmin * p->scale_x);
	Py[p->GM->nb_V + 3] = (INT) ((double) (b.ymin + b.ysize) * p->scale_y);

	for (i = 0; i < nb_nodes; i++) {
		user2vdev(vdev, Px[i], Py[i], 0, &x[i], &y[i], p->extrema);
		}
	
	for (i = 0; i < p->GM->nb_E; i++) {
		e1 = p->GM->e1[i];
		e2 = p->GM->e2[i];
		draw_line(vdev, x, y, e1, e2);
		}
#if 0
	draw_line(vdev, x, y, p->GM->nb_V + 0, p->GM->nb_V + 1);
	draw_line(vdev, x, y, p->GM->nb_V + 1, p->GM->nb_V + 2);
	draw_line(vdev, x, y, p->GM->nb_V + 2, p->GM->nb_V + 3);
	draw_line(vdev, x, y, p->GM->nb_V + 3, p->GM->nb_V + 0);
#endif

	if (p->f_dots) {
		for (i = 0; i < p->GM->nb_V; i++) {
			if (is_selected(p, i)) {
				draw_dot_black(vdev, x, y, i, rad);
				}
			else {
				draw_dot(vdev, x, y, i, rad);
				}
			}
		}
	else {
		if (p->f_select) {
			for (i = 0; i < p->select_n; i++) {
				draw_dot(vdev, x, y, p->selection[i], rad);
				}
			}
		}
	if (p->D3->f_labels) {
		for (i = 0; i < p->GM->nb_V; i++) {
			draw_label(vdev, p->D3->labels[i], "tr", x, y, i);
			}
		}
	
	my_free(Px);
	my_free(Py);
	my_free(Pz);
	my_free(x);
	my_free(y);
	my_free(pts);

	return OK;
}

INT is_selected(DRAW_LOCAL *p, INT a)
{
	INT i;
	
	if (!p->f_select)
		return FALSE;
	for (i = 0; i < p->select_n; i++) {
		if (p->selection[i] == a)
			return TRUE;
		}
	return FALSE;
}


INT draw_label(VDEVICE *vdev, LABEL *L, INT *x, INT *y)
{
	INT h_align = 1, v_align = 1;
	INT l, i, j;
	BYTE c;
	
	l = strlen(L->align);
	for (i = 0; i < l; i++) {
		c = L->align[i];
		if (c == 'r')
			h_align = 2;
		else if (c == 'l')
			h_align = 0;
		else if (c == 'b')
			v_align = 0;
		else if (c == 't')
			v_align = 2;
		else {
			printf("draw_label: unknown alignment character %c\n", c);
			}
		}
	Vst_alignment(vdev, h_align, v_align, &j, &j);
	V_gtext(vdev, x[L->point_index], y[L->point_index], L->label);
	return OK;
}


#include "t_graph.C"


