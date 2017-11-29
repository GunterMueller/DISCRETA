// t116.C
//
// 
// create_graph_from_lattice
// 
// reads a subgroup lattice file <fname>.sgl.dsc
// and produces .graph file <fname>_<lines>.graph containing the 3-d 
// placement of the lattice.
// here, lines stands for
// 0 = Asup, 1 = Ainf, 2 = all, 3 = poset
//
// (extension of t079.C)
//
// Anton Betten
// Bayreuth, 17.03.1999

#include <DISCRETA/discreta.h>
#include <DISCRETA/divs.h>
#include <DISCRETA/sgl.h>
#include <DISCRETA/perm.h> // for symmetric_generators_pp
#include <stdlib.h>


// t_graph.C:

#define STANDARD_SCALE_X 1500
#define STANDARD_SCALE_Y 1500
#define STANDARD_SCALE_Z 1500

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

INT get_vertex_koordinates(SGL_OP L, 
	VECTOR_OP Px, VECTOR_OP Py, VECTOR_OP Os, 
	INT *x, INT *y, INT *z, 
	INT orbit, INT r, INT f_2D)
{
	SGO_OP O;
	INT l, o, nb_groups;
	double angle, rad;
	
	L->orbit2lo(orbit, &l, &o);
	O = L->s_theOrbits_ij(l, o);
	nb_groups = O->s_o_len_i();
	if (nb_groups == 1) {
		*x = Px->s_ii(orbit);
		*y = 0;
		*z = Py->s_ii(orbit);
		return OK;
		}
	rad = ((double) Os->s_ii(orbit)) * 0.3;
	if (f_2D) {
		double dx;
		
		dx = 2 * rad / (nb_groups - 1);
		*x = Px->s_ii(orbit) + (INT)(- rad + r * dx);
		*y = 0;
		*z = Py->s_ii(orbit);
		}
	else {
		angle = 360. * r / (double) nb_groups;
		*x = Px->s_ii(orbit) + (INT)(rad * cos_grad(angle));
		*y = (INT)(rad * sin_grad(angle));
		*z = Py->s_ii(orbit);
		}
	return OK;
}

INT create_graph_from_lattice(GRAPH_MULTI **gm, DATA_3D **d3, 
	SGL_OP L, MATRIX_OP Acover, VECTOR_OP Px, VECTOR_OP Py, VECTOR_OP Os, 
	INT f_2D, INT lines)
{
	GRAPH_MULTI *GM;
	DATA_3D *D3;
	INT nb_V = 0, nb_E = 0, i, j, e, orbit, l, o;
	INT orbit1, l1, o1;
	INT orbit2, l2, o2;
	INT nb_groups, nb_groups1, nb_groups2;
	INT r, r1, r2, f_is_subgroup;
	INT nb_Layers, nb_Orbits;
	SGO_OP O, O1, O2;
	INT x, y, z, i0, i1, i2;
	INT r1_max, r2_max;
	
	nb_Layers = L->s_nb_Layers_i();
	nb_V = L->s_total_nb_subgroups_i();
	nb_Orbits = Px->s_li();
	if (nb_Orbits != L->s_total_nb_orbits_i())
		return error("create_graph_from_lattice() nb_Orbits != L->s_total_nb_orbits_i()");
	for (orbit = 0; orbit < nb_Orbits; orbit++) {
		L->orbit2lo(orbit, &l, &o);
		O = L->s_theOrbits_ij(l, o);
		nb_groups = O->s_o_len_i();
		if (nb_groups > 2)
			nb_E += nb_groups;
		else
			nb_E += nb_groups - 1;
		}
#if 1
	for (orbit1 = 0; orbit1 < nb_Orbits; orbit1++) {
		L->orbit2lo(orbit1, &l1, &o1);
		O1 = L->s_theOrbits_ij(l1, o1);
		nb_groups1 = O1->s_o_len_i();
		for (orbit2 = 0; orbit2 < nb_Orbits; orbit2++) {
			L->orbit2lo(orbit2, &l2, &o2);
			O2 = L->s_theOrbits_ij(l2, o2);
			nb_groups2 = O2->s_o_len_i();
			if (orbit2 <= orbit1)
				continue;
			if (Acover->s_iji(orbit1, orbit2)) {
				if (lines == 0) {
					r1_max = nb_groups1;
					r2_max = 1;
					}
				else if (lines == 1) {
					r1_max = 1;
					r2_max = nb_groups2;
					}
				else if (lines == 2) {
					r1_max = nb_groups1;
					r2_max = nb_groups2;
					}
				else {
					return error("lines must be 0, 1 or 2");
					}
				for (r1 = 0; r1 < r1_max; r1++) {
					for (r2 = 0; r2 < r2_max; r2++) {
						L->IsSubgroup(l1, o1, r1, 
							l2, o2, r2, &f_is_subgroup);
						if (f_is_subgroup) {
							nb_E++;
							}
						}
					}
				}
			}
		}
#endif
	printf("create_graph_from_lattice(): nb_V = %ld nb_E = %ld\n", nb_V, nb_E);
	D3 = open_data_3D(nb_V);
	GM = open_graph_multi(nb_V, nb_E, FALSE /* f_multi */);
	i = 0;
	e = 0;
	for (orbit = 0; orbit < nb_Orbits; orbit++) {
		L->orbit2lo(orbit, &l, &o);
		O = L->s_theOrbits_ij(l, o);
		nb_groups = O->s_o_len_i();
		for (r = 0; r < nb_groups; r++, i++) {
			get_vertex_koordinates(L, Px, Py, Os, 
				&x, &y, &z, orbit, r, f_2D);
			// printf("x=%ld y=%ld z=%ld\n", x, y, z); fflush(stdout);
			D3->x[i] = x;
			D3->y[i] = y;
			D3->z[i] = z;
			// printf("orbit=%ld l=%ld o=%ld r=%ld\n", orbit, l, o, r);
			if (r == 0) {
				i0 = i;
				}
			else {
				GM->e1[e] = i;
				GM->e2[e] = i - 1;
				GM->mult[e] = 1;
				// printf("edge (i, i-1)=(%ld,%ld)\n", i, i - 1);
				e++;
				}
			if (r == nb_groups - 1 && nb_groups > 2) {
				GM->e1[e] = i;
				GM->e2[e] = i0;
				GM->mult[e] = 1;
				// printf("edge (i, i0)=(%ld,%ld)\n", i, i0);
				e++;
				}
			}
		}
	if (i != nb_V)
		return error("i != nb_V");
	printf("vertices initialized\n"); fflush(stdout);

#if 1
	i1 = 0;
	for (orbit1 = 0; orbit1 < nb_Orbits; orbit1++, i1 += nb_groups1) {
		L->orbit2lo(orbit1, &l1, &o1);
		O1 = L->s_theOrbits_ij(l1, o1);
		nb_groups1 = O1->s_o_len_i();
		for (orbit2 = orbit1, i2 = i1; orbit2 < nb_Orbits; orbit2++, i2 += nb_groups2) {
			L->orbit2lo(orbit2, &l2, &o2);
			O2 = L->s_theOrbits_ij(l2, o2);
			nb_groups2 = O2->s_o_len_i();
			if (orbit2 <= orbit1)
				continue;
			if (Acover->s_iji(orbit1, orbit2)) {
				if (lines == 0) {
					r1_max = nb_groups1;
					r2_max = 1;
					}
				else if (lines == 1) {
					r1_max = 1;
					r2_max = nb_groups2;
					}
				else if (lines == 2) {
					r1_max = nb_groups1;
					r2_max = nb_groups2;
					}
				else {
					return error("lines must be 0, 1 or 2");
					}
				for (r1 = 0; r1 < r1_max; r1++) {
					for (r2 = 0; r2 < r2_max; r2++) {
						L->IsSubgroup(l1, o1, r1, 
							l2, o2, r2, &f_is_subgroup);
						if (f_is_subgroup) {
							GM->e1[e] = i1 + r1;
							GM->e2[e] = i2 + r2;
							GM->mult[e] = 1;
							e++;
							}
						}
					}
				}
			}
		}
#endif

	if (e != nb_E)
		return error("e != nb_E");
	*d3 = D3;
	*gm = GM;
	return OK;
}

INT do_sgl(BYTE *fname, INT f_2D, INT lines, INT f_verbose, INT size_x, INT size_y, INT f_upside_down)
{
	BYTE str[1024];
	VECTOR_OB V;
	VECTOR_OP gen, G, Inn_gen;
	STRING_OP g_label, g_label_tex;
	SGL_OP L;	
	// SG_LATTICE_LOCAL *sgll;
	// BYTE cmd[1000];
	// BYTE path[1000];
	MATRIX_OB Asup, Ainf, D, M, B, Acover;
	VECTOR_OB nl, orbit_size;
	SYM_OB d;
	VECTOR_OB Px, Py, Os;
	GRAPH_MULTI *GM = NULL;
	DATA_3D *D3 = NULL;

	sprintf(str, "%s.sgl.dsc", fname);
	V.load(str);
	if (V.s_li() < 6) 
		error("V too short");
	gen = (VECTOR_OP) V.s_i(0);
	G = (VECTOR_OP) V.s_i(1);
	Inn_gen = (VECTOR_OP) V.s_i(2);
	g_label = (STRING_OP) V.s_i(3);
	g_label_tex = (STRING_OP) V.s_i(4);
	L = (SGL_OP) V.s_i(5);
	printf("read a subgroup lattice of a group with ");
	printf("%ld generators\n", gen->s_li());
	printf("%ld elements\n", G->s_li());
	printf("%ld generators for Inn(G)\n", Inn_gen->s_li());
	printf("g_label = %s\n", g_label->s_str());
	printf("g_label_tex = %s\n", g_label_tex->s_str());
	
	// sprintf(path, "sgl_%ld_", lines);

	L->Burnside_info(&orbit_size, &Asup, &Ainf, &D, &M, &B, &d, f_verbose);
	Asup.Asup2Acover(&Acover);
	if (f_verbose) {
		printf("Acover =\n");
		Acover.Print();
		}
	Acover.Acover2nl(&nl);
	if (f_verbose) {
		printf("nl =\n");
		nl.println();
		}
	place_lattice(&nl, &orbit_size, size_x, size_y, 
		&Px, &Py, &Os, 
		f_upside_down, f_verbose);
	
	create_graph_from_lattice(&GM, &D3, L, &Acover, &Px, &Py, &Os, f_2D, lines);
	sprintf(str, "%s_%ld.graph", fname, lines);
	save_graph_multi_to_file(str, GM, D3);

#if 0
	vbp(nl, orbit_size, 
		SGL_VBP_X_PIX, SGL_VBP_Y_PIX, 
		vbp_plaz, vbp_o_dx, TRUE /* f_upside_down */);
	p->extrema[0] = 0.;
	p->extrema[1] = (double) SGL_VBP_X_PIX /* * 2. */;
	p->extrema[2] = 0.;
	p->extrema[3] = (double) SGL_VBP_Y_PIX /* * 2. */;
	p->extrema[4] = -1.;
	p->extrema[5] = 1.;

	sgll = open_sgll(L, 
		FALSE /* f_L_allocated */, 
		FALSE /* f_show_generators */, 
		FALSE /* f_with_perm */, 
		lines /* draw_lines_type */, 
		TRUE /* f_draw_sgl */, 
		TRUE /* f_burnside_tex */, 
		path /* path_name */, 
		g_label->s_str(), 
		g_label_tex->s_str(), 
		TRUE /* f_verbose */);
		/* draw_lines_type: 0 = Asup, 1 = Ainf, 2 = all, 3 = poset */
	free_sgll(sgll);

	sprintf(cmd, "mp %s%s.mp", path, g_label->s_str());
	printf("calling system: %s\n", cmd);
	fflush(stdout);
	call_system(cmd);
	
	sprintf(cmd, "latex %s%s.tex", path, g_label->s_str());
	printf("calling system: %s\n", cmd);
	fflush(stdout);
	call_system(cmd);
	
	sprintf(cmd, "dvips %s%s.dvi -o", path, g_label->s_str());
	printf("calling system: %s\n", cmd);
	fflush(stdout);
	call_system(cmd);
	
	sprintf(cmd, "ghostview %s%s.ps", path, g_label->s_str());
	printf("calling system: %s\n", cmd);
	fflush(stdout);
	call_system(cmd);
#endif
	
	return OK;
}

int main(int argc, char **argv)
{
	INT t0, t1, user_time;
	BYTE s[256];
	INT f_verbose = FALSE;
	INT size_x = 500;
	INT size_y = 500;
	INT f_upside_down = TRUE;

	if( argc < 3 ) {
		printf("usage: t116.out [options] fname lines\n");
		printf("reads a subgroup lattice file <fname>.sgl.dsc\n");
		printf("and produces a latex file <fname>_<lines>.graph\n");
		printf("the argument lines determines the drawing of the lattice:\n");
		printf("lines = 0   : draw the Asup picture\n");
		printf("lines = 1   : draw the Ainf picture\n");
		printf("lines = 2   : draw the full lattice\n");
		printf("lines = 3   : draw only the poset of orbits\n");
		printf("options:\n");
		printf("-2D         : 2D placement\n");
		return 1;
		}
	discreta_init();
	{
	// INT t0, t1, user_time;
	INT i, lines;
	BYTE s[256], *fname;
	INT f_2D = FALSE;
	
	t0 = os_ticks();
	{
	for (i = 1; i < argc - 2; i++) {
		if (strcmp(argv[i], "-2D") == 0) {
			f_2D = TRUE;
			}
		else
			break;
		}
	fname = argv[i++];
	sscanf(argv[i++], "%ld", &lines);
	
	do_sgl(fname, f_2D, lines, f_verbose, size_x, size_y, f_upside_down);
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

#include "t_graph.C"



