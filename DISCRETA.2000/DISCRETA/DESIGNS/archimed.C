#define MAX_V 1000
#define MAX_E 1000
#define MAX_F 1000

typedef struct archimed ARCHIMED;

struct archimed {
	INT nb_V;
	INT Px[MAX_V];
	INT Py[MAX_V];
	INT Pz[MAX_V];
	INT x[MAX_V];
	INT y[MAX_V];
	INT f_vertex_labels;
	INT vertex_label[MAX_V];
	
	INT nb_E;
	INT e1[MAX_E];
	INT e2[MAX_E];
	INT f1[MAX_E];
	INT f2[MAX_E];
	
	INT nb_F;
	INT nb_e[MAX_F];
	INT edge[MAX_F][MAX_E];
	INT neighbour[MAX_F][MAX_F];
	
};

#define EPSILON 400
#define EPSILON_GAUSS 0.00001


INT dode_simum1(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v);
INT archimed_snub_cube(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v);
INT archimed_russian_snub_cube(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v);
INT cubus_simus1(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v);
INT get_automorphism(double P[][3], 
	double Y[3][3], double y0[3], ARCHIMED *X, PERMUTATION_OP p);
INT apply_motion_to_archimed(double P[][3], 
	double Y[3][3], double y0[3], ARCHIMED *X, INT f_v);
INT find_point(double P[][3], INT len, INT i);
INT determine_motion(double P[][3], double Rot[3][3], double x0[3], 
	INT i1, INT i2, INT i3, INT i4, 
	INT j1, INT j2, INT j3, INT j4, 
	double s, INT f_v);
void apply_motion(double P[][3], double Rot[3][3], double x0[3], INT i, INT j);
INT write_graphfile(BYTE *fname, ARCHIMED *A);
INT add_face3(ARCHIMED *A, INT i1, INT i2, INT i3);
INT add_face4(ARCHIMED *A, INT i1, INT i2, INT i3, INT i4);
INT add_face5(ARCHIMED *A, INT i1, INT i2, INT i3, INT i4, INT i5);
INT add_face_n(ARCHIMED *A, INT n, INT *points);
INT find_and_add_edge(ARCHIMED *A, INT i1, INT i2);
INT find_edge(ARCHIMED *A, INT v1, INT v2);
INT find_face_by_two_edges(ARCHIMED *A, INT e1, INT e2);
INT find_faces_at_edge(ARCHIMED *A, INT e, INT *f1, INT *f2);
INT find_face(ARCHIMED *A, INT e, INT *f1, INT *j1, INT *f2, INT *j2);
INT add_edge(ARCHIMED *A, INT v1, INT v2);
void print_points(double P[][3], INT len);
void print_point_i(double P[][3], INT i);
double distance(double P[][3], INT i, INT j);
INT standard_labels(ARCHIMED *A);

#define MAX_N 1000

INT dode_simum1(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v)
{
	double phi, psi, r, h;
	double P[MAX_N][3];
	double Y[3][3], y0[3];
	double Z[3][3], z0[3];
	double x1, x2, x3, y1, y2, y3, t;
	INT i;
	double a;
	
	printf("s=%f\n", s);
	phi = 60. * 360. / 348.;
	psi = 108. * 360. / 348.;
	printf("phi=%f\n", phi);
	printf("psi=%f\n", psi);
	
	a = 2. * sin_grad(phi * .5);
	r = s / a;
	printf("r=%f\n", r);
	h = sqrt(s * s - r * r);
	printf("h=%f\n", h);
	P[0][0] = 0.;
	P[0][1] = 0.;
	P[0][2] = h;
	P[1][0] = r * cos_grad(2. * phi);
	P[1][1] = r * sin_grad(2. * phi);
	P[1][2] = 0.;
	P[4][0] = r * cos_grad(-2. * phi);
	P[4][1] = r * sin_grad(-2. * phi);
	P[4][2] = 0.;
	P[5][0] = r * cos_grad(-1. * phi);
	P[5][1] = r * sin_grad(-1. * phi);
	P[5][2] = 0.;
	P[6][0] = r * cos_grad(0. * phi);
	P[6][1] = r * sin_grad(0. * phi);
	P[6][2] = 0.;
	P[7][0] = r * cos_grad(1. * phi);
	P[7][1] = r * sin_grad(1. * phi);
	P[7][2] = 0.;
	
	x1 = P[1][0] - P[0][0];
	x2 = P[1][1] - P[0][1];
	x3 = P[1][2] - P[0][2];
	y1 = P[4][0] - P[0][0];
	y2 = P[4][1] - P[0][1];
	y3 = P[4][2] - P[0][2];
	t = 2 * (s * cos_grad(72)) / s + 1.;
	P[2][0] = P[0][0] + y1 + t * x1;
	P[2][1] = P[0][1] + y2 + t * x2;
	P[2][2] = P[0][2] + y3 + t * x3;
	P[3][0] = P[0][0] + x1 + t * y1;
	P[3][1] = P[0][1] + x2 + t * y2;
	P[3][2] = P[0][2] + x3 + t * y3;

	determine_motion(P, Y, y0, 
		0, 5, 4, 3, 
		1, 7, 0, 4, s, FALSE);
	apply_motion(P, Y, y0, 7, 9);
	apply_motion(P, Y, y0, 9, 11);
	apply_motion(P, Y, y0, 11, 13);
	apply_motion(P, Y, y0, 6, 8);
	apply_motion(P, Y, y0, 8, 10);
	apply_motion(P, Y, y0, 10, 12);
	apply_motion(P, Y, y0, 12, 14);
	
	X->f_vertex_labels = FALSE;
	for (i = 0; i <= 14; i++) {
		X->Px[i] = (INT)(P[i][0]);
		X->Py[i] = (INT)(P[i][1]);
		X->Pz[i] = (INT)(P[i][2]);
		}
	X->nb_V = 15;
	X->nb_E = 0;
	X->nb_F = 0;
	
	INT e, f1, j1, f2, j2;
	
	add_face5(X, 0, 1, 2, 3, 4);
	add_face3(X, 0, 6, 7);
	add_face3(X, 0, 5, 6);
	add_face3(X, 0, 4, 5);

	standard_labels(X);
	if (f_v) {
		write_graphfile(fname, X);
		}
	
	determine_motion(P, Z, z0, 
		0, 4, 5, 14, 
		9, 10, 2, 11, s, TRUE);

#if 0
	apply_motion(P, Z, z0, 6, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("6 -> %ld\n", idx);

	apply_motion(P, Z, z0, 5, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("5 -> %ld\n", idx);

	apply_motion(P, Z, z0, 0, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("0 -> %ld\n", idx);

	apply_motion(P, Z, z0, 4, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("4 -> %ld\n", idx);
#endif
	
	BYTE fname1[100000];
	BYTE fname2[100000], *pp;
	strcpy(fname1, fname);
	if ((pp = strrchr(fname1, '.')) != NULL)
		*pp = 0;
		
	while (TRUE) {
		INT f_change = FALSE;
		
		for (i = 0; i < 4; i++) {
			printf("apply Y:\n");
			apply_motion_to_archimed(P, Y, y0, X, f_v);	
			strcat(fname1, "_y");
			strcpy(fname2, fname1);
			strcat(fname2, ".graph");
			standard_labels(X);
			printf("nb_V=%ld nb_E=%ld nb_F=%ld\n", X->nb_V, X->nb_E, X->nb_F);
			fflush(stdout);
			if (f_v) {
				write_graphfile(fname2, X);
				}
			}
		
		printf("apply Z:\n");
		if (apply_motion_to_archimed(P, Z, z0, X, f_v))
			f_change = TRUE;
		strcat(fname1, "_z");
		strcpy(fname2, fname1);
		strcat(fname2, ".graph");
		standard_labels(X);
		printf("nb_V=%ld nb_E=%ld nb_F=%ld\n", X->nb_V, X->nb_E, X->nb_F);
		fflush(stdout);
		if (f_v) {
			write_graphfile(fname2, X);
			}


		if (!f_change)
			break;
		}


	for (e = 0; e < X->nb_E; e++) {
		find_face(X, e, &f1, &j1, &f2, &j2);
		X->f1[e] = f1;
		X->f2[e] = f2;
		X->neighbour[f1][j1] = f2;
		X->neighbour[f2][j2] = f1;
		}
	// write_graphfile(fname, X);

	{
	PERMUTATION_OB per;
	
	V->m_il(2);
	get_automorphism(P, Y, y0, X, &per);
	per.swap(V->s_i(0));
	get_automorphism(P, Z, z0, X, &per);
	per.swap(V->s_i(1));
	V->Print();
	}
	return OK;
}

INT archimed_snub_cube(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v)
{
	INT i;
	double r, s2, sr;
	double Y[3][3] = { { 0., 0., -1. }, { 0., 1., 0. }, { 1., 0., 0. } }, 
		y0[3] = { 0., 0., 0. };
	double Z[3][3] = { { 0., -1., 0. }, { 1., 0., 0. }, { 0., 0., 1. } }, 
		z0[3] = { 0., 0., 0. };
	double P[MAX_N][3];
	
	printf("s=%f\n", s);
	r = s * sin_grad(45);
	printf("r=%f\n", r);
	s2 = s * .5;
	sr = s2 + r;
	P[0][0] = s2;
	P[0][1] = -sr;
	P[0][2] = -s2;
	P[1][0] = s2;
	P[1][1] = -sr;
	P[1][2] = +s2;
	P[2][0] = -s2;
	P[2][1] = -sr;
	P[2][2] = +s2;
	P[3][0] = -s2;
	P[3][1] = -sr;
	P[3][2] = -s2;
	P[4][0] = sr;
	P[4][1] = -s2;
	P[4][2] = -s2;
	P[5][0] = sr;
	P[5][1] = -s2;
	P[5][2] = s2;
	P[6][0] = s2;
	P[6][1] = -s2;
	P[6][2] = sr;

	X->f_vertex_labels = FALSE;
	for (i = 0; i <= 6; i++) {
		X->Px[i] = (INT)(P[i][0]);
		X->Py[i] = (INT)(P[i][1]);
		X->Pz[i] = (INT)(P[i][2]);
		}
	X->nb_V = 7;
	X->nb_E = 0;
	X->nb_F = 0;
	
	add_face4(X, 0, 1, 2, 3);
	add_face4(X, 0, 4, 5, 1);
	add_face3(X, 1, 5, 6);

	while (TRUE) {
		INT f_change = FALSE;
		
		printf("apply Y:\n");
		if (apply_motion_to_archimed(P, Y, y0, X, f_v))
			f_change = TRUE;
		
		printf("apply Z:\n");
		if (apply_motion_to_archimed(P, Z, z0, X, f_v))
			f_change = TRUE;

		if (!f_change)
			break;
		}

	INT e, f1, f2, j1, j2;
	
	for (e = 0; e < X->nb_E; e++) {
		find_face(X, e, &f1, &j1, &f2, &j2);
		X->f1[e] = f1;
		X->f2[e] = f2;
		X->neighbour[f1][j1] = f2;
		X->neighbour[f2][j2] = f1;
		}

	{
	PERMUTATION_OB per;
	
	V->m_il(2);
	get_automorphism(P, Y, y0, X, &per);
	per.swap(V->s_i(0));
	get_automorphism(P, Z, z0, X, &per);
	per.swap(V->s_i(1));
	V->Print();
	}
	return OK;
}

INT archimed_russian_snub_cube(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v)
{
	INT i, j;
	double r, s2, sr, s45;
	double Z[3][3] = { { 0., -1., 0. }, { 1., 0., 0. }, { 0., 0., 1. } }, 
		z0[3] = { 0., 0., 0. };
	double W[3][3], w0[3] = { 0., 0., 0. };
	double ZZ[3][3] = { { 0., 0., 0. }, { 0., 0., 0. }, { 0., 0., 1. } }, 
		zz0[3] = { 0., 0., 0. };
	double P[MAX_N][3];
	
	printf("s=%f\n", s);
	s45 = sin_grad(45);
	r = s * s45;
	printf("r=%f\n", r);
	s2 = s * .5;
	sr = s2 + r;
	ZZ[0][0] = s45;
	ZZ[0][1] = s45;
	ZZ[1][0] = -s45;
	ZZ[1][1] = s45;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			W[i][j] = ZZ[i][j];
			}
		}
	W[2][2] = -1.;
	P[0][0] = s2;
	P[0][1] = -sr;
	P[0][2] = -s2;
	P[1][0] = s2;
	P[1][1] = -sr;
	P[1][2] = +s2;
	P[2][0] = -s2;
	P[2][1] = -sr;
	P[2][2] = +s2;
	P[3][0] = -s2;
	P[3][1] = -sr;
	P[3][2] = -s2;
	P[4][0] = sr;
	P[4][1] = -s2;
	P[4][2] = -s2;
	P[5][0] = sr;
	P[5][1] = -s2;
	P[5][2] = s2;
	P[6][0] = s2;
	P[6][1] = -s2;
	P[6][2] = sr;
	P[7][0] = -s2;
	P[7][1] = -s2;
	P[7][2] = sr;
	P[8][0] = s2;
	P[8][1] = s2;
	P[8][2] = sr;
	P[9][0] = -s2;
	P[9][1] = s2;
	P[9][2] = sr;
	
	P[14][0] = s2;
	P[14][1] = -s2;
	P[14][2] = -sr;
	P[15][0] = s2;
	P[15][1] = s2;
	P[15][2] = -sr;
	P[16][0] = -s2;
	P[16][1] = s2;
	P[16][2] = -sr;
	P[17][0] = -s2;
	P[17][1] = -s2;
	P[17][2] = -sr;
	
	apply_motion(P, ZZ, zz0, 14, 10);
	apply_motion(P, ZZ, zz0, 15, 11);
	apply_motion(P, ZZ, zz0, 16, 12);
	apply_motion(P, ZZ, zz0, 17, 13);

	X->f_vertex_labels = FALSE;
	for (i = 0; i <= 13; i++) {
		X->Px[i] = (INT)(P[i][0]);
		X->Py[i] = (INT)(P[i][1]);
		X->Pz[i] = (INT)(P[i][2]);
		}
	X->nb_V = 14;
	X->nb_E = 0;
	X->nb_F = 0;
	
	add_face4(X, 0, 1, 2, 3);
	add_face4(X, 0, 4, 5, 1);
	add_face3(X, 1, 5, 6);
	add_face4(X, 1, 6, 7, 2);
	add_face4(X, 6, 8, 9, 7);
	
	add_face4(X, 10, 11, 12, 13);
	add_face3(X, 0, 3, 10);
	add_face4(X, 0, 4, 11, 10);

	while (TRUE) {
		INT f_change = FALSE;
		
		printf("apply Z:\n");
		if (apply_motion_to_archimed(P, Z, z0, X, f_v))
			f_change = TRUE;

		if (!f_change)
			break;
		}

	INT e, f1, f2, j1, j2;
	
	for (e = 0; e < X->nb_E; e++) {
		find_face(X, e, &f1, &j1, &f2, &j2);
		X->f1[e] = f1;
		X->f2[e] = f2;
		X->neighbour[f1][j1] = f2;
		X->neighbour[f2][j2] = f1;
		}

	{
	PERMUTATION_OB per;
	
	V->m_il(2);
	get_automorphism(P, Z, z0, X, &per);
	per.swap(V->s_i(0));
	get_automorphism(P, W, w0, X, &per);
	per.swap(V->s_i(1));
	V->Print();
	}
	return OK;
}

INT cubus_simus1(double s, BYTE *fname, ARCHIMED *X, VECTOR_OP V, INT f_v)
{
	double phi, psi, r, h;
	double P[MAX_N][3];
	double Rot[3][3], x0[3];
	double Y[3][3], y0[3];
	double Z[3][3], z0[3];
	INT i;
	double a;
	
	printf("s=%f\n", s);
	phi = 20. * 36. / 11.;
	psi = 30. * 36. / 11.;
	printf("phi=%f\n", phi);
	printf("psi=%f\n", psi);
	
	a = 2. * sin_grad(phi * .5);
	r = s / a;
	printf("r=%f\n", r);
	h = sqrt(s * s - r * r);
	printf("h=%f\n", h);
	P[0][0] = 0.;
	P[0][1] = 0.;
	P[0][2] = h;
	P[1][0] = r * cos_grad(0. * phi);
	P[1][1] = r * sin_grad(0. * phi);
	P[1][2] = 0.;
	P[2][0] = r * cos_grad(1. * phi);
	P[2][1] = r * sin_grad(1. * phi);
	P[2][2] = 0.;
	P[3][0] = r * cos_grad(2. * phi);
	P[3][1] = r * sin_grad(2. * phi);
	P[3][2] = 0.;
	P[4][0] = r * cos_grad(-1. * phi);
	P[4][1] = r * sin_grad(-1. * phi);
	P[4][2] = 0.;
	P[5][0] = r * cos_grad(-2. * phi);
	P[5][1] = r * sin_grad(-2. * phi);
	P[5][2] = 0.;
	P[6][0] = P[3][0] + P[5][0];
	P[6][1] = P[3][1] + P[5][1];
	P[6][2] = -h;

	determine_motion(P, Rot, x0, 
		0, 2, 4, 1, 
		2, 1, 3, 0, s, FALSE);
	apply_motion(P, Rot, x0, 6, 9);
	apply_motion(P, Rot, x0, 5, 8);
	apply_motion(P, Rot, x0, 3, 7);
	
	X->f_vertex_labels = FALSE;
	for (i = 0; i <= 9; i++) {
		X->Px[i] = (INT)(P[i][0]);
		X->Py[i] = (INT)(P[i][1]);
		X->Pz[i] = (INT)(P[i][2]);
		}
	X->nb_V = 10;
	X->nb_E = 0;
	X->nb_F = 0;
	
	INT e, f1, j1, f2, j2;
	
	add_face4(X, 0, 3, 6, 5);
	add_face3(X, 0, 1, 2);
	add_face3(X, 0, 2, 3);
	add_face3(X, 0, 1, 4);
	add_face3(X, 0, 4, 5);
	add_face4(X, 2, 7, 9, 8);
	add_face3(X, 2, 3, 8);
	add_face3(X, 1, 2, 7);

	standard_labels(X);
	write_graphfile(fname, X);
	
	determine_motion(P, Y, y0, 
		6, 5, 0, 4, 
		5, 0, 3, 2, s, TRUE);
		
#if 0
	INT idx, nb_V = X->nb_V;

	apply_motion(P, Y, y0, 6, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("6 -> %ld\n", idx);

	apply_motion(P, Y, y0, 5, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("5 -> %ld\n", idx);

	apply_motion(P, Y, y0, 0, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("0 -> %ld\n", idx);

	apply_motion(P, Y, y0, 4, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("4 -> %ld\n", idx);

	apply_motion(P, Y, y0, 1, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("1 -> %ld\n", idx);

	print_points(P, 11);
	double d = distance(P, 8, 10);
	printf("d 8, 10=%f\n", d);
#endif

	determine_motion(P, Z, z0, 
		6, 5, 0, 4, 
		8, 2, 7, 1, s, TRUE);

#if 0
	apply_motion(P, Z, z0, 6, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("6 -> %ld\n", idx);

	apply_motion(P, Z, z0, 5, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("5 -> %ld\n", idx);

	apply_motion(P, Z, z0, 0, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("0 -> %ld\n", idx);

	apply_motion(P, Z, z0, 4, nb_V);
	idx = find_point(P, nb_V, nb_V);
	printf("4 -> %ld\n", idx);
#endif
	
	BYTE fname1[100000];
	BYTE fname2[100000], *pp;
	strcpy(fname1, fname);
	if ((pp = strrchr(fname1, '.')) != NULL)
		*pp = 0;
	
	while (TRUE) {
		INT f_change = FALSE;
		
		printf("apply Y:\n");
		if (apply_motion_to_archimed(P, Y, y0, X, f_v))
			f_change = TRUE;
		strcat(fname1, "_y");
		strcpy(fname2, fname1);
		strcat(fname2, ".graph");
		standard_labels(X);
		// write_graphfile(fname2, X);
		
		printf("apply Z:\n");
		if (apply_motion_to_archimed(P, Z, z0, X, f_v))
			f_change = TRUE;
		strcat(fname1, "_z");
		strcpy(fname2, fname1);
		strcat(fname2, ".graph");
		standard_labels(X);
		// write_graphfile(fname2, X);


		if (!f_change)
			break;
		}


	for (e = 0; e < X->nb_E; e++) {
		find_face(X, e, &f1, &j1, &f2, &j2);
		X->f1[e] = f1;
		X->f2[e] = f2;
		X->neighbour[f1][j1] = f2;
		X->neighbour[f2][j2] = f1;
		}
	// write_graphfile(fname, X);

	{
	PERMUTATION_OB per;
	
	V->m_il(2);
	get_automorphism(P, Y, y0, X, &per);
	per.swap(V->s_i(0));
	get_automorphism(P, Z, z0, X, &per);
	per.swap(V->s_i(1));
	V->Print();
	}
	return OK;
}

INT get_automorphism(double P[][3], 
	double Y[3][3], double y0[3], ARCHIMED *X, PERMUTATION_OP p)
{
	INT i, j, nb_V;
	
	nb_V = X->nb_V;
	p->m_il(nb_V);
	for (i = 0; i < nb_V; i++) {
		apply_motion(P, Y, y0, i, nb_V);
		j = find_point(P, nb_V, nb_V);
		if (j == -1)
			return error("get_automorphism() point not found");
		p->m_ii(i, j + 1);
		}
	return OK;
}

INT apply_motion_to_archimed(double P[][3], 
	double Y[3][3], double y0[3], ARCHIMED *X, INT f_v)
{
	INT i, j, nb_e, e, ee, e1, e2, v1, v2, v3, v4, idx1, idx2, last_vertex;
	INT nb_V, nb_F, old_vertices[MAX_N], new_vertices[MAX_N];
	INT f;
	INT face_added = FALSE;
	
	nb_V = X->nb_V;
	nb_F = X->nb_F;
	for (i = 0; i < nb_F; i++) {
		if (f_v) {
			printf("face %ld:\n", i);
			}
		nb_e = X->nb_e[i];
		for (j = 0; j < nb_e; j++) {
			e = X->edge[i][j];
			v1 = X->e1[e];
			v2 = X->e2[e];
			apply_motion(P, Y, y0, v1, nb_V);
			idx1 = find_point(P, nb_V, nb_V);
			if (idx1 == -1) {
				X->Px[nb_V] = (INT)(P[nb_V][0]);
				X->Py[nb_V] = (INT)(P[nb_V][1]);
				X->Pz[nb_V] = (INT)(P[nb_V][2]);
				idx1 = nb_V;
				if (f_v) {
					printf("point %ld maps to unknown point, this becomes new point %ld\n", v1, idx1);
					}
				nb_V++;
				}
			else {
				if (f_v) {
					printf("point %ld maps to %ld\n", v1, idx1);
					}
				}
			apply_motion(P, Y, y0, v2, nb_V);
			idx2 = find_point(P, nb_V, nb_V);
			if (idx2 == -1) {
				X->Px[nb_V] = (INT)(P[nb_V][0]);
				X->Py[nb_V] = (INT)(P[nb_V][1]);
				X->Pz[nb_V] = (INT)(P[nb_V][2]);
				idx2 = nb_V;
				if (f_v) {
					printf("point %ld maps to unknown point, this becomes new point %ld\n", v2, idx2);
					}
				nb_V++;
				}
			else {
				if (f_v) {
					printf("point %ld maps to %ld\n", v2, idx2);
					}
				}
			if (f_v) {
				printf("edge between %ld and %ld is ", idx1, idx2);
				}
			ee = find_and_add_edge(X, idx1, idx2);
			if (f_v) {
				printf("%ld\n", ee);
				}
			if (j == 0)
				e1 = ee;
			else if (j == 1)
				e2 = ee;

			}
		if (f_v) {
			printf("find face for edges %ld and %ld ", e1, e2);
			}
		f = find_face_by_two_edges(X, e1, e2);
		if (f_v) {
			printf("%ld\n", f);
			}
		if (f == -1) {
			if (f_v) {
				printf("old face: ");
				}
			e = X->edge[i][0];
			v1 = X->e1[e];
			v2 = X->e2[e];
			for (j = 1; j < nb_e; j++) {
				e = X->edge[i][j];
				v3 = X->e1[e];
				v4 = X->e2[e];
				if (j == 1) {
					if (v3 == v2) {
						old_vertices[0] = v1;
						old_vertices[1] = v2;
						last_vertex = v4;
						}
					else if (v3 == v1) {
						old_vertices[0] = v2;
						old_vertices[1] = v1;
						last_vertex = v4;
						}
					else if (v4 == v2) {
						old_vertices[0] = v1;
						old_vertices[1] = v2;
						last_vertex = v3;
						}
					else if (v4 == v1) {
						old_vertices[0] = v2;
						old_vertices[1] = v1;
						last_vertex = v3;
						}
					else {
						return error("error edges not adjacent! (j=1)");
						}
					}
				else {
					if (v3 == last_vertex) {
						old_vertices[j] = v3;
						last_vertex = v4;
						}
					else if (v4 == last_vertex) {
						old_vertices[j] = v4;
						last_vertex = v3;
						}
					else {
						return error("error edges not adjacent!");
						}
					}
				}
			if (f_v) {
				for (j = 0; j < nb_e; j++) {
					printf(" %ld", old_vertices[j]);
					}
				printf("\n");
				}
			for (j = 0; j < nb_e; j++) {
				apply_motion(P, Y, y0, old_vertices[j], nb_V);
				idx1 = find_point(P, nb_V, nb_V);
				if (idx1 == -1) {
					return error("error: point not found");
					}
				new_vertices[j] = idx1;
				}
			if (f_v) {
				printf("new face: ");
				for (j = 0; j < nb_e; j++) {
					printf(" %ld", new_vertices[j]);
					}
				printf("\n");
				}
			add_face_n(X, nb_e, new_vertices);
			face_added = TRUE;
			}
		}
	for (i = X->nb_V; i < nb_V; i++) {
		X->Px[i] = (INT)(P[i][0]);
		X->Py[i] = (INT)(P[i][1]);
		X->Pz[i] = (INT)(P[i][2]);
		}
	X->nb_V = nb_V;
	return face_added;
}

INT find_point(double P[][3], INT len, INT i)
{
	INT j;
	double d;
	
	for (j = 0; j < len; j++) {
		d = distance(P, i, j);
		if (d < EPSILON) {
			return j;
			}
		}
	return -1;
}


INT determine_motion(double P[][3], double Rot[3][3], double x0[3], 
	INT i1, INT i2, INT i3, INT i4, 
	INT j1, INT j2, INT j3, INT j4, 
	double s, INT f_v)
{
	double a, c, d, tmp, sv = 1. / s;
	INT i, j, k;
	SYM_OB ob1, ob2;
	double A[12][24], B[12], B1[12];
	
	for (i = 0; i < 12; i++) {
		for (j = 0; j < 12; j++) {
			A[i][j] = 0.;
			}
		}
	A[0][0] = sv * P[i1][0];
	A[0][1] = sv * P[i1][1];
	A[0][2] = sv * P[i1][2];
	A[1][3] = sv * P[i1][0];
	A[1][4] = sv * P[i1][1];
	A[1][5] = sv * P[i1][2];
	A[2][6] = sv * P[i1][0];
	A[2][7] = sv * P[i1][1];
	A[2][8] = sv * P[i1][2];
	A[0][9] = 1.;
	A[1][10] = 1.;
	A[2][11] = 1.;

	A[3][0] = sv * P[i2][0];
	A[3][1] = sv * P[i2][1];
	A[3][2] = sv * P[i2][2];
	A[4][3] = sv * P[i2][0];
	A[4][4] = sv * P[i2][1];
	A[4][5] = sv * P[i2][2];
	A[5][6] = sv * P[i2][0];
	A[5][7] = sv * P[i2][1];
	A[5][8] = sv * P[i2][2];
	A[3][9] = 1.;
	A[4][10] = 1.;
	A[5][11] = 1.;

	A[6][0] = sv * P[i3][0];
	A[6][1] = sv * P[i3][1];
	A[6][2] = sv * P[i3][2];
	A[7][3] = sv * P[i3][0];
	A[7][4] = sv * P[i3][1];
	A[7][5] = sv * P[i3][2];
	A[8][6] = sv * P[i3][0];
	A[8][7] = sv * P[i3][1];
	A[8][8] = sv * P[i3][2];
	A[6][9] = 1.;
	A[7][10] = 1.;
	A[8][11] = 1.;

	A[9][0] = sv * P[i4][0];
	A[9][1] = sv * P[i4][1];
	A[9][2] = sv * P[i4][2];
	A[10][3] = sv * P[i4][0];
	A[10][4] = sv * P[i4][1];
	A[10][5] = sv * P[i4][2];
	A[11][6] = sv * P[i4][0];
	A[11][7] = sv * P[i4][1];
	A[11][8] = sv * P[i4][2];
	A[9][9] = 1.;
	A[10][10] = 1.;
	A[11][11] = 1.;


	B[0] = 	sv * P[j1][0];
	B[1] = 	sv * P[j1][1];
	B[2] = 	sv * P[j1][2];
	B[3] = 	sv * P[j2][0];
	B[4] = 	sv * P[j2][1];
	B[5] = 	sv * P[j2][2];
	B[6] = 	sv * P[j3][0];
	B[7] = 	sv * P[j3][1];
	B[8] = 	sv * P[j3][2];
	B[9] = 	sv * P[j4][0];
	B[10] = sv * P[j4][1];
	B[11] = sv * P[j4][2];
	
	for (i = 0; i < 12; i++) {
		for (j = 0; j < 12; j++) {
			if (i == j)
				A[i][12 + j] = 1.;
			else 
				A[i][12 + j] = 0.;
			}
		}
#if 0
	for (i = 0; i < 12; i++) {
		for (j = 0; j < 12; j++) {
			printf("%f\t", A[i][j]);
			}
		printf("\n");
		}
#endif

	for (j = 0; j < 12; j++) {
		for (i = j; i < 12; i++) {
			if (ABS(A[i][j]) > EPSILON_GAUSS) {
				if (i != j) {
					for (k = j; k < 24; k++) {
						tmp = A[j][k];
						A[j][k] = A[i][k];
						A[i][k] = tmp;
						}
					}
				break;
				}
			}
		if (i == 12) {
			return error("no pivot element!");
			}
		a = 1. / A[j][j];
		for (k = j; k < 24; k++) {
			A[j][k] *= a;
			}
		for (i = j + 1; i < 12; i++) {
			c = A[i][j];
			for (k = 0; k < 24; k++) {
				A[i][k] = A[i][k] - c * A[j][k];
				}
			}
		}
	for (j = 11; j >= 0; j--) {
		for (i = j - 1; i >= 0; i--) {
			c = A[i][j];
			for (k = j; k < 24; k++) {
				A[i][k] = A[i][k] - c * A[j][k];
				}
			}
		}
	for (i = 0; i < 12; i++) {
		B1[i] = 0.;
		}
	for (i = 0; i < 12; i++) {
		for (j = 0; j < 12; j++) {
			B1[i] += A[i][12 + j] * B[j];
			}
		}
	Rot[0][0] = B1[0];
	Rot[0][1] = B1[1];
	Rot[0][2] = B1[2];
	Rot[1][0] = B1[3];
	Rot[1][1] = B1[4];
	Rot[1][2] = B1[5];
	Rot[2][0] = B1[6];
	Rot[2][1] = B1[7];
	Rot[2][2] = B1[8];
	x0[0] = s * B1[9];
	x0[1] = s * B1[10];
	x0[2] = s * B1[11];

#if 0	
	printf("Rot:\n");
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			printf("%f\t", Rot[i][j]);
			}
		printf("\n");
		}
	printf("x0:\n");
	for (j = 0; j < 3; j++) {
		printf("%f\t", x0[j]);
		}
	printf("\n");
#endif

	d = det3((double *) Rot);
	
	// printf("det = %f\n", d);
	return OK;	
}

void apply_motion(double P[][3], double Rot[3][3], double x0[3], INT i, INT j)
{
	INT ii, jj;
	double x[3] = {0., 0., 0.};
	
	for (ii = 0; ii < 3; ii++) {
		for (jj = 0; jj < 3; jj++) {
			x[ii] += Rot[ii][jj] * P[i][jj];
			}
		}
	for (ii = 0; ii < 3; ii++) {
		x[ii] += x0[ii];
		P[j][ii] = x[ii];
		}
}

INT write_graphfile(BYTE *fname, ARCHIMED *A)
{
	FILE *fp;
	INT i;
	
	fp = fopen(fname, "w");
	fprintf(fp, "NUMBER_OF_VERTICES %ld\n", A->nb_V);
	fprintf(fp, "NUMBER_OF_EDGES %ld\n", A->nb_E);
	fprintf(fp, "VERTICES\n");
	for (i = 0; i < A->nb_V; i++) {
		fprintf(fp, "%ld %ld %ld\n", 
			A->Px[i], 
			A->Py[i], 
			A->Pz[i]);
		}
	fprintf(fp, "EDGES\n");
	for (i = 0; i < A->nb_E; i++) {
		fprintf(fp, "%ld %ld\n", A->e1[i], A->e2[i]);
		}
	fprintf(fp, "\n");
	if (A->f_vertex_labels) {
		fprintf(fp, "VERTEX-LABELS\n");
		for (i = 0; i < A->nb_V; i++) 
		{
			fprintf(fp, "%ld\n", A->vertex_label[i]); 
		}
	}
	fclose(fp);
	return OK;
}

INT add_face3(ARCHIMED *A, INT i1, INT i2, INT i3)
{
	INT e;
	
	A->nb_e[A->nb_F] = 3;
	
	e = find_and_add_edge(A, i1, i2);
	A->edge[A->nb_F][0] = e;
	
	e = find_and_add_edge(A, i2, i3);
	A->edge[A->nb_F][1] = e;
	
	e = find_and_add_edge(A, i3, i1);
	A->edge[A->nb_F][2] = e;
	
	A->nb_F++;
	return OK;
}

INT add_face4(ARCHIMED *A, INT i1, INT i2, INT i3, INT i4)
{
	INT e;
	
	A->nb_e[A->nb_F] = 4;
	
	e = find_and_add_edge(A, i1, i2);
	A->edge[A->nb_F][0] = e;
	
	e = find_and_add_edge(A, i2, i3);
	A->edge[A->nb_F][1] = e;
	
	e = find_and_add_edge(A, i3, i4);
	A->edge[A->nb_F][2] = e;
	
	e = find_and_add_edge(A, i4, i1);
	A->edge[A->nb_F][3] = e;
	
	A->nb_F++;
	return OK;
}

INT add_face5(ARCHIMED *A, INT i1, INT i2, INT i3, INT i4, INT i5)
{
	INT e;
	
	A->nb_e[A->nb_F] = 5;
	
	e = find_and_add_edge(A, i1, i2);
	A->edge[A->nb_F][0] = e;
	
	e = find_and_add_edge(A, i2, i3);
	A->edge[A->nb_F][1] = e;
	
	e = find_and_add_edge(A, i3, i4);
	A->edge[A->nb_F][2] = e;
	
	e = find_and_add_edge(A, i4, i5);
	A->edge[A->nb_F][3] = e;
	
	e = find_and_add_edge(A, i5, i1);
	A->edge[A->nb_F][4] = e;
	
	A->nb_F++;
	return OK;
}

INT add_face_n(ARCHIMED *A, INT n, INT *points)
{
	INT e, j;
	
	A->nb_e[A->nb_F] = n;
	
	for (j = 0; j < n - 1; j++) {
		e = find_and_add_edge(A, points[j], points[j + 1]);
		A->edge[A->nb_F][j] = e;
		}
	e = find_and_add_edge(A, points[n - 1], points[0]);
	A->edge[A->nb_F][n - 1] = e;
	
	A->nb_F++;
	return OK;
}

INT find_and_add_edge(ARCHIMED *A, INT i1, INT i2)
{
	INT e;
	
	e = find_edge(A, i1, i2);
	if (e == -1)
		e = add_edge(A, i1, i2);
	return e;
}

INT find_edge(ARCHIMED *A, INT v1, INT v2)
{
	INT i, e1, e2;
	
	for (i = 0; i < A->nb_E; i++) {
		e1 = A->e1[i];
		e2 = A->e2[i];
		if ((e1 == v1 && e2 == v2) || (e1 == v2 && e2 == v1)) {
			// printf("found edge %ld (v1=%ld,v2=%ld) (e1=%ld,e2=%ld)\n", i, v1, v2, e1, e2);
			return i;
			}
		}
	return -1;
}

INT find_face_by_two_edges(ARCHIMED *A, INT e1, INT e2)
{
	INT f1, f2, f3, f4;
	
	find_faces_at_edge(A, e1, &f1, &f2);
	find_faces_at_edge(A, e2, &f3, &f4);
	if (f1 != -1) {
		if (f1 == f3 || f1 == f4) {
			return f1;
			}
		}
	if (f2 != -1) {
		if (f2 == f3 || f2 == f4) {
			return f2;
			}
		}
	return -1;
}

INT find_faces_at_edge(ARCHIMED *A, INT e, INT *f1, INT *f2)
{
	INT i, j, l, n = 0;
	
	*f1 = -1;
	*f2 = -1;
	for (i = 0; i < A->nb_F; i++) {
		l = A->nb_e[i];
		for (j = 0; j < l; j++) {
			if (A->edge[i][j] == e) {
				if (n == 0) {
					*f1 = i;
					n++;
					}
				else if (n == 1) {
					*f2 = i;
					n++;
					}
				else {
					printf("e=%ld f1=%ld f2=%ld i=%ld\n", e, *f1, *f2, i);
					return error("too many faces for this edge");
					}
				}
			}
		}
	return OK;
}

INT find_face(ARCHIMED *A, INT e, INT *f1, INT *j1, INT *f2, INT *j2)
{
	INT i, j, l, ff1 = -1, ff2 = -1;
	
	for (i = 0; i < A->nb_F; i++) {
		l = A->nb_e[i];
		for (j = 0; j < l; j++) {
			if (A->edge[i][j] == e) {
				if (ff1 != -1) {
					ff2 = i;
					*j2 = j;
					}
				else {
					ff1 = i;
					*j1 = j;
					}
				}
			}
		}
	if (ff1 == -1 || ff2 == -1) {
		printf("e = %ld\n", e);
		for (i = 0; i < A->nb_F; i++) {
			l = A->nb_e[i];
			for (j = 0; j < l; j++) {
				printf("%ld ", A->edge[i][j]);
				}
			printf("\n");
			}
		}
	if (ff1 == -1)
		return error("face not found");
	if (ff2 == -1)
		return error("face not found");
	*f1 = ff1;
	*f2 = ff2;
	return OK;
}

INT add_edge(ARCHIMED *A, INT v1, INT v2)
{
	INT i;
	
	A->e1[A->nb_E] = v1;
	A->e2[A->nb_E] = v2;
	A->f1[A->nb_E] = -1;
	A->f2[A->nb_E] = -1;
	i = A->nb_E;
	A->nb_E++;
	// printf("add_edge: edge %ld = (%ld, %ld) = (%ld, %ld)added\n", i, v1, v2, A->e1[i], A->e2[i]);
	return i;
}

void print_points(double P[][3], INT len)
{
	INT i;
	
	for (i = 0; i < len; i++) {
		print_point_i(P, i);
		}
}

void print_point_i(double P[][3], INT i)
{
	printf("%ld: %f\t%f\t%f\n", i, P[i][0], P[i][1], P[i][2]);
}

double distance(double P[][3], INT i, INT j)
{
	double a, d = 0.;
	INT ii;
	
	for (ii = 0; ii < 3; ii++) {
		a = P[j][ii] - P[i][ii];
		d += a * a;
		}
	d = sqrt(d);
	return d;
}

INT standard_labels(ARCHIMED *A)
{
	INT i;
	
	for (i = 0; i < A->nb_V; i++) {
		A->vertex_label[i] = i;
		}
	A->f_vertex_labels = TRUE;
	return OK;
}

