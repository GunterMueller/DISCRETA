/* solid.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten, Evi Haberberger, 1999
 */


#include <DISCRETA/discreta.h>
#include <DISCRETA/graphics.h>
#include <stdlib.h>

#ifdef SOLID_TRUE

#include <DISCRETA/perm.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/fga.h>
#include <DISCRETA/lb.h>
#include <DISCRETA/DESIGNS/solid.h>


#if TEXDOCU
INT SOLID_OB::init()
#endif
{
	INT erg = OK;
	
	erg += m_il(14);
	c_obj_k(SOLID_KIND);
	s_group_generators()->m_il(0);
	s_nb_V()->m_i(0);
	s_nb_E()->m_i(0);
	s_nb_F()->m_i(0);
	s_placement()->m_il(0);
	s_v1()->m_il(0);
	s_v2()->m_il(0);
	s_f1()->m_il(0);
	s_f2()->m_il(0);
	s_nb_e()->m_il(0);
	s_edge()->m_il(0);
	s_neighbour_faces()->m_il(0);
	s_f_vertex_labels()->m_i(FALSE);
	s_vertex_labels()->m_il(0);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::init_V(INT nb_V)
#endif
{
	s_nb_V()->m_i(nb_V);
	s_placement()->m_il(3);
	s_x()->m_il_n(nb_V);
	s_y()->m_il_n(nb_V);
	s_z()->m_il_n(nb_V);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::init_E(INT nb_E)
#endif
{
	s_nb_E()->m_i(nb_E);
	s_v1()->m_il_n(nb_E);
	s_v2()->m_il_n(nb_E);
	s_f1()->m_il_n(nb_E);
	s_f2()->m_il_n(nb_E);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::init_F(INT nb_F)
#endif
{
	s_nb_F()->m_i(nb_F);
	s_nb_e()->m_il_n(nb_F);
	s_edge()->m_il(nb_F);
	s_neighbour_faces()->m_il(nb_F);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::standard_vertex_labels(INT f_start_with_zero)
#endif
{
	INT i, l;
	BYTE str[1000];
	
	l = s_nb_V_i();
	s_vertex_labels()->m_il(l);
	for (i = 0; i < l; i++) {
		if (f_start_with_zero)
			sprintf(str, "%ld", i);
		else
			sprintf(str, "%ld", i + 1);
		s_vertex_labels_i(i)->init(str);
		}
	s_f_vertex_labels()->m_i(TRUE);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::sprint(BYTE *s)
#endif
{
	INT i, l;
	
	sprintf(s + strlen(s), "SOLID\n");
	sprintf(s + strlen(s), "nb_V = %ld\n", s_nb_V_i());
	sprintf(s + strlen(s), "nb_E = %ld\n", s_nb_E_i());
	sprintf(s + strlen(s), "nb_F = %ld\n", s_nb_V_i());
	l = s_group_generators()->s_li();
	sprintf(s + strlen(s), "nb of automorphisms = %ld\n", l);
	for (i = 0; i < l; i++) 
	{
		s_group_generators()->s_i(i)->sprint(s);
		sprintf(s + strlen(s), "\n");
	}
	if (s_f_vertex_labels_i()) {
		sprintf(s + strlen(s), "with vertex labels\n");	
		}
	return OK;
}

#include "archimed.C"

INT init_solid_from_archimed(SOLID_OP S, ARCHIMED *X)
{
	INT i, j, l, nb_V, nb_E, nb_F;

	nb_V = X->nb_V;
	nb_E = X->nb_E;
	nb_F = X->nb_F;
	
	S->init();

	S->s_f_vertex_labels()->m_i(FALSE);
	S->s_placement()->m_il(3);
	S->s_x()->m_il(nb_V);
	S->s_y()->m_il(nb_V);
	S->s_z()->m_il(nb_V);
	for (i = 0; i < nb_V; i++)
	{
		S->s_x_i(i)->m_i(X->Px[i]);
		S->s_y_i(i)->m_i(X->Py[i]);
		S->s_z_i(i)->m_i(X->Pz[i]);
	}

	S->s_nb_V()->m_i(nb_V);
	
	S->s_v1()->m_il(nb_E);
	S->s_v2()->m_il(nb_E);
	S->s_f1()->m_il(nb_E);
	S->s_f2()->m_il(nb_E);
	for (i = 0; i < nb_E; i++) {
		S->s_v1()->m_ii(i, X->e1[i]);
		S->s_v2()->m_ii(i, X->e2[i]);
		S->s_f1()->m_ii(i, X->f1[i]);
		S->s_f2()->m_ii(i, X->f2[i]);
		}
	S->s_nb_E()->m_i(nb_E);
	
	S->s_nb_e()->m_il(nb_F);
	S->s_edge()->m_il(nb_F);
	S->s_neighbour_faces()->m_il(nb_F);
	for (i = 0; i < nb_F; i++) {
		VECTOR_OB E, N;
		
		l = X->nb_e[i];
		S->s_nb_e()->m_ii(i, l);
		E.m_il(l);
		N.m_il(l);
		for (j = 0; j < l; j++) {
			E.m_ii(j, X->edge[i][j]);
			N.m_ii(j, X->neighbour[i][j]);
			}
		E.swap(S->s_edge_i(i));
		N.swap(S->s_neighbour_faces_i(i));
		}
	S->s_nb_F()->m_i(nb_F);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::cubus_simus(INT s)
// s should be 1000
#endif
{
	ARCHIMED *X;
	VECTOR_OB V;
	INT f_v = FALSE;
	
	X = (ARCHIMED *) my_malloc(sizeof(ARCHIMED), "ARCHIMED");
	cubus_simus1(s, "Cubus_simus.graph", X, &V, f_v);
	init_solid_from_archimed(this, X);
	V.copy(s_group_generators());
	
	standard_vertex_labels(FALSE /* f_start_with_zero */);

	my_free(X);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::dode_simum(INT s)
// s should be 1000
#endif
{
	ARCHIMED *X;
	VECTOR_OB V;
	INT f_v = FALSE;
	
	X = (ARCHIMED *) my_malloc(sizeof(ARCHIMED), "ARCHIMED");
	dode_simum1(s, "dode_simum.graph", X, &V, f_v);
	init_solid_from_archimed(this, X);
	V.copy(s_group_generators());
	
	standard_vertex_labels(FALSE /* f_start_with_zero */);

	my_free(X);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::snub_cube(INT s)
// s should be 1000
#endif
{
	ARCHIMED *X;
	VECTOR_OB V;
	INT f_v = FALSE;
	
	X = (ARCHIMED *) my_malloc(sizeof(ARCHIMED), "ARCHIMED");
	archimed_snub_cube(s, "Snub_cube.graph", X, &V, f_v);
	init_solid_from_archimed(this, X);
	V.copy(s_group_generators());
	
	standard_vertex_labels(FALSE /* f_start_with_zero */);

	my_free(X);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::russian_snub_cube(INT s)
// s should be 1000
#endif
{
	ARCHIMED *X;
	VECTOR_OB V;
	INT f_v = FALSE;
	
	X = (ARCHIMED *) my_malloc(sizeof(ARCHIMED), "ARCHIMED");
	archimed_russian_snub_cube(s, "Snub_cube_russian.graph", X, &V, f_v);
	init_solid_from_archimed(this, X);
	V.copy(s_group_generators());
	
	standard_vertex_labels(FALSE /* f_start_with_zero */);

	my_free(X);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::tetrahedron(INT r)
#endif
{
	double phi, h, c;
	INT i;
	PERMUTATION_OB P;
	
	init();
	s_group_generators()->m_il(2);
	P.m_il(4);
	P.one();
	P.Add3Cycle(1, 2, 3);
	P.copy((PERMUTATION_OP) s_group_generators_i(0));
	P.m_il(4);
	P.one();
	P.Add3Cycle(1, 2, 4);
	P.copy((PERMUTATION_OP) s_group_generators_i(1));
	s_f_vertex_labels()->m_i(FALSE);
	phi = 120.;
	h = 4. * r * 0.33333;
	c = sqrt(8. / 9.) * r;
	
	s_placement()->m_il(3);
	s_x()->m_il(4);
	s_y()->m_il(4);
	s_z()->m_il(4);
	for (i = 0; i < 3; i++)
	{
		On_circle_int(s_x(), s_y(), i, (INT) ((double) i * phi), (INT)c);
		s_z_i(i)->m_i((INT)(-0.333333 * r));
	}
	s_x_i(3)->m_i(0);
	s_y_i(3)->m_i(0);
	s_z_i(3)->m_i(r);
	
	s_nb_V()->m_i(4);
	add_edge(0, 1, 2, 3);	
	add_edge(1, 2, 0, 3);	
	add_edge(2, 0, 1, 3);	
	add_edge(0, 3, 1, 2);	
	add_edge(1, 3, 2, 0);	
	add_edge(2, 3, 0, 1);	

	add_face3(1, 5, 4, 3, 1, 2);
	add_face3(2, 3, 5, 3, 2, 0);
	add_face3(0, 4, 3, 3, 0, 1);
	add_face3(0, 1, 2, 2, 0, 1);
	
	determine_neighbours();
	standard_vertex_labels(FALSE /* f_start_with_zero */);
	return OK;
}

#if 0
#if TEXDOCU
INT SOLID_OB::cube(INT r)
#endif
{
	PERMUTATION_OB p;
	SOLID_OB A, B;
	INT l;
	
	A.octahedron(r);
	A.dual(&B);
	p.m_il(8);
	p.Add4Cycle(1, 2, 3, 4);
	p.Add4Cycle(5, 6, 7, 8);
	l = B.s_group_generators()->s_li();
	B.s_group_generators()->inc();
	p.swap(B.s_group_generators()->s_i(l));
	p.m_il(8);
	p.m_ii(0, 8);
	p.m_ii(1, 7);
	p.m_ii(2, 3);
	p.m_ii(3, 4);
	p.m_ii(4, 6);
	p.m_ii(5, 5);
	p.m_ii(6, 1);
	p.m_ii(7, 2);
	B.relabel_points(this, &p, TRUE /* f_relabel_vertex_labels */);
	standard_vertex_labels(FALSE /* f_start_with_zero */);
	return OK;
}
#endif

#if TEXDOCU
INT SOLID_OB::cube(INT r)
#endif
{
	PERMUTATION_OB p;
	INT r2 = r >> 1;
	
	init();
	s_nb_V()->m_i(8);
	s_placement()->m_il(3);
	s_x()->m_il(8);
	s_y()->m_il(8);
	s_z()->m_il(8);
	s_x_i(0)->m_i(-r2);
	s_y_i(0)->m_i(r2);
	s_z_i(0)->m_i(-r2);
	s_x_i(1)->m_i(r2);
	s_y_i(1)->m_i(r2);
	s_z_i(1)->m_i(-r2);
	s_x_i(2)->m_i(-r2);
	s_y_i(2)->m_i(-r2);
	s_z_i(2)->m_i(-r2);
	s_x_i(3)->m_i(r2);
	s_y_i(3)->m_i(-r2);
	s_z_i(3)->m_i(-r2);
	s_x_i(4)->m_i(-r2);
	s_y_i(4)->m_i(r2);
	s_z_i(4)->m_i(r2);
	s_x_i(5)->m_i(r2);
	s_y_i(5)->m_i(r2);
	s_z_i(5)->m_i(r2);
	s_x_i(6)->m_i(-r2);
	s_y_i(6)->m_i(-r2);
	s_z_i(6)->m_i(r2);
	s_x_i(7)->m_i(r2);
	s_y_i(7)->m_i(-r2);
	s_z_i(7)->m_i(r2);

	add_face4(0, 1, 5, 4);
	add_face4(1, 3, 7, 5);
	add_face4(3, 2, 6, 7);
	add_face4(2, 0, 4, 6);
	add_face4(6, 7, 5, 4);
	add_face4(2, 3, 1, 0);

	determine_neighbours();
	standard_vertex_labels(FALSE /* f_start_with_zero */);

	PERMUTATION_OB P;

	s_group_generators()->m_il(2);
	P.m_il(8);
	P.Add4Cycle(1, 2, 4, 3);
	P.Add4Cycle(5, 6, 8, 7);
	P.swap((PERMUTATION_OP) s_group_generators_i(0));
	P.m_il(8);
	P.Add4Cycle(1, 3, 7, 5);
	P.Add4Cycle(2, 4, 8, 6);
	P.swap((PERMUTATION_OP) s_group_generators_i(1));


	return OK;
}

#if TEXDOCU
INT SOLID_OB::octahedron(INT r)
#endif
{
	PERMUTATION_OB P;
	init();
	
	s_group_generators()->m_il(2);
	P.m_il(6);
	P.Add3Cycle(1, 3, 6);
	P.Add3Cycle(2, 4, 5);
	P.swap((PERMUTATION_OP) s_group_generators_i(0));
	P.m_il(6);
	P.Add3Cycle(1, 3, 5);
	P.Add3Cycle(2, 4, 6);
	P.swap((PERMUTATION_OP) s_group_generators_i(1));
	s_f_vertex_labels()->m_i(FALSE);
	s_nb_V()->m_i(6);
	s_placement()->m_il(3);
	s_x()->m_il(6);
	s_y()->m_il(6);
	s_z()->m_il(6);
	s_x_i(0)->m_i(r);
	s_y_i(0)->m_i(0);
	s_z_i(0)->m_i(0);
	s_x_i(1)->m_i(-r);
	s_y_i(1)->m_i(0);
	s_z_i(1)->m_i(0);
	s_x_i(2)->m_i(0);
	s_y_i(2)->m_i(r);
	s_z_i(2)->m_i(0);
	s_x_i(3)->m_i(0);
	s_y_i(3)->m_i(-r);
	s_z_i(3)->m_i(0);
	s_x_i(4)->m_i(0);
	s_y_i(4)->m_i(0);
	s_z_i(4)->m_i(r);
	s_x_i(5)->m_i(0);
	s_y_i(5)->m_i(0);
	s_z_i(5)->m_i(-r);
	add_edge(0, 4, 0, 4);	
	add_edge(4, 1, 1, 5);	
	add_edge(1, 5, 2, 6);	
	add_edge(5, 0, 3, 7);	
	add_edge(0, 2, 0, 3);	
	add_edge(4, 2, 0, 1);	
	add_edge(1, 2, 1, 2);	
	add_edge(5, 2, 2, 3);	
	add_edge(0, 3, 7, 4);	
	add_edge(4, 3, 4, 5);	
	add_edge(1, 3, 5, 6);	
	add_edge(5, 3, 6, 7);	

	add_face3(0, 5, 4, 4, 1, 3);
	add_face3(1, 6, 5, 5, 2, 0);
	add_face3(2, 7, 6, 6, 3, 1);
	add_face3(3, 4, 7, 7, 0, 2);
	add_face3(0, 8, 9, 0, 7, 5);
	add_face3(1, 9, 10, 1, 4, 6);
	add_face3(2, 10, 11, 2, 5, 7);
	add_face3(3, 11, 8, 3, 6, 4);

	determine_neighbours();
	standard_vertex_labels(FALSE /* f_start_with_zero */);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::dodecahedron(INT r)
#endif
{
	PERMUTATION_OB P;
	
	init();
	double phi, phi_2, sin_phi_2, R, RR, s, d, dr, h, hh, H, a, hH;
	INT i;
	
	s_group_generators()->m_il(2);
	P.m_il(20);
	P.Add5Cycle(1,2,3,4,5);
	P.Add5Cycle(6,7,8,9,10);
	P.Add5Cycle(11,12,13,14,15);
	P.Add5Cycle(16,17,18,19,20);
	P.copy((PERMUTATION_OP)s_group_generators_i(0));
	P.m_il(20);
	P.Add5Cycle(1,2,7,19,6);
	P.Add5Cycle(3,20,14,18,5);
	P.Add5Cycle(4,8,15,13,10);
	P.Add5Cycle(9,16,11,12,17);
	P.copy((PERMUTATION_OP)s_group_generators_i(1));
	s_f_vertex_labels()->m_i(FALSE);
	s_nb_V()->m_i(20);
	s_placement()->m_il(3);
	s_x()->m_il(20);
	s_y()->m_il(20);
	s_z()->m_il(20);
	phi = 72.;
	phi_2 = 36.;
	sin_phi_2 = sin_grad(phi * .5);
	s = 2. * (double) r * sin_phi_2;
	d = 2. * (double) r * sin_grad(phi);
	R = 0.5 * d / sin_phi_2;
	dr = R - r;
	hh = s * s - dr * dr;
	h = sqrt(hh);
	RR = R * R;
	H = .5 * (RR - hh - r * r) / h;
	a = sqrt(RR + H * H);
	hH = h + H;
	for (i = 0; i < 5; i++) 
	{
		On_circle_int(s_x(), s_y(), i, (INT)((double) i * phi), r);
		s_z_i(i)->m_i(- (INT) hH);
	}
	for (i = 0; i < 5; i++) 
	{
		On_circle_int(s_x(), s_y(), 5 + i, (INT)((double) i * phi), (INT) R);
		s_z_i(5 + i)->m_i(- (INT) H);
	}
	for (i = 0; i < 10; i++) 
	{
		s_x_i(10 + i)->m_i(-s_x_i(i)->s_i());
		s_y_i(10 + i)->m_i(-s_y_i(i)->s_i());
		s_z_i(10 + i)->m_i(-s_z_i(i)->s_i());
	}

	add_face5(0, 1, 2, 3, 4);
	add_face5(0, 1, 6, 18, 5);
	add_face5(1, 2, 7, 19, 6);
	add_face5(2, 3, 8, 15, 7);
	add_face5(3, 4, 9, 16, 8);
	add_face5(4, 0, 5, 17, 9);
	add_face5(5, 18, 13, 12, 17);
	add_face5(6, 19, 14, 13, 18);
	add_face5(7, 15, 10, 14, 19);
	add_face5(8, 16, 11, 10, 15);
	add_face5(9, 17, 12, 11, 16);
	add_face5(10, 11, 12, 13, 14);
	
	determine_neighbours();
	standard_vertex_labels(FALSE /* f_start_with_zero */);
	return OK;
	
}

#if TEXDOCU
INT SOLID_OB::icosahedron(INT r)
#endif
{
	SOLID_OB A;
	
	A.dodecahedron(r);
	A.standard_vertex_labels(FALSE /* f_start_with_zero */);
	A.dual(this);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::dual(SOLID_OP A)
#endif
{
	INT adj[1000], nb_adj;
	INT v, i, j, f, e, e1, e2, v1, v2, v3, v4, V1, V2, V3, V4, l;
	SOLID_OB B;
	PERMUTATION_OB P;
	PERMUTATION_OP Q;	
	INT nb_F = s_nb_F_i();
	INT nb_E = s_nb_E_i();
	INT nb_V = s_nb_V_i();
	
	printf("dual()\n");
	B.init();
	B.init_V(nb_F);
	B.init_E(nb_E);
	B.init_F(nb_V);
	l = s_group_generators()->s_li();
	B.s_group_generators()->m_il(l);
	for (i = 0; i < l; i++) 
	{
		printf("i=%ld\n", i);
		Q = (PERMUTATION_OP) s_group_generators_i(i);
		Q->println();
		P.m_il(nb_F);
		for (j = 0; j < nb_F; j++) 
		{
			e1 = s_edge_iji(j, 0);
			e2 = s_edge_iji(j, 1);
			v1 = s_v1_ii(e1);
			v2 = s_v2_ii(e1);
			v3 = s_v1_ii(e2);
			v4 = s_v2_ii(e2);
			V1 = Q->s_ii(v1) - 1;
			V2 = Q->s_ii(v2) - 1;
			V3 = Q->s_ii(v3) - 1;
			V4 = Q->s_ii(v4) - 1;
			e1 = find_edge(V1, V2);
			e2 = find_edge(V3, V4);
			if (e1 == -1) 
			{
				printf("v1=%ld\n", v1);
				printf("v2=%ld\n", v2);
				printf("v3=%ld\n", v3);
				printf("v4=%ld\n", v4);
				printf("V1=%ld\n", V1);
				printf("V2=%ld\n", V2);
				printf("V3=%ld\n", V3);
				printf("V4=%ld\n", V4);
				return error("dual: edge e1 not found");
			}
			if (e2 == -1) 
			{
				printf("v1=%ld\n", v1);
				printf("v2=%ld\n", v2);
				printf("v3=%ld\n", v3);
				printf("v4=%ld\n", v4);
				printf("V1=%ld\n", V1);
				printf("V2=%ld\n", V2);
				printf("V3=%ld\n", V3);
				printf("V4=%ld\n", V4);
				return error("dual: edge e2 not found");
			}
			f = find_face_2(e1, e2);
			P.m_ii(j, f + 1);
		}
		P.copy((PERMUTATION_OP) B.s_group_generators_i(i));
	}
	
	for (i = 0; i < nb_F; i++) 
	{
		center(i, B.s_x(), B.s_y(), B.s_z());
	}
	for (i = 0; i < nb_E; i++) 
	{
		B.s_v1_i(i)->m_i(s_f1_ii(i));
		B.s_v2_i(i)->m_i(s_f2_ii(i));
		B.s_f1_i(i)->m_i(s_v1_ii(i));
		B.s_f2_i(i)->m_i(s_v2_ii(i));
	}
	
	for (v = 0; v < nb_V; v++) 
	{
	
		adjacency_list(v, adj, &nb_adj);
		B.s_nb_e_i(v)->m_i(nb_adj);
		((VECTOR_OP)B.s_edge_i(v))->m_il_n(nb_adj);
		((VECTOR_OP)B.s_neighbour_faces_i(v))->m_il_n(nb_adj);
		
		for (i = 0; i < nb_adj; i++) 
		{
			e = adj[i];
			v1 = s_v1_ii(e);
			v2 = s_v2_ii(e);
			B.s_edge_ij(v, i)->m_i(e);
			if (v1 == v) 
			{
				B.s_neighbour_faces_ij(v, i)->m_i(v2);
			}
			else if (v2 == v) 
			{
				B.s_neighbour_faces_ij(v, i)->m_i(v1);
			}
			else
				return error("edge does not contain vertex v");
		}
	}
	B.swap(A);
	if (s_f_vertex_labels_i())
		A->standard_vertex_labels(FALSE /* f_start_with_zero */);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::center(INT f, VECTOR_OP Px, VECTOR_OP Py, VECTOR_OP Pz)
#endif
{
	INT i, nb_e, e, e1, e2;
	INT x = 0, y = 0, z = 0;
	
	nb_e = s_nb_e_ii(f);
	for (i = 0; i < nb_e; i++) 
	{
		e = s_edge_iji(f, i);
		e1 = s_v1_ii(e);
		e2 = s_v2_ii(e);
		x += s_x_ii(e1);
		y += s_y_ii(e1);
		z += s_z_ii(e1);
		x += s_x_ii(e2);
		y += s_y_ii(e2);
		z += s_z_ii(e2);
	}
	nb_e <<= 1;
	Px->m_ii(f, (INT)((double)x / (double)nb_e));
	Py->m_ii(f, (INT)((double)y / (double)nb_e));
	Pz->m_ii(f, (INT)((double)z / (double)nb_e));
	return OK;
}

#if TEXDOCU
INT SOLID_OB::adjacency_list(INT vertex, INT *adj, INT *nb_adj)
#endif
{
	INT i, j, l, ll, a, b, f_found;
	VECTOR_OB Adj;
	INT nb_E = s_nb_E_i();
	
	Adj.m_il(0);
	ll = 0;
	for (i = 0; i < nb_E; i++) {
		if (s_v1_ii(i) == vertex || s_v2_ii(i) == vertex) {
			Adj.inc();
			Adj.m_ii(ll++, i);
		}
	}
	l = 0;
	adj[l++] = a = Adj.s_ii(--ll);
	Adj.dec();
	while (ll) {
		f_found = FALSE;
		for (i = 0; i < ll; i++) {
			b = Adj.s_ii(i);
			if (find_face_by_two_edges(a, b) != -1) {
				adj[l++] = b;
				for (j = i + 1; j < ll; j++) {
					Adj.m_ii(j - 1, Adj.s_ii(j));
					}
				Adj.dec();
				ll--;
				a = b;
				f_found = TRUE;
				break;
				}
			}
		if (!f_found) {
			return error("found no adjacent edge");
			}
		}
	*nb_adj = l;
	return OK;
}

#if TEXDOCU
INT SOLID_OB::add_edge(INT v1, INT v2, INT f1, INT f2)
#endif
{
	INT nb_E = s_nb_E_i();
	
	s_v1()->inc();
	s_v2()->inc();
	s_f1()->inc();
	s_f2()->inc();
	s_v1()->m_ii(nb_E, v1);
	s_v2()->m_ii(nb_E, v2);
	s_f1()->m_ii(nb_E, f1);
	s_f2()->m_ii(nb_E, f2);
	s_nb_E()->inc();
	return OK;
}

#if TEXDOCU
INT SOLID_OB::add_face3(INT e1, INT e2, INT e3, INT n1, INT n2, INT n3)
#endif
{
	INT nb_F = s_nb_F_i();
	
	s_nb_e()->inc();
	s_edge()->inc();
	s_neighbour_faces()->inc();
	s_nb_e()->m_ii(nb_F, 3);
	s_edge_i(nb_F)->m_il(3);
	s_edge_ij(nb_F, 0)->m_i(e1);
	s_edge_ij(nb_F, 1)->m_i(e2);
	s_edge_ij(nb_F, 2)->m_i(e3);
	s_neighbour_faces_i(nb_F)->m_il(3);
	s_neighbour_faces_ij(nb_F, 0)->m_i(n1);
	s_neighbour_faces_ij(nb_F, 1)->m_i(n2);
	s_neighbour_faces_ij(nb_F, 2)->m_i(n3);
	s_nb_F()->inc();
	return OK;
}


#if TEXDOCU
INT SOLID_OB::add_face4(INT i1, INT i2, INT i3, INT i4)
#endif
{
	VECTOR_OB V;
	
	V.m_il(4);
	V.m_ii(0, i1);
	V.m_ii(1, i2);
	V.m_ii(2, i3);
	V.m_ii(3, i4);
	add_face_n(&V);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::add_face5(INT i1, INT i2, INT i3, INT i4, INT i5)
#endif
{
	INT e;
	INT nb_F = s_nb_F_i();
	INT f_v = FALSE;
	
	s_nb_e()->inc();
	s_edge()->inc();
	s_neighbour_faces()->inc();
	s_nb_e()->m_ii(nb_F, 5);
	s_edge_i(nb_F)->m_il(5);
	s_neighbour_faces_i(nb_F)->m_il(5);
	
	e = find_and_add_edge(i1, i2, f_v);
	s_edge_ij(nb_F, 0)->m_i(e);
	
	e = find_and_add_edge(i2, i3, f_v);
	s_edge_ij(nb_F, 1)->m_i(e);
	
	e = find_and_add_edge(i3, i4, f_v);
	s_edge_ij(nb_F, 2)->m_i(e);
	
	e = find_and_add_edge(i4, i5, f_v);
	s_edge_ij(nb_F, 3)->m_i(e);
	
	e = find_and_add_edge(i5, i1, f_v);
	s_edge_ij(nb_F, 4)->m_i(e);
		
	s_nb_F()->inc();
	return OK;
}

#if TEXDOCU
INT SOLID_OB::add_face_n(VECTOR_OP vertices)
#endif
{
	INT e, i, l;
	INT nb_F = s_nb_F_i();
	INT f_v = FALSE;
	
	l = vertices->s_li();
	s_nb_e()->inc();
	s_edge()->inc();
	s_neighbour_faces()->inc();
	s_nb_e()->m_ii(nb_F, l);
	s_edge_i(nb_F)->m_il(l);
	s_neighbour_faces_i(nb_F)->m_il(l);
	
	for (i = 0; i < l - 1; i++) {
		e = find_and_add_edge(vertices->s_ii(i), vertices->s_ii(i + 1), f_v);
		s_edge_ij(nb_F, i)->m_i(e);
		}
	e = find_and_add_edge(vertices->s_ii(l - 1), vertices->s_ii(0), f_v);
	s_edge_ij(nb_F, l - 1)->m_i(e);
	
	s_nb_F()->inc();
	return OK;
}


#if TEXDOCU
INT SOLID_OB::find_and_add_edge(INT i1, INT i2, INT f_v)
#endif
{
	INT e;
	
	e = find_edge(i1, i2);
	if (e == -1) 
	{
		if (f_v) 
		{
			printf("adding edge %ld %ld ", i1, i2); fflush(stdout);
		}
		e = add_edge(i1, i2);
		if (f_v) 
		{
			printf("- edge no %ld\n", e); fflush(stdout);
		}
		
	}
	return e;
}

#if TEXDOCU
INT SOLID_OB::find_edge(INT v1, INT v2)
#endif
{
	INT i, e1, e2;
	INT nb_E = s_nb_E_i();
	
	for (i = 0; i < nb_E; i++) 
	{
		e1 = s_v1_i(i)->s_i();
		e2 = s_v2_i(i)->s_i();
		if ((e1 == v1 && e2 == v2) || (e1 == v2 && e2 == v1)) 
		{
			return i;
		}
	}
	return -1;
}

#if TEXDOCU
INT SOLID_OB::add_edge(INT v1, INT v2)
#endif
{
	INT i;
	INT nb_E = s_nb_E_i();

	s_v1()->inc();
	s_v2()->inc();	
	s_f1()->inc();
	s_f2()->inc();	
	s_v1()->m_ii(nb_E, v1);
	s_v2()->m_ii(nb_E, v2);
	i = nb_E;
	s_nb_E()->inc();
	return i;
}

#if TEXDOCU
INT SOLID_OB::find_face_by_two_edges(INT e1, INT e2)
#endif
{
	INT f1, f2, f3, f4;
	
	find_faces_at_edge(e1, &f1, &f2);
	find_faces_at_edge(e2, &f3, &f4);
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

#if TEXDOCU
INT SOLID_OB::find_faces_at_edge(INT e, INT *f1, INT *f2)
#endif
{
	INT nb_F, i, j, l, n = 0;
	
	*f1 = -1;
	*f2 = -1;
	nb_F = s_nb_F_i();
	for (i = 0; i < nb_F; i++) {
		l = s_nb_e_ii(i);
		for (j = 0; j < l; j++) {
			if (s_edge_iji(i, j) == e) {
				if (n == 0) {
					*f1 = i;
					n++;
					}
				else if (n == 1) {
					*f2 = i;
					n++;
					}
				else {
					return error("SOLID_OB::find_faces_at_edge(): too many faces for this edge");
					}
				}
			}
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::find_face(INT e, INT *f1, INT *j1, INT *f2, INT *j2)
#endif
{
	INT i, j, l, ff1 = -1, ff2 = -1;
	INT nb_F = s_nb_F_i();
	
	for (i = 0; i < nb_F; i++) 
	{
		l = s_nb_e_ii(i);
		for (j = 0; j < l; j++) 
		{
			if (s_edge_iji(i, j) == e) 
			{
				if (ff1 != -1) 
				{
					ff2 = i;
					*j2 = j;
				}
				else 
				{
					ff1 = i;
					*j1 = j;
				}
			}
		}
	}
	if (ff1 == -1 || ff2 == -1) 
	{
		printf("SOLID::find_face() face not found\n");
		printf("edge e = %ld ff1=%ld ff2=%ld\n", e, ff1, ff2);
		printf("the faces and their edges are:\n");
		for (i = 0; i < nb_F; i++) 
		{
			printf("face %ld: ", i);
			l = s_nb_e_ii(i);
			for (j = 0; j < l; j++) 
			{
				printf("%ld ", s_edge_iji(i, j));
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

#if TEXDOCU
INT SOLID_OB::find_face_2(INT e1, INT e2)
#endif
{
	INT f1, f2, f3, f4, j1, j2, j3, j4;
	
	find_face(e1, &f1, &j1, &f2, &j2);
	find_face(e2, &f3, &j3, &f4, &j4);
	if (f1 == f3 || f1 == f4) 
		return f1;
	if (f2 == f3 || f2 == f4) 
		return f2;
	return error("find_face_2: face not found");
}

#if TEXDOCU
INT SOLID_OB::determine_neighbours()
#endif
{
	INT nb_E = s_nb_E_i();
	INT e, f1, j1, f2, j2;
	
	for (e = 0; e < nb_E; e++) {
		find_face(e, &f1, &j1, &f2, &j2);
		s_f1_i(e)->m_i(f1);
		s_f2_i(e)->m_i(f2);
		s_neighbour_faces_ij(f1, j1)->m_i(f2);
		s_neighbour_faces_ij(f2, j2)->m_i(f1);
		}
	return OK;
}


#if TEXDOCU
INT SOLID_OB::write_graphfile(BYTE *fname)
#endif
{
	BYTE name[1000], str[256], *p;
	FILE *fp, *fp1;
	INT i, l;
	INT nb_V = s_nb_V()->s_i();
	INT nb_E = s_nb_E()->s_i();
	
	printf("SOLID::write_graphfile() fname=%s nb_V=%ld nb_E=%ld\n", 
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
	fprintf(fp, "<GRAPH NAME=\"%s\" NUMBER_OF_VERTICES=%ld NUMBER_OF_EDGES=%ld>\n", 
		name, s_nb_V_i(), s_nb_E_i());
	fprintf(fp, "<!-- created by DISCRETA, %s -->\n", str);
	fprintf(fp, "<COORDS3D_INT>\n");
	for (i = 0; i < nb_V; i++) 
	{
		fprintf(fp, "%ld %ld %ld\n", 
			s_x_ii(i), 
			s_y_ii(i), 
			s_z_ii(i)); 
	}
	fprintf(fp, "</COORDS3D_INT>\n");
	fprintf(fp, "<EDGELIST>\n");
	for (i = 0; i < nb_E; i++) 
	{
		fprintf(fp, "%ld %ld\n", s_v1_ii(i), s_v2_ii(i));
	}
	fprintf(fp, "</EDGELIST>\n");
	if (s_f_vertex_labels_i()) {
		fprintf(fp, "<VERTEXLABELS>\n");
		for (i = 0; i < nb_V; i++) 
		{
			fprintf(fp, "\"%s\"\n", s_vertex_labels_is(i)); 
		}
		fprintf(fp, "</VERTEXLABELS>\n");
	}

	fprintf(fp, "</GRAPH>\n\n");
	fclose(fp);
	return OK;
}

#if 0
#if TEXDOCU
INT SOLID_OB::write_graphfile(BYTE *fname)
#endif
{
	FILE *fp;
	INT i;
	INT nb_V = s_nb_V()->s_i();
	INT nb_E = s_nb_E()->s_i();
	
	printf("SOLID::write_graphfile() fname=%s nb_V=%ld nb_E=%ld\n", 
		fname, nb_V, nb_E); fflush(stdout);
	fp = fopen(fname, "w");
	fprintf(fp, "NUMBER_OF_VERTICES %ld\n", s_nb_V_i());
	fprintf(fp, "NUMBER_OF_EDGES %ld\n", s_nb_E_i());
	fprintf(fp, "VERTICES\n");
	for (i = 0; i < nb_V; i++) 
	{
		fprintf(fp, "%ld %ld %ld\n", 
			s_x_ii(i), 
			s_y_ii(i), 
			s_z_ii(i)); 
	}
	fprintf(fp, "EDGES\n");
	for (i = 0; i < nb_E; i++) 
	{
		fprintf(fp, "%ld %ld\n", s_v1_ii(i), s_v2_ii(i));
	}
	if (s_f_vertex_labels_i()) {
		fprintf(fp, "VERTEX-LABELS\n");
		for (i = 0; i < nb_V; i++) 
		{
			fprintf(fp, "%s\n", s_vertex_labels_is(i)); 
		}
	}

	fprintf(fp, "\n");
	fclose(fp);
	return OK;
}
#endif

#if TEXDOCU
INT SOLID_OB::archimed_scale(double f)
#endif
{
	INT i, a;
	
	for (i = 0; i < s_nb_V_i(); i++) 
	{
		a = s_x_ii(i);
		a = (INT)((double) a * f);
		s_x_i(i)->m_i(a);
		a = s_y_ii(i);
		a = (INT)((double) a * f);
		s_y_i(i)->m_i(a);
		a = s_z_ii(i);
		a = (INT)((double) a * f);
		s_z_i(i)->m_i(a);
	}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::cut_vertices(double r, SOLID_OP A)
#endif
{
	INT adj[1000], adj1[1000], nb_adj, nb_adj1;
	INT nb_V, nb_E, nb_F, nb_new_V;
	INT v, v1;
	INT ei, ej, e, ee, e1, e2, f, i, j, l, first_new_vertex, k;
	INTEGER_OB Px, Py, Pz;
	SOLID_OB B;
	PERMUTATION_OB P, Pind;
	PERMUTATION_OP Q;	
	VECTOR_OB Pindv;
	VECTOR_OB first_new_vertex_number;
	
	
	printf("cut_vertices()\n");
	nb_V = s_nb_V_i();
	nb_E = s_nb_E_i();
	nb_F = s_nb_F_i();
	copy(&B);
	l = s_group_generators()->s_li();

	first_new_vertex_number.m_il(nb_V);
	first_new_vertex = nb_V;
	for (v = 0; v < s_nb_V_i(); v++) {
		first_new_vertex_number.m_ii(v, first_new_vertex);
		adjacency_list(v, adj, &nb_adj);
		first_new_vertex += nb_adj;
		} 
	
	for (v = 0; v < s_nb_V_i(); v++)
	{
	
		first_new_vertex = nb_V;
		
		adjacency_list(v, adj, &nb_adj);
		
		// shortens the edges from one side (near v)
		for (ei = 0; ei < nb_adj; ei++) 
		{
			e = adj[ei];
			e1 = s_v1_ii(e);
			B.s_x()->inc();
			B.s_y()->inc();
			B.s_z()->inc();
			B.s_nb_V()->inc();
			if (e1 == v)
				Ratio(e, r, B.s_x_i(nb_V), B.s_y_i(nb_V), B.s_z_i(nb_V));
			else
				Ratio(e, 1. - r, B.s_x_i(nb_V), B.s_y_i(nb_V), B.s_z_i(nb_V));
			nb_V++;
		}		
		
		// define a new face:
		B.s_nb_e()->inc();
		B.s_edge()->inc();
		B.s_neighbour_faces()->inc();
		B.s_nb_e_i(nb_F)->m_i(0);
		B.s_edge_i(nb_F)->m_il(nb_adj);
		B.s_neighbour_faces_i(nb_F)->m_il(nb_adj);
		
		// all the edges of the new face:
		for (ei = 0; ei < nb_adj; ei++) 
		{
			e = adj[ei];
			for (ej = ei + 1; ej < nb_adj; ej++) 
			{
				ee = adj[ej];
				if (find_common_face(e, ee, &f)) 
				{
					
					// new edge:
					B.s_v1()->inc();
					B.s_v1_i(nb_E)->m_i(first_new_vertex + ei);
					B.s_v2()->inc();
					B.s_v2_i(nb_E)->m_i(first_new_vertex + ej);
					B.s_f1()->inc();
					B.s_f1_i(nb_E)->m_i(f);
					B.s_f2()->inc();
					B.s_f2_i(nb_E)->m_i(nb_F);
					B.s_nb_E()->inc();
					
					// new face:
					
					k = B.s_nb_e_ii(nb_F);
					B.s_edge_ij(nb_F, k)->m_i(nb_E);
					B.s_neighbour_faces_ij(nb_F, k)->m_i(f);
					B.s_nb_e_i(nb_F)->inc();
					
					// old face:
					k = B.s_nb_e_ii(f);
					B.s_edge_i(f)->inc();
					B.s_neighbour_faces_i(f)->inc();
					B.s_edge_ij(f, k)->m_i(nb_E);
					B.s_neighbour_faces_ij(f, k)->m_i(nb_F);
					B.s_nb_e_i(f)->inc();
					
					nb_E++;
				}
			}
		}
		
		
		// face:
		B.s_nb_F()->inc();
		nb_F++;
	} // next v

	// induce the group:
	printf("inducing the group, number of generators = %ld\n", l);
	printf("nb_V = %ld s_nb_V_i()=%ld\n", nb_V, s_nb_V_i());
	fflush(stdout);
	for (INT z = 0; z < l; z++) {
		P.m_il(nb_V);
		Q = (PERMUTATION_OP) s_group_generators_i(z);
		printf("z=%ld\n", z);
		Q->println();
		fflush(stdout);
		for (v = 0; v < s_nb_V_i(); v++) 
		{
			INT k, v2;
			INT pii, pij;
			INT new_vertex, image_new_vertex;
			
			k = Q->s_ii(v) - 1;
			printf("v=%ld -> k=%ld\n", v, k); fflush(stdout);
			P.m_ii(v, k + 1);
			adjacency_list(v, adj, &nb_adj);
			for (ei = 0; ei < nb_adj; ei++) 
			{
				e = adj[ei];
				e1 = s_v1_ii(e);
				e2 = s_v2_ii(e);
				if (v == e1) 
				{
					i = e1;
					j = e2;
				}
				if (v == e2)
				{
					i = e2;
					j = e1;
				}
				pii = Q->s_ii(i) - 1;
				pij = Q->s_ii(j) - 1;
				printf("ei=%ld i=%ld ->pii=%ld j=%ld ->pij=%ld\n", 
					ei, i, pii, j, pij); fflush(stdout);
				adjacency_list(pii, adj1, &nb_adj1);
				k = -1;
				for (INT ii = 0; ii < nb_adj1; ii++)
				{
					e1 = adj1[ii];
					v1 = s_v1_ii(e1);
					v2 = s_v2_ii(e1);
					if (v1 == pij || v2 == pij) k = ii;
				}
				if (k == -1) 
					return error("ERROR: image of pij not found!");
				// new_vertex = s_nb_V_i() + v * nb_adj + ei;
				new_vertex = first_new_vertex_number.s_ii(v) + ei;
				// image_new_vertex = s_nb_V_i() + pii * nb_adj + k;
				image_new_vertex = first_new_vertex_number.s_ii(pii) + k;
				printf("new_vertex=%ld image_new_vertex=%ld\n", 
					new_vertex, image_new_vertex); fflush(stdout);
				P.m_ii(new_vertex, image_new_vertex + 1);
			}
		} // next v
		P.swap((PERMUTATION_OP) B.s_group_generators_i(z));
	} // next z
	printf("finished\n");
	fflush(stdout);
	B.s_group_generators()->println();
	fflush(stdout);
	first_new_vertex = s_nb_V_i();
	
	printf("updating the edges:\n"); fflush(stdout);
	for (v = 0; v < s_nb_V_i(); v++) 
	{

		adjacency_list(v, adj, &nb_adj);

		// update the edges (from one side):
		for (ei = 0; ei < nb_adj; ei++) 
		{
			e = adj[ei];
			e1 = s_v1_ii(e);
			e2 = s_v2_ii(e);
			if (e1 == v) 
			{
				B.s_v1_i(e)->m_i(first_new_vertex + ei);
			}
			else if (e2 == v) 
			{
				B.s_v2_i(e)->m_i(first_new_vertex + ei);
			}
			else
				return error("edge does not contain vertex v");
		}

		first_new_vertex += nb_adj;
	}
	if (first_new_vertex != nb_V)
		return error("first_new_vertex != nb_V");

	printf("eliminating the old vertices:\n"); fflush(stdout);
	nb_new_V = nb_V - s_nb_V_i();
	// eliminate the old vertices:
	for (i = 0; i < nb_new_V; i++) 
	{
		B.s_x_i(i)->m_i(B.s_x_ii(s_nb_V_i() + i));
		B.s_y_i(i)->m_i(B.s_y_ii(s_nb_V_i() + i));
		B.s_z_i(i)->m_i(B.s_z_ii(s_nb_V_i() + i));
	}
	for (i = 0; i < B.s_nb_E_i(); i++) 
	{
		if (B.s_v1_ii(i) < s_nb_V_i())
			return error("B.s_v1_ii(i) < s_nb_V_i()");
		if (B.s_v2_ii(i) < s_nb_V_i())
			return error("B.s_v2_ii(i) < s_nb_V_i()");
			
		B.s_v1_i(i)->m_i(B.s_v1_ii(i) - s_nb_V_i());
		B.s_v2_i(i)->m_i(B.s_v2_ii(i) - s_nb_V_i());
	}
	
	Pindv.m_il(l);
	for (INT zz = 0; zz < l; zz++) 
	{
		Pind.m_il(nb_new_V);
		for (v = 0; v < nb_new_V; v++)
		{
			Pind.m_ii(v, ((PERMUTATION_OP) B.s_group_generators_i(zz))
				->s_ii(s_nb_V_i()+v) - s_nb_V_i());
			Pind.copy((PERMUTATION_OP)Pindv.s_i(zz));
		}
	}
	printf("new generators:\n");
	for (INT zzz = 0; zzz < l; zzz++) 
	{
		Pindv.s_i(zzz)->println();
	}
	fflush(stdout);
	Pindv.swap(B.s_group_generators());
	B.s_nb_V()->m_i(nb_new_V);
	B.s_x()->realloc_z(nb_new_V);
	B.s_y()->realloc_z(nb_new_V);
	B.s_z()->realloc_z(nb_new_V);
	if (B.s_f_vertex_labels_i()) {
		B.standard_vertex_labels(FALSE /* f_start_with_zero */);
		}
	B.swap(A);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::Ratio(INT e, double r, SYM_OP Px, SYM_OP Py, SYM_OP Pz)
#endif
{
	INT dx, dy, dz;
	INT e1, e2;

	e1 = s_v1_ii(e);
	e2 = s_v2_ii(e);
	dx = (INT)((double)(s_x_ii(e2) - s_x_ii(e1)) * r);
	dy = (INT)((double)(s_y_ii(e2) - s_y_ii(e1)) * r);
	dz = (INT)((double)(s_z_ii(e2) - s_z_ii(e1)) * r);
	Px->m_i_i(s_x_ii(e1) + dx);
	Py->m_i_i(s_y_ii(e1) + dy);
	Pz->m_i_i(s_z_ii(e1) + dz);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::find_vertex(INTEGER_OB Px, INTEGER_OB Py, INTEGER_OB Pz, INT *vertex_image)
#endif
{
	INT v, vx, vy, vz; 
	INT px, py, pz;
	
	for (v = 0; v < s_nb_V_i(); v++) 
	{
		vx = s_x_ii(v);
		vy = s_y_ii(v);
		vz = s_z_ii(v);
		px = Px.s_i();
		py = Py.s_i();
		pz = Pz.s_i();
		if ((vx == px) && (vy == py) && (vz == pz)) 
		{
			*vertex_image = v;
			return 1;
		}
		return 0;
	}
	if (*vertex_image == -1) 
	{
		printf("find_vertex: no vertex found!\n");
	}
}

#if TEXDOCU
INT SOLID_OB::find_common_face(INT e1, INT e2, INT *f)
#endif
{
	INT n11, n12, n21, n22;
	
	n11 = s_f1_ii(e1);
	n12 = s_f2_ii(e1);
	n21 = s_f1_ii(e2);
	n22 = s_f2_ii(e2);
	if (n11 == n21 || n11 == n22) 
	{
		*f = n11;
		return TRUE;
	}
	if (n12 == n21 || n12 == n22) 
	{
		*f = n12;
		return TRUE;
	}
	return FALSE;
}

#if TEXDOCU
INT SOLID_OB::direct_product(VECTOR_OP gen, SOLID_OP J)
#endif
{
	double f;
	INT d;
	VECTOR_OP gen0;
	SOLID_OB A0, A1, A2;
	INT i;
	
	printf("SOLID::direct_product()\n"); fflush(stdout);
	copy(&A0);
	gen0 = s_group_generators();
	printf("gen=\n");
	gen->println();
	fflush(stdout);
	d = perm_vec_get_degree(gen);
	for (i = 1; i < d; i++){
		f = 1. + (double) i * 1.;
		copy(&A1);
		A1.archimed_scale(f);
		A0.join_disjoint(&A1, &A2);
		A2.swap(&A0);
		}
	printf("calling vec_generators_direct_product\n"); fflush(stdout);
	vec_generators_direct_product(gen, gen0, A0.s_group_generators());
	A0.swap(J);
	

	return OK;
}

#if TEXDOCU
INT SOLID_OB::direct_sum(SOLID_OP B, SOLID_OP J)
#endif
{
	VECTOR_OP gen1, gen2;
	SOLID_OB C;
	
	printf("SOLID::direct_sum()\n"); fflush(stdout);
	gen1 = s_group_generators();
	gen2 = B->s_group_generators();
	join_disjoint(B, &C);
	printf("calling vec_generators_diagonal_sum\n"); fflush(stdout);
	vec_generators_diagonal_sum(gen1, gen2, C.s_group_generators());
	C.swap(J);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::join_disjoint(SOLID_OP A, SOLID_OP J)
#endif
{
	VECTOR_OB gen;
	INT nb_F1, nb_F2, nb_E1, nb_E2, nb_V1, nb_V2, nb_F, nb_E, nb_V;
	INT i, j, v;
	INT ll;
	INT f_vertex_labels = FALSE;
	
	printf("join_disjoint()\n");
	nb_F1 = s_nb_F_i();
	nb_F2 = A->s_nb_F_i();
	nb_E1 = s_nb_E_i();
	nb_E2 = A->s_nb_E_i();
	nb_V1 = s_nb_V_i();
	nb_V2 = A->s_nb_V_i();
	printf("nb_V1 = %ld\n", nb_V1);
	printf("nb_V2 = %ld\n", nb_V2);
	printf("nb_E1 = %ld\n", nb_E1);
	printf("nb_E2 = %ld\n", nb_E2);
	printf("nb_F1 = %ld\n", nb_F1);
	printf("nb_F2 = %ld\n", nb_F2);
	nb_F = nb_F1 + nb_F2;
	nb_E = nb_E1 + nb_E2;
	nb_V = nb_V1 + nb_V2;
	J->init();
	J->init_V(nb_V);
	J->init_E(nb_E);
	J->init_F(nb_F);
	if (s_f_vertex_labels_i() && A->s_f_vertex_labels_i())
		f_vertex_labels = TRUE;
	if (f_vertex_labels) {
		J->s_f_vertex_labels()->m_i(TRUE);
		J->s_vertex_labels()->m_il(nb_V);
		}
	else {
		J->s_f_vertex_labels()->m_i(FALSE);
		}
	printf("1:\n"); fflush(stdout);
	for (i = 0; i < nb_V1; i++)
	{
		J->s_x_i(i)->m_i(s_x_ii(i));
		J->s_y_i(i)->m_i(s_y_ii(i));
		J->s_z_i(i)->m_i(s_z_ii(i));
		if (f_vertex_labels)
			J->s_vertex_labels_i(i)->init(s_vertex_labels_is(i));
	}
	printf("2:\n"); fflush(stdout);
	for (i = nb_V1; i < nb_V; i++)
	{
		j = i - nb_V1;
		J->s_x_i(i)->m_i(A->s_x_ii(j));
		J->s_y_i(i)->m_i(A->s_y_ii(j));
		J->s_z_i(i)->m_i(A->s_z_ii(j));
		if (f_vertex_labels) {
			BYTE str[1000];
			INT a;
			
			sscanf(A->s_vertex_labels_is(j), "%ld", &a);
			a += nb_V1;
			sprintf(str, "%ld", a);
			J->s_vertex_labels_i(i)->init(str);
			}
	}
	printf("3:\n"); fflush(stdout);
	for (i = 0; i < nb_E1; i++)
	{
		J->s_v1_i(i)->m_i(s_v1_ii(i));
		J->s_v2_i(i)->m_i(s_v2_ii(i));
		J->s_f1_i(i)->m_i(s_f1_ii(i));
		J->s_f2_i(i)->m_i(s_f2_ii(i));
	}
	printf("4:\n"); fflush(stdout);
	// J->s_v1()->println();
	// J->s_v2()->println();
	// J->s_f1()->println();
	// J->s_f2()->println();
	for (i = nb_E1; i < nb_E; i++)
	{
		j = i - nb_E1;
		J->s_v1_i(i)->m_i(A->s_v1_ii(j) + nb_V1);
		J->s_v2_i(i)->m_i(A->s_v2_ii(j) + nb_V1);
		J->s_f1_i(i)->m_i(A->s_f1_ii(j) + nb_F1);
		J->s_f2_i(i)->m_i(A->s_f2_ii(j) + nb_F1);
		
	}
	printf("faces 1\n"); fflush(stdout);
	for (i = 0; i < nb_F1; i++)
	{
		ll = s_nb_e_ii(i);
		J->s_nb_e_i(i)->m_i(ll);
		J->s_edge_i(i)->m_il(ll);
		J->s_neighbour_faces_i(i)->m_il(ll);
		for (v = 0; v < ll; v++)
		{
			J->s_edge_ij(i, v)->m_i(s_edge_iji(i, v));
			J->s_neighbour_faces_ij(i, v)->m_i(s_neighbour_faces_iji(i, v));
		}
	}
	printf("faces 2\n"); fflush(stdout);
	for (i = nb_F1; i < nb_F; i++)
	{
		j = i - nb_F1;
		ll = A->s_nb_e_ii(j);
		J->s_nb_e_i(i)->m_i(ll);
		J->s_edge_i(i)->m_il(ll);
		J->s_neighbour_faces_i(i)->m_il(ll);
		for (v = 0; v < ll; v++)
		{
			J->s_edge_ij(i, v)->m_i(A->s_edge_iji(j, v) + nb_E1);
			J->s_neighbour_faces_ij(i, v)->m_i(A->s_neighbour_faces_iji(j, v) + nb_F1);
		}
	}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::relabel_points(SOLID_OP A, PERMUTATION_OP p, 
	INT f_relabel_vertex_labels)
#endif
{
	INT i, l, nb_E, v1, v2;
	
	copy(A);
	l = p->s_li();
	if (l != s_nb_V_i())
		return error("SOLID_OB::relabel_points(): l != s_nb_V_i()");
	A->s_x()->apply_perm(p);
	A->s_y()->apply_perm(p);
	A->s_z()->apply_perm(p);
	if (f_relabel_vertex_labels) {
		if (A->s_f_vertex_labels_i())
			A->s_vertex_labels()->apply_perm(p);
		}
	nb_E = A->s_nb_E_i();
	for (i = 0; i < nb_E; i++) {
		v1 = A->s_v1_ii(i);
		v2 = A->s_v2_ii(i);
		v1 = p->s_ii(v1) - 1;
		v2 = p->s_ii(v2) - 1;
		A->s_v1_i(i)->m_i(v1);
		A->s_v2_i(i)->m_i(v2);
		}
	vec_conjugate(A->s_group_generators(), p);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::cube4D(INT r1, INT r2)
#endif
{
	SOLID_OB A, B, C;
	INT nb_V, nb_E1, nb_F;
	INT i;

	A.cube(r1);
	B.cube(r2);
	A.join_disjoint(&B, &C);
	nb_V = A.s_nb_V_i();
	nb_E1 = C.s_nb_E_i();
	nb_F = A.s_nb_F_i();
	C.s_v1()->realloc_z(nb_E1 + nb_V);
	C.s_v2()->realloc_z(nb_E1 + nb_V);
	C.s_f1()->realloc_z(nb_E1 + nb_V);
	C.s_f2()->realloc_z(nb_E1 + nb_V);
	for (i = 0; i < nb_V; i++) {
		C.s_v1_i(nb_E1 + i)->m_i(i);
		C.s_v2_i(nb_E1 + i)->m_i(nb_V + i);
		C.s_f1_i(nb_E1 + i)->m_i(0);
		C.s_f2_i(nb_E1 + i)->m_i(0);
		}
	C.s_nb_E()->m_i(nb_E1 + nb_V);
	vec_generators_aut_cube_nd(4, C.s_group_generators());
	swap(&C);
	return OK;
}

#if TEXDOCU
INT vec_generators_aut_cube_nd(INT n, VECTOR_OP gen)
#endif
{
	VECTOR_OB gen1, gen2;
	PERMUTATION_OB p;
	PERMUTATION_OP q;
	INT i, j, a, b, x, l, ll;
	INT *v, *w;
	
	// printf("vec_generators_aut_cube_nd()\n");
	l = i_power_j(2, n);
	v = (INT *) my_malloc(n * sizeof(INT), "vec_generators_aut_cube_nd v");
	w = (INT *) my_malloc(n * sizeof(INT), "vec_generators_aut_cube_nd w");
	symmetric_generators(&gen2, n);
	ll = gen2.s_li();
	gen1.m_il(n + ll);
	for (x = 0; x < n; x++) {
		p.m_il(l);
		for (i = 0; i < l; i++) {
			number_to_binary(i, v, n);
			if (v[x])
				v[x] = 0;
			else
				v[x] = 1;
			j = binary_to_number(v, n);
			p.m_ii(i, j + 1);
			}
		// p.println();
		p.swap((PERMUTATION_OP) gen1.s_i(x));
		}
	for (x = 0; x < ll; x++) {
		q = (PERMUTATION_OP) gen2.s_i(x);
		p.m_il(l);
		for (i = 0; i < l; i++) {
			number_to_binary(i, v, n);
			for (a = 0; a < n; a++) {
				b = q->s_ii(a) - 1;
				w[b] = v[a];
				}
			j = binary_to_number(w, n);
			p.m_ii(i, j + 1);
			}
		// p.println();
		p.swap((PERMUTATION_OP) gen1.s_i(n + x));
		}
	my_free(v);
	my_free(w);
	gen1.swap(gen);
	return OK;
}

#undef DEBUG_BINARY_CONVERSION

#if TEXDOCU
INT number_to_binary(INT n, INT *v, INT digits)
#endif
{
	INT i;
#ifdef DEBUG_BINARY_CONVERSION
	INT n1 = n;
#endif
	
	for (i = 0; i < digits; i++) {
		if (ODD(n))
			v[i] = 1;
		else
			v[i] = 0;
		n >>= 1;
		}
#ifdef DEBUG_BINARY_CONVERSION
	printf("number %ld to binary: ", n1);
	for (i = digits - 1; i >= 0; i--) {
		printf("%ld", v[i]);
		}
	printf("\n");
#endif
	return OK;
}

#if TEXDOCU
INT binary_to_number(INT *v, INT digits)
#endif
{
	INT i, n = 0;
	
#ifdef DEBUG_BINARY_CONVERSION
	printf("binary ");
	for (i = digits - 1; i >= 0; i--) {
		printf("%ld", v[i]);
		}
#endif
	for (i = digits - 1; i >= 0; i--) {
		n <<= 1;
		if (v[i])
			n++;
		}
#ifdef DEBUG_BINARY_CONVERSION
	printf(" to number %ld\n", n);
#endif
	return n;
}

#if TEXDOCU
INT SOLID_OB::induced_group_on_edges(VECTOR_OP gen, VECTOR_OP gen_e)
#endif
{
	PERMUTATION_OP p, q;
	INT nb_E, i, l;
	
	nb_E = s_nb_E_i();
	l = gen->s_li();
	gen_e->m_il(l);
	for (i = 0; i < l; i++) {
		p = (PERMUTATION_OP) gen->s_i(i);
		q = (PERMUTATION_OP) gen_e->s_i(i);
		induced_action_on_edges(p, q);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::induced_action_on_edges(PERMUTATION_OP p, PERMUTATION_OP q)
#endif
{
	INT j, nb_E, v1, v2, w1, w2, k;
	
	nb_E = s_nb_E_i();
	q->m_il(nb_E);
	for (j = 0; j < nb_E; j++) {
		v1 = s_v1_ii(j);
		v2 = s_v2_ii(j);
		// printf("edge %ld = (%ld,%ld)", j, v1, v2); fflush(stdout);
		w1 = p->s_ii(v1) - 1;
		w2 = p->s_ii(v2) - 1;
		// printf("maps to (%ld,%ld)\n", w1, w2); fflush(stdout);
		k = find_edge(w1, w2);
		if (k < 0) {
			printf("p=");
			p->println();
			printf("edge %ld = (%ld,%ld) maps to (%ld,%ld)\n", j, v1, v2, w1, w2);
			return error("SOLID_OB::induced_action_on_edges() error in find_edge");
			}
		q->m_ii(j, k + 1);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::induced_group_on_edges_and_faces(VECTOR_OP gen, 
	VECTOR_OP gen_e, VECTOR_OP gen_f)
#endif
{
	INT i, l;
	PERMUTATION_OP p, pe, pf;
	
	l = gen->s_li();
	gen_e->m_il(l);
	gen_f->m_il(l);
	for (i = 0; i < l; i++) {
		p = (PERMUTATION_OP) gen->s_i(i);
		pe = (PERMUTATION_OP) gen_e->s_i(i);
		pf = (PERMUTATION_OP) gen_f->s_i(i);
		induced_action_on_edges_and_faces(p, pe, pf);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::induced_action_on_edges_and_faces(PERMUTATION_OP p, 
	PERMUTATION_OP pe, PERMUTATION_OP pf)
#endif
{
	INT nb_F, j, e1, e2, ee1, ee2, f;
	
	induced_action_on_edges(p, pe);
	nb_F = s_nb_F_i();
	pf->m_il(nb_F);
	for (j = 0; j < nb_F; j++) {
		if (s_nb_e_ii(j) < 2)
			return error("SOLID_OB::induced_action_on_edges_and_faces() s_nb_e_ii(j) < 2");
		e1 = s_edge_iji(j, 0);
		e2 = s_edge_iji(j, 1);
		ee1 = pe->s_ii(e1) - 1;
		ee2 = pe->s_ii(e2) - 1;
		find_common_face(ee1, ee2, &f);
		if (f < 0) {
			return error("SOLID_OB::induced_action_on_edges_and_faces() face not found");
			}
		pf->m_ii(j, f + 1);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::add_central_point(SOLID_OP A)
#endif
{
	SOLID_OB J;
	INT nb_V;
	
	printf("add_central_point()\n");
	copy(&J);
	nb_V = s_nb_V_i();
	vec_add_fixpoint_at_end(J.s_group_generators());
	J.s_nb_V()->inc();
	J.s_x()->inc();
	J.s_y()->inc();
	J.s_z()->inc();
	J.s_x_i(nb_V)->m_i(0);
	J.s_y_i(nb_V)->m_i(0);
	J.s_z_i(nb_V)->m_i(0);
	if (J.s_f_vertex_labels_i()) {
		BYTE str[1000];
		
		sprintf(str, "%ld", nb_V + 1);
		J.s_vertex_labels()->inc();
		J.s_vertex_labels_i(nb_V)->init(str);
		}
	J.swap(A);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::edge_midpoints(SOLID_OP A)
#endif
{
	INT adj[1000], nb_adj;
	INT nb_V, nb_E, nb_E_old, nb_F;
	INT v;
	INT ei, ej, e, ee, f, i, j, l, k;
	INTEGER_OB Px, Py, Pz;
	SOLID_OB B;
	VECTOR_OB gen_new;
	
	
	printf("edge_midpoints()\n");
	induced_group_on_edges(s_group_generators(), &gen_new);
	nb_V = s_nb_V_i();
	nb_E = nb_E_old = s_nb_E_i();
	nb_F = s_nb_F_i();
	copy(&B);
	l = s_group_generators()->s_li();

	B.s_x()->realloc_z(nb_V + nb_E);
	B.s_y()->realloc_z(nb_V + nb_E);
	B.s_z()->realloc_z(nb_V + nb_E);
	for (i = 0; i < nb_E; i++) {
		Ratio(i, 0.5, B.s_x_i(nb_V + i), B.s_y_i(nb_V + i), B.s_z_i(nb_V + i));
		}
	// B.s_nb_F()->m_i(nb_F + nb_V);
	B.s_nb_e()->realloc_z(nb_F + nb_V);
	B.s_edge()->realloc_z(nb_F + nb_V);
	B.s_neighbour_faces()->realloc_z(nb_F + nb_V);
	
	// define new faces:
	for (v = 0; v < nb_V; v++) 
	{
		adjacency_list(v, adj, &nb_adj);
		
		B.s_nb_e_i(nb_F + v)->m_i(0);
		B.s_edge_i(nb_F + v)->m_il(nb_adj);
		B.s_neighbour_faces_i(nb_F + v)->m_il(nb_adj);
		
		// all the edges of the new face:
		for (ei = 0; ei < nb_adj; ei++) 
		{
			e = adj[ei];
			for (ej = ei + 1; ej < nb_adj; ej++) 
			{
				ee = adj[ej];
				if (find_common_face(e, ee, &f)) 
				{
					
					// new edge:
					B.s_v1()->inc();
					B.s_v1_i(nb_E)->m_i(nb_V + e);
					B.s_v2()->inc();
					B.s_v2_i(nb_E)->m_i(nb_V + ee);
					B.s_f1()->inc();
					B.s_f1_i(nb_E)->m_i(f);
					B.s_f2()->inc();
					B.s_f2_i(nb_E)->m_i(nb_F + v);
					B.s_nb_E()->inc();
					
					// new face:
					
					k = B.s_nb_e_ii(nb_F + v);
					B.s_edge_ij(nb_F + v, k)->m_i(nb_E);
					B.s_neighbour_faces_ij(nb_F + v, k)->m_i(f);
					B.s_nb_e_i(nb_F + v)->inc();
					
					// old face:
					k = B.s_nb_e_ii(f);
					B.s_edge_i(f)->inc();
					B.s_neighbour_faces_i(f)->inc();
					B.s_edge_ij(f, k)->m_i(nb_E);
					B.s_neighbour_faces_ij(f, k)->m_i(nb_F + v);
					B.s_nb_e_i(f)->inc();
					
					nb_E++;
				}
			}
		}
		
		
		// face:
		B.s_nb_F()->inc();
	} // next v
	
	// eliminate old edges in old faces:
	for (f = 0; f < nb_F; f++) {
		VECTOR_OB edge, neighbour_faces;
		INT n;
		
		k = B.s_nb_e_ii(f);
		edge.m_il(0);
		neighbour_faces.m_il(0);
		i = 0;
		for (j = 0; j < k; j++) {
			e = B.s_edge_iji(f, j);
			n = B.s_neighbour_faces_iji(f, j);
			if (e >= nb_E_old) {
				edge.inc();
				neighbour_faces.inc();
				edge.m_ii(i, e);
				neighbour_faces.m_ii(i, n);
				i++;
				}
			}
		edge.swap(B.s_edge_i(f));
		neighbour_faces.swap(B.s_neighbour_faces_i(f));
		B.s_nb_e_i(f)->m_i(i);
		}


	printf("eliminating the old vertices:\n"); fflush(stdout);
	// eliminate the old vertices:
	for (i = 0; i < nb_E_old; i++) 
	{
		B.s_x_i(i)->m_i(B.s_x_ii(nb_V + i));
		B.s_y_i(i)->m_i(B.s_y_ii(nb_V + i));
		B.s_z_i(i)->m_i(B.s_z_ii(nb_V + i));
	}
	B.s_nb_V()->m_i(nb_E_old);
	
	// eliminate old edges and update vertex labels of new edges:
	printf("eliminate old edges and update vertex labels of new edges:\n"); fflush(stdout);
	for (i = nb_E_old; i < B.s_nb_E_i(); i++) 
	{
		if (B.s_v1_ii(i) < nb_V)
			return error("B.s_v1_ii(i) < nb_V");
		if (B.s_v2_ii(i) < nb_V)
			return error("B.s_v2_ii(i) < nb_V");
			
		B.s_v1_i(i - nb_E_old)->m_i(B.s_v1_ii(i) - nb_V);
		B.s_v2_i(i - nb_E_old)->m_i(B.s_v2_ii(i) - nb_V);
		B.s_f1_i(i - nb_E_old)->m_i(B.s_f1_ii(i));
		B.s_f2_i(i - nb_E_old)->m_i(B.s_f2_ii(i));
	}
	B.s_nb_E()->m_i(B.s_nb_E_i() - nb_E_old);
	
	for (f = 0; f < B.s_nb_F_i(); f++) {
		INT n;
		
		k = B.s_nb_e_ii(f);
		for (j = 0; j < k; j++) {
			e = B.s_edge_iji(f, j);
			n = B.s_neighbour_faces_iji(f, j);
			if (e < nb_E_old) {
				printf("f=%ld j=%ld k=%ld e=%ld n=%ld\n", f, j, k, e, n);
				return error("SOLID::edge_midpoints: e < nb_E_old");
				}
			B.s_edge_ij(f, j)->m_i(e - nb_E_old);
			}
		}
	
	
	
	gen_new.swap(B.s_group_generators());

	B.s_x()->realloc_z(B.s_nb_V_i());
	B.s_y()->realloc_z(B.s_nb_V_i());
	B.s_z()->realloc_z(B.s_nb_V_i());
	if (B.s_f_vertex_labels_i()) {
		B.standard_vertex_labels(FALSE /* f_start_with_zero */);
		}
	B.s_v1()->realloc_z(B.s_nb_E_i());
	B.s_v2()->realloc_z(B.s_nb_E_i());
	B.s_f1()->realloc_z(B.s_nb_E_i());
	B.s_f2()->realloc_z(B.s_nb_E_i());
	
	B.swap(A);
	return OK;
}

static INT search_and_insert_orbit_below(VECTOR_OP orbit_below1, 
	VECTOR_OP orbit_below2, INT o)
{
	INTEGER_OB int_ob;
	INT idx, f_found, l, ii;
	
	int_ob.m_i(o);
	orbit_below1->search(orbit_below1->s_li(), TRUE, &int_ob, &idx, &f_found);
	if (f_found) {
		idx--;
		orbit_below2->s_i(idx)->inc();
		}
	else {
		l = orbit_below1->s_li();
		orbit_below1->inc();
		orbit_below2->inc();
		for (ii = l; ii > idx; ii--) {
			orbit_below1->s_i(ii)->swap(orbit_below1->s_i(ii - 1));
			orbit_below2->s_i(ii)->swap(orbit_below2->s_i(ii - 1));
			}
		orbit_below1->m_ii(idx, o);
		orbit_below2->m_ii(idx, 1);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::plesken(VECTOR_OP Orbits, VECTOR_OP Orbit_ago, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2, 
	INT check_involution, INT f_v)
#endif
{
	INT nb_V, nb_E, nb_F, nb_o_v, nb_o_e, nb_o_f, nb_o;
	INT i, j, v1, v2, o1, o2, l, ii, e;
	VECTOR_OP gen_v;
	VECTOR_OB gen_e, gen_f;
	VECTOR_OB orbit_v;
	VECTOR_OB orbit_e;
	VECTOR_OB orbit_f;
	VECTOR_OB SVlast_v, SVgen_v;
	VECTOR_OB SVlast_e, SVgen_e;
	VECTOR_OB SVlast_f, SVgen_f;
	VECTOR_OB orbit_v_first, orbit_v_size, orbit_v_ago;
	VECTOR_OB orbit_e_first, orbit_e_size, orbit_e_ago;
	VECTOR_OB orbit_f_first, orbit_f_size, orbit_f_ago;
	VECTOR_OB ol, obelow1, obelow2, ob1, ob2;
	VECTOR_OB tmp_v;
	INTEGER_OB tmp_a;
	SYM_OB go;
	LABRA_OB L;
	PERMUTATION_OB inv_v, inv_e, inv_f;
	
	nb_V = s_nb_V_i();
	nb_E = s_nb_E_i();
	nb_F = s_nb_F_i();

	gen_v = s_group_generators();
	L.build(gen_v);
	L.group_order(&go);

	induced_group_on_edges_and_faces(s_group_generators(), &gen_e, &gen_f);
	if (f_v) {
		printf("SOLID::plesken induced group on edges:\n");
		gen_e.Print();
		printf("SOLID::plesken induced group on faces:\n");
		gen_f.Print();
		printf("\n");
		fflush(stdout);
		}
	if (check_involution) {
		get_central_involution(&inv_v);
		induced_action_on_edges_and_faces(&inv_v, &inv_e, &inv_f);
		}
	
	orbits(gen_v, &orbit_v, &SVlast_v, &SVgen_v, &orbit_v_first, &orbit_v_size);
	orbits(&gen_e, &orbit_e, &SVlast_e, &SVgen_e, &orbit_e_first, &orbit_e_size);
	orbits(&gen_f, &orbit_f, &SVlast_f, &SVgen_f, &orbit_f_first, &orbit_f_size);
		
	SVlast_v.dec_all_entries();
	SVlast_e.dec_all_entries();
	SVlast_f.dec_all_entries();
	orbit_v_first.dec_all_entries();
	orbit_e_first.dec_all_entries();
	orbit_f_first.dec_all_entries();
	l = orbit_v_size.s_li();
	orbit_v_ago.m_il(l);
	for (i = 0; i < l; i++) {
		go.ganzdiv_integral(orbit_v_size.s_i(i), orbit_v_ago.s_i(i));
		}
	l = orbit_e_size.s_li();
	orbit_e_ago.m_il(l);
	for (i = 0; i < l; i++) {
		go.ganzdiv_integral(orbit_e_size.s_i(i), orbit_e_ago.s_i(i));
		}
	l = orbit_f_size.s_li();
	orbit_f_ago.m_il(l);
	for (i = 0; i < l; i++) {
		go.ganzdiv_integral(orbit_f_size.s_i(i), orbit_f_ago.s_i(i));
		}
	
	nb_o_v = orbit_v_first.s_li();
	nb_o_e = orbit_e_first.s_li();
	nb_o_f = orbit_f_first.s_li();
	nb_o = 1 + nb_o_v + nb_o_e + nb_o_f + 1;

	if (check_involution) {
		for (i = 0; i < nb_o_v; i++) {
			INT j, k, o;
			
			j = orbit_v_first.s_ii(i);
			k = inv_v.s_ii(j) - 1;
			o = orbit_v.s_ii(k);
			if (o != i) {
				printf("vertex orbit %ld and %ld fuse !\n", i, o);
				}
			}
		for (i = 0; i < nb_o_e; i++) {
			INT j, k, o;
			
			j = orbit_e_first.s_ii(i);
			k = inv_e.s_ii(j) - 1;
			o = orbit_e.s_ii(k);
			if (o != i) {
				printf("edge orbit %ld and %ld fuse !\n", i, o);
				}
			}
		for (i = 0; i < nb_o_f; i++) {
			INT j, k, o;
			
			j = orbit_f_first.s_ii(i);
			k = inv_f.s_ii(j) - 1;
			o = orbit_f.s_ii(k);
			if (o != i) {
				printf("face orbit %ld and %ld fuse !\n", i, o);
				}
			}
		}

	Orbits->m_il(5);
	tmp_v.m_il(1);
	tmp_v.m_ii(0, 0);
	tmp_v.copy((VECTOR_OP) Orbits->s_i(0));
	orbit_v_first.copy((VECTOR_OP) Orbits->s_i(1));
	orbit_e_first.copy((VECTOR_OP) Orbits->s_i(2));
	orbit_f_first.copy((VECTOR_OP) Orbits->s_i(3));
	tmp_v.copy((VECTOR_OP) Orbits->s_i(4));
	
	Orbit_ago->m_il(5);
	ol.m_il(1);
	go.copy(ol.s_i(0));
	ol.copy((VECTOR_OP) Orbit_ago->s_i(0));
	orbit_v_ago.copy((VECTOR_OP) Orbit_ago->s_i(1));
	orbit_e_ago.copy((VECTOR_OP) Orbit_ago->s_i(2));
	orbit_f_ago.copy((VECTOR_OP) Orbit_ago->s_i(3));
	ol.copy((VECTOR_OP) Orbit_ago->s_i(4));
	
	Orbits_below1->m_il(5);
	Orbits_below2->m_il(5);
	obelow1.m_il(nb_o_v);	
	obelow2.m_il(nb_o_v);
	for (i = 0; i < nb_o_v; i++) {
		ob1.m_il(1);
		ob2.m_il(1);
		ob1.m_ii(0, 0);
		ob2.m_ii(0, 1);
		ob1.swap((VECTOR_OP) obelow1.s_i(i));
		ob2.swap((VECTOR_OP) obelow2.s_i(i));
		}
	obelow1.swap((VECTOR_OP) Orbits_below1->s_i(1));
	obelow2.swap((VECTOR_OP) Orbits_below2->s_i(1));
	
	obelow1.m_il(nb_o_e);	
	obelow2.m_il(nb_o_e);
	for (i = 0; i < nb_o_e; i++) {
		ob1.m_il(0);
		ob2.m_il(0);
		j = orbit_e_first.s_ii(i);
		v1 = s_v1_ii(j);
		v2 = s_v2_ii(j);
		o1 = orbit_v.s_ii(v1);
		o2 = orbit_v.s_ii(v2);
		search_and_insert_orbit_below(&ob1, &ob2, o1);  
		search_and_insert_orbit_below(&ob1, &ob2, o2);  
		ob1.swap((VECTOR_OP) obelow1.s_i(i));
		ob2.swap((VECTOR_OP) obelow2.s_i(i));
		}
	obelow1.swap((VECTOR_OP) Orbits_below1->s_i(2));
	obelow2.swap((VECTOR_OP) Orbits_below2->s_i(2));
	
	obelow1.m_il(nb_o_f);	
	obelow2.m_il(nb_o_f);
	for (i = 0; i < nb_o_f; i++) {
		ob1.m_il(0);
		ob2.m_il(0);
		j = orbit_f_first.s_ii(i);
		l = s_nb_e_ii(j);
		for (ii = 0; ii < l; ii++) {
			e = s_edge_iji(j, ii);
			o1 = orbit_e.s_ii(e);
			search_and_insert_orbit_below(&ob1, &ob2, o1);  
			}
		ob1.swap((VECTOR_OP) obelow1.s_i(i));
		ob2.swap((VECTOR_OP) obelow2.s_i(i));
		}
	obelow1.swap((VECTOR_OP) Orbits_below1->s_i(3));
	obelow2.swap((VECTOR_OP) Orbits_below2->s_i(3));
	
	obelow1.m_il(1);	
	obelow2.m_il(1);
	ob1.m_il(nb_o_f);
	ob2.m_il(nb_o_f);
	for (i = 0; i < nb_o_f; i++) {
		ob1.m_ii(i, i);
		ob2.m_ii(i, orbit_f_size.s_ii(i));
		}
	ob1.swap((VECTOR_OP) obelow1.s_i(0));
	ob2.swap((VECTOR_OP) obelow2.s_i(0));
	obelow1.swap((VECTOR_OP) Orbits_below1->s_i(4));
	obelow2.swap((VECTOR_OP) Orbits_below2->s_i(4));
	
	return OK;
}

#if TEXDOCU
INT SOLID_OB::get_central_involution(PERMUTATION_OP p)
#endif
{
	double Y[3][3] = { { -1., 0., 0. }, { 0., -1., 0. }, { 0., 0., -1. } };
	double y0[3] = { 0., 0., 0. };
	
	get_automorphism(Y, y0, p);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::get_automorphism(double Y[3][3], double y0[3], PERMUTATION_OP p)
#endif
{
	INT i, j, nb_V, x, y, z;
	
	nb_V = s_nb_V_i();
	p->m_il(nb_V);
	for (i = 0; i < nb_V; i++) {
		apply_motion(Y, y0, i, &x, &y, &z);
		j = find_point(x, y, z);
		if (j == -1)
			return error("SOLID_OB::get_automorphism() point not found");
		p->m_ii(i, j + 1);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::find_point(INT x, INT y, INT z)
#endif
{
	INT nb_V, i;
	double d, a;
	
	nb_V = s_nb_V_i();
	for (i = 0; i < nb_V; i++) {
		a = (double) (s_x_ii(i) - x);
		d = a * a;
		a = (double) (s_y_ii(i) - y);
		d += a * a;
		a = (double) (s_z_ii(i) - z);
		d += a * a;
		d = sqrt(d);
		if (d < EPSILON) {
			return i;
			}
		}
	return -1;
}

#if TEXDOCU
INT SOLID_OB::identify_points(SOLID_OP A, VECTOR_OP map, VECTOR_OP new_points, 
	INT *nb_new_points)
#endif
{
	INT i, j, l, k = 0;
	
	l = A->s_nb_V_i();
	map->m_il(l);
	new_points->m_il_n(l);
	for (i = 0; i < l; i++) {
		j = find_point(A->s_x_ii(i), A->s_y_ii(i), A->s_z_ii(i));
		// printf("identify_points: i=%ld j=%ld\n", i, j);
		if (j >= 0)
			map->m_ii(i, j);
		else {
			map->m_ii(i, -1);
			new_points->m_ii(i, k++);
			}
		}
	*nb_new_points = k;
	return OK;
}

#if TEXDOCU
INT SOLID_OB::apply_motion(double Rot[3][3], double x0[3], INT i, 
	INT *px, INT *py, INT *pz)
#endif
{
	INT ii;
	double x[3] = {0., 0., 0.};
	
	for (ii = 0; ii < 3; ii++) {
		x[ii] += Rot[ii][0] * (double) s_x_ii(i);
		x[ii] += Rot[ii][1] * (double) s_y_ii(i);
		x[ii] += Rot[ii][2] * (double) s_z_ii(i);
		}
	for (ii = 0; ii < 3; ii++) {
		x[ii] += x0[ii];
		}
	*px = (INT) x[0];
	*py = (INT) x[1];
	*pz = (INT) x[2];
	return OK;
}

#if TEXDOCU
INT SOLID_OB::apply_motion_to_all_points(double Rot[3][3], double x0[3])
#endif
{
	INT i, l, x, y, z;
	
	l = s_nb_V_i();
	for (i = 0; i < l; i++) {
		apply_motion(Rot, x0, i, &x, &y, &z);
		s_x_i(i)->m_i(x);
		s_y_i(i)->m_i(y);
		s_z_i(i)->m_i(z);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::determine_motion(double Y[3][3], double y0[3], 
	INT i1, INT i2, INT i3, INT i4, 
	SOLID_OP S2, INT j1, INT j2, INT j3, INT j4, 
	INT f_v)
#endif
{
	double P[8][3];
	
	P[0][0] = (double) s_x_ii(i1);
	P[0][1] = (double) s_y_ii(i1);
	P[0][2] = (double) s_z_ii(i1);
	P[1][0] = (double) s_x_ii(i2);
	P[1][1] = (double) s_y_ii(i2);
	P[1][2] = (double) s_z_ii(i2);
	P[2][0] = (double) s_x_ii(i3);
	P[2][1] = (double) s_y_ii(i3);
	P[2][2] = (double) s_z_ii(i3);
	P[3][0] = (double) s_x_ii(i4);
	P[3][1] = (double) s_y_ii(i4);
	P[3][2] = (double) s_z_ii(i4);
	P[4][0] = (double) S2->s_x_ii(j1);
	P[4][1] = (double) S2->s_y_ii(j1);
	P[4][2] = (double) S2->s_z_ii(j1);
	P[5][0] = (double) S2->s_x_ii(j2);
	P[5][1] = (double) S2->s_y_ii(j2);
	P[5][2] = (double) S2->s_z_ii(j2);
	P[6][0] = (double) S2->s_x_ii(j3);
	P[6][1] = (double) S2->s_y_ii(j3);
	P[6][2] = (double) S2->s_z_ii(j3);
	P[7][0] = (double) S2->s_x_ii(j4);
	P[7][1] = (double) S2->s_y_ii(j4);
	P[7][2] = (double) S2->s_z_ii(j4);
	::determine_motion(P, Y, y0, 0, 1, 2, 3, 4, 5, 6, 7, 1., f_v);
	return OK;
}

#if TEXDOCU
INT SOLID_OB::get_vertices_of_face(INT i, VECTOR_OP V)
#endif
{
	INT e, nb_e, j, v1, v2, v3, v4, last_vertex;
	
	nb_e = s_nb_e_ii(i);
	V->m_il_n(nb_e);
	
	e = s_edge_iji(i, 0);
	v1 = s_v1_ii(e);
	v2 = s_v2_ii(e);
	printf("e=%ld v1=%ld v2=%ld\n", e, v1, v2); fflush(stdout);
	for (j = 1; j < nb_e; j++) {
		e = s_edge_iji(i, j);
		v3 = s_v1_ii(e);
		v4 = s_v2_ii(e);
		printf("j=%ld e=%ld v3=%ld v4=%ld\n", j, e, v3, v4); fflush(stdout);
		if (j == 1) {
			if (v3 == v2) {
				V->m_ii(0, v1);
				V->m_ii(1, v2);
				last_vertex = v4;
				}
			else if (v3 == v1) {
				V->m_ii(0, v2);
				V->m_ii(1, v1);
				last_vertex = v4;
				}
			else if (v4 == v2) {
				V->m_ii(0, v1);
				V->m_ii(1, v2);
				last_vertex = v3;
				}
			else if (v4 == v1) {
				V->m_ii(0, v2);
				V->m_ii(1, v1);
				last_vertex = v3;
				}
			else {
				printf("face %ld: edges=", i);
				s_edge_i(i)->println();
				printf("v1=%ld v2=%ld v3=%ld v4=%ld last_vertex=%ld\n", v1, v2, v3, v4, last_vertex);
				V->println();
				fflush(stdout);
				return error("SOLID_OB::get_vertices_of_face() error edges not adjacent! (j=1)");
				}
			}
		else {
			if (v3 == last_vertex) {
				V->m_ii(j, v3);
				last_vertex = v4;
				}
			else if (v4 == last_vertex) {
				V->m_ii(j, v4);
				last_vertex = v3;
				}
			else {
				printf("face %ld: edges=", i);
				s_edge_i(i)->println();
				printf("v1=%ld v2=%ld v3=%ld v4=%ld last_vertex=%ld\n", v1, v2, v3, v4, last_vertex);
				V->println();
				fflush(stdout);
				return error("SOLID_OB::get_vertices_of_face() error edges not adjacent!");
				}
			}
		printf("last_vertex=%ld\n", last_vertex);
		}
	return OK;
}

#if TEXDOCU
INT SOLID_OB::join_with(SOLID_OP B, INT f_v)
#endif
{
	VECTOR_OB map, new_points;
	INT nb_new_points;
	INT nb_V;
	INT i, j, e, nb_e, v1, v2, vv1, vv2, k, ee, e1, e2, f;
	VECTOR_OB old_vertices, new_vertices;
	INT face_added = FALSE;
	
	identify_points(B, &map, &new_points, &nb_new_points);
	if (f_v) {
		printf("identify points: map=");
		map.println();
		printf("new_points=");
		new_points.println();
		fflush(stdout);
		}
	nb_V = s_nb_V_i();
	s_x()->realloc_z(nb_V + nb_new_points);
	s_y()->realloc_z(nb_V + nb_new_points);
	s_z()->realloc_z(nb_V + nb_new_points);
	k = 0;
	for (i = 0; i < B->s_nb_V_i(); i++) {
		if (map.s_ii(i) == -1) {
			s_x_i(nb_V + k)->m_i(B->s_x_ii(i));
			s_y_i(nb_V + k)->m_i(B->s_y_ii(i));
			s_z_i(nb_V + k)->m_i(B->s_z_ii(i));
			map.m_ii(i, nb_V + k);
			k++;
			}
		}
	if (k != nb_new_points) {
		return error("SOLID::join_with() k != nb_new_points");
		}
	s_nb_V()->m_i(nb_V + nb_new_points);
	// s_vertex_labels()->realloc_z(nb_V + nb_new_points);
	standard_vertex_labels(FALSE /* f_start_with_zero */);
	
	for (i = 0; i < B->s_nb_F_i(); i++) {
		nb_e = B->s_nb_e_ii(i);
		if (f_v) {
			printf("face %ld:\n", i); fflush(stdout);
			}
		for (j = 0; j < nb_e; j++) {
			e = B->s_edge_iji(i, j);
			v1 = B->s_v1_ii(e);
			v2 = B->s_v2_ii(e);
			vv1 = map.s_ii(v1);
			vv2 = map.s_ii(v2);
			ee = find_and_add_edge(vv1, vv2, FALSE /* f_v */);
			if (f_v) {
				printf("new edge %ld (%ld,%ld)\n", ee, vv1, vv2);
				}
			if (j == 0)
				e1 = ee;
			else if (j == 1)
				e2 = ee;
			}
		f = find_face_by_two_edges(e1, e2);
		if (f == -1) {
			B->get_vertices_of_face(i, &old_vertices);
			if (f_v) {
				printf("old face: ");
				old_vertices.println();
				}
			new_vertices.m_il(nb_e);
			for (j = 0; j < nb_e; j++) {
				k = old_vertices.s_ii(j);
				new_vertices.m_ii(j, map.s_ii(k));
				}
			if (f_v) {
				printf("new face: ");
				new_vertices.println();
				}
			add_face_n(&new_vertices);
			if (f_v) {
				INT f;
				
				f = s_nb_F_i() - 1;
				printf("new face %ld", f);
				s_edge_i(f)->println();
				}
			face_added = TRUE;
			}
		} // next i
	printf("SOLID::join_with() calling determine_neighbours()\n");
	determine_neighbours();
	return face_added;
}

#endif /* SOLID_TRUE */
