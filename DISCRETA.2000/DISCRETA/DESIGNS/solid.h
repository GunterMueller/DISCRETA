/* solid.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten, Evi Haberberger 1999
 */


#ifndef SOLID_INCLUDED
#define SOLID_INCLUDED

#include <DISCRETA/divs.h>

/* solid.c */

class solid_ob : public VECTOR_OB {
public:
	VECTOR_OP s_group_generators() {
		return((VECTOR_OP)s_i(0));
	};
	SYM_OP s_group_generators_i(INT i) {
		return(s_group_generators()->s_i(i));
	}
	INTEGER_OP s_nb_V() {
		return((INTEGER_OP)s_i(1));
	};
	INT s_nb_V_i() {
		return(s_nb_V()->s_i());
	};
	VECTOR_OP s_placement() {
		return((VECTOR_OP)s_i(2));
	};
	VECTOR_OP s_x() {
		return((VECTOR_OP)s_placement()->s_i(0));
	};
	INTEGER_OP s_x_i(INT i) {
		return((INTEGER_OP)s_x()->s_i(i));
	};
	INT s_x_ii(INT i) {
		return(s_x_i(i)->s_i());
	};
	VECTOR_OP s_y() {
		return((VECTOR_OP)s_placement()->s_i(1));
	};
	INTEGER_OP s_y_i(INT i) {
		return((INTEGER_OP)s_y()->s_i(i));
	};
	INT s_y_ii(INT i) {
		return(s_y_i(i)->s_i());
	};
	VECTOR_OP s_z() {
		return((VECTOR_OP)s_placement()->s_i(2));
	};
	INTEGER_OP s_z_i(INT i) {
		return((INTEGER_OP)s_z()->s_i(i));
	};
	INT s_z_ii(INT i) {
		return(s_z_i(i)->s_i());
	};
	INTEGER_OP s_nb_E() {
		return((INTEGER_OP)s_i(3));
	};
	INT s_nb_E_i() {
		return(s_nb_E()->s_i());
	};
	VECTOR_OP s_v1() {
		return((VECTOR_OP)s_i(4));
	};
	INTEGER_OP s_v1_i(INT i) {
		return((INTEGER_OP)s_v1()->s_i(i));
	};
	INT s_v1_ii(INT i) {
		return(s_v1_i(i)->s_i());
	};
	VECTOR_OP s_v2() {
		return((VECTOR_OP)s_i(5));
	};
	INTEGER_OP s_v2_i(INT i) {
		return((INTEGER_OP)s_v2()->s_i(i));
	};
	INT s_v2_ii(INT i) {
		return(s_v2_i(i)->s_i());
	};
	VECTOR_OP s_f1() {
		return((VECTOR_OP)s_i(6));
	};
	INTEGER_OP s_f1_i(INT i) {
		return((INTEGER_OP)s_f1()->s_i(i));
	};
	INT s_f1_ii(INT i) {
		return(s_f1_i(i)->s_i());
	};
	VECTOR_OP s_f2() {
		return((VECTOR_OP)s_i(7));
	};
	INTEGER_OP s_f2_i(INT i) {
		return((INTEGER_OP)s_f2()->s_i(i));
	};
	INT s_f2_ii(INT i) {
		return(s_f2_i(i)->s_i());
	};
	INTEGER_OP s_nb_F() {
		return((INTEGER_OP)s_i(8));
	};
	INT s_nb_F_i() {
		return(s_nb_F()->s_i());
	};
	VECTOR_OP s_nb_e() {
		return((VECTOR_OP)s_i(9));
	};
	INTEGER_OP s_nb_e_i(INT i) {
		return((INTEGER_OP)s_nb_e()->s_i(i));
	};
	INT s_nb_e_ii(INT i) {
		return(s_nb_e_i(i)->s_i());
	};
	VECTOR_OP s_edge() {
		return((VECTOR_OP)s_i(10));
	};
	VECTOR_OP s_edge_i(INT i) {
		return((VECTOR_OP)s_edge()->s_i(i));
	};
	INTEGER_OP s_edge_ij(INT i, INT j) {
		return((INTEGER_OP)s_edge_i(i)->s_i(j));
	};
	INT s_edge_iji(INT i, INT j) {
		return(s_edge_ij(i, j)->s_i());
	};
	VECTOR_OP s_neighbour_faces() {
		return((VECTOR_OP)s_i(11));
	};
	VECTOR_OP s_neighbour_faces_i(INT i) {
		return((VECTOR_OP)s_neighbour_faces()->s_i(i));
	};
	INTEGER_OP s_neighbour_faces_ij(INT i, INT j) {
		return((INTEGER_OP)s_neighbour_faces_i(i)->s_i(j));
	};
	INT s_neighbour_faces_iji(INT i, INT j) {
		return(s_neighbour_faces_ij(i, j)->s_i());
	};
	INTEGER_OP s_f_vertex_labels() {
		return((INTEGER_OP)s_i(12));
	};
	INT s_f_vertex_labels_i() {
		return(s_f_vertex_labels()->s_i());
	};
	VECTOR_OP s_vertex_labels() {
		return((VECTOR_OP)s_i(13));
	};
	STRING_OP s_vertex_labels_i(INT i) {
		return((STRING_OP)s_vertex_labels()->s_i(i));
	};
	BYTE *s_vertex_labels_is(INT i) {
		return(s_vertex_labels_i(i)->s_str());
	};
INT init();
INT init_V(INT nb_V);
INT init_E(INT nb_E);
INT init_F(INT nb_F);
INT standard_vertex_labels(INT f_start_with_zero);
INT sprint(BYTE *s);
INT cubus_simus(INT s);
INT dode_simum(INT s);
INT snub_cube(INT s);
INT russian_snub_cube(INT s);
INT tetrahedron(INT r);
INT cube(INT r);
INT octahedron(INT r);
INT dodecahedron(INT r);
INT icosahedron(INT r);
INT dual(SOLID_OP A);
INT center(INT f, VECTOR_OP Px, VECTOR_OP Py, VECTOR_OP Pz);
INT adjacency_list(INT vertex, INT *adj, INT *nb_adj);
INT add_edge(INT v1, INT v2, INT f1, INT f2);
INT add_face3(INT e1, INT e2, INT e3, INT n1, INT n2, INT n3);
INT add_face4(INT i1, INT i2, INT i3, INT i4);
INT add_face5(INT i1, INT i2, INT i3, INT i4, INT i5);
INT add_face_n(VECTOR_OP vertices);
INT find_and_add_edge(INT i1, INT i2, INT f_v);
INT find_edge(INT v1, INT v2);
INT add_edge(INT v1, INT v2);
INT find_face_by_two_edges(INT e1, INT e2);
INT find_faces_at_edge(INT e, INT *f1, INT *f2);
INT find_face(INT e, INT *f1, INT *j1, INT *f2, INT *j2);
INT find_face_2(INT e1, INT e2);
INT determine_neighbours();
INT write_graphfile(BYTE *fname);
INT archimed_scale(double f);
INT cut_vertices(double r, SOLID_OP A);
INT Ratio(INT e, double r, SYM_OP Px, SYM_OP Py, SYM_OP Pz);
INT find_vertex(INTEGER_OB Px, INTEGER_OB Py, INTEGER_OB Pz, INT *vertex_image);
INT find_common_face(INT e1, INT e2, INT *f);
INT direct_product(VECTOR_OP gen, SOLID_OP J);
INT direct_sum(SOLID_OP B, SOLID_OP J);
INT join_disjoint(SOLID_OP A, SOLID_OP J);
INT relabel_points(SOLID_OP A, PERMUTATION_OP p, INT f_relabel_vertex_labels);
INT cube4D(INT r1, INT r2);
INT induced_group_on_edges(VECTOR_OP gen, VECTOR_OP gen_e);
INT induced_action_on_edges(PERMUTATION_OP p, PERMUTATION_OP q);
INT induced_group_on_edges_and_faces(VECTOR_OP gen, 
	VECTOR_OP gen_e, VECTOR_OP gen_f);
INT induced_action_on_edges_and_faces(PERMUTATION_OP p, 
	PERMUTATION_OP pe, PERMUTATION_OP pf);
INT add_central_point(SOLID_OP A);
INT edge_midpoints(SOLID_OP A);
INT plesken(VECTOR_OP Orbits, VECTOR_OP Orbit_ago, 
	VECTOR_OP Orbits_below1, VECTOR_OP Orbits_below2, 
	INT check_involution, INT f_v);
INT get_central_involution(PERMUTATION_OP p);
INT get_automorphism(double Y[3][3], double y0[3], PERMUTATION_OP p);
INT find_point(INT x, INT y, INT z);
INT identify_points(SOLID_OP A, VECTOR_OP map, VECTOR_OP new_points, 
	INT *nb_new_points);
INT apply_motion(double Rot[3][3], double x0[3], INT i, 
	INT *px, INT *py, INT *pz);
INT apply_motion_to_all_points(double Rot[3][3], double x0[3]);
INT determine_motion(double Y[3][3], double y0[3], 
	INT i1, INT i2, INT i3, INT i4, 
	SOLID_OP S2, INT j1, INT j2, INT j3, INT j4, 
	INT f_v);
INT get_vertices_of_face(INT i, VECTOR_OP V);
INT join_with(SOLID_OP B, INT f_v);
};

INT vec_generators_aut_cube_nd(INT n, VECTOR_OP gen);
INT number_to_binary(INT n, INT *v, INT digits);
INT binary_to_number(INT *v, INT digits);

#endif /* SOLID_INCLUDED */

