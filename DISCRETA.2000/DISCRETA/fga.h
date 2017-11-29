/* fga.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef FGA_INCLUDED
#define FGA_INCLUDED

#ifndef VDI_INCLUDED
#include <DISCRETA/graphics.h>
#endif
#ifndef DIVS_INCLUDED
#include <DISCRETA/divs.h>
#endif

/*
 * FGA - finite group action
 */

#define FGA_GROUP_SYM 0
#define FGA_GROUP_SYM_TRANSPOSITIONS 1
#define FGA_GROUP_ALT 2
#define FGA_GROUP_DIHEDRAL 3
#define FGA_GROUP_CYCLIC 4
/* #define FGA_GROUP_PGL_2_P 5 */
#define FGA_GROUP_M11 6
#define FGA_GROUP_M12 7
#define FGA_GROUP_M23 8
#define FGA_GROUP_M24 9
#define FGA_GROUP_CUBE 10
#define FGA_GROUP_PSL 11
#define FGA_GROUP_PGL 12
#define FGA_GROUP_PSSL 13
#define FGA_GROUP_PGGL 14
#define FGA_GROUP_PSU_3_Q2 15
#define FGA_GROUP_SZ_Q 16
#define FGA_GROUP_Zn_MULTIPLICATOR 17
#define FGA_GROUP_FROM_FILE 18
#define FGA_GROUP_ASL 19
#define FGA_GROUP_AGL 20
#define FGA_GROUP_ASSL 21
#define FGA_GROUP_AGGL 22
#define FGA_GROUP_YOUNG2 23
	// young groups of type [t,n-t] (2-rowed partitions)
#define FGA_GROUP_YOUNG3 24
	// not yet implemented
#define FGA_GROUP_SN_WREATH_SM 25
#define FGA_GROUP_HOLOMORPH_OF_CYCLIC_GROUP 26
#define FGA_GROUP_TRIVIAL 27
#define FGA_GROUP_BY_PERMUTATION_GENERATOR 28


#define FGA_MODE_ON_N 0
#define FGA_MODE_ON_N2 1
#define FGA_MODE_ON_N2_01 2
#define FGA_MODE_BY_CONJUGATION 3

#define FGA_GROUP_MATHIEU 30
#define FGA_GROUP_HIGMAN_SIMS_176 31
#define FGA_SOLVABLE_GROUP 32

#define FGA_DIRECT_SUM 40
#define FGA_DIRECT_PRODUCT 41
#define FGA_INDUCE_2_SETS 42
#define FGA_INDUCE_2_TUPLES 43
#define FGA_INDUCE_INJECTIVE_2_TUPLES 44
#define FGA_ADD_FIXPOINT 45
#define FGA_STABILIZE_POINT 46
#define FGA_WREATH_PRODUCT 47
#define FGA_EXPONENTIATION 48
#define FGA_HOLOMORPH 49
#define FGA_INDUCE_3_SETS 50
#define FGA_COMMA 51
#define FGA_ON_MAPPINGS 52
#define FGA_SELECT_ITH 53
#define FGA_SOLVABLE_GROUP_PERMUTATION_REPRESENTATION 54
#define FGA_ABSTRACT_PRESENTATION 55
#define FGA_SUBGROUP_OF_HOLOMORPH_OF_CYCLIC_GROUP 56

#define FGA_GROUP_SL 57
#define FGA_GROUP_GL 58
#define FGA_GROUP_SSL 59
#define FGA_GROUP_GGL 60

#define FGA_POWER 61
#define FGA_GROUP_AFFINE_TRANSLATIONS 62

#define FGA_TETRAHEDRON 63
#define FGA_CUBE 64
#define FGA_OCTAHEDRON 65
#define FGA_DODECAHEDRON 66
#define FGA_ICOSAHEDRON 67

#define FGA_SOLID_TRUNCATE 68
#define FGA_SOLID_DUAL 69
#define FGA_SOLID_TRUNCATE_DODE 70
#define FGA_SOLID_TRUNCATE_CUBE 71
#define FGA_CUBE_4D 72
#define FGA_SOLID_RELABEL_POINTS 73
#define FGA_SOLID_INDUCED_GROUP_ON_EDGES 74
#define FGA_SOLID_EDGE_MIDPOINTS 75
#define FGA_SOLID_ADD_CENTRAL_POINT 76
#define FGA_SOLID_ADD_CENTRAL_INVOLUTION 77
#define FGA_SOLID_CUBUS_SIMUS 78
#define FGA_SOLID_DODE_SIMUM 79
#define FGA_SOLID_CUBUS_EE 80
#define FGA_SOLID_CUBUS_EE_RUSSIAN 81



class group_selection_ob : public VECTOR_OB {
public:
	INTEGER_OP s_type() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_type_i() { 
		return(s_type()->s_i()); };
	INTEGER_OP s_val1() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_val1_i() { 
		return(s_val1()->s_i()); };
	INTEGER_OP s_val2() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_val2_i() { 
		return(s_val2()->s_i()); };
	STRING_OP s_s() {
		return (STRING_OP) s_i(3); };
	BYTE *s_s_s() {
		return s_s()->s_str(); };

INT field_name(INT i, INT j, BYTE *str);
INT init(INT type, INT val1, INT val2, BYTE *s);
INT sprint(BYTE *s);
};

INT orbits(VECTOR_OP G, 
	VECTOR_OP SVorbit, VECTOR_OP SVlast, 
	VECTOR_OP SVgen, 
	VECTOR_OP Ofirst, VECTOR_OP Osize);
INT trace_schreier_vectors(INT i, 
	VECTOR_OP SVlast, VECTOR_OP SVgen, 
	VECTOR_OP G, SYM_OP p);
INT schreier_stabilizer(VECTOR_OP G, 
	VECTOR_OP SVorbit, VECTOR_OP SVlast, 
	VECTOR_OP SVgen, 
	INT which_orbit, VECTOR_OP stab);
INT transitivity(VECTOR_OP G, INT *t, 
	INT *f_sharp, INT *f_point_five, INT f_verbose);



INT report_group_latex_stdout(VECTOR_OP G_gen);
INT perm_vec_get_degree(VECTOR_OP V);
void go(VECTOR_OP p);
INT gen_reduce(
	VECTOR_OP gen, SYM_OP group_order, 
	INT f_verbose, FILE *fp);
INT reduce_generators(VECTOR_OP p, SYM_OP group_order, INT f_verbose);
INT reduce_generators_labra(VECTOR_OP p, SYM_OP group_order, 
	INT f_verbose, LABRA_OP L);
INT group_elements(VECTOR_OP G_gen, VECTOR_OP G);
INT dimino_ex(VECTOR_OP G, VECTOR_OP G_gen, SYM_OP g);

INT km_get_group_generators1(INT type, INT val1, INT val2, BYTE *s, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len);
INT km_get_group_generators2(INT type, INT val1, INT val2, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len);
INT km_get_group_solid(INT type, INT val1, INT val2, BYTE *s, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len);
INT km_get_group_generators_from_file(INT type, INT val1, INT val2, BYTE *fname, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len);
INT km_binary_operations(INT type, INT val1, INT val2, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len);
INT km_unary_operations(INT type, INT val1, INT val2, 
	BYTE *g_label, BYTE *g_label_tex, VECTOR_OP S, INT *S_len);
VECTOR_OP km_get_generators(VECTOR_OP S, INT S_len);
INT km_get_deg(VECTOR_OP S, INT S_len);
INT print_deg(VECTOR_OP S, INT S_len);
INT km_get_group_from_selection(VECTOR_OP gsel, INT nb_g_sel, 
	VECTOR_OP generators, BYTE *g_label, BYTE *g_label_tex);
INT compose_gsel_from_strings(VECTOR_OP V, INT num_args, BYTE **args);



#endif /* FGA_INCLUDED */

