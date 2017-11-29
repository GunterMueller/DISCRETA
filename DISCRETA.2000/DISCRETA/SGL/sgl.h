/* sgl.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef SGL_INCLUDED
#define SGL_INCLUDED

#ifndef VDI_INCLUDED
#include <DISCRETA/graphics.h>
#endif
#ifndef DIMINO_INCLUDED
#include <DISCRETA/dimino.h>
#endif

/* type 0: ordinary SYMMETRICA mult()
 * 1: reverse mult (PERMUTATIONs)
 * 2: ENUM mult */
/* sym.h, used in nu.C: */
/* #define DO_TYPE_SYM 0 */
/* #define DO_TYPE_PERM 1 */
/* #define DO_TYPE_ENUM 2 */
/* #define DO_TYPE_IPERM 3 */
/* #define DO_TYPE_FG 4 */

#define SGL_VBP_X_PIX 10000
#define SGL_VBP_Y_PIX 10000

#define SGL_RADIUS 90


class sgo_ob : public VECTOR_OB {
public:
	/* o_len = number of groups in this orbit 
	 * (the total number of groups, the orbit according to 
	 * outer conjugation)
	 * so_len = the orbit length under conjugation by inner automorphisms
	 */
	INTEGER_OP s_f_has_aut() {
		return((INTEGER_OP) s_i(0)); };
	INT s_f_has_aut_i() {
		return(s_f_has_aut()->s_i()); };
	INTEGER_OP s_o_len() {
		return((INTEGER_OP) s_i(1)); };
	INT s_o_len_i() {
		return(s_o_len()->s_i()); };
	INTEGER_OP s_so_len() {
		return((INTEGER_OP) s_i(2)); };
	INT s_so_len_i() {
		return(s_so_len()->s_i()); };
	INTEGER_OP s_go() { 
		return((INTEGER_OP) s_i(3)); };
	INT s_go_i() { 
		return(s_go()->s_i()); };
	
	/* Erzeugendensystem aus Zuppo - Indices: 
	 * diese sind eine Teilmenge von 
	 * s_Zidx_i(s_first_i()). */
	VECTOR_OP s_gen_Zidx() { 
		return((VECTOR_OP) s_i(4)); };
	INTEGER_OP s_gen_Zidx_i(INT i) { 
		return((INTEGER_OP) s_gen_Zidx()->s_i(i)); };
	INT s_gen_Zidx_ii(INT i) { 
		return(s_gen_Zidx_i(i)->s_i()); };

	/* first: Position der urspruenglichen Gruppe H 
	 * (welche von gen erzeugt ist) in Zidx */
	INTEGER_OP s_first() { 
		return((INTEGER_OP) s_i(5)); };
	INT s_first_i() { 
		return(s_first()->s_i()); };

	/* alle Indices von Zuppos, 
	 * welche in der Gruppe liegen: 
	 * Zusammengefasst alle konjugierten 
	 * Untergruppen in einem Vektor.
	 * i laeuft ueber die Konjugationsklasse 
	 *   von Untergruppen, 
	 * j ueber die Zuppos einer Gruppe.  */
	VECTOR_OP s_Zidx() { 
		return((VECTOR_OP) s_i(6)); };
	VECTOR_OP s_Zidx_i(INT i) { 
		return((VECTOR_OP) s_Zidx()->s_i(i)); };
	INTEGER_OP s_Zidx_ij(INT i, INT j) { 
		return((INTEGER_OP) s_Zidx_i(i)->s_i(j)); };
	INT s_Zidx_iji(INT i, INT j) { 
		return(s_Zidx_ij(i, j)->s_i()); };
	
	/* Schreier Vektoren, welche die 
	 * Bahn der (Unter-) Gruppe
	 * unter Konjugation beschreiben. 
	 * Die Generator - Indices 
	 * beziehen sich auf den verwendeten 
	 * Vektor von Generatoren von G. 
	 * SVlast muss beim sortierten Einfuegen neuer 
	 * Gruppen laufend aktualisiert werden. */
	
	 /* SVlast / SVgen
	  * Schreier vectors for the orbit(s) of inner automorphisms 
	  * on the groups; 
	  * eventually, there are more that one entries -1, i.e. 
	  * orbit representatives in SVlast. 
	  * The correspond to the representatives of classes 
	  * modulo inner automorphisms.
	  * The first -1 is the representative of the orbit 
	  * under outer automorphisms. This position is 
	  * also stored in 'first' (see above). */
	VECTOR_OP s_SVlast() {
		return((VECTOR_OP) s_i(7)); };
	INTEGER_OP s_SVlast_i(INT i) {
		return((INTEGER_OP) s_SVlast()->s_i(i)); };
	INT s_SVlast_ii(INT i) {
		return((INT) s_SVlast_i(i)->s_i()); };
	VECTOR_OP s_SVgen() {
		return((VECTOR_OP) s_i(8)); };
	INTEGER_OP s_SVgen_i(INT i) {
		return((INTEGER_OP) s_SVgen()->s_i(i)); };
	INT s_SVgen_ii(INT i) {
		return((INT) s_SVgen_i(i)->s_i()); };

	/* Aut_SVlast / Aut_SVgen
	 * Schreier vectors (of length o_len) 
	 * for the orbit under outer automorphisms 
	 */
	VECTOR_OP s_Aut_SVlast() {
		return((VECTOR_OP) s_i(9)); };
	INTEGER_OP s_Aut_SVlast_i(INT i) {
		return((INTEGER_OP) s_Aut_SVlast()->s_i(i)); };
	INT s_Aut_SVlast_ii(INT i) {
		return((INT) s_Aut_SVlast_i(i)->s_i()); };
	VECTOR_OP s_Aut_SVgen() {
		return((VECTOR_OP) s_i(10)); };
	INTEGER_OP s_Aut_SVgen_i(INT i) {
		return((INTEGER_OP) s_Aut_SVgen()->s_i(i)); };
	INT s_Aut_SVgen_ii(INT i) {
		return((INT) s_Aut_SVgen_i(i)->s_i()); };
	
	PERMUTATION_OP s_p() {
		return((PERMUTATION_OP) s_i(11)); };
	PERMUTATION_OP s_pv() {
		return((PERMUTATION_OP) s_i(12)); };
	
	/* the generators of G acting 
	 * by conjugation
	 * on the groups of this orbit 
	 * (degree = orbit length): */
	VECTOR_OP s_G_gen_on_orbit() { 
		return((VECTOR_OP) s_i(13)); };
	PERMUTATION_OP s_G_gen_on_orbit_i(INT i) { 
		return((PERMUTATION_OP) 
			s_G_gen_on_orbit()->s_i(i)); };
	VECTOR_OP s_Aut_gen_on_orbit() { 
		return((VECTOR_OP) s_i(14)); };
	PERMUTATION_OP s_Aut_gen_on_orbit_i(INT i) { 
		return((PERMUTATION_OP) 
			s_Aut_gen_on_orbit()->s_i(i)); };

	INTEGER_OP s_Nlayer() { 
		return((INTEGER_OP) s_i(15)); };
	INT s_Nlayer_i() { 
		return(s_Nlayer()->s_i()); };
	INTEGER_OP s_Norbit() { 
		return((INTEGER_OP) s_i(16)); };
	INT s_Norbit_i() { 
		return(s_Norbit()->s_i()); };
		
	/* die Normalisatoren der Gruppen des Orbits 
	 * im Orbit (Nlayer / Norbit): */
	VECTOR_OP s_Nrep() {
		return((VECTOR_OP) s_i(17)); };
	INTEGER_OP s_Nrep_i(INT i) {
		return((INTEGER_OP) s_Nrep()->s_i(i)); };
	INT s_Nrep_ii(INT i) {
		return((INT) s_Nrep_i(i)->s_i()); };
	

	INT field_name(INT i, INT j, BYTE *str);
	INT Print(SGL_OP L, 
		INT f_zuppos_expanded, 
		INT type, void *data);
	INT Init(
		VECTOR_OP H_gen, VECTOR_OP H, 
		VECTOR_OP N_gen, VECTOR_OP N, 
		SGL_OP L, INT f_verbose, 
		INT type, void *data);
	INT calc_splitting_orbits(SGL_OP L, 
		VECTOR_OP N, VECTOR_OP N_gen, 
		INT f_verbose, INT f_very_verbose, 
		INT type, void *data);
	INT calc_orbit_Aut_SV(SGL_OP L, VECTOR_OP H, 
		INT f_verbose, INT f_very_verbose);
	INT calc_nreps(SGO_OP N_orbit, INT nrep0, INT f_verbose);
	INT recalc_Zidx(SGL_OP L, INT type, void *data, INT f_v);
	INT calc_group_table(SGL_OP L, INT rep, 
		MATRIX_OP T, VECTOR_OP generators, 
		VECTOR_OP embedding, VECTOR_OP embedding_inv, 
		INT type, void *data);
	INT calc_group_elements(SGL_OP L, INT rep, 
		VECTOR_OP ge, VECTOR_OP generators, INT type, void *data);
	INT calc_Zidx(SGL_OP L, INT rep, 
		VECTOR_OP Zidx, INT type, void *data);
	INT calc_gen_Zidx(SGL_OP L, INT rep, VECTOR_OP gen_Zidx);
	INT gen_on_orbit(SGL_OP L, 
		VECTOR_OP gen_on_zuppos, VECTOR_OP gen_on_orbit, INT f_verbose);
	INT find_group(VECTOR_OP Zidx, INT *idx, INT *f_found);
	INT trace_on(INT f_Aut_SV, INT i, SYM_OP g, 
	VECTOR_OP Gen, INT type, void *data);
	INT trace_G_gen(SGL_OP L, INT i, SYM_OP g, INT type, void *data);
	INT trace_G_gen_on_zuppos(SGL_OP L, INT i, PERMUTATION_OP g);
	INT trace_Aut_gen_on_zuppos(SGL_OP L, INT i, PERMUTATION_OP g);
	INT generators(SGL_OP L, VECTOR_OP gen, INT type, void *data);
};

class sgl_ob : public VECTOR_OB {
public:
	/* nb_zuppos, nb_PP entfallen */
	/* Z, Zidx, Zorder, Zprime, Zk, PP, PPidx 
	 * berechnet von zuppos() */
	VECTOR_OP s_Z() { 
		return((VECTOR_OP) s_i(0)); };
	SYM_OP s_Z_i(INT i) { 
		return(s_Z()->s_i(i)); };
	VECTOR_OP s_Zidx() { 
		return((VECTOR_OP) s_i(1)); };
	INTEGER_OP s_Zidx_i(INT i) { 
		return((INTEGER_OP) s_Zidx()->s_i(i)); };
	INT s_Zidx_ii(INT i) { 
		return(s_Zidx_i(i)->s_i()); };
	VECTOR_OP s_Zorder() { 
		return((VECTOR_OP) s_i(2)); };
	INTEGER_OP s_Zorder_i(INT i) { 
		return((INTEGER_OP) s_Zorder()->s_i(i)); };
	INT s_Zorder_ii(INT i) { 
		return(s_Zorder_i(i)->s_i()); };
	VECTOR_OP s_Zprime() { 
		return((VECTOR_OP) s_i(3)); };
	INTEGER_OP s_Zprime_i(INT i) { 
		return((INTEGER_OP) s_Zprime()->s_i(i)); };
	INT s_Zprime_ii(INT i) { 
		return(s_Zprime_i(i)->s_i()); };
	VECTOR_OP s_Zk() { 
		return((VECTOR_OP) s_i(4)); };
	INTEGER_OP s_Zk_i(INT i) { 
		return((INTEGER_OP) s_Zk()->s_i(i)); };
	INT s_Zk_ii(INT i) { 
		return(s_Zk_i(i)->s_i()); };
	VECTOR_OP s_PP() { 
		return((VECTOR_OP) s_i(5)); };
	SYM_OP s_PP_i(INT i) { 
		return(s_PP()->s_i(i)); };
	VECTOR_OP s_PPidx() { 
		return((VECTOR_OP) s_i(6)); };
	INTEGER_OP s_PPidx_i(INT i) { 
		return((INTEGER_OP) s_PPidx()->s_i(i)); };
	INT s_PPidx_ii(INT i) { 
		return(s_PPidx_i(i)->s_i()); };

	VECTOR_OP s_zuppos_on_zuppos() { 
		return((VECTOR_OP) s_i(7)); };
	PERMUTATION_OP s_zuppos_on_zuppos_i(INT i) { 
		return((PERMUTATION_OP) 
			s_zuppos_on_zuppos()->s_i(i)); };
	
	INTEGER_OP s_nb_Layers() { 
		return((INTEGER_OP) s_i(8)); };
	INT s_nb_Layers_i() { 
		return(s_nb_Layers()->s_i()); };
	VECTOR_OP s_theOrbits() { 
		return((VECTOR_OP) s_i(9)); };
	VECTOR_OP s_theOrbits_i(INT i) { 
		return((VECTOR_OP) 
			s_theOrbits()->s_i(i)); };
	SGO_OP s_theOrbits_ij(INT i, INT j) { 
		return((SGO_OP) 
			s_theOrbits_i(i)->s_i(j)); };
	/* nb_orbits VECTOR ueber INT entfaellt */
	VECTOR_OP s_nb_subgroups() { 
		return((VECTOR_OP) s_i(10)); };
	INTEGER_OP s_nb_subgroups_i(INT i) { 
		return((INTEGER_OP) 
			s_nb_subgroups()->s_i(i)); };
	INT s_nb_subgroups_ii(INT i) { 
		return(s_nb_subgroups_i(i)->s_i()); };
	
	INTEGER_OP s_total_nb_orbits() { 
		return((INTEGER_OP) s_i(11)); };
	INT s_total_nb_orbits_i() { 
		return(s_total_nb_orbits()->s_i()); };
	INTEGER_OP s_total_nb_subgroups() { 
		return((INTEGER_OP) s_i(12)); };
	INT s_total_nb_subgroups_i() { 
		return(s_total_nb_subgroups()->s_i()); };
	
	/* the generators of the group; 
	 * then as permutations on the zuppos */
	VECTOR_OP s_G_gen() { 
		return((VECTOR_OP) s_i(13)); };
	SYM_OP s_G_gen_i(INT i) { 
		return(s_G_gen()->s_i(i)); };
	VECTOR_OP s_G_gen_on_zuppos() { 
		return((VECTOR_OP) s_i(14)); };
	PERMUTATION_OP s_G_gen_on_zuppos_i(INT i) { 
		return((PERMUTATION_OP) 
			s_G_gen_on_zuppos()->s_i(i)); };
	
	/* generators for the automorphism group:
	 * as permutations of the group elements;
	 * the group elements have to be integers 
	 * in 0...group_order - 1
	 */
	INTEGER_OP s_f_has_aut_group() {
		return((INTEGER_OP) s_i(15)); };
	INT s_f_has_aut_group_i() {
		return(s_f_has_aut_group()->s_i()); };
	VECTOR_OP s_Aut_gen() { 
		return((VECTOR_OP) s_i(16)); };
	PERMUTATION_OP s_Aut_gen_i(INT i) { 
		return((PERMUTATION_OP) s_Aut_gen()->s_i(i)); };
	VECTOR_OP s_Aut_gen_on_zuppos() { 
		return((VECTOR_OP) s_i(17)); };
	PERMUTATION_OP s_Aut_gen_on_zuppos_i(INT i) { 
		return((PERMUTATION_OP) s_Aut_gen_on_zuppos()->s_i(i)); };

	/* group order: */
	INTEGER_OP s_go() { 
		return((INTEGER_OP) s_i(18)); };
	INT s_go_i() { 
		return(s_go()->s_i()); };
	
	INTEGER_OP s_f_verbose() { 
		return((INTEGER_OP) s_i(19)); };
	INT s_f_verbose_i() { 
		return(s_f_verbose()->s_i()); };
	INTEGER_OP s_f_very_verbose() { 
		return((INTEGER_OP) s_i(20)); };
	INT s_f_very_verbose_i() { 
		return(s_f_very_verbose()->s_i()); };

	/* sgl.C: */
	INT field_name(INT i, INT j, BYTE *str);
	INT Init_with_Aut_generators(
		VECTOR_OP G_gen, VECTOR_OP G, VECTOR_OP Aut_gen, 
		INT f_verbose, INT f_very_verbose, 
		INT type, void *data);
	INT Init(
		VECTOR_OP G_gen, VECTOR_OP G, 
		INT f_verbose, INT f_very_verbose, 
		INT type, void *data);
	INT Init1(VECTOR_OP G, INT f_verbose, INT f_very_verbose, INT type, void *data);
	INT Init2(VECTOR_OP G_gen, INT type, void *data);
	INT add_trivial_subgroups(VECTOR_OP G_gen, VECTOR_OP G, INT f_verbose, INT type, void *data);
	INT PrintOrbits(INT f_zuppos_expanded, INT type, void *data);
	INT PrintOrbit(INT layer, INT orbit, INT f_zuppos_expanded, INT type, void *data);
	INT PrintStatistics();
	INT sprint_statistics(BYTE *str1, BYTE *str2, BYTE *str3);
	INT sprint_statistics_tex(BYTE *str1, BYTE *str2, BYTE *str3);
	INT get_orbit_generators(VECTOR_OP Gen, INT type, void *data);
	INT add_group(VECTOR_OP U_gen, VECTOR_OP U, INT f_verbose, INT type, void *data);
	INT all_layers(INT f_verbose, INT type, void *data);
	INT extend_layer(INT layer, INT f_verbose, INT type, void *data);
	INT extend_group(INT layer, INT orbit, VECTOR_OP Gamma, INT f_verbose, INT type, void *data);
	INT extend_init_gamma(INT layer, INT orbit, VECTOR_OP Gamma);
	INT sylow_normalizing_zuppos( INT layer, INT orbit, 
		INT U_go, VECTOR_OP Gamma, INT f_verbose, INT type, void *data);
	INT find_normalizing_p_zuppo(VECTOR_OP U_Zidx, VECTOR_OP U, 
		VECTOR_OP G_Zidx, INT p, INT *idx, INT f_verbose, INT type, void *data);
	INT find_max_p_element(VECTOR_OP Z_idx, INT p, INT ny_p0, 
		INT *idx, INT *ny_p1, INT type, void *data);
	INT find_group(VECTOR_OP Zidx, INT layer, INT *orbit, INT *rep, INT *f_found);
	INT representation_on_zuppos(SYM_OP g, PERMUTATION_OP p, INT type, void *data);
	INT Aut_representation_on_zuppos(VECTOR_OP G, PERMUTATION_OP aut, PERMUTATION_OP p);
	INT conjugate_zuppos_by_perm(VECTOR_OP z_idx, PERMUTATION_OP g, VECTOR_OP z_idx1);
	INT conjugate_z_idx(VECTOR_OP z_idx, SYM_OP g, VECTOR_OP z_idx1, INT type, void *data);
	INT Aut_conjugate_z_idx(VECTOR_OP z_idx, PERMUTATION_OP aut, VECTOR_OP z_idx1);
	
	/* sgls2.C: */
	INT zuppos(VECTOR_OP G, INT f_verbose, INT type, void *data);
	INT zuppos_on_zuppos(INT type, void *data);
	INT calc_G_gen_on_zuppos(INT type, void *data);
	INT calc_aut_on_zuppos(VECTOR_OP G);
	INT latex_Z_info(FILE *fp, BYTE *tex_group_name, INT type, void *data);
	INT print_zuppo_vector(VECTOR_OP V, INT f_vertically, INT type, void *data);
	INT calc_zidx(VECTOR_OP H, VECTOR_OP H_zidx);
	INT Z2Zidx_repetitions_allowed(VECTOR_OP V, VECTOR_OP V_Zidx, INT type, void *data);
	INT Z2Zidx(VECTOR_OP V, VECTOR_OP V_Zidx, INT type, void *data);
	INT Zidx2Z(VECTOR_OP V_Zidx, VECTOR_OP V);
	INT IndexInNormalizer(INT layer, INT orbit, INT *idx_in_normalizer);
	INT orbit2lo(INT orbit, INT *l, INT *o);
	INT lo2orbit(INT l, INT o, INT *orbit);
	INT IsSubgroup(
		INT layer1, INT orbit1, INT rep1, 
		INT layer2, INT orbit2, INT rep2, 
		INT *f_is_subgroup);
	INT IsSubgroup_with_recalc(
		INT layer1, INT orbit1, INT rep1, 
		INT layer2, INT orbit2, INT rep2, 
		INT *f_is_subgroup, INT type, void *data);
	INT Asup(MATRIX_OP A);
	INT Ni(MATRIX_OP D);
	INT OrbitSizes(VECTOR_OP v);
	INT shrink();
	INT grow(INT type, void *data, INT f_v);
	INT Burnside_info(VECTOR_OP orbit_size, 
		MATRIX_OP As, MATRIX_OP Ai, MATRIX_OP D, 
		MATRIX_OP M, MATRIX_OP B, SYM_OP d, INT f_v);
};

// sgll:

typedef struct sg_lattice_local SG_LATTICE_LOCAL;

struct sg_lattice_local {
	MEM_OP plaz;
	MEM_OP o_dx;
	double extrema[6];
	SGL_OP L;
	INT f_L_allocated;
	MATRIX_OP Asup;
	MATRIX_OP Ainf;
	MATRIX_OP D;
	MATRIX_OP M;
	MATRIX_OP B;
	MATRIX_OP Acover;
	VECTOR_OP nl;
	VECTOR_OP orbit_size;
	SYM_OP d; // = det M
	INT f_show_generators;
	/* INT f_enumerated_group; (entfaellt) */
	INT f_with_perm;
	INT draw_lines_type;
		/* 0 = Asup, 1 = Ainf, 2 = all */
	BYTE path_name[1024];
	BYTE group_name[64];
	BYTE tex_group_name[64];
	VDEVICE vdev;
	
	/* only used inside fg_tape.C 
	 * as fg_tape_draw_data: */
	INT f_has_fg;
	FG_OP fg;
	SGL_OP sgl;
	INT lines3;
	VECTOR_OP bottom_text;

	VECTOR_OP subgroups; /* computed via fg_sgl.C */
	double extrema2[6];
	double extrema3[6];
	double extrema4[6];
	INT font_size3;
	INT font_size4;
	double line_height3;
	double line_height4;
	INT lines4;
	INT ago_n;
	INT ago_m;

	VECTOR_OP generators;
};


INT calc_sgl(VECTOR_OP G_gen, VECTOR_OP G, VECTOR_OP Inn_gen, SGL_OP L);
INT get_Inn_G(VECTOR_OP G_gen, VECTOR_OP G, VECTOR_OP Aut_gen);
SG_LATTICE_LOCAL *open_sgll(SGL_OP L, 
	INT f_L_allocated, 
	INT f_show_generators, 
	INT f_with_perm, INT draw_lines_type, 
	INT f_draw_sgl, INT f_burnside_tex, 
	BYTE *path_name, BYTE *group_name, 
	BYTE *tex_group_name, INT f_verbose);
INT free_sgll(SG_LATTICE_LOCAL *p);
INT sgll_output_layers(SG_LATTICE_LOCAL *sgll, FILE *fp);
INT sgll_init_generators(SG_LATTICE_LOCAL *sgll);

/* sgld.C: */

extern INT sgll_draw_sims;

INT sgll_draw_labra(void *sgll, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev);
INT sgll_draw_boxed_func(void *sgll, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev);
INT sgll_draw_func(void *sgll, 
	INT x, INT y, INT w, INT h, VDEVICE *vdev);
INT sgll_d_orbit(SG_LATTICE_LOCAL *sgll, INT orbit, INT rad, VDEVICE *vdev);
INT sgll_d_neighbour(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT neighbour, VDEVICE *vdev);
INT sgll_d_neighbour2(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT neighbour, 
	INT l1, INT o1, INT rep1, INT l2, INT o2, INT rep2, 
	VDEVICE *vdev);
INT sgll_get_ko_center(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT *pix_x1, INT *pix_y1, VDEVICE *vdev);
INT sgll_get_ko(SG_LATTICE_LOCAL *sgll, 
	INT orbit, INT rep, 
	INT *pix_x1, INT *pix_y1, VDEVICE *vdev);
INT sgll_draw_dot(VDEVICE *vdev, INT *x, INT *y, INT i, INT rad);

#endif /* SGL_INCLUDED */

