/* graphics.h */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#ifndef GRAPHICS_INCLUDED
#define GRAPHICS_INCLUDED

#ifdef GRAPHICS_TRUE

#ifndef MA_INCLUDED
#include <DISCRETA/ma.h>
#endif

typedef struct ged_local GED_LOCAL;

typedef struct vdevice VDEVICE;


#define GREY_INTERIOR 90
#define GREY_COLOR 0
#define GREY1_INTERIOR 85
#define GREY1_COLOR 0
#define GREY2_INTERIOR 90
#define GREY2_COLOR 0
#define GREY3_INTERIOR 95
#define GREY3_COLOR 0

#define WHITE_INTERIOR 100
#define WHITE_COLOR 0

#define BLACK_INTERIOR 0
#define BLACK_COLOR 0


/*
 * contrl[0]: main opcode
 * contrl[1]: number of points in ptsin
 * contrl[2]:
 * contrl[3]: number of entries in intin
 * contrl[4]:
 * contrl[5]: minor opcode
 */

#define Vi_intin(p, ptr)  (p)->vdiblk[1] = (INT *)(ptr)
#define Vi_intout(p, ptr) (p)->vdiblk[3] = (INT *)(ptr)
#define Vi_ptsin(p, ptr)  (p)->vdiblk[2] = (INT *)(ptr)
#define Vi_ptsout(p, ptr) (p)->vdiblk[4] = (INT *)(ptr)

#define Vi_ptr(p, addr) (*(INT *)(&(p)->contrl[7]) = (INT)addr)
#define Vi_ptr2(p, addr) (*(INT *)(&(p)->contrl[9]) = (INT)addr)
#define Vm_lptr2(p, addr) (*(INT *)addr = *(INT *)(&(p)->contrl[9]))

#ifndef M_EPS
#define M_EPS 0.000001
#endif
#ifndef M_PI
#define M_PI 3.1415927
#endif
#ifndef M_2PI
#define M_2PI 6.283185307
#endif
#ifndef M_PI_2
#define M_PI_2 (M_PI / 2.)
#endif
#ifndef M_SQRT2INV
#define M_SQRT2INV 0.707106781
#endif

struct vdevice {
	INT type;
		/* 0 = PsDraw() file: ps.cp */
	/* INT (*vdi_draw) (VDEVICE *vdev); */
	INT *contrl;
	INT *intin;
	INT *ptsin;
	INT *intout;
	INT *ptsout;
	INT *vdiblk[5];
	INT extra[5];
	/* UI_XWINDOWS: 
	 * extra[0] - Display pointer
	 * extra[1] - drawable 
	 * extra[2] - gc pointer
	 * alles gesetzt von MACLIB/WndDraw(), verwendet von x11_draw() 
	 * UI_SUNVIEW:
	 * extra[0] - Pixwin *pw
	 * gesetzt von MACLIB/draw_canvas(), verwendet von sv_draw() 
	 * VDEVICE - META:
	 * extra[0] - META pointer 
	 * gesetzt von MACLIB/opn_dev() geloescht von cls_dev(), 
	 * verwendet von VDEVICE/META/MetaDraw() */
	INT co[4]; /* llx, lly, urx, ury */
	double width;
	double height;
};

typedef struct lw_tapestry_contents LW_TAPESTRY_CONTENTS;
typedef INT (*TAPE_DBL_CLICK_FUNC)(LW_TAPESTRY_CONTENTS *p, INT active);

struct lw_tapestry_contents {
	void *data;
	INT len;
	INT width;
	INT picture_width;
	INT picture_height;
	INT (*calc_size) (LW_TAPESTRY_CONTENTS *p);
	INT (*draw_data_get) (LW_TAPESTRY_CONTENTS *p, 
		INT i, INT j, INT lines, INT columns, 
		INT first, void **draw_data);
	INT (*draw_data_free) (LW_TAPESTRY_CONTENTS *p, 
		INT i, INT j, INT lines, INT columns, void **draw_data);
	INT (*draw_picture) (LW_TAPESTRY_CONTENTS *p, 
		void *draw_data_this_pic, 
		INT m, INT n, INT x, INT y, INT w, INT h, VDEVICE *vdev);
	TAPE_DBL_CLICK_FUNC dbl_click;
	INT (*exit_data) (LW_TAPESTRY_CONTENTS *p);
	VDEVICE vdev;
};
/* Attribute definitions */

#define IP_HOLLOW 0
#define IP_1PATT 1
#define IP_2PATT 2
#define IP_3PATT 3
#define IP_4PATT 4
#define IP_5PATT 5
#define IP_6PATT 6
#define IP_SOLID 7

/* gsx modes */

#define MD_REPLACE 1
#define MD_TRANS 2
#define MD_XOR 3
#define MD_ERASE 4

/* gsx styles */

#define FIS_HOLLOW 0
#define FIS_SOLID 1
#define FIS_PATTERN 2
#define FIS_HATCH 3
#define FIS_USER 4

/* bit blt rules */

#define ALL_WHITE 0
#define S_AND_D 1
#define S_ONLY 3
#define NOTS_AND_D 4
#define S_XOR_D 6
#define S_OR_D 7
#define D_INVERT 10
#define NOTS_OR_D 13
#define ALL_BLACK 15

extern INT gl_vdi_contrl[12];
extern INT gl_vdi_intin[128];
extern INT gl_vdi_intout[128];
extern INT gl_vdi_ptsin[128];
extern INT gl_vdi_ptsout[128];

extern void *gl_vdi_draw_tbl[10];
/* extern INT (*gl_vdi_draw_tbl[10])(VDEVICE *p); */

void Vdi_draw(VDEVICE *p);
INT LwDTapestry2(LW_TAPESTRY_CONTENTS *p, VDEVICE *vdev, 
	INT x0, INT y0, INT w0, INT h0, 
	INT i, INT j, INT lines, INT columns, INT first, INT num);

INT Vswr_mode(VDEVICE *p, INT mode);
void Vs_color(VDEVICE *p, INT index, INT *rgb_in);
INT Vsl_type(VDEVICE *p, INT style);
void Vsl_udsty(VDEVICE *p, INT pattern);
INT Vsl_width(VDEVICE *p, INT width);
INT Vsl_color(VDEVICE *p, INT color_index);
void Vsl_ends(VDEVICE *p, INT beg_style, INT end_style);
INT Vsm_type(VDEVICE *p, INT symbol);
INT Vsm_height(VDEVICE *p, INT height);
INT Vsm_color(VDEVICE *p, INT color_index);
void Vst_alignment(VDEVICE *p, INT hor_in, 
	INT vert_in, INT *hor_out, INT *vert_out);
void Vst_height(VDEVICE *p, INT height, 
	INT *char_width, INT *char_height, 
	INT *cell_width, INT *cell_height);
INT Vst_point(VDEVICE *p, INT point, 
	INT *char_width, INT *char_height, 
	INT *cell_width, INT *cell_height); 
INT Vst_rotation(VDEVICE *p, INT angle);
INT Vst_font(VDEVICE *p, INT font);
INT Vst_color(VDEVICE *p, INT color_index);
INT Vst_effects(VDEVICE *p, INT effect);
INT Vsf_interior(VDEVICE *p, INT style);
INT Vsf_style(VDEVICE *p, INT style_index);
INT Vsf_color(VDEVICE *p, INT color_index);
INT Vsf_shape(VDEVICE *p, INT shape);
INT Vsf_outline(VDEVICE *p, INT outline);
INT Vsf_nofill(VDEVICE *p, INT nofill);
INT Vsf_perimeter(VDEVICE *p, INT per_vis);
// void Vsf_updat(VDEVICE *p, INT pfill_pat[], INT planes);
INT Vst_boxed(VDEVICE *p, INT boxed);
INT Vst_overwrite(VDEVICE *p, INT overwrite);

/* Output definitions */

void V_pline(VDEVICE *p, INT count, INT pxyarray[]);
void V_bezier(VDEVICE *p, INT count, INT pxyarray[]);
void V_pmarker(VDEVICE *p, INT count, INT pxyarray[]);
void V_gtext(VDEVICE *p, INT x, INT y, BYTE *string);
void V_fillarea(VDEVICE *p, INT count, INT pxyarray[]);
void V_cellarray(VDEVICE *p, INT pxyarray[4], 
	INT row_length, INT el_used, INT num_rows, 
	INT wrt_mode, INT colarray[]);
void V_contourfill(VDEVICE *p, INT x, INT y, 
	INT index);
void Vr_recfl(VDEVICE *p, INT *pxyarray);
void V_bar(VDEVICE *p, INT pxyarray[]);
void V_arc(VDEVICE *p, INT x, INT y, INT radius, 
	INT begang, INT endang);
void V_pie(VDEVICE *p, INT x, INT y, INT radius, 
	INT begang, INT endang);
void V_circle(VDEVICE *p, INT x, INT y, INT radius);
void V_ellipse(VDEVICE *p, INT x, INT y, 
	INT xradius, INT yradius);
void V_ellarc(VDEVICE *p, INT x, INT y, 
	INT xradius, INT yradius, INT begang, INT endang);
void V_ellpie(VDEVICE *p, INT x, INT y, 
	INT xradius, INT yradius, INT begang, INT endang);
void V_rbox(VDEVICE *p, INT xyarray[]);
void V_rfbox(VDEVICE *p, INT xyarray[]);
void V_justified(VDEVICE *p, INT x, INT y, BYTE *string, 
	INT length, INT word_space, INT char_space);
void V_pie_text(VDEVICE *p, INT x, INT y, INT radius, 
	INT begang, INT endang, BYTE *string);

INT VdevSExtra(VDEVICE *p, INT num, INT extra);
INT VdevGExtra(VDEVICE *p, INT num, INT *extra);
INT VdevI2(VDEVICE *p, INT *contrl, INT *intin, INT *ptsin, 
	INT *intout, INT *ptsout);
INT VdevGArrays(VDEVICE *p, INT **pcontrl, INT **pintin, 
	INT **pptsin, INT **pintout, INT **pptsout);
/* Es werden nur die arrays ausgelesen, deren Doppelpointer 
 * gesetzt (!=NIL) ist. */

INT user2vdev(VDEVICE *vdev, double dx, double dy, double dz, 
	INT *ix, INT *iy, double *extrema);
INT vdev2user(VDEVICE *vdev, INT ix, INT iy, double *dx, double *dy, 
	double *extrema);
INT user2xywh(double dx, double dy, double dz, 
	INT *ix, INT *iy, double *extrema, INT *xywh);
INT xywh2user(INT ix, INT iy, double *dx, double *dy, 
	double *extrema, INT *xywh);
INT place_lattice(VECTOR_OP nl, VECTOR_OP orbit_size, 
	INT size_x, INT size_y, 
	VECTOR_OP Px, VECTOR_OP Py, VECTOR_OP Os, 
	INT f_upside_down, INT f_v);
INT vbp(VECTOR_OP nl, VECTOR_OP orbit_size, 
	INT x_pix, INT y_pix, 
	MEM_OP plazierung, MEM_OP orbit_dx, 
	INT f_upside_down);
void draw_box(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1);
void draw_kreuz(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1);
void draw_kreuz2(VDEVICE *vdev, INT x, INT y, INT rad);
void draw_mol(VDEVICE *vdev, INT x, INT y, BYTE *text);
void bind_mol(VDEVICE *vdev, INT x1, INT y1, INT x2, INT y2, INT n);
void Bind_mol(VDEVICE *vdev, INT dist, INT x1, INT y1, INT x2, INT y2, INT n);
INT my_line(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1);
double cos_grad(double phi);
double sin_grad(double phi);
double tan_grad(double phi);
double atan_grad(double x);
INT draw_label(VDEVICE *vdev, BYTE *text, BYTE *align, INT *x, INT *y, INT idx);
INT draw_dot(VDEVICE *vdev, INT *x, INT *y, INT idx, INT rad);
INT draw_dot_white(VDEVICE *vdev, INT *x, INT *y, INT idx, INT rad);
INT draw_dot_black(VDEVICE *vdev, INT *x, INT *y, INT idx, INT rad);
INT draw_line(VDEVICE *vdev, INT *x, INT *y, INT idx1, INT idx2);
INT draw_bezier3(VDEVICE *vdev, INT *x, INT *y, INT i1, INT i2, INT i3);
INT draw_bezier4(VDEVICE *vdev, INT *x, INT *y, INT i1, INT i2, INT i3, INT i4);
INT draw_bezier5(VDEVICE *vdev, INT *x, INT *y, INT i1, INT i2, INT i3, INT i4, INT i5);
INT draw_bezier6(VDEVICE *vdev, INT *x, INT *y, 
	INT i1, INT i2, INT i3, INT i4, INT i5, INT i6);
INT draw_bezier7(VDEVICE *vdev, INT *x, INT *y, 
	INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7);
INT halve_way_int(INT *Px, INT *Py, INT i, INT j, INT k);
INT rotate_int(INT *Px, INT *Py, INT i, INT j, INT angle);
INT on_circle_int(INT *Px, INT *Py, INT idx, INT angle_in_degree, INT rad);
INT On_circle_int(VECTOR_OP Px, VECTOR_OP Py, INT idx, INT angle_in_degree, INT rad);
INT intersect_lines(INT *Px, INT *Py, 
	INT i1, INT i2, INT j1, INT j2, INT res);
INT intersect_lines_int(INT Px, INT Py, INT Qx, INT Qy, 
	INT Rx, INT Ry, INT Sx, INT Sy, INT *Ix, INT *Iy);
INT intersect_lines_P_u_Q_v(INT Px, INT Py, INT ux, INT uy, 
	INT Rx, INT Ry, INT vx, INT vy, INT *Ix, INT *Iy);
INT intersect_lines_P_u_Q_v_factor(INT Px, INT Py, INT ux, INT uy, 
	INT Rx, INT Ry, INT vx, INT vy, double *f);
INT intersect_line_area(VDEVICE *vdev, INT *x, INT *y, INT k, INT *points, 
	INT Px, INT Py, INT dx, INT dy, INT f_rim, INT rim_length);
INT line_with_whole(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1, INT rim_length);
INT lotfuss_int(INT *x, INT *y, INT i, INT a, INT b, INT *lot_x, INT *lot_y);
INT ratio_int(INT *x, INT *y, INT a, INT b, double r, INT *xx, INT *yy);
INT Ratio(INT *Px, INT *Py, INT *Pz, INT from, INT to, double r, INT idx);
INT stuetzgeraden_x(INT *x, INT *y, INT k, INT *points, INT phi, INT y0, 
	INT *p_xmin, INT *p_xmax);
INT stuetzgeraden_y(INT *x, INT *y, INT k, INT *points, INT phi, INT x0, 
	INT *p_ymin, INT *p_ymax);
INT build_path(VDEVICE *vdev, INT *x, INT *y, INT k, INT *points);
INT shade_area(VDEVICE *vdev, INT *x, INT *y, INT k, INT *points, INT phi, INT delta, 
	INT f_rim, INT rim_length);
INT translate(INT *Px, INT *Py, INT *Pz, INT dx, INT dy, INT dz, INT n);
INT rotate(INT *Px, INT *Py, INT *Pz, INT n);
INT random_isometry(double *A);
INT isometry(double *A, double nx, double ny, double nz, INT phi);
INT orthogonal(double nx, double ny, double nz, 
	double *ax, double *ay, double *az, 
	double *bx, double *by, double *bz);
double det3(double *A);
double det4(double a11, double a12, double a21, double a22);
double detAij(double *A, INT i, INT j);
INT transpose_mat(double *X, double *Y);
INT inverse_mat(double *X, double *Y);
INT mult_mat(double *A, double *B, double *C);
INT mult_mat_vec(double *A, double *v, double *v1);
INT print_mat33(double *A);
INT print_vec3(double *v);

// vbp.C:
int verband_plazieren(
	ULONG *nachfolgerliste, ULONG *orbit_size, 
	FLOAT **plazierung, FLOAT **orbit_dx);



/* ps.C: */
extern FILE *ps_draw_fp;
extern INT ps_draw_dev[4]; /* llx/lly/urx/ury */
extern INT ps_draw_ps[4]; /* llx/lly/urx/ury */
extern INT gl_ps_lines_per_page_tape;
extern INT gl_ps_cols_per_page_tape;
INT PsDraw(VDEVICE *vdev);
void PrintHeaderPS(FILE *fp, BYTE *fname);
void PrintFooterPS(FILE *fp);
INT ps_font_name(INT font_nr, BYTE str[256]);
INT draw_PS_file(void *data, VDEVICE *vdev, 
	INT pix_x_max, INT pix_y_max, BYTE *name, 
	INT (*draw_func)(void *gp_lattice, INT x, INT y, INT w, INT h, VDEVICE *vdev));
INT draw_tape_PS(LW_TAPESTRY_CONTENTS *p, BYTE *name);
INT draw_tape_PS_multiple_files(LW_TAPESTRY_CONTENTS *p, BYTE *name);

INT EpicDraw(VDEVICE *vdev);
INT draw_epic_file(void *data, VDEVICE *vdev, BYTE *name, 
	INT factor_1000, BYTE *caption, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev));
INT draw_tape_EPIC(LW_TAPESTRY_CONTENTS *p, 
	INT factor_1000, BYTE *name);

INT mp_draw(VDEVICE *vdev);
INT mp_set_output_coordinates(INT xmin, INT ymin, INT xmax, INT ymax);
INT mp_get_output_coordinates(INT *xmin, INT *ymin, INT *xmax, INT *ymax);
INT draw_mp_file(void *data, VDEVICE *vdev, BYTE *name, INT factor_1000, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev));
INT draw_mp_picture(FILE *fp, 
	void *data, VDEVICE *vdev, 
	INT xmin, INT ymin, INT xmax, INT ymax, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev));;
/* file output of a page;
 * output will be written into open file fp. 
 * mp_draw_dev[] has to be set 
 * before calling this routine. */
INT mp_header(FILE *fp, BYTE *fname);
INT mp_footer(FILE *fp);
INT mp_begin_figure(FILE *fp, INT factor_1000);
INT mp_end_figure(FILE *fp);
INT mp_begin_drawing(FILE *fp);
INT mp_end_drawing(FILE *fp);


/* tree.h */


class tree_ob : public VECTOR_OB {
public:
	INTEGER_OP s_depth() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_depth_i() { 
		return(s_depth()->s_i()); };
	
	INTEGER_OP s_nb_leaves() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_nb_leaves_i() { 
		return(s_nb_leaves()->s_i()); };
	
	INTEGER_OP s_nb_nodes() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_nb_nodes_i() { 
		return(s_nb_nodes()->s_i()); };


	// V is a vector of length nb_leaves
	// it is a vector of vectors, 
	// each vector containing a path to a leave. 
	// the leaves (or better: the paths leading 
	// to the leaves) are sorted. 
	// the last entry in each path-vector contains a 
	// string which labels that leave. 
	// so, the path length is v->s_li() - 1.
	VECTOR_OP s_V() { 
		return((VECTOR_OP)s_i(3)); };
	VECTOR_OP s_V_i(INT i) { 
		return((VECTOR_OP) s_V()->s_i(i)); };
	
	// # rows =  nb_nodes
	// first, len, for the suns of this node (indices into Tree),
	// first_V, len_V for the entries in V according to this subtree, 
	// node_label
	MATRIX_OP s_Tree() {
		return((MATRIX_OP)s_i(4)); };
	
	VECTOR_OP s_X() { 
		return((VECTOR_OP)s_i(5)); };
	INTEGER_OP s_X_i(INT i) { 
		return((INTEGER_OP)s_X()->s_i(i)); };
	INT s_X_ii(INT i) { 
		return(s_X_i(i)->s_i()); };
	VECTOR_OP s_Y() { 
		return((VECTOR_OP)s_i(6)); };
	INTEGER_OP s_Y_i(INT i) { 
		return((INTEGER_OP)s_Y()->s_i(i)); };
	INT s_Y_ii(INT i) { 
		return(s_Y_i(i)->s_i()); };

	INTEGER_OP s_max_x() { 
		return((INTEGER_OP)s_i(7)); };
	INT s_max_x_i() { 
		return(s_max_x()->s_i()); };
	INTEGER_OP s_max_y() { 
		return((INTEGER_OP)s_i(8)); };
	INT s_max_y_i() { 
		return(s_max_y()->s_i()); };

	INT field_name(INT i, INT j, BYTE *str);
	INT Init();
	INT sprint(BYTE *s);
	INT read_tree_from_file(INT f_verbose, BYTE *fname);
	INT build_tree(VECTOR_OP V, INT f_v);
	INT place(INT f_v);
	INT init_ged(GED_OP G, INT f_v);
	INT init_ged_labels(GED_OP G, INT f_v);
	INT init_ged_user_labels(GED_OP G, INT f_v);
	INT leaves_in_array();
};

/* ged.h */


#define KO_LENGTH 7
#define LINES_LENGTH 3

/* sel_idx, old_sel_idx: 
 * -1 means: no selection */

/* KO: */
/* j = 0: id, j = 1: x-KO, j = 2: y-KO 
 * 3: f_visible 4: string 
 * 5: font_nr 6: font_size */

// needed: 
// string_a, halign_a, valign_a
// string_b, halign_b, valign_b
// f_movable (f_fixed)




/* lines: */
/* j = 0: id, j = 1: 
 * id startpoint, j = 2: id endpoint */

// needed:
// groups (overlapping)
// a vector of vectors


class ged_ob : public VECTOR_OB {
public:
	INTEGER_OP s_max_x() { 
		return((INTEGER_OP)s_i(0)); };
	INT s_max_x_i() { 
		return(s_max_x()->s_i()); };
	INTEGER_OP s_max_y() { 
		return((INTEGER_OP)s_i(1)); };
	INT s_max_y_i() { 
		return(s_max_y()->s_i()); };
	INTEGER_OP s_cur_id() { 
		return((INTEGER_OP)s_i(2)); };
	INT s_cur_id_i() { 
		return(s_cur_id()->s_i()); };
	INTEGER_OP s_sel_idx() { 
		return((INTEGER_OP)s_i(3)); };
	INT s_sel_idx_i() { 
		return(s_sel_idx()->s_i()); };
	INTEGER_OP s_old_sel_idx() { 
		return((INTEGER_OP)s_i(4)); };
	INT s_old_sel_idx_i() { 
		return(s_old_sel_idx()->s_i()); };
	MATRIX_OP s_KO() { 
		return((MATRIX_OP)s_i(5)); };
	INTEGER_OP s_KO_ij(INT i, INT j) { 
		return((INTEGER_OP) 
			s_KO()->s_ij(i, j)); };
	INT s_KO_iji(INT i, INT j) { 
		return(s_KO_ij(i, j)->s_i()); };
		/* j = 0: id, j = 1: x-KO, j = 2: y-KO 
		 * 3: f_visible 4: string 
		 * 5: font_nr 6: font_size */
	MATRIX_OP s_lines() { 
		return((MATRIX_OP)s_i(6)); };
	INTEGER_OP s_lines_ij(INT i, INT j) { 
		return((INTEGER_OP) s_lines()->s_ij(i, j)); };
	INT s_lines_iji(INT i, INT j) { 
		return(s_lines_ij(i, j)->s_i()); };
		/* j = 0: id, j = 1: 
		 * id startpoint, j = 2: id endpoint */

INT read_graph_from_file(INT f_verbose, BYTE *fname);
INT field_name(INT i, INT j, BYTE *str);
INT Init();
INT sprint(BYTE *s);
INT click(GED_LOCAL *gl, VDEVICE *vdev, 
	INT x, INT y, INT click_count, 
	INT f_shift, INT f_control, 
	INT button, INT f_points, INT f_lines);
INT track(GED_LOCAL *gl, VDEVICE *vdev, 
	INT x0, INT y0, INT x1, INT y1, 
	INT f_shift, INT f_control, 
	INT button, INT f_points, INT f_lines, 
	INT *f_have_tracked);
INT track_smooth(
	GED_LOCAL *gl, VDEVICE *vdev, 
	INT x0, INT y0, INT x1, INT y1, 
	INT f_shift, INT f_control, 
	INT button, INT f_points, INT f_lines, 
	INT *f_have_tracked);
INT delete_KO_idx_full(INT i);
INT KO_edit_text(INT i);
INT new_selection(GED_LOCAL *gl, 
	VDEVICE *vdev, INT new_sel);
INT add_line_idx(INT idx1, INT idx2, INT *f_not_added, INT f_v);
INT add_line(INT id0, INT id1, INT *f_not_added, INT f_v);
INT del_line(INT i);
INT search_line_id0_id1(INT id0, INT id1, 
	INT *idx, INT *f_found);
INT add_KO(INT x, INT y, INT f_fisible);
INT init_KO(INT i, INT id);
INT init_KO_xy(INT i, INT x, INT y);
INT init_KO_string(INT i, BYTE *str);
INT del_KO(INT i);
INT search_KO_id(INT id, INT *idx, INT *f_found);
INT search_nearest(GED_LOCAL *gl, 
	INT x, INT y, INT *idx, INT *f_found);
INT draw(GED_LOCAL *gl, VDEVICE *vdev);
INT draw_KO_full(GED_LOCAL *gl, VDEVICE *vdev, INT i);
INT draw_ko(GED_LOCAL *gl, VDEVICE *vdev, INT i);
INT draw_ko_selection(GED_LOCAL *gl, 
	VDEVICE *vdev, INT i);
INT draw_line(GED_LOCAL *gl, VDEVICE *vdev, INT i);
};

struct ged_local {
	GED_OP theGED;
	double extrema[6];
};

#endif /* GRAPHICS_TRUE */

#endif /* GRAPHICS_INCLUDED */

