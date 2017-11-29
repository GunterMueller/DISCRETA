// disc.h


#include <stdio.h>
#include <stdlib.h>

#include <Xm/List.h>
#include <Xm/ToggleB.h>
#include <Xm/Text.h>
#include <Xm/TextF.h>


#include <DISCRETA/discreta.h>
#include <DISCRETA/fga.h> // for GROUP_SELECTION and FGA_GROUP_...
#include <DISCRETA/ladder.h>
#include <DISCRETA/gfq.h>
#include <DISCRETA/generators.h>

#ifdef SYM_GEO
#include <DISCRETA/geo.h>
#endif

#include <DISCRETA/db.h>

#include <DISCRETA/sgl.h>

#define BUFSIZE 16000

#define MAX_GROUP 50
/* #define RB_GROUP_DEFAULT 12 */
/* #define RB_M_DEFAULT 3
#define RB_LIN_DEFAULT 3
#define RB_SPORADIC_DEFAULT 0
#define RB_OTHER_DEFAULT 1 */

#define LOC_DESPAR_OUT "discreta_despar"
#define LOC_DB_CREATE "discreta_dbcreate.out"
#define LOC_DB_INSERT "discreta_dbinsert.out"
#define LOC_DB_SELECT "discreta_dbselect.out"
#define LOC_DB_DELETE "discreta_dbremove.out"

#define MAX_PG 20
#define MAX_PD 20
#define MAX_RB 20
#define DIAL_FRACTION_BASE 100

#define MAX_ITEMS 500
	// max nb of items in one listbox


// options for calc.C:
#define DESIGNS_LOG_FNAME "designs.txt"
#define DATA_PATH "DATA"

#undef WRITE_INTO_DATA_PATH




typedef struct discreta_form_data DISCRETA_FORM_DATA;

extern DISCRETA_FORM_DATA *form_data;
extern Widget Discreta_scrolled_text;

#define STRING_SIZE  1024

struct discreta_form_data {
	
	// group dialog:
	BYTE group_label[STRING_SIZE];
	INT v;
	BYTE group_order[STRING_SIZE];
	
	// parameter dialog:
	INT lambda, p;
	
	// compute KM-matrices dialog:
	INT t, k;
	BYTE KM_fname[STRING_SIZE];
	INT f_strong_generators;
	INT f_orderly_generation;
	INT f_TDO;
	INT f_extension_construction;
	INT f_canonical_representatives;
	INT f_k_k2_design;
	INT k2;
	INT f_graphical_design;
	INT f_graphical_design_Sn;
	INT graphical_design_n;
	INT graphical_nrow;
	INT graphical_ncol;
	
	INT toggle_lines;

	// gsel dialog:
	INT f_choosing_normalizer;
	BYTE gsel_file[STRING_SIZE];
	BYTE gsel_generator[STRING_SIZE];
	
	INT i_th;
	
	// Widget g_form;
	// Widget g_box;
	XmString g_Items[MAX_ITEMS];
	VECTOR_OP g_sel;
	// GROUP_SELECTION g_sel[MAX_ITEMS];
	INT g_nb_items;
	
	// linear group:
	INT linear_n, linear_q;
	INT linear_toggle;

	// well known group;
	INT well_n, well_m;
	INT well_rb_state;

	// solvable group:
	INT solvable_n;
	INT solvable_m;
	INT solvable_toggle;
	
	// regular solid:
	INT solid_toggle;
	
	// parameter:
	INT LLL_c0, LLL_beta, LLL_p, silence_level, nb_restart;
	INT largeset_max_rounds, largeset_N, largeset_max_solutions, largeset_max_loops;

	// report options:
	INT report_select_lambda;
	INT report_select_first;
	INT report_select_len;
	
	// sporadic group:
	INT sporadic_toggle;

	// database:
	INT t_min, t_max;
	INT v_min, v_max;
	INT k_min, k_max;
	INT lambda_min, lambda_max;
	
};

// discreta_stubs.C:

void set_group_label(BYTE *s);
void set_v(INT v);
void set_group_order(BYTE *s);
void text_field_set_string(Widget w, char *p);
void text_field_set_int(Widget w, INT v);
void text_field_get_string(Widget w, char *p, INT p_size);
INT text_field_get_int(Widget w);
INT get_toggle_state(Widget w);
INT set_toggle_state(Widget w, INT v);
void get_form_data(DISCRETA_FORM_DATA *p);
void alloc_form();
void dealloc_form();
INT decode_form();
void init_discreta_text();
void out_str(char *s);
void program_init();

Widget widget_sgl_type(INT i);



// dial_g.C:

INT g_box_free(DISCRETA_FORM_DATA *p);
INT g_box_refresh(DISCRETA_FORM_DATA *p);
INT get_group_from_selection(VECTOR_OP gsel, INT nb_g_sel, 
	VECTOR_OP generators, BYTE *g_label, BYTE *g_label_tex);
void group_selection_clear(DISCRETA_FORM_DATA *p);
void group_selection_delete(DISCRETA_FORM_DATA *p);
void group_selection_add(DISCRETA_FORM_DATA *p, INT type, INT val1, INT val2, BYTE *s);

// dial_well.C:
Widget widget_well_known(INT i);
void choose_well_known_group();

// dial_solvable.C:
Widget widget_solvable(INT i);
void choose_solvable_group();

// dial_linear.C:
Widget widget_linear(INT i);
void choose_linear_group();

// dial_sporadic.C:
Widget widget_sporadic(INT i);
void choose_sporadic_group();

// dial_solid.C:
Widget widget_solid(INT i);
void choose_regular_solid();

// cb.C
void calc_lambda_i(INT t, INT v, INT k, INT lambda);
INT km_fname_available();
INT read_km_fname(BYTE *KM_fname);
INT solid_fname_available();
INT read_solid_fname(BYTE *fname);
void do_mendelsohn(INT t, INT v, INT k, INT lambda);
void do_koehler(INT t, INT v, INT k, INT lambda, INT m);
void show_KM_matrix(BYTE *KM_fname);
void cb_do_LLL(INT f_silent, BYTE *KM_fname, 
	INT c0_factor, INT beta, INT p, INT lambda, INT f_iterate, INT nb_iterate);
void cb_do_iso_test(BYTE *KM_fname, INT lambda);
void cb_do_iso_test_sylow(BYTE *KM_fname, INT p, INT lambda, INT fv);

// db.C

#define DP_PATH "."

typedef struct despar DESPAR;

struct despar {
	INT t, v, k, lambda, id;
};

void db_add(INT t, INT v, INT k, INT lambda);
void db_add_vec_data(VECTOR_OP vec_data);
INT despar_add(MATRIX_OP data);
INT despar_buf_flush(char *db_fname);
INT despar_db_add(char *db_fname, INT t, INT v, INT k, INT lambda, INT id);
INT db_add(DESPAR *dp, INT len, char *db_fname);
INT despar_search(
	INT t_min, INT t_max, 
	INT v_min, INT v_max, 
	INT k_min, INT k_max, 
	INT lambda_min, INT lambda_max);
INT despar_delete(
	INT t_min, INT t_max, 
	INT v_min, INT v_max, 
	INT k_min, INT k_max, 
	INT lambda_min, INT lambda_max);
INT despar_db_get(char *sdf_fname, char *out_fname);
INT despar_db_delete(char *fname);
INT despar_read(BYTE *file_name, DESPAR **p_dp, INT *len);
INT despar2dbf(DESPAR *dp, INT len, BYTE **p_data, INT *p_data_len);
INT dbf2despar(char *pp, DESPAR *dp, INT *rec_len);
void despar_db_create(char *fname);
void despar_print(DESPAR *dp);


// tvt.C
void db_dp_tvt_closure();


