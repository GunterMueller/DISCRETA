// view.h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>

#include <Xm/Xm.h>
#include <X11/StringDefs.h>
#include <Xm/CascadeB.h>
#include <Xm/DialogS.h>
#include <Xm/FileSB.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/LabelG.h>
#include <Xm/List.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/ScrollBar.h>
#include <Xm/Separator.h>
#include <Xm/Text.h>
#include <Xm/TextF.h>
#include <Xm/ToggleB.h>

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include <GL/GLwMDrawA.h> 

// #include <GL/glut.h>

#undef VERBOSE

#define MINIMUM(x, y)   ( ((x) < (y)) ?  (x) : (y) )
#define MAXIMUM(x, y)   ( ((x) > (y)) ?  (x) : (y) )

// #define glColor3ubv(a) glColor3ub(a[0], a[1], a[2])
 
typedef struct shell1_widgets SHELL1_WIDGETS;
typedef struct glxarea_data GLXAREA_DATA;
typedef struct data_3D DATA_3D;
typedef struct data_3D_bounds DATA_3D_BOUNDS;
typedef struct graph_multi GRAPH_MULTI;
typedef struct appearance APPEARANCE;
typedef struct file_select_dialog FILE_SELECT_DIALOG;


// main.C:
extern XtAppContext app_context;
extern Display *display;       /*  Display             */
extern XVisualInfo             *vi;
extern int has_double_buffer;
extern GLXContext *share_context;
extern int display_lists_initialised;
extern float curquat[4];
extern double fovy;
extern int new_display_list;
extern APPEARANCE appearance;
extern int dsplnr, dsplnr_hi;
extern GLubyte black[3], brown[3], edge_col[3], ball_col[3], ball2_col[3];
extern GLfloat rot[3][3];
extern int redisplayPending;
extern XtWorkProcId redisplayID;
extern int rotate_beginx, rotate_beginy;
extern int spinning;
extern XtWorkProcId animateID;
extern int sphereVersion;
extern int pendingAutoHiRes;
extern XtWorkProcId hiResID;
extern float lastquat[4];
extern char *graph_fname;
extern FILE_SELECT_DIALOG *fsel_load_save;
extern int nb_select, selected_vertices[1000];

extern SHELL1_WIDGETS *s1;
extern GLXAREA_DATA *glx_data;

extern int graph_loaded;
extern GRAPH_MULTI *GM;
extern DATA_3D *D3;
extern DATA_3D_BOUNDS D3_bounds;

void create_shell1 (Display *display, char *app_name, int app_argc, char **app_argv);
void cb_empty(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_quit(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_curquat(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_toggle_print_show_labels(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_toggle_print_show_balls(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_toggle_show_balls(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_save(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_save_screen(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_reset(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_zoom_in(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_zoom_out(Widget w, XtPointer client_data, XtPointer xt_call_data);
void cb_zoom_reset(Widget w, XtPointer client_data, XtPointer xt_call_data);
void load_save_graph_func(FILE_SELECT_DIALOG *fsel);
int display_graph_latex(char *fname_graph, int f_print_show_balls);


//callbacks.C:
void glxarea_init(Widget w, XtPointer data, XtPointer callData);
void glxarea_resize(Widget w, XtPointer data, XtPointer callData);
void glxarea_draw(Widget w, XtPointer data, XtPointer callData);
void postRedisplay();
Boolean handleRedisplay(XtPointer callData);


// file_select.C:
FILE_SELECT_DIALOG *fsel_open(
	Display *display, char *app_name, int app_argc, char **app_argv, 
	char *filter_text, 
	char *list_text, 
	char *select_text, 
	char *pattern, 
	int f_name_used,
	char *f_name, 
	void (* do_func)(FILE_SELECT_DIALOG *fsel));
void fsel_realize(FILE_SELECT_DIALOG *fsel, int mode);
void fsel_unrealize(FILE_SELECT_DIALOG *fsel);
void cb_file_select(Widget w, XtPointer client_data, XtPointer xt_call_data);


// render.C:

#define HI_RES_SPHERE  1
#define LO_RES_SPHERE  2
#define DISK           3
#define RAND           4
#define ZYLINDER       5

#define STANDARD_SCALE_X 1500
#define STANDARD_SCALE_Y 1500
#define STANDARD_SCALE_Z 1500

extern GLfloat winWidth, winHeight;

void save_graph(char *fname, int f_show_labels);
void save_graph_screen(char *fname, int f_show_labels);
void renderScene();
void swap();
void init_viewport();
void init_projection();
void init_modelview();
void render_graph(GRAPH_MULTI *gm, DATA_3D *D3);
int is_selected(int v);
void compute_rotation(GLfloat x1, GLfloat y1, GLfloat z1, 
	GLfloat x2, GLfloat y2, GLfloat z2, 
	GLfloat *phi, GLfloat *nx, GLfloat *ny, GLfloat *nz, GLfloat *length, int f_debug);
void render_init(int share);
void render_reshape(GLXAREA_DATA *glx_data, int width, int height);
GLdouble compute_fovy();
GLdouble calc_edge_radius();
GLdouble calc_ball_radius();
void MyPushMatrix(int line);
void MyPopMatrix(int line);

// rotate.C:
void startRotation(Widget w, XEvent * event, 
	String * params, Cardinal * num_params);
void rotation(Widget w, XEvent * event, 
	String * params, Cardinal * num_params);
void animate(XtPointer closure, XtIntervalId *id);
void stopSpinning();
void makeHiRes();
void makeLoRes();
void stopAutoHiRes();
void postRedisplaySpin();
void hiresTimeout(XtPointer closure, XtIntervalId *id);
Boolean handleRedisplaySpin(XtPointer closure);

// graph.C:
DATA_3D *open_data_3D(int nb_V);
void free_data_3D(DATA_3D *p);
GRAPH_MULTI *open_graph_multi(int nb_V, int nb_E, int f_multi);
void free_graph_multi(GRAPH_MULTI *p);
int save_graph_multi_to_file(char *fname, GRAPH_MULTI *gm, DATA_3D *d3);
int read_graph_multi_from_file(char *fname, GRAPH_MULTI **gm, DATA_3D **d3);
void graph_standard_coordinates(GRAPH_MULTI *gm, DATA_3D *d3);
void graph_center(GRAPH_MULTI *gm, DATA_3D *d3, float *x0, float *y0, float *z0);
void data_3d_compute_bounds(DATA_3D *p, DATA_3D_BOUNDS *b);
void bounds_3D_print(DATA_3D_BOUNDS *b);
void bounds_3D_center(DATA_3D_BOUNDS *b, int *x0, int *y0, int *z0);
void data_3d_translate(DATA_3D *p, float dx, float dy, float dz);
void data_3d_scale(DATA_3D *p, float sx, float sy, float sz);


// quat.C
void vzero(float *v);
void vset(float *v, float x, float y, float z);
void vsub(const float *src1, const float *src2, float *dst);
void vcopy(const float *v1, float *v2);
void vcross(const float *v1, const float *v2, float *cross);
float vlength(const float *v);
void vscale(float *v, float div);
void vnormal(float *v);
float vdot(const float *v1, const float *v2);
void vadd(const float *src1, const float *src2, float *dst);
void trackball(float q[4], float p1x, float p1y, float p2x, float p2y);
/* Given an axis and angle, compute quaternion. */
void axis_to_quat(float a[3], float phi, float q[4]);
/* Project an x,y pair onto a sphere of radius r OR a
   hyperbolic sheet if we are away from the center of the
   sphere. */
float tb_project_to_sphere(float r, float x, float y);
/* Given two rotations, e1 and e2, expressed as quaternion
   rotations, figure out the equivalent single rotation and
   stuff it into dest.  This routine also normalizes the result
   every RENORMCOUNT times it is called, to keep error from
   creeping in.  NOTE: This routine is written so that q1 or q2
   may be the same as dest (or each other). */
void add_quats(float q1[4], float q2[4], float dest[4]);
/* Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0 If
   they don't add up to 1.0, dividing by their magnitude will
   renormalize them. */
void normalize_quat(float q[4]);
/* Build a rotation matrix, given a quaternion rotation. */
void build_rotmatrix(GLfloat m[4][4], float q[4]);


struct glxarea_data {
	Widget glxarea;
        Window glxwin;
	Dimension viewWidth, viewHeight;
		// Dimension of glxarea
	int GLXContext_initialised;
        GLXContext cx;
};

struct shell1_widgets {
  Widget shell1;
    Widget form57;
      // Widget frame4;
        // Widget Back_form;
          Widget Menu_bar;
            Widget File_menu;
              Widget menu3;
                Widget Save_mb;
                Widget Save_screen_mb;
                Widget separator_file1;
                Widget Toggle_print_show_labels;
                Widget Toggle_print_show_balls;
                Widget separator_file2;
                Widget Quit_mb;
	    Widget Control_menu;
	      Widget menu_control;
	        Widget curquat;
	        Widget reset;
	        Widget zoom_in;
	        Widget zoom_out;
	        Widget zoom_reset;
	        Widget separator_control;
		Widget Toggle_show_balls;
          Widget Discreta_title;
          Widget Main_form;
            Widget glxarea;
};

struct data_3D {
	int nb_V;
        float *x, *y, *z;
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

struct appearance {
	int show_balls;
	float viewpoint0[3];
	float viewpoint[3];
	float quat[4];
};



struct file_select_dialog {

	int f_realized;
	int mode;
	void (*do_func)(FILE_SELECT_DIALOG *fsel);

	/* labels in the dialog: */
	char filter_text[10000];
	char list_text[10000];
	char select_text[10000];

	char directory[10000];
	char pattern[10000];
	
	int f_name_used;
	char f_name[10000];
	
	Widget shell;
	Widget fsel_box;
};


