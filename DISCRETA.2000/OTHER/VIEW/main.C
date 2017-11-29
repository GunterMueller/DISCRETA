// main.C



#include "view.h"

// int glwidth = 300;
int glheight = 500;

XtAppContext app_context;
Display *display;       /*  Display             */
XVisualInfo             *vi;
int has_double_buffer = TRUE;
GLXContext *share_context = NULL;
int display_lists_initialised = FALSE;
float curquat[4];
double fovy = 30.;
int new_display_list = TRUE;
APPEARANCE appearance;
int dsplnr = 10, dsplnr_hi = 11;
GLubyte black[3], brown[3], edge_col[3], ball_col[3], ball2_col[3];
GLfloat rot[3][3];
int redisplayPending = FALSE;
XtWorkProcId redisplayID;
int rotate_beginx, rotate_beginy;
int spinning = FALSE;
XtWorkProcId animateID;
int sphereVersion = HI_RES_SPHERE;
int pendingAutoHiRes = FALSE;
XtWorkProcId hiResID;
float lastquat[4];
char *graph_fname = NULL;
char graph_fname_base[1024];
FILE_SELECT_DIALOG *fsel_load_save = NULL;
int nb_select = 0, selected_vertices[1000];

// for save / latex display:
int f_print_show_labels = FALSE;
int f_print_show_balls = FALSE;


SHELL1_WIDGETS s1_data, *s1 = &s1_data;
GLXAREA_DATA glx_data_, *glx_data = &glx_data_;

int graph_loaded = FALSE;
GRAPH_MULTI gm_, *GM = &gm_;
DATA_3D d3_, *D3 = &d3_;
DATA_3D_BOUNDS D3_bounds;

void create_shell1 (Display *display, char *app_name, int app_argc, char **app_argv);
static void initialise_objects( Widget parent );
void cb_empty(Widget w, XtPointer client_data, XtPointer xt_call_data);

static int _xd_initialised = 0;

typedef struct FontResources_s { 
	XmFontList fontlist1;
	XmFontList bfont;
	XmFontList courier;
} FontResources_t, *FontResources_p;

static FontResources_t font_resources;

typedef struct PixelResources_s { 
	Pixel color1;
	Pixel color2;
	Pixel color3;
} PixelResources_t, *PixelResources_p;

static PixelResources_t pixel_resources;






XtActionsRec actionsTable[] =
{
  {"startRotation", startRotation},
  {"rotation", rotation},
};

char *glxareaTranslations =
  "#override\n\
  <Btn1Down>:startRotation()\n\
  <Btn1Motion>:rotation()\n";
XtTranslations xt_translations_glxarea;

static int config[] =
{
  None, None,
  GLX_DOUBLEBUFFER, GLX_RGBA, GLX_DEPTH_SIZE, 12,
  GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1,
  None
};


static int no_depth_config[] =
{ None, None, GLX_RGBA, GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1,
GLX_BLUE_SIZE, 1, None };


static int *dblBuf = &config[2];
static int *snglBuf = &config[3];
static int *no_depthBuf = no_depth_config;











int main (int argc, char **argv)
{
	int i, j;
	char *p;
	
	XtSetLanguageProc ( (XtAppContext) NULL, (XtLanguageProc) NULL, (XtPointer) NULL );
	XtToolkitInitialize ();
	app_context = XtCreateApplicationContext ();
	display = XtOpenDisplay (app_context, NULL, argv[0], "XApplication",
	                         NULL, 0, &argc, argv);
	if (!display)
	{
	    printf("%s: can't open display, exiting...\n", argv[0]);
	    exit (-1);
	}
	appearance.show_balls = FALSE;
	appearance.viewpoint0[0] = 0.;
	appearance.viewpoint0[1] = 0;
	appearance.viewpoint0[2] = -4.5;
	for (i = 0; i < 3; i++) {
		appearance.viewpoint[i] = appearance.viewpoint0[i];
		}
	appearance.quat[0] = -0.115618;
	appearance.quat[1] = 0.001670;
	appearance.quat[2] = 0.003509;
	appearance.quat[3] = 0.993286;
	{
	FILE *fp;
	fp = fopen("default_quat", "r");
	if (fp != NULL) {
		fscanf(fp, "%f %f %f %f", 
			&appearance.quat[0], 
			&appearance.quat[1], 
			&appearance.quat[2], 
			&appearance.quat[3]);
		fclose(fp);
		}
	}

	black[0] = 0;
	black[1] = 0;
	black[2] = 0;
	brown[0] = 0;
	brown[1] = 100;
	brown[2] = 100;
	edge_col[0] = 255;
	edge_col[1] = 0;
	edge_col[2] = 0;
	ball_col[0] = 255;
	ball_col[1] = 255;
	ball_col[2] = 0;
	ball2_col[0] = 255;
	ball2_col[1] = 0;
	ball2_col[2] = 255;
	if (argc < 2) {
		printf("usage: view [options] graph_file\n");
		printf("where options can be:\n");
		printf("--balls                  : show balls\n");
		printf("--quat f1 f2 f3 f4       : initial rotation,\n");
		printf("                         : f1 , ... , f4 are floating point numbers\n");
		printf("--select n v1 v2 ... vn  : show vertices v1,..,vn in a different color\n");
		XtAppError(app_context, "wrong number of arguments");
		}
#ifdef VERBOSE
	printf("argc=%d\n", argc);
	for (i = 1; i < argc - 1; i++) {
		printf("arg %d: %s\n", i, argv[i]);
		}
#endif
	for (i = 1; i < argc - 1; i++) {
#ifdef VERBOSE
		printf("parsing arg %d: %s\n", i, argv[i]);
#endif
		if (strcmp(argv[i], "--balls") == 0) {
			appearance.show_balls = TRUE;
			}
		else if (strcmp(argv[i], "--quat") == 0) {
			if (argc < i + 4)
				XtAppError(app_context, "wrong number of arguments");
			sscanf(argv[++i], "%f", &appearance.quat[0]);
			sscanf(argv[++i], "%f", &appearance.quat[1]);
			sscanf(argv[++i], "%f", &appearance.quat[2]);
			sscanf(argv[++i], "%f", &appearance.quat[3]);
			}
		else if (strcmp(argv[i], "--select") == 0) {
			sscanf(argv[++i], "%d", &nb_select);
			if (argc < i + nb_select)
				XtAppError(app_context, "wrong number of arguments");
			for (j = 0; j < nb_select; j++) {
				sscanf(argv[++i], "%d", &selected_vertices[j]);
				}
#ifdef VERBOSE
			printf("-select %d ", nb_select);
			for (j = 0; j < nb_select; j++) {
				printf("%d ", selected_vertices[j]);
				}
			printf("\n");
#endif
			}
		}
	graph_fname = argv[argc - 1];
	strcpy(graph_fname_base, graph_fname);
	if ((p = strrchr(graph_fname_base, '.')))
		*p = 0;
	// glutInit(&argc, argv);

	vi = glXChooseVisual(display, DefaultScreen(display), dblBuf);
	has_double_buffer = TRUE;
	if (vi == NULL) {
		has_double_buffer = FALSE;
		vi = glXChooseVisual(display, DefaultScreen(display), snglBuf);
		if (vi == NULL) {
			vi = glXChooseVisual(display, DefaultScreen(display), no_depthBuf);
			if (vi)
				printf("no depth buffer\n");
			else
 				XtAppError(app_context, "no RGB visual");
			}
		}

	printf("opening file %s for reading\n", graph_fname); fflush(stdout);
	
	read_graph_multi_from_file(graph_fname, &GM, &D3);
	graph_standard_coordinates(GM, D3);
	graph_standard_coordinates(GM, D3);
		// we do it twice because the first time errors 
		// may come in from rounding
	save_graph_multi_to_file("normalized.graph", GM, D3);
#ifdef VERBOSE
	bounds_3D_print(&D3_bounds);
#endif
	graph_loaded = TRUE;
	
	
	
	xt_translations_glxarea = XtParseTranslationTable(glxareaTranslations);
#ifdef VERBOSE
	printf("translations parsed\n"); fflush(stdout);
#endif
	
	XtAppAddActions(app_context, actionsTable, 2);
#ifdef VERBOSE
	printf("xt actions added\n"); fflush(stdout);
#endif
	
	// discreta_init();
	create_shell1 ( display, argv[0], argc, argv );

	printf("calling fsel_open\n"); fflush(stdout);
	fsel_load_save = fsel_open(display, argv[0], argc, argv, 
		"file filter", 
		"graph file", 
		"selection", 
		"*.graph", 
		TRUE /* f_name_used */,
		graph_fname, load_save_graph_func);
	printf("finished with fsel_open\n"); fflush(stdout);



	// program_init();
	XtRealizeWidget (s1->shell1);
	XtAppMainLoop (app_context);
	exit (0);
}

void create_shell1 (Display *display, char *app_name, int app_argc, char **app_argv)
{
	Widget children[12];      /* Children to manage */
	Arg al[64];                    /* Arg List */
	register int ac = 0;           /* Arg Count */
	// XrmValue from_value, to_value; /* For resource conversion */
	// Pixel fg, bg;                    /* colour values for pixmaps */ 
	XmString xmstrings[16];    /* temporary storage for XmStrings */

	XtSetArg(al[ac], XmNallowShellResize, TRUE); ac++;
	XtSetArg(al[ac], XmNargc, app_argc); ac++;
	XtSetArg(al[ac], XmNargv, app_argv); ac++;
	s1->shell1 = XtAppCreateShell ( app_name, "XApplication", applicationShellWidgetClass, display, al, ac );
	ac = 0;
	initialise_objects ( s1->shell1 );
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	XtSetArg(al[ac], XmNautoUnmanage, FALSE); ac++;
	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	s1->form57 = XmCreateForm ( s1->shell1, "form57", al, ac );
#if 0
	ac = 0;
	s1->frame4 = XmCreateFrame ( s1->form57, "frame4", al, ac );
	XtSetArg(al[ac], XmNborderColor, pixel_resources.color1); ac++;
	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	s1->Back_form = XmCreateForm ( s1->frame4, "back_form", al, ac );
#endif
	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	s1->Menu_bar = XmCreateMenuBar ( s1->form57 /* s1->Back_form */, "menu_bar", al, ac );
	
	
	
	// File menu:
	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "File", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->File_menu = XmCreateCascadeButton ( s1->Menu_bar, "file_menu", al, ac );
	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	s1->menu3 = XmCreatePulldownMenu ( s1->Menu_bar, "menu3", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "Save", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->Save_mb = XmCreatePushButton ( s1->menu3, "save_mb", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "Save from screen", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->Save_screen_mb = XmCreatePushButton ( s1->menu3, "save_from_screen_mb", al, ac );

	ac = 0;
	s1->separator_file1 = XmCreateSeparator ( s1->menu3, "separator_file1", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNstringDirection, XmSTRING_DIRECTION_L_TO_R); ac++;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "show labels", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNmarginTop, 1); ac++;
	XtSetArg(al[ac], XmNmarginBottom, 1); ac++;
	XtSetArg(al[ac], XmNmarginLeft, 19); ac++;
	XtSetArg(al[ac], XmNalignment, XmALIGNMENT_BEGINNING); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	XtSetArg(al[ac], XmNfillOnSelect, TRUE); ac++;
	XtSetArg(al[ac], XmNindicatorType, XmN_OF_MANY); ac++;
	if (f_print_show_labels) {
		XtSetArg(al[ac], XmNset, TRUE); ac++;
		}
	else {
		XtSetArg(al[ac], XmNset, FALSE); ac++;
		}
	XtSetArg(al[ac], XmNvisibleWhenOff, TRUE); ac++;
	XtSetArg(al[ac], XmNindicatorSize, 14); ac++;
	s1->Toggle_print_show_labels = XmCreateToggleButton ( s1->menu3, "toggle_print_show_labels", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNstringDirection, XmSTRING_DIRECTION_L_TO_R); ac++;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "show balls", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNmarginTop, 1); ac++;
	XtSetArg(al[ac], XmNmarginBottom, 1); ac++;
	XtSetArg(al[ac], XmNmarginLeft, 19); ac++;
	XtSetArg(al[ac], XmNalignment, XmALIGNMENT_BEGINNING); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	XtSetArg(al[ac], XmNfillOnSelect, TRUE); ac++;
	XtSetArg(al[ac], XmNindicatorType, XmN_OF_MANY); ac++;
	if (f_print_show_balls) {
		XtSetArg(al[ac], XmNset, TRUE); ac++;
		}
	else {
		XtSetArg(al[ac], XmNset, FALSE); ac++;
		}
	XtSetArg(al[ac], XmNvisibleWhenOff, TRUE); ac++;
	XtSetArg(al[ac], XmNindicatorSize, 14); ac++;
	s1->Toggle_print_show_balls = XmCreateToggleButton ( s1->menu3, "toggle_print_show_balls", al, ac );

	ac = 0;
	s1->separator_file2 = XmCreateSeparator ( s1->menu3, "separator_file2", al, ac );




	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "Quit", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->Quit_mb = XmCreatePushButton ( s1->menu3, "quit_mb", al, ac );

	// Control menu
	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "Control", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->Control_menu = XmCreateCascadeButton ( s1->Menu_bar, "control_menu", al, ac );
	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	s1->menu_control = XmCreatePulldownMenu ( s1->Menu_bar, "menu_control", al, ac );
	
	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "show current quat", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->curquat = XmCreatePushButton ( s1->menu_control, "show_curquat", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "reset rotation", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->reset = XmCreatePushButton ( s1->menu_control, "reset_rotation", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "zoom in", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->zoom_in = XmCreatePushButton ( s1->menu_control, "zoom_in", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "zoom out", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->zoom_out = XmCreatePushButton ( s1->menu_control, "zoom_out", al, ac );

	ac = 0;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "reset zoom", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	s1->zoom_reset = XmCreatePushButton ( s1->menu_control, "reset_zoom", al, ac );

	ac = 0;
	s1->separator_control = XmCreateSeparator ( s1->menu_control, "separator_control", al, ac );
	XtSetArg(al[ac], XmNstringDirection, XmSTRING_DIRECTION_L_TO_R); ac++;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color3); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "show balls", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNmarginTop, 1); ac++;
	XtSetArg(al[ac], XmNmarginBottom, 1); ac++;
	XtSetArg(al[ac], XmNmarginLeft, 19); ac++;
	XtSetArg(al[ac], XmNalignment, XmALIGNMENT_BEGINNING); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	XtSetArg(al[ac], XmNfillOnSelect, TRUE); ac++;
	XtSetArg(al[ac], XmNindicatorType, XmN_OF_MANY); ac++;
	if (appearance.show_balls) {
		XtSetArg(al[ac], XmNset, TRUE); ac++;
		}
	else {
		XtSetArg(al[ac], XmNset, FALSE); ac++;
		}
	XtSetArg(al[ac], XmNvisibleWhenOff, TRUE); ac++;
	XtSetArg(al[ac], XmNindicatorSize, 14); ac++;
	s1->Toggle_show_balls = XmCreateToggleButton ( s1->menu_control, "toggle_show_balls", al, ac );




	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	XtSetArg(al[ac], XmNstringDirection, XmSTRING_DIRECTION_L_TO_R); ac++;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color1); ac++;
	xmstrings[0] = XmStringCreateLtoR ( "*** DISCRETA ***", (XmStringCharSet)XmFONTLIST_DEFAULT_TAG );
	XtSetArg(al[ac], XmNlabelString, xmstrings[0]); ac++;
	XtSetArg(al[ac], XmNfontList, font_resources.bfont); ac++;
	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, s1->Menu_bar); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	s1->Discreta_title = XmCreateLabel ( s1->form57 /* s1->Back_form */, "title", al, ac );
	ac = 0;
	XmStringFree ( xmstrings [ 0 ] );
	XtSetArg(al[ac], XmNwidth, 700); ac++;
	XtSetArg(al[ac], XmNstringDirection, XmSTRING_DIRECTION_L_TO_R); ac++;
	XtSetArg(al[ac], XmNbackground, pixel_resources.color1); ac++;
	XtSetArg(al[ac], XmNmarginWidth, 0); ac++;
	XtSetArg(al[ac], XmNmarginHeight, 0); ac++;
	XtSetArg(al[ac], XmNautoUnmanage, FALSE); ac++;
	XtSetArg(al[ac], XmNbuttonFontList, font_resources.fontlist1); ac++;
	XtSetArg(al[ac], XmNlabelFontList, font_resources.fontlist1); ac++;
	XtSetArg(al[ac], XmNtextFontList, font_resources.fontlist1); ac++;
	XtSetArg(al[ac], XmNtopAttachment, XmATTACH_WIDGET); ac++;
	XtSetArg(al[ac], XmNtopWidget, s1->Discreta_title); ac++;
	XtSetArg(al[ac], XmNleftAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNrightAttachment, XmATTACH_FORM); ac++;
	XtSetArg(al[ac], XmNbottomAttachment, XmATTACH_FORM); ac++;
	s1->Main_form = XmCreateForm ( s1->form57 /* s1->Back_form */, "main_form", al, ac );



	s1->glxarea = 
		XtVaCreateManagedWidget("glxarea", glwMDrawingAreaWidgetClass,
					s1->Main_form, 
					GLwNvisualInfo, vi, 
					// XmNwidth, glwidth,
					XmNheight, glheight-15,
					XmNtopAttachment, XmATTACH_FORM,
					// XmNtopOffset, 1,
					XmNleftAttachment, XmATTACH_FORM,
					XmNrightAttachment, XmATTACH_FORM,
					XmNbottomAttachment, XmATTACH_FORM,
					NULL);



	ac = 0;
	XtAddCallback( s1->Save_mb, XmNactivateCallback, cb_save, (XtPointer) 0 );
	children[ac++] = s1->Save_mb;
	XtManageChildren(children, ac);
	ac = 0;
	XtAddCallback( s1->Save_screen_mb, XmNactivateCallback, cb_save_screen, (XtPointer) 0 );
	children[ac++] = s1->Save_screen_mb;
	XtManageChildren(children, ac);
	ac = 0;
	children[ac++] = s1->separator_file1;
	children[ac++] = s1->Toggle_print_show_labels;
	children[ac++] = s1->Toggle_print_show_balls;
	children[ac++] = s1->separator_file2;
	XtManageChildren(children, ac);

	ac = 0;
	XtAddCallback( s1->Quit_mb, XmNactivateCallback, cb_quit, (XtPointer) 0 );
	children[ac++] = s1->Quit_mb;
	XtManageChildren(children, ac);
	ac = 0;
	XtSetArg(al[ac], XmNsubMenuId, s1->menu3); ac++;
	XtSetValues ( s1->File_menu, al, ac );

	XtAddCallback( s1->Toggle_print_show_labels, XmNvalueChangedCallback, cb_toggle_print_show_labels, (XtPointer) 0 );
	XtAddCallback( s1->Toggle_print_show_balls, XmNvalueChangedCallback, cb_toggle_print_show_balls, (XtPointer) 0 );

	XtAddCallback( s1->curquat, XmNactivateCallback, cb_curquat, (XtPointer) 0 );
	XtAddCallback( s1->reset, XmNactivateCallback, cb_reset, (XtPointer) 0 );
	XtAddCallback( s1->zoom_in, XmNactivateCallback, cb_zoom_in, (XtPointer) 0 );
	XtAddCallback( s1->zoom_out, XmNactivateCallback, cb_zoom_out, (XtPointer) 0 );
	XtAddCallback( s1->zoom_reset, XmNactivateCallback, cb_zoom_reset, (XtPointer) 0 );
	XtAddCallback( s1->Toggle_show_balls, XmNvalueChangedCallback, cb_toggle_show_balls, (XtPointer) 0 );

	ac = 0;
	children[ac++] = s1->curquat;
	children[ac++] = s1->reset;
	children[ac++] = s1->zoom_in;
	children[ac++] = s1->zoom_out;
	children[ac++] = s1->zoom_reset;
	children[ac++] = s1->separator_control;
	children[ac++] = s1->Toggle_show_balls;
	XtManageChildren(children, ac);
	ac = 0;
	XtSetArg(al[ac], XmNsubMenuId, s1->menu_control); ac++;
	XtSetValues ( s1->Control_menu, al, ac );


	ac = 0;
	children[ac++] = s1->File_menu;
	children[ac++] = s1->Control_menu;
	XtManageChildren(children, ac);

	ac = 0;
	children[ac++] = s1->glxarea;
	XtManageChildren(children, ac);

	ac = 0;
	// XtManageChild(xdref_11);
	children[ac++] = s1->Menu_bar;
	children[ac++] = s1->Discreta_title;
	children[ac++] = s1->Main_form;
	XtManageChildren(children, ac);
	ac = 0;
	XtManageChild ( s1->form57);


	glx_data->glxarea = s1->glxarea;
	
	XtVaGetValues(s1->glxarea, 
		XtNwidth, &glx_data->viewWidth, 
		XtNheight, &glx_data->viewHeight, 
		NULL);
#ifdef VERBOSE
	printf("viewWidth=%ld viewHeight=%ld\n", 
		(long int) glx_data->viewWidth, 
		(long int) glx_data->viewHeight);
#endif
	glx_data->GLXContext_initialised = FALSE;
	
	XtAddCallback(s1->glxarea, XmNexposeCallback, glxarea_draw, (void*)NULL);
	XtAddCallback(s1->glxarea, XmNresizeCallback, glxarea_resize, (void*)NULL);
	XtAddCallback(s1->glxarea, GLwNginitCallback, glxarea_init, (void*)NULL);
	
	XtOverrideTranslations( s1->glxarea, xt_translations_glxarea );
}

static void initialise_objects( Widget parent )
{
	XrmValue from_value, to_value; /* For resource conversions */
	if ( _xd_initialised ) return;
	_xd_initialised = 1;
	while ( XtParent ( parent ) )
		parent = XtParent ( parent );
	from_value.addr = "-adobe-helvetica-bold-r-normal--20-140-100-100-p-105-iso8859-1";
	from_value.size = strlen(from_value.addr)+1;
	to_value.addr = NULL;
	XtConvertAndStore( parent, XmRString, &from_value, XmRFontList, &to_value);
	if ( to_value.addr )
		font_resources.fontlist1 = *(XmFontList*)to_value.addr;
	from_value.addr = "-adobe-helvetica-bold-r-normal--20-140-100-100-p-105-iso8859-1";
	from_value.size = strlen(from_value.addr)+1;
	to_value.addr = NULL;
	XtConvertAndStore( parent, XmRString, &from_value, XmRFontList, &to_value);
	if ( to_value.addr )
		font_resources.bfont = *(XmFontList*)to_value.addr;
	from_value.addr = "-*-courier-medium-r-*-*-14-*-*-*-*-*-*-*";
	from_value.size = strlen(from_value.addr)+1;
	to_value.addr = NULL;
	XtConvertAndStore( parent, XmRString, &from_value, XmRFontList, &to_value);
	if ( to_value.addr )
		font_resources.courier = *(XmFontList*)to_value.addr;
	if (DefaultDepthOfScreen(DefaultScreenOfDisplay(XtDisplay(parent))) != 1) {
	from_value.addr = "wheat1";
	from_value.size = strlen(from_value.addr)+1;
	to_value.addr = NULL;
	XtConvertAndStore( parent, XmRString, &from_value, XmRPixel, &to_value);
	if ( to_value.addr )
		pixel_resources.color1 = *(Pixel*)to_value.addr;
	from_value.addr = "royalblue4";
	from_value.size = strlen(from_value.addr)+1;
	to_value.addr = NULL;
	XtConvertAndStore( parent, XmRString, &from_value, XmRPixel, &to_value);
	if ( to_value.addr )
		pixel_resources.color2 = *(Pixel*)to_value.addr;
	from_value.addr = "wheat2";
	from_value.size = strlen(from_value.addr)+1;
	to_value.addr = NULL;
	XtConvertAndStore( parent, XmRString, &from_value, XmRPixel, &to_value);
	if ( to_value.addr )
		pixel_resources.color3 = *(Pixel*)to_value.addr;
	}
}

void cb_empty(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
}

void cb_quit(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	exit(0);
}

void cb_curquat(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	
	normalize_quat(curquat);
	printf("current quat: %f %f %f %f\n", 
		curquat[0], curquat[1], curquat[2], curquat[3]);
	
}

void cb_toggle_print_show_labels(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	XmToggleButtonCallbackStruct *call_data = (XmToggleButtonCallbackStruct *) xt_call_data ;
	
	printf("cb_toggle_print_show_labels() value changed set=%d\n", call_data->set); fflush(stdout);
	if (call_data->set)
		f_print_show_labels = TRUE;
	else
		f_print_show_labels = FALSE;
}

void cb_toggle_print_show_balls(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	XmToggleButtonCallbackStruct *call_data = (XmToggleButtonCallbackStruct *) xt_call_data ;
	
	printf("cb_toggle_print_show_balls() value changed set=%d\n", call_data->set); fflush(stdout);
	if (call_data->set)
		f_print_show_balls = TRUE;
	else
		f_print_show_balls = FALSE;
}

void cb_toggle_show_balls(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	XmToggleButtonCallbackStruct *call_data = (XmToggleButtonCallbackStruct *) xt_call_data ;
	
	printf("cb_toggle_show_balls() value changed set=%d\n", call_data->set); fflush(stdout);
	if (call_data->set)
		appearance.show_balls = TRUE;
	else
		appearance.show_balls = FALSE;
	new_display_list = TRUE;
	postRedisplay();
}

void cb_save(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	printf("cb_save()\n"); fflush(stdout);
	// save_graph();
	fsel_realize(fsel_load_save, 1);
}

void cb_save_screen(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	printf("cb_save_screen()\n"); fflush(stdout);
	// save_graph_screen();
	fsel_realize(fsel_load_save, 2);
}

void cb_reset(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	printf("cb_reset()\n"); fflush(stdout);
	// trackball(curquat, 0.0, 0.0, 0.0, 0.0);
	curquat[0] = appearance.quat[0];
	curquat[1] = appearance.quat[1];
	curquat[2] = appearance.quat[2];
	curquat[3] = appearance.quat[3];
	postRedisplay();
}

void cb_zoom_in(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	appearance.viewpoint[2] += 0.2;
	// new_display_list = TRUE;
	postRedisplay();
}

void cb_zoom_out(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	appearance.viewpoint[2] -= 0.2;
	// new_display_list = TRUE;
	postRedisplay();
}

void cb_zoom_reset(Widget w, XtPointer client_data, XtPointer xt_call_data)
{
	int i;
	// XmPushButtonCallbackStruct *call_data = (XmPushButtonCallbackStruct *) xt_call_data ;
	for (i = 0; i < 3; i++) {
		appearance.viewpoint[i] = appearance.viewpoint0[i];
		}
	// new_display_list = TRUE;
	postRedisplay();
}

void load_save_graph_func(FILE_SELECT_DIALOG *fsel)
{

	if (fsel->mode == 0) { // load
		printf("load\n"); fflush(stdout);
		}
	else if (fsel->mode == 1) { // save
		printf("save\n"); fflush(stdout);
		save_graph(fsel->f_name, f_print_show_labels);
		// save_graph_screen(fsel->f_name);
		display_graph_latex(fsel->f_name, f_print_show_balls);
		}
	else if (fsel->mode == 2) { // save_screen
		printf("save_screen\n"); fflush(stdout);
		save_graph_screen(fsel->f_name, f_print_show_labels);
		display_graph_latex(fsel->f_name, f_print_show_balls);
		}
	fsel_unrealize(fsel);
}

int display_graph_latex(char *fname_graph, int f_print_show_balls)
{
	char fname_base[1024];
	char str[1024];
	char *fname_latex_base = "view";
	char cmd[1024];
	char *p;
	FILE *fp;
	int i;
	
	strcpy(fname_base, fname_graph);
	if ((p = strrchr(fname_base, '.')) != NULL)
		*p = 0;
	sprintf(cmd, "t117.out -scale 0.3 0.3");
	if (nb_select) {
		sprintf(cmd + strlen(cmd), " -select %d", nb_select);
		for (i = 0; i < nb_select; i++) {
			sprintf(cmd + strlen(cmd), " %d", selected_vertices[i]);
			}
		}
	if (f_print_show_balls)
		sprintf(cmd + strlen(cmd), " -dots");
	sprintf(cmd + strlen(cmd), " %s %s", fname_base, fname_base);
	printf("executing: %s\n", cmd);
	system(cmd);
	sprintf(cmd, "mpost %s.mp", fname_base);
	system(cmd);
	sprintf(str, "%s.tex", fname_latex_base);
	fp = fopen(str, "w");
	fprintf(fp, "\\documentclass[]{article}\n");
	fprintf(fp, "\\usepackage{epsfig}\n");
	fprintf(fp, "\\evensidemargin 0in\n");
	fprintf(fp, "\\oddsidemargin 0in\n");
	fprintf(fp, "\\marginparwidth 0pt\n");
	fprintf(fp, "\\marginparsep 0pt\n");
	fprintf(fp, "\\topmargin -1in\n");
	fprintf(fp, "\\headheight 0.7cm\n");
	fprintf(fp, "\\headsep 1.8cm\n");
	fprintf(fp, "%%\\footheight 0.7cm\n");
	fprintf(fp, "\\footskip 2cm\n");
	fprintf(fp, "\\textheight 22cm\n");
	fprintf(fp, "\\textwidth 6.2in\n");
	fprintf(fp, "\\marginparpush 0pt\n");
	fprintf(fp, "%%\\title{}\n");
	fprintf(fp, "%%\\author{{\\sc }}\n");
	fprintf(fp, "%%\\date{}\n");
	fprintf(fp, "\n");
	fprintf(fp, "\\begin{document}\n");
	fprintf(fp, "\n");
	fprintf(fp, "\\begin{center}\n");
	fprintf(fp, "\\epsfig{file=%s.1,width=150mm}\n", fname_base);
	fprintf(fp, "\\end{center}\n");
	fprintf(fp, "\\end{document}\n");
	fprintf(fp, "\n");
	fclose(fp);
	sprintf(cmd, "latex %s.tex", fname_latex_base);
	system(cmd);
	sprintf(cmd, "dvips %s.dvi -o", fname_latex_base);
	system(cmd);
	sprintf(cmd, "ghostview %s.ps &", fname_latex_base);
	system(cmd);
	return TRUE;
}



