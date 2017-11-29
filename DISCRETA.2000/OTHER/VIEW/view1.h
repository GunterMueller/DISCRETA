#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <gl.h>
#include <glx.h>

class MolGLZeichenWindowData
{
public:
  Bool                   GLXContext_initialised;
  GLXContext             cx;
  
  int                    spinning;
  int                    pendingAutoHiRes;
  int                    beginx, beginy;
  int                    redisplayPending ;

  int                    picked_atom;
  int                    pick_x, pick_y;

  float                  lastquat[4];

  MOLBOOL                mol_loaded;

  MOLECULE_DATA          moldat;

  APPEAR_OPTIONS         appearence;

  int                    displaylist_nr, displaylist_nr_hi;
  MOLBOOL                neue_displayliste;

};




class MolGLZeichenWindow : public MolGLZeichenWindowData
{
public:
  Widget                 glxarea;
  Window                 glxwin;
  Dimension              viewWidth, viewHeight;

  XtWorkProcId           animateID, hiResID, redisplayID;

  Widget                 *formula_labels;
  int                    formula_labels_anz;
  int                    formula_labels_mapped;
};





typedef struct
{
  // Hauptfenster:
  Widget toplevel;
    Widget main_form;
      Widget menu;
        Widget File;
        Widget file_menu;
          Widget Open;
          Widget Reread;
          Widget Export;
          Widget Print;
          Widget Exit;
        Widget Options;
        Widget opt_menu;
          Widget Grid;
          Widget Representation;
          Widget Repr_pulldown;
            Widget toggle2D_3D;
            Widget toggle_recomp_coords;
            Widget spin_toggle;
            Widget anzeige; 
          Widget Select;
          Widget GeheZu;
          Widget Perspective;
        Widget  Help;
        Widget  help_menu;
          Widget  help_help;
          Widget  help_about;
      Widget scroll; 
      Widget gl_form;
        Widget *main_frame; 
           Widget *frame_forms;
              Widget *frame_headings;
              Widget *heading_number;
              Widget *heading_name;

  // Popup auf rechte Maustaste:
  Widget popup;
     Widget  popup_export_pulldown;
     Widget  popup_export;
        Widget  popup_exp_file;
        Widget  popup_moled;
     Widget  popup_select;
     Widget  popup_name;
     Widget  popup_representation;
     Widget  popup_representation_2D_pulldown;
        Widget  popup_2DTo3D;
        Widget  popup_arrange;
        Widget  popup_appear_2D;
     Widget  popup_representation_3D_pulldown;
        Widget  popup_3DTo2D;
        Widget  popup_momentum;
        Widget  popup_optimize_pulldown;
        Widget  popup_optimize;
           Widget  popup_mm2;
           Widget  popup_eval;
        Widget  popup_appear_3D;
     

  // Dialoge:
  MV_GRID_STR          grid_str;
  MV_GOTO_STR          goto_str;
  MV_EVAL_STR          eval_str;
  MV_MOLNAME_STR       name_str;
  MV_SELECT_STR        sel_str;
  MV_EXPORT_STR        exp_str;
  MV_EXPORTCOUNT_STR   exp_count_str;
  MV_PRINT_STR         print_str;
  MV_PRINTCOUNT_STR    print_count_str;
  MV_OPTI_STR          opti_str;
  MV_PERSPEKTIVE_STR   persp_str;

  // 
  APPEAR_OPTIONS          appearence;
  APPEAR_STRUCT           app_struct;
  LIZENZ_WINDOW           liz_str;
  MOLGEN_FONT_STRUCT      m_fstr;
  MOLGEN_COLOR_STRUCT     m_cstr;

  // globale Optionen
  // render:
  double                  fovy;
  float                   scale_factor;
  VEK_1D<float>           atom_radius;
  GLubyte                 atomcolor[20][3]; 
  GLubyte                 schwarz[3], brown[3], edge_col[3];
  float                   default_radius;
  GLubyte                 default_color[3];
  VEK_2D<double>          M, S;
  MOLBOOL                 male_atomnummern;


  // molfile:
  MOLBOOL                 use_stored_coords;



  // gui_run:
  STRING_TG               mol_in_file;
  int                     mols_in_file;
  char                    file_format;
  VEK_1D<int>             molfile_read_idx;
  MOLBOOL                 file_opened;
  GLXContext              *share_context;
  int                     last_spin_context;




  



  int                     el_anzg;
  ARRAY<STRING_TG>        el_names;
  VEK_1D<int>             el_anz;



  int                     popup_over_idx;
  MOLBOOL                 mit_filearg;
  STRING_TG               start_file;
  int                     contexts_init_anz;
  char                    is_drawn[300];
  MOLBOOL                 start_file_loaded;
  GLfloat                 rot[3][3], c[3][3], rotinv[3][3];
  ARRAY<STRING_TG>        new_molnames;
  VEK_1D<int>             mol_name_idx;

  // Aendern der 2D Koordinaten:
  int                     actual_alter_idx;
  VEK_1D<int>             altered_numbers;
  VEK_2D<double>          x2_altered, y2_altered;

  // Aendern der 3D Koordinaten:
  VEK_1D<int>             altered_numbers_3D;
  VEK_2D<double>          x3_altered, y3_altered, z3_altered;

  VEK_1D<int>             altered_numbers_quat;
  VEK_2D<float>           altered_quats;

  // Drucken in 3D:
  GLXContext              *cx3D;
  Pixmap                  *pmap;
  GLXPixmap               *glxpmap;
  MOLECULE_DATA           *moldat3D;

  DIR_ENTRY_TG            drag_export;
  MOLBOOL                 is_dragging;
  int                     pageup_to, pagedown_to;
  Widget                  PageWidget;
  STRING_TG               fontfile;
  MOLBOOL                 mit_double_buffer;
  MOLBOOL                 direct_rendering;

  MOLBOOL                 display_lists_initialised;

  Visual                  *overlayVisual;
  int                     overlayDepth;
  Colormap                overlayColormap;
  XtTranslations          trans;
  XtWorkProcId            load_id;
  XtAppContext            *app;
  Display                 *dpy;
  XVisualInfo             *vi;
  
  STRING_TG               exe_dir;

} MolGLGlobals;

MOLBOOL LeseMolecule(MOLECULE_DATA &moldat, APPEAR_OPTIONS &appearence, int molnr, MOLBOOL mrt);
extern void ErrorMessageTimeOutGL(XtPointer closure, XtIntervalId *id);
extern double BerechneFovy();
extern MolGLZeichenWindow *zeichen_windows;
extern MolGLGlobals       *MGLG;
extern int zeichen_window_anz, zeichen_window_used;
extern void ExportiereMols( Widget widget, XtPointer client_data, XtPointer call_data );
extern void ErstellePlazierung(MOLECULE_DATA &moldat, APPEAR_OPTIONS &appearence, int w, int h);
extern void ScrollCB( Widget widget, XtPointer client_data, XtPointer call_data );

extern void Export2DClicked( Widget widget, XtPointer client_data, XtPointer call_data );
extern void CopyGlobalAppearFlags(APPEAR_OPTIONS &appearence);
extern void ExportCB( Widget widget, XtPointer client_data, XtPointer call_data );
extern void PrintCB( Widget widget, XtPointer client_data, XtPointer call_data );
extern void PlaziereMol2D(MOLECULE_DATA &moldat);
extern void InitFarbe(MOLECULE_DATA &moldat );
extern MOLBOOL IsSelected(int molnr);
extern void PopupSelectCB(Widget w, XtPointer clientData, XtPointer callData);
extern void MenuSelectCB(Widget w, XtPointer clientData, XtPointer callData);
extern void SelectListClicked(Widget w, XtPointer clientData, XtPointer callData);
extern void select_selectbuttonCB(Widget w, XtPointer clientData, XtPointer callData);
extern void select_unselectbuttonCB(Widget w, XtPointer clientData, XtPointer callData);

extern XtActionsRec actionsTable[];

extern void startRotation(Widget w, XEvent * event, String * params, Cardinal * num_params);
extern void rotation(Widget w, XEvent * event, String * params, Cardinal * num_params);
extern void doPick(Widget w, XEvent * event, String * params, Cardinal * num_params);

extern void mapStateChanged(Widget w, XtPointer data, XEvent *event, Boolean *cont);

extern void ensurePulldownColormapInstalled(Widget w, XtPointer clientData, XtPointer callData);
extern void HideFormula(int idx);
extern void ZeigeFormula(int idx);
extern void openMolecule(Widget w, XtPointer data, XtPointer callData);
extern void quit(Widget w, XtPointer data, XtPointer callData);
extern void draw(Widget w, XtPointer data, XtPointer callData);
extern void resize(Widget w, XtPointer data, XtPointer callData);
extern void init(Widget w, XtPointer data, XtPointer callData);
extern void processMenuUse(Widget w, XtPointer clientData, XtPointer callData);
extern void PrintHeadingNumber(MolGLZeichenWindow *ZW, int j1);
extern void PrintHeadingName(MolGLZeichenWindow *ZW, int j1);
extern void ChangeAppearence(Widget w, XtPointer clientData, XtPointer callData);
extern void AppearOkCB(Widget w, XtPointer clientData, XtPointer callData);
extern void AppearCancelCB(Widget w, XtPointer clientData, XtPointer callData);
extern void OkMolnameCB(Widget w, XtPointer clientData, XtPointer callData);
extern void CancelMolnameCB(Widget w, XtPointer clientData, XtPointer callData);
extern void EnterMoleculeName(Widget w, XtPointer clientData, XtPointer callData);
extern void GeheZuMolNr(Widget w, XtPointer clientData, XtPointer callData);
extern void activateMenutg(Widget w, XEvent * event, String * params, Cardinal * num_params);
extern void activateMenu(Widget w, XtPointer clientData, XEvent *event, Boolean *cont);
extern void ChangeGrid(Widget w, XtPointer data, XtPointer callData);
extern void swap(MolGLZeichenWindow *ZW);
extern void ScrBAnpassen(int j1);
extern void SelectOkCB(Widget w, XtPointer data, XtPointer callData);
extern void SelectCancelCB(Widget w, XtPointer data, XtPointer callData);
extern void MolGLGlobals_Init(int argc, char *argv[], STRING_TG &fontfile);
extern void MolGLZeichenWindow_Init( MolGLZeichenWindow *zw );
extern void CreateWidgets();
extern int pickScene(MolGLZeichenWindow *ZW, int x, int y);

/*-- render.h --*/

#define HI_RES_SPHERE  1
#define LO_RES_SPHERE  2
#define DISK           3
#define RAND           4
#define ZYLINDER       5



extern void renderInit(MOLBOOL share);
extern void renderReshape(MolGLZeichenWindow *ZW, int width, int height);
extern void renderMolecule(APPEAR_OPTIONS &appearence, MOLECULE &mol40, MOLBOOL is_selected, float *curquat, 
			   VEK_1D<double> &x3,  VEK_1D<double> &y3,  VEK_1D<double> &z3, 
			   VEK_1D<double> &x32, VEK_1D<double> &y32, VEK_1D<double> &z32,
			   double xsize3, double xmin3, double ysize3, double ymin3, double maxdim3,
			   VEK_1D<short> &color_rad_idx, int sphereVersion,
			   VEK_1D<double> &x2, VEK_1D<double> &y2,
			   MOLBOOL mit_displayliste, int dsplnr, int dsplnr_hi,MOLBOOL &neue_dspl_liste );
extern void renderScene(MolGLZeichenWindow *ZW);
extern void TransfMolKoord(MOLECULE_DATA &moldat, int width, int height);

/*-- trackball.h --*/

void trackball(float q[4], float p1x, float p1y, float p2x, float p2y);
void add_quats(float *q1, float *q2, float *dest);
void build_rotmatrix(float m[4][4], float q[4]);
void axis_to_quat(float a[3], float phi, float q[4]);



