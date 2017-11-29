// render.C


#include "view.h"
#include <math.h>

// #include "glutstroke.h"
// #include "glut_roman.c"
// #include "glut_stroke.c"

GLfloat winWidth, winHeight;

#define ABS(x)      ( ((x) <  0 ) ? (-(x)) : (x) )


static float mat_specular[] = {.72, .8, .93, 1.0};
static float mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};
static float mat_shininess[] = {128.0};

static float light_ambient[] = {0.1, 0.1, 0.1, 1.0};
static float light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
static float light_specular[] = {1.0, 1.0, 1.0, 1.0};
static float light_position[] = {1.0, 1.0, 1.5, 0.0};
static float light0_position[] = {-1.0, -1.0, 1.5, 0.0};

#if 0
GLubyte rasters[24] = {
0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0xc0, 0x00, 
0xc0, 0x00, 0xff, 0x00, 0xff, 0x00, 0xc0, 0x00, 
0xc0, 0x00, 0xc0, 0x00, 0xff, 0xc0, 0xff, 0xc0 };
#endif

void compute_rotation(GLfloat x1, GLfloat y1, GLfloat z1, 
	GLfloat x2, GLfloat y2, GLfloat z2, 
	GLfloat *phi, GLfloat *nx, GLfloat *ny, GLfloat *nz, GLfloat *length, int f_debug);

void save_graph(char *fname, int f_show_labels)
{
	int i, j, k;
	GLfloat units_inv = 1. / 75.;
	GRAPH_MULTI *gm = NULL;
	DATA_3D *d3 = NULL;
	char **s_labels;

	if (!graph_loaded)
		return;
	
	d3 = open_data_3D(GM->nb_V);
	gm = open_graph_multi(GM->nb_V, GM->nb_E, GM->f_multi);
	for (i = 0; i < GM->nb_E; i++) {
		gm->e1[i] = GM->e1[i];
		gm->e2[i] = GM->e2[i];
		gm->mult[i] = GM->mult[i];
		}
	GLfloat m[4][4];

	build_rotmatrix(m, curquat);
	glMultMatrixf(&m[0][0]);

	for(i = 0; i < gm->nb_V; i++) {
		GLfloat x[4], y[4];
		GLboolean b;
		char str[256];
		
		x[0] = (GLfloat) D3->x[i];
		x[1] = (GLfloat) D3->y[i];
		x[2] = (GLfloat) D3->z[i];
		x[3] = 1.;
		// printf("vertex %d: %f %f %f ->", i, x[0], x[1], x[2]);
		
		y[0] = 0.;
		y[1] = 0.;
		y[2] = 0.;
		y[3] = 0.;
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				y[j] += m[k][j] * x[k];
				}
			}
		// printf(" %f %f %f %f\n", y[0], y[1], y[2], y[3]);
		
		d3->x[i] = (float) y[0];
		d3->z[i] = (float) y[1];
		d3->y[i] = (float) y[2];
		
		}
	d3->f_labels = D3->f_labels;
	s_labels = d3->labels;
	d3->labels = D3->labels;
	if (!f_show_labels)
		d3->f_labels = FALSE;
	save_graph_multi_to_file(fname, gm, d3);
	free_graph_multi(gm);
	d3->f_labels = FALSE;
	d3->labels = s_labels;
	free_data_3D(d3);
}

void save_graph_screen(char *fname, int f_show_labels)
{
	int i, j;
	GLfloat units_inv = 2. / STANDARD_SCALE_X;
	GRAPH_MULTI *gm = NULL;
	DATA_3D *d3 = NULL;
	char **s_labels;

	if (!graph_loaded)
		return;
	
	d3 = open_data_3D(GM->nb_V);
	gm = open_graph_multi(GM->nb_V, GM->nb_E, GM->f_multi);
	for (i = 0; i < GM->nb_E; i++) {
		gm->e1[i] = GM->e1[i];
		gm->e2[i] = GM->e2[i];
		gm->mult[i] = GM->mult[i];
		}
	
	glXMakeCurrent(display, glx_data->glxwin, glx_data->cx);
	glXWaitX();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
	
	// render_graph(GM, D3);
	
	init_viewport();
	init_projection();
	init_modelview();

	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
	glEnable(GL_COLOR_MATERIAL);
	glLineWidth(1.5);
	
	GLfloat m[4][4];

	build_rotmatrix(m, curquat);
	glMultMatrixf(&m[0][0]);

	for(i = 0; i < gm->nb_V; i++) {
		GLfloat x1, y1, z1, f[4], fv;
		GLboolean b;
		char str[256];
		
		x1 = (GLfloat) D3->x[i] * units_inv;
		y1 = (GLfloat) D3->y[i] * units_inv;
		z1 = (GLfloat) D3->z[i] * units_inv;
		// printf("vertex %d: %f %f %f ->", i, x1, y1, z1);
		
		glRasterPos3f(x1, y1, z1);
		glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID, &b);
		if (!b) {
			printf("raster position invalid, unable to save!\n");
			return;
			}
		glGetFloatv(GL_CURRENT_RASTER_POSITION, f);
		if (f[3] < .000001) {
			printf("4-th coordinate is zero, unable to save!\n");
			return;
			}
		// fv = 1. / f[3];
		// for (j = 0; j < 3; j++) {
		// 	f[j] *= fv;
		//	}
		// printf("vertex %d: %f %f %f %f\n", i, f[0], f[1], f[2], f[3]);
		// d3->x[i] = (int) f[0];
		// d3->y[i] = (int) f[1];
		// d3->z[i] = (int) f[2];
		
		d3->x[i] = (float) (f[0] * (double) STANDARD_SCALE_X);
		d3->z[i] = (float) (f[1] * (double) STANDARD_SCALE_Y);
		d3->y[i] = 0.0; // (float) (f[2] * (double) STANDARD_SCALE_Z);
		
		// sprintf(str, "%d", i);
		// glutStrokeCharacter(glutBitmapHelvetica10, str[0]);
		// glBitmap(10, 12, 0., 0., 12., 0., rasters);
		}
	d3->f_labels = D3->f_labels;
	s_labels = d3->labels;
	d3->labels = D3->labels;
	if (!f_show_labels)
		d3->f_labels = FALSE;
	save_graph_multi_to_file(fname, gm, d3);
	free_graph_multi(gm);
	d3->f_labels = FALSE;
	d3->labels = s_labels;
	free_data_3D(d3);
	// swap();
}

void renderScene()
{
	glXMakeCurrent(display, glx_data->glxwin, glx_data->cx);
	glXWaitX();
	if (graph_loaded) {
		render_graph(GM, D3);
		}
	else {
		glClearColor(0.6, 0.85, 1.0, 0.0);
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
		}
	swap();
}

void swap()
{
	if (has_double_buffer)
		glXSwapBuffers(display, glx_data->glxwin);
}

void init_viewport()
{
	glViewport(0, 0, (int)winWidth, (int)winHeight);
}

void init_projection()
{
	// GLdouble aspect, left, right, top, bottom, near, far;
	

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	GLdouble aspect = ((GLdouble) winWidth) / ((GLdouble) winHeight);
	gluPerspective( compute_fovy(), aspect, 1.0, 90000000.0);

#if 0
	if (winWidth > winHeight) {
		aspect = ((float) winWidth) / ((float) winHeight);
		left = -1.2;
		right = 1.2;
		bottom = left * aspect;
		top = right * aspect;
		}
	else {
		aspect = ((float) winHeight) / ((float) winWidth);
		bottom = -1.2;
		top = 1.2;
		left = bottom * aspect;
		right = top * aspect;
		}
	near = 3.;
	far = 7;
	printf("left=%f right=%f bottom=%f top=%f\n", 
		left, right, bottom, top);
	glFrustum(left, right, bottom, top, near, far);
#endif

}

void init_modelview()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(appearance.viewpoint[0], 
		appearance.viewpoint[1], 
		appearance.viewpoint[2]);
}

void render_graph(GRAPH_MULTI *gm, DATA_3D *D3)
{
	int e1, e2, i;
	GLfloat units_inv = 2. / STANDARD_SCALE_X;

	// printf("in render_graph()\n"); fflush(stdout);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
	
	init_viewport();
	init_projection();
	init_modelview();
#if 0
	glViewport(0, 0, (int)winWidth, (int)winHeight);
#endif
	// glLoadIdentity();
	
	if (new_display_list) {
		new_display_list = FALSE;
#ifdef VERBOSE
		printf("new_display_list, dsplnr=%d, dsplnr_hi=%d\n", dsplnr, dsplnr_hi);
#endif
		for (int u = 1; u <= 2; u++) {
			if (u == 1)
				glNewList(dsplnr, GL_COMPILE);
			else
				glNewList(dsplnr_hi, GL_COMPILE);	      
			
			glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
			//MyPushMatrix(__LINE__);
			GLfloat radius;
			// int name;
			glEnable(GL_COLOR_MATERIAL);
			glLineWidth(1.5);

			// zeichne die Kanten zwischen Heteroatomen
			for (i = 0; i < gm->nb_E; i++) {
				e1 = gm->e1[i];
				e2 = gm->e2[i];
				GLfloat x1, y1, z1, x2, y2, z2;
				x1 = (GLfloat) D3->x[e1] * units_inv;
				y1 = (GLfloat) D3->y[e1] * units_inv;
				z1 = (GLfloat) D3->z[e1] * units_inv;
				x2 = (GLfloat) D3->x[e2] * units_inv;
				y2 = (GLfloat) D3->y[e2] * units_inv;
				z2 = (GLfloat) D3->z[e2] * units_inv;
				// printf("edge no %d: e1=%d e2=%d\n", i, e1, e2);
				// printf("x1=%f y1=%f z1=%f\n", x1, y1, z1);
				// printf("x2=%f y2=%f z2=%f\n", x2, y2, z2);
				MyPushMatrix(__LINE__);
				// d1 = x1 - x2;
				// d2 = y1 - y2;
				// d3 = z1 - z2;
				// d4 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
				
				if (appearance.show_balls) {
					GLfloat d4;
					int f_debug = 0;
					float phi, nx, ny, nz;
					
					glTranslatef(x2, y2, z2);
					//if (e1 == 7 && e2 == 4)
						f_debug = 0;
					if (f_debug) {
						printf("e1 = %d e2 = %d\n", e1, e2);
						}
					compute_rotation(x1, y1, z1, x2, y2, z2, 
						&phi, &nx, &ny, &nz, &d4, f_debug);
					if (f_debug) {
						printf("nx=%f ny=%f nz=%f\n", nx, ny, nz);
						}
					if (ABS(phi) > 0.001) {
						glRotatef((GLfloat)phi, 
							(GLfloat)nx, 
							(GLfloat)ny, 
							(GLfloat)nz); 
						}
					glColor3ubv(edge_col);
					GLdouble edge_radius = calc_edge_radius();
					GLUquadricObj *quadObj1 = gluNewQuadric();
					gluCylinder(quadObj1, 
						(GLdouble)edge_radius, 
						(GLdouble)edge_radius, 
						(GLdouble)d4, 
						(GLint)10, 
						(GLint)10);
					gluDeleteQuadric(quadObj1);
					}
				else {
					glColor3ubv(black);
					glBegin(GL_LINES);
					glVertex3f((GLfloat)x1, (GLfloat)y1, (GLfloat)z1);
					glVertex3f((GLfloat)x2, (GLfloat)y2, (GLfloat)z2);
					glEnd();
					}
				MyPopMatrix(__LINE__);
				} // next i (next edge) 
	      
			// zeichne die Atome
			for(i = 0; i < gm->nb_V; i++) {
				GLfloat x1, y1, z1;
				x1 = (GLfloat) D3->x[i] * units_inv;
				y1 = (GLfloat) D3->y[i] * units_inv;
				z1 = (GLfloat) D3->z[i] * units_inv;
				

				if (appearance.show_balls) {
					if (is_selected(i))
						glColor3ubv(ball2_col);
					else
						glColor3ubv(ball_col);
					// i2 = color_rad_idx(XX,i1);
					// glLoadName(name);
					MyPushMatrix(__LINE__);
					glTranslatef(x1, y1, z1);
					// printf("ball at x1=%f y1=%f z1=%f\n", x1, y1, z1);
		      
					radius = calc_ball_radius();
					// calc_radius( color_rad_idx(XX,i1), appearence, mol40.atom_anz);
					glScalef( (GLfloat)radius, 
						(GLfloat)radius, 
						(GLfloat)radius);
#if 0
					if((i2==-1)||( ! appearence.show_color))
						glColor3ubv(MGLG->default_color);
					else
						glColor3ubv(MGLG->atomcolor[i2]);
#endif
					if (u == 1)
						glCallList(LO_RES_SPHERE);
					else
						glCallList(HI_RES_SPHERE);
					MyPopMatrix(__LINE__);
					}

				} // next i
			//MyPopMatrix(__LINE__);
			glEndList();
			} // next u
		} // if (new_display_list)
	

	GLfloat m[4][4];

	MyPushMatrix(__LINE__);
	build_rotmatrix(m, curquat);
	glMultMatrixf(&m[0][0]);

	if (sphereVersion == LO_RES_SPHERE)
		glCallList(dsplnr);
	else
		glCallList(dsplnr_hi);
	MyPopMatrix(__LINE__);

}

int is_selected(int v)
{
	int i;
	
	for (i = 0; i < nb_select; i++) {
		if (selected_vertices[i] == v)
			return TRUE;
		}
	return FALSE;
}

void compute_rotation(GLfloat x1, GLfloat y1, GLfloat z1, 
	GLfloat x2, GLfloat y2, GLfloat z2, 
	GLfloat *phi, GLfloat *nx, GLfloat *ny, GLfloat *nz, GLfloat *length, int f_debug)
{
	GLfloat d1, d2, d3, d4, d5, dd;
	
	d1 = x1 - x2;
	d2 = y1 - y2;
	d3 = z1 - z2;
	d4 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
	*length = d4;
	if (f_debug) {
		printf("x1 = %f y1 = %f z1 = %f\n", x1, y1, z1);
		printf("x2 = %f y2 = %f z2 = %f\n", x2, y2, z2);
		printf("%f %f %f\n", d1, d2, d3);
		}
	dd = d1 * d1 + d2 * d2;
	if (dd < 0.001) {
#ifdef VERBOSE
		printf("warning: dd small\n");
#endif
		if (d3 < 0) {
			*phi = 180;
			*nx = 0;
			*ny = 1;
			*nz = 0;
			}
		else {
			*phi = 0;
			*nx = 0;
			*ny = 1;
			*nz = 0;
			}
		}
	else {

#if 0
		d5 = d3 / d4;
		if (d5 > 1.) d5 = 1.;
		if (d5 < -1.) d5 = -1.;
		d5 = acos(d5);
		d5 *= 360.;
		d5 /= 6.28;
		if (d3 < 0)
			d5 = 180. - d5;
		*phi = -d5;
		*nx = d2;
		*ny = -d1;
		*nz = 0;
#else
		d5 = sqrt(d1 * d1 + d2 * d2);
		d5 /= d4;
		if (d5 > 1.) d5 = 1.;
		if (d5 < -1.) d5 = -1.;
		d5 = asin(d5);
		d5 *= 360.;
		d5 /= 6.28;
		if (d3 < 0.0)
			d5 = 180.0 - d5;
		*phi = -d5;
		*nx = d2;
		*ny = -d1;
		*nz = 0;

#endif
		}

	if (f_debug) {
		printf("phi=%f\n", *phi);
		}
	if (*nx * *nx + *ny * *ny + *nz * *nz < 0.001)
		*phi = 0.;
}

void render_init(int share)
{
	if(! share) {
		display_lists_initialised = TRUE; 
		GLUquadricObj *quadObj;

		quadObj = gluNewQuadric();
		gluQuadricDrawStyle(quadObj, (GLenum)GLU_FILL);
		gluQuadricOrientation(quadObj, (GLenum)GLU_OUTSIDE);
		gluQuadricNormals(quadObj, (GLenum)GLU_SMOOTH);
      
		/* hi-detail sphere */
		glNewList(HI_RES_SPHERE, GL_COMPILE);
		gluSphere(quadObj, 1.0, 32, 32);
		glEndList();
      
		/* lo-detail sphere */
		glNewList(LO_RES_SPHERE, GL_COMPILE);
		gluSphere(quadObj, 1.0, 12, 12);
		glEndList();
      
		glNewList(DISK, GL_COMPILE);
		gluDisk(quadObj, 0.0, 1.0, 10, 1);
		glEndList();
      
		glNewList(RAND, GL_COMPILE);
		gluDisk(quadObj, (GLdouble)1.0, (GLdouble)1.2, (GLint)20, (GLint)1);
		//gluDisk(quadObj, 1.0, 1.2, 20, 1);
		glEndList();

		gluDeleteQuadric(quadObj);
		}

	// MOLGEN_FONT_STRUCT &mfs = MGLG->m_fstr;


#if 0
	glXUseXFont(mfs.MyGLFont1->fid, 32, 96, MYGLFONT1+32);
	glXUseXFont(mfs.MyGLFont2->fid, 32, 96, MYGLFONT2+32);
	glXUseXFont(mfs.MyGLFont3->fid, 32, 96, MYGLFONT3+32);
	glXUseXFont(mfs.MyGLFont4->fid, 32, 96, MYGLFONT4+32);
	glXUseXFont(mfs.MyGLFont5->fid, 32, 96, MYGLFONT5+32);
	glXUseXFont(mfs.MyGLFont6->fid, 32, 96, MYGLFONT6+32);
#endif

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.6, 0.85, 1.0, 0.0);
	glClearDepth(1.0);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, light_position);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);




	// glMatrixMode(GL_MODELVIEW);
	// glLoadIdentity();
	// glTranslatef(0, 0, 0);
}

void render_reshape(GLXAREA_DATA *glx_data, int width, int height)
{
	// GLdouble aspect, left, right, top, bottom, d;
	glXMakeCurrent(display, glx_data->glxwin, glx_data->cx);

#ifdef VERBOSE
	printf("in render_reshape()\n"); fflush(stdout);
#endif
	winWidth = width;
	winHeight = height;

#if 0
	glViewport(0, 0, (GLsizei) winWidth, (GLsizei) winHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	aspect = ((GLdouble) winHeight) / (GLdouble) winWidth;
	printf("aspect=%f\n", aspect);
	
	// gluPerspective( compute_fovy(), aspect, 5.0, 90000000.0);
	left = -.9;
	right = .9;
	bottom = left * aspect;
	top = right * aspect;
	printf("left=%f right=%f bottom=%f top=%f\n", 
		left, right, bottom, top);
	if (top < .5) {
		d = 1. / top;
		left *= d;
		right *= d;
		top *= d;
		bottom *= d;
		}
	printf("left=%f right=%f bottom=%f top=%f\n", 
		left, right, bottom, top);
	glFrustum(left, right, bottom, top, 2, 7);
	glMatrixMode(GL_MODELVIEW);
#endif

#if 0
	if (glx_data->mol_loaded)
		TransfMolKoord(glx_data->moldat, width, height);
	glx_data->neue_displayliste = WAHR;
#endif
}



GLdouble compute_fovy()
{
	// zwischen 20 und 50
	return((GLdouble)fovy);
}

GLdouble calc_edge_radius()
{
	GLdouble radius;
	int n, s;

	// printf("in calc_edge_radius()\n"); fflush(stdout);
	if(winHeight <= 100) radius = 2.0;
	else if(winHeight <= 200) radius = 1.9;
	else if(winHeight <= 300) radius = 1.9;
	else if(winHeight <= 400) radius = 1.8;
	else if(winHeight <= 500) radius = 1.7;
	else if(winHeight <= 600) radius = 1.65;
	else if(winHeight <= 700) radius = 1.6;
	else if(winHeight <= 800) radius = 1.5;
	else radius = 1.4;

	radius *= 0.85;
	radius *= 0.01;
	n = GM->nb_V;
	s = 20;
	if (n > 10) {
		radius *= 0.7;
		n -= 10;
		}
	if (n > 15) {
		radius *= 0.7;
		n -= 15;
		}
	if (n > 15) {
		radius *= 0.7;
		n -= 15;
		}
	while (n > 21) {
		if (n > 21) {
			radius *= 0.8;
			n -= s;
			s *= 2;
			// printf("n=%d s=%d\n", n, s); fflush(stdout);
			if (s > 1000)
				s = 1000;
			}
		}
	// printf("radius=%f\n", radius);
	return(radius);
}

GLdouble calc_ball_radius()
{
	GLdouble radius;
	int n, s;

	// printf("in calc_ball_radius()\n"); fflush(stdout);
	if(winHeight <= 100) radius = 2.0;
	else if(winHeight <= 200) radius = 1.9;
	else if(winHeight <= 300) radius = 1.9;
	else if(winHeight <= 400) radius = 1.8;
	else if(winHeight <= 500) radius = 1.7;
	else if(winHeight <= 600) radius = 1.65;
	else if(winHeight <= 700) radius = 1.6;
	else if(winHeight <= 800) radius = 1.5;
	else radius = 1.4;

	radius *= 0.85;
	radius *= 0.09;
	n = GM->nb_V;
	s = 20;
	if (n > 10) {
		radius *= 0.7;
		n -= 10;
		}
	if (n > 15) {
		radius *= 0.7;
		n -= 15;
		}
	if (n > 15) {
		radius *= 0.7;
		n -= 15;
		}
	while (n > 21) {
		if (n > 21) {
			radius *= 0.8;
			n -= s;
			s *= 2;
			// printf("n=%d s=%d\n", n, s); fflush(stdout);
			if (s > 1000)
				s = 1000;
			}
		}
	// printf("radius=%f\n", radius);
	return(radius);
}

void MyPushMatrix(int line)
{
	glPushMatrix();
}

void MyPopMatrix(int line)
{
	glPopMatrix();
}




