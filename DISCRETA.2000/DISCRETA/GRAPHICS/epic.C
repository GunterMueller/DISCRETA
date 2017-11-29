/* epic.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <DISCRETA/graphics.h>

#define ASPECT_RATIO

/*#define DEBUG */

/* die draw-routinen gehen davon aus, 
 * dass epic_draw_fp / epic_draw_dev / epic_draw_epic gesetzt sind. */

/* globals: */
FILE *epic_draw_fp = NIL;
INT epic_draw_dev[4]; /* llx/lly/urx/ury */
INT epic_draw_epic[4]; /* llx/lly/urx/ury */

static INT x_min, x_max, y_min, y_max, f_min_max_set;

static INT transform_llur(INT *ko_dev, INT *ko_out, 
	INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out);
static INT transform_dist(INT *ko_dev, INT *ko_out, 
	INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out);
static INT epic_dev2epic(INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out);
static INT epic_epic2dev(INT ix_out, INT iy_out, INT *ix_dev, INT *iy_dev);
static INT epic_dev2epic_dist(INT dx_dev, INT dy_dev, INT *dx_out, INT *dy_out);
static INT epic_epic2dev_dist(INT dx_out, INT dy_out, INT *dx_dev, INT *dy_dev);
static void ko_min_max(INT x, INT y);
static INT pline(INT *contrl, INT *ptsin);
static INT gtext(INT *contrl, INT *intin, INT *ptsin);
static INT gtext2(BYTE *text, INT x_dev, INT y_dev);
static INT arc(INT *contrl, INT *ptsin, INT *intin);
static INT pie(INT *contrl, INT *ptsin, INT *intin);
static INT pie2(INT x_dev, INT y_dev, INT radius_dev, INT begang, INT endang);
static INT pie_text(INT *contrl, INT *ptsin, INT *intin);
static INT epic_header(FILE *fp, INT w, INT h, INT factor_1000);
static INT epic_footer(FILE *fp, BYTE *caption);
static INT draw_epic_page_file(FILE *fp, BYTE *caption, 
	void *data, VDEVICE *vdev, INT factor_1000, 
	INT f_tape /* TRUE: LwDTapestry2, FALSE: (*draw_func) */ , 
	LW_TAPESTRY_CONTENTS *p, /* only for tape */
	INT x0, INT y0, INT w0, INT h0,  /* only for tape */
	INT i, INT j, INT lines, INT columns,  /* only for tape */
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev));
/* file output of a page;
 * output will be written into open file fp. 
 * can be draw_func or a tapestry.
 * epic_draw_dev[] has to be set 
 * before calling this routine. */
static INT PrintTapestryEPIC(
	LW_TAPESTRY_CONTENTS *lw, INT factor_1000, FILE *fp);

INT EpicDraw(VDEVICE *vdev)
{
	INT *contrl = NIL;
	INT *intin = NIL;
	INT *ptsin = NIL;
	
	if (!VdevGArrays(vdev, &contrl, &intin, &ptsin, NIL, NIL)) {
		Srff("EpicDraw", "VdevGArrays");
		return(FALSE);
		}
	if (contrl[0] == 6) { /* v_pline */
		pline(contrl, ptsin);
		}
	else if (contrl[0] == 8) { /* v_gtext */
		gtext(contrl, intin, ptsin);
		}
	else if (contrl[0] == 9) { /* v_fillarea */
		/* fillarea(contrl, ptsin); */;
		}
	else if (contrl[0] == 11) { /*  */
		if (contrl[5] == 1) /* v_bar(): #ptsin = 2, #intin = 0 */
			/*bar(contrl, ptsin)*/;
		else if (contrl[5] == 2) /* v_arc(): #ptsin = 4, #intin = 2 */
			arc(contrl, ptsin, intin);
		else if (contrl[5] == 3) /* v_pie(): #ptsin = 4, #intin = 2 */
			pie(contrl, ptsin, intin);
		else if (contrl[5] == 4) /* v_circle(): #ptsin = 3, #intin = 0 */
			/*circle(contrl, ptsin)*/;
		else if (contrl[5] == 5) /* v_ellipse(): #ptsin = 2, #intin = 0 */
			/*ellipse(contrl, ptsin)*/;
		else if (contrl[5] == 6) /* v_ellarc(): #ptsin = 2, #intin = 2 */
			/*ellarc(contrl, ptsin, intin)*/;
		else if (contrl[5] == 7) /* v_ellpie(): #ptsin = 2, #intin = 2 */
			/*ellpie(contrl, ptsin, intin)*/;
		else if (contrl[5] == 8) /* v_rbox(): #ptsin = 2, #intin = 0 */
			/*rbox(contrl, ptsin)*/;
		else if (contrl[5] == 9) /* v_rfbox(): #ptsin = 2, #intin = 0 */
			/*rfbox(contrl, ptsin, intin)*/;
		else if (contrl[5] == 11) /* v_pie_text(): #ptsin = 2, #intin = 2 + strlen */
			pie_text(contrl, ptsin, intin);
		}
	else if (contrl[0] == 15) { /* vsl_type */
		/* sl_type(intin) */;
		}
	else if (contrl[0] == 21) { /* vst_font */
		/* st_font(intin) */;
		}
	else if (contrl[0] == 23) { /* vsf_interior */
		/* sf_interior(intin) */;
		}
	else if (contrl[0] == 24) { /* vsf_style */
		/* sf_interior(intin) */;
		}
	else if (contrl[0] == 32) { /* vswr_mode */
		/* swr_mode(intin) */;
		}
	else if (contrl[0] == 39) { /* vst_alignment */
		/* st_alignment(intin) */;
		}
	else if (contrl[0] == 107) { /* vst_point */
		/* st_point(intin) */;
		}
	return(TRUE);
}

static INT transform_llur(INT *ko_dev, INT *ko_out, 
	INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out)
{
	INT dx, dy;
	double a, b;

	dx = ix_dev - ko_dev[0];
	dy = iy_dev - ko_dev[1];
	a = (double) dx / (double)(ko_dev[2] - ko_dev[0]);
	b = (double) dy / (double)(ko_dev[3] - ko_dev[1]);
	dx = (INT)(a * (double)(ko_out[2] - ko_out[0]));
	dy = (INT)(b * (double)(ko_out[3] - ko_out[1]));
	*ix_out = dx + ko_out[0];
	*iy_out = dy + ko_out[1];
	return(TRUE);
}

static INT transform_dist(INT *ko_dev, INT *ko_out, 
	INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out)
{
	INT dx, dy;
	double a, b;

	a = (double) ix_dev / (double)(ko_dev[2] - ko_dev[0]);
	b = (double) iy_dev / (double)(ko_dev[3] - ko_dev[1]);
	dx = (INT)(a * (double) (ko_out[2] - ko_out[0]));
	dy = (INT)(b * (double) (ko_out[3] - ko_out[1]));
	*ix_out = dx;
	*iy_out = dy;
	return(TRUE);
}

#undef DEBUG_DEV2EPIC

static INT f_first_time = TRUE;

void print_4INT(INT *p)
{
	printf("%ld %ld %ld %ld\n", p[0], p[1], p[2], p[3]);
}

static INT epic_dev2epic(INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out)
{
	transform_llur(epic_draw_dev, epic_draw_epic, 
		ix_dev, iy_dev, ix_out, iy_out);
#ifdef DEBUG_DEV2EPIC
	if (f_first_time) {
		print_4INT(epic_draw_dev);
		print_4INT(epic_draw_epic);
		f_first_time = FALSE;
		}
	printf("ix_dev = %ld iy_dev = %ld ix_out = %ld iy_out = %ld\n", 
		ix_dev, iy_dev, *ix_out, *iy_out);
#endif
	return TRUE;
}

static INT epic_epic2dev(INT ix_out, INT iy_out, INT *ix_dev, INT *iy_dev)
{
	transform_llur(epic_draw_epic, epic_draw_dev, 
		ix_out, iy_out, ix_dev, iy_dev);
	return TRUE;
}

static INT epic_dev2epic_dist(INT dx_dev, INT dy_dev, INT *dx_out, INT *dy_out)
{
	transform_llur(epic_draw_dev, epic_draw_epic, 
		dx_dev, dy_dev, dx_out, dy_out);
	return TRUE;
}

static INT epic_epic2dev_dist(INT dx_out, INT dy_out, INT *dx_dev, INT *dy_dev)
{
	transform_dist(epic_draw_epic, epic_draw_dev, 
		dx_out, dy_out, dx_dev, dy_dev);
	return TRUE;
}

static void ko_min_max(INT x, INT y)
{
	if (!f_min_max_set) {
		x_min = x_max = x;
		y_min = y_max = y;
		}
	else {
		x_min = MINIMUM(x_min, x);
		y_min = MINIMUM(y_min, y);
		x_max = MAXIMUM(x_max, x);
		y_max = MAXIMUM(y_max, y);
		}
	f_min_max_set = TRUE;
}

static INT pline(INT *contrl, INT *ptsin)
{
	INT x_dev, y_dev, x_out, y_out;
	INT i, i2;
	BYTE str[256];

	fputs("\\put(0, 0){\\drawline", epic_draw_fp);
	for (i = 0; i < contrl[1]; i++) {
		i2 = i << 1;
		x_dev = (INT) ptsin[i2];
		y_dev = (INT) ptsin[i2 + 1];
		ko_min_max(x_dev, y_dev);
		epic_dev2epic(x_dev, y_dev, &x_out, &y_out);
		sprintf(str, "(%ld, %ld)", x_out, y_out);
		fputs(str, epic_draw_fp);
		}
	fputs("}\n", epic_draw_fp);
	return(TRUE);
}

static INT gtext(INT *contrl, INT *intin, INT *ptsin)
{
	INT str_len, i, len;
	INT x_dev, y_dev;
	BYTE str[10000], c;
	
	x_dev = (INT)ptsin[0];
	y_dev = (INT)ptsin[1];
	len = (INT) contrl[3]; /* Anzahl Zeichen ohne Nullbyte */
	for (i = 0; i < len; i++) {
		c = (BYTE)intin[i];
		str[i] = c;
		}
	str[i] = 0;
	gtext2(str, x_dev, y_dev);
	return OK;
}

static INT gtext2(BYTE *text, INT x_dev, INT y_dev)
{
	INT str_len, i, len;
	INT x_out, y_out;
	BYTE str[10000], str1[10000], c;
	
	ko_min_max(x_dev, y_dev);
	epic_dev2epic(x_dev, y_dev, &x_out, &y_out);
	
	len = strlen(text);
	if (len <= 0)
		return(TRUE);
	str_len = 0;
	if (len) {
		c = text[0];
#if 0
		if (c != '$') {
			for (i = 0; i < len; i++) {
				c = (BYTE)intin[i];
				if (c == '#') {
					str[str_len++] = '\\';
					str[str_len++] = c;
					}
				else if (c == '^') {
					str[str_len++] = '\\';
					str[str_len++] = c;
					}
				else if (c == '\\') {
					str[str_len++] = '\\';
					str[str_len++] = c;
					}
				else
					str[str_len++] = c;
				}
			}
		else {
#endif
			for (i = 0; i < len; i++) {
				c = (BYTE)text[i];
				str[str_len++] = c;
				}
#if 0
			}
#endif
		}
	str[str_len] = 0;
	sprintf(str1, "\\put(%ld,%ld){%s}\n", x_out, y_out, str);
	fputs(str1, epic_draw_fp);

	return(TRUE);
}

static INT arc(INT *contrl, INT *ptsin, INT *intin)
{
	INT x_dev, y_dev, x_out, y_out;
	INT radius_dev, radius_out, tmp, begang, endang;
	BYTE str[256];
	
	x_dev = (INT) ptsin[0];
	y_dev = (INT) ptsin[1];
	ko_min_max(x_dev, y_dev);
	epic_dev2epic(x_dev, y_dev, &x_out, &y_out);
	
	radius_dev = (INT)ptsin[6];
	epic_dev2epic_dist(radius_dev, 0, &radius_out, &tmp);
	if (radius_out <= 0)
		radius_out = 1;

	begang = (INT) intin[0];
	endang = (INT) intin[1];
	begang /= 10;
	endang /= 10;
	sprintf(str, "\\put(%ld,%ld){\\circle{%ld}}\n", x_out, y_out, radius_out);
	fputs(str, epic_draw_fp);
	return(TRUE);
}

static INT pie(INT *contrl, INT *ptsin, INT *intin)
{
	INT x_dev, y_dev;
	INT radius_dev, begang, endang;
	
	x_dev = (INT)ptsin[0];
	y_dev = (INT)ptsin[1];
	radius_dev = (INT)ptsin[6];
	begang = (INT)intin[0];
	endang = (INT)intin[1];
	pie2(x_dev, y_dev, radius_dev, begang, endang);
	
	return(TRUE);
}

static INT pie2(INT x_dev, INT y_dev, INT radius_dev, INT begang, INT endang)
{
	INT x_out, y_out;
	INT radius_out, tmp;
	BYTE str[256];
	
	ko_min_max(x_dev, y_dev);
	epic_dev2epic(x_dev, y_dev, &x_out, &y_out);
	
	epic_dev2epic_dist(radius_dev, 0, &radius_out, &tmp);
	if (radius_out <= 0)
		radius_out = 1;

	begang /= 10;
	endang /= 10;
	/* !!! filled ellipse statt pie */
	sprintf(str, "\\put(%ld,%ld){\\circle*{%ld}}\n", x_out, y_out, radius_out);
	fputs(str, epic_draw_fp);
	return(TRUE);
}

static INT pie_text(INT *contrl, INT *ptsin, INT *intin)
{
	BYTE text[10000], c;
	INT len, i;
	INT x_dev, y_dev, radius_dev, begang, endang;
	
	len = (INT)contrl[3] - 2; /* strlen = Anzahl Zeichen ohne Nullbyte */
	for (i = 0; i < len; i++) {
		text[i] = (BYTE) intin[2 + i];
		}
	text[i] = 0;
	// printf("epic: pie_text, text = %s\n", text);
	// fflush(stdout);
	x_dev = (INT) ptsin[0];
	y_dev = (INT) ptsin[1];
	gtext2(text, x_dev, y_dev);
	
	radius_dev = (INT) ptsin[6];
	begang = (INT) intin[0];
	endang = (INT) intin[1];
	pie2(x_dev, y_dev, radius_dev, begang, endang);
	return OK;
}

static INT epic_header(FILE *fp, INT w, INT h, INT factor_1000)
{
	BYTE str[256];
	double d, d_w;
	
	f_min_max_set = FALSE;
	d = (double) factor_1000 * 0.001  * 0.1;
	d_w = 0.1 * d * (double) w + 200.;
		// 1600 * 0.1 = 160 mm width
	fprintf(fp, "\\begin{figure}[ht]\n");
	fprintf(fp, "\\begin{center}\n");
	fprintf(fp, "{\n");
	// fprintf(fp, "\\fbox{\n");
	// fprintf(fp, "\\begin{minipage}[t]{%lfmm}\n", d_w);
	fprintf(fp, "\\unitlength%lfmm\n", d);
	sprintf(str, "\\begin{picture}(%ld, %ld)\n", w, h);
	fputs(str, fp);
	return OK;
}

#if 0
\begin{figure}[ht]
\begin{center}
\fbox{
\begin{minipage}[t]{8cm}
#endif

static INT epic_footer(FILE *fp, BYTE *caption)
{
	fputs("\\end{picture}\n", fp);
	fprintf(fp, "\\caption{%s}\n", caption);
	// fprintf(fp, "\\end{minipage}\n");
	// fprintf(fp, "}\n");
	fprintf(fp, "}\n");
	fprintf(fp, "\\end{center}\n");
	fprintf(fp, "\\end{figure}\n");
	return OK;
}

#if 0
\caption{\label{g4} the graphs on 4 points}
\end{minipage}
}
\end{center}
\end{figure}
#endif

static INT draw_epic_page_file(FILE *fp, BYTE *caption, 
	void *data, VDEVICE *vdev, INT factor_1000, 
	INT f_tape /* TRUE: LwDTapestry2, FALSE: (*draw_func) */ , 
	LW_TAPESTRY_CONTENTS *p, /* only for tape */
	INT x0, INT y0, INT w0, INT h0,  /* only for tape */
	INT i, INT j, INT lines, INT columns,  /* only for tape */
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev))
/* file output of a page;
 * output will be written into open file fp. 
 * can be draw_func or a tapestry.
 * epic_draw_dev[] has to be set 
 * before calling this routine. */
{
	INT s_co[4], s_type;
	INT ix_dev, iy_dev;
	BYTE *tmp_name = "/tmp/tmp.tmp";
	FILE *tmp_fp;
	INT dx, dy, Dx, Dy;
	double ar;
	
	VdevI2(vdev, gl_vdi_contrl, gl_vdi_intin, gl_vdi_ptsin, 
		gl_vdi_intout, gl_vdi_ptsout);
	s_type = vdev->type;
	s_co[0] = vdev->co[0];
	s_co[1] = vdev->co[1];
	s_co[2] = vdev->co[2];
	s_co[3] = vdev->co[3];
	vdev->type = 3; /* Epic output */
	vdev->co[0] = 0;
	vdev->co[1] = 0;
	vdev->co[2] = 1000000; // pix_x_max
	vdev->co[3] = 1000000; // pix_y_max

	/* 
	 * the first drawing goes into a temporary file; 
	 * we only want to know the extrema 
	 */
	if ((tmp_fp = fopen(tmp_name, "w")) == NIL) {
		Srff("draw_epic_page_file", "fopen");
		return ERROR;
		}
	epic_header(tmp_fp, 1400, 1600, factor_1000);
	
	epic_draw_fp = tmp_fp;
	epic_draw_epic[0] = 0; // x_min (llx)
	epic_draw_epic[1] = 0; // y_min (lly)
	epic_draw_epic[2] = 1400; // x_max (urx)
	epic_draw_epic[3] = 1600; // y_max (ury)
	Vswr_mode(vdev, 1);
	Vsf_color(vdev, 0);
	Vsf_interior(vdev, 0);
	Vsf_perimeter(vdev, 0);
	if (f_tape) {
		if (!LwDTapestry2(p, vdev, x0, y0, w0, h0, i, j, 
			lines, columns, 0, lines * columns)) {
			Srff("draw_epic_page_file", "LwDTapestry2");
			return FALSE;
			}
		}
	else {
		(*draw_func)(data, 0, 0, 0, 0, vdev);
		}
	Vswr_mode(vdev, 1);
	epic_footer(tmp_fp, caption);
	if (fflush(tmp_fp) != NIL) {
		Srff("draw_epic_file", "fflush");
		return ERROR;
		}
	if (fclose(tmp_fp) != NIL) {
		Srff("draw_epic_file", "fclose");
		return ERROR;
		}

	printf("epic_draw_dev = [%ld,%ld,%ld,%ld]\n", 
		epic_draw_dev[0], 
		epic_draw_dev[1], 
		epic_draw_dev[2], 
		epic_draw_dev[3]);
	printf("x_min = %ld x_max = %ld y_min = %ld y_max = %ld\n", x_min, x_max, y_min, y_max);

#ifdef ASPECT_RATIO
	dx = x_max - x_min;
	dy = y_max - y_min;
	ar = (double) dy / (double) dx;
	Dx = 1600;
	Dy = (INT)(((double) Dx) * ar);
	printf("epic: ASPECT_RATIO dx = %ld dy = %ld, ar = %lf, Dx = %ld, Dy = %ld\n", 
		dx, dy, ar, Dx, Dy);
	epic_draw_epic[0] = 0; // x_min (llx)
	epic_draw_epic[1] = 0; // y_min (lly)
	epic_draw_epic[2] = Dx; // x_max (urx)
	epic_draw_epic[3] = Dy; // y_max (ury)
#else
	Dx = 1600;
	Dy = 1400;
#endif
	epic_draw_dev[0] = x_min; // llx
	epic_draw_dev[1] = y_max; // lly
	epic_draw_dev[2] = x_max; // urx
	epic_draw_dev[3] = y_min; // ury

	printf("epic_draw_dev = [%ld,%ld,%ld,%ld]\n", 
		epic_draw_dev[0], 
		epic_draw_dev[1], 
		epic_draw_dev[2], 
		epic_draw_dev[3]);
	printf("epic_draw_epic = [%ld,%ld,%ld,%ld]\n", 
		epic_draw_epic[0], 
		epic_draw_epic[1], 
		epic_draw_epic[2], 
		epic_draw_epic[3]);
	/* 
	 * now we draw into the given file fp 
	 */
	// epic_header(fp, x_max - x_min, y_max - y_min, factor_1000);
	epic_header(fp, Dx, Dy, factor_1000);
	epic_draw_fp = fp;
	Vswr_mode(vdev, 1);
	Vsf_color(vdev, 0);
	Vsf_interior(vdev, 0);
	Vsf_perimeter(vdev, 0);
	if (f_tape) {
		if (!LwDTapestry2(p, vdev, x0, y0, w0, h0, i, j, 
			lines, columns, 0, lines * columns)) {
			Srff("draw_epic_page_file", "LwDTapestry2");
			return FALSE;
			}
		}
	else {
		(*draw_func)(data, 0, 0, 0, 0, vdev);
		}
	Vswr_mode(vdev, 1);
	epic_footer(fp, caption);

	vdev->type = s_type;
	vdev->co[0] = s_co[0];
	vdev->co[1] = s_co[1];
	vdev->co[2] = s_co[2];
	vdev->co[3] = s_co[3];
	return OK;
}

INT draw_epic_file(void *data, VDEVICE *vdev, BYTE *name, 
	INT factor_1000, BYTE *caption, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev))
{
	FILE *fp;
	
	if ((fp = fopen(name, "w")) == NIL) {
		Srff("draw_epic_file", "fopen");
		return ERROR;
		}

	// INT epic_draw_dev[4]; /* llx/lly/urx/ury */
	// the drawing is in typical X-windows coordinates,
	// the origin (0,0) lies at the {\em top left} corner of the window !
	epic_draw_dev[0] = 0;
	epic_draw_dev[1] = 1000000;
	epic_draw_dev[2] = 1000000;
	epic_draw_dev[3] = 0;
	draw_epic_page_file(fp, caption, 
		data, vdev, factor_1000, 
		FALSE /* f_tape */, 
		NIL, /* only for tape */
		0, 0, 0, 0,  /* only for tape */
		0, 0, 0, 0,  /* only for tape */
		draw_func);
	
	if (fflush(fp) != NIL) {
		Srff("draw_epic_file", "fflush");
		return ERROR;
		}
	if (fclose(fp) != NIL) {
		Srff("draw_epic_file", "fclose");
		return ERROR;
		}

	printf("draw_epic_file()|written file %s of size %ld\n", 
		name, file_size(name));
	return OK;
}

#define TAPESTRY_PRINT_PAGE_COLUMNS 2
#define TAPESTRY_PRINT_PAGE_LINES 3
/* same as in ps.C */

INT draw_tape_EPIC(LW_TAPESTRY_CONTENTS *p, 
	INT factor_1000, BYTE *name)
{
	FILE *fp;

	if ((fp = fopen(name, "w")) == NIL) {
		Srff("draw_tape_EPIC", "fopen");
		return(FALSE);
		}
	PrintTapestryEPIC(p, factor_1000, fp);
	if (fflush(fp) != NIL) {
		Srff("draw_tape_EPIC", "fflush");
		return(FALSE);
		}
	if (fclose(fp) != NIL) {
		Srff("draw_tape_EPIC", "fclose");
		return(FALSE);
		}
	printf("draw_tape_EPIC()|written file %s of size %ld\n", 
		name, file_size(name));
	return TRUE;
}

static INT PrintTapestryEPIC(
	LW_TAPESTRY_CONTENTS *lw, INT factor_1000, FILE *fp)
{
	VDEVICE vdev_, *vdev = &vdev_;
	BYTE str[1024];
	INT pics_per_page, page, pages, cur_line, old_type;
	INT pageHeight, pageWidth;
	INT x0, y0, w0, h0;
	BYTE caption[256];
	
	VdevI2(vdev, gl_vdi_contrl, gl_vdi_intin, gl_vdi_ptsin, 
		gl_vdi_intout, gl_vdi_ptsout);
	old_type = vdev->type;
	vdev->type = 0; /* PS output */
	vdev->co[0] = 0;
	vdev->co[1] = 0;
	vdev->co[2] = lw->picture_width * TAPESTRY_PRINT_PAGE_COLUMNS;
	vdev->co[3] = lw->picture_height * TAPESTRY_PRINT_PAGE_LINES;

	epic_draw_dev[0] = 0;
	epic_draw_dev[1] = 0;
	epic_draw_dev[2] = lw->picture_width * TAPESTRY_PRINT_PAGE_COLUMNS;
	epic_draw_dev[3] = lw->picture_height * TAPESTRY_PRINT_PAGE_LINES;

	pics_per_page = TAPESTRY_PRINT_PAGE_LINES * TAPESTRY_PRINT_PAGE_COLUMNS;
	pages = lw->len / pics_per_page;
	if (pages * pics_per_page < lw->len) 
		pages++;
	for (page = 0; page < pages; page++) {
		printf("page %ld / %ld\n", page, pages);
		cur_line = page * TAPESTRY_PRINT_PAGE_LINES;
		x0 = 0;
		y0 = 0;
		w0 = lw->picture_width;
		h0 = lw->picture_height;
		sprintf(caption, "%ld / %ld", page + 1, pages);
		draw_epic_page_file(fp, caption, 
			NIL /* data */, vdev, factor_1000, 
			TRUE /* f_tape */, 
			lw, 
			x0, y0, w0, h0, cur_line, 0L, 
			TAPESTRY_PRINT_PAGE_LINES, 
			TAPESTRY_PRINT_PAGE_COLUMNS, 
			NIL /* draw_func */);
		fputs("\n\\clearpage\n\n", fp);
		}
	vdev->type = old_type;
	return TRUE;
}

#endif /* GRAPHICS_TRUE */


