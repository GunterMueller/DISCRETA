/* mp.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1997
 */


#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <stdlib.h>

#include <DISCRETA/graphics.h>
#include <DISCRETA/divs.h> // for STRING_OP

#define ASPECT_RATIO

/*#define DEBUG */

/* die draw-routinen gehen davon aus, 
 * dass mp_draw_fp / mp_draw_dev / mp_draw_mp gesetzt sind. */

/* globals: */
FILE *mp_draw_fp = NIL;
INT mp_draw_dev[4]; /* llx/lly/urx/ury */
INT mp_draw_mp[4]; /* llx/lly/urx/ury */

static INT x_min, x_max, y_min, y_max, f_min_max_set;

static INT txt_halign = 0; 
/* 0: links, 1: zentriert, 2: rechts */
static INT txt_valign = 0;
/* 0: unten, 1:mittig, 2: oben (haengender Text) */
static INT txt_boxed = 0;
static INT txt_overwrite = 0;

static INT line_beg_style = 0;
static INT line_end_style = 0; // 1: arrow

static INT fill_interior = 0; // in 1/100th, 0= none (used for pie-drawing)
static INT fill_color = 0; // 0 = white, 1 = black
static INT fill_shape = 1; // 0 =  .., 1 = -- 
static INT fill_outline = 0; 
static INT fill_nofill = 0; 


static INT line_dashing = 0; // 0 = no dashing, otherwise scaling factor 1/100th evenly

static INT transform_llur(INT *ko_dev, INT *ko_out, 
	INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out);
static INT transform_dist(INT *ko_dev, INT *ko_out, 
	INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out);
static INT mp_dev2mp(INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out);
static INT mp_mp2dev(INT ix_out, INT iy_out, INT *ix_dev, INT *iy_dev);
static INT mp_dev2mp_dist(INT dx_dev, INT dy_dev, INT *dx_out, INT *dy_out);
static INT mp_mp2dev_dist(INT dx_out, INT dy_out, INT *dx_dev, INT *dy_dev);
static void ko_min_max(INT x, INT y);
static INT pline(INT *contrl, INT *ptsin);
static INT bezier(INT *contrl, INT *ptsin);
static INT pline2(INT *contrl, INT *ptsin, BYTE *symbol);
static INT gtext(INT *contrl, INT *intin, INT *ptsin);
static INT fillarea(INT *contrl, INT *ptsin);
static INT arc(INT *contrl, INT *ptsin, INT *intin);
static INT pie(INT *contrl, INT *ptsin, INT *intin);
static INT pie_text(INT *contrl, INT *ptsin, INT *intin);
static INT st_alignment(INT *intin);
static void get_alignment(BYTE *align);
static INT sl_udsty(INT *intin);
static INT sl_ends(INT *intin);
static INT sf_interior(INT *intin);
static INT sf_color(INT *intin);
static INT sf_shape(INT *intin);
static INT sf_outline(INT *intin);
static INT sf_nofill(INT *intin);
static INT st_boxed(INT *intin);
static INT st_overwrite(INT *intin);
static BYTE *get_label(INT x, INT y);
static void mp_free_local_data();
static INT draw_boxes();

INT mp_draw(VDEVICE *vdev)
{
	INT *contrl = NIL;
	INT *intin = NIL;
	INT *ptsin = NIL;
	
	if (!VdevGArrays(vdev, &contrl, &intin, &ptsin, NIL, NIL)) {
		error("mp_draw() error in VdevGArrays");
		return(FALSE);
		}
	if (contrl[0] == 6) { /* v_pline */
		pline(contrl, ptsin);
		}
	else if (contrl[0] == 7) { /* v_pline */
		bezier(contrl, ptsin);
		}
	else if (contrl[0] == 8) { /* v_gtext */
		gtext(contrl, intin, ptsin);
		}
	else if (contrl[0] == 9) { /* v_fillarea */
		fillarea(contrl, ptsin);
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
		sf_interior(intin);
		}
	else if (contrl[0] == 24) { /* vsf_style */
		/* sf_interior(intin) */;
		}
	else if (contrl[0] == 25) { /* vsf_color */
		sf_color(intin);
		}
	else if (contrl[0] == 26) { /* vsf_shape */
		sf_shape(intin);
		}
	else if (contrl[0] == 27) { /* vsf_outline */
		sf_outline(intin);
		}
	else if (contrl[0] == 28) { /* vsf_nofill */
		sf_nofill(intin);
		}
	else if (contrl[0] == 32) { /* vswr_mode */
		/* swr_mode(intin) */;
		}
	else if (contrl[0] == 39) { /* vst_alignment */
		st_alignment(intin);
		}
	else if (contrl[0] == 107) { /* vst_point */
		/* st_point(intin) */;
		}
	else if (contrl[0] == 108) { /* vsl_ends */
		sl_ends(intin);
		}
	else if (contrl[0] == 110) { /* vst_boxed */
		st_boxed(intin);
		}
	else if (contrl[0] == 111) { /* vst_overwrite */
		st_overwrite(intin);
		}
	else if (contrl[0] == 113) { /* vsl_udsty */
		sl_udsty(intin);
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

static INT mp_dev2mp(INT ix_dev, INT iy_dev, INT *ix_out, INT *iy_out)
{
	transform_llur(mp_draw_dev, mp_draw_mp, 
		ix_dev, iy_dev, ix_out, iy_out);
	return TRUE;
}

static INT mp_mp2dev(INT ix_out, INT iy_out, INT *ix_dev, INT *iy_dev)
{
	transform_llur(mp_draw_mp, mp_draw_dev, 
		ix_out, iy_out, ix_dev, iy_dev);
	return TRUE;
}

static INT mp_dev2mp_dist(INT dx_dev, INT dy_dev, INT *dx_out, INT *dy_out)
{
	transform_llur(mp_draw_dev, mp_draw_mp, 
		dx_dev, dy_dev, dx_out, dy_out);
	return TRUE;
}

static INT mp_mp2dev_dist(INT dx_out, INT dy_out, INT *dx_dev, INT *dy_dev)
{
	transform_dist(mp_draw_mp, mp_draw_dev, 
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
	return pline2(contrl, ptsin, " -- ");
}

static INT bezier(INT *contrl, INT *ptsin)
{
	return pline2(contrl, ptsin, " .. ");
}

static INT pline2(INT *contrl, INT *ptsin, BYTE *symbol)
{
	INT x_dev, y_dev, x_out, y_out;
	INT i, i2;

	if (line_end_style == 1)
		fprintf(mp_draw_fp, "drawarrow ");
	else
		fprintf(mp_draw_fp, "draw ");
	for (i = 0; i < contrl[1]; i++) {
		i2 = i << 1;
		x_dev = (INT) ptsin[i2];
		y_dev = (INT) ptsin[i2 + 1];
		ko_min_max(x_dev, y_dev);
		mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
		// printf("%ld %ld -> %ld %ld\n", x_dev, y_dev, x_out, y_out);
		// fflush(stdout);
		if (i) {
			fprintf(mp_draw_fp, symbol);
			}
		fprintf(mp_draw_fp, "(%ldu,%ldu)", x_out, y_out);

#if 0		
		p = get_label(x_dev, y_dev);
		strcpy(label2, p);
		ko_min_max(x_dev, y_dev);
		mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
		fprintf(mp_draw_fp, "%s.c = (%ldu,%ldu);\n", label2, x_out, y_out);
		if (i) {
			fprintf(mp_draw_fp, "cuta(%s, %s) %s.c--%s.c;\n", 
				label1, label2, label1, label2);
			}
		strcpy(label1, label2);
#endif
		}
	if (line_dashing) {
		fprintf(mp_draw_fp, " dashed evenly");
		if (line_dashing != 100) {
			fprintf(mp_draw_fp, " scaled %g", (double) line_dashing / 100.);
			}
		}
	fprintf(mp_draw_fp, ";\n");
	// fflush(mp_draw_fp);
	return(TRUE);
}

static INT gtext(INT *contrl, INT *intin, INT *ptsin)
{
	INT str_len, i, len;
	INT x_dev, y_dev, x_out, y_out;
	BYTE str[10000], c;
	BYTE align[64], *lab;
	
	x_dev = (INT)ptsin[0];
	y_dev = (INT)ptsin[1];
	mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
	
	len = (INT)contrl[3]; /* Anzahl Zeichen ohne Nullbyte */
	if (len <= 0)
		return(TRUE);
	ko_min_max(x_dev, y_dev);
	str_len = 0;
	if (len) {
		for (i = 0; i < len; i++) {
			c = (BYTE)intin[i];
			str[str_len++] = c;
			}
		}
	str[str_len] = 0;
	get_alignment(align);
	if (txt_boxed) {
		lab = get_label(x_out, y_out);
		fprintf(mp_draw_fp, "boxit.%s(btex %s etex);\n", lab, str);
		fprintf(mp_draw_fp, "%s.c=(%ldu,%ldu);\n", lab, x_out, y_out);
		if (txt_overwrite) {
			fprintf(mp_draw_fp, "unfill bpath %s;\n", lab);
			}
		fprintf(mp_draw_fp, "drawboxed(%s);\n", lab);
		
		}
	else {
		fprintf(mp_draw_fp, "label%s(btex %s etex, (%ldu,%ldu));\n", 
			align, str, x_out, y_out);
		}

	return(TRUE);
}

static INT fillarea(INT *contrl, INT *ptsin)
{
	INT x_dev, y_dev, x_out, y_out;
	INT i, i2;
	BYTE *shape;
	
	if (fill_shape == 0)
		shape = " .. ";
	else
		shape = " -- ";
	fprintf(mp_draw_fp, "path pp;\n");
	fprintf(mp_draw_fp, "pp = ");
	for (i = 0; i < contrl[1]; i++) {
		i2 = i << 1;
		x_dev = (INT) ptsin[i2];
		y_dev = (INT) ptsin[i2 + 1];
		ko_min_max(x_dev, y_dev);
		mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
		if (i) {
			fprintf(mp_draw_fp, shape);
			}
		fprintf(mp_draw_fp, "(%ldu,%ldu)", x_out, y_out);
		}
	fprintf(mp_draw_fp, shape);
	fprintf(mp_draw_fp, " cycle;\n");
	if (!fill_nofill) {
		fprintf(mp_draw_fp, "fill pp withcolor ");
		if (fill_interior > 99) {
			fprintf(mp_draw_fp, "1 ");
			}
		else
			fprintf(mp_draw_fp, ".%02ld ", fill_interior);
		if (fill_color == 1)
			fprintf(mp_draw_fp, "black");
		else
			fprintf(mp_draw_fp, "white");
		fprintf(mp_draw_fp, ";\n");
		}
	if (fill_outline) {
		fprintf(mp_draw_fp, "draw pp ");
		if (line_dashing) {
			fprintf(mp_draw_fp, " dashed evenly ");
			if (line_dashing != 100) {
				fprintf(mp_draw_fp, " scaled %g", (double) line_dashing / 100.);
				}
			}
		fprintf(mp_draw_fp, ";\n");
		}
	return(TRUE);
}


static INT arc(INT *contrl, INT *ptsin, INT *intin)
{
	INT x_dev, y_dev, x_out, y_out, i;
	INT x[10], y[10];
	INT radius_dev, radius_out, tmp, begang, endang;
	double t1, t2;
	
	x_dev = (INT) ptsin[0];
	y_dev = (INT) ptsin[1];
	ko_min_max(x_dev, y_dev);
	mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
	
	radius_dev = (INT)ptsin[6];
	mp_dev2mp_dist(radius_dev, 0, &radius_out, &tmp);
	if (radius_out <= 0)
		radius_out = 1;

	begang = (INT) intin[0];
	endang = (INT) intin[1];

	x[0] = x_out +  radius_out;
	y[0] = y_out;
	x[1] = x_out;
	y[1] = y_out +  radius_out;
	x[2] = x_out -  radius_out;
	y[2] = y_out;
	x[3] = x_out;
	y[3] = y_out -  radius_out;
	x[4] = x_out +  radius_out;
	y[4] = y_out;
#if 0
	fprintf(mp_draw_fp, "draw ");
	for (i = 0; i < 5; i++) {
		if (i) {
			fprintf(mp_draw_fp, " .. ");
			}
		fprintf(mp_draw_fp, "(%ldu,%ldu)", x[i], y[i]);
		}
	if (line_dashing) {
		fprintf(mp_draw_fp, " dashed evenly");
		if (line_dashing != 100) {
			fprintf(mp_draw_fp, " scaled %lf", (double) line_dashing / 100.);
			}
		}
	fprintf(mp_draw_fp, ";\n");
#endif
	fprintf(mp_draw_fp, "path p; p = ");
	for (i = 0; i < 5; i++) {
		if (i) {
			fprintf(mp_draw_fp, " .. ");
			}
		fprintf(mp_draw_fp, "(%ldu,%ldu)", x[i], y[i]);
		}
	fprintf(mp_draw_fp, ";\n");
	t1 = ((double) begang * 4) / 3600.;
	t2 = ((double) endang * 4) / 3600.;
	// printf("mp: arc begang = %ld endang = %ld t1 = %lf t2 = %lf\n", begang, endang, t1, t2);
	fprintf(mp_draw_fp, "path q; q = subpath (%g,%g) of p;\n", t1, t2);
	fprintf(mp_draw_fp, "draw q ");
	if (line_dashing) {
		fprintf(mp_draw_fp, " dashed evenly");
		if (line_dashing != 100) {
			fprintf(mp_draw_fp, " scaled %g", (double) line_dashing / 100.);
			}
		}
	fprintf(mp_draw_fp, ";\n");
	
	
#if 0
	begang = (INT) intin[0];
	endang = (INT) intin[1];
	begang /= 10;
	endang /= 10;
	sprintf(str, "\\put(%ld,%ld){\\circle{%ld}}\n", x_out, y_out, radius_out);
	fputs(str, mp_draw_fp);
#endif

	return(TRUE);
}

static INT pie(INT *contrl, INT *ptsin, INT *intin)
{
	INT x_dev, y_dev, x_out, y_out, i;
	INT x[10], y[10];
	INT radius_dev, radius_out, tmp;
	
	x_dev = (INT) ptsin[0];
	y_dev = (INT) ptsin[1];
	ko_min_max(x_dev, y_dev);
	mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
	
	radius_dev = (INT)ptsin[6];
	mp_dev2mp_dist(radius_dev, 0, &radius_out, &tmp);
	if (radius_out <= 0)
		radius_out = 1;
	x[0] = x_out +  radius_out;
	y[0] = y_out;
	x[1] = x_out;
	y[1] = y_out +  radius_out;
	x[2] = x_out -  radius_out;
	y[2] = y_out;
	x[3] = x_out;
	y[3] = y_out -  radius_out;
	x[4] = x_out +  radius_out;
	y[4] = y_out;
	fprintf(mp_draw_fp, "path pp;\n");
	fprintf(mp_draw_fp, "pp = ");
	for (i = 0; i < 5; i++) {
		if (i) {
			fprintf(mp_draw_fp, " .. ");
			}
		fprintf(mp_draw_fp, "(%ldu,%ldu)", x[i], y[i]);
		}
	fprintf(mp_draw_fp, " .. cycle;\n");
	fprintf(mp_draw_fp, "fill pp withcolor ");
	if (fill_interior > 99) {
		fprintf(mp_draw_fp, "1 ");
		}
	else
		fprintf(mp_draw_fp, ".%02ld ", fill_interior);
	if (fill_color == 1)
		fprintf(mp_draw_fp, "black");
	else
		fprintf(mp_draw_fp, "white");
	fprintf(mp_draw_fp, ";\n");
#if 0
	begang = (INT) intin[0];
	endang = (INT) intin[1];
	begang /= 10;
	endang /= 10;
	sprintf(str, "\\put(%ld,%ld){\\circle{%ld}}\n", x_out, y_out, radius_out);
	fputs(str, mp_draw_fp);
#endif

	return(TRUE);
#if 0
	INT x_dev, y_dev, x_out, y_out;
	INT radius_dev, radius_out, tmp, begang, endang;
	BYTE str[256], *s;
	
	x_dev = (INT) ptsin[0];
	y_dev = (INT) ptsin[1];
	s = get_label(x_dev, y_dev);
	strcpy(str, s);
	printf("pie %s\n", str);
	
	ko_min_max(x_dev, y_dev);
	mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
	
	radius_dev = (INT) ptsin[6];
	mp_dev2mp_dist(radius_dev, 0, &radius_out, &tmp);
	if (radius_out <= 0)
		radius_out = 1;

	begang = (INT) intin[0];
	endang = (INT) intin[1];
	begang /= 10;
	endang /= 10;
	/* !!! filled ellipse statt pie */

	fprintf(mp_draw_fp, "circleit.%s(" ");\n", str);
	fprintf(mp_draw_fp, "%s.c = (%ldu,%ldu);\n", str, x_out, y_out);
	fprintf(mp_draw_fp, "drawboxed(%s);\n", str);
	return(TRUE);
#endif
}

static INT pie_text(INT *contrl, INT *ptsin, INT *intin)
{
	INT x_dev, y_dev, x_out, y_out;
	INT radius_dev, radius_out, tmp, begang, endang;
	BYTE str[256], *s;
	BYTE text[10000], c;
	INT str_len, len, i;
	
	len = (INT)contrl[3] - 2; /* strlen = Anzahl Zeichen ohne Nullbyte */
#if 1
	if (len <= 0)
		return(TRUE);
#endif
	str_len = 0;
	for (i = 0; i < len; i++) {
		c = (BYTE) intin[2 + i];
		text[str_len++] = c;
		}
	text[str_len] = 0;
	// printf("mp: pie_text, text = %s\n", text);
	// fflush(stdout);
	
	x_dev = (INT) ptsin[0];
	y_dev = (INT) ptsin[1];
	s = get_label(x_dev, y_dev);
	strcpy(str, s);
	// printf("pietext %s at %ld %ld\n", str, x_dev, y_dev);
	
	ko_min_max(x_dev, y_dev);
	mp_dev2mp(x_dev, y_dev, &x_out, &y_out);
	
	radius_dev = (INT) ptsin[6];
	mp_dev2mp_dist(radius_dev, 0, &radius_out, &tmp);
	if (radius_out <= 0)
		radius_out = 1;

	begang = (INT) intin[0];
	endang = (INT) intin[1];
	begang /= 10;
	endang /= 10;
	/* !!! filled ellipse statt pie */

	fprintf(mp_draw_fp, "circleit.%s(btex %s etex);\n", str, text);
	fprintf(mp_draw_fp, "%s.c = (%ldu,%ldu);\n", str, x_out, y_out);
	fprintf(mp_draw_fp, "unfill bpath %s;\n", str);
	fprintf(mp_draw_fp, "drawboxed(%s);\n", str);
	return(TRUE);
}

static INT st_alignment(INT *intin)
{
	txt_halign = intin[0];
	txt_valign = intin[1];
#ifdef DEBUG
	printf("vst_alignment():|txt_halign = %d txt_valign = %d\n", 
		txt_halign, txt_valign);
#endif
	return(TRUE);
}

static void get_alignment(BYTE *align)
{
	if (txt_halign == 2) { // right aligned, text to the left of the current position
		if (txt_valign == 2) 
			strcpy(align, ".llft");
		else if (txt_valign == 1) 
			strcpy(align, ".lft");
		else if (txt_valign == 0) 
			strcpy(align, ".ulft");
		}
	else if (txt_halign == 1) { // horizontally centered
		if (txt_valign == 2) 
			strcpy(align, ".bot");
		else if (txt_valign == 1) 
			strcpy(align, "");
		else if (txt_valign == 0) 
			strcpy(align, ".top");
		}
	else if (txt_halign == 0) {
		if (txt_valign == 2) 
			strcpy(align, ".lrt");
		else if (txt_valign == 1) 
			strcpy(align, ".rt");
		else if (txt_valign == 0) 
			strcpy(align, ".urt");
		}
}

static INT sl_udsty(INT *intin)
{
	line_dashing = intin[0];
	printf("mp: sl_udsty: line_dashing = %ld\n", line_dashing);
	return OK;
}

static INT sl_ends(INT *intin)
{
	line_beg_style = intin[0];
	line_end_style = intin[1];
	return OK;	
}

static INT sf_interior(INT *intin)
{
	fill_interior = intin[0];
#ifdef DEBUG
	printf("mp: sf_interior():|fill_interior = %d\n", fill_interior);
#endif
	return(TRUE);
}

static INT sf_color(INT *intin)
{
	fill_color = intin[0];
#ifdef DEBUG
	printf("mp: sf_color():|fill_color = %d\n", fill_color);
#endif
	return(TRUE);
}

static INT sf_shape(INT *intin)
{
	fill_shape = intin[0];
#ifdef DEBUG
	printf("mp: sf_shape():|shape = %d\n", fill_shape);
#endif
	return(TRUE);
}

static INT sf_outline(INT *intin)
{
	fill_outline = intin[0];
#ifdef DEBUG
	printf("mp: sf_outline():|outline = %d\n", fill_outline);
#endif
	return(TRUE);
}

static INT sf_nofill(INT *intin)
{
	fill_nofill = intin[0];
#ifdef DEBUG
	printf("mp: sf_nofill():|nofill = %d\n", fill_nofill);
#endif
	return(TRUE);
}

static INT st_boxed(INT *intin)
{
	txt_boxed = intin[0];
#ifdef DEBUG
	printf("mp: st_boxed():|boxed = %d\n", txt_boxed);
#endif
	return(TRUE);
}

static INT st_overwrite(INT *intin)
{
	txt_overwrite = intin[0];
#ifdef DEBUG
	printf("mp: st_overwrite():|overwrite = %d\n", txt_overwrite);
#endif
	return(TRUE);
}


static INTEGER_OP co_table_len = NULL;
static VECTOR_OP co_table = NULL;
static VECTOR_OP co_label = NULL;

static BYTE *get_label(INT x, INT y)
{
	VECTOR_OB v;
	STRING_OP s;
	STRING_OB ss;
	BYTE str[1024];
	INT idx, f_found;
	INT l;
	
	v.m_il(2);
	v.m_ii(0, x);
	v.m_ii(1, y);
	co_table->search(co_table_len->s_i(), TRUE, &v, &idx, &f_found);
	if (f_found) {
		s = (STRING_OP) co_label->s_i(idx - 1);
		return s->s_str();
		}
	co_table->insert_at(co_table_len, idx, &v);
	co_table_len->dec();
	l = co_table_len->s_i();
	sprintf(str, "l%ld", l);
	ss.init(str);
	// printf("mp: new label %s at idx=%ld\n", str, idx);
	// fflush(stdout);
	co_label->insert_at(co_table_len, idx, &ss);
	s = (STRING_OP) co_label->s_i(idx);
	return s->s_str();
}

static void mp_free_local_data()
{
	freeall(co_table_len);
	freeall(co_table);
	freeall(co_label);
	co_table_len = NULL;
	co_table = NULL;
	co_label = NULL;
}

static INT draw_boxes()
{
	INT i, l;
	BYTE *str;
	STRING_OP s;
	
	l = co_table_len->s_i();
	for (i = 0; i < l; i++) {
		s = (STRING_OP) co_label->s_i(i);
		str = s->s_str();
		fprintf(mp_draw_fp, "unfill bpath %s;\n", str);
		fprintf(mp_draw_fp, "drawboxed(%s);\n", str);
		}
	return OK;
}

static INT mp_xmin = 0, mp_ymin = 0, mp_xmax = 1500, mp_ymax = 1500;

INT mp_set_output_coordinates(INT xmin, INT ymin, INT xmax, INT ymax)
{
	mp_xmin = xmin;
	mp_ymin = ymin;
	mp_xmax = xmax;
	mp_ymax = ymax;
	return OK;
}

INT mp_get_output_coordinates(INT *xmin, INT *ymin, INT *xmax, INT *ymax)
{
	*xmin = mp_xmin;
	*ymin = mp_ymin;
	*xmax = mp_xmax;
	*ymax = mp_ymax;
	return OK;
}

INT draw_mp_file(void *data, VDEVICE *vdev, BYTE *name, INT factor_1000, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev))
{
	FILE *fp;
	INT xmin = mp_xmin, ymin = mp_ymin, xmax = mp_xmax, ymax = mp_ymax;
	
	if ((fp = fopen(name, "w")) == NIL) {
		return error("draw_mp_file() error in fopen()");
		}

	mp_header(fp, name);
	mp_begin_figure(fp, factor_1000);
	mp_begin_drawing(fp);
	draw_mp_picture(fp, data, vdev, 
		xmin, ymin, xmax, ymax, draw_func);
	mp_end_drawing(fp);
	mp_end_figure(fp);
	mp_footer(fp);
	
	fclose(fp);
	printf("draw_mp_file()|written file %s of size %ld\n", 
		name, file_size(name));
	fflush(stdout);
	return OK;
}

INT draw_mp_picture(FILE *fp, 
	void *data, VDEVICE *vdev, 
	INT xmin, INT ymin, INT xmax, INT ymax, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev))
/* file output of a page;
 * output will be written into open file fp. 
 * can be draw_func or a tapestry.
 * mp_draw_dev[] has to be set 
 * before calling this routine. */
{
	INT s_co[4], s_type;
	
	VdevI2(vdev, gl_vdi_contrl, gl_vdi_intin, gl_vdi_ptsin, 
		gl_vdi_intout, gl_vdi_ptsout);
	s_type = vdev->type;
	s_co[0] = vdev->co[0];
	s_co[1] = vdev->co[1];
	s_co[2] = vdev->co[2];
	s_co[3] = vdev->co[3];
	vdev->type = 4; /* MetaPost output */
	vdev->co[0] = 0;
	vdev->co[1] = 0;
	vdev->co[2] = 1000000; // pix_x_max
	vdev->co[3] = 1000000; // pix_y_max

	mp_draw_fp = fp;
	
	/* llx/lly/urx/ury */
	mp_draw_dev[0] = vdev->co[0];
	mp_draw_dev[1] = vdev->co[1];
	mp_draw_dev[2] = vdev->co[2];
	mp_draw_dev[3] = vdev->co[3];
	
	mp_draw_mp[0] = xmin; // (llx)
	mp_draw_mp[1] = ymin; // (lly)
	mp_draw_mp[2] = xmax; // (urx)
	mp_draw_mp[3] = ymax; // (ury)


	printf("mp_draw_dev = [%ld,%ld,%ld,%ld]\n", 
		mp_draw_dev[0], 
		mp_draw_dev[1], 
		mp_draw_dev[2], 
		mp_draw_dev[3]);
	printf("mp_draw_mp = [%ld,%ld,%ld,%ld]\n", 
		mp_draw_mp[0], 
		mp_draw_mp[1], 
		mp_draw_mp[2], 
		mp_draw_mp[3]);


	Vswr_mode(vdev, 1);
	Vsf_color(vdev, 0);
	Vsf_interior(vdev, 0);
	Vsf_perimeter(vdev, 0);
	(*draw_func)(data, 0, 0, 0, 0, vdev);
	Vswr_mode(vdev, 1);
	
	
	printf("mp_draw_dev = [%ld,%ld,%ld,%ld]\n", 
		mp_draw_dev[0], 
		mp_draw_dev[1], 
		mp_draw_dev[2], 
		mp_draw_dev[3]);
	printf("x_min = %ld x_max = %ld y_min = %ld y_max = %ld\n", x_min, x_max, y_min, y_max);

	vdev->type = s_type;
	vdev->co[0] = s_co[0];
	vdev->co[1] = s_co[1];
	vdev->co[2] = s_co[2];
	vdev->co[3] = s_co[3];
	return OK;
}

INT mp_header(FILE *fp, BYTE *fname)
{
	BYTE str[256];
	FILE *fp1;
	
	f_min_max_set = FALSE;
	fprintf(fp, "%% file: %s\n", fname);
	fprintf(fp, "%% created by DISCRETA MetaPost interface\n");
	system("rm a");
	system("date >a");
	fp1 = fopen("a", "r");
	fgets(str, 256, fp1);
	fclose(fp1);
	fprintf(fp, "%% creation date: %s\n\n", str);
	
	fprintf(fp, "input boxes\n\n");
#if 0
	fprintf(fp, "vardef  cuta(suffix a,b) expr p =\n");
	fprintf(fp, "  draw p cutbefore bpath.a cutafter bpath.b;\n");
	fprintf(fp, "enddef;\n\n");
#endif
	

	return OK;
}

INT mp_footer(FILE *fp)
{

	fprintf(fp, "end\n\n");
	return OK;
}

INT mp_begin_figure(FILE *fp, INT factor_1000)
{
	double d;
	BYTE str[1000];
	INT i, l;
	
	d = (double) factor_1000 * 0.001  * 0.1;
	fprintf(fp, "defaultfont:=\"cmr7\";\n");
	sprintf(str, "u=%gmm;", d);
	l = strlen(str);
	for (i = 0; i < l; i++) {
		if (str[i] == ',')
			str[i] = '.';
		}
	fprintf(fp, "%s\n", str);
	fprintf(fp, "beginfig(1);\n");

	return OK;
}

INT mp_end_figure(FILE *fp)
{
	fprintf(fp, "endfig;\n\n");
	return OK;
}

INT mp_begin_drawing(FILE *fp)
{
	if (co_table_len == NULL) {
		co_table_len = (INTEGER_OP) callocobject("co_table_len");
		co_table = (VECTOR_OP) callocobject("co_table");
		co_label = (VECTOR_OP) callocobject("co_label");
		co_table_len->m_i(0);
		co_table->m_il(VECTOR_OVERSIZE);
		co_label->m_il(VECTOR_OVERSIZE);
		}

	// INT mp_draw_dev[4]; /* llx/lly/urx/ury */
	// the drawing is in typical X-windows coordinates,
	// the origin (0,0) lies at the {\em top left} corner of the window !
	mp_draw_dev[0] = 0;
	mp_draw_dev[1] = 1000000;
	mp_draw_dev[2] = 1000000;
	mp_draw_dev[3] = 0;
	return OK;
}

INT mp_end_drawing(FILE *fp)
{
	draw_boxes();
	mp_free_local_data();
	return OK;
}


#endif /* GRAPHICS_TRUE */


