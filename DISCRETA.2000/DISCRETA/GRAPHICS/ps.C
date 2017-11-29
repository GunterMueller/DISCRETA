/* ps.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <DISCRETA/graphics.h>
#include <DISCRETA/divs.h>

#include <stdlib.h>

/*#define DEBUG */

#define POWERSHOW
#define ASPECT_RATIO

/* die draw-routinen gehen davon aus, 
 * dass ps_draw_fp / ps_draw_dev / ps_draw_ps gesetzt sind. */

/* globals: */
FILE *ps_draw_fp = NIL;
INT ps_draw_dev[4]; /* llx/lly/urx/ury */
INT ps_draw_ps[4]; /* llx/lly/urx/ury */

// static INT plot_turn;
// static INT angle; 
/* in 1/10 Altgrad (0 bis 3600) plot_turn schon beruecksichtigt */
// static INT text_height; /* in 1/360 Zoll */
static INT txt_halign = 0; 
/* 0: links, 1: zentriert, 2: rechts */
static INT txt_valign = 0;
/* 0: unten, 1:mittig, 2: oben (haengender Text) */
/* the magic formula:
 * halign + 3 * valign
 */
/* alt: */
/* 0: Basislinie, 1: Halblinie, 2: Zeichenoberkante
 * 3: Zeichenzellenunterkante 4: Zeichenunterkante
 * 5: Zeichenzellenoberkante */
static INT fill_interior;
/* 0: leer
 * 1: voll (einfarbig)
 * 2: Muster
 * 3: Schraffur
 * 4:frei definierbar */
static INT fill_style;
static INT line_type = -1; /* vsl_type() */

#define MAX_PS_FONTS 32

static INT font_nr, font_size;
static INT selected_font_nr, selected_font_size;
static INT line_width = 1000, selected_line_width;
/* 1000 = 1.0 = default */
/* initialization of selected_* with -1 in PrintHeaderPS() ! */
static BYTE f_font_used[MAX_PS_FONTS];
static INT x_min, x_max, y_min, y_max, f_min_max_set;

static INT ps_dev2ps(INT ix_dev, INT iy_dev, INT *ix_ps, INT *iy_ps);
static INT ps_ps2dev(INT ix_ps, INT iy_ps, INT *ix_dev, INT *iy_dev);
static INT ps_dev2ps_dist(INT dx_dev, INT dy_dev, INT *dx_ps, INT *dy_ps);
static INT ps_ps2dev_dist(INT dx_ps, INT dy_ps, INT *dx_dev, INT *dy_dev);
static void ko_min_max(INT x, INT y);
static INT pline(INT *contrl, INT *ptsin);
static INT gtext(INT *contrl, INT *intin, INT *ptsin);
static INT gtext2(BYTE *text, INT x_dev, INT y_dev);
static INT fillarea(INT *contrl, INT *ptsin);
static INT arc(INT *contrl, INT *ptsin, INT *intin);
static INT pie(INT *contrl, INT *ptsin, INT *intin);
static INT pie2(INT x_dev, INT y_dev, INT radius_dev, INT begang, INT endang);
static INT pie_text(INT *contrl, INT *ptsin, INT *intin);
static INT sl_width(INT *ptsin);
static INT sl_type(INT *intin);
static INT swr_mode(INT *intin);
static INT st_alignment(INT *intin);
static INT st_font(INT *intin);
static INT st_point(INT *intin);
static INT sf_interior(INT *intin);
static INT sf_style(INT *intin);
static INT draw_ps_page_file(FILE *fp, BYTE *fname, 
	void *data, VDEVICE *vdev, 
	INT pix_x_max, INT pix_y_max, 
	INT f_tape /* TRUE: LwDTapestry2, FALSE: (*draw_func) */ , 
	LW_TAPESTRY_CONTENTS *p, /* only for tape */
	INT x0, INT y0, INT w0, INT h0,  /* only for tape */
	INT i, INT j, INT lines, INT columns, INT first,  INT num, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev));
static INT PrintTapestryPS(LW_TAPESTRY_CONTENTS *lw, 
	INT f_multiple_files, 
	BYTE *name /* base_name in case multiple files */, INT first, INT len);

INT PsDraw(VDEVICE *vdev)
{
	INT *contrl = NIL;
	INT *intin = NIL;
	INT *ptsin = NIL;
	
	if (!VdevGArrays(vdev, &contrl, &intin, &ptsin, NIL, NIL)) {
		Srff("PsDraw", "VdevGArrays");
		return(FALSE);
		}
	if (contrl[0] == 6) { /* v_pline */
		pline(contrl, ptsin);
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
		sl_type(intin);
		}
	else if (contrl[0] == 16) { /* vsl_width */
		sl_width(ptsin);
		}
	else if (contrl[0] == 21) { /* vst_font */
		st_font(intin);
		}
	else if (contrl[0] == 23) { /* vsf_interior */
		sf_interior(intin);
		}
	else if (contrl[0] == 24) { /* vsf_style */
		sf_interior(intin);
		}
	else if (contrl[0] == 32) { /* vswr_mode */
		swr_mode(intin);
		}
	else if (contrl[0] == 39) { /* vst_alignment */
		st_alignment(intin);
		}
	else if (contrl[0] == 107) { /* vst_point */
		st_point(intin);
		}
	return(TRUE);
}

static INT ps_dev2ps(INT ix_dev, INT iy_dev, INT *ix_ps, INT *iy_ps)
{
	INT dx, dy;
	double a, b;

	dx = ix_dev - ps_draw_dev[0];
	dy = iy_dev - ps_draw_dev[1];
	a = (double)dx / (double)(ps_draw_dev[2] - ps_draw_dev[0]);
	b = (double)dy / (double)(ps_draw_dev[3] - ps_draw_dev[1]);
	dx = (INT)(a * (double)(ps_draw_ps[2] - ps_draw_ps[0]));
	dy = (INT)(b * (double)(ps_draw_ps[3] - ps_draw_ps[1]));
	*ix_ps = dx + ps_draw_ps[0];
	*iy_ps = dy + ps_draw_ps[1];
	return(TRUE);
}

static INT ps_ps2dev(INT ix_ps, INT iy_ps, INT *ix_dev, INT *iy_dev)
{
	INT dx, dy;
	double a, b;

	dx = ix_ps - ps_draw_ps[0];
	dy = iy_ps - ps_draw_ps[1];
	a = (double)dx / (double)(ps_draw_ps[2] - ps_draw_ps[0]);
	b = (double)dy / (double)(ps_draw_ps[3] - ps_draw_ps[1]);
	dx = (INT)(a * (double)(ps_draw_dev[2] - ps_draw_dev[0]));
	dy = (INT)(b * (double)(ps_draw_dev[3] - ps_draw_dev[1]));
	*ix_dev = dx + ps_draw_dev[0];
	*iy_dev = dy + ps_draw_dev[1];
	return(TRUE);
}

static INT ps_dev2ps_dist(INT dx_dev, INT dy_dev, INT *dx_ps, INT *dy_ps)
{
	INT dx, dy;
	double a, b;

	a = (double)dx_dev / (double)(ps_draw_dev[2] - ps_draw_dev[0]);
	b = (double)dy_dev / (double)(ps_draw_dev[3] - ps_draw_dev[1]);
	dx = (INT)(a * (double)(ps_draw_ps[2] - ps_draw_ps[0]));
	dy = (INT)(b * (double)(ps_draw_ps[3] - ps_draw_ps[1]));
	*dx_ps = dx;
	*dy_ps = dy;
	return(TRUE);
}

static INT ps_ps2dev_dist(INT dx_ps, INT dy_ps, INT *dx_dev, INT *dy_dev)
{
	INT dx, dy;
	double a, b;

	a = (double)dx_ps / (double)(ps_draw_ps[2] - ps_draw_ps[0]);
	b = (double)dy_ps / (double)(ps_draw_ps[3] - ps_draw_ps[1]);
	dx = (INT)(a * (double)(ps_draw_dev[2] - ps_draw_dev[0]));
	dy = (INT)(b * (double)(ps_draw_dev[3] - ps_draw_dev[1]));
	*dx_dev = dx;
	*dy_dev = dy;
	return(TRUE);
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
	INT x_dev, y_dev, x_ps, y_ps;
	INT i, i2;
	BYTE str[256];

	/* printf("selected_line_width = %ld\n", selected_line_width);
	printf("line_width = %ld\n", line_width); */
	fflush(stdout);
	if (selected_line_width != line_width) {
		sprintf(str, "%ld 1000 div setlinewidth\n", line_width);
		fputs(str, ps_draw_fp);
		/* printf("%s", str);
		fflush(stdout); */
		selected_line_width = line_width;
		}
	for (i = 0; i < contrl[1]; i++) {
		i2 = i << 1;
		x_dev = (INT)ptsin[i2];
		y_dev = (INT)ptsin[i2 + 1];
		ko_min_max(x_dev, y_dev);
		ps_dev2ps(x_dev, y_dev, &x_ps, &y_ps);
		if (i == 0) {
			fputs("newpath\n", ps_draw_fp);
			sprintf(str, "%ld %ld moveto\n", x_ps, y_ps);
			fputs(str, ps_draw_fp);
#ifdef DEBUG
			printf(str);
#endif
			}
		else {
			sprintf(str, "%ld %ld lineto\n", x_ps, y_ps);
			fputs(str, ps_draw_fp);
#ifdef DEBUG
			printf(str);
#endif
			}
		}
	fputs("stroke\n", ps_draw_fp);
	return(TRUE);
}

static INT gtext(INT *contrl, INT *intin, INT *ptsin)
{
	INT i, len;
	INT x_dev, y_dev;
	BYTE text[10000];
	
	x_dev = (INT)ptsin[0];
	y_dev = (INT)ptsin[1];
	len = (INT)contrl[3]; /* Anzahl Zeichen ohne Nullbyte */
	for (i = 0; i < len; i++) {
		text[i] = (BYTE) intin[i];
		}
	text[i] = 0;
	gtext2(text, x_dev, y_dev);

	return(TRUE);
}

static INT gtext2(BYTE *text, INT x_dev, INT y_dev)
{
	INT str_len, i, len, align;
	INT x_ps, y_ps;
	BYTE str[10000], c;
	
	if ((selected_font_nr != font_nr) || (selected_font_size != font_size)) {
		strcpy(str, "/");
		if (!ps_font_name(font_nr, str))
			return FALSE;
		sprintf(str + strlen(str), " findfont %ld scalefont setfont\n", font_size);
		/* sprintf(str, "/Times-Roman findfont %ld scalefont setfont\n", text_height); */
		fputs(str, ps_draw_fp);
		fputs("/y_height newpath 0 0 moveto\n", ps_draw_fp);
		fputs("  (XgjQ) true charpath flattenpath pathbbox\n", ps_draw_fp);
		fputs("  3 2 roll sub 3 1 roll pop pop def\n", ps_draw_fp);
		selected_font_nr = font_nr;
		selected_font_size = font_size;
		f_font_used[font_nr] = TRUE;
		}
	ko_min_max(x_dev, y_dev);
	ps_dev2ps(x_dev, y_dev, &x_ps, &y_ps);
	fputs("newpath\n", ps_draw_fp);
	fprintf(ps_draw_fp, "%ld %ld moveto\n", x_ps, y_ps);
	// printf("ps: gtext2 text = %s\n", text);
	
	len = strlen(text);
	if (len <= 0)
		return(TRUE);
	str_len = 0;
	str[str_len++] = '(';
	for (i = 0; i < len; i++) {
		c = text[i];
		if (c == '(') {
			str[str_len++] = '\\';
			str[str_len++] = c;
			}
		else if (c == ')') {
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
	str[str_len] = 0;
#ifdef POWERSHOW
	align = txt_halign + 3 * txt_valign;
	sprintf(str + strlen(str), ") %ld powershow\n", align);
#else
	strcat(str, ") show\n");
#endif
	fputs(str, ps_draw_fp);

	return(TRUE);
}

static INT fillarea(INT *contrl, INT *ptsin)
{
	INT x_dev, y_dev, x_ps, y_ps;
	INT i, i2;
	BYTE str[256];

	for (i = 0; i < contrl[1]; i++) {
		i2 = i << 1;
		x_dev = (INT)ptsin[i2];
		y_dev = (INT)ptsin[i2 + 1];
		ko_min_max(x_dev, y_dev);
		ps_dev2ps(x_dev, y_dev, &x_ps, &y_ps);
		
		if (i == 0) {
			fputs("newpath\n", ps_draw_fp);
			sprintf(str, "%ld %ld moveto\n", x_ps, y_ps);
			fputs(str, ps_draw_fp);
#ifdef DEBUG
			printf(str);
#endif
			}
		else {
			sprintf(str, "%ld %ld lineto\n", x_ps, y_ps);
			fputs(str, ps_draw_fp);
#ifdef DEBUG
			printf(str);
#endif
			}
		}
	if (fill_interior == 0) { /* fill it with white: */
		fputs("gsave 1 setgray fill grestore stroke\n", ps_draw_fp);
		}
	else if (fill_interior == 1) { /* fill it with black: */
		fputs("gsave 0 setgray fill stroke grestore\n", ps_draw_fp);
		}
	return(TRUE);
}

static INT arc(INT *contrl, INT *ptsin, INT *intin)
{
	INT x_dev, y_dev, x_ps, y_ps;
	INT radius_dev, radius_ps, tmp, begang, endang;
	BYTE str[256];
	
	x_dev = (INT)ptsin[0];
	y_dev = (INT)ptsin[1];
	ko_min_max(x_dev, y_dev);
	ps_dev2ps(x_dev, y_dev, &x_ps, &y_ps);
	
	radius_dev = (INT)ptsin[6];
	ps_dev2ps_dist(radius_dev, 0L, &radius_ps, &tmp);

	begang = (INT)intin[0];
	endang = (INT)intin[1];
	begang /= 10;
	endang /= 10;
	sprintf(str, "%ld %ld %ld %ld %ld arc stroke\n", 
		x_ps, y_ps, radius_ps, begang, endang);
	fputs("newpath\n", ps_draw_fp);
	fputs(str, ps_draw_fp);

#ifdef DEBUG
	printf(str);
#endif
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
	INT x_ps, y_ps;
	INT radius_ps, tmp;
	BYTE str[256];
	
	ko_min_max(x_dev, y_dev);
	ps_dev2ps(x_dev, y_dev, &x_ps, &y_ps);
	
	ps_dev2ps_dist(radius_dev, 0, &radius_ps, &tmp);
	/* printf("radius_dev = %ld radius_ps = %ld\n", radius_dev, radius_ps);
	fflush(stdout); */

	begang /= 10;
	endang /= 10;
	/* !!! filled ellipse statt pie */
	/* following change was made by wolfb
	/* (date: 11 march 1996): */
	sprintf(str, "%ld %ld %ld %ld %ld arc fill stroke\n", 
		x_ps, y_ps, radius_ps, begang, endang);
	fputs("newpath\n", ps_draw_fp);
	fputs(str, ps_draw_fp);

#ifdef DEBUG
	printf(str);
#endif
	return(TRUE);
}

static INT pie_text(INT *contrl, INT *ptsin, INT *intin)
{
	BYTE text[10000];
	INT len, i;
	INT x_dev, y_dev, radius_dev, begang, endang;
	
	len = (INT)contrl[3] - 2; /* strlen = Anzahl Zeichen ohne Nullbyte */
	for (i = 0; i < len; i++) {
		text[i] = (BYTE) intin[2 + i];
		}
	text[i] = 0;
	// printf("ps: pie_text, text = %s\n", text);
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

static INT sl_width(INT *ptsin)
{
	line_width = ptsin[0];
#ifdef DEBUG
	printf("vsl_width():line_width = %d\n", line_width);
	fflush(stdout);
#endif
	return(TRUE);
}

static INT sl_type(INT *intin)
{
	line_type = intin[0];
#if FALSE
#ifdef DEBUG
	printf("vsl_type():|line_type = %d\n", line_type);
#endif
	switch (line_type) {
		case 1:
			PenPat(black);
			break;
		case 2:
			PenPat(dkGray);
			break;
/*		case 3:
			PenPat(gray);
			break;*/
		case 3:
		case 4:
			PenPat(ltGray);
			break;
		case 5:
			PenPat(white);
			break;
		}
#endif
	return(TRUE);
}

static INT swr_mode(INT *intin)
{
	INT mode;
	
	mode = intin[0];
#if FALSE
	switch (mode) {
		case 1:
			PenMode(patCopy);
#ifdef DEBUG
			printf("PenMode(patCopy)\n");
#endif
			break;
		case 2:
			PenMode(patXor);
#ifdef DEBUG
			printf("PenMode(patXor)\n");
#endif
			break;
		}
#endif
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

static INT st_font(INT *intin)
{
	font_nr = intin[0];
#if 0
	INT i;
	BYTE str[256];
	/* Texthoehe in Druckerpunkten von 1/360 bzw. 1/72 Zoll */
	text_height = intin[0];
	i = (INT)((double)text_height / 3 + 0.5);
	sprintf(str, "/Times-Roman findfont %ld scalefont setfont\n", text_height);
	fputs(str, ps_draw_fp);
#ifdef DEBUG
	printf("vst_point():|text_height = %d i = %d\n", text_height, i);
#endif
#endif
	return(TRUE);
}

static INT st_point(INT *intin)
{
	font_size = intin[0];
#if 0
	INT i;
	BYTE str[256];
	
	/* Texthoehe in Druckerpunkten von 1/360 bzw. 1/72 Zoll */
	text_height = intin[0];
	i = (INT)((double)text_height / 3 + 0.5);
	sprintf(str, "/Times-Roman findfont %ld scalefont setfont\n", text_height);
	fputs(str, ps_draw_fp);
#ifdef DEBUG
	printf("vst_point():|text_height = %d i = %d\n", text_height, i);
#endif
#endif
	return(TRUE);
}

static INT sf_interior(INT *intin)
{
	fill_interior = intin[0];
#ifdef DEBUG
	printf("sf_interior():|fill_interior = %d\n", fill_interior);
#endif
	return(TRUE);
}

static INT sf_style(INT *intin)
{
	fill_style = intin[0];
#ifdef DEBUG
	printf("sf_style():|fill_style = %d\n", fill_style);
#endif
	return(TRUE);
}

void PrintHeaderPS(FILE *fp, BYTE *fname)
{
	BYTE str[256];
	FILE *fp1;
	INT i;
	
	fputs("%!PS-Adobe-1.0\n", fp);
	fputs("%%Creator: DISCRETA PostScript output interface\n", fp);
	fprintf(fp, "%%Title: %s\n", fname);
	
	system("rm a");
	system("date >a");
	fp1 = fopen("a", "r");
	fgets(str, 256, fp1);
	fclose(fp1);
	
	fprintf(fp, "%%CreationDate: %s\n", str);
	fputs("%%Pages: 1\n", fp);
#if 0
	fputs("%%DocumentFonts: (atend)\n", fp);
	fputs("%%BoundingBox: (atend)\n", fp);
#endif
	fputs("%%DocumentFonts: ", fp);
	for (i = 0; i < MAX_PS_FONTS; i++) {
		if (f_font_used[i]) {
			str[0] = 0;
			ps_font_name(i, str);
			strcat(str, " ");
			fputs(str, fp);
			}	
		}
	fputs("\n", fp);
	/* sprintf(str, "%%%%BoundingBox: %ld %ld %ld %ld\n", 
		x_min, y_min, x_max, y_max); */
	sprintf(str, "%%%%BoundingBox: %ld %ld %ld %ld\n", 
		x_min /* before: 0 */, 
		ps_draw_ps[1] /* y_min (before: 0) */, 
		x_max /* before: 612 */, 
		ps_draw_ps[1] + y_max - y_min /* y_max */ );
	fputs(str, fp);
	/* fputs("%%BoundingBox: 0 0 612 826\n", fp); */
	/* fputs("%%BoundingBox: 0 0 612 826\n", fp); */
	fputs("%%EndComments\n", fp);
	fputs("%%EndProlog\n", fp);
	fputs("%%Page: 0 1\n", fp);
#ifdef POWERSHOW
	fputs("/powershow\n", fp);
	fputs("{ dup 3 mod -2 div\n", fp);
	fputs("  exch 3 idiv -2 div\n", fp);
	fputs("  3 2 roll y_height\n", fp);
	fputs("  exch dup stringwidth pop\n", fp);
	fputs("  3 1 roll 5 1 roll\n", fp);
	fputs("  3 2 roll mul\n", fp);
	fputs("  3 1 roll exch mul exch\n", fp);
	fputs("  rmoveto\n", fp);
	fputs("  show } def\n", fp);
#endif
	for (i = 0; i < MAX_PS_FONTS; i++)
		f_font_used[i] = FALSE;
	font_nr = 0;
	font_size = 12;
	selected_font_nr = -1;
	selected_font_size = -1;
	selected_line_width = -1;
	f_min_max_set = FALSE;
}

void PrintFooterPS(FILE *fp)
{
	fputs("%%Trailer\n", fp);
	fputs("%%EOF\n", fp);
}

INT ps_font_name(INT font_nr, BYTE str[256])
{
	if (font_nr == 0)
		strcat(str, "Times-Roman");
	else if (font_nr == 1)
		strcat(str, "Times-Italic");
	else if (font_nr == 2)
		strcat(str, "Times-Bold");
	else if (font_nr == 3)
		strcat(str, "Times-BoldItalic");
	else if (font_nr == 4)
		strcat(str, "Helvetica");
	else if (font_nr == 5)
		strcat(str, "Helvetica-Oblique");
	else if (font_nr == 6)
		strcat(str, "Helvetica-Bold");
	else if (font_nr == 7)
		strcat(str, "Helvetica-BoldOblique");
	else if (font_nr == 8)
		strcat(str, "Courier");
	else if (font_nr == 9)
		strcat(str, "Courier-Oblique");
	else if (font_nr == 10)
		strcat(str, "Courier-Bold");
	else if (font_nr == 11)
		strcat(str, "Courier-BoldOblique");
	else if (font_nr == 12)
		strcat(str, "Symbol");
	else
		return FALSE;
	return TRUE;
}

INT draw_PS_file(void *data, VDEVICE *vdev, 
	INT pix_x_max, INT pix_y_max, BYTE *name, 
	INT (*draw_func)(void *gp_lattice, INT x, INT y, INT w, INT h, VDEVICE *vdev))
{
#if 0
	INT s_co[4], s_type;
	INT ix_dev, iy_dev;
#endif
	FILE *fp;
	
	if ((fp = fopen(name, "w")) == NIL) {
		Srff("draw_PS_file", "fopen");
		return ERROR;
		}

	// the drawing is in typical X-windows coordinates,
	// the origin (0,0) lies at the {\em top left} corner of the window !
	ps_draw_dev[0] = 0;
	ps_draw_dev[1] = pix_y_max;
	ps_draw_dev[2] = pix_x_max;
	ps_draw_dev[3] = 0;
	draw_ps_page_file(fp, name, 
		data, vdev, 
		pix_x_max, pix_y_max, 
		FALSE /* f_tape */, 
		NIL, /* only for tape */
		0, 0, 0, 0, 
		0, 0, 0, 0, 
		0, 0,  
		draw_func);

	if (fflush(fp) != NIL) {
		Srff("draw_PS_file", "fflush");
		return ERROR;
		}
	if (fclose(fp) != NIL) {
		Srff("draw_PS_file", "fclose");
		return ERROR;
		}
	printf("draw_PS_file()|written file %s of size %ld\n", 
		name, file_size(name));

#if 0
	vdev->type = s_type;
	vdev->co[0] = s_co[0];
	vdev->co[1] = s_co[1];
	vdev->co[2] = s_co[2];
	vdev->co[3] = s_co[3];
#endif
	return OK;
}

static INT draw_ps_page_file(FILE *fp, BYTE *fname, 
	void *data, VDEVICE *vdev, 
	INT pix_x_max, INT pix_y_max, 
	INT f_tape /* TRUE: LwDTapestry2, FALSE: (*draw_func) */ , 
	LW_TAPESTRY_CONTENTS *p, /* only for tape */
	INT x0, INT y0, INT w0, INT h0,  /* only for tape */
	INT i, INT j, INT lines, INT columns, INT first, INT num, 
	INT (*draw_func)(void *data, INT x, INT y, INT w, INT h, VDEVICE *vdev))
/* file output of a page;
 * output will be written into open file fp. 
 * can be draw_func or a tapestry.
 * ps_draw_dev[] has to be set 
 * before calling this routine. */
{
	INT s_co[4], s_type;
	BYTE *tmp_name = "/tmp/tmp.tmp";
	FILE *tmp_fp;
	
	/* printf("draw_ps_page_file: first = %ld num = %ld\n", first, num); */
	VdevI2(vdev, gl_vdi_contrl, gl_vdi_intin, gl_vdi_ptsin, 
		gl_vdi_intout, gl_vdi_ptsout);
	s_type = vdev->type;
	s_co[0] = vdev->co[0];
	s_co[1] = vdev->co[1];
	s_co[2] = vdev->co[2];
	s_co[3] = vdev->co[3];
	vdev->type = 0; /* PS output */

	ps_draw_dev[0] = 0;
	ps_draw_dev[1] = pix_y_max;
	ps_draw_dev[2] = pix_x_max;
	ps_draw_dev[3] = 0;
	vdev->co[0] = 0;
	vdev->co[1] = 0;
	vdev->co[2] = pix_x_max;
	vdev->co[3] = pix_y_max;
#if 0
	vdev->co[0] = 0;
	vdev->co[1] = pix_x_max;
	vdev->co[2] = pix_y_max;
	vdev->co[3] = 0;
#endif

	/* 
	 * the first drawing goes into a temporary file; 
	 * we only want to know the extrema 
	 */
	if ((tmp_fp = fopen(tmp_name, "w")) == NIL) {
		printf("draw_ps_page_file()| error in fopen('%s')\n", tmp_name);
		Srff("draw_ps_page_file", "fopen");
		return ERROR;
		}
	PrintHeaderPS(tmp_fp, tmp_name);
	
	ps_draw_fp = tmp_fp;
	ps_draw_ps[0] = 36;
	ps_draw_ps[1] = 36;
	ps_draw_ps[2] = (INT)((21. / 2.54) * 72. - 36. );
	ps_draw_ps[3] = (INT)((29.7 / 2.54) * 72. - 72.);
	/* fputs("/epson finddevice setdevice\n\n", fp); */
	Vswr_mode(vdev, 1);
	Vsf_color(vdev, 0);
	Vsf_interior(vdev, 0);
	Vsf_perimeter(vdev, 0);
	if (f_tape) {
		if (!LwDTapestry2(p, vdev, x0, y0, w0, h0, i, j, 
			lines, columns, first, num)) {
			Srff("draw_ps_page_file", "LwDTapestry2");
			return FALSE;
			}
		}
	else {
		(*draw_func)(data, 0, 0, pix_x_max, pix_y_max, vdev);
		}
	fputs("\nshowpage\n\n", tmp_fp);
	Vswr_mode(vdev, 1);
	PrintFooterPS(tmp_fp);
	if (fflush(tmp_fp) != NIL) {
		Srff("draw_ps_page_file", "fflush");
		return ERROR;
		}
	if (fclose(tmp_fp) != NIL) {
		Srff("draw_ps_page_file", "fclose");
		return ERROR;
		}


#ifdef ASPECT_RATIO
	{
	INT dx, dy, Dx, Dy;
	double ar;
	
	dx = x_max - x_min;
	dy = y_max - y_min;
	ar = (double) dy / (double) dx;
	Dx = (INT)((21. / 2.54) * 72. - 36. -36.);
	Dy = (INT)(((double) Dx) * ar);
	printf("ps: ASPECT_RATIO dx = %ld dy = %ld, ar = %lf, Dx = %ld, Dy = %ld\n", 
		dx, dy, ar, Dx, Dy);
	ps_draw_ps[0] = 36;
	ps_draw_ps[1] = 36;
	ps_draw_ps[2] = 36 + Dx;
	ps_draw_ps[3] = 36 + Dy;
	}
#endif
#if 1
	printf("draw_ps_page_file: x_min = %ld x_max = %ld\n", x_min, x_max);
	printf("draw_ps_page_file: y_min = %ld y_max = %ld\n", y_min, y_max);
	// printf("draw_ps_page_file: iy_dev = %ld\n", iy_dev);
#endif
	
	ps_draw_dev[0] = x_min;
	ps_draw_dev[1] = y_max;
	ps_draw_dev[2] = x_max;
	ps_draw_dev[3] = y_min;
	
	
	/* 
	 * now we draw into the given file fp 
	 */
	PrintHeaderPS(fp, fname);
	ps_draw_fp = fp;
	Vswr_mode(vdev, 1);
	Vsf_color(vdev, 0);
	Vsf_interior(vdev, 0);
	Vsf_perimeter(vdev, 0);
	if (f_tape) {
		if (!LwDTapestry2(p, vdev, x0, y0, w0, h0, i, j, 
			lines, columns, first, num)) {
			Srff("draw_ps_page_file", "LwDTapestry2");
			return FALSE;
			}
		}
	else {
		(*draw_func)(data, 0, 0, 30000, 30000, vdev);
		}
	fputs("\nshowpage\n\n", fp);
	Vswr_mode(vdev, 1);
	PrintFooterPS(fp);

	vdev->type = s_type;
	vdev->co[0] = s_co[0];
	vdev->co[1] = s_co[1];
	vdev->co[2] = s_co[2];
	vdev->co[3] = s_co[3];
	return OK;
}

INT draw_tape_PS(LW_TAPESTRY_CONTENTS *p, BYTE *name)
{
#if 0
	FILE *fp;

	if ((fp = fopen(name, "w")) == NIL) {
		Srff("draw_tape_PS", "fopen");
		return(FALSE);
		}
	PrintHeaderPS(fp);
	PrintTapestryPS(p, fp);
	PrintFooterPS(fp);
	if (fflush(fp) != NIL) {
		Srff("draw_tape_PS", "fflush");
		return(FALSE);
		}
	if (fclose(fp) != NIL) {
		Srff("draw_tape_PS", "fclose");
		return(FALSE);
		}
	printf("draw_tape_PS()|written file %s of size %ld\n", 
		name, file_size(name));
#endif
	PrintTapestryPS(p, FALSE /* f_multiple_files */, name, -1, -1);
	return TRUE;
}

INT draw_tape_PS_multiple_files(LW_TAPESTRY_CONTENTS *p, BYTE *name)
{
#if 0
	PrintTapestryPS(p, TRUE, "upto63_.ps", 0, 316);
	PrintTapestryPS(p, TRUE, "65_95_.ps", 584, 207);
	PrintTapestryPS(p, TRUE, "97_127_.ps", 1022, 241);
	/* PrintTapestryPS(p, TRUE, "64_.ps", 317, 267);
	PrintTapestryPS(p, TRUE, "96_.ps", 791, 231); */
#else
	VECTOR_OB args;
	STRING_OP s;
	INTEGER_OP first, len;

	parse_args(&args);
	if (args.s_li()) {
		if (args.s_li() != 3)
			return error("draw_tape_PS_multiple_files()" 
				"3 args expected: file name, first, len");
		s = (STRING_OP) args.s_i(0);
		first = (INTEGER_OP) args.s_i(1);
		len = (INTEGER_OP) args.s_i(2);
		PrintTapestryPS(p, TRUE, s->s_str(), first->s_i(), len->s_i());
		}
	else
		PrintTapestryPS(p, TRUE, name, -1, -1);
#endif
	return TRUE;
}

INT gl_ps_lines_per_page_tape = 5;
INT gl_ps_cols_per_page_tape = 2;

static INT PrintTapestryPS(LW_TAPESTRY_CONTENTS *lw, 
	INT f_multiple_files, 
	BYTE *name /* base_name in case multiple files */, INT first, INT len)
{
	VDEVICE vdev_, *vdev = &vdev_;
	BYTE str[1024];
	INT pics_per_page, page, pages, old_type;
	INT pageHeight, pageWidth;
	INT x0, y0, w0, h0;
	FILE *fp = NIL;
	BYTE base_name[256], ext[256];
	INT i, at, num;
	
	if (f_multiple_files) {
		ext[0] = 0;
		strcpy(base_name, name);
		for (i = strlen(base_name) - 1; i >= 0; i--) {
			if (base_name[i] == '.') {
				strcpy(ext, base_name + i + 1);
				base_name[i] = 0;
				break;
				}
			}
		printf("base_name = '%s'\n", base_name);
		printf("ext = '%s'\n", ext);
		}
	if (!f_multiple_files) {
		if ((fp = fopen(name, "w")) == NIL) {
			Srff("PrintTapestryPS", "fopen");
			return(FALSE);
			}
		}
	VdevI2(vdev, gl_vdi_contrl, gl_vdi_intin, gl_vdi_ptsin, 
		gl_vdi_intout, gl_vdi_ptsout);
	old_type = vdev->type;
	vdev->type = 0; /* PS output */
	vdev->co[0] = 0;
	vdev->co[1] = 0;
	vdev->co[2] = lw->picture_width * gl_ps_cols_per_page_tape;
	vdev->co[3] = lw->picture_height * gl_ps_lines_per_page_tape;

	/* Ausmasse des bedruckten Bereichs: */
	pageWidth = (INT)((21. / 2.54) * 72. - 72. - 36.);
	pageHeight = (INT)((29.7 / 2.54) * 72. - 72. - 36.);
	ps_draw_fp = fp;
	ps_draw_dev[0] = 0;
	ps_draw_dev[1] = 0;
	ps_draw_dev[2] = lw->picture_width  * 10 * gl_ps_cols_per_page_tape;
	ps_draw_dev[3] = lw->picture_height  * 10 * gl_ps_lines_per_page_tape;
	ps_draw_ps[0] = 72;
	ps_draw_ps[1] = 36;
	ps_draw_ps[2] = ps_draw_ps[0] + pageWidth;
	ps_draw_ps[3] = ps_draw_ps[1] + pageHeight;
	/* fputs("/epson finddevice setdevice\n\n", fp); */
	/* fputs("/Times-Roman findfont 15 scalefont setfont\n", fp); */
	pics_per_page = gl_ps_lines_per_page_tape * 
			gl_ps_cols_per_page_tape;
	if (first == -1) {
		first = 0;
		len = lw->len;
		}
	printf("first = %ld len = %ld\n", first, len);
	fflush(stdout);
	pages = len / pics_per_page;
	if (pages * pics_per_page < len) 
		pages++;
	for (page = 0; page < pages; page++) {
		printf("page %ld / %ld ", page + 1, pages);
		if (f_multiple_files) {
			fflush(stdout);
			sprintf(str, "%s%03ld.%s", base_name, (INT) (page + 1), ext);
			if ((fp = fopen(str, "w")) == NIL) {
				Srff("PrintTapestryPS", "fopen");
				printf("str = '%s'\n", str);
				return(FALSE);
				}
			ps_draw_fp = fp;
			}
		else {
			printf("\n");
			fflush(stdout);
			}
		at = page * gl_ps_lines_per_page_tape;
		num = MINIMUM(gl_ps_lines_per_page_tape * 
			gl_ps_cols_per_page_tape, 
			len - at * gl_ps_cols_per_page_tape);
		x0 = 0;
		y0 = 0;
		w0 = lw->picture_width * 10; /* for more exact placement ! */
		h0 = lw->picture_height * 10;
		draw_ps_page_file(fp, str, 
			NIL /* data */, vdev, 
			lw->picture_width * gl_ps_cols_per_page_tape, 
			lw->picture_height * gl_ps_lines_per_page_tape, 
			TRUE /* f_tape */, 
			lw, 
			x0, y0, w0, h0, at, 0, 
			gl_ps_lines_per_page_tape, 
			gl_ps_cols_per_page_tape, first, num, 
			NIL /* draw_func */);
		fputs("\nshowpage\n\n", fp);
		if (f_multiple_files) {
			if (fflush(fp) != NIL) {
				Srff("PrintTapestryPS", "fflush");
				return(FALSE);
				}
			if (fclose(fp) != NIL) {
				Srff("PrintTapestryPS", "fclose");
				return(FALSE);
				}
			printf("written file '%s' of size %ld\n", 
				str, file_size(str));
			{
			FILE *f;
			f = fopen("tape_include.tex", "a");
			fprintf(f, "\\epsfxsize=15.35cm\n");
			fprintf(f, "\\begin{center}\n");
			fprintf(f, "\\mbox{ \\epsffile{PIC/PS/GROUPS/%s} }\n", str);
			fprintf(f, "\\end{center}\n\n");
			fclose(f);
			}
			}
		}
	vdev->type = old_type;
	if (!f_multiple_files) {
		if (fflush(fp) != NIL) {
			Srff("PrintTapestryPS", "fflush");
			return(FALSE);
			}
		if (fclose(fp) != NIL) {
			Srff("PrintTapestryPS", "fclose");
			return(FALSE);
			}
		printf("PrintTapestryPS()|written file %s of size %ld\n", 
			name, file_size(name));
		}
	return TRUE;
}

#endif /* GRAPHICS_TRUE */


