/* vdi.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <DISCRETA/graphics.h>
#include <DISCRETA/divs.h>
#include <stdlib.h>		// srand, rand
#include <math.h>		// sin, cos etc.

#ifndef PI
#define PI 3.1415927
#endif

INT gl_vdi_contrl[12];
INT gl_vdi_intin[128];
INT gl_vdi_intout[128];
INT gl_vdi_ptsin[128];
INT gl_vdi_ptsout[128];

void *gl_vdi_draw_tbl[10];
/* INT (*gl_vdi_draw_tbl[10])(VDEVICE *p); */
/* 0 PsDraw
 * 1 x11_draw
 * 2 MqdDraw
 * 3 EpicDraw
 * 4 mp_draw
 */

#if TEXDOCU
void Vdi_draw(VDEVICE *p)
#endif
{
	INT (* draw)(VDEVICE *);
	
	if (p->type == 0) {
		PsDraw(p);
		return;
		}
	if (p->type == 3) {
		EpicDraw(p);
		return;
		}
	if (p->type == 4) {
		mp_draw(p);
		return;
		}
	if (p->type < 10 && gl_vdi_draw_tbl[p->type]) {
		draw = (INT (* )(VDEVICE *)) gl_vdi_draw_tbl[p->type];
		(*draw)(p);
		return;
		}
	printf("Vdi_draw() unknown type\n");
}

#if TEXDOCU
INT LwDTapestry2(LW_TAPESTRY_CONTENTS *p, VDEVICE *vdev, 
	INT x0, INT y0, INT w0, INT h0, 
	INT i, INT j, INT lines, INT columns, INT first, INT num)
#endif
/* before: LISTWINDOW *p */ 
{
	INT m, n, len, l, x0_1, y0_1, w0_1, h0_1;
	INT size;
	void **draw_data = NIL;
	INT ret = FALSE;
	INT f_normal;
	
	if (vdev->type == 3 || vdev->type == 0) /* 3: epic 0: .ps*/
		f_normal = TRUE;
	else
		f_normal = FALSE; /* screen */
	len = lines * columns;
	size = len * sizeof(BYTE *);
	if (p->draw_picture == NIL) {
		printf("LwDTapestry2()|draw_picture == NIL\n");
		goto l_exit;
		}
	draw_data = (void **) my_malloc(size, "vdi.C: LwDTapestry2");
	if (draw_data == NIL) {
		printf("lines = %ld columns = %ld", lines, columns);
		return error("LwDTapestry2() no memory");
		}
	for (l = 0; l < len; l++) {
		draw_data[l] = NIL;
		}
	if (p->draw_data_get) {
		if (!(*p->draw_data_get)(p, i, j, lines, columns, first, 
			draw_data)) {
			return error("LwDTapestry2() error in draw_data_get");
			}
		}
	Vswr_mode(vdev, 1);
	Vsf_color(vdev, 0);
	Vsf_interior(vdev, 0);
	Vsf_perimeter(vdev, 0);
	/* wnd_f_clip_interior = TRUE;
	wnd_interior = p->interior; */
#ifdef SYSTEMMAC
	/* wnd_redraw_region = p->wp->visRgn; */ /* Nur Verweis; keine Kopie ! */
#endif
	l = 0;
	for (m = 0; m < lines; m++) {
		for (n = 0; n < columns; n++) {
			if (f_normal) { /* epic */
				x0_1 = x0 + w0 * n;
				y0_1 = y0 + lines * h0 - h0 * m;
				w0_1 = w0;
				h0_1 = -h0;
				}
			else { /* screen */
				x0_1 = x0 + w0 * n;
				y0_1 = y0 + h0 * m;
				w0_1 = w0;
				h0_1 = h0;
				}
			if (l < num) {
				(*p->draw_picture)(p, draw_data[l], 
					i + m, j + n, 
					x0_1, y0_1, w0_1, h0_1, vdev);
				}
			l++;
			}
		}
	/* wnd_f_clip_interior = FALSE; */
	Vswr_mode(vdev, 1);
	ret = TRUE;
l_exit:
	if (draw_data) {
		if (p->draw_data_free) {
			if (!(*p->draw_data_free)(p, 
				i, j, lines, columns, draw_data)) {
				return error("LwDTapestry() error in draw_data_free");
				}
			}
		my_free(draw_data);
		draw_data = NIL;
		}
	return (ret);
}

/*
 * VDEVICE - ATTRIB.C
 */

#if TEXDOCU
INT Vswr_mode(VDEVICE *p, INT mode)
#else
V\_set\_writing\_mode: unused
#endif
{
	p->intin[0] = mode;
	p->contrl[0] = 32;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
void Vs_color(VDEVICE *p, INT index, INT *rgb_in)
#else
unused
#endif
{
	INT i;
	
	p->intin[0] = index;
	for(i = 1; i<4; i++) 
		p->intin[i] = *rgb_in++;
	p->contrl[0] = 14;
	p->contrl[1] = 0;
	p->contrl[3] = 4;
	Vdi_draw(p);
}

#if TEXDOCU
INT Vsl_type(VDEVICE *p, INT style)
#else
V\_set\_line\_type
#endif
{
	p->intin[0] = style;
	p->contrl[0] = 15;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
void Vsl_udsty(VDEVICE *p, INT pattern)
#else
V\_set\_line\_user\_defined\_style
#endif
{
	p->intin[0] = pattern;
	p->contrl[0] = 113;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
}

#if TEXDOCU
INT Vsl_width(VDEVICE *p, INT width)
#else
V\_set\_line\_width
#endif
{
	p->ptsin[0] = width;
	p->ptsin[1] = 0;
	p->contrl[0] = 16;
	p->contrl[1] = 1,
	p->contrl[3] = 0;
	Vdi_draw(p);
	return(p->ptsout[0]);
}

#if TEXDOCU
INT Vsl_color(VDEVICE *p, INT color_index)
#else
V\_set\_line\_color
#endif
{
	p->intin[0] = color_index;
	p->contrl[0] = 17;
	p->contrl[1] = 3;
	Vdi_draw(p);
	return (p->intout[0]);
}

#if TEXDOCU
void Vsl_ends(VDEVICE *p, INT beg_style, INT end_style)
#else
V\_set\_line\_ends
#endif
{
	p->intin[0] = beg_style;
	p->intin[1] = end_style;
	p->contrl[0] = 108;
	p->contrl[1] = 0;
	p->contrl[3] = 2;
	Vdi_draw(p);
}

#if TEXDOCU
INT Vsm_type(VDEVICE *p, INT symbol)
#else
V\_set\_marker\_type
#endif
{
	p->intin[0] = symbol;
	p->contrl[0] = 18;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsm_height(VDEVICE *p, INT height)
#else
V\_set\_marker\_height
#endif
{
	p->ptsin[0] = 0;
	p->ptsin[1] = height;
	p->contrl[0] = 19;
	p->contrl[1] = 1;
	p->contrl[3] = 0;
	Vdi_draw(p);
	return(p->ptsout[1]);
}

#if TEXDOCU
INT Vsm_color(VDEVICE *p, INT color_index)
#else
V\_set\_marker\_color
#endif
{
	p->intin[0] = color_index;
	p->contrl[0] = 20;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
void Vst_alignment(VDEVICE *p, INT hor_in, 
	INT vert_in, INT *hor_out, INT *vert_out)
#else
V\_set\_text\_alignment:
vetrical and horizontal text alignment. 
hor\_in: 0= left aligned, 1 = centered, 2 = right aligned.
vert\_in: 0= bottom aligned, 1 = middle aligned, 2 = top aligned.
#endif
{
	p->intin[0] = hor_in;
	p->intin[1] = vert_in;
	p->contrl[0] = 39;
	p->contrl[1] = 0;
	p->contrl[3] = 2;
	Vdi_draw(p);
	*hor_out = p->intout[0];
	*vert_out = p->intout[1];
}

#if TEXDOCU
void Vst_height(VDEVICE *p, INT height, 
	INT *char_width, INT *char_height, 
	INT *cell_width, INT *cell_height)
#else
V\_set\_text\_height: unused
#endif
{
	p->ptsin[0] = 0;
	p->ptsin[1] = height;
	p->contrl[0] = 12;
	p->contrl[1] = 1;
	p->contrl[3] = 0;
	Vdi_draw(p);
	*char_width = p->ptsout[0];
	*char_height = p->ptsout[1];
	*cell_width = p->ptsout[2];
	*cell_height = p->ptsout[3];
}

#if TEXDOCU
INT Vst_point(VDEVICE *p, INT point, 
	INT *char_width, INT *char_height, 
	INT *cell_width, INT *cell_height)
#else
V\_set\_text\_height\_in\_points: unused
#endif
{
	p->intin[0] = point;
	p->contrl[0] = 107;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	*char_width = p->ptsout[0];
	*char_height = p->ptsout[1];
	*cell_width = p->ptsout[2];
	*cell_height = p->ptsout[3];
	return(p->intout[0]);
}

#if TEXDOCU
INT Vst_rotation(VDEVICE *p, INT angle)
#else
V\_set\_text\_rotation: unused
#endif
{
	p->intin[0] = angle;
	p->contrl[0] = 13;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vst_font(VDEVICE *p, INT font)
#else
V\_set\_text\_font: unused
#endif
{
	p->intin[0] = font;
	p->contrl[0] = 21;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vst_color(VDEVICE *p, INT color_index)
#else
V\_set\_text\_color: unused
#endif
{
	p->intin[0] = color_index;
	p->contrl[0] = 22;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vst_effects(VDEVICE *p, INT effect)
#else
V\_set\_text\_effects: unused
#endif
{
	p->intin[0] = effect;
	p->contrl[0] = 106;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_interior(VDEVICE *p, INT style)
#else
V\_set\_fill\_interior:
in 1/100th, 0= none (used for pie-drawing)
#endif
{
	p->intin[0] = style;
	p->contrl[0] = 23;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_style(VDEVICE *p, INT style_index)
#else
V\_set\_fill\_style:
#endif
{
	p->intin[0] = style_index;
	p->contrl[0] = 24;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_color(VDEVICE *p, INT color_index)
#else
V\_set\_fill\_color: 
0 = white, 1 = black
#endif
{
	p->intin[0] = color_index;
	p->contrl[0] = 25;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_shape(VDEVICE *p, INT shape)
#else
V\_set\_fill\_shape: unused
#endif
{
	p->intin[0] = shape;
	p->contrl[0] = 26;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_outline(VDEVICE *p, INT outline)
#else
V\_set\_fill\_outline: 
#endif
{
	p->intin[0] = outline;
	p->contrl[0] = 27;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_nofill(VDEVICE *p, INT nofill)
#else
V\_set\_fill\_nofill: 
#endif
{
	p->intin[0] = nofill;
	p->contrl[0] = 28;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vsf_perimeter(VDEVICE *p, INT per_vis)
#else
V\_set\_fill\_perimeter:
#endif
{
	p->intin[0] = per_vis;
	p->contrl[0] = 104;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if 0
void Vsf_updat(VDEVICE *p, INT pfill_pat[], INT planes)
{
	Vi_intin(p, pfill_pat);
	p->contrl[0] = 112;
	p->contrl[1] = 0;
	p->contrl[3] = 16 * planes;
	Vdi_draw(p);
	Vi_intin(p, p->intin);
}
#endif

#if TEXDOCU
INT Vst_boxed(VDEVICE *p, INT boxed)
#else
V\_set\_text\_boxed: draw text in a box
#endif
{
	p->intin[0] = boxed;
	p->contrl[0] = 110;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

#if TEXDOCU
INT Vst_overwrite(VDEVICE *p, INT overwrite)
#else
V\_set\_text\_overwrite: unused
#endif
{
	p->intin[0] = overwrite;
	p->contrl[0] = 111;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	Vdi_draw(p);
	return(p->intout[0]);
}

/*
 * VDEVICE - CONTROL.C
 */

#if FALSE
void v_opnwk(INT work_in[], INT *handle, INT work_out[])
{
	i_intin(work_in);
	i_intout(work_out);
	i_ptsout(work_out + 45);
	contrl[0] = 1;
	contrl[1] = 0;
	contrl[3] = 11;
	vdi();
	*handle = contrl[6];
	i_intin(intin);
	i_intout(intout);
	i_ptsout(ptsout);
	i_ptsin(ptsin);
}

void v_clswk(INT handle)
{
	contrl[0] = 2;
	contrl[1] = 0;
	contrl[3] = 0;
	contrl[6] = handle;
	vdi();
}

void v_opnvwk(INT work_in[], INT *handle, INT work_out[])
{
	i_intin(work_in);
	i_intout(work_out);
	i_ptsout(work_out + 45);
	contrl[0] = 100;
	contrl[1] = 0;
	contrl[2] = 11;
	contrl[6] = *handle;
	vdi();
	*handle = contrl[6];
	i_intin(intin);
	i_intout(intout);
	i_ptsout(ptsout);
	i_ptsin(ptsin);
}

void v_clsvwk(INT handle)
{
	contrl[0] = 101;
	contrl[1] = 0;
	contrl[3] = 0;
	contrl[6] = handle;
	vdi();
}

void V_clrwk(VDEVICE *p)
{
	p->contrl[0] = 3;
	p->contrl[1] = 0;
	p->contrl[2] = 0;
	Vdi_draw(p);
}

void V_updwk(VDEVICE *p)
{
	p->contrl[0] = 4;
	p->contrl[1] = 0;
	p->contrl[3] = 0;
	Vdi_draw(p);
}

INT Vst_load_fonts(VDEVICE *p, INT select)
{
	p->contrl[0] = 119;
	p->contrl[1] = 0;
	p->contrl[3] = 1;
	p->intin[0] = select;
	Vdi_draw(p);
	return(p->intout[0]);
}

void Vst_unload_fonts(VDEVICE *p, INT select)
{
	p->contrl[0] = 120;
	p->contrl[1] = 0;
	p->contrl[3] = 0;
	Vdi_draw(p);
}

void Vs_clip(VDEVICE *p, INT clip_flag, INT pxyarray[])
{
	Vi_ptsin(p, pxyarray);
	p->intin[0] = clip_flag;
	p->contrl[0] = 129;
	p->contrl[1] = 2;
	p->contrl[3] = 1;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}
#endif

/*
 * VDEVICE - OUTPUT.C
 */

#if TEXDOCU
void V_pline(VDEVICE *p, INT count, INT pxyarray[])
#else
V\_poly\_line
#endif
{
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 6;
	p->contrl[1] = count;
	p->contrl[3] = 0;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_bezier(VDEVICE *p, INT count, INT pxyarray[])
#else
V\_bezier
#endif
{
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 7;
	p->contrl[1] = count;
	p->contrl[3] = 0;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_pmarker(VDEVICE *p, INT count, INT pxyarray[])
#else
V\_put\_marker: unused
#endif
{
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 7;
	p->contrl[1] = count;
	p->contrl[3] = 0;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_gtext(VDEVICE *p, INT x, INT y, BYTE *string)
#else
V\_graphical\_text: text at $(x,y)$
#endif
{
	INT i;
	
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	i = 0;
	while (*string)
		p->intin[i++] = *string++;
	p->contrl[0] = 8;
	p->contrl[1] = 1;
	p->contrl[3] = i;
	Vdi_draw(p);
}

#if TEXDOCU
void V_fillarea(VDEVICE *p, INT count, INT pxyarray[])
#else
V\_fill\_area
#endif
{
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 9;
	p->contrl[1] = count;
	p->contrl[3] = 0;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_cellarray(VDEVICE *p, INT pxyarray[4], 
	INT row_length, INT el_used, INT num_rows, 
	INT wrt_mode, INT colarray[])
#else
unused
#endif
{
	Vi_intin(p, colarray);
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 10;
	p->contrl[1] = 2;
	p->contrl[3] = row_length * num_rows;
	p->contrl[7] = row_length;
	p->contrl[8] = el_used;
	p->contrl[9] = num_rows;
	p->contrl[10] = wrt_mode;
	Vdi_draw(p);
	Vi_intin(p, p->intin);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_contourfill(VDEVICE *p, INT x, INT y, INT index)
#else
unused
#endif
{
	p->intin[0] = index;
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->contrl[0] = 103;
	p->contrl[1] = 1;
	p->contrl[3] = 1;
	Vdi_draw(p);
}

#if TEXDOCU
void Vr_recfl(VDEVICE *p, INT *pxyarray)
#else
unused
#endif
{
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 114;
	p->contrl[1] = 2;
	p->contrl[3] = 0;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_bar(VDEVICE *p, INT pxyarray[])
#else
unused
#endif
{
	Vi_ptsin(p, pxyarray);
	p->contrl[0] = 11;
	p->contrl[1] = 2;
	p->contrl[3] = 0;
	p->contrl[5] = 1;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_arc(VDEVICE *p, INT x, INT y, INT radius, 
	INT begang, INT endang)
#else
V\_arc: arc around $(x,y)$ of radius form begang to endang, 
where begang and endang are in 1/10-th degree.
#endif
{
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = 0;
	p->ptsin[3] = 0;
	p->ptsin[4] = 0;
	p->ptsin[5] = 0;
	p->ptsin[6] = radius;
	p->ptsin[7] = 0;
	p->intin[0] = begang;
	p->intin[1] = endang;
	p->contrl[0] = 11;
	p->contrl[1] = 4;
	p->contrl[3] = 2;
	p->contrl[5] = 2;
	Vdi_draw(p);
}

#if TEXDOCU
void V_pie(VDEVICE *p, INT x, INT y, INT radius, INT begang, INT endang)
#else
V\_pie: pie around $(x,y)$ of radius form begang to endang, 
where begang and endang are in 1/10-th degree.
#endif
{
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = 0;
	p->ptsin[3] = 0;
	p->ptsin[4] = 0;
	p->ptsin[5] = 0;
	p->ptsin[6] = radius;
	p->ptsin[7] = 0;
	p->intin[0] = begang;
	p->intin[1] = endang;
	p->contrl[0] = 11;
	p->contrl[1] = 4;
	p->contrl[3] = 2;
	p->contrl[5] = 3;
	Vdi_draw(p);
}

#if TEXDOCU
void V_circle(VDEVICE *p, INT x, INT y, INT radius)
#else
V\_circle: full circle around $(x,y)$ of radius.
#endif
{
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = 0;
	p->ptsin[3] = 0;
	p->ptsin[4] = radius;
	p->ptsin[5] = 0;
	p->contrl[0] = 11;
	p->contrl[1] = 3;
	p->contrl[3] = 0;
	p->contrl[5] = 4;
	Vdi_draw(p);
}

#if TEXDOCU
void V_ellipse(VDEVICE *p, INT x, INT y, 
	INT xradius, INT yradius)
#else
unused
#endif
{
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = xradius;
	p->ptsin[3] = yradius;
	p->contrl[0] = 11;
	p->contrl[1] = 2;
	p->contrl[3] = 0;
	p->contrl[5] = 5;
	Vdi_draw(p);
}

#if TEXDOCU
void V_ellarc(VDEVICE *p, INT x, INT y, 
	INT xradius, INT yradius, INT begang, INT endang)
#else
unused
#endif
{
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = xradius;
	p->ptsin[3] = yradius;
	p->intin[0] = begang;
	p->intin[1] = endang;
	p->contrl[0] = 11;
	p->contrl[1] = 2;
	p->contrl[3] = 2;
	p->contrl[5] = 6;
	Vdi_draw(p);
}

#if TEXDOCU
void V_ellpie(VDEVICE *p, INT x, INT y, 
	INT xradius, INT yradius, INT begang, INT endang)
#else
unused
#endif
{
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = xradius;
	p->ptsin[3] = yradius;
	p->intin[0] = begang;
	p->intin[1] = endang;
	p->contrl[0] = 11;
	p->contrl[1] = 2;
	p->contrl[3] = 2;
	p->contrl[5] = 7;
	Vdi_draw(p);
}

#if TEXDOCU
void V_rbox(VDEVICE *p, INT xyarray[])
#else
unused
#endif
{
	Vi_ptsin(p, xyarray);
	p->contrl[0] = 11;
	p->contrl[1] = 2;
	p->contrl[3] = 0;
	p->contrl[5] = 8;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_rfbox(VDEVICE *p, INT xyarray[])
#else
unused
#endif
{
	Vi_ptsin(p, xyarray);
	p->contrl[0] = 11;
	p->contrl[1] = 2;
	p->contrl[3] = 0;
	p->contrl[5] = 9;
	Vdi_draw(p);
	Vi_ptsin(p, p->ptsin);
}

#if TEXDOCU
void V_justified(VDEVICE *p, INT x, INT y, BYTE *string, 
	INT length, INT word_space, INT char_space)
#else
unused
#endif
{
	INT *tmp;
	
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = length;
	p->ptsin[3] = 0;
	p->intin[0] = word_space;
	p->intin[1] = char_space;
	tmp = &(p->intin[2]);
	while (*string)
		*tmp++ = *string++;
	p->contrl[0] = 1;
	p->contrl[1] = 2;
	p->contrl[3] = (INT)(tmp - p->intin[2]);
	p->contrl[5] = 10;
	Vdi_draw(p);
}

#if TEXDOCU
void V_pie_text(VDEVICE *p, INT x, INT y, INT radius, 
	INT begang, INT endang, BYTE *string)
#else
#endif
{
	INT i;
	
	p->ptsin[0] = x;
	p->ptsin[1] = y;
	p->ptsin[2] = 0;
	p->ptsin[3] = 0;
	p->ptsin[4] = 0;
	p->ptsin[5] = 0;
	p->ptsin[6] = radius;
	p->ptsin[7] = 0;
	p->intin[0] = begang;
	p->intin[1] = endang;
	i = 0;
	while (*string)
		p->intin[2 + i++] = *string++;
	p->contrl[0] = 11;
	p->contrl[1] = 4;
	p->contrl[3] = 2 + i;
	p->contrl[5] = 11;
	Vdi_draw(p);
}

#if TEXDOCU
INT VdevSExtra(VDEVICE *p, INT num, INT extra)
#endif
{
	if (p == NIL) {
		return error("VdevSExtra() p == NIL");
		}
	if (num >= 0 && num < 5)
		p->extra[num] = extra;
	return(TRUE);
}

#if TEXDOCU
INT VdevGExtra(VDEVICE *p, INT num, INT *extra)
#endif
{
	if (p == NIL) {
		return error("VdevGExtra() p == NIL");
		}
	if (num >= 0 && num < 5)
		*extra = p->extra[num];
	return(TRUE);
}

#if TEXDOCU
INT VdevI2(VDEVICE *p, INT *contrl, INT *intin, INT *ptsin, 
	INT *intout, INT *ptsout)
#endif
{
	if (p == NIL) {
		return error("VdevI2() p == NIL");
		}
	p->contrl = contrl;
	p->intin = intin;
	p->ptsin = ptsin;
	p->intout = intout;
	p->ptsout = ptsout;
	p->vdiblk[0] = contrl;
	p->vdiblk[1] = intin;
	p->vdiblk[2] = ptsin;
	p->vdiblk[3] = intout;
	p->vdiblk[4] = ptsout;
	return(TRUE);
}

#if TEXDOCU
INT VdevGArrays(VDEVICE *p, INT **pcontrl, INT **pintin, 
	INT **pptsin, INT **pintout, INT **pptsout)
#endif
/* Es werden nur die arrays ausgelesen, deren Doppelpointer 
 * gesetzt (!=NIL) ist. */
{
	if (p == NIL) {
		return error("VdevGArrays() p == NIL");
		}
	if (pcontrl)
		*pcontrl = p->vdiblk[0];
	if (pintin)
		*pintin = p->vdiblk[1];
	if (pptsin)
		*pptsin = p->vdiblk[2];
	if (pintout)
		*pintout = p->vdiblk[3];
	if (pptsout)
		*pptsout = p->vdiblk[4];
	return(TRUE);
}

#undef DEBUG_USER2VDEV

#if TEXDOCU
INT user2vdev(VDEVICE *vdev, double dx, double dy, double dz, 
	INT *ix, INT *iy, double *extrema)
#endif
/* In: extrema[0] <= dx <= extrema[1]
 *     extrema[2] <= dy <= extrema[3]
 *     extrema[4] <= dz <= extrema[5]
 * Out: Pixel - Koordinaten ins vdev */
{
	INT size_x, size_y, size_z, ix1, iy1, iz, iz0;
	INT f_use_z = TRUE;
	double a, b, c;

	if (vdev == NIL || ix == NIL || iy == NIL) {
		return error("user2vdev() args NIL");
		}
	dx -= extrema[0];
	dy -= extrema[2];
	dz -= extrema[4];
	a = extrema[1] - extrema[0];
	b = extrema[3] - extrema[2];
	c = extrema[5] - extrema[4];
	if (ABS(a) > M_EPS)
		dx /= a;
	else
		return error("user2vdev() ABS(a) <= M_EPS");
	if (ABS(b) > M_EPS)
		dy /= b;
	else
		return error("user2vdev() ABS(b) <= M_EPS");
	if (ABS(c) > M_EPS)
		dz /= c;
	else
		f_use_z = FALSE;
	size_x = vdev->co[2] - vdev->co[0];
	size_y = vdev->co[3] - vdev->co[1];
	dx *= (double) size_x;
	dy *= (double) size_y;
	ix1 = (INT) dx;
	iy1 = (INT) dy;
	if (f_use_z) {
		size_z = (size_x + size_y) >> 1;
		dz *= (double) size_z;
		dz *= M_SQRT2INV;
		iz0 = (INT) ((double)size_z * (0. - extrema[4]) / c * M_SQRT2INV);
		iz = -(INT) dz + iz0;
		ix1 += iz;
		iy1 += iz;
		}
	iy1 = /* size_y - */ iy1;
	ix1 += vdev->co[0]; /*(lw->interior.x - lw->hwork.p) */
	iy1 += vdev->co[1]; /*(lw->interior.y - lw->vwork.p) */
	*ix = /* lw_grid_cont->border_width + */ ix1;
	*iy = /* lw_grid_cont->border_height + */ iy1;
#ifdef DEBUG_USER2VDEV
	printf("user2vdev(): dx = %lf dy = %lf ix = %ld iy = %ld\n", dx, dy, *ix, *iy);
#endif
	return TRUE;
}

#if TEXDOCU
INT vdev2user(VDEVICE *vdev, INT ix, INT iy, double *dx, double *dy, 
	double *extrema)
#endif
{
	double dx1, dy1;
	INT size_x, size_y, ix1, iy1;

	if (vdev == NIL || dx == NIL || dy == NIL) {
		return error("vdev2user() args NIL");
		}
	ix1 = ix; /* - lw_grid_cont->border_width */
	iy1 = iy; /* - lw_grid_cont->border_height */
	ix1 -= vdev->co[0]; /* (lw->interior.x - lw->hwork.p) */
	iy1 -= vdev->co[1]; /* (lw->interior.y - lw->vwork.p) */
	size_x = vdev->co[2] - vdev->co[0];
	size_y = vdev->co[3] - vdev->co[1];
	iy1 = /* size_y - */ iy1;
	dx1 = (double) ix1 / (double)size_x;
	dy1 = (double) iy1 / (double)size_y;
	dx1 *= (extrema[1] - extrema[0]);
	dy1 *= (extrema[3] - extrema[2]);
	dx1 += extrema[0];
	dy1 += extrema[2];
	*dx = dx1;
	*dy = dy1;
	return TRUE;
}

#if TEXDOCU
INT user2xywh(double dx, double dy, double dz, 
	INT *ix, INT *iy, double *extrema, INT *xywh)
#endif
/* In: extrema[0] <= dx <= extrema[1]
 *     extrema[2] <= dy <= extrema[3]
 *     extrema[4] <= dz <= extrema[5]
 * Out: Pixel - Koordinaten im Bereich von xywh[]. */
{
	INT size_x, size_y, size_z, ix1, iy1, iz, iz0;
	INT f_use_z = TRUE;
	double a, b, c;

	if (ix == NIL || iy == NIL || xywh == NIL) {
		return error("user2xywh() args NIL");
		}
	dx -= extrema[0];
	dy -= extrema[2];
	dz -= extrema[4];
	a = extrema[1] - extrema[0];
	b = extrema[3] - extrema[2];
	c = extrema[5] - extrema[4];
	if (ABS(a) > M_EPS)
		dx /= a;
	else
		return error("user2xywh() ABS(a) <= M_EPS");
	if (ABS(b) > M_EPS)
		dy /= b;
	else
		return error("user2xywh() ABS(b) <= M_EPS");
	if (ABS(c) > M_EPS)
		dz /= c;
	else
		f_use_z = FALSE;
	size_x = xywh[2] /* vdev->co[2] - vdev->co[0] */;
	size_y = xywh[3] /* vdev->co[3] - vdev->co[1] */;
	dx *= (double) size_x;
	dy *= (double) size_y;
	ix1 = (INT) dx;
	iy1 = (INT) dy;
	if (f_use_z) {
		size_z = (size_x + size_y) >> 1;
		dz *= (double) size_z;
		dz *= M_SQRT2INV;
		iz0 = (INT) ((double)size_z * (0. - extrema[4]) / c * M_SQRT2INV);
		iz = -(INT) dz + iz0;
		ix1 += iz;
		iy1 += iz;
		}
	iy1 = /* size_y - */ iy1;
	ix1 += xywh[0]; /*(lw->interior.x - lw->hwork.p) */
	iy1 += xywh[1]; /*(lw->interior.y - lw->vwork.p) */
	*ix = /* lw_grid_cont->border_width + */ ix1;
	*iy = /* lw_grid_cont->border_height + */ iy1;
	return TRUE;
}

#if TEXDOCU
INT xywh2user(INT ix, INT iy, double *dx, double *dy, 
	double *extrema, INT *xywh)
#endif
{
	double dx1, dy1;
	INT size_x, size_y, ix1, iy1;

	if (dx == NIL || dy == NIL || xywh == NIL) {
		return error("xywh2user() args NIL");
		}
	ix1 = ix; /* - lw_grid_cont->border_width */
	iy1 = iy; /* - lw_grid_cont->border_height */
	ix1 -= xywh[0]; /* (lw->interior.x - lw->hwork.p) */
	iy1 -= xywh[1]; /* (lw->interior.y - lw->vwork.p) */
	size_x = xywh[2];
	size_y = xywh[3];
	iy1 = /* size_y -  */ iy1;
	dx1 = (double) ix1 / (double)size_x;
	dy1 = (double) iy1 / (double)size_y;
	dx1 *= (extrema[1] - extrema[0]);
	dy1 *= (extrema[3] - extrema[2]);
	dx1 += extrema[0];
	dy1 += extrema[2];
	*dx = dx1;
	*dy = dy1;
	return TRUE;
}

#if FALSE
INT VdevSType(VDEVICE *p, INT type)
{
	if (p == NIL) {
		Srfs("VdevSType", "p == NIL");
		return(FALSE);
		}
	p->type = type;
	return(TRUE);
}

INT VdevGType(VDEVICE *p, INT *type)
{
	if (p == NIL) {
		Srfs("VdevGType", "p == NIL");
		return(FALSE);
		}
	*type = p->type;
	return(TRUE);
}

INT VdevSvdi_draw(VDEVICE *p, INT (*vdi_draw)(VDEVICE *vdev))
{
	if (p == NIL) {
		Srfs("VdevSvdi_draw", "p == NIL");
		return(FALSE);
		}
	p->vdi_draw = vdi_draw;
	return(TRUE);
}

INT VdevGvdi_draw(VDEVICE *p, INT (**vdi_draw)(VDEVICE *vdev))
{
	if (p == NIL || vdi_draw == NIL) {
		Srfs("VdevGvdi_draw", "args NIL");
		return(FALSE);
		}
	*vdi_draw = p->vdi_draw;
	return(TRUE);
}
#endif

#if TEXDOCU
INT place_lattice(VECTOR_OP nl, VECTOR_OP orbit_size, 
	INT size_x, INT size_y, 
	VECTOR_OP Px, VECTOR_OP Py, VECTOR_OP Os, 
	INT f_upside_down, INT f_v)
#endif
{
	MEM_OB vbp_plaz, vbp_o_dx;
	INT *plazierung;
	INT *o_dx;
	INT n, i, x, y, o;
	INT l, s0 = 0, s1 = 0, S0;
	INTEGER_OB S;
	
	l = orbit_size->s_li();
	if (l > 4) {
		s0 = orbit_size->s_ii(0);
		s1 = orbit_size->s_ii(l - 1);
		orbit_size->sum(&S);
		S0 = S.s_i() >> 2;
		orbit_size->m_ii(0, S0);
		orbit_size->m_ii(l - 1, S0);
		}
	vbp(nl, orbit_size, 
		size_x, size_y, 
		&vbp_plaz, &vbp_o_dx, f_upside_down);
	if (l > 4) {
		orbit_size->m_ii(0, s0);
		orbit_size->m_ii(l - 1, s1);
		}
	n = nl->s_ii(0) - 1;
	plazierung = (INT *) vbp_plaz.ob_self.ob_charpointer;
	o_dx = (INT *) vbp_o_dx.ob_self.ob_charpointer;
	Px->m_il_n(n);
	Py->m_il_n(n);
	Os->m_il_n(n);
#if 0
	INT max_dy;
	
	max_dy = 0;
	for (i = 1; i < n - 1; i++) {
		max_dy = MAXIMUM(max_dy, 
			plazierung[2 * (i + 1) + 1] - plazierung[2 * i + 1]);
		}
	printf("max_dy = %ld\n", max_dy);
	plazierung[2 * 0 + 1] -= max_dy;
	plazierung[2 * (n - 1) + 1] += max_dy;
#endif
	for (i = 0; i < n; i++) {
		x = plazierung[2 * i];
		y = plazierung[2 * i + 1];
		o = o_dx[i];
		Px->m_ii(i, x);
		Py->m_ii(i, y);
		Os->m_ii(i, o);
		}
	if (f_v) {
		printf("place_lattice: placement:\n");
		for (i = 0; i < n; i++) {
			printf("%ld: %ld %ld %ld\n", 
				i, Px->s_ii(i), Py->s_ii(i), Os->s_ii(i));
			}
		}
	return OK;
}

#if TEXDOCU
INT vbp(VECTOR_OP nl, VECTOR_OP orbit_size, 
	INT x_pix, INT y_pix, 
	MEM_OP plazierung, MEM_OP orbit_dx, 
	INT f_upside_down)
#endif
{
	INT i, len_plaz, len_o_dx;
	ULONG *adj_list = NIL;
	ULONG *os = NIL;
	float *plazierung1 = NIL;
	float *orbit_dx1 = NIL;
	double extrema[6];
	double x, y, dx;
	INT *plaz = NIL;
	INT *o_dx = NIL;
	
	if (nl->s_ii(0) - 1 != orbit_size->s_li())
		return error("nl->s_ii(0) - 1 != orbit_size->s_li()\n");
	adj_list = (ULONG *) my_malloc(sizeof(ULONG) * nl->s_li(), "vbp");
	os = (ULONG *) my_malloc(sizeof(ULONG) * orbit_size->s_li(), "vbp");
	if (adj_list == NIL || os == NIL) {
		return error("vbp() no memory");
		}
	for (i = 0; i < nl->s_li(); i++)
		adj_list[i] = nl->s_ii(i);
	for (i = 0; i < orbit_size->s_li(); i++)
		os[i] = orbit_size->s_ii(i);
	len_plaz = sizeof(INT) * (nl->s_ii(0) - 1) * 2;
	len_o_dx = sizeof(INT) * orbit_size->s_li();
	plaz = (INT *) my_malloc(len_plaz, "vbp");
	o_dx = (INT *) my_malloc(len_o_dx, "vbp");
	if (plaz == NIL || o_dx == NIL) {
		return error("vbp() no memory");
		}
	if (verband_plazieren(adj_list, os, 
		&plazierung1, &orbit_dx1) != OK)
		return error("error in verband_plazieren()\n");
	
	if (f_upside_down) {
		/* Spiegelung: oben und unten vertauschen */
		for (i = 0; i < (INT) adj_list[0] - 1; i++) {
			plazierung1[2 * i + 1] = 1. - plazierung1[2 * i + 1];
			}
		}
	/*x_streckung(n, plazierung1);*/
	/* for (i = 0; i < (INT) adj_list[0] - 1; i++) {
		printf("%ld: %f %f\n", i, 
			plazierung1[2 * i], plazierung1[2 * i + 1]);
		} */
	/* t_gerichtet_ungerichtet(nl[i], &erg); */
	/* free(plazierung1);
	free(erg); */
	extrema[0] = 0.;
	extrema[1] = 1.;
	extrema[2] = 0.;
	extrema[3] = 1.;
	extrema[4] = -1.;
	extrema[5] = 1.;
	for (i = 0; i < (INT) adj_list[0] - 1; i++) {
		x = plazierung1[2 * i];
		y = plazierung1[2 * i + 1];
		dx = orbit_dx1[i];
		plaz[2 * i] = (INT)( (double) x_pix * x );
		plaz[2 * i + 1] = (INT)( (double) y_pix * y );
		o_dx[i] = (INT)( (double) x_pix * dx );
		}
	my_free(plazierung1);
	my_free(orbit_dx1);
	plazierung->init(len_plaz, (BYTE *) plaz);
	orbit_dx->init(len_o_dx, (BYTE *) o_dx);
	
	if (adj_list) {
		my_free(adj_list);
		adj_list = NIL;
		}
	if (os) {
		my_free(os);
		os = NIL;
		}
	if (plaz) {
		my_free(plaz);
		plaz = NIL;
		}
	if (o_dx) {
		my_free(o_dx);
		o_dx = NIL;
		}
	return OK;
}

#if TEXDOCU
void draw_box(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1)
#endif
{
	INT pxy[32];

	pxy[0] = x0;
	pxy[1] = y0;
	pxy[2] = x1;
	pxy[3] = y0;
	pxy[4] = x1;
	pxy[5] = y1;
	pxy[6] = x0;
	pxy[7] = y1;
	pxy[8] = x0;
	pxy[9] = y0;
	V_pline(vdev, 5, pxy);
}

#if TEXDOCU
void draw_kreuz(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1)
#endif
{
	INT pxy[32];

	pxy[0] = x0;
	pxy[1] = y0;
	pxy[2] = x1;
	pxy[3] = y1;
	V_pline(vdev, 2, pxy);
	pxy[0] = x0;
	pxy[1] = y1;
	pxy[2] = x1;
	pxy[3] = y0;
	V_pline(vdev, 2, pxy);
}

#if TEXDOCU
void draw_kreuz2(VDEVICE *vdev, INT x, INT y, INT rad)
#endif
{
	INT pxy[32];

	pxy[0] = x - rad;
	pxy[1] = y;
	pxy[2] = x + rad;
	pxy[3] = y;
	V_pline(vdev, 2, pxy);
	pxy[0] = x;
	pxy[1] = y - rad;
	pxy[2] = x;
	pxy[3] = y + rad;
	V_pline(vdev, 2, pxy);
}

#if TEXDOCU
void draw_mol(VDEVICE *vdev, INT x, INT y, BYTE *text)
#endif
{
	INT rad, off_x, off_y;
	
#ifdef SYSTEMUNIX
	rad = 5;
	off_x = 15;
	off_y = 20;

#if 0
	if (vdev->type == 0) { /* ps output */
		rad = 40;
		}
#endif

#else
	rad = 4;
	off_x = 5;
	off_y = 8;
#endif
	V_pie(vdev, x, y, rad, 0, 3600);
	V_gtext(vdev, x + off_x, y + off_y, text);
}

#if TEXDOCU
void bind_mol(VDEVICE *vdev, INT x1, INT y1, INT x2, INT y2, INT n)
#endif
{
	INT k;
	float dx, dy, winkel, laenge;

	if (n == 1) /* Einfachbindung */
		my_line(vdev, x1, y1, x2, y2);
	else { /* Mehrfachbindung */
		dx = x2 - x1;
		dy = -y2 + y1;

		if (dx != 0)
			winkel = atan(dy / dx);
		else
			if (dy > 0)
				winkel = M_PI_2;
			else
				winkel = 3. * M_PI_2;

		winkel += M_PI_2;
		laenge = (n - 1) * 1.5;
		dx = cos(winkel) * laenge;
		dy = sin(winkel) * laenge;

		x1 += (INT) dx;
		y1 -= (INT) dy;
		x2 += (INT) dx;
		y2 -= (INT) dy;
		my_line(vdev, x1, y1, x2, y2);

		winkel += M_PI;
		dx = cos(winkel) * 3.;
		dy = sin(winkel) * 3.;

		for (k = 1; k < n; k++) {
			x1 += (INT) dx;
			y1 -= (INT) dy;
			x2 += (INT) dx;
			y2 -= (INT) dy;
			my_line(vdev, x1, y1, x2, y2);
			}
		}
}

#if TEXDOCU
void Bind_mol(VDEVICE *vdev, INT dist, INT x1, INT y1, INT x2, INT y2, INT n)
#endif
{
	INT k;
	float dx, dy, winkel, laenge;

	if (n == 1) /* Einfachbindung */
		my_line(vdev, x1, y1, x2, y2);
	else { /* Mehrfachbindung */
		dx = x2 - x1;
		dy = -y2 + y1;

		if (dx != 0)
			winkel = atan(dy / dx);
		else
			if (dy > 0)
				winkel = M_PI_2;
			else
				winkel = 3. * M_PI_2;

		winkel += M_PI_2;
		laenge = ((n - 1) * (double) dist) / 2.;
		dx = cos(winkel) * laenge;
		dy = sin(winkel) * laenge;

		x1 += (INT) dx;
		y1 -= (INT) dy;
		x2 += (INT) dx;
		y2 -= (INT) dy;
		my_line(vdev, x1, y1, x2, y2);

		winkel += M_PI;
		dx = cos(winkel) * (double) dist;
		dy = sin(winkel) * (double) dist;

		for (k = 1; k < n; k++) {
			x1 += (INT) dx;
			y1 -= (INT) dy;
			x2 += (INT) dx;
			y2 -= (INT) dy;
			my_line(vdev, x1, y1, x2, y2);
			}
		}
}

#if TEXDOCU
INT my_line(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1)
#endif
{
	INT pxy[4];
	
	pxy[0] = x0;
	pxy[1] = y0;
	pxy[2] = x1;
	pxy[3] = y1;
	V_pline(vdev, 2, pxy);
	return TRUE;
}

#if TEXDOCU
double cos_grad(double phi)
#endif
{
	double x;

	x = (phi * PI) / 180.;
	return cos(x);
}

#if TEXDOCU
double sin_grad(double phi)
#endif
{
	double x;

	x = (phi * PI) / 180.;
	return sin(x);
}

#if TEXDOCU
double tan_grad(double phi)
#endif
{
	double x;

	x = (phi * PI) / 180.;
	return tan(x);
}

#if TEXDOCU
double atan_grad(double x)
#endif
{
	double y, phi;

	y = atan(x);
	phi = (y * 180.) / PI;
	return phi;
}

#if TEXDOCU
INT draw_label(VDEVICE *vdev, BYTE *text, BYTE *align, INT *x, INT *y, INT idx)
#endif
{
	INT h_align = 1, v_align = 1;
	INT l, i, j;
	BYTE c;
	
	l = strlen(align);
	for (i = 0; i < l; i++) {
		c = align[i];
		if (c == 'r')
			h_align = 2;
		else if (c == 'l')
			h_align = 0;
		else if (c == 'b')
			v_align = 0;
		else if (c == 't')
			v_align = 2;
		else {
			printf("vdi.C: draw_label: unknown alignment character %c\n", c);
			}
		}
	Vst_alignment(vdev, h_align, v_align, &j, &j);
	V_gtext(vdev, x[idx], y[idx], text);
	return OK;
}

#if TEXDOCU
INT draw_dot(VDEVICE *vdev, INT *x, INT *y, INT idx, INT rad)
#endif
{
	Vsf_interior(vdev, GREY_INTERIOR);
	Vsf_color(vdev, GREY_COLOR);
	V_pie(vdev, x[idx], y[idx], rad, 0, 0);
	Vsf_interior(vdev, 0);
	V_arc(vdev, x[idx], y[idx], rad, 0 /* begang */, 3600 /* endang */);
	return OK;
}

#if TEXDOCU
INT draw_dot_white(VDEVICE *vdev, INT *x, INT *y, INT idx, INT rad)
#endif
{
	Vsf_interior(vdev, WHITE_INTERIOR);
	Vsf_color(vdev, WHITE_COLOR);
	V_pie(vdev, x[idx], y[idx], rad, 0, 0);
	Vsf_interior(vdev, 0);
	V_arc(vdev, x[idx], y[idx], rad, 0 /* begang */, 3600 /* endang */);
	return OK;
}

#if TEXDOCU
INT draw_dot_black(VDEVICE *vdev, INT *x, INT *y, INT idx, INT rad)
#endif
{
	Vsf_interior(vdev, BLACK_INTERIOR);
	Vsf_color(vdev, BLACK_COLOR);
	V_pie(vdev, x[idx], y[idx], rad, 0, 0);
	Vsf_interior(vdev, 0);
	V_arc(vdev, x[idx], y[idx], rad, 0 /* begang */, 3600 /* endang */);
	return OK;
}

#if TEXDOCU
INT draw_line(VDEVICE *vdev, INT *x, INT *y, INT idx1, INT idx2)
#endif
{
	Bind_mol(vdev, 1, x[idx1], y[idx1], x[idx2], y[idx2], 1);
	return OK;
}

#if TEXDOCU
INT draw_bezier3(VDEVICE *vdev, INT *x, INT *y, INT i1, INT i2, INT i3)
#endif
{
	INT pts[50];

	pts[0] = x[i1];
	pts[1] = y[i1];
	pts[2] = x[i2];
	pts[3] = y[i2];
	pts[4] = x[i3];
	pts[5] = y[i3];
	V_bezier(vdev, 3, pts);
	return OK;
}

#if TEXDOCU
INT draw_bezier4(VDEVICE *vdev, INT *x, INT *y, INT i1, INT i2, INT i3, INT i4)
#endif
{
	INT pts[50];

	pts[0] = x[i1];
	pts[1] = y[i1];
	pts[2] = x[i2];
	pts[3] = y[i2];
	pts[4] = x[i3];
	pts[5] = y[i3];
	pts[6] = x[i4];
	pts[7] = y[i4];
	V_bezier(vdev, 4, pts);
	return OK;
}

#if TEXDOCU
INT draw_bezier5(VDEVICE *vdev, INT *x, INT *y, INT i1, INT i2, INT i3, INT i4, INT i5)
#endif
{
	INT pts[50];

	pts[0] = x[i1];
	pts[1] = y[i1];
	pts[2] = x[i2];
	pts[3] = y[i2];
	pts[4] = x[i3];
	pts[5] = y[i3];
	pts[6] = x[i4];
	pts[7] = y[i4];
	pts[8] = x[i5];
	pts[9] = y[i5];
	V_bezier(vdev, 5, pts);
	return OK;
}

#if TEXDOCU
INT draw_bezier6(VDEVICE *vdev, INT *x, INT *y, 
	INT i1, INT i2, INT i3, INT i4, INT i5, INT i6)
#endif
{
	INT pts[50];

	pts[0] = x[i1];
	pts[1] = y[i1];
	pts[2] = x[i2];
	pts[3] = y[i2];
	pts[4] = x[i3];
	pts[5] = y[i3];
	pts[6] = x[i4];
	pts[7] = y[i4];
	pts[8] = x[i5];
	pts[9] = y[i5];
	pts[10] = x[i6];
	pts[11] = y[i6];
	V_bezier(vdev, 6, pts);
	return OK;
}

#if TEXDOCU
INT draw_bezier7(VDEVICE *vdev, INT *x, INT *y, 
	INT i1, INT i2, INT i3, INT i4, INT i5, INT i6, INT i7)
#endif
{
	INT pts[50];

	pts[0] = x[i1];
	pts[1] = y[i1];
	pts[2] = x[i2];
	pts[3] = y[i2];
	pts[4] = x[i3];
	pts[5] = y[i3];
	pts[6] = x[i4];
	pts[7] = y[i4];
	pts[8] = x[i5];
	pts[9] = y[i5];
	pts[10] = x[i6];
	pts[11] = y[i6];
	pts[12] = x[i7];
	pts[13] = y[i7];
	V_bezier(vdev, 7, pts);
	return OK;
}

#if TEXDOCU
INT halve_way_int(INT *Px, INT *Py, INT i, INT j, INT k)
#endif
{
	Px[k] = (Px[i] + Px[j]) >> 1;
	Py[k] = (Py[i] + Py[j]) >> 1;
	return OK;
}

#if TEXDOCU
INT rotate_int(INT *Px, INT *Py, INT i, INT j, INT angle)
#endif
{
	double c, s, ms;
	double a;
	INT x, y;

	a = ((double) angle * PI) / 180.;
	c = cos(a);
	s = sin(a);
	ms = -1. * sin(a);
	x = (INT)(c * (double) Px[i] + s * (double) Py[i]);
	y = (INT)(ms * (double) Px[i] + c * (double) Py[i]);
	Px[j] = x;
	Py[j] = y;
	return OK;
}

#if TEXDOCU
INT on_circle_int(INT *Px, INT *Py, INT idx, INT angle_in_degree, INT rad)
#endif
{
	
	Px[idx] = (INT)(cos_grad(angle_in_degree) * (double) rad);
	Py[idx] = (INT)(sin_grad(angle_in_degree) * (double) rad);
	return OK;
}

#if TEXDOCU
INT On_circle_int(VECTOR_OP Px, VECTOR_OP Py, INT idx, INT angle_in_degree, INT rad)
#endif
{
	INT x, y;
	
	x = (INT)(cos_grad(angle_in_degree) * (double) rad);
	y = (INT)(sin_grad(angle_in_degree) * (double) rad);
	Px->m_ii(idx, x);
	Py->m_ii(idx, y);
	return OK;
}

#ifndef EPSILON
#define EPSILON 0.000001
#endif

#if TEXDOCU
INT intersect_lines(INT *Px, INT *Py, 
	INT i1, INT i2, INT j1, INT j2, INT res)
#endif
{
	intersect_lines_int(Px[i1], Py[i1], Px[i2], Py[i2], 
		Px[j1], Py[j1], Px[j2], Py[j2], &Px[res], &Py[res]);
	return OK;
}

#if TEXDOCU
INT intersect_lines_int(INT Px, INT Py, INT Qx, INT Qy, 
	INT Rx, INT Ry, INT Sx, INT Sy, INT *Ix, INT *Iy)
#endif
{
	INT u1, u2, v1, v2, b1, b2;
	double det, lambda;
	
	u1 = Qx - Px;
	u2 = Qy - Py;
	v1 = -(Sx - Rx);
	v2 = -(Sy - Ry);
	b1 = Rx - Px;
	b2 = Ry - Py;
	det = u1 * v2 - u2 * v1;
	if (ABS(det) < EPSILON)
		return error("intersect_lines det = 0");
	lambda = v2 * b1 - v1 * b2;
	// mu = -u2 * b1 + u1 * b2;
	lambda = lambda / det;
	// mu = mu / det;
	*Ix = Px + (INT)(lambda * u1);
	*Iy = Py + (INT)(lambda * u2);
	return OK;
}

INT intersect_lines_P_u_Q_v(INT Px, INT Py, INT ux, INT uy, 
	INT Rx, INT Ry, INT vx, INT vy, INT *Ix, INT *Iy)
{
	double lambda;
	
	if (!intersect_lines_P_u_Q_v_factor(Px, Py, ux, uy, Rx, Ry, vx, vy, &lambda))
		return FALSE;
	*Ix = Px + (INT)(lambda * (double) ux);
	*Iy = Py + (INT)(lambda * (double) uy);
	return TRUE;
}

INT intersect_lines_P_u_Q_v_factor(INT Px, INT Py, INT ux, INT uy, 
	INT Rx, INT Ry, INT vx, INT vy, double *f)
{
	INT v1, v2, b1, b2, nu, nv;
	double det, lambda;
	
	nu = (INT) sqrt((double) ux * (double) ux + (double) uy * (double) uy);
	if (nu < EPSILON) {
		printf("intersect_lines_P_u_Q_v(): nu < EPSILON\n");
		printf("ux = %ld uy=%ld nu=%ld\n", ux, uy, nu);
		return FALSE;
		}
	nv = (INT) sqrt((double) vx * (double) vx + (double) vy * (double) vy);
	nv = vx * vx + vy * vy;
	if (nv < EPSILON) {
		printf("intersect_lines_P_u_Q_v(): nv < EPSILON\n");
		return FALSE;
		}
	
	v1 = -vx;
	v2 = -vy;
	b1 = Rx - Px;
	b2 = Ry - Py;
	det = ux * v2 - uy * v1;
	if (ABS(det) < EPSILON) {
		printf("intersect_lines_P_u_Q_v det = 0, cannot determine point of intersection\n");
		return FALSE;
		}
	lambda = v2 * b1 - v1 * b2;
	// mu = -u2 * b1 + u1 * b2;
	lambda = lambda / det;
	// mu = mu / det;
	*f = lambda;
	return TRUE;
}

INT intersect_line_area(VDEVICE *vdev, INT *x, INT *y, INT k, INT *points, 
	INT Px, INT Py, INT dx, INT dy, INT f_rim, INT rim_length)
{
	INT x0, y0, x1, y1;
	double lambda;
	INT *intersect_idx = NULL;
	double *intersect_lambda = NULL;
	INT i, nb_i = 0;
	
	intersect_idx = (INT *) my_malloc(sizeof(INT) * k, "intersect_idx");
	intersect_lambda = (double *) my_malloc(sizeof(double) * k, "intersect_idx");
	for (i = 0; i < k; i++) {
		if (!intersect_lines_P_u_Q_v_factor(
			x[points[i]], y[points[i]], 
			x[points[(i + 1) % k]] - x[points[i]], 
			y[points[(i + 1) % k]] - y[points[i]], 
			Px, Py, dx, dy, &lambda)) {
			printf("intersect_line_area() i=%ld, problem with points "
				"p[i]=%ld p[i+1]=%ld\n", i, points[i], points[(i + 1) % k]);
			continue;
			}
		if (lambda >= 0. && lambda <= 1.) {
			intersect_idx[nb_i] = i;
			intersect_lambda[nb_i] = lambda;
			nb_i++;
			}
		}
	if (nb_i != 2) {
		printf("intersect_line_area() no two points of intersection, nb_i = %ld\n", nb_i);
		return OK;
		}
	x0 = x[points[intersect_idx[0]]] + (INT)(intersect_lambda[0] * 
		(double)(x[points[(intersect_idx[0] + 1) % k]] - x[points[intersect_idx[0]]]));
	y0 = y[points[intersect_idx[0]]] + (INT)(intersect_lambda[0] * 
		(double)(y[points[(intersect_idx[0] + 1) % k]] - y[points[intersect_idx[0]]]));
	x1 = x[points[intersect_idx[1]]] + (INT)(intersect_lambda[1] * 
		(double)(x[points[(intersect_idx[1] + 1) % k]] - x[points[intersect_idx[1]]]));
	y1 = y[points[intersect_idx[1]]] + (INT)(intersect_lambda[1] * 
		(double)(y[points[(intersect_idx[1] + 1) % k]] - y[points[intersect_idx[1]]]));

	if (f_rim) {
		line_with_whole(vdev, x0, y0, x1, y1, rim_length);
		}
	else {
		my_line(vdev, x0, y0, x1, y1);
		}
	return OK;
}

INT line_with_whole(VDEVICE *vdev, INT x0, INT y0, INT x1, INT y1, INT rim_length)
{
	double dx, dy, dx0, dy0, l;
	INT xx, yy;
	
	dx = (double) (x1 - x0);
	dy = (double) (y1 - y0);
	l = sqrt(dx * dx + dy * dy);
	if ((INT) l < rim_length) {
		my_line(vdev, x0, y0, x1, y1);
		return OK;
		}
	dx0 = dx / l;
	dy0 = dy / l;
	xx = x0 + (INT)(dx0 * (double) rim_length);
	yy = y0 + (INT)(dy0 * (double) rim_length);
	my_line(vdev, x0, y0, xx, yy);
	xx = x1 - (INT)(dx0 * (double) rim_length);
	yy = y1 - (INT)(dy0 * (double) rim_length);
	my_line(vdev, xx, yy, x1, y1);
	return OK;
}

#if TEXDOCU
INT lotfuss_int(INT *x, INT *y, INT i, INT a, INT b, INT *lot_x, INT *lot_y)
#endif
{
	double ux, uy, vx, vy;
	INT Px, Py, Sx, Sy;
	
	ux = (double) (x[b] - x[a]);
	uy = (double) (y[b] - y[a]);
	if (ABS(ux) > EPSILON) {
		vx = - uy / ux;
		vy = 1.;
		}
	else {
		if (ABS(uy) < EPSILON) {
			return error("ux < EPSILON and uy < EPSILON");
			}
		vy = - ux / uy;
		vx = 1.;
		}
	Px = x[i] + (INT) (1000. * vx);
	Py = y[i] + (INT) (1000. * vy);
	intersect_lines_int(Px, Py, x[i], y[i], x[a], y[a], x[b], y[b], &Sx, &Sy);
	*lot_x = Sx;
	*lot_y = Sy;
	return OK;
}

#if TEXDOCU
INT ratio_int(INT *x, INT *y, INT a, INT b, double r, INT *xx, INT *yy)
#endif
{
	INT dx, dy;

	dx = (INT)((double)(x[b] - x[a]) * r);
	dy = (INT)((double)(y[b] - y[a]) * r);
	*xx = x[a] + dx;
	*yy = y[a] + dy;
	return OK;
}

INT Ratio(INT *Px, INT *Py, INT *Pz, INT from, INT to, double r, INT idx)
{
	INT dx, dy, dz;

	dx = (INT)((double)(Px[to] - Px[from]) * r);
	dy = (INT)((double)(Py[to] - Py[from]) * r);
	dz = (INT)((double)(Pz[to] - Pz[from]) * r);
	Px[idx] = Px[from] + dx;
	Py[idx] = Py[from] + dy;
	Pz[idx] = Pz[from] + dz;
	return OK;
}

INT center3(INT *Px, INT *Py, INT *Pz, INT i1, INT i2, INT i3, INT idx)
{
	INT x, y, z;
	
	x = Px[i1] + Px[i2] + Px[i3];
	y = Py[i1] + Py[i2] + Py[i3];
	z = Pz[i1] + Pz[i2] + Pz[i3];
	Px[idx] = (INT)((double) x * 0.3333333);
	Py[idx] = (INT)((double) y * 0.3333333);
	Pz[idx] = (INT)((double) z * 0.3333333);
	return OK;
}

INT center4(INT *Px, INT *Py, INT *Pz, INT i1, INT i2, INT i3, INT i4, INT idx)
{
	INT x, y, z;
	
	x = Px[i1] + Px[i2] + Px[i3] + Px[i4];
	y = Py[i1] + Py[i2] + Py[i3] + Py[i4];
	z = Pz[i1] + Pz[i2] + Pz[i3] + Pz[i4];
	Px[idx] = (INT)((double) x * 0.25);
	Py[idx] = (INT)((double) y * 0.25);
	Pz[idx] = (INT)((double) z * 0.25);
	return OK;
}

INT stuetzgeraden_x(INT *x, INT *y, INT k, INT *points, INT phi, INT y0, 
	INT *p_xmin, INT *p_xmax)
{
	INT i, x00, y00, x_min = 0, x_max = 0, f_first = TRUE;
	INT dx, dy;
	
	dx = (INT)(100. * cos_grad(phi * 0.1));
	dy = (INT)(100. * sin_grad(phi * 0.1));
	for (i = 0; i < k; i++) {
		if (intersect_lines_P_u_Q_v(x[points[i]], y[points[i]], dx, dy, 
			0, y0, 100, 0, 
			&x00, &y00)) {
			if (f_first) {
				x_min = x00;
				x_max = x00;
				f_first = FALSE;
				}
			else {
				x_min = MIN(x_min, x00);
				x_max = MAX(x_max, x00);
				}
			}
		else {
			printf("stuetzgeraden_x() problem in "
				"intersect_lines_P_u_Q_v for i=%ld\n", i);
			}
		}
	*p_xmin = x_min;
	*p_xmax = x_max;
	return OK;
}

INT stuetzgeraden_y(INT *x, INT *y, INT k, INT *points, INT phi, INT x0, 
	INT *p_ymin, INT *p_ymax)
{
	INT i, x00, y00, y_min = 0, y_max = 0, f_first = TRUE;
	INT dx, dy;
	
	dx = (INT)(100. * cos_grad(phi * 0.1));
	dy = (INT)(100. * sin_grad(phi * 0.1));	
	for (i = 0; i < k; i++) {
		if (intersect_lines_P_u_Q_v(x[points[i]], y[points[i]], dx, dy, 
			x0, 0, 0, 100, 
			&x00, &y00)) {
			if (f_first) {
				y_min = y00;
				y_max = y00;
				f_first = FALSE;
				}
			else {
				y_min = MIN(y_min, y00);
				y_max = MAX(y_max, y00);
				}
			}
		else {
			printf("stuetzgeraden_y() problem in "
				"intersect_lines_P_u_Q_v for i=%ld\n", i);
			}
		}
	*p_ymin = y_min;
	*p_ymax = y_max;
	return OK;
}

INT build_path(VDEVICE *vdev, INT *x, INT *y, INT k, INT *points)
{
	INT i;
	
	for (i = 0; i < k; i++) {
		my_line(vdev, 
			x[points[i]], y[points[i]], 
			x[points[(i + 1) % k]], 
			y[points[(i + 1) % k]]);
		}
	return OK;
}
 
INT shade_area(VDEVICE *vdev, INT *x, INT *y, INT k, INT *points, INT phi, INT delta, 
	INT f_rim, INT rim_length)
// 0 \le phi \le 3600
{
	INT dia = 0, f_x, dx, dy, xmin = 0, xmax = 0, ymin = 0, ymax = 0, x0, x1, y0, y1;
	INT xc, yc, i, j;
	
	if (delta <= 0) {
		printf("shade_area() waring: delta <= 0, setting delta to 3");
		delta = 3;
		}
	if (k < 2) {
		printf("shade_area() no area");
		return OK;
		}
	while (phi < 0)
		phi += 3600;
	phi = phi % 3600;
	
	if ((phi > 450 && phi < 1350) || (phi > 2250 && phi < 3150)) {
		f_x = TRUE;
		}
	else {
		f_x = FALSE;
		}
	printf("shade_area() f_x = %ld\n", f_x);
	for (i = 0; i < k; i++) {
		for (j = i + 1; j < k; j++) {
			dx = ABS(x[points[j]] - x[points[i]]);
			dy = ABS(y[points[j]] - y[points[i]]);
			dia = MAX(dia, dx);
			dia = MAX(dia, dy);
			}
		if (i == 0) {
			xmin = xmax = x[points[i]];
			ymin = ymax = y[points[i]];
			}
		else {
			xmin = MIN(xmin, x[points[i]]);
			xmax = MAX(xmax, x[points[i]]);
			ymin = MIN(ymin, y[points[i]]);
			ymax = MAX(ymax, y[points[i]]);
			}
		}
	printf("diameter=%ld\n", dia);
	printf("xmin=%ld xmax=%ld ymin=%ld ymax=%ld\n", xmin, xmax, ymin, ymax);
	dx = (INT)(100. * cos_grad(phi * 0.1));
	dy = (INT)(100. * sin_grad(phi * 0.1));
	if (f_x) {
		y0 = ymin;
		stuetzgeraden_x(x, y, k, points, phi, y0, &x0, &x1);
		yc = y0;
		xc = x0;
		while (TRUE) {
			intersect_line_area(vdev, x, y, k, points, 
				xc, yc, dx, dy, f_rim, rim_length);
			
			xc += delta;
			if (xc > x1)
				break;
			}
		}
	else {
		x0 = xmin;
		stuetzgeraden_y(x, y, k, points, phi, x0, &y0, &y1);
		yc = y0;
		xc = x0;
		while (TRUE) {
			intersect_line_area(vdev, x, y, k, points, 
				xc, yc, dx, dy, f_rim, rim_length);
			
			yc += delta;
			if (yc > y1)
				break;
			}
		}
	return OK;
}

INT translate(INT *Px, INT *Py, INT *Pz, INT dx, INT dy, INT dz, INT n)
{
	INT i;
	
	for (i = 0; i < n; i++) {
		Px[i] += dx;
		Py[i] += dy;
		Pz[i] += dz;
		}
	return OK;
}

INT rotate(INT *Px, INT *Py, INT *Pz, INT n)
{
	double A[9];
	double x[3], y[3];
	INT i;
	
	random_isometry(A);
	for (i = 0; i < n; i++) {
		x[0] = (double) Px[i];
		x[1] = (double) Py[i];
		x[2] = (double) Pz[i];
		mult_mat_vec(A, x, y);
		Px[i] = (INT) y[0];
		Py[i] = (INT) y[1];
		Pz[i] = (INT) y[2];
		}
	return OK;
}

INT random_isometry(double *A)
{
	double nx, ny, nz, nv;
	INT phi;
	
	nx = (double) rand();
	ny = (double) rand();
	nz = (double) rand();
	nv = 1. / sqrt(nx * nx + ny * ny + nz * nz);
	nx *= nv;
	ny *= nv;
	nz *= nv;
	phi = (INT)((double)rand() * 360. / (double)RAND_MAX);
	isometry(A, nx, ny, nz, phi);
	return OK;
}

INT isometry(double *A, double nx, double ny, double nz, INT phi)
{
	double ax, ay, az;
	double bx, by, bz;
	double B[9], Bv[9], Phi[9], C[9], At[9], D[9];
	double c, s, d;
	
	orthogonal(nx, ny, nz, &ax, &ay, &az, &bx, &by, &bz);
	B[0 * 3 + 0] = nx;
	B[1 * 3 + 0] = ny;
	B[2 * 3 + 0] = nz;
	B[0 * 3 + 1] = ax;
	B[1 * 3 + 1] = ay;
	B[2 * 3 + 1] = az;
	B[0 * 3 + 2] = bx;
	B[1 * 3 + 2] = by;
	B[2 * 3 + 2] = bz;
	inverse_mat(B, Bv);
	c = cos_grad(phi);
	s = sin_grad(phi);
	Phi[0 * 3 + 0] = 1.;
	Phi[1 * 3 + 0] = 0.;
	Phi[2 * 3 + 0] = 0.;
	Phi[0 * 3 + 1] = 0.;
	Phi[1 * 3 + 1] = c;
	Phi[2 * 3 + 1] = -s;
	Phi[0 * 3 + 2] = 0.;
	Phi[1 * 3 + 2] = s;
	Phi[2 * 3 + 2] = c;
	
	transpose_mat(Phi, At);
	mult_mat(Phi, At, D);
	printf("Phi * Phi^t:\n");
	print_mat33(D);

	// A = B * Phi * Bv
	mult_mat(B, Phi, C);
	mult_mat(C, Bv, A);
	
	printf("isometry:\n");
	print_mat33(A);
	d = det3(A);
	printf("det = %lf\n", d);
	transpose_mat(A, At);
	mult_mat(A, At, D);
	printf("A * A^t:\n");
	print_mat33(D);
	return OK;
}

INT orthogonal(double nx, double ny, double nz, 
	double *ax, double *ay, double *az, 
	double *bx, double *by, double *bz)
{
	double u, v, nuv, nvv, sp, nbv;
	double a[3], b[3];
	
	if (ABS(nx) < EPSILON) {
		if (ABS(ny) < EPSILON) {
			orthogonal(nz, nx, ny, az, ax, ay, bz, bx, by);
			}
		else {
			orthogonal(ny, nx, nz, ay, ax, az, by, bx, bz);
			}
		}
	u = ny / nx;
	v = nz / nx;
	nuv = 1. / sqrt(1. + u * u);
	nvv = 1. / sqrt(1. + v * v);
	a[0] = u * nuv;
	a[1] = -1. * nuv;
	a[2] = 0.;
	b[0] = v * nvv;
	b[1] = 0.;
	b[2] = -1. * nvv;
	sp = a[0] * b[0];
	b[0] = b[0] - sp * a[0];
	b[1] = b[1] - sp * a[1];
	b[2] = b[2] - sp * a[2];
	nbv = 1. / sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
	b[0] *= nbv;
	b[1] *= nbv;
	b[2] *= nbv;
	*ax = a[0];
	*ay = a[1];
	*az = a[2];
	*bx = b[0];
	*by = b[1];
	*bz = b[2];
	return OK;
}

double det3(double *A)
// sarrus
{
	double d, a;
	
	d = A[0 * 3 + 0] * A[1 * 3 + 1] * A[2 * 3 + 2];
	a = A[1 * 3 + 0] * A[2 * 3 + 1] * A[0 * 3 + 2]; d += a;
	a = A[2 * 3 + 0] * A[0 * 3 + 1] * A[1 * 3 + 2]; d += a;
	a = A[2 * 3 + 0] * A[1 * 3 + 1] * A[0 * 3 + 2]; d -= a;
	a = A[1 * 3 + 0] * A[0 * 3 + 1] * A[2 * 3 + 2]; d -= a;
	a = A[0 * 3 + 0] * A[2 * 3 + 1] * A[1 * 3 + 2]; d -= a;
	return d;
}

double det4(double a11, double a12, double a21, double a22)
{
	return (a11 * a22 - a21 * a12);
}

double detAij(double *A, INT i, INT j)
{
	INT i1 = 0, i2 = 0, j1 = 0, j2 = 0;
	
	if (i == 0) {
		i1 = 1;
		i2 = 2;
		}
	else if (i == 1) {
		i1 = 0;
		i2 = 2;
		}
	else if (i == 2) {
		i1 = 0;
		i2 = 1;
		}
	if (j == 0) {
		j1 = 1;
		j2 = 2;
		}
	else if (j == 1) {
		j1 = 0;
		j2 = 2;
		}
	else if (j == 2) {
		j1 = 0;
		j2 = 1;
		}
	return det4(A[i1 * 3 + j1], A[i1 * 3 + j2], A[i2 * 3 + j1], A[i2 * 3 + j2]);
}

INT transpose_mat(double *X, double *Y)
{
	Y[0 * 3 + 0] = X[0 * 3 + 0];
	Y[0 * 3 + 1] = X[1 * 3 + 0];
	Y[0 * 3 + 2] = X[2 * 3 + 0];
	Y[1 * 3 + 0] = X[0 * 3 + 1];
	Y[1 * 3 + 1] = X[1 * 3 + 1];
	Y[1 * 3 + 2] = X[2 * 3 + 1];
	Y[2 * 3 + 0] = X[0 * 3 + 2];
	Y[2 * 3 + 1] = X[1 * 3 + 2];
	Y[2 * 3 + 2] = X[2 * 3 + 2];
	return OK;
}

INT inverse_mat(double *X, double *Y)
{
	double detX, dxv, s;
	INT i, j;
	
	detX = det3(X);
	dxv = 1. / detX;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			if (EVEN(i + j))
				s = 1.;
			else
				s = -1.;
			Y[i * 3 + j] = s * dxv * detAij(X, j, i);
			}
		}
	// mult_mat(Y, Yinv, A);
	return OK;
}

INT mult_mat(double *A, double *B, double *C)
{
	double a, b;
	INT i, j, k;
	
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			a = 0;
			for (k = 0; k < 3; k++) {
				b = A[i * 3 + k] * B[k * 3 + j];
				a += b;
				}
			C[i * 3 + j] = a;
			}
		}
	return OK;
}

INT mult_mat_vec(double *A, double *v, double *v1)
// v1 := A * v
{
	double a, b;
	INT i, k;
	
	for (i = 0; i < 3; i++) {
		a = 0;
		for (k = 0; k < 3; k++) {
			b = A[i * 3 + k] * v[k];
			a += b;
			}
		v1[i] = a;
		}
	return OK;
}

INT print_mat33(double *A)
{
	INT i, j;
	
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			printf("%5.3lf ", A[i * 3 + j]);
			}
		printf("\n");
		}
	return OK;
}

INT print_vec3(double *v)
{
	INT i;
	
	for (i = 0; i < 3; i++) {
		printf("%2.1lf ", v[i]);
		printf("\n");
		}
	return OK;
}



#endif /* GRAPHICS_TRUE */

