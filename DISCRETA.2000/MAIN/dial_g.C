// dial_g.C


#if 0
#include <stdlib.h>
#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>

#include <Xm/Xm.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>

#include "discreta.h"
#include "disc.h"
#endif


INT g_box_free(DISCRETA_FORM_DATA *p)
{
	INT i;
	
	for (i = 0; i < p->g_nb_items; i++) {
		XmStringFree(p->g_Items[i]);
		p->g_sel->s_i(i)->freeself();
		}
	return OK;
}

INT g_box_refresh(DISCRETA_FORM_DATA *p)
{
	Arg args[64];
	BYTE str[10000];
	INT i;
	GROUP_SELECTION_OP gs;
	
	for (i = 0; i < p->g_nb_items; i++) {
		
		str[0] = 0;
		gs = (GROUP_SELECTION_OP) p->g_sel->s_i(i);
		gs->sprint(str);
		p->g_Items[i] = (XmString) XmStringCreateLtoR(str, 
			(XmStringCharSet)XmFONTLIST_DEFAULT_TAG);
		}

	i = 0;
	XtSetArg(args[i], XmNitems, p->g_Items); i++;
	XtSetArg(args[i], XmNitemCount, p->g_nb_items); i++;
	XtSetValues(Gsel_list, args, i);
	return OK;
}

INT get_group_from_selection(VECTOR_OP gsel, INT nb_g_sel, 
	VECTOR_OP generators, BYTE *g_label, BYTE *g_label_tex)
{
	PERMUTATION_OP p;
	SYM_OB go;
	INT deg;
	BYTE *q;
	BYTE s_go[STRING_SIZE];
	INT i;

	km_get_group_from_selection(gsel, 
		nb_g_sel, generators, g_label, g_label_tex);

	if (generators->s_li() <= 0)
		return error("get_group_from_selection(): no generators");
	p = (PERMUTATION_OP) generators->s_i(0);
	i = p->s_li();
	set_v(i);
	set_group_label(g_label);
	
	vec_generators_group_order(generators, &go);
	s_go[0] = 0;
	go.sprint(s_go);
	set_group_order(s_go);
	
	{ 
		BYTE fname[1000];
		FILE *fp;
		
		sprintf(fname, "gen_%s.txt", g_label);
		fp = fopen(fname, "w");
		write_generators(generators, fp);
		fprintf(fp, "group: %s ($%s$), order=%s\n\n", g_label, g_label_tex, s_go);
		fclose(fp);
		sprintf(fname, "gen_%s.g", g_label);
		fp = fopen(fname, "w");
		generators->fprint_GAP(fp);
		fprintf(fp, "# group: %s ($%s$), order=%s\n\n", g_label, g_label_tex, s_go);
		fclose(fp);
	}
	
	return OK;
}

void group_selection_clear(DISCRETA_FORM_DATA *p)
{
	g_box_free(p);
	p->g_nb_items = 0;
	g_box_refresh(p);
}

void group_selection_delete(DISCRETA_FORM_DATA *p)
{
	Arg args[64];
	int *pos_list;
	int pos_count;
	INT i, j, k, n;

	if (!XmListGetSelectedPos(Gsel_list, &pos_list, &pos_count))
		return; /* no selected items */
	if (pos_count) {
		for (i = 0; i < pos_count; i++) {
			j = pos_list[i] - 1;
			XmStringFree(p->g_Items[j]);
			p->g_sel->s_i(j)->freeself();
			for (k = j + 1; k < p->g_nb_items; k++) {
				p->g_Items[k - 1] = p->g_Items[k];
				p->g_sel->s_i(k - 1)->swap(p->g_sel->s_i(k));
				}
			p->g_nb_items--;
			}
		XtFree((char *) pos_list);

		n = 0;
		XtSetArg(args[n], XmNitems, p->g_Items); n++;
		XtSetArg(args[n], XmNitemCount, p->g_nb_items); n++;
		XtSetValues(Gsel_list, args, n);
		}
}

void group_selection_add(DISCRETA_FORM_DATA *p, INT type, INT val1, INT val2, BYTE *s)
{
	Arg args[64];
	GROUP_SELECTION_OP gs;
	BYTE str[10000];
	INT i, n;
	
	if (p->g_nb_items >= MAX_ITEMS) {
		printf("group_selection_add() too many items");
		return;
		}
	i = p->g_nb_items;
	gs = (GROUP_SELECTION_OP) p->g_sel->s_i(i);
	gs->init(type, val1, val2, s);
	
	str[0] = 0;
	// group_selection_print(gs, str);
	gs->sprint(str);
	printf("added item %ld: %s\n", i, str);
	fflush(stdout);
	p->g_Items[i] = (XmString) XmStringCreateLtoR(str, 
		(XmStringCharSet)XmFONTLIST_DEFAULT_TAG);
	p->g_nb_items++;
	
	n = 0;
	XtSetArg(args[n], XmNitems, p->g_Items); n++;
	XtSetArg(args[n], XmNitemCount, p->g_nb_items); n++;
	XtSetValues(Gsel_list, args, n);
}

