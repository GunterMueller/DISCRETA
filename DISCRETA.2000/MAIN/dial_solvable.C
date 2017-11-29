// dial_solvable.C

Widget widget_solvable(INT i)
{
	Widget w;
	
	switch (i) {
		case 0:
			w = To_solvable_n_m;
			break;
		default:
			printf("widget_solvable() unknown toggle state !\n");
			break;
		}
	return w;
}

void choose_solvable_group()
{
	DISCRETA_FORM_DATA *p = form_data;
	INT type, val1, val2;
	BYTE *q;
	
	get_form_data(p);
	switch (p->solvable_toggle) {
		case 0:
			type = FGA_SOLVABLE_GROUP;
			break;
		default:
			printf("choose_solvable_group() unknown toggle state !\n");
			break;
		}
	val1 = p->solvable_n;
	val2 = p->solvable_m;

	group_selection_add(p, type, val1, val2, NIL);
}

