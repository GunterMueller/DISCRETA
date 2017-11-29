// dial_sporadic.C

Widget widget_sporadic(INT i)
{
	Widget w;
	
	switch (i) {
		case 0:
			w = To_M11;
			break;
		case 1:
			w = To_M12;
			break;
		case 2:
			w = To_M23;
			break;
		case 3:
			w = To_M24;
			break;
		case 4:
			w = To_HS;
			break;
		default:
			printf("widget_sporadic() unknown toggle state !\n");
			break;
		}
	return w;
}

void choose_sporadic_group()
{
	DISCRETA_FORM_DATA *p = form_data;
	INT type, val1 = 0, val2 = 0;
	BYTE *q;
	
	get_form_data(p);
	switch (p->sporadic_toggle) {
		case 0:
			type = FGA_GROUP_MATHIEU;
			val1 = 11;
			break;
		case 1:
			type = FGA_GROUP_MATHIEU;
			val1 = 12;
			break;
		case 2:
			type = FGA_GROUP_MATHIEU;
			val1 = 23;
			break;
		case 3:
			type = FGA_GROUP_MATHIEU;
			val1 = 24;
			break;
		case 4:
			type = FGA_GROUP_HIGMAN_SIMS_176;
			break;
		default:
			printf("choose_sporadic_group() unknown toggle state !\n");
			break;
		}
	group_selection_add(p, type, val1, val2, NIL);
}

