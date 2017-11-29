// dial_linear.C

Widget widget_linear(INT i)
{
	Widget w;
	
	switch (i) {
		case 0:
			w = To_PSL;
			break;
		case 1:
			w = To_PGL;
			break;
		case 2:
			w = To_PSSL;
			break;
		case 3:
			w = To_PGGL;
			break;
		case 4:
			w = To_SL;
			break;
		case 5:
			w = To_GL;
			break;
		case 6:
			w = To_SSL;
			break;
		case 7:
			w = To_GGL;
			break;
		case 8:
			w = To_ASL;
			break;
		case 9:
			w = To_AGL;
			break;
		case 10:
			w = To_ASSL;
			break;
		case 11:
			w = To_AGGL;
			break;
		case 12:
			w = To_T;
			break;
		case 13:
			w = To_PSU_3_q2;
			break;
		case 14:
			w = To_Sz_q;
			break;
		default:
			printf("widget_linear() unknown toggle state !\n");
			break;
		}
	return w;
}

void choose_linear_group()
{
	DISCRETA_FORM_DATA *p = form_data;
	INT type, val1, val2;
	BYTE *q;
	
	get_form_data(p);
	switch (p->linear_toggle) {
		case 0:
			type = FGA_GROUP_PSL;
			break;
		case 1:
			type = FGA_GROUP_PGL;
			break;
		case 2:
			type = FGA_GROUP_PSSL;
			break;
		case 3:
			type = FGA_GROUP_PGGL;
			break;
		case 4:
			type = FGA_GROUP_SL;
			break;
		case 5:
			type = FGA_GROUP_GL;
			break;
		case 6:
			type = FGA_GROUP_SSL;
			break;
		case 7:
			type = FGA_GROUP_GGL;
			break;
		case 8:
			type = FGA_GROUP_ASL;
			break;
		case 9:
			type = FGA_GROUP_AGL;
			break;
		case 10:
			type = FGA_GROUP_ASSL;
			break;
		case 11:
			type = FGA_GROUP_AGGL;
			break;
		case 12:
			type = FGA_GROUP_AFFINE_TRANSLATIONS;
			break;
		case 13:
			type = FGA_GROUP_PSU_3_Q2;
			break;
		case 14:
			type = FGA_GROUP_SZ_Q;
			break;
		default:
			printf("choose_solvable_group() unknown toggle state !\n");
			break;
		}
	val1 = p->linear_n;
	val2 = p->linear_q;

	group_selection_add(p, type, val1, val2, NIL);
}

