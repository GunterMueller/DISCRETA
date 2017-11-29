// dial_well.C

Widget widget_well_known(INT i)
{
	Widget w;
	
	switch (i) {
		case 0:
			w = To_Sn_short;
			break;
		case 1:
			w = To_Sn_trans;
			break;
		case 2:
			w = To_An;
			break;
		case 3:
			w = To_Dn;
			break;
		case 4:
			w = To_Cn;
			break;
		case 5:
			w = To_HolCn;
			break;
		case 6:
			w = To_Idn;
			break;
		case 7:
			w = To_subgroup_of_hol;
			break;
		case 8:
			w = To_SnwrSm;
			break;
		default:
			printf("unknown toggle state !\n");
			break;
		}
	return w;
}

void choose_well_known_group()
{
	DISCRETA_FORM_DATA *p = form_data;
	INT type, val1, val2;
	BYTE *q;
	
	get_form_data(p);
	switch (p->well_rb_state) {
		case 0:
			type = FGA_GROUP_SYM;
			break;
		case 1:
			type = FGA_GROUP_SYM_TRANSPOSITIONS;
			break;
		case 2:
			type = FGA_GROUP_ALT;
			break;
		case 3:
			type = FGA_GROUP_DIHEDRAL;
			break;
		case 4:
			type = FGA_GROUP_CYCLIC;
			break;
		case 5:
			type = FGA_GROUP_HOLOMORPH_OF_CYCLIC_GROUP;
			break;
		case 6:
			type = FGA_GROUP_TRIVIAL;
			break;
		case 7:
			type = FGA_SUBGROUP_OF_HOLOMORPH_OF_CYCLIC_GROUP;
			// type = FGA_GROUP_Zn_MULTIPLICATOR;
			break;
		case 8:
			type = FGA_GROUP_SN_WREATH_SM;
			break;
		default:
			printf("unknown toggle state !\n");
			break;
		}
	val1 = p->well_n;
	val2 = p->well_m;

	group_selection_add(p, type, val1, val2, NIL);
}

