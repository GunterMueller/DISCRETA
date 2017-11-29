// dial_solid.C

Widget widget_solid(INT i)
{
	Widget w;
	
	switch (i) {
		case 0:
			w = To_tetra;
			break;
		case 1:
			w = To_cube;
			break;
		case 2:
			w = To_octa;
			break;
		case 3:
			w = To_dode;
			break;
		case 4:
			w = To_ico;
			break;
		case 5:
			w = To_cube4d;
			break;
		case 6:
			w = To_cubus_simus;
			break;
		case 7:
			w = To_dode_simum;
			break;
		case 8:
			w = To_cube_ee;
			break;
		case 9:
			w = To_cube_ee_russian;
			break;
		case 10:
			w = To_dual;
			break;
		case 11:
			w = To_trunc;
			break;
		case 12:
			w = To_truncdode;
			break;
		case 13:
			w = To_trunccube;
			break;
		case 14:
			w = To_midpoints_of_edges;
			break;
		case 15:
			w = To_central_point;
			break;
		case 16:
			w = To_central_involution;
			break;
		case 17:
			w = To_group_on_edges;
			break;
		case 18:
			w = To_relabel_points;
			break;
		case 19:
			w = To_ext1;
			break;
		case 20:
			w = To_ext2;
			break;
		case 21:
			w = To_ext3;
			break;
		default:
			printf("widget_solid() unknown toggle state !\n");
			break;
		}
	return w;
}

void choose_regular_solid()
{
	DISCRETA_FORM_DATA *p = form_data;
	INT type, val1, val2;
	BYTE *q;
	
	get_form_data(p);
	switch (p->solid_toggle) {
		case 0:
			type = FGA_TETRAHEDRON;
			break;
		case 1:
			type = FGA_CUBE;
			break;
		case 2:
			type = FGA_OCTAHEDRON;
			break;
		case 3:
			type = FGA_DODECAHEDRON;
			break;
		case 4:
			type = FGA_ICOSAHEDRON;
			break;
		case 5:
			type = FGA_CUBE_4D;
			break;
		case 6:
			type = FGA_SOLID_CUBUS_SIMUS;
			break;
		case 7:
			type = FGA_SOLID_DODE_SIMUM;
			break;
		case 8:
			type = FGA_SOLID_CUBUS_EE;
			break;
		case 9:
			type = FGA_SOLID_CUBUS_EE_RUSSIAN;
			break;
		case 10:
			type = FGA_SOLID_DUAL;
			break;
		case 11:
			type = FGA_SOLID_TRUNCATE;
			break;
		case 12:
			type = FGA_SOLID_TRUNCATE_DODE;
			break;
		case 13:
			type = FGA_SOLID_TRUNCATE_CUBE;
			break;
		case 14:
			type = FGA_SOLID_EDGE_MIDPOINTS;
			break;
		case 15:
			type = FGA_SOLID_ADD_CENTRAL_POINT;
			break;
		case 16:
			type = FGA_SOLID_ADD_CENTRAL_INVOLUTION;
			break;
		case 17:
			type = FGA_SOLID_INDUCED_GROUP_ON_EDGES;
			break;
		case 18:
			type = FGA_SOLID_RELABEL_POINTS;
			break;
		case 19:
			type = FGA_TETRAHEDRON; // unused 
			break;
		case 20:
			type = FGA_TETRAHEDRON; // unused 
			break;
		case 21:
			type = FGA_TETRAHEDRON; // unused 
			break;
		default:
			printf("choose_regular_solid() unknown toggle state %d!\n", 
				p->solid_toggle);
			break;
		}
	val1 = 0;
	val2 = 0;

	group_selection_add(p, type, val1, val2, NIL);
}

