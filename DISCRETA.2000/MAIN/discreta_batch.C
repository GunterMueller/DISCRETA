/* discreta_batch.C 
 * Anton Betten
 */

#include <stdio.h>
#include <stdlib.h>

#include <DISCRETA/discreta.h>
#include <DISCRETA/ladder.h>
#include <DISCRETA/geo.h> // for db_geo_by_base_blocks_add()

INT Argc;
BYTE **Argv;



void check_next_arg(INT i);
INT get_INT_arg(INT i);
BYTE *get_charp_arg(INT i);
INT group_label_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *g_label, BYTE *g_label_tex);
INT KM_fname_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, BYTE *KM_fname);
INT vtk_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, BYTE *KM_fname, 
	INT v, INT t, INT k);
INT KM_matrix_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT tt, INT kk);
INT stab_order_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable);
INT plesken_matrix_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT k_min, INT kk);
INT plesken_matrix_with_inverse_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT k_min, INT kk);
INT generators_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable);
INT orbit_representatives_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT f_with_stabilizer_orders);
INT stabilizer_orders_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable);
INT orbit_representatives_to_ASCII(BYTE *KM_fname, 
	BYTE *export_fname, INT k, INT f_with_stabilizer_orders);
INT nb_of_solutions_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT lambda);
INT read_until_lambdaend(FILE *fp);
INT solutions_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT lambda, INT from, INT len);
INT normalizing_permutation_on_orbits(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, BYTE *from_GAP, INT on_k_sets);
INT batch_fuse_orbits(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, BYTE *generators_fname, INT on_k_sets);
INT batch_fuse_orbits_by_representatives(
	BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *reps_fname, BYTE *generators_fname);
INT id_fname_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *db_prefix, BYTE *id_fname);
INT left_transversal_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *fname_G, BYTE *fname_U);


void check_next_arg(INT i)
{
	if (i >= Argc) {
		error("argument missing");
		exit(1);
		}
}

INT get_INT_arg(INT i)
{
	char *p;
	INT k;
	
	if (i >= Argc) {
		error("argument missing");
		exit(1);
		}
	p = Argv[i];
	sscanf(p, "%ld", &k);
	return k;
}

BYTE *get_charp_arg(INT i)
{
	char *p;
	
	if (i >= Argc) {
		error("argument missing");
		exit(1);
		}
	p = Argv[i];
	return p;
}

int main(int argc, char **argv)
// argc counts also the program name, so a call to 
// discreta_batch cmd
// gives argc = 2, argv[0] = discreta_batch, argv[1] = cmd
{
	INT t0, t1, user_time;
	BYTE s[256];

	Argc = argc;
	Argv = argv;
	if (argc < 2) {
		fprintf(stderr, "usage: %s command\n", argv[0]);
		exit(1);
		}
	
	discreta_init();
	{
	BYTE *command;
	INT i;
	BYTE *KM_fname;
	
	t0 = os_ticks();
	
	command = argv[1];
	printf("in discreta_batch: command = %s argc=%d\n", command, argc);
	for (i = 0; i < argc; i++) {
		printf("%ld: %s\n", i, argv[i]);
		}
	fflush(stdout);
	i = 2; // next argument
	if (strcmp(command, "calc_delta_lambda") == 0) {
		INT v, t, k;
		
		printf("calc_delta_lambda v t k\n"); fflush(stdout);
		v = get_INT_arg(i++);
		t = get_INT_arg(i++);
		k = get_INT_arg(i++);
		calc_delta_lambda(v, t, k);
		}
	
	else if (strcmp(command, "mendelsohn") == 0) {
		INT v, t, k, lambda, s_max;
		MATRIX_OB M;
		MATRIX_OB LAmbda;
		VECTOR_OB RHS;
		
		printf("mendelsohn v t k lambda s_max\n"); fflush(stdout);
		v = get_INT_arg(i++);
		t = get_INT_arg(i++);
		k = get_INT_arg(i++);
		lambda = get_INT_arg(i++);
		s_max = get_INT_arg(i++);
		calc_Lambda(v, t, k, lambda, &RHS);
		Mendelsohn(v, t, k, s_max, lambda, &M, &RHS, &LAmbda, TRUE /* f_v */);
		}
	
	
	else if (strcmp(command, "koehler") == 0) {
		INT v, t, k, lambda, m, s_max;
		VECTOR_OB RHS;
		VECTOR_OB Coeff, Constant_term;
		
		printf("koehler v t k lambda m s_max\n"); fflush(stdout);
		v = get_INT_arg(i++);
		t = get_INT_arg(i++);
		k = get_INT_arg(i++);
		lambda = get_INT_arg(i++);
		m = get_INT_arg(i++);
		s_max = get_INT_arg(i++);
		calc_Lambda(v, t, k, lambda, &RHS);
		// Koehler(v, t, m, s_max, &RHS, &Coeff, &Constant_term, TRUE /* f_v */);
		}
	
	else if (strcmp(command, "compose_group") == 0) {
		INT nb_args;
		BYTE **args;
		VECTOR_OB V, gen;
		BYTE fname[1024];
		BYTE g_label[1024];
		BYTE g_label_tex[1024];
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		FILE *fp;
		
		printf("compose_group GAP_fname GAP_variable args\n"); fflush(stdout);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		args = argv + i;
		nb_args = argc - i;
		compose_gsel_from_strings(&V, nb_args, args);
		km_get_group_from_selection(&V, 
			V.s_li(), &gen, g_label, g_label_tex);
		sprintf(fname, "%s.txt", g_label);
		printf("writing generators file %s\n", fname); fflush(stdout);
		printf("g_label_tex = %s\n", g_label_tex); fflush(stdout);
		write_file_of_generators(&gen, fname);
		group_label_to_GAP(GAP_fname, GAP_variable, g_label, g_label_tex);
		fp = fopen("group_label", "w");
		fprintf(fp, "%s\n", g_label);
		fprintf(fp, "%s\n", g_label_tex);
		fclose(fp);
		}
	
	else if (strcmp(command, "compute_KM") == 0) {
		INT t, k, f_TDO;
		BYTE *group_name;
		BYTE *group_name_tex;
		INT f_strong_generators = FALSE;
		INT f_orderly_generation = FALSE;
		INT f_extension_construction = FALSE;
		INT f_canonical_reps = TRUE;
		INT f_k_k2_design = FALSE, k2 = 0;
		VECTOR_OB gen;
		BYTE fname[1024];
		BYTE KM_fname[1024];
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		
		printf("compute_KM GAP_fname GAP_variable group_name group_name_tex t k f_TDO\n"); fflush(stdout);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		group_name = get_charp_arg(i++);
		group_name_tex = get_charp_arg(i++);
		t = get_INT_arg(i++);
		k = get_INT_arg(i++);
		f_TDO = get_INT_arg(i++);
		sprintf(fname, "%s.txt", group_name);
		sprintf(KM_fname, "KM_%s_t%ld_k%ld.txt", group_name, t, k);
		printf("trying to read generators from file %s\n", fname); fflush(stdout);
		read_file_of_generators(&gen, fname);
		compute_KM(&gen, group_name, group_name_tex, t, k, 
			f_strong_generators, f_orderly_generation, 
			f_TDO, f_extension_construction, 
			f_canonical_reps, f_k_k2_design, k2);
		KM_fname_to_GAP(GAP_fname, GAP_variable, KM_fname);
		}
	
	else if (strcmp(command, "get_vtk") == 0) {
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		BYTE *KM_fname;
		INT v, t, k;
		
		printf("get_vtk GAP_fname GAP_variable KM_fname\n"); fflush(stdout);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		KM_fname = get_charp_arg(i++);
		printf("KM_fname = %s\n", KM_fname); fflush(stdout);
		km_read_ascii_vtk(KM_fname, &v, &t, &k);
		vtk_to_GAP(GAP_fname, GAP_variable, KM_fname, v, t, k);
		}
	
	else if (strcmp(command, "show_KM_matrix") == 0) {
		BYTE *KM_fname;
		
		printf("show_KM_matrix KM_fname\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		printf("KM_fname = %s\n", KM_fname); fflush(stdout);
		show_km_matrix(KM_fname);
		}
	
	else if (strcmp(command, "do_LLL") == 0) {
		INT f_with, c0, beta, p, lambda;
		BYTE *KM_fname;
		
		printf("do_LLL f_with KM_fname c0 beta p lambda\n"); fflush(stdout);
		f_with = get_INT_arg(i++);
		KM_fname = get_charp_arg(i++);
		c0 = get_INT_arg(i++);
		beta = get_INT_arg(i++);
		p = get_INT_arg(i++);
		lambda = get_INT_arg(i++);
		do_LLL(!f_with, KM_fname, c0, beta, p, lambda, FALSE /* f_restart */, 0);
		}
	
	else if (strcmp(command, "do_McKay") == 0) {
		INT lambda;
		BYTE *KM_fname;
		
		printf("do_McKay KM_fname lambda\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		lambda = get_INT_arg(i++);
		do_mckay(KM_fname, lambda);
		}
	
	else if (strcmp(command, "get_solutions_from_solver") == 0) {
		INT lambda;
		BYTE *KM_fname;
		
		printf("get_solutions KM_fname lambda\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		lambda = get_INT_arg(i++);
		get_solutions(KM_fname, lambda);
		}
	
	else if (strcmp(command, "check_solutions") == 0) {
		INT lambda;
		BYTE *KM_fname;
		
		printf("check_solutions KM_fname lambda\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		lambda = get_INT_arg(i++);
		check_solutions(KM_fname, lambda);
		}
	
	else if (strcmp(command, "report") == 0) {
		INT f_select = FALSE;
		INT design_select_lambda = 0, select_first = 0, select_length = 0;
		
		printf("report KM_fname\n"); fflush(stdout);
#if 0
		if (strcmp(argv[i], "-select") == 0) {
			f_select = TRUE;
			sscanf(argv[i + 1], "%ld", &design_select_lambda);
			i++;
			sscanf(argv[i + 1], "%ld", &select_first);
			i++;
			sscanf(argv[i + 1], "%ld", &select_length);
			i++;
			printf("-select %ld %ld %ld\n", design_select_lambda, 
				select_first, select_length);
			}
#endif
		KM_fname = get_charp_arg(i++);
		printf("KM_fname = %s\n", KM_fname); fflush(stdout);

		do_report(KM_fname, TRUE /* f_html */, f_select, 
			design_select_lambda, select_first, select_length);
#if 0
		// the internal function:
		design_report(KM_fname, f_select, 
			design_select_lambda, select_first, select_length);
#endif
		}
	
	else if (strcmp(command, "report_plesken") == 0) {
		BYTE *KM_fname;
		
		printf("plesken KM_fname\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		printf("KM_fname = %s\n", KM_fname); fflush(stdout);
		design_report_plesken(KM_fname);
		}
	
	// retreive the data from DISCRETA and put it into GAP:

	else if (strcmp(command, "get_KM_matrix") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		INT t, k;
		
		printf("get_KM_matrix KM_fname GAP_fname GAP_variable t k\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		t = get_INT_arg(i++);
		k = get_INT_arg(i++);
		KM_matrix_to_GAP(KM_fname, GAP_fname, GAP_variable, t, k);
		}
	
#if 0
	else if (strcmp(command, "get_stab_order") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		
		printf("get_stab_order KM_fname GAP_fname GAP_variable\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		stab_order_to_GAP(KM_fname, GAP_fname, GAP_variable);
		}
#endif
	
	else if (strcmp(command, "get_plesken_matrix") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		INT t, k_min, k;
		
		printf("get_plesken_matrix KM_fname GAP_fname GAP_variable k_min k\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		k_min = get_INT_arg(i++);
		k = get_INT_arg(i++);
		plesken_matrix_to_GAP(KM_fname, GAP_fname, GAP_variable, k_min, k);
		}
	
	else if (strcmp(command, "get_plesken_matrix_with_inverse") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		INT t, k_min, k;
		
		printf("get_plesken_matrix KM_fname GAP_fname GAP_variable k\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		k_min = get_INT_arg(i++);
		k = get_INT_arg(i++);
		plesken_matrix_with_inverse_to_GAP(KM_fname, GAP_fname, GAP_variable, k_min, k);
		}
	
	else if (strcmp(command, "get_generators") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		INT t, k;
		
		printf("get_generators KM_fname GAP_fname GAP_variable\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		generators_to_GAP(KM_fname, GAP_fname, GAP_variable);
		}
	
	else if (strcmp(command, "get_orbit_representatives") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		
		printf("get_orbit_representatives KM_fname GAP_fname GAP_variable\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		orbit_representatives_to_GAP(KM_fname, GAP_fname, GAP_variable, 
			FALSE /* f_with_stabilizer_orders */);
		}
	
	else if (strcmp(command, "get_stabilizer_orders") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		
		printf("get_stabilizer_orders KM_fname GAP_fname GAP_variable\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		stabilizer_orders_to_GAP(KM_fname, GAP_fname, GAP_variable);
		}
	
	else if (strcmp(command, "get_orbit_representatives_ASCII") == 0) {
		BYTE *KM_fname;
		BYTE *export_fname;
		INT k;
		
		printf("get_orbit_representatives KM_fname export_fname k\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		export_fname = get_charp_arg(i++);
		k = get_INT_arg(i++);
		orbit_representatives_to_ASCII(KM_fname, export_fname, k, TRUE);
		}
	
	else if (strcmp(command, "get_number_of_solutions") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		INT lambda;
		
		printf("get_number_of_solutions KM_fname GAP_fname GAP_variable lambda\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		lambda = get_INT_arg(i++);
		nb_of_solutions_to_GAP(KM_fname, GAP_fname, GAP_variable, lambda);
		}
	
	else if (strcmp(command, "get_solutions") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		INT lambda, from, len;
		
		printf("get_solutions KM_fname GAP_fname GAP_variable lambda from len\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		lambda = get_INT_arg(i++);
		from = get_INT_arg(i++);
		len = get_INT_arg(i++);
		solutions_to_GAP(KM_fname, GAP_fname, GAP_variable, lambda, from, len);
		}
	
	else if (strcmp(command, "normalizing_permutation_on_orbits") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		BYTE *from_GAP;
		INT k;
		
		printf("normalizing_permutation_on_orbits KM_fname GAP_fname GAP_variable from_GAP k\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		from_GAP = get_charp_arg(i++);
		k = get_INT_arg(i++);
		normalizing_permutation_on_orbits(KM_fname, GAP_fname, GAP_variable, from_GAP, k);
		}
	
	else if (strcmp(command, "fuse_orbits") == 0) {
		BYTE *KM_fname;
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		BYTE *generators_fname;
		INT k;
		
		printf("fuse_orbits KM_fname GAP_fname GAP_variable generators_fname on_k_sets\n"); fflush(stdout);
		KM_fname = get_charp_arg(i++);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		generators_fname = get_charp_arg(i++);
		k = get_INT_arg(i++);
		batch_fuse_orbits(KM_fname, GAP_fname, GAP_variable, generators_fname, k);
		}
	
	else if (strcmp(command, "fuse_orbits_by_representatives") == 0) {
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		BYTE *reps_fname;
		BYTE *generators_fname;
		
		printf("fuse_orbits_by_representatives GAP_fname GAP_variable reps_fname generators_fname\n"); fflush(stdout);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		reps_fname = get_charp_arg(i++);
		generators_fname = get_charp_arg(i++);
		batch_fuse_orbits_by_representatives(GAP_fname, GAP_variable, reps_fname, generators_fname);
		}
	
	else if (strcmp(command, "geo_db_build_from_bb") == 0) {
		BYTE *generators_fname;
		BYTE *bb_fname;
		BYTE *db_prefix;
		INT f_create;
		INT f_0 = FALSE;
		INT f_tdo = TRUE;
		INT f_tdo2 = TRUE;
		INT f_tda = FALSE;
		INT f_range = FALSE;
		INT f = 0, l = 0;
		INT f_V = FALSE;
		INT f_B = FALSE;
		VECTOR_OB V, B;
		INT f_test2design = FALSE;
		INT f_do_not_add = FALSE;
		INT f_v = FALSE;
		
		printf("geo_db_build_from_bb generators_fname bb_fname db_prefix f_create\n"); fflush(stdout);
		generators_fname = get_charp_arg(i++);
		bb_fname = get_charp_arg(i++);
		db_prefix = get_charp_arg(i++);
		f_create = get_INT_arg(i++);
		
		db_geo_by_base_blocks_add(db_prefix, f_0, f_create, 
			generators_fname, bb_fname, 
			FALSE /* f_inc_file */, NIL /* inc_file_name */, 
			f_tdo, f_tdo2, 
			f_tda, f_range, f, l, f_test2design, 
			f_V, &V, f_B, &B, 
			f_do_not_add, 
			f_v, FALSE /* f_vv */);
		
		}
	
	else if (strcmp(command, "geo_db_export_id") == 0) {
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		BYTE *db_prefix;
		BYTE *id_fname;
		
		printf("geo_db_export_id GAP_fname GAP_variable db_prefix id_fname\n"); fflush(stdout);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		db_prefix = get_charp_arg(i++);
		id_fname = get_charp_arg(i++);

		db_geo_by_base_blocks_export_id(db_prefix, id_fname);
		id_fname_to_GAP(GAP_fname, GAP_variable, db_prefix, id_fname);
		}
	
	else if (strcmp(command, "left_transversal") == 0) {
		BYTE *GAP_fname;
		BYTE *GAP_variable;
		BYTE *fname_G;
		BYTE *fname_H;
		INT k;
		
		printf("left_transversal GAP_fname GAP_variable fname_G fname_H\n"); fflush(stdout);
		GAP_fname = get_charp_arg(i++);
		GAP_variable = get_charp_arg(i++);
		fname_G = get_charp_arg(i++);
		fname_H = get_charp_arg(i++);
		left_transversal_to_GAP(GAP_fname, GAP_variable, fname_G, fname_H);
		}
	
	else {
		printf("discreta_batch: unrecognized command\n");
		}
	
	
	
	printf("discreta_batch finished\n"); fflush(stdout);
	
	t1 = os_ticks();
	user_time = t1 - t0;
	s[0] = 0;
	print_delta_time(user_time, s);
	printf("total computing time: %s\n", s);
	fflush(stdout);
	}
	discreta_exit();
	return 0;
}

INT group_label_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *g_label, BYTE *g_label_tex)
{
	FILE *fp;

	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# group_label_to_GAP(%s %s %s %s);\n\n", 
		GAP_fname, GAP_variable, g_label, g_label_tex);
	fprintf(fp, "%s := [\"%s\", \"%s\"];\n", 
		GAP_variable, g_label, g_label_tex);
	fclose(fp);
	fflush(stdout);
	return OK;
}

INT KM_fname_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, BYTE *KM_fname)
{
	FILE *fp;

	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# KM_fname_to_GAP(%s %s %s);\n\n", 
		GAP_fname, GAP_variable, KM_fname);
	fprintf(fp, "%s := \"%s\";\n", GAP_variable, KM_fname);
	fclose(fp);
	fflush(stdout);
	return OK;
}

INT vtk_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, BYTE *KM_fname, 
	INT v, INT t, INT k)
{
	FILE *fp;

	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# get_vtk(%s %s %s);\n\n", 
		GAP_fname, GAP_variable, KM_fname);
	fprintf(fp, "%s := [%ld,%ld,%ld];\n", GAP_variable, v, t, k);
	fclose(fp);
	fflush(stdout);
	return OK;
}

INT KM_matrix_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT tt, INT kk)
{
	MATRIX_OB Mtk, Mttkk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	dc_Mtk_via_MM(tt, kk, &Mttkk, &MM, FALSE /* f_v */);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# KM_matrix_to_GAP(%s %s %s %ld %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, tt, kk);
	fprintf(fp, "%s := ", GAP_variable);
	Mttkk.fprint_GAP(fp);
	fprintf(fp, ";\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

#if 0
INT stab_order_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go;
	MATRIX_OB I, Ik2;
	VECTOR_OB K_first, K_len;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	INT i, f, l;
	MATRIX_OP pM;
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &I, &Ik2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	K_first.m_il(k + 1);
	K_len.m_il(k + 1);
	f = 0;
	for (i = 0; i < k; i++) {
		pM = (MATRIX_OP) MM.s_i(i);
		l = pM->s_hi();
		K_first.m_ii(i, f);
		K_len.m_ii(i, l);
		f += l;
		if (i == k - 1) {
			l = pM->s_li();
			K_first.m_ii(k, f);
			K_len.m_ii(k, l);
			}
		}
	

	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# stab_order_to_GAP(%s %s %s);\n\n", 
		KM_fname, GAP_fname, GAP_variable);
	fprintf(fp, "%s := [\n", GAP_variable);
	stab_go.fprint_GAP(fp);
	fprintf(fp, ",\n");
	K_first.fprint_GAP(fp);
	fprintf(fp, ",\n");
	K_len.fprint_GAP(fp);
	fprintf(fp, "];\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}
#endif

INT plesken_matrix_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT k_min, INT kk)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	MATRIX_OB Ainf, Ainf_inv, Ainf_block_wise;
	VECTOR_OB K_first, K_len;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	dc_plesken_matrices_prepare(k_min, kk, &MM, 
		&Ainf_block_wise, &Ainf, &Ainf_inv, 
		&K_first, &K_len, f_v, f_vv);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# plesken_matrix_to_GAP(%s %s %s %ld %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, kk);
	fprintf(fp, "%s := ", GAP_variable);
	Ainf.fprint_GAP(fp);
	fprintf(fp, ";\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT plesken_matrix_with_inverse_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT k_min, INT kk)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	MATRIX_OB Ainf, Ainf_inv, Ainf_block_wise;
	VECTOR_OB K_first, K_len;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	dc_plesken_matrices_prepare(k_min, kk, &MM, 
		&Ainf_block_wise, &Ainf, &Ainf_inv, 
		&K_first, &K_len, f_v, f_vv);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# plesken_matrix_to_GAP(%s %s %s %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, kk);
	fprintf(fp, "%s := [\n", GAP_variable);
	Ainf.fprint_GAP(fp);
	fprintf(fp, ", \n");
	Ainf_inv.fprint_GAP(fp);
	fprintf(fp, ", \n");
	K_first.fprint_GAP(fp);
	fprintf(fp, ", \n");
	K_len.fprint_GAP(fp);
	fprintf(fp, "];\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT generators_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# generators_to_GAP(%s %s %s);\n\n", 
		KM_fname, GAP_fname, GAP_variable);
	fprintf(fp, "%s := ", GAP_variable);
	G_gen.fprint_GAP(fp);
	fprintf(fp, ";\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT orbit_representatives_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT f_with_stabilizer_orders)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# orbit_representatives_to_GAP(%s %s %s);\n\n", 
		KM_fname, GAP_fname, GAP_variable);
	if (f_with_stabilizer_orders) {
		fprintf(fp, "%s := [\n", GAP_variable);
		RR.fprint_GAP(fp);
		fprintf(fp, ",\n");
		stab_go.fprint_GAP(fp);
		fprintf(fp, "]");
		}
	else {
		fprintf(fp, "%s := ", GAP_variable);
		RR.fprint_GAP(fp);
		}
	fprintf(fp, ";\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT stabilizer_orders_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2; //, S, s;
	FILE *fp;
	// INT i, l, j, ll, i0;
	VECTOR_OP R;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# orbit_representatives_to_GAP(%s %s %s);\n\n", 
		KM_fname, GAP_fname, GAP_variable);
#if 0
	l = RR.s_li();
	S.m_il(l);
	i0 = 0;
	for (i = 0; i < l; i++) {
		R = (VECTOR_OP) RR.s_i(i);
		ll = R->s_li();
		s.m_il(ll);
		for (j = 0; j < ll; j++) {
			if (i0 >= stab_go.s_li()) {
				printf("i=%ld l=%ld j=%ld ll=%ld i0=%ld stab_go.s_li()=%ld\n", 
					i, l, j, ll, i0, stab_go.s_li());
				stab_go.println();
				return error("stabilizer_orders_to_GAP() i0 >= stab_go.s_li()");
				}
			stab_go.s_i(i0)->swap(s.s_i(j));
			i0++;
			}
		s.swap((VECTOR_OP) S.s_i(i));
		}
	fprintf(fp, "%s := ", GAP_variable);
	S.fprint_GAP(fp);
#endif
	stab_go.fprint_GAP(fp);
	
	fprintf(fp, ";\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT orbit_representatives_to_ASCII(BYTE *KM_fname, 
	BYTE *export_fname, INT kk, INT f_with_stabilizer_orders)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	INT i, j, l, ll;
	VECTOR_OP R, S;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	FILE *fp;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	fp = fopen(export_fname, "w");
	// fprintf(fp, "# output from discreta_batch\n");
	// fprintf(fp, "# due to a call to\n");
	// fprintf(fp, "# orbit_representatives_to_GAP(%s %s %s);\n\n", 
	// 	KM_fname, GAP_fname, GAP_variable);
#if 0
	if (f_with_stabilizer_orders) {
		fprintf(fp, "%s := [\n", GAP_variable);
		RR.fprint_GAP(fp);
		fprintf(fp, ",\n");
		stab_go.fprint_GAP(fp);
		fprintf(fp, "]");
		}
	else {
		fprintf(fp, "%s := ", GAP_variable);
		RR.fprint_GAP(fp);
		}
#endif
	if (kk > RR.s_li())
		return error("orbit_representatives_to_ASCII() kk > RR.s_li()");
	R = (VECTOR_OP) RR.s_i(kk);
	fprintf(fp, "%ld %ld\n", R->s_li(), kk);
	l = R->s_li();
	for (i = 0; i < l; i++) {
		S = (VECTOR_OP) R->s_i(i);
		ll = S->s_li();
		for (j = 0; j < ll; j++) {
			fprintf(fp, "%ld ", S->s_ii(j));
			}
		fprintf(fp, "\n");
		}
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

#ifndef BUFSIZE
#define BUFSIZE 10000
#endif

INT nb_of_solutions_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT lambda)
{
	FILE *fp;
	INT nb_sol;
	
	nb_sol = km_nb_of_solutions(KM_fname, lambda);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# nb_of_solutions_to_GAP(%s %s %s %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, lambda);
	
	if (nb_sol == -1) {
		printf("no soutions for lambda=%ld found in file %s\n", 
			lambda, KM_fname);
		}
	
	fprintf(fp, "%s := %ld;\n", GAP_variable, nb_sol);
	fclose(fp);
	fflush(stdout);
	return OK;
	
}


INT read_until_lambdaend(FILE *fp)
{
	BYTE buf[BUFSIZE];
	INT l;
	
	// search for LAMBDAEND:
	while (TRUE) {
		if (fgets(buf, BUFSIZE, fp) == NULL)
			break;
		l = strlen(buf);
		if (l)
			buf[l - 1] = 0;
		if (strncmp(buf, "LAMBDAEND", 9) == 0)
			break;
		}
	return OK;
}

INT solutions_to_GAP(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, INT lambda, INT from, INT len)
{
	FILE *fp;
	VECTOR_OB S;
	VECTOR_OP s;
	INT i, j, a, l;
	BYTE buf[BUFSIZE];
	
	km_get_solutions(KM_fname, lambda, from, len, &S);
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# nb_of_solutions_to_GAP(%s %s %s %ld %ld %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, lambda, from, len);
	
	fprintf(fp, "%s := [\n", GAP_variable);
	for (i = 0; i < len; i++) {
		s = (VECTOR_OP) S.s_i(i);
		l = s->s_li();
		for (j = 0; j < l; j++) {
			if (s->s_ii(j))
				buf[j] = '1';
			else
				buf[j] = '0';
			}
		fprintf(fp, "\"%s\"", buf);
		if (i < len - 1)
			fprintf(fp, ",\n");
		}
	fprintf(fp, "\n];\n");
		
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT normalizing_permutation_on_orbits(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, BYTE *from_GAP, INT on_k_sets)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	VECTOR_OP pR;
	PERMUTATION_OB p, q;
	SYM_OB go;
	LABRA_OB lab_G;
	FILE *fp, *fp_in;
	INT i, l, a;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	fp_in = fopen(from_GAP, "r");
	fscanf(fp_in, "%ld", &l);
	if (l > v)
		return error("normalizing_permutation_on_orbits() l > v");
	p.m_il(v);
	for (i = 0; i < v; i++) {
		p.m_ii(i, i + 1);
		}
	for (i = 0; i < l; i++) {
		fscanf(fp_in, "%ld", &a);
		p.m_ii(i, a);
		}
	fclose(fp_in);
	printf("read permutation from file %s:", from_GAP);
	p.println();
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# normalizing_permutation_on_orbits(%s %s %s %s %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, from_GAP, on_k_sets);
	if (on_k_sets > RR.s_li()) {
		return error("normalizing_permutation_on_orbits() on_k_sets > RR");
		}
	pR = (VECTOR_OP) RR.s_i(on_k_sets);
	reduce_generators_labra(&G_gen, &go, 
		FALSE /* f_verbose */, &lab_G);
	km_normalizer_action_on_orbits(&lab_G, pR, &p, &q, FALSE /* f_v */),
	fprintf(fp, "%s := \n", GAP_variable);
	q.fprint_GAP(fp);
	fprintf(fp, ";\n");
	fclose(fp);
	fflush(stdout);
	return OK;
	
}

INT batch_fuse_orbits(BYTE *KM_fname, 
	BYTE *GAP_fname, BYTE *GAP_variable, BYTE *generators_fname, INT on_k_sets)
{
	MATRIX_OB Mtk;
	INT v, t, k, m, n;
	VECTOR_OB G_gen, RR, MM, stab_go, Orbits_below1, Orbits_below2;
	VECTOR_OP pR;
	SYM_OB H_go;
	VECTOR_OB H_gen;
	LABRA_OB lab_H;
	FILE *fp;
	INT i, l, a;
	VECTOR_OB new_rep_idx, new_reps;
	
	km_read_ascii(KM_fname, &Mtk, &v, &t, &k, 
		&G_gen, &RR, &MM, &stab_go, &Orbits_below1, &Orbits_below2);
	printf("read KM file, v=%ld t=%ld k=%ld\n", v, t, k);
	m = Mtk.s_hi();
	n = Mtk.s_li();
	printf("the KM matrix has size %ld x %ld\n", m, n);
	
	read_file_of_generators(&H_gen, generators_fname);
	reduce_generators_labra(&H_gen, &H_go, 
		FALSE /* f_verbose */, &lab_H);
	
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# fuse_orbits(%s %s %s %s %ld);\n\n", 
		KM_fname, GAP_fname, GAP_variable, generators_fname, on_k_sets);
	if (on_k_sets > RR.s_li()) {
		return error("normalizing_permutation_on_orbits() on_k_sets > RR");
		}
	pR = (VECTOR_OP) RR.s_i(on_k_sets);
	fuse_orbits(&lab_H, pR, &new_rep_idx, &new_reps, TRUE /* f_v */),
	fprintf(fp, "%s := \n", GAP_variable);
	fprintf(fp, "[ ");
	new_reps.fprint_GAP(fp);
	fprintf(fp, ",  ");
	l = new_rep_idx.s_li();
	for (i = 0; i < l; i++) {
		new_rep_idx.s_i(i)->inc();
		}
	new_rep_idx.fprint_GAP(fp);
	fprintf(fp, "];\n");
	fclose(fp);
	fflush(stdout);
	return OK;
}

INT batch_fuse_orbits_by_representatives(
	BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *reps_fname, BYTE *generators_fname)
{
	INT v, k, nb_reps, i, j, a, l;
	VECTOR_OB RR, R;
	SYM_OB H_go;
	VECTOR_OB H_gen;
	LABRA_OB lab_H;
	FILE *fp;
	VECTOR_OB new_rep_idx, new_reps;
	
	fp = fopen(reps_fname, "r");
	fscanf(fp, "%ld %ld\n", &nb_reps, &k);
	RR.m_il(nb_reps);
	R.m_il(k);
	for (i = 0; i < nb_reps; i++) {
		for (j = 0; j < k; j++) {
			fscanf(fp, "%ld", &a);
			R.m_ii(j, a);
			}
		R.copy((VECTOR_OP) RR.s_i(i));
		}
	fclose(fp);
	printf("read file %s with %ld representatives, k=%ld\n", reps_fname, nb_reps, k);

	read_file_of_generators(&H_gen, generators_fname);
	reduce_generators_labra(&H_gen, &H_go, 
		FALSE /* f_verbose */, &lab_H);
	
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# fuse_orbits_by_representatives(%s %s %s %s);\n\n", 
		GAP_fname, GAP_variable, reps_fname, generators_fname);
	fuse_orbits(&lab_H, &RR, &new_rep_idx, &new_reps, TRUE /* f_v */),
	fprintf(fp, "%s := \n", GAP_variable);
	fprintf(fp, "[ ");
	new_reps.fprint_GAP(fp);
	fprintf(fp, ",  ");
	l = new_rep_idx.s_li();
	for (i = 0; i < l; i++) {
		new_rep_idx.s_i(i)->inc();
		}
	new_rep_idx.fprint_GAP(fp);
	fprintf(fp, "];\n");
	fclose(fp);
	fflush(stdout);
	return OK;
}

INT id_fname_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *db_prefix, BYTE *id_fname)
{
	FILE *fp, *fp_in;
	INT l, i, id;

	fp_in = fopen(id_fname, "r");
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# geo_db_export_id(%s %s %s %s);\n\n", 
		GAP_fname, GAP_variable, db_prefix, id_fname);
	fprintf(fp, "%s := [", GAP_variable);
	fscanf(fp_in, "%ld", &l);
	for (i = 0; i < l; i++) {
		fscanf(fp_in, "%ld", &id);
		fprintf(fp, "%ld", id);
		if (i < l - 1) {
			fprintf(fp, ", \n");
			}
		}
	fprintf(fp, "];\n");
	fclose(fp_in);
	fclose(fp);
	fflush(stdout);
	return OK;
}


INT left_transversal_to_GAP(BYTE *GAP_fname, BYTE *GAP_variable, 
	BYTE *fname_G, BYTE *fname_H)
{
	SYM_OB go_G, go_H, ll;
	VECTOR_OB G_gen, H_gen;
	LABRA_OB lab_G, lab_H;
	VECTOR_OB R;
	PERMUTATION_OB rep;
	INT no, l, i;
	SINGLE_COSET_WORK *scw;
	FILE *fp;
	
	read_file_of_generators(&G_gen, fname_G);
	reduce_generators_labra(&G_gen, &go_G, FALSE /* f_verbose */, &lab_G);
	read_file_of_generators(&H_gen, fname_H);
	reduce_generators_labra(&H_gen, &go_H, FALSE /* f_verbose */, &lab_H);
	printf("left_transversal_to_GAP():");
	printf("go_G = ");
	go_G.println();
	printf("go_H = ");
	go_H.println();
	go_G.ganzdiv(&go_H, &ll);
	if (ll.s_obj_k() != INTEGER)
		return error("left_transversal_to_GAP() too many cosets");
	l = ll.s_i_i();
	printf("#cosets = %ld\n", l);
	
	scw = single_coset_open(&lab_G, &lab_H, TRUE);
	R.m_il(l);
	no = 0;
	single_coset_first(scw, &rep);
	do {
		if (no >= l)
			return error("left_transversal_to_GAP() no >= l, too many cosets");
		rep.swap(R.s_i(no));
		no++;
		printf("%ld ", no);
		if (no % 10 == 0)
			printf("\n");
		fflush(stdout);
		;
	} while (single_coset_next(scw, &rep));
	printf("\n");
	if (no < l)
		return error("left_transversal_to_GAP() no < l, too few cosets");
	
	single_coset_free(scw);
	
	fp = fopen(GAP_fname, "w");
	fprintf(fp, "# output from discreta_batch\n");
	fprintf(fp, "# due to a call to\n");
	fprintf(fp, "# left_transversal(%s %s %s %s);\n\n", 
		GAP_fname, GAP_variable, fname_G, fname_H);
	fprintf(fp, "%s := [\n", GAP_variable);
	for (i = 0; i < l; i++) {
		R.s_i(i)->fprint_GAP(fp);
		if (i < l - 1) {
			fprintf(fp, ", \n");
			}
		}
	fprintf(fp, "];\n");
	fclose(fp);
	fflush(stdout);
	return OK;
}









