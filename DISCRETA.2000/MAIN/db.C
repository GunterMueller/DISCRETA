// db.C

void db_add(INT t, INT v, INT k, INT lambda)
{
	FILE *fp_tex;
	BYTE *comment = NIL;
	INT f_v = FALSE;
	DATABASE_OP db;
	MATRIX_OB data;
	
	i_db_dp(&db, DP_PATH);
	
	printf("t = %ld v = %ld k = %ld lambda = %ld\n", t, v, k, lambda);
	fflush(stdout);
	if (t == 0 || v == 0 || k == 0 || lambda == 0)
		return;
		
	// fp_tex = fopen("parameter_sets.tex", "w");
	fp_tex = stdout;
	fprintf(fp_tex, "\\noindent\n{\\bf Parameter sets obtained from %ld-(%ld,%ld,%ld)}\n", 
		t, v, k, lambda);
	fprintf(fp_tex, "\\[\n\\begin{array}{rll}\n");
	db_dp_design_parameter_produce(db, t, v, k, lambda, comment, f_v, fp_tex, &data);
	fprintf(fp_tex, "\\end{array}\n\\]\n");
	// fclose(fp_tex);
	
	e_db_dp(db);
	
	printf("adding %ld entries into multi-dimensional db\n", data.s_hi());
	fflush(stdout);
	despar_add(&data);
}

#define DB_FNAME "despar"

void db_add_vec_data(VECTOR_OP vec_data)
{
	BYTE *db_fname = DB_FNAME;
	INT i, l;
	MATRIX_OP data;
	// INT f_create, f_termination;

	l = vec_data->s_li();
	printf("***** db_add_vec_data() adding l = %ld sets of parameter sets\n", l);
	for (i = 0; i < l; i++) {
		printf("***** db_add_vec_data() i = %ld l = %ld\n", i, l);
		data = (MATRIX_OP) vec_data->s_i(i);
		if (data->s_obj_k() != MATRIX) {
			printf("i = %ld l = %ld\n", i, l);
			fflush(stdout);
			error("db_add_vec_data() s_obj_k() != MATRIX !");
			return;
			}
		printf("adding %ld entries into multi-dimensional db\n", data->s_hi());
		fflush(stdout);
#if 0
		if (i == 0)
			f_create = TRUE;
		else
			f_create = FALSE;
		if (i == l - 1)
			f_termination = TRUE;
		else
			f_termination = FALSE;
#endif
		despar_add(data);
		// despar_add(data, f_create, f_termination);
		}
	// despar_add2();
	despar_buf_flush(db_fname);
}

INT despar_add(MATRIX_OP data)
{
	BYTE *db_fname = DB_FNAME;
	INT i, l, t, v, k, lambda, id;
	
	l = data->s_hi();
	printf("despar_add() adding %ld data sets\n", l);
	for (i = 0; i < l; i++) {
		t = data->s_iji(i, 0);
		v = data->s_iji(i, 1);
		k = data->s_iji(i, 2);
		lambda = data->s_iji(i, 3);
		id = data->s_iji(i, 4);
		printf("%ld %ld %ld %ld %ld\n", t, v, k, lambda, id);
		despar_db_add(db_fname, t, v, k, lambda, id);
		}
	return OK;
}

#define DESPAR_BUF_SIZE 50

static DESPAR *despar_buf = NIL;
static INT despar_buf_len = 0;

INT despar_buf_flush(char *db_fname)
{
	db_add(despar_buf, despar_buf_len, db_fname);
	despar_buf_len = 0;
	return OK;
}

INT despar_db_add(char *db_fname, INT t, INT v, INT k, INT lambda, INT id)
{
	DESPAR dp;
	
	dp.t = t;
	dp.v = v;
	dp.k = k;
	dp.lambda = lambda;
	dp.id = id;
	if (despar_buf == NIL) {
		despar_buf = (DESPAR *) my_malloc(sizeof(DESPAR) * DESPAR_BUF_SIZE, "DESPAR buffer");
		despar_buf_len = 0;
		}
	if (despar_buf_len >= DESPAR_BUF_SIZE)
		return error("despar_db_add(): despar_buf_len >= DESPAR_BUF_SIZE");
	printf("despar_db_add() adding %ld-(%ld,%ld,%ld) id=%ld\n",dp.t, dp.v, dp.k, dp.lambda, dp.id);
	despar_buf[despar_buf_len++] = dp;
	if (despar_buf_len >= DESPAR_BUF_SIZE) {
		despar_buf_flush(db_fname);
		}
	// db_add(&dp, db_fname);
	return OK;
}

INT db_add(DESPAR *dp, INT len, char *db_fname)
{
	char *pp = NIL;
	INT size = 0, l, i;
	FILE *f;
	INT f_verbose = TRUE;
	BYTE in_file_name[256];
	BYTE str[256];
	
	sprintf(in_file_name, "%s.in.bin", db_fname);

	printf("adding %ld design parameter sets: ", len);
	despar_print(dp);
	printf("\n");
	fflush(stdout);
	despar2dbf(dp, len, &pp, &size);

	f = fopen(in_file_name, "wb+");
	fwrite(pp, 1 /* size */, size /* items */, f);
	fclose(f);
	if (f_verbose) {
		printf("db_add()|wrote file %s of size %ld\n", 
			in_file_name, file_size(in_file_name));
		fflush(stdout);
		}
	
	sprintf(str, "%s from:%s to:%s", 
		LOC_DB_INSERT, in_file_name, db_fname);
	call_system(str);
	
	if (pp) {
		my_free(pp);
		pp = NIL;
		}
	return OK;
}

#ifndef BUFSIZE
#define BUFSIZE 1024
#endif

INT despar_search(
	INT t_min, INT t_max, 
	INT v_min, INT v_max, 
	INT k_min, INT k_max, 
	INT lambda_min, INT lambda_max)
{
	BYTE buf[BUFSIZE];
	BYTE str[1024];
	FILE *fp;
	INT i, t, v, k, lambda, id;
	DESPAR *dp = NIL;
	INT dp_len;
	DATABASE_OP db;
	DESIGN_PARAMETER_OB dp_obj;
	INT prev;
	BYTE s0[10000];
	BYTE s1[1024];
	BYTE s2[1024];
	BYTE s3[1024];
	
	fp = fopen("despar.sdf", "w");
	fprintf(fp, "from despar\n");
	fprintf(fp, "to despar.out.bin\n");
	fprintf(fp, "where t\n");
	fprintf(fp, ">= %ld\n", t_min);
	fprintf(fp, "<= %ld\n", t_max);
	fprintf(fp, "where v\n");
	fprintf(fp, ">= %ld\n", v_min);
	fprintf(fp, "<= %ld\n", v_max);
	fprintf(fp, "where k\n");
	fprintf(fp, ">= %ld\n", k_min);
	fprintf(fp, "<= %ld\n", k_max);
	fprintf(fp, "where lambda\n");
	fprintf(fp, ">= %ld\n", lambda_min);
	fprintf(fp, "<= %ld\n", lambda_max);
	fclose(fp);
	sprintf(str, "rm despar_searchresult.txt");
	call_system(str);
	sprintf(str, LOC_DESPAR_OUT " get despar.sdf despar.out.bin >despar_searchresult.txt", DISCRETA_ARCH);
	call_system(str);
	despar_db_get("despar.sdf", "despar.out.bin");
	
	i_db_dp(&db, DP_PATH);
	despar_read("despar.out.bin", &dp, &dp_len);
	for (i = 0; i < dp_len; i++) {
		t = dp[i].t;
		v = dp[i].v;
		k = dp[i].k;
		lambda = dp[i].lambda;
		id = dp[i].id;
		// printf("%ld-(%ld,%ld,%ld) id = %ld\n", t, v, k, lambda, id);
		// fflush(stdout);
		
		db_dp_load_id(db, id, &dp_obj);
#if 0
		printf("%ld-(%ld,%ld,%ld) (id=%ld)\n", 
			dp_obj.s_t_i(),
			dp_obj.s_v_i(),
			dp_obj.s_k_i(),
			dp_obj.s_lambda_i(),
			dp_obj.s_id_i());
#endif
		s0[0] = 0;
		dp_obj.sprint(s0);
		printf("%ld: %s (id=%ld)\n", i, s0, dp_obj.s_id_i());

#if 0
		prev = dp_obj.s_prev_i();
		s0[0] = 0;
		s1[0] = 0;
		s2[0] = 0;
		dp_obj.sprint_text012(s0, s1, s2);
		printf("%ld: %s %s ", i, s0, s1);
		if (prev >= 0)
			printf("%ld ", prev);
		printf(" %s (id=%ld)\n", s2, dp_obj.s_id_i());
		// dp_obj.println();
#endif
		}
	if (dp) {
		my_free(dp);
		dp = NIL;
		}
	e_db_dp(db);

#if 0
	fp = fopen("despar_searchresult.txt", "r");
	while(fgets(buf, BUFSIZE, fp) != NULL) {
		printf(buf);
#if 0
		sscanf(buf, "%ld %ld %ld %ld %ld", &t, &v, &k, &lambda, &id);
		printf("%ld-(%ld,%ld,%ld) id = %ld\n", t, v, k, lambda, id);
#endif
		}
	fclose(fp);
#endif
	return OK;
}

INT despar_delete(
	INT t_min, INT t_max, 
	INT v_min, INT v_max, 
	INT k_min, INT k_max, 
	INT lambda_min, INT lambda_max)
{
	BYTE buf[BUFSIZE];
	BYTE str[1024];
	FILE *fp;
	INT i, t, v, k, lambda, id;
	DESPAR *dp = NIL;
	INT dp_len;
	DATABASE_OP db;
	DESIGN_PARAMETER_OB dp_obj;
	DATATYPE data;
	INT prev;
	BYTE s0[10000];
	BYTE s1[1024];
	BYTE s2[1024];
	BYTE s3[1024];
	
	fp = fopen("despar.sdf", "w");
	fprintf(fp, "from despar\n");
	fprintf(fp, "to despar.out.bin\n");
	fprintf(fp, "where t\n");
	fprintf(fp, ">= %ld\n", t_min);
	fprintf(fp, "<= %ld\n", t_max);
	fprintf(fp, "where v\n");
	fprintf(fp, ">= %ld\n", v_min);
	fprintf(fp, "<= %ld\n", v_max);
	fprintf(fp, "where k\n");
	fprintf(fp, ">= %ld\n", k_min);
	fprintf(fp, "<= %ld\n", k_max);
	fprintf(fp, "where lambda\n");
	fprintf(fp, ">= %ld\n", lambda_min);
	fprintf(fp, "<= %ld\n", lambda_max);
	fclose(fp);
	// sprintf(str, "rm despar_searchresult.txt");
	// call_system(str);
	// sprintf(str, LOC_DESPAR_OUT " get despar.sdf despar.out.bin >despar_searchresult.txt", DISCRETA_ARCH);
	// call_system(str);
	despar_db_get("despar.sdf", "despar.out.bin");
	
	i_db_dp(&db, DP_PATH);
	despar_read("despar.out.bin", &dp, &dp_len);
	for (i = 0; i < dp_len; i++) {
		t = dp[i].t;
		v = dp[i].v;
		k = dp[i].k;
		lambda = dp[i].lambda;
		id = dp[i].id;
		// printf("%ld-(%ld,%ld,%ld) id = %ld\n", t, v, k, lambda, id);
		// fflush(stdout);
		db_dp_load_id_data(db, id, &dp_obj, &data);
		s0[0] = 0;
		dp_obj.sprint(s0);
		printf("%ld: %s (id=%ld)...", i, s0, dp_obj.s_id_i());

		db->open();
		db->del_op(&dp_obj, data.datref, data.data_size);
		db->close();
		printf("deleted\n");

		}
	despar_db_delete("despar.out.bin");
	printf("deleted in gridfile\n");
	

	if (dp) {
		my_free(dp);
		dp = NIL;
		}
	e_db_dp(db);

	return OK;
}

INT despar_db_get(char *sdf_fname, char *out_fname)
{
	/* BYTE out_file_name[256]; */
	BYTE str[256];
	DESPAR *dp = NIL;
	INT i, len;

	sprintf(str, "%s %s", LOC_DB_SELECT, sdf_fname);
	call_system(str);
	despar_read(out_fname, &dp, &len);
#if 0
	for (i = 0; i < len; i++) {
		printf("%ld: ", i);
		despar_print(dp + i);
		printf("\n");
		}
	fflush(stdout);
#endif
	if (dp) {
		my_free(dp);
		dp = NIL;
		}
	return OK;
}

INT despar_db_delete(char *fname)
{
	BYTE str[256];
	DESPAR *dp = NIL;
	INT i, len;

	sprintf(str, "%s %s despar", LOC_DB_DELETE, fname);
	call_system(str);
	return OK;
}

INT despar_read(BYTE *file_name, DESPAR **p_dp, INT *len)
{
	FILE *f;
	INT size, rec_len, rec_len1, i, l;
	BYTE *pp = NIL;
/*	INT f_verbose = TRUE; */
	INT f_verbose = FALSE;
	DESPAR dp1, *dp;
	
	size = file_size(file_name);
	if (size == 0) {
		*len = 0;
		return OK;
		}
	if (f_verbose) {
		printf("despar_read()|reading file %s of size %ld\n", 
			file_name, size);
		fflush(stdout);
		}
	pp = (BYTE *) my_malloc(size, "despar_read");
	f = fopen(file_name, "rb");
	fread(pp, 1 /* size */, size /* items */, f);
	fclose(f);
	l = 0;
	rec_len1 = 0;
	while (rec_len1 < size) {
		dbf2despar(pp + rec_len1, &dp1, &rec_len);
		l++;
		rec_len1 += rec_len;
		}
	if (f_verbose) {
		printf("despar_read()|reading %ld structs\n", l);
		fflush(stdout);
		}
	dp = (DESPAR *) my_malloc(l * sizeof(DESPAR), "despar_read");
	*p_dp = dp;
	rec_len1 = 0;
	for (i = 0; i < l; i++) {
		dbf2despar(pp + rec_len1, dp + i, &rec_len);
		rec_len1 += rec_len;
		}
	*len = l;
	if (pp) {
		my_free(pp);
		pp = NIL;
		}
	return OK;
}

INT despar2dbf(DESPAR *dp, INT len, BYTE **p_data, INT *p_data_len)
{
	unsigned int *pi;
	char *pp = NIL;
	INT i, h;
	INT nb_index = 5; /* t, v, k, lambda, id */
	INT size_pos, size_key, size_info;
	INT size_head; /* = size_pos + size_key */
	INT size_rec; /* = size_head + size_info */
	
	size_pos = (nb_index + 2) * sizeof(int);
	size_key = nb_index * sizeof(int);
	size_head = size_pos + size_key;
	size_info = 0; /* no info part */
	size_rec = size_head + size_info;

	pp = (char *) my_malloc(size_rec * len, "despar2dbf");
	
	for (h = 0; h < len; h++) {
		pi = (unsigned int *) (pp + size_rec * h);

		pi[0] = size_rec;
		for (i = 0; i < nb_index; i++) {
			pi[1 + i] = size_pos + i * sizeof(int);
			}
		pi[1 + nb_index] = size_head;
		pi[2 + nb_index + 0] = (int) dp[h].t;
		pi[2 + nb_index + 1] = (int) dp[h].v;
		pi[2 + nb_index + 2] = (int) dp[h].k;
		pi[2 + nb_index + 3] = (int) dp[h].lambda;
		pi[2 + nb_index + 4] = (int) dp[h].id;
		/* no info part ! */
		}

	*p_data = pp;
	*p_data_len = size_rec * len;
	return OK;
}

INT dbf2despar(char *pp, DESPAR *dp, INT *rec_len)
{
	unsigned int *pi = (unsigned int *) pp;
	INT i;
	INT nb_index = 5;
	
	if (pp == NIL) {
		printf("dbf2despar() pp == NIL\n");
		return ERROR;
		}
	if (rec_len)
		*rec_len = pi[0];
	dp->t = pi[2 + nb_index + 0];
	dp->v = pi[2 + nb_index + 1];
	dp->k = pi[2 + nb_index + 2];
	dp->lambda = pi[2 + nb_index + 3];
	dp->id = pi[2 + nb_index + 4];
	return OK;
}

void despar_db_create(char *fname)
{
	BYTE str[256];
	BYTE tbl_fname[256];
	BYTE tdf_fname[256];
	
	sprintf(tbl_fname, "%s.tbl", fname);
	sprintf(tdf_fname, "%s.tdf", fname);
	if (file_size(tbl_fname) >0) {
		sprintf(str, "rm %s", tbl_fname);
		call_system(str);
		}
	if (file_size(tdf_fname) >0) {
		sprintf(str, "rm %s", tdf_fname);
		call_system(str);
		}
	sprintf(str, "%s %s.def", LOC_DB_CREATE, fname);
	call_system(str);
}

void despar_print(DESPAR *dp)
{
	printf("%ld-(%ld, %ld, %ld) id=%ld", dp->t, dp->v, dp->k, dp->lambda, dp->id);
}


