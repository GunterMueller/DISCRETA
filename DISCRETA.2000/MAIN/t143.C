// t143.C
//
// 02.12.1999
//
// transposes inc-files 
//
// 
//
// Anton Betten
// Bayreuth

#include <DISCRETA/discreta.h>

INT do_it(char *geo_fname1, char *geo_fname2);

int main(int argc, char **argv)
{
	INT t0, t1, user_time;
	BYTE s[256];

	if( argc != 3 ) {
		printf("usage: t143.out inc_file inc_file_transposed\n");
		return 1;
		}
	discreta_init();
	{
	// INT t0, t1, user_time;
	BYTE s[256];
	INT i, m, n;
	INT f_v = FALSE;
	INT f_vv = FALSE;
	INT f_vvv = FALSE;
	INT f_noimages = FALSE;
	BYTE *geo_fname1, *geo_fname2;
	
	t0 = os_ticks();
	
	geo_fname1 = argv[argc - 2];
	geo_fname2 = argv[argc - 1];

#if 0
	for (i = 1; i < argc - 2; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			f_v = TRUE;
			}
		else if (strcmp(argv[i], "-vv") == 0) {
			f_vv = TRUE;
			}
		else if (strcmp(argv[i], "-vvv") == 0) {
			f_vvv = TRUE;
			}
		else if (strcmp(argv[i], "-noimages") == 0) {
			f_noimages = TRUE;
			}
#if 0
		if (strcmp(argv[i], "-n") == 0) {
			sscanf(argv[i + 1], "%ld", &n);
			i++;
			}
#endif
		}
	sscanf(argv[argc - 2], "%ld", &m);
	sscanf(argv[argc - 1], "%ld", &n);
	do_it(m, n, f_v, f_vv, f_vvv, f_noimages);
#endif
	do_it(geo_fname1, geo_fname2);


	t1 = os_ticks();
	user_time = t1 - t0;
	s[0] = 0;
	print_delta_time_100(user_time, s);
	printf("total computing time: %s\n", s);
	fflush(stdout);
	
	}
	discreta_exit();
	return 0;
}

INT do_it(char *geo_fname1, char *geo_fname2)
{
	FILE *in_fp;
	FILE *out_fp;
	INT nrow, ncol, nb_X, *X1, *X2, a, i, j, k, l = 0;

	in_fp = fopen(geo_fname1, "r");
	if (in_fp == NIL) {
		printf("geo_fname1 = %s\n", geo_fname1);
		return error("can't open geo-file1 for reading");
		}
	out_fp = fopen(geo_fname2, "w");
	if (out_fp == NIL) {
		printf("geo_fname2 = %s\n", geo_fname2);
		return error("can't open geo-file2 for writing");
		}
	
	fscanf(in_fp, "%ld %ld %ld", &nrow, &ncol, &nb_X);
	fprintf(out_fp, "%ld %ld %ld\n", ncol, nrow, nb_X);
	
	printf("reading geometries of size %ld %ld with nb_X=%ld\n", 
		nrow, ncol, nb_X);
	X1 = (INT *) my_malloc(nrow * ncol * sizeof(INT), "X1");
	X2 = (INT *) my_malloc(nrow * ncol * sizeof(INT), "X2");
	while (TRUE) {

		for (i = 0; i < nrow * ncol; i++) {
			X1[i] = 0;
			X2[i] = 0;
			}

		fscanf(in_fp, "%ld", &a);
		if (a == -1)
			break;
		i = a / ncol;
		j = a % ncol;
		X1[i * ncol + j] = 1;
		X2[j * nrow + i] = 1;
		for (k = 1; k < nb_X; k++) {
			fscanf(in_fp, "%ld", &a);
			i = a / ncol;
			j = a % ncol;
			X1[i * ncol + j] = 1;
			X2[j * nrow + i] = 1;
			}
		for (i = 0; i < ncol; i++) {
			for (j = 0; j < nrow; j++) {
				if (X2[i * nrow + j])
					fprintf(out_fp, "%ld ", i * nrow + j);
				}
			}
		fprintf(out_fp, "\n");
		l++;
		}
	fprintf(out_fp, "-1 %ld geometries\n", l);
	return OK;
}





