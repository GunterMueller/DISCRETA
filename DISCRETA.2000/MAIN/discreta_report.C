/* discreta_report.C 
 * Anton Betten
 */

#include <stdio.h>
#include <stdlib.h>

#include <DISCRETA/discreta.h>
#include <DISCRETA/ladder.h>

int main(int argc, char **argv)
{
	INT t0, t1, user_time;
	BYTE s[256];

	if (argc < 2) {
		fprintf(stderr, "usage: %s KM_file\n", argv[0]);
		exit(1);
		}
	
	discreta_init();
	{
	BYTE *KM_fname;
	INT f_select = FALSE;
	INT i, design_select_lambda = 0, select_first = 0, select_length = 0;
	INT f_plesken = FALSE;
	INT f_html = FALSE;
	
	t0 = os_ticks();

	for (i = 1; i < argc - 1; i++) {
		if (strcmp(argv[i], "-plesken") == 0) {
			printf("-plesken\n");
			f_plesken = TRUE;
			}
		else if (strcmp(argv[i], "-html") == 0) {
			printf("-html\n");
			f_html = TRUE;
			}
#if 1
		else if (strcmp(argv[i], "-select") == 0) {
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
		}
	KM_fname = argv[i];
	if (f_plesken) {
		design_report_plesken(KM_fname);
		}
	else {
		design_report(KM_fname, f_html, f_select, 
			design_select_lambda, select_first, select_length);
		}

	
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



