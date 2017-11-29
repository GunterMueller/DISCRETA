#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "llll.h"

int main(int argc, char *argv[])
{
	int c0,c1,beta,p,sel;
	char *filename;
	char commandline[1024];

	if ((argc<6)||(argc>7)) {
		printf("Wrong number of parameters in command line input!\n");
		printf("discreta_lllpacking_shell c0 c1 beta p filename [SEL]\n");
		return 1;
	}
	
	c0 = atoi(argv[1]);
	c1 = atoi(argv[2]);
	beta = atoi(argv[3]);
	p = atoi(argv[4]);
	filename = argv[5];
	sel = 0;
	if (argc>=7) sel = atoi(argv[6]);
	
	if (argc>=7) {
		sprintf(commandline,"discreta_deletelines %d %d %d %d %s %d >infile",c0,c1,beta,p,filename,sel);
	} else {
		sprintf(commandline,"discreta_deletelines %d %d %d %d %s >infile",c0,c1,beta,p,filename);
	}
	system(commandline);
	
	if (argc>=7) {
		sprintf(commandline,"discreta_lllpacking %d %d %d %d infile %d",c0,c1,beta,p,sel);
	} else {
		sprintf(commandline,"discreta_lllpacking %d %d %d %d infile",c0,c1,beta,p);
	}
	system(commandline);

	system("rm infile");
	printf("finished.\n"); fflush(stdout);
	
	return 0;
}

/* \endc */
