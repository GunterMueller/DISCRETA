/*-------------------------------------------------
	Control program for computing t-designs
	Calls
		stripcolumns
		diophant
-------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define zlength 16000
#define MAX_INT 2147483647

char *filename;
long lambda;
int selection;

char commandline[256];
int iterate, iterate_no;
int bkz_beta,bkz_p;
long bkz_fac;
int silent;

int main(int argc, char *argv[])
{
	int i;
/*------------------------------------------------------
	Read the commandline parameters
------------------------------------------------------*/
	if ((argc<6)&&(argc>9)) {
		printf("Wrong number of parameters in command line input!\n");
		printf("solvedesign factor lambda [silent] bkz beta p filename [SEL]\n");
		printf("solvedesign 100 6 bkz 80 18 KM.dat 1\n");
		printf("solvedesign factor lambda [silent] iterate number filename [SEL]\n");
		printf("solvedesign 100 6 iterate 5 KM.dat 1\n");
		return 1;
	}
	bkz_fac = atol(argv[1]);
	lambda = atol(argv[2]);
	if (lambda<=0) lambda = MAX_INT;
	
	iterate = silent = 0;
	for (i=3;i<argc;i++) {
		if (strcmp(argv[i],"silent")==0) {
			silent = 1;
			printf("%s\n",argv[i]);
		}
		if (strcmp(argv[i],"iterate")==0) {
			iterate = 1;
			printf("%s\n",argv[i]);
			iterate_no = atoi(argv[i+1]);	
		}
		if (strcmp(argv[i],"bkz")==0) {
			iterate = 0;
			printf("%s\n",argv[i]);
			bkz_beta = atoi(argv[i+1]);	
			bkz_p = atoi(argv[i+2]);	
		}
	}	
	
	selection = 0;
	if (argc==7+silent) 
		selection = atoi(argv[argc-1]);
	
	if (selection>0)
		filename = argv[argc-2];
	else
		filename = argv[argc-1];

/*------------------------------------------------------
	Start stripcolumns
------------------------------------------------------*/
	if (selection==0) {
		sprintf(commandline,"stripcolumns %ld %s >diophant_infile.tmp",lambda,filename);
	} else {
		sprintf(commandline,"stripcolumns %ld %s %d>diophant_infile.tmp",
		                                           lambda,filename,selection);
	}
	printf("%s\n",commandline); fflush(stdout);
	system(commandline);

/*------------------------------------------------------
	Start diophant
------------------------------------------------------*/
#if defined(MPREC)
	sprintf(commandline,"solvediophant_mp %ld %ld ", bkz_fac,lambda);
#else
	sprintf(commandline,"solvediophant %ld %ld ", bkz_fac,lambda);
#endif															  
	if (silent)	sprintf(commandline,"%s silent ",commandline);
	if (iterate)	
		sprintf(commandline,"%s iterate %d diophant_infile.tmp",commandline,iterate_no);
	else 
		sprintf(commandline,"%s bkz %d %d diophant_infile.tmp",commandline,bkz_beta,bkz_p);

	printf("\n%s\n",commandline); fflush(stdout);
	system(commandline);

/*------------------------------------------------------
	Remove the temp. file diophant_infile.tmp
------------------------------------------------------*/
/*
	sprintf(commandline,"rm -f diophant_infile.tmp");
	printf("\n%s\n",commandline); fflush(stdout);
	system(commandline);
*/
	return 0;
}

