#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define zlength 16000
#define MAX_INT 2147483647

char *filename;
char *zp;
long *stabs;

char zeile[zlength];
FILE *txt;
char detectstring[100];
long grouporder;

int no_stab;
long stab_size[1024];
int no_orbit[1024];

int main(int argc, char *argv[])
{
	int i,j,flag;
	int rows,columns;
	int stab_found;
	long stab;
	
	if (argc!=2) {
		printf("Wrong number of parameters in command line input!\n");
		printf("orbithint filename\n");
		printf("orbithint KM.dat\n");
		return 1;
	}

	filename = argv[1];

/*-------------------------------------------------
	Read the preamble of the file containing the
	linear system.
-------------------------------------------------*/
	
   txt = fopen(filename, "r");
	
	do {
		fgets(zeile,zlength,txt); 
      if (strstr(zeile,"% order:")!=NULL) {
			sscanf(zeile,"%% order: %ld",&grouporder);
		}
	}
	while (zeile[0]=='%');

/*-------------------------------------------------
	Read the matrix size.

	If a 3. parameter is given and positiv then
	the last column of the matrix is the RHS.
	
	RHSdirect is set in this case.
	The number of columns is given without the RHS !
-------------------------------------------------*/
	flag = 0;
	sscanf(zeile,"%d%d%d",&(rows),&(columns),&flag); 
	if (flag!=1) sscanf(zeile,"%d%d",&(rows),&(columns));

/*-------------------------------------------------
	Memory allocation of the matrix and the 
	selection vector.
-------------------------------------------------*/
	stabs = (long*)calloc(columns,sizeof(long));
	
/*-------------------------------------------------
	Search for the stabilizer order
-------------------------------------------------*/
	sprintf(detectstring,"STABILIZER-ORDER-K-SETS");
	do {
		zp=fgets(zeile,zlength,txt);  
	} while ((zp!=NULL)&&(strstr(zeile,detectstring)==NULL));

	if (zp!=NULL) {
		for (i=0;i<columns;i++) fscanf(txt,"%ld",&(stabs[i]));
		stab_found = 1;
	} else {
		printf("%s not found\n",detectstring); fflush(stdout); 
		stab_found = 0;
	}

/*-------------------------------------------------
	Close the input file.
-------------------------------------------------*/
	fclose(txt); 

	for(i=0;i<1024;i++) {
		stab_size[i] = 0;
		no_orbit[i] = 0;
	}
	
	no_stab = 0;
	for(i=0;i<columns;i++) {
		stab = stabs[i];
		flag = 0;
		for (j=0;j<no_stab;j++) {
			if (stab==stab_size[j]) {
				no_orbit[j]++;
				flag = 1;
				break;
			}
		}
		if (!flag) {
			stab_size[no_stab] = stab;
			no_orbit[no_stab] = 1;
			no_stab++;
		}
	}

	printf("There are the following types of orbits: \n");
	printf("len.\tnumber \n");
	for(i=0;i<no_stab;i++) {
		printf("%ld\t",grouporder/stab_size[i]);
		printf("%d",no_orbit[i]);	
		printf("\n");
	}
	printf("\n");

	return 0;
}

/* \endc */
