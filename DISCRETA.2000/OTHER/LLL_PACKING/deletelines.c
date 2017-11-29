#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "llll.h"
#define zlength 16000

int selection;
int *permutation;
int inequalities;
int slackvar;
int c1_global;
int s_sel;
long stopafter;
long stoploops;
DOUBLE DELTA;
char *filename;

long lambda;
int KMs_org;
int glob_v, glob_t, glob_k;
long grouporder;

int design_eingabe(int ***b, long **x, int *s, int *z, int *KMs, int *KMz, long nom, long denom, int selection, int **permutation, int *s_sel, char *filename);
/* \section{Main program} */
int main(int argc, char *argv[])
{
	int i,j;
	int **b;
	long *x;
	int s,z;
	int KMs,KMz;

	long c0,c1;
	long nom, denom;
	
	int flag;

/* Here is the input of the command line parameters tested. We need exactly 9
   parameters. 
*/
	if ((argc<6)||(argc>7)) {
		printf("Wrong number of parameters in command line input!\n");
		printf("design c0 c1 beta p filename [SEL]\n");
		printf("design 1000 lambda 120 14 KM.dat 1\n");
		return 1;
	}

/* Test if they are well chosen. There a more sophisticated test will be 
   welcome. */
	for (i=1;i<=4;i++) {
		if ((atoi(argv[i])<=0)) {
			printf("Something is wrong with input parameter %d in the command line!\n",i);
			printf("design c0 c1 beta p filename [SEL]\n");
			printf("design 1000 lambda 120 14 KM.dat 1\n");
			return 1;
		}
	}
	
/*
The constants $c_0$, $c_1$, $c_2$. These constants are controlling
very sensitive the behaviour of the \|bkz|-algorithm.
*/
	c0 = atoi(argv[1]);
	c1 = atoi(argv[2]);
	lambda = c1;
	
	nom = 1;
	denom = 2;
	filename = argv[5];
	if (argc>6) {
		selection = atoi(argv[6]);
		if (selection<0) selection = 0;
	} else {
		selection = 0;
	}
	s_sel = 0;
	DELTA = LLLCONST_LOW;
	
/*   Read $b$ and $x$. $b$ has $s-1$ columns plus one zero column
     for LLL reduction.
*/
	flag = design_eingabe(&b,&x,&s,&z,&KMs,&KMz,nom,denom,selection,&permutation,&s_sel,filename);
	
	printf("%%\n");
	printf("%% %s\n",filename);
	printf("%%\n");
	printf("%% degree: %d\n",glob_v);
	printf("%% t, k:\n");
	printf("%% %d %d\n",glob_t,glob_k);
	printf("%% order: %ld\n",grouporder);
	printf("%%\n");
	
	printf("%d %d\n",KMz,KMs);
	for (i=0;i<KMz;i++) {
		for (j=0;j<KMs_org;j++) if (b[j][z-2]==1) printf("%d ",b[j][i+1]);
		printf("\n");
	}
	printf("\n");
	printf("STABILIZER-ORDER-K-SETS\n");
	for (j=0;j<KMs_org;j++) if (b[j][z-2]==1) printf("%d ",b[j][z-1]);
	
	
	return 0;
	/* finito!!! */
}

/* \section{Special functions}
\subsection{Scaling function}
The Kramer-Mesner matrix is multiplied by $c_0$.
The identity matrix which is appended to the Kramer-Mesner matrix is
scaled with $c_1$. The last two rows are multiplied by $c_2$.
*/
int design_eingabe(int ***b,long **x, int *s, int *z, int *KMs, int *KMz, long nom, long denom, int selection, int **permutation, int *s_sel, char *filename)
{
	int i,i0,j,l,flag;
	char zeile[zlength];
	char select[100];
	int *swap;
	int swap2;
	char *zp;
	FILE *txt;

   txt = fopen(filename, "r");
	sprintf(select,"%% with RHS");
	inequalities = 0;
	stopafter = 0;
	stoploops = 0;
	do {
		fgets(zeile,zlength,txt); 
		if (strstr(zeile,"% degree:")!=NULL) {
			sscanf(zeile,"%% degree: %d",&glob_v);
		}
		if (strstr(zeile,"% t, k:")!=NULL) {
			fgets(zeile,zlength,txt);
			sscanf(zeile,"%% %d %d",&glob_t,&glob_k);
		}
      if (strstr(zeile,"% order:")!=NULL) {
			sscanf(zeile,"%% order: %ld",&grouporder);
		}
      if (strstr(zeile,select)!=NULL) inequalities = 1;
      if (strstr(zeile,"% stopafter")!=NULL) {
			sscanf(zeile,"%% stopafter %ld",&stopafter);
		}
      if (strstr(zeile,"% stoploops")!=NULL) {
			sscanf(zeile,"%% stoploops %ld",&stoploops);
		}
	}
	while (zeile[0]=='%');

/*printf("%ld\n",stoploops);*/
	
   /* \|flag|$=1$ means arbitrary right hand side. */
	flag = 0;       
	sscanf(zeile,"%d%d%d",&(*KMz),&(*KMs),&flag); 
	if (flag!=1) sscanf(zeile,"%d%d",&(*KMz),&(*KMs));

	KMs_org = (*KMs);
	if (inequalities==1) {
		(*KMs) += (*KMz);
	}
	(*z) = (*KMz)+5;
	/* Two additional columns: the $\lambda$-column and a zero-column. */
	(*s) = (*KMs)+3;


/* allocation */
	(*x) = (long*)calloc((*s),sizeof(long));
	for (i=0;i<(*s);i++) (*x)[i] = 0;

	(*b) = (int**)calloc((*s),sizeof(int*));
	for(j=0;j<(*s);j++) {
		(*b)[j] = (int*)calloc((*z)+1,sizeof(int));
		for (i=0;i<=(*z);i++) (*b)[j][i] = 0; 
	}
	
/* Read the Kramer-Mesner-matrix */
	for (j=0;j<(*KMz);j++) {
		for (i=0;i<KMs_org;i++) fscanf(txt,"%d",&((*b)[i][j+1]));
		if (flag==1) {
			fscanf(txt,"%d",&((*b)[(*KMs)][j+1]));
		} else {
			(*b)[(*KMs)][j+1] = 1; 
		}
	}

	sprintf(select,"STABILIZER-ORDER-K-SETS");
	do {
		zp=fgets(zeile,zlength,txt);  
	} while ((zp!=NULL)&&(strstr(zeile,select)==NULL));

	if (zp==NULL) {
		printf("%s not found\n",select); fflush(stdout);
		fclose(txt); 
		exit(0);
	}
	for (i=0;i<(*KMs);i++) {
		fscanf(txt,"%d",&((*b)[i][(*z)-1]));
	}

	(*s_sel) = (*s);
	(*permutation) = (int*)calloc((*s),sizeof(int));
	i0=0;
	for (i=0;i<(*s_sel)-3;i++) {
		l = 1;
		for (j=0;j<(*KMz);j++) {
			if ((*b)[i][j+1] > lambda) {
				l = 0;
				break;
			}
		}
		if (l==0) {	
			(*b)[i][(*z)-2] = 0;
		} else {
			(*b)[i][(*z)-2] = 1;
			i0++;
		}
	}
	(*KMs) = i0;

	if (selection>0) {
		sprintf(select,"SELECTION %d",selection);
		do {
			zp=fgets(zeile,zlength,txt);  
		}
		while ((zp!=NULL)&&(strstr(zeile,select)==NULL));

		if (zp==NULL) {
			printf("%s not found\n",select);
			fclose(txt); 
			exit(0);
		}
		zp=fgets(zeile,zlength,txt); 
		fclose(txt); 
		if (zp==NULL) {
			printf("Selection vector not found\n");
			exit(0);
		}
		printf("%s:\n",select);
		printf("%s\n",zeile);
		
		(*permutation) = (int*)calloc((*s),sizeof(int));
		i0=0;
		zp = zeile;
		for (i=0;i<(*s_sel);i++) (*permutation)[i]=i;
		for (i=0;i<(*s_sel)-3;i++) {
			sscanf(zp,"%d",&l);
			zp=&(zp[2]);
			/* columns corresponding to 0 are deleted */
			if (l==0) {
				swap = (*b)[i0];
				swap2 = (*permutation)[i0];
				for (j=i0+1;j<(*s_sel);j++) {
					(*b)[j-1] = (*b)[j];
					(*permutation)[j-1] = (*permutation)[j];
				}
				(*b)[(*s_sel)-1] = swap;
				(*permutation)[(*s_sel)-1] = swap2;
				(*s)--;
			} else {
				i0++;
			}
		}
		(*KMs) = (*s)-3;
	}

	
	return flag;
}

/* \endc */
