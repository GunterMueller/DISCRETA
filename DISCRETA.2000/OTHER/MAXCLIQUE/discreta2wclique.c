#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define zlength 16000
#define MAXINC 10000

FILE *infile;
FILE *outfile;	
long stopafter, stoploops;
int lambda;
int numeqn, numvar;
int selection;
int **matrix;
int *weights;
long grouporder;

int numbervertices;
int numberedges;
int *vertices;

long bin (long n, long k) {
	int i;
	long b;
	
	if ((n<k)||(n<0)||(k<0)) {
		printf("n=%ld and k=%ld not correct for bin\n", n,k);
		exit(0);
	}
	b=1;
	for(i=0;i<k;i++) {
		b *= (n-i);
		b /= (i+1);
	}
	return b;
}

long scalarproduct(int *v, int *w, int len) {
	int i;
	long b;
	b = 0;
	for (i=0;i<len;i++) b += v[i]*w[i];
	return b;
}

long maxnorm(int *v, int len) {
	int i;
	int b;
	b = 0;
	for (i=0;i<len;i++) if (b<v[i]) b = v[i];
	return b;
/*	return 0;*/
}
	
/******************************************************************
*  The main program
*******************************************************************/
int main(int argc, char *argv[]) {
	int i,j,p,s,b;
	long v,t,k,blocks;
	char zeile[zlength];
	char *filename;
	char select[100];
	char *zp;
	
/******************************************************************
*  The command line parameters are read
*******************************************************************/
/*
	if ((argc!=3)&&(argc!=4)) {
		printf("Wrong parameters! Correct are: lambda filename [SEL]\n");
		exit(0) ;
	}
*/	
 	lambda = atoi(argv[1]);
	filename = argv[2];
	if (argc>3) {
		selection = atoi(argv[3]);
		if (selection<0) selection = 0;
	} else {
		selection = 0;
	}

/******************************************************************
* Look for maximal number of solutions or loops
*******************************************************************/
	stopafter = -1;
	stoploops = -1;
	infile = fopen(filename,"r");
	do {
		fgets(zeile,zlength,infile);
		if (strstr(zeile,"% degree:")!=NULL) {
			sscanf(zeile,"%% degree: %ld",&v);
		}
		if (strstr(zeile,"% t, k:")!=NULL) {
			fgets(zeile,zlength,infile);
			sscanf(zeile,"%% %ld %ld",&t,&k);
		}
      if (strstr(zeile,"% stopafter")!=NULL) {
			sscanf(zeile,"%% stopafter %ld",&stopafter);
		}
      if (strstr(zeile,"% stoploops")!=NULL) {
			sscanf(zeile,"%% stoploops %ld",&stoploops);
		}
      if (strstr(zeile,"% order:")!=NULL) {
			sscanf(zeile,"%% order: %ld",&grouporder);
		}
	} while (zeile[0]=='%');
	sscanf(zeile,"%d%d",&numeqn,&numvar);

/******************************************************************
* Allocate memory 
*******************************************************************/
	matrix = (int**)calloc(numvar,sizeof(int));
	for (i=0;i<numvar;i++) matrix[i] = (int*)calloc(numeqn,sizeof(int));
	weights = (int*)calloc(numvar,sizeof(int));
	vertices = (int*)calloc(numvar,sizeof(int));

/******************************************************************
* Read the matrix
*******************************************************************/
	for (i=0;i<numeqn;i++) {
		for (j=0;j<numvar;j++) fscanf(infile,"%d",&(matrix[j][i]));
	}
/*
	for (i=0;i<numeqn;i++) {
		for (j=0;j<numvar;j++) printf("%d ",matrix[j][i]);
		printf("\n");
	}
	printf("\n");
*/

/******************************************************************
* find the weights
*******************************************************************/
	sprintf(select,"STABILIZER-ORDER-K-SETS");
	do {
		zp=fgets(zeile,zlength,infile);  
	}
	while ((zp!=NULL)&&(strstr(zeile,select)==NULL));

	if (zp==NULL) {
		printf("%s not found\n",select);
		fclose(infile); 
		exit(0);
	}
	for (i=0;i<numvar;i++) {
		fscanf(infile,"%d",&s);
		weights[i] = grouporder/s;
/*		printf("%d ",weights[i]);*/
	}
/*	printf("\n");*/
	fclose(infile);
	
/******************************************************************
* delete colums with entries > lambda
*******************************************************************/
	for (i=0;i<numvar;i++) {
		if (maxnorm(matrix[i],numeqn) > lambda) {
			vertices[i] = 0;
/*			printf("delete column %d\n",i); */
		} else {
			vertices[i] = 1;
		}
	}

/******************************************************************
* count remaining columns
*******************************************************************/
	numbervertices = 0;
	for (i=0;i<numvar;i++) if (vertices[i] == 1) numbervertices++;
	
	outfile = fopen("infile","w");
	fprintf(outfile,"%d %d \n",numbervertices,1);          /* 1 is a dummy number of edges. */
	printf("Number of vertices %d \n",numbervertices); 

/******************************************************************
* for each vertice find the adjacent vertices
*******************************************************************/
	numberedges = 0;
	for (i=0;i<numvar;i++) if (vertices[i]==1) {
		fprintf(outfile,"%d ",weights[i]);
		b = 0;
		for (j=0;j<numvar;j++) if (i!=j && vertices[j]==1) {
			if (scalarproduct(matrix[i],matrix[j],numeqn)==0) b++;
		}
		fprintf(outfile,"%d  ",b);
		if (b>0) {
			p = 0;
			for (j=0;j<numvar;j++) if (vertices[j]==1) {
				if (i!=j && scalarproduct(matrix[i],matrix[j],numeqn)==0) {
					fprintf(outfile,"%d ",p);
					if (i>j) numberedges++;
				}
				p++;
			}
		} 
		fprintf(outfile,"\n");
	}
	fclose(outfile);
	printf("Number of edges %d \n",numberedges); 

	blocks = bin(v,t)/bin(k,t);	
	printf("Number of blocks in a %ld-(%ld,%ld,%d) design: %ld \n",t,v,k,lambda,blocks); 

	return 1;
}
