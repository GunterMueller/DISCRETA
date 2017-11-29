/* 
  $Header: /usr/local/cvsroot/largesets/largeset.c,v 1.6 1997/12/08 15:08:25 alfred Exp $ 
  largeset.c
*/  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define zlength 16000

FILE *KMfile;
FILE *solfile;
FILE *lsfile;
FILE *mataus;

int nosolutions;
int KMs, KMz, N;
long lambdamax, lambda;
char *filename;
char zeile[zlength];
char zeile2[zlength];
char cmdline[200];

/*-----------------------------------------------------------------------*/
/* bin computes {n over k} with n >= k >= 0 */
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

/*-----------------------------------------------------------------------*/
void main(int argc, char *argv[]) {
	long v,k,t;
	int found;
	int i,j,l, n,m;
	char x;
	int *mattrans;
	long C0;
	int BETA, BKZP;
		
	if (argc!=6) {
		printf("Wrong parameters! Correct are: filename N cofactor beta p\n");
		exit(0) ;
	}
	filename = argv[1];
	N = atoi(argv[2]);

	C0 = atol(argv[3]);
	BETA = atol(argv[4]);
	BKZP = atol(argv[5]);

	/* Find v,t,k, and size of the matrix in KM-file */
	KMfile = fopen(filename,"r");
	found = 0;
	do {
		fgets(zeile,zlength,KMfile);
		if (strstr(zeile,"% degree:")!=NULL) {
			found++;
			sscanf(zeile,"%% degree: %ld",&v);
		}
		if (strstr(zeile,"% t, k:")!=NULL) {
			found++;
			fgets(zeile,zlength,KMfile);
			sscanf(zeile,"%% %ld %ld",&t,&k);
		}
	} while (zeile[0]=='%');
	sscanf(zeile,"%d%d",&KMz,&KMs);
	fclose(KMfile);
	if (found<2) {
		printf("v,t,k not found in KM-file -> exit\n");
		exit(0);
	}
	printf("v=%ld t=%ld k=%ld\n",v,t,k);
	lambdamax = bin(v-t,k-t);
	lambda = lambdamax / N;
	if (lambdamax%N!=0) {
		printf("No largesets possible: lambdamax=%ld N=%d -> exit\n",lambdamax,N);
		exit(0);
	}
	printf("%ld\n", lambdamax);
	
	/* Start discreta_lll_with */
	sprintf(cmdline,"discreta_lll %ld %ld bkz %d %d %s",C0*lambda,lambda,BETA,BKZP,filename);
	printf("%s\n",cmdline);
	i = system(cmdline);
	i = system("cp solutions ls_solutions");
	i = system("wc ls_solutions >ls_tmp");
	lsfile = fopen("ls_tmp","r");
	fscanf(lsfile,"%d",&nosolutions);
	fclose(lsfile);
	if (nosolutions==0) {
		printf("No designs found -> exit\n");
		exit(0);
	} else {
		printf("%d solutions found\n",nosolutions);
	}

	/* Transponation of solutions and writing the new
     matrix as KM-file */
	m = nosolutions;
	n = KMs;
	solfile = fopen("ls_solutions","r");
	mattrans = (int*)calloc (m*n,sizeof(int));
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++) {
			fscanf(solfile,"%c",&x);
			if(x == '1') mattrans[j*m+i] = 1;
			if(x == '0') mattrans[j*m+i] = 0;
		}
		fscanf(solfile,"%c",&x);
	}
	fclose(solfile);
	mataus = fopen("KMtmp","w"); 
	fprintf(mataus,"%% stopafter 1\n"); 
	fprintf(mataus,"%%\n%%\n%d %d \n",n+1,m);

	for(j=0;j<m;j++) fprintf(mataus,"%d ",1);
	fprintf(mataus,"\n");
	for(i=0;i<n;i++) {
		for(j=0;j<m;j++) fprintf(mataus,"%d ",N*mattrans[i*m+j]);
		fprintf(mataus,"\n");
	}
	fclose(mataus);

	/* Start discreta_mckay for large sets*/
/*	sprintf(cmdline,"discreta_lll_with %d %d %d %d KMtmp",C0*N,N,BETA,BKZP);*/
/*	sprintf(cmdline,"discreta_mckay mckay.log %d <KMtmp",N);*/
	/* Start discreta_spread for large sets*/
	sprintf(cmdline,"discreta_spread %d KMtmp",N);
	printf("%s\n",cmdline);
	i = system(cmdline);
	i = system("wc solutions >ls_tmp");
	lsfile = fopen("ls_tmp","r");
	fscanf(lsfile,"%d",&nosolutions);
	fclose(lsfile);
	if (nosolutions==0) {
		printf("No large sets found -> exit\n");
		exit(0);
	} else {
		printf("%d large set(s) found\n",nosolutions);
	}

	/* Output of solutions */
	solfile = fopen("solutions","r");
	printf("\nLarge set(s) found:\n");
	for(i=0;i<nosolutions;i++) {
		fgets(zeile,zlength,solfile);
		for(j=0;j<m;j++) {
			if(zeile[j]=='1') {
				for (l=0;l<n;l++) printf("%d",mattrans[l*m+j]);
				printf("\n");
			}
		}
		printf("\n");
	}
	free(mattrans);
	
	exit(0);
}
