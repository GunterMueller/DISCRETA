/* 
  $Header: /usr/local/cvsroot/LARGESET/largeset_r.c,v 1.1.1.1 1998/01/15 08:51:13 alfred Exp $ 
  largeset_iterativ.c
*/  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define zlength 16000

FILE *KMfile;
FILE *solfile;
FILE *lsfile;
FILE *mataus;

long nosolutions;
int KMs, KMz, N;
long lambdamax, lambda;
char *filename;
char zeile[zlength];
char zeileold[zlength];
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
	long rmax;
	int i,j,iterat;
	int found;
	int m,n,l,n_new;
	long maxsolutions;
	long maxloops;
	int pick;
	int runde;
	int successful;
	int iterat0done;
	long C0,MAX_ROUND;
	int BETA,BKZP;

	if (argc!=9) {
		printf("Wrong parameters! Correct are: filename N maxsolutions maxloops maxrounds c0factor beta p\n");
		exit(0) ;
	}
	filename = argv[1];
	N = atoi(argv[2]);
	maxsolutions = atol(argv[3]);
	maxloops = atol(argv[4]);

	MAX_ROUND = atol(argv[5]);
	C0  = atol(argv[6]);
	BETA = atoi(argv[7]);
	BKZP = atoi(argv[8]);

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
	
	/* Compute lambdamax and lambda */
	lambdamax = bin(v-t,k-t);
	lambda = lambdamax / N;
	if (lambdamax%N!=0) {
		printf("No largesets possible: lambdamax=%ld N=%d -> exit\n",lambdamax,N);
		exit(0);
	}
	printf("lambdamax = %ld\n", lambdamax);
	
	/*-----------------------------------------------------------------*/
iterat0done = 0;
for (runde=0;runde<MAX_ROUND;runde++) {
	i = system("rm -f solutions.ls");
	for (j=0;j<KMs;j++) zeileold[j] = '0';
	KMfile = fopen("KMtmp","w");
	fprintf(KMfile,"%% stopafter %ld\n",maxsolutions);
	fprintf(KMfile,"%% stoploops %ld\n",maxloops);
	fclose(KMfile);
	sprintf(cmdline,"cat %s >>KMtmp",filename);
	i = system(cmdline);
	successful = 0;
	
	n = KMs;
	m = KMz;
/*	srand(time());*/
	for (iterat=0;iterat<N;iterat++) { 
		if ((iterat==0)&&(iterat0done == 1)) {
			i = system("cp solutions.1 solutions");
			printf("%d. step: already done\n",iterat+1);
		} else {
			sprintf(cmdline,"discreta_lll %ld %ld bkz %d %d KMtmp >>/dev/null",C0*lambda,lambda,BETA,BKZP);
			printf("%d. step: %s\n",iterat+1,cmdline);
			i = system(cmdline);
		}
		i = system("cp solutions ls_solutions");
		if ((iterat==0)&&(iterat0done == 0)) {
			 i = system("cp solutions solutions.1");
			 iterat0done = 1;
		}
		i = system("wc ls_solutions >ls_tmp");

		lsfile = fopen("ls_tmp","r");
		fscanf(lsfile,"%ld",&nosolutions);
		fclose(lsfile);
	
		if (nosolutions==0) {
			printf("No designs found in step %d -> exit\n", iterat+1);
			printf("Covered so far:\n%s\n",zeileold);			
			break;
		} else {
			printf("%ld solutions found in step %d \n",nosolutions, iterat+1);
		}

		rmax = (nosolutions<maxsolutions)?nosolutions:maxsolutions;
		pick = 0;
		if (rmax>RAND_MAX) pick = rand();
		pick = (int)((pick*RAND_MAX + rand()) % rmax);		
		printf("Pick No: %d\n",pick);
		
		/* Read from file solutions one design */
		solfile = fopen("ls_solutions","r");
		for (i=0;i<pick;i++) fgets(zeile, zlength/*n+1*/, solfile);
		fgets(zeile, n+1, solfile);
		fclose(solfile);

		n_new = 0;
		for (i=0;i<n;i++) {
			if(zeile[i]=='0') n_new++;
		}
		
		/*-----------------------------------------------------------------*/
		/* The reduced KM-file is written */
		i = system("cp KMtmp ls_tmp");
		KMfile = fopen("ls_tmp","r");
		mataus = fopen("KMtmp","w");
		fprintf(mataus,"%% stopafter %ld\n",maxsolutions);
		fprintf(mataus,"%% stoploops %ld\n",maxloops);
		fprintf(mataus,"%% degree: %ld\n",v);
		fprintf(mataus,"%% t, k:\n");
		fprintf(mataus,"%% %ld %ld\n",t,k);
		fprintf(mataus,"%%\n%%\n %d %d\n",KMz,n_new);
		do {
			fgets(zeile2,zlength,KMfile);			
		} while (zeile2[0]=='%');
		sscanf(zeile2,"%d%d",&m,&n);  /* m = rows, n = columns */
		for (i=0;i<m;i++) {
			for (j=0;j<n;j++) {
				fscanf(KMfile,"%d",&l);
				if (zeile[j]=='0') fprintf(mataus,"%d ",l);
			}
			fprintf(mataus,"\n");
		}
		fclose(mataus);
		fclose(KMfile);
		
		lsfile = fopen("ls_tmp","w");
		if (iterat>0) {
			l = 0;
			for (j=0;j<KMs;j++) {
				if (zeileold[j] == '0') {
					fprintf(lsfile,"%c",zeile[l]); 
					if (zeile[l]=='1') zeileold[j] = '1';
					l++;
				} else 
					if (zeileold[j] == '1') fprintf(lsfile,"0");
			}		
		} else {
			fprintf(lsfile,"%s",zeile);
			strcpy(zeileold,zeile);
		}
		fprintf(lsfile,"\n");
		if (iterat==N-1) fprintf(lsfile,"\nLarge set found!!!\n");
		fclose(lsfile);
		fflush(stdout);
		i = system("cat ls_tmp >>solutions.ls"); 
		
		if (iterat==N-1) {
			successful = 1;
			goto ausstieg;
		}
		n = n_new;

	}  
}

ausstieg:

	/* successful exit */
	if (successful == 1) {
		printf("\nLS[%d](%ld,%ld,%ld) lambda=%ld:\n",N,t,k,v,lambda);
		fflush(stdout);
		i = system("cat solutions.ls");
	}
	exit(0);
}
