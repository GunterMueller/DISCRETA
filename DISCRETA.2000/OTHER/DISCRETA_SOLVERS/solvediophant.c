/*5:*/
#line 92 "solvediophant.w"

/*6:*/
#line 116 "solvediophant.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include "diophant.h"

/*:6*/
#line 93 "solvediophant.w"
;
/*7:*/
#line 123 "solvediophant.w"

#if defined(MPREC)
verylong matrix_fact= 0;
verylong max_nrm= 0;
verylong*upperbounds;
verylong**A,*b;
#else
long matrix_fact;
long max_nrm;
long*upperbounds;
long**A,*b;
#endif
int bkz_beta,bkz_p;
int system_rows,system_columns;

long stop_after_solutions;
long stop_after_loops;
int free_RHS;
FILE*txt;
char*inputfile_name,*rowp;
#define zlength 16000
char zeile[zlength];
char detectstring[100];

int*original_columns;
int no_original_columns;

/*:7*/
#line 94 "solvediophant.w"
;
int main(int argc,char*argv[]){
int i,j,flag;
int silent;
int iterate;
int iterate_no;

/*8:*/
#line 150 "solvediophant.w"

if(argc<5){
printf("Wrong number of parameters in command line input!\n");
printf("solvediophant factor max_norm bkz beta p filename\n");
printf("solvediophant 1000 6 bkz 120 14 KM.in\n");
printf(" or \n");
printf("solvediophant factor max_norm iterate number filename\n");
printf("solvediophant 1000 6 iterate 10 KM.in\n");
fflush(stdout);
return 1;
}
#if 0
for(i= 1;i<=4;i++){
if((atoi(argv[i])<=0)){
printf("solvediophant 1000 6 120 14 KM.in\n");
fflush(stdout);
return 1;
}
}
#endif 
#if defined(MPREC)
zstrtoz(argv[1],&matrix_fact);
zstrtoz(argv[2],&max_nrm);
#else
matrix_fact= atol(argv[1]);
max_nrm= atol(argv[2]);
#endif

iterate= silent= 0;
for(i= 3;i<argc;i++){
if(strcmp(argv[i],"silent")==0){
silent= 1;
printf("%s\n",argv[i]);
}
if(strcmp(argv[i],"iterate")==0){
iterate= 1;
printf("%s\n",argv[i]);
iterate_no= atoi(argv[i+1]);
}
if(strcmp(argv[i],"bkz")==0){
iterate= 0;
printf("%s\n",argv[i]);
bkz_beta= atoi(argv[i+1]);
bkz_p= atoi(argv[i+2]);
}
}
inputfile_name= argv[argc-1];

/*:8*/
#line 101 "solvediophant.w"
;
/*9:*/
#line 200 "solvediophant.w"

txt= fopen(inputfile_name,"r");
flag= 0;
free_RHS= 0;
stop_after_loops= 0;
stop_after_solutions= 0;
do{
fgets(zeile,zlength,txt);
if(strstr(zeile,"% stopafter")!=NULL){
sscanf(zeile,"%% stopafter %ld",&stop_after_solutions);
}
if(strstr(zeile,"% stoploops")!=NULL){
sscanf(zeile,"%% stoploops %ld",&stop_after_loops);
}
if(strstr(zeile,"% FREERHS")!=NULL){
free_RHS= 1;
}
}
while(zeile[0]=='%');
sscanf(zeile,"%d%d%d",&system_rows,&system_columns,&flag);

/*:9*/
#line 102 "solvediophant.w"
;
/*10:*/
#line 221 "solvediophant.w"

#if defined(MPREC)
A= (verylong**)calloc(system_rows,sizeof(verylong*));
for(j= 0;j<system_rows;j++){
A[j]= (verylong*)calloc(system_columns,sizeof(verylong));
for(i= 0;i<system_columns;i++)A[j][i]= 0;
for(i= 0;i<system_columns;i++)zzero(&(A[j][i]));
}
b= (verylong*)calloc(system_rows,sizeof(verylong));
for(i= 0;i<system_rows;i++)b[i]= 0;
for(i= 0;i<system_rows;i++)zzero(&(b[i]));
#else
A= (long**)calloc(system_rows,sizeof(long*));
for(j= 0;j<system_rows;j++){
A[j]= (long*)calloc(system_columns,sizeof(long));
for(i= 0;i<system_columns;i++)A[j][i]= 0;
}
b= (long*)calloc(system_rows,sizeof(long));
for(i= 0;i<system_rows;i++)b[i]= 0;
#endif  

/*:10*/
#line 103 "solvediophant.w"
;
/*11:*/
#line 242 "solvediophant.w"

#if defined(MPREC)
for(j= 0;j<system_rows;j++){
for(i= 0;i<system_columns;i++)
zfread(txt,&(A[j][i]));
zfread(txt,&(b[j]));
}
#else
for(j= 0;j<system_rows;j++){
for(i= 0;i<system_columns;i++)
fscanf(txt,"%ld",&(A[j][i]));
fscanf(txt,"%ld",&(b[j]));
}
#endif
/*:11*/
#line 104 "solvediophant.w"
;
/*12:*/
#line 258 "solvediophant.w"

sprintf(detectstring,"BOUNDS");
do{
rowp= fgets(zeile,zlength,txt);
}while((rowp!=NULL)&&(strstr(zeile,detectstring)==NULL));

if(rowp==NULL){
upperbounds= NULL;
printf("No %s \n",detectstring);fflush(stdout);
}else{
#if defined(MPREC)
upperbounds= (verylong*)calloc(system_columns,sizeof(verylong));
for(i= 0;i<system_columns;i++){
upperbounds[i]= 0;
zfread(txt,&(upperbounds[i]));
}
#else
upperbounds= (long*)calloc(system_columns,sizeof(long));
for(i= 0;i<system_columns;i++)
fscanf(txt,"%ld",&(upperbounds[i]));
#endif
}
fclose(txt);
txt= fopen(inputfile_name,"r");

/*:12*/
#line 105 "solvediophant.w"
;
/*13:*/
#line 284 "solvediophant.w"

sprintf(detectstring,"SELECTEDCOLUMNS");
do{
rowp= fgets(zeile,zlength,txt);
}while((rowp!=NULL)&&(strstr(zeile,detectstring)==NULL));

if(rowp!=NULL)fscanf(txt,"%d",&(no_original_columns));
else no_original_columns= system_columns;

original_columns= (int*)calloc(no_original_columns,sizeof(int));

if(rowp!=NULL){
for(i= 0;i<no_original_columns;i++)fscanf(txt,"%d",&(original_columns[i]));
}else{
for(i= 0;i<no_original_columns;i++)original_columns[i]= 1;
}
fclose(txt);

/*:13*/
#line 106 "solvediophant.w"
;
printf("Input finished\n");fflush(stdout);

diophant(A,b,upperbounds,system_columns,system_rows,
matrix_fact,max_nrm,silent,iterate,iterate_no,bkz_beta,bkz_p,
stop_after_solutions,stop_after_loops,
free_RHS,original_columns,no_original_columns);

return 1;
}
/*:5*/
