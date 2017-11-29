#define DEEPINSERT 1
#define DEEPINSERT_CONST 2000
#define VERBOSE 0 \

#define GIVENS 1
#define LASTLINESFACTOR 100
#define EPSILON 0.000001
#define LLLCONST_LOW 0.75
#define LLLCONST_HIGH 0.99
#define LLLCONST_HIGHER 0.999 \

#define SQRT sqrt
#define DOUBLE double
#define COEFF struct coe \

#define round(r) ceil(r-0.5)  \

#define FINCKEPOHST 1 \

/*2:*/
#line 43 "diophant.w"

/*5:*/
#line 127 "diophant.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <malloc.h> 
#include <math.h> 
#if defined(MPREC)
#include "freelip/lip.h"
#undef BLAS
#endif
#if defined(BLAS)
#include "blas/blas1.h"
#else
#if defined(BLAS2)
#include "blas_intel/cblas.h"
#endif
#endif

/*:5*/
#line 44 "diophant.w"
;
/*6:*/
#line 146 "diophant.w"

#if defined(MPREC)
struct coe{
verylong c;
int p;
};
#else
struct coe{
long c;
int p;
};
#endif

/*:6*/
#line 45 "diophant.w"
;
/*7:*/
#line 163 "diophant.w"

#if defined(MPREC)
verylong matrix_factor= 0;
verylong max_norm= 0;
verylong max_norm_initial= 0;
verylong max_up= 0;
verylong dummy= 0;
#else
long matrix_factor;
long max_norm,max_norm_initial,max_up;
#endif
long nom,denom;
long lastlines_factor;

/*:7*//*8:*/
#line 178 "diophant.w"

int system_rows,system_columns;
int lattice_rows,lattice_columns;
COEFF**lattice;
int free_RHS;
int iszeroone;
#if defined(MPREC)
verylong*upperbounds;
verylong upperbounds_max= 0;
verylong upfac= 0;
verylong upfac_h= 0;
#else
long*upperbounds;
long upperbounds_max;
long upfac;
#endif

/*:8*//*9:*/
#line 196 "diophant.w"

int*original_columns;
int no_original_columns;
long stop_after_solutions;
long stop_after_loops;
long nosolutions;
int iterate;
int no_iterates;
int bkz_beta,bkz_p;
int SILENT;

/*:9*/
#line 46 "diophant.w"
;
/*10:*/
#line 211 "diophant.w"

#if defined(MPREC)
#define put_to(i,j,val) zcopy(val,&(lattice[i][j+1].c))
#define smult_lattice(i,j,factor) zmulin(factor,&(lattice[i][j+1].c))
#else
#define put_to(i,j,val) lattice[i][j+1].c =  val
#define smult_lattice(i,j,factor) lattice[i][j+1].c *=  factor
#define sdiv_lattice(i,j,factor) lattice[i][j+1].c /=  factor
#endif
#define get_entry(i,j) lattice[i][j+1].c

/*:10*/
#line 47 "diophant.w"
;

/*27:*/
#line 531 "diophant.w"

/*28:*/
#line 539 "diophant.w"

void debug_print(char*m,int l){
if(VERBOSE>=l){
printf("debug>> %s\n",m);fflush(stdout);
}
return;
}

/*:28*/
#line 532 "diophant.w"
;
/*29:*/
#line 549 "diophant.w"

void print_lattice(){
int i,j;
for(i= 0;i<lattice_columns;i++){
for(j= 0;j<lattice_rows;j++){
#if defined(MPREC)
zwrite(get_entry(i,j));
printf(" ");
#else
printf("%ld ",get_entry(i,j));
#endif   
}
printf("\n");
}
printf("\n");fflush(stdout);
return;
}
/*:29*/
#line 533 "diophant.w"
;
/*30:*/
#line 568 "diophant.w"

long gcd(long n1,long n2){
long a,b,c;

if(n1> n2){
a= n1;b= n2;
}else{
a= n2;b= n1;
}

while((c= a%b)> 0){
a= b;b= c;
}
return b;
}
/*:30*/
#line 534 "diophant.w"
;
/*31:*/
#line 584 "diophant.w"

void coeffinit(COEFF*v,int z)
{
short r= 0;
short i;
for(i= z;i>=0;i--){
v[i].p= r;
#if defined(MPREC)
if(!ziszero(v[i].c))r= i;
#else  
if(v[i].c!=0)r= i;
#endif  
}
return;
}

/*:31*/
#line 535 "diophant.w"
;
/*32:*/
#line 605 "diophant.w"

int cutlattice(){
int j,i,flag;

/*33:*/
#line 623 "diophant.w"

j= 0;
do{
if(lattice[j][0].p> system_rows)
j++;
else{
for(i= j+1;i<lattice_columns;i++)lattice[i-1]= lattice[i];
lattice_columns--;
}
}while(j<lattice_columns-1);

/*:33*/
#line 609 "diophant.w"
;
/*34:*/
#line 637 "diophant.w"

flag= 0;
#if defined(MPREC)
for(i= 0;i<lattice_columns;i++)if(!ziszero(get_entry(i,lattice_rows-1))){
#else
for(i= 0;i<lattice_columns;i++)if(get_entry(i,lattice_rows-1)!=0){
#endif
flag= 1;
break;
}
if(flag==0){
printf("Nonhomogenous solution not possible.\n");fflush(stdout);
return 0;
}
/*:34*/
#line 610 "diophant.w"
;

for(j= 0;j<lattice_columns;j++)
for(i= system_rows;i<lattice_rows;i++)
put_to(j,i-system_rows,get_entry(j,i));
lattice_rows-= system_rows;

for(j= 0;j<lattice_columns;j++)coeffinit(lattice[j],lattice_rows);

return 1;
}

/*:32*/
#line 536 "diophant.w"
;
/*35:*/
#line 652 "diophant.w"

int solutiontest(int position){
int i,j;
int low,up;

#if defined(MPREC)
static verylong u= 0;
static verylong s= 0;
#else
long u,s;
#endif

/*36:*/
#line 681 "diophant.w"

#if defined(MPREC)
zcopy(get_entry(position,lattice_rows-1),&dummy);
zabs(&dummy);
if(zcompare(dummy,max_norm)!=0)return 0;
if(ziszero(get_entry(position,lattice_rows-1-free_RHS)))return 0;
#else
if(labs(get_entry(position,lattice_rows-1))!=max_norm)return 0;
if(get_entry(position,lattice_rows-1-free_RHS)==0)return 0;
#endif
/*:36*/
#line 664 "diophant.w"
;
/*37:*/
#line 698 "diophant.w"

low= 0;
up= lattice_rows-1-free_RHS;
#if defined(MPREC)
if(lattice_columns==system_columns+2+free_RHS){
for(i= 0;i<system_rows;i++)
if(!ziszero(get_entry(position,i)))return 0;
low= system_rows;
}

if(iszeroone){
for(i= low;i<up;i++){
zcopy(get_entry(position,i),&dummy);
zabs(&dummy);
if(zcompare(dummy,max_norm)!=0)return 0;
}
}else{
for(i= low;i<up;i++){
zcopy(get_entry(position,i),&dummy);
zabs(&dummy);
if(zcompare(dummy,max_norm)==1)return 0;
}
}
#else
if(lattice_columns==system_columns+2+free_RHS){
for(i= 0;i<system_rows;i++)
if(get_entry(position,i)!=0)return 0;
low= system_rows;
}

if(iszeroone){
for(i= low;i<up;i++)
if(labs(get_entry(position,i))!=max_norm)return 0;
}else{
for(i= low;i<up;i++)
if(labs(get_entry(position,i))> max_norm)return 0;
}
#endif
/*:37*/
#line 665 "diophant.w"
;

#if defined(MPREC)
zone(&upfac);
zsdiv(get_entry(position,lattice_rows-1),lastlines_factor,&s);
#else
upfac= 1;
s= get_entry(position,lattice_rows-1)/lastlines_factor;
#endif 

/*38:*/
#line 737 "diophant.w"

i= low;
for(j= 0;j<no_original_columns;j++){
#if defined(MPREC)
if(original_columns[j]==0){
zzero(&u);
}else{
if(!iszeroone)
zdiv(upperbounds_max,upperbounds[i-low],&upfac,&upfac_h);
zsub(get_entry(position,i),s,&u);
zdiv(u,max_norm_initial,&u,&dummy);
zdiv(u,upfac,&u,&dummy);
zsdiv(u,denom,&u);
zabs(&u);
i++;
}
zwrite(u);
printf(" ");
#else
if(original_columns[j]==0){
u= 0;
}else{
if(!iszeroone)
upfac= upperbounds_max/upperbounds[i-low];
u= labs(get_entry(position,i)-s)/(denom*max_norm_initial*upfac);
i++;
}
printf("%ld ",u);
#endif 
}

if(free_RHS){
#if defined(MPREC)
zdiv(get_entry(position,up),max_up,&u,&dummy);
zsdiv(u,lastlines_factor,&u);
zabs(&u);
printf(" L = ");
zwrite(u);
#else
u= labs(get_entry(position,up))/(lastlines_factor*max_up);
printf(" L = %ld",u);
#endif 
}
printf("\n");fflush(stdout);

/*:38*/
#line 675 "diophant.w"
;
/*39:*/
#line 783 "diophant.w"

if(stop_after_solutions==1){
printf("Stopped in phase 1 after finding a random solution\n");
return(1);
}
/*:39*/
#line 676 "diophant.w"
;

return 1;
}
/*:35*/
#line 537 "diophant.w"
;

/*:27*/
#line 49 "diophant.w"
;
/*40:*/
#line 792 "diophant.w"

/*58:*/
#line 1176 "diophant.w"

/*59:*/
#line 1191 "diophant.w"

DOUBLE scalarproductlfp(COEFF*v,COEFF*w)
{
DOUBLE erg;
long t1,t2;
COEFF*vv,*ww;

erg= 0.0;
t1= v[0].p;
t2= w[0].p;
if((t1==0)||(t2==0))return 0;

do{
if(t2> t1){
t1= v[t2-1].p;
if(t2!=t1){
if(t1==0)break;
t2= w[t2].p;
if(t2==0)break;
}
else goto gleich;
}
else if(t2<t1){
t2= w[t1-1].p;
if(t2!=t1){
if(t2==0)break;
t1= v[t1].p;
if(t1==0)break;
}
else goto gleich;
}
else{
gleich:vv= &(v[t1]);
ww= &(w[t2]);
#if defined(MPREC)   
erg+= (DOUBLE)zdoub(vv->c)*(DOUBLE)zdoub(ww->c);
#else
erg+= (DOUBLE)(vv->c)*(DOUBLE)(ww->c);
#endif   
t1= vv->p;
if(t1==0)break;
t2= ww->p;
if(t2==0)break;
}
}
while(1);

return(erg);
}
/*:59*/
#line 1177 "diophant.w"
;
/*60:*/
#line 1240 "diophant.w"

DOUBLE scalarproductfp(DOUBLE*v,DOUBLE*w,int n)
{
DOUBLE r;
#if defined(BLAS)
r= ddot(n,v,1,w,1);
#else
#if defined(BLAS2)
r= cblas_ddot(n,v,1,w,1);
#else
int i;

r= 0.0;
for(i= n-1;i>=0;i--)r+= v[i]*w[i];
#endif 
#endif 
return r;
}

/*:60*/
#line 1178 "diophant.w"
;
/*61:*/
#line 1261 "diophant.w"

#if !defined(MPREC)
DOUBLE normfp(COEFF*v){
DOUBLE erg;
long t;
COEFF*vv;

erg= 0.0;

t= v[0].p;
if(t==0)return 0;

do{
vv= &(v[t]);
erg+= (DOUBLE)(vv->c)*(DOUBLE)(vv->c);
t= vv->p;
}
while(t!=0);

return erg;
}
#endif


/*:61*/
#line 1179 "diophant.w"
;
/*62:*/
#line 1288 "diophant.w"

#if defined(MPREC)     
int lllalloc(DOUBLE***mu,DOUBLE**c,DOUBLE**N,DOUBLE***bs,int s,int z)
#else 
int lllalloc(DOUBLE***mu,DOUBLE**c,DOUBLE**N,int s,int z)
#endif 
{
int i,m;

if((z<1)||(s<1))return 0;

(*c)= (DOUBLE*)calloc(s,sizeof(DOUBLE));
(*N)= (DOUBLE*)calloc(s,sizeof(DOUBLE));
(*mu)= (DOUBLE**)calloc(s,sizeof(DOUBLE*));
for(i= 0;i<s;i++)(*mu)[i]= (DOUBLE*)calloc(z,sizeof(DOUBLE));

m= (z> s)?z:s;
#if defined(MPREC)     
(*bs)= (DOUBLE**)calloc(m,sizeof(DOUBLE*));

for(i= 0;i<m;i++)
(*bs)[i]= (DOUBLE*)calloc(m,sizeof(DOUBLE));
#endif  

return 1;
}
/*:62*/
#line 1180 "diophant.w"
;
/*63:*/
#line 1314 "diophant.w"

#if defined(MPREC)
int lllfree(DOUBLE**mu,DOUBLE*c,DOUBLE*N,DOUBLE**bs,int s)
#else
int lllfree(DOUBLE**mu,DOUBLE*c,DOUBLE*N,int s)
#endif
{
int i;

#if defined(MPREC)
for(i= 0;i<s;++i)free(bs[i]);
free(bs);
#endif 
for(i= 0;i<s;++i)free(mu[i]);
free(mu);
free(N);
free(c);

return 1;
}
/*:63*/
#line 1181 "diophant.w"
;

/*:58*/
#line 793 "diophant.w"
;
/*41:*/
#line 805 "diophant.w"

#define TWOTAUHALF 67108864.0 
#if defined(MPREC)
int lllfp(COEFF**b,DOUBLE**mu,DOUBLE*c,DOUBLE*N,DOUBLE**bs,
int start,int s,int z,DOUBLE delta)
#else
int lllfp(COEFF**b,DOUBLE**mu,DOUBLE*c,DOUBLE*N,
int start,int s,int z,DOUBLE delta)
#endif
{
/*43:*/
#line 864 "diophant.w"

int i,j,k;
DOUBLE ss;
int counter;
/*:43*//*46:*/
#line 919 "diophant.w"

int Fc,Fr;
DOUBLE mus,cc;
#if defined(MPREC)
verylong musvl= 0;
verylong hv= 0;
char musch[100];
DOUBLE*swapd;
#else
long musvl;
#endif

/*:46*//*48:*/
#line 948 "diophant.w"

int ii,iii;
COEFF*swapvl;
COEFF*bb;
/*:48*//*55:*/
#line 1100 "diophant.w"

int Fi;

/*:55*/
#line 815 "diophant.w"
;
if((z<=1)||(s<=1)){
printf("Wrong dimensions in lllfp\n");fflush(stdout);
return(0);
}

k= (start> 1)?start:1;

/*42:*/
#line 850 "diophant.w"

if(k<1)k= 1;
for(i= k-1;i<s;++i){
#if defined(MPREC)
ss= 0.0;
for(j= 0;j<z;++j){
bs[i][j]= (DOUBLE)zdoub(b[i][j+1].c);
ss+= bs[i][j]*bs[i][j];
}
#else   
ss= normfp(b[i]);
#endif
N[i]= SQRT(ss);
}
/*:42*/
#line 823 "diophant.w"
;
counter= 0;
while(k<s){
#if VERBOSE >  3
if((counter%500)==0){printf("LLL: %d k:%d\n",counter,k);fflush(stdout);}
counter++;
#endif

/*44:*/
#line 871 "diophant.w"

if(k==1)c[0]= N[0]*N[0];
c[k]= N[k]*N[k];
for(j= 0;j<k;j++){
#if defined(MPREC)
ss= scalarproductfp(bs[k],bs[j],z);
if(fabs(ss)<N[k]*N[j]/TWOTAUHALF){
ss= (DOUBLE)scalarproductlfp(b[k],b[j]);
}
#else
ss= scalarproductlfp(b[k],b[j]);
#endif 
for(i= 0;i<j;i++)ss-= mu[j][i]*mu[k][i]*c[i];
mu[k][j]= ss/c[j];
c[k]-= ss*mu[k][j];
}
/*:44*/
#line 831 "diophant.w"
;
/*45:*/
#line 901 "diophant.w"

Fc= Fr= 0;
for(j= k-1;j>=0;j--){
if(fabs(mu[k][j])> 0.5){
/*47:*/
#line 934 "diophant.w"

mus= round(mu[k][j]);
#if defined(MPREC)
sprintf(musch,"%.0f",mus);
zsread(musch,&musvl);
#else
musvl= (long)mus;
#endif    
if(fabs(mus)> TWOTAUHALF){
Fc= 1;
printf("correct possible rounding errors\n");fflush(stdout);
}


/*:47*/
#line 905 "diophant.w"
;
Fr= 1;
/*49:*/
#line 961 "diophant.w"

#if defined(MPREC)    
switch(ztoint(musvl)){
#else
switch(musvl){
#endif    
case 1:
/*50:*/
#line 980 "diophant.w"

i= b[j][0].p;
while(i!=0){
bb= &(b[k][i]);
#if defined(MPREC)    
zsub(bb->c,b[j][i].c,&(bb->c));
iii= bb->p;
if((b[k][i-1].p!=i)&&(ziszero(bb->c)!=1))
for(ii= i-1;(b[k][ii].p==iii)&&(ii>=0);ii--)b[k][ii].p= i;
else if(ziszero(bb->c)==1){
#else
bb->c-= b[j][i].c;
iii= bb->p;
if((b[k][i-1].p!=i)&&(bb->c!=0))
for(ii= i-1;(b[k][ii].p==iii)&&(ii>=0);ii--)b[k][ii].p= i;
else if(bb->c==0){
#endif
for(ii= i-1;(b[k][ii].p==i)&&(ii>=0);ii--)b[k][ii].p= iii;
}
i= b[j][i].p;
}
for(i= 0;i<j;i++)mu[k][i]-= mu[j][i];
/*:50*/
#line 968 "diophant.w"
;
break;

case-1:
/*51:*/
#line 1004 "diophant.w"

i= b[j][0].p;
while(i!=0){
bb= &(b[k][i]);
#if defined(MPREC)    
zadd(bb->c,b[j][i].c,&(bb->c));
iii= bb->p;
if((b[k][i-1].p!=i)&&(ziszero(bb->c)!=1))
for(ii= i-1;(b[k][ii].p==iii)&&(ii>=0);ii--)b[k][ii].p= i;
else if(ziszero(bb->c)==1){
#else
bb->c+= b[j][i].c;
iii= bb->p;
if((b[k][i-1].p!=i)&&(bb->c!=0))
for(ii= i-1;(b[k][ii].p==iii)&&(ii>=0);ii--)b[k][ii].p= i;
else if(bb->c==0){
#endif
for(ii= i-1;(b[k][ii].p==i)&&(ii>=0);ii--)b[k][ii].p= iii;
}
i= b[j][i].p;
}
for(i= 0;i<j;i++)mu[k][i]+= mu[j][i];
/*:51*/
#line 972 "diophant.w"
;
break;

default:
/*52:*/
#line 1028 "diophant.w"

i= b[j][0].p;
while(i!=0){
bb= &(b[k][i]);
#if defined(MPREC)    
zmul(b[j][i].c,musvl,&hv);
zsub(bb->c,hv,&(bb->c));
iii= bb->p;
if((b[k][i-1].p!=i)&&(ziszero(bb->c)!=1))
for(ii= i-1;(b[k][ii].p==iii)&&(ii>=0);ii--)b[k][ii].p= i;
else if(ziszero(bb->c)==1){
#else
bb->c-= b[j][i].c*musvl;
iii= bb->p;
if((b[k][i-1].p!=i)&&(bb->c!=0))
for(ii= i-1;(b[k][ii].p==iii)&&(ii>=0);ii--)b[k][ii].p= i;
else if(bb->c==0){
#endif
for(ii= i-1;(b[k][ii].p==i)&&(ii>=0);ii--)b[k][ii].p= iii;
}
i= b[j][i].p;
}
#if 0
daxpy(j,-mus,mu[k],1,mu[j],1);
#endif
for(i= 0;i<j;i++)mu[k][i]-= mu[j][i]*mus;
/*:52*/
#line 976 "diophant.w"
;
}
/*:49*/
#line 907 "diophant.w"
;
mu[k][j]-= mus;
solutiontest(k);
}
}
/*53:*/
#line 1054 "diophant.w"

{
#if defined(MPREC)
N[k]= 0.0;
for(i= 0;i<z;i++){
bs[k][i]= (DOUBLE)zdoub(b[k][i+1].c);
N[k]+= bs[k][i]*bs[k][i];
}
#else
N[k]= normfp(b[k]);
#endif
N[k]= SQRT(N[k]);
}
/*:53*/
#line 912 "diophant.w"
;
if(Fc==1){
k= (k-1> 1)?k-1:1;
}else{
/*54:*/
#line 1072 "diophant.w"

if(N[k]<-EPSILON){
printf("Nk negativ! contact the author.\n");fflush(stdout);
exit(0);
}
if(N[k]<0.5){
swapvl= b[k];
ss= N[k];
#if defined(MPREC)     
swapd= bs[k];
#endif     
for(i= k+1;i<s;i++){
b[i-1]= b[i];
N[i-1]= N[i];
#if defined(MPREC)     
bs[i-1]= bs[i];
#endif     
}
b[s-1]= swapvl;
N[s-1]= ss;
#if defined(MPREC)     
bs[s-1]= swapd;
#endif  
s= s-1;
k= 1;
continue;
}

/*:54*/
#line 916 "diophant.w"
;
}

/*:45*/
#line 832 "diophant.w"
;
#if defined(DEEPINSERT)
/*56:*/
#line 1111 "diophant.w"

cc= N[k]*N[k];
j= 0;
Fi= 0;
while(j<k){
#if 0
if((j> DEEPINSERT_CONST&&j<k-DEEPINSERT_CONST)||delta*c[j]<=cc){
#endif    
if(delta*c[j]<=cc){
cc-= mu[k][j]*mu[k][j]*c[j];
j++;
}else{
swapvl= b[k];
ss= N[k];
#if defined(MPREC)     
swapd= bs[k];
#endif     
for(i= k-1;i>=j;i--){
b[i+1]= b[i];
N[i+1]= N[i];
#if defined(MPREC)     
bs[i+1]= bs[i];
#endif     
}
b[j]= swapvl;
N[j]= ss;
#if defined(MPREC)     
bs[j]= swapd;
#endif     

Fi= 1;
break;
}
}
if(Fi==1)k= (j-1> 1)?j-1:1;
else{
k++;
}

/*:56*/
#line 834 "diophant.w"
;
#else  
/*57:*/
#line 1157 "diophant.w"

if(delta*c[k-1]> c[k]+mu[k][k-1]*mu[k][k-1]*c[k-1]){
swapvl= b[k];
b[k]= b[k-1];
b[k-1]= swapvl;
ss= N[k];
N[k]= N[k-1];
N[k-1]= ss;
#if defined(MPREC)     
swapd= bs[k];
bs[k]= bs[k-1];
bs[k-1]= swapd;
#endif    

k= (k-1> 1)?k-1:1;
}else k++;

/*:57*/
#line 836 "diophant.w"
;
#endif  
}
return(1);

}
/*:41*/
#line 794 "diophant.w"
;
/*64:*/
#line 1337 "diophant.w"

double orthogonal_defect(COEFF**lattice,DOUBLE*c,int s,int z){
double defect= 0.0;

#if !defined(MPREC)
int i;
for(i= 0;i<s;i++)defect+= log((double)normfp(lattice[i]))
-log((double)c[i]);
#endif 
defect/= 2.0;
return defect;
}

/*:64*/
#line 795 "diophant.w"
;
/*65:*/
#line 1351 "diophant.w"

void lll(COEFF**b,int s,int z,DOUBLE quality)
{
DOUBLE**mu;
DOUBLE*c;
DOUBLE*N;
#if defined(MPREC) 
DOUBLE**bs;
#endif 
int r;

#if defined(MPREC) 
lllalloc(&mu,&c,&N,&bs,s,z);
r= lllfp(b,mu,c,N,bs,1,s,z,quality);
lllfree(mu,c,N,bs,s);
#else
lllalloc(&mu,&c,&N,s,z);
r= lllfp(b,mu,c,N,1,s,z,quality);
printf("Orthogonal defect: %f\n",orthogonal_defect(b,c,s,z));
fflush(stdout);
lllfree(mu,c,N,s);
#endif 

return;
}

/*:65*/
#line 796 "diophant.w"
;
/*66:*/
#line 1381 "diophant.w"

void iteratedlll(COEFF**b,int s,int z,int no_iterates,DOUBLE quality)
{
DOUBLE**mu;
DOUBLE*c;
DOUBLE*N;
#if defined(MPREC) 
DOUBLE**bs;
#endif 
int r,l,i,j,runs;
COEFF*swapvl;

#if defined(MPREC) 
lllalloc(&mu,&c,&N,&bs,s,z);
r= lllfp(b,mu,c,N,bs,1,s,z,quality);
lllfree(mu,c,N,bs,s);
#else
lllalloc(&mu,&c,&N,s,z);
r= lllfp(b,mu,c,N,1,s,z,quality);
#endif

printf("   Orthogonal defect: %f\n",orthogonal_defect(b,c,s,z));
fflush(stdout);


for(runs= 1;runs<no_iterates;runs++){


for(j= s-1;j> 0;j--){
for(l= j-1;l>=0;l--){
if(N[l]<N[j]){
swapvl= b[l];
for(i= l+1;i<=j;i++)b[i-1]= b[i];
b[j]= swapvl;
}
}
}

#if defined(MPREC) 
r= lllfp(b,mu,c,N,bs,1,s,z,quality);
#else
r= lllfp(b,mu,c,N,1,s,z,quality);
#endif  
printf("%d: Orthogonal defect: %f\n",runs,orthogonal_defect(b,c,s,z));
fflush(stdout);
}

#if defined(MPREC) 
lllfree(mu,c,N,bs,s);
#else
lllfree(mu,c,N,s);
#endif 

return;
}

/*:66*/
#line 797 "diophant.w"
;
/*67:*/
#line 1438 "diophant.w"

/*78:*/
#line 1797 "diophant.w"

DOUBLE laurin(DOUBLE x){
static DOUBLE K1= 0.9181938533204672741780329736405620;
static DOUBLE K2= 0.0833333333333333333333333333333333333;
static DOUBLE K3= 0.0027777777777777777777777777777777777;
static DOUBLE K4= 0.000793650793650793650793650793650793650;
static DOUBLE K5= 0.0005952380952380952380952380952380952380;
static DOUBLE K6= 0.000841750841750841750841750841750841750;
static DOUBLE K7= 0.001917526917526917526917526917526917526;
static DOUBLE K8= 0.00641025641025641025641025641025641025;
static DOUBLE K9= 0.0295506529510021209716796875;
static DOUBLE K10= 0.17968122661113739013671875;
static DOUBLE K11= 1.39243221282958984375;
DOUBLE y;

y= 1.0/(x*x);
y= (x-0.5)*log(x)-x+K1+
(1.0/x)*(K2-y*(K3-y*(K4-y*(K5-y*(K6-y*(K7-y*(K8-y*(K9-y*(K10-y*K11)))))))));
return y;
}

DOUBLE log_gamma(DOUBLE x){
DOUBLE y;
int i,n;
static int MM= 13;

if(x<=0.0)return-1.0;

if(x> 100000000.0){
y= x*(log(x)-1.0);
}else{
if(x>=MM){
y= laurin(x);
}else{
n= MM-(int)(floor(x));
y= x-floor(x)+MM;
y= laurin(y);
for(i= 0;i<n;i++)y-= log(y+i);
}
}
return y;
}


/*:78*/
#line 1439 "diophant.w"
;
/*75:*/
#line 1674 "diophant.w"

DOUBLE enumerate(DOUBLE**mu,DOUBLE*c,long*u,int s,
int start_block,int end_block,int p)
{
DOUBLE cd,dum;
DOUBLE*y,*cs,*eta;
long*us,*delta,*d,*v;
int t,i,t_up;
int tmax;
static DOUBLE pi= 3.141592653589793238462643383;
int SCHNITT= 10;

if(c[start_block]<=EPSILON){
printf("Hier ist was faul! start_block=%d %f\n",start_block,
(double)c[start_block]);fflush(stdout);exit(0);
}

us= (long*)calloc(s+1,sizeof(long));
cs= (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
y= (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
delta= (long*)calloc(s+1,sizeof(long));
d= (long*)calloc(s+1,sizeof(long));
eta= (DOUBLE*)calloc(s+1,sizeof(DOUBLE));
v= (long*)calloc(s+1,sizeof(long));

for(i= start_block;i<=end_block+1;i++){
cs[i]= y[i]= 0.0;
u[i]= us[i]= v[i]= delta[i]= 0;
d[i]= 1;
}
us[start_block]= u[start_block]= 1;
cd= c[start_block];

t= tmax= start_block;


/*77:*/
#line 1784 "diophant.w"

eta[start_block]= 0.0;
if(end_block-start_block<SCHNITT){
for(i= start_block+1;i<=end_block;i++)eta[i]= 0.0;
}
else{
dum= log(c[start_block]);
for(i= start_block+1;i<=end_block;i++){
t_up= i-start_block;
eta[i]= exp((log_gamma(t_up/2.0+1)-p*log(2.0))*2.0/t_up+dum/t_up)/pi;
if(i<end_block)dum+= log(c[i]);
}
}
/*:77*/
#line 1710 "diophant.w"
;
while(t<=end_block){
/*76:*/
#line 1726 "diophant.w"

{
dum= us[t]+y[t];
cs[t]= cs[t+1]+dum*dum*c[t];

if(cs[t]<cd-eta[t]){
if(t> start_block){
t--;
delta[t]= 0;
dum= 0.0;
for(i= t+1;i<=tmax;i++)dum+= us[i]*mu[i][t];
y[t]= dum;
us[t]= v[t]= (long)(round(-y[t]));
d[t]= (v[t]> -y[t])?-1:1;
}else{
cd= cs[start_block];
for(i= start_block;i<=end_block;i++)u[i]= us[i];
goto nextstep;
}
}else{
t++;
nextstep:
if(tmax<t)tmax= t;
if(t<tmax)delta[t]= -delta[t];
if(delta[t]*d[t]>=0)delta[t]+= d[t];
us[t]= v[t]+delta[t];
}
}

/*:76*/
#line 1712 "diophant.w"
;
}

free(us);
free(cs);
free(y);
free(delta);
free(d);
free(eta);
free(v);
return(cd);
}

/*:75*/
#line 1440 "diophant.w"
;
/*68:*/
#line 1459 "diophant.w"

int bkz(COEFF**b,int s,int z,DOUBLE delta,int beta,int p)
{
/*69:*/
#line 1526 "diophant.w"

DOUBLE**mu,*c,*N;
#if defined(MPREC) 
DOUBLE**bs;
static verylong hv= 0;
#endif 
int zaehler;
int h,i,last;
int start_block,end_block;
long*u;
DOUBLE new_cj;


/*:69*//*71:*/
#line 1570 "diophant.w"

int g,ui,q,swapi,j;
COEFF*swapvl;

/*:71*/
#line 1462 "diophant.w"
;

last= s-2;
if(last<2){
printf("BKZ: s too small.\n");
printf("Please increase c0 (the first parameter)!\n");
return 0;
}

u= (long*)calloc(s,sizeof(long));
for(i= 0;i<s;i++)u[i]= 0;

#if defined(MPREC)
lllalloc(&mu,&c,&N,&bs,s,z);
#else
lllalloc(&mu,&c,&N,s,z);
#endif 


#if defined(MPREC)
lllfp(b,mu,c,N,bs,1,s,z,delta);
#else 
lllfp(b,mu,c,N,1,s,z,delta);
#endif 

start_block= zaehler= -1;
while(zaehler<last){

start_block++;
if(start_block==last)start_block= 0;
end_block= (start_block+beta-1<last)?start_block+beta-1:last;


new_cj= enumerate(mu,c,u,s,start_block,end_block,p);

h= (end_block+1<last)?end_block+1:last;


if(delta*c[start_block]> new_cj){
/*70:*/
#line 1549 "diophant.w"

#if defined(ORIGINAL_SCHNORR_EUCHNER)
/*73:*/
#line 1639 "diophant.w"

#if defined(MPREC)
for(l= 1;l<=z;l++)zzero(&(b[last+1][l].c));
for(i= start_block;i<=end_block;i++)
for(l= 1;l<=z;l++){
zsmul(b[i][l].c,u[i],&hv);
zadd(b[last+1][l].c,hv,&(b[last+1][l].c));
}
#else
for(l= 1;l<=z;l++)b[last+1][l].c= 0;
for(i= start_block;i<=end_block;i++)
for(l= 1;l<=z;l++)b[last+1][l].c+= b[i][l].c*u[i];
#endif  
coeffinit(b[last+1],z);
solutiontest(last+1);
swapvl= b[last+1];
for(i= last;i>=start_block;i--)b[i+1]= b[i];
b[start_block]= swapvl;

/*:73*/
#line 1551 "diophant.w"
;
#else
/*72:*/
#line 1581 "diophant.w"

#if defined(MPREC)
for(j= 1;j<=z;j++)zzero(&(b[last+1][j].c));
for(i= start_block;i<=end_block;i++){
if(u[i]!=0)for(j= 1;j<=z;j++){
zsmul(b[i][j].c,u[i],&hv);
zadd(b[last+1][j].c,hv,&(b[last+1][j].c));
}
}
#else
for(j= 1;j<=z;j++)b[last+1][j].c= 0;
for(i= start_block;i<=end_block;i++){
if(u[i]!=0)for(j= 1;j<=z;j++)
b[last+1][j].c+= u[i]*b[i][j].c;
}
#endif 
g= end_block;
while(u[g]==0)g--;

i= g-1;
while(labs(u[g])> 1){
while(u[i]==0)i--;
q= (int)round((1.0*u[g])/u[i]);
ui= u[i];
u[i]= u[g]-q*u[i];
u[g]= ui;

#if defined(MPREC)
for(j= 1;j<=z;j++){
zcopy(b[g][j].c,&hv);
zsmul(b[g][j].c,(long)q,&(b[g][j].c));
zadd(b[g][j].c,b[i][j].c,&(b[g][j].c));
zcopy(hv,&(b[i][j].c));
}
#else
for(j= 1;j<=z;j++){
swapi= b[g][j].c;
b[g][j].c= q*b[g][j].c+b[i][j].c;
b[i][j].c= swapi;
}
#endif  
coeffinit(b[g],z);
coeffinit(b[i],z);
}

swapvl= b[g];
for(i= g;i> start_block;i--)b[i]= b[i-1];
b[start_block]= b[last+1];
coeffinit(b[start_block],z);

b[last+1]= swapvl;
for(j= 1;j<=z;j++)b[last+1][j].c= 0;
coeffinit(b[last+1],z);

/*:72*/
#line 1553 "diophant.w"
;
#endif

#if defined(MPREC)
lllfp(b,mu,c,N,bs,start_block-1,h+1,z,delta);
#else
lllfp(b,mu,c,N,start_block-1,h+1,z,delta);
#endif  

if(N[h]<-EPSILON){
printf("NN negativ\n");fflush(stdout);
exit(0);
}
#if defined(ORIGINAL_SCHNORR_EUCHNER)
/*74:*/
#line 1658 "diophant.w"

if(N[h]<0.5){
swapvl= b[h];
for(i= h+1;i<=last+1;i++)b[i-1]= b[i];
b[last+1]= swapvl;
}else{
printf("Not linear dependent; %f\n",(double)(N[h-1]));fflush(stdout);
exit(0);
}

/*:74*/
#line 1567 "diophant.w"
;
#endif

/*:70*/
#line 1501 "diophant.w"
;
zaehler= -1;
}
else{
#if defined(MPREC)
lllfp(b,mu,c,N,bs,h-2,h+1,z,delta);
#else 
lllfp(b,mu,c,N,h-2,h+1,z,delta);
#endif 
zaehler++;
}
}

printf("bkz: Orthogonal defect: %f\n",orthogonal_defect(b,c,s-1,z));
fflush(stdout);

#if defined(MPREC)
lllfree(mu,c,N,bs,s);
#else
lllfree(mu,c,N,s);
#endif 
free(u);

return 1;
}
/*:68*/
#line 1441 "diophant.w"
;

/*:67*/
#line 798 "diophant.w"
;
/*79:*/
#line 1845 "diophant.w"

/*80:*/
#line 1909 "diophant.w"

static FILE*fp;


/*:80*//*83:*/
#line 1961 "diophant.w"

#if defined(FINCKEPOHST)
DOUBLE**muinv;
#endif 
int*upb,*lowb;
long fipo_success;

/*:83*//*85:*/
#line 1977 "diophant.w"

#if defined(MPREC)
verylong dummy1= 0;
verylong dummy2= 0;
#endif 
/*:85*//*89:*/
#line 2042 "diophant.w"

DOUBLE dum1;

/*:89*//*92:*/
#line 2098 "diophant.w"

long only_zeros_no,only_zeros_success,hoelder_no,hoelder_success;
long cs_success;
DOUBLE dums;

/*:92*/
#line 1846 "diophant.w"
;
/*101:*/
#line 2258 "diophant.w"

void open_solution_file(){
fp= fopen("solutions","w");
}
void close_solution_file(){
fclose(fp);
}

/*:101*/
#line 1847 "diophant.w"
;
/*102:*/
#line 2266 "diophant.w"

DOUBLE compute_y(DOUBLE**mu,DOUBLE*us,int level,int level_max){
int i;
DOUBLE dum;

i= level_max;
dum= 0.0;
while(i>=level+1){
dum+= mu[i][level]*us[i];
i--;
}

return dum;
}

void compute_w(DOUBLE**w,DOUBLE**bd,DOUBLE dum,int level,int rows){
#if defined(BLAS)
dcopy(rows,w[level+1],1,w[level],1);
daxpy(rows,dum,bd[level],1,w[level],1);
#else
#if defined(BLAS2)
cblas_dcopy(rows,w[level+1],1,w[level],1);
cblas_daxpy(rows,dum,bd[level],1,w[level],1);
#else
int l;

l= rows-1;
do{
w[level][l]= w[level+1][l]+dum*bd[level][l];
l--;
}while(l>=0);
#endif
#endif
return;
}

/*:102*/
#line 1848 "diophant.w"
;
/*103:*/
#line 2304 "diophant.w"

/*104:*/
#line 2309 "diophant.w"

void gramschmidt(COEFF**lattice,int columns,int rows,DOUBLE**mu,
DOUBLE**bd,DOUBLE*c,DOUBLE*N,DOUBLE Fq){
int i,l,j;
DOUBLE dum;

for(i= 0;i<columns;i++){
#if defined(MPREC)
for(l= 0;l<rows;l++)bd[i][l]= zdoub(get_entry(i,l));
#else
for(l= 0;l<rows;l++)bd[i][l]= (DOUBLE)get_entry(i,l);
#endif
N[i]= 0.0;
for(j= 0;j<i;j++){
dum= 0.0;
#if defined(MPREC)
for(l= 0;l<rows;l++)dum+= zdoub(get_entry(i,l))*bd[j][l];
#else
for(l= 0;l<rows;l++)dum+= (DOUBLE)get_entry(i,l)*bd[j][l];
#endif
mu[i][j]= dum/c[j];
for(l= 0;l<rows;l++)bd[i][l]-= mu[i][j]*bd[j][l];
}

c[i]= scalarproductfp(bd[i],bd[i],rows);
for(l= 0;l<rows;l++)N[i]+= fabs(bd[i][l]);
N[i]/= c[i];
N[i]*= Fq;
#if VERBOSE >  0  
printf("%f ",(double)c[i]);
#endif
N[i]*= (1.0+EPSILON);
}
#if VERBOSE >  0  
printf("\n\n");
#endif
return;
}

/*:104*/
#line 2305 "diophant.w"
;
/*105:*/
#line 2349 "diophant.w"

void givens(COEFF**lattice,int columns,int rows,DOUBLE**mu,
DOUBLE**bd,DOUBLE*c,DOUBLE*N,DOUBLE Fq){
int i,l,j;
int mm;
DOUBLE d1,d2;
DOUBLE gc,gs;
DOUBLE t;





for(i= 0;i<columns;i++){
for(l= 0;l<rows;l++){
#if defined(MPREC)
mu[i][l]= zdoub(get_entry(i,l));
#else
mu[i][l]= (DOUBLE)get_entry(i,l);
#endif
}
}

for(i= 0;i<rows;i++){
for(l= 0;l<rows;l++)bd[i][l]= 0.0;
bd[i][i]= 1.0;
}

for(j= 1;j<rows;j++){
mm= (j<columns)?j:columns;
for(i= 0;i<mm;i++){
if(mu[i][j]==0.0){

gc= 1.0;
gs= 0.0;
}else{



if(fabs(mu[i][j])>=fabs(mu[i][i])){
t= mu[i][i]/mu[i][j];
gs= 1.0/SQRT(1.0+t*t);
gc= gs*t;
}else{
t= mu[i][j]/mu[i][i];
gc= 1.0/SQRT(1.0+t*t);
gs= gc*t;
}

for(l= i;l<columns;l++){
d1= mu[l][i];
d2= mu[l][j];
mu[l][i]= gc*d1+gs*d2;
mu[l][j]= -gs*d1+gc*d2;
}

for(l= 0;l<rows;l++){
d1= bd[i][l];
d2= bd[j][l];
bd[i][l]= gc*d1+gs*d2;
bd[j][l]= -gs*d1+gc*d2;
}
}
}
}


for(i= 0;i<columns;i++){
c[i]= mu[i][i]*mu[i][i];
N[i]= 0.0;
for(j= 0;j<rows;j++){
bd[i][j]*= mu[i][i];
N[i]+= fabs(bd[i][j]);
}
N[i]/= c[i];
N[i]*= Fq;
N[i]*= 1.0+EPSILON;

for(j= columns-1;j>=i;j--){
mu[j][i]/= mu[i][i];
}
#if VERBOSE >  0  
printf("%f ",(double)c[i]);
#endif
}
#if VERBOSE >  0
printf("\n\n");fflush(stdout);
#endif

return;
}
/*:105*/
#line 2306 "diophant.w"
;

/*:103*/
#line 1849 "diophant.w"
;
/*106:*/
#line 2441 "diophant.w"

#if defined(FINCKEPOHST)
void inverse(DOUBLE**mu,DOUBLE**muinv,int columns){
int i,j,k;
DOUBLE sum;
for(j= 0;j<columns;j++)
for(i= j;i>=0;i--){
sum= 0.0;
for(k= i+1;k<columns;k++)sum+= mu[k][i]*muinv[k][j];
if(i==j)
muinv[i][j]= 1.0-sum;
else
muinv[i][j]= -sum;
}
return;
}
#endif

/*:106*/
#line 1850 "diophant.w"
;
/*107:*/
#line 2460 "diophant.w"

/*108:*/
#line 2468 "diophant.w"

int exacttest(DOUBLE*v,int rows,DOUBLE Fq){
int i;

i= rows-1;
do{
if(fabs(v[i])> Fq*(1.0+EPSILON))return 0;
i--;
}while(i>=0);
return 1;
}

/*:108*/
#line 2461 "diophant.w"
;
/*109:*/
#line 2481 "diophant.w"

int prune0(DOUBLE li,DOUBLE re){
if(li> re){
return 1;
}else{
return 0;
}
}

/*:109*/
#line 2462 "diophant.w"
;
/*110:*/
#line 2491 "diophant.w"

int prune(DOUBLE*w,DOUBLE cs,int rows,DOUBLE Fq){
DOUBLE reseite;
#if !defined(BLAS)
int i;
#endif 

hoelder_no++;
reseite= -cs*(1.0-EPSILON)/Fq;

#if defined(BLAS)
reseite+= dasum(rows,w,1);
#else
#if defined(BLAS2)
reseite+= cblas_dasum(rows,w,1);
#else
i= rows-1;
do{
reseite+= fabs(w[i]);
i--;
}while(i>=0);
#endif 
#endif 
if(0.0<reseite)return 0;

hoelder_success++;
return 1;
}

/*:110*/
#line 2463 "diophant.w"
;
/*111:*/
#line 2522 "diophant.w"

int prune_only_zeros(DOUBLE*w,int level,int rows,DOUBLE Fq,
int*first_nonzero_in_column,int*firstp){
int i;
int f;

only_zeros_no++;

if(iszeroone){
for(i= 0;i<first_nonzero_in_column[firstp[level]];i++){
f= first_nonzero_in_column[firstp[level]+1+i];
if(fabs(fabs(w[f])-Fq)> EPSILON){
only_zeros_success++;
return 1;
}
}
}else{
for(i= 0;i<first_nonzero_in_column[firstp[level]];i++){
f= first_nonzero_in_column[firstp[level]+1+i];
if(fabs(w[f])> Fq*(1+EPSILON)){
only_zeros_success++;
return 1;
}
}
}
return 0;
}
/*:111*/
#line 2464 "diophant.w"
;
/*112:*/
#line 2550 "diophant.w"

int print_solution(DOUBLE*w,int rows,DOUBLE Fq){
int i,j;
long u,s;
int upper;
#if defined(MPREC)
static verylong dummy= 0;
static verylong upfac= 0;
static verylong upfac_h= 0;
#else
long upfac;
#endif

if(fabs(fabs(w[rows-1])-Fq)> EPSILON)return 0;
upper= rows-1-free_RHS;
if(free_RHS&&fabs(fabs(w[upper])-Fq)> EPSILON)return 0;

nosolutions++;

if(!SILENT){
#if defined(MPREC)
zone(&upfac);
#else
upfac= 1;
#endif
s= (long)round(w[rows-1]);
i= 0;
for(j= 0;j<no_original_columns;j++){
if(original_columns[j]==0){
u= 0;
}else{
#if defined(MPREC)
if(!iszeroone)
zdiv(upperbounds_max,upperbounds[i],&upfac,&upfac_h);
zmul(max_norm_initial,upfac,&dummy);
zsmul(dummy,denom,&dummy);
u= labs((long)round(w[i])-s)/(long)zdoub(dummy);
#else
if(!iszeroone)
upfac= upperbounds_max/upperbounds[i];
u= labs((long)round(w[i])-s)/(denom*max_norm_initial*upfac);
#endif
i++;
}
printf("%ld",u);
fprintf(fp,"%ld",u);
if(!iszeroone){
printf(" ");
fprintf(fp," ");
}
}

if(free_RHS){
#if defined(MPREC)
u= labs((long)round(w[upper]))/(long)zdoub(max_up);
#else
u= labs((long)round(w[upper]))/max_up;
#endif  
printf(" L = %ld",u);
}

printf("\n");fflush(stdout);
fprintf(fp,"\n");fflush(fp);
}
if(nosolutions%10000==0){
printf("%ld\n",nosolutions);fflush(stdout);
}
return 1;
}

/*:112*/
#line 2465 "diophant.w"
;

/*:107*/
#line 1851 "diophant.w"
;

DOUBLE explicit_enumeration(COEFF**lattice,int columns,int rows)
{
int level,level_max;
int i,j,l,m;
long loops;

DOUBLE*y,*cs,*us;

long*delta,*d,*eta,*v;
int*first_nonzero,*first_nonzero_in_column,*firstp;

DOUBLE*N,**mu,*c,**w,**bd;

DOUBLE Fd,Fq;
DOUBLE dum;
COEFF*swap_vec;

#if defined(FINCKEPOHST)
DOUBLE*fipo;
int fiposort;
#endif


/*81:*/
#line 1914 "diophant.w"

printf("Dimension of solution space (k): %d compared to s-z+2: %d\n",
columns,system_columns-system_rows+1+free_RHS);
fflush(stdout);

if(columns<system_columns-system_rows+1+free_RHS){
printf("LLL didn't succeed in computing a basis of the kernel.\n");
printf("Please increase c0 (the first parameter)!\n");
return 0;
}


/*:81*/
#line 1876 "diophant.w"
;
/*82:*/
#line 1928 "diophant.w"

#if defined(MPREC)
lllalloc(&mu,&c,&N,&bd,columns,rows);
#else
lllalloc(&mu,&c,&N,columns,rows);
m= (rows> columns)?rows:columns;
bd= (DOUBLE**)calloc(m,sizeof(DOUBLE*));
for(i= 0;i<rows;i++)bd[i]= (DOUBLE*)calloc(rows,sizeof(DOUBLE));
#endif 

us= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
cs= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
y= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
delta= (long*)calloc(columns+1,sizeof(long));
d= (long*)calloc(columns+1,sizeof(long));
first_nonzero= (int*)calloc(rows,sizeof(int));
first_nonzero_in_column= (int*)calloc(columns+rows+1,sizeof(int));
if(first_nonzero_in_column==NULL)return(0);
firstp= (int*)calloc(columns+1,sizeof(int));

eta= (long*)calloc(columns+1,sizeof(long));
v= (long*)calloc(columns+1,sizeof(long));
w= (DOUBLE**)calloc(columns+1,sizeof(DOUBLE*));
for(i= 0;i<=columns;i++)w[i]= (DOUBLE*)calloc(rows,sizeof(DOUBLE));

#if defined(FINCKEPOHST)
fipo= (DOUBLE*)calloc(columns+1,sizeof(DOUBLE));
muinv= (DOUBLE**)calloc(columns,sizeof(DOUBLE*));
for(i= 0;i<columns;++i)muinv[i]= (DOUBLE*)calloc(rows,sizeof(DOUBLE));
#endif
upb= (int*)calloc(columns+1,sizeof(int));
lowb= (int*)calloc(columns+1,sizeof(int));

/*:82*/
#line 1877 "diophant.w"
;
/*84:*/
#line 1969 "diophant.w"

for(i= 0;i<=columns;i++){
cs[i]= y[i]= us[i]= 0.0;
v[i]= delta[i]= 0;
eta[i]= d[i]= 1;
for(l= 0;l<rows;l++)w[i][l]= 0.0;
}

/*:84*/
#line 1878 "diophant.w"
;
/*87:*/
#line 2004 "diophant.w"

if(free_RHS){
i= 0;
#if defined(MPREC)
for(j= columns-1;j>=0;j--)if(!ziszero(get_entry(j,rows-2)))i++;
#else
for(j= columns-1;j>=0;j--)if(get_entry(j,rows-2)!=0)i++;
#endif
printf("Number of nonzero entries in the second last row: %d\n",i);
fflush(stdout);
}

i= 0;
#if defined(MPREC)
for(j= columns-1;j>=0;j--)if(!ziszero(get_entry(j,rows-1)))i++;
#else
for(j= columns-1;j>=0;j--)if(get_entry(j,rows-1)!=0)i++;
#endif
printf("Number of nonzero entries in the last row: %d\n",i);
fflush(stdout);

/*:87*/
#line 1879 "diophant.w"
;
/*86:*/
#line 1985 "diophant.w"

for(j= columns-1;j> 0;j--){
for(l= j-1;l>=0;l--){
#if defined(MPREC)
zcopy(get_entry(l,rows-1),&dummy1);
zcopy(get_entry(j,rows-1),&dummy2);
zabs(&dummy1);
zabs(&dummy2);
if(zcompare(dummy1,dummy2)==1){
#else
if(labs(get_entry(l,rows-1))> labs(get_entry(j,rows-1))){
#endif
swap_vec= lattice[l];
for(i= l+1;i<=j;i++)lattice[i-1]= lattice[i];
lattice[j]= swap_vec;
}
}
}

/*:86*/
#line 1880 "diophant.w"
;
/*88:*/
#line 2030 "diophant.w"

#if defined(MPREC)
Fq= zdoub(max_norm);
#else
Fq= (DOUBLE)max_norm;
#endif 
Fd= (rows*Fq*Fq)*(1.0+EPSILON);
#if VERBOSE >  0
printf("Fq: %f\n",(double)Fq);
printf("Fd: %f\n",(double)Fd);
#endif

/*:88*/
#line 1881 "diophant.w"
;

#if defined(GIVENS)
givens(lattice,columns,rows,mu,bd,c,N,Fq);
#else 
gramschmidt(lattice,columns,rows,mu,bd,c,N,Fq);
#endif 

#if defined(FINCKEPOHST)
/*90:*/
#line 2046 "diophant.w"

fiposort= 0;
fipo_success= 0;
inverse(mu,muinv,columns);

for(i= 0;i<columns;i++){
fipo[i]= 0.0;
dum1= 0.0;
for(j= 0;j<rows;j++){
dum= 0.0;
for(l= i;l<columns;l++)dum+= muinv[i][l]*bd[l][j]/c[l];
fipo[i]+= dum*dum;
dum1+= fabs(dum);
}
fipo[i]= SQRT(fipo[i]*Fd);
dum1= fabs(dum1*Fq);
if(dum1<fipo[i])fipo[i]= dum1;
fipo[i]*= (1.0+EPSILON);
}

/*:90*/
#line 1890 "diophant.w"
;
#endif

/*91:*/
#line 2070 "diophant.w"

for(l= 0;l<rows;l++){
#if defined(MPREC)
for(i= 0;i<columns;i++)if(!ziszero(get_entry(i,l))){
#else
for(i= 0;i<columns;i++)if(get_entry(i,l)!=0){
#endif 
first_nonzero[l]= i;
break;
}
}
j= 0;
for(l= 0;l<columns;l++){
firstp[l]= j;
first_nonzero_in_column[j]= 0;
j++;
for(i= 0;i<rows;i++){
if(first_nonzero[i]==l){
first_nonzero_in_column[j]= i;
first_nonzero_in_column[firstp[l]]++;
j++;
}
}
}
firstp[columns]= j;
first_nonzero_in_column[j]= 0;


/*:91*/
#line 1893 "diophant.w"
;

/*93:*/
#line 2103 "diophant.w"

level= first_nonzero[rows-1];
if(level<0)level= 0;
level_max= level;
us[level]= v[level]= 1;

only_zeros_no= only_zeros_success= 0;
hoelder_no= hoelder_success= 0;
cs_success= nosolutions= loops= 0;

/*:93*/
#line 1895 "diophant.w"
;

if(!SILENT)open_solution_file();

/*95:*/
#line 2131 "diophant.w"

do{
/*94:*/
#line 2115 "diophant.w"

loops++;
if((stop_after_loops> 0)&&(stop_after_loops<=loops))goto afterloop;
#if VERBOSE >  -1
if(loops%1000000==0){
#if defined(FINCKEPOHST)
printf("%ld loops, solutions: %ld, fipo: %ld\n",loops,nosolutions,fipo_success);
#else   
printf("%ld loops, solutions: %ld\n",loops,nosolutions);
#endif   
fflush(stdout);
}
#endif


/*:94*/
#line 2133 "diophant.w"
;
/*96:*/
#line 2185 "diophant.w"

dum= us[level]+y[level];
cs[level]= cs[level+1]+dum*dum*c[level];
/*:96*/
#line 2134 "diophant.w"
;

#if 0  
if(cs[level]>=Fd&&level<columns-1){
dums= dum-(us[level+1]+y[level+1])*mu[level+1][level];
if(dums*dums*cs[level]+cs[level+2]<Fd){
printf("step_back: %d\n",level);fflush(stdout);
level++;
goto step_back;
}
}
#endif

if((cs[level]<Fd)&&(!prune0(fabs(dum),N[level]))){
#if defined(FINCKEPOHST)
if(fabs(us[level])> fipo[level]){
fipo_success++;
goto side_step;
}
#endif
compute_w(w,bd,dum,level,rows);

if(level> 0){
/*97:*/
#line 2191 "diophant.w"

if(prune_only_zeros(w[level],level,rows,Fq,first_nonzero_in_column,firstp))
goto side_step;

if(prune(w[level],cs[level],rows,Fq)){
if(eta[level]==1)goto step_back;
eta[level]= 1;
delta[level]*= -1;
if(delta[level]*d[level]>=0)delta[level]+= d[level];
us[level]= v[level]+delta[level];
}else{
level--;
eta[level]= 0;
delta[level]= 0;
dum= compute_y(mu,us,level,level_max);
y[level]= dum;
v[level]= (long)round(-dum);
us[level]= v[level];
d[level]= (v[level]> -y[level])?-1:1;
}
/*:97*/
#line 2157 "diophant.w"
;
}else{
/*98:*/
#line 2212 "diophant.w"

if(exacttest(w[0],rows,Fq)==1){
print_solution(w[level],rows,Fq);
if((stop_after_solutions> 0)&&(stop_after_solutions<=nosolutions))
goto afterloop;
}
goto side_step;
/*:98*/
#line 2159 "diophant.w"
;
}
}else{
cs_success++;
step_back:

level++;
if(level_max<level)level_max= level;
side_step:




if(eta[level]==0){
delta[level]*= -1;
if(delta[level]*d[level]>=0)delta[level]+= d[level];
}else{
delta[level]+= (delta[level]*d[level]>=0)?d[level]:-d[level];
}
us[level]= v[level]+delta[level];
}
}while(level<columns);
afterloop:

/*:95*/
#line 1899 "diophant.w"
;

if(!SILENT)close_solution_file();

/*99:*/
#line 2220 "diophant.w"

printf("Prune_cs: %ld\n",cs_success);
printf("Prune_only_zeros: %ld of %ld\n",only_zeros_success,only_zeros_no);
printf("Prune_hoelder: %ld of %ld\n",hoelder_success,hoelder_no);
#if defined(FINCKEPOHST)
printf("Fincke-Pohst: %ld\n",fipo_success);
#endif 
printf("Loops: %ld\n",loops);
if((stop_after_solutions<=nosolutions&&stop_after_solutions> 0)||
(stop_after_loops<=loops&&stop_after_loops> 0)){
printf("Stopped after number of solutions: %ld\n",nosolutions);
}else{
printf("Total number of solutions: %ld\n",nosolutions);
}
printf("\n");fflush(stdout);

/*:99*/
#line 1903 "diophant.w"
;
/*100:*/
#line 2236 "diophant.w"

free(us);
free(cs);
free(y);
free(delta);
free(d);
free(first_nonzero);
free(first_nonzero_in_column);
free(firstp);
free(eta);
free(v);
for(l= 0;l<=columns;l++)free(w[l]);
free(w);
#if defined(FINCKEPOHST)
free(fipo);
#endif 
#if defined(MPREC)
lllfree(mu,c,N,bd,columns);
#else
lllfree(mu,c,N,columns);
#endif 

/*:100*/
#line 1904 "diophant.w"
;

return 1;
}

/*:79*/
#line 799 "diophant.w"
;

/*:40*/
#line 50 "diophant.w"
;

/*:2*//*3:*/
#line 53 "diophant.w"

#if defined(MPREC)
long diophant(verylong**a_input,verylong*b_input,verylong*upperbounds_input,
int no_columns,int no_rows,
verylong factor_input,verylong norm_input,
int silent,int iterate,int iterate_no,
int bkz_beta_input,int bkz_p_input,
int stop_after_sol_input,int stop_after_loops_input,
int free_RHS_input,int*org_col_input,int no_org_col_input)
#else
long diophant(long**a_input,long*b_input,long*upperbounds_input,
int no_columns,int no_rows,
long factor_input,long norm_input,
int silent,int iterate,int iterate_no,
int bkz_beta_input,int bkz_p_input,
int stop_after_sol_input,int stop_after_loops_input,
int free_RHS_input,int*org_col_input,int no_org_col_input)
#endif
{
int i,j;
COEFF*swap_vec;

/*11:*/
#line 224 "diophant.w"

#if defined(MPREC)
zcopy(factor_input,&matrix_factor);
zcopy(norm_input,&max_norm);
#else
matrix_factor= factor_input;
max_norm= norm_input;
#endif
if(iterate){
no_iterates= iterate_no;
}else{
bkz_beta= bkz_beta_input;
bkz_p= bkz_p_input;
}
SILENT= silent;
stop_after_solutions= stop_after_sol_input;
stop_after_loops= stop_after_loops_input;
free_RHS= free_RHS_input;
nom= 1;
denom= 2;

system_rows= no_rows;
system_columns= no_columns;
/*:11*/
#line 75 "diophant.w"
;
/*12:*/
#line 252 "diophant.w"

lattice_rows= system_rows+system_columns+1;
lattice_columns= system_columns+2;

if(free_RHS){
lattice_rows++;
lattice_columns++;
}else{
printf("The RHS is fixed !\n");fflush(stdout);
}
/*:12*/
#line 76 "diophant.w"
;
/*13:*/
#line 263 "diophant.w"

lattice= (COEFF**)calloc(lattice_columns,sizeof(COEFF*));
for(j= 0;j<lattice_columns;j++){
lattice[j]= (COEFF*)calloc(lattice_rows+1,sizeof(COEFF));
#if defined(MPREC)
for(i= 0;i<=lattice_rows;i++)zzero(&(lattice[j][i].c));
#else
for(i= 0;i<=lattice_rows;i++)lattice[j][i].c= 0;
#endif  
}
/*:13*/
#line 77 "diophant.w"
;
/*14:*/
#line 275 "diophant.w"

#if defined(MPREC)
for(j= 0;j<system_rows;j++){
for(i= 0;i<system_columns;i++){
zcopy(a_input[j][i],&(lattice[i][j+1].c));
smult_lattice(i,j,matrix_factor);
}
zcopy(b_input[j],&(lattice[system_columns][j+1].c));
smult_lattice(system_columns,j,matrix_factor);
}
#else
for(j= 0;j<system_rows;j++){
for(i= 0;i<system_columns;i++){
lattice[i][j+1].c= a_input[j][i];
smult_lattice(i,j,matrix_factor);
}
lattice[system_columns][j+1].c= b_input[j];
smult_lattice(system_columns,j,matrix_factor);
}
#endif 
/*:14*/
#line 78 "diophant.w"
;
/*15:*/
#line 308 "diophant.w"

#if defined(MPREC)
zone(&upperbounds_max);
iszeroone= 1;
if(upperbounds_input==NULL){
printf("No upper bounds: 0/1 variables are assumed \n");fflush(stdout);
}else{
upperbounds= (verylong*)calloc(system_columns,sizeof(verylong));
for(i= 0;i<system_columns;i++){
upperbounds[i]= 0;
zcopy(upperbounds_input[i],&(upperbounds[i]));
zgcd(upperbounds_max,upperbounds[i],&dummy);
zmulin(upperbounds[i],&upperbounds_max);
zdiv(upperbounds_max,dummy,&upperbounds_max,&dummy);
}
if(zscompare(upperbounds_max,1)==1)iszeroone= 0;
printf("upper bounds found. Max=");
zwriteln(upperbounds_max);fflush(stdout);
}
#else
upperbounds_max= 1;
iszeroone= 1;
if(upperbounds_input==NULL){
printf("No upper bounds: 0/1 variables are assumed \n");fflush(stdout);
}else{
upperbounds= (long*)calloc(system_columns,sizeof(long));
for(i= 0;i<system_columns;i++){
upperbounds[i]= upperbounds_input[i];
upperbounds_max*= upperbounds[i]/gcd(upperbounds[i],upperbounds_max);
}
if(upperbounds_max> 1)iszeroone= 0;
printf("upper bounds found. Max=%ld\n",upperbounds_max);fflush(stdout);
}
#endif
/*:15*/
#line 79 "diophant.w"
;
/*16:*/
#line 346 "diophant.w"

if(org_col_input!=NULL)no_original_columns= no_org_col_input;
else no_original_columns= system_columns;


original_columns= (int*)calloc(no_original_columns,sizeof(int));

if(org_col_input!=NULL)
for(i= 0;i<no_original_columns;i++)original_columns[i]= org_col_input[i];
else{
for(i= 0;i<no_original_columns;i++)original_columns[i]= 1;
printf("No preselected columns \n");fflush(stdout);
}

/*:16*/
#line 80 "diophant.w"
;
/*17:*/
#line 368 "diophant.w"

#if defined(MPREC)
for(j= system_rows;j<lattice_rows;j++){
zsmul(max_norm,denom,&(lattice[j-system_rows][j+1].c));
zsmul(max_norm,nom,&(lattice[lattice_columns-2][j+1].c));
}

zcopy(max_norm,&(lattice[system_columns+free_RHS][lattice_rows].c));

for(j= system_rows;j<lattice_rows;j++){
zsmul(max_norm,denom,&(lattice[j-system_rows][j+1].c));
zsmul(max_norm,nom,&(lattice[lattice_columns-2][j+1].c));
}

if(free_RHS){
zone(&(lattice[system_columns][lattice_rows-1].c));
zzero(&(lattice[system_columns+1][lattice_rows-1].c));
}
zcopy(max_norm,&(lattice[system_columns+free_RHS][lattice_rows].c));
#else
for(j= system_rows;j<lattice_rows;j++){
put_to(j-system_rows,j,denom*max_norm);
put_to(lattice_columns-2,j,nom*max_norm);
}

put_to(system_columns+free_RHS,lattice_rows-1,max_norm);

for(j= system_rows;j<lattice_rows;j++){
put_to(j-system_rows,j,denom*max_norm);
put_to(lattice_columns-2,j,nom*max_norm);
}

if(free_RHS){
put_to(system_columns,lattice_rows-2,1);
put_to(system_columns+1,lattice_rows-2,0);
}
put_to(system_columns+free_RHS,lattice_rows-1,max_norm);
#endif
for(i= 0;i<lattice_columns-1;i++)coeffinit(lattice[i],lattice_rows);
/*:17*/
#line 81 "diophant.w"
;
/*18:*/
#line 415 "diophant.w"

#if defined(MPREC)
zcopy(max_norm,&max_norm_initial);
zone(&max_up);
if(!iszeroone){
for(j= 0;j<system_columns;j++){
zdiv(upperbounds_max,upperbounds[j],&upfac,&upfac_h);
smult_lattice(j,j+system_rows,upfac);
smult_lattice(system_columns+free_RHS,j+system_rows,upperbounds_max);
}
zcopy(upperbounds_max,&max_up);
zmulin(max_up,&max_norm);
#else
max_norm_initial= max_norm;
max_up= 1;
if(!iszeroone){
for(j= 0;j<system_columns;j++){
upfac= upperbounds_max/upperbounds[j];
smult_lattice(j,j+system_rows,(long)upfac);
smult_lattice(system_columns+free_RHS,j+system_rows,upperbounds_max);
}
max_up= upperbounds_max;
max_norm*= max_up;
#endif
if(free_RHS)
smult_lattice(system_columns,lattice_rows-2,max_up);

smult_lattice(system_columns+free_RHS,lattice_rows-1,max_up);
}
/*:18*/
#line 82 "diophant.w"
;
/*19:*/
#line 447 "diophant.w"

swap_vec= lattice[lattice_columns-2];
for(i= lattice_columns-2;i> 0;i--)lattice[i]= lattice[i-1];
lattice[0]= swap_vec;
/*:19*/
#line 83 "diophant.w"
;
/*20:*/
#line 453 "diophant.w"

lastlines_factor= 1;
printf("\n");fflush(stdout);
lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_LOW);

/*:20*/
#line 84 "diophant.w"
;
#if 0
print_lattice();
exit(0);
#endif
/*22:*/
#line 473 "diophant.w"

if(cutlattice()){
printf("First reduction successful\n");fflush(stdout);
}else{
printf("First reduction not successful\n");fflush(stdout);
return 0;
}

for(j= 0;j<lattice_columns-1;j++)solutiontest(j);
/*:22*/
#line 89 "diophant.w"
;
/*21:*/
#line 460 "diophant.w"

lastlines_factor= 1;
lll(lattice,lattice_columns-1,lattice_rows,LLLCONST_HIGH);
printf("Second reduction successful\n");fflush(stdout);

/*:21*/
#line 90 "diophant.w"
;
/*23:*/
#line 488 "diophant.w"

lastlines_factor= LASTLINESFACTOR;
#if defined(MPREC)
for(i= 0;i<lattice_columns;i++)
zsmul(lattice[i][lattice_rows].c,lastlines_factor,&(lattice[i][lattice_rows].c));
if(free_RHS)
for(i= 0;i<lattice_columns;i++)
zsmul(lattice[i][lattice_rows-1].c,lastlines_factor,&(lattice[i][lattice_rows-1].c));
#else
for(i= 0;i<lattice_columns;i++)smult_lattice(i,lattice_rows-1,lastlines_factor);
if(free_RHS)
for(i= 0;i<lattice_columns;i++)smult_lattice(i,lattice_rows-2,lastlines_factor);
#endif
/*:23*/
#line 91 "diophant.w"
;
/*24:*/
#line 503 "diophant.w"

printf("\n");fflush(stdout);
if(iterate)
iteratedlll(lattice,lattice_columns-1,lattice_rows,no_iterates,LLLCONST_HIGH);
else
bkz(lattice,lattice_columns,lattice_rows,LLLCONST_HIGH,bkz_beta,bkz_p);
printf("Third reduction successful\n");fflush(stdout);

/*:24*/
#line 92 "diophant.w"
;
/*25:*/
#line 513 "diophant.w"

#if defined(MPREC)
for(i= 0;i<lattice_columns;i++)
zsdiv(lattice[i][lattice_rows].c,lastlines_factor,&(lattice[i][lattice_rows].c));
if(free_RHS)
for(i= 0;i<lattice_columns;i++)
zsdiv(lattice[i][lattice_rows-1].c,lastlines_factor,&(lattice[i][lattice_rows-1].c));
#else
for(i= 0;i<lattice_columns;i++)sdiv_lattice(i,lattice_rows-1,lastlines_factor);
if(free_RHS)
for(i= 0;i<lattice_columns;i++)sdiv_lattice(i,lattice_rows-2,lastlines_factor);
#endif
/*:25*/
#line 93 "diophant.w"
;
/*26:*/
#line 526 "diophant.w"

printf("\n");fflush(stdout);
nosolutions= explicit_enumeration(lattice,lattice_columns-1,lattice_rows);

/*:26*/
#line 94 "diophant.w"
;

return nosolutions;
}
/*:3*/
