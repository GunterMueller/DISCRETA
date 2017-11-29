/* file:lo.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


#include <DISCRETA/discreta.h>

#ifdef LONGINTTRUE

#include <DISCRETA/lo.h>
#include <DISCRETA/ma.h>
#include <DISCRETA/bruch.h>
#include <DISCRETA/poly.h>
#include <DISCRETA/divs.h>


#define EXP 15
#define B 32768           /*              1000000000000000*/
#define BMINUSEINS 32767  /*               111111111111111*/
#define B1 2147450880L    /* 111111111111111000000000000000*/
#define B2MINUSEINS 2147483647L    /*1000000000000000000000000000000 - 1*/
#define Basis 45
#define MSB  16384 
#define FREE_LOC(a) my_free(a)

struct ganzdaten gd;

/* alles.c */

static INT locadd(struct loc *lx, struct loc *ly, INT cy);
static INT locdiv(struct loc *qu, struct loc *rest,
	struct loc *dd, struct loc *dv);
static INT locint(struct loc *lx, INT i);
static INT locms1(struct loc *lx);
static INT locmul(struct loc *ly, struct loc *lx,
	struct loc *la, struct loc *lb);
static INT locneg(struct loc *lx, INT cy);
static INT locnull(struct loc *lx);
static INT locodd(struct loc *lx);
static INT locpsl(struct loc *lx, struct loc *ly, INT a);
static INT locpsr(struct loc *lx, struct loc *ly, INT a);
static INT locsadd(struct loc *lx, INT i);
static INT locsdiv(struct loc *qu, INT di, struct loc *dd, INT dv);
static INT locsgn(struct loc *lx);
static INT locsmul(struct loc *lx, INT i, INT ue);
static INT locssub(struct loc *lx, INT i);
static INT locsub(struct loc *lx, struct loc *ly, INT cy);
static INT locvgl(struct loc *lx, struct loc *ly);
static INT ganzadd(struct longint *x, struct longint *y);
static INT ganzsquores(struct longint *x, INT *rest, INT y);
static INT ganzquores(struct longint *x, struct longint *rest, struct longint *y);
static INT ganzganzdiv(struct longint *x, struct longint *y);
static INT ganzmod(struct longint *x, struct longint *rest, struct longint *y);
static INT ganzein(FILE *fp, struct longint *x);
static INT holeziffer(struct zahldaten *zd);
static INT ganzfziffer(struct zahldaten *zd, struct ganzdaten *d);
static INT retteziffer(INT z, struct zahldaten *zd);
static INT ganz1ziffer(struct zahldaten *zd, struct longint	*x);
static INT ganzaus(FILE *fp, struct longint *x);
static INT ganzmul(struct longint *x, struct longint *y);
static INT ganzsmul(struct longint *x, INT a);
static INT ganzsadd(struct longint *x, INT y);
static INT ganzvergleich(struct longint *x, struct longint *y);
static INT intganz(struct longint *x);
static INT ganzint(struct longint *x, INT i);
static INT ganzsignum(struct longint *x);
static INT ganzeven(struct longint *x);
static INT ganzodd(struct longint *x);
static INT ganzkopiere(struct longint *x, struct longint *a);
static INT ganzneg(struct longint *x);
static INT lochole(struct loc **aloc);
static INT loclisterette(struct loc **aloc);
static INT locrette(struct loc **aloc);
static INT locrezliste(struct loc **aloc);
static struct longint * calloclongint(void);

static INT ganzparam(INT basis, INT auspos, INT auslaenge, char folgezeichen,
	struct ganzdaten *d);
static INT ganzanfang(struct ganzdaten *d);
static INT ganzdefiniere(struct longint *x);
static INT ganzloesche(struct longint *x);
static INT ganzein_str(BYTE *str, struct longint *x);
static INT lo_s_scan_int(BYTE **s, INT *i);
static INT lo_s_scan_token(BYTE **s, BYTE *str);
static INT ganzaus_str(struct longint *x, BYTE *str);

/* dieser Teil wurde von Peter Hain in Karlsruhe
entworfen. Er schrieb diese Langzahl arithmetik in Pascal und
Assembler. In Bayreuth wurden in Form eines Seminars die
Assemblerteile in C geschrieben und spaeter wurde von
Axel Kohnert die restlichen Pascalteile in C uebersetzt.
Die Ein und Ausgabe routinen wurden vollstaendig in Bayreuth
entworfen */

static INT locadd(struct loc *lx, struct loc *ly, INT cy)
	{
	static INT hh;
	hh=ly->w0+cy+lx->w0;
	lx->w0=(hh&BMINUSEINS);
	cy = hh >>EXP;
	hh=ly->w1+cy+lx->w1;
	lx->w1=(hh&BMINUSEINS);
	cy = hh >>EXP;
	hh=ly->w2+cy+lx->w2;
	lx->w2=(hh&BMINUSEINS);
	cy = hh >>EXP;
	return((INT)cy);
	}

#define LOCADD(lx,ly,cy)\
	hh=(ly)->w0+cy+(lx)->w0, (lx)->w0=(hh&BMINUSEINS), cy = hh >>EXP,\
	hh=(ly)->w1+cy+(lx)->w1, (lx)->w1=(hh&BMINUSEINS), cy = hh >>EXP,\
	hh=(ly)->w2+cy+(lx)->w2, (lx)->w2=(hh&BMINUSEINS), cy = hh >>EXP,cy


#define LOCBAS2() Basis
    

#define LOCASS(lx,ly) ((lx)->w2=(ly)->w2,(lx)->w1=(ly)->w1,(lx)->w0=(ly)->w0)

static INT locdiv(struct loc *qu, struct loc *rest,
	struct loc *dd, struct loc *dv)
/* Division. Bei Eingabe muss gelten: rest<dv.            */
/* (qu,rest) := ((rest*B+dd) DIV dv, (rest*B+dd) MOD dv). */
	{
	INT d6,d5,d4,d3,d2,d1,
	     h6,h5,h4,h3,h2,h1, 
	     m2,m1,m0;


 /* d=rest*B+dd */
 d6=rest->w2; d5=rest->w1; d4=rest->w0; d3=dd->w2; d2=dd->w1; d1=dd->w0;

 /* h=dv */
 h6=0L; h5=0L; h4=0L; h3=dv->w2; h2=dv->w1; h1=dv->w0; 

 /* qu=0 */
 qu->w2=0L; qu->w1=0L; qu->w0=0L;

 /* m=1 */
 m2=0L; m1=0L; m0=1L;

 while /* h<=d */
 (h6 <d6 ||
  (h6==d6 && h5<d5 ) ||
  (h6==d6 && h5==d5 && h4<d4 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3<d3 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2<d2 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2==d2 && h1<d1 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2==d2 && h1==d1 )
 )
 {
  /* h=h*2 */
  h6=h6<<1; h5=h5<<1; h4=h4<<1; h3=h3<<1; h2=h2<<1; h1=h1<<1; 
  if (h1&B1) {h1=h1&BMINUSEINS; h2++; };
  if (h2&B1) {h2=h2&BMINUSEINS; h3++; };
  if (h3&B1) {h3=h3&BMINUSEINS; h4++; };
  if (h4&B1) {h4=h4&BMINUSEINS; h5++; };
  if (h5&B1) {h5=h5&BMINUSEINS; h6++; };

  /* m=m*2 */
  m2=m2<<1; m1=m1<<1; m0=m0<<1;
  if (m0&B1) {m0=m0&BMINUSEINS; m1++; };
  if (m1&B1) {m1=m1&BMINUSEINS; m2++; };
 }
  while /* d>=dv */
  (d6 >0 || d5>0  || d4>0  ||
   (d6==0 && d5==0 && d4==0 && d3>dv->w2 ) || 
   (d6==0 && d5==0 && d4==0 && d3==dv->w2 && d2>dv->w1 ) || 
   (d6==0 && d5==0 && d4==0 && d3==dv->w2 && d2==dv->w1 && d1>dv->w0 ) ||
   (d6==0 && d5==0 && d4==0 && d3==dv->w2 && d2==dv->w1 && d1==dv->w0 )
  )
  {
   while /*h>d */
   (h6 >d6 ||
    (h6==d6 && h5>d5) ||
    (h6==d6 && h5==d5 && h4>d4) ||
    (h6==d6 && h5==d5 && h4==d4 && h3>d3) ||
    (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2>d2) ||
    (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2==d2 && h1>d1 )
   )
   {
    /* h=h/2 */
    if (h6&1) { h6--; h5=h5|B; };  h6=h6>>1;
    if (h5&1) { h5--; h4=h4|B; };  h5=h5>>1;
    if (h4&1) { h4--; h3=h3|B; };  h4=h4>>1;
    if (h3&1) { h3--; h2=h2|B; };  h3=h3>>1;
    if (h2&1) { h2--; h1=h1|B; };  h2=h2>>1;
                                   h1=h1>>1;

    /* m=m/2 */
    if (m2&1) { m2--; m1=m1|B; };  m2=m2>>1;
    if (m1&1) { m1--; m0=m0|B; };  m1=m1>>1;
                                   m0=m0>>1;
   }
   
   /* d=d-h */
   if (h1>d1) { d1=d1+B; d2--; };  d1=d1-h1;
   if (h2>d2) { d2=d2+B; d3--; };  d2=d2-h2;
   if (h3>d3) { d3=d3+B; d4--; };  d3=d3-h3;
   if (h4>d4) { d4=d4+B; d5--; };  d4=d4-h4;
   if (h5>d5) { d5=d5+B; d6--; };  d5=d5-h5;
                                   d6=d6-h6;

   /* qu=qu+m */
   qu->w0=qu->w0|m0; qu->w1=qu->w1|m1; qu->w2=qu->w2|m2;
  }
	rest->w2=d3; rest->w1=d2; rest->w0=d1;
	return(OK);
}

static INT locint(struct loc *lx, INT i)
/* Umwandlung Integer in loc: lx:=abs(i); locint:=sgn(i)  */
	{
	INT s;

	if (i<0L) 		{s=0-1L; i=0-i;}
	else if (i>0L) 		{s=1L;}
	else 			{s=0L;}

	lx->w0=i & BMINUSEINS;
	lx->w1=(i & B1)>>EXP; 
	lx->w2=(lx->w1 & B1)>>EXP;
	lx->w1=lx->w1 & BMINUSEINS;
	return(s);
	}


#define LOCMAX(lx) ((lx)->w0=BMINUSEINS,\
	(lx)->w1=BMINUSEINS,(lx)->w2=BMINUSEINS)

static INT locms1(struct loc *lx)
	{
	INT j,c,cc;

	c=Basis; cc=1L;
    for (j=14L; (j >=0L) && cc ; j--)
    {
      if ( lx->w2 & ( 1L << j )) 
         cc=0L;
      c--;
    } 
  if (cc) 
    for (j=14L; (j >=0L) && cc ; j--)
    {
      if ( lx->w1 & ( 1L << j )) 
         cc=0L;
      c--;
    }
  if (cc) 
    for (j=14L; (j >=0L) && cc ; j--)
    {
      if ( lx->w0 & ( 1L << j )) 
         cc=0L;
      c--;
    }
  if (cc) error("locms1:");
  return(c);
 }

 
static INT locmul(struct loc *ly, struct loc *lx,
	struct loc *la, struct loc *lb)
{
static INT hh;
#define teile(z) (hh=(z)>>EXP, (z) &= BMINUSEINS, hh)
/* AK 260390 */

 lx->w0 =                          la->w0 * lb->w0;
 lx->w1 =          teile(lx->w0) + la->w1 * lb->w0;
 lx->w2 =          teile(lx->w1) + la->w2 * lb->w0;
 ly->w0 =          teile(lx->w2)                  ;
 lx->w1 +=                  la->w0 * lb->w1;
 lx->w2 +=  teile(lx->w1) + la->w1 * lb->w1;
 ly->w0 +=  teile(lx->w2) + la->w2 * lb->w1;
 ly->w1 =          teile(ly->w0)                  ;
 lx->w2 +=                  la->w0 * lb->w2;
 ly->w0 +=  teile(lx->w2) + la->w1 * lb->w2;
 ly->w1 +=  teile(ly->w0) + la->w2 * lb->w2;
 ly->w2 =          teile(ly->w1)                  ;
 return OK;
}

#define LOCMUL(ly,lx,la,lb) /* hh ist noetig */\
 (lx)->w0 =                          (la)->w0 * (lb)->w0,\
 (lx)->w1 =          teile((lx)->w0) + (la)->w1 * (lb)->w0,\
 (lx)->w2 =          teile((lx)->w1) + (la)->w2 * (lb)->w0,\
 (ly)->w0 =          teile((lx)->w2)                  ,\
 (lx)->w1 +=                  (la)->w0 * (lb)->w1,\
 (lx)->w2 +=  teile((lx)->w1) + (la)->w1 * (lb)->w1,\
 (ly)->w0 +=  teile((lx)->w2) + (la)->w2 * (lb)->w1,\
 (ly)->w1 =          teile((ly)->w0)                  ,\
 (lx)->w2 +=                  (la)->w0 * (lb)->w2,\
 (ly)->w0 +=  teile((lx)->w2) + (la)->w1 * (lb)->w2,\
 (ly)->w1 +=  teile((ly)->w0) + (la)->w2 * (lb)->w2,\
 (ly)->w2 =          teile((ly)->w1)                  

static INT locneg(struct loc *lx, INT cy)
/* AK 130789 V1.0 */ /* AK 050790 V1.1 */ /* AK 210891 V1.3 */
	{
if ((cy==0L)&&(lx->w0==0L)&&(lx->w1==0L)&&(lx->w2==0L)) {
		return(0L); }
else

	{
	lx->w0 ^= BMINUSEINS;
	lx->w1 ^= BMINUSEINS;
	lx->w2 ^= BMINUSEINS;
	if (cy == 0L )
	  { ++lx->w0;
	  if (lx->w0 & B)
	    {
	    ++lx->w1;
	    lx->w0 &= BMINUSEINS;
	    if (lx->w1 & B)
	      { ++lx->w2; lx->w1 &= BMINUSEINS; }
            }
          }
       return(1L);
       }
  } /* locneg */

static INT locnull(struct loc *lx)
	{
	lx->w2 = 0L; lx->w1 = 0L; lx->w0 = 0L;
	return OK;
	}  /* Ende von locnull */

#define LOCNULL(lx) ((lx)->w0=0L,(lx)->w1=0L,(lx)->w2=0L)

static INT locodd(struct loc *lx)
/*locodd:=lx ist ungerade */
	{
	return (INT) (lx->w0 & 1);
	}


static INT locpsl(struct loc *lx, struct loc *ly, INT a)
	{
   INT s1,s2,s3,s4,s5,i;
   struct loc lyy;


   lyy.w0 = ly->w0; lyy.w1 = ly->w1; lyy.w2 = ly->w2;

   if ( a >= Basis) error("locpsl:");
   for (i=1L; i <= a;i++)
   {
     s1= (lyy.w0 & MSB) >> 14;
     s2= (lyy.w1 & MSB) >> 14;
     s3= (lyy.w2 & MSB) >> 14;
     s4= (lx->w0 & MSB) >> 14;
     s5= (lx->w1 & MSB) >> 14;
     lyy.w0 = lyy.w0 << 1;
     lyy.w1 = (lyy.w1 << 1) | s1;
     lyy.w2 = (lyy.w2 << 1) | s2;
     lx->w0 = (lx->w0 << 1) | s3;
     lx->w1 = (lx->w1 << 1) | s4;
     lx->w2 = (lx->w2 << 1) | s5;
   }
   lx->w0 = lx->w0 & BMINUSEINS;
   lx->w1 = lx->w1 & BMINUSEINS;
   lx->w2 = lx->w2 & BMINUSEINS;
   return OK;
 } /* Ende von locpsl */
     
static INT locpsr(struct loc *lx, struct loc *ly, INT a)
 {
   INT s1,s2,s3,s4,s5,i;
   struct loc lxx;


   lxx.w0 = lx->w0; lxx.w1 = lx->w1; lxx.w2 = lx->w2;


   if ( a >= Basis) error("locpsr:");
   for (i=1L; i <= a;i++)
   {
     s1= (ly->w1 & 1) << 14;
     s2= (ly->w2 & 1) << 14;
     s3= (lxx.w0 & 1) << 14;
     s4= (lxx.w1 & 1) << 14;
     s5= (lxx.w2 & 1) << 14;
     ly->w0 = (ly->w0 >> 1) | s1;
     ly->w1 = (ly->w1 >> 1) | s2;
     ly->w2 = (ly->w2 >> 1) | s3;
     lxx.w0 = (lxx.w0 >> 1) | s4;
     lxx.w1 = (lxx.w1 >> 1) | s5;
     lxx.w2 = (lxx.w2 >> 1);
   }
   return OK;
 } /* Ende von locpsr */

static INT locsadd(struct loc *lx, INT i)
{
	INT cy,hh;
	if (i<0L) i=(-i);
	hh=lx->w0+(i%B);
	lx->w0=(hh & BMINUSEINS);
	cy = hh >>EXP;
	hh=lx->w1+(i/B)+cy;
	lx->w1=(hh & BMINUSEINS);
	cy = hh >>EXP;
	hh=lx->w2+cy;
	lx->w2=(hh & BMINUSEINS);
	cy = hh >>EXP;
	return(cy);
	}


static INT locsdiv(struct loc *qu, INT di, struct loc *dd, INT dv)
/* Division. Bei Eingabe muss gelten: di<dv.             */
/* (locsdiv,qu) := ((di*B+dd) MOD dv, (di*B+dd) DIV dv). */
	{
 INT d6,d5,d4,d3,d2,d1,
     h6,h5,h4,h3,h2,h1, 
     m2,m1,m0,
     dv2,dv1,dv0;

 /* di umwandeln */
 if (di<0L) {return error("locsdiv:di<0");};
 d4=di & BMINUSEINS;
 d5=(di & B1)>>EXP; 
 d6=(d5 & B1)>>EXP;
 d5=d5 & BMINUSEINS;

 /* dv umwandeln */
 if (dv<0L) {return error("locsdiv:dv<0");};
 dv0=dv & BMINUSEINS;
 dv1=(dv & B1)>>EXP; 
 dv2=(dv1 & B1)>>EXP;
 dv1=dv1 & BMINUSEINS;

 /* d=di*B+dd */
 d3=dd->w2; d2=dd->w1; d1=dd->w0;

 /* h=dv */
 h6=0L; h5=0L; h4=0L; h3=dv2; h2=dv1; h1=dv0; 

 /* qu=0 */
 qu->w2=0L; qu->w1=0L; qu->w0=0L;

 /* m=1 */
 m2=0L; m1=0L; m0=1L;

 while /* h<=d */
 (h6 <d6 ||
  (h6==d6 && h5<d5 ) ||
  (h6==d6 && h5==d5 && h4<d4 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3<d3 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2<d2 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2==d2 && h1<d1 ) ||
  (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2==d2 && h1==d1) 
 )
 {
	  /* h=h*2 */
	  h6=h6<<1; h5=h5<<1; h4=h4<<1; h3=h3<<1; h2=h2<<1; h1=h1<<1; 
	  if (h1&B1) {h1&=BMINUSEINS; h2++; };
	  if (h2&B1) {h2&=BMINUSEINS; h3++; };
	  if (h3&B1) {h3&=BMINUSEINS; h4++; };
	  if (h4&B1) {h4&=BMINUSEINS; h5++; };
	  if (h5&B1) {h5&=BMINUSEINS; h6++; };

	  /* m=m*2 */
	  m2=m2<<1; m1=m1<<1; m0=m0<<1;
	  if (m0&B1) {m0&=BMINUSEINS; m1++; };
	  if (m1&B1) {m1&=BMINUSEINS; m2++; };
 }
  while /* d>=dv */
  (d6 >0L || d5>0L  || d4>0L  ||
   (d6==0L && d5==0L && d4==0L && d3>dv2 ) || 
   (d6==0L && d5==0L && d4==0L && d3==dv2 && d2>dv1 ) || 
   (d6==0L && d5==0L && d4==0L && d3==dv2 && d2==dv1 && d1>dv0 ) ||
   (d6==0L && d5==0L && d4==0L && d3==dv2 && d2==dv1 && d1==dv0 )
  )
  {
		   while /*h>d */
		   (h6 >d6 ||
		    (h6==d6 && h5>d5) ||
		    (h6==d6 && h5==d5 && h4>d4) ||
		    (h6==d6 && h5==d5 && h4==d4 && h3>d3) ||
		    (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2>d2) ||
		    (h6==d6 && h5==d5 && h4==d4 && h3==d3 && h2==d2 && h1>d1 )
		   )
		   {
			    /* h=h/2 */
			    if (h6&1) { h6--; h5|=B; };  h6=h6>>1;
			    if (h5&1) { h5--; h4|=B; };  h5=h5>>1;
			    if (h4&1) { h4--; h3|=B; };  h4=h4>>1;
			    if (h3&1) { h3--; h2|=B; };  h3=h3>>1;
			    if (h2&1) { h2--; h1|=B; };  h2=h2>>1;
							   h1=h1>>1;

			    /* m=m/2 */
			    if (m2&1) { m2--; m1=m1|B; };  m2=m2>>1;
			    if (m1&1) { m1--; m0=m0|B; };  m1=m1>>1;
							   m0=m0>>1;
		   }
		   
		   /* d=d-h */
		   if (h1>d1) { d1+=B; d2--; };  d1-=h1;
		   if (h2>d2) { d2+=B; d3--; };  d2-=h2;
		   if (h3>d3) { d3+=B; d4--; };  d3-=h3;
		   if (h4>d4) { d4+=B; d5--; };  d4-=h4;
		   if (h5>d5) { d5+=B; d6--; };  d5-=h5;
						 d6-=h6;

		   /* qu=qu+m */
		   qu->w0|=m0; qu->w1|=m1; qu->w2|=m2;
  }
  d3=d3<<EXP|d2;
  d3=d3<<EXP|d1;
  return(d3);
	}

static INT locsgn(struct loc *lx)
	{
	  if (lx->w2 || lx->w1 || lx->w0 ) 
	     return(1L);
	  else return(0L);
	} /* Ende locsgn */


static INT locsmul(struct loc *lx, INT i, INT ue)
	{
	INT cy,h0,h1,h2,i0,i1,i2,u0,u1,u2;
	if (i<0)  {i=~i;++i;}
	if (ue<0) {ue=~ue;++ue;} 
	i0=i&BMINUSEINS; i1=((i&B1)>>15); i2=((i&(B2MINUSEINS+1))>>30);
	u0=ue&BMINUSEINS; u1=((ue&B1)>>15); u2=((ue&(B2MINUSEINS+1))>>30);
	h0=(lx->w0)*i0+u0;
		cy = (h0 >> 15);
		h0 &= BMINUSEINS;
	h1=(lx->w0)*i1+(lx->w1)*i0+cy+u1;
		cy = (h1 >> 15);
		h1 &= BMINUSEINS;
	h2=(lx->w0)*i2+(lx->w1)*i1+(lx->w2)*i0+cy+u2;
	cy = (h2 >> 15);
	h2 &= BMINUSEINS;
	cy += (lx->w1)*i2+(lx->w2)*i1+(((lx->w2)*i2)<<15);
	lx->w0 = h0; lx->w1 = h1; lx->w2 = h2;
	return(cy);
	}

static INT locssub(struct loc *lx, INT i)
	{
	INT cy;
	if (i<0L) i=(-i);
	lx->w0=lx->w0-(i%B);
	if (lx->w0 < 0) { lx->w0 += B;
			   cy = 1L; }
	else cy = 0L;
	lx->w1=lx->w1-(i/B)-cy;
	if (lx->w1 < 0) { lx->w1 += B;
			   cy = 1L; }
	else cy = 0L;
	lx->w2=lx->w2-cy;
	if (lx->w2 < 0) { lx->w2 += B;
			   cy = 1L; }
	else cy = 0L;
	return(cy);
	}

static INT locsub(struct loc *lx, struct loc *ly, INT cy)
{
	lx->w0=lx->w0- ly->w0 -cy;
	if (lx->w0 < 0L) { lx->w0 += B;
			   cy = 1L; }
	else cy = 0L;
	lx->w1=lx->w1- ly->w1- cy;
	if (lx->w1 < 0L) { lx->w1 += B;
			   cy = 1L; }
	else cy = 0L;
	lx->w2=lx->w2- ly->w2- cy;
	if (lx->w2 < 0L) { lx->w2 += B;
			   cy = 1L; }
	else cy = 0L;
	return(cy);
}

static INT locvgl(struct loc *lx, struct loc *ly)
	{
	  if (lx->w2 > ly->w2) return(1L);
	    else if (lx->w2 < ly->w2 ) return(-1L);
	    else if (lx->w1 > ly->w1) return(1L);
	    else if (lx->w1 < ly->w1) return(-1L);
	    else if (lx->w0 > ly->w0) return(1L);
	    else if (lx->w0 < ly->w0) return(-1L);
	    else return(0L); 
	}  /* Ende locvgl */



static INT ganzadd(struct longint *x, struct longint *y)
	{
	struct loc *alocx, *alocy, *lloc, *plocx, *plocy;
	INT	cy,xl,ll;
	signed char xs,ys;

	xs = x->signum; ys = y->signum; xl = x->laenge;
	if (((xs>=(signed char)0) && (ys>=(signed char)0)) || 
			((xs<(signed char)0) && (ys<(signed char)0)))
		{ alocx = x->floc; alocy = y->floc; cy = 0;
		do	{ cy = locadd(alocx,alocy,cy);
			plocx = alocx; plocy = alocy; alocx = alocx->nloc;
			alocy = alocy->nloc; }
		while ((alocx != NULL) && (alocy != NULL));

		/* fuege rest an */
		if (alocy != NULL)
			{ do 	{ lochole(&alocx); plocx->nloc = alocx;
				xl++; cy = locadd(alocx,alocy,cy);
				plocx = alocx; alocx = NULL; plocy = alocy;
				alocy = alocy->nloc; }
			while (alocy != NULL);
			}
		else	{ while ((alocx != NULL) && (cy != 0))
				{ cy = locsadd(alocx,cy); plocx = alocx;
				alocx = alocx->nloc; }
			}

		/* noch ein cy? */
		if (cy != 0)
			{ lochole(&alocx);
			plocx->nloc = alocx; locint(alocx,cy); xl++; }
		if (xs == 0) xs = ys;
		} /* end of first if */
	else	{
		alocx = x->floc; alocy = y->floc; cy = 0;
		/* subtract y from x */
		do	{ cy = locsub(alocx,alocy,cy); plocx = alocx;
			alocx = alocx->nloc; plocy = alocy; alocy = alocy->nloc;
			}
		while ((alocx != NULL) && (alocy != NULL));

		/* append the remaining part */
		if (alocy != NULL)
			{
			do 	{ lochole(&alocx); plocx->nloc = alocx; xl++;
				cy = locsub(alocx,alocy,cy); plocx = alocx;
				alocx = NULL;plocy = alocy;alocy = alocy->nloc;
				}
			while (alocy != NULL);
			}
		else 	{
			while	((alocx != NULL) && (cy != 0))
				{ cy = locssub(alocx,cy);
				plocx = alocx; alocx = alocx->nloc; }
			};

		/* normieren  von x */
		if (cy != 0)
			{
			alocx = x->floc; lloc = NULL; ll = 1; cy = 0;
			do 	{ cy = locneg(alocx,cy);
				if (locsgn(alocx) != 0)
					{ lloc = alocx; xl = ll; }
				alocx = alocx->nloc; ll++;
				}
			while (alocx != NULL);
			loclisterette(&(lloc->nloc)); xs = -xs;
			if (xs == 0) xs = -1;
			}
		else	{
			alocx = x->floc; lloc = NULL; ll = 1;
			do	{
				if (locsgn(alocx) != 0)
					{ lloc = alocx; xl = ll; }
				alocx = alocx->nloc; ll++;
				}
			while (alocx != NULL);
			if (lloc == NULL)
				/* das ergebnis der addition ist null */
				{ loclisterette(&(x->floc->nloc));
				xs = 0; xl =1; }
			else	loclisterette(&(lloc->nloc));
			}
		}
	x->laenge = xl; x->signum = xs;
	return(OK);			
	}


static INT ganzsquores(struct longint *x, INT *rest, INT y)
	{

	struct loc *alocx, *blocx, *slocx;
	INT	r;
	signed char sx,sy;

	sx = x->signum;
	if (y>0L) sy=(signed char)1; 
	else if (y<0L) sy = (signed char)-1; 
	else sy=(signed char)0;
	if (y<0L) y = -y;
	blocx = x->floc; x->floc = NULL; locrezliste(&blocx);
	alocx = blocx; slocx = alocx->nloc; r=0L;
	while (slocx != NULL)
		{ r = locsdiv(alocx,r,alocx,y); alocx = slocx;
		slocx = alocx->nloc; }


	r = locsdiv(alocx,r,alocx,y);
	*rest = r * sx;
	if (locsgn(blocx) != 0L) x->signum = sx*sy;
	else if (x->laenge == 1L) x->signum = (signed char)0;
	else	{
		alocx = blocx; blocx = blocx->nloc; alocx->nloc =  NULL;
		locrette(&alocx); x->laenge --; x->signum = sx*sy;
		};

	locrezliste(&blocx);
	x->floc = blocx;
	return(OK);
	}


static INT ganzquores(struct longint *x, struct longint *rest, struct longint *y)
{
	INT		vgl,cy,cyn,a,i,rl=0L,ql;
	signed char 		sx,sy;
	struct loc 	*alocx, *plocx,*slocx,*blocx,*rlocx,*llocx,
			*alocy,*plocy,*blocy,*blocq,
			*locx2,*locx1,*locx0,*locy1,*locy0;

	struct loc	null,q,r,ov,hi,lo;
	INT		fertig;


if ((x->floc == y->floc) || (x->floc == rest->floc) || (y->floc == rest->floc))
	error("ganzquores: (1) equal variables");

loclisterette(&rest->floc);
sx = x->signum;
sy = y->signum;
if (y->laenge == 1L)	/* einfache divison */
	{
	locnull(&null); LOCASS(&lo,y->floc); blocx = x->floc;
	x->floc = NULL; locrezliste(&blocx); alocx = blocx; slocx = alocx->nloc;
	LOCASS(&r,&null);
	while (slocx != NULL)
		{
		locdiv(alocx,&r,alocx,&lo);
		alocx = slocx;
		slocx = slocx->nloc;
		}
	locdiv(alocx,&r,alocx,&lo);
	/* Seite 6 von test.p */
	if (locsgn(&r) == 0L) rest->signum = (signed char)0; 
	else rest->signum = sx;
	lochole(&rest->floc);
	LOCASS(rest->floc,&r);
	rest->laenge = 1L;
	if (locsgn(blocx) != 0L) x->signum = sx * sy;
	else if (x->laenge == 1L) x->signum = (signed char)0;
	else	{ alocx = blocx; blocx = blocx->nloc; alocx->nloc = NULL;
		locrette(&alocx); x->laenge --; x->signum = sx * sy;
		}
	locrezliste(&blocx);
	x->floc = blocx;
	} /* ende der einfachen division */
else	if (x->laenge < y->laenge)	/* trivial */
	{
	*rest = *x; x->floc = NULL; lochole(&x->floc);
	x->signum = (signed char)0; x->laenge = 1L;
	}	/* ende des trivialfalles */
else	{	/* normalfall x->laenge >= y->laenge >= 2 */
		/* lange division */
	locnull(&null); blocy = y->floc; y->floc = NULL;
	locrezliste(&blocy); locy1 = blocy; locy0 = blocy->nloc;
	a = LOCBAS2() - 1L - locms1(locy1);
	locx1 = x->floc; x->floc = NULL; locrezliste(&locx1);
	locx2 = NULL; lochole(&locx2); locx2->nloc = locx1;
	locx0 = locx1->nloc;

	/* dividend und divisor normieren. dividend zerlegen */
	locpsl(locx2,locx1,a);
	alocy = locy0; plocy = locy1; alocx = locx0; plocx = locx1;
	do	{
		locpsl(plocy,alocy,a); locpsl(plocx,alocx,a);
		plocy = alocy; alocy = alocy->nloc;
		plocx = alocx; alocx = alocx->nloc;
		}
	while (alocy != NULL);
	locpsl(plocy,&null,a);

	llocx = plocx;
	rlocx = alocx;

	while (alocx != NULL)	/* rest des dividenden normieren */
		{
		locpsl(plocx,alocx,a);
		plocx = alocx; alocx = alocx->nloc;
		}
	locpsl(plocx,&null,a);

	llocx->nloc = NULL; 	/* dividend getrennt */

	/* listen fuer teildividend und divisor umkehren */
	blocx = locx2; locrezliste(&blocx); locrezliste(&blocy);

	/* quotientenliste mit laenge */
	blocq = NULL; ql = 0L;

	do	{	/* divisionsschritt */
		if (locvgl(locx2,locy1) == 0L) LOCMAX(&q);
		else	{
			LOCASS(&r,locx2);
			locdiv(&q,&r,locx1,locy1);
			locmul(&hi,&lo,&q,locy0);
			/* falls (hi,lo) <= (r,locx0):fertig */
			vgl = locvgl(&hi,&r);
			if ((vgl >0) || ((vgl == 0L) && 
						(locvgl(&lo,locx0) > 0L)))
				{
				locssub(&q,1L);
				cy = locadd(&r,locy1,0L);
				if (cy == 0L)
					{
					cy = locsub(&lo,locy0,0L);
					if (cy == 1L) cy = locssub(&hi,1L);
			/* seite 7 von test.p */
					vgl = locvgl(&hi,&r);
					if (
						(vgl > 0L) ||
					 ((vgl == 0L) 
				&& 
	/* bug 050790 */	  (locvgl(&lo,locx0) > 0L )))
						cy = locssub(&q,1L);
					}
				}
			};

	/* subtrahiere q*divisor von teildivdend llocx = vorgaenger locx0 */
		alocy = blocy; alocx = blocx; cy = 0; cyn = 0; locnull(&ov);
		llocx = NULL; plocx = NULL;
		do	{
			locmul(&hi,&lo,alocy,&q);
			cy = locadd(&lo,&ov,cy);
			LOCASS(&ov,&hi);
			cyn = locsub(alocx,&lo,cyn);
			plocx = alocx; alocx = alocx->nloc; alocy = alocy->nloc;
			if (alocx == locx0) llocx = plocx;
			}
		while (alocy != NULL);
		cy = locsadd(&ov,cy); cyn = locsub(alocx,&ov,cyn);
		if (cy != 0L) error("ganzquores:(2) cy != 0");

		/* falls differenz negativ, q war um 1zu gross. korrektur */

		if (cyn == 1L)
			{
			cyn = locssub(&q,1L);
			alocx = blocx; alocy = blocy; cy = 0L;
			do	{
				cy = locadd(alocx,alocy,cy);
				alocx = alocx->nloc; alocy = alocy->nloc;
				}
			while (alocy != NULL);
			cy = locsadd(alocx,cy); 
			if (cy != 1L)  error("ganzquores:(3) cy != 1");
			}

		/* quotientenziffer q abspeichern . locx2 ist frei dafuer */
		locx1->nloc = NULL;
		if ((blocq == NULL) && (locsgn(&q) == 0L)) locrette(&locx2);
		else	{
			locx2->nloc = blocq; blocq = locx2; locx2 = NULL;
			LOCASS(blocq,&q); ql ++;
			};

		/* neuer teildividend */
		fertig = (rlocx == NULL);
		if (! fertig)
			{
			alocx = blocx; blocx = rlocx; rlocx = rlocx->nloc;
			blocx->nloc = alocx;
			locx2 = locx1; locx1 = locx0; locx0 = llocx;
			if (locx0 == NULL) locx0 = blocx; 
			}
		}
	while (! fertig);	/* ende divisionsschritt */

	/* quotient */
	if (blocq == NULL)
		{
		lochole (&x->floc); x->signum = (signed char)0;
		x->laenge = 1L;
		}
	else	{
		x->floc = blocq; blocq = NULL;
		x->signum = sx * sy; x->laenge = ql;
		}

	/* rest normierung rueckgaengig machen fuehrende nullen entfernen */
	i = 0L; llocx = NULL;
	plocx = blocx; alocx = plocx->nloc; 
	do	{
		i++;
	/* Seite 8 von test.p */
		locpsr(alocx,plocx,a);
		if (locsgn(plocx) != 0L) 
			{
			llocx = plocx; rl = i;
			}
		plocx = alocx; alocx = alocx->nloc;
		}
	while (alocx != NULL);
	locpsr(&null,plocx,a);
	if (locsgn(plocx) != 0L) { llocx = plocx; rl = i+ 1L; }

	if (llocx == NULL)	/* rest 0 */
		{
		loclisterette(&blocx->nloc); rest->floc = blocx;
		blocx = NULL; rest->signum = (signed char)0; rest->laenge = 1L;
		}
	else	{
		loclisterette(&llocx->nloc); rest->floc = blocx;
		blocx = NULL; rest->signum = sx; rest->laenge = rl;
		}

	/* divisor. normierung rueckgaengig machen */
	plocy = blocy; alocy = plocy->nloc;
	do	{ locpsr(alocy,plocy,a); plocy = alocy; alocy = alocy->nloc; }
	while (alocy != NULL);
	locpsr(&null,plocy,a); y->floc = blocy; blocy = NULL;
	} /* lange divison */
return(OK);
} /* ende ganzquores */

	
static INT ganzganzdiv(struct longint *x, struct longint *y)
	{
	struct longint rest;

	rest.floc = NULL; ganzquores(x,&rest,y); 
	return ganzloesche(&rest);
	}

static INT ganzmod(struct longint *x, struct longint *rest, struct longint *y)
	{
	return ganzquores(x,rest,y);
	}

static INT ganzein(FILE *fp, struct longint *x)
	{
	INT i;
	signed char sgn=(signed char)1;
	char c;
	

	fscanf(fp,"%ld",&i);
	if (i < 0L) 
		{
		sgn = (signed char)-1;
		i = i * -1L;
		}
	ganzint(x,  i % gd.basis);
	while ((c=getc(fp)) == gd.folgezeichen)
		{
		fscanf(fp,"%ld",&i);
		if (i < 0L) 
			{
			return error("ganzein:i < 0");
			}
		ganzsmul(x,gd.basis);
		ganzsadd(x,i % gd.basis);
		}

	x->signum = sgn;
	return OK;
	}


static INT holeziffer(struct zahldaten *zd)
	{
	struct loc *adez;
	INT zzmod3,erg = OK;

	zd->ziffernzahl --;
	zzmod3 = zd->ziffernzahl % 3L;

	if (zzmod3 == 0L) erg=zd->fdez->w0;
	if (zzmod3 == 1L) erg=zd->fdez->w1;
	if (zzmod3 == 2L) erg=zd->fdez->w2;

	if (zzmod3 == 0L) 
		{ adez = zd->fdez; zd->fdez = zd->fdez->nloc;
		adez->nloc = NULL; locrette(&adez); }
	
	return(erg);
	}

static INT ganzfziffer(struct zahldaten *zd, struct ganzdaten *d)
	{
	INT z,f0;
	char buffer[200];

	if (zd->ziffernzahl == 0)
		{ zd->mehr = FALSE; strcpy(zd->ziffer," "); }
	else	{
		z = holeziffer(zd /* ,d */ ); /* AB 161193 */
		if (zd->ziffernzahl > 0) zd->mehr=TRUE; else zd->mehr=FALSE;
		sprintf(buffer,"%ld",z);
		f0 = d->basislaenge-strlen(buffer);
		sprintf(zd->ziffer,"%s","000000000000");
			/* max. 12 Nullen */
		sprintf(zd->ziffer + f0,"%ld",z);
		if (zd->mehr == TRUE)
			sprintf(zd->ziffer,"%s%c",
				zd->ziffer,d->folgezeichen);
		else
			sprintf(zd->ziffer,"%s%c",zd->ziffer,' ');
		}
	return(OK);
	}

static INT retteziffer(INT z, struct zahldaten *zd)
	{
	struct loc *adez;
	INT zzmod3;

	zzmod3 = zd->ziffernzahl % 3L;

	if (zzmod3 == 0L) {
		adez = NULL; lochole(&adez);
		adez ->nloc = zd->fdez; zd->fdez = adez; }
	if (zzmod3 == 0L) zd->fdez->w0 = z;
	if (zzmod3 == 1L) zd->fdez->w1 = z;
	if (zzmod3 == 2L) zd->fdez->w2 = z;

	zd->ziffernzahl ++; return(OK);
	}

static INT ganz1ziffer(struct zahldaten *zd, struct longint	*x)
	{
	INT z;
	signed char sgn;
	struct longint xx;

	zd->fdez =  NULL; zd->ziffernzahl = 0L; xx.floc = NULL;
	lochole(&xx.floc); ganzkopiere(&xx,x); sgn = xx.signum;
	if (xx.signum < (signed char)0) xx.signum = -xx.signum;

	while (xx.signum > (signed char)0)
		{ ganzsquores(&xx,&z,gd.basis); retteziffer(z,zd /* ,&gd */ ); } /* AB 161193 */
	if (zd->ziffernzahl == 0L)
		{ zd->mehr = FALSE; strcpy(zd->ziffer," "); }
	else	{
		z = holeziffer(zd /* ,&gd */ ); z = sgn * z; /* AB 161193 */
		zd->mehr = (zd->ziffernzahl > 0L);
		if (zd->mehr == TRUE)
		  sprintf(zd->ziffer,"%s%ld%c",zd->ziffer,z,gd.folgezeichen);
		else	
		  sprintf(zd->ziffer,"%s%ld",zd->ziffer,z);
		}
	locrette(& xx.floc);
	return(OK);
	}

static INT ganzaus(FILE *fp, struct longint *x)
	{
	struct zahldaten zd;
	char *blanks = (char *) my_calloc(201 * sizeof(char), "lo.C: ganzaus");
	char *zeile  = (char *) my_calloc(201 * sizeof(char), "lo.C: ganzaus");
	INT	i;

	for (i=1;i<gd.auspos;i++) blanks[i-1]=' ';
	blanks[(gd.auspos)-1]='\0';

	zd.ziffer[0] = '\0';

	zeile[0]='\0';
	gd.auszz = 0;

	ganz1ziffer(&zd,x);
	strcat(zeile,zd.ziffer);

	while (zd.mehr == TRUE)
		{
		ganzfziffer(&zd,&gd);
		if ((INT)strlen(zeile) + (INT)strlen(zd.ziffer) > gd.auslaenge)
			{ fprintf(fp,"%s%s\n",blanks,zeile);
			strcpy(zeile,zd.ziffer); gd.auszz++; }
		else	strcat(zeile,zd.ziffer);

		}
	

	fprintf(fp,"%s%s",blanks,zeile);
		
	gd.auszz++; my_free(blanks); my_free(zeile); return(OK);
	}

static INT ganzmul(struct longint *x, struct longint *y)
/* x = x * y */
	{
	struct loc *alocx, *alocy,   *ploca, *floca, *bloca, *aloca;
	struct loc hi,lo,ov;
	INT	cy,cya;
	INT 	hh; /* fuer LOCADD,LOCMUL */



	x->signum = x->signum * y ->signum;
	if (x->signum == (signed char)0)
		{ 
		loclisterette(& (x->floc->nloc));
		locnull(x->floc); x->laenge = 1L; return OK;
		/* das ergebnis ist null */ }

	/* das ergebnis ist nicht null */
	x->laenge = x->laenge + y->laenge;
	floca = NULL; lochole(&floca); bloca = floca; alocx = x->floc->nloc;
	ploca = floca; aloca = NULL;
	while (alocx != NULL)
		{ lochole(&aloca); ploca->nloc = aloca;
		ploca = aloca; aloca = NULL; alocx = alocx->nloc; }

	alocy = y->floc;

	do	{
		cya = 0L; locnull(&ov); cy = 0L; alocx = x->floc; aloca = bloca;
		do	{ LOCMUL(&hi,&lo,alocx,alocy); 
			cy = LOCADD(&lo,&ov, cy);
			ov = hi; cya = LOCADD(aloca,&lo,cya);
			alocx=alocx->nloc; ploca=aloca; aloca = aloca->nloc; }
		while (alocx !=  NULL);
		cy = locsadd(&ov, cy+cya);
			/* cy ist jetzt 0 */
			if (cy != 0) return error("ganzmul:cy <> 0");
		lochole(&aloca);
		ploca->nloc = aloca;
		LOCASS(aloca,&ov);
		bloca = bloca->nloc;
		alocy = alocy->nloc;
		}
	while (alocy != NULL);

	if (locsgn(aloca ) == 0L)
		{
		locrette(&(ploca->nloc));
		x->laenge --;
		}
	loclisterette(&x->floc);
	x->floc = floca;
	return OK;
	}

static INT ganzsmul(struct longint *x, INT a)
{
	struct loc *alocx, *plocx;
	INT	ue;

if (a==0L)
	{
	loclisterette(&(x->floc->nloc));
	x->signum = (signed char)locint(x->floc,0L);
	x->laenge = 1L;
	}
else if (a == 1L) ;
else if (a == -1L) x->signum = - x->signum;
else	{
	if (a<0L) x->signum = - x->signum;
	alocx = x->floc; plocx = NULL; if (a<0L) a = -a;
	ue = 0L;
	do	{ ue = locsmul(alocx,a,ue);
		plocx = alocx; alocx = alocx ->nloc; }
	while (alocx != NULL);
	if (ue != 0L)
		{ lochole(&alocx); plocx->nloc = alocx; x->laenge ++;
		ue = locint(alocx,ue); }
	}
return OK;
}


static INT ganzsadd(struct longint *x, INT y)
	{
	INT cy,xl,ll;
	signed char xs,ys;
	struct loc *lloc,*alocx,*plocx;

	xl = x->laenge; 
	xs = x->signum;
	if (y>0L) 
		ys=(signed char)1; 
	else if (y<0L) 
		ys = (signed char)-1; 
	else 
		ys=(signed char)0;
	if (y<0L) y = -y;

	if (	((xs>=(signed char)0)&&(ys>=(signed char)0))
		||
		((xs<(signed char)0)&&(ys < (signed char)0))
	   )
		{
		alocx = x->floc; 
		cy = y;
		while ((alocx != NULL)&&(cy != 0L))
			{
			cy = locsadd(alocx,cy); 
			plocx = alocx;
			alocx = alocx->nloc;
			}
		if (cy != 0L) 
			{
			lochole(&alocx); 
			plocx->nloc = alocx;
			locint(alocx,cy);
			xl ++;
			}
		if (xs == (signed char)0) 
			xs = ys;
		}
	else	{
		alocx = x->floc;
		cy = y;
		while ((alocx != NULL) && (cy != 0L))
			{ 
			cy = locssub(alocx,cy);
			plocx = alocx; 
			alocx = alocx->nloc; 
			}
		if (cy != 0L)
			{
			alocx = x->floc; 
			lloc = NULL; 
			ll = 1L; 
			cy = 0L;
			do	{ 
				cy = locneg(alocx,cy);
				if (locsgn(alocx) != 0L) /* changed to 0L from NULL AB 051093 */
					{ lloc = alocx; xl = ll; }
				alocx = alocx->nloc;
				ll ++;
				}
			while (alocx != NULL);
			loclisterette(&(lloc->nloc));
			xs = -xs;
			if (xs == (signed char)0) 
				xs = (signed char)-1;
			}
		else 	{ 
			alocx = x->floc; 
			lloc = NULL; 
			ll = 1L;
			do	{
				if (locsgn(alocx) != 0L)
					{ 
					lloc = alocx; 
					xl = ll; 
					}
				alocx = alocx->nloc;
				ll ++;
				}
			while (alocx != NULL);
			if (lloc == NULL)
				{
				loclisterette(&(x->floc->nloc));
				xs = (signed char)0; 
				xl = 1L;
				}
			else 	loclisterette(&(lloc->nloc));
			}
		}
	x->laenge = xl; 
	x->signum = xs;
	return(OK);
	}
			


static INT ganzvergleich(struct longint *x, struct longint *y)
	{
	struct loc *alocx, *alocy;
	INT	av,lv;
	signed char sx,sy;

	sx = x->signum; sy = y->signum;
	if (sx>sy) 
		return(1L);
	if (sx<sy) 
		return(-1L);
	/* es gilt nun gleiches vorzeichen */
	if (sx==(signed char)0) 
		return(0L);
	if (x->laenge > y->laenge) 
		return((INT)sx);
	if (x->laenge < y->laenge) 
		return((INT)-sy);
	/* es gilt nicht nur gleiches vorzeichen sondern auch
	gleiche laenge */
	alocx = x->floc;
	alocy = y->floc;
	lv = 0;
	do	{ av = locvgl(alocx,alocy);
		if (av != 0) lv = av;
		alocx = alocx->nloc; alocy = alocy->nloc; }
	while (alocx != NULL);
	if (sx>(signed char)0) 
		return(lv); 
	else 
		return(-lv);
	}

static INT intganz(struct longint *x)
/* umwandlung longint to int falls moeglich sonst fehler */
	{
	INT wert;
	if (x->laenge != 1L) {
		error("intganz: no integer value ");
		}
	wert = x->floc->w0;
	wert += x->floc->w1 * B;
	return( (INT)x->signum * wert );
	}

static INT ganzint(struct longint *x, INT i)
	{
	if (x->floc->nloc != NULL) loclisterette(& x->floc->nloc );
	x->laenge = 1L;
	x->signum = (signed char)locint(x->floc,i);
	return(OK);
	}


static INT ganzsignum(struct longint *x)
	{
	return (INT)x->signum;
	}

static INT ganzeven(struct longint *x)
	{
	return ! locodd(x->floc);
	}

static INT ganzodd(struct longint *x)
	{
	return locodd(x->floc);
	}

static INT ganzkopiere(struct longint *x, struct longint *a)
/* x:= a AK umgeschrieben in C */
	{
	struct loc *alocx, *plocx, *aloca;

	if (x->floc == a->floc)
		 return error("ganzkopiere: identic lists ");

	x->signum = a->signum; x->laenge = a->laenge;
	aloca = a->floc; alocx = x->floc; plocx = NULL;

	do	{ if (alocx == NULL)
			{ lochole(&alocx); plocx->nloc = alocx; }
		LOCASS(alocx,aloca);
		aloca = aloca->nloc; plocx = alocx; alocx = alocx->nloc;
		}
	while (aloca != NULL);
	loclisterette(&(plocx->nloc));return(OK);
	}


static INT ganzneg(struct longint *x)
	{
	x->signum = -x->signum;return(OK);
	}


static INT lochole(struct loc **aloc)
	{
	*aloc = (struct loc *) my_malloc(sizeof(struct loc), "lo.C: lochole");
	if (*aloc == NULL) error("lochole:no memory");
	locnull(*aloc);
	(*aloc)->nloc = NULL;
	return(OK);
	}

static INT loclisterette(struct loc **aloc)
/* AK 150290 ohne struct ganzdaten */
	{
	struct loc *aloc1, *ploc1;
	if (*aloc != NULL)
		{
		aloc1=   (*aloc);
		do { ploc1 = aloc1->nloc; FREE_LOC(aloc1); aloc1 = ploc1; }
		while 	(aloc1 != NULL);
		*aloc = NULL; 
		}
	return(OK);
	}


static INT locrette(struct loc **aloc)
	{
	if (*aloc != NULL)
		{
		/* kommentar AK 100190 */
		/*
		if ((*aloc)->nloc != NULL) error("locrette:next != NULL");
		*/
		FREE_LOC(*aloc);
		*aloc = NULL;
		}
	return(OK);
	}

static INT locrezliste(struct loc **aloc)
/* dreht liste um */
	{
	struct loc *lloc,*rloc,*hloc;
	if (*aloc != NULL)
		{
		lloc = NULL; rloc = *aloc;
		while (rloc != NULL)
			{
			hloc = rloc->nloc;
			rloc->nloc=lloc;
			lloc=rloc;
			rloc=hloc;
			}
		*aloc = lloc;
		}
	return(OK);
	}

INT start_longint(void)
{
	ganzanfang(&gd);
	ganzparam(1000000L,2L,78L,'.',&gd);
	return(OK);
}

static struct longint * calloclongint(void)
{
	struct longint *ergebnis = 
		(struct longint *) my_malloc( sizeof(struct longint), "lo.C: calloclongint");
	if (ergebnis == NULL) error("calloclongint: no mem");
	return ergebnis;
}

static INT ganzparam(INT basis, INT auspos, INT auslaenge, char folgezeichen,
	struct ganzdaten *d) 
	{
	if (basis>1L) d->basis=basis; else error("ganzparam:basis negativ");

	if (auspos>1L) d->auspos=auspos; else d->auspos = 2;
	if (basis <= 10L) d->basislaenge = 1; else
	if (basis <= 100L) d->basislaenge = 2; else
	if (basis <= 1000L) d->basislaenge = 3; else
	if (basis <= 10000L) d->basislaenge = 4; else
	if (basis <= 100000L) d->basislaenge = 5; else
	if (basis <= 1000000L) d->basislaenge = 6; else
	if (basis <= 10000000L) d->basislaenge = 7; else
	if (basis <= 100000000L) d->basislaenge = 8; else
	if (basis <= 1000000000L) d->basislaenge = 9; else
				  d->basislaenge = 10;
	if (auslaenge > d->basislaenge) d->auslaenge = auslaenge;
	else d->auslaenge = d->basislaenge+1;
	d->folgezeichen = folgezeichen;return(OK);
	}

static INT ganzanfang(struct ganzdaten *d)
	{
	d->vorrat = 0L; d->speicher = 0L; d->floc = NULL; d->auszz = 0L;
	d->basis = 1000000L; d->basislaenge = 6L; d->folgezeichen = '.';
	d->auspos = 2L; d->auslaenge = 78L; return(OK);
	}

static INT ganzdefiniere(struct longint *x)
	{
	x->signum = (signed char)0; x->laenge = 1L; x->floc = NULL;
	lochole(&x->floc);
	return(OK);
	}

static INT ganzloesche(struct longint *x)

	{
	loclisterette(&x->floc);
	x->laenge = 0L; x->signum = (signed char)0;
	return(OK);
	}

INT LONGINT_OB::init()
/* init_longint */
{
	OBJECTSELF c;

	c.ob_longint = calloclongint();
	b_ks(LONGINT, c);
	c = s_obj_s();
	ganzdefiniere(c.ob_longint);
	return(OK);
}

INT LONGINT_OB::fprint(FILE *f)
/* fprint_longint */
{
	OBJECTSELF c;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::fprint() no LO obj");
		return ((SYM_OP)this)->fprint(f);
		}
	c = s_obj_s();
	ganzaus(f, c.ob_longint /* ,&gd */); /* AB 161193 */
	return(OK);
}

INT LONGINT_OB::copy(LONGINT_OP b)
/* copy_longint */
{
	OBJECTSELF bs, as;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::copy() no LO obj");
		return ((SYM_OP)this)->copy(b);
		}
	if (!b->emptyp())
		((SYM_OP)b)->freeself();
	((SYM_OP)b)->init(LONGINT);
	bs = b->s_obj_s();
	as = s_obj_s();
	return ganzkopiere(bs.ob_longint, as.ob_longint);
	/* erstes element wird gleich dem zweiten */
}

INT LONGINT_OB::freeself()
/* freeself_longint */
{
	OBJECTSELF c;
	
	if (s_obj_k() != LONGINT)
		return error("LO::freeself: no LONGINT");
	c = s_obj_s();
	ganzloesche(c.ob_longint);
	my_free(c.ob_longint);
	return(OK);
}

INT t_longint_int(LONGINT_OP a)
/* umwandlung in INTEGER falls moeglich */
{
	OBJECTSELF cs;
	INT wert;
	
	if (a->s_obj_k() == INTEGER)
		return OK;
	if (a->s_obj_k() != LONGINT)
		return error("t_longint_int: no LONGINT");
	cs = a->s_obj_s();
	if (cs.ob_longint ->laenge == 1L)
		if (cs.ob_longint ->floc ->w2 == 0L)
		if (cs.ob_longint ->floc ->w1 < 100L)
		{
		wert = intganz(cs.ob_longint);
		((INTEGER_OP)a)->m_i(wert);
		}
	return OK;
}

INT t_int_longint(INTEGER_OP a, LONGINT_OP c)
/* umwandeln von INTEGER -> LONGINT */
{
	INT erg = OK;
	OBJECTSELF cs;
	INT av;

	if (a->s_obj_k() != INTEGER)
		return error("LO::t_int_longint(): wrong obj_k");
	av = a->s_i();
	erg += ((SYM_OP)c)->init(LONGINT); 
	cs = c->s_obj_s();
	erg += ganzint(cs.ob_longint, av);
	return erg;
}

INT LONGINT_OB::add(SYM_OP b, SYM_OP c)
/* add_longint */
{
	INT erg = OK;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::add() no LO obj");
		return ((SYM_OP)this)->add(b, c);
		}
	
	switch (b->s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH:
			erg += ((BRUCH_OP)b)->add(this, c);
			if (c->s_obj_k() == LONGINT)
				erg += t_longint_int((LONGINT_OP) c);
			break;
#endif /* BRUCHTRUE */
		case INTEGER: 
			erg += add_integer((INTEGER_OP) b, (LONGINT_OP) c);
			break;
		case LONGINT:{
			OBJECTSELF bs, cs;

			erg += copy((LONGINT_OP) c);
			cs = c->s_obj_s();
			bs = b->s_obj_s();
			erg += ganzadd(cs.ob_longint, bs.ob_longint);
				/* longinteger-addition ist x:= x+y */
			erg += t_longint_int((LONGINT_OP) c);
			break;
			}
		default:{
			c->printobjectkind();
			error("LO::add: wrong second type");
			return(ERROR);
			}
		};
	if (erg != OK)
		{
		error("LO::add: error during computation");
		}
	return erg;
}

INT LONGINT_OB::add_integer(INTEGER_OP b, LONGINT_OP c)
/* add_longint_integer */
{
	OBJECTSELF cs;
	
	if (s_obj_k() != LONGINT || b->s_obj_k() != INTEGER)
		return error("LO::add_integer(): wrong obj_k");
	copy(c);
	cs = c->s_obj_s();
	ganzsadd(cs.ob_longint, b->s_i());
		/* longinteger-addition ist x:= x+y */
	t_longint_int(c);
	return(OK);
}

INT LONGINT_OB::mult(SYM_OP b, SYM_OP c)
/* mult_longint */
{
	INT erg = OK;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::mult() no LO obj");
		return ((SYM_OP)this)->mult((SYM_OP) b, (SYM_OP) c);
		}
	
	switch (b->s_obj_k()) {
#ifdef BRUCHTRUE
		case BRUCH: 
			erg += ((BRUCH_OP)b)->mult(this, c); break;
#endif
#ifdef INTEGERTRUE
		case INTEGER: 
			erg += mult_integer((INTEGER_OP) b, (LONGINT_OP) c); break;
#endif
#ifdef MATRIXTRUE
		case MATRIX: 
			erg += ((MATRIX_OP)b)->mult_scalar(this, (MATRIX_OP) c); break;
#endif
#ifdef MONOMTRUE
		case MONOM: 
			erg += ((MONOM_OP)b)->mult_scalar(this, (MONOM_OP) c); break;
#endif
		case LONGINT: 
			erg += mult_longint(
				(LONGINT_OP) b, (LONGINT_OP) c); break;
		case POLYNOM: 
			erg += ((POLYNOM_OP)b)->mult_scalar(
				this, (POLYNOM_OP) c); break;
#ifdef SCHURTRUE
		case SCHUR: erg += mult_scalar_schur(this, b, c); break;
#endif
#ifdef CHARTRUE
		case SYMCHAR: erg += mult_scalar_symchar(this, b, c); break;
#endif
		default:
			{
			b->printobjectkind();
			return error("LO::mult:wrong second type");
			}
		};
	if (erg != OK)
		return error("LO::mult: error in computation");
	return erg;
}

INT LONGINT_OB::mult_longint(LONGINT_OP b, LONGINT_OP c)
/* mult_longint_longint */
{
	OBJECTSELF cs, bs;
	
	if (s_obj_k() != LONGINT || b->s_obj_k() != LONGINT)
		return error("LO::mult_longint(): wrong obj_k");
	copy(c);
	cs = c->s_obj_s();
	bs = b->s_obj_s();
	return ganzmul(cs.ob_longint, bs.ob_longint);
		/* longinteger-multiplikation ist x:= x*y */
}

INT LONGINT_OB::mult_integer(INTEGER_OP b, LONGINT_OP c)
/* mult_longint_integer */
{
	INT erg = OK;
	OBJECTSELF cs;

	if (s_obj_k() != LONGINT || b->s_obj_k() != INTEGER)
		return error("LO::mult_integer(): wrong obj_k");
	erg += copy(c);
	cs = c->s_obj_s();
	erg += ganzsmul(cs.ob_longint, b->s_i());
		/* longinteger-multiplikation ist x:= x*y */
	erg += t_longint_int(c);
	return erg;
}

INT LONGINT_OB::addinvers_apply()
/* addinvers_apply_longint */
{
	OBJECTSELF c;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::addinvers_apply() no LO obj");
		return ((SYM_OP)this)->addinvers_apply();
		}
	
	c = s_obj_s();
	return(ganzneg(c.ob_longint));
}

INT LONGINT_OB::addinvers(LONGINT_OP b)
/* addinvers_longint */
{
	OBJECTSELF c;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::addinvers() no LO obj");
		return ((SYM_OP)this)->addinvers(b);
		}
	copy(b);
	c = b->s_obj_s();
	ganzneg(c.ob_longint);
		/* longinteger-addinvers ist x:= -x */
	return(OK);
}

INT LONGINT_OB::add_apply(SYM_OP b)
/* add_apply_longint */
{
	INT erg = OK;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::add_apply() no LO obj");
		return ((SYM_OP)this)->add_apply(b);
		}
	switch (b->s_obj_k()) {
#if 0
#ifdef BRUCHTRUE
		case BRUCH:
			erg += add_apply_scalar_bruch(this, b); break;
#endif
#endif
		case INTEGER: erg += add_apply_integer((INTEGER_OP) b); break;
		case LONGINT: erg += add_apply_longint((LONGINT_OP) b); break;
		default: /* AK 190291 */
			{
			SYM_OB c;
			
			c = *b;
			b->c_obj_k(EMPTY);
			erg += ((SYM_OP) this)->add(&c, b);
			}
		}
	if (erg != OK) {
		error("LO::add_apply: error during computation");
		}
	return erg;
}

INT LONGINT_OB::add_apply_integer(INTEGER_OP b)
/* add_apply_longint_integer */
{
	INT erg = OK;
	SYM_OB c;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::add_apply_integer() no LO obj");
		return ((SYM_OP)this)->add_apply(b);
		}
	c = *b;
	b->c_obj_k(EMPTY);
	erg += ((SYM_OP)this)->add(&c, b);
	erg += t_longint_int((LONGINT_OP) b);
	return erg;
}

INT LONGINT_OB::add_apply_longint(LONGINT_OP b)
/* add_apply_longint_longint */
{
	INT erg = OK;
	OBJECTSELF as, bs;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::add_apply_longint() no LO obj");
		return ((SYM_OP)this)->add_apply(b);
		}
	as = s_obj_s();
	bs = b->s_obj_s();
	erg += ganzadd(bs.ob_longint, as.ob_longint);
	erg += t_longint_int(b);
	if (erg != OK)
		error("LO::add_apply_longint: error during computation");
	return erg;
}

#if 0
INT add_apply_integer_longint(OP a, OP b)
/* b = a + b */ /* b ist LONGINT, a ist INTEGER */
	{
	OBJECTSELF ls;
	ls = S_O_S(b);
	return(ganzsadd(ls.ob_longint,S_I_I(a)));
	}
#endif

INT LONGINT_OB::mult_apply(SYM_OP b)
/* mult_apply_longint */
{
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::mult_apply() no LO obj");
		return ((SYM_OP)this)->mult_apply(b);
		}
	switch (b->s_obj_k()) {
#if 0
#ifdef BRUCHTRUE
		case BRUCH: 
				mult_apply(S_B_O(b));
				kuerzen(b); 
				return(OK);
#endif
#endif
		case INTEGER: 
			return mult_apply_integer((INTEGER_OP) b);
		case LONGINT: 
			return mult_apply_longint((LONGINT_OP) b);
#ifdef MATRIXTRUE
		case MATRIX: 
			return mult_apply_matrix((MATRIX_OP) b);
#endif
		default: /* AK 190291 */
			{
			SYM_OB c;
			INT erg = OK;
			
			c = *b;
			b->c_obj_k(EMPTY);
			erg += ((SYM_OP)this)->mult(&c, b);
			return erg;
			}
		}
}

INT LONGINT_OB::mult_apply_longint(LONGINT_OP b)
/* mult_apply_longint_longint */
{
	OBJECTSELF as,bs;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::mult_apply_longint() no LO obj");
		return ((SYM_OP)this)->mult_apply(b);
		}
	if (b->s_obj_k() != LONGINT)
		return error("LO::mult_apply_longint() b not a LO obj");
	as = s_obj_s();
	bs = b->s_obj_s();
	return(ganzmul(bs.ob_longint, as.ob_longint));
}

INT LONGINT_OB::mult_apply_integer(INTEGER_OP b)
/* mult_apply_longint_integer */
{
	SYM_OB c;
	INT erg = OK;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::mult_apply_integer() no LO obj");
		return ((SYM_OP)this)->mult_apply(b);
		}
	if (b->s_obj_k() != INTEGER)
		return error("LO::mult_apply_integer() b not a INT obj");
	c = *b;
	b->c_obj_k(EMPTY);
	erg += ((SYM_OP)this)->mult(&c, b);
	return erg;
}

#if 0
INT mult_apply_integer_longint(OP a, OP b)
/* b = a * b */ /* b ist LONGINT, a ist INTEGER */
	{
	OBJECTSELF ls;
	INT erg = OK;
	ls = S_O_S(b);
	erg += ganzsmul(ls.ob_longint,S_I_I(a));
	erg += t_longint_int(b);
	return erg;
	}

#endif

#ifdef MATRIXTRUE
INT LONGINT_OB::mult_apply_matrix(MATRIX_OP b)
/* mult_apply_longint_matrix */
{
	SYM_OP z;
	INT i, erg=OK;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::mult_apply_matrix() no LO obj");
		return ((SYM_OP)this)->mult_apply(b);
		}
	z = b->s_s();
	i = b->s_hi() * b->s_li();
	for(; i > 0; i--, z++)
		erg += mult_apply(z);
	return erg;
}
#endif /* MATRIXTRUE */

INT LONGINT_OB::dec()
/* dec_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::dec() no LO obj");
		return ((SYM_OP)this)->dec();
		}
	as = s_obj_s();
	ganzsadd(as.ob_longint, -1L);
	return(OK);
}

INT LONGINT_OB::inc()
/* inc_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::inc() no LO obj");
		return ((SYM_OP)this)->inc();
		}
	as = s_obj_s();
	return(ganzsadd(as.ob_longint, 1L));
}

#if 0
INT scan_longint(OP a)
	{
	OBJECTSELF as;
	
	printeingabe("longint:");
	init(LONGINT,a);as=S_O_S(a);
	ganzein(stdin,as.ob_longint /* ,&gd */ ); /* AB 161193 */
	if (nullp(a) ) { /* AK 020889 V1.0 */
		freeself(a);
		m_i_i(0L,a); 
		}
	return(OK);
	}
#endif

INT LONGINT_OB::posp()
/* posp_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT)
		return error("LO::posp() no LO obj");
	as = s_obj_s();
	return (ganzsignum(as.ob_longint) == 1L);
}

INT LONGINT_OB::negp()
/* negp_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT)
		return error("LO::negp() no LO obj");
	as = s_obj_s();
	return(ganzsignum(as.ob_longint) == -1);
}

INT LONGINT_OB::odd()
/* odd_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT)
		return error("LO::odd() no LO obj");
	as = s_obj_s();
	return ganzodd(as.ob_longint);
}

INT LONGINT_OB::even()
/* even_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT)
		return error("LO::even() no LO obj");
	as = s_obj_s();
	return ganzeven(as.ob_longint);
}

INT LONGINT_OB::nullp()
/* nullp_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::nullp() no LO obj");
		return ((SYM_OP)this)->nullp();
		}
	as = s_obj_s();
	return (ganzsignum(as.ob_longint) == 0L);
}

INT LONGINT_OB::einsp()
/* einsp_longint */
{
	OBJECTSELF as;
	
	if (s_obj_k() != LONGINT) {
		printf("warning: LO::einsp() no LO obj");
		return ((SYM_OP)this)->einsp();
		}
	as = s_obj_s();
	if (as.ob_longint ->laenge == 1L)
	if (as.ob_longint ->signum == 1L)
	if (as.ob_longint ->floc ->w2 == 0L)
	if (as.ob_longint ->floc ->w1 == 0L)
	if (as.ob_longint ->floc ->w0 == 1L)
			return TRUE;
	return FALSE;
}

INT LONGINT_OB::m_i(INT i)
{
	INTEGER_OB int_ob;
	
	int_ob.m_i(i);
	t_int_longint(&int_ob, this);
	return OK;
}

static INT ganzaus_str(struct longint *x, BYTE *str);

INT LONGINT_OB::lo2string(STRING_OP p)
{
	BYTE s[1024];
	
	s[0] = 0;
	sprint(s);
	p->init(s);
	return OK;
}

INT LONGINT_OB::string2lo(STRING_OP p)
{
	sscan(p->s_str());
	return OK;
}

INT LONGINT_OB::sprint(BYTE *str)
/* sprint_longint */
/* appends to str. writes to maximal strlength of 200. */
{
	BYTE str1[256];
	OBJECTSELF c;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::sprint() no LO obj");
		return ((SYM_OP)this)->sprint(str);
		}
	str1[0] = 0;
	if (nullp()) {
		strcat(str, "0");
		return OK;
		}
	c = s_obj_s();
	ganzaus_str(c.ob_longint, str1);
	strcat(str, str1);
	return(OK);
}

INT LONGINT_OB::sprint_latex(BYTE *str)
{
	sprint(str);
	return OK;
}

INT LONGINT_OB::sscan(BYTE *str)
{
	INTEGER_OB int_ob;
	OBJECTSELF c;

	int_ob.m_i(0);
	t_int_longint(&int_ob, this);
	c = s_obj_s();
	ganzein_str(str, c.ob_longint);
	return(OK);
}

static INT ganzein_str(BYTE *str, struct longint *x)
	{
	INT i;
	signed char sgn=(signed char)1;
	char c;
	BYTE *ps;
	
	ps = str;
	lo_s_scan_int(&ps, &i);
	/* fscanf(fp,"%ld",&i); */
	if (i < 0L) 
		{
		sgn = (signed char)-1;
		i = i * -1L;
		}
	ganzint(x,  i % gd.basis);
	while ((c = *ps++ /* c=getc(fp) */) == gd.folgezeichen)
		{
		lo_s_scan_int(&ps, &i);
	/* fscanf(fp,"%ld",&i); */
		if (i < 0L) 
			{
			return error("ganzein:i < 0");
			}
		ganzsmul(x,gd.basis);
		ganzsadd(x,i % gd.basis);
		}

	x->signum = sgn;
	return OK;
	}


static INT lo_s_scan_int(BYTE **s, INT *i)
{
	BYTE str1[512];
	
	if (!lo_s_scan_token(s, str1))
		return FALSE;
	sscanf(str1, "%ld", i);
	return TRUE;
}

static INT lo_s_scan_token(BYTE **s, BYTE *str)
{
	BYTE c;
	INT len;
	
	while (TRUE) {
		c = **s;
		if (c == 0) {
			return(FALSE);
			}
		if (c == ' ' || c == '\t' || 
			c == '\r' || c == 10 || c == 13) {
			(*s)++;
			continue;
			}
		break;
		}
	len = 0;
	while (TRUE) {
		c = **s;
		if (c == ' ' || c == '\t' || 
			c == '\r' || c == 10 || c == 13 || c == '.') {
			break;
			}
		if (c == 0) {
			break;
			}
		if (c == '\\') {
			(*s)++;
			c = **s;
			str[len] = c;
			len++;
			}
		else {
			str[len] = c;
			len++;
			}
		(*s)++;
		}
	str[len] = 0;
	return TRUE;
}


static INT ganzaus_str(struct longint *x, BYTE *str)
	{
	struct zahldaten zd;
	char *blanks = (char *) my_calloc(201 * sizeof(char), "lo.C: ganzaus_str");
	char *zeile  = (char *) my_calloc(201 * sizeof(char), "lo.C: ganzaus_str");
	INT	i;
	
	str[0] = 0;
	for (i=1;i<gd.auspos;i++) blanks[i-1]=' ';
	blanks[(gd.auspos)-1]='\0';

	zd.ziffer[0] = '\0';

	zeile[0]='\0';
	gd.auszz = 0;

	ganz1ziffer(&zd,x);
	strcat(zeile,zd.ziffer);

	while (zd.mehr == TRUE)
		{
		ganzfziffer(&zd,&gd);
		if ((INT)strlen(zeile) + (INT)strlen(zd.ziffer) > gd.auslaenge)
			{ sprintf(eostr(str), "%s%s\n",blanks,zeile);
			strcpy(zeile,zd.ziffer); gd.auszz++; }
		else	strcat(zeile,zd.ziffer);
		
		if (strlen(str) > 200)
			goto fertig;
		}
	

	sprintf(eostr(str), "%s%s",blanks,zeile);
	if (strlen(str) > 200)
		goto fertig;

		
	gd.auszz++; 
fertig:
	my_free(blanks); my_free(zeile); return(OK);
	}


INT LONGINT_OB::compare(SYM_OP c)
/* comp_longint(a, b) */
{
	INT erg;

	if (s_obj_k() != LONGINT) {
		printf("warning: LO::compare() no LO obj");
		return ((SYM_OP)this)->sym_comp(c);
		}
	switch (c->s_obj_k())
		{
		case INTEGER:
			{
			LONGINT_OB d;
			OBJECTSELF as, ds;
			
			t_int_longint((INTEGER_OP)c, &d);
			as = s_obj_s();
			ds = d.s_obj_s();
			erg = ganzvergleich(as.ob_longint, ds.ob_longint);
			return(erg);
			}
		case LONGINT:
			{
			OBJECTSELF as, cs;
			as = s_obj_s();
			cs = c->s_obj_s();
			erg = ganzvergleich(as.ob_longint, cs.ob_longint);
			return(erg);
			}
#ifdef BRUCHTRUE
		case BRUCH:
			{
			BRUCH_OB d; 
			d.m_scalar(this);
			erg = d.sym_comp(c);
			return erg;
			}
#endif
		default: {
			c->printobjectkind();
			error("LO::compare: wrong second type");
			return(ERROR);
			}
		};
	}

INT LONGINT_OB::mod(SYM_OP e, SYM_OP c)
/* mod_longint */
/* c = a mod e */
	{
	LONGINT_OB d; 

	if (s_obj_k() != LONGINT)
		return error("LO::mod() no LO obj");
	copy(&d);
		/* noetig da sonst in a das divisionsergebnis steht */
	switch (e->s_obj_k()) {
		case INTEGER: {
			OBJECTSELF cs, ds;
	
			((INTEGER_OP)c)->m_i(0);
			cs = c->s_obj_s();
			ds = d.s_obj_s();
			ganzsquores(ds.ob_longint, &cs.ob_INT, ((INTEGER_OP)e)->s_i());
			c->c_obj_s(cs);
			return(OK);
			}
		case LONGINT: {
			OBJECTSELF es, ds, cs;

			c->init(LONGINT);
			ds = d.s_obj_s();
			cs = c->s_obj_s();
			es = e->s_obj_s();
			ganzmod(ds.ob_longint, cs.ob_longint, es.ob_longint);
			return(OK);
			}
		default: {
			c->freeself();
			e->printobjectkind();
			error("LO::mod:wrong second type");
			return(ERROR);
			}
		}
	}

INT LONGINT_OB::ganzdiv(SYM_OP e, SYM_OP c)
/* ganzdiv_longint */
/* c = a / e */
{
	INT erg = OK;
	OBJECTSELF es, cs;
	INT	rest;
	
	if (s_obj_k() != LONGINT)
		return error("LO::ganzdiv() no LO obj");
	switch (e->s_obj_k()) {
		case INTEGER:
			erg += copy((LONGINT_OP) c);
			cs = c->s_obj_s(); 
			erg += ganzsquores(cs.ob_longint, &rest, ((INTEGER_OP)e)->s_i());
			erg += t_longint_int((LONGINT_OP) c);
/* AK 201092 */		if ((((INTEGER_OP)e)->s_i() < 0) && (posp())) c->dec();
/* AK 201092 */		if ((((INTEGER_OP)e)->s_i() > 0) && (negp())) c->dec();
			break;
		case LONGINT:
			erg += copy((LONGINT_OP) c); 
			cs = c->s_obj_s(); 
			es = e->s_obj_s(); 
			erg += ganzganzdiv(cs.ob_longint, es.ob_longint);
			erg += t_longint_int((LONGINT_OP) c);
/* AK 201092 */		if ((((LONGINT_OP)e)->negp()) && (posp())) c->dec();
/* AK 201092 */		if ((((LONGINT_OP)e)->posp()) && (negp())) c->dec();
			break;
#ifdef BRUCHTRUE
		case BRUCH:
			error("LO::ganzdiv: BRUCH");
#if 0
			if (einsp(S_B_U(e)))
				{
				erg += ganzdiv_longint(S_B_O(e), c);
				erg += t_longint_int(c);
				}
			else
				erg += error("LO::ganzdiv: wrong bruch");
#endif
			break;
#endif /* BRUCHTRUE */
		default:
			{
			e->printobjectkind();
			erg += error("LO::ganzdiv: wrong second type");
			break;
			}
		};
	if (erg != OK) {
		printobjectkind();
		e->printobjectkind();
		error("LO::ganzdiv: error during computation");
		return erg;
		}
	return erg;
}

INT LONGINT_OB::quores(SYM_OP e, SYM_OP c, SYM_OP d)
/* quores_longint */
/* this = c * e + d */
/* c = a / e */
{
	if (s_obj_k() != LONGINT)
		return error("LO::quores() no LO obj");
	switch (e->s_obj_k()) {
		case INTEGER:
			{
			OBJECTSELF cs;
			INT	rest;
			
			copy((LONGINT_OP) c);
			cs = c->s_obj_s(); 
			ganzsquores(cs.ob_longint, &rest, ((INTEGER_OP)e)->s_i());
			t_longint_int((LONGINT_OP) c);
			((INTEGER_OP)d)->m_i(rest);
			return(OK);
			}
		case LONGINT:
			{
			OBJECTSELF es, cs, ds;
			copy((LONGINT_OP) c);
			copy((LONGINT_OP) d);
			cs = c->s_obj_s();
			es = e->s_obj_s();
			ds = d->s_obj_s();
			ganzquores(cs.ob_longint, ds.ob_longint, es.ob_longint);
			t_longint_int((LONGINT_OP) c);
			t_longint_int((LONGINT_OP) d);
			return(OK);
			}
		default:
			{
			e->printobjectkind();
			return error("LO::quores: wrong second type");
			}
		}
	return OK;
}

INT LONGINT_OB::invers(LONGINT_OP b)
{
	if (s_obj_k() != LONGINT) {
		printf("warning: LONGINT::invers() not a LONGINT_OB\n");
		return ((SYM_OP)this)->invers(b);
		}

	if (einsp())
		return(copy(b));
	if (negeinsp())
		return(copy(b));
	
#ifdef BRUCHTRUE
	{
	BRUCH_OB c;
	c.m_ioiu(1, 1);
	copy((LONGINT_OP) c.s_u());
	c.swap(b);
	return OK;
	}
#else
	error("LONGINT::invers: BRUCH not available");
	printf("%ld\n", s_i());
	return(ERROR);
#endif
}


#endif /* LONGINTTRUE */

