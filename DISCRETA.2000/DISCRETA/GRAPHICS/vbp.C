/* vbp.C */

/* this file is part of the DISCRETA project
 * University of Bayreuth, Bavaria, Germany
 * Copyright Anton Betten 1995
 */


/* Plazieren eines Verbandes */

/* Alexander Thron und Martin Kellner '1992 */

#include <DISCRETA/discreta.h>

#ifdef GRAPHICS_TRUE

#include <DISCRETA/graphics.h>

/*
#define KOMMENTAR    1
#define XKDEB        1 
#define KOMMENTAR    1
#define DEBUG_ST     1
#define DEB_ABSTSCH  1
#define XKDEB        1 
#define DEB_MANABS   1
#define FORTGANG     1
#define X_REORG      1
*/


#define NOMEM -2             /* kein Speicherplatz vorhanden */

#define MAXKNOTEN 20
#define MAXKANTEN 100
#define MAXKLASSEN 100
#define MAXFARBEN 100
#define maxpfeile MAXKANTEN*2



typedef struct knoten KNOTEN;
typedef struct st_schicht SCHICHT;

struct knoten {
	LONG nummer; /* laufende Knotennummer in nl */
	LONG orbit_size; /* AB 080993 */
	LONG gewicht;
	struct knoten *nachfolger;
};

struct pfeile {
	int anfang, ende, nummer, zuruecknummer, richtung;
};

struct inhalt {
	int kapazitaet, fluss;
};

struct fluss_knoten {
	int vorgaenger, pfeilnummer, fluss, ausgelastet;
};

extern  int gf[MAXKNOTEN];
extern  struct pfeile gpf[maxpfeile];
extern  struct inhalt ginh[maxpfeile];

/* Rueckgabewerte */
#define GT_OK            0
#define GT_MAXFLUSS   -100

/* Deklarationen fuer die Verband-Plazierung */

/* Alexander Thron und Martin Kellner '1992 */

/* Strukturen fuer die topologische Sortierung */

struct st_schicht
{
	ULONG nr;
	ULONG knzahl;
	ULONG total_orbit_size; /* AB 080993 */ 
		/* the sum of all orbit sizes of knoten in this schicht */
	KNOTEN *knoten;
	struct st_schicht *naechste;
};


/* Anfang der Listen fuer GRAPH_TEILEN() */
#define LSTART   1
/* Minimale Knotennummer in GRAPH_TEILEN() */
#define MINKNR   1



/* Deklarationen der Rueckgabeparameter */
#define VP_QSERR     -11      /* Quelle und Senke waren unbestimmbar! */
#define VP_GRAPH     -12      /* Der Graph konnte nicht geteilt werden! */
#define VP_PARAM     -13      /* Unbelegte Parameter! */
#define VB_PLAZ_H



/* VB_PLAZ.C: */
int verband_plazieren(ULONG *nachfolgerliste, ULONG *orbit_size, 
	FLOAT **plazierung, FLOAT **orbit_dx);
int verband_plazieren_qs(ULONG *nachfolgerliste, ULONG *orbit_size, 
	FLOAT **plazierung, FLOAT **orbit_dx, 
	int (*qsfkt)(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten));
int x_streckung(ULONG *nachfolgerliste, ULONG *orbit_size, FLOAT *plazierung);
int t_gerichtet_ungerichtet(ULONG *ein, ULONG **aus);
int qs_abst(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten);
int qs_abstsch(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten);
int qs_manabs(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten);
int qs_manach(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten);
int qs_klsch(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten);
/* int main(void); */

/* GR_TEIL.C: */
int GRAPH_TEILEN(KNOTEN liste[], KNOTEN liste1[], KNOTEN liste2[], 
	int quelle, int senke, int korrektur, int kommentar);
int markiere_knoten(int gf[], struct pfeile gpf[], struct inhalt ginh[], 
	int kanten, int index, int nr);
int sch_korrektur(struct pfeile gpf[], int gf[], int knoten, int kanten, int start, 
	int qs_nicht, int eb, int rich, int index);
void pfeilliste_erst(KNOTEN list[MAXKNOTEN], struct pfeile pfl[maxpfeile], 
	struct inhalt inhl[maxpfeile], int *kantenzahl, int *knotenzahl);
INT sort_nachfolger(int untergr, int obergr, KNOTEN list[MAXKNOTEN]);
int anz_nachfolger(KNOTEN *lauf_z);
void quelle_senke_bestimmen(KNOTEN list[MAXKNOTEN], int *q, int *s, int knotenzahl);
void liste_trennen(KNOTEN list[MAXKNOTEN], KNOTEN list1[MAXKNOTEN],
	KNOTEN list2[MAXKNOTEN], int f[MAXKNOTEN], int knotenzahl);
void schnittkanten_zu_liste1_loeschen(int f[MAXKNOTEN], KNOTEN *zeiger, int knotenzahl);
void schnittkanten_zu_liste2_loeschen(int f[MAXKNOTEN], KNOTEN *zeiger, int knotenzahl);
void listenausgabe(KNOTEN list[MAXKNOTEN]);
int in_vektor(int zu_suchen, int feld[MAXKNOTEN],int von,int bis);


/* FF.C: */
int max_fluss(struct pfeile pf[maxpfeile], struct inhalt inh[maxpfeile], 
	int kantenzahl, int knotenzahl, int quelle, int senke, int f[MAXKNOTEN], 
	int korrektur);

extern INT (*Qsort_cmp_func)(void *p1, void *p2, LONG *res);
/* result > 0: p1 > p2
 *        = 0: p1 = p2
 *        < 0: p1 < p2 
 */
extern INT (*Qsort_swap_func)(void *p1, void *p2);
extern INT Qsort_f_ascending;
extern INT Qsort_ElementSizeof;

INT Q2sort(BYTE *arr, LONG left, LONG right);

/* end of vb_plaz.h */

#ifdef SYSTEMMAC
#include <stdlib.h> 
	/* for calloc() */
#endif

/*#define ULONG_MAX 429496729*/
#define VBP_ULONG_MAX (2L << 31)


/* Prototypen */
static int x_plaz(KNOTEN *gtliste, SCHICHT *schichten, FLOAT a, FLOAT b, 
	int (*qs)(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten),                  /* Funktion zur Quellen-Senken-Bestimmung */
	FLOAT *plazierung, FLOAT *orbit_dx);
static int knoten_plaz(SCHICHT *schichten, FLOAT a, FLOAT b, 
	FLOAT *plazierung, FLOAT *orbit_dx);
static int schicht_trennen(SCHICHT *schichten, SCHICHT **schichten1, 
	SCHICHT **schichten2, KNOTEN *gtliste1, KNOTEN *gtliste2, 
	ULONG *anz1, ULONG *anz2);
static void gtliste_free(KNOTEN *gtliste);
static void schichten_free(SCHICHT *schichten);
static SCHICHT *topol_sort(ULONG *ein, ULONG *orbit_size);
static void y_plaz(SCHICHT *schichten, FLOAT a, FLOAT b, FLOAT *plazierung);
static KNOTEN *gtliste(ULONG *ein);
static ULONG *abstandsmatrix(ULONG **umsetz_out, KNOTEN *gtliste, 
	ULONG knzahl, ULONG maxnr);
static int sch_null(ULONG *pneu, ULONG knzahl, ULONG *umsetz, SCHICHT *schichten);
static ULONG nachz(KNOTEN *kn);
static void print_schichten(SCHICHT *schichten, FLOAT *plazierung);
static void print_kn_aus_sc(SCHICHT *schichten);
static void print_liste(ULONG *ein);
static void print_gtliste(KNOTEN *ein);

/* Plazieren eines Verbandes */
int verband_plazieren(ULONG *nachfolgerliste, ULONG *orbit_size, 
	FLOAT **plazierung, FLOAT **orbit_dx)
{
  return(verband_plazieren_qs(nachfolgerliste, orbit_size, 
  	plazierung, orbit_dx, qs_abst));
}

int verband_plazieren_qs(ULONG *nachfolgerliste, ULONG *orbit_size, 
	FLOAT **plazierung, FLOAT **orbit_dx, 
	int (*qsfkt)(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten))
{
  SCHICHT *schicht;
  KNOTEN *gtl;
  int rueck;

  /* Test des Zeigers auf die Nachfolgerliste */
  if (nachfolgerliste == NULL) {
  	printf("verband_plazieren_qs()|nachfolgerliste == NULL\n");
    return(VP_PARAM);
    }

  /* Speicher fuer Plazierung holen */
  if (plazierung == NULL) {
   printf("verband_plazieren_qs()|plazierung == NULL\n");
   return(VP_PARAM);
   }
  *plazierung = (FLOAT *) my_calloc(nachfolgerliste[0] * 2 * sizeof(FLOAT), "vbp");
  if (*plazierung == NULL) {
    printf("verband_plazieren_qs()|no memory for plazierung\n");
    return(NOMEM);
    }
	if (orbit_size != NULL) { /* AB 080993 */
		if (orbit_dx == NULL) {
			printf("verband_plazieren_qs()|orbit_dx == NULL\n");
			return(VP_PARAM);
			}
		*orbit_dx = (FLOAT *) my_calloc(nachfolgerliste[0] * sizeof(FLOAT), "vbp");
		if (*orbit_dx == NULL) {
			printf("verband_plazieren_qs()|no memory for orbit_dx\n");
			return(VP_PARAM);
			}
		}

#ifdef KOMMENTAR
  printf("Nachbarschaftsliste:\n");
  print_liste(nachfolgerliste);
#endif
#ifdef FORTGANG
  printf("Y-Plazierung:\n");
#endif

  /* Verband topologisch sortieren und in der y-Koordinate plazieren */
  schicht = topol_sort(nachfolgerliste, orbit_size);
  if(schicht == NULL)
  {
    printf("verband_plazieren_qs()|error in topol_sort()\n");
    my_free(*plazierung);
    return(NOMEM);
  }
  y_plaz(schicht, 0.0, 1.0, *plazierung);

#ifdef KOMMENTAR
  printf("Y-plaziert:\n");
  print_schichten(schicht, *plazierung);
#endif
#ifdef FORTGANG
  printf("X-Plazierung:\n");
#endif

  /* Verband mittels Ford-Fulkerson in der x-Koordinate plazieren */
  /* qs_ stellt dabei eine Funktion zur Quelle-Senke-Bestimmung dar */
  gtl = gtliste(nachfolgerliste);
  if (gtl == NULL)
  {
     printf("verband_plazieren_qs()|error in gtliste()\n");
   schichten_free(schicht);
    my_free(*plazierung);
    return(NOMEM);
  }
#ifdef KOMMENTAR
	printf("Adjazenzliste:\n");
	print_gtliste(gtl);
#endif
	if (orbit_size != NULL) /* AB 080993 */
		rueck = x_plaz(gtl, schicht, (FLOAT)0.0, (FLOAT)1.0, qsfkt, *plazierung, *orbit_dx);
	else
		rueck = x_plaz(gtl, schicht, (FLOAT)0.0, (FLOAT)1.0, qsfkt, *plazierung, NULL);
	if(rueck != OK)   
  {
     printf("verband_plazieren_qs()|error in x_plaz()\n");
    my_free(*plazierung);
    return(rueck);
  }

#ifdef X_REORG
  x_streckung(nachfolgerliste, *plazierung);
#endif
#ifdef FORTGANG
  printf("Plazierung fertig!\n");
#endif

  return(0);
}
 

int x_streckung(ULONG *nachfolgerliste, ULONG *orbit_size, FLOAT *plazierung)
{
  SCHICHT *schicht, *sl;
  KNOTEN *kn;
  FLOAT min, max, neu, eps, lalt, lneu;
  
  /* Test der Zeiger */
  if(nachfolgerliste == NULL)
    return(VP_PARAM);
  if(plazierung == NULL)
    return(VP_PARAM);

  /* Schichten errechnen */
  schicht = topol_sort(nachfolgerliste, orbit_size);
  if (schicht == NULL)
  {
    return(NOMEM);
  }

  sl = schicht;
  while(sl != NULL)
  {
    if(sl->knzahl <= 1)
    {
      sl = sl->naechste;
      continue;
    }

    /* Maximum und Minimum der x-Koordinate suchen */
    min = max = 0.0;
    kn = sl->knoten;
    if(kn != NULL)
      min = max = plazierung[kn->nummer *2];
    while(kn != NULL)
    {
      neu = plazierung[kn->nummer * 2];
      if(neu > max)
	max = neu;
      else if(neu < min)
	min = neu;
      kn = kn->nachfolger;
    }

    /* Streckung */
    eps = 1.0 / (2.0 * sl->knzahl);
    lneu = 1.0 - 1.0 / sl->knzahl;
    lalt = max - min;
    kn = sl->knoten;
    while(kn != NULL)
    {
      plazierung[kn->nummer * 2] = 
      	(plazierung[kn->nummer * 2] - min) * lneu / lalt + eps;
      kn = kn->nachfolger;
    }

    sl = sl->naechste;
  }
  return(OK);
}



/* Rueckgabeparameter der Knotenplazierungsfunktion */
#define FERTIG           0
#define NICHT_FERTIG     1


/* Bestimmen der x-Koordinate eines topologisch sort. Graphen */
static int x_plaz(KNOTEN *gtliste, SCHICHT *schichten, FLOAT a, FLOAT b, 
	int (*qs)(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten), 
	FLOAT *plazierung, FLOAT *orbit_dx)
#if FALSE
KNOTEN * gtliste;     /* Teilgraph gemaess graphteilen() */
SCHICHT * schichten;         /* Schichtendarstellung */
FLOAT a,b;                   /* Rechte und Linke Koordinatengrenze */
int (*qs)(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten);                 /* Funktion zur Quellen-Senken-Bestimmung */
FLOAT * plazierung;          /* Vektor fuer Plazierungskoordinaten */
#endif
{
  KNOTEN * gtliste1, * gtliste2;
  SCHICHT * schichten1 = NULL, * schichten2 = NULL;
  ULONG quelle, senke, anz1, anz2;
  FLOAT t;
  int rueck;
  
#ifdef XKDEB
printf("----- x_plaz() in [%f,%f] -----\n", a, b);
print_schichten(schichten,plazierung);
/*printf("knoten_plaz():\n");*/
#endif
#ifdef FORTGANG
  printf("-");
  fflush(stdout);
#endif

  /* Plazierung einzelner Knoten; falls fertig mit allen Knoten: Ruecksprung */
  if (knoten_plaz(schichten, a, b, plazierung, orbit_dx)==FERTIG) 
  {
#ifdef XKDEB
printf("FERTIG!\n");
#endif
    /* Freigabe des Speicherplatzes von gtliste samt den Nachfolgern */
    gtliste_free(gtliste);
    /* ... und von schichten */
    schichten_free(schichten);
    return(OK);
  }
  /* Quelle und Senke bestimmen */
  if((rueck = (*qs)(&quelle, &senke, gtliste, schichten)) != OK)
  {
    /* Freigabe des Speicherplatzes von gtliste samt den Nachfolgern */
    gtliste_free(gtliste);
    /* ... und von schichten */
    schichten_free(schichten);
    return(rueck);
  }
#ifdef XKDEB
printf("qs: q=%ld s=%ld\n", quelle, senke);
#endif
  { /* Speicherplatz fuer Teillisten reservieren */
    int laenge = 0;
   
    gtliste1 = gtliste;
    while(gtliste1[laenge + LSTART].nummer > 0)
      laenge++;
    gtliste1 = (KNOTEN *) my_calloc((laenge + 1 + LSTART) * sizeof(KNOTEN), "vbp");
    gtliste2 = (KNOTEN *) my_calloc((laenge + 1 + LSTART) * sizeof(KNOTEN), "vbp");
    if(gtliste1 == NULL || gtliste2 == NULL)
    {
      /* Freigabe des Speicherplatzes von gtliste samt den Nachfolgern */
      gtliste_free(gtliste);
      /* ... und von schichten */
      schichten_free(schichten);
      if (gtliste1 != NULL) my_free(gtliste1);
      if (gtliste2 != NULL) my_free(gtliste2);
      return(NOMEM);  
    }
  }

#ifdef FORTGANG
  printf(">");
  fflush(stdout);
#endif
#ifdef XKDEB
printf("GRAPH_TEILEN():\n");
#endif

  /* Teilen des Graphen mit Ford-Fulkerson */
  /* schnittkorrektur = 1 funktioniert bei GRAPH_TEILEN leider nicht! */
  if(GRAPH_TEILEN(gtliste, gtliste1, gtliste2, 
  	(int)quelle + MINKNR, (int)senke + MINKNR, 1, 0) != GT_OK)
  {
    gtliste_free(gtliste);
    schichten_free(schichten);
    gtliste_free(gtliste1);
    gtliste_free(gtliste2);
    return(VP_GRAPH);
  }
  /* Freigabe des Speicherplatzes von gtliste */
  my_free(gtliste);

#ifdef FORTGANG
  printf("<");
  fflush(stdout);
#endif

  /* Trennen des Schichtenmodells */
#ifdef XKDEB
printf("schicht_trennen():\n");
#endif
  if((rueck = schicht_trennen(schichten, &schichten1, &schichten2, 
  	gtliste1, gtliste2, &anz1, &anz2)) != OK)
  {
    return(rueck);
  }
  /* schichten ist nach dem Trennen ungueltig */
  schichten = NULL;
  if(anz1 == 0L || anz2 == 0L)
  {
    schichten_free(schichten1);
    schichten_free(schichten2);
    gtliste_free(gtliste1);
    gtliste_free(gtliste2);
    return(VP_GRAPH);
  }
#ifdef FORTGANG
  printf("%ld:%ld", anz1, anz2);
  fflush(stdout);
#endif

#ifdef XKDEB
printf("Halbgraph 1:\n");
print_schichten(schichten1, plazierung);
printf("Halbgraph 2:\n");
print_schichten(schichten2, plazierung);
printf("a=%f b=%f, anz1=%ld, anz2=%ld\n", a, b, anz1, anz2);
#endif

  /* Berechnung des Intervallschnittes a < t < b */
  t = a + (FLOAT) anz1 / (anz1 + anz2) * (b - a);

#ifdef XKDEB
if(anz1 == 0L || anz2 == 0L)
exit(-1);
printf("rekursiv x_plaz(): [%f,%f] - [%f,%f]\n", a, t, t, b);
#endif

  /* Rekursiver Aufruf der Funktion mit den Teilgraphen */
  rueck = x_plaz(gtliste1, schichten1, a, t, qs, plazierung, orbit_dx);
  if(rueck != OK)
  {
    gtliste_free(gtliste2);
    schichten_free(schichten2);
    return(rueck);
  }
  rueck = x_plaz(gtliste2, schichten2, t, b, qs, plazierung, orbit_dx);
  if(rueck != OK)
    return(rueck);
  return(OK);
}


static int knoten_plaz(SCHICHT *schichten, FLOAT a, FLOAT b, 
	FLOAT *plazierung, FLOAT *orbit_dx)
/* Plaziert die Knoten aus schichten, die allein in der Schicht stehen */
/* Rueckgabe: NICHT_FERTIG = Noch nicht alle Knoten plaziert */
/*            FERTIG = Alle Knoten in schichten konnten plaziert werden */
{
  SCHICHT * sl = schichten;
  KNOTEN * kl;
  int rueck = FERTIG;
  int nb_knoten, i;
  FLOAT t, d, d1;
  INT f_use_orbit_size;

	f_use_orbit_size = (schichten->total_orbit_size > 0);
	if (f_use_orbit_size && orbit_dx == NULL) {
		printf("knoten_plaz()|orbit_dx == NULL\n");
		return (rueck);
		}
	/* Durchlaufe die Schichten */
	while (sl != NULL) {
		/* Betrachte die Knoten der Schicht */
		if (f_use_orbit_size) {
			nb_knoten = sl->total_orbit_size;
			}
		else {
			kl = sl->knoten;
			nb_knoten = 1;
			if (kl != NULL) {
				while (kl->nachfolger) {
					nb_knoten++;
					kl = kl->nachfolger;
					}
				}
			}
		if (nb_knoten == 0) {
			printf("knoten_plaz()|nb_knoten == 0\n");
			return (rueck);
			}
		t = b - a;
		i = 0;
		kl = sl->knoten;
		while (TRUE) {
			if (f_use_orbit_size)
				d = (FLOAT)kl->orbit_size / (FLOAT)nb_knoten;
			else
				d = (FLOAT)1. / (FLOAT)nb_knoten;
			
			d1 = (FLOAT)i / (FLOAT)nb_knoten;
			
			if (plazierung[kl->nummer * 2] == 0.0) {
				plazierung[kl->nummer * 2] = a + (d1 + d * .5) * t;
				}
			
			if (f_use_orbit_size) {
				i += kl->orbit_size;
				orbit_dx[kl->nummer] = d * t;
				}
			else
				i++;
			
			if (kl->nachfolger == NULL)
				break;
			kl = kl->nachfolger;
			}
		sl = sl->naechste;
		}
	return (rueck);
	
#if FALSE
		/* we do not use ford fulkerson AB 080993 */
      if(kl->nachfolger == NULL)
      {
	/* Es ist nur ein Knoten in der Schicht: Er wird mittig plaziert */
        if(plazierung[kl->nummer * 2] == 0.0)
        {
          plazierung[kl->nummer * 2] = (a + b) / 2;
        }
      }
      else   /* sonst sind mehrere Knoten in der Schicht => Ford-Fulkerson! */
        rueck = NICHT_FERTIG;
#endif
} 


static int schicht_trennen(SCHICHT *schichten, SCHICHT **schichten1, 
	SCHICHT **schichten2, KNOTEN *gtliste1, KNOTEN *gtliste2, 
	ULONG *anz1, ULONG *anz2)
/* Trennen der schichten gemaess gtliste1 und gtliste2 */
/* Liefert in anz_ die Anzahl der Knoten in schichten_ */
/* Die Knoten aus schichten werden dabei aufgeteilt auf die beiden 
   Teilschichten. Der Graph liegt danach nicht mehr als gesamtes
   Schichtenmodell in schichten vor */
{
  SCHICHT * sl = schichten, * sl1, * neu, * akt;
  KNOTEN * kl, * kl_d, * kl1;
  int i;
  
  /* Initialisieren */
  *anz1 = *anz2 = (ULONG) 0;

  /* Speicher fuer die Teilliste schichten1 holen */
  while(sl != NULL)
  {
    if((neu = (SCHICHT *) my_calloc(sizeof(SCHICHT), "vbp")) == NULL)
    {
      schichten_free(schichten);
      schichten_free(*schichten1);
      gtliste_free(gtliste1);
      gtliste_free(gtliste2);
      return(NOMEM);
    }
    if((*schichten1) != NULL)
      akt->naechste = neu;
    else
      *schichten1 = neu;
    neu->nr = sl->nr;
    neu->knzahl = 0L;
    neu->knoten = NULL;
    neu->naechste = NULL;
    akt = neu;
    sl = sl->naechste;
  } 

#ifdef DEBUG_ST
print_gtliste(gtliste1);
print_kn_aus_sc(schichten);
print_kn_aus_sc(*schichten1);
#endif

  /* Durchlaufen der Schichten und Ausschneiden der Knoten aus gtliste1. */
  /* Haenge diese Knoten (aus gtliste1) in die Struktur schichten1. */
  /* Die verbliebenen Knoten bilden die Struktur schichten2. */
  sl = schichten;
  sl1 = *schichten1;
  while(sl != NULL && sl1 != NULL)
  {
    /* Knoten der aktuellen Schicht absuchen */
    kl = sl->knoten;
    kl_d = NULL;
    kl1 = NULL;
    while(kl != NULL)
    {
#ifdef DEBUG_ST
printf("Teilgrafen [0...]:\n");
print_kn_aus_sc(schichten);
print_kn_aus_sc(*schichten1);
printf("Knoten %lu betrachten:\n", kl->nummer+MINKNR);
#endif
      /* Knoten kl untersuchen */
      for(i = LSTART; gtliste1[i].nummer > 0; ++i)
      { 
#ifdef DEBUG_ST
printf("sch_trenn : gtliste1[%lu] = %d\n",i,gtliste1[i].nummer);
#endif
        if(gtliste1[i].nummer == kl->nummer + MINKNR)
        {
	  /* Der Knoten kl ist in der Liste gtliste1. */
#ifdef DEBUG_ST
printf("%lu aus schichten nach schichten1\n",kl->nummer+MINKNR);
#endif
          /* Knoten in die Liste schichten1 bewegen */
          if(kl1 == NULL)
            sl1->knoten = kl;
          else
            kl1->nachfolger = kl;
          kl1 = kl;
          if(kl_d == NULL)
            sl->knoten = kl->nachfolger;
          else
            kl_d->nachfolger = kl->nachfolger;
          kl = kl->nachfolger;
          kl1->nachfolger = NULL;
	  sl1->knzahl++;
	  sl->knzahl--;
          (*anz1)++;
	  break;
        }  /* if() */
      }    /* for() */
      if(gtliste1[i].nummer <= 0)
      {
        /* Knoten bleibt wo er ist */
#ifdef DEBUG_ST
printf("Knoten nicht bewegen\n");
#endif
        kl_d = kl;
        kl = kl->nachfolger;
        (*anz2)++;
      }
    }   /* while(kl != NULL) */
    /* Mit der naechsten Schicht weiter */
    sl = sl->naechste;
    sl1 = sl1->naechste;
  }

  /* Der Rest bildet schichten2 */
  *schichten2 = schichten;
  return(OK);
}


static void gtliste_free(KNOTEN *gtliste)
/* Freigabe des Speicherplatzes von gtliste samt den Nachfolgern */
{
  KNOTEN *gtlistez;
  int i = 0;
  KNOTEN *alt, *neu;

  /* Durchlauf der Liste */
  gtlistez = gtliste;
  while (gtlistez[i + LSTART].nummer > 0)
  {
    /* Durchlauf der Nachfolger */
    alt = gtlistez[i + LSTART].nachfolger;
    while (alt != NULL)
    {
      neu = alt->nachfolger;
      my_free(alt);
      alt = neu;
    } 
    i++;
  }
  my_free(gtliste);
}


static void schichten_free(SCHICHT *schichten)
/* Freigabe des Speichers fuer den Graphen in schichten */
{
  SCHICHT *sa = schichten, *sn;
  KNOTEN *ka, *kn;

  /* Durchlauf der Schichten */
  while(sa != NULL)
  {
    sn = sa->naechste;
    ka = sa->knoten;
    /* Durchlauf der Knoten der Schicht */
    while(ka != NULL)
    {
      kn = ka->nachfolger;
      my_free(ka);
      ka = kn;
    }
    my_free(sa);
    sa = sn;
  }
}



/* Bestimmen der y-Koordinate eines topologisch sort. Graphen */

#define BEARBEITET    9999L


static SCHICHT *topol_sort(ULONG *ein, ULONG *orbit_size)
/* Topologisches Sortieren einer Nachbarschaftsliste */
/* Das Ergebnis ist das Schichtenmodell */
{
  ULONG i, knoten, nachbar, *va, ausgegeben, schichtentiefe;
  ULONG weg = ein[0] - 1; /* Anzahl der knoten noch zu bearbeiten */
  SCHICHT *schichten = NULL, *arbeit = NULL, *dummy;
  KNOTEN *aknoten = NULL, *kdummy;

  /* Speicher fuer die n Knoten zum Merken der Vorgaengerzahl holen */
  if (NULL == (va = (ULONG *) my_calloc((ein[0] - 1) * sizeof(ULONG), "vbp"))) 
    return(NULL);
  /* Initialisierung */
  for(i = 0; i < ein[0] - 1; ++i)
    va[i] = (ULONG) 0;

  /* Nachbarschaftsliste ein durchlaufen, um die Vorgaenger zu zaehlen */
  for (knoten = 0; knoten < ein[0] - 1; ++knoten)
    for (nachbar = ein[knoten]; nachbar < ein[knoten + 1]; ++nachbar) {
      /* ein[nachbar] hat noch einen Vorgaenger! */
      if (ein[nachbar] < 0 || ein[nachbar] >= ein[0] - 1) { 
		my_free(va);
		return(NULL);
      }
      va[ein[nachbar]]++;
    }
  /* Durchlaufe, bis alle Knoten behandelt wurden. */
  schichtentiefe = (ULONG) 0;
  do { 
    ausgegeben = (ULONG) 0;
 
    /* Schichtstruktur reservieren und initialisieren */
    if (NULL == (dummy = (SCHICHT *) my_calloc(sizeof(SCHICHT), "vbp"))) {
      my_free(va);
      schichten_free(schichten);
      return(NULL);
    }
    else {
      if (schichten == NULL)
        schichten = dummy;
      else
        arbeit->naechste = dummy;
      arbeit = dummy;
      arbeit->nr = schichtentiefe++;
      arbeit->knoten = NULL;
      arbeit->naechste = NULL;
      arbeit->knzahl = (ULONG) 0;
      arbeit->total_orbit_size = (ULONG) 0;

      /* Knoten der Schicht x herausfinden */
      for (knoten = 0; knoten < ein[0] - 1; ++knoten) {
        if (va[knoten] == 0) {  /* != BEARBEITET */
          /* Knotenstruktur reservieren */
          if (NULL == (kdummy = (KNOTEN *) my_calloc(sizeof(KNOTEN), "vbp"))) {
	    	my_free(va);
	    	schichten_free(schichten);
	    	return(NULL);
	  		}
          else {
          	/*printf("Schicht %ld knoten %ld\n", schichtentiefe - 1L, knoten);*/
	    	/* und den Knoten darin speichern */
            if (arbeit->knoten == NULL)
              arbeit->knoten = kdummy;
            else
              aknoten->nachfolger = kdummy;
            aknoten = kdummy;
            aknoten->nummer = knoten;
            aknoten->orbit_size = 0;
            if (orbit_size) {
            	aknoten->orbit_size = orbit_size[knoten];
            	arbeit->total_orbit_size += aknoten->orbit_size;
            	}
            aknoten->nachfolger = NULL;
	    	(arbeit->knzahl)++;
            va[knoten] = BEARBEITET;
            ausgegeben = (ULONG) 1;
            --weg;
          }
        }
      }

      /* Vorgaengeranzahlen der Nachfolgeknoten dekrementieren */
      aknoten = arbeit->knoten;
      while(aknoten) {
        knoten = aknoten->nummer;
        for (nachbar = ein[knoten]; nachbar < ein[knoten + 1]; ++nachbar)
          va[ein[nachbar]]--;
        aknoten = aknoten->nachfolger;
      }
    }
  } while (ausgegeben && weg > 0L);
  
  my_free(va);
  return(schichten);
}


static void y_plaz(SCHICHT *schichten, FLOAT a, FLOAT b, FLOAT *plazierung)
{
  ULONG summe = 0;
  SCHICHT * sl;
  KNOTEN * kl;
  FLOAT t, y, l, d;
  INT f_use_orbit_size;

	/* Anzahl der Knoten in jeder Schicht errechnen */
	sl = schichten;
	f_use_orbit_size = (schichten->total_orbit_size > 0);
	while (sl != NULL) {
  		if (f_use_orbit_size)
			summe += sl->total_orbit_size;
		else
  			summe += sl->knzahl;
		sl = sl->naechste;
		}

	/* Knoten plazieren in y */
	l = b - a;
	sl = schichten;
	/* Schichten durchlaufen */
	while (sl != NULL) {
  		if (f_use_orbit_size)
  			d = (FLOAT) sl->total_orbit_size / (FLOAT) summe;
  		else
  			d = (FLOAT) sl->knzahl / (FLOAT) summe;
		t = b - d * l;
		y = (b + t) / 2;
		kl = sl->knoten;
		/* Knoten dieser durchlaufen und plazieren */
		while (kl != NULL) {
			plazierung[kl->nummer * 2 + 1] = y;
			kl = kl->nachfolger;
			}
		b = t;
		sl = sl->naechste;
		}
}


static KNOTEN *gtliste(ULONG *ein)
/* Erstellen der Datenstruktur fuers Graphteilen */
/* Wandelt eine Nachfolgerliste in die Datenorganisation [KNOTEN] um */
/* Die Nachbarschaftsliste muss dabei gerichtet sein. Das Ergebnis ist eine */
/* Struktur mit ungerichteten Kanten. */
{
  KNOTEN *aus, *neu;
  ULONG i, knot, nach;

  /* Speicher fuer die Liste holen */
  if(NULL == (aus = (KNOTEN *) my_calloc((ein[0] + LSTART) * sizeof(KNOTEN), "vbp")))
    return(NULL);
  
  /* Liste vorbelegen */
  for(i = 0; i < ein[0]; ++i)
  {
    aus[i + LSTART].nummer = i + MINKNR;
    aus[i + LSTART].gewicht = 0;
    aus[i + LSTART].nachfolger = NULL;
  }
  aus[ein[0] + LSTART - 1].nummer = -1;

  /* Durchlaufen der Nachbarschaftsliste und Aufbau der neuen Liste */
  /* die ungerichtet ist */
  for(knot = 0; knot < ein[0] - 1; ++knot)
    for(nach = ein[knot]; nach < ein[knot +1]; ++nach)
    {
      /* Pfeil knot->nach in die Liste einhaengen */
      if(NULL == (neu = (KNOTEN *) my_calloc(sizeof(KNOTEN), "vbp")))
      {
	gtliste_free(aus);
	return(NULL);
      }
      neu->nachfolger = aus[knot + LSTART].nachfolger;
      neu->nummer = ein[nach] + MINKNR;
      neu->gewicht = 1;
      aus[knot + LSTART].nachfolger = neu;
      
      /* Pfeil nach->knot in die Liste einhaengen */
      if(NULL == (neu = (KNOTEN *) my_calloc(sizeof(KNOTEN), "vbp")))
      {
	gtliste_free(aus);
	return(NULL);
      }
      neu->nachfolger = aus[ein[nach] + LSTART].nachfolger;
      neu->nummer = knot + MINKNR;
      neu->gewicht = 1;
      aus[ein[nach] + LSTART].nachfolger = neu;
    }

  return(aus);  
}


int t_gerichtet_ungerichtet(ULONG *ein, ULONG **aus)
/* Umwandlung eines gerichteten in einen ungerichteten Graphen */
{
  ULONG i, j, knoten, nachbar, **nl, zahl;

  /* Speicher fuer die n Knoten und die moeglichen n Nachbarn jedes holen */
  if(NULL == (nl = (ULONG **) my_calloc((ein[0] - 1) * sizeof(ULONG *), "vbp"))) 
    return(NOMEM);

  /* Initialisierung */
  for(i = 0; i < ein[0] - 1; ++i)
  {
    if(NULL == (nl[i] = (ULONG *) my_calloc((ein[0] - 1) * sizeof(ULONG), "vbp")))
    {
      for(j = 0; j < i; ++j)
	my_free(nl[j]);
      my_free(nl);
      return(NOMEM);
    }
    *nl[i] = (ULONG) 0;
  }

  /* Gerichtete Nachbarschaftsliste ein durchlaufen */
  for(knoten = 0; knoten < ein[0] - 1; ++knoten)
    for(nachbar = ein[knoten]; nachbar < ein[knoten + 1]; ++nachbar)
    {
      /* knoten-->ein[nachbar] und ein[nachbar]-->knoten in nl vermerken */
      if(ein[nachbar] < 0 || ein[nachbar] >= ein[0] - 1)
      { 
        for(i = 0L; i < ein[0] - 1; ++i) 
          my_free(nl[i]);
        my_free(nl);
	return(NOMEM);
      }
      for(i = 1; i <= *nl[knoten]; ++i)
        if(nl[knoten][i] == ein[nachbar])
          break;
      if(i == *nl[knoten] + 1)
      {
        nl[knoten][i] = ein[nachbar];        /* knoten-->nachbar */
        (*nl[knoten])++;
      }
      for(i = 1; i <= *nl[ein[nachbar]]; ++i)
        if(nl[ein[nachbar]][i] == knoten)
          break;
      if(i == *nl[ein[nachbar]] + 1)
      {
        nl[ein[nachbar]][i] = knoten;        /* knoten<--nachbar */
        (*nl[ein[nachbar]])++;
      }
    }

  /* Ungerichtete Nachbarschaftsliste zusammensetzen */
  for(knoten = 0, zahl = ein[0]; knoten < ein[0] - 1; ++knoten)
    zahl += *nl[knoten];
  if(NULL == (*aus = (ULONG *) my_calloc(zahl * sizeof(ULONG), "vbp")))
  {
    for(i = 0L; i < ein[0] - 1; ++i) 
      my_free(nl[i]);
    my_free(nl);
    return(NOMEM);
  }
  else
  {
    /* Fuellen des Feldes mit den Listendaten */
    (*aus)[0]=ein[0];
    i=ein[0];
    for(knoten=0; knoten<ein[0]-1; ++knoten)
    {
      (*aus)[knoten+1]=i+(*nl[knoten]);
      for(nachbar=0; nachbar<*nl[knoten]; ++nachbar)
      {
        (*aus)[i++]=nl[knoten][nachbar+1];
      }
    }
  }

  /* Speicher freigeben */
  for(i = 0L; i < ein[0] - 1; ++i) 
    my_free(nl[i]);
  my_free(nl);
  return(OK);
}


/* Quelle und Senke bestimmen: Verschiedene Funktionen */


static ULONG *abstandsmatrix(ULONG **umsetz_out, KNOTEN *gtliste, 
	ULONG knzahl, ULONG maxnr)
/* Abstandsmatrix berechnen */
{
  ULONG i, j, d, k, min, t;
  ULONG *adja = NULL, *pneu = NULL, *palt = NULL;
  KNOTEN *kn;
  int fertig;
  ULONG *umsetz = NULL;

  /* Umsetzungsfeld erstellen */
  if((umsetz = (ULONG *) my_calloc((maxnr + 1) * sizeof(ULONG), "vbp")) == NULL)
    return(NULL);
  for(i = LSTART; gtliste[i].nummer > 0; ++i)
    umsetz[gtliste[i].nummer-MINKNR] = i - LSTART;

  /* Adjazenzmatrix aufstellen */
  if((adja = (ULONG *) my_calloc(knzahl * knzahl * sizeof(ULONG), "vbp")) == NULL)
  {
    my_free(umsetz);
    return(NULL);
  }
  for(i = 0L; i < knzahl; ++i) 
  {
    kn = gtliste[i+LSTART].nachfolger;
    while(kn != NULL)
    {
      adja[i*knzahl+umsetz[kn->nummer-MINKNR]] = 1L;
      kn = kn->nachfolger;
    }
  }

  /* Abstandsberechnung durch "Potenzierung" */
  if((pneu = (ULONG *) my_calloc(knzahl * knzahl * sizeof(ULONG), "vbp")) == NULL)
  {
    my_free(umsetz);
    my_free(adja);
    return(NULL);
  }
  if((palt = (ULONG *) my_calloc(knzahl * knzahl * sizeof(ULONG), "vbp")) == NULL)
  {
    my_free(umsetz);
    my_free(adja);
    my_free(pneu);
    return(NULL);
  }
  for(i = 0L; i < knzahl*knzahl; ++i)
    pneu[i] = palt[i] = adja[i];
  for(d = 0L; d < knzahl; ++d)
  {
    fertig = 1;
    for(i = 0L; i < knzahl; ++i)
      for(j = 0L; j < knzahl; ++j)
	if(pneu[i*knzahl+j] == 0L)
        {
	  min = 0L;
	  for(k = 0L; k < knzahl; ++k)
	    if(palt[i*knzahl+k] != 0L && adja[k*knzahl+j] != 0L)
	      if((t = palt[i*knzahl+k] + adja[k*knzahl+j]) < min || min == 0L)
	      {
	        min = t;
	        fertig = 0;
	      }
	  pneu[i*knzahl+j] = min;
        }
    for(i = 0L; i < knzahl*knzahl; ++i)
      palt[i] = pneu[i];

    if(fertig != 0)
      break;
  }

  /* Diagonale loeschen */
  for(i = 0L; i < knzahl; ++i)
    pneu[i*knzahl+i] = 0L;

  my_free(adja);
  my_free(palt);
  *umsetz_out = umsetz;
  return(pneu);
}


int qs_abst(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten)
/* Groesster Abstand mittels Abstandsmatrix */
{
  ULONG knzahl = 0L, i, j, min, t;
  ULONG *pneu = NULL, *umsetz = NULL;
  ULONG maxnr = 0L;

  /* Initialisierung von Quelle und Senke */
  if(gtliste[LSTART].nummer <= 0L || gtliste[LSTART + 1].nummer <= 0L)
    return(VP_QSERR);
  *quelle = (ULONG) gtliste[LSTART].nummer - MINKNR;
  *senke = (ULONG) gtliste[LSTART+1].nummer - MINKNR;

#ifdef DEBUG
printf("qs_abst():\n");
print_kn_aus_sc(schichten);
print_gtliste(gtliste);
printf("Knotenzahl bestimmen:\n");
fflush(stdout);
#endif

  /* Maximale Knotennummer und Knotenanzahl bestimmen */
  for(i = LSTART; gtliste[i].nummer > 0; ++i, ++knzahl)
    if(maxnr < gtliste[i].nummer)
      maxnr = gtliste[i].nummer;

#ifdef DEBUG
printf("knzahl=%lu\n", knzahl);
fflush(stdout);
#endif

  /* Abstandsmatrix berechnen */
  pneu = abstandsmatrix(&umsetz, gtliste, knzahl, maxnr);
  if(pneu == NULL)
    return(NOMEM);

  /* Quelle und Senke aussuchen */
  min = 0L;
  for(i = 0L; i < knzahl; ++i)
    for(j = 0L; j < knzahl; ++j)
      if((t = pneu[i*knzahl+j]) > min)
      {
	min = t;
	*quelle = gtliste[i+LSTART].nummer - MINKNR;
	*senke = gtliste[j+LSTART].nummer - MINKNR;
      }

  my_free(umsetz);
  my_free(pneu);
  return(OK);
}


static int sch_null(ULONG *pneu, ULONG knzahl, ULONG *umsetz, SCHICHT *schichten)
/* Knoten aus verschiedenen Schichten loeschen */
{
  KNOTEN *kn1, *kn2;
  SCHICHT *sl, *sl2;

  /* Schichten durchlaufen */
#ifdef DEB_ABSTSCH
printf("Knoten aus verschiedenen Schichten loeschen:\n");
#endif
  sl = schichten;
  while(sl != NULL)
  {
#ifdef DEB_ABSTSCH
printf("\nSchicht %lu \n", sl->nummer);
#endif
    /* Knoten der Schicht durchlaufen */
    kn1 = sl->knoten;
    while(kn1 != NULL)
    {
#ifdef DEB_ABSTSCH
printf("\nKnoten %lu \n", kn1->nummer);
#endif
      /* Alle anderen Schichten durchlaufen */
      sl2 = schichten;
      while(sl2 != NULL)
      {
	if(sl2 != sl)
	{
#ifdef DEB_ABSTSCH
printf("\ns%lu: ", sl2->nummer);
#endif
	  /* Die Knoten dieser anderen Schichten durchlaufen */
	  kn2 = sl2->knoten;
	  while(kn2 != NULL)
	  {
#ifdef DEB_ABSTSCH
printf("%lu ", kn2->nummer);
if(umsetz[kn1->nummer] * knzahl + umsetz[kn2->nummer] >= knzahl*knzahl)
printf("\nMatrixzugriffsfehler!\n");
#endif
	    /* Loeschen des Abstandes zwischen kn1 und kn2 */
            pneu[umsetz[kn1->nummer] * knzahl + umsetz[kn2->nummer]] = 0L;
	    kn2 = kn2->nachfolger;
	  }
	}
        sl2 = sl2->naechste;
      }
      kn1 = kn1->nachfolger;
    }
    sl = sl->naechste;
  }
  return(OK);
}


int qs_abstsch(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten)
/* Groesster Abstand mittels Abstandsmatrix in einer Schicht */
{
  ULONG knzahl = 0L, i, j, min, t;
  ULONG *pneu = NULL, *umsetz = NULL;
  KNOTEN *kn1, *kn2;
  SCHICHT *sl, *sl2;
  ULONG maxnr = 0L;

  /* Initialisierung von Quelle und Senke */
  if(gtliste[LSTART].nummer <= 0L || gtliste[LSTART + 1].nummer <= 0L)
    return(VP_PARAM);
  *quelle = (ULONG) gtliste[LSTART].nummer - MINKNR;
  *senke = (ULONG) gtliste[LSTART+1].nummer - MINKNR;

#ifdef DEB_ABSTSCH
printf("qs_abstsch():\n");
print_kn_aus_sc(schichten);
/*print_gtliste(gtliste);*/
printf("Knotenzahl bestimmen:\n");
fflush(stdout);
#endif

  /* Maximale Knotennummer und Knotenanzahl bestimmen */
  for(i = LSTART; gtliste[i].nummer > 0; ++i, ++knzahl)
    if(maxnr < gtliste[i].nummer)
      maxnr = gtliste[i].nummer;

#ifdef DEB_ABSTSCH
printf("knzahl=%lu\n", knzahl);
fflush(stdout);
#endif

  /* Abstandsmatrix berechnen */
  pneu = abstandsmatrix(&umsetz, gtliste, knzahl, maxnr);
  if(pneu == NULL)
    return(NOMEM);
  if(OK != sch_null(pneu, knzahl, umsetz, schichten))
  {
    my_free(pneu);
    my_free(umsetz);
    return(NOMEM);
  }

  /* Quelle und Senke aussuchen */
  min = 0L;
  for(i = 0L; i < knzahl; ++i)
    for(j = 0L; j < knzahl; ++j)
      if((t = pneu[i*knzahl+j]) > min)
      {
	min = t;
	*quelle = gtliste[i+LSTART].nummer - MINKNR;
	*senke = gtliste[j+LSTART].nummer - MINKNR;
      }

#ifdef DEB_ABSTSCH
printf("\nQuelle und Sebke: %lu %lu\n", *quelle, *senke);
#endif
  my_free(pneu);
  my_free(umsetz);
  return(OK);
}


int qs_manabs(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten)
/* Maximale Nachbarn mit maximalem Abstand */
{
  ULONG i, knzahl = 0L, maxnr = 0L, *pneu = NULL, *umsetz = NULL;
  KNOTEN *qu = NULL, *se = NULL;
  ULONG quz = 0L, sez = 0L, naz;
  ULONG *erste, *zweite, *copy;
  ULONG erste_zahl = 0L, zweite_zahl = 0L, erste_nachz = 0L, zweite_nachz = 0L;
  int rueck = OK;

  if (gtliste[LSTART].nummer <= 0L || gtliste[LSTART + 1].nummer <= 0L)
    return(VP_PARAM);

  /* Maximale Knotennummer und Knotenanzahl bestimmen */
  for (i = LSTART; gtliste[i].nummer > 0; ++i, ++knzahl)
    if(maxnr < gtliste[i].nummer)
      maxnr = gtliste[i].nummer;

  /* Abstandsmatrix berechnen */
  pneu = abstandsmatrix(&umsetz, gtliste, knzahl, maxnr);
  if (pneu == NULL)
    return(NOMEM);
  /*sch_null(pneu, knzahl, umsetz, schichten);*/   /* Bringt nichts! */

  /* Felder fuer maximale Knoten holen */
  if ((erste = (ULONG *) my_calloc(knzahl * sizeof(ULONG), "vbp")) == NULL)
  {
    my_free(pneu);
    my_free(umsetz);
    return(NOMEM);
  }
  if((zweite = (ULONG *) my_calloc(knzahl * sizeof(ULONG), "vbp")) == NULL)
  {
    my_free(pneu);
    my_free(umsetz);
    my_free(erste);
    return(NOMEM);
  }

  /* Initialisieren */
  erste[erste_zahl++] = gtliste[LSTART].nummer - MINKNR;
  naz = erste_nachz = nachz(&gtliste[LSTART]);

#ifdef DEB_MANABS
printf("\n %lu mit %lu \n", gtliste[LSTART].nummer - MINKNR, naz);
#endif
  /* Suche Knoten mit mehr Nachbarn */
  i = LSTART + 1;
  while(gtliste[i].nummer > 0L)
  {
    naz = nachz(&gtliste[i]);
#ifdef DEB_MANABS
printf(" %lu mit %lu \n", gtliste[i].nummer - MINKNR, naz);
#endif
    if(erste_nachz == naz)   /* Nachfolgerzahl wie bisher ersten */
    {
      erste[erste_zahl++] = gtliste[i].nummer - MINKNR;
    }
    else if(erste_nachz < naz)   /* Mehr als bisher ersten */
    {
      copy = erste; erste = zweite; zweite = copy;
      zweite_zahl = erste_zahl; zweite_nachz = erste_nachz;
      erste_zahl = 0L; 
      erste_nachz = naz;
      erste[erste_zahl++] = gtliste[i].nummer - MINKNR;
    }
    else if(zweite_nachz == naz)   /* wie bisher zweiten */
    {
      zweite[zweite_zahl++] = gtliste[i].nummer - MINKNR;
      zweite_nachz = naz;
    }
    else if(zweite_nachz < naz)   /* Zwischen ersten und zweiten */
    {
      zweite_nachz = naz;
      zweite_zahl = 0L;
      zweite[zweite_zahl++] = gtliste[i].nummer - MINKNR;
    }
    i++;
  }

#ifdef DEB_MANABS
  printf("\nErste:");
  for(i=0;i<erste_zahl;++i)
    printf(" erste[%d]=%lu ", i, erste[i]);
  printf("\nZweite:");
  for(i=0;i<zweite_zahl;++i)
    printf(" zweite[%d]=%lu ", i, zweite[i]);
  printf("\n");
#endif

  if(erste_zahl == 2L)
  {
    /* Es gibt genau zwei Knoten mit maximaler Nachfolgerzahl */
    *quelle = erste[0L];
    *senke = erste[1L];
  }
  else if(erste_zahl > 2L)
  {
    /* Es gibt mehr als zwei Knoten mit maximaler Nachfolgerzahl */
    ULONG j, quellind, neu_abst, abst;
    *quelle = erste[0L];
    *senke = erste[1L];
    abst = pneu[umsetz[*quelle]*knzahl + umsetz[*senke]];
    for(i=0L; i<erste_zahl-1L; i++)
    {
      quellind = umsetz[erste[i]]*knzahl;
      for(j=0L; j<erste_zahl; j++)
      {
        neu_abst = pneu[quellind + umsetz[erste[j]]];
        if (neu_abst > abst)  
        {
          *quelle = erste[i];
          *senke = erste[j];
          abst = neu_abst;
        }
      }
    }
  }
  else if(erste_zahl == 1L)
  {
    /* Es gibt nur einen Knoten mit maximaler Nachfolgerzahl */
    *quelle = erste[0L];
    if(zweite_zahl == 1L)
    {
      /* Es gibt nur einen Knoten mit zweitmaximaler Nachfolgerzahl */
      *senke = zweite[0L];
    }
    else if(zweite_zahl > 1L)
    {
      /* Es gibt mehrere Knoten mit zweitmaximaler Nachfolgerzahl */
      ULONG neu, quellind = umsetz[*quelle] * knzahl;
      ULONG sn_nr = zweite[0L];
      /*ULONG sn_ab = pneu[umsetz[*quelle] * knzahl + umsetz[sn_nr]];*/
      ULONG sn_ab = pneu[quellind + umsetz[sn_nr]];

#ifdef DEB_MANABS
printf(" (%lu ,%lu) ", sn_nr, sn_ab);
#endif

      for(i = 1L; i < zweite_zahl; ++i)   
      {
	/* Suche den Knoten mit maximalem Abstand heraus */
	neu = pneu[quellind + umsetz[zweite[i]]];
#ifdef DEB_MANABS
printf(" (%lu ,%lu) ", zweite[i], neu);
#endif
	if(neu > sn_ab)
	{
	  sn_nr = zweite[i];
	  sn_ab = neu;
	}
      }
      *senke = sn_nr;
    }
    else
      rueck = VP_QSERR;
  }
  else
    rueck = VP_QSERR;


#ifdef DEB_MANABS
printf("Quelle %lu und Senke %lu!\n", *quelle, *senke);
#endif

  my_free(umsetz);
  my_free(pneu);
  my_free(erste); 
  my_free(zweite);
  return(rueck);
}


int qs_manach(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten)
/* Maximale Nachbarn */
{
  ULONG i;
  KNOTEN *qu = NULL, *se = NULL;
  ULONG quz = 0L, sez = 0L, naz;

  if(gtliste[LSTART].nummer <= 0L || gtliste[LSTART + 1].nummer <= 0L)
    return(VP_PARAM);

  /* Initialisieren */
  qu = &gtliste[LSTART];
  se = &gtliste[LSTART + 1];
  quz = nachz(qu);
  sez = nachz(se);

  /* Suche Knoten mit mehr Nachbarn als qu oder se */
  i = LSTART + 2;
  while(gtliste[i].nummer > 0L)
  {
    naz = nachz(&gtliste[i]);
    if(quz < naz)
    {
      qu = &gtliste[i];
      quz = naz;
    }
    else if(sez < naz)
    {
      se = &gtliste[i];
      sez = naz;
    }
    i++;
  }

  /* Quelle und Senke belegen */
  *quelle = qu->nummer - MINKNR;
  *senke = se->nummer - MINKNR;

  return(OK);
}


static ULONG nachz(KNOTEN *kn)
/* Berechnung der Anzahl der Nachfolger */
{
  ULONG naz = 0L;
  
  while(kn->nachfolger != NULL)
  {
    ++naz;
    kn = kn->nachfolger;
  }
  return(naz);
}


int qs_klsch(ULONG *quelle, ULONG *senke, KNOTEN *gtliste, SCHICHT *schichten)
/* Quelle und Senke bestimmen: Kleinste Schicht */
{
  ULONG min_zahl = VBP_ULONG_MAX;
  SCHICHT *sl, *min_schicht = NULL;
  KNOTEN *kn;

  /* Schicht mit minimaler Knotenanzahl errechnen */
  sl = schichten;
  while(sl != NULL)
  {
    if(min_zahl > sl->knzahl && sl->knzahl > 1)
    {
      min_zahl = sl->knzahl;
      min_schicht = sl;
    }
    sl = sl->naechste;
  }
  if(min_schicht == NULL) 
    return(VP_PARAM);

  /* Quelle und Senke belegen */
  kn = min_schicht->knoten;
  *quelle = kn->nummer;
  kn = kn->nachfolger;
  *senke = kn->nummer;

  return(OK);
}


/* Ausgabe- und Fehlerfunktionen */
#ifdef KOMMENTAR    
#ifndef DEB_AUSGABE
#define DEB_AUSGABE   1
#endif
#endif
#ifdef XKDEB   
#ifndef DEB_AUSGABE
#define DEB_AUSGABE   1
#endif
#endif
#ifdef DEB_ABSTSCH  
#ifndef DEB_AUSGABE
#define DEB_AUSGABE   1
#endif
#endif
#ifdef DEBUG_ST     
#ifndef DEB_AUSGABE
#define DEB_AUSGABE   1
#endif
#endif
#ifdef FORTGANG     
#ifndef DEB_AUSGABE
#define DEB_AUSGABE   1
#endif
#endif

#ifdef DEB_AUSGABE

static void print_schichten(SCHICHT *schichten, FLOAT *plazierung)
{
  SCHICHT *arbeit = schichten;
  KNOTEN *knoten;

  while(arbeit != NULL)
  {
    printf("S %ld(%ld): ", arbeit->nr, arbeit->knzahl);
    knoten = arbeit->knoten;
    while(knoten != NULL)
    {
      printf("%ld (%.3f,%.3f)  ", knoten->nummer, plazierung[knoten->nummer * 2], plazierung[knoten->nummer * 2 + 1]);
      knoten = knoten->nachfolger;
    }
    printf("\n");
    arbeit = arbeit->naechste;
  }
}


static void print_kn_aus_sc(SCHICHT *schichten)
{
  SCHICHT *arbeit = schichten;
  KNOTEN *knoten;

  printf("Schichtenmodell:\n");
  if(schichten == NULL) printf("Leere Schichten!");
  else
  while(arbeit != NULL)
  {
    printf("S %ld (%ld): ", arbeit->nr, arbeit->knzahl);
    knoten = arbeit->knoten;
    while(knoten != NULL)
    {
      printf("%lu  ", knoten->nummer);
      knoten = knoten->nachfolger;
    }
    arbeit = arbeit->naechste;
    printf("\n");
  }
  /*printf("\n");*/
}

static void print_liste(ULONG *ein)
{
  ULONG i, knoten, nachbar, **nl, zahl;

  /* Gerichtete Nachbarschaftsliste ein durchlaufen */
  for(knoten = 0; knoten < ein[0] - 1; ++knoten)
  {
    printf("\n%lu: ", knoten);
    for(nachbar = ein[knoten]; nachbar < ein[knoten + 1]; ++nachbar)
    {
      printf("%lu ", ein[nachbar]);
    }
  }
  printf("\n");
}


static void print_gtliste(KNOTEN *ein)
{
  ULONG i, knoten, nachbar, **nl, zahl;
  KNOTEN *k;

  /* gtliste ein durchlaufen */
  for(i = LSTART; ein[i].nummer > 0; ++i)
  {
    printf("\ngtliste[%lu] = %ld: ", i, ein[i].nummer);
    k = ein[i].nachfolger;
    while (k) {
    	printf("%ld/%ld, ", k->nummer, k->gewicht);
    	k = k->nachfolger;
    	}
  }
  printf("\ngtliste[%lu] = %ld", i, ein[i].nummer);
  printf("\n");
}

#endif

/* Beispielverband */
/* ULONG nl[] = {7, 9, 9, 11, 14, 15, 17, 2, 4, 5, 1, 1, 2, 4, 1, 1, 4}; */
/* ULONG nl[] = {5, 7, 8, 9, 9, 1, 2, 3, 3}; */
ULONG nl[] = {9, 11, 14, 15, 16, 17, 18, 19, 19, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7};

#if FALSE
int main()
{
  ULONG *erg = NULL;
  FLOAT *plazierung = NULL;
  LONG i;
  
  if(verband_plazieren_qs(nl, &plazierung, qs_abst) == OK)
  {
    /*x_streckung(n, plazierung);*/
    for (i = 0; i < nl[0] - 1; i++) {
    	printf("%ld: %f %f\n", i, plazierung[2 * i], plazierung[2 * i + 1]);
    	}
    /* t_gerichtet_ungerichtet(nl[i], &erg); */
    my_free(plazierung);
    my_free(erg);
  }
  return(0);
}
#endif

/* gr_teil.c: */
/*
#define GTDEBUG     1
#define FORTGANG    1
*/


#define KORR_VOR    0                     /* Vorwaerts: Hinzunahme von Knoten */
#define KORR_ZUR    1                         /* Zurueck: Wegnahme von Knoten */


int gf[MAXKNOTEN];
struct pfeile gpf[maxpfeile];
struct inhalt ginh[maxpfeile];
 
int GRAPH_TEILEN(KNOTEN liste[], KNOTEN liste1[], KNOTEN liste2[], int quelle, 
	int senke, int korrektur, int kommentar)
/* liste1: fuer markierbare Knoten 
 * liste2: fuer nicht markierbare Knoten */
{
  int kant, knot, ebenen = 1, index;
	int f[MAXKNOTEN];
	
  if(kommentar == 1)                            /* Ausgabe vor der Aufteilung */
  {
    printf("Knoten-nachfolger-liste : vor der Aufteilung\n");
    listenausgabe(liste);
  }

#ifdef FORTGANG
  printf("p");
  fflush(stdout);
#endif

  pfeilliste_erst(liste, gpf, ginh, &kant, &knot);     /* erstellt Pfeilliste */
  if(quelle == 0)                           /* Quelle und Senke nicht bekannt */
  {
    sort_nachfolger(1, knot, liste);
    quelle_senke_bestimmen(liste, &quelle, &senke, knot);
  }
  if(kant != 0)                                      /* Kanten sind vorhanden */
  {                                              /* Maximalen Fluss bestimmen */
#ifdef FORTGANG
  printf("m");
  fflush(stdout);
#endif
    if(!max_fluss(gpf, ginh, kant, knot, quelle, senke, f, kommentar)) 
    	/* f, kommentar hinzugefuegt !!! */
    {
      printf("Fehler: Der Maximale Fluss konnte nicht bestimmt werden!\n");
      return(GT_MAXFLUSS);
    }
#ifdef FORTGANG
  printf("k");
  fflush(stdout);
#endif
    gf[1] = quelle;                                /* Quelle immer markierbar */
    index = markiere_knoten(gf, gpf, ginh, kant, 1, quelle);
    gf[index+1] = 0;
    if(korrektur) 
      if(gf[2] == 0)                      /* Schnitt direkt hinter der Quelle */
        index = sch_korrektur(gpf, gf, knot, kant, quelle, senke, ebenen, KORR_VOR, index);
      else if(gf[knot - 1] != 0)              /* Schnitt direkt vor der Senke */
        index = sch_korrektur(gpf, gf, knot, kant, senke, quelle, ebenen, KORR_ZUR, index);
    gf[index+1] = 0;
#ifdef GTDEBUG
{
int i;
printf("Knoten in gf[]:");
for(i = 1; gf[i] > 0; ++i)
printf(" %d ", gf[i]);
printf("\n");
}
#endif
  }
  else                                              /* Keine Kanten vorhanden */
  {
    gf[1] = quelle;
    gf[2] = 0;
  }
#ifdef FORTGANG
  printf("l");
  fflush(stdout);
#endif
  liste_trennen(liste, liste1, liste2, gf, knot);    /* Trennt die Pfeilliste */
  if(kommentar == 1)
  {
    printf("Gefundener minimaler Schnitt :\n");
    printf("*******************************\n");
    printf("liste1 : markierbare knoten\n");
    listenausgabe(liste1);
    printf("liste2 : nicht markierbare knoten\n");
    listenausgabe(liste2);
  }
#ifdef FORTGANG
  printf("r");
  fflush(stdout);
#endif
  return(GT_OK);
}                                                         /* ende graphteilen */


                    /* markiere_knoten - Markierbare Knoten in gf abspeichern */
int markiere_knoten(int gf[], struct pfeile gpf[], struct inhalt ginh[], 
	int kanten, int index, int nr)
            /* (Knoten ist markierbar <=> die Kante zu diesem Knoten ist noch */
              /* nicht erschoepft) und bestimmt die moegliche Flussteigerung. */
/* int kanten: Kantenzahl 
 * int index: Letzter besetzter Platz in gf 
 * int nr: Nummer des aktuellen Knotens */
{
  int l, i = 1; 
 
  while(gpf[i].anfang != nr)        /* Sucht erste Kante mit Anfangsknoten nr */
    ++i;
  while(i <= kanten && gpf[i].anfang == nr)      /* Die von nr ausgehendenden */
  {                                               /* Kanten werden betrachtet */
    l = gpf[i].nummer;
    if(ginh[l].kapazitaet > ginh[l].fluss)         /* Kante nicht ausgelastet */
      if(!in_vektor(gpf[i].ende, gf, 1, index))       /* Endknoten markierbar */
      {                /* wenn dieser noch nicht in gf ist, wird er angefuegt */
#ifdef GTDEBUG
printf("Neu: %d da %d > %d \n", gpf[i].ende, ginh[l].kapazitaet, ginh[l].fluss);
#endif
        gf[++index] = gpf[i].ende;
        index = markiere_knoten(gf, gpf, ginh, kanten, index, gpf[i].ende);
      }
    ++i;
  }
  gf[index+1] = 0;

#ifdef GTDEBUG
{
int i;
printf("markiere_knoten: in gf[]:");
for(i = 1; gf[i] > 0; ++i)
printf(" %d ", gf[i]);
printf("\n");
}
#endif
  return(index);
}                                                 /* Ende von markiere_knoten */


            /* sch_korrektur - Schnittkorrektur um Ebenengroesse durchfuehren */
	    /* Aenderung: sch_korrektur fuegt weder Quelle noch Senke an (JJ) */
int sch_korrektur(struct pfeile gpf[], int gf[], int knoten, int kanten, 
	int start, int qs_nicht, int eb, int rich, int index)
#if FALSE
int knoten;                                   /* Anzahl der Knoten im Graphen */
int kanten;                                   /* Anzahl der Kanten im Graphen */
int start;                              /* Startknoten zur Schnitterweiterung */
int qs_nicht;                       /* Quelle bzw Senke, wird nicht angefuegt */
int eb;                                  /* Anzahl der Ebenen zur Erweiterung */
int rich;                                         /* Richtung der Erweiterung */
int index;
#endif
{
  int j, e, i = 1;

  if(knoten > 3 && eb > 0)                              /* Korrektur sinnvoll */
  {
    if(rich == KORR_VOR)                             /* Hinzunahme von Knoten */
    {
      while(gpf[i].anfang != start)   /* Die erste vom Startknoten ausgehende */
        ++i;                                            /* Kante wird gesucht */
      while(i <= kanten && gpf[i].anfang == start)     /* Alle Kanten mit An= */
      {                                                  /* fangsknoten start */
        e = gpf[i].ende;
	if(e != qs_nicht)                                               /* JJ */
        if(!in_vektor(e, gf, 1, index))   /* ist der Endknoten der Kante noch */
          gf[++index] = e;                  /* nicht in gf, wird er angefuegt */
        index = sch_korrektur(gpf, gf, knoten, kanten, e, qs_nicht, eb-1, rich, index);
                                                   /* naechste Korrekturebene */
        ++i;
      }
    }
    else                                               /* Wegnahme von Knoten */
    {
      while(i <= kanten)
      {
        if(gpf[i].ende == start)            /* Kante zum Startknoten gefunden */
        {
          e = gpf[i].anfang;
	  if(e != qs_nicht)                                               /* JJ */
	  {
          for(j = 1; j <= index; j++)    /* Sucht den Anfangsknoten der Kante */
            if(gf[j] == e)                                           /* in gf */
              break;
          if(j <= index)                            /* Wenn Suche erfolgreich */
          {
            if(j < index)
              gf[j] = gf[index]; /* Knoten aus gf loeschen, index verkleinern */
            gf[index--] = 0; 
          }
	  }
          index = sch_korrektur(gpf, gf, knoten, kanten, e, qs_nicht, eb-1, rich, index);
                                                  /* Naechste Korrekturebene */
        }
        ++i;
      }
    }
  }
  return(index);
  }                                                 /* Ende von sch_korrektur */


void pfeilliste_erst(KNOTEN l[MAXKNOTEN],struct pfeile pfl[maxpfeile],
	struct inhalt inhl[maxpfeile], int *kantenzahl, int *knotenzahl)
/**  Diese Funktion erstellt aus der Nachfolgerliste 'l' eine    */
/**  Pfeilliste 'pfl' (Infos ueber die Kanten: Anfangs-, Endknoten, */
/**  Kantennr. und Richtung) und eine Inhaltsliste 'inhl' (Kanten-  */
/**  gewicht = Kapazitaet,Fluss), die zur Berechnung des minimalen  */
/**  Schnitts benoetigt werden.                                     */
{
	KNOTEN *lauf_z, *knoten_z;
	int i, j = 0;
	int schon_da[MAXKNOTEN];

	for (i = 1; l[i].nummer > 0; i++) {
		lauf_z = &l[i];
		knoten_z = &l[i];
		schon_da[i] = l[i].nummer;
		schon_da[i + 1] = 0; 
		while (lauf_z->nachfolger != NULL) {
			lauf_z = lauf_z->nachfolger;
			if (!in_vektor(lauf_z->nummer, schon_da, 1, i)) {
				pfl[++j].anfang = knoten_z->nummer;
				pfl[j].ende = lauf_z->nummer;
				pfl[j].nummer = j;
				pfl[j].zuruecknummer = j + 1;
				pfl[j].richtung = 1;
				inhl[j].kapazitaet = lauf_z->gewicht;
				inhl[j].fluss = 0;
				pfl[++j].anfang = lauf_z->nummer;
				pfl[j].ende = knoten_z->nummer;
				pfl[j].nummer = j;
				pfl[j].zuruecknummer = j - 1;
				pfl[j].richtung = 1;
				inhl[j].kapazitaet = lauf_z->gewicht;
				inhl[j].fluss = 0;
				}
			}
		}
	*kantenzahl = j;
	*knotenzahl = i - 1;
}

static INT kn_cmp(void *p1, void *p2, LONG *res);
static INT kn_swap(void *p1, void *p2);

static INT kn_cmp(void *p1, void *p2, LONG *res)
{
	LONG n1, n2;
	KNOTEN *pk1 = (KNOTEN *) p1;
	KNOTEN *pk2 = (KNOTEN *) p2;
	
	n1 = anz_nachfolger(pk1);
	n2 = anz_nachfolger(pk2);
	if (n1 < n2) {
		*res = -1L;
		return(TRUE);
		}
	if (n1 > n2) {
		*res = 1L;
		return(TRUE);
		}
	*res = 0L;
	return(TRUE);
}

static INT kn_swap(void *p1, void *p2)
{
	KNOTEN tmp;
	KNOTEN *pk1 = (KNOTEN *) p1;
	KNOTEN *pk2 = (KNOTEN *) p2;
	
	tmp = *pk1;
	*pk1 = *pk2;
	*pk2 = tmp;
	return(TRUE);
}

INT sort_nachfolger(int untergr, int obergr, KNOTEN list[MAXKNOTEN])
/* Diese Funktion sortiert die Nachfolgerliste 'list' nach der      */
/* Anzahl der Nachfolger eines jeden Knotens.(rekursives Quicksort) */
{   
	Qsort_cmp_func = kn_cmp;
	Qsort_swap_func = kn_swap;
	Qsort_f_ascending = TRUE;
	Qsort_ElementSizeof = sizeof(KNOTEN);
	if (!Q2sort((BYTE *)list, untergr, obergr)) {
		printf("sort_nachfolger()|error in Q2sort\n");
		return(FALSE);
		}
	return(TRUE);
#if FALSE
	int x,i,m,i1,i2,i3;
    KNOTEN test[MAXKNOTEN];
    KNOTEN l1[MAXKNOTEN],l2[MAXKNOTEN],l3[MAXKNOTEN];
    
    m = (untergr + obergr) / 2;
    i1 = 0;
    i2 = 0;
    i3 = 0;
    test[1] = list[m];
    for (i = untergr; i <= obergr; ++i)
    {   x = anz_nachfolger(&list[i]);
        if (x > anz_nachfolger(&test[1]))
        {   ++i1;
            l1[i1] = list[i];
        }
        else if ( x < anz_nachfolger(&test[1]) )
        {   ++i3;
            l3[i3] = list[i];
        }
        else if (x == anz_nachfolger(&test[1]))
        {   ++i2;
            l2[i2] = list[i];
        }
    }
    for (i = 1; i <= i1; ++i)
        list[untergr - 1 + i] = l1[i];
    for (i = 1; i <= i2; ++i)
        list[untergr - 1 + i1 + i] = l2[i];
    for (i = 1; i <= i3; ++i)
        list[untergr - 1 + i1 + i2 + i] = l3[i];
    if (i1 > 1)
        sort_nachfolger(untergr,untergr + i1 - 1,list);
    if (i3 > 1)
        sort_nachfolger(untergr + i1 + i2, obergr,list);
#endif
}

int anz_nachfolger(KNOTEN *lauf_z)
/**  Diese Funktion bestimmt die Anzahl der Nachfolger des Knotens, */
/**  auf den der Zeiger 'lauf_z' zeigt.                             */
{
	int i;
	for (i=0; lauf_z->nachfolger != NULL; i++, lauf_z = lauf_z->nachfolger) ;
	return(i);
}

void quelle_senke_bestimmen(KNOTEN l[MAXKNOTEN],int *q, int *s, int knotenzahl)
/* Diese Funktion bestimmt als Quelle und Senke zwei moeglichst     */
/* nicht benachbarte Knoten mit moeglichst vielen Nachfolgern.      */
{
int i,j;
KNOTEN *lauf_z;
int benachbart=1,whileabbruch;

for(i=1;i<=knotenzahl&&benachbart;++i)
   {
   *q=l[i].nummer;
   for(j=i+1;j<=knotenzahl&&benachbart;j++)
      {
      lauf_z= &l[j];
      *s=l[j].nummer;
      whileabbruch=0;
      while(lauf_z->nachfolger != NULL &&!whileabbruch)
           {
           lauf_z=lauf_z->nachfolger;
           if(lauf_z->nummer == *q)
             {
             benachbart=1;
             whileabbruch=1;
             }
           else
             {
             benachbart=0;
             }
           }
      }
   }
if(benachbart)
   {
   *q=l[1].nummer;
   *s=l[2].nummer;
   }
}

void liste_trennen(KNOTEN li[MAXKNOTEN],KNOTEN list1[MAXKNOTEN],
	KNOTEN list2[MAXKNOTEN],int f[MAXKNOTEN],int knotenzahl)
/* Diese Funktion teilt die Nachfolgerliste 'li' in zwei Teillisten */
/* 'list1' und 'list2'. Dies bedeutet, dass der Graph bzgl. des      */
/* minimalen Schnitts in zwei Teilgraphen aufgeteilt wird.           */
{
int i,j,k,l,l1=0,l2=0;

k = 1;
for(j=1;j<=knotenzahl;j++)
   {
   i=li[j].nummer;
   if(in_vektor(i,f,1,knotenzahl))
     {
     l1++;
     for(l=1;l<=knotenzahl&&li[l].nummer >0;l++)  
      {
      if(i==li[l].nummer)
         {
         list1[l1]=li[l];
         schnittkanten_zu_liste2_loeschen(f,&list1[l1],knotenzahl);
         }
     }
     }
   else 
     {
     l2++;
     for(l=1;l<=knotenzahl&&li[l].nummer >0;l++)
      {
      if(i==li[l].nummer)
         {
         list2[l2]=li[l];
         schnittkanten_zu_liste1_loeschen(f,&list2[l2],knotenzahl);
         }
      }
     }
   k++;
   }

#ifdef GTDEBUG
printf("liste_trennen(): %d - %d Knoten\n", l1, l2);
#endif

list1[l1+1].nummer= -1;
list2[l2+1].nummer= -1;
}

void schnittkanten_zu_liste1_loeschen(int f[MAXKNOTEN],KNOTEN *zeiger,
	int knotenzahl)
/* Diese Funktion loescht die Nachfolger, die nicht mehr zu diesem  */
/* Teilgraphen gehoeren.                                            */
{
while(zeiger->nachfolger != NULL)
     {
     if( in_vektor(zeiger->nachfolger->nummer,f,1,knotenzahl) )
       {
       zeiger->nachfolger=zeiger->nachfolger->nachfolger;
       }
     else
       {
       zeiger=zeiger->nachfolger;
       }
     }
}

void schnittkanten_zu_liste2_loeschen(int f[MAXKNOTEN],KNOTEN *zeiger,
	int knotenzahl)
/* Diese Funktion loescht die Nachfolger, die nicht mehr zu diesem  */
/* Teilgraphen gehoeren.                                            */
{
while(zeiger->nachfolger != NULL)
     {
     if( !in_vektor(zeiger->nachfolger->nummer,f,1,knotenzahl) )
       {
       zeiger->nachfolger=zeiger->nachfolger->nachfolger;
       }
     else
       {
       zeiger=zeiger->nachfolger;
       }
     }
}

void listenausgabe(KNOTEN l[MAXKNOTEN])
/* Diese Funktion gibt die Nachfolgeliste 'l' aus.               */
{
KNOTEN *lauf_z;
int i;

for(i=1;l[i].nummer>0;i++)
   {
   lauf_z= &l[i];
   printf("knoten: %d nachfolger: \n",l[i].nummer);
   while(lauf_z->nachfolger != NULL)
     {
     lauf_z=lauf_z->nachfolger;
     printf("          ---------->%d    mit kapazitaet: %d\n",
     	lauf_z->nummer,lauf_z->gewicht);
     }
   }
}

int in_vektor(int zu_suchen,int feld[MAXKNOTEN], int von, int bis)
/**  Dieses ist eine Hilfsfunktion,die bestimmt ob das Element          **/
/**  'zu_suchen' in 'feld' enthalten ist.                               **/
{
int i,enthalten=0;
for( i=von ; i<=bis && !enthalten && feld[i] != 0; i++ )
   {
   if(feld[i]==zu_suchen)
     {
     enthalten=1;
     return(1);
     }
   }
return(0);
}
/* end gr_teil.c */

/* ff.c: */
static INT sort_liste(int untergr,int obergr, struct pfeile pf[maxpfeile]);
static void intervallvektor(int kantenzahl, struct pfeile pf[maxpfeile],
	int vektor[maxpfeile], int merke[MAXKNOTEN]);
static void flussteigern(int knotenzahl, int kantenzahl, int senke, int quelle,
	int f[maxpfeile], int vektor[maxpfeile], struct pfeile pf[maxpfeile], 
	struct fluss_knoten kn[maxpfeile], struct inhalt inh[maxpfeile], 
	int merke[MAXKNOTEN]);
static void markiere(int diff, int j, int k, int x, int *y, 
	int f[MAXKNOTEN], struct pfeile pf[maxpfeile], 
	struct fluss_knoten kn[maxpfeile], struct inhalt inh[maxpfeile]);

int max_fluss(struct pfeile pf[maxpfeile], struct inhalt inh[maxpfeile], 
	int kantenzahl, int knotenzahl, int quelle, int senke, int f[MAXKNOTEN], 
	int korrektur)
/**  Diese Funktion bestimmt nach dem Algorithmus von Ford & Fulkerson  **/
/**  in 'flussteigern' den maximalen Fluss durch den Graphen und be-    **/
/**  stimmt den minimalen Schnitt.                                      **/
{
	struct fluss_knoten kn[maxpfeile];
	int vektor[maxpfeile];
	int merke[MAXKNOTEN];

	sort_liste(1,kantenzahl,pf);
	intervallvektor(kantenzahl,pf,vektor,merke);
	flussteigern(knotenzahl,kantenzahl,senke,quelle,f,vektor,pf,kn,inh,merke);
	return(TRUE);
}

INT (*Qsort_cmp_func)(void *p1, void *p2, LONG *res) = NIL;
/* result > 0: p1 > p2
 *        = 0: p1 = p2
 *        < 0: p1 < p2 
 */
INT (*Qsort_swap_func)(void *p1, void *p2) = NIL;
INT Qsort_f_ascending = TRUE;
INT Qsort_ElementSizeof = 0L;

static INT pf_cmp(void *p1, void *p2, LONG *res);
static INT pf_swap(void *p1, void *p2);
static INT partition(BYTE *arr, LONG left, LONG right, LONG *middle);

static INT pf_cmp(void *p1, void *p2, LONG *res)
{
	struct pfeile *pf1 = (struct pfeile *) p1;
	struct pfeile *pf2 = (struct pfeile *) p2;
	
	if (pf1->anfang < pf2->anfang) {
		*res = -1L;
		return(TRUE);
		}
	if (pf1->anfang > pf2->anfang) {
		*res = 1L;
		return(TRUE);
		}
	*res = 0L;
	return(TRUE);
}

static INT pf_swap(void *p1, void *p2)
{
	struct pfeile tmp;
	struct pfeile *pf1 = (struct pfeile *) p1;
	struct pfeile *pf2 = (struct pfeile *) p2;
	
	tmp = *pf1;
	*pf1 = *pf2;
	*pf2 = tmp;
	return(TRUE);
}

static INT sort_liste(int untergr,int obergr, struct pfeile pf[maxpfeile])
/**  Diese Funktion sortiert die Pfeilliste 'pf' nach der Nummer der    **/
/**  Anfangsknoten.(rekursives Quicksort)                               **/
{
	Qsort_cmp_func = pf_cmp;
	Qsort_swap_func = pf_swap;
	Qsort_f_ascending = TRUE;
	Qsort_ElementSizeof = sizeof(struct pfeile);
	if (!Q2sort((BYTE *)pf, untergr, obergr)) {
		printf("sort_liste()|error in Q2sort\n");
		return(FALSE);
		}
	return(TRUE);
#if FALSE
   int x,i,m,i1,i2,i3;
    struct pfeile test[maxpfeile];
    struct pfeile l1[maxpfeile],l2[maxpfeile],l3[maxpfeile];
    
    m = (untergr + obergr) / 2;
    i1 = 0;
    i2 = 0;
    i3 = 0;
    test[1] = pf[m];
    for (i = untergr; i <= obergr; ++i)
    {   x = pf[i].anfang;
     if (x < test[1].anfang)
        {   ++i1;
            l1[i1] = pf[i];
        }
        else if (x > test[1].anfang)
        {   ++i3;
            l3[i3] = pf[i];
        }
        else if (x == test[1].anfang)
        {   ++i2;
            l2[i2] = pf[i];
        }
    }
    for (i = 1; i <= i1; ++i)
        pf[untergr - 1 + i] = l1[i];
    for (i = 1; i <= i2; ++i)
        pf[untergr - 1 + i1 + i] = l2[i];
    for (i = 1; i <= i3; ++i)
        pf[untergr - 1 + i1 + i2 + i] = l3[i];
    if (i1 > 1)
        sort_liste(untergr,untergr + i1 - 1,pf);
    if (i3 > 1)
        sort_liste(untergr + i1 + i2, obergr,pf);
#endif
}

INT Q2sort(BYTE *arr, LONG left, LONG right)
{
	LONG middle;
	
	if (left < right) {
		if (!partition(arr, left, right, &middle)) {
			printf("Q2sort()|partition()\n");
			return(FALSE);
			}
		if (!Q2sort(arr, left, middle - 1L) || 
			!Q2sort(arr, middle + 1L, right) ) {
			printf("Q2sort()|Q2sort()\n");
			return(FALSE);
			}
		}
	return(TRUE);
}

static INT partition(BYTE *arr, LONG left, LONG right, LONG *middle)
{
	void *pivot;
	LONG l, r, m;
	LONG res;
	
	/* Pivot Strategie: nimm' linkes Element: */
	pivot = arr + left * Qsort_ElementSizeof;
	l = left;
	r = right;
	while (l < r) {
		while (TRUE) {
			if (l > right)
				break;
			if (!(*Qsort_cmp_func)(arr + l * Qsort_ElementSizeof, pivot, &res)) {
				printf("partition()|Qsort_cmp_func()\n");
				return(FALSE);
				}
			if (!Qsort_f_ascending)
				res *= -1L;
			if (res > 0L)
				break;
			l++;
			}
		while (TRUE) {
			if (r < left)
				break;
			if (!(*Qsort_cmp_func)(arr + r * Qsort_ElementSizeof, pivot, &res)) {
				printf("partition()|Qsort_cmp_func()\n");
				return(FALSE);
				}
			if (!Qsort_f_ascending)
				res *= -1L;
			if (res <= 0L)
				break;
			r--;
			}
		if (l < r) {
			if (!(*Qsort_swap_func)(arr + l * Qsort_ElementSizeof, 
					arr + r * Qsort_ElementSizeof)) {
				printf("partition()|Qsort_swap_func()\n");
				return(FALSE);
				}
			}
		}
	m = r;
	if (!(*Qsort_swap_func)(arr + left * Qsort_ElementSizeof, 
			arr + m * Qsort_ElementSizeof)) {
		printf("partition()|Qsort_swap_func()\n");
		return(FALSE);
		}
	*middle = m;
	return(TRUE);
}

static void intervallvektor(int kantenzahl, struct pfeile pf[maxpfeile],
	int vektor[maxpfeile], int merke[MAXKNOTEN])
/**  Diese Funktion speichert in 'vektor' ab, in welchem Intervall      **/
/**  in der Pfeilliste 'pf' die Pfeile der jeweiligen Anfangsknoten     **/
/**  stehen.                                                            **/
{   int x,i,j,k;
 
    j = 1;
    x = pf[1].anfang;
    merke[x] = j ;
    for (i = 2; i <= kantenzahl; ++i)
    {   vektor[i] = 0;
        vektor[j] = i - 1;
        if (pf[i].anfang > x)
        {   x = pf[i].anfang;
            ++j;
            merke[x] = j ;
        }
    }
        vektor[j] = i - 1;
}

static void flussteigern(int knotenzahl, int kantenzahl, int senke, int quelle,
	int f[maxpfeile], int vektor[maxpfeile], struct pfeile pf[maxpfeile], 
	struct fluss_knoten kn[maxpfeile], struct inhalt inh[maxpfeile], 
	int merke[MAXKNOTEN])
/**  Diese Funktion versucht iterativ den Fluss durch den Graphen zu    **/
/**  erhoehen, bis der maximale Fluss erreicht ist.                     **/
{   int i,j,l,y,u,fl,abbruch1,abbruch2,x,k,diff;
    int hilfdiff = 0;
    diff = 0;
    abbruch1 = 0;
    fl = 0;
    while( !abbruch1 )
    {
        abbruch2 = 0;
        for( i = 1; i <= MAXKNOTEN; ++i )
        {   f[i] = 0;
            kn[i].vorgaenger = i;
            kn[i].pfeilnummer = 0;
            kn[i].fluss = 0;
            kn[i].ausgelastet = 1;
        }
        kn[quelle].fluss = 100000;
        kn[quelle].ausgelastet = 0;
        f[1] = quelle;
        y = 1;
        for( i = 1; i <= kantenzahl && !abbruch1 && !abbruch2; ++i)
        {    
            j = f[i];
            if( j == 0 )
            {  
               abbruch1 = 1;
               
            }
            if( !abbruch1 )
            {   if(merke[j] == 1 )
                   u = 1;
                else
                   u = vektor[merke[j] - 1] + 1;
                for( x = u; x <= vektor[merke[j]] && !abbruch2; ++x)
                {   k = pf[x].ende;
                    if( kn[k].ausgelastet )
                    { 
                        markiere(diff,j,k,x,&y,f,pf,kn,inh);                               
                        if( k == senke )
                           if( kn[k].ausgelastet == 0 )
                              abbruch2 = 1;
                    }
                }
            }
        }
        if( !abbruch1 && quelle != k )  
           {   
           /* Fuss um diff von quelle nach k erhoehen: 
            * Der Fluss zwischen zwei Knoten wird in inh[].fluss eingetragen. 
            * Ein evtl. Fluss in der Rueckrichtung wird verringert. 
            * Die Route ergibt sich aus kn[].vorgaenger. */
           diff = kn[k].fluss;
           fl = fl + diff;
           i = k;
           do
              {   x = kn[i].pfeilnummer;
                  j = pf[x].nummer;
                  l = pf[x].zuruecknummer;
                  if(inh[l].fluss>0)
                    {
                    hilfdiff=diff;
                    if(hilfdiff<=inh[l].fluss)
                      {
                      inh[l].fluss -= hilfdiff;
                      }
                    else
                      {
                      hilfdiff-=inh[l].fluss;
                      inh[l].fluss = 0;
                      inh[j].fluss += hilfdiff;
                      }
                    }
                  else
                    {
                    inh[j].fluss += diff;
                    }
                  i = kn[i].vorgaenger;
              } while( i != quelle );
           }
    else
       {
       abbruch1 = 1;
       }
    }
}

static void markiere(int diff, int j, int k, int x, int *y, 
	int f[MAXKNOTEN], struct pfeile pf[maxpfeile], 
	struct fluss_knoten kn[maxpfeile], struct inhalt inh[maxpfeile])
/**  Diese Funktion speichert in 'f' die markierbaren Knoten ab.        **/
/**  (Knoten ist markierbar <=> die Kante zu diesem Knoten ist noch     **/
/**  nicht erschoepft) und bestimmt die moegliche Flussteigerung.       **/
{   int i,l,wert;
    
    i = pf[x].nummer;
    l = pf[x].zuruecknummer;
    wert = kn[j].fluss;
    diff = inh[l].fluss+inh[i].kapazitaet -inh[i].fluss;
    if(diff<wert)
      {
      wert = diff;
      }
    if( wert > 0 )
        {
        kn[k].vorgaenger = j;
        kn[k].pfeilnummer = x;
        kn[k].fluss = wert;
        kn[k].ausgelastet = 0;
        (*y)++;
        f[*y] = k;
    }
}

/* end ff.c */

#endif /* GRAPHICS_TRUE */

