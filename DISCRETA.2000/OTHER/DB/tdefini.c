/* ----------------------------------------------------------------------- */
/* TDEFINI.C                                                               */
/*                                                                         */
/* Diese Datei enthaelt Funktionen zum Initialisieren der Indexstruktur und */
/* zur Behandlung der Indizes sowie die Vergleichs- und Lesefunktionen fuer */
/* die vordefinierten Standardtypen.                                       */
/*                                                                         */
/* Funktionen:                                                             */
/* get_definition                                                          */
/* read_txt                                                                */
/* cmpinfo                                                                 */
/* Standardtypenfunktionen                                                 */
/* ----------------------------------------------------------------------- */


#include <stdio.h>
#include <string.h>

#include "dbms.h"           /* Benutzerstrukturen/operationen              */
#include "dbmsdef.h"        /* Definitionen                                */
#include "dbmstype.h"       /* Interne Verwaltungsstrukturen               */
#include "dbmsprot.h"       /* Prototypen der internen Funktionen          */

/* ----------------------------------------------------------------------- */
/* Prototypen der Vergleichs- und Lesefunktionen der Standardtypen         */
/* ----------------------------------------------------------------------- */

int   ubyte_cmp(  void *, void *, unsigned );
void *ubyte_read( char *line, unsigned *elem, unsigned *size );
int   sbyte_cmp(  void *, void *, unsigned );
void *sbyte_read( char *line, unsigned *elem, unsigned *size );
int   uint_cmp(   void *, void *, unsigned );
void *uint_read( char *line, unsigned *elem, unsigned *size );
int   sint_cmp(   void *, void *, unsigned );
void *sint_read( char *line, unsigned *elem, unsigned *size );
int   ulong_cmp(  void *, void *, unsigned );
void *ulong_read( char *line, unsigned *elem, unsigned *size );
int   slong_cmp(  void *, void *, unsigned );
void *slong_read( char *line, unsigned *elem, unsigned *size );
int   float_cmp(  void *, void *, unsigned );
void *float_read( char *line, unsigned *elem, unsigned *size );
int   double_cmp( void *, void *, unsigned );
void *double_read( char *line, unsigned *elem, unsigned *size );
int   string_cmp( void *, void *, unsigned );
void *string_read( char *line, unsigned *elem, unsigned *size );

ULONG    page_size;
UINT     data_length;
UINT     index_count;
INDEX    index_db[INDEXC+1];
CHAR     read_access;
CHAR     change_access;
UINT     index_end;
char    *adminpw;

/* ----------------------------------------------------------------------- */
/* Name       get_dafinition                                               */
/* Definition int get_definition( char *name );                            */
/* Funktion   Liest die TDF-Datei der Tabelle (name) und setzt sie         */
/*            obigen Strukturen.                                           */
/* Ergebnis   0 OK, sonst Fehlermeldung                                    */
/* ----------------------------------------------------------------------- */
int get_definition( char *name ) {

   FILE *file;
   TDF tdf;
   int n;
   int error = DBMS_ALL_OK;

   file = fopen( name, "rb" );

   if ( file == NULL ) {
      error = DBMS_TDF_READ_ERROR;
   } else {
      if ( fread( &tdf, sizeof(TDF), 1, file ) != 1 ) {
         error = DBMS_TDF_READ_ERROR;
      }
      fclose( file );
   }

   if ( ! error ) {

      read_access   = tdf.read_access;
      change_access = tdf.change_access;

      page_size   = tdf.psize;
      data_length = tdf.reclen;
      index_count = tdf.icnt;
      index_end   = tdf.rpos;
      adminpw     = NULL;

      for ( n=0; n<tdf.icnt; n++ ) {
         index_db[n].special  = 0;
         index_db[n].number   = 0;
         index_db[n].position = tdf.pos[n];
         index_db[n].length   = tdf.elem[n];
         if ( tdf.type[n] > 9 && tdf.type[n] < 30 ) {
            index_db[n].length   = 1;
            index_db[n].cmpfunc  = slong_cmp;
            if ( tdf.type[n] < 20 ) {
               index_db[n].special = AUTOINCREMENT;
               index_db[n].number = tdf.type[n] - 10;
            } else {
               index_db[n].special = AUTODECREMENT;
               index_db[n].number = tdf.type[n] - 20;
            }
         } else {
            switch( tdf.type[n] ) {

               case 1: index_db[n].cmpfunc = ubyte_cmp; break;
               case 2: index_db[n].cmpfunc = sbyte_cmp; break;
               case 3: index_db[n].cmpfunc = uint_cmp; break;
               case 4: index_db[n].cmpfunc = sint_cmp; break;
               case 5: index_db[n].cmpfunc = ulong_cmp; break;
               case 6: index_db[n].cmpfunc = slong_cmp; break;
               case 7: index_db[n].cmpfunc = float_cmp; break;
               case 8: index_db[n].cmpfunc = double_cmp; break;
               case 9: index_db[n].cmpfunc = string_cmp; break;

               default:   error = DBMS_NO_DEFINITION ;
            }
         }
      }

      index_db[index_count].position = index_end;
      index_db[index_count].length   = tdf.rsize;
   }

   return( error );
}

/* ----------------------------------------------------------------------- */
/* Name        read_txt                                                    */
/* Definition  void *read_txt( unsigned type, char *line, unsigned *elem   */
/*                             unsigned *size );                           */
/* Funktion    Einlesen der Standardtypen im Textformat                    */
/* ----------------------------------------------------------------------- */
void *read_txt( unsigned type, char *line, unsigned *elem, unsigned *size ) {

   if ( type > 9 && type < 30 ) type = 6;

   switch ( type ) {

      case 1: return( ubyte_read( line, elem, size ) );
      case 2: return( sbyte_read( line, elem, size ) );
      case 3: return( uint_read( line, elem, size ) );
      case 4: return( sint_read( line, elem, size ) );
      case 5: return( ulong_read( line, elem, size ) );
      case 6: return( slong_read( line, elem, size ) );
      case 7: return( float_read( line, elem, size ) );
      case 8: return( double_read( line, elem, size ) );
      case 9: return( string_read( line, elem, size ) );

      default: return( NULL );
   }
}

/* ----------------------------------------------------------------------- */
/* Name       cmpinfo                                                      */
/* Definition void cmpinfo();                                              */
/* Funktion   Ausgabe einer Uebersicht mit den Standardtypen                */
/* ----------------------------------------------------------------------- */

void cmpinfo( void ) {

   printf("   ss                                                                     \n");
   printf("     DBMS   V.1.1       Compare Functions      (c) Juergen Mueller 1992  \n");
   printf("                                                                         \n");
   printf("     Number   Compares                                       Type        \n");
   printf("                                                                         \n");
   printf("        1     Arrays of unsigned byte                        Lex.        \n");
   printf("                                                                         \n");
   printf("        2     Arrays of signed byte                          Lex.        \n");
   printf("                                                                         \n");
   printf("        3     Arrays of unsigned int                         Lex.        \n");
   printf("                                                                         \n");
   printf("        4     Arrays of signed int                           Lex.        \n");
   printf("                                                                         \n");
   printf("        5     Arrays of unsigned long                        Lex.        \n");
   printf("                                                                         \n");
   printf("        6     Arrays of signed long                          Lex.        \n");
   printf("                                                                         \n");
   printf("        7     Arrays of float                                Lex.        \n");
   printf("                                                                         \n");
   printf("        8     Arrays of double                               Lex.        \n");
   printf("                                                                         \n");
   printf("        9     Strings ('\\0' terminated)                      strcmp      \n");
   printf("                                                                        ss\n");

   return;
}


/* ----------------------------------------------------------------------- */
/*                                                                         */
/* Standard-Vergleichsfunktionen                                           */
/*                                                                         */
/* uchar_cmp  - unsigned char arrays                                       */
/* schar_cmp  - signed char arrays                                         */
/* uint_cmp   - unsigned int arrays                                        */
/* sint_cmp   - signed int arrays                                          */
/* ulong_cmp  - unsigned long arrays                                       */
/* slong_cmp  - unsigned long arrays                                       */
/* floar_cmp  - float arrays                                               */
/* double_cmp - double arrays                                              */
/* string_cmp - strings terminated by '\0'                                 */
/*                                                                         */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/* ubyte                                                                   */
/* ----------------------------------------------------------------------- */

int ubyte_cmp( void *data1, void *data2, unsigned length ) {

   int       res = 0;
   unsigned  n, count;
   unsigned  d1len, d2len;
   ubyte    *d1, *d2;

   if ( length > 0 )   {
      for ( n=0; n<length; n++ ) {
         if( ((ubyte*)data1)[n] > ((ubyte*)data2)[n] ) return(  1 );
         if( ((ubyte*)data1)[n] < ((ubyte*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (ubyte*)data1 + sizeof(int); d1len = *(unsigned*)data1;
      d2 = (ubyte*)data2 + sizeof(int); d2len = *(unsigned*)data2;

      count = (d1len>d2len)?(d2len):(d1len);

      for ( n=0; n<count; n++ )   {
         if ( d1[n] < d2[n] ) { res = -1; break; }
         if ( d1[n] > d2[n] ) { res =  1; break; }
      }
      if ( res == 0 )   {
         if ( d1len > d2len ) res = 1;
         if ( d1len < d2len ) res = -1;
      }
   }
   return( res );
}

void *ubyte_read( char *line, unsigned *elem, unsigned *size ) {

   char   buffer[10];
   ubyte *a = NULL;
   char  *h = NULL;
   int    l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(char);
      h = (char*)malloc(*size);
      if ( h != NULL ) {
         a = (ubyte*)(h+sizeof(int));
         *(unsigned*)h = *elem;
      }
   } else {
      *size = *elem * sizeof(char);
      h = (char*)malloc(*size);
      if ( h != NULL ) a = (ubyte*)h ;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf( line+p, "%s ", buffer );
         l = atoi(buffer);
         if ( l < 0 || l > 255 ) { free( h ); h = NULL; break; }
         a[n] = (ubyte)l;
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* sbyte                                                                   */
/* ----------------------------------------------------------------------- */

int sbyte_cmp( void *data1, void *data2, unsigned length ) {

   int       res = 0;
   unsigned  n, count;
   unsigned  d1len, d2len;
   sbyte    *d1, *d2;

   if ( length > 0 ) {
      for( n = 0; n < length; n++ )
      {
         if( ((sbyte*)data1)[n] > ((sbyte*)data2)[n] ) return(  1 );
         if( ((sbyte*)data1)[n] < ((sbyte*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (sbyte*)data1 + sizeof(int); d1len = *(unsigned*)data1;
      d2 = (sbyte*)data2 + sizeof(int); d2len = *(unsigned*)data2;

      count = (d1len>d2len)?(d2len):(d1len);
      for( n=0; n < count; n++ )   {
         if ( d1[n] < d2[n] ) { res = -1; break; }
         if ( d1[n] > d2[n] ) { res =  1; break; }
      }
      if ( res == 0 ) {
         if ( d1len > d2len ) res = 1;
         if ( d1len < d2len ) res = -1;
      }
   }
   return( res );
}

void *sbyte_read( char *line, unsigned *elem, unsigned *size ) {

   char   buffer[10];
   sbyte *a = NULL;
   char  *h = NULL;
   int    l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(char);
      h = (char*)malloc(*size);
      if ( h != NULL ) {
         a = (sbyte*)(h+sizeof(int));
         *(unsigned*)h = *elem;
      }
   } else {
      *size = *elem * sizeof(char);
      h = (char*)malloc(*size);
      if ( h != NULL ) a = (sbyte*)h ;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf( line+p, "%s ", buffer );
         l = atoi(buffer);
         if ( l < 0 || l > 255 ) { free( h ); h = NULL; break; }
         a[n] = (sbyte)l;
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* uint                                                                    */
/* ----------------------------------------------------------------------- */

int uint_cmp( void *data1, void *data2, unsigned length ) {

   unsigned  n, count;
   unsigned  d1len, d2len;
   unsigned *d1, *d2;

   if ( length > 0 )   {
      for( n = 0; n < length; n++ )   {
         if( ((unsigned*)data1)[n] > ((unsigned*)data2)[n] ) return(  1 );
         if( ((unsigned*)data1)[n] < ((unsigned*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (unsigned*)((char*)data1 + sizeof(int)); d1len = *(unsigned*)data1;
      d2 = (unsigned*)((char*)data2 + sizeof(int)); d2len = *(unsigned*)data2;

      count = (d1len > d2len) ? d2len : d1len;
      for( n = 0; n < count; n++ ) {
         if ( d1[n] > d2[n] ) return(  1 );
         if ( d1[n] < d2[n] ) return( -1 );
      }
      if ( d1len > d2len ) return(  1 );
      if ( d1len < d2len ) return( -1 );
   }
   return( 0 );
}

void *uint_read( char *line, unsigned *elem, unsigned *size ) {

   char      buffer[15];
   char     *h = NULL;
   unsigned *a = NULL;
   int       l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(int);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) {
         *(unsigned*)h = *elem;
         a = (unsigned*)(h+sizeof(int));
      }
   } else {
      *size = *elem * sizeof(int);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) a = (unsigned*)h;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf(line+p, "%s ", buffer );
         if ( sscanf(buffer,"%u", a+n ) != 1 ) { free( h ); h = NULL; break; }
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* sint                                                                    */
/* ----------------------------------------------------------------------- */

int sint_cmp( void *data1, void *data2, unsigned length ) {
   unsigned  n, count;
   unsigned  d1len, d2len;
   int      *d1, *d2;

   if ( length > 0 ) {
      for( n = 0; n < length; n++ )   {
         if( ((int*)data1)[n] > ((int*)data2)[n] ) return(  1 );
         if( ((int*)data1)[n] < ((int*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (int*)((char*)data1 + sizeof(int)); d1len = *(unsigned*)data1;
      d2 = (int*)((char*)data2 + sizeof(int)); d2len = *(unsigned*)data2;

      count = (d1len > d2len) ? d2len : d1len;
      for( n = 0; n < count; n++ ) {
         if ( d1[n] > d2[n] ) return(  1 );
         if ( d1[n] < d2[n] ) return( -1 );
      }
      if ( d1len > d2len ) return(  1 );
      if ( d1len < d2len ) return( -1 );
   }
   return( 0 );
}

void *sint_read( char *line, unsigned *elem, unsigned *size ) {

   char    buffer[15];
   char   *h = NULL;
   int    *a = NULL;
   int     l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(int);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) {
         *(unsigned*)h = *elem;
         a = (int*)(h+sizeof(int));
      }
   } else {
      *size = *elem * sizeof(int);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) a = (int*)h;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf(line+p, "%s ", buffer );
         if ( sscanf(buffer,"%d", a+n ) != 1 ) { free( h ); h = NULL; break; }
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* ulong                                                                   */
/* ----------------------------------------------------------------------- */

int ulong_cmp( void *data1, void *data2, unsigned length ) {

   unsigned  n, count;
   unsigned  d1len, d2len;
   ulong    *d1, *d2;

   if ( length > 0 ) {
      for( n = 0; n < length; n++ )
      {
         if( ((ulong*)data1)[n] > ((ulong*)data2)[n] ) return(  1 );
         if( ((ulong*)data1)[n] < ((ulong*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (ulong*)((char*)data1 + sizeof(int)); d1len = *(unsigned*)data1;
      d2 = (ulong*)((char*)data2 + sizeof(int)); d2len = *(unsigned*)data2;

      count = (d1len > d2len) ? d2len : d1len;
      for( n = 0; n < count; n++ ) {
         if ( d1[n] > d2[n] ) return(  1 );
         if ( d1[n] < d2[n] ) return( -1 );
      }
      if ( d1len > d2len ) return(  1 );
      if ( d1len < d2len ) return( -1 );
   }
   return( 0 );
}

void *ulong_read( char *line, unsigned *elem, unsigned *size ) {

   char   buffer[15];
   char  *h = NULL;
   ulong *a = NULL;
   int    l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(ulong);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) {
         *(unsigned*)h = *elem;
         a = (ulong*)(h+sizeof(int));
      }
   } else {
      *size = *elem * sizeof(long);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) a = (ulong*)h;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf(line+p, "%s ", buffer );
         if ( sscanf(buffer,"%lu", a+n ) != 1 ) { free( h ); h = NULL; break; }
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* slong                                                                   */
/* ----------------------------------------------------------------------- */

int slong_cmp( void *data1, void *data2, unsigned length ) {

   unsigned  n, count;
   unsigned  d1len, d2len;
   long     *d1, *d2;

   if ( length > 0 )   {
      for ( n = 0; n < length; n++ ) {
         if( ((long*)data1)[n] > ((long*)data2)[n] ) return(  1 );
         if( ((long*)data1)[n] < ((long*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (long*)((char*)data1 + sizeof(int)); d1len = *(unsigned*)data1;
      d2 = (long*)((char*)data2 + sizeof(int)); d2len = *(unsigned*)data2;

      count = (d1len > d2len) ? d2len : d1len;
      for ( n = 0; n < count; n++ )   {
         if ( d1[n] > d2[n] ) return(  1 );
         if ( d1[n] < d2[n] ) return( -1 );
      }
      if ( d1len > d2len ) return(  1 );
      if ( d1len < d2len ) return( -1 );
   }
   return( 0 );
}

void *slong_read( char *line, unsigned *elem, unsigned *size ) {

   char    buffer[20];
   char   *h = NULL;
   long   *a = NULL;
   int     l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(long);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) {
         *(unsigned*)h = *elem;
         a = (long*)(h + sizeof(int));
      }
   } else {
      *size = *elem * sizeof(long);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) a = (long*)h;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf(line+p, "%s ", buffer );
         if ( sscanf(buffer,"%ld", a+n ) != 1 ) { free( h ); h = NULL; break; }
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* float                                                                   */
/* ----------------------------------------------------------------------- */

int float_cmp( void *data1, void *data2, unsigned length ) {

   unsigned  n, count;
   unsigned  d1len, d2len;
   float    *d1, *d2;

   if ( length > 0 ) {
      for ( n = 0; n < length; n++ ) {
         if( ((float*)data1)[n] > ((float*)data2)[n] ) return(  1 );
         if( ((float*)data1)[n] < ((float*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (float*)((char*)data1 + sizeof(int)); d1len = *(unsigned*)data1;
      d2 = (float*)((char*)data2 + sizeof(int)); d2len = *(unsigned*)data2;
      count = (d1len > d2len) ? d2len : d1len;
      for( n = 0; n < count; n++ ) {
         if ( d1[n] > d2[n] ) return(  1 );
         if ( d1[n] < d2[n] ) return( -1 );
      }
      if ( d1len > d2len ) return(  1 );
      if ( d1len < d2len ) return( -1 );
   }
   return( 0 );
}

void *float_read( char *line, unsigned *elem, unsigned *size ) {

   char   buffer[25];
   char  *h = NULL;
   float *a = NULL;
   int    l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(float);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) {
         *(unsigned*)h = *elem;
         a = (float*)(h+sizeof(int));
      }
   } else {
      *size = *elem * sizeof(float);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) a = (float*)h;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf(line+p, "%s ", buffer );
         if ( sscanf(buffer,"%f", a+n ) != 1 ) { free( h ); h = NULL; break; }
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* double                                                                  */
/* ----------------------------------------------------------------------- */

int double_cmp( void *data1, void *data2, unsigned length ) {

   unsigned  n, count;
   unsigned  d1len, d2len;
   double   *d1, *d2;

   if ( length > 0 ) {
      for ( n = 0; n < length; n++ ) {
         if( ((double*)data1)[n] > ((double*)data2)[n] ) return(  1 );
         if( ((double*)data1)[n] < ((double*)data2)[n] ) return( -1 );
      }
   } else {
      d1 = (double*)((char*)data1 + sizeof(int)); d1len = *(unsigned*)data1;
      d2 = (double*)((char*)data2 + sizeof(int)); d2len = *(unsigned*)data2;
      count = (d1len > d2len) ? d2len : d1len;
      for( n = 0; n < count; n++ ) {
         if ( d1[n] > d2[n] ) return(  1 );
         if ( d1[n] < d2[n] ) return( -1 );
      }
      if ( d1len > d2len ) return(  1 );
      if ( d1len < d2len ) return( -1 );
   }
   return( 0 );
}

void *double_read( char *line, unsigned *elem, unsigned *size ) {

   char    buffer[25];
   char   *h = NULL;
   double *a = NULL;
   int     l,n,p;

   if ( *elem == 0 ) {
      p = 0; l = strlen( line );
      while ( p < l ) {
         if ( sscanf(line+p, "%s ", buffer ) == EOF ) break;
         (*elem)++; p += strlen(buffer)+1;
      }
      *size = sizeof(int) + *elem * sizeof(double);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) {
         *(unsigned*)h = *elem;
         a = (double*)(h+sizeof(int));
      }
   } else {
      *size = *elem * sizeof(long);
      h = (char*)memory_alloc(*size);
      if ( h != NULL ) a = (double*)h;
   }

   if ( a != NULL ) {
      p = 0;
      for ( n=0; n<*elem; n++ ) {
         sscanf(line+p, "%s ", buffer );
         if ( sscanf(buffer,"%lf", a+n ) != 1 ) { free( h ); h = NULL; break; }
         p += strlen(buffer)+1;
      }
   }

   return( h );
}

/* ----------------------------------------------------------------------- */
/* String-Typ                                                              */
/* ----------------------------------------------------------------------- */

int string_cmp( void *data1, void *data2, unsigned length ) {
   int res;

   res = strcmp( (char*)data1, (char*)data2 );

   if ( res < 0 ) { res = -1; }
   if ( res > 0 ) { res =  1; }

   return( res );
}

void *string_read( char *line, unsigned *elem, unsigned *size ) {
   char *h = NULL;

   if ( *elem == 0 ) { *elem = strlen(line)+1; }
   *size = *elem;
   h = (char*)memory_alloc(*size);
   if ( h != NULL ) memory_copy( h, line, *size );
   return( h );
}

/* ----------------------------------------------------------------------- */
/* Ende von TDEFINI.C                                                      */
/* ----------------------------------------------------------------------- */
