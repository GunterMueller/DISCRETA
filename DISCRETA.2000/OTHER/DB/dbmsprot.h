/* ----------------------------------------------------------------------- */
/* DBMSPROT.H                                                              */
/*                                                                         */
/* Prototypen der im Projekt verfuegbaren, internen Funktionen.             */
/* ----------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>

/* ----------------------------------------------------------------------- */
/* memory_... Speichermanagement                                           */
/* ----------------------------------------------------------------------- */
/* #define memory_alloc malloc */
#define memory_alloc(a) calloc(a, 1)
#define memory_free  free
/* #define memory_copy  memmove */
#define memory_copy Memmove
#define memory_zero Memory_zero

/* ----------------------------------------------------------------------- */
/* INIT.C  Initialisierung                                                 */
/* ----------------------------------------------------------------------- */
int  init_table( TABLE **table, char *name );
void drop_table( TABLE *table );

/* ----------------------------------------------------------------------- */
/* INSERTH.C - Hilfsfunktionen zum Einfuegen                                */
/* INSERTO.C - Einfuegen auf Overlayseiten                                  */
/* SPLIT.C   - Spalten von Datenseiten                                     */
/* SPLITI.C  - Spalten von Indexseiten                                     */
/* SPLITIH.C - Hilfsfunktionen zum Indexseitenspalten                      */
/* ----------------------------------------------------------------------- */
int   create_new_overlay( TABLE *table, UINT num );
int   insert_overlay( TABLE *table, UINT num, ULONG *pages_to_del, UINT *del_cnt );
void  put_R_comp( UCHAR *overlay, UCHAR *r_comp, UINT r_size, DINFO *di );
ULONG get_next_page( TABLE *table, UINT *next );
int   insert_data( TABLE *table, PAGEBUF *pagebuf, UINT *last, ULONG *pages_to_del, UINT *del_c, ULONG *offset );
int   split_data_page( TABLE *table, DPS *dps, UINT in );
int   d_new_root( TABLE *table, ULONG *offset, DPS *dps );
int   d_insert_index( TABLE *table, ULONG *offset, DPS *dps, UINT onum );

/*--------------------------------------------------------------------------*/
/* DELETE                                                                   */
/*--------------------------------------------------------------------------*/
UINT delete_data( TABLE *table, ULONG *offset, ULONG *pages_to_del, UINT *del_c, USER_RIGHTS *user );

/* ----------------------------------------------------------------------- */
/* BASEROUT.C - Basisfunktionen, werden in allen Modulen benoetigt          */
/* ----------------------------------------------------------------------- */
void  newPage( TABLE *table, UINT type );
ULONG getFreeFilePosition( TABLE *table );
void  getIndexOrder( TABLE *table, UINT *I );
UINT  compSize( TABLE *table, UINT comp, UCHAR *data );
void  copyData( TABLE *table, UCHAR *data, DINFO *dinfo );
void  makefname( char *fname, char *path, char *filename, char *extension );

/* ----------------------------------------------------------------------- */
/* FILEROUT.C - Funktionen fuer saemtliche Dateioperationen                  */
/* ----------------------------------------------------------------------- */
int file_exists( char *filename );
int create_table( TABLE *table );
int open_table( TABLE *table );
int close_table( TABLE *table );
int delete_table( TABLE *table );
int read_page( TABLE *table );
int write_page( TABLE *table );
int read_base_page( TABLE *table );
int write_base_page( TABLE *table );

/*---< VIEW.C - Anzeige wichtiger Strukturen >----------------------------------------------------*/

void view_table( TABLE *table );
void view_page( TABLE *table );

/* ----------------------------------------------------------------------- */
/* TDEFINI.C - Initialisierung Index, Indexfunktionen                      */
/* ----------------------------------------------------------------------- */
int   get_definition( char *name );
void *read_txt( unsigned type, char *line, unsigned *elem, unsigned *size );
void  cmpinfo( void );

/* ----------------------------------------------------------------------- */
/* Makro fuer die Projektion eines Datensatzes auf die i-te Komponente      */
/* ----------------------------------------------------------------------- */
#define proj(t,i,d) ( (((BASE_PAGE*)((t)->base))->variable)?((d)+*(((UINT*)(d))+i+1)):((d)+(t)->index_db[i].position))

/* ----------------------------------------------------------------------- */
/* Ende von DBMSPROT.H                                                     */
/* ----------------------------------------------------------------------- */
