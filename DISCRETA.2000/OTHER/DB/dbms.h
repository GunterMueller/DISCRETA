/*


    D B M S . H                                              Version 1.1      


   This header files contains the dbms function prototypes and their          
   parameters. This header is required in all modules using the dbms          
   functions.                                                                 

   (c) Juergen Mueller - January 1993                                         

*/
#ifndef __DBMS__
#define __DBMS__

#define byte  unsigned char
#define ubyte unsigned char
#define sbyte signed char
#define uint  unsigned int
#define sint  signed int

#if 0
#define ulong unsigned long
#define slong long
#endif
#define ulong unsigned int
#define slong int
/* AB 950118 */

/*

  CREATE_PARM, Parameter zum Erzeugen von Tabellen                            

*/
typedef struct _create_parm {
   char         *filename;       /* Name der zu erzeugenden Tabelle */
   char          read_access;    /* '0',...,'9'                     */
   char          change_access;  /* '0',...,'9' oder '-'            */
   unsigned      page_size;      /* size of all pages in bytes */
   unsigned      record_length;  /* size of all records, if variable 0 */
   unsigned      index_count;    /* count of comparable components */
   unsigned      info_pos;       /* begin of uncomparable record part */
   unsigned      info_size;      /* size of uncomparable part */
   unsigned      pos[32];        /* component positions, if variable 0 */
   unsigned      elem[32];       /* count of elements in component, variable 0 */
   unsigned      type[32];       /* types of each component */
   unsigned      flags[32];      /* unused */
   char          name[32][10];   /* component names */
} CREATE_PARM;

#define DENIED 255

/*

  INSERT_PARM, parameter for function dbms_insert                             

*/
        typedef struct _insert_parm {

                char *name;
                char *data;
                char  dataaccess;
                int   usernumber;

        } INSERT_PARM;

#define READ0     0
#define READ1     1
#define READ2     2
#define READ3     3
#define READ4     4
#define READ5     5
#define READ6     6
#define READ7     7
#define CHANGE0   0
#define CHANGE1   8
#define CHANGE2  16
#define CHANGE3  24
#define CHANGE4  32
#define CHANGE5  40
#define CHANGE6  48
#define CHANGE7  56


/*

  FIO_PARM, parameter for function dbms_fio                                   

  name       Path+filename of the table, without extension                    
  data       Buffer for records to insert                                     
  readb      Function to read records to insert int data buffer               
  usernumber Identifcation of user                                            

*/
        typedef struct _fio_parm {

                char           *name;
                char           *data;
                int           (*readb)( struct _fio_parm * );
                char            dataaccess;
                unsigned long   count;
                unsigned int    special;
                int             usernumber;

        } FIO_PARM;

#define NO_BREAK_ON_IDENTICALS  1


/*

  SELECT_PARM, parameter for function dbms_select                             

*/
typedef struct _select_parm {
   char            *tblname;
   char            *low;
   char            *up;
   char             mode[32];
   int            (*outfunc)( char * );
   long             count;
} SELECT_PARM;

#define LOW_GREATER   0
#define UP_SMALLER    0
#define LOW_EQUAL     1
#define UP_EQUAL      2
#define LOW_UNLIMITED 4
#define UP_UNLIMITED  8

/*

  UPDATE_PARM, parameter for function dbms_update                             

  name       Path+filename of the table, without extension                    
  old_data   Record with old data to update                                   
  new_data   Record with new data                                             
  usernumber Identifcation of user                                            

*/
        typedef struct _update_parm {

                char *name;
                char *old_data;
                char *new_data;
                int   usernumber;

        } UPDATE_PARM;


/*

  DELETE_PARM, parameter for function db_delete                               

  name       Path+filename of the table, without extension                    
  data       Record to delete                                                 
  usernumber Identifcation of user                                            

*/
        typedef struct _delete_parm {

                char *name;
                char *data;
                int   usernumber;

        } DELETE_PARM;


/*

  Prototypes of the main dbms functions                                       

*/

int dbmsCreate( CREATE_PARM *cp );
int dbmsCreateF( char *defname );

int dbms_insert( INSERT_PARM *ip ); int dbms_insert_f( char *tblname, char *in, int msg );
int dbms_select( SELECT_PARM *sp ); int dbms_select_f( char *defname, int msg );
int dbms_update( UPDATE_PARM *up );
int dbms_delete( DELETE_PARM *dp ); int dbms_remove_f( char *tblname, char *in, int msg );

int dbms_fio( FIO_PARM *fp );

/*

  Prototypes of the help dbms functions                                       

*/
        int dbms_set_counter( char *name, int counter, long value, char *password );
        int dbms_get_counter( char *name, int counter, long *value, char *password );
        int dbms_statistics(  char *name, char *password );

#endif
