/*
 *      FILE: inout
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *  DISCLAIMER:   This program is free software; you can redistribute
 *              it and/or modify it under the terms of the GNU General
 *              Public License as published by the Free Software
 *              Foundation; either version 2 of the License, or
 *              (at your option) any later version.
 *                This program is distributed in the hope that it will
 *              be useful, but WITHOUT ANY WARRANTY; without even the
 *              implied warranty of MERCHANTABILITY or FITNESS FOR A
 *              PARTICULAR PURPOSE. See the GNU General Public License
 *              for more details.
 *                You should have received a copy of the GNU General
 *              Public License along with this program; if not, write
 *              to the Free Software Foundation, Inc., 675 Mass Ave,
 *              Cambridge, MA 02139, USA.
 *
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch>
 *            Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            16-09-1997  --  added parsing functions
 *            
 *
 */
#ifndef GGL_INOUT
#define GGL_INOUT


/* 
 * function for scanning 1 line from a data file 
 * (used in read_colfile)
 */
typedef void * (*readlinefunction)(char *, void *);

/*
 *  Struct to hold char columns
 */
struct coordcol {
        int count;
        char *col[7];
} ;

struct scoordcol {
	int count;
	int slen;
	char **col;
} ;

int read_comment (FILE *fpin, const int idchar)
;

/*
 *  Overread comment lines, assuming next char read is first of a line.
 *  RETURN VALUE: zero on success, EOF upon end of file.
 */
int overread_comment (FILE *fpin,	/* input stream */
		      FILE *fpout,	/* output: if NULL, lines are discarded */
		      const int commentchar);

char * read_option (FILE *fpin, const int idchar);

char * read_soption (FILE *fpin, const int idchar);

/*
 *  The following functions will have to be reviewed when porting
 *  the library to another plattform or using another compiler 
 *  than the GNU C Compiler (gcc)!
 */
#if defined(OS_LINUX) || defined(OS_AIX) || defined(OS_CYGWIN)
int read_delim (FILE *fpin, int delimiter, char **lineptr);

int read_line (FILE *fpin, char **lineptr);
int read_line_static (FILE *fpin, char *s, const int size);

/* 
 * read column file 
 */
void ** read_colfile (FILE *fp,       int *count, void *param, readlinefunction readline);
void ** read_colfile2(char *filename, int *count, void *param, readlinefunction readline);


/* for use with returned int of %n on sscanf */
int strpos_withblanks(const char *s, const int count);

int erase_file (const char *filename);

/*
 *  uses the access function from unistd.h
 */
int check_file_name (const char *filename,  int mode);
FILE * open_text_file (const char *fname, const char *mode);
FILE * open_binary_file (const char *fname, const char *mode);

/*
 * check, whether the file fname is a zipped file by looking at its magic
 * - if the file is zipped, the function closes it, 
 *   and a pipe on the unzipped file is passed back.
 * - if the file is unzipped, the filepointer is left untouched
 * ispipe: TRUE is zipped, FALSE if unzipped, EOF on error
 */
FILE * openzipped (const char *fname, int *ispipe);

#endif /* OS_LINUX */ 

/*
 * returns TRUE or FALSE
 * lines starting with # and blank lines are considered as comments
 */
int iscomment(char *line);

/*
 * skipping col values in string, returning offset 
 */
int col_offset(char * s, int col);


/*
 *  this function modifies the string passed to it!
 *  the extension including its leading point will
 *  be removed (in fact: the point is replaced by a
 *  null character)
 */
void remove_extension (char *filename);

/*
 *  this function returns a new allocated string
 *  containing the extension including the leading point
 */
char * get_extension (const char *filename);

/*
 *  this function returns a new allocated string
 *  containing the filename with a new extension
 */
char * add_extension (
	const char *oldname,
	const char *ext );

/*
 *  this function returns a new allocated string
 *  containing the filename excluding the leading path
 */
char * get_fname_wopath (const char *filename);

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void uppercase (char *s);

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void lowercase (char *s);

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void replacechar (char *s, int oldc, int newc);

/*
 *  extract columns from an input line
 */
int extract_values (
	const char *line,
	char **prefix,
	const unsigned int prefix_size,
	struct coordcol *infix,
	char **suffix,
	const unsigned long linecount);

/*
 *  extract columns from an input line
 */
int extract_values_static (
	const char *line,
	char prefix[],
	const int prefix_size,
	struct scoordcol *infix,
	char suffix[],
	const int max_suffix_size,
	const unsigned long linecount);

/*
 *  strip directory and suffix from filenames
 *  returns a string 
 */ 
char * base_name(char *filename, char *ext);



/*
************************************************************
* parse.c
************************************************************
*/
#define ASCII_LINE 4
#define ASCII_BLOCK 5
#define BINARY_BLOCK 1

/* byte of type meaning */

#define BLOCK 1
/* unused 2 */
#define ASCII 4
/* unused 8 */

#define SKIP_NONE 0
#define SKIP_CHAR 1
#define SKIP_FIELD 2

int parse_line (
  FILE *fpin,
  const char *comment_delimiters,
  const char *column_delimiters,
  const int skip,
  const int skip_type,
  char **irrelevant,
  char **datastring)
;

char * get_fileformat (
  FILE *fpin,
  const int type,
  const char *comment_delimiters,
  const char *column_delimiters)
;

#endif


