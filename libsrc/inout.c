/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/inout.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: file input and output functions including
 *            error handling
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: inout.c,v 1.24 2001/08/07 11:02:16 cocard Exp favey $
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
 *   HISTORY: 
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/inout.c $
 * 
 * 2     10.12.04 16:14 Nobs
 * make compile under darwin
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.24  2001/08/07 11:02:16  cocard
 * *** empty log message ***
 *
 * Revision 1.23  2001/07/03 09:16:06  cocard
 * *** empty log message ***
 *
 * Revision 1.22  2001/07/03 08:47:13  cocard
 * *** empty log message ***
 *
 * Revision 1.21  2001/02/21 14:26:58  cocard
 * *** empty log message ***
 *
 * Revision 1.20  2001/02/21 14:00:19  cocard
 * *** empty log message ***
 *
 * Revision 1.19  2000/08/15 16:38:13  favey
 * extended openzipped to look for .gz
 *
 * Revision 1.18  2000/07/25 13:18:46  favey
 * introducing OS_AIX
 *
 * Revision 1.17  2000/03/28 17:08:54  cocard
 * *** empty log message ***
 *
 * Revision 1.16  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.15  1999/08/17 12:59:14  cocard
 * *** empty log message ***
 *
 * Revision 1.14  1999/06/30 13:21:50  favey
 * bugfixes with openzipped
 *
 * Revision 1.13  1999/06/30 12:17:30  favey
 * new function iszipped()
 *
 * Revision 1.12  1999/06/21 12:22:44  favey
 * bugfix in get_fname_wopath
 *
 * Revision 1.11  1999/06/02 16:13:19  cocard
 * *** empty log message ***
 *
 * Revision 1.10  1999/02/10 18:03:59  favey
 * changes to make commonlib compile on Solaris
 *
 * Revision 1.9  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.8  1998/07/21 15:18:10  favey
 * added new function overread_comment()
 *
 * Revision 1.7  1998/07/17 16:24:11  favey
 * changed various function calls and cleaned up
 *
 * Revision 1.6  1998/07/16 11:44:31  favey
 * *** empty log message ***
 *
 * Revision 1.5  1998/07/16 11:00:33  mace
 * added new function called add_extension
 *
 * Revision 1.4  1998/07/03 11:09:11  cocard
 *  strpos_withblank added
 *
 * Revision 1.3  1998/05/06 09:43:19  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1998/05/04 10:05:36  favey
 * extract_values_static and read_line_static added, please use those now
 *
 * Revision 1.1  1998/04/22 12:31:28  favey
 * Initial revision
 *
 *
 * DOCUMENTATION: pipe this file into
 *
 *            grep " \*% " | cut -c 5- | nroff -man | less -is
 *
 *            to view the manpages describing the functions
 *            (That's how I do it under Linux)
 *
 */

/*
 *  standard C libraries
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>


#include "global.h"
#include "chain.h"
/* we need this up here, because it defines OS_LINUX */

#if defined(OS_LINUX) || defined(OS_AIX) || defined(OS_CYGWIN)

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#endif


#ifdef WITH_ERRORS
#include "errno.h"
#include "inout.h"
#endif

#include "errproc.h"


/*
 *  Overread all characters up to a comment specifier.
 */
int read_comment (FILE *fpin,
		  const int idchar)
{
  int c;

#ifdef WITH_ERRORS
  clearerr (fpin);
#endif

  while ((c=getc(fpin)) != EOF)
    if (c==idchar) return 0;
#ifdef WITH_ERRORS
  if (ferror (fpin) != 0)
    {
      fprintf (stderr, "GGL Common Library -- inout\n");
      fprintf (stderr, "read_comment: could not read file: %s\n",
	       strerror(errno));
    }
#endif
  return EOF;
} 

/*
 *  Overread comment lines, assuming next char read is first of a line.
 *  RETURN VALUE: zero on success, EOF upon end of file.
 */
int overread_comment (FILE *fpin,	/* input stream */
		      FILE *fpout,	/* output: if NULL, lines are discarded */
		      const int commentchar)
{
  int c, d;

  do {
    if ((c=getc(fpin)) == EOF) return EOF;
    if (c==commentchar) {
      /* overread entire line, and write to output if avail. */
      if (fpout) putc(c, fpout);
      while ((d=getc(fpin)) != EOF) {
	if (fpout) putc(d, fpout);
	if (d=='\n') break;
      }
      if (d==EOF) return EOF;
    }
  } while (c==commentchar);
  ungetc (c, fpin);
  return 0;
}

/*
 *  read option entry following a specific character
 */
char * read_option (FILE *fpin,
		 const int idchar)
{
  char *option_line;

  if (read_comment (fpin, idchar) != 0) return NULL;
  if (read_line (fpin, &option_line) == EOF) return NULL;
  return option_line;
}

/*
 *  read string option entry following a specific character
 */
char * read_soption (FILE *fpin,
		 const int idchar)
{
  char *option_line;
  char soption[256];

  if (read_comment (fpin, idchar) != 0) return NULL;
  if (read_line (fpin, &option_line) == EOF) return NULL;
  if (strlen (option_line) > 255) option_line[255] = '\0';
  if (sscanf (option_line, " %s", soption) < 1) return NULL;
  free (option_line);
  return strdup (soption);
}

/*
 * read_line, read_delim \- read an array of char upto a newline or delimiter
 */
int read_delim (FILE *fpin,
		int delimiter,
		char **lineptr)
{
#if defined(OS_SOLARIS) || defined(OS_AIX) || defined(OS_CYGWIN) || defined(__APPLE__)
  /*
   *  Sadly enough, Solaris does not support `getdelim',
   *  so we need to use fgets instead. Note that in this case, the allocation
   *  space is set to 512 characters, this is not really as it should be.
   *  Beware, the delimiter is not looked at under Solaris, so read_delim()
   *  really converts to read_line().
   *
   *  Consider it a quick fix to get it running under Solaris.
   */
  int n = 511;
  if (delimiter != '\n')
    {
      fprintf (stderr, "GGL Common Library --- inout:\n");
      fprintf (stderr, "read_delim: Solaris version: ");
      fprintf (stderr, "delimiter ignored! Using newline instead.\n");
    }
  *lineptr = (char *) malloc ((n+1)*sizeof(char));
  fgets (*lineptr, n, fpin);
  return (feof (fpin) ? EOF : 0);
#else
  /*
   *  We make use of the GNU extension `getdelim', this may or may not
   *  be supported, when you port to other platforms.
   */
  size_t n = 0;
  *lineptr = NULL;
  return (int) getdelim (lineptr, &n, delimiter, fpin);
#endif
}


int read_line (FILE *fpin,
	       char **lineptr)
{
  return read_delim (fpin, '\n', lineptr);
}

/*
 *  Read entire line into static string buffer
 */
int read_line_static (FILE *fpin, char *s, const int size)
{
	char *tmp;
	tmp = fgets (s, size, fpin);
	if (tmp == NULL) return EOF;       
	return strlen(s);
}


int strpos_withblanks(const char *s, const int count)
{
  int i=count;
  while (i>0 && isspace(s[i-1])) i--;
  return i;
}

/*
 *  For the following former inout funcions:
 *  FileExists, CheckFileName
 *  the standard library function access from unistd.h may be used:
 *  int access (const char *filename, int how)
 *   where how may be R_OK, W_OK, X_OK, F_OK
 *
 *  I've still included them here, but they are obsolete and should
 *  not be used anymore. That's why I did not bother to write any
 *  documentation for it.
 *
 *  They may get removed in later versions!
 */


/* DEBUG XXX still used in gpskin */
int check_file_name (const char *filename,
		     int mode)
{
  char atype;
  int err;
  switch (mode)
    {
    case 'w': 
    case W_OK:
      mode = W_OK;
      atype = 'w';
      break;
    case 'r':
    case R_OK:
      mode = R_OK;
      atype = 'r';
      break;
    case 'x':
    case X_OK:
      mode = X_OK;
      atype = 'x';
      break;
    default:
#ifdef WITH_ERRORS
      fprintf (stderr, 
        "GGL Common Library -- inout\ncheck_file_name: invalid option\n");
#endif
      return -1;
    }
  err = access (filename, mode);

#ifdef WITH_ERRORS
  if (err != 0)
    {
      if (errno != ENOENT)
	{
	  fprintf (stderr, "GGL Common Library -- inout\n");
	  fprintf (stderr,
	  "check_file_name: could not access file %s in mode %c: %s\n",
		   filename,
		   atype,
		   strerror(errno));
	}
    }
#endif
  return err;  
}


/* 
 * read column file 
 */
void ** read_colfile(FILE *fp, int *count, void *param, readlinefunction readline)
{
  char line[MAXLINELENGTH];
  struct Node *chain;
  void *data;
  void **array;
  int linenr = 0;

  chain=chain_init(count);
  while (!feof(fp))
    {
       linenr++;
       if (read_line_static(fp, line, MAXLINELENGTH) <= 0) continue;
 
      /* we overread comment and jump to the next line */
      if  (iscomment(line)) continue;
      
      /* read data from line */
      if ((data = readline(line, param)) == NULL) 
	{
	  fprintf(stderr,"GGL Common Library -- inout\n");
	  fprintf(stderr,"read_colfile: error on line %d\n", linenr);
	  return NULL;
	}
      /* include data in list */
      chain=chain_append(chain, chain_newelement(data), count);
    }
  
  array=(void **) chain_to_array(chain, *count);
  chain_free(chain, FALSE, *count);
  return array;     
 
}

/* 
 * read file of column type 
*/
void ** read_colfile2(char *filename, int *count, void *param, readlinefunction readline)
{
  FILE *fp;
  void **array;

  if ((fp=open_text_file(filename, "r")) == NULL) return NULL;
  if ((array = read_colfile(fp, count, param, readline)) == NULL)
    {
      fprintf(stderr,"GGL Common Library -- inout\n");
      fprintf(stderr,"read_colfile2: error in file %s\n", filename);
    }
  fclose(fp);
  return array;
}

/* 
 * check on comment 
*/
int iscomment(char *line)
{
  unsigned int i;
  /* lines starting with #  */
  if (line[0]=='#') return 1;

  /* blank lines are also considered as comments */
  for (i=0;i<strlen(line); i++)
    if (!isspace(line[i])) return 0;
  return 1;
}

/*
 *  enhanced file removal
 */
int erase_file (const char *filename)
{
  int err = (unlink (filename));
#ifdef WITH_ERRORS
  if (err != 0)
    {
      fprintf (stderr, "GGL Common Library -- inout\n");
      fprintf (stderr, "erase_file: could not erase file %s: %s\n",
	       filename,
	       strerror(errno));
    }
#endif
  return err;
}

/*
 *  enhanced file opening
 */
FILE * open_text_file (const char *fname,
		       const char *mode)
{
  FILE *fp;
  fp = fopen (fname, mode);
#ifdef WITH_ERRORS
  if (!fp)
    {
      fprintf (stderr, "GGL Common Library -- inout\n");
      fprintf (stderr, "open_text_file: could not open file %s : %s\n",
	       fname,
	       strerror(errno));
    }
#endif
  return fp;
}

FILE * open_binary_file (const char *fname,
			 const char *mode)
{
  FILE *fp;
  fp = fopen (fname, mode);
#ifdef WITH_ERRORS
  if (!fp)
    {
      fprintf (stderr, "GGL Common Library -- inout\n");
      fprintf (stderr, "open_binary_file: could not open file %s : %s\n",
	       fname,
	       strerror(errno));
    }
#endif
  return fp;
}


/*
 * check, whether the file *fp is a zipped file by looking at its magic
 * *fp must already be open for reading!
 * - if the file is zipped, the function closes it, 
 *   and a pipe on the unzipped file is passed back.
 * - if the file is unzipped, the filepointer is left untouched
 * RETURN VALUES: TRUE is zipped, FALSE if unzipped, EOF on error
 */
FILE * openzipped (const char *fname, int *ispipe)
{
	int gza, gzb;
	char *fname2 = NULL;
	FILE *fp;
	fp = fopen (fname, "r");
	if (!fp) {
	  fname2 = (char *) malloc ((strlen(fname)+5)*sizeof(char));
	  if (fname2==NULL) { *ispipe = EOF; return NULL; }
	  sprintf (fname2,"%s.gz",fname);
	  fp = fopen (fname2, "r");
	  if (!fp) { 
	    puterr ("inout","openzipped");
	    fprintf (stderr, "could not open neither %s nor %s : %s\n",
		     fname, fname2, strerror(errno));
	    free (fname2); 
	    *ispipe = EOF; 
	    return NULL; }
	} else {
	  fname2 = strdup (fname);
	}

	gza = getc(fp); gzb = getc(fp);
	
	if (gza == 0x1f && gzb == 0x8b) {
		/* gzipped data determined by magic, use pipe */
		char cmd[1024];
		fclose (fp);
		sprintf (cmd, "/bin/zcat %s", fname2);
		free (fname2);
		fp = popen (cmd, "r");
		if (fp == NULL) {
			fprintf (stderr, "GGL Common Library -- inout\n");
			fprintf (stderr,
				 "openzipped: could not open pipe %s : %s\n",
				 cmd, strerror(errno));
			*ispipe = EOF;
			return NULL;
		}
		*ispipe = TRUE;
		return fp;
	}
	free (fname2);
	if (isprint(gza) && isprint(gzb)) {
		rewind (fp);
		*ispipe = FALSE;
		return fp;
	}
	fprintf (stderr, "GGL Common Library -- inout\n");
	fprintf (stderr, "openzipped: apparently neither zipped nor ASCII\n");
	*ispipe = EOF;
	return NULL;
}


/*
 *  filename extension operations
 */


/*
 *  this function modifies the string passed to it!
 *  the extension including its leading point will
 *  be removed (in fact: the point is replaced by a
 *  null character)
 */
void remove_extension (char *filename)
{
  char * tmp = strrchr (filename, '.');
  if (tmp == NULL) return;     /* there was no extension */
  tmp[0] = '\0';           /* this seems to be a kind of hack ;-) */
}

/*
 *  this function returns a new allocated string
 *  containing the extension including the leading point
 */
char * get_extension (const char *filename)
{
  char * tmp = strrchr (filename, '.');
  if (tmp == NULL) return NULL;
  return strdup(tmp);
}


/*
 *  this function returns a new allocated string
 *  containing the filename with a new extension
 */
char * add_extension (
	const char *oldname,
	const char *ext )
{
	char *tempname;

	if ((tempname =(char *)malloc(strlen(oldname)+strlen(ext)+2))==NULL)
		return NULL;
	strcpy(tempname, oldname);
	strcat(tempname, ".");
	strcat(tempname, ext);
	return tempname;
}


/*
 *  this function returns a new allocated string
 *  containing the filename excluding the leading path
 */
char * get_fname_wopath (const char *filename)
{
  char * tmp = strrchr (filename, '/');
  if (tmp == NULL) return strdup(filename);
  return strdup(&(tmp[1]));
}



/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void uppercase (char *s)
{
  char *tmp = s;
  while (*tmp!='\0')
    {
      *tmp = (char) toupper (*tmp); 
      tmp++;
    }
}

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void lowercase (char *s)
{
  char *tmp = s;
  while (*tmp!='\0')
    {
      *tmp = (char) tolower (*tmp); 
      tmp++;
    }
}

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void replacechar (char *s, int oldc, int newc)
{
  char *tmp = s;
  while ((tmp = strchr(tmp, oldc)) != NULL)
    tmp[0] = (char) newc;
}

/*
 * returns offset  in string
 */
int col_offset(char * s, int col)
{
  char * tmp = s;
  int i;
  
  if (col==0) return 0;
  while (strchr(" \t,\n", (int)tmp[0])) tmp++;   /* find start of first column */

  for (i=0; i<col; i++)   /* overread col columns */
    {
      while (!strchr (" ,\t\r\f\v\n", (int)tmp[0])) tmp++; /* overread column */
      while ( strchr (" ,\t\r\f\v",   (int)tmp[0])) tmp++; /* find beginnig of next column */
      if (tmp[0] == '\n')
	{
	  fprintf (stderr, "GGL Common Lib -- col_offset\n");
	  fprintf (stderr, "not enough columns on line\n");
	  fprintf (stderr, "line:%s\n", s);
	  return -1;
	}
    }
  return (int) (tmp - s);
}

/*
 *  extract columns from an input line
 */
int extract_values (
	const char *line,
	char **prefix,
	const unsigned int prefix_size,
	struct coordcol *infix,
	char **suffix,
	const unsigned long linecount)
{
	const char delimiters[] = " \f\t\v\r\n";
	char *ank, *element, *token;
	int i;

	if (strlen(line) < prefix_size) {
#ifdef WITH_ERRORS
		fprintf (stderr, "GGL Common Lib -- parse\n");
		fprintf (stderr, "extract_values: %s%ld\n",
			 "prefix size bigger than length of line: line ",
			 linecount);
#endif
		return -1;
	}
	if (infix->count <= 0 || infix->count > 7) {
#ifdef WITH_ERRORS
		fprintf (stderr, "GGL Common Lib -- parse\nextract_values: %s\n",
			 "invalid column count");
#endif
	}
	if (prefix_size > 0) {
               char *tmppre;
		tmppre = (char *) malloc ((prefix_size+1) * sizeof (char));
		tmppre[prefix_size] = '\0';
		strncpy (tmppre, line, prefix_size);
		*prefix = tmppre;
	} else {
		*prefix = NULL;
	}
	ank = element = strdup (line + prefix_size);
	token = strtok (element, delimiters);
	i=0;
	while (token != NULL && i < infix->count) {
		infix->col[i] = strdup (token);
		i++;
		if (i == infix->count)
			token = strtok (NULL, "\r\n");
		else
			token = strtok (NULL, delimiters);
	}
	free (ank);
	if (i < infix->count) {
#ifdef WITH_ERRORS
		fprintf (stderr, "GGL Common Lib -- parse\nextract_values: ");
		fprintf (stderr,
			 "not enough valid tokens in line %ld, required: %d, found %d\n",
			 linecount, infix->count, i);
#endif
		return -1;
	}
	*suffix = (token != NULL) ? strdup(token):(char*)NULL;
	return 0;
}


/*
 *  extract columns from an input line
 *  this static version is for faster processing without
 *  new memory allocation
 */
int extract_values_static (
	const char *line,
	char prefix[],
	const int prefix_size,
	struct scoordcol *infix,
	char suffix[],
	const int max_suffix_size,
	const unsigned long linecount)
{
	int i, j, ccount, state, linelen;

	if ((linelen = strlen(line)) < prefix_size) {
#ifdef WITH_ERRORS
		fprintf (stderr, "GGL Common Lib -- parse\n");
		fprintf (stderr, "extract_values_static: %s%ld\n",
			 "prefix size bigger than length of line: line ",
			 linecount);
#endif
		return -1;
	}
	if (prefix_size > 0) {
		strncpy (prefix, line, prefix_size);
		prefix[prefix_size] = '\0';
	}
	ccount = state = j = 0;
	i = prefix_size;
	while (i<linelen && ccount < infix->count && j < infix->slen) {
		if (isspace(line[i])) {
			if (state) {
				state = j = 0;
				ccount++;
			}
			i++;
			continue;
		} else {
			state = 1;
			infix->col[ccount][j++] = line[i++];
		}
	}


	if (ccount < infix->count) {
#ifdef WITH_ERRORS
		fprintf (stderr, 
			 "GGL Common Lib -- parse\nextract_values_static: ");
		fprintf (stderr,
			 "not enough valid tokens in line %ld, required: %d, found %d\n",
			 linecount, infix->count, i);
#endif
		return -1;
	}
	if (i < linelen) {
		strncpy (suffix, line+i, max_suffix_size);
		suffix[max_suffix_size] = '\0';
	} else {
		for (i=0;i<max_suffix_size;i++) suffix[i] = '\0';
	}
	if (suffix[strlen(suffix)-1] == '\n') {
		suffix[strlen(suffix)-1] = '\0';
	}

 	return 0;
}

/*
 *  strip directory and suffix from filenames
 *  returns a string 
 */ 

char * base_name(char *filename, char *ext)
{
  char *s;
  int i, j;

  for (i=strlen(filename); i>0; i--) 
    if (filename[i-1] == '/') break;
  if ((s=strdup(filename+i)) == NULL) return NULL;

  if (ext != NULL)
    {
      i=strlen(s);
      j=strlen(ext);
      while (i*j!=0 && s[--i] == ext[--j]) s[i]='\0';
    }
  
  return s;
}

