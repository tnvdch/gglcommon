/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/parse.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: file input and output functions including
 *            error handling
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: parse.c,v 1.6 2000/07/25 13:18:46 favey Exp favey $
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
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/parse.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.6  2000/07/25 13:18:46  favey
 * introducing OS_AIX
 *
 * Revision 1.5  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.4  1999/06/02 16:13:19  cocard
 * *** empty log message ***
 *
 * Revision 1.3  1998/11/16 17:58:06  favey
 * changed syntax of puterr
 *
 * Revision 1.2  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.1  1998/04/22 12:31:43  favey
 * Initial revision
 *
 * Revision 0.3  1997/12/24 12:11:53  favey
 * followup of V. 0.29
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
#include <stddef.h>
#include <ctype.h>

/*
 *  GGL Common Libs Headers
 */
#include "global.h"

#ifdef WITH_ERRORS
#include <errno.h>
#endif

#if defined(OS_LINUX) || defined(OS_AIX) || defined(OS_CYGWIN)

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#endif

#include "errproc.h"
#include "inout.h"

/*
 * Function declarations
 */

/*
 *% .TH PARSE_LINE 3ggl "26 September 1997" "GGL/IGP ETHZ" "GGL Common Lib"
 *% .SH NAME
 *% parse_line \- extract data from line into given structure
 *% .SH SYNOPSIS
 *% .B #include <common/inout.h>
 *% 
 *% .RB "int " parse_line
 *% .RI "(FILE *" fpin ", const char *" comment_delimiters ", const char *" column_delimiters ", const int " skip ", const int " skip_type ", char **" irrelevant ", char **" datastring ");"
 *% .SH DESCRIPTION
 *% (still TODO)
 *% .I skip_type
 *% may be SKIP_CHAR or SKIP_FIELD.
 *% .SH "RETURN VALUES"
 *% .IP -2
 *% fatal, wrong parameters
 *% .IP EOF
 *% end of file occurred
 *% .IP 0
 *% no valid datastring, irrelevant may contain string.
 *% .IP >0
 *% length of string in datastring
 *% .SH CAVEATS
 *% You have to compile your program using
 *% .I -lcommon
 *% as a compiler flag to properly link
 *% your function to the library.
 *% .SH AUTHORS
 *%   (c) 1997 Geodesy and Geodynamics Lab, ETH Zurich:
 *% .RS
 *% Etienne Favey <favey@geod.ethz.ch>,
 *% .RE
 *% .SH "SEE ALSO"
 *% .BR read_line (3ggl),
 *% .BR get_fileformat (3ggl),
 */
int parse_line (
  FILE *fpin,
  const char *comment_delimiters,
  const char *column_delimiters,
  const  int skip,
  const int skip_type,
  char **irrelevant,
  char **datastring)
{
  char *s;
  int skip_size, cols_count;
  int inwhite, delim_found;
  int ssize;
  int numdelim;
  int empty_line;
  int i;

  /* check parameters */
  if (column_delimiters == NULL)
    column_delimiters = strdup(" \f\n\r\t\v");
  if ((skip_type != SKIP_CHAR &&
       skip_type != SKIP_FIELD  ) ||
      datastring == NULL) {
#ifdef WITH_ERRORS
	  puterr("parse","parse_line");
	  fprintf (stderr, "invalid skip type\n");
#endif
    return -2;
  }

  ssize = read_line (fpin, &s);
  numdelim = strlen(comment_delimiters);
  empty_line = TRUE;

  /* check for comment lines */
  while (empty_line != FALSE && ssize != EOF) {
    empty_line = FALSE;
    for (i=0;i<numdelim;i++) {
      if (s[0] == comment_delimiters[i]) empty_line = TRUE;
    }
    if (s[0] == '\n' || s[0] == '\r') empty_line = TRUE;
    if (empty_line != FALSE) {
      free (s);
      ssize = read_line (fpin, &s);
    }
  }
  if (ssize == EOF) return EOF;

  /* apply skip of irrelevant data */
  switch (skip_type) {
  case SKIP_CHAR:
    if (skip == 0){
      *datastring = strdup (s);
      free (s);
      return strlen(*datastring);
    }
    if (strlen(s) > skip) {
      if (irrelevant != NULL) {
	*irrelevant = (char *) malloc ((skip+1)*sizeof(char));
	memset (*irrelevant, 0, skip+1);
	strncat (*irrelevant, s, skip);
      }
      *datastring = (char *) malloc ((strlen(s)-skip+1)*sizeof(char));
      memset (*datastring, 0, strlen(s)-skip+1);
      strncat (*datastring, s+skip, strlen(s)-skip);
      free (s);
      return strlen(*datastring);
    } else {
      if (irrelevant != NULL) {
	*irrelevant = strdup (s);
      }
      *datastring = NULL;
      free (s);
      return 0;
    }
    break;
  case SKIP_FIELD:
    if (skip == 0){
      *datastring = strdup (s);
      free (s);
      return strlen(*datastring);
    }
    
    cols_count=-1; inwhite=TRUE;
    skip_size=0;
    while (cols_count < skip && ssize > skip_size) {
      delim_found = FALSE;
      
      if (strchr(column_delimiters, s[skip_size]) != NULL) delim_found = TRUE;
      if (s[skip_size] == '\n' || s[skip_size] == '\r') delim_found = TRUE;
      if (inwhite != delim_found) {
	inwhite = delim_found;
	if (delim_found == FALSE) {
	  cols_count++;
	}
      }
      skip_size++;
    }
    /* at hits point, skip_size should be index of first valid char */
    if (strlen(s) > skip_size) {
      skip_size--;
      if (irrelevant != NULL) {
	*irrelevant = (char *) malloc ((skip_size+1)*sizeof(char));
	memset (*irrelevant, 0, skip_size+1);
	strncat (*irrelevant, s, skip_size);
      }
      *datastring = (char *) malloc ((strlen(s)-skip_size+1)*sizeof(char));
      memset (*datastring, 0, strlen(s)-skip_size+1);
      strncat (*datastring, s+skip_size, strlen(s)-skip_size);
      free (s);
      return strlen(*datastring);
    } else {
      if (irrelevant != NULL) {
	*irrelevant = strdup (s);
      }
      *datastring = NULL;
      free (s);
      return 0;
    }
    break;
  default:
#ifdef WITH_ERRORS
	  puterr("parse","parse_line");
	  fprintf (stderr, "invalid value\n");
#endif
    break;
  }
  return -2;
}



/*
 *% .TH GET_FILEFORMAT 3ggl "16 September 1997" "GGL/IGP ETHZ" "GGL Common Lib"
 *% .SH NAME
 *% get_fileformat \- check file format and give clue for template
 *% .SH SYNOPSIS
 *% .B #include <common/inout.h>
 *% 
 *% .RB "char * " get_fileformat
 *% .RI "(FILE *" fpin ", int " type ", const char *" comment_delimiters ", const char *" column_delimiters ");"
 *% .SH DESCRIPTION
 *% This function reads the first line or segment of a file
 *% .I fpin
 *% and attempts to find out, which template string has to be passed
 *% to a scanf or an analogous reading function to input the file
 *% properly.
 *% .B get_fileformat()
 *% support up to now only the ASCII file format, but other formats,
 *% like binary or ASCII block oriented format are planned to be
 *% supperted in later releases.
 *% 
 *% The
 *% .I comment_delimiters
 *% are ignored at the moment, so when the beginning of the file
 *% consists of comment lines,
 *% .B get_fileformat()
 *% will return a template to read strings.
 *% .SH "RETURN VALUES"
 *% The template returned will either contain a %li
 *% for a long integer, an octal number or a hex number
 *% which had been written in the 0x00 format,
 *% a %lf for float or double numbers,
 *% a %x for a string that only contains hexadecimal digits,
 *% or a %[^ \\n\\r\\t\\v\\f] for ordinary strings.
 *% This latter exclusion string is planned to be built
 *% using the
 *% .I column_delimiters
 *% passed in the function call.
 *% .SH CAVEATS
 *% You have to compile your program using
 *% .I -lcommon
 *% as a compiler flag to properly link
 *% your function to the library.
 *% .SH AUTHORS
 *%   (c) 1997 Geodesy and Geodynamics Lab, ETH Zurich:
 *% .RS
 *% Etienne Favey <favey@geod.ethz.ch>,
 *% .RE
 *% .SH "SEE ALSO"
 *% .BR parse_line (3ggl),
 *% .BR read_line (3ggl),
 *% .BR isdigit (3),
 *% .BR isxdigit (3),
 *% .BR fscanf (3),
 *% .BR strtok (3)
 */

char * get_fileformat (
  FILE *fpin,
  const int type,
  const char *comment_delimiters,
  const char *column_delimiters)
/*
 *  DEBUG: extend funtionality to dealing with comments and empty lines
 */
{
  char *s, *stmp;
  char *tmpstr = NULL;
  int length;
  char *token;
  int l_tok;
  char *format;
  int format_pos; /* only for debugging purposes */
  int err ,i;
  int type_found;
  /* we take white space as delimiters and currently ignore the passed arg */
  const char delimiters[] = " \f\n\r\t\v";
  char *output_delim;

  output_delim = (char *) malloc (strlen(delimiters)+4*sizeof(char));
  sprintf (output_delim, " %%[^%s]", delimiters);

  switch (type) {
  case ASCII_LINE:
    /*
     *  process one line of ASCII file, try to find out the data types
     */
    err = read_line (fpin, &s);
    if (err == EOF) {
#ifdef WITH_ERRORS
      fprintf (stderr, "GGL Common Library -- parse\n");
      fprintf (stderr,
	       "get_fileformat: unexpected end of file, check file\n");
#endif
      free (s);
      free(output_delim);
      return NULL;
    }

    /*
     *  check whether it really _is_ ASCII
     */
    length = strlen (s);
    if (length == 0) {
#ifdef WITH_ERRORS
      fprintf (stderr, "GGL Common Library -- parse\n");
      fprintf (stderr,
	       "get_fileformat: empty line read, empty file?\n");
#endif
      free (s);
      free(output_delim);
      return NULL;
    }


    for (i=0;i<length-1;i++) if (!isprint(s[i])) {
#ifdef WITH_ERRORS
      fprintf (stderr, "GGL Common Library -- parse\n");
      fprintf (stderr,
	       "get_fileformat: binary type, ascii was desired\n");
#endif
      free (s);
      free(output_delim);
      return NULL;
    }

    /*
     *  extract blocks delimited by white space, count blocks
     *  DEBUG: allow blocks delimited by arbitrary character?
     *  DEBUG: allow array of delimiters to be exchanged?
     */
    format = (char *) malloc (length * (strlen(delimiters)+3) * sizeof(char));
    memset (format, 0, sizeof(format));
    format_pos = 0;
    stmp = strdup (s);
    token = strtok (stmp, delimiters);
    while (token != NULL) {

      type_found = FALSE;
      l_tok = strlen (token);
#ifdef LIB_DEBUG
      printf ("token: %s | length %d\n", token, l_tok);
#endif
      /*
       *  is current block long integer (only digits or sign)
       */
      if (FALSE == type_found) type_found = 2;
      for (i=0;(i<l_tok && (2 == type_found || 6 == type_found));i++) {
	switch (i) {
	case 0:
	  if ( (isdigit(token[i]) || (token[i]=='+') ||
		(token[i]=='-')) ) {
	    continue;
	  }
	  break;
	case 1:
	  if (token[0]=='0' && tolower(token[1])=='x') {
	    /* may be hex */
	    type_found = 6;
	    continue;
	  }
	  if ( isdigit(token[i]) ) {
	    continue;
	  }
	  break;
	default:
	  if ( 2 == type_found && isdigit(token[i]) ) {
	    continue;
	  }
	  if ( 6 == type_found && isxdigit(token[i]) ) {
	    continue;
	  }
	  break;
	}
	/* not an integer */
	type_found = FALSE;
      }
#ifdef LIB_DEBUG
      printf ("int i: %d | type_found %d\n", i, type_found);
#endif
      if (2 == type_found || 6 == type_found) { 
	/* integer or octal or hex found */
	if (length * (strlen(delimiters)+3) * sizeof(char) 
	    - format_pos - 4 < 0)
	  {
	    /* this traps my programming mistakes, check malloc statement */
	    fprintf (stderr, "GGL Common Library -- parse\n");
	    fprintf (stderr,
		     "get_fileformat: DEBUG: format size too small\n");
	    fflush (stderr);
	    abort();
	  }
	strcat (format, " %li");
	format_pos = strlen(format);
/* 	type_found = TRUE; */
      } 
      /*
       *  is current block double (sign, digits, point, letter e)
       */
      if (FALSE == type_found) type_found = 3;
      for (i=0;(i<l_tok && 3 == type_found);i++) {
	switch (i) {
	case 0:
	  if ( isdigit(token[i]) || (token[i]=='+') ||
	       (token[i]=='-') || (token[i]=='.') ) {
	    continue;
	  }
	  break;
	default:
	  if (isdigit(token[i]) || (token[i]=='.') ||
		       (tolower(token[i])=='e') ) {
	    continue;
	  }
	  break;
	}
	/* not a double */
	type_found = FALSE;
      }
#ifdef LIB_DEBUG
      printf ("double i: %d | type_found %d\n", i, type_found);
#endif
      if (3 == type_found) { 
	/* float or double found */
	if (length * (strlen(delimiters)+3) * sizeof(char)
	    - format_pos - 4 < 0)
	  {
	    /* this traps my programming mistakes, check malloc statement */
	    fprintf (stderr, "GGL Common Library -- parse\n");
	    fprintf (stderr,
		     "get_fileformat: DEBUG: format size too small\n");
	    fflush (stderr);
	    abort();
	  }
	strcat (format, " %lf");
	format_pos = strlen(format);
/* 	type_found = TRUE; */
      }
      /*
       *  is current block hex (only xdigits)
       */
      if (FALSE == type_found) type_found = 4;
      for (i=0;(i<l_tok && 4 == type_found);i++) {
	if ( isxdigit(token[i]) ) {
	  continue;
	}
	/* not a hex */
	type_found = FALSE;
      }
#ifdef LIB_DEBUG
      printf ("hex i: %d | type_found %d\n", i, type_found);
#endif
      if (4 == type_found) { 
	if (length * (strlen(delimiters)+3) * sizeof(char)
	    - format_pos - 3 < 0)
	  {
	    /* this traps my programming mistakes, check malloc statement */
	    fprintf (stderr, "GGL Common Library -- parse\n");
	    fprintf (stderr,
		     "get_fileformat: DEBUG: format size too small\n");
	    fflush (stderr);
	    abort();
	  }
	strcat (format, " %x");
	format_pos = strlen(format);
/* 	type_found = TRUE; */
      }
      /*
       *  is current block our common lib angular format
       *  (like double, but with letters d or o @ 2nd or 3rd spot,
       *   and m or ')
       */
      if (type_found == FALSE) type_found = 5;
      for (i=0;(i<l_tok && type_found == 5);i++) { 
	switch (i) {
	case 0:
	  if ((token[i]=='-') ||
	      (token[i]=='+') || isdigit(token[i])
	      ) continue; /* for loop */
	  break; /* switch statement */
	case 1:
	  /* fall through */
	case 2:
	  if (
	      ( (token[i]=='o') || (token[i]=='d')  ) ||
	      isdigit(token[i])
	      ) continue; /* for loop */
	  break; /* switch statement */
	case 3:
	  /* fall through */
	case 4:
	  if (
	      ( (token[i]=='o') || (token[i]=='d')  ) ||
	      ( (token[i]=='m') || (token[i]=='\'') ) ||
	      isdigit(token[i])
	      ) continue; /* for loop */
	  break; /* switch statement */
	case 5:
	  /* fall through */
	case 6:
	  /* fall through */
	case 7:
	  if (
	      ( (token[i]=='m') || (token[i]=='\'') ) ||
	      (token[i]=='.') || 
	      isdigit(token[i])
	      ) continue; /* for loop */
	  break; /* switch statement */
	default:
	  if (
	      (token[i]=='.') || 
	      isdigit(token[i])
	      ) continue; /* for loop */
	  break; /* switch statement */
	}
	/* not an angle so we leave it a string */
	type_found = FALSE;
      }
      if (5 == type_found) {
	/* this is our common lib angle format */
	if (length * (strlen(delimiters)+3) * sizeof(char) - format_pos - 3
	    - strlen(delimiters) < 0)
	  {
	    /* this traps my programming mistakes, check malloc statement */
	    fprintf (stderr, "GGL Common Library -- parse\n");
	    fprintf (stderr,
		     "get_fileformat: DEBUG: format size too small\n");
	    fflush (stderr);
	    abort();
	  }
	strcat (format, output_delim);
	format_pos = strlen(format);
/* 	type_found = TRUE; */
#ifdef LIB_DEBUG
	fprintf (stderr, "GGL Common Library -- parse\n");
	fprintf (stderr,
		 "get_fileformat: recognized angle format: %s\n",token);
#endif
      }

      if (type_found == FALSE) {
	if (length * (strlen(delimiters)+3) * sizeof(char) - format_pos - 3
	    - strlen(delimiters) < 0)
	  {
	    /* this traps my programming mistakes, check malloc statement */
	    fprintf (stderr, "GGL Common Library -- parse\n");
	    fprintf (stderr,
		     "get_fileformat: DEBUG: format size too small\n");
	    fflush (stderr);
	    abort();
	  }
	strcat (format, output_delim);
	format_pos = strlen(format);
	type_found = TRUE;
      }
      /*
       *  read next token
       */
      token = strtok (NULL, delimiters);
    }
    tmpstr = strdup (format);
	free(format);
    free(output_delim);
    free(stmp);
    free (s);
	return tmpstr;
  case ASCII_BLOCK:
    /* fall through */
  case BINARY_BLOCK:
    /* fall through */
  default:
    /* error -- type not implemented */
#ifdef WITH_ERRORS
    fprintf (stderr, "GGL Common Library -- parse\n");
    fprintf (stderr,
	     "get_fileformat: desired file type not implemented\n");
#endif
    free(output_delim);
    return NULL;
  }
}
