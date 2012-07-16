/* gtools.c : Common routines for gtools programs. */
/* Version 1.4, August 2002. */

/* Todo: size check if MAXN>0
         sparse graph format */

#include "gtools.h"

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

long ogf_linelen;
boolean is_pipe;

#if !FTELL_DEC
extern long ftell(FILE*);
extern int fseek(FILE*,long,int);
#endif
#if !POPEN_DEC
extern FILE *popen(const char*,const char*);
#endif

/*
  Version 1.1: Fixed sparse6 input for powers of 2.  May 9, 1998
  Version 1.2: Added "cmd: ..." option for opengraphfile(). 
	       Fixed readg() bug (could not be invisible).  Oct 5, 1998
  Version 1.3: Added "is_pipe".  June 20, 2002
  Version 1.4: Stuff for autoconf.  August 30, 2002
*/

/*********************************************************************
opengraphfile(filename,codetype,assumefixed,position) 
          opens and positions a file for reading graphs.

  filename = the name of the file to open 
		(NULL means stdin, assumed already open)
             If filename starts with "cmd:", the remainder is taken
             to be a command to open a subshell for, using a pipe.
  codetype   = returns a code for the format.
		This is a combination of SPARSE6, GRAPH6,
		UNKNOWN_TYPE and HAS_HEADER.  If a header is
		present, that overrides the data.  If there is
                no header, the first graph is examined.
  assumefixed = nonzero if files other than stdin or pipes should be
		assumed to be seekable and have equal record sizes.
		Ignored if there is a sparse6 header or the first
		graph has sparse6 format.
  position = the number of the record to position to
		(the first is number 1; 0 and -NOLIMIT also mean
                 to position at start)

  If the file starts with ">", there must be a header, either
  GRAPH6_HEAD or SPARSE6_HEAD.  Otherwise opengraphfile() fails.

  The value returned is a file pointer or NULL.  
  If assumedfixed is not zero and position > 1, the global variable
  ogf_linelen is set to the length (including \n) of the length of the 
  first record.

  The global variable is_pipe is set to whether the input file is a pipe.

**********************************************************************/

FILE*
opengraphfile(char *filename, int *codetype, int assumefixed, long position)
{
	FILE *f;
	int c,firstc;
	long i,l,pos,pos1,pos2;
	boolean bad_header;

	is_pipe = FALSE;

	if (filename == NULL)
	{
	    f = stdin;
	    assumefixed = FALSE;
	}
	else
	{
	    if (filename[0] == 'c' && filename[1] == 'm'
		&& filename[2] == 'd' && filename[3] == ':')
	    {
#if !HAVE_POPEN
		gt_abort
		   (">E The "cmd:" option is not available in this version");
#else
		filename += 4;
		while (*filename == ' ') ++filename;
		f = popen(filename,"r");
#endif
		assumefixed = FALSE;
		is_pipe = TRUE;
	    }
	    else
	        f = fopen(filename,"r");

	    if (f == NULL)
	    {
	 	fprintf(stderr,">E opengraphfile: can't open %s\n",filename);
		return NULL;
	    }
	}

	firstc = c = getc(f);
	if (c == EOF)
	{
	    *codetype = GRAPH6;
	    return f;
	}

	if (c != '>')
	{
	    *codetype = firstc == ':' ? SPARSE6 : GRAPH6;
	    ungetc(c,f);
	}
	else
	{
	    bad_header = FALSE;
	    if ((c = getc(f)) == EOF || c != '>')
	        bad_header = TRUE;
	    if (!bad_header && ((c = getc(f)) == EOF || c != 'g' && c != 's'))
		bad_header = TRUE;	
	    if (!bad_header && c == 'g')
	    {
		if ((c = getc(f)) == EOF || c != 'r' ||
		    (c = getc(f)) == EOF || c != 'a' ||
		    (c = getc(f)) == EOF || c != 'p' ||
		    (c = getc(f)) == EOF || c != 'h' ||
		    (c = getc(f)) == EOF || c != '6' ||
		    (c = getc(f)) == EOF || c != '<' ||
		    (c = getc(f)) == EOF || c != '<')
			bad_header = TRUE;
		else
		    *codetype = GRAPH6 | HAS_HEADER;
	    }
	    else if (!bad_header && c == 's')
	    {
		if ((c = getc(f)) == EOF || c != 'p' ||
		    (c = getc(f)) == EOF || c != 'a' ||
		    (c = getc(f)) == EOF || c != 'r' ||
		    (c = getc(f)) == EOF || c != 's' ||
		    (c = getc(f)) == EOF || c != 'e' ||
		    (c = getc(f)) == EOF || c != '6' ||
		    (c = getc(f)) == EOF || c != '<' ||
		    (c = getc(f)) == EOF || c != '<')
			bad_header = TRUE;
		else
		    *codetype = SPARSE6 | HAS_HEADER;
	    }
	    if (bad_header)
	    {
		fprintf(stderr,">E opengraphfile: illegal header in %s\n",
			filename == NULL ? "stdin" : filename);
		*codetype = UNKNOWN_TYPE | HAS_HEADER;
		return NULL;
	    }
	}

	if (position <= 1) return f;

	if (!assumefixed || (*codetype&SPARSE6) || firstc == ':')
	{
	    l = 1;
	    while ((c = getc(f)) != EOF)
	    {
	        if (c == '\n')
		{
		    ++l;
		    if (l == position) break;
		}
	    }
	    if (l == position) return f;

	    fprintf(stderr,
               ">E opengraphfile: can't find line %ld in %s\n",position,
		filename == NULL ? "stdin" : filename);
	    return NULL;
	}
	else
	{
	    pos1 = ftell(f);
	    if (pos1 < 0)
	    {
		fprintf(stderr,">E opengraphfile: error on first ftell\n");
		return NULL;
	    }

	    for (i = 1; (c = getc(f)) != EOF && c != '\n'; ++i) {}
	    ogf_linelen = i;

	    if (c == EOF)
	    {
		fprintf(stderr,
		        ">E opengraphfile: required record no present\n");
		return NULL;
	    }
	    
	    pos2 = ftell(f);
	    if (pos2 < 0)
            {
                fprintf(stderr,">E opengraphfile: error on second ftell\n");
                return NULL;
            }

	    pos = pos1 + (position-1)*(pos2-pos1);
	    if (fseek(f,pos,SEEK_SET) < 0)
	    {
		fprintf(stderr,">E opengraphfile: seek failed\n");
		return NULL;
	    }
	}

	return f;
}

/*********************************************************************/

void
writeline(FILE *f, char *s)
/* write a line with error checking */
/* \n is not appended automatically */
{
	if (fputs(s,f) == EOF || ferror(f))
	{
	    fprintf(stderr,">E writeline : error on writing\n");
	    ABORT(">E writeline");
	}
}

/*********************************************************************/

char*
getline(FILE *f)     /* read a line with error checking */
/* includes \n (if present) and \0.  Immediate EOF causes NULL return. */
{
	DYNALLSTAT(char,s,s_sz);
	int c;
	long i;

	DYNALLOC1(char,s,s_sz,500,"getline");

	i = 0;
	while ((c = getc(f)) != EOF && c != '\n')
	{
	    if (i == s_sz-2) DYNREALLOC(char,s,s_sz,s_sz+1000,"getline");
	    s[i++] = c;
	}

	if (i == 0 && c == EOF) return NULL;

	if (c == '\n') s[i++] = '\n';
	s[i] = '\0';
	return s;
}

/****************************************************************************/

int
graphsize(char *s)
/* Get size of graph out of graph6 or sparse6 string. */
{
	char *p;
	int n;

	if (s[0] == ':') p = s+1;
	else             p = s;
	n = *p++ - BIAS6;

	if (n > SMALLN) 
	{
	    n = *p++ - BIAS6;
	    n = (n << 6) | (*p++ - BIAS6);
	    n = (n << 6) | (*p++ - BIAS6);
        }
	return n;
}

/****************************************************************************/

void
stringtograph(char *s, graph *g, int m)
/* Convert string (graph6 or sparse6 format) to graph. */
/* Assumes g is big enough to hold it.                 */
{
	char *p;
	int n,i,j,k,v,x,nb;
	long ii;
	set *gi,*gj;

	n = graphsize(s);

	p = s + 1 + (s[0] == ':');
	if (n > SMALLN) p += 3;

	if (TIMESWORDSIZE(m) < n)
	    gt_abort(">E stringtograph: impossible m value\n");

	for (ii = (long)m*n; --ii >= 0;)
	    g[ii] = 0;

	if (s[0] != ':')       /* graph6 format */
	{
	    k = 1;
	    for (j = 1; j < n; ++j)
	    {
	        gj = GRAPHROW(g,j,m);
    
	        for (i = 0; i < j; ++i)
	        {
	            if (--k == 0)
	            {
		        k = 6;
		        x = *(p++) - BIAS6;
	            }
	    
	            if (x & TOPBIT6)
	            {
		        gi = GRAPHROW(g,i,m);
		        ADDELEMENT(gi,j);
		        ADDELEMENT(gj,i);
	            }
	            x <<= 1;
	        }
	    }
	}
	else    /* sparse6 format */
	{
	    for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb)
	    {}

	    k = 1;
	    v = 0;
	    for (;;)
	    {
		if (--k == 0)
		{
		    k = 6;
		    if (*p == '\n' || *p == '\0') break;
		    else x = *p - BIAS6;
		    ++p;
		}
		else
		    x <<= 1;

		if (x & TOPBIT6) ++v;
		j = 0;
		for (i = 0; i < nb; ++i)
		{
		    if (--k == 0)
		    {
		 	k = 6;
			if (*p == '\n' || *p == '\0') break;
			else x = *p - BIAS6;
			++p;
		    }
		    else
			x <<= 1;
		    if (x & TOPBIT6) j = (j << 1) | 1;
		    else             j <<= 1;
	 	}
		if (i < nb) break;
		if (j > v)
		    v = j;
		else if (v < n)
		{
		    ADDELEMENT(GRAPHROW(g,v,m),j);
		    ADDELEMENT(GRAPHROW(g,j,m),v);
		}
	    }
	}
}

/***********************************************************************/

graph*                 /* read graph into nauty format */
readg(FILE *f, graph *g, int reqm, int *pm, int *pn) 
/* graph6 and sparse6 formats are supported */
/* f = an open file */
/* g = place to put the answer (NULL for dynamic allocation) */
/* reqm = the requested value of m (0 => compute from n) */
/* *pm = the actual value of m */
/* *pn = the value of n */
{
	char *s,*p;
	int m,n;

	if ((readg_line = getline(f)) == NULL) return NULL;

	s = readg_line;
	if (s[0] == ':')
	{
	    readg_code = SPARSE6;
	    p = s + 1;
	}
	else
	{
	    readg_code = GRAPH6;
            p = s;
	}

	while (*p >= BIAS6 && *p <= MAXBYTE) 
	    ++p;
	if (*p == '\0')
	    gt_abort(">E readg: missing newline\n");
	else if (*p != '\n')
	    gt_abort(">E readg: illegal character\n");

	n = graphsize(s);
	if (readg_code == GRAPH6 && p - s != G6LEN(n))
	    gt_abort(">E readg: truncated graph6 line\n");

	if (reqm > 0 && TIMESWORDSIZE(reqm) < n)
	    gt_abort(">E readg: reqm too small\n");
	else if (reqm > 0)
	    m = reqm;
	else
	    m = (n + WORDSIZE - 1) / WORDSIZE;

	if (g == NULL)
	{
	    if ((g = (graph*)ALLOCS(n,m*sizeof(graph))) == NULL)
		gt_abort(">E readg: malloc failed\n");
	}

	*pn = n;
	*pm = m;

	stringtograph(s,g,m);
	return g;
}

/****************************************************************************/

DYNALLSTAT(char,gcode,gcode_sz);  /* Used by ntog6 and ntos6 */
int readg_code;
char *readg_line;

/****************************************************************************/

char*
ntog6(graph *g, int m, int n)
/* convert nauty graph to graph6 string, including \n and \0 */
{
        register int i,j,k;
        register char *p,x;
	register set *gj;
	long ii;

	ii = G6LEN(n)+3;

	DYNALLOC1(char,gcode,gcode_sz,ii,"ntog6");

	p = gcode;
	if (n <= SMALLN) 
            *p++ = BIAS6 + n;
	else
	{
	    *p++ = MAXBYTE;
	    *p++ = BIAS6 + (n >> 12);
	    *p++ = BIAS6 + ((n >> 6) & C6MASK);
	    *p++ = BIAS6 + (n & C6MASK);
	}

	k = 6;
	x = 0;

	for (j = 1; j < n; ++j)
	{
	    gj = GRAPHROW(g,j,m);
	    for (i = 0; i < j; ++i)
	    {
		x <<= 1;
		if (ISELEMENT(gj,i)) x |= 1;
		if (--k == 0)
		{
		    *p++ = BIAS6 + x;
		    k = 6;
		    x = 0;
		}
	    }
	}

	if (k != 6) *p++ = BIAS6 + (x << k);

	*p++ = '\n';
	*p = '\0';

	return gcode;
}

/****************************************************************************/

char*
ntos6(graph *g, int m, int n)
/* convert nauty graph to sparse6 string, including \n and \0 */
{
        register int i,j,k;
        register char *p,x;
	register set *gj;
	long ii;
	int r,rr,topbit,nb,lastj;
	char *plim;

	DYNALLOC1(char,gcode,gcode_sz,500,"ntos6");

	plim = gcode + gcode_sz - 20;

	gcode[0] = ':';
	p = gcode+1;
	if (n <= SMALLN) 
            *p++ = BIAS6 + n;
	else
	{
	    *p++ = MAXBYTE;
	    *p++ = BIAS6 + (n >> 12);
	    *p++ = BIAS6 + ((n >> 6) & C6MASK);
	    *p++ = BIAS6 + (n & C6MASK);
	}

	for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb)
	{}
	topbit = 1 << (nb-1);
	k = 6;
	x = 0;

	lastj = 0;
	for (j = 0; j < n; ++j)
	{
	    gj = GRAPHROW(g,j,m);
	    for (i = 0; i <= j; ++i)
	    {
		if (ISELEMENT(gj,i))
		{
	       	    if (p >= plim)
		    {
			ii = p - gcode;
			DYNREALLOC(char,gcode,gcode_sz,
				           gcode_sz+1000,"ntos6");
			p = gcode + ii;
			plim = gcode + gcode_sz - 20;
		    }
		    if (j == lastj)
		    {
		        x <<= 1;
		        if (--k == 0)
		        {
		            *p++ = BIAS6 + x;
		            k = 6;
		            x = 0;
		        }
		    }
		    else
		    {
			x = (x << 1) | 1;
			if (--k == 0)
			{
		            *p++ = BIAS6 + x;
			    k = 6;
			    x = 0;
			}
			if (j > lastj+1)
			{
			    for (r = 0, rr = j; r < nb; ++r, rr <<= 1)
			    {
			        if (rr & topbit) x = (x << 1) | 1;
			        else             x <<= 1;
			        if (--k == 0)
			        {
				    *p++ = BIAS6 + x;
				    k = 6;
				    x = 0;
			        }
			    }
			    x <<= 1;
			    if (--k == 0)
			    {
				*p++ = BIAS6 + x;
				k = 6;
				x = 0;
			    }
			}
			lastj = j;
		    }
		    for (r = 0, rr = i; r < nb; ++r, rr <<= 1)
		    {
			if (rr & topbit) x = (x << 1) | 1;
			else             x <<= 1;
			if (--k == 0)
			{
			    *p++ = BIAS6 + x;
			    k = 6;
			    x = 0;
			}
		    }
		}
	    }
	}

	if (k != 6) *p++ = BIAS6 + ((x << k) | ((1 << k) - 1));

	*p++ = '\n';
	*p = '\0';
	return gcode;
}

/**************************************************************************/

void
writeg6(FILE *f, graph *g, int m, int n)
/* write graph to file in graph6 format */
{
        writeline(f,ntog6(g,m,n));
}

/**************************************************************************/

void
writes6(FILE *f, graph *g, int m, int n)
/* write graph to file in sparse6 format */
{
        writeline(f,ntos6(g,m,n));
}

/**************************************************************************/

void
writelast(FILE *f)
/* write last graph read by readg() assuming no intervening getline() */
{
	writeline(f,readg_line);
}

/**************************************************************************/

int
longvalue(char **ps, long *l)
{
	boolean neg,pos;
	long sofar,last;
	char *s;

	s = *ps;
	pos = neg = FALSE;
	if (*s == '-')
	{
	    neg = TRUE;
	    ++s;
	}
	else if (*s == '+')
	{
	    pos = TRUE;
	    ++s;
	}

	if (*s < '0' || *s > '9') 
	{
	    *ps = s;
	    return (pos || neg) ? ARG_ILLEGAL : ARG_MISSING;
	}

	sofar = 0;

	for (; *s >= '0' && *s <= '9'; ++s)
	{
	    last = sofar;
	    sofar = sofar * 10 + (*s - '0');
	    if (sofar < last || sofar > MAXARG)
	    {
		*ps = s;
		return ARG_TOOBIG;
	    }
	}
	*ps = s;
	*l = neg ? -sofar : sofar;
	return ARG_OK;
}
	
/*************************************************************************/

void
arg_long(char **ps, long *val, char *id)
{
	int code;

	code = longvalue(ps,val);
	if (code == ARG_MISSING || code == ARG_ILLEGAL)
	{
	    fprintf(stderr,">E %s: missing argument value\n",id);
	    gt_abort(NULL);
	}
	else if (code == ARG_TOOBIG)
	{
	    fprintf(stderr,">E %s: argument value too large\n",id);
	    gt_abort(NULL);
	}
}

/*************************************************************************/

void
arg_int(char **ps, int *val, char *id)
{
	int code;
	long longval;

	code = longvalue(ps,&longval);
	*val = longval;
	if (code == ARG_MISSING || code == ARG_ILLEGAL)
	{
	    fprintf(stderr,">E %s: missing argument value\n",id);
	    gt_abort(NULL);
	}
	else if (code == ARG_TOOBIG || *val != longval)
	{
	    fprintf(stderr,">E %s: argument value too large\n",id);
	    gt_abort(NULL);
	}
}

/************************************************************************/

boolean
strhaschar(char *s, int c)
/* Check if s contains c.  Saves the bother of figuring out whether
  strchr() is available, or index() or whatever.  */
{
	int i;

	for (i = 0; s[i] != '\0'; ++i)
	    if (s[i] == c) return TRUE;

	return FALSE;
}

/************************************************************************/

void
arg_range(char **ps, char *sep, long *val1, long *val2, char *id)
{
	int code,i;
	char *s;

	s = *ps;
	code = longvalue(&s,val1);
	if (code != ARG_MISSING)
	{
	    if (code == ARG_ILLEGAL)
	    {
		fprintf(stderr,">E %s: bad range\n",id);
		gt_abort(NULL);
	    }
	    else if (code == ARG_TOOBIG)
	    {
		fprintf(stderr,">E %s: value too big\n",id);
		gt_abort(NULL);
	    }
	}
	else if (*s == '\0' || !strhaschar(sep,*s))
	{
	    fprintf(stderr,">E %s: missing value\n",id);
	    gt_abort(NULL);
	}
	else
	    *val1 = -NOLIMIT;

	if (*s != '\0' && strhaschar(sep,*s))
	{
	    ++s;
	    code = longvalue(&s,val2);
	    if (code == ARG_MISSING)
		*val2 = NOLIMIT;
	    else if (code == ARG_TOOBIG)
	    {
		fprintf(stderr,">E %s: value too big\n",id);
		gt_abort(NULL);
	    }
	    else if (code == ARG_ILLEGAL)
	    {
		fprintf(stderr,">E %s: illegal range\n",id);
		gt_abort(NULL);
	    }
	}
	else
	    *val2 = *val1;

	*ps = s;
}

/***********************************************************************/

void
writerange(FILE *f, int c, long lo, long hi)    /* Write a range. */
{
	if (c != '\0') fprintf(f,"%c",c);
	if (lo != -NOLIMIT) fprintf(f,"%ld",lo);
	if (lo != hi)
	{
	    fprintf(stderr,":");
	    if (hi != NOLIMIT) fprintf(f,"%ld",hi);
	}
}

/************************************************************************/

void
gt_abort(char *msg)     /* Write message and halt. */
{
	if (msg) fprintf(stderr,msg);
	ABORT(">E gtools");
}

/************************************************************************/

char*
stringcopy(char *s)   /* duplicate string */
{
	char *scopy;
	size_t i,len;

	for (len = 0; s[len] != '\0'; ++len)
	{}

	if ((scopy = (char*)ALLOCS(len+1,1)) == NULL)
	    gt_abort(">E stringcopy: malloc failed\n");

	for (i = 0; i <= len; ++i)
	    scopy[i] = s[i];

	return scopy;
}
