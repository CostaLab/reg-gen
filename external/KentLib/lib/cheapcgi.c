/* Routines for getting variables passed in from web page
 * forms via CGI.
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "hash.h"
#include "cheapcgi.h"
#include "portable.h"
#include "linefile.h"
#include "errabort.h"
#ifndef GBROWSE
#include "mime.h"
#endif /* GBROWSE */
#include <signal.h>

static char const rcsid[] = "$Id: cheapcgi.c,v 1.137 2010/06/11 22:34:53 tdreszer Exp $";

/* These three variables hold the parsed version of cgi variables. */
static char *inputString = NULL;
static unsigned long inputSize;
static struct hash *inputHash = NULL;
static struct cgiVar *inputList = NULL;

static boolean haveCookiesHash = FALSE;
static struct hash *cookieHash = NULL;
static struct cgiVar *cookieList = NULL;

/* should cheapcgi use temp files to store uploaded files */
static boolean doUseTempFile = FALSE;

void dumpCookieList()
/* Print out the cookie list. */
{
struct cgiVar *v;
for (v=cookieList; v != NULL; v = v->next)
    printf("%s=%s (%d)\n", v->name, v->val, v->saved);
}

void useTempFile()
/* tell cheapcgi to use temp files */
{
doUseTempFile = TRUE;
}

boolean cgiIsOnWeb()
/* Return TRUE if looks like we're being run as a CGI. */
{
return getenv("REQUEST_METHOD") != NULL;
}

char *cgiRequestMethod()
/* Return CGI REQUEST_METHOD (such as 'GET/POST/PUT/DELETE/HEAD') */
{
return getenv("REQUEST_METHOD");
}

char *cgiRequestUri()
/* Return CGI REQUEST_URI */
{
return getenv("REQUEST_URI");
}

char *cgiRequestContentLength()
/* Return HTTP REQUEST CONTENT_LENGTH if available*/
{
return getenv("CONTENT_LENGTH");
}

char *cgiScriptName()
/* Return name of script so libs can do context-sensitive stuff. */
{
return getenv("SCRIPT_NAME");
}

char *cgiServerName()
/* Return name of server, better to use cgiServerNamePort() for
   actual URL construction */
{
return getenv("SERVER_NAME");
}

char *cgiServerPort()
/* Return port number of server, default 80 if not found */
{
char *port = getenv("SERVER_PORT");
if (port)
    return port;
else
    return "80";
}

char *cgiServerNamePort()
/* Return name of server with port if different than 80 */
{
char *port = cgiServerPort();
char *namePort = cgiServerName();
struct dyString *result = newDyString(256);
if (namePort)
    {
    dyStringPrintf(result,"%s",namePort);
    if (differentString(port, "80"))
	dyStringPrintf(result,":%s",port);
    return dyStringCannibalize(&result);
    }
else
    return NULL;
}

char *cgiRemoteAddr()
/* Return IP address of client (or "unknown"). */
{
static char *dunno = "unknown";
char *remoteAddr = getenv("REMOTE_ADDR");
if (remoteAddr == NULL)
    remoteAddr = dunno;
return remoteAddr;
}

char *cgiUserAgent()
/* Return remote user agent (HTTP_USER_AGENT) or NULL if remote user agent is not known */
{
return getenv("HTTP_USER_AGENT");
}

enum browserType cgiClientBrowser(char **browserQualifier, enum osType *clientOs, char **clientOsQualifier)
/* Return client browser type determined from (HTTP_USER_AGENT)
   Optionally requuest the additional info about the client */
{
// WARNING: The specifics of the HTTP_USER_AGENT vary widely.
//          This has only been tested on a few cases.
static enum browserType clientBrowser = btUnknown;
static enum browserType clientOsType  = osUnknown;
static char *clientBrowserExtra       = NULL;
static char *clientOsExtra            = NULL;

if (clientBrowser == btUnknown)
    {
    char *userAgent = cgiUserAgent();
    if (userAgent != NULL)
        {
        //warn(userAgent);  // Use this to investigate other cases
        char *ptr=NULL;

        // Determine the browser
        if ((ptr = stringIn("Opera",userAgent)) != NULL) // Must be before IE
            {
            clientBrowser = btOpera;
            }
        else if ((ptr = stringIn("MSIE ",userAgent)) != NULL)
            {
            clientBrowser = btIE;
            ptr += strlen("MSIE ");
            clientBrowserExtra = cloneFirstWordByDelimiter(ptr,';');
            }
        else if ((ptr = stringIn("Firefox",userAgent)) != NULL)
            {
            clientBrowser = btFF;
            ptr += strlen("(Firefox/");
            clientBrowserExtra = cloneFirstWordByDelimiter(ptr,' ');
            }
        else if ((ptr = stringIn("Chrome",userAgent)) != NULL)  // Must be before Safari
            {
            clientBrowser = btChrome;
            ptr += strlen("Chrome/");
            clientBrowserExtra = cloneFirstWordByDelimiter(ptr,' ');
            }
        else if ((ptr = stringIn("Safari",userAgent)) != NULL)
            {
            clientBrowser = btSafari;
            ptr += strlen("Safari/");
            clientBrowserExtra = cloneFirstWordByDelimiter(ptr,' ');
            }
        else
            {
            clientBrowser = btOther;
            }

        // Determine the OS
        if ((ptr = stringIn("Windows",userAgent)) != NULL)
            {
            clientOsType = osWindows;
            ptr += strlen("Windows ");
            clientOsExtra = cloneFirstWordByDelimiter(ptr,';');
            }
        else if ((ptr = stringIn("Linux",userAgent)) != NULL)
            {
            clientOsType = osLinux;
            ptr += strlen("Linux ");
            clientOsExtra = cloneFirstWordByDelimiter(ptr,';');
            }
        else if ((ptr = stringIn("Mac ",userAgent)) != NULL)
            {
            clientOsType = osMac;
            ptr += strlen("Mac ");
            clientOsExtra = cloneFirstWordByDelimiter(ptr,';');
            }
        else
            {
            clientOsType = osOther;
            }
        }
    }
if (browserQualifier != NULL)
    {
    if (clientBrowserExtra != NULL)
        *browserQualifier = cloneString(clientBrowserExtra);
    else
        *browserQualifier = NULL;
    }
if (clientOs != NULL)
    *clientOs = clientOsType;
if (clientOsQualifier != NULL)
    {
    if (clientOsExtra != NULL)
        *clientOsQualifier = cloneString(clientOsExtra);
    else
        *clientOsQualifier = NULL;
    }

return clientBrowser;
}

char *_cgiRawInput()
/* For debugging get the unprocessed input. */
{
return inputString;
}

static void getQueryInputExt(boolean abortOnErr)
/* Get query string from environment if they've used GET method. */
{
inputString = getenv("QUERY_STRING");
if (inputString == NULL)
    {
    if (abortOnErr)
	errAbort("No QUERY_STRING in environment.");
    inputString = cloneString("");
    return;
    }
inputString = cloneString(inputString);
}

static void getQueryInput()
/* Get query string from environment if they've used GET method. */
{
getQueryInputExt(TRUE);
}

static void getPostInput()
/* Get input from file if they've used POST method.
 * Grab any GET QUERY_STRING input first. */
{
char *s;
long i;
int r;

getQueryInputExt(FALSE);
int getSize = strlen(inputString);

s = getenv("CONTENT_LENGTH");
if (s == NULL)
    errAbort("No CONTENT_LENGTH in environment.");
if (sscanf(s, "%lu", &inputSize) != 1)
    errAbort("CONTENT_LENGTH isn't a number.");
s = getenv("CONTENT_TYPE");
if (s != NULL && startsWith("multipart/form-data", s))
    {
    /* use MIME parse on input stream instead, can handle large uploads */
    /* inputString must not be NULL so it knows it was set */
    return;
    }
int len = getSize + inputSize;
if (getSize > 0)
    ++len;
char *temp = needMem((size_t)len+1);
for (i=0; i<inputSize; ++i)
    {
    r = getc(stdin);
    if (r == EOF)
	errAbort("Short POST input.");
    temp[i] = r;
    }
if (getSize > 0)
  temp[i++] = '&';
strncpy(temp+i, inputString, getSize);
temp[len] = 0;
freeMem(inputString);
inputString = temp;
}

#define memmem(hay, haySize, needle, needleSize) \
    memMatch(needle, needleSize, hay, haySize)

#ifndef GBROWSE
static void cgiParseMultipart(struct hash **retHash, struct cgiVar **retList)
/* process a multipart form */
{
char h[1024];  /* hold mime header line */
char *s = NULL, *ct = NULL;
struct dyString *dy = newDyString(256);
struct mimeBuf *mb = NULL;
struct mimePart *mp = NULL;
char **env = NULL;
struct hash *hash = newHash(6);
struct cgiVar *list = NULL, *el;
extern char **environ;


//debug
//fprintf(stderr,"GALT: top of cgiParseMultipart()\n");
//fflush(stderr);

/* find the CONTENT_ environment strings, use to make Alternate Header string for MIME */
for(env=environ; *env; env++)
    if (startsWith("CONTENT_",*env))
	{
	//debug
    	//fprintf(stderr,"%s\n",*env);  //debug
	safef(h,sizeof(h),"%s",*env);
	s = strchr(h,'_');    /* change env syntax to MIME style header, from _= to -: */
	if (!s)
	    errAbort("expecting '_' parsing env var %s for MIME alt header", *env);
	*s = '-';
	s = strchr(h,'=');
	if (!s)
	    errAbort("expecting '=' parsing env var %s for MIME alt header", *env);
	*s = ':';
	dyStringPrintf(dy,"%s\r\n",h);
	}
dyStringAppend(dy,"\r\n");  /* blank line at end means end of headers */

//debug
//fprintf(stderr,"Alternate Header Text:\n%s",dy->string);
//fflush(stderr);
mb = initMimeBuf(STDIN_FILENO);
//debug
//fprintf(stderr,"got past initMimeBuf(STDIN_FILENO)\n");
//fflush(stderr);
mp = parseMultiParts(mb, cloneString(dy->string)); /* The Alternate Header will get freed */
freeDyString(&dy);
if(!mp->multi) /* expecting multipart child parts */
    errAbort("Malformatted multipart-form.");

//debug
//fprintf(stderr,"GALT: Wow got past parse of MIME!\n");
//fflush(stderr);

ct = hashFindVal(mp->hdr,"content-type");
//debug
//fprintf(stderr,"GALT: main content-type: %s\n",ct);
//fflush(stderr);
if (!startsWith("multipart/form-data",ct))
    errAbort("main content-type expected starts with [multipart/form-data], found [%s]",ct);

for(mp=mp->multi;mp;mp=mp->next)
    {
    char *cd = NULL, *cdMain = NULL, *cdName = NULL, *cdFileName = NULL, *ct = NULL;
    cd = hashFindVal(mp->hdr,"content-disposition");
    ct = hashFindVal(mp->hdr,"content-type");
    //debug
    //fprintf(stderr,"GALT: content-disposition: %s\n",cd);
    //fprintf(stderr,"GALT: content-type: %s\n",ct);
    //fflush(stderr);
    cdMain=getMimeHeaderMainVal(cd);
    cdName=getMimeHeaderFieldVal(cd,"name");
    cdFileName=getMimeHeaderFieldVal(cd,"filename");
    //debug
    //fprintf(stderr,"cgiParseMultipart: main:[%s], name:[%s], filename:[%s]\n",cdMain,cdName,cdFileName);
    //fflush(stderr);
    if (!sameString(cdMain,"form-data"))
	errAbort("main content-type expected [form-data], found [%s]",cdMain);

    //debug
    //fprintf(stderr,"GALT: mp->size[%llu], mp->binary=[%d], mp->fileName=[%s], mp=>data:[%s]\n",
	//(unsigned long long) mp->size, mp->binary, mp->fileName,
	//mp->binary && mp->data ? "<binary data not safe to print>" : mp->data);
    //fflush(stderr);

    /* filename if there is one */
    /* Internet Explorer on Windows is sending full path names, strip
     * directory name from those.  Using \ and / and : as potential
     * path separator characters, e.g.:
     *	 C:\Documents and Settings\tmp\file.txt.gz
     */
    if(cdFileName)
	{
	char *lastPathSep = strrchr(cdFileName, (int) '\\');
	if (!lastPathSep)
		lastPathSep = strrchr(cdFileName, (int) '/');
	if (!lastPathSep)
		lastPathSep = strrchr(cdFileName, (int) ':');
	char varNameFilename[256];
	safef(varNameFilename, sizeof(varNameFilename), "%s__filename", cdName);
	AllocVar(el);
	if (lastPathSep)
	    el->val = cloneString(lastPathSep+1);
	else
	    el->val = cloneString(cdFileName);
	slAddHead(&list, el);
	hashAddSaveName(hash, varNameFilename, el, &el->name);
    	}

    if (mp->data)
	{
	if (mp->binary)
	    {
	    char varNameBinary[256];
	    char addrSizeBuf[40];
	    safef(varNameBinary,sizeof(varNameBinary),"%s__binary",cdName);
	    safef(addrSizeBuf,sizeof(addrSizeBuf),"%lu %llu",
		(unsigned long)mp->data,
		(unsigned long long)mp->size);
	    AllocVar(el);
	    el->val = cloneString(addrSizeBuf);
	    slAddHead(&list, el);
	    hashAddSaveName(hash, varNameBinary, el, &el->name);
	    }
	else  /* normal variable, not too big, does not contain zeros */
	    {
	    AllocVar(el);
	    el->val = mp->data;
	    slAddHead(&list, el);
	    hashAddSaveName(hash, cdName, el, &el->name);
	    }
	}
    else if (mp->fileName)
	{
	char varNameData[256];
	safef(varNameData, sizeof(varNameData), "%s__data", cdName);
	AllocVar(el);
	el->val = mp->fileName;
	slAddHead(&list, el);
	hashAddSaveName(hash, varNameData, el, &el->name);
	//debug
    	//fprintf(stderr,"GALT special: saved varNameData:[%s], mp=>fileName:[%s]\n",el->name,el->val);
       	//fflush(stderr);
	}
    else if (mp->multi)
	{
	warn("unexpected nested MIME structures");
	}
    else
	{
	errAbort("mp-> type not data,fileName, or multi - unexpected MIME structure");
	}

    freez(&cdMain);
    freez(&cdName);
    freez(&cdFileName);
    }

slReverse(&list);
*retList = list;
*retHash = hash;
}
#endif /* GBROWSE */



static void parseCookies(struct hash **retHash, struct cgiVar **retList)
/* parses any cookies and puts them into the given hash and list */
{
char* str;
char *namePt, *dataPt, *nextNamePt;
struct hash *hash;
struct cgiVar *list = NULL, *el;

/* don't build the hash table again */
if(haveCookiesHash == TRUE)
	return;

str = cloneString(getenv("HTTP_COOKIE"));
if(str == NULL) /* don't have a cookie */
	return;

hash = newHash(6);

namePt = str;
while (isNotEmpty(namePt))
    {
    dataPt = strchr(namePt, '=');
    if (dataPt == NULL)
	errAbort("Mangled Cookie input string: no = in '%s' (offset %d in complete cookie string: '%s')",
		 namePt, (int)(namePt - str), getenv("HTTP_COOKIE"));
    *dataPt++ = 0;
    nextNamePt = strchr(dataPt, ';');
    if (nextNamePt != NULL)
	{
         *nextNamePt++ = 0;
	 if (*nextNamePt == ' ')
	     nextNamePt++;
	}
    cgiDecode(dataPt,dataPt,strlen(dataPt));
    AllocVar(el);
    el->val = dataPt;
    slAddHead(&list, el);
    hashAddSaveName(hash, namePt, el, &el->name);
    namePt = nextNamePt;
    }

haveCookiesHash = TRUE;

slReverse(&list);
*retList = list;
*retHash = hash;
}

char *findCookieData(char *varName)
/* Get the string associated with varName from the cookie string. */
{
struct hashEl *hel;
char *firstResult;

/* make sure that the cookie hash table has been created */
parseCookies(&cookieHash, &cookieList);
if (cookieHash == NULL)
    return NULL;
/* Watch out for multiple cookies with the same name (hel is a list) --
 * warn if we find them. */
hel = hashLookup(cookieHash, varName);
if (hel == NULL)
    return NULL;
else
    firstResult = ((struct cgiVar *)hel->val)->val;
hel = hel->next;
while (hel != NULL)
    {
    char *val = ((struct cgiVar *)(hel->val))->val;
    if (sameString(varName, hel->name) && !sameString(firstResult, val))
	{
	/* This is too early to call warn -- it will mess up html output. */
	fprintf(stderr,
		"findCookieData: Duplicate cookie value from IP=%s: "
		"%s has both %s and %s\n",
		cgiRemoteAddr(),
		varName, firstResult, val);
	}
    hel = hel->next;
    }
return firstResult;
}

static char *cgiInputSource(char *s)
/* For NULL sources make a guess as to real source. */
{
char *qs;
if (s != NULL)
    return s;
qs = getenv("QUERY_STRING");
if (qs == NULL)
    return "POST";
char *cl = getenv("CONTENT_LENGTH");
if (cl != NULL && atoi(cl) > 0)
    return "POST";
return "QUERY";
}

static void _cgiFindInput(char *method)
/* Get raw CGI input into inputString.  Method can be "POST", "QUERY", "GET" or NULL
 * for unknown. */
{
if (inputString == NULL)
    {
    method = cgiInputSource(method);
    if (sameWord(method, "POST"))
        getPostInput();
    else if (sameWord(method, "QUERY") || sameWord(method, "GET"))
        getQueryInput();
    else
        errAbort("Unknown form method");
    }
}

static void cgiParseInputAbort(char *input, struct hash **retHash,
	struct cgiVar **retList)
/* Parse cgi-style input into a hash table and list.  This will alter
 * the input data.  The hash table will contain references back
 * into input, so please don't free input until you're done with
 * the hash. Prints message aborts if there's an error.*/
{
char *namePt, *dataPt, *nextNamePt;
struct hash *hash = *retHash;
struct cgiVar *list = *retList, *el;

if (!hash)
  hash = newHash(6);
slReverse(&list);

namePt = input;
while (namePt != NULL && namePt[0] != 0)
    {
    dataPt = strchr(namePt, '=');
    if (dataPt == NULL)
	{
	errAbort("Mangled CGI input string %s", namePt);
	}
    *dataPt++ = 0;
    nextNamePt = strchr(dataPt, '&');
    if (nextNamePt == NULL)
	nextNamePt = strchr(dataPt, ';');	/* Accomodate DAS. */
    if (nextNamePt != NULL)
         *nextNamePt++ = 0;
    cgiDecode(namePt,namePt,strlen(namePt));	/* for unusual ct names */
    cgiDecode(dataPt,dataPt,strlen(dataPt));
    AllocVar(el);
    el->val = dataPt;
    slAddHead(&list, el);
    hashAddSaveName(hash, namePt, el, &el->name);
    namePt = nextNamePt;
    }
slReverse(&list);
*retList = list;
*retHash = hash;
}

static jmp_buf cgiParseRecover;

static void cgiParseAbort()
/* Abort cgi parsing. */
{
longjmp(cgiParseRecover, -1);
}

boolean cgiParseInput(char *input, struct hash **retHash,
	struct cgiVar **retList)
/* Parse cgi-style input into a hash table and list.  This will alter
 * the input data.  The hash table will contain references back
 * into input, so please don't free input until you're done with
 * the hash. Prints message and returns FALSE if there's an error.*/
{
boolean ok = TRUE;
int status = setjmp(cgiParseRecover);
if (status == 0)    /* Always true except after long jump. */
    {
    pushAbortHandler(cgiParseAbort);
    cgiParseInputAbort(input, retHash, retList);
    }
else    /* They long jumped here because of an error. */
    {
    ok = FALSE;
    }
popAbortHandler();
return ok;
}


static boolean dumpStackOnSignal = FALSE;  // should a stack dump be generated?

static void catchSignal(int sigNum)
/* handler for various terminal signals for logging purposes */
{
char *sig = NULL;
switch (sigNum)
    {
    case SIGABRT:
      sig = "SIGABRT";
      break;
    case SIGSEGV:
      sig = "SIGSEGV";
      break;
    case SIGFPE:
      sig = "SIGFPE";
      break;
    case SIGBUS:
      sig = "SIGBUS";
      break;
    }
    logCgiToStderr();
    fprintf(stderr, "Received signal %s\n", sig);
    if (dumpStackOnSignal)
        dumpStack("Stack for signal %s\n", sig);
    raise(SIGKILL);
}

void initSigHandlers(boolean dumpStack)
/* set handler for various terminal signals for logging purposes.
 * if dumpStack is TRUE, attempt to dump the stack. */
{
if (cgiIsOnWeb())
    {
    signal(SIGABRT, catchSignal);
    signal(SIGSEGV, catchSignal);
    signal(SIGFPE, catchSignal);
    signal(SIGBUS, catchSignal);
    dumpStackOnSignal = dumpStack;
    }
}


static void initCgiInput()
/* Initialize CGI input stuff.  After this CGI vars are
 * stored in an internal hash/list regardless of how they
 * were passed to the program. */
{
char* s;

if (inputString != NULL)
    return;

_cgiFindInput(NULL);

#ifndef GBROWSE
/* check to see if the input is a multipart form */
s = getenv("CONTENT_TYPE");
if (s != NULL && startsWith("multipart/form-data", s))
    {
    cgiParseMultipart(&inputHash, &inputList);
    }
#endif /* GBROWSE */

cgiParseInputAbort(inputString, &inputHash, &inputList);

/* now parse the cookies */
parseCookies(&cookieHash, &cookieList);

/* Set enviroment variables CGIs to enable sql tracing and/or profiling */
s = cgiOptionalString("JKSQL_TRACE");
if (s != NULL)
    envUpdate("JKSQL_TRACE", s);
s = cgiOptionalString("JKSQL_PROF");
if (s != NULL)
    envUpdate("JKSQL_PROF", s);

}

struct cgiVar *cgiVarList()
/* return the list of cgiVar's */
{
initCgiInput();
return inputList;
}

static char *findVarData(char *varName)
/* Get the string associated with varName from the query string. */
{
struct cgiVar *var;

initCgiInput();
if ((var = hashFindVal(inputHash, varName)) == NULL)
    return NULL;
return var->val;
}

void cgiBadVar(char *varName)
/* Complain about a variable that's not there. */
{
if (varName == NULL) varName = "";
errAbort("Sorry, didn't find input variable %s\n"
        "Probably the web page didn't mean to call this program.",
        varName);
}

static char *mustFindVarData(char *varName)
/* Find variable and associated data or die trying. */
{
char *res = findVarData(varName);
if (res == NULL)
    cgiBadVar(varName);
return res;
}

char *javaScriptLiteralEncode(char *inString)
/* Use backslash escaping on newline
 * and quote chars, backslash and others.
 * Intended that the encoded string will be
 * put between quotes at a higher level and
 * then interpreted by Javascript. */
{
char c;
int outSize = 0;
char *outString, *out, *in;

if (inString == NULL)
    return(cloneString(""));

/* Count up how long it will be */
in = inString;
while ((c = *in++) != 0)
    {
    if (c == '\''
     || c == '\"'
     || c == '&'
     || c == '\\'
     || c == '\n'
     || c == '\r'
     || c == '\t'
     || c == '\b'
     || c == '\f'
	)
        outSize += 2;
    else
        outSize += 1;
    }
outString = needMem(outSize+1);

/* Encode string */
in = inString;
out = outString;
while ((c = *in++) != 0)
    {
    if (c == '\''
     || c == '\"'
     || c == '&'
     || c == '\\'
     || c == '\n'
     || c == '\r'
     || c == '\t'
     || c == '\b'
     || c == '\f'
	)
        *out++ = '\\';
    *out++ = c;
    }
*out++ = 0;
return outString;

}


void cgiDecode(char *in, char *out, int inLength)
/* Decode from cgi pluses-for-spaces format to normal.
 * Out will be a little shorter than in typically, and
 * can be the same buffer. */
{
char c;
int i;
for (i=0; i<inLength;++i)
    {
    c = *in++;
    if (c == '+')
	*out++ = ' ';
    else if (c == '%')
	{
	int code;
	if (sscanf(in, "%2x", &code) != 1)
	    code = '?';
	in += 2;
	i += 2;
	*out++ = code;
	}
    else
	*out++ = c;
    }
*out++ = 0;
}

char *cgiEncode(char *inString)
/* Return a cgi-encoded version of inString.
 * Alphanumerics kept as is, space translated to plus,
 * and all other characters translated to %hexVal. */
{
char c;
int outSize = 0;
char *outString, *out, *in;

if (inString == NULL)
    return(cloneString(""));

/* Count up how long it will be */
in = inString;
while ((c = *in++) != 0)
    {
    if (isalnum(c) || c == ' ' || c == '.' || c == '_')
        outSize += 1;
    else
        outSize += 3;
    }
outString = needMem(outSize+1);

/* Encode string */
in = inString;
out = outString;
while ((c = *in++) != 0)
    {
    if (isalnum(c) || c == '.' || c == '_')
        *out++ = c;
    else if (c == ' ')
        *out++ = '+';
    else
        {
        unsigned char uc = c;
        char buf[4];
        *out++ = '%';
        safef(buf, sizeof(buf), "%02X", uc);
        *out++ = buf[0];
        *out++ = buf[1];
        }
    }
*out++ = 0;
return outString;
}

char *cgiEncodeFull(char *inString)
/* Return a cgi-encoded version of inString (no + for space!).
 * Alphanumerics/./_ kept as is and all other characters translated to
 * %hexVal. */
{
char c;
int outSize = 0;
char *outString, *out, *in;

if (inString == NULL)
    return(cloneString(""));

/* Count up how long it will be */
in = inString;
while ((c = *in++) != 0)
    {
    if (isalnum(c) || c == '.' || c == '_')
        outSize += 1;
    else
        outSize += 3;
    }
outString = needMem(outSize+1);

/* Encode string */
in = inString;
out = outString;
while ((c = *in++) != 0)
    {
    if (isalnum(c) || c == '.' || c == '_')
        *out++ = c;
    else
        {
        unsigned char uc = c;
        char buf[4];
        *out++ = '%';
        safef(buf, sizeof(buf), "%02X", uc);
        *out++ = buf[0];
        *out++ = buf[1];
        }
    }
*out++ = 0;
return outString;
}

char *cgiOptionalString(char *varName)
/* Return value of string if it exists in cgi environment, else NULL */
{
return findVarData(varName);
}


char *cgiString(char *varName)
/* Return string value of cgi variable. */
{
return mustFindVarData(varName);
}

char *cgiUsualString(char *varName, char *usual)
/* Return value of string if it exists in cgi environment.
 * Otherwise return 'usual' */
{
char *pt;
pt = findVarData(varName);
if (pt == NULL)
    pt = usual;
return pt;
}

struct slName *cgiStringList(char *varName)
/* Find list of cgi variables with given name.  This
 * may be empty.  Free result with slFreeList(). */
{
struct hashEl *hel;
struct slName *stringList = NULL, *string;

initCgiInput();
for (hel = hashLookup(inputHash, varName); hel != NULL; hel = hel->next)
    {
    if (sameString(hel->name, varName))
        {
	struct cgiVar *var = hel->val;
	string = newSlName(var->val);
	slAddHead(&stringList, string);
	}
    }
return stringList;
}


int cgiInt(char *varName)
/* Return int value of cgi variable. */
{
char *data;
char c;

data = mustFindVarData(varName);
data = skipLeadingSpaces(data);
c = data[0];
if (!(isdigit(c) || (c == '-' && isdigit(data[1]))))
     errAbort("Expecting number in %s, got \"%s\"\n", varName, data);
return atoi(data);
}

int cgiIntExp(char *varName)
/* Evaluate an integer expression in varName and
 * return value. */
{
return intExp(cgiString(varName));
}

int cgiOptionalInt(char *varName, int defaultVal)
/* This returns integer value of varName if it exists in cgi environment
 * and it's not just the empty string otherwise it returns defaultVal. */
{
char *s = cgiOptionalString(varName);
s = skipLeadingSpaces(s);
if (isEmpty(s))
    return defaultVal;
return cgiInt(varName);
}

double cgiDouble(char *varName)
/* Returns double value. */
{
char *data;
double x;

data = mustFindVarData(varName);
if (sscanf(data, "%lf", &x)<1)
     errAbort("Expecting real number in %s, got \"%s\"\n", varName, data);
return x;
}

double cgiOptionalDouble(char *varName, double defaultVal)
/* Returns double value. */
{
if (!cgiVarExists(varName))
    return defaultVal;
return cgiDouble(varName);
}

boolean cgiVarExists(char *varName)
/* Returns TRUE if the variable was passed in. */
{
initCgiInput();
return hashLookup(inputHash, varName) != NULL;
}

boolean cgiBoolean(char *varName)
{
return cgiVarExists(varName);
}

int cgiOneChoice(char *varName, struct cgiChoice *choices, int choiceSize)
/* Returns value associated with string variable in choice table. */
{
char *key = cgiString(varName);
int i;
int val = -1;

for (i=0; i<choiceSize; ++i)
    {
    if (sameWord(choices[i].name, key))
        {
        val =  choices[i].value;
        return val;
        }
    }
errAbort("Unknown key %s for variable %s\n", key, varName);
return val;
}

void cgiMakeSubmitButton()
/* Make 'submit' type button. */
{
cgiMakeButton("Submit", "Submit");
}

void cgiMakeResetButton()
/* Make 'reset' type button. */
{
printf("<INPUT TYPE=RESET NAME=\"Reset\" VALUE=\" Reset \">");
}

void cgiMakeClearButton(char *form, char *field)
/* Make button to clear a text field. */
{
char javascript[1024];

safef(javascript, sizeof(javascript),
    "document.%s.%s.value = ''; document.%s.submit();", form, field, form);
cgiMakeOnClickButton(javascript, " Clear  ");
}

void cgiMakeButtonWithMsg(char *name, char *value, char *msg)
/* Make 'submit' type button. Display msg on mouseover, if present*/
{
printf("<INPUT TYPE=SUBMIT NAME=\"%s\" VALUE=\"%s\" %s%s%s>",
        name, value,
        (msg ? " TITLE=\"" : ""), (msg ? msg : ""), (msg ? "\"" : "" ));
}

void cgiMakeButtonWithOnClick(char *name, char *value, char *msg, char *onClick)
/* Make 'submit' type button, with onclick javascript */
{
printf("<input type=\"submit\" name=\"%s\" value=\"%s\" onclick=\"%s\" %s%s%s>",
        name, value, onClick,
        (msg ? " TITLE=\"" : ""), (msg ? msg : ""), (msg ? "\"" : "" ));
}

void cgiMakeButton(char *name, char *value)
/* Make 'submit' type button */
{
cgiMakeButtonWithMsg(name, value, NULL);
}

void cgiMakeOnClickButton(char *command, char *value)
/* Make 'push' type button with client side onClick (java)script. */
{
printf("<INPUT TYPE=\"button\" VALUE=\"%s\" onClick=\"%s\">", value, command);
}

void cgiMakeOnClickSubmitButton(char *command, char *name, char *value)
/* Make submit button with both variable name and value with client side
 * onClick (java)script. */
{
printf("<INPUT TYPE=SUBMIT NAME=\"%s\" VALUE=\"%s\" onClick=\"%s\">",
       name, value, command);
}

void cgiMakeOptionalButton(char *name, char *value, boolean disabled)
/* Make 'submit' type button that can be disabled. */
{
printf("<INPUT TYPE=SUBMIT NAME=\"%s\" VALUE=\"%s\"", name, value);
if (disabled)
    printf(" DISABLED");
printf(">");
}

void cgiMakeFileEntry(char *name)
/* Make file entry box/browser */
{
    printf("<INPUT TYPE=FILE NAME=\"%s\">", name);
}

void cgiSimpleTableStart()
/* start HTML table  -- no customization. Leaves room
 * for a fancier implementation */
{
printf("<TABLE>\n");
}

void cgiTableEnd()
/* end HTML table */
{
printf("</TABLE>\n");
}

void cgiSimpleTableRowStart()
/* Start table row */
{
printf("<TR>\n");
}

void cgiTableRowEnd()
/* End table row */
{
printf("</TR>\n");
}

void cgiSimpleTableFieldStart()
/* Start table field */
{
printf("<TD>");
}

void cgiTableFieldStartAlignRight()
/* Start table field and align right*/
{
printf("<TD ALIGN = RIGHT>");
}

void cgiTableFieldEnd()
/* End table field */
{
printf("</TD>\n");
}

void cgiTableField(char *text)
/* Make table field entry */
{
printf("<TD> %s </TD>\n", text);
}

void cgiTableFieldWithMsg(char *text, char *msg)
/* Make table field entry with mouseover */
{
printf("<TD %s%s%s> %s </TD>\n",
        (msg ? " TITLE=\"" : ""), (msg ? msg : ""), (msg ? "\"" : "" ),
        text);
}

void cgiParagraph(char *text)
/* Make text paragraph */
{
printf("<P> %s\n", text);
}

void cgiMakeRadioButton(char *name, char *value, boolean checked)
/* Make radio type button.  A group of radio buttons should have the
 * same name but different values.   The default selection should be
 * sent with checked on. */
{
printf("<INPUT TYPE=RADIO NAME=\"%s\" VALUE=\"%s\" %s>", name, value,
   (checked ? "CHECKED" : ""));
}

void cgiMakeOnClickRadioButton(char *name, char *value, boolean checked,
                                        char *command)
/* Make radio type button with onClick command.
 *  A group of radio buttons should have the
 * same name but different values.   The default selection should be
 * sent with checked on. */
{
printf("<INPUT TYPE=RADIO NAME=\"%s\" VALUE=\"%s\" %s %s>",
        name, value, command, (checked ? "CHECKED" : ""));
}

char *cgiBooleanShadowPrefix()
/* Prefix for shadow variable set with boolean variables. */
{
return "boolshad.";
}

boolean cgiBooleanDefined(char *name)
/* Return TRUE if boolean variable is defined (by
 * checking for shadow. */
{
char buf[256];
safef(buf, sizeof(buf), "%s%s", cgiBooleanShadowPrefix(), name);
return cgiVarExists(buf);
}

static void cgiMakeCheckBox2Bool(char *name, boolean checked, boolean enabled, char *id, char *moreHtml)
/* Make check box - designed to be called by the variously overloaded
 * cgiMakeCheckBox functions, and should NOT be called directly.
 * moreHtml: optional additional html like javascript call or mouseover msg (may be NULL)
 * id: button id (may be NULL)
 * Also make a shadow hidden variable and support 2 boolean states:
 *    checked/unchecked and enabled/disabled. */
{
char buf[256], idBuf[256];

if(id)
    safef(idBuf, sizeof(idBuf), " id=\"%s\"", id);
else
    idBuf[0] = 0;

printf("<INPUT TYPE=CHECKBOX NAME=\"%s\"%s VALUE=on %s%s%s>", name, idBuf,
    (moreHtml ? moreHtml : ""),
    (checked ? " CHECKED" : ""),
    (enabled ? "" : " DISABLED"));
safef(buf, sizeof(buf), "%s%s", cgiBooleanShadowPrefix(), name);
cgiMakeHiddenVar(buf, ( enabled ? "0" : (checked ? "-1" : "-2")));
}

void cgiMakeCheckBoxUtil(char *name, boolean checked, char *msg, char *id)
/* Make check box - can be called directly, though it was originally meant
 * as the common code for all lower level checkbox routines.
 * However, it's util functionality has been taken over by
 * cgiMakeCheckBoxWithIdAndOptionalHtml() */
{
char buf[256];

if(msg)
    safef(buf, sizeof(buf), "TITLE=\"%s\"", msg);
else
    buf[0] = 0;

cgiMakeCheckBox2Bool(name, checked, TRUE, id, buf);
}

void cgiMakeCheckBoxWithMsg(char *name, boolean checked, char *msg)
{
cgiMakeCheckBox2Bool(name, checked, TRUE, NULL, msg);
}

void cgiMakeCheckBoxWithId(char *name, boolean checked, char *id)
/* Make check box, which includes an ID. */
{
cgiMakeCheckBox2Bool(name, checked, TRUE, id, NULL);
}

void cgiMakeCheckBox(char *name, boolean checked)
/* Make check box. */
{
cgiMakeCheckBox2Bool(name, checked, TRUE, NULL, NULL);
}

void cgiMakeCheckBoxJS(char *name, boolean checked, char *javascript)
/* Make check box with javascript. */
{
cgiMakeCheckBox2Bool(name,checked,TRUE,NULL,javascript);
}

void cgiMakeCheckBoxIdAndJS(char *name, boolean checked, char *id, char *javascript)
/* Make check box with ID and javascript. */
{
cgiMakeCheckBox2Bool(name,checked,TRUE,id,javascript);
}

void cgiMakeCheckBoxFourWay(char *name, boolean checked, boolean enabled, char *id, char *classes, char *moreHtml)
/* Make check box - with fourWay functionality (checked/unchecked by enabled/disabled)
 * Also makes a shadow hidden variable that supports the 2 boolean states. */
{
char shadName[256], extra[256];

printf("<INPUT TYPE=CHECKBOX NAME='%s'", name);
if(id)
    printf(" id='%s'", id);
if(checked)
    printf(" CHECKED");
if(!enabled)
    printf(" DISABLED");
if(classes)
    printf(" class='%s'",classes);
if(moreHtml)
    printf(" %s",moreHtml);
printf(">");

// The hidden var needs to hold the 4way state
safef(shadName, sizeof(shadName), "%s%s", cgiBooleanShadowPrefix(), name);
safef(extra, sizeof(extra), "id='%s_4way'",name);
cgiMakeHiddenVarWithExtra(shadName, ( enabled ? "0" : (checked ? "-1" : "-2")),extra); // Doesn't need enabled/checked!
}


void cgiMakeHiddenBoolean(char *name, boolean on)
/* Make hidden boolean variable. Also make a shadow hidden variable so we
 * can distinguish between variable not present and
 * variable set to false. */
{
char buf[256];
cgiMakeHiddenVar(name, on ? "on" : "off");
safef(buf, sizeof(buf), "%s%s", cgiBooleanShadowPrefix(), name);
cgiMakeHiddenVar(buf, "1");
}

void cgiMakeTextArea(char *varName, char *initialVal, int rowCount, int columnCount)
/* Make a text area with area rowCount X columnCount and with text: intialVal */
{
cgiMakeTextAreaDisableable(varName, initialVal, rowCount, columnCount, FALSE);
}

void cgiMakeTextAreaDisableable(char *varName, char *initialVal, int rowCount, int columnCount, boolean disabled)
/* Make a text area that can be disabled. The rea has rowCount X
 * columnCount and with text: intialVal */
{
printf("<TEXTAREA NAME=\"%s\" ROWS=%d COLS=%d %s>%s</TEXTAREA>", varName,
       rowCount, columnCount, disabled ? "DISABLED" : "",
       (initialVal != NULL ? initialVal : ""));
}

void cgiMakeOnKeypressTextVar(char *varName, char *initialVal, int charSize,
			      char *script)
/* Make a text control filled with initial value, with a (java)script
 * to execute every time a key is pressed.  If charSize is zero it's
 * calculated from initialVal size. */
{
if (initialVal == NULL)
    initialVal = "";
if (charSize == 0) charSize = strlen(initialVal);
if (charSize == 0) charSize = 8;

printf("<INPUT TYPE=TEXT NAME=\"%s\" SIZE=%d VALUE=\"%s\"", varName,
	charSize, initialVal);
if (isNotEmpty(script))
    printf(" onkeypress=\"%s\"", script);
printf(">\n");
}

void cgiMakeTextVar(char *varName, char *initialVal, int charSize)
/* Make a text control filled with initial value.  If charSize
 * is zero it's calculated from initialVal size. */
{
cgiMakeOnKeypressTextVar(varName, initialVal, charSize, NULL);
}

void cgiMakeTextVarWithExtraHtml(char *varName, char *initialVal, int width, char *extra)
/* Make a text control filled with initial value. */
{
if (initialVal == NULL)
    initialVal = "";
if(width==0)
    width=strlen(initialVal)*10;
if(width==0)
    width = 100;

printf("<INPUT TYPE=TEXT class='inputBox' NAME=\"%s\" style='width: %dpx' VALUE=\"%s\"", varName,width, initialVal);
if (isNotEmpty(extra))
    printf(" %s",extra);
printf(">\n");
}


void cgiMakeIntVar(char *varName, int initialVal, int maxDigits)
/* Make a text control filled with initial value.  */
{
if (maxDigits == 0) maxDigits = 4;

printf("<INPUT TYPE=TEXT NAME=\"%s\" SIZE=%d VALUE=%d>", varName,
	maxDigits, initialVal);
}

void cgiMakeIntVarInRange(char *varName, int initialVal, char *title, int width, char *min, char *max)
/* Make a integer control filled with initial value.
   If min and/or max are non-NULL will enforce range
   Requires utils.js jQuery.js and inputBox class */
{
if(width==0)
    {
    if(max)
        width=strlen(max)*10;
    else
        {
        int sz=initialVal+1000;
        if(min)
            sz=atoi(min) + 1000;
        width = 10;
        while(sz/=10)
            width+=10;
        }
    }
if (width < 65)
    width = 65;

printf("<INPUT TYPE=TEXT class='inputBox' name=\"%s\" style='width: %dpx' value=%d",varName,width,initialVal);
printf(" onChange='return validateInt(this,%s,%s);'",(min?min:"\"null\""),(max?max:"\"null\""));
if(title)
    printf(" title='%s'",title);
printf(">\n");
}
void cgiMakeIntVarWithLimits(char *varName, int initialVal, char *title, int width, int min, int max)
{
char minLimit[20];
char maxLimit[20];
char *minStr=NULL;
char *maxStr=NULL;
if(min != NO_VALUE)
    {
    safef(minLimit,sizeof(minLimit),"%d",min);
    minStr = minLimit;
    }
if(max != NO_VALUE)
    {
    safef(maxLimit,sizeof(maxLimit),"%d",max);
    maxStr = maxLimit;
    }
cgiMakeIntVarInRange(varName,initialVal,title,width,minStr,maxStr);
}
void cgiMakeIntVarWithMin(char *varName, int initialVal, char *title, int width, int min)
{
char minLimit[20];
char *minStr=NULL;
if(min != NO_VALUE)
    {
    safef(minLimit,sizeof(minLimit),"%d",min);
    minStr = minLimit;
    }
cgiMakeIntVarInRange(varName,initialVal,title,width,minStr,NULL);
}
void cgiMakeIntVarWithMax(char *varName, int initialVal, char *title, int width, int max)
{
char maxLimit[20];
char *maxStr=NULL;
if(max != NO_VALUE)
    {
    safef(maxLimit,sizeof(maxLimit),"%d",max);
    maxStr = maxLimit;
    }
cgiMakeIntVarInRange(varName,initialVal,title,width,NULL,maxStr);
}

void cgiMakeDoubleVar(char *varName, double initialVal, int maxDigits)
/* Make a text control filled with initial floating-point value.  */
{
if (maxDigits == 0) maxDigits = 4;

printf("<INPUT TYPE=TEXT NAME=\"%s\" SIZE=%d VALUE=%g>", varName,
	maxDigits, initialVal);
}

void cgiMakeDoubleVarInRange(char *varName, double initialVal, char *title, int width, char *min, char *max)
/* Make a floating point control filled with initial value.
   If min and/or max are non-NULL will enforce range
   Requires utils.js jQuery.js and inputBox class */
{
if(width==0)
    {
    if(max)
        width=strlen(max)*10;
    }
if (width < 65)
    width = 65;

printf("<INPUT TYPE=TEXT class='inputBox' name=\"%s\" style='width: %dpx' value=%g",varName,width,initialVal);
printf(" onChange='return validateFloat(this,%s,%s);'",(min?min:"\"null\""),(max?max:"\"null\""));
if(title)
    printf(" title='%s'",title);
printf(">\n");
}

void cgiMakeDoubleVarWithLimits(char *varName, double initialVal, char *title, int width, double min, double max)
{
char minLimit[20];
char maxLimit[20];
char *minStr=NULL;
char *maxStr=NULL;
if((int)min != NO_VALUE)
    {
    safef(minLimit,sizeof(minLimit),"%g",min);
    minStr = minLimit;
    }
if((int)max != NO_VALUE)
    {
    safef(maxLimit,sizeof(maxLimit),"%g",max);
    maxStr = maxLimit;
    }
cgiMakeDoubleVarInRange(varName,initialVal,title,width,minStr,maxStr);
}

void cgiMakeDoubleVarWithMin(char *varName, double initialVal, char *title, int width, double min)
{
char minLimit[20];
char *minStr=NULL;
if((int)min != NO_VALUE)
    {
    safef(minLimit,sizeof(minLimit),"%g",min);
    minStr = minLimit;
    }
cgiMakeDoubleVarInRange(varName,initialVal,title,width,minStr,NULL);
}
void cgiMakeDoubleVarWithMax(char *varName, double initialVal, char *title, int width, double max)
{
char maxLimit[20];
char *maxStr=NULL;
if((int)max != NO_VALUE)
    {
    safef(maxLimit,sizeof(maxLimit),"%g",max);
    maxStr = maxLimit;
    }
cgiMakeDoubleVarInRange(varName,initialVal,title,width,NULL,maxStr);
}

void cgiMakeDropListClassWithStyleAndJavascript(char *name, char *menu[],
    int menuSize, char *checked, char *class, char *style,char *javascript)
/* Make a drop-down list with names, text class, style and javascript. */
{
int i;
char *selString;
if (checked == NULL) checked = menu[0];
printf("<SELECT");
if (name)
    printf(" NAME='%s'", name);
if (class)
    printf(" class='%s'", class);
if (style)
    printf(" style='%s'", style);
if (javascript)
    printf(" %s", javascript);
printf(">\n");
for (i=0; i<menuSize; ++i)
    {
    if (sameWord(menu[i], checked))
        selString = " SELECTED";
    else
        selString = "";
    printf("<OPTION%s>%s</OPTION>\n", selString, menu[i]);
    }
printf("</SELECT>\n");
}

void cgiMakeDropListClassWithStyle(char *name, char *menu[],
    int menuSize, char *checked, char *class, char *style)
/* Make a drop-down list with names, text class and style. */
{
    cgiMakeDropListClassWithStyleAndJavascript(name,menu,menuSize,checked,class,style,"");
}

void cgiMakeDropListClass(char *name, char *menu[],
	int menuSize, char *checked, char *class)
/* Make a drop-down list with names. */
{
    cgiMakeDropListClassWithStyle(name, menu, menuSize, checked,
                                        class, NULL);
}

void cgiMakeDropList(char *name, char *menu[], int menuSize, char *checked)
/* Make a drop-down list with names.
 * uses style "normalText" */
{
    cgiMakeDropListClass(name, menu, menuSize, checked, "normalText");
}

char *cgiMultListShadowPrefix()
/* Prefix for shadow variable set with multi-select inputs. */
{
return "multishad.";
}

void cgiMakeMultList(char *name, char *menu[], int menuSize, struct slName *checked, int length)
/* Make a list of names with window height equalt to length,
 * which can have multiple selections. Same as drop-down list
 * except "multiple" is added to select tag. */
{
int i;
char *selString;
if (checked == NULL) checked = slNameNew(menu[0]);
printf("<SELECT MULTIPLE SIZE=%d ALIGN=CENTER NAME=\"%s\">\n", length, name);
for (i=0; i<menuSize; ++i)
    {
    if (slNameInList(checked, menu[i]))
        selString = " SELECTED";
    else
        selString = "";
    printf("<OPTION%s>%s</OPTION>\n", selString, menu[i]);
    }
printf("</SELECT>\n");
char buf[512];
safef(buf, sizeof(buf), "%s%s", cgiMultListShadowPrefix(), name);
cgiMakeHiddenVar(buf, "1");
}

void cgiMakeCheckboxGroupWithVals(char *name, char *menu[], char *values[], int menuSize,
				  struct slName *checked, int tableColumns)
/* Make a table of checkboxes that have the same variable name but different
 * values (same behavior as a multi-select input), with nice labels in menu[]. */
{
int i;
if (values == NULL) values = menu;
if (menu == NULL) menu = values;
puts("<TABLE BORDERWIDTH=0><TR>");
for (i = 0;  i < menuSize;  i++)
    {
    if (i > 0 && (i % tableColumns) == 0)
	printf("</TR><TR>");
    printf("<TD><INPUT TYPE=CHECKBOX NAME=\"%s\" VALUE=\"%s\" %s> %s</TD>\n", name, values[i],
	   (slNameInList(checked, values[i]) ? "CHECKED" : ""), menu[i]);
    }
if ((i % tableColumns) != 0)
    while ((i++ % tableColumns) != 0)
	printf("<TD></TD>");
puts("</TR></TABLE>");
char buf[512];
safef(buf, sizeof(buf), "%s%s", cgiMultListShadowPrefix(), name);
cgiMakeHiddenVar(buf, "0");
}

void cgiMakeCheckboxGroup(char *name, char *menu[], int menuSize, struct slName *checked,
			  int tableColumns)
/* Make a table of checkboxes that have the same variable name but different
 * values (same behavior as a multi-select input). */
{
cgiMakeCheckboxGroupWithVals(name, menu, NULL, menuSize, checked, tableColumns);
}

void cgiMakeDropListFull(char *name, char *menu[], char *values[],
                         int menuSize, char *checked, char *extraAttribs)
/* Make a drop-down list with names and values. */
{
int i;
char *selString;
if (checked == NULL) checked = menu[0];

if (NULL != extraAttribs)
    {
    printf("<SELECT NAME=\"%s\" %s>\n", name, extraAttribs);
    }
else
    {
    printf("<SELECT NAME=\"%s\">\n", name);
    }

for (i=0; i<menuSize; ++i)
    {
    if (sameWord(values[i], checked))
        selString = " SELECTED";
    else
        selString = "";
    printf("<OPTION%s VALUE=\"%s\">%s</OPTION>\n", selString, values[i], menu[i]);
    }
printf("</SELECT>\n");
}

char *cgiMakeSelectDropList(boolean multiple, char *name, struct slPair *valsAndLabels,char *selected, char *anyAll,char *extraClasses, char *extraHtml)
// Returns allocated string of HTML defining a drop-down select (if multiple, REQUIRES ui-dropdownchecklist.js)
// In valsAndLabels, val (pair->name) must be filled in but label (pair->val) may be NULL.
// selected, if not NULL is a val found in the valsAndLabels (multiple then comma delimited list).  If null and anyAll not NULL, that will be selected
// anyAll, if not NULL is the string for an initial option.  It can contain val and label, delimited by a comma
// extraHtml, if not NULL contains id, javascript calls and style.  It does NOT contain class definitions
{
struct dyString *output = dyStringNew(1024);
boolean checked = FALSE;

dyStringPrintf(output,"<SELECT name='%s'",name);
if (multiple)
    dyStringAppend(output," MULTIPLE");
if (extraClasses != NULL)
    dyStringPrintf(output," class='%s%s'",extraClasses,(multiple?" filterBy":""));
else if (multiple)
    dyStringAppend(output," class='filterBy'");

if (extraHtml != NULL)
    dyStringPrintf(output," %s",extraHtml);
dyStringAppend(output,">\n");

// Handle initial option "Any" or "All"
if (anyAll != NULL)
    {
    char *val = anyAll;  // Could contain a label after the value
    char *label = strchr(val,',');  // Could contain a label after the value
    if (label != NULL)
        {
        val = cloneString(anyAll);
        label = strchr(val,',');  // again because this is new mem
        *label = '\0';
        label = label+1;
        }
    else
        label = val;
    checked = TRUE; // The default case
    if (selected != NULL)
        {
        if (multiple)
            checked = (findWordByDelimiter(val,',', selected) != NULL);
        else
            checked = sameString(val,selected);
        }
    dyStringPrintf(output, "<OPTION%s VALUE='%s'>%s</OPTION>\n",(checked?" SELECTED":""),
                   val, javaScriptLiteralEncode(label));
    if (label != val)
        freeMem(val);
    }

// All other options
struct slPair *valPair = valsAndLabels;
for (; valPair != NULL; valPair = valPair->next)
    {
    checked = FALSE;
    if (selected != NULL)
        {
        if (multiple)
            checked = (findWordByDelimiter(valPair->name,',', selected) != NULL);
        else
            checked = sameString(valPair->name,selected);
        }
    char *label = valPair->name;
    if (valPair->val != NULL)
        label = valPair->val;
    dyStringPrintf(output, "<OPTION%s VALUE='%s'>%s</OPTION>\n",(checked?" SELECTED":""),
                   (char *)valPair->name, javaScriptLiteralEncode(label));
    }

dyStringPrintf(output,"</SELECT>\n");

return dyStringCannibalize(&output);
}

void cgiMakeDropListWithVals(char *name, char *menu[], char *values[],
                         int menuSize, char *checked)
/* Make a drop-down list with names and values. In this case checked
 * corresponds to a value, not a menu. */
{
int i;
char *selString;
if (checked == NULL) checked = values[0];

printf("<SELECT NAME=\"%s\">\n", name);
for (i=0; i<menuSize; ++i)
    {
    if (sameWord(values[i], checked))
        selString = " SELECTED";
    else
        selString = "";
    printf("<OPTION%s VALUE=\"%s\">%s</OPTION>\n", selString, values[i], menu[i]);
    }
printf("</SELECT>\n");
}

void cgiDropDownWithTextValsAndExtra(char *name, char *text[], char *values[],
    int count, char *selected, char *extra)
/* Make a drop-down list with both text and values. */
{
int i;
char *selString;
assert(values != NULL && text != NULL);
if (selected == NULL)
    selected = values[0];
printf("<SELECT");
if (name)
    printf(" NAME='%s'", name);
if (extra)
    printf("%s", extra);
printf(">\n");
for (i=0; i<count; ++i)
    {
    if (sameWord(values[i], selected))
        selString = " SELECTED";
    else
        selString = "";
    printf("<OPTION%s value='%s'>%s</OPTION>\n", selString, values[i], text[i]);
    }
printf("</SELECT>\n");
}


void cgiMakeHiddenVarWithExtra(char *varName, char *string,char *extra)
/* Store string in hidden input for next time around. */
{
printf("<INPUT TYPE=HIDDEN NAME='%s' VALUE='%s'", varName, string);
if(extra)
    printf(" %s>\n",extra);
else
    puts(">");
}

void cgiContinueHiddenVar(char *varName)
/* Write CGI var back to hidden input for next time around. */
{
if (cgiVarExists(varName))
    cgiMakeHiddenVar(varName, cgiString(varName));
}

void cgiVarExclude(char *varName)
/* If varName exists, remove it. */
{
if (cgiVarExists(varName))
    {
    struct cgiVar *cv = hashRemove(inputHash, varName);
    slRemoveEl(&inputList, cv);
    }
}

void cgiVarExcludeExcept(char **varNames)
/* Exclude all variables except for those in NULL
 * terminated array varNames.  varNames may be NULL
 * in which case nothing is excluded. */
{
struct hashEl *list, *el;
struct hash *exclude = newHash(8);
char *s;

/* Build up hash of things to exclude */
if (varNames != NULL)
   {
   while ((s = *varNames++) != NULL)
       hashAdd(exclude, s, NULL);
   }

/* Step through variable list and remove them if not
 * excluded. */
initCgiInput();
list = hashElListHash(inputHash);
for (el = list; el != NULL; el = el->next)
    {
    if (!hashLookup(exclude, el->name))
        cgiVarExclude(el->name);
    }
hashElFreeList(&list);
freeHash(&exclude);
}

void cgiVarSet(char *varName, char *val)
/* Set a cgi variable to a particular value. */
{
struct cgiVar *var;
initCgiInput();
AllocVar(var);
var->val = cloneString(val);
hashAddSaveName(inputHash, varName, var, &var->name);
}

struct dyString *cgiUrlString()
/* Get URL-formatted that expresses current CGI variable state. */
{
struct dyString *dy = newDyString(0);
struct cgiVar *cv;
char *e;


for (cv = inputList; cv != NULL; cv = cv->next)
    {
    if (cv != inputList)
       dyStringAppend(dy, "&");
    e = cgiEncode(cv->val);
    dyStringPrintf(dy, "%s=", cv->name);
    dyStringAppend(dy, e);
    freez(&e);
    }
return dy;
}

void cgiContinueAllVars()
/* Write back all CGI vars as hidden input for next time around. */
{
struct cgiVar *cv;
for (cv = inputList; cv != NULL; cv = cv->next)
    cgiMakeHiddenVar(cv->name, cv->val);
}


boolean cgiFromCommandLine(int *pArgc, char *argv[], boolean preferWeb)
/* Use the command line to set up things as if we were a CGI program.
 * User types in command line (assuming your program called cgiScript)
 * like:
 *        cgiScript nonCgiArg1 var1=value1 var2=value2 var3=value3 nonCgiArg2
 * or like
 *        cgiScript nonCgiArg1 var1=value1&var2=value2&var3=value3 nonCgiArg2
 * or even like
 *        cgiScript nonCgiArg1 -x -y=bogus z=really
 * (The non-cgi arguments can occur anywhere.  The cgi arguments (all containing
 * the character '=' or starting with '-') are erased from argc/argv.  Normally
 * you call this cgiSpoof(&argc, argv);
 */
{
int argc = *pArgc;
int i;
int argcLeft = argc;
char *name;
static char queryString[16384];
char *q = queryString;
boolean needAnd = TRUE;
boolean gotAny = FALSE;
boolean startDash;
boolean gotEq;
static char hostLine[512];

if (preferWeb && cgiIsOnWeb())
    return TRUE;	/* No spoofing required! */
q += safef(q, queryString + sizeof(queryString) - q,
	   "%s", "QUERY_STRING=cgiSpoof=on");
for (i=0; i<argcLeft; )
    {
    name = argv[i];
    if ((startDash = (name[0] == '-')))
       ++name;
    gotEq = (strchr(name, '=') != NULL);
    if (gotEq || startDash)
        {
        if (needAnd)
            *q++ = '&';
        q += safef(q, queryString + sizeof(queryString) - q, "%s", name);
        if (!gotEq || strchr(name, '&') == NULL)
            needAnd = TRUE;
	if (!gotEq)
	    q += safef(q, queryString + sizeof(queryString) - q, "=on");
        memcpy(&argv[i], &argv[i+1], sizeof(argv[i]) * (argcLeft-i-1));
        argcLeft -= 1;
        gotAny = TRUE;
        }
    else
        i++;
    }
if (gotAny)
    {
    *pArgc = argcLeft;
    }
putenv("REQUEST_METHOD=GET");
putenv(queryString);
char *host = getenv("HOST");
if (host == NULL)
    host = "unknown";
safef(hostLine, sizeof(hostLine), "SERVER_NAME=%s", host);
putenv(hostLine);
initCgiInput();
return gotAny;
}

boolean cgiSpoof(int *pArgc, char *argv[])
/* If run from web line set up input
 * variables from web line, otherwise
 * set up from command line. */
{
return cgiFromCommandLine(pArgc, argv, TRUE);
}

boolean cgiFromFile(char *fileName)
/* Set up a cgi environment using parameters stored in a file.
 * Takes file with arguments in the form:
 *       argument1=someVal
 *       # This is a comment
 *       argument2=someOtherVal
 *       ...
 * and puts them into the cgi environment so that the usual
 * cgiGetVar() commands can be used. Useful when a program
 * has a lot of possible parameters.
 */
{
char **argv = NULL;
int argc = 0;
int maxArgc = 10;
int i;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *line;
boolean spoof= FALSE;
AllocArray(argv, maxArgc);
/* Remember that first arg is program name.
   Put filename there instead. */
argc = 1;
argv[0] = cloneString(fileName);
for(;;)
    {
    /* If we are at the end we're done. */
    if(!lineFileNext(lf, &line, NULL))
	break;
    /* If it is a comment or blank line skip it. */
    if (line[0] == '#' || sameString(line, ""))
        continue;
    /* If our argv array is full expand it. */
    if((argc+1) >= maxArgc)
	{
	ExpandArray(argv, maxArgc, 2*maxArgc);
	maxArgc *= 2;
	}
    /* Fill in another argument to our psuedo arguments. */
    argv[argc++] = cloneString(line);
    }
spoof = cgiSpoof(&argc, argv);
/* Cleanup. */
lineFileClose(&lf);
for(i=0; i<argc; i++)
    freez(&argv[i]);
freez(&argv);
return spoof;
}

void logCgiToStderr()
/* Log useful CGI info to stderr */
{
char *ip = getenv("REMOTE_ADDR");
char *cgiBinary = getenv("SCRIPT_FILENAME");
char *requestUri = getenv("REQUEST_URI");
char *hgsid = cgiOptionalString("hgsid");
char *cgiFileName = NULL;
time_t nowTime = time(NULL);
struct tm *tm;
tm = localtime(&nowTime);
char *ascTime = asctime(tm);
size_t len = strlen(ascTime);
if (cgiBinary)
    cgiFileName = basename(cgiBinary);
else
    cgiFileName = "cgi-bin";
if (len > 3) ascTime[len-2] = '\0';
if (!ip)
    ip = "unknown";
if (!hgsid)
    hgsid = "unknown";
if (!requestUri)
    requestUri = "unknown";
fprintf(stderr, "[%s] [%s] [client %s] [hgsid=%.24s] [%.1024s] ", ascTime, cgiFileName, ip,
	hgsid, requestUri);
}

void cgiResetState()
/* This is for reloading CGI settings multiple times in the same program
 * execution.  No effect if state has not yet been initialized. */
{
freez(&inputString);
inputString = NULL;
if (inputHash != NULL)
    hashFree(&inputHash);
inputHash = NULL;
inputList = NULL;

haveCookiesHash = FALSE;
if (cookieHash != NULL)
    hashFree(&cookieHash);
cookieHash = NULL;
cookieList = NULL;
}

void cgiDown(float lines)
// Drop down a certain number of lines (may be fractional)
{
printf("<div style='height:%fem;'></div>\n",lines);
}

char *commonCssStyles()
/* Returns a string of common CSS styles */
{
// Contents currently is OBSOLETE as these have been moved to HGStyle.css
// However, don't loose this function call yet, as it may have future uses.
return "";
#ifdef OMIT
static boolean commonStylesWritten = FALSE;
if(commonStylesWritten)
    return "";
commonStylesWritten = TRUE;

struct dyString *style = newDyString(256);

dyStringPrintf(style,"\n<style type='text/css'>\n");
dyStringPrintf(style,".ghost {background-color:%s;}\n",COLOR_BG_GHOST);
dyStringPrintf(style,".pale {background-color:%s;}\n",COLOR_BG_PALE);

// These are for dreagReorder: both in imageV2 and in hgTrackUi subtrack list
dyStringPrintf(style,".trDrag {background-color:%s;}\n",COLOR_LTGREEN);
dyStringPrintf(style,".dragHandle {cursor: s-resize;}\n");

// These are for imageV2 sideButtons:
dyStringPrintf(style,".btn  {border-style:outset; background-color:%s; border-color:%s;}\n",COLOR_LTGREY,COLOR_DARKGREY);
dyStringPrintf(style,".btnN {border-width:1px 1px 1px 1px; margin:1px 1px 0px 1px;}\n"); // connect none
dyStringPrintf(style,".btnU {border-width:0px 1px 1px 1px; margin:0px 1px 0px 1px;}\n"); // connect up
dyStringPrintf(style,".btnD {border-width:1px 1px 0px 1px; margin:1px 1px 0px 1px;}\n"); // connect down
dyStringPrintf(style,".btnL {border-width:0px 1px 0px 1px; margin:0px 1px 0px 1px;}\n"); // connect linear
dyStringPrintf(style,".btnBlue {background-color:#91B3E6; border-color:%s;}\n",COLOR_BLUE_BUTTON);

// Common boxes
dyStringPrintf(style,".inputBox {border: 2px inset %s;}\n",COLOR_LTGREY);
dyStringPrintf(style,".greenRoof {border-top: 3px groove %s;}\n",COLOR_DARKGREEN);
//dyStringPrintf(style,".greenFloor {border-bottom: 3px ridge %s;}\n",COLOR_DARKGREEN);      // Unused
//dyStringPrintf(style,".hiddenRoof {border-top: 0px solid %s;}\n",COLOR_BG_ALTDEFAULT);     // Doesn't work
//dyStringPrintf(style,".hiddenFloor {border-bottom: 0px solid %s;}\n",COLOR_BG_ALTDEFAULT); // Doesn't work
dyStringPrintf(style,".greenBox {border: 5px outset %s;}\n",COLOR_DARKGREEN); // Matrix
dyStringPrintf(style,".blueBox {border: 4px inset %s;}\n",COLOR_DARKBLUE);    // cfg box
dyStringPrintf(style,".redBox {border: 3px ridge %s; background:Beige; padding:10px; margin:10px; text-align:left;}\n",COLOR_RED); // Special alert

// Experiments with squeezing giant matrices
//dyStringPrintf(style,".slantUp {-moz-transform:rotate(-75deg); -moz-transform-origin: bottom left; -webkit-transform:rotate(-75deg); -webkit-transform-origin: bottom left; white-space:nowrap; position:relative;left: 16px; filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=3)}\n");
//dyStringPrintf(style,".slantDn {-moz-transform:rotate( 75deg); -moz-transform-origin: top left;    -webkit-transform:rotate( 75deg); -webkit-transform-origin: top left;    white-space:nowrap; position:relative;left: 16px; filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=3)}\n");
//dyStringPrintf(style,".rotate90 {-webkit-transform: rotate(-90deg); -moz-transform: rotate(-90deg);}\n");

dyStringPrintf(style,".hintCol {font-size:70%%; line-height:80%%; border-style: hidden; background-color:%s;}\n",COLOR_BG_ALTDEFAULT);
dyStringPrintf(style,".hintRow {font-size:70%%; line-height:80%%; border-style: hidden; background-color:%s;}\n",COLOR_BG_ALTDEFAULT);
//dyStringPrintf(style,".halfVis {opacity: 0.5; filters.alpha.opacity=50;}\n");   // not ready for prime time because ff and ie can't agree

// waitMask allows waiting on long running javascript using utils.js::waitOnFunction
dyStringPrintf(style,".waitMask {display: none; cursor: wait; z-index: 9999; position: absolute; top: 0; left: 0; height: 100%%; width: 100%%; background-color: #fff; opacity: 0;}\n");
dyStringPrintf(style,".inOutButton {height:24px; width:24px; border-style: outset;}\n"); // A [+][-] button can be toggled by waitOnFunction during long running scripts

dyStringPrintf(style,"</style>\n");
return dyStringCannibalize(&style);
#endif///def OMIT
}

