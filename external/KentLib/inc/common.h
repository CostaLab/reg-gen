/* Common.h - functions that are commonly used.  Includes
 * routines for managing singly linked lists, some basic
 * string manipulation stuff, and other stuff of the
 * short but useful nature.
 *
 * This file is copyright 2002-2005 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef COMMON_H	/* Wrapper to avoid including this twice. */
#define COMMON_H

/* Some stuff to support large files in Linux. */
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE 1
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

/* Some stuff for safer pthreads. */
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>
#include <fcntl.h>
#include <assert.h>
#include <setjmp.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>

#if defined(MACHTYPE_ppc)
#include <sys/wait.h>
#endif

#if defined(__APPLE__)
#if defined(__i686__)
/* The i686 apple math library defines warn. */
#define warn jkWarn
#endif
#endif

#ifdef __CYGWIN32__
#include <mingw/math.h>
#endif

#ifndef NAN
#define NAN (0.0 / 0.0)
#endif

#ifndef WIFEXITED
#define WIFEXITED(stat)  (((*((int *) &(stat))) & 0xff) == 0)
#endif

#ifndef WEXITSTATUS
#define WEXITSTATUS(stat) (short)(((*((int *) &(stat))) >> 8) & 0xffff)
#endif

#ifndef WIFSIGNALED
#define WIFSIGNALED(stat) (((*((int *) &(stat)))) && ((*((int *) &(stat))) == ((*((int *) &(stat))) &0x00ff)))
#endif

#ifndef WTERMSIG
#define WTERMSIG(stat)    ((*((int *) &(stat))) & 0x7f)
#endif

#ifndef WIFSTOPPED
#define WIFSTOPPED(stat)  (((*((int *) &(stat))) & 0xff) == 0177)
#endif

#ifndef WSTOPSIG
#define WSTOPSIG(stat)    (((*((int *) &(stat))) >> 8) & 0xff)
#endif

#ifndef HUGE
#define HUGE MAXFLOAT
#endif


/* Let's pretend C has a boolean type. */
#define TRUE 1
#define FALSE 0
#define boolean int
#ifndef	__cplusplus
#ifndef bool
#define bool char
#endif
#endif

/* Some other type synonyms */
#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short	      /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int	      /* Wants to be signed 32 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */

#define BIGNUM 0x3fffffff	/* A really big number */
#define BIGDOUBLE 1.7E+308	/* Close to biggest double-precision number */

#define LIMIT_2or8GB (2147483647 * ((sizeof(size_t)/4)*(sizeof(size_t)/4)))
/*      == 2 Gb for 32 bit machines, 8 Gb for 64 bit machines */
#define LIMIT_2or6GB (2147483647 + (2147483647 * ((sizeof(size_t)/4)-1)) + \
	(2147483647 * ((sizeof(size_t)/4)-1)))
/*      == 2 Gb for 32 bit machines, 6 Gb for 64 bit machines */

/* Default size of directory path, file name and extension string buffers */
#define PATH_LEN 512
#define FILENAME_LEN 128
#define FILEEXT_LEN 64

/* inline functions: To declare a function inline, place the entire function
 * in a header file and prefix it with the INLINE macro.  If used with a
 * compiler that doesn't support inline, change the INLINE marco to be simply
 * `static'.
 */
#ifndef INLINE
#define INLINE static inline
#endif

/* stdargs compatibility: in a unix universe a long time ago, types of va_list
 * were passed by value.  It was assume one could do things like:
 *
 *     va_start(args);
 *     vfprintf(fh1, fmt, args);
 *     vfprintf(fh2, fmt, args);
 *     va_end(args);
 *
 * and life would good.  However this is not true on some modern systems (for
 * instance gcc/glibc on x86_64), where va_args can be a pointer to some type
 * of object).  The second call to vfprintf() would then crash, since the
 * first call modified the object that va_args was pointing to. C99 adds a
 * va_copy macro that to address this issue.  Many non-C99 system include this
 * macro, sometimes called __va_copy.  Here we ensure that va_copy is defined.
 * If if doesn't exist, we try to define it in terms of __va_copy.  If that is
 * not available, we make the assumption that va_list can be copied by value
 * and create our own.  Our implementation is the same as used on Solaris.
 */
#if defined(__va_copy) && !defined(va_copy)
#   define va_copy __va_copy
#endif
#if !defined(va_copy)
#   define va_copy(to, from) ((to) = (from))
#endif

/* Cast a pointer to a long long. Use to printf format points as long-longs
 * in a 32/64bit portable manner.  Format should use %llx for the result.
 * Needed because casting a pointer to a different sized number cause a
 * warning with gcc */
#define ptrToLL(p) ((long long)((size_t)p))

/* How big is this array? */
#define ArraySize(a) (sizeof(a)/sizeof((a)[0]))

#define uglyf printf  /* debugging printf */
#define uglyAbort errAbort /* debugging error abort. */
#define uglyOut stdout /* debugging fprintf target. */

void *needMem(size_t size);
/* Need mem calls abort if the memory allocation fails. The memory
 * is initialized to zero. */

void *needLargeMem(size_t size);
/* This calls abort if the memory allocation fails. The memory is
 * not initialized to zero. */

void *needLargeZeroedMem(size_t size);
/* Request a large block of memory and zero it. */

void *needLargeMemResize(void* vp, size_t size);
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL,
 * a new memory block is allocated.  Memory not initted. */

void *needLargeZeroedMemResize(void* vp, size_t oldSize, size_t newSize);
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL, a
 * new memory block is allocated.  If block is grown, new memory is zeroed. */

void *needHugeMem(size_t size);
/* No checking on size.  Memory not initted to 0. */

void *needHugeZeroedMem(size_t size);
/* Request a large block of memory and zero it. */

void *needHugeMemResize(void* vp, size_t size);
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL,
 * a new memory block is allocated.  No checking on size.  Memory not
 * initted. */

void *needHugeZeroedMemResize(void* vp, size_t oldSize, size_t newSize);
/* Adjust memory size on a block, possibly relocating it.  If vp is NULL, a
 * new memory block is allocated.  No checking on size.  If block is grown,
 * new memory is zeroed. */

void *needMoreMem(void *old, size_t copySize, size_t newSize);
/* Allocate a new buffer, copy old buffer to it, free old buffer. */

void *cloneMem(void *pt, size_t size);
/* Allocate a new buffer of given size, and copy pt to it. */

#define CloneVar(pt) cloneMem(pt, sizeof((pt)[0]))
/* Allocate copy of a structure. */

void *wantMem(size_t size);
/* Want mem just calls malloc - no zeroing of memory, no
 * aborting if request fails. */

void freeMem(void *pt);
/* Free memory will check for null before freeing. */

void freez(void *ppt);
/* Pass address of pointer.  Will free pointer and set it
 * to NULL. Typical use:
 *     s = needMem(1024);
 *          ...
 *     freez(&s); */

#define AllocVar(pt) (pt = needMem(sizeof(*pt)))
/* Shortcut to allocating a single variable on the heap and
 * assigning pointer to it. */

#define AllocArray(pt, size) (pt = needLargeZeroedMem(sizeof(*pt) * (size)))

#define AllocA(type) needMem(sizeof(type))
/* Shortcut to allocating a variable on heap of a specific type. */

#define AllocN(type,count) ((type*)needLargeZeroedMem(sizeof(type) * (count)))
/* Shortcut to allocating an array on the heap of a specific type. */

#define ExpandArray(array, oldCount, newCount) \
  (array = needMoreMem((array), (oldCount)*sizeof((array)[0]), (newCount)*sizeof((array)[0])))
/* Expand size of dynamically allocated array. */

#define CopyArray(source, dest,count) memcpy(dest,source,(count)*sizeof(dest[0]))
/* Copy count elements of array from source to dest. */

#define CloneArray(a, count) cloneMem(a, (count)*sizeof(a[0]))
/* Make new dynamic array initialized with  count elements of a */

void errAbort(char *format, ...)
/* Abort function, with optional (printf formatted) error message. */
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

void errnoAbort(char *format, ...)
/* Prints error message from UNIX errno first, then does errAbort. */
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

#define internalErr()  errAbort("Internal error %s %d", __FILE__, __LINE__)
/* Generic internal error message */

void warn(char *format, ...)
/* Issue a warning message. */
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

void warnWithBackTrace(char *format, ...)
/* Issue a warning message and append backtrace. */
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

void verbose(int verbosity, char *format, ...)
/* Write printf formatted message to log (which by
 * default is stdout) if global verbose variable
 * is set to verbosity or higher.  Default level is 1. */
#if defined(__GNUC__)
__attribute__((format(printf, 2, 3)))
#endif
    ;

void verboseTimeInit(void);
/* Initialize or reinitialize the previous time for use by verboseTime. */

void verboseTime(int verbosity, char *label, ...)
/* Print label and how long it's been since last call.  Start time can be
 * initialized with verboseTimeInit, otherwise the elapsed time will be
 * zero. */
#if defined(__GNUC__)
__attribute__((format(printf, 2, 3)))
#endif
    ;

void verboseDot();
/* Write I'm alive dot (at verbosity level 1) */

int verboseLevel();
/* Get verbosity level. */

void verboseSetLevel(int verbosity);
/* Set verbosity level in log.  0 for no logging,
 * higher number for increasing verbosity. */

INLINE void zeroBytes(void *vpt, int count)
/* fill a specified area of memory with zeroes */
{
memset(vpt, '\0', count);
}

#define ZeroVar(v) zeroBytes(v, sizeof(*v))

void reverseBytes(char *bytes, long length);
/* Reverse the order of the bytes. */

void reverseInts(int *a, int length);
/* Reverse the order of the integer array. */

void reverseUnsigned(unsigned *a, int length);
/* Reverse the order of the unsigned array. */

void reverseDoubles(double *a, int length);
/* Reverse the order of the double array. */

void reverseStrings(char **a, int length);
/* Reverse the order of the char* array. */

void swapBytes(char *a, char *b, int length);
/* Swap buffers a and b. */

/* Some things to manage simple lists - structures that begin
 * with a pointer to the next element in the list. */
struct slList
    {
    struct slList *next;
    };

int slCount(const void *list);
/* Return # of elements in list.  */

void *slElementFromIx(void *list, int ix);
/* Return the ix'th element in list.  Returns NULL
 * if no such element. */

int slIxFromElement(void *list, void *el);
/* Return index of el in list.  Returns -1 if not on list. */

INLINE void slAddHead(void *listPt, void *node)
/* Add new node to start of list.
 * Usage:
 *    slAddHead(&list, node);
 * where list and nodes are both pointers to structure
 * that begin with a next pointer.
 */
{
struct slList **ppt = (struct slList **)listPt;
struct slList *n = (struct slList *)node;
n->next = *ppt;
*ppt = n;
}

INLINE void slSafeAddHead(void *listPt, void *node)
/* Add new node to start of list.  Now that slAddHead is an inline instead of
 * a macro, this function is obsolete.
 */
{
slAddHead(listPt, node);
}

void slAddTail(void *listPt, void *node);
/* Add new node to tail of list.
 * Usage:
 *    slAddTail(&list, node);
 * where list and nodes are both pointers to structure
 * that begin with a next pointer. This is sometimes
 * convenient but relatively slow.  For longer lists
 * it's better to slAddHead, and slReverse when done.
 */

void *slPopHead(void *listPt);
/* Return head of list and remove it from list. (Fast) */

void *slPopTail(void *listPt);
/* Return tail of list and remove it from list. (Not so fast) */

void *slCat(void *a, void *b);
/* Return concatenation of lists a and b.
 * Example Usage:
 *   struct slName *a = getNames("a");
 *   struct slName *b = getNames("b");
 *   struct slName *ab = slCat(a,b)
 * After this it is no longer safe to use a or b.
 */

void *slLastEl(void *list);
/* Returns last element in list or NULL if none. */

void slReverse(void *listPt);
/* Reverse order of a list.
 * Usage:
 *    slReverse(&list);
 */

typedef int CmpFunction(const void *elem1, const void *elem2);

void slSort(void *pList, CmpFunction *compare);
/* Sort a singly linked list with Qsort and a temporary array.
 * The arguments to the compare function in real, non-void, life
 * are pointers to pointers. */

void slUniqify(void *pList, CmpFunction *compare, void (*free)());
/* Return sorted list with duplicates removed.
 * Compare should be same type of function as slSort's compare (taking
 * pointers to pointers to elements.  Free should take a simple
 * pointer to dispose of duplicate element, and can be NULL. */

boolean slRemoveEl(void *vpList, void *vToRemove);
/* Remove element from doubly linked list.  Usage:
 *    slRemove(&list, el);
 * Returns TRUE if element in list.  */

void slFreeList(void *listPt);
/* Free all elements in list and set list pointer to null.
 * Usage:
 *    slFreeList(&list);
 */

struct slInt
/* List of integers. */
    {
    struct slInt *next;	/* Next in list. */
    int val;		/* Integer value. */
    };

struct slInt *slIntNew(int x);
#define newSlInt slIntNew
/* Return a new double. */

int slIntCmp(const void *va, const void *vb);
/* Compare two slInts. */

int slIntCmpRev(const void *va, const void *vb);
/* Compare two slInts in reverse direction. */

struct slInt * slIntFind(struct slInt *list, int target);
/* Find target in slInt list or return NULL */

void doubleSort(int count, double *array);
/* Sort an array of doubles. */

double doubleMedian(int count, double *array);
/* Return median value in array.  This will sort
 * the array as a side effect. */

void doubleBoxWhiskerCalc(int count, double *array, double *retMin,
	double *retQ1, double *retMedian, double *retQ3, double *retMax);
/* Calculate what you need to draw a box and whiskers plot from an array of doubles. */

struct slDouble
/* List of double-precision numbers. */
    {
    struct slDouble *next;	/* Next in list. */
    double val;			/* Double-precision value. */
    };

struct slDouble *slDoubleNew(double x);
#define newSlDouble slDoubleNew
/* Return a new int. */

int slDoubleCmp(const void *va, const void *vb);
/* Compare two slDoubles. */

double slDoubleMedian(struct slDouble *list);
/* Return median value on list. */

void slDoubleBoxWhiskerCalc(struct slDouble *list, double *retMin,
	double *retQ1, double *retMedian, double *retQ3, double *retMax);
/* Calculate what you need to draw a box and whiskers plot from a list of slDoubles. */

void intSort(int count, int *array);
/* Sort an array of ints. */

int intMedian(int count, int *array);
/* Return median value in array.  This will sort
 * the array as a side effect. */

struct slName
/* List of names. The name array is allocated to accommodate full name
 */
    {
    struct slName *next;	/* Next in list. */
    char name[1];               /* Allocated at run time to length of string. */
    };

struct slName *newSlName(char *name);

#define slNameNew newSlName
/* Return a new slName. */

#define slNameFree freez
/* Free a single slName */

#define slNameFreeList slFreeList
/* Free a list of slNames */

struct slName *slNameNewN(char *name, int size);
/* Return new slName of given size. */

int slNameCmpCase(const void *va, const void *vb);
/* Compare two slNames, ignore case. */

int slNameCmp(const void *va, const void *vb);
/* Compare two slNames. */

int slNameCmpStringsWithEmbeddedNumbers(const void *va, const void *vb);
/* Compare strings such as gene names that may have embedded numbers,
 * so that bmp4a comes before bmp14a */

void slNameSortCase(struct slName **pList);
/* Sort slName list, ignore case. */

void slNameSort(struct slName **pList);
/* Sort slName list. */

boolean slNameInList(struct slName *list, char *string);
/* Return true if string is in name list -- case insensitive. */

boolean slNameInListUseCase(struct slName *list, char *string);
/* Return true if string is in name list -- case sensitive. */

void *slNameFind(void *list, char *string);
/* Return first element of slName list (or any other list starting
 * with next/name fields) that matches string. */

int slNameFindIx(struct slName *list, char *string);
/* Return index of first element of slName list (or any other
 * list starting with next/name fields) that matches string.
 * ... Return -1 if not found. */

char *slNameStore(struct slName **pList, char *string);
/* Put string into list if it's not there already.
 * Return the version of string stored in list. */

struct slName *slNameAddHead(struct slName **pList, char *name);
/* Add name to start of list and return it. */

struct slName *slNameAddTail(struct slName **pList, char *name);
/* Add name to end of list (not efficient for long lists),
 * and return it. */

struct slName *slNameCloneList(struct slName *list);
/* Return clone of list. */

struct slName *slNameListFromString(char *s, char delimiter);
/* Return list of slNames gotten from parsing delimited string.
 * The final delimiter is optional. a,b,c  and a,b,c, are equivalent
 * for comma-delimited lists. */

#define slNameListFromComma(s) slNameListFromString(s, ',')
/* Parse out comma-separated list. */

struct slName *slNameListOfUniqueWords(char *text,boolean respectQuotes);
// Return list of unique words found by parsing string delimited by whitespace.
// If respectQuotes then ["Lucy and Ricky" 'Fred and Ethyl'] will yield 2 slNames no quotes

struct slName *slNameListFromStringArray(char *stringArray[], int arraySize);
/* Return list of slNames from an array of strings of length arraySize.
 * If a string in the array is NULL, the array will be treated as
 * NULL-terminated. */

char *slNameListToString(struct slName *list, char delimiter);
/* Return string created by joining all names with the delimiter. */

struct slName *slNameLoadReal(char *fileName);
/* load file lines that are not blank or start with a '#' into a slName
 * list */

struct slName *slNameIntersection(struct slName *a, struct slName *b);
/* return intersection of two slName lists.  */

struct slRef
/* Singly linked list of generic references. */
    {
    struct slRef *next;	/* Next in list. */
    void *val;		/* A reference to something. */
    };

struct slRef *slRefNew(void *val);
/* Create new slRef element. */

struct slRef *refOnList(struct slRef *refList, void *val);
/* Return ref if val is already on list, otherwise NULL. */

void refAdd(struct slRef **pRefList, void *val);
/* Add reference to list. */

void refAddUnique(struct slRef **pRefList, void *val);
/* Add reference to list if not already on list. */

struct slRef *refListFromSlList(void *list);
/* Make a reference list that mirrors a singly-linked list. */

struct slPair
/* A name/value pair. */
    {
    struct slPair *next;	/* Next in list. */
    char *name;			/* Name of item. */
    void *val;			/* Pointer to item data. */
    };

struct slPair *slPairNew(char *name, void *val);
/* Allocate new name/value pair. */

void slPairAdd(struct slPair **pList, char *name, void *val);
/* Add new slPair to head of list. */

void slPairFree(struct slPair **pEl);
/* Free up struct and name.  (Don't free up values.) */

void slPairFreeList(struct slPair **pList);
/* Free up list.  (Don't free up values.) */

void slPairFreeVals(struct slPair *list);
/* Free up all values on list. */

void slPairFreeValsAndList(struct slPair **pList);
/* Free up all values on list and list itself */

struct slPair *slPairFind(struct slPair *list, char *name);
/* Return list element of given name, or NULL if not found. */

void *slPairFindVal(struct slPair *list, char *name);
/* Return value associated with name in list, or NULL if not found. */

struct slPair *slPairListFromString(char *str,boolean respectQuotes);
// Return slPair list parsed from list in string like:  [name1=val1 name2=val2 ...]
// if respectQuotes then string can have double quotes: [name1="val 1" "name 2"=val2 ...]
//    resulting pair strips quotes: {name1}={val 1},{name 2}={val2}
// Returns NULL if parse error.  Free this up with slPairFreeValsAndList.
#define slPairFromString(s) slPairListFromString(s,FALSE)

char *slPairListToString(struct slPair *list,boolean quoteIfSpaces);
// Returns an allocated string of pairs in form of [name1=val1 name2=val2 ...]
// If requested, will wrap name or val in quotes if contain spaces: [name1="val 1" "name 2"=val2]

char *slPairNameToString(struct slPair *list, char delimiter,boolean quoteIfSpaces);
// Return string created by joining all names (ignoring vals) with the delimiter.
// If requested, will wrap name in quotes if contain spaces: [name1,"name 2" ...]

int slPairCmpCase(const void *va, const void *vb);
/* Compare two slPairs, ignore case. */

void slPairSortCase(struct slPair **pList);
/* Sort slPair list, ignore case. */

int slPairCmp(const void *va, const void *vb);
/* Compare two slPairs. */

int slPairValCmpCase(const void *va, const void *vb);
/* Case insensitive compare two slPairs on their values (must be string). */

int slPairValCmp(const void *va, const void *vb);
/* Compare two slPairs on their values (must be string). */

void slPairValSortCase(struct slPair **pList);
/* Sort slPair list on values (must be string), ignore case. */

void slPairValSort(struct slPair **pList);
/* Sort slPair list on values (must be string). */

int slPairIntCmp(const void *va, const void *vb);
// Compare two slPairs on their integer values.

void slPairIntSort(struct slPair **pList);
// Sort slPair list on integer values.

int slPairAtoiCmp(const void *va, const void *vb);
// Compare two slPairs on their strings interpreted as integer values.

void slPairValAtoiSort(struct slPair **pList);
// Sort slPair list on string values interpreted as integers.

void gentleFree(void *pt);
/* check pointer for NULL before freeing.
 * (Actually plain old freeMem does that these days.) */

/*******  Some stuff for processing strings. *******/

char *cloneStringZ(const char *s, int size);
/* Make a zero terminated copy of string in memory */

char *cloneString(const char *s);
/* Make copy of string in dynamic memory */

char *cloneLongString(char *s);
/* Make clone of long string. */

char *catTwoStrings(char *a, char *b);
/* Allocate new string that is a concatenation of two strings. */

int differentWord(char *s1, char *s2);
/* strcmp ignoring case - returns zero if strings are
 * the same (ignoring case) otherwise returns difference
 * between first non-matching characters. */

#define sameWord(a,b) (!differentWord(a,b))
/* Return TRUE if two strings are same ignoring case */

#define differentString(a,b) (strcmp(a,b))
/* Returns FALSE if two strings same. */

int differentStringNullOk(char *a, char *b);
/* Returns 0 if two strings (either of which may be NULL)
 * are the same.  Otherwise it returns a positive or negative
 * number depending on the alphabetical order of the two
 * strings.
 * This is basically a strcmp that can handle NULLs in
 * the input.  If used in a sort the NULLs will end
 * up before any of the cases with data.   */

#define sameOk(a,b) (differentStringNullOk(a,b) == 0)
/* returns TRUE if two strings same, NULLs OK */

#define sameString(a,b) (strcmp(a,b)==0)
/* Returns TRUE if two strings same. */

#define sameStringN(a,b,c) (strncmp(a,b,c)==0)
/* Returns TRUE if two strings start with the same c characters. */

#define isEmpty(string) ((string) == NULL || (string)[0] == 0)
#define isNotEmpty(string) (! isEmpty(string))

int cmpStringsWithEmbeddedNumbers(const char *a, const char *b);
/* Compare strings such as gene names that may have embedded numbers,
 * so that bmp4a comes before bmp14a */

int cmpWordsWithEmbeddedNumbers(const char *a, const char *b);
/* Case insensitive version of cmpStringsWithEmbeddedNumbers. */

boolean startsWith(const char *start, const char *string);
/* Returns TRUE if string begins with start. */

boolean startsWithWord(char *firstWord, char *line);
/* Return TRUE if first white-space-delimited word in line
 * is same as firstWord.  Comparison is case sensitive. */

boolean startsWithWordByDelimiter(char *firstWord,char delimit, char *line);
/* Return TRUE if first word in line is same as firstWord as delimited by delimit.
   Comparison is case sensitive. Delimit of ' ' uses isspace() */

#define stringIn(needle, haystack) strstr(haystack, needle)
/* Returns position of needle in haystack or NULL if it's not there. */
/*        char *stringIn(char *needle, char *haystack);      */

char *rStringIn(char *needle, char *haystack);
/* Return last position of needle in haystack, or NULL if it's not there. */

char *stringBetween(char *start, char *end, char *haystack);
/* Return string between start and end strings, or NULL if
 * none found.  The first such instance is returned.
 * String must be freed by caller. */

char * findWordByDelimiter(char *word,char delimit, char *line);
/* Return pointer to first word in line matching 'word' and delimited
   by 'delimit'. Comparison is case sensitive. Delimit of ' ' uses isspace() */

boolean endsWith(char *string, char *end);
/* Returns TRUE if string ends with end. */

char lastChar(char *s);
/* Return last character in string. */

char *lastNonwhitespaceChar(char *s);
// Return pointer to last character in string that is not whitespace.

char *matchingCharBeforeInLimits(char *limit, char *s, char c);
/* Look for character c sometime before s, but going no further than limit.
 * Return NULL if not found. */

boolean wildMatch(const char *wildCard, const char *string);
/* does a case insensitive wild card match with a string.
 * * matches any string or no character.
 * ? matches any single character.
 * anything else etc must match the character exactly. */

boolean sqlMatchLike(char *wildCard, char *string);
/* Match using % and _ wildcards. */

boolean anyWild(const char *string);
/* Return TRUE if any wild card characters in string. */

char *memMatch(char *needle, int nLen, char *haystack, int hLen);
/* Returns first place where needle (of nLen chars) matches
 * haystack (of hLen chars) */

void toUpperN(char *s, int n);
/* Convert a section of memory to upper case. */

void toLowerN(char *s, int n);
/* Convert a section of memory to lower case. */

void toggleCase(char *s, int size);
/* toggle upper and lower case chars in string. */

char *strUpper(char *s);
#define touppers(s) (void)strUpper(s)
/* Convert entire string to upper case. */

char *strLower(char *s);
#define tolowers(s) (void)strLower(s)
/* Convert entire string to lower case */

char *replaceChars(char *string, char *oldStr, char *newStr);
/*
  Replaces the old with the new.
 The old and new string need not be of equal size
 Can take any length string.
 Return value needs to be freeMem'd.
*/

int strSwapStrs(char *string, int sz,char *oldStr, char *newStr);
/* Swaps all occurrences of the oldStr with the newStr in string. Need not be same size
   Swaps in place but restricted by sz.  Returns count of swaps or -1 for sz failure.*/

char * memSwapChar(char *s, int len, char oldChar, char newChar);
/* Substitute newChar for oldChar throughout memory of given length.
   old or new may be null */
#define strSwapChar(s,old,new)   memSwapChar((s),strlen(s),(old),(new))
#define subChar(s,old,new) (void)memSwapChar((s),strlen(s),(old),(new))
/* Substitute newChar for oldChar throughout string s. */

void stripChar(char *s, char c);
/* Remove all occurences of c from s. */

char *stripEnclosingChar(char *inout,char encloser);
// Removes enclosing char if found at both beg and end, preserving pointer
// Note: handles brackets '(','{' and '[' by complement at end
#define stripEnclosingDoubleQuotes(inout) stripEnclosingChar((inout),'"')
#define stripEnclosingSingleQuotes(inout) stripEnclosingChar((inout),'\'')

void stripString(char *s, char *strip);
/* Remove all occurences of strip from s. */

int countChars(char *s, char c);
/* Return number of characters c in string s. */

int countCharsN(char *s, char c, int size);
/* Return number of characters c in string s of given size. */

int countLeadingChars(char *s, char c);
/* Count number of characters c at start of string. */

int countLeadingDigits(const char *s);
/* Return number of leading digits in s */

int countLeadingNondigits(const char *s);
/* Count number of leading non-digit characters in s. */

int countSame(char *a, char *b);
/* Count number of characters that from start in a,b that are same. */

int chopString(char *in, char *sep, char *outArray[], int outSize);
/* int chopString(in, sep, outArray, outSize); */
/* This chops up the input string (cannabilizing it)
 * into an array of zero terminated strings in
 * outArray.  It returns the number of strings.
 * If you pass in NULL for outArray, it will just
 * return the number of strings that it *would*
 * chop. */

extern char crLfChopper[];
extern char whiteSpaceChopper[];
/* Some handy predefined separators. */

int chopByWhite(char *in, char *outArray[], int outSize);
/* Like chopString, but specialized for white space separators. */

#define chopLine(line, words) chopByWhite(line, words, ArraySize(words))
/* Chop line by white space. */

int chopByWhiteRespectDoubleQuotes(char *in, char *outArray[], int outSize);
/* Like chopString, but specialized for white space separators.
 * Further, any doubleQuotes (") are respected.
 * If doubleQuote encloses whole string, then they are removed:
 *   "Fred and Ethyl" results in word [Fred and Ethyl]
 * If doubleQuotes exist inside string they are retained:
 *   Fred "and Ethyl" results in word [Fred "and Ethyl"]
 * Special note "" is a valid, though empty word. */

int chopByChar(char *in, char chopper, char *outArray[], int outSize);
/* Chop based on a single character. */

#define chopTabs(string, words) chopByChar(string, '\t', words, ArraySize(words))
/* Chop string by tabs. */

#define chopCommas(string, words) chopByChar(string, ',', words, ArraySize(words))
/* Chop string by commas. */


char *skipBeyondDelimit(char *s,char delimit);
/* Returns NULL or pointer to first char beyond one (or more contiguous) delimit char.
   If delimit is ' ' then skips beyond first patch of whitespace. */

char *skipLeadingSpaces(char *s);
/* Return first non-white space */

char *skipToSpaces(char *s);
/* Return first white space. */

void eraseTrailingSpaces(char *s);
/* Replace trailing white space with zeroes. */

void eraseWhiteSpace(char *s);
/* Remove white space from a string */

void eraseNonAlphaNum(char *s);
/* Remove non-alphanumeric chars from string */

char *trimSpaces(char *s);
/* Remove leading and trailing white space. */

void repeatCharOut(FILE *f, char c, int count);
/* Write character to file repeatedly. */

void spaceOut(FILE *f, int count);
/* Put out some spaces to file. */

void starOut(FILE *f, int count);
/* Put out some asterisks to file. */

boolean hasWhiteSpace(char *s);
/* Return TRUE if there is white space in string. */

char *firstWordInLine(char *line);
/* Returns first word in line if any (white space separated).
 * Puts 0 in place of white space after word. */

char *lastWordInLine(char *line);
/* Returns last word in line if any (white space separated).
 * Returns NULL if string is empty.  Removes any terminating white space
 * from line. */

char *nextWord(char **pLine);
/* Return next word in *pLine and advance *pLine to next
 * word. Returns NULL when no more words. */

char *nextWordRespectingQuotes(char **pLine);
// return next word but respects single or double quotes surrounding sets of words.

char *cloneFirstWord(char *line);
/* Clone first word in line */

char *nextTabWord(char **pLine);
/* Return next tab-separated word. */

char *cloneFirstWordByDelimiter(char *line,char delimit);
/* Returns a cloned first word, not harming the memory passed in
   Delimiter of ' ' will delimit by isspace() */
#define cloneFirstWordInLine(line) cloneFirstWordByDelimiter((line),' ')
#define cloneFirstWordByTab(line)  cloneFirstWordByDelimiter((line),'\t')
/* Returns a cloned first word, not harming the memory passed in */

char *cloneNextWordByDelimiter(char **line,char delimit);
/* Returns a cloned first word, advancing the line pointer
   but not harming memory passed in. Delimit ' ' uses isspace() */
#define cloneNextWord(line)      cloneNextWordByDelimiter((line),' ')
#define cloneNextWordByTab(line) cloneNextWordByDelimiter((line),'\t')

char *nextStringInList(char **pStrings);
/* returns pointer to the first string and advances pointer to next in
   list of strings dilimited by 1 null and terminated by 2 nulls. */

int cntStringsInList(char *pStrings);
/* returns count of strings in a
   list of strings dilimited by 1 null and terminated by 2 nulls. */

int stringArrayIx(char *string, char *array[], int arraySize);
/* Return index of string in array or -1 if not there. */

int ptArrayIx(void *pt, void *array, int arraySize);
/* Return index of pt in array or -1 if not there. */

#define stringIx(string, array) stringArrayIx( (string), (array), ArraySize(array))

/* Some stuff that is left out of GNU .h files!? */
#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif

#ifndef SEEK_END
#define SEEK_END 2
#endif

#ifndef FILEPATH_H
void splitPath(char *path, char dir[PATH_LEN], char name[FILENAME_LEN],
	       char extension[FILEEXT_LEN]);
/* Split a full path into components.  The dir component will include the
 * trailing / if any.  The extension component will include the starting
 * . if any.   Pass in NULL for dir, name, or extension if you don't care about
 * that part. */
#endif /* FILEPATH_H */

char *addSuffix(char *head, char *suffix);
/* Return a needMem'd string containing "headsuffix". Should be free'd
 when finished. */

void chopSuffix(char *s);
/* Remove suffix (last . in string and beyond) if any. */

void chopSuffixAt(char *s, char c);
/* Remove end of string from last occurrence of char c.
 * chopSuffixAt(s, '.') is equivalent to regular chopSuffix. */

char *chopPrefix(char *s);
/* This will replace the first '.' in a string with
 * 0, and return the character after this.  If there
 * is no '.' in the string this will just return the
 * unchanged s passed in. */

char *chopPrefixAt(char *s, char c);
/* Like chopPrefix, but can chop on any character, not just '.' */

FILE *mustOpen(char *fileName, char *mode);
/* Open a file - or squawk and die. */

void mustWrite(FILE *file, void *buf, size_t size);
/* Write to file or squawk and die. */

#define writeOne(file, var) mustWrite((file), &(var), sizeof(var))
/* Write out one variable to file. */

void mustRead(FILE *file, void *buf, size_t size);
/* Read size bytes from a file or squawk and die. */

#define mustReadOne(file, var) mustRead((file), &(var), sizeof(var))
/* Read one variable from file or die. */

void mustGetLine(FILE *file, char *buf, int charCount);
/* Read at most charCount-1 bytes from file, but stop after newline if one is
 * encountered.  The string in buf is '\0'-terminated.  (See man 3 fgets.)
 * Die if there is an error. */

int mustOpenFd(char *fileName, int flags);
/* Open a file descriptor (see man 2 open) or squawk and die. */

void mustReadFd(int fd, void *buf, size_t size);
/* Read size bytes from a file descriptor or squawk and die. */

void mustWriteFd(int fd, void *buf, size_t size);
/* Write size bytes to file descriptor fd or die.  (See man 2 write.) */

off_t mustLseek(int fd, off_t offset, int whence);
/* Seek to given offset, relative to whence (see man lseek) in file descriptor fd or errAbort.
 * Return final offset (e.g. if this is just an (fd, 0, SEEK_CUR) query for current position). */

void mustCloseFd(int *pFd);
/* Close file descriptor *pFd if >= 0, abort if there's an error, set *pFd = -1. */

#define writeOneFd(fd, var) mustWriteFd((fd), &(var), sizeof(var))
/* Write out one variable to file descriptor fd. */

#define readOne(file, var) (fread(&(var), sizeof(var), 1, (file)) == 1)
/* Read one variable from file. Returns FALSE if can't do it. */

#define readOneFd(fd, var) (read((fd), &(var), sizeof(var)) == sizeof(var))
/* Read one variable from file. Returns FALSE if can't do it. */

#define mustReadOneFd(fd, var) mustReadFd((fd), &(var), sizeof(var))
/* Read one variable from file or die. */

#define memReadOne(pPt, var) memRead((pPt), &(var), sizeof(var))
/* Read one variable from memory. */

void writeString(FILE *f, char *s);
/* Write a 255 or less character string to a file.
 * This will write the length of the string in the first
 * byte then the string itself. */

char *readString(FILE *f);
/* Read a string (written with writeString) into
 * memory.  freeMem the result when done. Returns
 * NULL at EOF. */

char *mustReadString(FILE *f);
/* Read a string.  Squawk and die at EOF or if any problem. */

boolean fastReadString(FILE *f, char buf[256]);
/* Read a string into buffer, which must be long enough
 * to hold it.  String is in 'writeString' format.
 * Returns FALSE at EOF. */

void msbFirstWriteBits64(FILE *f, bits64 x);
/* Write out 64 bit number in manner that is portable across architectures */

bits64 msbFirstReadBits64(FILE *f);
/* Write out 64 bit number in manner that is portable across architectures */

void carefulClose(FILE **pFile);
/* Close file if open and null out handle to it. */

boolean carefulCloseWarn(FILE **pFile);
/* Close file if open and null out handle to it.
 * Return FALSE and print a warning message if there
 * is a problem.*/

char *firstWordInFile(char *fileName, char *wordBuf, int wordBufSize);
/* Read the first word in file into wordBuf. */

struct fileOffsetSize
/* A piece of a file. */
   {
   struct fileOffsetSize *next;	/* Next in list. */
   bits64	offset;		/* Start offset of block. */
   bits64	size;		/* Size of block. */
   };

int fileOffsetSizeCmp(const void *va, const void *vb);
/* Help sort fileOffsetSize by offset. */

struct fileOffsetSize *fileOffsetSizeMerge(struct fileOffsetSize *inList);
/* Returns a new list which is inList transformed to have adjacent blocks
 * merged.  Best to use this with a sorted list. */

void fileOffsetSizeFindGap(struct fileOffsetSize *list,
	struct fileOffsetSize **pBeforeGap, struct fileOffsetSize **pAfterGap);
/* Starting at list, find all items that don't have a gap between them and the previous item.
 * Return at gap, or at end of list, returning pointers to the items before and after the gap. */

void mustSystem(char *cmd);
/* Execute cmd using "sh -c" or die.  (See man 3 system.) fail on errors */

int roundingScale(int a, int p, int q);
/* returns rounded a*p/q */

int intAbs(int a);
/* Return integer absolute value */

#define logBase2(x)(log(x)/log(2))
/* return log base two of number */

#define round(a) ((int)((a)+0.5))
/* Round floating point val to nearest integer. */

#define roundll(a) ((long long)((a)+0.5))
/* Round floating point val to nearest long long. */

#ifndef min
#define min(a,b) ( (a) < (b) ? (a) : (b) )
/* Return min of a and b. */
#endif

#ifndef max
#define max(a,b) ( (a) > (b) ? (a) : (b) )
/* Return max of a and b. */
#endif

int  rangeIntersection(int start1, int end1, int start2, int end2);
/* Return amount of bases two ranges intersect over, 0 or negative if no
 * intersection. */

int  positiveRangeIntersection(int start1, int end1, int start2, int end2);
/* Return amount of bases two ranges intersect over, 0 if no
 * intersection. */

INLINE void memRead(char **pPt, void *buf, int size)
/* Copy memory from *pPt to buf, and advance *pPt by size. */
{
memcpy(buf, *pPt, size);
*pPt += size;
}

INLINE void memWrite(char **pPt, void *buf, int size)
/* Copy memory from buf to *pPt and advance *pPt by size. */
{
memcpy(*pPt, buf, size);
*pPt += size;
}

#define memWriteOne(pPt, var) memWrite((pPt), &(var), sizeof(var))
/* Write out one variable to memory stream. */

INLINE void memWriteFloat(char **pPt, float val)
/* Write out floating point val to file.  Mostly to convert from double... */
{
memWriteOne(pPt, val);
}

bits64 byteSwap64(bits64 a);
/* Swap from intel to sparc order of a 64 bit quantity. */

bits64 readBits64(FILE *f, boolean isSwapped);
/* Read and optionally byte-swap 64 bit entity. */

bits64 fdReadBits64(int fd, boolean isSwapped);
/* Read and optionally byte-swap 64 bit entity. */

bits64 memReadBits64(char **pPt, boolean isSwapped);
/* Read and optionally byte-swap 64 bit entity from memory buffer pointed to by
 * *pPt, and advance *pPt past read area. */

bits32 byteSwap32(bits32 a);
/* Swap from intel to sparc order of a 32 bit quantity. */

bits32 readBits32(FILE *f, boolean isSwapped);
/* Read and optionally byte-swap 32 bit entity. */

bits32 fdReadBits32(int fd, boolean isSwapped);
/* Read and optionally byte-swap 32 bit entity. */

bits32 memReadBits32(char **pPt, boolean isSwapped);
/* Read and optionally byte-swap 32 bit entity from memory buffer pointed to by
 * *pPt, and advance *pPt past read area. */

bits16 byteSwap16(bits16 a);
/* Swap from intel to sparc order of a 16 bit quantity. */

bits16 readBits16(FILE *f, boolean isSwapped);
/* Read and optionally byte-swap 16 bit entity. */

bits16 fdReadBits16(int fd, boolean isSwapped);
/* Read and optionally byte-swap 16 bit entity. */

bits16 memReadBits16(char **pPt, boolean isSwapped);
/* Read and optionally byte-swap 32 bit entity from memory buffer pointed to by
 * *pPt, and advance *pPt past read area. */

double byteSwapDouble(double a);
/* Return byte-swapped version of a */

double readDouble(FILE *f, boolean isSwapped);
/* Read and optionally byte-swap double-precision floating point entity. */

double memReadDouble(char **pPt, boolean isSwapped);
/* Read and optionally byte-swap double-precision floating point entity
 * from memory buffer pointed to by *pPt, and advance *pPt past read area. */

float byteSwapFloat(float a);
/* Return byte-swapped version of a */

float readFloat(FILE *f, boolean isSwapped);
/* Read and optionally byte-swap single-precision floating point entity. */

float memReadFloat(char **pPt, boolean isSwapped);
/* Read and optionally byte-swap single-precision floating point entity
 * from memory buffer pointed to by *pPt, and advance *pPt past read area. */

void removeReturns(char* dest, char* src);
/* Removes the '\r' character from a string.
 * the source and destination strings can be the same,
 * if there are no threads */

int intExp(char *text);
/* Convert text to integer expression and evaluate.
 * Throws if it finds a non-number. */

double doubleExp(char *text);
/* Convert text to floating point expression and
 * evaluate. */

char* readLine(FILE* fh);
/* Read a line of any size into dynamic memory, return null on EOF */

off_t fileSize(char *fileName);
/* The size of a file. */

boolean fileExists(char *fileName);
/* Does a file exist? */

/*
 Friendly name for strstrNoCase
*/
char *containsStringNoCase(char *haystack, char *needle);

char *strstrNoCase(char *haystack, char *needle);
/* A case-insensitive strstr */

int vasafef(char* buffer, int bufSize, char *format, va_list args);
/* Format string to buffer, vsprintf style, only with buffer overflow
 * checking.  The resulting string is always terminated with zero byte. */

int safef(char* buffer, int bufSize, char *format, ...)
/* Format string to buffer, vsprintf style, only with buffer overflow
 * checking.  The resulting string is always terminated with zero byte. */
#ifdef __GNUC__
__attribute__((format(printf, 3, 4)))
#endif
;

void safecpy(char *buf, size_t bufSize, const char *src);
/* copy a string to a buffer, with bounds checking.*/

void safencpy(char *buf, size_t bufSize, const char *src, size_t n);
/* copy n characters from a string to a buffer, with bounds checking.
 * Unlike strncpy, always null terminates the result */

void safecat(char *buf, size_t bufSize, const char *src);
/* Append  a string to a buffer, with bounds checking.*/

void safencat(char *buf, size_t bufSize, const char *src, size_t n);
/* append n characters from a string to a buffer, with bounds checking. */

char *naForNull(char *s);
/* Return 'n/a' if s is NULL, otherwise s. */

char *naForEmpty(char *s);
/* Return n/a if s is "" or NULL, otherwise s. */

char *emptyForNull(char *s);
/* Return "" if s is NULL, otherwise s. */

char *nullIfAllSpace(char *s);
/* Return NULL if s is all spaces, otherwise s. */

char *trueFalseString(boolean b);
/* Return "true" or "false" */

void uglyTime(char *label, ...);
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */

/*	In case the development environment does not supply INFINITY	*/
#if !defined(INFINITY)
#define INFINITY (1.0/0.0)
#endif

void makeDirs(char* path);
/* make a directory, including parent directories */

char *skipNumeric(char *s);
/* Return first char of s that's not a digit */

char *skipToNumeric(char *s);
/* skip up to where numeric digits appear */

char *splitOffNonNumeric(char *s);
/* Split off non-numeric part, e.g. mm of mm8. Result should be freed when done */

char *splitOffNumber(char *db);
/* Split off number part, e.g. 8 of mm8. Result should be freed when done */

void vaDumpStack(char *format, va_list args);
/* debugging function to run the pstack program on the current process. In
 * prints a message, following by a new line, and then the stack track.  Just
 * prints errors to stderr rather than aborts. For debugging purposes
 * only.  */

void dumpStack(char *format, ...)
/* debugging function to run the pstack program on the current process. In
 * prints a message, following by a new line, and then the stack track.  Just
 * prints errors to stderr rather than aborts. For debugging purposes
 * only.  */
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))
#endif
;

// SETTING_ON set of macros are frequently used comparisons of string values for boolean questions.
// Notice the subtle difference between NOT_ON and IS_OFF.  NOT_ON could be NULL but IS_OFF must be explicitly set
#define SETTING_IS_ON(setting)    (setting && (sameWord(setting,"on") || sameWord(setting,"true") || sameWord(setting,"yes") || sameWord(setting,"enabled") || atoi(setting) != 0))
#define SETTING_NOT_ON(setting)   (!SETTING_IS_ON(setting))
#define SETTING_IS_OFF(setting)   (setting && (sameWord(setting,"off") || sameWord(setting,"false") || sameWord(setting,"no") || sameWord(setting,"disabled") || sameWord(setting,"0")))

// Standard bit mask macros
#define BITS_ADD(    flags,bits) ((flags) = ((flags) |  (bits)))
#define BITS_REMOVE( flags,bits) ((flags) = ((flags) & ~(bits)))
#define BITS_ARE_ON( flags,bits) (((flags) & (bits)) == (bits))
#define BITS_ARE_OFF(flags,bits) (((flags) & (bits)) == 0)

// It is sometimes useful to distinguish between 3 "boolean" states: TRUE, FALSE and UNKNOWN
enum enumBool
    {
    beUnknown=0,              // Not yet set
    ebYes=1,                  // already set to TRUE
    ebNo=-1                   // already set to FALSE
    };
#define SET_TO_YES(ebool) { (ebool) = ebYes; }
#define SET_TO_NO(ebool)  { (ebool) = ebNo; }
#define IS_YES(ebool)     ((ebool) == ebYes)
#define IS_NO(ebool)      ((ebool) == ebNo)
#define IS_KNOWN(ebool)   (IS_YES(ebool) || IS_NO(ebool))
#define IS_TRUE           IS_YES
#define IS_FALSE          IS_NO

time_t mktimeFromUtc (struct tm *t);
/* Return time_t for tm in UTC (GMT)
 * Useful for stuff like converting to time_t the
 * last-modified HTTP response header
 * which is always GMT. Returns -1 on failure of mktime */


time_t dateToSeconds(const char *date,const char*format);
// Convert a string date to time_t

boolean dateIsOld(const char *date,const char*format);
// Is this string date older than now?

char *dateAddTo(char *date,char *format,int addYears,int addMonths,int addDays);
/* Add years,months,days to a formatted date and returns the new date as a cloned string
*  format is a strptime/strftime format: %F = yyyy-mm-dd */

#endif /* COMMON_H */
