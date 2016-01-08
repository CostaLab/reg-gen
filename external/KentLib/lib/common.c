/* Commonly used routines in a wide range of applications.
 * Strings, singly-linked lists, and a little file i/o.
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "errabort.h"
#include "portable.h"
#include "linefile.h"
#include "hash.h"

static char const rcsid[] = "$Id: common.c,v 1.151 2010/06/02 19:06:41 tdreszer Exp $";

void *cloneMem(void *pt, size_t size)
/* Allocate a new buffer of given size, and copy pt to it. */
{
void *newPt = needLargeMem(size);
memcpy(newPt, pt, size);
return newPt;
}

static char *cloneStringZExt(const char *s, int size, int copySize)
/* Make a zero terminated copy of string in memory */
{
char *d = needMem(copySize+1);
copySize = min(size,copySize);
memcpy(d, s, copySize);
d[copySize] = 0;
return d;
}

char *cloneStringZ(const char *s, int size)
/* Make a zero terminated copy of string in memory */
{
return cloneStringZExt(s, strlen(s), size);
}

char *cloneString(const char *s)
/* Make copy of string in dynamic memory */
{
int size = 0;
if (s == NULL)
    return NULL;
size = strlen(s);
return cloneStringZExt(s, size, size);
}

char *cloneLongString(char *s)
/* Make clone of long string. */
{
size_t size = strlen(s);
return cloneMem(s, size+1);
}

char *catTwoStrings(char *a, char *b)
/* Allocate new string that is a concatenation of two strings. */
{
int aLen = strlen(a), bLen = strlen(b);
int len = aLen + bLen;
char *newBuf = needLargeMem(len+1);
memcpy(newBuf, a, aLen);
memcpy(newBuf+aLen, b, bLen);
newBuf[len] = 0;
return newBuf;
}

/* Reverse the order of the bytes. */
void reverseBytes(char *bytes, long length)
{
long halfLen = (length>>1);
char *end = bytes+length;
char c;
while (--halfLen >= 0)
    {
    c = *bytes;
    *bytes++ = *--end;
    *end = c;
    }
}

void reverseInts(int *a, int length)
/* Reverse the order of the integer array. */
{
int halfLen = (length>>1);
int *end = a+length;
int c;
while (--halfLen >= 0)
    {
    c = *a;
    *a++ = *--end;
    *end = c;
    }
}

void reverseUnsigned(unsigned *a, int length)
/* Reverse the order of the unsigned array. */
{
int halfLen = (length>>1);
unsigned *end = a+length;
unsigned c;
while (--halfLen >= 0)
    {
    c = *a;
    *a++ = *--end;
    *end = c;
    }
}

void reverseDoubles(double *a, int length)
/* Reverse the order of the double array. */
{
int halfLen = (length>>1);
double *end = a+length;
double c;
while (--halfLen >= 0)
    {
    c = *a;
    *a++ = *--end;
    *end = c;
    }
}

void reverseStrings(char **a, int length)
/* Reverse the order of the char* array. */
{
int halfLen = (length>>1);
char **end = a+length;
char *c;
while (--halfLen >= 0)
    {
    c = *a;
    *a++ = *--end;
    *end = c;
    }
}

/* Swap buffers a and b. */
void swapBytes(char *a, char *b, int length)
{
char c;
int i;

for (i=0; i<length; ++i)
    {
    c = a[i];
    a[i] = b[i];
    b[i] = c;
    }
}


/** List managing routines. */

/* Count up elements in list. */
int slCount(const void *list)
{
struct slList *pt = (struct slList *)list;
int len = 0;

while (pt != NULL)
    {
    len += 1;
    pt = pt->next;
    }
return len;
}

void *slElementFromIx(void *list, int ix)
/* Return the ix'th element in list.  Returns NULL
 * if no such element. */
{
struct slList *pt = (struct slList *)list;
int i;
for (i=0;i<ix;i++)
    {
    if (pt == NULL) return NULL;
    pt = pt->next;
    }
return pt;
}

int slIxFromElement(void *list, void *el)
/* Return index of el in list.  Returns -1 if not on list. */
{
struct slList *pt;
int ix = 0;

for (pt = list, ix=0; pt != NULL; pt = pt->next, ++ix)
    if (el == (void*)pt)
	return ix;
return -1;
}

void *slLastEl(void *list)
/* Returns last element in list or NULL if none. */
{
struct slList *next, *el;
if ((el = list) == NULL)
    return NULL;
while ((next = el->next) != NULL)
    el = next;
return el;
}

/* Add new node to tail of list.
 * Usage:
 *    slAddTail(&list, node);
 * where list and nodes are both pointers to structure
 * that begin with a next pointer.
 */
void slAddTail(void *listPt, void *node)
{
struct slList **ppt = (struct slList **)listPt;
struct slList *n = (struct slList *)node;

while (*ppt != NULL)
    {
    ppt = &((*ppt)->next);
    }
n->next = NULL;
*ppt = n;
}

void *slPopHead(void *vListPt)
/* Return head of list and remove it from list. (Fast) */
{
struct slList **listPt = (struct slList **)vListPt;
struct slList *el = *listPt;
if (el != NULL)
    {
    *listPt = el->next;
    el->next = NULL;
    }
return el;
}

void *slPopTail(void *vListPt)
/* Return tail of list and remove it from list. (Not so fast) */
{
struct slList **listPt = (struct slList **)vListPt;
struct slList *el = *listPt;
if (el != NULL)
    {
    for (;;)
        {
        if (el->next == NULL)
            {
            *listPt = NULL;
            break;
            }
        listPt = &el->next;
        el = el->next;
        }
    }
return el;
}



void *slCat(void *va, void *vb)
/* Return concatenation of lists a and b.
 * Example Usage:
 *   struct slName *a = getNames("a");
 *   struct slName *b = getNames("b");
 *   struct slName *ab = slCat(a,b)
 */
{
struct slList *a = va;
struct slList *b = vb;
struct slList *end;
if (a == NULL)
    return b;
for (end = a; end->next != NULL; end = end->next)
    ;
end->next = b;
return a;
}

void slReverse(void *listPt)
/* Reverse order of a list.
 * Usage:
 *    slReverse(&list);
 */
{
struct slList **ppt = (struct slList **)listPt;
struct slList *newList = NULL;
struct slList *el, *next;

next = *ppt;
while (next != NULL)
    {
    el = next;
    next = el->next;
    el->next = newList;
    newList = el;
    }
*ppt = newList;
}

void slFreeList(void *listPt)
/* Free list */
{
struct slList **ppt = (struct slList**)listPt;
struct slList *next = *ppt;
struct slList *el;

while (next != NULL)
    {
    el = next;
    next = el->next;
    freeMem((char*)el);
    }
*ppt = NULL;
}

void slSort(void *pList, int (*compare )(const void *elem1,  const void *elem2))
/* Sort a singly linked list with Qsort and a temporary array. */
{
struct slList **pL = (struct slList **)pList;
struct slList *list = *pL;
int count;
count = slCount(list);
if (count > 1)
    {
    struct slList *el;
    struct slList **array;
    int i;
    array = needLargeMem(count * sizeof(*array));
    for (el = list, i=0; el != NULL; el = el->next, i++)
        array[i] = el;
    qsort(array, count, sizeof(array[0]), compare);
    list = NULL;
    for (i=0; i<count; ++i)
        {
        array[i]->next = list;
        list = array[i];
        }
    freeMem(array);
    slReverse(&list);
    *pL = list;
    }
}

void slUniqify(void *pList, int (*compare )(const void *elem1,  const void *elem2), void (*free)())
/* Return sorted list with duplicates removed.
 * Compare should be same type of function as slSort's compare (taking
 * pointers to pointers to elements.  Free should take a simple
 * pointer to dispose of duplicate element, and can be NULL. */
{
struct slList **pSlList = (struct slList **)pList;
struct slList *oldList = *pSlList;
struct slList *newList = NULL, *el;

slSort(&oldList, compare);
while ((el = slPopHead(&oldList)) != NULL)
    {
    if ((newList == NULL) || (compare(&newList, &el) != 0))
        slAddHead(&newList, el);
    else if (free != NULL)
        free(el);
    }
slReverse(&newList);
*pSlList = newList;
}

boolean slRemoveEl(void *vpList, void *vToRemove)
/* Remove element from doubly linked list.  Usage:
 *    slRemove(&list, el);
 * Returns TRUE if element in list.  */
{
struct slList **pList = vpList;
struct slList *toRemove = vToRemove;
struct slList *el, *next, *newList = NULL;
boolean didRemove = FALSE;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    if (el != toRemove)
	{
	slAddHead(&newList, el);
	}
    else
        didRemove = TRUE;
    }
slReverse(&newList);
*pList = newList;
return didRemove;
}

struct slInt *slIntNew(int x)
/* Return a new int. */
{
struct slInt *a;
AllocVar(a);
a->val = x;
return a;
}

int slIntCmp(const void *va, const void *vb)
/* Compare two slInts. */
{
const struct slInt *a = *((struct slInt **)va);
const struct slInt *b = *((struct slInt **)vb);
return a->val - b->val;
}

int slIntCmpRev(const void *va, const void *vb)
/* Compare two slInts in reverse direction. */
{
const struct slInt *a = *((struct slInt **)va);
const struct slInt *b = *((struct slInt **)vb);
return b->val - a->val;
}

struct slInt * slIntFind(struct slInt *list, int target)
/* Find target in slInt list or return NULL */
{
struct slInt *i;
for (i=list;i;i=i->next)
    if (i->val == target)
	return i;
return NULL;
}

static int doubleCmp(const void *va, const void *vb)
/* Compare function to sort array of doubles. */
{
const double *a = va;
const double *b = vb;
double diff = *a - *b;
if (diff < 0)
    return -1;
else if (diff > 0)
    return 1;
else
    return 0;
}

void doubleSort(int count, double *array)
/* Sort an array of doubles. */
{
if (count > 1)
qsort(array, count, sizeof(array[0]), doubleCmp);
}

double doubleMedian(int count, double *array)
/* Return median value in array.  This will sort
 * the array as a side effect. */
{
double median;
doubleSort(count, array);
if ((count&1) == 1)
    median = array[count>>1];
else
    {
    count >>= 1;
    median = (array[count] + array[count-1]) * 0.5;
    }
return median;
}

void doubleBoxWhiskerCalc(int count, double *array, double *retMin,
	double *retQ1, double *retMedian, double *retQ3, double *retMax)
/* Calculate what you need to draw a box and whiskers plot from an array of doubles. */
{
doubleSort(count, array);
*retMin = array[0];
*retQ1 = array[(count+2)/4];
int halfCount = count>>1;
if ((count&1) == 1)
    *retMedian = array[halfCount];
else
    {
    *retMedian = (array[halfCount] + array[halfCount-1]) * 0.5;
    }
*retQ3 = array[(3*count+2)/4];
*retMax = array[count-1];
}

struct slDouble *slDoubleNew(double x)
/* Return a new double. */
{
struct slDouble *a;
AllocVar(a);
a->val = x;
return a;
}

int slDoubleCmp(const void *va, const void *vb)
/* Compare two slDoubles. */
{
const struct slDouble *a = *((struct slDouble **)va);
const struct slDouble *b = *((struct slDouble **)vb);
double diff = a->val - b->val;
if (diff < 0)
    return -1;
else if (diff > 0)
    return 1;
else
    return 0;
}

double slDoubleMedian(struct slDouble *list)
/* Return median value on list. */
{
int i,count = slCount(list);
struct slDouble *el;
double *array, median;
if (count == 0)
    errAbort("Can't take median of empty list");
AllocArray(array,count);
for (i=0, el=list; i<count; ++i, el=el->next)
    array[i] = el->val;
median = doubleMedian(count, array);
freeMem(array);
return median;
}

void slDoubleBoxWhiskerCalc(struct slDouble *list, double *retMin,
	double *retQ1, double *retMedian, double *retQ3, double *retMax)
/* Calculate what you need to draw a box and whiskers plot from a list of slDoubles. */
{
int i,count = slCount(list);
struct slDouble *el;
double *array;
if (count == 0)
    errAbort("Can't take do slDoubleBoxWhiskerCalc of empty list");
AllocArray(array,count);
for (i=0, el=list; i<count; ++i, el=el->next)
    array[i] = el->val;
doubleBoxWhiskerCalc(count, array, retMin, retQ1, retMedian, retQ3, retMax);
freeMem(array);
}

static int intCmp(const void *va, const void *vb)
/* Compare function to sort array of ints. */
{
const int *a = va;
const int *b = vb;
int diff = *a - *b;
if (diff < 0)
    return -1;
else if (diff > 0)
    return 1;
else
    return 0;
}

void intSort(int count, int *array)
/* Sort an array of ints. */
{
if (count > 1)
qsort(array, count, sizeof(array[0]), intCmp);
}

int intMedian(int count, int *array)
/* Return median value in array.  This will sort
 * the array as a side effect. */
{
int median;
intSort(count, array);
if ((count&1) == 1)
    median = array[count>>1];
else
    {
    count >>= 1;
    median = (array[count] + array[count-1]) * 0.5;
    }
return median;
}


struct slName *newSlName(char *name)
/* Return a new name. */
{
struct slName *sn;
if (name != NULL)
    {
    int len = strlen(name);
    sn = needMem(sizeof(*sn)+len);
    strcpy(sn->name, name);
    return sn;
    }
else
    {
    AllocVar(sn);
    }
return sn;
}

struct slName *slNameNewN(char *name, int size)
/* Return new slName of given size. */
{
struct slName *sn = needMem(sizeof(*sn) + size);
memcpy(sn->name, name, size);
return sn;
}

int slNameCmpCase(const void *va, const void *vb)
/* Compare two slNames, ignore case. */
{
const struct slName *a = *((struct slName **)va);
const struct slName *b = *((struct slName **)vb);
return strcasecmp(a->name, b->name);
}

void slNameSortCase(struct slName **pList)
/* Sort slName list, ignore case. */
{
slSort(pList, slNameCmpCase);
}

int slNameCmp(const void *va, const void *vb)
/* Compare two slNames. */
{
const struct slName *a = *((struct slName **)va);
const struct slName *b = *((struct slName **)vb);
return strcmp(a->name, b->name);
}

int slNameCmpStringsWithEmbeddedNumbers(const void *va, const void *vb)
/* Compare strings such as gene names that may have embedded numbers,
 * so that bmp4a comes before bmp14a */
{
const struct slName *a = *((struct slName **)va);
const struct slName *b = *((struct slName **)vb);
return cmpStringsWithEmbeddedNumbers(a->name, b->name);
}



void slNameSort(struct slName **pList)
/* Sort slName list. */
{
slSort(pList, slNameCmp);
}

boolean slNameInList(struct slName *list, char *string)
/* Return true if string is in name list -- case insensitive. */
{
struct slName *el;
for (el = list; el != NULL; el = el->next)
    if (sameWord(string, el->name))
        return TRUE;
return FALSE;
}

boolean slNameInListUseCase(struct slName *list, char *string)
/* Return true if string is in name list -- case sensitive. */
{
struct slName *el;
for (el = list; el != NULL; el = el->next)
    if (string != NULL && !strcmp(string, el->name))
        return TRUE;
return FALSE;
}

void *slNameFind(void *list, char *string)
/* Return first element of slName list (or any other list starting
 * with next/name fields) that matches string. */
{
struct slName *el;
for (el = list; el != NULL; el = el->next)
    if (sameWord(string, el->name))
        return el;
return NULL;
}

int slNameFindIx(struct slName *list, char *string)
/* Return index of first element of slName list (or any other
 * list starting with next/name fields) that matches string.
 * Return -1 if not found. */
{
struct slName *el;
int ix = 0;
for (el = list; el != NULL; el = el->next, ix++)
    if (sameString(string, el->name))
        return ix;
return -1;
}

char *slNameStore(struct slName **pList, char *string)
/* Put string into list if it's not there already.
 * Return the version of string stored in list. */
{
struct slName *el;
for (el = *pList; el != NULL; el = el->next)
    {
    if (sameString(string, el->name))
	return el->name;
    }
el = newSlName(string);
slAddHead(pList, el);
return el->name;
}

struct slName *slNameAddHead(struct slName **pList, char *name)
/* Add name to start of list and return it. */
{
struct slName *el = slNameNew(name);
slAddHead(pList, el);
return el;
}

struct slName *slNameAddTail(struct slName **pList, char *name)
/* Add name to end of list (not efficient for long lists),
 * and return it. */
{
struct slName *el = slNameNew(name);
slAddTail(pList, el);
return el;
}

struct slName *slNameCloneList(struct slName *list)
/* Return clone of list. */
{
struct slName *el, *newEl, *newList = NULL;
for (el = list; el != NULL; el = el->next)
    {
    newEl = slNameNew(el->name);
    slAddHead(&newList, newEl);
    }
slReverse(&newList);
return newList;
}


struct slName *slNameListFromString(char *s, char delimiter)
/* Return list of slNames gotten from parsing delimited string.
 * The final delimiter is optional. a,b,c  and a,b,c, are equivalent
 * for comma-delimited lists. */
{
char *e;
struct slName *list = NULL, *el;
while (s != NULL && s[0] != 0)
    {
    e = strchr(s, delimiter);
    if (e == NULL)
	el = slNameNew(s);
    else
	{
        el = slNameNewN(s, e-s);
	e += 1;
	}
    slAddHead(&list, el);
    s = e;
    }
slReverse(&list);
return list;
}

struct slName *slNameListOfUniqueWords(char *text,boolean respectQuotes)
// Return list of unique words found by parsing string delimited by whitespace.
// If respectQuotes then ["Lucy and Ricky" 'Fred and Ethyl'] will yield 2 slNames no quotes
{
struct slName *list = NULL;
char *word = NULL;
while (text != NULL)
    {
    if (respectQuotes)
        {
        word = nextWordRespectingQuotes(&text);
        if (word != NULL)
            {
            if (word[0] == '"')
                stripChar(word, '"');
            else if (word[0] == '\'')
                stripChar(word, '\'');
            }
        }
    else
        word = nextWord(&text);
    if (word)
        slNameStore(&list, word);
    else
        break;
    }

slReverse(&list);
return list;
}

struct slName *slNameListFromStringArray(char *stringArray[], int arraySize)
/* Return list of slNames from an array of strings of length arraySize.
 * If a string in the array is NULL, the array will be treated as
 * NULL-terminated (shorter than arraySize). */
{
char *s;
struct slName *list = NULL, *el;
int i;
if (stringArray == NULL)
    return NULL;
for (i = 0;  i < arraySize;  i++)
    {
    s = stringArray[i];
    if (s == NULL)
	break;
    el = slNameNew(s);
    slAddHead(&list, el);
    }
slReverse(&list);
return list;
}

char *slNameListToString(struct slName *list, char delimiter)
/* Return string created by joining all names with the delimiter. */
{
struct slName *el;
int elCount = 0;
int len = 0;
char del[2];
char *s;

del[0] = delimiter;
del[1] = '\0';

for (el = list; el != NULL; el = el->next, elCount++)
	len += strlen(el->name);
len += elCount;

AllocArray(s, len);

for (el = list; el != NULL; el = el->next)
	{
	strcat(s, el->name);
	if (el->next != NULL)
		strcat(s, del);
	}
return s;
}

struct slName *slNameLoadReal(char *fileName)
/* load file lines that are not blank or start with a '#' into a slName
 * list */
{
struct slName *lines = NULL;
char *line;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
while (lineFileNextReal(lf, &line))
    slSafeAddHead(&lines, slNameNew(line));
lineFileClose(&lf);
slReverse(&lines);
return lines;
}

struct slName *slNameIntersection(struct slName *a, struct slName *b)
/* return intersection of two slName lists.  */
{
struct hash *hashA = newHash(0);
struct slName *el, *retval = NULL;

for (el = a; el != NULL; el = el->next)
    hashAddInt(hashA, el->name, 1);
for (el = b; el != NULL; el = el->next)
    if(hashLookup(hashA, el->name) != NULL)
        slNameAddHead(&retval, el->name);
hashFree(&hashA);
return retval;
}

struct slRef *refOnList(struct slRef *refList, void *val)
/* Return ref if val is already on list, otherwise NULL. */
{
struct slRef *ref;
for (ref = refList; ref != NULL; ref = ref->next)
    if (ref->val == val)
        return ref;
return NULL;
}

struct slRef *slRefNew(void *val)
/* Create new slRef element. */
{
struct slRef *ref;
AllocVar(ref);
ref->val = val;
return ref;
}

void refAdd(struct slRef **pRefList, void *val)
/* Add reference to list. */
{
struct slRef *ref;
AllocVar(ref);
ref->val = val;
slAddHead(pRefList, ref);
}

void refAddUnique(struct slRef **pRefList, void *val)
/* Add reference to list if not already on list. */
{
if (refOnList(*pRefList, val) == NULL)
    {
    refAdd(pRefList, val);
    }
}

struct slRef *refListFromSlList(void *list)
/* Make a reference list that mirrors a singly-linked list. */
{
struct slList *el;
struct slRef *refList = NULL, *ref;
for (el= list; el != NULL; el = el->next)
    {
    ref = slRefNew(el);
    slAddHead(&refList, ref);
    }
slReverse(&refList);
return refList;
}


struct slPair *slPairNew(char *name, void *val)
/* Allocate new name/value pair. */
{
struct slPair *el;
AllocVar(el);
el->name = cloneString(name);
el->val = val;
return el;
}

void slPairAdd(struct slPair **pList, char *name, void *val)
/* Add new slPair to head of list. */
{
struct slPair *el = slPairNew(name, val);
slAddHead(pList, el);
}

void slPairFree(struct slPair **pEl)
/* Free up struct and name.  (Don't free up values.) */
{
struct slPair *el = *pEl;
if (el != NULL)
    {
    freeMem(el->name);
    freez(pEl);
    }
}

void slPairFreeList(struct slPair **pList)
/* Free up list.  (Don't free up values.) */
{
struct slPair *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    slPairFree(&el);
    }
*pList = NULL;
}

void slPairFreeVals(struct slPair *list)
/* Free up all values on list. */
{
struct slPair *el;
for (el = list; el != NULL; el = el->next)
    freez(&el->val);
}

void slPairFreeValsAndList(struct slPair **pList)
/* Free up all values on list and list itself */
{
slPairFreeVals(*pList);
slPairFreeList(pList);
}

struct slPair *slPairFind(struct slPair *list, char *name)
/* Return list element of given name, or NULL if not found. */
{
struct slPair *el;
for (el = list; el != NULL; el = el->next)
    if (sameString(name, el->name))
        break;
return el;
}

void *slPairFindVal(struct slPair *list, char *name)
/* Return value associated with name in list, or NULL if not found. */
{
struct slPair *el = slPairFind(list, name);
if (el == NULL)
    return NULL;
return el->val;
}

struct slPair *slPairListFromString(char *str,boolean respectQuotes)
// Return slPair list parsed from list in string like:  [name1=val1 name2=val2 ...]
// if respectQuotes then string can have double quotes: [name1="val 1" "name 2"=val2 ...]
//    resulting pair strips quotes: {name1}={val 1},{name 2}={val2}
// Returns NULL if parse error.  Free this up with slPairFreeValsAndList.
{
char *s = skipLeadingSpaces(str);  // Would like to remove this and tighten up the standard someday.
if (isEmpty(s))
    return NULL;

struct slPair *list = NULL;
char name[1024];
char val[1024];
char buf[1024];
bool inQuote = FALSE;
char *b = buf;
char sep = '=';
char c = ' ';
int mode = 0;
while(1)
    {
    c = *s++;
    if (mode == 0 || mode == 2) // reading name or val
	{
	boolean term = FALSE;
	if (respectQuotes && b == buf && !inQuote && c == '"')
	    inQuote = TRUE;
	else if (inQuote && c == '"')
	    term = TRUE;
	else if ((c == sep || c == 0) && !inQuote)
	    {
	    term = TRUE;
	    --s;  // rewind
	    }
	else if (c == ' ' && !inQuote)
	    {
	    warn("slPairListFromString: Unexpected whitespace in %s", str);
	    return NULL;
	    }
	else if (c == 0 && inQuote)
	    {
	    warn("slPairListFromString: Unterminated quote in %s", str);
	    return NULL;
	    }
	else
	    {
	    *b++ = c;
	    if ((b - buf) > sizeof buf)
		{
		warn("slPairListFromString: pair name or value too long in %s", str);
		return NULL;
		}
	    }
	if (term)
	    {
	    inQuote = FALSE;
	    *b = 0;
	    if (mode == 0)
		{
		safecpy(name, sizeof name, buf);
		if (strlen(name)<1)
		    {
		    warn("slPairListFromString: Pair name cannot be empty in %s", str);
		    return NULL;
		    }
		// Shall we check for name being alphanumeric, at least for the respectQuotes=FALSE case?
		}
	    else // mode == 2
		{
		safecpy(val, sizeof val, buf);
		if (!respectQuotes && (hasWhiteSpace(name) || hasWhiteSpace(val))) // should never happen
		    {
		    warn("slPairListFromString() Unexpected white space in name=value pair: [%s]=[%s] in string=[%s]\n", name, val, str);
		    break;
		    }
		slPairAdd(&list, name, cloneString(val));
		}
	    ++mode;
	    }
	}
    else if (mode == 1) // read required "=" sign
	{
	if (c != '=')
	    {
	    warn("slPairListFromString: Expected character = after name in %s", str);
	    return NULL;
	    }
	++mode;
	sep = ' ';
	b = buf;
	}
    else // (mode == 3) reading optional separating space
	{
	if (c == 0)
	    break;
	if (c != ' ')
	    {
	    mode = 0;
	    --s;
	    b = buf;
	    sep = '=';
	    }
	}
    }
slReverse(&list);
return list;
}

char *slPairListToString(struct slPair *list,boolean quoteIfSpaces)
// Returns an allocated string of pairs in form of [name1=val1 name2=val2 ...]
// If requested, will wrap name or val in quotes if contain spaces: [name1="val 1" "name 2"=val2]
{
// Don't rely on dyString.  We should do the accounting ourselves and not create extra dependencies.
int count = 0;
struct slPair *pair = list;
for(;pair != NULL; pair = pair->next)
    {
    assert(pair->name != NULL && pair->val != NULL); // Better assert and get this over with, complete with stack
    count += strlen(pair->name);
    count += strlen((char *)(pair->val));
    count += 2; // = and ' ' delimit
    if (quoteIfSpaces)
        {
        if (hasWhiteSpace(pair->name))
            count += 2; // " and "
        if (hasWhiteSpace((char *)(pair->val)))
            count += 2; // " and "
        }
    }
if (count == 0)
    return NULL;

char *str = needMem(count+5); // A bit of slop

char *strPtr = str;
for(pair = list; pair != NULL; pair = pair->next, strPtr += strlen(strPtr))
    {
    if (pair != list) // Not first cycle
        *strPtr++ = ' ';
    if (hasWhiteSpace(pair->name))
        {
        if (quoteIfSpaces)
            sprintf(strPtr,"\"%s\"=",pair->name);
        else
            {
            warn("slPairListToString() Unexpected white space in name: [%s]\n", pair->name);
            sprintf(strPtr,"%s=",pair->name); // warn but still make string
            }
        }
    else
        sprintf(strPtr,"%s=",pair->name);
    strPtr += strlen(strPtr);
    if (hasWhiteSpace((char *)(pair->val)))
        {
        if (quoteIfSpaces)
            sprintf(strPtr,"\"%s\"",(char *)(pair->val));
        else
            {
            warn("slPairListToString() Unexpected white space in val: [%s]\n", (char *)(pair->val));
            sprintf(strPtr,"%s",(char *)(pair->val)); // warn but still make string
            }
        }
    else
        sprintf(strPtr,"%s",(char *)(pair->val));
    }
return str;
}

char *slPairNameToString(struct slPair *list, char delimiter,boolean quoteIfSpaces)
// Return string created by joining all names (ignoring vals) with the delimiter.
// If requested, will wrap name in quotes if contain spaces: [name1,"name 2" ...]
{
int elCount = 0;
int count = 0;
struct slPair *pair = list;
for (; pair != NULL; pair = pair->next, elCount++)
    {
    assert(pair->name != NULL);
    count += strlen(pair->name);
    if (quoteIfSpaces && hasWhiteSpace(pair->name))
        count += 2;
    }
count += elCount;
if (count == 0)
    return NULL;

char *str = needMem(count+5); // A bit of slop

char *strPtr = str;
for(pair = list; pair != NULL; pair = pair->next, strPtr += strlen(strPtr))
    {
    if (pair != list)
        *strPtr++ = delimiter;
    if (hasWhiteSpace(pair->name))
        {
        if (quoteIfSpaces)
            sprintf(strPtr,"\"%s\"",pair->name);
        else
            {
            if (delimiter == ' ')  // if delimied by commas, this is entirely okay!
                warn("slPairListToString() Unexpected white space in name delimied by space: [%s]\n", pair->name);
            sprintf(strPtr,"%s",pair->name); // warn but still make string
            }
        }
    else
        sprintf(strPtr,"%s",pair->name);
    }
return str;
}

int slPairCmpCase(const void *va, const void *vb)
/* Compare two slPairs, ignore case. */
{
const struct slPair *a = *((struct slPair **)va);
const struct slPair *b = *((struct slPair **)vb);
return strcasecmp(a->name, b->name);
}

void slPairSortCase(struct slPair **pList)
/* Sort slPair list, ignore case. */
{
slSort(pList, slPairCmpCase);
}

int slPairCmp(const void *va, const void *vb)
/* Compare two slPairs. */
{
const struct slPair *a = *((struct slPair **)va);
const struct slPair *b = *((struct slPair **)vb);
return strcmp(a->name, b->name);
}


int slPairValCmpCase(const void *va, const void *vb)
/* Case insensitive compare two slPairs on their values (must be string). */
{
const struct slPair *a = *((struct slPair **)va);
const struct slPair *b = *((struct slPair **)vb);
return strcasecmp((char *)(a->val), (char *)(b->val));
}

int slPairValCmp(const void *va, const void *vb)
/* Compare two slPairs on their values (must be string). */
{
const struct slPair *a = *((struct slPair **)va);
const struct slPair *b = *((struct slPair **)vb);
return strcmp((char *)(a->val), (char *)(b->val));
}

void slPairValSortCase(struct slPair **pList)
/* Sort slPair list on values (must be string), ignore case. */
{
slSort(pList, slPairValCmpCase);
}

void slPairValSort(struct slPair **pList)
/* Sort slPair list on values (must be string). */
{
slSort(pList, slPairValCmp);
}

int slPairIntCmp(const void *va, const void *vb)
// Compare two slPairs on their integer values.
{
const struct slPair *a = *((struct slPair **)va);
const struct slPair *b = *((struct slPair **)vb);
return ((char *)(a->val) - (char *)(b->val)); // cast works and val is 0 vased integer
}

void slPairIntSort(struct slPair **pList)
// Sort slPair list on integer values.
{
slSort(pList, slPairIntCmp);
}

int slPairAtoiCmp(const void *va, const void *vb)
// Compare two slPairs on their strings interpreted as integer values.
{
const struct slPair *a = *((struct slPair **)va);
const struct slPair *b = *((struct slPair **)vb);
return (atoi((char *)(a->val)) - atoi((char *)(b->val)));
}

void slPairValAtoiSort(struct slPair **pList)
// Sort slPair list on string values interpreted as integers.
{
slSort(pList, slPairAtoiCmp);
}

void gentleFree(void *pt)
{
if (pt != NULL) freeMem((char*)pt);
}

int differentWord(char *s1, char *s2)
/* strcmp ignoring case - returns zero if strings are
 * the same (ignoring case) otherwise returns difference
 * between first non-matching characters. */
{
char c1, c2;
for (;;)
    {
    c1 = toupper(*s1++);
    c2 = toupper(*s2++);
    if (c1 != c2) /* Takes care of end of string in one but not the other too */
	return c2-c1;
    if (c1 == 0)  /* Take care of end of string in both. */
	return 0;
    }
}

int differentStringNullOk(char *a, char *b)
/* Returns 0 if two strings (either of which may be NULL)
 * are the same.  Otherwise it returns a positive or negative
 * number depending on the alphabetical order of the two
 * strings.
 * This is basically a strcmp that can handle NULLs in
 * the input.  If used in a sort the NULLs will end
 * up before any of the cases with data.   */
{
if (a == b)
    return FALSE;
else if (a == NULL)
    return -1;
else if (b == NULL)
    return 1;
else
    return strcmp(a,b) != 0;
}

boolean startsWith(const char *start, const char *string)
/* Returns TRUE if string begins with start. */
{
char c;
int i;

for (i=0; ;i += 1)
    {
    if ((c = start[i]) == 0)
        return TRUE;
    if (string[i] != c)
        return FALSE;
    }
}

boolean startsWithWord(char *firstWord, char *line)
/* Return TRUE if first white-space-delimited word in line
 * is same as firstWord.  Comparison is case sensitive. */
{
int len = strlen(firstWord);
int i;
for (i=0; i<len; ++i)
   if (firstWord[i] != line[i])
       return FALSE;
char c = line[len];
return c == 0 || isspace(c);
}

boolean startsWithWordByDelimiter(char *firstWord,char delimit, char *line)
/* Return TRUE if first word in line is same as firstWord as delimited by delimit.
   Comparison is case sensitive. Delimit of ' ' uses isspace() */
{
if(delimit == ' ')
    return startsWithWord(firstWord,line);
if (!startsWith(firstWord,line))
    return FALSE;
char c = line[strlen(firstWord)];
return (c == '\0' || c == delimit);
}

char * findWordByDelimiter(char *word,char delimit, char *line)
/* Return pointer to first occurance of word in line broken by 'delimit' char
   Comparison is case sensitive. Delimit of ' ' uses isspace() */
{
int ix;
char *p=line;
while(p!=NULL && *p!='\0')
    {
    for (ix = 0;
         word[ix] != '\0' && word[ix] == *p;
         ix++,p++); // advance as long as they match
    if(ix == strlen(word))
        {
        if(*p=='\0'
        || *p==delimit
        || (delimit == ' ' && isspace(*p)))
            return p - ix; // matched and delimited
        }
    for(;   *p!='\0'
         && *p!=delimit
         && (delimit != ' ' || !isspace(*p));
            p++); // advance to next delimit
    if(*p!='\0')
        {
        p++;
        continue;  // delimited so start again after delimit
        }
    }
return NULL;
}

char *rStringIn(char *needle, char *haystack)
/* Return last position of needle in haystack, or NULL if it's not there. */
{
int nSize = strlen(needle);
char *pos;
for (pos = haystack + strlen(haystack) - nSize; pos >= haystack; pos -= 1)
    {
    if (memcmp(needle, pos, nSize) == 0)
        return pos;
    }
return NULL;
}

char *stringBetween(char *start, char *end, char *haystack)
/* Return string between start and end strings, or NULL if
 * none found.  The first such instance is returned.
 * String must be freed by caller. */
{
char *pos, *p;
int len;
if ((p = stringIn(start, haystack)) != NULL)
    {
    pos = p + strlen(start);
    if ((p = stringIn(end, pos)) != NULL)
        {
        len = p - pos;
        pos = cloneMem(pos, len + 1);
        pos[len] = 0;
        return pos;
        }
    }
return NULL;
}

boolean endsWith(char *string, char *end)
/* Returns TRUE if string ends with end. */
{
int sLen, eLen, offset;
sLen = strlen(string);
eLen = strlen(end);
offset = sLen - eLen;
if (offset < 0)
    return FALSE;
return sameString(string+offset, end);
}

char lastChar(char *s)
/* Return last character in string. */
{
if (s == NULL || s[0] == 0)
    return 0;
return s[strlen(s)-1];
}

char *lastNonwhitespaceChar(char *s)
// Return pointer to last character in string that is not whitespace.
{
if (s == NULL || s[0] == 0)
    return NULL;

char *sPos = s + (strlen(s) - 1);
for (;sPos >= s;sPos--)
    {
    if (!isspace(*sPos))
        return sPos;
    }
return NULL;
}

char *matchingCharBeforeInLimits(char *limit, char *s, char c)
/* Look for character c sometime before s, but going no further than limit.
 * Return NULL if not found. */
{
while (--s >= limit)
    if (*s == c)
        return s;
return NULL;
}

char *memMatch(char *needle, int nLen, char *haystack, int hLen)
/* Returns first place where needle (of nLen chars) matches
 * haystack (of hLen chars) */
{
char c = *needle++;
nLen -= 1;
hLen -= nLen;
while (--hLen >= 0)
    {
    if (*haystack++ == c && memcmp(needle, haystack, nLen) == 0)
        {
        return haystack-1;
        }
    }
return NULL;
}

void toUpperN(char *s, int n)
/* Convert a section of memory to upper case. */
{
int i;
for (i=0; i<n; ++i)
    s[i] = toupper(s[i]);
}

void toLowerN(char *s, int n)
/* Convert a section of memory to upper case. */
{
int i;
for (i=0; i<n; ++i)
    s[i] = tolower(s[i]);
}

void toggleCase(char *s, int size)
/* toggle upper and lower case chars in string. */
{
char c;
int i;
for (i=0; i<size; ++i)
    {
    c = s[i];
    if (isupper(c))
        c = tolower(c);
    else if (islower(c))
        c = toupper(c);
    s[i] = c;
    }
}


char *strUpper(char *s)
/* Convert entire string to upper case. */
{
char c;
char *ss=s;
for (;;)
    {
    if ((c = *ss) == 0) break;
    *ss++ = toupper(c);
    }
return s;
}

char *replaceChars(char *string, char *old, char *new)
/*
  Replaces the old with the new. The old and new string need not be of equal size
 Can take any length string.
 Return value needs to be freeMem'd.
*/
{
int numTimes = 0;
int oldLen = strlen(old);
int newLen = strlen(new);
int strLen = 0;
char *result = NULL;
char *ptr = strstr(string, old);
char *resultPtr = NULL;

while(NULL != ptr)
    {
    numTimes++;
    ptr += oldLen;
    ptr = strstr(ptr, old);
    }
strLen = max(strlen(string) + (numTimes * (newLen - oldLen)), strlen(string));
result = needMem(strLen + 1);

ptr = strstr(string, old);
resultPtr = result;
while(NULL != ptr)
    {
    strLen = ptr - string;
    strcpy(resultPtr, string);
    string = ptr + oldLen;

    resultPtr += strLen;
    strcpy(resultPtr, new);
    resultPtr += newLen;
    ptr = strstr(string, old);
    }

strcpy(resultPtr, string);
return result;
}

int strSwapStrs(char *string, int sz,char *oldStr, char *newStr)
/* Swaps all occurrences of the old with the new in string. Need not be same size
   Swaps in place but restricted by sz.  Returns count of swaps or -1 for sz failure. */
{
// WARNING: called at low level, so no errors allowed.
int count = 0;
char *p=NULL;
for(p=strstr(string,oldStr);p!=NULL;p=strstr(p+strlen(oldStr),oldStr))
    count++;
if(count == 0)
    return 0;
if((strlen(string)+(count*(strlen(newStr) - strlen(oldStr))))>=sz)
    return -1;
for(p=strstr(string,oldStr);p!=NULL;p=strstr(p+strlen(newStr),oldStr))
    {
    memmove(p+strlen(newStr),p+strlen(oldStr),strlen(p+strlen(oldStr))+1); // NULL at end is also moved!
    memcpy(p,newStr,strlen(newStr));
    }
return count;
}

char *strLower(char *s)
/* Convert entire string to lower case */
{
char c;
char *ss=s;
for (;;)
    {
    if ((c = *ss) == 0) break;
    *ss++ = tolower(c);
    }
return s;
}

char * memSwapChar(char *s, int len, char oldChar, char newChar)
/* Substitute newChar for oldChar throughout memory of given length. */
{
int ix=0;
for (;ix<len;ix++)
    {
    if (s[ix] == oldChar)
        s[ix] =  newChar;
    }
    return s;
}

void stripChar(char *s, char c)
/* Remove all occurences of c from s. */
{
char *in = s, *out = s;
char b;

for (;;)
    {
    b = *out = *in++;
    if (b == 0)
       break;
    if (b != c)
       ++out;
    }
}

char *stripEnclosingChar(char *inout,char encloser)
// Removes enclosing char if found at both beg and end, preserving pointer
// Note: handles brackets '(','{' and '[' by complement at end
{
if(inout == NULL || strlen(inout) < 2 || *inout != encloser)
    return inout;

char *end = inout + (strlen(inout) - 1);
char closer = encloser;
switch (closer)
    {
    case '(': closer = ')'; break;
    case '{': closer = '}'; break;
    case '[': closer = ']'; break;
    default: break;
    }
if(*end  != closer)
    return inout;
*end = '\0';
return memmove(inout,inout+1,strlen(inout));  // use memmove to safely copy in place
}

void stripString(char *s, char *strip)
/* Remove all occurences of strip from s. */
{
char c, *in = s, *out = s;
int stripSize = strlen(strip);
char stripFirst = strip[0];

while ((c = *in) != 0)
    {
    c = *in;
    if (c == stripFirst)
        {
	if (startsWith(strip, in))
	    {
	    in += stripSize;
	    continue;
	    }
	}
    *out = c;
    ++out;
    ++in;
    }
*out = 0;
}

int countChars(char *s, char c)
/* Return number of characters c in string s. */
{
char a;
int count = 0;
while ((a = *s++) != 0)
    if (a == c)
        ++count;
return count;
}

int countCharsN(char *s, char c, int size)
/* Return number of characters c in string s of given size. */
{
int i;
int count = 0;
for (i=0; i<size; ++i)
    if (s[i] == c)
        ++count;
return count;
}

int countLeadingChars(char *s, char c)
/* Count number of characters c at start of string. */
{
int count = 0;
while (*s++ == c)
   ++count;
return count;
}

int countLeadingDigits(const char *s)
/* Return number of leading digits in s */
{
int count = 0;
while (isdigit(*s))
   {
   ++count;
   ++s;
   }
return count;
}

int countLeadingNondigits(const char *s)
/* Count number of leading non-digit characters in s. */
{
int count = 0;
char c;
while ((c = *s++) != 0)
   {
   if (isdigit(c))
       break;
   ++count;
   }
return count;
}

int cmpStringsWithEmbeddedNumbers(const char *a, const char *b)
/* Compare strings such as gene names that may have embedded numbers,
 * so that bmp4a comes before bmp14a */
{
for (;;)
   {
   /* Figure out number of digits at start, and do numerical comparison if there
    * are any.  If numbers agree step over numerical part, otherwise return difference. */
   int aNum = countLeadingDigits(a);
   int bNum = countLeadingDigits(b);
   if (aNum >= 0 && bNum >= 0)
       {
       int diff = atoi(a) - atoi(b);
       if (diff != 0)
           return diff;
       a += aNum;
       b += bNum;
       }

   /* Count number of non-digits at start. */
   int aNonNum = countLeadingNondigits(a);
   int bNonNum = countLeadingNondigits(b);

   // If different sizes of non-numerical part, then don't match, let strcmp sort out how
   if (aNonNum != bNonNum)
       return strcmp(a,b);
   // If no characters left then they are the same!
   else if (aNonNum == 0)
       return 0;
   // Non-numerical part is the same length and non-zero.  See if it is identical.  Return if not.
   else
       {
       int diff = memcmp(a,b,aNonNum);
       if (diff != 0)
            return diff;
       a += aNonNum;
       b += bNonNum;
       }
   }
}

int cmpWordsWithEmbeddedNumbers(const char *a, const char *b)
/* Case insensitive version of cmpStringsWithEmbeddedNumbers. */
{
char *A = cloneString(a);
char *B = cloneString(b);
int diff = cmpStringsWithEmbeddedNumbers(strUpper(A), strUpper(B));
freeMem(A);
freeMem(B);
return diff;
}

int countSame(char *a, char *b)
/* Count number of characters that from start in a,b that are same. */
{
char c;
int i;
int count = 0;
for (i=0; ; ++i)
   {
   c = a[i];
   if (b[i] != c)
       break;
   if (c == 0)
       break;
   ++count;
   }
return count;
}


/* int chopString(in, sep, outArray, outSize); */
/* This chops up the input string (cannabilizing it)
 * into an array of zero terminated strings in
 * outArray.  It returns the number of strings.
 * If you pass in NULL for outArray, it will just
 * return the number of strings that it *would*
 * chop. */
int chopString(char *in, char *sep, char *outArray[], int outSize)
{
int recordCount = 0;

for (;;)
    {
    if (outArray != NULL && recordCount >= outSize)
	break;
    /* Skip initial separators. */
    in += strspn(in, sep);
    if (*in == 0)
	break;
    if (outArray != NULL)
	outArray[recordCount] = in;
    recordCount += 1;
    in += strcspn(in, sep);
    if (*in == 0)
	break;
    if (outArray != NULL)
	*in = 0;
    in += 1;
    }
return recordCount;
}

int chopByWhite(char *in, char *outArray[], int outSize)
/* Like chopString, but specialized for white space separators. */
{
int recordCount = 0;
char c;
for (;;)
    {
    if (outArray != NULL && recordCount >= outSize)
	break;

    /* Skip initial separators. */
    while (isspace(*in)) ++in;
    if (*in == 0)
	break;

    /* Store start of word and look for end of word. */
    if (outArray != NULL)
	outArray[recordCount] = in;
    recordCount += 1;
    for (;;)
        {
        if ((c = *in) == 0)
            break;
        if (isspace(c))
            break;
        ++in;
        }
    if (*in == 0)
	break;

    /* Tag end of word with zero. */
    if (outArray != NULL)
	*in = 0;
    /* And skip over the zero. */
    in += 1;
    }
return recordCount;
}

int chopByWhiteRespectDoubleQuotes(char *in, char *outArray[], int outSize)
/* Like chopString, but specialized for white space separators.
 * Further, any doubleQuotes (") are respected.
 * If doubleQuote is encloses whole string, then they are removed:
 *   "Fred and Ethyl" results in word [Fred and Ethyl]
 * If doubleQuotes exist inside string they are retained:
 *   Fred" and Ethyl" results in word [Fred" and Ethyl"]
 * Special note "" is a valid, though empty word. */
{
int recordCount = 0;
char c;
char *quoteBegins = NULL;
boolean quoting = FALSE;
for (;;)
    {
    if (outArray != NULL && recordCount >= outSize)
        break;

    /* Skip initial separators. */
    while (isspace(*in)) ++in;
    if (*in == 0)
        break;

    /* Store start of word and look for end of word. */
    if (outArray != NULL)
        {
        outArray[recordCount] = in;
        if((*in == '"'))
            quoteBegins = (in+1);
        else
            quoteBegins = NULL;
        }
    recordCount += 1;
    quoting = FALSE;
    for (;;)
        {
        if ((c = *in) == 0)
            break;
        if(quoting)
            {
            if(c == '"')
                {
                quoting = FALSE;
                if(quoteBegins != NULL) // implies out array
                    {
                    if((c = *(in+1) == 0 )|| isspace(c)) // whole word is quoted.
                        {
                        outArray[recordCount-1] = quoteBegins; // Fix beginning of word
                        quoteBegins = NULL;
                        break;
                        }
                    }
                }
            }
        else
            {
            quoting = (c == '"');
            if (isspace(c))
                break;
            }
        ++in;
        }
    if (*in == 0)
        break;

    /* Tag end of word with zero. */
    if (outArray != NULL)
        *in = 0;
    /* And skip over the zero. */
    in += 1;
    }
    return recordCount;
}

int chopByChar(char *in, char chopper, char *outArray[], int outSize)
/* Chop based on a single character. */
{
int i;
char c;
if (*in == 0)
    return 0;
for (i=0; (i<outSize) || (outArray==NULL); ++i)
    {
    if (outArray != NULL)
        outArray[i] = in;
    for (;;)
	{
	if ((c = *in++) == 0)
	    return i+1;
	else if (c == chopper)
	    {
            if (outArray != NULL)
                in[-1] = 0;
	    break;
	    }
	}
    }
return i;
}

char crLfChopper[] = "\n\r";
char whiteSpaceChopper[] = " \t\n\r";


char *skipBeyondDelimit(char *s,char delimit)
/* Returns NULL or pointer to first char beyond one (or more contiguous) delimit char.
   If delimit is ' ' then skips beyond first patch of whitespace. */
{
if(s != NULL)
    {
    char *beyond = NULL;
    if(delimit == ' ')
        return skipLeadingSpaces(skipToSpaces(s));
    else
        beyond = strchr(s,delimit);
    if(beyond != NULL)
        {
        for(beyond++;*beyond == delimit;beyond++);
        if(*beyond != '\0')
            return beyond;
        }
    }
return NULL;
}

char *skipLeadingSpaces(char *s)
/* Return first non-white space. */
{
char c;
if (s == NULL) return NULL;
for (;;)
    {
    c = *s;
    if (!isspace(c))
	return s;
    ++s;
    }
}

/* Return first white space or NULL if none.. */
char *skipToSpaces(char *s)
{
char c;
if (s == NULL)
    return NULL;
for (;;)
    {
    c = *s;
    if (c == 0)
        return NULL;
    if (isspace(c))
	return s;
    ++s;
    }
}



void eraseTrailingSpaces(char *s)
/* Replace trailing white space with zeroes. */
{
int len = strlen(s);
int i;
char c;

for (i=len-1; i>=0; --i)
    {
    c = s[i];
    if (isspace(c))
	s[i] = 0;
    else
	break;
    }
}

/* Remove white space from a string */
void eraseWhiteSpace(char *s)
{
char *in, *out;
char c;

in = out = s;
for (;;)
    {
    c = *in++;
    if (c == 0)
	break;
    if (!isspace(c))
	*out++ = c;
    }
*out++ = 0;
}

/* Remove non-alphanumeric chars from string */
void eraseNonAlphaNum(char *s)
{
char *in, *out;
char c;

in = out = s;
for (;;)
    {
    c = *in++;
    if (c == 0)
        break;
    if (isalnum(c))
        *out++ = c;
    }
*out = 0;
}

char *trimSpaces(char *s)
/* Remove leading and trailing white space. */
{
if (s != NULL)
    {
    s = skipLeadingSpaces(s);
    eraseTrailingSpaces(s);
    }
return s;
}

void repeatCharOut(FILE *f, char c, int count)
/* Write character to file repeatedly. */
{
while (--count >= 0)
    fputc(c, f);
}

void spaceOut(FILE *f, int count)
/* Put out some spaces to file. */
{
repeatCharOut(f, ' ', count);
}

void starOut(FILE *f, int count)
/* Put out some asterisks to file. */
{
repeatCharOut(f, '*', count);
}

boolean hasWhiteSpace(char *s)
/* Return TRUE if there is white space in string. */
{
char c;
while ((c = *s++) != 0)
    if (isspace(c))
        return TRUE;
return FALSE;
}

char *firstWordInLine(char *line)
/* Returns first word in line if any (white space separated).
 * Puts 0 in place of white space after word. */
{
char *e;
line = skipLeadingSpaces(line);
if ((e = skipToSpaces(line)) != NULL)
    *e = 0;
return line;
}

char *cloneFirstWord(char *line)
/* Clone first word in line */
{
char *startFirstWord = skipLeadingSpaces(line);
if (startFirstWord == NULL)
    return NULL;
char *endFirstWord = skipToSpaces(startFirstWord);
if (endFirstWord == NULL)
    return cloneString(startFirstWord);
else
    return cloneStringZ(startFirstWord, endFirstWord - startFirstWord);
}

char *lastWordInLine(char *line)
/* Returns last word in line if any (white space separated).
 * Returns NULL if string is empty.  Removes any terminating white space
 * from line. */
{
char *s = line;
char *word = NULL, *wordEnd = NULL;
for (;;)
    {
    s = skipLeadingSpaces(s);
    if (s == NULL || s[0] == 0)
	break;
    word = s;
    s = wordEnd = skipToSpaces(s);
    if (s == NULL)
        break;
    }
if (wordEnd != NULL)
    *wordEnd = 0;
return word;
}

char *nextWord(char **pLine)
/* Return next word in *pLine and advance *pLine to next
 * word. */
{
char *s = *pLine, *e;
if (s == NULL || s[0] == 0)
    return NULL;
s = skipLeadingSpaces(s);
if (s[0] == 0)
    return NULL;
e = skipToSpaces(s);
if (e != NULL)
    *e++ = 0;
*pLine = e;
return s;
}

char *nextWordRespectingQuotes(char **pLine)
// return next word but respects single or double quotes surrounding sets of words.
{
char *s = *pLine, *e;
if (s == NULL || s[0] == 0)
    return NULL;
s = skipLeadingSpaces(s);
if (s[0] == 0)
    return NULL;
if (s[0] == '"')
    {
    e = skipBeyondDelimit(s+1,'"');
    if (e != NULL && !isspace(e[0]))
        e = skipToSpaces(s);
    }
else if (s[0] == '\'')
    {
    e = skipBeyondDelimit(s+1,'\'');
    if (e != NULL && !isspace(e[0]))
        e = skipToSpaces(s);
    }
else
    e = skipToSpaces(s);
if (e != NULL)
    *e++ = 0;
*pLine = e;
return s;
}

char *nextTabWord(char **pLine)
/* Return next tab-separated word. */
{
char *s = *pLine;
char *e;
if (s == NULL || *s == '\n' || *s == 0)
    {
    *pLine = NULL;
    return NULL;
    }
e = strchr(s, '\t');
if (e == NULL)
    {
    e = strchr(s, '\n');
    if (e != NULL)
        *e = 0;
    *pLine = NULL;
    }
else
    {
    *e++ = 0;
    *pLine = e;
    }
return s;
}

char *cloneFirstWordByDelimiter(char *line,char delimit)
/* Returns a cloned first word, not harming the memory passed in */
{
if(line == NULL || *line == 0)
    return NULL;
line = skipLeadingSpaces(line);
if(*line == 0)
    return NULL;
int size=0;
char *e;
for(e=line;*e!=0;e++)
    {
    if(*e==delimit)
        break;
    else if(delimit == ' ' && isspace(*e))
        break;
    size++;
    }
if(size == 0)
    return NULL;
char *new = needMem(size + 2); // Null terminated by 2
memcpy(new, line, size);
return new;
}

char *cloneNextWordByDelimiter(char **line,char delimit)
/* Returns a cloned first word, advancing the line pointer but not harming memory passed in */
{
char *new = cloneFirstWordByDelimiter(*line,delimit);
if(new != NULL)
    {
    *line = skipLeadingSpaces(*line);
    *line += strlen(new);
    if( **line != 0)
        (*line)++;
    }
return new;
}

char *nextStringInList(char **pStrings)
/* returns pointer to the first string and advances pointer to next in
   list of strings dilimited by 1 null and terminated by 2 nulls. */
{
if(pStrings == NULL || *pStrings == NULL || **pStrings == 0)
    return NULL;
char *p=*pStrings;
*pStrings += strlen(p)+1;
return p;
}

int cntStringsInList(char *pStrings)
/* returns count of strings in a
   list of strings dilimited by 1 null and terminated by 2 nulls. */
{
int cnt=0;
char *p = pStrings;
while(nextStringInList(&p) != NULL)
    cnt++;
return cnt;
}

int stringArrayIx(char *string, char *array[], int arraySize)

/* Return index of string in array or -1 if not there. */
{
int i;
for (i=0; i<arraySize; ++i)
    if (!differentWord(array[i], string))
        return i;
return -1;
}

int ptArrayIx(void *pt, void *array, int arraySize)
/* Return index of pt in array or -1 if not there. */
{
int i;
void **a = array;
for (i=0; i<arraySize; ++i)
    {
    if (pt == a[i])
        return i;
    }
return -1;
}



FILE *mustOpen(char *fileName, char *mode)
/* Open a file - or squawk and die. */
{
FILE *f;

if (sameString(fileName, "stdin"))
    return stdin;
if (sameString(fileName, "stdout"))
    return stdout;
if ((f = fopen(fileName, mode)) == NULL)
    {
    char *modeName = "";
    if (mode)
        {
        if (mode[0] == 'r')
            modeName = " to read";
        else if (mode[0] == 'w')
            modeName = " to write";
        else if (mode[0] == 'a')
            modeName = " to append";
        }
    errAbort("Can't open %s%s: %s", fileName, modeName, strerror(errno));
    }
return f;
}

void mustWrite(FILE *file, void *buf, size_t size)
/* Write to a file or squawk and die. */
{
if (size != 0 && fwrite(buf, size, 1, file) != 1)
    {
    errAbort("Error writing %lld bytes: %s\n", (long long)size, strerror(ferror(file)));
    }
}


void mustRead(FILE *file, void *buf, size_t size)
/* Read size bytes from a file or squawk and die. */
{
if (size != 0 && fread(buf, size, 1, file) != 1)
    {
    if (ferror(file))
	errAbort("Error reading %lld bytes: %s", (long long)size, strerror(ferror(file)));
    else
	errAbort("End of file reading %lld bytes", (long long)size);
    }
}

void writeString(FILE *f, char *s)
/* Write a 255 or less character string to a file.
 * This will write the length of the string in the first
 * byte then the string itself. */
{
UBYTE bLen;
int len = strlen(s);

if (len > 255)
    {
    warn("String too long in writeString (%d chars):\n%s", len, s);
    len = 255;
    }
bLen = len;
writeOne(f, bLen);
mustWrite(f, s, len);
}

char *readString(FILE *f)
/* Read a string (written with writeString) into
 * memory.  freeMem the result when done. */
{
UBYTE bLen;
int len;
char *s;

if (!readOne(f, bLen))
    return NULL;
len = bLen;
s = needMem(len+1);
if (len > 0)
    mustRead(f, s, len);
return s;
}

char *mustReadString(FILE *f)
/* Read a string.  Squawk and die at EOF or if any problem. */
{
char *s = readString(f);
if (s == NULL)
    errAbort("Couldn't read string");
return s;
}


boolean fastReadString(FILE *f, char buf[256])
/* Read a string into buffer, which must be long enough
 * to hold it.  String is in 'writeString' format. */
{
UBYTE bLen;
int len;
if (!readOne(f, bLen))
    return FALSE;
if ((len = bLen)> 0)
    mustRead(f, buf, len);
buf[len] = 0;
return TRUE;
}

void msbFirstWriteBits64(FILE *f, bits64 x)
/* Write out 64 bit number in manner that is portable across architectures */
{
int i;
UBYTE buf[8];
for (i=7; i>=0; --i)
    {
    buf[i] = (UBYTE)(x&0xff);
    x >>= 8;
    }
mustWrite(f, buf, 8);
}

bits64 msbFirstReadBits64(FILE *f)
/* Write out 64 bit number in manner that is portable across architectures */
{
int i;
UBYTE buf[8];
bits64 x = 0;
mustRead(f, buf, 8);
for (i=0; i<8; ++i)
    {
    x <<= 8;
    x |= buf[i];
    }
return x;
}

void mustGetLine(FILE *file, char *buf, int charCount)
/* Read at most charCount-1 bytes from file, but stop after newline if one is
 * encountered.  The string in buf is '\0'-terminated.  (See man 3 fgets.)
 * Die if there is an error. */
{
char *success = fgets(buf, charCount, file);
if (success == NULL && charCount > 0)
    buf[0] = '\0';
if (ferror(file))
    errAbort("mustGetLine: fgets failed: %s", strerror(ferror(file)));
}

int mustOpenFd(char *fileName, int flags)
/* Open a file descriptor (see man 2 open) or squawk and die. */
{
if (sameString(fileName, "stdin"))
    return STDIN_FILENO;
if (sameString(fileName, "stdout"))
    return STDOUT_FILENO;
// mode is necessary when O_CREAT is given, ignored otherwise
int mode = 00664;
int fd = open(fileName, flags, mode);
if (fd < 0)
    {
    char *modeName = "";
    if ((flags & (O_WRONLY | O_CREAT | O_TRUNC)) == (O_WRONLY | O_CREAT | O_TRUNC))
	modeName = " to create and truncate";
    else if ((flags & (O_WRONLY | O_CREAT)) == (O_WRONLY | O_CREAT))
	modeName = " to create";
    else if ((flags & O_WRONLY) == O_WRONLY)
	modeName = " to write";
    else if ((flags & O_RDWR) == O_RDWR)
	modeName = " to append";
    else
	modeName = " to read";
    errnoAbort("Can't open %s%s", fileName, modeName);
    }
return fd;
}

void mustReadFd(int fd, void *buf, size_t size)
/* Read size bytes from a file or squawk and die. */
{
ssize_t actualSize;
char *cbuf = buf;
// using a loop because linux was not returning all data in a single request when request size exceeded 2GB.
while (size > 0) 
    {
    actualSize = read(fd, cbuf, size);
    if (actualSize < 0)
	errnoAbort("Error reading %lld bytes", (long long)size);
    if (actualSize == 0)
	errAbort("End of file reading %llu bytes (got %lld)", (unsigned long long)size, (long long)actualSize);
    cbuf += actualSize;
    size -= actualSize;
    }
}

void mustWriteFd(int fd, void *buf, size_t size)
/* Write size bytes to file descriptor fd or die.  (See man 2 write.) */
{
ssize_t result = write(fd, buf, size);
if (result < size)
    errAbort("mustWriteFd: write failed: %s", strerror(errno));
}

off_t mustLseek(int fd, off_t offset, int whence)
/* Seek to given offset, relative to whence (see man lseek) in file descriptor fd or errAbort.
 * Return final offset (e.g. if this is just an (fd, 0, SEEK_CUR) query for current position). */
{
off_t ret = lseek(fd, offset, whence);
if (ret < 0)
    errnoAbort("lseek(%d, %lld, %s (%d)) failed", fd, (long long)offset,
	       ((whence == SEEK_SET) ? "SEEK_SET" : (whence == SEEK_CUR) ? "SEEK_CUR" :
		(whence == SEEK_END) ? "SEEK_END" : "invalid 'whence' value"), whence);
return ret;
}

void mustCloseFd(int *pFd)
/* Close file descriptor *pFd if >= 0, abort if there's an error, set *pFd = -1. */
{
if (pFd != NULL && *pFd >= 0)
    {
    if (close(*pFd) < 0)
	errnoAbort("close failed");
    *pFd = -1;
    }
}

char *addSuffix(char *head, char *suffix)
/* Return a needMem'd string containing "headsuffix". Should be free'd
 when finished. */
{
char *ret = NULL;
int size = strlen(head) + strlen(suffix) +1;
ret = needMem(sizeof(char)*size);
snprintf(ret, size, "%s%s", head, suffix);
return ret;
}

void chopSuffix(char *s)
/* Remove suffix (last . in string and beyond) if any. */
{
char *e = strrchr(s, '.');
if (e != NULL)
    *e = 0;
}

void chopSuffixAt(char *s, char c)
/* Remove end of string from first occurrence of char c.
 * chopSuffixAt(s, '.') is equivalent to regular chopSuffix. */
{
char *e = strrchr(s, c);
if (e != NULL)
    *e = 0;
}

char *chopPrefixAt(char *s, char c)
/* Like chopPrefix, but can chop on any character, not just '.' */
{
char *e = strchr(s, c);
if (e == NULL) return s;
*e++ = 0;
return e;
}

char *chopPrefix(char *s)
/* This will replace the first '.' in a string with
 * 0, and return the character after this.  If there
 * is no '.' in the string this will just return the
 * unchanged s passed in. */
{
return chopPrefixAt(s, '.');
}



boolean carefulCloseWarn(FILE **pFile)
/* Close file if open and null out handle to it.
 * Return FALSE and print a warning message if there
 * is a problem.*/
{
FILE *f;
boolean ok = TRUE;
if ((pFile != NULL) && ((f = *pFile) != NULL))
    {
    if (f != stdin && f != stdout)
        {
        if (fclose(f) != 0)
	    {
            errnoWarn("fclose failed");
	    ok = FALSE;
	    }
        }
    *pFile = NULL;
    }
return ok;
}

void carefulClose(FILE **pFile)
/* Close file if open and null out handle to it.
 * Warn and abort if there's a problem. */
{
if (!carefulCloseWarn(pFile))
    noWarnAbort();
}

char *firstWordInFile(char *fileName, char *wordBuf, int wordBufSize)
/* Read the first word in file into wordBuf. */
{
FILE *f = mustOpen(fileName, "r");
mustGetLine(f, wordBuf, wordBufSize);
fclose(f);
return trimSpaces(wordBuf);
}

int fileOffsetSizeCmp(const void *va, const void *vb)
/* Help sort fileOffsetSize by offset. */
{
const struct fileOffsetSize *a = *((struct fileOffsetSize **)va);
const struct fileOffsetSize *b = *((struct fileOffsetSize **)vb);
if (a->offset > b->offset)
    return 1;
else if (a->offset == b->offset)
    return 0;
else
    return -1;
}

struct fileOffsetSize *fileOffsetSizeMerge(struct fileOffsetSize *inList)
/* Returns a new list which is inList transformed to have adjacent blocks
 * merged.  Best to use this with a sorted list. */
{
struct fileOffsetSize *newList = NULL, *newEl = NULL, *oldEl, *nextOld;

for (oldEl = inList; oldEl != NULL; oldEl = nextOld)
    {
    nextOld = oldEl->next;
    if (nextOld != NULL && nextOld->offset < oldEl->offset)
        errAbort("Unsorted inList in fileOffsetSizeMerge %llu %llu", oldEl->offset, nextOld->offset);
    if (newEl == NULL || newEl->offset + newEl->size < oldEl->offset)
        {
	newEl = CloneVar(oldEl);
	slAddHead(&newList, newEl);
	}
    else
        {
	newEl->size = oldEl->offset + oldEl->size - newEl->offset;
	}
    }
slReverse(&newList);
return newList;
}

void fileOffsetSizeFindGap(struct fileOffsetSize *list,
	struct fileOffsetSize **pBeforeGap, struct fileOffsetSize **pAfterGap)
/* Starting at list, find all items that don't have a gap between them and the previous item.
 * Return at gap, or at end of list, returning pointers to the items before and after the gap. */
{
struct fileOffsetSize *pt, *next;
for (pt = list; ; pt = next)
    {
    next = pt->next;
    if (next == NULL || next->offset != pt->offset + pt->size)
	{
	*pBeforeGap = pt;
	*pAfterGap = next;
	return;
	}
    }
}


void mustSystem(char *cmd)
/* Execute cmd using "sh -c" or die.  (See man 3 system.) fail on errors */
{
if (cmd == NULL) // don't allow (system() supports testing for shell this way)
    errAbort("mustSystem: called with NULL command.");
int status = system(cmd);
if (status == -1)
    errnoAbort("error starting command: %s", cmd);
else if (WIFSIGNALED(status))
    errAbort("command terminated by signal %d: %s", WTERMSIG(status), cmd);
else if (WIFEXITED(status))
    {
    if (WEXITSTATUS(status) != 0)
        errAbort("command exited with %d: %s", WEXITSTATUS(status), cmd);
    }
else
    errAbort("bug: invalid exit status for command: %s", cmd);
}

int roundingScale(int a, int p, int q)
/* returns rounded a*p/q */
{
if (a > 100000 || p > 100000)
    {
    double x = a;
    x *= p;
    x /= q;
    return round(x);
    }
else
    return (a*p + q/2)/q;
}

int intAbs(int a)
/* Return integer absolute value */
{
return (a >= 0 ? a : -a);
}

int  rangeIntersection(int start1, int end1, int start2, int end2)
/* Return amount of bases two ranges intersect over, 0 or negative if no
 * intersection. */
{
int s = max(start1,start2);
int e = min(end1,end2);
return e-s;
}

int positiveRangeIntersection(int start1, int end1, int start2, int end2)
/* Return number of bases in intersection of two ranges, or
 * zero if they don't intersect. */
{
int ret = rangeIntersection(start1,end1,start2,end2);
if (ret < 0)
    ret = 0;
return ret;
}

bits64 byteSwap64(bits64 a)
/* Return byte-swapped version of a */
{
union {bits64 whole; UBYTE bytes[4];} u,v;
u.whole = a;
v.bytes[0] = u.bytes[7];
v.bytes[1] = u.bytes[6];
v.bytes[2] = u.bytes[5];
v.bytes[3] = u.bytes[4];
v.bytes[4] = u.bytes[3];
v.bytes[5] = u.bytes[2];
v.bytes[6] = u.bytes[1];
v.bytes[7] = u.bytes[0];
return v.whole;
}

bits64 readBits64(FILE *f, boolean isSwapped)
/* Read and optionally byte-swap 64 bit entity. */
{
bits64 val;
mustReadOne(f, val);
if (isSwapped)
    val = byteSwap64(val);
return val;
}

bits64 fdReadBits64(int fd, boolean isSwapped)
/* Read and optionally byte-swap 64 bit entity. */
{
bits64 val;
mustReadOneFd(fd, val);
if (isSwapped)
    val = byteSwap64(val);
return val;
}

bits64 memReadBits64(char **pPt, boolean isSwapped)
/* Read and optionally byte-swap 64 bit entity from memory buffer pointed to by
 * *pPt, and advance *pPt past read area. */
{
bits64 val;
memcpy(&val, *pPt, sizeof(val));
if (isSwapped)
    val = byteSwap64(val);
*pPt += sizeof(val);
return val;
}

bits32 byteSwap32(bits32 a)
/* Return byte-swapped version of a */
{
union {bits32 whole; UBYTE bytes[4];} u,v;
u.whole = a;
v.bytes[0] = u.bytes[3];
v.bytes[1] = u.bytes[2];
v.bytes[2] = u.bytes[1];
v.bytes[3] = u.bytes[0];
return v.whole;
}

bits32 readBits32(FILE *f, boolean isSwapped)
/* Read and optionally byte-swap 32 bit entity. */
{
bits32 val;
mustReadOne(f, val);
if (isSwapped)
    val = byteSwap32(val);
return val;
}

bits32 fdReadBits32(int fd, boolean isSwapped)
/* Read and optionally byte-swap 32 bit entity. */
{
bits32 val;
mustReadOneFd(fd, val);
if (isSwapped)
    val = byteSwap32(val);
return val;
}

bits32 memReadBits32(char **pPt, boolean isSwapped)
/* Read and optionally byte-swap 32 bit entity from memory buffer pointed to by
 * *pPt, and advance *pPt past read area. */
{
bits32 val;
memcpy(&val, *pPt, sizeof(val));
if (isSwapped)
    val = byteSwap32(val);
*pPt += sizeof(val);
return val;
}

bits16 byteSwap16(bits16 a)
/* Return byte-swapped version of a */
{
union {bits16 whole; UBYTE bytes[2];} u,v;
u.whole = a;
v.bytes[0] = u.bytes[1];
v.bytes[1] = u.bytes[0];
return v.whole;
}

bits16 readBits16(FILE *f, boolean isSwapped)
/* Read and optionally byte-swap 16 bit entity. */
{
bits16 val;
mustReadOne(f, val);
if (isSwapped)
    val = byteSwap16(val);
return val;
}

bits16 fdReadBits16(int fd, boolean isSwapped)
/* Read and optionally byte-swap 16 bit entity. */
{
bits16 val;
mustReadOneFd(fd, val);
if (isSwapped)
    val = byteSwap16(val);
return val;
}

bits16 memReadBits16(char **pPt, boolean isSwapped)
/* Read and optionally byte-swap 16 bit entity from memory buffer pointed to by
 * *pPt, and advance *pPt past read area. */
{
bits16 val;
memcpy(&val, *pPt, sizeof(val));
if (isSwapped)
    val = byteSwap16(val);
*pPt += sizeof(val);
return val;
}

double byteSwapDouble(double a)
/* Return byte-swapped version of a */
{
union {double whole; UBYTE bytes[4];} u,v;
u.whole = a;
v.bytes[0] = u.bytes[7];
v.bytes[1] = u.bytes[6];
v.bytes[2] = u.bytes[5];
v.bytes[3] = u.bytes[4];
v.bytes[4] = u.bytes[3];
v.bytes[5] = u.bytes[2];
v.bytes[6] = u.bytes[1];
v.bytes[7] = u.bytes[0];
return v.whole;
}


double readDouble(FILE *f, boolean isSwapped)
/* Read and optionally byte-swap double-precision floating point entity. */
{
double val;
mustReadOne(f, val);
if (isSwapped)
    val = byteSwapDouble(val);
return val;
}

double memReadDouble(char **pPt, boolean isSwapped)
/* Read and optionally byte-swap double-precision floating point entity
 * from memory buffer pointed to by *pPt, and advance *pPt past read area. */
{
double val;
memcpy(&val, *pPt, sizeof(val));
if (isSwapped)
    val = byteSwapDouble(val);
*pPt += sizeof(val);
return val;
}

float byteSwapFloat(float a)
/* Return byte-swapped version of a */
{
union {float whole; UBYTE bytes[4];} u,v;
u.whole = a;
v.bytes[0] = u.bytes[3];
v.bytes[1] = u.bytes[2];
v.bytes[2] = u.bytes[1];
v.bytes[3] = u.bytes[0];
return v.whole;
}


float readFloat(FILE *f, boolean isSwapped)
/* Read and optionally byte-swap single-precision floating point entity. */
{
float val;
mustReadOne(f, val);
if (isSwapped)
    val = byteSwapFloat(val);
return val;
}

float memReadFloat(char **pPt, boolean isSwapped)
/* Read and optionally byte-swap single-precision floating point entity
 * from memory buffer pointed to by *pPt, and advance *pPt past read area. */
{
float val;
memcpy(&val, *pPt, sizeof(val));
if (isSwapped)
    val = byteSwapFloat(val);
*pPt += sizeof(val);
return val;
}


void removeReturns(char *dest, char *src)
/* Removes the '\r' character from a string.
 * The source and destination strings can be the same, if there are
 * no other threads */
{
int i = 0;
int j = 0;

/* until the end of the string */
for (;;)
    {
    /* skip the returns */
    while(src[j] == '\r')
	j++;

    /* copy the characters */
    dest[i] = src[j];

    /* check to see if done */
    if(src[j] == '\0')
	break;

    /* advance the counters */
    i++;
    j++;
    }
}

char* readLine(FILE* fh)
/* Read a line of any size into dynamic memory, return null on EOF */
{
int bufCapacity = 256;
int bufSize = 0;
char* buf = needMem(bufCapacity);
int ch;

/* loop until EOF of EOLN */
while (((ch = getc(fh)) != EOF) && (ch != '\n'))
    {
    /* expand if almost full, always keep one extra char for
     * zero termination */
    if (bufSize >= bufCapacity-2)
        {
        bufCapacity *= 2;
        buf = realloc(buf, bufCapacity);
        if (buf == NULL)
            {
            errAbort("Out of memory in readline - request size %d bytes", bufCapacity);
            }
        }
    buf[bufSize++] = ch;
    }

/* only return EOF if no data was read */
if ((ch == EOF) && (bufSize == 0))
    {
    freeMem(buf);
    return NULL;
    }
buf[bufSize] = '\0';
return buf;
}

boolean fileExists(char *fileName)
/* Return TRUE if file exists (may replace this with non-
 * portable faster way some day). */
{
/* To make piping easier stdin and stdout always exist. */
if (sameString(fileName, "stdin")) return TRUE;
if (sameString(fileName, "stdout")) return TRUE;

return fileSize(fileName) != -1;
}

/*
 Friendly name for strstrNoCase
*/
char *containsStringNoCase(char *haystack, char *needle)
{
return strstrNoCase(haystack, needle);
}

char *strstrNoCase(char *haystack, char *needle)
/*
  A case-insensitive strstr function
Will also robustly handle null strings
param haystack - The string to be searched
param needle - The string to look for in the haystack string

return - The position of the first occurence of the desired substring
or -1 if it is not found
 */
{
char *haystackCopy = NULL;
char *needleCopy = NULL;
int index = 0;
int haystackLen = 0;
int needleLen = 0;
char *p, *q;

if (NULL == haystack || NULL == needle)
    {
    return NULL;
    }

haystackLen = strlen(haystack);
needleLen = strlen(needle);

haystackCopy = (char*) needMem(haystackLen + 1);
needleCopy = (char*) needMem(needleLen + 1);

for(index = 0; index < haystackLen;  index++)
    {
    haystackCopy[index] = tolower(haystack[index]);
    }
haystackCopy[haystackLen] = 0; /* Null terminate */

for(index = 0; index < needleLen;  index++)
    {
    needleCopy[index] = tolower(needle[index]);
    }
needleCopy[needleLen] = 0; /* Null terminate */

p=strstr(haystackCopy, needleCopy);
q=haystackCopy;

freeMem(haystackCopy);
freeMem(needleCopy);

if(p==NULL) return NULL;

return p-q+haystack;
}

int vasafef(char* buffer, int bufSize, char *format, va_list args)
/* Format string to buffer, vsprintf style, only with buffer overflow
 * checking.  The resulting string is always terminated with zero byte. */
{
int sz = vsnprintf(buffer, bufSize, format, args);
/* note that some version return -1 if too small */
if ((sz < 0) || (sz >= bufSize))
    {
    buffer[bufSize-1] = (char) 0;
    errAbort("buffer overflow, size %d, format: %s, buffer: '%s'", bufSize, format, buffer);
    }
return sz;
}

int safef(char* buffer, int bufSize, char *format, ...)
/* Format string to buffer, vsprintf style, only with buffer overflow
 * checking.  The resulting string is always terminated with zero byte. */
{
int sz;
va_list args;
va_start(args, format);
sz = vasafef(buffer, bufSize, format, args);
va_end(args);
return sz;
}

void safecpy(char *buf, size_t bufSize, const char *src)
/* copy a string to a buffer, with bounds checking.*/
{
size_t slen = strlen(src);
if (slen > bufSize-1)
    errAbort("buffer overflow, size %lld, string size: %lld", (long long)bufSize, (long long)slen);
strcpy(buf, src);
}

void safencpy(char *buf, size_t bufSize, const char *src, size_t n)
/* copy n characters from a string to a buffer, with bounds checking.
 * Unlike strncpy, always null terminates the result */
{
if (n > bufSize-1)
    errAbort("buffer overflow, size %lld, substring size: %lld", (long long)bufSize, (long long)n);
size_t slen = strlen(src);
if (slen > n)
    slen = n;
strncpy(buf, src, n);
buf[slen] = '\0';
}

void safecat(char *buf, size_t bufSize, const char *src)
/* Append  a string to a buffer, with bounds checking.*/
{
size_t blen = strlen(buf);
size_t slen = strlen(src);
if (blen+slen > bufSize-1)
    errAbort("buffer overflow, size %lld, new string size: %lld", (long long)bufSize, (long long)(blen+slen));
strcat(buf, src);
}

void safencat(char *buf, size_t bufSize, const char *src, size_t n)
/* append n characters from a string to a buffer, with bounds checking. */
{
size_t blen = strlen(buf);
if (blen+n > bufSize-1)
    errAbort("buffer overflow, size %lld, new string size: %lld", (long long)bufSize, (long long)(blen+n));
size_t slen = strlen(src);
if (slen > n)
    slen = n;
strncat(buf, src, n);
buf[blen+slen] = '\0';
}


static char *naStr = "n/a";
static char *emptyStr = "";

char *naForNull(char *s)
/* Return 'n/a' if s is NULL, otherwise s. */
{
if (s == NULL)
   s = naStr;
return s;
}

char *naForEmpty(char *s)
/* Return n/a if s is "" or NULL, otherwise s. */
{
if (s == NULL || s[0] == 0)
    s = naStr;
return s;
}

char *emptyForNull(char *s)
/* Return "" if s is NULL, otherwise s. */
{
if (s == NULL)
   s = emptyStr;
return s;
}

char *nullIfAllSpace(char *s)
/* Return NULL if s is all spaces, otherwise s. */
{
s = skipLeadingSpaces(s);
if (s != NULL)
    if (s[0] == 0)
        s = NULL;
return s;
}

char *trueFalseString(boolean b)
/* Return "true" or "false" */
{
return (b ? "true" : "false");
}

void uglyTime(char *label, ...)
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */
{
static long lastTime = 0;
long time = clock1000();
va_list args;
va_start(args, label);
if (label != NULL)
    {
    fprintf(stdout, "<span class='timing'>");
    vfprintf(stdout, label, args);
    fprintf(stdout, ": %ld millis<BR></span>\n", time - lastTime);
    }
lastTime = time;
va_end(args);
}

void makeDirs(char* path)
/* make a directory, including parent directories */
{
char pathBuf[PATH_LEN];
char* next = pathBuf;

strcpy(pathBuf, path);
if (*next == '/')
    next++;

while((*next != '\0')
      && (next = strchr(next, '/')) != NULL)
    {
    *next = '\0';
    makeDir(pathBuf);
    *next = '/';
    next++;
    }
makeDir(pathBuf);
}

char *skipNumeric(char *s)
/* Return first char of s that's not a digit */
{
while (isdigit(*s))
   ++s;
return s;
}

char *skipToNumeric(char *s)
/* skip up to where numeric digits appear */
{
while (*s != 0 && !isdigit(*s))
    ++s;
return s;
}

char *splitOffNonNumeric(char *s)
/* Split off non-numeric part, e.g. mm of mm8. Result should be freed when done */
{
return cloneStringZ(s,skipToNumeric(s)-s);
}

char *splitOffNumber(char *db)
/* Split off number part, e.g. 8 of mm8. Result should be freed when done */
{
return cloneString(skipToNumeric(db));
}


time_t mktimeFromUtc (struct tm *t)
/* Return time_t for tm in UTC (GMT)
 * Useful for stuff like converting to time_t the
 * last-modified HTTP response header
 * which is always GMT. Returns -1 on failure of mktime */
{
    time_t time;
    char *tz;
    char save_tz[100];
    tz=getenv("TZ");
    if (tz)
        safecpy(save_tz, sizeof(save_tz), tz);
    setenv("TZ", "GMT0", 1);
    tzset();
    t->tm_isdst = 0;
    time=mktime(t);
    if (tz)
        setenv("TZ", save_tz, 1);
    else
        unsetenv("TZ");
    tzset();
    return (time);
}


time_t dateToSeconds(const char *date,const char*format)
// Convert a string date to time_t
{
    struct tm storage={0,0,0,0,0,0,0,0,0};
    if(strptime(date,format,&storage)==NULL)
        return 0;
    else
        return mktime(&storage);
}

boolean dateIsOld(const char *date,const char*format)
// Is this string date older than now?
{
time_t test = dateToSeconds(date,format);
time_t now = clock1();
return (test < now);
}

static int daysOfMonth(struct tm *tp)
/* Returns the days of the month given the year */
{
int days=0;
switch(tp->tm_mon)
    {
    case 3:
    case 5:
    case 8:
    case 10:    days = 30;   break;
    case 1:     days = 28;
                if( (tp->tm_year % 4) == 0
                && ((tp->tm_year % 20) != 0 || (tp->tm_year % 100) == 0) )
                    days = 29;
                break;
    default:    days = 31;   break;
    }
return days;
}

static void dateAdd(struct tm *tp,int addYears,int addMonths,int addDays)
/* Add years,months,days to a date */
{
tp->tm_mday  += addDays;
tp->tm_mon   += addMonths;
tp->tm_year  += addYears;
int dom=28;
while( (tp->tm_mon >11  || tp->tm_mon <0)
    || (tp->tm_mday>dom || tp->tm_mday<1) )
    {
    if(tp->tm_mon>11)   // First month: tm.tm_mon is 0-11 range
        {
        tp->tm_year += (tp->tm_mon / 12);
        tp->tm_mon  = (tp->tm_mon % 12);
        }
    else if(tp->tm_mon<0)
        {
        tp->tm_year += (tp->tm_mon / 12) - 1;
        tp->tm_mon  =  (tp->tm_mon % 12) + 12;
        }
    else
        {
        dom = daysOfMonth(tp);
        if(tp->tm_mday>dom)
            {
            tp->tm_mday -= dom;
            tp->tm_mon  += 1;
            dom = daysOfMonth(tp);
            }
        else if(tp->tm_mday < 1)
            {
            tp->tm_mon  -= 1;
            dom = daysOfMonth(tp);
            tp->tm_mday += dom;
            }
        }
    }
}

char *dateAddTo(char *date,char *format,int addYears,int addMonths,int addDays)
/* Add years,months,days to a formatted date and returns the new date as a cloned string
*  format is a strptime/strftime format: %F = yyyy-mm-dd */
{
char *newDate = needMem(12);
struct tm tp;
if(strptime (date,format, &tp))
    {
    dateAdd(&tp,addYears,addMonths,addDays); // tp.tm_year only contains years since 1900
    strftime(newDate,12,format,&tp);
    }
return cloneString(newDate);  // newDate is never freed!
}

