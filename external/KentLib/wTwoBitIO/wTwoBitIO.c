// TwoBit file manipulation
#include <Python.h>
#include <string.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dnaseq.h"
#include "fa.h"
#include "twoBit.h"

static void unknownToN(char *s, int size)
/* Convert non ACGT characters to N. */
{
char c;
int i;
for (i=0; i<size; ++i)
    {
    c = s[i];
    if (ntChars[(int)c] == 0)
        {
    if (isupper(c))
        s[i] = 'N';
    else
        s[i] = 'n';
    }
    }
}



static PyObject *faToTwoBit(PyObject *self, PyObject *args)
{
	char *faname;
	char *tbfname;
	struct twoBit *twoBitList = NULL, *twoBit;
	if(!PyArg_ParseTuple(args,"ss",&faname, &tbfname)) return Py_BuildValue("i",1);
	// faToTwoBit
	struct hash *uniqHash = newHash(18);
	FILE *f;
	struct lineFile *lf = lineFileOpen(faname, TRUE);
	struct dnaSeq seq;
	ZeroVar(&seq);
	while (faMixedSpeedReadNext(lf, &seq.dna, &seq.size, &seq.name))
	{
		// seq size check
		if( seq.size ==0 )
		{
			warn("Skipping item %s which has no sequence.\n",seq.name);
			continue;
		}
		// Duplication check
        if (hashLookup(uniqHash, seq.name))
        {
            errAbort("Duplicate sequence name %s", seq.name);
        }
	    hashAdd(uniqHash, seq.name, NULL);
        unknownToN(seq.dna, seq.size);
    	twoBit = twoBitFromDnaSeq(&seq, TRUE);
    	slAddHead(&twoBitList, twoBit);
	}
	lineFileClose(&lf);
	slReverse(&twoBitList);
	f = mustOpen(tbfname, "wb");
	twoBitWriteHeader(twoBitList, f);
	for (twoBit = twoBitList; twoBit != NULL; twoBit = twoBit->next)
    {
    	twoBitWriteOne(twoBit, f);
    }
	carefulClose(&f);
	return Py_BuildValue("i",0);
}


char *cgetSeq(char *tbfname, char *seqID, 
        int start, int end)
/* Output sequence. */
{
    struct twoBitFile *tbf;
    tbf=twoBitOpen(tbfname);
    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, seqID, start, end);
    /*printf(seq->name);
    printf(seq->dna);*/
    twoBitClose(&tbf);
    return seq->dna;
};
static PyObject *getTwoBitSeqNames(PyObject *self,PyObject *args)
{
    char *tbfname;
    struct slName *s;
    struct twoBitIndex *index;
    PyObject *list=Py_BuildValue("[]");
    struct twoBitFile *tbf;
    if(!PyArg_ParseTuple(args,"s",&tbfname)) return NULL;
    tbf=twoBitOpen(tbfname);
    for (index = tbf->indexList; index != NULL; index = index->next)
    {
        PyObject *name=Py_BuildValue("s",index->name);
        PyList_Append(list,name);
    }
    twoBitClose(&tbf);
    return list;
};

static PyObject *getSeq(PyObject * self, PyObject *args)
{
    char *tbfname,*seqid;
    int start,end;
    char *seq;
    if(!PyArg_ParseTuple(args,"ssii",&tbfname,&seqid,&start,&end)) return NULL;
    seq=cgetSeq(tbfname,seqid,start,end);
    return Py_BuildValue("s",seq);
};
static PyObject *getSeqSizes(PyObject *self, PyObject *args)
{
    char *tbfname;
    struct twoBitFile *tbf;
    struct twoBitIndex *index;
    int size;
    if(!PyArg_ParseTuple(args,"s",&tbfname)) return NULL;
    tbf=twoBitOpen(tbfname);
    PyObject *hash=PyDict_New();
    for (index = tbf->indexList; index != NULL; index = index->next)
    {
        PyObject *name=Py_BuildValue("s",index->name);
        size=twoBitSeqSize(tbf,index->name);
        PyObject *psize=Py_BuildValue("i",size);
        PyDict_SetItem(hash,name,psize);
    }
    twoBitClose(&tbf);
    return hash;
    
}
static struct PyMethodDef wTwoBitIOMethods[]=
{
    {"getSeq",getSeq,1},
    {"getTwoBitSeqNames",getTwoBitSeqNames,1},
    {"getSeqSizes",getSeqSizes,1},
	{"faToTwoBit", faToTwoBit, 1},
    {NULL,NULL}
};

void initwTwoBitIO()
{
    PyObject * m;
    m=Py_InitModule("wTwoBitIO",wTwoBitIOMethods);
};



