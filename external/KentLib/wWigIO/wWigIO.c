#include <Python.h>
#include <stdlib.h>
//#include <stdio.h> // for test output
#include <string.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "localmem.h"
#include "udc.h"
#include "bigWig.h"
#include "obscure.h"
#include "bwgInternal.h"
#include "zlibFace.h"

char *clChrom = NULL;
int clStart = -1; 
int clEnd = -1;

// wigToBigWig convension
static PyObject *wigToBigWig(PyObject *self,PyObject *args)
{
	char* wigFile;
	char* sizeFile;
	char* bwFile;
	if(!PyArg_ParseTuple(args,"sss", &wigFile, &sizeFile, &bwFile))
		return Py_BuildValue("");
	bigWigFileCreate(wigFile, sizeFile, 256, 1024, FALSE, FALSE, bwFile);
	return Py_BuildValue("");
};

// bigWigToWig convension
static PyObject *bigWigToWig(PyObject *self, PyObject *args)
{
	char *inFile;
	char *outFile;
	if(!PyArg_ParseTuple(args,"ss", &inFile, &outFile))
		return Py_BuildValue("");
	
	// 
	FILE *f = mustOpen(outFile, "w");
	struct bbiFile *bwf = bigWigFileOpen(inFile);
	struct bbiChromInfo *chrom, *chromList = bbiChromList(bwf);
	for (chrom = chromList; chrom != NULL; chrom = chrom->next)
	{
		bigWigIntervalDump(bwf, chrom->name, 0, chrom->size, 0, f);
	}
	bbiChromInfoFreeList(&chromList);
	carefulClose(&f);
	bbiFileClose(&bwf);
	return Py_BuildValue("");
};



// struct BigWigFile 
struct BigWigFile
{
	char inFile[255]; // File name.
	short  count; // number of handles for this file
	struct bbiFile *bwf; // File handle
	struct bbiChromInfo *chromList; // chrom list
	struct BigWigFile *next; // next BigWigFile
};

// BigWigList. Save file handels of bigwig files to avoid frequent opening and closing. Allow to open multiple files.
static struct BigWigFile *BigWigList=NULL;

// Cleanup BigWigFile struct
static void BigWigFileCleanup( char* inFile)
{
	// Find position of file
	struct BigWigFile *p,*q;
	p=BigWigList;
	// First node
	if (sameString(inFile, p->inFile)) 
	{
		if (p->count == 1)
		{
			BigWigList=p->next;
			bbiFileClose(&(p->bwf));
			bbiChromInfoFreeList(&(p->chromList));
			free(p);
			// printf("close file(first node):%s\n", inFile);
		}
		else
		{
			p->count--; // count of file handles - 1
			// printf("number of %s -1\n",inFile);
		}
		return;
	}
	// Other nodes
	for (;p->next!=NULL;p=p->next)
		if (sameString(inFile, p->next->inFile))
			break;
	// File name not in list
	if (p->next==NULL) return;
	// file name found
	if (p->next->count == 1)
	{
		q=p->next;
		p->next=p->next->next; 
		bbiFileClose(&(q->bwf));
		bbiChromInfoFreeList(&(q->chromList));
		free(q);
		// printf("close file(not first node):%s\n", inFile);
	}
	else
	{
		p->next->count--;
		// printf("number of %s -1\n",inFile);
	}
	return ;
}

// close wig file 
static   PyObject *closeWig(PyObject *self,PyObject *args)
{
	char* inFile;
	if (BigWigList==NULL) // empty node list
		return Py_BuildValue("");
	if(!PyArg_ParseTuple(args,"s", &inFile)) return Py_BuildValue("");
	BigWigFileCleanup(inFile);
	return Py_BuildValue("");
}

// Build BigWigFile struct
static void BigWigFileBuild( char *inFile)
{
	struct bbiFile *bwf = NULL;
	struct bbiChromInfo *chromList = NULL;
	struct BigWigFile *p;
	p = BigWigList;

	// if file exists, count ++
	for(;p!=NULL;p=p->next)
	{
		if( sameString(inFile, p->inFile))
		{
			p->count++;
			// printf("number of %s +1\n",inFile);
			return;
		}
	}

	// Open bigwig file
	bwf = bigWigFileOpen(inFile);
	chromList = bbiChromList(bwf);
	// Create new BigWigFile node
	struct BigWigFile *cur;
	cur=(struct BigWigFile *)malloc(sizeof(struct BigWigFile));
	strcpy(cur->inFile,inFile);
	cur->bwf=bwf;
	cur->count=1;
	cur->chromList=chromList;
	cur->next=BigWigList; // Add node to top of BigWigList;
	BigWigList=cur;
	// printf("open file:%s\n",inFile);
	return ;
}

// Open bigwig file
static PyObject *openWig(PyObject *self,PyObject *args)
{
	char *inFile;

	if(!PyArg_ParseTuple(args,"s",&inFile)) return Py_BuildValue("");
	BigWigFileBuild(inFile);
	return Py_BuildValue("");
}

static PyObject *getChromSize(PyObject *self,PyObject *args)
{
	char* inFile;
	char *clChrom;
	struct bbiChromInfo *chrom=NULL;
	PyObject *chroms = Py_BuildValue("[]");
	PyObject *sizes = Py_BuildValue("[]");

	if(!PyArg_ParseTuple(args,"s",&inFile)) return Py_BuildValue("");
	struct BigWigFile *p;
	for(p=BigWigList;p!=NULL;p=p->next) // Find inFile name
		if(sameString(inFile,p->inFile))
		{
			chrom=p->chromList;
			break;
		}
	for(; chrom !=NULL; chrom = chrom->next)
	{
		PyList_Append(chroms,Py_BuildValue("s",chrom->name));
		PyList_Append(sizes,Py_BuildValue("i",chrom->size));
	}
	PyObject *rslt = PyTuple_New(2);
	PyTuple_SetItem(rslt, 0, chroms);
	PyTuple_SetItem(rslt, 1, sizes);
	return rslt;
}

// Get intervals from BigWig file. Coordinates are 0 based.
static PyObject *getIntervals(PyObject *self,PyObject *args) 
{
	//Parse the args
	char *clChrom;
	int clStart,clEnd;
	char *inFile;
	PyObject *intervals=Py_BuildValue("[]");
	if(!PyArg_ParseTuple(args,"ssii",&inFile,&clChrom,&clStart,&clEnd)) return Py_BuildValue("i",1);
	struct bbiFile *bwf=NULL;
	struct bbiChromInfo *chrom=NULL;
	struct BigWigFile *p=NULL;

	// Find inFile
	for(p=BigWigList;p!=NULL;p=p->next)
	{
		if(sameString(inFile,p->inFile))
		{
			bwf=p->bwf;
			chrom=p->chromList;
			break;
		}
	}
	
	// if not found, return "2"
	if (bwf==NULL) return Py_BuildValue("i",2);

	for (; chrom != NULL; chrom = chrom->next)
	{
		if (clChrom != NULL && !sameString(clChrom, chrom->name))
			continue;
		char *chromName = chrom->name;
		struct lm *lm = lmInit(0);
		int start = 0, end = chrom->size;
		if (clStart > 0)
			start = clStart;
		if (clEnd > 0)
			end = clEnd;
		struct bbiInterval *interval, *intervalList = bigWigIntervalQuery(bwf, chromName, start, end, lm);
		for (interval = intervalList; interval != NULL; interval = interval->next)
		{
			PyObject *rslt = PyTuple_New(3);
			PyTuple_SetItem(rslt, 0, Py_BuildValue("i",interval->start));
			PyTuple_SetItem(rslt, 1, Py_BuildValue("i",interval->end));
			PyTuple_SetItem(rslt, 2, Py_BuildValue("f",interval->val));
			PyList_Append(intervals,rslt);
		}
		lmCleanup(&lm);
	}
	return intervals;
}

static struct PyMethodDef wWigIOMethods[]=
{
		{"close",closeWig,1},
		{"open",openWig,1},
		{"getIntervals",getIntervals,1},
		{"getChromSize",getChromSize,1},
		{"wigToBigWig",wigToBigWig,1},
		{"bigWigToWig",bigWigToWig,1},
		{NULL,NULL}
};  


void initwWigIO()
{
	    PyObject * m;
		    m=Py_InitModule("wWigIO",wWigIOMethods);
}


