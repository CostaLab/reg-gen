#ifndef _LIBRGT_H_
#define _LIBRGT_H_

bool overlap(const char *chromA, const int initialA, const int finalA, const char *chromB, const int initialB, const int finalB);

int compareGenomicRegions(const char *chromA, const int initialA, const int finalA, const char *chromB, const int initialB, const int finalB);

void intersectGenomicRegionSetsOverlap (
    const char **chromsA,
    const int *initialsA,
    const int *finalsA,
    const int sizeA,
    const char **chromsB,
    const int *initialsB,
    const int *finalsB,
    const int sizeB,
    int **indicesR,
    int **initialsR,
    int **finalsR,
    int *sizeR
);

void intersectGenomicRegionSetsOriginal (
    const char **chromsA,
    const int *initialsA,
    const int *finalsA,
    const int sizeA,
    const char **chromsB,
    const int *initialsB,
    const int *finalsB,
    const int sizeB,
    int **indicesR,
    int **initialsR,
    int **finalsR,
    int *sizeR
);


void intersectGenomicRegionSetsCompletelyIncluded (
    const char **chromsA,
    const int *initialsA,
    const int *finalsA,
    const int sizeA,
    const char **chromsB,
    const int *initialsB,
    const int *finalsB,
    const int sizeB,
    int **indicesR,
    int **initialsR,
    int **finalsR,
    int *sizeR
);

void intersectGenomicRegionSets (
    const int overlapType,
    const char **chromsA,
    const int *initialsA,
    const int *finalsA,
    const int sizeA,
    const char **chromsB,
    const int *initialsB,
    const int *finalsB,
    const int sizeB,
    int **indicesR,
    int **initialsR,
    int **finalsR,
    int *sizeR
);

double jaccard (
    const char **chromsA,
    const int *initialsA,
    const int *finalsA,
    const int sizeA,
    const char **chromsB,
    const int *initialsB,
    const int *finalsB,
    const int sizeB
);

#endif // _LIBRGT_H_
