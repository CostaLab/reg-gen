/* Fixed colors - always in color map. */
#include "common.h"
#include "memgfx.h"

static char const rcsid[] = "$Id: fixColor.c,v 1.3 2003/05/06 07:33:42 kate Exp $";

struct rgbColor mgFixedColors[9] = {
/* These correspond to MG_WHITE, MG_BLACK, etc. */
    { 255, 255, 255},
    { 0, 0, 0},
    { 255, 0, 0},
    { 0, 255, 0},
    { 0, 0, 255},
    { 0, 255, 255},
    { 255, 0, 255},
    { 255, 255, 0},
    { 140, 140, 140},
};



