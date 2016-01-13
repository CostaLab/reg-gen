/* rainbow - stuff to generate rainbow colors. */

#ifndef RAINBOW_H
#define RAINBOW_H

#ifndef MEMGFX_H
#include "memgfx.h"
#endif

struct rgbColor saturatedRainbowAtPos(double pos);
/* Given pos, a number between 0 and 1, return a saturated rainbow rgbColor
 * where 0 maps to red,  0.1 is orange, and 0.9 is violet and 1.0 is back to red */

struct rgbColor lightRainbowAtPos(double pos);
/* Given pos, a number between 0 and 1, return a lightish rainbow rgbColor
 * where 0 maps to red,  0.1 is orange, and 0.9 is violet and 1.0 is back to red */

#endif /* RAINBOW_H */
