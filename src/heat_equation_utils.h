#ifndef HEAT_EQUATION_UTILS_H
#define HEAT_EQUATION_UTILS_H

#include <math.h>

double maxd(double v1, double v2)
{
  return v1 < v2 ? v2 : v1;
}

float maxf(float v1, float v2)
{
  return v1 < v2 ? v2 : v1;
}

long double maxld(long double v1, long double v2)
{
  return v1 < v2 ? v2 : v1;
}

#define fabs(x) _Generic((x),        \
                         long double \
                         : fabsl,    \
                           float     \
                         : fabsf,    \
                           default   \
                         : fabs)(x)

#define fmax(x, y) _Generic((x) + (y),  \
                            long double \
                            : maxld,    \
                              float     \
                            : maxf,     \
                              default   \
                            : maxd)(x, y)

#define xstr(s) str(s)
#define str(s) #s

#endif //  HEAT_EQUATION_UTILS_H
