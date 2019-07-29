/* util.h - declarations for "util.c"
 *
 * Copyright (C) 2013  Jochen Voss, Andreas Voss.
 */

#ifndef FILE_UTIL_H_SEEN
#define FILE_UTIL_H_SEEN

#include <stdio.h>

extern void params_write(FILE *fd, double para[p_count],
						 double zr, double precision, int n);
extern void params_check(double para[p_count], double zr);

#endif /* FILE_UTIL_H_SEEN */
