/* util.c - auxiliary function for the command line tools
 *
 * Copyright (C) 2013  Jochen Voss, Andreas Voss.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301 USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fast-dm.h"
#include "util.h"

void
params_write(FILE *fd, double para[p_count],
			 double zr, double precision, int n)
{
	fprintf(fd, "a = %g\n", para[p_a]);
	fprintf(fd, "zr = %g\n", zr);
	fprintf(fd, "v = %g\n", para[p_v]);
	fprintf(fd, "t0 = %g\n", para[p_t0]);
	fprintf(fd, "d = %g\n", para[p_d]);
	fprintf(fd, "szr = %g\n", para[p_szr]);
	fprintf(fd, "sv = %g\n", para[p_sv]);
	fprintf(fd, "st0 = %g\n", para[p_st0]);
	fprintf(fd, "precision = %g\n", precision);
	if (n > 0) {
		fprintf(fd, "n = %d\n", n);
	}
}

void
params_check(double para[p_count], double zr)
{
	if (para[p_a] <= 0) {
		fprintf(stderr, "error: invalid parameter a=%g\n",
				para[p_a]);
		exit(1);
	}
	if (para[p_szr] < 0 || para[p_szr] > 1) {
		fprintf(stderr, "error: invalid parameter szr=%g\n",
				para[p_szr]);
		exit(1);
	}
	if (para[p_st0] < 0) {
		fprintf(stderr, "error: invalid parameter st0=%g\n",
				para[p_st0]);
		exit(1);
	}
	if (para[p_sv] < 0) {
		fprintf(stderr, "error: invalid parameter sv=%g\n",
				para[p_sv]);
		exit(1);
	}

	if (para[p_t0] - fabs(0.5*para[p_d]) - 0.5*para[p_st0] < 0) {
		fprintf(stderr,
				"error: invalid parameter combination t0=%g, d=%g, st0=%g\n",
				para[p_t0], para[p_d], para[p_st0]);
		exit(1);
	}
	if (zr - 0.5*para[p_szr] <= 0) {
		fprintf(stderr,
				"error: invalid parameter combination zr=%g, szr=%g\n",
				zr, para[p_szr]);
		exit(1);
	}
	if (zr + 0.5*para[p_szr] >= 1) {
		fprintf(stderr,
				"error: invalid parameter combination zr=%g, szr=%g\n",
				zr, para[p_szr]);
		exit(1);
	}
}
