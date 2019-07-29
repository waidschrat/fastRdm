/* parameters.c - functions related to checking/setting parameter values
 *
 * Copyright (C) 2006-2012  Jochen Voss, Andreas Voss.
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

#include "fast-dm.h"

#include <math.h>


double
check_bounds(const double *para, const struct samples *samples, int check_zr)
/* Check whether the parameters 'para' are valid.
 *
 * If 'samples' is non-null, also check if t0, st0 and d are
 * compatible with the smallest reposonse times.  If 'check_zr' is
 * non-zero, also check whether the combination of zr and szr is
 * valid.
 *
 * If the parameter set is valid, 0 is returned.  In case of invalid
 * parameters, a value >1 is returned, the magnitude gives the
 * 'badness' of the violation.
 */
{
	double  t_min, penalty = 0;
	int  bad = 0;

	if (para[p_a] <= 0) {
		bad = 1;
		penalty += -para[p_a];
	} else if (para[p_a] > LIMIT_a) {
		bad = 1;
		penalty += para[p_a] - LIMIT_a;
	}
	if ((para[p_t0] - 0.5*fabs(para[p_d]) - 0.5*para[p_st0]) < 0) {
		bad = 1;
		penalty += -(para[p_t0] - 0.5*fabs(para[p_d]) - 0.5*para[p_st0]);
	}
	if (check_zr) {
		if (para[p_count] - 0.5*para[p_szr] < 0) {
			bad = 1;
			penalty += 0.5*para[p_szr] - para[p_count];
		}
		if (para[p_count] + 0.5*para[p_szr] > 1) {
			bad = 1;
			penalty += para[p_count] + 0.5*para[p_szr] - 1;
		}
	}
	if (para[p_szr] < 0) {
		bad = 1;
		penalty += -para[p_szr];
	} else if (para[p_szr] > 1) {
		bad = 1;
		penalty += para[p_szr]-1;
	}
	if (para[p_sv] < 0) {
		bad = 1;
		penalty += -para[p_sv];
	}
	if (para[p_st0] < 0) {
		bad = 1;
		penalty += -para[p_st0];
	}
	if (fabs(para[p_v]) > LIMIT_v) {
		bad = 1;
		penalty += fabs(para[p_v]) - LIMIT_v;
	}
	if (para[p_p] > 1) {
		bad = 1;
		penalty += para[p_p] - 1;
	}
	if (para[p_p] < 0) {
		bad = 1;
		penalty -= para[p_p];
	}

	if (samples) {
		if (samples->plus_used > 0) {
			t_min = para[p_t0] - 0.5*para[p_st0] - 0.5*para[p_d];
			if (samples->plus_data[0] < t_min) {
				bad = 1;
				penalty += t_min - samples->plus_data[0];
			}
		}
		if (samples->minus_used > 0) {
			t_min = para[p_t0]- 0.5*para[p_st0] + 0.5*para[p_d] ;
			if (samples->minus_data[0] < t_min) {
				bad = 1;
				penalty += t_min - samples->minus_data[0];
			}
		}
	}

	/* avoid problems caused by rounding errors */
	return  bad ? penalty+1.0 : 0.0;
}

void
initialise_parameters(const struct dataset *ds, double *x, double *eps)
{
	double  def_x [p_count], def_eps [p_count];
	int  i;

	def_x[p_a] = 1;  def_eps[p_a] = 0.5;
	def_x[p_v] = 0;  def_eps[p_v] = 1;
	def_x[p_t0] = 0.3;  def_eps[p_t0] = 0.5;
	def_x[p_d] = 0.0;  def_eps[p_d] = 0.05;
	def_x[p_szr] = 0.3;  def_eps[p_szr] = 0.2;
	def_x[p_sv] = 0.5;  def_eps[p_sv] = 0.2;
	def_x[p_st0] = 0.2;  def_eps[p_st0] = 0.1;
	def_x[p_p] = 0.1;  def_eps[p_p] = 0.05;

	// trick: we run through the commands list backwards, so that
	// c_run is seen before c_copy_param.  This way, we can have c_run
	// modify the default values and c_copy_param then moves the
	// adjusted values into the correct slots.
	for (i=ds->cmds_used-1; i>=0; --i) {
		int  arg1 = ds->cmds[i].arg1;
		int  arg2 = ds->cmds[i].arg2;
		switch (ds->cmds[i].cmd) {
		case c_run:
			EZ_par(ds->samples[arg1],
				   def_x+p_a, def_x+p_v, def_x+p_t0);
			if (def_x[p_t0]-def_x[p_st0]/2 < 0.01) {
				/* avoid invalid starting point for t0 */
				def_x[p_t0] = def_x[p_st0] / 2 + 0.01;
			}
			break;
		case c_copy_param:
			if (arg1 < 0)  break;
			x[arg2] = def_x[arg1];
			eps[arg2] = def_eps[arg1];
			break;
		default:
			break;
		}
	}
}
