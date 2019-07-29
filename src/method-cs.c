/* method-cs.c - Chi-Square method for parameter estimation
 *
 * Copyright (C) 2007-2012  Jochen Voss, Andreas Voss.
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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fast-dm.h"


struct distance_cs  {
	int  N;
	double *CS;              /* array of distances (length N+1) */
	int *df;
	double  zmin, zmax;     /* absolute z-values for CS[0] and CS[N] */
	double  a;              /* interval width used when computing CS */
};


static struct distance_cs *
new_distance(int N, double zmin, double zmax, double a)
/* Allocate a new structure to hold N+1 values in [zmin,zmax] */
{
	struct distance_cs *d;

	assert(N >= 0);

	d = xnew(struct distance_cs, 1);
	d->N = N;
	d->CS = xnew(double, N+1);
	d->df = xnew(int, N+1);
	d->zmin = zmin;
	d->zmax = zmax;
	d->a = a;

	return d;
}

static void
delete_distance(struct distance_cs *d)
{
	xfree(d->CS);
	xfree(d->df);
	xfree(d);
}

static double
distance_interpolate(const struct distance_cs *d, double z)
/* Use linear interpolation to get the distance for a given
 * absolute value of z.
 */
{
	int  N = d->N;
	double  dz, ratio;
	int  step_before;

	if (N == 0)
		return d->CS[0];

	dz = (d->zmax - d->zmin) / N;
	step_before = (int)((z - d->zmin) / dz);
	if (step_before >= N)  step_before = N-1;
	ratio = (z - (d->zmin + dz*step_before)) / dz;
	return  (1-ratio) * d->CS[step_before]
		+ ratio * d->CS[step_before+1];
}

static void
find_best_CS(struct distance_cs *const*result, const struct dataset *ds,
			 const int *z_used, double dist_ret[2], double *zr_ret)
/* Combine the T-values from 'result' into a common p-value.
 *
 * The function computes p-values as the product of the probabilities
 * for all Ts from 'result' and optimises over 'z'.  This takes the
 * different experimental conditions into account.
 *
 * On entry 'zr_ret' must point to an array of size 'ds->z->used'.
 *
 * The distance is returned in 'dist_ret'.  If the parameters are
 * invalid, dist_ret[0]>0 on exit.  Otherwise dist_ret[0]=0 and
 * dist_ret[1] = -log(p).  The relative z parameter values
 * corresponding to this (minimal) distance are returned in '*zr_ret'.
 */
{
	double  total_CS;
	double  penalty;
	int  k;

	total_CS = 0;
	penalty = 0;
	for (k=0; k<ds->z->used; k++) {
		double  dz;
		double  best_CS, best_zr;
		double  zrmin, zrmax;
		int  i, j, n;

		zrmin = -DBL_MAX;
		zrmax = DBL_MAX;
		assert(ds->samples_used > 0);
		for (j=0; j<ds->samples_used; ++j) {
			if (z_used[j]!=k) continue;
			if (zrmin < result[j]->zmin / result[j]->a)
				zrmin = result[j]->zmin / result[j]->a;
			if (zrmax > result[j]->zmax / result[j]->a)
				zrmax = result[j]->zmax / result[j]->a;
		}
		if (zrmax < zrmin) {
			/* no common z is available: abort the
			 * computation and return a penalty value >1.  */
			zr_ret[k] = 0.5;
			penalty = 1 + (zrmin-zrmax);
			break;
		}

		/* Use sub-sampling of the interval zmin..zmax
		 * to find the best z.  */
		n = (int)((zrmax - zrmin) / 0.0001 + 1.5);
		dz = (zrmax - zrmin) / n;
		best_CS = DBL_MAX;
		best_zr = 0.5;
		for (i=0; i<=n; ++i) {
			double  zr = zrmin + i*dz;
			double  CS = 0;

			for (j=0; j<ds->samples_used; ++j) {
				if (z_used[j] != k)  continue;
				CS += distance_interpolate(result[j], zr * result[j]->a);
			}
			if (CS < best_CS) {
				best_CS = CS;
				best_zr = zr;
			}
		}
		total_CS += best_CS;
		zr_ret[k] = best_zr;
	}

	if (penalty > 0) {
		dist_ret[0] = penalty;
		dist_ret[1] = 0;
	} else {
		dist_ret[0] = 0;
		dist_ret[1] = total_CS;
	}
}

static void
find_fixed_CS(struct distance_cs *const*result, const struct dataset *ds,
			  double dist_ret[2], double zr)
/* Combine the T-values from 'result' into a common p-value.
 *
 * 'z' gives the fixed z-value.
 *
 * If the parameters are invalid, dist_ret[0]>0 on exit.  Otherwise
 * dist_ret[0]=0 and dist_ret[1] = -log(p).
 */
{
	double  CS = 0;
	double  zrmin, zrmax;
	int  j;

	zrmin = -DBL_MAX;
	zrmax = DBL_MAX;
	assert(ds->samples_used > 0);

	for (j=1; j<ds->samples_used; ++j) {
		if (zrmin < result[j]->zmin / result[j]->a)
			zrmin = result[j]->zmin / result[j]->a;
		if (zrmax > result[j]->zmax / result[j]->a)
			zrmax = result[j]->zmax / result[j]->a;
	}
	if (zr < zrmin) {
		dist_ret[0] = 1 + (zrmin-zr);
		dist_ret[1] = 0;
		return;
	}
	if (zr > zrmax) {
		dist_ret[0] = 1 + (zr-zrmax);
		dist_ret[1] = 0;
		return;
	}


	for (j=0; j<ds->samples_used; ++j) {
		CS += distance_interpolate(result[j], zr * result[j]->a);
	}
	dist_ret[0] = 0;
	dist_ret[1] = CS;
}


struct distance_cs *
CS_get_distance(const struct samples *samples, const double *para)
/* Compute the CS test statistic
 *
 * The function returns a 'struct distance_cs ' which the caller must free
 * after use.
 */
{
	struct F_calculator *fc;
	struct distance_cs *result;
	struct quantiles *quantiles = samples->quantiles;
	const double *F;
	double *CS;
	int *df;
	const double CORRECTION = 1e-5;
	double *last_F = NULL, last_level;
	int plus_used, minus_used, S;
	int i, j, N;

	fc = F_new(para);
	N = F_get_N(fc);
	last_F = xnew(double, N+1);

	result = new_distance(N, F_get_z(fc, 0), F_get_z(fc, N), para[p_a]);
	CS = result->CS;
	df = result->df;
	for (i=0; i<=N; ++i)  {
		CS[i] = 0;
		df[i] = -1;
	}

	plus_used = samples->plus_used;
	minus_used = samples->minus_used;
	S = plus_used + minus_used;


	if (plus_used>11) {
		F_start(fc, b_upper);
		F = F_get_F(fc, 0);
		for (i=0; i<=N; ++i) last_F[i] = F[i];
		last_level = 0;
		for (j=0; j<quantiles->m; ++j) {
			double level;
			double emp_q;

			if (quantiles->plus_bounds[j] > 0) {
				F = F_get_F(fc, quantiles->plus_bounds[j]);
			} else {
				F = NULL;
			}

			/* empirical size of bin (% responses) */
			level = 0.001 * quantiles->permille_levels[j];
			emp_q = level - last_level;
			last_level = level;

			for (i=0; i<=N; i++) {
				double pred_q, Ne, Np;
				if (F) {
					pred_q = F[i] - last_F[i];
					last_F[i] = F[i];
				} else {
					pred_q=0;
				}

				Ne = emp_q * plus_used;
				Np = pred_q * S;
				if (Np<CORRECTION) Np = CORRECTION;

				CS[i] += (Np-Ne) * (Np-Ne) / Np;
				df[i]++;
			}
		}
		for (i=0; i<=N; i++) {
			double pred_q, Ne, Np;

			pred_q = 1 - last_F[i];

			Ne = (1-last_level) * plus_used;
			Np = pred_q * S;
			if (Np<CORRECTION) Np = CORRECTION;

			CS[i] += (Np-Ne) * (Np-Ne) / Np;
			df[i]++;
		}
	}

	if (minus_used>11) {
		F_start(fc, b_lower);
		F = F_get_F(fc, 0);
		for (i=0; i<=N; ++i) last_F[i] = F[i];
		last_level = 0;
		for (j=0; j<quantiles->m; ++j) {
			double level;
			double emp_q;

			if (quantiles->minus_bounds[j] > 0) {
				F = F_get_F(fc, quantiles->minus_bounds[j]);
			} else {
				F = NULL;
			}

			/* empirical size of bin (% responses) */
			level = 0.001 * quantiles->permille_levels[j];
			emp_q = level - last_level;
			last_level = level;

			for (i=0; i<=N; i++) {
				double pred_q, Ne, Np;

				if (F) {
					pred_q = last_F[i] - F[i];
					last_F[i] = F[i];
				} else {
					pred_q = 0;
				}

				Ne = emp_q * minus_used;
				Np = pred_q * S;
				if (Np<CORRECTION) Np = CORRECTION;

				CS[i] += (Np-Ne) * (Np-Ne) / Np;
				df[i]++;
			}
		}
		for (i=0; i<=N; i++) {
			double pred_q, Ne, Np;

			pred_q = last_F[i];

			Ne = (1-last_level) * minus_used;
			Np = pred_q * S;
			if (Np<CORRECTION) Np = CORRECTION;

			CS[i] += (Np-Ne) * (Np-Ne) / Np;
			df[i]++;
		}
	}

	xfree(last_F);
	F_delete(fc);

	return  result;
}

void
distance_cs(const struct dataset *ds, const double *x,
			double cs_ret[2], double *zr_ret)
/* Get the 'distance' between theoretical and target distribution.
 *
 * The observed distribution is described by the dataset 'ds', the
 * predicted distribution is described by the parameters 'x'.  The
 * correnspondence between entries of 'x' and the model parameters is
 * encoded in the 'ds->cmds' field.
 *
 * If the parameter 'z' is being optimised, 'zr_ret' must on entry
 * point to an array of size 'ds->z->used'.  If the parameter 'z' is
 * fixed, the value 'zr_ret' must be 'NULL'.
 *
 * If the parameters 'x' are invalid, a positive penalty value is
 * returned in cs_ret[0].  Otherwise cs_ret[0]=0 and the distance
 * between empirical and theoretically predicted distribution function
 * is returned in 'cs_ret[1]' as -log(p).  If 'zr_ret' is non-null,
 * the z-parameter values corresponding to this (minimal) distance are
 * returned in '*zr_ret'.
 */
{
	double  para[p_count], para_zr, penalty;
	struct samples *samples;
	struct distance_cs **result;
	int *z_used;
	int  i, z_idx = -1;

	result = xnew(struct distance_cs *, ds->samples_used);
	for (i=0; i<ds->samples_used; ++i) result[i] = NULL;
	z_used = xnew(int, ds->samples_used);

	para_zr = -1;

	penalty = 0;
	for (i=0; i<ds->cmds_used; ++i) {
		int  arg1 = ds->cmds[i].arg1;
		int  arg2 = ds->cmds[i].arg2;
		switch (ds->cmds[i].cmd) {
		case c_copy_param:
			if (arg1 >= 0) {
				para[arg1] = x[arg2];
			} else {
				assert(zr_ret);
				z_idx = arg2;
			}
			break;
		case c_copy_const:
			if (arg1 >= 0) {
				para[arg1] = ds->consts[arg2];
			} else {
				assert(! zr_ret);
				para_zr = ds->consts[arg2];
			}
			break;
		case c_run:
			samples = ds->samples[arg1];
			penalty += check_bounds(para, NULL, 0);
			if (! zr_ret) {
				double zr = para_zr;
				assert(zr >= 0);
				if (zr - 0.5*para[p_szr] < 0)
					penalty += 1 + 0.5*para[p_szr] - zr;
				if (zr + 0.5*para[p_szr] > 1)
					penalty += 1 + zr + 0.5*para[p_szr] - 1;
			}
			if (penalty > 0) continue;

			result[arg1] = CS_get_distance(samples, para);
			z_used[arg1] = z_idx;
			break;
		}
	}

	if (penalty>0) {
		cs_ret[0] = penalty;
		cs_ret[1] = 0;
	} else if (zr_ret) {
		find_best_CS(result, ds, z_used, cs_ret, zr_ret);
	} else {
		assert(para_zr >= 0);
		find_fixed_CS(result, ds, cs_ret, para_zr);
	}

	xfree(z_used);
	for (i=0; i<ds->samples_used; ++i) {
		if (result[i]) delete_distance(result[i]);
	}
	xfree(result);
}

static void
minimiser(const double *x, double cs_ret[2], void *data)
/* Wrapper to call 'distance' from inside the 'simplex2' function.  */
{
	struct dataset *ds = data;
	double *zr;

	zr = xnew(double, ds->z->used);
	distance_cs(ds, x, cs_ret, zr);
	xfree(zr);
}

static void
print_cs_dist(const double *dist)
{
	if (dist[0] > 0)
		printf("  ... penalty %g\n", dist[0]);
	else
		printf("  ... CS = %g\n", dist[1]);
}

struct dataset_results *
analyse_dataset_cs(struct dataset *ds)
{
	double  *x, *eps, dist[2], *z;
	clock_t  start, stop;
	double cpu_time_used;
	int i;

	for(i=0; i<ds->samples_used; i++) {
		if ((ds->samples[i]->plus_used > 11)
			|| (ds->samples[i]->minus_used > 11))
			continue;
		printf("  Not enough trials in at least one condition for the CS method!\n");
		printf("  (Required: 12 trials for one response in each condition)\n");
		return NULL;
	}

	dataset_find_quantiles(ds, M_QUANTILES, quantile_tab);

	start = clock();
	x = xnew(double, ds->param->used);
	eps = xnew(double, ds->param->used);
	set_precision(ds->precision);
	initialise_parameters(ds, x, eps);
	if (ds->z->used > 0)
		z = xnew(double, ds->z->used);
	else
		z = NULL;

	if (ds->param->used == 0) {
		printf("  nothing to do!\n");
		distance_cs(ds, x, dist, z);
		print_cs_dist(dist);
	} else {
		simplex2(ds->param->used, eps, 0.01, minimiser, x, dist, ds, NULL);
		print_cs_dist(dist);
		simplex2(ds->param->used, eps, 0.001, minimiser, x, dist, ds, NULL);
		print_cs_dist(dist);
		simplex2(ds->param->used, eps, 0.0001, minimiser, x, dist, ds, NULL);
		print_cs_dist(dist);
	}

	stop = clock();
	cpu_time_used = (double)(stop - start) / CLOCKS_PER_SEC;
	xfree(z);

	if (ds->z->used > 0)
		z = xnew(double, ds->z->used);
	else
		z = NULL;
	distance_cs(ds, x, dist, z);

	for (i=0; i<ds->z->used; ++i) {
		printf("  -> %s = %f\n", ds->z->entry[i], z[i]);
	}
	for (i=0; i<ds->param->used; ++i) {
		printf("  -> %s = %f\n", ds->param->entry[i], x[i]);
	}

	//dataset_save_result(ds, x, z, dist[0], dist[1], cpu_time_used);
	//experiment_log(ex, ds, x, z, dist[0], dist[1], cpu_time_used);

	struct dataset_results *resu = xnew(struct dataset_results, 1);
	resu->ds = ds;
	resu->x = x;
	resu->z = z;
	resu->penalty = dist[0];
	resu->fit = exp(-dist[1]);
	resu->time = cpu_time_used;

	//xfree(z);
	//xfree(x);
	xfree(eps);
	return resu;
}
