/* method-ks.c - Kolmogorov-Smirnov method for parameter estimation
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


struct distance_ks {
	int  N;
	double *T;              /* array of distances (length N+1) */
	double  zmin, zmax;     /* absolute z-values for T[0] and T[N] */
	double  a;              /* interval width used when computing T */
};

static double
Q_KS(double lambda)
/* Implement equation (14.3.7) from numerical recipes.  */
{
	double accuracy = 0.001;
	double  c, d, limit, sum;
	int  n, nn;

	/* We want to avoid the random fluctuations (caused by
	 * truncating the sum) for small values of lambda.  The exact
	 * value at lambda=0 is 1.  The code below returns
	 * 1-1.02072e-05 for lambda=0.3 with accuracy=0.0001.  It
	 * returns 1-0.000314775 for lambda=0.35 with accuracy=0.001.
	 * We use linear interpolation between these points.  */
	if (lambda <= 0.3) {
		double a = 1.02072e-05;
		return 1.0 - a*lambda/0.3;
	} else if (lambda <= 0.35) {
		double a = 1.02072e-05;
		double b = 0.000314775;
		return 1-a - (b-a)*(lambda-0.3)/(0.35-0.3);
	}

	c = 2*lambda*lambda;
	d = 1/c;
	limit = d*log(d/(6*accuracy));
	n = 0;
	sum = 0;
	do {
		n += 1;
		nn = n*n;
		sum -= exp(-c*nn);
		n += 1;
		nn = n*n;
		sum += exp(-c*nn);
	} while (n<6 || nn<=limit);
	return  -2*sum;
}

static double
KS_T_to_p(double T, const struct samples *s)
/* Implement equation (14.3.9) from numerical recipes.  */
{
	int n = s->plus_used + s->minus_used;
	double sqrtn = sqrt(n);

	return  Q_KS((sqrtn + 0.12 + 0.11/sqrtn)*T);
}

static struct distance_ks *
new_distance(int N, double zmin, double zmax, double a)
/* Allocate a new structure to hold N+1 values in [zmin,zmax] */
{
	struct distance_ks *d;

	assert(N >= 0);

	d = xnew(struct distance_ks, 1);
	d->N = N;
	d->T = xnew(double, N+1);
	d->zmin = zmin;
	d->zmax = zmax;
	d->a = a;

	return d;
}

static void
delete_distance(struct distance_ks *d)
{
	xfree(d->T);
	xfree(d);
}

static double
distance_interpolate(const struct distance_ks *d, double z)
/* Use linear interpolation to get the KS distance for a given
 * (absolute) value of z.  */
{
	int  N = d->N;
	double  dz, ratio;
	int  step_before;

	if (N == 0)
		return d->T[0];

	dz = (d->zmax - d->zmin) / N;
	step_before = (int)((z - d->zmin) / dz);
	if (step_before >= N)  step_before = N-1;
	ratio = (z - (d->zmin + dz*step_before)) / dz;
	return  (1-ratio) * d->T[step_before]
		+ ratio * d->T[step_before+1];
}

static void
find_best_log_p(struct distance_ks *const*result, const struct dataset *ds,
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
	double  total_log;
	double  penalty;
	int  k;

	total_log = 0;
	penalty = 0;
	for (k=0; k<ds->z->used; k++) {
		double  dz;
		double  best_logp, best_zr;
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

		/* Use sub-sampling of the interval zrmin..zrmax
		 * to find the best z.  */
		n = (int)((zrmax - zrmin) / 0.0001 + 1.5);
		dz = (zrmax - zrmin) / n;
		best_logp = -DBL_MAX;
		best_zr = 0.5;
		for (i=0; i<=n; ++i) {
			double  zr = zrmin + i*dz;
			double  logp = 0;

			for (j=0; j<ds->samples_used; ++j) {
				double  T, p;

				if (z_used[j] != k)  continue;
				T = distance_interpolate(result[j],
										 zr * result[j]->a);
				p = KS_T_to_p(T, ds->samples[j]);
				if (p > 0)
					logp += log(p);
				else logp = -DBL_MAX;
			}
			if (logp > best_logp) {
				best_logp = logp;
				best_zr = zr;
			}
		}
		total_log += best_logp;
		zr_ret[k] = best_zr;
	}

	if (penalty > 0) {
		dist_ret[0] = penalty;
		dist_ret[1] = 0;
	} else {
		dist_ret[0] = 0;
		dist_ret[1] = -total_log;
	}
}

static void
find_fixed_log_p(struct distance_ks *const*result, const struct dataset *ds,
				 double dist_ret[2], double zr)
/* Combine the T-values from 'result' into a common p-value.
 *
 * 'zr' gives the fixed relative z-value.
 *
 * If the parameters are invalid, dist_ret[0]>0 on exit.  Otherwise
 * dist_ret[0]=0 and dist_ret[1] = -log(p).
 */
{
	double  logp;
	double  zrmin, zrmax;
	int  j;

	zrmin = -DBL_MAX;
	zrmax = DBL_MAX;
	assert(ds->samples_used > 0);
	for (j=0; j<ds->samples_used; ++j) {
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

	logp = 0;
	for (j=0; j<ds->samples_used; ++j) {
		double T;

		T = distance_interpolate(result[j], zr * result[j]->a);
		logp += log(KS_T_to_p(T, ds->samples[j]));
	}
	dist_ret[0] = 0;
	dist_ret[1] = -logp;
}

static struct distance_ks *
KS_get_distance(const struct samples *samples, const double *para)
/* Compute the KS test statistic
 *
 * The function returns a 'struct distance_ks' which the caller must free
 * after use.
 */
{
	struct F_calculator *fc;
	struct distance_ks *result;
	double  p_mid, dp;
	double *T;
	int  i, j, N;

	fc = F_new(para);
	N = F_get_N(fc);

	result = new_distance(N, F_get_z(fc, 0), F_get_z(fc, N), para[p_a]);
	T = result->T;
	for (i=0; i<=N; ++i)  T[i] = 0;

	dp = 1.0/(samples->plus_used+samples->minus_used);
	p_mid = samples->minus_used*dp;

	F_start (fc, b_upper);
	for (j=0; j<samples->plus_used; ++j) {
		const double *F = F_get_F(fc, samples->plus_data[j]);
		for (i=0; i<=N; ++i) {
			double p_theo = F[i];
			double dist;

			dist = fabs(p_mid+j*dp-p_theo);
			if (dist > T[i])  T[i] = dist;
			dist = fabs(p_mid+(j+1)*dp-p_theo);
			if (dist > T[i])  T[i] = dist;
		}
	}

	F_start(fc, b_lower);
	for (j=0; j<samples->minus_used; ++j) {
		const double *F = F_get_F(fc, samples->minus_data[j]);
		for (i=0; i<=N; ++i) {
			double p_theo = F[i];
			double dist;

			dist = fabs(p_mid-j*dp-p_theo);
			if (dist > T[i])  T[i] = dist;
			dist = fabs(p_mid-(j+1)*dp-p_theo);
			if (dist > T[i])  T[i] = dist;
		}
	}

	F_delete(fc);
	return  result;
}

static void
distance(const struct dataset *ds, const double *x,
		 double ks_ret[2], double *zr_ret)
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
 * returned in ks_ret[0].  Otherwise ks_ret[0]=0 and the distance
 * between empirical and theoretically predicted distribution function
 * is returned in 'ks_ret[1]' as -log(p).  If 'zr_ret' is non-null,
 * the z-parameter values corresponding to this (minimal) distance are
 * returned in '*zr_ret'.
 */
{
	double  para[p_count], para_zr, penalty;
	struct samples *samples;
	struct distance_ks **result;
	int *z_used;
	int  i, z_idx = -1;

	result = xnew(struct distance_ks *, ds->samples_used);
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

			result[arg1] = KS_get_distance(samples, para);
			z_used[arg1] = z_idx;
			break;
		}
	}

	if (penalty>0) {
		ks_ret[0] = penalty;
		ks_ret[1] = 0;
	} else if (zr_ret) {
		find_best_log_p(result, ds, z_used, ks_ret, zr_ret);
	} else {
		assert(para_zr >= 0);
		find_fixed_log_p(result, ds, ks_ret, para_zr);
	}

	xfree(z_used);
	for (i=0; i<ds->samples_used; ++i) {
		if (result[i]) delete_distance(result[i]);
	}
	xfree(result);
}

static void
minimiser(const double *x, double ks_ret[2], void *data)
/* Wrapper to call 'distance' from inside the 'simplex2' function.  */
{
	struct dataset *ds = data;
	double *zr;

	zr = xnew(double, ds->z->used);
	distance(ds, x, ks_ret, zr);
	xfree(zr);
}

static void
print_ks_dist(const double *dist)
{
	if (dist[0] > 0)
		printf("  ... penalty %g\n", dist[0]);
	else
		printf("  ... p = %g\n", exp(-dist[1]));
}

struct dataset_results *
analyse_dataset_ks(struct dataset *ds)
{
	double  *x, *eps, dist[2], *z;
	clock_t  start, stop;
	double cpu_time_used;
	int i;

	for(i=0; i<ds->samples_used; i++) {
		if (ds->samples[i]->plus_used + ds->samples[i]->minus_used >= 10)
			continue;
		printf("  Not enough trials for the KS method!\n");
		printf("  (Required: 10 trials for each condition)\n");
		return NULL;
	}

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
		distance(ds, x, dist, z);
		print_ks_dist(dist);
	} else {
		simplex2(ds->param->used, eps, 0.01, minimiser, x, dist, ds, NULL);
		print_ks_dist(dist);
		simplex2(ds->param->used, eps, 0.001, minimiser, x, dist, ds, NULL);
		print_ks_dist(dist);
		simplex2(ds->param->used, eps, 0.0001, minimiser, x, dist, ds, NULL);
		print_ks_dist(dist);
	}

	stop = clock();
	cpu_time_used = (double)(stop - start) / CLOCKS_PER_SEC;

	distance(ds, x, dist, z);

	for (i=0; i<ds->z->used; ++i) {
		printf("  -> %s = %f\n", ds->z->entry[i], z[i]);
	}
	for (i=0; i<ds->param->used; ++i) {
		printf("  -> %s = %f\n", ds->param->entry[i], x[i]);
	}

	//dataset_save_result(ds, x, z, dist[0], exp(-dist[1]), cpu_time_used);
	//experiment_log(ex, ds, x, z, dist[0], exp(-dist[1]), cpu_time_used);

	struct dataset_results *resu = xnew(struct dataset_results, 1);
	resu->ds = ds;
	resu->x = x;
	resu->z = z;
	resu->penalty = dist[0];
	resu->fit = exp(-dist[1]);
	resu->time = cpu_time_used;

	//xfree(x);
	//xfree(z);
	xfree(eps);
	return resu;
}
