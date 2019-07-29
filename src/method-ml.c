/* method-ml.c - maximum likelihood method for parameter estimation
 *
 * Copyright (C) 2012  Andreas Voss, Jochen Voss.
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fast-dm.h"


static double
log_likelihood(struct samples *s, const double *para, double *penalty)
/* Returns the sum of log-likelihood across all responses */
{
	double max = s->max, min = s->min;
	double res = 0;
	int i;

	for (i=0; i<s->plus_used; ++i) {
		double dens = (1-para[p_p]) * g_plus(s->plus_data[i], para) + para[p_p] / (max-min);
		if (dens <= 0)
			penalty += 1;
		else
			res += log(dens);
	}
	for (i=0; i<s->minus_used; ++i) {
		double dens = (1-para[p_p]) * g_minus(s->minus_data[i], para)  + para[p_p] / (max-min);
		if (dens <= 0)
			penalty += 1;
		else
			res += log(dens);
	}
	return res;
}

static void
distance_ml(const struct dataset *ds, const double *x,
			double ml_ret[2])
/* Get the negative log-likelihood function for use in the simplex
 * algorithm.
 *
 * The observed distribution is described by the dataset 'ds', the
 * predicted distribution is described by the parameters 'x'.  The
 * correnspondence between entries of 'x' and the model parameters is
 * encoded in the 'ds->cmds' field.
 *
 * If the parameters 'x' are invalid, a positive penalty value is
 * returned in ml_ret[0].  Otherwise ml_ret[0]=0 and the negative
 * log-likelihood is returned in 'ml_ret[1]'.
 */
{
	double  para[p_count + 1], penalty;
	struct samples *samples;
	double  log_lik = 0;
	int  j;

	penalty = 0;
	for (j=0; j<ds->cmds_used; ++j) {
		int  arg1 = ds->cmds[j].arg1;
		int  arg2 = ds->cmds[j].arg2;
		switch (ds->cmds[j].cmd) {
		case c_copy_param:
			if (arg1 >= 0) {
				para[arg1] = x[arg2];
			} else {
				para[p_count] = x[ds->param->used+arg2];
			}
			break;
		case c_copy_const:
			if (arg1 >= 0) {
				para[arg1] = ds->consts[arg2];
			} else {
				para[p_count] = ds->consts[arg2];
			}
			break;
		case c_run:
			samples = ds->samples[arg1];
			penalty += check_bounds(para, samples, 1);
			if (penalty > 0) continue;

			log_lik += log_likelihood(samples, para, &penalty);
			break;
		}
	}

	ml_ret[0] = penalty;
	ml_ret[1] = - log_lik; // we minimise the _negative_ likelihood
}


static void
minimiser(const double *x, double ml_ret[2], void *data)
{
	struct dataset *ds = data;
	distance_ml(ds, x, ml_ret);
}

static void
print_ml_dist(const double *dist)
{
	if (dist[0] > 0)
		printf("  ... penalty %g\n", dist[0]);
	else
		printf("  ... -LL = %g\n", dist[1]);
}

struct dataset_results *
analyse_dataset_ml(struct dataset *ds)
{
	clock_t  start, stop;
	double cpu_time_used;
	double  *x, *eps, dist[2], *z;
	int i;

	for (i=0; i<ds->samples_used; i++) {
		if (ds->samples[i]->plus_used + ds->samples[i]->minus_used >= 10)
			continue;
		printf("  Not enough trials for the ML method!\n");
		printf("  (Required: 10 trials for each condition)\n");
		return NULL;
	}

	start = clock();
	x = xnew(double, ds->param->used + ds->z->used);
	eps = xnew(double, ds->param->used + ds->z->used);
	set_precision(ds->precision);
	initialise_parameters(ds, x, eps);
	for (i=ds->param->used; i<ds->param->used+ds->z->used; ++i) {
		x[i] = 0.5;
		eps[i] = 0.1;
	}

	if (ds->param->used + ds->z->used == 0) {
		minimiser(x, dist, ds);
		printf("  nothing to do!\n");
		print_ml_dist(dist);
	} else {
		simplex2(ds->param->used + ds->z->used, eps, 0.50, minimiser,
				 x, dist, ds, NULL);
		print_ml_dist(dist);
		simplex2(ds->param->used + ds->z->used, eps, 0.01, minimiser,
				 x, dist, ds, NULL);
		print_ml_dist(dist);
		simplex2(ds->param->used + ds->z->used, eps, 0.001, minimiser,
				 x, dist, ds, NULL);
		print_ml_dist(dist);
	}

	stop = clock();
	cpu_time_used = (double)(stop - start) / CLOCKS_PER_SEC;

	z = x + ds->param->used;

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

	//xfree(x);
	xfree(eps);
	return resu;
}
