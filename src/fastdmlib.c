#include <stdio.h>
#include <stdlib.h>

#include <float.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#include "fast-dm.h"
#include "fastdmlib.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) < (Y) ? (Y) : (X))


static int
  string_to_double (const char *str, double *res)
  /* return 1 on success and 0 on error */
  {
    char *tail;
    double  x;

    x = strtod (str, &tail);
    if (tail == str)  return 0;
    if (*tail != '\0')  return 0;
    *res = x;
    return 1;
  }

static void
  str_add (char **s1, const char *s2, size_t *used, size_t *alloc)
  {
    size_t  l = strlen(s2);

    if (*used+l+1 >= *alloc) {
      *alloc += 64;
      *s1 = xrenew(char, *s1, *alloc);
    }
    memcpy (*s1+*used, s2, l);
    *used += l;
    (*s1)[*used] = '\0';
  }


struct dataset_results *process(double precision, enum method meth, struct set *cond, struct param_info oparams[], size_t oparam_count,
                                int trialcount, double ts[], int resps[], const char *conds[]) {

  struct param_info* params;
  //struct experiment *ex;
  struct dataset *ds;

  params = params_override_default(oparams, oparam_count);
  //ex = new_experiment_2(precision, METHOD_KS, cond, params);
  //ds = experiment_get_dataset_2 (ex, trialcount, ts, resps, conds);

  ds = build_dataset( precision, meth, cond, params, trialcount, ts, resps, conds);

  struct dataset_results *resu;
  switch (ds->method) {
  case METHOD_KS:
    resu = analyse_dataset_ks(ds);
    break;
  case METHOD_ML:
    resu = analyse_dataset_ml(ds);
    break;
  case METHOD_CS:
    resu = analyse_dataset_cs(ds);
    break;
  default:
    return NULL;
  }

  //delete_dataset(ds); // TODO
  return resu;
}

struct param_info *params_override_default(struct param_info oparams[], size_t oparam_count) {
  int i, n;

  struct param_info *params;

  params = xmalloc(sizeof (default_params));
  memcpy (params, default_params, sizeof (default_params));
  for (i = 0; params[i].name; ++i) {
    int overridden = 0;
    for (n = 0; n < oparam_count; n++) {
      const char* oname = oparams[n].name;

      if (!oname)
        continue;

      if (strcmp(params[i].name, oname) == 0) {
        params[i].name = xstrdup(oname);
        //printf("override param: [%s]\n",params[i].name);
        params[i].value = oparams[n].value != NULL ? xstrdup(oparams[n].value) : NULL;
        params[i].depends = oparams[n].depends ? cpy_set(oparams[n].depends) : new_set();
        overridden = 1;
      }
    }
    if (!overridden)
      params[i].depends = new_set();
  }
  return params;
}




/* Add the necessary commands to 'ds' to initalise parameter 'param'.
 *
 * The dictionary 'condv' maps parameter names to experimental conditions
 * where required.
 */
void
  dataset_init_param (struct dataset *ds, struct param_info *param,
                      const struct dict *condv)
  {
    size_t  name_used, name_alloc;
    char *name;
    int  j;

    name_alloc = 64;
    name = xnew(char, name_alloc);
    name_used = 0;
    *name = '\0';

    str_add (&name, param->name, &name_used, &name_alloc);
    for (j=0; j<param->depends->used; ++j) {
      const char *cv;

      cv = dict_lookup (condv, param->depends->item[j]);
      str_add (&name, "_", &name_used, &name_alloc);
      str_add (&name, cv, &name_used, &name_alloc);
    }

    if (param->value) {
      double  x;
      if (! string_to_double (param->value, &x)) {
        fprintf (stderr, "invalid value '%s' for '%s'\n",
                 param->value, name);
        exit (1);
      }
      dataset_add_cmd (ds, c_copy_const, param->idx,
                       dataset_add_const (ds, x));
    } else {
      if (param->idx >= 0) {
        dataset_add_cmd (ds, c_copy_param, param->idx,
                         dataset_add_param (ds, name));
      } else {
        dataset_add_cmd (ds, c_copy_param, param->idx,
                         dataset_add_z (ds, name));
      }
    }
    xfree(name);
  }

/* Build dataset / samples directly.
 */
struct dataset *
  build_dataset (double precision, enum method meth, struct set *cond, struct param_info *params, int trialcount, double ts[], int resps[], const char *conds[])
  {
    struct dataset *ds;
    int  i, g;

    ds = new_dataset ();
    ds->precision = precision;
    ds->method = meth;

    // not used:
    ds->fname = NULL;//"[NONAME]";
    ds->logname = NULL;//"[NORESULTS]"; // results file (NULL = dont write)
    ds->key = NULL;


    struct dict *condv;
    size_t  sample_name_used, sample_name_alloc;
    char *sample_name;

    condv = new_dict ();

    sample_name_alloc = 80;
    sample_name = xnew(char, sample_name_alloc);
    sample_name_used = 1;
    sample_name[0] = '\0';

    for (g = 0; g < trialcount; g++) {
      struct samples *samples;
      double  t = -1;
      long int  resp = -1;
      int  idx;

      dict_clear (condv);

      t = ts[g];
      resp = resps[g]; // this has to be 0 or 1

      const int condcount = cond->used;
      for (i = 0; i < cond->used; i++) {
        dict_add (condv, cond->item[i], conds[g + i * condcount] );
      }

      /* get the sample set name */
      sample_name_used = 0;
      sample_name[0] = '\0';
      str_add (&sample_name, "",
                 &sample_name_used, &sample_name_alloc);
      for (i=0; i<cond->used; ++i) {
        const char *val;

        val = dict_lookup (condv, cond->item[i]);
        str_add (&sample_name, "_",
                   &sample_name_used, &sample_name_alloc);
        str_add (&sample_name, val,
                   &sample_name_used, &sample_name_alloc);
      }

      /* create new sample sets as needed */
      idx = dataset_samples_idx (ds, sample_name, 0);
      if (idx < 0) {
        idx = dataset_samples_idx (ds, sample_name, 1);
        for (i=0; params[i].name; ++i) {
          dataset_init_param (ds, params+i, condv);
        }
        dataset_add_cmd (ds, c_run, idx, 0);
      }
      samples = ds->samples[idx];

      /* register the sample data */
      assert(resp >= 0);
      samples_add_sample (samples, t, (int)resp);
    }

    xfree(sample_name);
    delete_dict(condv);

    for (i=0; i<ds->samples_used; ++i) {
      struct samples *samples = ds->samples[i];
      samples_sort(samples);
      int n_plus = samples->plus_used;
      double first_plus = n_plus > 0 ? samples->plus_data[0] : DBL_MAX;
      double last_plus = n_plus > 0 ? samples->plus_data[n_plus-1] : 0.0;
      int n_minus = samples->minus_used;
      double first_minus = n_minus > 0 ? samples->minus_data[0] : DBL_MAX;
      double last_minus = n_minus > 0 ? samples->minus_data[n_minus-1] : 0.0;

      samples->min = MIN(first_plus, first_minus);
      samples->max = MAX(last_plus, last_minus);
    }

    return  ds;
  }


void print_dataset_results(struct dataset_results *resu) {
  int i;
  for (i=0; i<resu->ds->z->used; ++i)
    printf ("%s = %f\n",
            resu->ds->z->entry[i], resu->z[i]);

  for (i=0; i<resu->ds->param->used; ++i)
    printf ("%s = %f\n",
            resu->ds->param->entry[i], resu->x[i]);

  printf ("precision = %f\n", resu->ds->precision);
  printf ("method = %s\n", method_name(resu->ds->method));
  printf ("penalty = %f\n", resu->penalty);
  printf ("fit index = %f\n", resu->fit);
  printf ("time = %f\n", resu->time);
}

void print_params(struct param_info *params) {
  printf("name: %s; value: %s; depends: ",params->name,params->value);
  print_set(params->depends);
}

void delete_dataset_results(struct dataset_results *resu) {
  delete_dataset(resu->ds);
  xfree(resu->x);
  xfree(resu->z);
  xfree(resu);
}
