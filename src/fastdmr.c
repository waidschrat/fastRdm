#include <stdio.h>
#include <stdlib.h>

#include <float.h>
#include <time.h>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "fast-dm.h"
#include "fastdmlib.h"
#include "fastdmr.h"

/**
 * prec = REALSXP precision
 * meth = STRSXP method
 * cond = STRSXP condition headers
 * oparams = VECSXP (list of lists)
 *           > STRSXP name
 *           > STRSXP depends (vector of condition headers)
 *           > REALSXP value (NULL = dont overwrite)
 * ts = REALSXP time array
 * resps = LGLSXP response array
 * conds = STRSXP conditions array
 */
SEXP process_r(SEXP _prec, SEXP _meth, SEXP _cond, SEXP _oparams, SEXP _ts, SEXP _resps, SEXP _conds) {
    int i, n, g;
    struct dataset_results *resu;

    double prec = asReal(_prec); // precision
    const char* methstr = CHAR(asChar(_meth)); // method
    enum method meth = METHOD_KS;
    if (strcmp(methstr, "ks") == 0) {
				meth = METHOD_KS;
    } else if (strcmp(methstr, "ml") == 0) {
				meth = METHOD_ML;
		} else if (strcmp(methstr, "cs") == 0) {
				meth = METHOD_CS;
		}

    // condotion headers
    int cond_length = length(_cond);
    struct set *cond = new_set();

    for (i = 0; i < cond_length; i++) {
        const char *entry = CHAR(STRING_ELT(_cond, i));
        set_item(cond, entry, 1);
    }

    // trial count
    int trialcount = length(_ts);

    double *ts = REAL(_ts); // ts
    int *resps = LOGICAL(_resps); // resps

    // conditions
    const char **conds;
    int conds_length = trialcount * cond_length;
    conds = xnew(const char*, conds_length);

    for (i = 0; i < conds_length; i++) {
        conds[i] = CHAR(STRING_ELT(_conds, i));
    }

    // oparams
    size_t oparam_count = length(_oparams);
    struct param_info *oparams = xnew(struct param_info, oparam_count);

    for (i = 0; i < oparam_count; i++) {
        SEXP op = VECTOR_ELT(_oparams, i);
        int op_l = length(op);

        const char *op_name = NULL;
        struct set *op_depends = new_set();
        double op_value = 0.0; //
        int op_value_set = 0;
        for (n = 0; n<op_l; n++) {
            SEXP x = VECTOR_ELT(op, n);
            if (isNull(x)) continue;
            const char *x_name = CHAR(STRING_ELT(getAttrib(op, R_NamesSymbol), n));

            if (strcmp("name",x_name)==0) {
                op_name = CHAR(asChar(x));

            } else if (strcmp("depends",x_name)==0) {

                for (g = 0; g < cond_length; g++) {
                    const char *entry = CHAR(STRING_ELT(_cond, i));
                    set_item(op_depends, entry, 1);
                }

            } else if (strcmp("value",x_name)==0) {
                op_value = asReal(x);
                op_value_set = 1;
            }
        }
        oparams[i].name = op_name;
        oparams[i].depends = op_depends;
        if (op_value_set) {
            char numbuf[50];
            snprintf(numbuf, 50, "%f", op_value);
            oparams[i].value = _strdup(numbuf); //double to string
        } else {
            oparams[i].value = NULL;
        }
        //printf("name: %s; value: %f; depends: ",op_name, op_value);
        //print_set(op_depends);
    }


    resu = process(prec, meth, cond, oparams, oparam_count,
                   trialcount, ts, resps, conds);


    const int resu_length = 4;
    SEXP result = PROTECT(allocVector(VECSXP, resu_length));
    SEXP list_names = PROTECT(allocVector(STRSXP, resu_length));

    const int resu_length_par = resu->ds->z->used + resu->ds->param->used;
    SEXP result_par = PROTECT(allocVector(VECSXP, resu_length_par));
    SEXP list_names_par = PROTECT(allocVector(STRSXP, resu_length_par));

    int off = 0;
    SEXP myreal;
    for (i=0; i<resu->ds->z->used; ++i) {
		 //Rprintf ("%s = %f\n",
		 //		 resu->ds->z->entry[i], resu->z[i]);

        SET_STRING_ELT(list_names_par, off, mkChar(resu->ds->z->entry[i]));
        myreal = PROTECT(allocVector(REALSXP, 1));
        REAL(myreal)[0] = resu->z[i];
        SET_VECTOR_ELT(result_par, off++, myreal);
    }
	for (i=0; i<resu->ds->param->used; ++i) {
		// Rprintf ("%s = %f\n",
		// 		 resu->ds->param->entry[i], resu->x[i]);
        SET_STRING_ELT(list_names_par, off, mkChar(resu->ds->param->entry[i]));
        myreal = PROTECT(allocVector(REALSXP, 1));
        REAL(myreal)[0] = resu->x[i];
        SET_VECTOR_ELT(result_par, off++, myreal);
    }

	SET_STRING_ELT(list_names, 0, mkChar("par") );
  SET_VECTOR_ELT(result, 0, result_par);

	// printf("precision = %f\n", resu->ds->precision);
  //  SET_STRING_ELT(list_names, off, mkChar("precision"));
  //  myreal = PROTECT(allocVector(REALSXP, 1));
  //  REAL(myreal)[0] = resu->ds->precision;
  //  SET_VECTOR_ELT(result, off++, myreal);
	// printf("method = %s\n", method_name(resu->ds->method));
	// printf("penalty = %f\n", resu->penalty);
    SET_STRING_ELT(list_names, 1, mkChar("penalty") );
    myreal = PROTECT(allocVector(REALSXP, 1));
    REAL(myreal)[0] = resu->penalty;
    SET_VECTOR_ELT(result, 1, myreal);
	// printf("fit index = %f\n", resu->fit);
    SET_STRING_ELT(list_names, 2, mkChar("fit") );
    myreal = PROTECT(allocVector(REALSXP, 1));
    REAL(myreal)[0] = resu->fit;
    SET_VECTOR_ELT(result, 2, myreal);
	// printf("time = %f\n", resu->time);
    SET_STRING_ELT(list_names, 3, mkChar("time") );
    myreal = PROTECT(allocVector(REALSXP, 1));
    REAL(myreal)[0] = resu->time;
    SET_VECTOR_ELT(result, 3, myreal);

    setAttrib(result, R_NamesSymbol, list_names);
    setAttrib(result_par, R_NamesSymbol, list_names_par);

    //delete_dataset_results(resu);

    UNPROTECT(4 + off);

    return result;
}
