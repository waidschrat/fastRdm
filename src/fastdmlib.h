
struct dataset_results *process(double precision, enum method meth, struct set *cond, struct param_info *oparams, size_t oparam_count,
                    int trialcount, double ts[], int resps[], const char *conds[]);

//struct experiment *new_experiment_2(double precision, enum method meth, struct set *cond, struct param_info *params);
//struct dataset *build_dataset (double precision, enum method meth, struct set *cond, struct param_info *params, struct experiment *ex, int trialcount, double ts[], int resps[], const char *conds[]);

struct dataset *
  build_dataset (double precision, enum method meth, struct set *cond, struct param_info *params, int trialcount, double ts[], int resps[], const char *conds[]);


//void experiment_print_2 (const struct experiment *ex);

struct param_info *params_override_default(struct param_info oparams[], size_t oparam_count);

void print_dataset_results(struct dataset_results *resu);
void print_params(struct param_info *params);
