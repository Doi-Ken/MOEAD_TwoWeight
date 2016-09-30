#ifndef MOEADDEF
#define MOEADDEF

void t_neighborhood(int **B, double **lambda, int N, int T, int ob);

double a_g_t_e(double *fit, double *lambda, double *z, int ob);
double g_t_e(double *fit, double *lambda, double *z, int ob);
double weight_sum(double *fit, double *lambda, int ob);

double p_b_i_d2(double *fit, double *lambda, double *z, int ob, double PENALTY);
double p_b_i_min_d2(double *fit, double *lambda, double *z, int ob, double PENALTY);


double p_b_i(double *fit, double *lambda, double *z, int ob, double PENALTY);
double p_b_i_double(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita);
double p_b_i_para0(double *fit, double *lambda, double *z, int ob, double shita);
double quadratic_p_b_i(double *fit, double *lambda, double *z, int ob, double PENALTY);
double p_b_i_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita);

double p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY);
double p_b_i_min_double(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita);
double p_b_i_min_para0(double *fit, double *lambda, double *z, int ob, double shita);
double quadratic_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY);
double p_b_i_min_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita);


double inverted_p_b_i(double *fit, double *lambda, double *z, int ob, double PENALTY);

double inverted_p_b_i_max(double *fit, double *lambda, double *z, int ob, double PENALTY);
double quadratic_inverted_p_b_i_max(double *fit, double *lambda, double *z, int ob, double PENALTY);
double inverted_p_b_i_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita);

double inverted_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY);
double quadratic_inverted_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY);
double inverted_p_b_i_min_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita);



void max_z(double *fit, double *z, int ob, double ALPHA, double *fitmax);
void min_z(double *fit, double *z, int ob, double ALPHA, double *fitmin);

bool function_valuation(char *problem,char *scalarizing_function,
	double *original_objective, double *children_objective,
	bool improve_lambda2, double *lambda, double *lambda2,
	double *maxz, double *minz, 
	int ob, 
	double PENALTY, double PENALTY2, double shita, char *max_min);

#endif