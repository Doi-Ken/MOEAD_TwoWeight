#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<vector>
#include<list>
#include<queue>
#include<stack>
#include<fstream>
#include<iomanip>
#include<ctime>
#include<cmath>
#include<direct.h>
//#include"mersenne.h"
#include "MT.h"
#include"wfg.h"
#include"moead.h"


using namespace std;

#define X_MIN 0.0
#define X_MAX 1.0
#define EPS 0.00000000001



void SBXcrossover(double *par_x1, double *par_x2, double **off_x, int off_p_n, int item, double CO, double SBX_ETA){
	int off_p_n1 = off_p_n;
	int off_p_n2 = off_p_n + 1;

	double x1 = 0.0;
	double x2 = 0.0;
	double beta1 = 1.0;
	double beta2 = 1.0;
	double alpha1 = 2.0;
	double alpha2 = 2.0;
	double beta_q_1 = 1.0;
	double beta_q_2 = 1.0;
	double c1, c2;
	double u = nextDoubleIE();

	if (nextDoubleIE() < CO){
		for (int n = 0; n < item; n++){
			if (nextDoubleIE() < 0.5){
				if (abs(par_x2[n] - par_x1[n]) > EPS){
					x1 = 0.0;
					x2 = 0.0;
					beta1 = 1.0;
					beta2 = 1.0;
					alpha1 = 2.0;
					alpha2 = 2.0;
					beta_q_1 = 1.0;
					beta_q_2 = 1.0;

					u = nextDoubleIE();
					if (par_x1[n] > par_x2[n]){
						x2 = par_x1[n];
						x1 = par_x2[n];
					}
					else if (par_x1[n] < par_x2[n]){
						x2 = par_x2[n];
						x1 = par_x1[n];
					}


					beta1 += 2.0 * (x1 - X_MIN) / (x2 - x1);
					beta2 += 2.0 * (X_MAX - x2) / (x2 - x1);

					alpha1 -= pow(beta1, -1.0 * (SBX_ETA + 1.0));
					alpha2 -= pow(beta2, -1.0 * (SBX_ETA + 1.0));

					if (u <= 1.0 / alpha1){
						beta_q_1 = pow(u * alpha1, (1.0 / (SBX_ETA + 1.0)));
					}
					else{
						beta_q_1 = pow(1.0 / (2.0 - u * alpha1), (1.0 / (SBX_ETA + 1.0)));
					}

					if (u <= 1.0 / alpha2){
						beta_q_2 = pow(u * alpha2, (1.0 / (SBX_ETA + 1.0)));
					}
					else{
						beta_q_2 = pow(1.0 / (2.0 - u * alpha2), (1.0 / (SBX_ETA + 1.0)));
					}

					c1 = 0.5 *((x1 + x2) - beta_q_1 * (x2 - x1));
					c2 = 0.5 *((x1 + x2) + beta_q_2 * (x2 - x1));

					if (c1 < X_MIN) c1 = X_MIN;
					if (c1 > X_MAX) c1 = X_MAX;

					if (c2 < X_MIN) c2 = X_MIN;
					if (c2 > X_MAX) c2 = X_MAX;
					
					
					//improvement
					/*off_x[off_p_n1][n] = c1;
					off_x[off_p_n2][n] = c2;
					*/
					
				
					if (nextDoubleIE() < 0.5){
						off_x[off_p_n1][n] = c1;
						off_x[off_p_n2][n] = c2;

					}
					else{
						off_x[off_p_n1][n] = c2;
						off_x[off_p_n2][n] = c1;
					}
					//improvement

					
				}
				else{
					off_x[off_p_n1][n] = par_x1[n];
					off_x[off_p_n2][n] = par_x2[n];
				}
			}
			else{
				off_x[off_p_n1][n] = par_x1[n];
				off_x[off_p_n2][n] = par_x2[n];
			}
		}

	}
	else{
		for (int n = 0; n < item; n++){
			off_x[off_p_n1][n] = par_x1[n];
			off_x[off_p_n2][n] = par_x2[n];
		}

	}


}

void polynomial_mutation(double *off_x, int item, double MUT, double MUT_ETA){

	double r = 0.0;
	double shita1 = 0.0;
	double shita2 = 0.0;
	double shitaq = 0.0;
	for (int i = 0; i < item; i++){
		r = nextDoubleIE();
		if (r <= MUT){
			shita1 = (off_x[i] - X_MIN) / (X_MAX - X_MIN);
			shita2 = (X_MAX - off_x[i]) / (X_MAX - X_MIN);

			r = nextDoubleIE();

			if (r <= 0.5){
				shitaq = pow(((2.0 * r) + (1.0 - 2.0 * r) * pow(1.0 - shita1, MUT_ETA + 1.0)), 1.0 / (MUT_ETA + 1.0)) - 1.0;

				off_x[i] += shitaq * (X_MAX - X_MIN);

			}
			else{
				shitaq = 1.0 - pow(2.0 * (1.0 - r) + 2.0 * (r - 0.5) * pow(1.0 - shita2, MUT_ETA + 1.0), 1.0 / (MUT_ETA + 1.0));
				off_x[i] += shitaq * (X_MAX - X_MIN);
			}
			if (off_x[i] < X_MIN){
				off_x[i] = X_MIN;
			}
			if (off_x[i] > X_MAX){
				off_x[i] = X_MAX;
			}
		}

	}
}


//continuous optimization
void make_new_pop(double *x1, double *x2, double off_x[], int item, double MUT, double CO,double MUT_ETA, double SBX_ETA){

	double **temp_off_x = new double*[2];

	temp_off_x[0] = new double[item];
	temp_off_x[1] = new double[item];

	SBXcrossover(x1, x2, temp_off_x, 0, item, CO, SBX_ETA);
	if (nextDoubleIE() < 0.5){
		for (int i = 0; i < item; i++){
			off_x[i] = temp_off_x[0][i];
		}
	}
	else{
		for (int i = 0; i < item; i++){
			off_x[i] = temp_off_x[1][i];
		}
	}

	polynomial_mutation(off_x, item, MUT, MUT_ETA);

}


void make_new_pop(int *x1, int *x2, int off_x[], int item, double MUT, double CO){

	int *temp = new int[item];
	int **x = new int *[2];
	x[0] = new int[item];
	x[1] = new int[item];

	for (int i = 0; i < item; i++){
		x[0][i] = x1[i];
	}
	for (int i = 0; i < item; i++){
		x[1][i] = x2[i];
	}

	if (nextDoubleIE() < CO){
		for (int i = 0; i < item; i++){
			temp[i] = genrand_int31() % 2;
		}
		/*
		for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() <= MUT){
		off_x[i] = (off_x[i] + 1) % 2;
		}
		}*/
	}
	else{
		int k = genrand_int31() % 2;
		for (int i = 0; i < item; i++){
			temp[i] = k;
			/*	off_x[i] = x[temp[i]][i]; */
		}
	}

	for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() < MUT){
			off_x[i] = (off_x[i] + 1) % 2;
		}
	}
	delete[] temp;
	for (int n = 0; n < 2; n++){
		delete[] x[n];
	}
	delete[] x;
}

void make_new_pop_20160301(int *x1, int *x2, int off_x[], int item, double MUT, double CO, int ob){

	int *temp = new int[item];
	int **x = new int *[2];
	x[0] = new int[item];
	x[1] = new int[item];

	for (int i = 0; i < item; i++){
		x[0][i] = x1[i];
	}
	for (int i = 0; i < item; i++){
		x[1][i] = x2[i];
	}

	if (nextDoubleIE() < CO){
		for (int i = 0; i < item; i++){
			temp[i] = genrand_int31() % 2;
		}
		/*
		for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() <= MUT){
		off_x[i] = (off_x[i] + 1) % 2;
		}
		}*/
	}
	else{
		int k = genrand_int31() % 2;
		for (int i = 0; i < item; i++){
			temp[i] = k;
			/*	off_x[i] = x[temp[i]][i]; */
		}
	}

	for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() < MUT){
			int t = off_x[i];
			do{
				off_x[i] = genrand_int31() % ob + 1;
			} while (t == off_x[i]);
		}
	}
	delete[] temp;
	for (int n = 0; n < 2; n++){
		delete[] x[n];
	}
	delete[] x;
}


bool dominate(double *fita, double *fitb, int ob, int SIGN){

	for (int o = 0; o < ob; o++){
		if (SIGN * fita[o] < SIGN * fitb[o]){
			return false;
		}
	}
	for (int o = 0; o < ob; o++){
		if (SIGN * fita[o] > SIGN * fitb[o]){
			return true;
		}
	}
	return false;
}

void non_dominated_set(double **f, double **non_dom_set, int &ndsize, int M, int NPOP, int SIGN){
	int *same_solution = new int[NPOP];
	int *flag_dominated = new int[NPOP];
	int *same_group = new int[NPOP];
	for (int n = 0; n < NPOP; n++){
		same_solution[n] = 0;
		flag_dominated[n] = 0;
		same_group[n] = 0;
	}

	int same_group_count = 1;
	for (int b = 0; b < NPOP; b++){
		int flag = 0;
		for (int a = 0; a < NPOP; a++){
			if (a != b){
				if (dominate(f[a], f[b], M, SIGN)){
					flag_dominated[b] = 1;
				}
				if (same_solution[b] == 0 || same_solution[b] == same_group_count){
					int o = 0;
					for (; o < M; o++){
						if (f[a][o] != f[b][o]){
							break;
						}
					}
					if (o == M){
						same_solution[b] = same_group_count;
						same_solution[a] = same_group_count;
						flag = 1;
					}
				}
			}
		}
		if (flag == 1){
			same_group_count++;
		}
	}

	int nds = 0;
	for (int n = 0; n < NPOP; n++){
		if (flag_dominated[n] == 0 && same_solution[n] == 0){
			for (int o = 0; o < M; o++){
				non_dom_set[nds][o] = f[n][o];
			}
			nds++;

		}
		else if (flag_dominated[n] == 0 && same_solution[n] > 0){
			if (same_group[same_solution[n]] == 0){
				for (int o = 0; o < M; o++){
					non_dom_set[nds][o] = f[n][o];
				}
				same_group[same_solution[n]] = 1;
				nds++;
			}

		}
	}
	ndsize = nds;

	delete [] same_solution;
	delete [] flag_dominated;
	delete [] same_group;

}


int main(int argc, char *argv[]){


	bool combination = false;
	
	if (strcmp(argv[1], "-num") != 0){
		cout << "file num " << endl;
		exit(1);
	}
	
	if (argv[2][0] == '-'){
		cout << "experimental number" << endl;
		exit(1);
	}
	init_genrand(atoi(argv[2]));

	//knapsack
	//knapsack_file_read(); //file_read
	//sorting_profit_per_weight(); //q_j sorting
	//repair_output(); // Step1: repair_output
	//check_input_file_output(x_check); //Step2: knapsack_check
	//for (int o = 0; o < ob; o++){
	//	cout << max_weight[o] << endl;
	//}

	//input
	//ifstream fl("lambda_data_m2_h199_n200.txt");
	int EPN = 6000;

	for (int i = 1; i < argc; i++){
		if (i % 2){
			if (argv[i][0] != '-'){
				cout << "format is disavailable" << endl;
				cout << "-... value -... value -... value ........" << endl;
				exit(1);
			}
		}
	}

	char problem[50] = "wfg1";

	int read_ind = 3;
	if (strcmp(argv[read_ind], "-pro") != 0){
		cout << "default wfg1" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "knapsack") == 0 ||
			strcmp(argv[read_ind], "wfg1") == 0 || strcmp(argv[read_ind], "wfg2") == 0 ||
			strcmp(argv[read_ind], "wfg3") == 0 || strcmp(argv[read_ind], "wfg4") == 0 ||
			strcmp(argv[read_ind], "wfg5") == 0 || strcmp(argv[read_ind], "wfg6") == 0 ||
			strcmp(argv[read_ind], "wfg7") == 0 || strcmp(argv[read_ind], "wfg8") == 0 ||
			strcmp(argv[read_ind], "wfg9") == 0 ||
			strcmp(argv[read_ind], "dtlz1") == 0 || strcmp(argv[read_ind], "dtlz2") == 0 ||
			strcmp(argv[read_ind], "dtlz3") == 0 || strcmp(argv[read_ind], "dtlz4") == 0 || 
			strcmp(argv[read_ind], "dtlz7") == 0 ||
			strcmp(argv[read_ind], "zdt1") == 0 || strcmp(argv[read_ind], "zdt2") == 0 ||
			strcmp(argv[read_ind], "zdt3") == 0 || strcmp(argv[read_ind], "zdt4") == 0 ||
			strcmp(argv[read_ind], "mtest_dtlz") == 0 || strcmp(argv[read_ind], "mtest_wfg") == 0 ||
			strcmp(argv[read_ind], "test20160301") == 0 || strcmp(argv[read_ind], "test20160303") == 0){
			strcpy(problem, argv[read_ind]);
		}
		else{
			cout << "there is no problem." << endl;
			exit(1);
		}

	}


	int ob = 10;

	read_ind++;
	if (strcmp(argv[read_ind], "-obj") != 0){
		cout << "default 10obj" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		ob = atoi(argv[read_ind]);
	}

	int POPSIZE1 = 220;
	int MOEAD_H1 = 3;

	read_ind++;
	if (strcmp(argv[read_ind], "-pop") != 0){
		cout << "default pop size 220" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		POPSIZE1 = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-dev") != 0){
		cout << "default devision 3" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		MOEAD_H1 = atoi(argv[read_ind]);
	}


	int T = POPSIZE1 / 10;

	read_ind++;
	if (strcmp(argv[read_ind], "-nei") != 0){
		cout << "default neighborsize popsize / 10" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		T = atoi(argv[read_ind]);
	}

	char function[100] = "tchebycheff";
	

	read_ind++;
	if (strcmp(argv[read_ind], "-fun") != 0){
		cout << "default function tchebycheff" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		strcpy(function, argv[read_ind]);

	}



	double ALPHA_MAX = 1.1;
	double ALPHA_MIN = 1.0;

	read_ind++;
	if (strcmp(argv[read_ind], "-amax") != 0){
		cout << "default maxalpha 1.1" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		ALPHA_MAX = atof(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-amin") != 0){
		cout << "default minalpha 1.0" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		ALPHA_MIN = atof(argv[read_ind]);
	}

	double PENALTY = 0.1;
	double PENALTY2 = 50.0;

	read_ind++;
	if (strcmp(argv[read_ind], "-pena") != 0){
		cout << "default penalty 0.1" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		PENALTY = atof(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-pena2") != 0){
		cout << "default penalty2 50.0" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		PENALTY2 = atof(argv[read_ind]);
	}

	double shita = M_PI / (2.0 * MOEAD_H1);
	read_ind++;
	if (strcmp(argv[read_ind], "-shita") != 0){
		cout << "default shita pi / (2.0 * div)" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		shita = atof(argv[read_ind]);
	}

	int GEN = 10000;
	int VALNUM = 40000;

	read_ind++;
	if (strcmp(argv[read_ind], "-valnum") != 0){
		cout << "default valnum 40000" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		VALNUM = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-gen") != 0){
		cout << "default generation 10000" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		GEN = atoi(argv[read_ind]);
	}


	double SBX_ETA = 15.0;
	double MUT_ETA = 20.0;
	double CO = 0.8;
	double mut = 1.0;

	read_ind++;
	if (strcmp(argv[read_ind], "-sbxeta") != 0){
		cout << "default sbxeta 15.0" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		SBX_ETA = atof(argv[read_ind]);
	}
	
	read_ind++;
	if (strcmp(argv[read_ind], "-muteta") != 0){
		cout << "default muteta 20.0" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		MUT_ETA = atof(argv[read_ind]);
	}
	
	read_ind++;
	if (strcmp(argv[read_ind], "-co") != 0){
		cout << "default crossover probability" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		CO = atof(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-mut") != 0){
		cout << "default mutation probability mut(1.0) / item" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		mut = atof(argv[read_ind]);
	}

	int k_dtlz = 10;
	int k_fac = 2;
	int l_fac = 10;

	read_ind++;
	if (strcmp(argv[read_ind], "-k_fac") != 0){
		cout << "default k_factor 2" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		k_fac = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-l_fac") != 0){
		cout << "default l_factor 10" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		l_fac = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-k_dtlz") != 0){
		cout << "default k_dtlz 10" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		k_dtlz = atoi(argv[read_ind]);
	}

	bool improve_lambda = false;
	read_ind++;
	if (strcmp(argv[read_ind], "-lam") != 0){
		cout << "default improve_lambda false" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "true") == 0){
			improve_lambda = true;
		}

	}

	char max_min[10] = "min";
	read_ind++;
	if (strcmp(argv[read_ind], "-maxmin") != 0){
		cout << "default improve_lambda false" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "max") == 0){
			strcpy(max_min, "max");
		}
	}

	bool ref_non_dom_flag = false;
	read_ind++;
	if (strcmp(argv[read_ind], "-refnondom") != 0){
		cout << "default improve_lambda false" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "true") == 0){
			ref_non_dom_flag = true;
		}
	}
	
	bool normalize_flag = false;
	read_ind++;
	if (strcmp(argv[read_ind], "-normalize") != 0){
		cout << "default improve_lambda false" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "true") == 0){
			normalize_flag = true;
		}
	}

	int POPSIZE2 = 0;
	int MOEAD_H2 = 1;

	read_ind++;
	if (strcmp(argv[read_ind], "-pop2") != 0){
		cout << "default pop size 220" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		POPSIZE2 = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-dev2") != 0){
		cout << "default devision 3" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		MOEAD_H2 = atoi(argv[read_ind]);
	}

	int point_num = 5;
	double **P = new double *[point_num];
	for (int p = 0; p < point_num; p++){
		P[p] = new double[ob];
	}

	

	//bool switch_pbi = false;
	//read_ind++;
	//if (strcmp(argv[read_ind], "-sw") != 0){
	//	cout << "default switch_pbi false" << endl;
	//	read_ind--;
	//}
	//else{
	//	read_ind++;
	//	if (strcmp(argv[read_ind], "true") == 0){
	//		switch_pbi = true;
	//	}
	//}

	//double switch_rate = 0.7;
	//read_ind++;
	//if (strcmp(argv[read_ind], "-swr") != 0){
	//	cout << "default switch_rate 0.7" << endl;
	//	read_ind--;
	//}
	//else{
	//	read_ind++;
	//	switch_rate = atof(argv[read_ind]);
	//}


	int k_val = k_fac;
	int l_val = l_fac;
	int item = k_val + l_val;


	if (strcmp("dtlz1", problem) == 0 || strcmp("dtlz2", problem) == 0 ||
		strcmp("dtlz3", problem) == 0 || strcmp("dtlz4", problem) == 0 || 
		strcmp("dtlz7", problem) == 0){
		item = ob + k_dtlz - 1;
	}
	else if (strcmp("knapsack", problem) == 0){
		combination = true;

		char file_name[50];
		strcpy(file_name, problem);
		strcat(file_name, "_2_500");
		if (ob != 2){
			char ch[5];
			sprintf(ch, "%d", ob);
			strcat(file_name, "to");
			strcat(file_name, ch);
		}
		strcat(file_name, ".txt");
		knapsack_file_read(file_name, item, ob); //file_read
		sorting_profit_per_weight(item, ob); //q_j sorting
		repair_output(item, ob);
		int *cheeck = new int[item];
		for (int i = 0; i < item; i++){
			cheeck[i] = 1;
		}
		check_input_file_output(cheeck, item, ob);

	}
	else if (strcmp("zdt1", problem) == 0 || strcmp("zdt2", problem) == 0 ||
		strcmp("zdt3", problem) == 0 || strcmp("zdt4", problem) == 0){
		k_dtlz = 30;
		if (strcmp("zdt4", problem) == 0){
			k_dtlz = 10;
		}
		item = k_dtlz;
	}
	else if (strcmp("mtest_wfg", problem) == 0 || strcmp("mtest_dtlz", problem) == 0){
		point_num = 5;
		*P = new double[point_num];
		for (int i = 0; i < point_num; i++){
			P[i] = new double[ob];
		}
		ifstream fcp("cp8.txt");
		for (int p = 0; p < point_num; p++){
			for (int o = 0; o < ob; o++){
				fcp >> P[p][o];
			}
		}
		fcp.close();
		item = k_dtlz;
	}
	else if (strcmp("test20160301", problem) == 0){
		item = k_dtlz;
		combination = true;

	}
	else if (strcmp("test20160303", problem) == 0){
		item = k_dtlz;
		combination = true;

	}
	else{
		*P = new double[point_num];
		for (int i = 0; i < point_num; i++){
			P[i] = new double[ob];
			for (int o = 0; o < ob; o++){
				P[i][o] = 0.0;
			}
		}

	}
	
	cout << "number: " << atoi(argv[2]) << endl;
	cout << "problem: " << problem << ", " << ob << "ob" << endl;
	cout << "popsize: " << POPSIZE1 << ", devision: " << MOEAD_H1 << ", neigborsize: " << T << endl;
	cout << "function: " << function << ", maxalpha: " << ALPHA_MAX << ", minalpha: " << ALPHA_MIN << endl;
	cout << "penalty1: " << PENALTY << ", penaly2: " << PENALTY2 << endl;
	cout << "shita: " << shita << endl;
	cout << "valnum: " << VALNUM << ", generation: " << GEN << endl << endl;
	cout << "sbx_eta: " << SBX_ETA << ", mut_eta: " << MUT_ETA << endl;
	cout << "crossover probabilty: " << CO << ", mutation probability: " << mut << " / " << item << endl;
	cout << "k_fac: " << k_fac << ", l_fac: " << l_fac << endl;
	cout << "k_dtlz: " << k_dtlz << endl;
	cout << "improve_lambda: " << improve_lambda << endl;
	//cout << "switch_pbi: " << switch_pbi << endl;
	//cout << "switch_rate: " << switch_rate << endl;
	double MUT = mut / item;

	char datafile[50];
	char outputfile[30] = "graph";
	char ccc[5] = "00";
	sprintf(ccc, "%d", atoi(argv[2]));
	strcat(outputfile, ccc);
	strcpy(datafile, outputfile);
	strcat(outputfile, ".txt");
	
	char outputfile_second[30] = "graph_second";
	strcat(outputfile_second, ccc);
	strcat(outputfile_second, ".txt");

	//progress_output
	//_mkdir(datafile);


	int N1 = POPSIZE1;
	int H1 = MOEAD_H1;

	int N2 = POPSIZE2;
	int H2 = MOEAD_H2;

	int N = N1 + N2;


	//for (int n = 0; n < N; n++){
	//	lambda[n] = new double[ob];
	//	for (int o = 0; o < ob; o++){
	//		fl >> lambda[n][o];
	//		/*		if (lambda[n][o] == 0){
	//		lambda[n][o] = 0.0000001;
	//		}*/
	//	}
	//}
	//fl.close();

	/*N = (int)pow(2, ob) - 1 + ob;
	T = (int) N / 10;
	if (T == 0){
		T = 1;
	}*/
	double **lambda_h1 = new double*[N1];
	for (int i = 0; i < N1; i++) lambda_h1[i] = new double[ob];

	//simplex-centroid design
	/*double *first_lambda = new double[ob];
	for (int o = 0; o < ob; o++){
		first_lambda[o] = 0.0;
	}

	int nnnn = 0;
	for (int o = 1; o <= ob; o++){
		for (int oo = 0; oo < ob; oo++){
			first_lambda[oo] = 0.0;
		}
		for (int oo = 0; oo < o; oo++){
			first_lambda[oo] = 1.0 / o;
		}
		do{
			for (int oo = 0; oo < ob; oo++){
				lambda[nnnn][oo] = first_lambda[oo];
			}
			nnnn++;
		} while (prev_permutation(first_lambda, first_lambda + ob));
	}

	for (int o = 0; o < ob; o++){
		for (int oo = 0; oo < ob; oo++){
			lambda[nnnn][oo] = 1.0 / (2.0 * ob);
		}
		lambda[nnnn][o] = (ob - 1) / (2.0 * ob);
		nnnn++;
	}*/

	/////////////////////////////////////////////////////////////


	//simplex-lattice design
	int *t_lambda = new int[ob];
	for (int i = 0; i < ob - 1; i++) t_lambda[i] = 0;

	
	bool end_flag = false;
	for (int i = 0; i < N1; i++){
		int temp = 0;
		for (int j = 0; j < ob - 1; j++) temp += t_lambda[j];
		t_lambda[ob - 1] = H1 - temp;

		bool assign_flag = true;
		for (int j = 0; j < ob; j++){
			assign_flag &= (t_lambda[ob - 1] >= 0);
		}

		if (assign_flag){
			for (int j = 0; j < ob; j++){
				lambda_h1[i][j] = (double)(t_lambda[j] / (double)H1);
				//cout << t_lambda[j] << " ";
			}
			//cout << endl;
		}

		for (int j = ob - 1; j >= 1; j--){
			if (t_lambda[j] <= 0){ //continue affirmation
				if (j == 1) end_flag = true;
				int sum = 0;
				for (int k = 0; k < j; k++){
					sum += t_lambda[k];
					if (sum > H1){
						t_lambda[k] = 0;
					}
				}
			}
			else{
				t_lambda[j - 1] += 1;
				for (int k = j; k < ob - 1; k++){
					t_lambda[k] = 0;
				}
				break;
			}
		}
		if (end_flag) break;
	}




	double **lambda_h2 = new double*[N2];
	for (int i = 0; i < N2; i++) lambda_h2[i] = new double[ob];

	int *tt_lambda = new int[ob];
	for (int i = 0; i < ob - 1; i++) tt_lambda[i] = 0;


	bool end_flag2 = false;
	for (int i = 0; i < N2; i++){
		int temp = 0;
		for (int j = 0; j < ob - 1; j++) temp += tt_lambda[j];
		tt_lambda[ob - 1] = H2 - temp;

		bool assign_flag = true;
		for (int j = 0; j < ob; j++){
			assign_flag &= (tt_lambda[ob - 1] >= 0);
		}

		if (assign_flag){
			for (int j = 0; j < ob; j++){
				lambda_h2[i][j] = (double)(tt_lambda[j] / (double)H2);
				//cout << tt_lambda[j] << " ";
			}
			//cout << endl;
		}

		for (int j = ob - 1; j >= 1; j--){
			if (tt_lambda[j] <= 0){ //continue affirmation
				if (j == 1) end_flag2 = true;
				int sum = 0;
				for (int k = 0; k < j; k++){
					sum += tt_lambda[k];
					if (sum > H2){
						tt_lambda[k] = 0;
					}
				}
			}
			else{
				tt_lambda[j - 1] += 1;
				for (int k = j; k < ob - 1; k++){
					tt_lambda[k] = 0;
				}
				break;
			}
		}
		if (end_flag2) break;
	}

	for (int i = 0; i < N2; i++){
		for (int o = 0; o < ob; o++){
			lambda_h2[i][o] = ((lambda_h2[i][o] + (1.0 / (double)ob)) / 2.0);
		}
	}

	//if you use two layer this is deleted from here.
	double **lambda = new double *[N];
	for (int i = 0; i < N; i++) lambda[i] = new double[ob];

	int n12 = 0;
	for (n12 = 0; n12 < N1; n12++){
		for (int o = 0; o < ob; o++){
			lambda[n12][o] = lambda_h1[n12][o];
		}
	}
	for (; n12 < N; n12++){
		for (int o = 0; o < ob; o++){
			lambda[n12][o] = lambda_h2[n12 - N1][o];
		}
	}


	double **lambda2 = new double*[N];
	for (int i = 0; i < N; i++) lambda2[i] = new double[ob];


	/*for (int n = 0; n < N; n++){
		for (int o = 0; o < ob; o++){
			lambda[n][o] = 1.0 - lambda[n][o];
		}
	}*/

	for (int n = 0; n < N; n++){
		for (int o = 0; o < ob; o++){
			lambda2[n][o] = 1.0 - lambda[n][o];
		}
		double sum_lam = 0.0;
		for (int o = 0; o < ob; o++){
			sum_lam += lambda2[n][o];
		}
		for (int o = 0; o < ob; o++){
			lambda2[n][o] = lambda2[n][o] / sum_lam;
		}
	}
	//to here.


	int **B = new int*[N];
	for (int i = 0; i < N; i++){
		B[i] = new int[T];
	}
	int **B_second = new int*[N];
	for (int i = 0; i < N; i++){
		B_second[i] = new int[T];
	}

	double *zmin = new double[ob];
	double *zmax = new double[ob];

	double *fmin = new double[ob];
	double *fmax = new double[ob];


	double **x = new double*[N];
	for (int n = 0; n < N; n++){
		x[n] = new double[item];
	}

	int **x_com = new int *[N];
	for (int n = 0; n < N; n++){
		x_com[n] = new int[item];
	}

	double **FV = new double*[N];
	for (int n = 0; n < N; n++){
		FV[n] = new double[ob];
	}
	double *y = new double[item];
	int *y_com = new int[item];
	double *y_fit = new double[ob];


	//archive
	double **x_second = new double*[N];
	for (int n = 0; n < N; n++){
		x_second[n] = new double[item];
	}

	int **x_com_second = new int *[N];
	for (int n = 0; n < N; n++){
		x_com_second[n] = new int[item];
	}

	double **FV_second = new double*[N];
	for (int n = 0; n < N; n++){
		FV_second[n] = new double[ob];
	}
	double *y_second = new double[item];
	int *y_com_second = new int[item];
	double *y_fit_second = new double[ob];
	//archive



	bool vector_chase_flag = false;
	if (atoi(argv[2]) == 0){
		vector_chase_flag = true;
	}

	double **non_dom_set = new double *[N];
	for (int i = 0; i < N; i++){
		non_dom_set[i] = new double[ob];
	}


	ofstream fout_ref_max;
	ofstream fout_ref_min;

	

	if (vector_chase_flag == true){
		fout_ref_max.open("ref_max0.txt");
		fout_ref_min.open("ref_min0.txt");

	}

	if (atoi(argv[2]) != 0){
		for (int i = 0; i < N; i++){
			int j = genrand_int31() % N;
			for (int o = 0; o < ob; o++){
				double t = lambda[i][o];
				lambda[i][o] = lambda[j][o];
				lambda[j][o] = t;
				t = lambda2[i][o];
				lambda2[i][o] = lambda2[j][o];
				lambda2[j][o] = t;
			}
		}
	/*	ofstream vec("vec1.txt");
		for (int i = 0; i < N; i++){
			for (int o = 0; o < ob; o++){
				vec << lambda[i][o] << "\t";
			}
			vec << endl;
		}*/
	}
	else{
		if (improve_lambda == false){
			ofstream vec("vec0.txt");
			for (int i = 0; i < N; i++){
				for (int o = 0; o < ob; o++){
					vec << lambda[i][o] << "\t";
				}
				vec << endl;
			}
		}
		else{
			ofstream vec("vec0.txt");
			for (int i = 0; i < N; i++){
				for (int o = 0; o < ob; o++){
					vec << lambda2[i][o] << "\t";
				}
				vec << endl;
			}
		}
	}

	cout << N << " " << T << endl;
	for (int o = 0; o < ob; o++){
		cout << lambda[N - 1][o] << " ";
	}
	cout << endl;
	///************************************Step1 Initialization***********************************/
	//Step1.2

	t_neighborhood(B, lambda, N, T, ob);
	t_neighborhood(B_second, lambda2, N, T, ob);

	//Step1.3

	if (combination == false){
		for (int n = 0; n < N; n++){
			for (int i = 0; i < item; i++){
				x[n][i] = nextDoubleIE();
				x_second[n][i] = nextDoubleIE();
			}
		}
	}
	else{

		if (strcmp("test20160301", problem) == 0 || strcmp("test20160303", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					x_com[n][i] = genrand_int31() % ob + 1;
					x_com_second[n][i] = genrand_int31() % ob + 1;
				}
			}

		}
		else if (strcmp("test20160303", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					x_com[n][i] = genrand_int31() % k_fac + 1;
					x_com_second[n][i] = genrand_int31() % k_fac + 1;
				}
			}

		}
		else{
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					x_com[n][i] = genrand_int31() % 2;
					x_com_second[n][i] = genrand_int31() % 2;
				}
			}
		}
	}


	if (combination == false){
		for (int n = 0; n < N; n++){
			fitness(x[n], FV[n], problem, ob, k_val, l_val, k_dtlz, P, point_num);
			fitness(x_second[n], FV_second[n], problem, ob, k_val, l_val, k_dtlz, P, point_num);
			/*for (int o = 0; o < ob; o++){
				FV[n][o] = -1.0 * FV[n][o];
			}*/
		}
	}
	else{
		if (strcmp(problem, "knapsack") == 0){
			for (int n = 0; n < N; n++){
				fitness_knap(x_com[n], FV[n], item, ob);
				fitness_knap(x_com_second[n], FV_second[n], item, ob);
			}
		}
		else if (strcmp("test20160301", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					test20160301(x_com[n], FV[n], ob, item);
					test20160301(x_com_second[n], FV_second[n], ob, item);
				}
			}
		}
		else if (strcmp("test20160303", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					test20160303(x_com[n], FV[n], ob, item, k_fac);
					test20160303(x_com_second[n], FV_second[n], ob, item, k_fac);
				}
			}
		}
		else{
			exit(1);
		}
	}

	
	ofstream fffout("dtlz-ini.dat");

	for (int i = 0; i < N; i++){
	for (int o = 0; o < ob; o++){
	fffout << FV[i][o] << "\t";
	}
	fffout << endl;

	}

	fffout.close();
	
	//Step1.4
	for (int o = 0; o < ob; o++){
		zmin[o] = 100000.0;
		fmin[o] = 1000000.0;
	}
	for (int o = 0; o < ob; o++){
		zmax[o] = 0.0;
		fmax[o] = 0.0;
	}
	for (int n = 0; n < N; n++){
		max_z(FV[n], zmax, ob, ALPHA_MAX, fmax);
	}

	for (int n = 0; n < N; n++){
		min_z(FV[n], zmin, ob, ALPHA_MIN, fmin);
	}
	int plus_minus = -1;

	if (ref_non_dom_flag == true){
		int ndsize = 0;
		int SIGN = 1;
		if (strcmp(max_min, "min") == 0){
			SIGN = -1;
		}
		non_dominated_set(FV, non_dom_set, ndsize, ob, N, SIGN);
		for (int n = 0; n < ndsize; n++){
			max_z(non_dom_set[n], zmax, ob, ALPHA_MAX, fmax);
		}

		for (int n = 0; n < ndsize; n++){
			min_z(non_dom_set[n], zmin, ob, ALPHA_MIN, fmin);
		}
	}


	/************************************Step2 Update***********************************/

	//progress_output
	char detaildata[100] = "./";
	strcat(detaildata, datafile);
	strcat(detaildata, "/");
	strcat(detaildata, "eva");
	char gennumch[8] = "00";
	sprintf(gennumch, "%d", 0);
	strcat(detaildata, gennumch);
	strcat(detaildata, ".txt");

	double alpha = ALPHA_MAX;
	double alpha_max = ALPHA_MAX;
	double alpha_min = ALPHA_MIN;
	int divide = VALNUM / N;
	int finiflag = 0;
	for (int g = 0; g < GEN; g++){
		for (int i = 0; i < N; i++){
			if (g * N + i < VALNUM){
			
				//Step2.1 Reproduction
				int k = genrand_int31() % T;
				int l = genrand_int31() % T;
				
				if (combination == false){
					make_new_pop(x[B[i][k]], x[B[i][l]], y, item, MUT, CO, MUT_ETA, SBX_ETA);
				}
				else{
					if (strcmp(problem, "knapsack") == 0){
						make_new_pop(x_com[B[i][k]], x_com[B[i][l]], y_com, item, MUT, CO);
					}
					else if (strcmp(problem, "test20160301") == 0){
						make_new_pop_20160301(x_com[B[i][k]], x_com[B[i][l]], y_com, item, MUT, CO, ob);
					}
					else if (strcmp(problem, "test20160303") == 0){
						make_new_pop_20160301(x_com[B[i][k]], x_com[B[i][l]], y_com, item, MUT, CO, k_fac);
					}
					else{
						exit(1);
					}
				}


				int k_sec = genrand_int31() % T;
				int l_sec = genrand_int31() % T;

				if (combination == false){
					make_new_pop(x_second[B[i][k_sec]], x_second[B_second[i][l_sec]], y_second, item, MUT, CO, MUT_ETA, SBX_ETA);
				}
				else{
					if (strcmp(problem, "knapsack") == 0){
						make_new_pop(x_com_second[B[i][k_sec]], x_com_second[B[i][l_sec]], y_com_second, item, MUT, CO);
					}
					else if (strcmp(problem, "test20160301") == 0){
						make_new_pop_20160301(x_com_second[B[i][k_sec]], x_com_second[B[i][l_sec]], y_com_second, item, MUT, CO, ob);
					}
					else if (strcmp(problem, "test20160303") == 0){
						make_new_pop_20160301(x_com_second[B[i][k_sec]], x_com_second[B[i][l_sec]], y_com_second, item, MUT, CO, k_fac);
					}
					else{
						exit(1);
					}
				}


				//Step2.2 Improvement
				if (combination == false){
					fitness(y, y_fit, problem, ob, k_val, l_val, k_dtlz, P, point_num);
					fitness(y_second, y_fit_second, problem, ob, k_val, l_val, k_dtlz, P, point_num);

				
				}
				else{
					if (strcmp(problem, "knapsack") == 0){
						fitness_knap(y_com, y_fit, item, ob);
						fitness_knap(y_com_second, y_fit_second, item, ob);
					}
					else if (strcmp(problem, "test20160301") == 0){
						test20160301(y_com, y_fit, ob, item);
						test20160301(y_com_second, y_fit_second, ob, item);
					}
					else if (strcmp(problem, "test20160303") == 0){
						test20160303(y_com, y_fit, ob, item, k_fac);
						test20160303(y_com_second, y_fit_second, ob, item, k_fac);
					}
					else{
						exit(1);
					}
					
				}

				//Step2.3 Update of z
				if (strcmp(max_min, "min") == 0){
					min_z(y_fit, zmin, ob, ALPHA_MIN, fmin);
					min_z(y_fit_second, zmin, ob, ALPHA_MIN, fmin);

					for (int o = 0; o < ob; o++){
						zmax[o] = 0.0;
						fmax[o] = 0.0;
					}
					for (int n = 0; n < N; n++){
						max_z(FV[n], zmax, ob, ALPHA_MAX, fmax);
						max_z(FV_second[n], zmax, ob, ALPHA_MAX, fmax);
					}

				//	max_z(y_fit, zmax, ob, ALPHA_MAX, fmax);

				}
				else{
					max_z(y_fit, zmax, ob, ALPHA_MAX, fmax);
					max_z(y_fit_second, zmax, ob, ALPHA_MAX, fmax);

					for (int o = 0; o < ob; o++){
						zmin[o] = 100000.0;
						fmin[o] = 100000.0;
					}

					for (int n = 0; n < N; n++){
						min_z(FV[n], zmin, ob, ALPHA_MIN, fmin);
						min_z(FV_second[n], zmin, ob, ALPHA_MIN, fmin);

					}
			
				}

				for (int j = 0; j < T; j++){
					int jejeje = B[i][j];
					double shita_temp = shita;
					if (strcmp(function, "fix_double_pbi") == 0 || strcmp(function, "fix_double_inverted_pbi") == 0
						|| strcmp(function, "quadratic_pbi") == 0 || strcmp(function, "quadratic_inverted_pbi") == 0){
						double diagonal = 0.0;
						for (int o = 0; o < ob; o++){
							diagonal += (fmax[o] - fmin[o]);
						}
						shita = /*(1.0 + 100.0 * (1.0 - ((double)g + 1.0) / ((double)(VALNUM / N) + 1.0))) */ shita_temp * diagonal / (double)ob / (double)H1;
						
					}
		
					if (function_valuation(problem, function,
						FV[jejeje], y_fit, 
						improve_lambda,lambda[jejeje],lambda2[jejeje], 
						zmax, zmin, 
						ob, 
						PENALTY, PENALTY2, shita, max_min)){
						
						if (combination == false){
							for (int it = 0; it < item; it++){
								x[jejeje][it] = y[it];
							}
						}
						else{
							for (int it = 0; it < item; it++){
								x_com[jejeje][it] = y_com[it];
							}
						}
						for (int o = 0; o < ob; o++){
							FV[jejeje][o] = y_fit[o];
						}
						
					}
					shita = shita_temp;
				}

				for (int j = 0; j < T; j++){
					int jejeje = B_second[i][j];
					double shita_temp = shita;
					if (strcmp(function, "fix_double_pbi") == 0 || strcmp(function, "fix_double_inverted_pbi") == 0
						|| strcmp(function, "quadratic_pbi") == 0 || strcmp(function, "quadratic_inverted_pbi") == 0){
						double diagonal = 0.0;
						for (int o = 0; o < ob; o++){
							diagonal += (fmax[o] - fmin[o]);
						}
						shita = /*(1.0 + 100.0 * (1.0 - ((double)g + 1.0) / ((double)(VALNUM / N) + 1.0))) */ shita_temp * diagonal / (double)ob / (double)H1;

					}

					if (function_valuation(problem, function,
						FV_second[jejeje], y_fit_second,
						improve_lambda, lambda2[jejeje], lambda[jejeje],
						zmax, zmin,
						ob,
						PENALTY, PENALTY2, shita, max_min)){

						if (combination == false){
							for (int it = 0; it < item; it++){
								x_second[jejeje][it] = y_second[it];
							}
						}
						else{
							for (int it = 0; it < item; it++){
								x_com_second[jejeje][it] = y_com_second[it];
							}
						}
						for (int o = 0; o < ob; o++){
							FV_second[jejeje][o] = y_fit_second[o];
						}

					}
					shita = shita_temp;
				}

				if (((i % N) == 0) && (g != 0)){

					double **sec_temp = new double *[N];
					double **fis_temp = new double *[N];
					for (int je = 0; je < N; je++){
						sec_temp[je] = new double[item];
						fis_temp[je] = new double[item];

					}
					int **sec_com_temp = new int *[N];
					int **fis_com_temp = new int *[N];
					for (int je = 0; je < N; je++){
						sec_com_temp[je] = new int[item];
						fis_com_temp[je] = new int[item];

					}
					double **sec_fv_temp = new double *[N];
					double **fis_fv_temp = new double *[N];
					for (int je = 0; je < N; je++){
						sec_fv_temp[je] = new double[item];
						fis_fv_temp[je] = new double[item];

					}

					for (int je = 0; je < N; je++){
						for (int it = 0; it < item; it++){
							sec_temp[je][it] = x_second[je][it];
							fis_temp[je][it] = x[je][it];

							sec_com_temp[je][it] = x_com_second[je][it];
							fis_com_temp[je][it] = x_com[je][it];

						}
						for (int o = 0; o < ob; o++){
							sec_fv_temp[je][o] = FV_second[je][o];
							fis_fv_temp[je][o] = FV[je][o];

						}
					}

					for (int jejeje = 0; jejeje < N; jejeje++){
						for (int je = 0; je < N; je++){
							if (function_valuation(problem, function,
								FV_second[jejeje], FV[je],
								improve_lambda, lambda2[jejeje], lambda[jejeje],
								zmax, zmin,
								ob,
								PENALTY, PENALTY2, shita, max_min)){

								if (combination == false){
									for (int it = 0; it < item; it++){
										x_second[jejeje][it] = x[je][it];
									}
								}
								else{
									for (int it = 0; it < item; it++){
										x_com_second[jejeje][it] = x_com[je][it];
									}
								}
								for (int o = 0; o < ob; o++){
									FV_second[jejeje][o] = FV[je][o];
								}

							}
						}

					}
					for (int jejeje = 0; jejeje < N; jejeje++){
						for (int je = 0; je < N; je++){
							if (function_valuation(problem, function,
								FV[jejeje], sec_fv_temp[je],
								improve_lambda, lambda[jejeje], lambda2[jejeje],
								zmax, zmin,
								ob,
								PENALTY, PENALTY2, shita, max_min)){

								if (combination == false){
									for (int it = 0; it < item; it++){
										x[jejeje][it] = sec_temp[je][it];
									}
								}
								else{
									for (int it = 0; it < item; it++){
										x_com[jejeje][it] = sec_com_temp[je][it];
									}
								}
								for (int o = 0; o < ob; o++){
									FV[jejeje][o] = sec_fv_temp[je][o];
								}

							}
						}

					}


					for (int je = 0; je < N; je++){
						delete[] sec_temp[je];
						delete[] fis_temp[je];
					}
					delete[] sec_temp;
					delete[] fis_temp;

					for (int je = 0; je < N; je++){
						delete[] sec_com_temp[je];
						delete[] fis_com_temp[je];
					}
					delete[] sec_com_temp;
					delete[] fis_com_temp;
					for (int je = 0; je < N; je++){
						delete[] sec_fv_temp[je];
						delete[] fis_fv_temp[je];
					}
					delete[] sec_fv_temp;
					delete[] fis_fv_temp;
				}

			


				if (vector_chase_flag == true){
					if (!((g * N + i) % 105)){
						for (int o = 0; o < ob - 1; o++){
							fout_ref_max << zmax[o] << "\t";
						}
						fout_ref_max << zmax[ob - 1];
						fout_ref_max << endl;
						for (int o = 0; o < ob - 1; o++){
							fout_ref_min << zmin[o] << "\t";
						}
						fout_ref_min << zmin[ob - 1];
						fout_ref_min << endl;
					}
				}
				if (!((g * N + i) % 100000)){
					cout << g * N + i << endl;
					
				}

				//progress_output
				if (!((g * N + i) % 100)){
					strcpy(detaildata, "./");
					strcat(detaildata, datafile);
					strcat(detaildata, "/");
					strcat(detaildata, "eva");

					sprintf(gennumch, "%d", g * N + i);
					strcat(detaildata, gennumch);
					strcat(detaildata, ".txt");


				}
			}
			else{

				finiflag = 1;
				break;
			}

		}

		if (finiflag == 1){
			break;
		}
	}





	
	if (vector_chase_flag = true){
		for (int o = 0; o < ob - 1; o++){
			fout_ref_max << zmax[o] << "\t";
		}
		fout_ref_max << zmax[ob - 1];
		fout_ref_max << endl;
		for (int o = 0; o < ob - 1; o++){
			fout_ref_min << zmin[o] << "\t";
		}
		fout_ref_min << zmin[ob - 1];
		fout_ref_min << endl;
		fout_ref_max.close();
		fout_ref_min.close();
	}

	ofstream fout(outputfile);
	fout << setprecision(20);
	for (int i = 0; i < N; i++){
		for (int o = 0; o < ob - 1; o++){
			fout << FV[i][o] << "\t";
		}
		fout << FV[i][ob - 1];
		fout << endl;
	}

	fout.close();

	ofstream fout_second(outputfile_second);
	fout_second << setprecision(20);
	for (int i = 0; i < N; i++){
		for (int o = 0; o < ob - 1; o++){
			fout_second << FV_second[i][o] << "\t";
		}
		fout_second << FV_second[i][ob - 1];
		fout_second << endl;
	}

	fout_second.close();











	return 0;
}

