#include <array>
#include <iostream>
#include <vector>
#include <math.h>
#include "model.h"

Models::Models() {
	 models_ublas_vec = {&Models::model_0, &Models::model_1, &Models::model_2, &Models::model_3, &Models::model_4, &Models::model_5, &Models::model_6};
	 models_jac_vec = {&Models::jac_0, &Models::jac_1, &Models::jac_2, &Models::jac_3, &Models::jac_4, &Models::jac_5, &Models::jac_6};
};

void Models::run_model_ublas(const ublas_vec_t &y , ublas_vec_t &dxdt , double t, std::vector <double> &part_params, int &model_ref)
{
	(this->*models_ublas_vec[model_ref])(y, dxdt, t, part_params);
}

void Models::run_jac(const ublas_vec_t & x , ublas_mat_t &J , const double & t , ublas_vec_t &dfdt, std::vector <double> &part_params, int &model_ref)
{
	(this->*models_jac_vec[model_ref])(x, J, t, dfdt, part_params);
}

void Models::model_0(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	//Species order is: N_x N_c S_glu B_2 B_1 A_2 A_1 
	dydt[0] = -std::pow(y[4], n_omega_1)*y[0]*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_2)*y[1]*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = std::pow(y[5], nB_2)*y[1]*kBmax_2/(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2)) - y[3]*D;
	dydt[4] = std::pow(y[6], nB_1)*y[0]*kBmax_1/(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1)) - y[4]*D;
	dydt[5] = -y[5]*D + y[1]*kA_2;
	dydt[6] = -y[6]*D + y[0]*kA_1;

}
void Models::jac_0(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	J( 0 , 0 ) = -std::pow(y[4], n_omega_1)*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = std::pow(y[4], 2*n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*std::pow(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[4], n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 0 , 5 ) = 0;
	J( 0 , 6 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_2)*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*std::pow(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2), 2)) - std::pow(y[3], n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)));
	J( 1 , 4 ) = 0;
	J( 1 , 5 ) = 0;
	J( 1 , 6 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 2 , 5 ) = 0;
	J( 2 , 6 ) = 0;
	J( 3 , 0 ) = 0;
	J( 3 , 1 ) = std::pow(y[5], nB_2)*kBmax_2/(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2));
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = 0;
	J( 3 , 5 ) = -std::pow(y[5], 2*nB_2)*y[1]*kBmax_2*nB_2/(y[5]*std::pow(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2), 2)) + std::pow(y[5], nB_2)*y[1]*kBmax_2*nB_2/(y[5]*(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2)));
	J( 3 , 6 ) = 0;
	J( 4 , 0 ) = std::pow(y[6], nB_1)*kBmax_1/(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1));
	J( 4 , 1 ) = 0;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;
	J( 4 , 5 ) = 0;
	J( 4 , 6 ) = -std::pow(y[6], 2*nB_1)*y[0]*kBmax_1*nB_1/(y[6]*std::pow(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1), 2)) + std::pow(y[6], nB_1)*y[0]*kBmax_1*nB_1/(y[6]*(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1)));
	J( 5 , 0 ) = 0;
	J( 5 , 1 ) = kA_2;
	J( 5 , 2 ) = 0;
	J( 5 , 3 ) = 0;
	J( 5 , 4 ) = 0;
	J( 5 , 5 ) = -D;
	J( 5 , 6 ) = 0;
	J( 6 , 0 ) = kA_1;
	J( 6 , 1 ) = 0;
	J( 6 , 2 ) = 0;
	J( 6 , 3 ) = 0;
	J( 6 , 4 ) = 0;
	J( 6 , 5 ) = 0;
	J( 6 , 6 ) = -D;

}
void Models::model_1(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	//Species order is: N_x N_c S_glu B_2 B_1 A_2 A_1 
	dydt[0] = -std::pow(y[4], n_omega_1)*y[0]*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_2)*y[1]*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = std::pow(y[6], nB_2)*y[1]*kBmax_2/(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2)) - y[3]*D;
	dydt[4] = -y[4]*D + std::pow(KB_1, nB_1)*y[0]*kBmax_1/(std::pow(y[5], nB_1) + std::pow(KB_1, nB_1));
	dydt[5] = -y[5]*D + y[1]*kA_2;
	dydt[6] = -y[6]*D + y[0]*kA_1;

}
void Models::jac_1(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	J( 0 , 0 ) = -std::pow(y[4], n_omega_1)*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = std::pow(y[4], 2*n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*std::pow(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[4], n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 0 , 5 ) = 0;
	J( 0 , 6 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_2)*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*std::pow(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2), 2)) - std::pow(y[3], n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)));
	J( 1 , 4 ) = 0;
	J( 1 , 5 ) = 0;
	J( 1 , 6 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 2 , 5 ) = 0;
	J( 2 , 6 ) = 0;
	J( 3 , 0 ) = 0;
	J( 3 , 1 ) = std::pow(y[6], nB_2)*kBmax_2/(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2));
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = 0;
	J( 3 , 5 ) = 0;
	J( 3 , 6 ) = -std::pow(y[6], 2*nB_2)*y[1]*kBmax_2*nB_2/(y[6]*std::pow(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2), 2)) + std::pow(y[6], nB_2)*y[1]*kBmax_2*nB_2/(y[6]*(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2)));
	J( 4 , 0 ) = std::pow(KB_1, nB_1)*kBmax_1/(std::pow(y[5], nB_1) + std::pow(KB_1, nB_1));
	J( 4 , 1 ) = 0;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;
	J( 4 , 5 ) = -std::pow(y[5], nB_1)*std::pow(KB_1, nB_1)*y[0]*kBmax_1*nB_1/(y[5]*std::pow(std::pow(y[5], nB_1) + std::pow(KB_1, nB_1), 2));
	J( 4 , 6 ) = 0;
	J( 5 , 0 ) = 0;
	J( 5 , 1 ) = kA_2;
	J( 5 , 2 ) = 0;
	J( 5 , 3 ) = 0;
	J( 5 , 4 ) = 0;
	J( 5 , 5 ) = -D;
	J( 5 , 6 ) = 0;
	J( 6 , 0 ) = kA_1;
	J( 6 , 1 ) = 0;
	J( 6 , 2 ) = 0;
	J( 6 , 3 ) = 0;
	J( 6 , 4 ) = 0;
	J( 6 , 5 ) = 0;
	J( 6 , 6 ) = -D;

}
void Models::model_2(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double K_c = part_params[2];
	const double K_omega_1 = part_params[3];
	const double K_x = part_params[4];
	const double S0_glu = part_params[5];
	const double g_c = part_params[6];
	const double g_x = part_params[7];
	const double kA_1 = part_params[8];
	const double kBmax_1 = part_params[9];
	const double mu_max_c = part_params[10];
	const double mu_max_x = part_params[11];
	const double nB_1 = part_params[12];
	const double n_omega_1 = part_params[13];
	const double omega_max_1 = part_params[14];

	//Species order is: N_x N_c S_glu B_1 A_1 
	dydt[0] = -D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_1)*y[1]*omega_max_1/(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = -y[3]*D + std::pow(KB_1, nB_1)*y[0]*kBmax_1/(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1));
	dydt[4] = -y[4]*D + y[0]*kA_1;

}
void Models::jac_2(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double K_c = part_params[2];
	const double K_omega_1 = part_params[3];
	const double K_x = part_params[4];
	const double S0_glu = part_params[5];
	const double g_c = part_params[6];
	const double g_x = part_params[7];
	const double kA_1 = part_params[8];
	const double kBmax_1 = part_params[9];
	const double mu_max_c = part_params[10];
	const double mu_max_x = part_params[11];
	const double nB_1 = part_params[12];
	const double n_omega_1 = part_params[13];
	const double omega_max_1 = part_params[14];

	J( 0 , 0 ) = -D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_1)*omega_max_1/(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_1)*y[1]*n_omega_1*omega_max_1/(y[3]*std::pow(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[3], n_omega_1)*y[1]*n_omega_1*omega_max_1/(y[3]*(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 1 , 4 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 3 , 0 ) = std::pow(KB_1, nB_1)*kBmax_1/(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1));
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = -std::pow(y[4], nB_1)*std::pow(KB_1, nB_1)*y[0]*kBmax_1*nB_1/(y[4]*std::pow(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1), 2));
	J( 4 , 0 ) = kA_1;
	J( 4 , 1 ) = 0;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;

}
void Models::model_3(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double K_c = part_params[2];
	const double K_omega_1 = part_params[3];
	const double K_x = part_params[4];
	const double S0_glu = part_params[5];
	const double g_c = part_params[6];
	const double g_x = part_params[7];
	const double kA_1 = part_params[8];
	const double kBmax_1 = part_params[9];
	const double mu_max_c = part_params[10];
	const double mu_max_x = part_params[11];
	const double nB_1 = part_params[12];
	const double n_omega_1 = part_params[13];
	const double omega_max_1 = part_params[14];

	//Species order is: N_x N_c S_glu B_1 A_1 
	dydt[0] = -D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_1)*y[1]*omega_max_1/(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = std::pow(y[4], nB_1)*y[0]*kBmax_1/(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1)) - y[3]*D;
	dydt[4] = -y[4]*D + y[1]*kA_1;

}
void Models::jac_3(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double K_c = part_params[2];
	const double K_omega_1 = part_params[3];
	const double K_x = part_params[4];
	const double S0_glu = part_params[5];
	const double g_c = part_params[6];
	const double g_x = part_params[7];
	const double kA_1 = part_params[8];
	const double kBmax_1 = part_params[9];
	const double mu_max_c = part_params[10];
	const double mu_max_x = part_params[11];
	const double nB_1 = part_params[12];
	const double n_omega_1 = part_params[13];
	const double omega_max_1 = part_params[14];

	J( 0 , 0 ) = -D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_1)*omega_max_1/(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_1)*y[1]*n_omega_1*omega_max_1/(y[3]*std::pow(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[3], n_omega_1)*y[1]*n_omega_1*omega_max_1/(y[3]*(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 1 , 4 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 3 , 0 ) = std::pow(y[4], nB_1)*kBmax_1/(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1));
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = -std::pow(y[4], 2*nB_1)*y[0]*kBmax_1*nB_1/(y[4]*std::pow(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1), 2)) + std::pow(y[4], nB_1)*y[0]*kBmax_1*nB_1/(y[4]*(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1)));
	J( 4 , 0 ) = 0;
	J( 4 , 1 ) = kA_1;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;

}
void Models::model_4(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double K_c = part_params[2];
	const double K_omega_1 = part_params[3];
	const double K_x = part_params[4];
	const double S0_glu = part_params[5];
	const double g_c = part_params[6];
	const double g_x = part_params[7];
	const double kA_1 = part_params[8];
	const double kA_2 = part_params[9];
	const double kBmax_1 = part_params[10];
	const double mu_max_c = part_params[11];
	const double mu_max_x = part_params[12];
	const double nB_1 = part_params[13];
	const double n_omega_1 = part_params[14];
	const double omega_max_1 = part_params[15];

	//Species order is: N_x N_c S_glu B_1 A_2 A_1 
	dydt[0] = -D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_1)*y[1]*omega_max_1/(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = std::pow(y[4], nB_1)*std::pow(KB_1, nB_1)*y[0]*kBmax_1/((std::pow(y[5], nB_1) + std::pow(KB_1, nB_1))*(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1))) - y[3]*D;
	dydt[4] = -y[4]*D + y[1]*kA_2;
	dydt[5] = -y[5]*D + y[0]*kA_1;

}
void Models::jac_4(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double K_c = part_params[2];
	const double K_omega_1 = part_params[3];
	const double K_x = part_params[4];
	const double S0_glu = part_params[5];
	const double g_c = part_params[6];
	const double g_x = part_params[7];
	const double kA_1 = part_params[8];
	const double kA_2 = part_params[9];
	const double kBmax_1 = part_params[10];
	const double mu_max_c = part_params[11];
	const double mu_max_x = part_params[12];
	const double nB_1 = part_params[13];
	const double n_omega_1 = part_params[14];
	const double omega_max_1 = part_params[15];

	J( 0 , 0 ) = -D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = 0;
	J( 0 , 5 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_1)*omega_max_1/(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_1)*y[1]*n_omega_1*omega_max_1/(y[3]*std::pow(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[3], n_omega_1)*y[1]*n_omega_1*omega_max_1/(y[3]*(std::pow(y[3], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 1 , 4 ) = 0;
	J( 1 , 5 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 2 , 5 ) = 0;
	J( 3 , 0 ) = std::pow(y[4], nB_1)*std::pow(KB_1, nB_1)*kBmax_1/((std::pow(y[5], nB_1) + std::pow(KB_1, nB_1))*(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1)));
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = -std::pow(y[4], 2*nB_1)*std::pow(KB_1, nB_1)*y[0]*kBmax_1*nB_1/(y[4]*(std::pow(y[5], nB_1) + std::pow(KB_1, nB_1))*std::pow(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1), 2)) + std::pow(y[4], nB_1)*std::pow(KB_1, nB_1)*y[0]*kBmax_1*nB_1/(y[4]*(std::pow(y[5], nB_1) + std::pow(KB_1, nB_1))*(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1)));
	J( 3 , 5 ) = -std::pow(y[5], nB_1)*std::pow(y[4], nB_1)*std::pow(KB_1, nB_1)*y[0]*kBmax_1*nB_1/(y[5]*std::pow(std::pow(y[5], nB_1) + std::pow(KB_1, nB_1), 2)*(std::pow(y[4], nB_1) + std::pow(KB_1, nB_1)));
	J( 4 , 0 ) = 0;
	J( 4 , 1 ) = kA_2;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;
	J( 4 , 5 ) = 0;
	J( 5 , 0 ) = kA_1;
	J( 5 , 1 ) = 0;
	J( 5 , 2 ) = 0;
	J( 5 , 3 ) = 0;
	J( 5 , 4 ) = 0;
	J( 5 , 5 ) = -D;

}
void Models::model_5(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	//Species order is: N_x N_c S_glu B_2 B_1 A_2 A_1 
	dydt[0] = -std::pow(y[4], n_omega_1)*y[0]*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_2)*y[1]*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = -y[3]*D + std::pow(KB_2, nB_2)*y[0]*kBmax_2/(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2));
	dydt[4] = std::pow(y[6], nB_1)*y[0]*kBmax_1/(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1)) - y[4]*D;
	dydt[5] = -y[5]*D + y[0]*kA_2;
	dydt[6] = -y[6]*D + y[0]*kA_1;

}
void Models::jac_5(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	J( 0 , 0 ) = -std::pow(y[4], n_omega_1)*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = std::pow(y[4], 2*n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*std::pow(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[4], n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 0 , 5 ) = 0;
	J( 0 , 6 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_2)*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*std::pow(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2), 2)) - std::pow(y[3], n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)));
	J( 1 , 4 ) = 0;
	J( 1 , 5 ) = 0;
	J( 1 , 6 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 2 , 5 ) = 0;
	J( 2 , 6 ) = 0;
	J( 3 , 0 ) = std::pow(KB_2, nB_2)*kBmax_2/(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2));
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = 0;
	J( 3 , 5 ) = 0;
	J( 3 , 6 ) = -std::pow(y[6], nB_2)*std::pow(KB_2, nB_2)*y[0]*kBmax_2*nB_2/(y[6]*std::pow(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2), 2));
	J( 4 , 0 ) = std::pow(y[6], nB_1)*kBmax_1/(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1));
	J( 4 , 1 ) = 0;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;
	J( 4 , 5 ) = 0;
	J( 4 , 6 ) = -std::pow(y[6], 2*nB_1)*y[0]*kBmax_1*nB_1/(y[6]*std::pow(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1), 2)) + std::pow(y[6], nB_1)*y[0]*kBmax_1*nB_1/(y[6]*(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1)));
	J( 5 , 0 ) = kA_2;
	J( 5 , 1 ) = 0;
	J( 5 , 2 ) = 0;
	J( 5 , 3 ) = 0;
	J( 5 , 4 ) = 0;
	J( 5 , 5 ) = -D;
	J( 5 , 6 ) = 0;
	J( 6 , 0 ) = kA_1;
	J( 6 , 1 ) = 0;
	J( 6 , 2 ) = 0;
	J( 6 , 3 ) = 0;
	J( 6 , 4 ) = 0;
	J( 6 , 5 ) = 0;
	J( 6 , 6 ) = -D;

}
void Models::model_6(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	//Species order is: N_x N_c S_glu B_2 B_1 A_2 A_1 
	dydt[0] = -std::pow(y[4], n_omega_1)*y[0]*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D*y[0] + y[0]*y[2]*mu_max_x/(K_x + y[2]);
	dydt[1] = -std::pow(y[3], n_omega_2)*y[1]*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D*y[1] + y[1]*y[2]*mu_max_c/(K_c + y[2]);
	dydt[2] = D*(S0_glu - y[2]) - y[1]*y[2]*mu_max_c/(g_c*(K_c + y[2])) - y[0]*y[2]*mu_max_x/(g_x*(K_x + y[2]));
	dydt[3] = std::pow(y[5], nB_2)*std::pow(KB_2, nB_2)*y[0]*kBmax_2/((std::pow(y[6], nB_2) + std::pow(KB_2, nB_2))*(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2))) - y[3]*D;
	dydt[4] = std::pow(y[6], nB_1)*y[0]*kBmax_1/(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1)) - y[4]*D;
	dydt[5] = -y[5]*D + y[1]*kA_2;
	dydt[6] = -y[6]*D + y[0]*kA_1;

}
void Models::jac_6(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double D = part_params[0];
	const double KB_1 = part_params[1];
	const double KB_2 = part_params[2];
	const double K_c = part_params[3];
	const double K_omega_1 = part_params[4];
	const double K_omega_2 = part_params[5];
	const double K_x = part_params[6];
	const double S0_glu = part_params[7];
	const double g_c = part_params[8];
	const double g_x = part_params[9];
	const double kA_1 = part_params[10];
	const double kA_2 = part_params[11];
	const double kBmax_1 = part_params[12];
	const double kBmax_2 = part_params[13];
	const double mu_max_c = part_params[14];
	const double mu_max_x = part_params[15];
	const double nB_1 = part_params[16];
	const double nB_2 = part_params[17];
	const double n_omega_1 = part_params[18];
	const double n_omega_2 = part_params[19];
	const double omega_max_1 = part_params[20];
	const double omega_max_2 = part_params[21];

	J( 0 , 0 ) = -std::pow(y[4], n_omega_1)*omega_max_1/(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)) - D + y[2]*mu_max_x/(K_x + y[2]);
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*y[2]*mu_max_x/std::pow(K_x + y[2], 2) + y[0]*mu_max_x/(K_x + y[2]);
	J( 0 , 3 ) = 0;
	J( 0 , 4 ) = std::pow(y[4], 2*n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*std::pow(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1), 2)) - std::pow(y[4], n_omega_1)*y[0]*n_omega_1*omega_max_1/(y[4]*(std::pow(y[4], n_omega_1) + std::pow(K_omega_1, n_omega_1)));
	J( 0 , 5 ) = 0;
	J( 0 , 6 ) = 0;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -std::pow(y[3], n_omega_2)*omega_max_2/(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)) - D + y[2]*mu_max_c/(K_c + y[2]);
	J( 1 , 2 ) = -y[1]*y[2]*mu_max_c/std::pow(K_c + y[2], 2) + y[1]*mu_max_c/(K_c + y[2]);
	J( 1 , 3 ) = std::pow(y[3], 2*n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*std::pow(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2), 2)) - std::pow(y[3], n_omega_2)*y[1]*n_omega_2*omega_max_2/(y[3]*(std::pow(y[3], n_omega_2) + std::pow(K_omega_2, n_omega_2)));
	J( 1 , 4 ) = 0;
	J( 1 , 5 ) = 0;
	J( 1 , 6 ) = 0;
	J( 2 , 0 ) = -y[2]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 1 ) = -y[2]*mu_max_c/(g_c*(K_c + y[2]));
	J( 2 , 2 ) = -D + y[1]*y[2]*mu_max_c/(g_c*std::pow(K_c + y[2], 2)) - y[1]*mu_max_c/(g_c*(K_c + y[2])) + y[0]*y[2]*mu_max_x/(g_x*std::pow(K_x + y[2], 2)) - y[0]*mu_max_x/(g_x*(K_x + y[2]));
	J( 2 , 3 ) = 0;
	J( 2 , 4 ) = 0;
	J( 2 , 5 ) = 0;
	J( 2 , 6 ) = 0;
	J( 3 , 0 ) = std::pow(y[5], nB_2)*std::pow(KB_2, nB_2)*kBmax_2/((std::pow(y[6], nB_2) + std::pow(KB_2, nB_2))*(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2)));
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;
	J( 3 , 4 ) = 0;
	J( 3 , 5 ) = -std::pow(y[5], 2*nB_2)*std::pow(KB_2, nB_2)*y[0]*kBmax_2*nB_2/(y[5]*(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2))*std::pow(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2), 2)) + std::pow(y[5], nB_2)*std::pow(KB_2, nB_2)*y[0]*kBmax_2*nB_2/(y[5]*(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2))*(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2)));
	J( 3 , 6 ) = -std::pow(y[6], nB_2)*std::pow(y[5], nB_2)*std::pow(KB_2, nB_2)*y[0]*kBmax_2*nB_2/(y[6]*std::pow(std::pow(y[6], nB_2) + std::pow(KB_2, nB_2), 2)*(std::pow(y[5], nB_2) + std::pow(KB_2, nB_2)));
	J( 4 , 0 ) = std::pow(y[6], nB_1)*kBmax_1/(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1));
	J( 4 , 1 ) = 0;
	J( 4 , 2 ) = 0;
	J( 4 , 3 ) = 0;
	J( 4 , 4 ) = -D;
	J( 4 , 5 ) = 0;
	J( 4 , 6 ) = -std::pow(y[6], 2*nB_1)*y[0]*kBmax_1*nB_1/(y[6]*std::pow(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1), 2)) + std::pow(y[6], nB_1)*y[0]*kBmax_1*nB_1/(y[6]*(std::pow(y[6], nB_1) + std::pow(KB_1, nB_1)));
	J( 5 , 0 ) = 0;
	J( 5 , 1 ) = kA_2;
	J( 5 , 2 ) = 0;
	J( 5 , 3 ) = 0;
	J( 5 , 4 ) = 0;
	J( 5 , 5 ) = -D;
	J( 5 , 6 ) = 0;
	J( 6 , 0 ) = kA_1;
	J( 6 , 1 ) = 0;
	J( 6 , 2 ) = 0;
	J( 6 , 3 ) = 0;
	J( 6 , 4 ) = 0;
	J( 6 , 5 ) = 0;
	J( 6 , 6 ) = -D;

}
