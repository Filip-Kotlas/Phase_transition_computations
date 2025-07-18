#include "ACEProblem.hpp"

#define COMPUTE_PHASE
#define COMPUTE_CONCENTRATION

#define C_INIT 0
/*
*  0 - Linear by parts in circle around the middle
*  1 - Constant on the whole domain
*  2 - Fourier along x axis
*  3 - Fourier along y axis
*/

#define P_INIT 1
/*
*  All based on radius
*  0 - Hyperbolic tangent
*  1 - Linear by parts
*  2 - Constant by parts
*/

#define C_BOUND 3
/*
*  0 - Dirichlet everywhere
*  1 - Neumann along x, Dirichlet along y
*  2 - Dirichlet along x, Neumann along y
*  3 - Neumann everywhere
*/

#define P_BOUND 1
/*
*  0 - Dirichlet
*  1 - Neumann
*/

//#define C_TEST_X

#ifdef C_TEST_X
#define C_INIT 2
#define C_BOUND 1
#endif

#ifdef C_TEST_Y
#define C_INIT 3
#define C_BOUND 2
#endif

#define FORCE 0
/*
*  0 - Force equal 40
*  1 - Force inversely proportional to the distance from the middle
*/

ACEProblem::ACEProblem(int sizeX,
                       int sizeY,
                       Domain domain,
                       double alpha,
                       double beta,
                       double par_a,
                       double ksi,
                       MODEL model,
                       std::string output_folder)
:  sizeX(sizeX),
   sizeY(sizeY),
   domain(domain),
   hx((domain.x_right - domain.x_left)/(sizeX-1)),
   hy((domain.y_right - domain.y_left)/(sizeY-1)),
   alpha(alpha),
   beta(beta),
   par_a(par_a),
   ksi(ksi),
   model(model),
   output_folder(output_folder)
{
}

int ACEProblem::getDegreesOfFreedom()
{
    return 2 * this->sizeX * this->sizeY;
}

void ACEProblem::getRightHandSide(const double &t, double *u, double *fu)
{
   apply_concentration_physical_condition(u);
   
   #ifdef COMPUTE_PHASE
   for(int i = 1; i < this->sizeX-1; i++)
   {
      for(int j = 1; j < this->sizeY-1; j++)
      {
         fu[j*sizeX + i] = get_rhs_phase_at(u, i, j);
      }
   }
   #endif
   
   #ifdef COMPUTE_CONCENTRATION

   for(int i = 1; i < this->sizeX-1; i++)
   {
      for(int j = 1; j < this->sizeY-1; j++)
      {
         fu[sizeX*sizeY + j*sizeX + i] = get_rhs_concentration_at(t, u, i, j);
      }
   }
   #endif
   apply_boundary_condition(u, fu);
}

bool ACEProblem::writeSolution(const double &t, int step, const double *u)
{
   /****
    * Filename with step index
    */
   std::stringstream str;
   str << output_folder << "\\ACE-equation-" << std::setw( 5 ) << std::setfill( '0' ) << step << ".txt";
   
   /****
    * Open file
    */
   std::fstream file;
   file.open( str.str(), std::fstream::out | std::fstream::trunc );
   if( ! file )
   {
      std::cerr << "Unable to open the file " << str.str() << std::endl;
      return false;
   }

   /****
    * Write solution
    */
   for( int j = 0; j < sizeY; j++ )
   {
      for( int i = 0; i < sizeX; i++ )
      {
         file << domain.x_left + i * hx << " " << domain.y_left + j * hy << " "
              << phase_at(u, i, j) << " " << conc_at(u, i, j);
         file << std::endl;
      }
      file << std::endl;
   }
   return true;
}

void ACEProblem::setInitialCondition(double *u)
{
   set_phase_initial_condition(u);
   set_concentration_initial_condition(u);
}

void ACEProblem::set_phase_initial_condition(double *u)
{
   double r1 = (domain.x_right - domain.x_left)/4 - 0.5*ksi;
   double r2 = r1 + ksi;
   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));

         #if P_INIT == 0
         //Hyperbolic tangent
         u[j*sizeX + i] = 1.0/2 * tanh(-3/ksi*(radius - r1)) + 1.0/2; //needs to be changed according to the correct phases

         #elif P_INIT == 1
         //Linear by parts
         if( radius < r1 )
         {
            u[j*sizeX + i] = constants::p_alpha;
         }
         else if( radius < r2 )
         {
            u[j*sizeX + i] = constants::p_alpha - (constants::p_alpha - constants::p_beta)*(radius - r1) / (r2 - r1);
         }
         else
         {
            u[j*sizeX + i] = constants::p_beta;
         }
         
         #elif P_INIT == 2
         //Constant by parts
         if( radius < r1 )
         {
            u[j*sizeX + i] = constants::p_alpha;
         }
         else
         {
            u[j*sizeX + i] = constants::p_beta;
         }

         #endif
      }
   }
}

void ACEProblem::set_concentration_initial_condition(double *u)
{
   int offset = sizeX * sizeY;
   double r1 = (domain.x_right - domain.x_left)/4 - 0.5*ksi;
   double r2 = r1 + ksi;

   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {         
         #if C_INIT == 0
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
         if(radius < r1)
            u[offset + j*sizeX + i] = constants::c_init_alpha;
         else if(radius < r2)
            u[offset + j*sizeX + i] = constants::c_init_alpha - (constants::c_init_alpha - constants::c_init_beta) * (radius - r1) / (r2 - r1);
         else
            u[offset + j*sizeX + i] = constants::c_init_beta;

         #elif C_INIT == 1
         u[offset + j*sizeX + i] = constants::c_init_alpha;

         #elif C_INIT == 2
         u[offset + j*sizeX + i] = 0;
         for(int n = 1; n < 5; n++)
         {
            double C_n = pow(1.0/2, n);
            double lambda_n = pow(n*M_PI/(domain.x_right-domain.x_left), 2);
            u[offset + j*sizeX + i] += C_n * sin(sqrt(lambda_n) * (i*hx));
         }

         #elif C_INIT == 3
         u[offset + j*sizeX + i] = 0;
         for(int n = 1; n < 5; n++)
         {
            double C_n = pow(1.0/2, n);
            double lambda_n = pow(n*M_PI/(domain.y_right-domain.y_left), 2);
            u[offset + j*sizeX + i] += C_n * sin(sqrt(lambda_n) * (j*hy));
         }
         #endif
      }
   }
}

void ACEProblem::apply_boundary_condition(double *u, double *fu)
{
   apply_phase_boundary_condition(u, fu);
   apply_concentration_boundary_condition(u, fu);
}

void ACEProblem::apply_phase_boundary_condition(double *u, double *fu)
{
   for(int i = 0; i < this->sizeX; i++)
   {
      //Dirichlet
      #if P_BOUND == 0
      fu[i] = 0;
      fu[(sizeY-1)*sizeX + i] = 0;

      //Neumann
      #elif P_BOUND == 1
      u[i] = u[i + sizeX];
      u[(sizeY-1)*sizeX + i] = u[(sizeY-2)*sizeX + i];
      fu[i] = 0;
      fu[(sizeY-1)*sizeX + i] = 0;
      #endif
   }
   for(int j = 1; j < this->sizeY-1; j++)
   {
      //Dirichlet
      #if P_BOUND == 0
      fu[j*sizeX] = 0;
      fu[(j+1)*sizeX - 1] = 0;

      //Neumann
      #elif P_BOUND == 1
      u[j*sizeX] = u[j*sizeX + 1];
      u[(j+1)*sizeX - 1] = u[(j+1)*sizeX - 2];
      fu[j*sizeX] = 0;
      fu[(j+1)*sizeX - 1] = 0;
      #endif
   }
}

void ACEProblem::apply_concentration_boundary_condition(double *u, double *fu)
{
   int offset = sizeY * sizeX;

   //Boundary conditions along x direction
   for(int i = 0; i < this->sizeX; i++)
   {
      //Dirichlet
      #if C_BOUND == 0 || C_BOUND == 2
      fu[offset + i] = 0;
      fu[offset + (sizeY-1)*sizeX + i] = 0;
      
      //Neumann
      #elif C_BOUND == 1 || C_BOUND == 3
      u[offset + i] = u[offset + i + sizeX];
      u[offset + (sizeY-1)*sizeX + i] = u[offset + (sizeY-2)*sizeX + i];
      fu[offset + i] = 0;
      fu[offset + (sizeY-1)*sizeX + i] = 0;

      #endif
   }

   //Boundary conditions along y direction
   for(int j = 0; j < this->sizeY; j++)
   {
      //Dirichlet
      #if C_BOUND == 0 || C_BOUND == 1
      fu[offset + j*sizeX] = 0;
      fu[offset + (j+1)*sizeX - 1] = 0;

      //Neumann
      #elif C_BOUND == 2 || C_BOUND == 3
      u[offset + j*sizeX] = u[offset + j*sizeX + 1];
      u[offset + (j+1)*sizeX - 1] = u[offset + (j+1)*sizeX - 2];
      fu[offset + j*sizeX] = 0;
      fu[offset + (j+1)*sizeX - 1] = 0;
      
      #endif
   }

}

void ACEProblem::apply_concentration_physical_condition(double *u)
{
   for(int i = 0; i < this->sizeX; i++)
   {
      for(int j = 0; j < this->sizeY; j++)
      {
         // Check if the concentration is in allowed range of (c_min, c_max).
         if(u[sizeX*sizeY + j*sizeX + i] < constants::c_min)
         {
            u[sizeX*sizeY + j*sizeX + i] = constants::c_min;
         }
         else if (u[sizeX*sizeY + j*sizeX + i] > constants::c_max)
         {
            u[sizeX*sizeY + j*sizeX + i] = constants::c_max;
         }
      }
   }
}

double ACEProblem::get_rhs_phase_at(const double* u, int i, int j)
{
   double rhs = 0.0;
   double c = conc_at(u, i, j);
   double D = 100000;

   if(model == MODEL::MODEL_3)
      rhs = laplace(u, i, j) + f_0(u, i , j) / ksi / ksi + grad_norm(u, i, j)*F(u, i, j);
   
   else if(model == MODEL::MODEL_4)
      /*
      if(i == sizeX/2 && j < sizeY/2)
      {
         std::cout << j << ", "
                   << "Laplace: " << constants::M_phi_tilde(T) * pow(constants::epsilon_tilde(T), 2)
                   << ", f_0: " << 4*constants::M_phi_tilde(T)*constants::w_tilde(T)
                   << ", F: " << constants::M_phi_tilde(T)*(constants::G_m_alpha(c, T)/constants::R/T - constants::G_m_beta(c, T)/constants::R/T)
                   << std::endl;
      }
      rhs = constants::M_phi_tilde(T)*(pow(constants::epsilon_tilde(T), 2) * laplace(u, i, j)
            + der_polynom_p(u, i, j) * (constants::G_m_alpha(c, T) - constants::G_m_beta(c, T)) / constants::R / T
            + der_polynom_q(u, i, j) * constants::w_tilde(T));
      rhs *= D;
      */

      rhs = laplace(u, i, j)
            + der_polynom_p(u, i, j) * (constants::G_m_beta(c, T) - constants::G_m_alpha(c, T)) / constants::R / T * D
            + der_polynom_q(u, i, j) * 500;
      /*
      if(i == sizeX/2 && j < sizeY/2)
      {
         std::cout << i << ", " << j << ": "
                   << "f_0: " << 2500
                   << ", F: " << (constants::G_m_alpha(c, T) - constants::G_m_beta(c, T)) / constants::R / T * D
                   << ", c: " << c
                   << std::endl;
      }
      */
   return rhs;
}

double ACEProblem::get_rhs_concentration_at(const double &t, double *u, int i, int j)
{  
   return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
}

double ACEProblem::laplace(const double *u, int i, int j)
{
   return (phase_at(u, i - 1, j) - 2*phase_at(u, i, j) + phase_at(u, i + 1, j))/hx/hx +
          (phase_at(u, i, j - 1) - 2*phase_at(u, i, j) + phase_at(u, i, j + 1))/hy/hy;
}

double ACEProblem::grad_norm(const double *u, int i, int j)
{
   double derivative_x = (phase_at(u, i + 1, j) - phase_at(u, i - 1, j))/(2*hx);
   double derivative_y = (phase_at(u, i, j + 1) - phase_at(u, i, j - 1))/(2*hy);
   return sqrt(pow(derivative_x,2) + pow(derivative_y,2));
}

double ACEProblem::div_D_grad_concentration(const double *u, int i, int j)
{
   //if(i == sizeX/2)
      //std::cout <<"c: " << conc_at(u, i, j) << ", p: " << phase_at(u, i, j) << ", 2g w c: " << sec_deriv_of_g_w_resp_to_c(u, i, j) << std::endl;

	double coeff_plus_half = (get_conc_diff_coef(u, i + 1, j) + get_conc_diff_coef(u, i, j)) / 2;
	double coeff_minus_half = (get_conc_diff_coef(u, i, j) + get_conc_diff_coef(u, i - 1, j)) / 2;
	double x_direction = coeff_plus_half * (conc_at(u, i+1, j) - conc_at(u, i, j)) / hx
				         - coeff_minus_half * (conc_at(u, i, j) - conc_at(u, i-1, j)) / hx;

	coeff_plus_half = (get_conc_diff_coef(u, i, j+1) + get_conc_diff_coef(u, i, j)) / 2;
	coeff_minus_half = (get_conc_diff_coef(u, i, j) + get_conc_diff_coef(u, i, j-1)) / 2;
	double y_direction = coeff_plus_half * (conc_at(u, i, j+1) - conc_at(u, i, j)) / hy
                      	 - coeff_minus_half * (conc_at(u, i, j) - conc_at(u, i, j-1)) / hy;

	return x_direction / hx + y_direction / hy;
}

double ACEProblem::div_D_grad_phase(const double *u, int i, int j)
{
   //if(i == sizeX/2)
      //std::cout <<"c: " << conc_at(u, i, j) << ", p: " << phase_at(u, i, j) <<  ", g w p c: " << deriv_of_g_w_resp_to_c_and_p(u, i, j) << std::endl;
	double coeff_plus_half = (get_phas_diff_coef(u, i + 1, j) + get_phas_diff_coef(u, i, j)) / 2;
	double coeff_minus_half = (get_phas_diff_coef(u, i, j) + get_phas_diff_coef(u, i - 1, j)) / 2;
	double x_direction = coeff_plus_half * (phase_at(u, i+1, j) - phase_at(u, i, j)) / hx
				         - coeff_minus_half * (phase_at(u, i, j) - phase_at(u, i-1, j)) / hx;

	coeff_plus_half = (get_phas_diff_coef(u, i, j+1) + get_phas_diff_coef(u, i, j)) / 2;
	coeff_minus_half = (get_phas_diff_coef(u, i, j) + get_phas_diff_coef(u, i, j-1)) / 2;
	double y_direction = coeff_plus_half * (phase_at(u, i, j+1) - phase_at(u, i, j)) / hy
                      	 - coeff_minus_half * (phase_at(u, i, j) - phase_at(u, i, j-1)) / hy;

	return x_direction / hx + y_direction / hy;
}

double ACEProblem::get_conc_diff_coef(const double *u, int i, int j)
{
   double D = 0.0000125;
   if(!blowout && sec_deriv_of_g_w_resp_to_c(u, i, j) < -0.001)
   {
      std::cout << "Blow out conc: (" << i << ", " <<  j << std::fixed << std::setprecision(8) << "), c: "
                << conc_at(u, i, j) << ", p: " << phase_at(u, i, j) << ", "
                << pow(constants::M_Nb_beta(T), 1 - polynom_p(u, i, j)) / pow(constants::M_Nb_alpha(T), 1 - polynom_p(u, i, j)) << ", " 
		          << sec_deriv_of_g_w_resp_to_c(u, i, j) << std::endl;
      blowout = true;
   }
   return D
          * conc_at(u, i, j)
		    * (1 - conc_at(u, i, j))
          * pow(constants::M_Nb_beta(T), 1 - polynom_p(u, i, j))
		    / pow(constants::M_Nb_alpha(T), 1 - polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

double ACEProblem::get_phas_diff_coef(const double *u, int i, int j)
{
   double D = 0.00000900;
   if(!blowout && deriv_of_g_w_resp_to_c_and_p(u, i, j) < - 0.001)
   {
      std::cout << "Blow out phase: (" << i << ", " <<  j << "), c: "
                << conc_at(u, i, j) << ", p: " << phase_at(u, i, j) << ", "
                << pow(constants::M_Nb_beta(T), 1 - polynom_p(u, i, j)) / pow(constants::M_Nb_alpha(T), 1 - polynom_p(u, i, j)) << ", " 
		          << deriv_of_g_w_resp_to_c_and_p(u, i, j) << std::endl;
      blowout = true;
   }
   return D * conc_at(u, i, j)
		    * (1 - conc_at(u, i, j))
		    * pow(constants::M_Nb_beta(T), 1 - polynom_p(u, i, j))
		    / pow(constants::M_Nb_alpha(T), 1 - polynom_p(u, i, j))
		    * deriv_of_g_w_resp_to_c_and_p(u, i, j);
}

double ACEProblem::f_0(const double *u, int i, int j)
{
   return par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - 1.0/2.0);
}

double ACEProblem::F(const double *u, int i, int j)
{
   #if FORCE == 0
   return 40;

   #elif FORCE == 1
   double mid_x = (domain.x_right - domain.x_left)/2;
   double mid_y = (domain.y_right - domain.y_left)/2;
   double r = sqrt(pow(i*hx - mid_x, 2) + pow(j*hy - mid_y, 2));
   return 2/std::max(r, 0.1);
   
   #endif
}

double ACEProblem::G(const double &t, double *u, int i, int j)
{
   double coefficients[5] = {0, 1.0/2, 0, 1.0/4, 0};
   double G = 0; 
   for(int n = 1; n < 5; n++)
   {
      double F_n = coefficients[n];
      
      #ifdef C_TEST_X
      double lambda_n = pow(n*M_PI/(domain.x_right-domain.x_left), 2);
      G += F_n * sin(sqrt(lambda_n) * (i*hx));
      #endif

      #ifdef C_TEST_Y
      double lambda_n = pow(n*M_PI/(domain.y_right-domain.y_left), 2);
      G += F_n * sin(sqrt(lambda_n) * (j*hy));
      #endif
   }
   return G;
}

double ACEProblem::grade_4_polynom(const double *u, int i, int j)
{
   return pow(phase_at(u, i, j), 2) * pow(phase_at(u, i, j) - 1.0, 2);
}

double ACEProblem::polynom_p(const double *u, int i, int j)
{
   return 6*pow(phase_at(u, i, j), 5) - 15*pow(phase_at(u, i, j), 4) + 10*pow(phase_at(u, i, j), 3);
}

double ACEProblem::der_polynom_p(const double *u, int i, int j)
{
   return 30*(pow(phase_at(u, i, j), 4) - 2*pow(phase_at(u, i, j), 3) + pow(phase_at(u, i, j), 2));
}

double ACEProblem::der_polynom_q(const double *u, int i, int j)
{
   return 4 * phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - 1.0/2.0);
}

double ACEProblem::sec_deriv_of_g_w_resp_to_c(const double* u, int i, int j)
{
   double c = conc_at(u, i, j);
   double d2_G_alpha_wrt_c = constants::R * T / (c*(1-c))
   							 - 2 * constants::L_0_alpha;
   double d2_G_beta_wrt_c = constants::R * T / (c*(1-c))
   							 - 2 * constants::L_0_beta(T)
							 + (6 - 12*c) * constants::L_0_i_beta(T);
   return polynom_p(u, i, j)*d2_G_alpha_wrt_c + (1 - polynom_p(u, i, j))*d2_G_beta_wrt_c;
}

double ACEProblem::deriv_of_g_w_resp_to_c_and_p(const double* u, int i, int j)
{
   double c = conc_at(u, i, j);
   double d_G_alpha_wrt_c = constants::G_Nb_alpha_0(T)
                            - constants::G_Zr_alpha_0(T)
                            + constants::R * T * (log(c) - log(1-c))
                            + (1 - 2*c) * constants::L_0_alpha;
   double d_G_beta_wrt_c = constants::G_Nb_beta_0(T)
						   - constants::G_Zr_beta_0(T)
						   + constants::R * T * (log(c) - log(1-c))
						   + (1 - 2 * c) * constants::L_0_beta(T)
						   + (6*c - 6*pow(c, 2) - 1) * constants::L_0_i_beta(T);

   return der_polynom_p(u, i, j)*(d_G_alpha_wrt_c - d_G_beta_wrt_c);
}