#include "ACEProblem.hpp"

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
   //#define PHASE
   #ifdef PHASE
   for(int i = 1; i < this->sizeX-1; i++)
   {
      for(int j = 1; j < this->sizeY-1; j++)
      {
         fu[j*sizeX + i] = get_rhs_phase_at(u, i, j);
      }
   }
   #endif
   
   #define CONCENTRATION
   #ifdef CONCENTRATION
   for(int i = 1; i < this->sizeX-1; i++)
   {
      for(int j = 1; j < this->sizeY-1; j++)
      {
         fu[sizeX*sizeY + j*sizeX + i] = get_rhs_concentration_at(u, i, j);
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
              << u[ j * sizeX + i ] << " " << u[ sizeX * sizeY + j * sizeX + i];
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
   #define VERSION_P 1

   double r1 = 0.5 - 0.5*ksi;
   double r2 = r1 + ksi;
   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));

         #if VERSION_P == 0
         //Hyperbolic tangent
         u[j*sizeX + i] = 1.0/2 * tanh(-3/ksi*(radius - r1)) + 1.0/2;

         #elif VERSION_P == 1
         //Linear by parts
         if( radius < r1 )
         {
            u[j*sizeX + i] = 1;
         }
         else if( radius < r2 )
         {
            u[j*sizeX + i] = 1 - (radius - r1) / (r2 - r1);
         }
         else
         {
            u[j*sizeX + i] = 0;
         }
         
         #elif VERSION_P == 2
         //Constant by parts
         if( radius < r1 )
         {
            u[j*sizeX + i] = 1;
         }
         else
         {
            u[j*sizeX + i] = 0;
         }

         #endif
      }
   }
}

void ACEProblem::set_concentration_initial_condition(double *u)
{
   #define VERSION_C 2
   int offset = sizeX * sizeY;
   double r1 = 0.5 - 0.5*ksi;
   double r2 = r1 + ksi;
   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {
         double init_conc = 0.025;
         
         #if VERSION_C == 0
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
         if(radius < r1)
            u[offset + j*sizeX + i] = init_conc;
         else if(radius < r2)
            u[offset + j*sizeX + i] = init_conc - init_conc * (radius - r1) / (r2 - r1);
         else
            u[offset + j*sizeX + i] = 0;

         #elif VERSION_C == 1
         double x_norm = abs(i*hx - (domain.x_right - domain.x_left)/2);
         if(x_norm < r1)
            u[offset + j*sizeX + i] = init_conc;
         else if(x_norm < r2)
            u[offset + j*sizeX + i] = init_conc - init_conc * (x_norm - r1) / (r2 - r1);
         else
            u[offset + j*sizeX + i] = 0;

         #elif VERSION_C == 2
         u[offset + j*sizeX + i] = 0;
         for(int n = 1; n < 5; n++)
         {
            double C_n = pow(1.0/2, n);
            double lambda_n = pow(n*M_PI/(domain.x_right-domain.x_left)*2, 2);
            u[offset + j*sizeX + i] += C_n * sin(sqrt(lambda_n) * (i*hx - (domain.x_right - domain.x_left)/2));
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
      fu[i] = 0;
      fu[(sizeY-1)*sizeX + i] = 0;
   }
   for(int j = 1; j < this->sizeY-1; j++)
   {
      fu[j*sizeX] = 0;
      fu[(j+1)*sizeX - 1] = 0;
   }
}

void ACEProblem::apply_concentration_boundary_condition(double *u, double *fu)
{
   int offset = sizeY * sizeX;

   #define VERSION_C_BOUND 1
   for(int i = 0; i < this->sizeX; i++)
   {
      #if VERSION_C_BOUND == 0
      u[offset + i] = 0;
      u[offset + (sizeY-1)*sizeX + i] = 0;
      fu[offset + i] = 0;
      fu[offset + (sizeY-1)*sizeX + i] = 0;
      #elif VERSION_C_BOUND == 1
      u[offset + i] = u[offset + i + sizeX];
      u[offset + (sizeY-1)*sizeX + i] = u[offset + (sizeY-2)*sizeX + i];
      fu[offset + i] = 0;
      fu[offset + (sizeY-1)*sizeX + i] = 0;
      #endif
   }
   for(int j = 0; j < this->sizeY; j++)
   {
      u[offset + j*sizeX] = 0;
      u[offset + (j+1)*sizeX - 1] = 0;
      fu[offset + j*sizeX] = 0;
      fu[offset + (j+1)*sizeX - 1] = 0;
   }
}

double ACEProblem::get_rhs_phase_at(double* u, int i, int j)
{
   double rhs = laplace(u, i, j) + f_0(u, i , j) / ksi / ksi + grad_norm(u, i, j)*F(u, i, j);
   /*
   double rhs = get_M_phi_tilde()*(get_epsilon_phi_tilde() * laplace(_u, i, j)
                                   - get_w_tilde()*get_q_prime(_u, i, j)
                                   + (get_G_alpha_tilde(_u, i, j) - get_G_beta_tilde(_u, i, j))*get_p_prime(_u, i, j));
   */
   return rhs;
}

double ACEProblem::get_rhs_concentration_at(double *u, int i, int j)
{
   double rhs = div_D_grad_concentration(u, i, j) + G(u, i, j);
   return rhs;
}

double ACEProblem::laplace(double *u, int i, int j)
{
   return (u[j*sizeX + i - 1] - 2*u[j*sizeX + i] + u[j*sizeX + i + 1])/hx/hx +
          (u[(j-1)*sizeX + i] - 2*u[j*sizeX + i] + u[(j+1)*sizeX + i])/hy/hy;
}

double ACEProblem::grad_norm(double *u, int i, int j)
{
   double derivative_x = (u[j*sizeX + i + 1] - u[j*sizeX + i - 1])/(2*hx);
   double derivative_y = (u[(j+1)*sizeX + i] - u[(j-1)*sizeX + i])/(2*hy);
   return sqrt(pow(derivative_x,2)+pow(derivative_y,2));
}

double ACEProblem::div_D_grad_concentration(double *u, int i, int j)
{
   int offset = sizeX * sizeY;
   double x_direction = get_diffusion_coef(u, i + 1, j) * (u[offset + j*sizeX + i + 1] - u[offset + j*sizeX + i]) / hx
                        - get_diffusion_coef(u, i, j) * (u[offset + j*sizeX + i] - u[offset + j*sizeX + i - 1]) / hx;
   double y_direction = get_diffusion_coef(u, i, j + 1) * (u[offset + (j+1)*sizeX + i] - u[offset + j*sizeX + i]) / hy
                        - get_diffusion_coef(u, i, j) * (u[offset + j*sizeX + i] - u[offset + (j-1)*sizeX + i]) / hy;
   return x_direction / hx + y_direction / hy;
}

double ACEProblem::get_diffusion_coef(double *u, int i, int j)
{
    return 1;
}

double ACEProblem::f_0(double *_u, int i, int j)
{
   return par_a*_u[j*sizeX + i]*(1 - _u[j*sizeX + i])*(_u[j*sizeX + i] - 1.0/2.0);
}

double ACEProblem::F(double *u, int i, int j)
{
   double r = sqrt(pow(i*hx + domain.x_left, 2) + pow(j*hy + domain.y_left, 2));
   return 2/(std::max(r, 0.1));
}

double ACEProblem::G(double *u, int i, int j)
{
   return 0;
}

double ACEProblem::get_M_phi_tilde()
{
   double l = 0.9*5e-9;
   double b = 3.23e-10;
   double D_nb_alpha = 6.6e-10*exp(-15851.4/T);
   double D_nb_beta = 9e-9*pow(T/1136, 18.1)*exp(-(25100+35.5*(T-1136))/(1.98*T));
   double D_eff = pow(sqrt(D_nb_alpha) + sqrt(D_nb_beta), 2);
   double M_phi = l*l*0.0235*D_eff/(D_nb_alpha*b*b);
   return M_phi;
}

double ACEProblem::get_epsilon_phi_tilde()
{
   double delta_0 = 5e-9;
   double l = 0.9 * delta_0;
   double sigma_0 = 0.3;
   double R = 8.31446261815324;
   double V_m = 1.4060e-5;
   double epsilon = sqrt(6*sigma_0*delta_0);
   double epsilon_tilde = epsilon/(l*sqrt(R*T/V_m));
   return epsilon_tilde;
}

double ACEProblem::get_G_alpha_tilde(const double *_u, int i, int j)
{  
   double c = _u[sizeX*sizeY + j*sizeX + i];
   double R = 8.31446261815324;
   double G_0_zr_alpha = -7827.595 + 125.64905*T - 24.1618*T*log(T) - 4.37791e-3*T*T + 34971/T;
   double G_0_nb_alpha = 1480.647 + 144.445475*T - 26.4711*T*log(T) + 2.03475e-4*T*T - 3.5012e-7*T*T*T + 93399/T;
   double L_0_alpha = 24411;

   double G_alpha = c*G_0_nb_alpha + (1 - c)*G_0_zr_alpha + R*T*c*log(c) + R*T*(1-c)*log(1-c) + c*(1-c)*L_0_alpha;
   return G_alpha/T/R;
}

double ACEProblem::get_G_beta_tilde(const double *_u, int i, int j)
{
   double c = _u[sizeX*sizeY + j*sizeX + i];
   double R = 8.31446261815324;
   double G_0_zr_beta = -525.539 + 124.9457*T - 25.607406*T*log(T) - 3.40084E-4*T*T - 9.729e-9*T*T*T + 25233/T - 7.6143E-11*T*T*T*T;
   double G_0_nb_beta = -8519.353 + 142.045475*T - 26.4711*T*log(T) + 2.03475e-4*T*T - 3.5012e-7*T*T*T + 93399/T;
   double L_0_beta = 15911 + 3.35*T;
   double L_0_i_beta = 3919 - 1.091*T;

   double G_beta = c*G_0_nb_beta + (1 - c)*G_0_zr_beta + R*T*c*log(c) + R*T*(1-c)*log(1-c) + c*(1-c)*(L_0_beta + L_0_i_beta*(2*c-1));
   return G_beta/T/R;
}

double ACEProblem::get_w_tilde()
{
   double sigma_0 = 0.3;
   double delta_0 = 5e-9;
   double R = 8.31446261815324;
   double w = 3*sigma_0/delta_0;
   return w/R/T;
}

double ACEProblem::get_p_prime(double *_u, int i, int j)
{  
   double q = pow(_u[j*sizeX + i], 2) - 2*pow(_u[j*sizeX + i], 3) + pow(_u[j*sizeX + i], 4);
   return 30*q;
}

double ACEProblem::get_q_prime(double *_u, int i, int j)
{
   return 2*_u[j*sizeX + i] - 6*pow(_u[j*sizeX + i], 2) + 4*pow(_u[j*sizeX + i],3);
}
