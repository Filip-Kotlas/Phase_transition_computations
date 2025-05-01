#include "ACEProblem.hpp"

#define COMPUTE_PHASE
#define COMPUTE_CONCENTRATION

#define C_INIT 0
/*
*  0 - Linear by parts in circle around the middle
*  1 -
*  2 - Fourier along x axis
*  3 - Fourier along y axis.
*/

#define P_INIT 1
/*
*  All based on radius
*  0 - Hyperbolic tangent
*  1 - Linear by parts
*  2 - Constant by parts
*/

#define C_BOUND 0
/*
*  0 - Dirichlet everywhere
*  1 - Neumann along x, Dirichlet along y
*  2 - Dirichlet along x, Neumann along y
*  3 - Neumann everywhere
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

#define FORCE 1
/*
*  0 - Force equal 1
*  1 - Force inversly proportionate to the distance from the middle
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

   double r1 = (domain.x_right - domain.x_left)/4 - 0.5*ksi;
   double r2 = r1 + ksi;
   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));

         #if P_INIT == 0
         //Hyperbolic tangent
         u[j*sizeX + i] = 1.0/2 * tanh(-3/ksi*(radius - r1)) + 1.0/2;

         #elif P_INIT == 1
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
         
         #elif P_INIT == 2
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
   int offset = sizeX * sizeY;
   double r1 = (domain.x_right - domain.x_left)/4 - 0.5*ksi;
   double r2 = r1 + ksi;

   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {
         double init_conc = 0.025;
         
         #if C_INIT == 0
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
         if(radius < r1)
            u[offset + j*sizeX + i] = init_conc;
         else if(radius < r2)
            u[offset + j*sizeX + i] = init_conc - init_conc * (radius - r1) / (r2 - r1);
         else
            u[offset + j*sizeX + i] = 0;

         #elif C_INIT == 1
         double x_norm = abs(i*hx - (domain.x_right - domain.x_left)/2);
         if(x_norm < r1)
            u[offset + j*sizeX + i] = init_conc;
         else if(x_norm < r2)
            u[offset + j*sizeX + i] = init_conc - init_conc * (x_norm - r1) / (r2 - r1);
         else
            u[offset + j*sizeX + i] = 0;

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

   //Boundary coniditons along x direction
   for(int i = 0; i < this->sizeX; i++)
   {
      //Dirichlet
      #if C_BOUND == 0 || C_BOUND == 2
      u[offset + i] = 0;
      u[offset + (sizeY-1)*sizeX + i] = 0;
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

   //Boundary conditions along y conditions
   for(int j = 0; j < this->sizeY; j++)
   {
      //Dirichlet
      #if C_BOUND == 0 || C_BOUND == 1
      u[offset + j*sizeX] = 0;
      u[offset + (j+1)*sizeX - 1] = 0;
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

double ACEProblem::get_rhs_phase_at(double* u, int i, int j)
{
   double rhs = 0.0;

   if(model == MODEL::MODEL_3)
      rhs = laplace(u, i, j) + f_0(u, i , j) / ksi / ksi + grad_norm(u, i, j)*F(u, i, j);
   
   else if(model == MODEL::MODEL_4)
      rhs = laplace(u, i, j) + f_0(u, i , j) / ksi / ksi + 10/sqrt(8)*sqrt(par_a)*1.0/ksi*grade_4_polynom(u, i, j)*F(u, i, j);

   return rhs;
}

double ACEProblem::get_rhs_concentration_at(const double &t, double *u, int i, int j)
{
   return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
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
   double x_direction = get_conc_diff_coef(u, i + 1, j) * (u[offset + j*sizeX + i + 1] - u[offset + j*sizeX + i]) / hx
                        - get_conc_diff_coef(u, i, j) * (u[offset + j*sizeX + i] - u[offset + j*sizeX + i - 1]) / hx;
   double y_direction = get_conc_diff_coef(u, i, j + 1) * (u[offset + (j+1)*sizeX + i] - u[offset + j*sizeX + i]) / hy
                        - get_conc_diff_coef(u, i, j) * (u[offset + j*sizeX + i] - u[offset + (j-1)*sizeX + i]) / hy;
   return x_direction / hx + y_direction / hy;
}

double ACEProblem::div_D_grad_phase(double *u, int i, int j)
{
   double x_direction = get_phas_diff_coef(u, i + 1, j) * (u[j*sizeX + i + 1] - u[j*sizeX + i]) / hx
                        - get_phas_diff_coef(u, i, j) * (u[j*sizeX + i] - u[j*sizeX + i - 1]) / hx;
   double y_direction = get_phas_diff_coef(u, i, j + 1) * (u[(j+1)*sizeX + i] - u[j*sizeX + i]) / hy
                        - get_phas_diff_coef(u, i, j) * (u[j*sizeX + i] - u[+ (j-1)*sizeX + i]) / hy;
   return x_direction / hx + y_direction / hy;
}

double ACEProblem::get_conc_diff_coef(double *u, int i, int j)
{
   return 1;
}

double ACEProblem::get_phas_diff_coef(double *u, int i, int j)
{
   return 0.01;
}

double ACEProblem::f_0(double *_u, int i, int j)
{
   return par_a*_u[j*sizeX + i]*(1 - _u[j*sizeX + i])*(_u[j*sizeX + i] - 1.0/2.0);
}

double ACEProblem::F(double *u, int i, int j)
{
   #if FORCE == 0
   return 1;

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

double ACEProblem::grade_4_polynom(double *u, int i, int j)
{
   return pow(u[j*sizeX + i], 2) * pow(u[j*sizeX + i] - 1.0, 2);
}