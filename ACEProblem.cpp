#include "ACEProblem.h"

ACEProblem::ACEProblem(int sizeX,
                       int sizeY,
                       Domain domain,
                       double alpha,
                       double beta,
                       double par_a,
                       double ksi,
                       MODEL model)
:  sizeX(sizeX),
   sizeY(sizeY),
   domain(domain),
   hx((domain.x_right - domain.x_left)/(sizeX-1)),
   hy((domain.y_right - domain.y_left)/(sizeY-1)),
   alpha(alpha),
   beta(beta),
   par_a(par_a),
   ksi(ksi),
   model(model)
{
}

int ACEProblem::getDegreesOfFreedom()
{
    return this->sizeX * this->sizeY;
}

void ACEProblem::getRightHandSide(const double &t, double *_u, double *fu)
{
   for(int i = 1; i < this->sizeX-1; i++)
   {
      for(int j = 1; j < this->sizeY-1; j++)
      {
         fu[j*sizeX + i] = right_hand_side_at(_u, i, j);
      }
   }
   set_dirichlet_boundary(_u, fu);
}

bool ACEProblem::writeSolution(const double &t, int step, const double *u)
{
   /****
    * Filename with step index
    */
   std::stringstream str;
   str << "Results\\ACE-equation-" << std::setw( 5 ) << std::setfill( '0' ) << step << ".txt";
   
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
         file << domain.x_left + i * hx << " " << domain.y_left + j * hy << " " << u[ j * sizeX + i ];
         file << std::endl;
      }
      file << std::endl;
   }
   return true;
}

void ACEProblem::setInitialCondition(double *u)
{
   #define VERSION 1

   double r1 = 0.5 - 0.5*ksi;
   double r2 = r1 + ksi;
   for(int i = 0; i < sizeX; i++)
   {
      for(int j = 0; j < sizeY; j++)
      {
         double radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));

         #if VERSION == 0
         //Hyperbolic tangent
         u[j*sizeX + i] = 1.0/2 * tanh(-3/ksi*(radius - r1)) + 1.0/2;

         #elif VERSION == 1
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
         
         #elif VERSION == 2
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

void ACEProblem::set_dirichlet_boundary(double* _u, double* fu)
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

double ACEProblem::right_hand_side_at(double* _u, int i, int j)
{
   try
   {
      double rhs = 0;
      if(model == MODEL::MODEL_1)
      {
         double b = 1.0/6*sqrt(par_a/2);
         rhs = 1.0/alpha*laplace(_u, i, j)
               + 1.0/ksi/ksi/alpha*f_0(_u, i, j)
               + b*beta/ksi/alpha*F(_u, i, j); // TODO: Add physical condition.
      }
      else if(model == MODEL::MODEL_2)
      {
         rhs = 1.0/alpha*laplace(_u, i, j)
               + 1.0/ksi/ksi/alpha*f_0(_u, i, j)
               + 1.0/alpha*grad_norm(_u, i, j)*F(_u, i, j); // TODO:: Fix this
      }
      else if(model == MODEL::MODEL_3)
      {
         rhs = 1.0/alpha*laplace(_u, i, j)
               + 1.0/ksi/ksi/alpha*f_0(_u, i, j)
               + beta/alpha*grad_norm(_u, i, j)*F(_u, i, j);
      }
      else
      {
         throw model;
      }
      return rhs;
   }
   catch (MODEL model)
   {
      std::cout << "Unrecognized model. The model: " << int(model) << std::endl;
      assert(true);
      return 0;
   }
}

double ACEProblem::laplace(double *_u, int i, int j)
{
   return (_u[j*sizeX + i - 1] - 2*_u[j*sizeX + i] + _u[j*sizeX + i + 1])/hx/hx +
          (_u[(j-1)*sizeX + i] - 2*_u[j*sizeX + i] + _u[(j+1)*sizeX + i])/hy/hy;
}

double ACEProblem::grad_norm(double *_u, int i, int j)
{
   double derivative_x = (_u[j*sizeX + i + 1] - _u[j*sizeX + i - 1])/(2*hx);
   double derivative_y = (_u[(j+1)*sizeX + i] - _u[(j-1)*sizeX + i])/(2*hy);
   return sqrt(pow(derivative_x,2)+pow(derivative_y,2));
}

double ACEProblem::f_0(double *_u, int i, int j)
{
   return par_a*_u[j*sizeX + i]*(1 - _u[j*sizeX + i])*(_u[j*sizeX + i] - 1.0/2.0);
}

double ACEProblem::F(double *_u, int i, int j)
{
   return 2/(sqrt(pow(i*hx + domain.x_left, 2) + pow(j*hy + domain.y_left, 2)));
}
