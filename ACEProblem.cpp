#include "ACEProblem.h"

ACEProblem::ACEProblem(int sizeX, int sizeY, double alpha, double sigma, double ksi, double T)
:  sizeX(sizeX),
   sizeY(sizeY),
   hx(1.0/(sizeX-1)),
   hy(1.0/(sizeY-1)),
   alpha(alpha),
   sigma(sigma),
   ksi(ksi),
   T(T)
{
}

int ACEProblem::getDegreesOfFreedom()
{
    return this->sizeX * this->sizeY;
}

void ACEProblem::getRightHandSide(const double &t, double *_u, double *fu)
{
   for(int i = 0; i < this->sizeX; i++)
   {
      for(int j = 0; j < this->sizeY; j++)
      {
         fu[j*sizeX + i] = sigma/alpha*laplaceD(_u, i, j) + 1/ksi/ksi/alpha*f_0(_u, i, j);
      }
   }
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
         file << i * hx << " " << j * hy << " " << u[ j * sizeX + i ];
         file << std::endl;
      }
      file << std::endl;
   }
   return true;
}

void ACEProblem::setInitialCondition(double *u)
{
   double r1 = 0.20;
   double r2 = 0.30;
   for(int i = 0; i < this->sizeX; i++)
   {
      for(int j = 0; j < this->sizeY; j++)
      {
         // Tangenc hyperbolický
         double radius = sqrt( pow(i - this->sizeX/2, 2) + pow(j - this->sizeY/2, 2)) / std::min(sizeX, sizeY);
         u[j*this->sizeX + i] = 1.0/2 * tanh(- 20 * (radius - 0.25)) + 1.0/2;
         //std::cout << i << ", " << j << "= " << tanh(- 20 * (radius - 0.25)) << std::endl;


         /*
         //Lineární po částech.
         double radius = sqrt( pow(i - this->sizeX/2, 2) + pow(j - this->sizeY/2, 2)) / std::min(sizeX, sizeY);
         if( radius < r1 )
         {
            u[j*this->sizeX + i] = 1;
         }
         else if( radius < r2 )
         {
            u[j*this->sizeX + i] = 1 - (radius - r1) / (r2 - r1);
         }
         else
         {
            u[j*this->sizeX + i] = 0;
         }
         */
      }
   }
}

double ACEProblem::laplaceD(double *_u, int i, int j)
{
   //shift to accomodate boundary, boundary is free
   int I_shift = 0;
   int J_shift = 0;

   if(i == 0)
      I_shift = 1;
   if(i == sizeX-1)
      I_shift = -1;
   if(j == 0)
      J_shift = 1;
   if(j == sizeY-1)
      J_shift = -1;

   return (_u[j*sizeX + i - 1 + I_shift] - 2*_u[j*sizeX + i + I_shift] + _u[j*sizeX + i + 1 + I_shift])/hx/hx +
          (_u[(j-1+J_shift)*sizeX + i] - 2*_u[(j+J_shift)*sizeX + i] + _u[(j+1+J_shift)*sizeX + i])/hy/hy;
}

double ACEProblem::f_0(double *_u, int i, int j)
{
   return _u[j*sizeX + i]*(1 - _u[j*sizeX + i])*(_u[j*sizeX + i] - 1.0/2.0);
}

double ACEProblem::F(double *_u, int i, int j)
{
   return 1/(sqrt(pow(i-0.5, 2) + pow(j-0.5, 2)));
}

double get_M()
{
   
}