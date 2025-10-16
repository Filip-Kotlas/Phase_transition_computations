#include "Problem.hpp"

#define COMPUTE_PHASE
#define COMPUTE_CONCENTRATION

#define INIT 7

#define P_INIT INIT


#define C_INIT INIT
/*
*  20 - Fourier along x axis
*  30 - Fourier along y axis
*/

#define C_BOUND 3
/*
*  0 - Dirichlet everywhere
*  1 - Neumann along x, Dirichlet along y
*  2 - Dirichlet along x, Neumann along y
*  3 - Neumann everywhere
*  4 - Dirichlet everywhere with c = 0.007
*/

#define P_BOUND 1
/*
*  0 - Dirichlet
*  1 - Neumann
*/

//#define C_TEST_X

#ifdef C_TEST_X
#define C_INIT 20
#define C_BOUND 1
#endif

#ifdef C_TEST_Y
#define C_INIT 30
#define C_BOUND 2
#endif

#define FORCE 2
/*
*  0 - Force equal 40
*  1 - Force inversely proportional to the distance from the middle
*  2 - Force for zirconium model
*/

Problem::Problem(Parameters param)
:  sizeX(param.sizeX),
   sizeY(param.sizeY),
   domain(param.domain),
   hx((domain.x_right - domain.x_left)/(sizeX-1)),
   hy((domain.y_right - domain.y_left)/(sizeY-1)),
   alpha(param.alpha),
   par_a(param.par_a),
   par_b(param.par_b),
   par_d(param.par_d),
   T(param.T),
   ksi(param.ksi),
   model(param.model)
{
}

Index Problem::getDegreesOfFreedom()
{
    return 2 * this->sizeX * this->sizeY;
}

__cuda_callable__
void Problem::set_rhs_at(const VectorView& u, VectorView& fu, Index i, Index j)
{
    Index offset = this->sizeX * this->sizeY;
    fu[j*this->sizeX + i] = get_rhs_phase_at(u, i, j);
    fu[offset + j*this->sizeX + i] = get_rhs_concentration_at(u, i, j);
}

bool Problem::writeSolution(const Real &t, Index step, const VectorView& u, const std::filesystem::path& output_folder)
{
   /****
    * Filename with step index
    */
   std::stringstream str;
   str << "Result-" << std::setw( 5 ) << std::setfill( '0' ) << step << ".txt";
   std::filesystem::path file_path = output_folder / str.str();

   /****
    * Open file
    */
   std::fstream file;
   file.open( file_path, std::fstream::out | std::fstream::trunc );
   if( ! file )
   {
      std::cerr << "Unable to open the file " << str.str() << std::endl;
      return false;
   }

   /****
    * Write solution
    */
   file << std::scientific << std::setprecision(15);
   for( Index j = 0; j < sizeY; j++ )
   {
      for( Index i = 0; i < sizeX; i++ )
      {
         file << domain.x_left + i * hx << " " << domain.y_left + j * hy << " "
              << phase_at(u, i, j) << " " << conc_at(u, i, j);
         file << std::endl;
      }
      file << std::endl;
   }
   return true;
}

void Problem::set_init_cond_manually(Vector& u, InitialCondition ic)
{
   set_phase_initial_condition(u, ic);
   set_concentration_initial_condition(u, ic);
}

void Problem::set_phase_initial_condition(Vector& u, InitialCondition ic)
{
    Real r = (domain.x_right - domain.x_left)/6;
    Real r1 = r - 0.5*ksi;
    Real r2 = r1 + ksi;

    switch(ic) {
        case InitialCondition::HyperbolicTangent:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    value = 1.0/2 * tanh(-3/ksi*(radius - r1)) + 1.0/2; //needs to be changed according to the correct phases
                } );
            break;

        case InitialCondition::LinearByParts:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    if( radius < r1 )
                    {
                        value = constants::p_alpha;
                    }
                    else if( radius < r2 )
                    {
                        value = constants::p_alpha - (constants::p_alpha - constants::p_beta)*(radius - r1) / (r2 - r1);
                    }
                    else
                    {
                        value = constants::p_beta;
                    }
                } );
            break;
        
        case InitialCondition::ConstantCircle:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    if( radius < r1 )
                    {
                        value = constants::p_alpha;
                    }
                    else
                    {
                        value = constants::p_beta;
                    }
                } );
            break;
        
        case InitialCondition::ConstantHalves:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    //Index j = ind / sizeX;
                    if( i < sizeX/2 )
                    {
                        value = constants::p_alpha;
                    }
                    else
                    {
                        value = constants::p_beta;
                    }
                } );
            break;

        case InitialCondition::Stripe:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX;
                    if( i*hx < 0.2 )
                    {
                        value = constants::p_alpha;
                    }
                    else
                    {
                        value = constants::p_beta;
                    }
                } );
            break;

        case InitialCondition::TwoBumps:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX;
                    Real y = 2*j*hy;
                    if( (y < 1 && i*hx < y*y*(1-2*y+y*y)/0.625+0.1) ||
                        (y >= 1 && i*hx < (y-1)*(y-1)*(1-2*(y-1)+(y-1)*(y-1))/0.625+0.1))
                    {
                        value = constants::p_alpha;
                    }
                    else
                    {
                        value = constants::p_beta;
                    }
                } );
            break;
        
        case InitialCondition::Star:
            u.forElements( 0, sizeX*sizeY,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    Real phi = atan( (j*hy - (domain.y_right-domain.y_left)/2) / (i*hx - (domain.x_right - domain.x_left)/2));
                    if ((i*hx - (domain.x_right - domain.x_left)/2) < 0 )
                        phi = phi + M_PI;
                    else if ( (i*hx - (domain.x_right - domain.x_left)/2) > 0 && (j*hy - (domain.y_right-domain.y_left)/2) < 0)
                        phi = phi + 2 * M_PI;
                    
                    if ( radius < 0.15 + 0.1 * sin(6 * phi) )
                    {
                        value = constants::p_alpha;
                    }
                    else
                    {
                        value = constants::p_beta;
                    }
                } );
            break;

            case InitialCondition::Perpendicular_Stripes:
                u.forElements( 0, sizeX*sizeY,
                    [ = ] __cuda_callable__( Index ind, Real & value )
                    {
                        Index i = ind % sizeX;
                        Index j = ind / sizeX;
                        if ( (i*hx < 0.1) || (j*hy < 0.1 ) )
                        {
                            value = constants::p_alpha;
                        }
                        else
                        {
                            value = constants::p_beta;
                        }
                    } );
                break;
    }   
}

void Problem::set_concentration_initial_condition(Vector& u, InitialCondition ic)
{
    Index offset = sizeX * sizeY;
    Real r = (domain.x_right - domain.x_left)/6;
    Real r1 = r - 0.5*ksi;
    Real r2 = r1 + ksi;

    switch(ic) {
        case InitialCondition::HyperbolicTangent:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    value = constants::c_init_alpha; // TODO: Nutno předělat je to špatně.
                });
            break;
        
        case InitialCondition::LinearByParts:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    if( radius < r1 )
                    {
                        value = constants::c_init_alpha;
                    }
                    else if( radius < r2 )
                    {
                        value = constants::c_init_alpha - (constants::c_init_alpha - constants::c_init_beta)*(radius - r1) / (r2 - r1);
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;

        case InitialCondition::ConstantCircle:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    if( radius < r )
                    {
                        value = constants::c_init_alpha;
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;
        
        case InitialCondition::ConstantHalves:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    //Index j = ind / sizeX - sizeY;
                    if( i < sizeX/2 )
                    {
                        value = constants::c_init_alpha;
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;
        
        case InitialCondition::Stripe:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    if( i*hx < 0.2 )
                    {
                        value = constants::c_init_alpha;
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;
        
        case InitialCondition::TwoBumps:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    Real y = 2*j*hy;
                    if( (y < 1 && i*hx < y*y*(1-2*y+y*y)/0.625+0.1) ||
                        (y >= 1 && i*hx < (y-1)*(y-1)*(1-2*(y-1)+(y-1)*(y-1))/0.625+0.1))
                    {
                        value = constants::c_init_alpha;
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;

        case InitialCondition::Star:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    Real radius = sqrt(pow(i*hx - (domain.x_right - domain.x_left)/2, 2) + pow(j*hy - (domain.y_right-domain.y_left)/2, 2));
                    Real phi = atan( (j*hy - (domain.y_right-domain.y_left)/2) / (i*hx - (domain.x_right - domain.x_left)/2));
                    if ((i*hx - (domain.x_right - domain.x_left)/2) < 0 )
                        phi = phi + M_PI;
                    else if ( (i*hx - (domain.x_right - domain.x_left)/2) > 0 && (j*hy - (domain.y_right-domain.y_left)/2) < 0)
                        phi = phi + 2 * M_PI;
                    
                    if ( radius < 0.15 + 0.1 * sin(6 * phi) )
                    {
                        value = constants::c_init_alpha;
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;
        
        case InitialCondition::FourierX:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    value = 0;
                    for(Index n = 1; n < 5; n++)
                    {
                        Real C_n = pow(1.0/2, n);
                        Real lambda_n = pow(n*M_PI/(domain.x_right-domain.x_left), 2);
                        value += C_n * sin(sqrt(lambda_n) * (i*hx));
                    }
                } );
            break;

        case InitialCondition::FourierY:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    value = 0;
                    for(Index n = 1; n < 5; n++)
                    {
                        Real C_n = pow(1.0/2, n);
                        Real lambda_n = pow(n*M_PI/(domain.y_right-domain.y_left), 2);
                        value += C_n * sin(sqrt(lambda_n) * (j*hy));
                    }
                } );
            break;

        case InitialCondition::Perpendicular_Stripes:
            u.forElements( sizeX*sizeY, sizeX*sizeY*2,
                [ = ] __cuda_callable__( Index ind, Real & value )
                {
                    Index i = ind % sizeX;
                    Index j = ind / sizeX - sizeY;
                    if ( (i*hx < 0.1) || (j*hy < 0.1 ) )
                    {
                        value = constants::c_init_alpha;
                    }
                    else
                    {
                        value = constants::c_init_beta;
                    }
                } );
            break;
   }
}

bool Problem::set_init_cond_from_file(Vector& u, const std::filesystem::path& filename)
{
    std::ifstream file(filename);
    if (!file)
    {
        std::cerr << "Nelze otevřít soubor: " << filename << '\n';
        return false;
    }

    Vector host_u(getDegreesOfFreedom());

    std::string line;
    Index count = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue; // přeskoč prázdné řádky

        std::istringstream iss(line);
        Real x, y, p, c;
        if (!(iss >> x >> y >> p >> c)) {
            std::cout << "Chyba při čtení řádku: " << line << '\n';
            return false;
        }

        if (count >= sizeX * sizeY) {
            std::cout << "Soubor obsahuje více dat než očekáváno.\n";
            return false;
        }

        host_u[count] = p;
        host_u[sizeX * sizeY + count] = c;
        ++count;
    }

    if (count < sizeX * sizeY) {
        std::cout << "Soubor obsahuje méně dat než očekáváno.\n";
        return false;
    }

    u = host_u;

   return true;
}

__cuda_callable__
void Problem::apply_boundary_condition_xdir(Index ind, VectorView u, VectorView fu)
{
    apply_phase_boundary_condition_xdir(ind, u, fu);
    apply_concentration_boundary_condition_xdir(ind, u, fu);
}

__cuda_callable__
void Problem::apply_boundary_condition_ydir(Index ind, VectorView u, VectorView fu)
{
    apply_phase_boundary_condition_ydir(ind, u, fu);
    apply_concentration_boundary_condition_ydir(ind, u, fu);
}

__cuda_callable__
void Problem::apply_phase_boundary_condition_xdir(Index ind, VectorView u, VectorView fu)
{
    // Horní a dolní hrana
    fu[ind] = 0;
    fu[(sizeY - 1) * sizeX + ind] = 0;

    #if P_BOUND == 0  // Dirichlet
    u[ind] = 0;
    u[(sizeY - 1) * sizeX + ind] = 0;

    #elif P_BOUND == 1  // Neumann
    u[ind] = u[ind + sizeX];
    u[(sizeY - 1) * sizeX + ind] = u[(sizeY - 2) * sizeX + ind];
    #endif
}

__cuda_callable__
void Problem::apply_phase_boundary_condition_ydir(Index ind, VectorView u, VectorView fu)
{
    // Levá a pravá hrana
    fu[ind * sizeX] = 0;
    fu[(ind + 1) * sizeX - 1] = 0;

    #if P_BOUND == 0  // Dirichlet
    u[ind * sizeX] = 0;
    u[(ind + 1) * sizeX - 1] = 0;

    #elif P_BOUND == 1  // Neumann
    u[ind * sizeX] = u[ind * sizeX + 1];
    u[(ind + 1) * sizeX - 1] = u[(ind + 1) * sizeX - 2];
    #endif
}

__cuda_callable__
void Problem::apply_concentration_boundary_condition_xdir(Index ind, VectorView u, VectorView fu)
{
    Index offset = sizeX * sizeY;

    fu[offset + ind] = 0;
    fu[offset + (sizeY - 1) * sizeX + ind] = 0;

    #if C_BOUND == 0 || C_BOUND == 2
    u[offset + ind] = 0;
    u[offset + (sizeY - 1) * sizeX + ind] = 0;

    #elif C_BOUND == 1 || C_BOUND == 3
    u[offset + ind] = u[offset + ind + sizeX];
    u[offset + (sizeY - 1) * sizeX + ind] = u[offset + (sizeY - 2) * sizeX + ind];

    #elif C_BOUND == 4
    u[offset + ind] = constants::c_init_beta;
    u[offset + (sizeY - 1) * sizeX + ind] = constants::c_init_beta;
    #endif
}

__cuda_callable__
void Problem::apply_concentration_boundary_condition_ydir(Index ind, VectorView u, VectorView fu)
{
    Index offset = sizeX * sizeY;

    fu[offset + ind * sizeX] = 0;
    fu[offset + (ind + 1) * sizeX - 1] = 0;

    #if C_BOUND == 0 || C_BOUND == 1
    u[offset + ind * sizeX] = 0;
    u[offset + (ind + 1) * sizeX - 1] = 0;

    #elif C_BOUND == 2 || C_BOUND == 3
    u[offset + ind * sizeX] = u[offset + ind * sizeX + 1];
    u[offset + (ind + 1) * sizeX - 1] = u[offset + (ind + 1) * sizeX - 2];

    #elif C_BOUND == 4
    u[offset + ind * sizeX] = constants::c_init_beta;
    u[offset + (ind + 1) * sizeX - 1] = constants::c_init_beta;
    #endif
}


/*
void Problem::apply_boundary_condition(double *u, double *fu)
{
   apply_phase_boundary_condition(u, fu);
   apply_concentration_boundary_condition(u, fu);
}

void Problem::apply_phase_boundary_condition(double *u, double *fu)
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

void Problem::apply_concentration_boundary_condition(double *u, double *fu)
{
   int offset = sizeY * sizeX;

   //Boundary conditions along x direction
   for(int i = 0; i < this->sizeX; i++)
   {
      //Dirichlet 0
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

      //Dirichlet with c_init_beta
      #elif C_BOUND == 4
      u[offset + i] = constants::c_init_beta;
      u[offset + (sizeY-1)*sizeX + i] = constants::c_init_beta;
      fu[offset + i] = constants::c_init_beta;
      fu[offset + (sizeY-1)*sizeX + i] = constants::c_init_beta;
      #endif
   }

   //Boundary conditions along y direction
   for(int j = 0; j < this->sizeY; j++)
   {
      //Dirichlet 0
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
      
      //Dirichlet with c_init_beta
      #elif C_BOUND == 4
      u[offset + j*sizeX] = constants::c_init_beta;
      u[offset + (j+1)*sizeX - 1] = constants::c_init_beta;
      fu[offset + j*sizeX] = constants::c_init_beta;
      fu[offset + (j+1)*sizeX - 1] = constants::c_init_beta;
      #endif
   }

}

void Problem::apply_concentration_physical_condition(double *u)
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
*/

__cuda_callable__
Real Problem::get_rhs_phase_at(const VectorView& u, Index i, Index j)
{
    return 1.0/alpha * (laplace(u, i, j) + f_0(u, i , j) / ksi / ksi - par_b/ksi*grade_4_polynom(u, i, j)*F(u, i, j));
}

__cuda_callable__
Real Problem::get_rhs_concentration_at(const VectorView& u, Index i, Index j)
{
    return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
}

__cuda_callable__
Real Problem::laplace(const VectorView& u, Index i, Index j)
{
   return (phase_at(u, i - 1, j) - 2*phase_at(u, i, j) + phase_at(u, i + 1, j))/hx/hx +
          (phase_at(u, i, j - 1) - 2*phase_at(u, i, j) + phase_at(u, i, j + 1))/hy/hy;
}

__cuda_callable__
Real Problem::div_D_grad_concentration(const VectorView& u, Index i, Index j)
{
	Real coeff_plus_half = (get_conc_diff_coef(u, i + 1, j) + get_conc_diff_coef(u, i, j)) / 2;
	Real coeff_minus_half = (get_conc_diff_coef(u, i, j) + get_conc_diff_coef(u, i - 1, j)) / 2;
	Real x_direction = coeff_plus_half * (conc_at(u, i+1, j) - conc_at(u, i, j)) / hx
				         - coeff_minus_half * (conc_at(u, i, j) - conc_at(u, i-1, j)) / hx;

	coeff_plus_half = (get_conc_diff_coef(u, i, j+1) + get_conc_diff_coef(u, i, j)) / 2;
	coeff_minus_half = (get_conc_diff_coef(u, i, j) + get_conc_diff_coef(u, i, j-1)) / 2;
	Real y_direction = coeff_plus_half * (conc_at(u, i, j+1) - conc_at(u, i, j)) / hy
                      	 - coeff_minus_half * (conc_at(u, i, j) - conc_at(u, i, j-1)) / hy;

	return x_direction / hx + y_direction / hy;
}

__cuda_callable__
Real Problem::div_D_grad_phase(const VectorView& u, Index i, Index j)
{
	Real coeff_plus_half = (get_phas_diff_coef(u, i + 1, j) + get_phas_diff_coef(u, i, j)) / 2;
	Real coeff_minus_half = (get_phas_diff_coef(u, i, j) + get_phas_diff_coef(u, i - 1, j)) / 2;
	Real x_direction = coeff_plus_half * (phase_at(u, i+1, j) - phase_at(u, i, j)) / hx
				         - coeff_minus_half * (phase_at(u, i, j) - phase_at(u, i-1, j)) / hx;

	coeff_plus_half = (get_phas_diff_coef(u, i, j+1) + get_phas_diff_coef(u, i, j)) / 2;
	coeff_minus_half = (get_phas_diff_coef(u, i, j) + get_phas_diff_coef(u, i, j-1)) / 2;
	Real y_direction = coeff_plus_half * (phase_at(u, i, j+1) - phase_at(u, i, j)) / hy
                      	 - coeff_minus_half * (phase_at(u, i, j) - phase_at(u, i, j-1)) / hy;

	return x_direction / hx + y_direction / hy;
}

__cuda_callable__
Real Problem::get_conc_diff_coef(const VectorView& u, Index i, Index j)
{
    return par_b*par_d/30.0
            * conc_at(u, i, j)
            * (1 - conc_at(u, i, j))
            * constants::M_Nb_beta(T)
		    * pow(constants::M_Nb_alpha(T)/constants::M_Nb_beta(T), polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

__cuda_callable__
Real Problem::get_phas_diff_coef(const VectorView& u, Index i, Index j)
{
   return par_b*par_d/30.0
          * conc_at(u, i, j)
		    * (1 - conc_at(u, i, j))
          * constants::M_Nb_beta(T)
		    * pow(constants::M_Nb_alpha(T)/constants::M_Nb_beta(T), polynom_p(u, i, j))
		    * deriv_of_g_w_resp_to_c_and_p(u, i, j);
}

__cuda_callable__
Real Problem::f_0(const VectorView& u, Index i, Index j)
{
   return par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - 1.0/2.0);
}

__cuda_callable__
Real Problem::F(const VectorView& u, Index i, Index j)
{
   #if FORCE == 0
   return 20;

   #elif FORCE == 1
   Real mid_x = (domain.x_right - domain.x_left)/2;
   Real mid_y = (domain.y_right - domain.y_left)/2;
   Real r = sqrt(pow(i*hx - mid_x, 2) + pow(j*hy - mid_y, 2));
   return -2/std::max(r, 0.1);

   #elif FORCE == 2
   return constants::G_m_alpha(conc_at(u, i, j), T) - constants::G_m_beta(conc_at(u, i, j), T);
   
   #endif
}

/*
double Problem::G(const double &t, double *u, int i, int j)
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
*/

__cuda_callable__
Real Problem::grade_4_polynom(const VectorView& u, Index i, Index j)
{
   return pow(phase_at(u, i, j), 2) * pow(phase_at(u, i, j) - 1.0, 2);
}

__cuda_callable__
Real Problem::polynom_p(const VectorView& u, Index i, Index j)
{
   return 6*pow(phase_at(u, i, j), 5) - 15*pow(phase_at(u, i, j), 4) + 10*pow(phase_at(u, i, j), 3);
}

__cuda_callable__
Real Problem::der_polynom_p(const VectorView& u, Index i, Index j)
{
   return 30*(pow(phase_at(u, i, j), 4) - 2*pow(phase_at(u, i, j), 3) + pow(phase_at(u, i, j), 2));
}

__cuda_callable__
Real Problem::sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j)
{
   Real c = conc_at(u, i, j);
   Real d2_G_alpha_wrt_c = constants::R * T / (c*(1-c))
   							 - 2 * constants::L_0_alpha;
   Real d2_G_beta_wrt_c = constants::R * T / (c*(1-c))
   							 - 2 * constants::L_0_beta(T)
							 + (6 - 12*c) * constants::L_0_i_beta(T);
   return polynom_p(u, i, j)*d2_G_alpha_wrt_c + (1 - polynom_p(u, i, j))*d2_G_beta_wrt_c;
}

__cuda_callable__
Real Problem::deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j)
{
   Real c = conc_at(u, i, j);
   Real d_G_alpha_wrt_c = constants::G_Nb_alpha_0(T)
                            - constants::G_Zr_alpha_0(T)
                            + constants::R * T * (log(c) - log(1-c))
                            + (1 - 2*c) * constants::L_0_alpha;
   Real d_G_beta_wrt_c = constants::G_Nb_beta_0(T)
						   - constants::G_Zr_beta_0(T)
						   + constants::R * T * (log(c) - log(1-c))
						   + (1 - 2 * c) * constants::L_0_beta(T)
						   + (6*c - 6*pow(c, 2) - 1) * constants::L_0_i_beta(T);

   return der_polynom_p(u, i, j)*(d_G_alpha_wrt_c - d_G_beta_wrt_c);
}