#include "Problem.hpp"
#include "npy.hpp"

#define COMPUTE_PHASE
#define COMPUTE_CONCENTRATION

#define INIT 7

#define P_INIT INIT


#define C_INIT INIT
/*
*  20 - Fourier along x axis
*  30 - Fourier along y axis
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

#define FORCE 1
/*
*  0 - Force equal 0
*  1 - Force inversely proportional to the distance from the middle
*  2 - Force for zirconium model
*/

Problem::Problem(Parameters param)
:  sizeX(param.sizeX),
   sizeY(param.sizeY),
   domain(param.domain),
   hx((domain.x_right - domain.x_left)/(sizeX-1)),
   hy((domain.y_right - domain.y_left)/(sizeY-1)),
   param(param)
{}

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
   const bool binary = param.output_format == OutputFormat::Npy;

   /****
    * Filename with step index
    */
   std::stringstream str;
   str << "Result-" << std::setw( 5 ) << std::setfill( '0' ) << step << ( binary ? ".npy" : ".txt" );
   std::filesystem::path file_path = output_folder / str.str();

   /****
    * Open file
    */
   std::ofstream file( file_path, std::ios::out | std::ios::trunc | ( binary ? std::ios::binary : std::ios::openmode() ) );
   if( ! file )
   {
      std::cerr << "Unable to open the file " << str.str() << std::endl;
      return false;
   }

   /****
    * The whole solution is copied from the GPU at once. Reading it element by element
    * (u.getElement) would mean one small device-to-host transfer per element.
    */
   TNL::Containers::Vector< Real, TNL::Devices::Host, Index > host_u( u.getSize() );
   host_u = u;
   const Index offset = sizeX * sizeY;

   /****
    * Write solution
    */
   if( binary )
   {
      // Array of the shape (2, sizeY, sizeX): the phase block followed by the concentration block,
      // which is exactly the layout of the state vector (index j*sizeX + i).
      npy::write_header( file, npy::descriptor< Real >(),
                         { 2, static_cast< std::size_t >( sizeY ), static_cast< std::size_t >( sizeX ) } );
      file.write( reinterpret_cast< const char* >( host_u.getData() ),
                  static_cast< std::streamsize >( host_u.getSize() ) * static_cast< std::streamsize >( sizeof( Real ) ) );
   }
   else
   {
      file << std::scientific << std::setprecision(15);
      for( Index j = 0; j < sizeY; j++ )
      {
         for( Index i = 0; i < sizeX; i++ )
         {
            file << domain.x_left + i * hx << " " << domain.y_left + j * hy << " "
                 << host_u[ j * sizeX + i ] << " " << host_u[ offset + j * sizeX + i ];
            file << '\n';
         }
         file << '\n';
      }
   }

   file.close();
   if( ! file )
   {
      std::cerr << "Error while writing the file " << str.str() << std::endl;
      return false;
   }
   return true;
}

void Problem::set_init_cond_manually(Vector& u, InitialCondition& init_cond)
{
    u.forElements( 0, sizeX*sizeY,
        [ = ] __cuda_callable__( Index ind, Real& value )
        {
            value = init_cond.get_phase(ind);
        } );

    u.forElements( sizeX*sizeY, sizeX*sizeY*2,
        [ = ] __cuda_callable__( Index ind, Real& value )
        {
            value = init_cond.get_concentration(ind - sizeX*sizeY);
        } );
}

bool Problem::set_init_cond_from_file(Vector& u, const std::filesystem::path& filename, const bool scaling)
{
    // The format is given by the extension: .npy is the binary array written by writeSolution,
    // anything else is the text format "x y phase conc".
    const bool binary = filename.extension() == ".npy";
    std::ifstream file(filename, binary ? std::ios::in | std::ios::binary : std::ios::in);
    if (!file)
    {
        std::cerr << "Unable to open the file: " << filename << '\n';
        return false;
    }

    const Index N = sizeX * sizeY;

    // Values from the file, the point k = j*sizeX + i is on the position k.
    std::vector<Real> file_phase(N);
    std::vector<Real> file_conc(N);

    if (binary)
    {
        try {
            std::string descr;
            bool fortran_order = false;
            std::vector<std::size_t> shape;
            npy::read_header(file, descr, fortran_order, shape);

            const std::vector<std::size_t> expected_shape = { 2, static_cast<std::size_t>(sizeY), static_cast<std::size_t>(sizeX) };
            if (descr != npy::descriptor<Real>() || fortran_order || shape != expected_shape) {
                std::cout << "The .npy file has to be a C-ordered array of the type " << npy::descriptor<Real>()
                          << " and the shape (2, " << sizeY << ", " << sizeX << ").\n";
                return false;
            }
        }
        catch (const std::runtime_error& e) {
            std::cout << "Error while reading the .npy file: " << e.what() << '\n';
            return false;
        }

        const std::streamsize bytes = static_cast<std::streamsize>(N) * static_cast<std::streamsize>(sizeof(Real));
        file.read(reinterpret_cast<char*>(file_phase.data()), bytes);
        file.read(reinterpret_cast<char*>(file_conc.data()), bytes);
        if (!file) {
            std::cout << "The file contains less data than expected.\n";
            return false;
        }
        if (file.peek() != std::ifstream::traits_type::eof()) {
            std::cout << "The file contains more data than expected.\n";
            return false;
        }
    }
    else
    {
        std::string line;
        Index count = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue; // skip empty lines

            std::istringstream iss(line);
            Real x, y, p, c;
            if (!(iss >> x >> y >> p >> c)) {
                std::cout << "Error while reading line: " << line << '\n';
                return false;
            }

            if (count >= N) {
                std::cout << "The file contains more data than expected.\n";
                return false;
            }

            file_phase[count] = p;
            file_conc[count] = c;
            count++;
        }

        if (count < N) {
            std::cout << "The file contains less data than expected.\n";
            return false;
        }
    }

    TNL::Containers::Vector<Real, TNL::Devices::Host, Index> host_u(getDegreesOfFreedom());

    for (int k = 0; k < N; ++k)
    {
        host_u[k] = static_cast<int>(constants::Phase::beta); // phase initialization
        host_u[N + k] = constants::c_init_beta; // concentration initialization
    }

    const Index start_i = sizeX / 4;
    const Index start_j = sizeY / 4;

    for (Index count = 0; count < N; ++count) {
        if (scaling) {
            int i = count % sizeX;
            int j = count / sizeX;
            if (i % 2 == 0 && j % 2 == 0) {
                Index scale_index = start_i + (i / 2) + (start_j + (j / 2)) * sizeX;
                host_u[scale_index] = file_phase[count];
                host_u[N + scale_index] = file_conc[count];
            }
        }
        else {
            host_u[count] = file_phase[count];
            host_u[N + count] = file_conc[count];
        }
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
    // Top and bottom edges
    fu[ind] = 0;
    fu[(sizeY - 1) * sizeX + ind] = 0;

    if (param.bc_phase_x == BCType::Dirichlet) {
        u[ind] = 0;
        u[(sizeY - 1) * sizeX + ind] = 0;
    }
    else if (param.bc_phase_x == BCType::Neumann) {
        u[ind] = u[ind + sizeX];
        u[(sizeY - 1) * sizeX + ind] = u[(sizeY - 2) * sizeX + ind];
    }
}

__cuda_callable__
void Problem::apply_phase_boundary_condition_ydir(Index ind, VectorView u, VectorView fu)
{
    // Left and right edges
    fu[ind * sizeX] = 0;
    fu[(ind + 1) * sizeX - 1] = 0;

    if (param.bc_phase_y == BCType::Dirichlet) {
        u[ind * sizeX] = 0;
        u[(ind + 1) * sizeX - 1] = 0;
    }
    else if (param.bc_phase_y == BCType::Neumann) {
        u[ind * sizeX] = u[ind * sizeX + 1];
        u[(ind + 1) * sizeX - 1] = u[(ind + 1) * sizeX - 2];
    }
}

__cuda_callable__
void Problem::apply_concentration_boundary_condition_xdir(Index ind, VectorView u, VectorView fu)
{
    Index offset = sizeX * sizeY;

    fu[offset + ind] = 0;
    fu[offset + (sizeY - 1) * sizeX + ind] = 0;

    if (param.bc_conc_x == BCType::Dirichlet) {
        u[offset + ind] = constants::c_init_beta;
        u[offset + (sizeY - 1) * sizeX + ind] = constants::c_init_beta;
    }
    else if (param.bc_conc_x == BCType::Neumann) {
        u[offset + ind] = u[offset + ind + sizeX];
        u[offset + (sizeY - 1) * sizeX + ind] = u[offset + (sizeY - 2) * sizeX + ind];
    }
}

__cuda_callable__
void Problem::apply_concentration_boundary_condition_ydir(Index ind, VectorView u, VectorView fu)
{
    Index offset = sizeX * sizeY;

    fu[offset + ind * sizeX] = 0;
    fu[offset + (ind + 1) * sizeX - 1] = 0;

    if (param.bc_conc_y == BCType::Dirichlet) {
        u[offset + ind * sizeX] = constants::c_init_beta;
        u[offset + (ind + 1) * sizeX - 1] = constants::c_init_beta;
    }
    else if (param.bc_conc_y == BCType::Neumann) {
        u[offset + ind * sizeX] = u[offset + ind * sizeX + 1];
        u[offset + (ind + 1) * sizeX - 1] = u[offset + (ind + 1) * sizeX - 2];
    }
}

__cuda_callable__
Real Problem::get_rhs_phase_at(const VectorView& u, Index i, Index j)
{
    return Real(1.0)/param.alpha * (div_T0(u, i, j) + f_0(u, i , j) / param.ksi / param.ksi - param.par_b/param.ksi*grade_4_polynom(u, i, j)*F(u, i, j));
}

__cuda_callable__
Real Problem::get_rhs_concentration_at(const VectorView& u, Index i, Index j)
{
    if (param.force_term_type == FTType::Zirconium)
        return div_D_grad_concentration(u, i, j) + div_D_grad_phase(u, i, j);
    else
        return 0;
}

__cuda_callable__
Real Problem::laplace(const VectorView& u, Index i, Index j)
{
   return (phase_at(u, i - 1, j) - 2*phase_at(u, i, j) + phase_at(u, i + 1, j))/hx/hx +
          (phase_at(u, i, j - 1) - 2*phase_at(u, i, j) + phase_at(u, i, j + 1))/hy/hy;
}

__cuda_callable__
Real Problem::grad_p_1_central(const VectorView& u, Index i, Index j)
{
    return (phase_at(u, i+1, j) - phase_at(u, i-1, j))/(2*hx);
}

__cuda_callable__
Real Problem::grad_p_2_central(const VectorView& u, Index i, Index j)
{
    return (phase_at(u, i, j+1) - phase_at(u, i, j-1))/(2*hy);
}

__cuda_callable__
Real Problem::grad_p_1_forward(const VectorView& u, Index i, Index j)
{
    return (phase_at(u, i+1, j) - phase_at(u, i, j))/hx;
}

__cuda_callable__
Real Problem::grad_p_2_forward(const VectorView& u, Index i, Index j)
{
   return (phase_at(u, i, j+1) - phase_at(u, i, j))/hy;
}

__cuda_callable__
Real Problem::grad_p_1_backward(const VectorView& u, Index i, Index j)
{
    return (phase_at(u, i, j) - phase_at(u, i-1, j))/hx;
}

__cuda_callable__
Real Problem::grad_p_2_backward(const VectorView& u, Index i, Index j)
{
   return (phase_at(u, i, j) - phase_at(u, i, j-1))/hy;
}

__cuda_callable__
Real Problem::div_T0(const VectorView& u, Index i, Index j)
{
    //if (i == 1 || i == sizeX-2 || j == 1 || j == sizeY-2)
        //forward
        return (T0_1(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_1(grad_p_1_forward(u, i-1, j), grad_p_2_forward(u, i-1, j)))/hx
               + (T0_2(grad_p_1_forward(u, i, j), grad_p_2_forward(u, i, j)) - T0_2(grad_p_1_forward(u, i, j-1), grad_p_2_forward(u, i, j-1)))/hy;
        //backward
        //return (T0_1(grad_p_1_backward(u, i+1, j), grad_p_2_backward(u, i+1, j)) - T0_1(grad_p_1_backward(u, i, j), grad_p_2_backward(u, i, j)))/hx
        //       + (T0_2(grad_p_1_backward(u, i, j+1), grad_p_2_backward(u, i, j+1)) - T0_2(grad_p_1_backward(u, i, j), grad_p_2_backward(u, i, j)))/hy;

    // both central, has to fulfill conditions on i and j
    //else
        //return (T0_1(grad_p_1_central(u, i+1, j), grad_p_2_central(u, i+1, j)) - T0_1(grad_p_1_central(u, i-1, j), grad_p_2_central(u, i-1, j)))/(2*hx)
        //       + (T0_2(grad_p_1_central(u, i, j+1), grad_p_2_central(u, i, j+1)) - T0_2(grad_p_1_central(u, i, j-1), grad_p_2_central(u, i, j-1)))/(2*hy);

    
    
}

__cuda_callable__
Real Problem::T0_1(const Real grad_p_1, const Real grad_p_2)
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    return psi(theta) * psi(theta) * grad_p_1 - psi(theta) * der_psi(theta) * grad_p_2;
}

__cuda_callable__
Real Problem::T0_2(const Real grad_p_1, const Real grad_p_2)
{
    Real theta = atan2(-grad_p_2, -grad_p_1);
    return psi(theta) * psi(theta) * grad_p_2 + psi(theta) * der_psi(theta) * grad_p_1;
}

__cuda_callable__
Real Problem::psi(const Real theta)
{
    return 1 + param.A*cos(param.m*(theta - param.theta_0));
}

__cuda_callable__
Real Problem::der_psi(const Real theta)
{
    return -param.A*param.m*sin(param.m*(theta - param.theta_0));
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
    return param.par_b*param.par_d/Real(30.0)
            * conc_at(u, i, j)
            * (1 - conc_at(u, i, j))
            * constants::M_Nb_beta(param.T)
		    * pow(constants::M_Nb_alpha(param.T)/constants::M_Nb_beta(param.T), polynom_p(u, i, j))
		    * sec_deriv_of_g_w_resp_to_c(u, i, j);
}

__cuda_callable__
Real Problem::get_phas_diff_coef(const VectorView& u, Index i, Index j)
{
   return param.par_b*param.par_d/Real(30.0)
          * conc_at(u, i, j)
		  * (1 - conc_at(u, i, j))
          * constants::M_Nb_beta(param.T)
		  * pow(constants::M_Nb_alpha(param.T)/constants::M_Nb_beta(param.T), polynom_p(u, i, j))
		  * deriv_of_g_w_resp_to_c_and_p(u, i, j);
}

__cuda_callable__
Real Problem::f_0(const VectorView& u, Index i, Index j)
{
   return param.par_a*phase_at(u, i, j)*(1 - phase_at(u, i, j))*(phase_at(u, i, j) - Real(1.0)/Real(2.0));
}

__cuda_callable__
Real Problem::F(const VectorView& u, Index i, Index j)
{
    if (param.force_term_type == FTType::Constant)
        return param.force_term_size;
    else if (param.force_term_type == FTType::Radial)
    {
        Real mid_x = (domain.x_right - domain.x_left)/2;
        Real mid_y = (domain.y_right - domain.y_left)/2;
        Real r = sqrt(pow(i*hx - mid_x, Real(2)) + pow(j*hy - mid_y, Real(2)));
        return -2/std::max(r, Real(0.1));
    }
    else if (param.force_term_type == FTType::Zirconium)
        return constants::G_m_alpha(conc_at(u, i, j), param.T) - constants::G_m_beta(conc_at(u, i, j), param.T);
    else
        return 0;

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
   return pow(phase_at(u, i, j), Real(2)) * pow(phase_at(u, i, j) - Real(1.0), Real(2));
}

__cuda_callable__
Real Problem::polynom_p(const VectorView& u, Index i, Index j)
{
   return 6*pow(phase_at(u, i, j), Real(5)) - 15*pow(phase_at(u, i, j), Real(4)) + 10*pow(phase_at(u, i, j), Real(3));
}

__cuda_callable__
Real Problem::der_polynom_p(const VectorView& u, Index i, Index j)
{
   return 30*(pow(phase_at(u, i, j), Real(4)) - 2*pow(phase_at(u, i, j), Real(3)) + pow(phase_at(u, i, j), Real(2)));
}

__cuda_callable__
Real Problem::sec_deriv_of_g_w_resp_to_c(const VectorView& u, Index i, Index j)
{
   Real c = conc_at(u, i, j);
   Real d2_G_alpha_wrt_c = constants::R * param.T / (c*(1-c))
   							 - 2 * constants::L_0_alpha;
   Real d2_G_beta_wrt_c = constants::R * param.T / (c*(1-c))
   							 - 2 * constants::L_0_beta(param.T)
							 + (6 - 12*c) * constants::L_0_i_beta(param.T);
   return polynom_p(u, i, j)*d2_G_alpha_wrt_c + (1 - polynom_p(u, i, j))*d2_G_beta_wrt_c;
}

__cuda_callable__
Real Problem::deriv_of_g_w_resp_to_c_and_p(const VectorView& u, Index i, Index j)
{
   Real c = conc_at(u, i, j);
   Real d_G_alpha_wrt_c = constants::G_Nb_alpha_0(param.T)
                            - constants::G_Zr_alpha_0(param.T)
                            + constants::R * param.T * (log(c) - log(1-c))
                            + (1 - 2*c) * constants::L_0_alpha;
   Real d_G_beta_wrt_c = constants::G_Nb_beta_0(param.T)
						   - constants::G_Zr_beta_0(param.T)
						   + constants::R * param.T * (log(c) - log(1-c))
						   + (1 - 2 * c) * constants::L_0_beta(param.T)
						   + (6*c - 6*pow(c, Real(2)) - 1) * constants::L_0_i_beta(param.T);

   return der_polynom_p(u, i, j)*(d_G_alpha_wrt_c - d_G_beta_wrt_c);
}