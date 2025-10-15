#pragma once

#include <TNL/Containers/Vector.h>
#include <TNL/Solvers/ODE/ODESolver.h>
#include <TNL/Solvers/ODE/Methods/OriginalRungeKutta.h>
#include <TNL/Solvers/IterativeSolverMonitor.h>

using Real = double;
using Index = int;
using Device = TNL::Devices::Cuda;

using Vector = TNL::Containers::Vector< Real, Device, Index >;
using VectorView = typename Vector::ViewType;
using Method = TNL::Solvers::ODE::Methods::OriginalRungeKutta< Real >;
using ODESolver = TNL::Solvers::ODE::ODESolver< Method, Vector >;

