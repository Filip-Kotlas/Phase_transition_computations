#include <fstream>

template< typename Vector, typename Real = typename Vector::RealType, typename Index = typename Vector::IndexType >
void
write( std::fstream& file, const Vector& u, const Index n, const Real& h, const Real& time )
{
   file << "# time = " << time << std::endl;
   for( Index i = 0; i < n; i++ )
      file << i * h << " " << u.getElement( i ) << std::endl;
   file << std::endl;
}

template< typename Vector, typename Real = typename Vector::RealType >
void
write( std::fstream& file, const Vector& u, const Real& h, const Real& time )
{
   file << "# time = " << time << std::endl;
   const auto localRange = u.getLocalRange();
   for( auto i = localRange.getBegin(); i < localRange.getEnd(); i++ )
      file << i * h << " " << u.getElement( i ) << std::endl;
   file << std::endl;
}