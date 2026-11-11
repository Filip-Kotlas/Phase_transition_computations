#include "ResultIO.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "npy.hpp"

namespace {

// Reads values.size() numbers of the type T from the stream and converts them to Real.
template< typename T >
bool read_values(std::istream& in, std::vector<Real>& values)
{
    std::vector<T> buffer(values.size());
    in.read(reinterpret_cast<char*>(buffer.data()),
            static_cast<std::streamsize>(buffer.size()) * static_cast<std::streamsize>(sizeof(T)));
    if (!in)
        return false;
    std::transform(buffer.begin(), buffer.end(), values.begin(), [](T value) { return static_cast<Real>(value); });
    return true;
}

}  // namespace

bool write_frame(const std::filesystem::path& folder, Index step, const VectorView& u, Index sizeX, Index sizeY)
{
    std::stringstream name;
    name << "Result-" << std::setw(5) << std::setfill('0') << step << ".npy";
    const std::filesystem::path file_path = folder / name.str();

    std::ofstream file(file_path, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!file)
    {
        std::cerr << "Unable to open the file " << name.str() << std::endl;
        return false;
    }

    // The whole solution is copied from the GPU at once. Reading it element by element
    // (u.getElement) would mean one small device-to-host transfer per element.
    TNL::Containers::Vector<Real, TNL::Devices::Host, Index> host_u(u.getSize());
    host_u = u;

    npy::write_header(file, npy::descriptor<Real>(),
                      { 2, static_cast<std::size_t>(sizeY), static_cast<std::size_t>(sizeX) });
    file.write(reinterpret_cast<const char*>(host_u.getData()),
               static_cast<std::streamsize>(host_u.getSize()) * static_cast<std::streamsize>(sizeof(Real)));

    file.close();
    if (!file)
    {
        std::cerr << "Error while writing the file " << name.str() << std::endl;
        return false;
    }
    return true;
}

bool read_frame(const std::filesystem::path& file_path, Index sizeX, Index sizeY,
                std::vector<Real>& phase, std::vector<Real>& conc)
{
    std::ifstream file(file_path, std::ios::in | std::ios::binary);
    if (!file)
    {
        std::cerr << "Unable to open the file: " << file_path << '\n';
        return false;
    }

    std::string descr;
    try {
        bool fortran_order = false;
        std::vector<std::size_t> shape;
        npy::read_header(file, descr, fortran_order, shape);

        const std::vector<std::size_t> expected_shape = { 2, static_cast<std::size_t>(sizeY), static_cast<std::size_t>(sizeX) };
        const bool known_type = descr == npy::descriptor<float>() || descr == npy::descriptor<double>();
        if (!known_type || fortran_order || shape != expected_shape) {
            std::cerr << "The .npy file has to be a C-ordered array of the type " << npy::descriptor<float>()
                      << " or " << npy::descriptor<double>() << " and the shape (2, " << sizeY << ", " << sizeX << ").\n";
            return false;
        }
    }
    catch (const std::runtime_error& e) {
        std::cerr << "Error while reading the .npy file: " << e.what() << '\n';
        return false;
    }

    const std::size_t N = static_cast<std::size_t>(sizeX) * static_cast<std::size_t>(sizeY);
    phase.resize(N);
    conc.resize(N);
    const bool complete = descr == npy::descriptor<float>()
                          ? read_values<float>(file, phase) && read_values<float>(file, conc)
                          : read_values<double>(file, phase) && read_values<double>(file, conc);
    if (!complete) {
        std::cerr << "The file contains less data than expected.\n";
        return false;
    }
    if (file.peek() != std::ifstream::traits_type::eof()) {
        std::cerr << "The file contains more data than expected.\n";
        return false;
    }
    return true;
}
