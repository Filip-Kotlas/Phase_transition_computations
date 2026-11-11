#pragma once

// Reading and writing of the frames Result-NNNNN.npy in the calculations/ folder of a run (host code).
// A frame is a C-ordered NumPy array of the shape (2, sizeY, sizeX): the phase block followed by the
// concentration block, which is exactly the layout of the state vector u (index j*sizeX + i).
// The coordinates are not stored; the Python scripts compute them from the config.

#include <filesystem>
#include <vector>

#include "types.hpp"

// Writes u as folder/Result-NNNNN.npy, where NNNNN is step. The type of the values is Real.
// Prints the reason to std::cerr and returns false on failure.
bool write_frame(const std::filesystem::path& folder, Index step, const VectorView& u, Index sizeX, Index sizeY);

// Reads a frame of a sizeX x sizeY grid. The values may be stored as float or double; they are
// converted to Real. On success phase[k] and conc[k] hold the point k = j*sizeX + i.
// Prints the reason to std::cerr and returns false on failure.
bool read_frame(const std::filesystem::path& file_path, Index sizeX, Index sizeY,
                std::vector<Real>& phase, std::vector<Real>& conc);
