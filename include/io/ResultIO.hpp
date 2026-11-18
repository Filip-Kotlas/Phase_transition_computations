#pragma once

// Reading and writing of the frames Result-NNNNN.npy in the calculations/ folder of a run (host code).
// A frame is a C-ordered NumPy array of the shape (2, sizeY, sizeX): the phase block followed by the
// concentration block, which is exactly the layout of the state vector u (index j*sizeX + i).
// The coordinates are not stored; the Python scripts compute them from the config.

#include <filesystem>
#include <vector>

#include "core/types.hpp"
#include "core/Grid.hpp"

// Writes u as folder/Result-NNNNN.npy, where NNNNN is step. The type of the values is Real.
// Prints the reason to std::cerr and returns false on failure.
bool write_frame(const std::filesystem::path& folder, Index step, const VectorView& u, const Grid& grid);

// Writes u in the frame format to file_path (any name; write_frame uses it for the frames, the
// time loop for info/stuck_state.npy). Prints the reason to std::cerr and returns false on failure.
bool write_state(const std::filesystem::path& file_path, const VectorView& u, const Grid& grid);

// Reads a frame of the given grid (its shape has to be (2, grid.sizeY, grid.sizeX)). The values may be
// stored as float or double; they are converted to Real. On success phase[k] and conc[k] hold the point
// k = grid.index(i, j).
// Prints the reason to std::cerr and returns false on failure.
bool read_frame(const std::filesystem::path& file_path, const Grid& grid,
                std::vector<Real>& phase, std::vector<Real>& conc);
