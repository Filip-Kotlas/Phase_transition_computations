#pragma once

// Warnings of the solver (host code). warn() prints the warning to std::cerr, as before, and also
// keeps it, so that the run records its warnings in parameters.txt (RunFolder::append_warnings),
// also those that come before the run folder exists (e.g. while the config is read).

#include <string>
#include <vector>

// Prints "Warning: <message>" to std::cerr and keeps the message.
void warn(const std::string& message);

// The messages of all warn() calls so far, in order.
const std::vector<std::string>& warnings();
