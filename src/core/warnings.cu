#include "core/warnings.hpp"

#include <iostream>

namespace {

std::vector<std::string>& collected()
{
    static std::vector<std::string> messages;
    return messages;
}

}  // namespace

void warn(const std::string& message)
{
    std::cerr << "Warning: " << message << std::endl;
    collected().push_back(message);
}

const std::vector<std::string>& warnings()
{
    return collected();
}
