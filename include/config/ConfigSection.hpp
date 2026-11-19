#pragma once

// Strict reading of one section of the config (a JSON object), used for every section the solver reads
// (solver, solver.domain, init_cond, bound_cond and the model's model.<type>): every key has to be
// present (there are no defaults), so a typo ("ksl" instead of "ksi") is an error ("missing ksi"),
// and its value has to have the right type. check_all_read() only WARNS about the keys nobody read
// (a stray key, a parameter of another model, or a key that a newer version of the solver added):
// an error there would make every commit refuse the configs of later commits, so `cluster.sh test`
// and `bench` against an older commit would fail whenever a key is added. Host code only.

#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include "external/json.hpp"
#include "core/warnings.hpp"

class ConfigSection {
public:
    // data is the JSON object at the path name in the config (e.g. "model.allen_cahn"); the name
    // is used in the error messages.
    ConfigSection(const nlohmann::json& data, std::string name)
    :  data(data), name(std::move(name))
    {
        if (!this->data.is_object())
            throw std::runtime_error("The config entry " + this->name + " has to be a section {...}, it is "
                                     + this->data.dump());
    }

    bool contains(const std::string& key) const { return data.contains(key); }

    // The value of the key, converted to T. Throws when the key is missing or the value has the wrong
    // type (a string for a number, or a number with a decimal point for an integer: nlohmann would
    // silently truncate 4.5 to 4).
    template< class T >
    T get(const std::string& key)
    {
        const nlohmann::json& value = take(key);
        if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) {
            if (!value.is_number_integer())
                throw std::runtime_error("The config value " + path(key) + " has to be an integer, it is " + value.dump() + ".");
        }
        try {
            return value.get<T>();
        }
        catch (const nlohmann::json::exception&) {
            throw std::runtime_error("The config value " + path(key) + " has the wrong type: " + value.dump() + ".");
        }
    }

    // Warns (warn(): on std::cerr and in parameters.txt) about every key of the section that get()
    // did not read. Call it after reading the section.
    void check_all_read() const
    {
        for (const auto& item : data.items()) {
            if (read_keys.count(item.key()) == 0)
                warn("the config entry " + path(item.key()) + " is not used (a stray key,"
                     " a parameter of another model, or a key of a newer version of the solver?).");
        }
    }

private:
    std::string path(const std::string& key) const { return name + "." + key; }

    // The value of the key, marked as read; throws when it is missing.
    const nlohmann::json& take(const std::string& key)
    {
        if (!data.contains(key))
            throw std::runtime_error("The config is missing " + path(key) + ".");
        read_keys.insert(key);
        return data.at(key);
    }

    nlohmann::json data;
    std::string name;
    std::set<std::string> read_keys;
};
