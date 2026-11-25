#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "io/TdbDatabase.hpp"

// A node of a parsed expression.
struct TdbDatabase::Node {
    enum class Kind { Number, T, P, Function, Negate, Add, Subtract, Multiply, Divide, Power, Ln, Exp };
    Kind kind;
    double number = 0;          // Number
    std::string name;           // Function: the name of the FUNCTION
    std::shared_ptr<const Node> a, b;   // the operands
};

namespace {

using Node = TdbDatabase::Node;
using NodePtr = std::shared_ptr<const Node>;
using Kind = Node::Kind;

// A ,, limit of a range: the database's TEMPERATURE_LIMITS, resolved when the range is evaluated
// (the statement may come after the FUNCTIONs).
constexpr double default_limit = std::numeric_limits<double>::quiet_NaN();

std::string upper(std::string text)
{
    for (char& ch : text)
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    return text;
}

bool is_space(char ch) { return std::isspace(static_cast<unsigned char>(ch)) != 0; }

std::string trim(const std::string& text)
{
    size_t begin = 0, end = text.size();
    while (begin < end && is_space(text[begin])) ++begin;
    while (end > begin && is_space(text[end - 1])) --end;
    return text.substr(begin, end - begin);
}

std::string without_whitespace(const std::string& text)
{
    std::string result;
    for (char ch : text)
        if (!is_space(ch))
            result += ch;
    return result;
}

std::vector<std::string> split_whitespace(const std::string& text)
{
    std::istringstream stream(text);
    std::vector<std::string> words;
    for (std::string word; stream >> word; )
        words.push_back(word);
    return words;
}

// Splits at sep, keeping empty parts.
std::vector<std::string> split(const std::string& text, char sep)
{
    std::vector<std::string> parts(1);
    for (char ch : text) {
        if (ch == sep)
            parts.emplace_back();
        else
            parts.back() += ch;
    }
    return parts;
}

// A constituent or function name without the markers % (major constituent) and # (Thermo-Calc's mark of
// a function name in exported files).
std::string bare_name(std::string name)
{
    name.erase(std::remove_if(name.begin(), name.end(), [](char ch) { return ch == '%' || ch == '#'; }), name.end());
    return name;
}

// A phase name without the suffix :L, :G, ... (LIQUID:L).
std::string phase_name(const std::string& name)
{
    return name.substr(0, name.find(':'));
}

// Whether word abbreviates keyword as Thermo-Calc allows: the same number of parts separated by _ or -,
// each a prefix of the keyword's part (TEMP-LIM for TEMPERATURE_LIMITS, PAR for PARAMETER).
bool abbreviates(std::string word, const std::string& keyword)
{
    std::replace(word.begin(), word.end(), '-', '_');
    const std::vector<std::string> word_parts = split(word, '_');
    const std::vector<std::string> keyword_parts = split(keyword, '_');
    if (word_parts.size() != keyword_parts.size())
        return false;
    for (size_t n = 0; n < word_parts.size(); ++n) {
        if (word_parts[n].empty() || keyword_parts[n].compare(0, word_parts[n].size(), word_parts[n]) != 0)
            return false;
    }
    return true;
}

NodePtr make(Kind kind, NodePtr a = nullptr, NodePtr b = nullptr)
{
    auto node = std::make_shared<Node>();
    node->kind = kind;
    node->a = std::move(a);
    node->b = std::move(b);
    return node;
}

// A recursive descent parser of an expression without whitespace:
//   sum     = product {(+|-) product}
//   product = unary {(*|/) unary}
//   unary   = (+|-) unary | power
//   power   = primary [** unary]          (right associative, T**-1 is allowed)
//   primary = number | T | P | NAME | (LN|LOG|EXP) ( sum ) | ( sum )
class ExpressionParser {
public:
    explicit ExpressionParser(const std::string& text) : text(text) {}

    NodePtr parse()
    {
        NodePtr node = sum();
        if (pos != text.size())
            fail("unexpected '" + text.substr(pos) + "'");
        return node;
    }

private:
    bool at(char ch) const { return pos < text.size() && text[pos] == ch; }
    bool at_power() const { return text.compare(pos, 2, "**") == 0; }

    [[noreturn]] void fail(const std::string& message) const
    {
        throw std::runtime_error("cannot read the expression " + text + ": " + message);
    }

    NodePtr sum()
    {
        NodePtr left = product();
        while (at('+') || at('-')) {
            const Kind kind = text[pos++] == '+' ? Kind::Add : Kind::Subtract;
            left = make(kind, left, product());
        }
        return left;
    }

    NodePtr product()
    {
        NodePtr left = unary();
        while (at('/') || (at('*') && !at_power())) {
            const Kind kind = text[pos++] == '*' ? Kind::Multiply : Kind::Divide;
            left = make(kind, left, unary());
        }
        return left;
    }

    NodePtr unary()
    {
        if (at('+')) {
            ++pos;
            return unary();
        }
        if (at('-')) {
            ++pos;
            return make(Kind::Negate, unary());
        }
        return power();
    }

    NodePtr power()
    {
        NodePtr base = primary();
        if (at_power()) {
            pos += 2;
            return make(Kind::Power, base, unary());
        }
        return base;
    }

    NodePtr primary()
    {
        if (pos >= text.size())
            fail("it ends too early");
        const char ch = text[pos];
        if (ch == '(') {
            ++pos;
            NodePtr inner = sum();
            if (!at(')'))
                fail("missing )");
            ++pos;
            return inner;
        }
        if (std::isdigit(static_cast<unsigned char>(ch)) || ch == '.') {
            const char* begin = text.c_str() + pos;
            char* end = nullptr;
            const double value = std::strtod(begin, &end);
            if (end == begin)
                fail("bad number at '" + text.substr(pos) + "'");
            pos += static_cast<size_t>(end - begin);
            NodePtr node = make(Kind::Number);
            std::const_pointer_cast<Node>(node)->number = value;
            return node;
        }
        if (std::isalpha(static_cast<unsigned char>(ch)) || ch == '_') {
            const size_t begin = pos;
            while (pos < text.size()
                   && (std::isalnum(static_cast<unsigned char>(text[pos])) || text[pos] == '_' || text[pos] == '#'))
                ++pos;
            const std::string name = bare_name(text.substr(begin, pos - begin));
            if (at('(')) {
                Kind kind;
                if (name == "LN" || name == "LOG")
                    kind = Kind::Ln;
                else if (name == "EXP")
                    kind = Kind::Exp;
                else
                    fail("unknown function " + name + "()");
                ++pos;
                NodePtr argument = sum();
                if (!at(')'))
                    fail("missing ) after " + name + "(");
                ++pos;
                return make(kind, argument);
            }
            if (name == "T")
                return make(Kind::T);
            if (name == "P")
                return make(Kind::P);
            NodePtr node = make(Kind::Function);
            std::const_pointer_cast<Node>(node)->name = name;
            return node;
        }
        fail("unexpected '" + text.substr(pos) + "'");
    }

    const std::string text;
    size_t pos = 0;
};

}  // namespace

TdbDatabase TdbDatabase::load(const std::filesystem::path& file)
{
    TdbDatabase database;
    database.path = file;
    std::ifstream stream(file);
    if (!stream)
        throw std::runtime_error("Unable to open the thermodynamic database " + file.string() + ".");

    // The text without comments ($ to the end of the line), split into the statements ending with !.
    std::string text;
    for (std::string line; std::getline(stream, line); )
        text += line.substr(0, line.find('$')) + '\n';
    std::vector<std::string> statements = split(text, '!');
    if (!trim(statements.back()).empty())
        database.fail("the text at its end is not closed by !: " + trim(statements.back()));
    statements.pop_back();
    for (const std::string& statement : statements) {
        if (!trim(statement).empty())
            database.parse_statement(upper(trim(statement)));
    }
    return database;
}

void TdbDatabase::fail(const std::string& message) const
{
    throw std::runtime_error("The thermodynamic database " + path.string() + ": " + message);
}

void TdbDatabase::parse_statement(const std::string& statement)
{
    std::istringstream stream(statement);
    std::string keyword;
    stream >> keyword;
    const std::string rest = trim(statement.substr(keyword.size()));
    const std::vector<std::string> words = split_whitespace(rest);
    const auto need_words = [&](size_t count) {
        if (words.size() < count)
            fail("incomplete statement: " + statement);
    };

    if (abbreviates(keyword, "ELEMENT")) {
        need_words(1);
        elements.insert(words[0]);
    }
    else if (abbreviates(keyword, "PHASE")) {
        need_words(3);
        const std::string name = phase_name(words[0]);
        if (phases.count(name))
            fail("the phase " + name + " is defined twice");
        Phase phase;
        phase.type_codes = words[1];
        const int sublattices = std::atoi(words[2].c_str());
        if (sublattices < 1 || words[2] != std::to_string(sublattices))
            fail("bad number of sublattices in: " + statement);
        need_words(3 + static_cast<size_t>(sublattices));
        for (int s = 0; s < sublattices; ++s) {
            char* end = nullptr;
            const double sites = std::strtod(words[3 + s].c_str(), &end);
            if (*end != '\0' || !(sites > 0))
                fail("bad number of sites in: " + statement);
            phase.sites.push_back(sites);
        }
        phase.constituents.resize(static_cast<size_t>(sublattices));
        phases[name] = phase;
    }
    else if (abbreviates(keyword, "CONSTITUENT")) {
        // CONSTITUENT NAME : A,B : C : (the name is followed by whitespace, it may contain : itself)
        need_words(1);
        const std::string name = phase_name(words[0]);
        auto phase = phases.find(name);
        if (phase == phases.end())
            fail("CONSTITUENT of the undefined phase " + name);
        std::vector<std::string> parts = split(trim(rest.substr(words[0].size())), ':');
        if (parts.size() < 3 || !trim(parts.front()).empty() || !trim(parts.back()).empty())
            fail("bad CONSTITUENT statement: " + statement);
        parts.erase(parts.begin());
        parts.pop_back();
        if (parts.size() != phase->second.sites.size())
            fail("the CONSTITUENT statement of " + name + " does not have " + std::to_string(phase->second.sites.size())
                 + " sublattices: " + statement);
        for (size_t s = 0; s < parts.size(); ++s) {
            std::string list = parts[s];
            std::replace(list.begin(), list.end(), ',', ' ');
            for (const std::string& constituent : split_whitespace(list))
                phase->second.constituents[s].push_back(bare_name(constituent));
            if (phase->second.constituents[s].empty())
                fail("an empty sublattice in: " + statement);
        }
    }
    else if (abbreviates(keyword, "TYPE_DEFINITION")) {
        need_words(1);
        type_definitions[words[0]] = statement;
    }
    else if (abbreviates(keyword, "FUNCTION")) {
        need_words(2);
        const std::string name = bare_name(words[0]);
        if (functions.count(name))
            fail("the function " + name + " is defined twice");
        functions[name] = parse_piecewise(rest.substr(words[0].size()), "FUNCTION " + name);
    }
    else if (abbreviates(keyword, "PARAMETER")) {
        // PARAMETER TYPE(PHASE,A,B:C;ORDER) LOW EXPRESSION; HIGH N REFERENCE
        const size_t open = rest.find('('), close = rest.find(')');
        if (open == std::string::npos || close == std::string::npos || close < open)
            fail("bad PARAMETER statement: " + statement);
        Parameter parameter;
        parameter.text = without_whitespace(rest.substr(0, close + 1));
        parameter.type = trim(rest.substr(0, open));
        std::string inside = without_whitespace(rest.substr(open + 1, close - open - 1));
        parameter.order = 0;
        const size_t semicolon = inside.find(';');
        if (semicolon != std::string::npos) {
            const std::string order = inside.substr(semicolon + 1);
            parameter.order = std::atoi(order.c_str());
            if (order != std::to_string(parameter.order))
                fail("bad order in " + parameter.text);
            inside = inside.substr(0, semicolon);
        }
        const size_t comma = inside.find(',');
        if (comma == std::string::npos)
            fail("no constituents in " + parameter.text);
        parameter.phase = phase_name(inside.substr(0, comma));
        for (const std::string& sublattice : split(inside.substr(comma + 1), ':')) {
            std::vector<std::string> constituents;
            for (const std::string& constituent : split(sublattice, ','))
                constituents.push_back(bare_name(constituent));
            parameter.constituents.push_back(constituents);
        }
        parameter.value = parse_piecewise(rest.substr(close + 1), parameter.text);
        parameters.push_back(parameter);
    }
    else if (abbreviates(keyword, "TEMPERATURE_LIMITS")) {
        need_words(2);
        T_low = std::strtod(words[0].c_str(), nullptr);
        T_high = std::strtod(words[1].c_str(), nullptr);
        if (!(T_low > 0 && T_high > T_low))
            fail("bad TEMPERATURE_LIMITS: " + statement);
    }
    else if (abbreviates(keyword, "SPECIES") || abbreviates(keyword, "DEFINE_SYSTEM_DEFAULT")
             || abbreviates(keyword, "DEFAULT_COMMAND") || abbreviates(keyword, "ASSESSED_SYSTEMS")
             || abbreviates(keyword, "LIST_OF_REFERENCES") || abbreviates(keyword, "ADD_REFERENCES")
             || abbreviates(keyword, "DATABASE_INFORMATION") || abbreviates(keyword, "VERSION_DATE")
             || abbreviates(keyword, "VERSION_DATA") || abbreviates(keyword, "REFERENCE_FILE")) {
        // bookkeeping, not needed for the Gibbs energies
    }
    else
        fail("unknown statement " + keyword + " (only a subset of the TDB format is read, see TdbDatabase.hpp)");
}

TdbDatabase::Piecewise TdbDatabase::parse_piecewise(const std::string& text, const std::string& what) const
{
    // LOW EXPRESSION; HIGH Y EXPRESSION; HIGH Y ... EXPRESSION; HIGH N [REFERENCE]
    size_t pos = 0;
    const auto skip_space = [&]() { while (pos < text.size() && is_space(text[pos])) ++pos; };
    const auto limit = [&]() {
        skip_space();
        if (text.compare(pos, 2, ",,") == 0) {
            pos += 2;
            return default_limit;
        }
        const char* begin = text.c_str() + pos;
        char* end = nullptr;
        const double value = std::strtod(begin, &end);
        if (end == begin)
            fail("no temperature limit at '" + trim(text.substr(pos)) + "' in " + what);
        pos += static_cast<size_t>(end - begin);
        return value;
    };

    Piecewise ranges;
    double low = limit();
    while (true) {
        const size_t semicolon = text.find(';', pos);
        if (semicolon == std::string::npos)
            fail("no ; after the expression in " + what);
        const std::string expression = without_whitespace(text.substr(pos, semicolon - pos));
        pos = semicolon + 1;
        const double high = limit();
        try {
            ranges.push_back({low, high, ExpressionParser(expression).parse()});
        }
        catch (const std::runtime_error& e) {
            fail(std::string(e.what()) + " (in " + what + ")");
        }
        skip_space();
        const char flag = pos < text.size() ? text[pos++] : '\0';
        if (flag == 'N')
            return ranges;
        if (flag != 'Y')
            fail("expected Y or N after the upper limit in " + what);
        low = high;
    }
}

double TdbDatabase::evaluate(const Piecewise& piecewise, double T, const std::string& what, int depth) const
{
    const auto resolve = [](double limit, double fallback) { return std::isnan(limit) ? fallback : limit; };
    for (size_t n = 0; n < piecewise.size(); ++n) {
        const double low = resolve(piecewise[n].low, T_low);
        const double high = resolve(piecewise[n].high, T_high);
        const bool last = n + 1 == piecewise.size();
        if (T >= low && (T < high || (last && T <= high)))
            return evaluate(*piecewise[n].expression, T, depth);
    }
    fail("T = " + std::to_string(T) + " K is outside the range ["
         + std::to_string(resolve(piecewise.front().low, T_low)) + ", "
         + std::to_string(resolve(piecewise.back().high, T_high)) + "] K of " + what);
}

double TdbDatabase::evaluate(const Node& node, double T, int depth) const
{
    switch (node.kind) {
        case Kind::Number:   return node.number;
        case Kind::T:        return T;
        case Kind::P:        return P;
        case Kind::Negate:   return -evaluate(*node.a, T, depth);
        case Kind::Add:      return evaluate(*node.a, T, depth) + evaluate(*node.b, T, depth);
        case Kind::Subtract: return evaluate(*node.a, T, depth) - evaluate(*node.b, T, depth);
        case Kind::Multiply: return evaluate(*node.a, T, depth) * evaluate(*node.b, T, depth);
        case Kind::Divide:   return evaluate(*node.a, T, depth) / evaluate(*node.b, T, depth);
        case Kind::Power:    return std::pow(evaluate(*node.a, T, depth), evaluate(*node.b, T, depth));
        case Kind::Exp:      return std::exp(evaluate(*node.a, T, depth));
        case Kind::Ln: {
            const double argument = evaluate(*node.a, T, depth);
            if (!(argument > 0))
                fail("LN of " + std::to_string(argument) + " at T = " + std::to_string(T) + " K");
            return std::log(argument);
        }
        case Kind::Function: {
            const auto function = functions.find(node.name);
            if (function == functions.end())
                fail("the function " + node.name + " is not defined");
            if (depth > 100)
                fail("the functions refer to each other in a cycle (at " + node.name + ")");
            return evaluate(function->second, T, "FUNCTION " + node.name, depth + 1);
        }
    }
    fail("an unknown expression node");
}

double TdbDatabase::function_value(const std::string& name, double T) const
{
    const auto function = functions.find(upper(name));
    if (function == functions.end())
        fail("the function " + upper(name) + " is not defined");
    return evaluate(function->second, T, "FUNCTION " + upper(name), 0);
}

double TdbDatabase::reference_energy(const std::string& element, double T) const
{
    // The SGTE name of the Gibbs energy of the element in its reference state (SER): GHSERZR, GHSERNB,
    // GHSERHH (one-letter symbols doubled).
    const std::string name = "GHSER" + (element.size() == 1 ? element + element : element);
    if (!functions.count(name))
        fail("there is no function " + name + ", the reference state of " + element
             + " (G is given relative to the references of the elements)");
    return function_value(name, T);
}

TdbDatabase::PhaseValues TdbDatabase::phase_values(const std::string& phase_requested, const std::string& solute_requested,
                                                   const std::string& host_requested, double T,
                                                   SoluteSite site, int terms) const
{
    const std::string name = upper(phase_requested);
    const std::string solute = upper(solute_requested);
    const std::string host = upper(host_requested);
    const auto found = phases.find(name);
    if (found == phases.end()) {
        std::string known;
        for (const auto& item : phases)
            known += (known.empty() ? "" : ", ") + item.first;
        fail("there is no phase " + name + " (it has " + known + ")");
    }
    const Phase& phase = found->second;
    for (const std::string& element : {solute, host}) {
        if (element == "VA" || !elements.count(element))
            fail(element + " is not an element of the database");
    }

    // The type definitions: only SEQ and MAGNETIC (whose parameters TC, BMAGN are refused below).
    for (const char code : phase.type_codes) {
        if (code == '%')
            continue;
        const auto definition = type_definitions.find(std::string(1, code));
        if (definition == type_definitions.end())
            fail("the type code " + std::string(1, code) + " of the phase " + name + " has no TYPE_DEFINITION");
        const std::vector<std::string> words = split_whitespace(definition->second);
        const bool known = std::find(words.begin(), words.end(), "SEQ") != words.end()
                           || std::find(words.begin(), words.end(), "MAGNETIC") != words.end();
        if (!known)
            fail("the phase " + name + " uses an unsupported type definition: " + definition->second);
    }

    // The sublattices: exactly one where two constituents mix.
    const size_t sublattices = phase.sites.size();
    int mixing = -1;
    for (size_t s = 0; s < sublattices; ++s) {
        if (phase.constituents[s].empty())
            fail("the phase " + name + " has no CONSTITUENT statement");
        if (phase.constituents[s].size() > 1) {
            if (mixing != -1 || phase.constituents[s].size() != 2)
                fail("the phase " + name + " is not a binary phase with mixing on one sublattice");
            mixing = static_cast<int>(s);
        }
    }
    if (mixing == -1)
        fail("the phase " + name + " has no sublattice where two constituents mix (a line compound)");
    const std::vector<std::string>& mixed = phase.constituents[mixing];
    if (std::find(mixed.begin(), mixed.end(), solute) == mixed.end())
        fail("the solute " + solute + " does not mix in the phase " + name);
    const std::string other = mixed[0] == solute ? mixed[1] : mixed[0];

    // The host atoms per formula unit and y = k*u.
    SoluteSite found_site;
    double host_atoms, k;
    if (other == host) {
        found_site = SoluteSite::substitutional;
        // substitutional: (solute, host)_a (VA)_c ...
        for (size_t s = 0; s < sublattices; ++s) {
            if (static_cast<int>(s) != mixing && phase.constituents[s] != std::vector<std::string>{"VA"})
                fail("the phase " + name + " has other constituents than VA outside the sublattice of "
                     + solute + " and " + host);
        }
        host_atoms = phase.sites[mixing];
        k = 1;
    }
    else if (other == "VA") {
        // interstitial: (host)_a (solute, VA)_c
        found_site = SoluteSite::interstitial;
        const size_t h = 1 - static_cast<size_t>(mixing);
        if (sublattices != 2 || phase.constituents[h] != std::vector<std::string>{host})
            fail("the phase " + name + " is not (" + host + ")(" + solute + ",VA)");
        host_atoms = phase.sites[h];
        k = phase.sites[h] / phase.sites[mixing];
    }
    else
        fail("in the phase " + name + " the solute " + solute + " mixes with " + other + ", not with "
             + host + " or VA");
    const auto kind = [](SoluteSite s) { return s == SoluteSite::substitutional ? "substitutional" : "interstitial"; };
    if (found_site != site)
        fail("the phase " + name + " has " + solute + " " + kind(found_site) + ", the model needs it " + kind(site));
    if (terms < 1)
        fail("a phase needs at least one Redlich-Kister term");

    // The parameters of the phase at T.
    bool have_solute = false, have_other = false;
    double G_solute = 0, G_other = 0;
    std::vector<double> L(static_cast<size_t>(terms), 0.0);
    std::vector<bool> have_L(static_cast<size_t>(terms), false);
    for (const Parameter& parameter : parameters) {
        if (parameter.phase != name)
            continue;
        if (parameter.type != "G" && parameter.type != "L")
            fail("the parameter " + parameter.text + " is of a type that is not supported (only G and L)");
        if (parameter.constituents.size() != sublattices)
            fail("the parameter " + parameter.text + " does not have " + std::to_string(sublattices) + " sublattices");
        for (size_t s = 0; s < sublattices; ++s) {
            if (static_cast<int>(s) != mixing && parameter.constituents[s] != phase.constituents[s])
                fail("the parameter " + parameter.text + " has constituents the phase does not have");
        }
        const std::vector<std::string>& at_mixing = parameter.constituents[mixing];
        const double value = evaluate(parameter.value, T, parameter.text, 0);
        if (parameter.type == "G") {
            if (at_mixing.size() != 1 || parameter.order != 0)
                fail("the parameter " + parameter.text + " is not an end member");
            bool& have = at_mixing[0] == solute ? have_solute : have_other;
            if (at_mixing[0] != solute && at_mixing[0] != other)
                fail("the parameter " + parameter.text + " has constituents the phase does not have");
            if (have)
                fail("the parameter " + parameter.text + " is given twice");
            have = true;
            (at_mixing[0] == solute ? G_solute : G_other) = value;
        }
        else {
            // L(phase, A,B; n) multiplies (y_A - y_B)^n: written as (other, solute) it changes sign for odd n.
            double sign;
            if (at_mixing == std::vector<std::string>{solute, other})
                sign = 1;
            else if (at_mixing == std::vector<std::string>{other, solute})
                sign = parameter.order % 2 == 0 ? 1 : -1;
            else
                fail("the parameter " + parameter.text + " is not an interaction of " + solute + " and " + other);
            if (parameter.order < 0 || parameter.order >= terms)
                fail("the parameter " + parameter.text + " has the order " + std::to_string(parameter.order)
                     + ", but the model reads the Redlich-Kister terms only up to the order " + std::to_string(terms - 1)
                     + " (raise the number of terms of the phase in the model)");
            if (have_L[parameter.order])
                fail("the parameter " + parameter.text + " is given twice");
            have_L[parameter.order] = true;
            L[parameter.order] = sign*value;
        }
    }
    if (!have_solute || !have_other)
        fail("the phase " + name + " lacks the end member with " + (have_solute ? other : solute)
             + " on its mixing sublattice");

    // The references of the elements, per formula unit of each end member.
    const auto reference = [&](const std::string& at_mixing) {
        double sum = 0;
        for (size_t s = 0; s < sublattices; ++s) {
            const std::string& constituent = static_cast<int>(s) == mixing ? at_mixing : phase.constituents[s][0];
            if (constituent != "VA")
                sum += phase.sites[s]*reference_energy(constituent, T);
        }
        return sum;
    };

    PhaseValues result;
    result.k = k;
    result.G_solute = (G_solute - reference(solute)) / host_atoms;
    result.G_other = (G_other - reference(other)) / host_atoms;
    result.mixing_RT = R*T*phase.sites[mixing] / host_atoms;
    result.L = L;
    for (double& term : result.L)
        term /= host_atoms;
    return result;
}
