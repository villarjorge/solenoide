// ISO C++20:
#include <algorithm>
#include <cmath>   // contiene std::sqrt
#include <fstream>
#include <numbers>
#include <string>
#include <vector>

// ISO C++23:
#include <cppexten/generator>
#include <cppexten/print>
#include <cppexten/ranges>

namespace stdx = std::experimental;

// biblioteca de terceros:
#include <boost/tokenizer.hpp> // tokenizador para extraer los valores en una línea del CSV
#include <dlib/optimization.h> // contiene el método LM de regresión no-lineal

#include <units/format.h>
#include <units/math.h>
#include <units/isq/si/length.h>
#include <units/isq/si/electric_field_strength.h>
#include <units/isq/si/magnetic_induction.h>
#include <units/isq/si/permeability.h>

namespace un = units;
namespace usi = units::isq::si;

using len_t = usi::length<usi::millimetre>; // Tipo de longitud en milimetros
using curr_t = usi::electric_current<usi::ampere>; // Tipo de corriente en Amperios
using mag_ind_t = usi::magnetic_induction<usi::millitesla>; // Tipo de intensidad del campo magnético en militeslas
using perm_t = usi::permeability<usi::henry_per_metre>; // Tipo de permeabilidad del vacío en Henry · metro^(-1)

struct Data {
    len_t position;
    mag_ind_t magnetic_induction;
};

auto main() -> int {
    stdx::println("Lorem ipsum, dolor sit amet");
}