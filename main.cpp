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
#include <units/generic/dimensionless.h>
//#include <units/isq/si/electric_current.h>
#include <units/isq/si/electric_field_strength.h> // Por que ampere está aquí y no en electric_current?
#include <units/isq/si/magnetic_induction.h>
#include <units/isq/si/permeability.h>

namespace un = units;
namespace usi = units::isq::si;

using len_t = usi::length<usi::millimetre>; // Tipo de longitud en milimetros
using curr_t = usi::electric_current<usi::ampere>; // Tipo de corriente en Amperios
using mag_ind_t = usi::magnetic_induction<usi::millitesla>; // Tipo de intensidad del campo magnético en militeslas
using perm_t = usi::permeability<usi::henry_per_metre>; // Tipo de permeabilidad del vacío en Henry · metro^(-1)
using dimless_t = un::quantity<un::dim_one, un::one, float>;

struct Data {
    len_t position;
    mag_ind_t magnetic_induction;
};

auto get_data(std::string path) -> stdx::generator<Data> {
    auto ifs = std::ifstream{path, std::ios::binary};
    if (!ifs) throw std::ios::failure{"unable to open the CSV file"};
    
    auto ln = std::string{}; // string auxiliar para almacenar cada línea del CSV
    while (std::getline(ifs, ln)) { // procesamos cada línea del CSV
        auto values = boost::tokenizer{ln, boost::char_separator{","}}
                  | std::views::transform([](std::string v){ return std::stod(v); })
                  | stdx::ranges::to<std::vector<double>>();

        co_yield Data{
            .position = len_t{values[0]},
            .magnetic_induction = mag_ind_t{values[1]}
        };
   }
}

auto theoretic_magnetic_induction(len_t z, len_t z_0, perm_t mu_0) -> mag_ind_t {
    // Inducción magnética en un solenoide, z + z_0 es la distancia al centro del solenoide
    // Por ahora estas constantes no se corresponden a la realidad
    auto const N = dimless_t{100.0}; // número de espiras, adimensional
    auto const I = curr_t{1.0}; // Corriente que pasa por el solenoide en Amperios
    auto const L = len_t{70.0}; // Longitud del solenoide en milímetros
    auto const D = len_t{10.0}; // Diámetro del solenoide
    auto const R = D/2;

    auto const a = z - z_0 + L/2;
    auto const b = z - z_0 - L/2;

    auto const sqrt_one = un::sqrt(un::pow<2>(a) + un::pow<2>(R));
    auto const sqrt_two = un::sqrt(un::pow<2>(b) + un::pow<2>(R));

    return (mu_0/2)*(N*I/L)*(a/sqrt_one - b/sqrt_two);
}

auto main() -> int {
    auto data_vec = get_data("../../magnetic_data_no_z_correction.csv") 
                    | stdx::ranges::to<std::vector<Data>>();

    for (Data d : data_vec) {
        stdx::println("position: {}, field: {}", d.position, d.magnetic_induction);
    }
}