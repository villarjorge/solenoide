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
using dimless_t = un::quantity<un::dim_one, un::one, float>; // Tipo sin unidades

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
    auto const N = dimless_t{200.0}; // número de espiras, adimensional
    auto const I = curr_t{.9}; // Corriente que pasa por el solenoide en Amperios
    auto const L = len_t{105.0}; // Longitud del solenoide en milímetros
    auto const D = len_t{41.0}; // Diámetro del solenoide
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

    // dlib ---------------------------------------------------------------------

    // tenemos dos parametros, por lo que necesitamos una matriz de 2 por 1
    using parameter_t = dlib::matrix<double, 2, 1>; // Tipo de parámetro a optimizar (primer elemento es z_0 y segundo es mu_0)

    // Función de cálculo de residuos (diferencia entre el campo teórico y el real)
    auto residual = [](Data const& d, parameter_t const& param) -> double {
        // Desglosa los datos 
        auto const [pos, induc] = d;
        // Desglosa el parámetro y envuelvelo en unidades
        auto const z_0_param = len_t{param(0)};
        auto const mu_0_param = perm_t{param(1)};

        return (theoretic_magnetic_induction(pos, z_0_param, mu_0_param) - induc).number();
    };
    // Estimación inicial de los parámetros z_0 y mu_0
    auto parameters_array_without_units = parameter_t{5.0, 1e-6};

    // Ejecutamos el método LM con los datos experimentales, optimizando el valor de z_0 y mu_0. 
    // El resultado obtenido multiplicado por 2 proporciona el RSS (suma de los residuos al cuadrado) medido implícitamente en Teslas^2:
    // Ejemplo de uso de esta función: http://dlib.net/least_squares_ex.cpp.html

    auto rss_without_units = 2.0*dlib::solve_least_squares_lm(
        dlib::objective_delta_stop_strategy{1.0e-7}, // precisión
        residual,
        dlib::derivative(residual), // aproximación numérica de la derivada del residuo 
        data_vec, // todos los primeros argumentos que ir mandando a la función residual
        parameters_array_without_units // objeto parámetro a optimizar
    );

    // dlib ---------------------------------------------------------------------

    auto const z_0 = len_t{parameters_array_without_units(0)}; // Parámetro de corrección optimizado
    auto const mu_0 = perm_t{parameters_array_without_units(1)}; // Permeavilidad del vacío optimizada

    auto const rse = mag_ind_t{std::sqrt(rss_without_units/(data_vec.size() - 1))}; // Raiz cuadrada de la media de los residuos al cuadrado

    //stdx::println("z_0 = {:%.2Q %q} | mu_0 = {} | RSE = {:%.2Q %q}", z_0, mu_0, rse);
    stdx::println("Factor de corrección z_0 = {}", z_0);
    stdx::println("Permeabilidad magnética mu_0 = {} | mu_0/(4 * pi * 10^-7) = {}", mu_0, mu_0/(4 * std::numbers::pi * 1e-7));
    stdx::println("RSE = {}", rse);

    auto tod = [](auto q)->double {return q.number();}; // Elimina las unidades de los valores por que matplot no trabaja con unidades
    // Vector de Inducción magnética (sin unidades)
    auto const p = data_vec | std::views::transform(&Data::magnetic_induction)
                            | std::views::transform(tod) 
                            | stdx::ranges::to<std::vector<double>>(); 
    // Vector de posiciones (sin unidades)
    auto const Z = data_vec | std::views::transform(&Data::position)
                            | std::views::transform(tod)
                            | stdx::ranges::to<std::vector<double>>();
}