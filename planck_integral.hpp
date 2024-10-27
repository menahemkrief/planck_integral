#ifndef PLANCK_INTEGRAL_HPP
#define PLANCK_INTEGRAL_HPP

#include <boost/math/special_functions/pow.hpp>
#include <cmath>
#include <cassert>

#include "../units/units.hpp"

/**
 * @brief An accurate method for the calculation of integrals of the dimensionless Planck function in a given range [a,b].
 * Based on a variation of B.A. Clark "Computing multigroup radiation integrals using polylogarithm-based methods" JCP 70(2), 311-329 (1987)
 */

namespace planck_integral {

inline double Clark_Taylor(double const x){
    // Clark 1987 eq. 32 without the pre-factor
    // Analytic integral of the Taylor approximation at small x
    double const x2 = x*x;
    double const x3 = x2*x;

    return x3*(1./3. + x*(-1./8.+x*(1./60.+x2*(-1./5040.+x2*(1./272160.+x2*(-1./13305600.+x2/622702080.))))));
}

// // For N=21 terms, x=1 is a good stitching point giving accuracy better than ~1e-10 for all x. (see Clark 1987 Fig. 3). However, it is substantially slower. 
// static int constexpr N_clark = 21;
// static double constexpr x_clark = 2.;

// For N=5 terms, x=2 is a good stitching point giving accuracy better than ~1e-5 for all x. (see Clark 1987 Fig. 3)
static int constexpr N_clark = 5;
static double constexpr x_clark = 2.;


inline double Clark_series(double const x){
    // Clark 1987 eq. 38 without pre-factor of the brackets and without the constant - to prevent loss of accuracy for [a,b] with b>a>10.
    double const x2 = x*x;
    double const x3 = x2*x;
    double sum=0.;

    double const exp = std::exp(-x);
    double expn = exp;
    for(int n=1; n <= N_clark; ++n){
        double const in = 1. / static_cast<double>(n);
        sum += in*(x3+in*(3.*x2+6.*in*(x+in)))*expn;
        expn *= exp;
    }
    return -sum;
}

/**
 * @brief The definite integral of the dimensionless normalized Planck function b(x)=15/PI^4 * x^3/(e^x-1) on the range [a,b].
 * @param a lower integration range
 * @param b upper integration range
 * @return double the resulting integral on [a,b]
 */
inline double planck_integral(double const a, double const b){
    assert(a<b);
    using boost::math::pow;
    static double constexpr coeff = 15./pow<4>(M_PI);

    // based on Clark 1987 Fig. 3, we stitch the two approximations at x=x_clark which depends
    // on the order of the exponential expansion, N_clark.
    if(a > x_clark) return coeff*(Clark_series(b) - Clark_series(a));
    if(b < x_clark) return coeff*(Clark_Taylor(b) - Clark_Taylor(a));
    
    return 1.0 + coeff*(Clark_series(b) - Clark_Taylor(a));
}

/**
 * @brief Calculates the integral of the dimensional Planck function at temperature T, in the energy
 * range [E_low, E_high]. for E_low<<kT and E_high>>kT the result should approach aT^4.
 * 
 * @param E_low the lower energy group boundary [erg]
 * @param E_high the upper energy group boundary [erg]
 * @param T the temperature [K]
 * @return double the integral of the Planck function in the given temperature and energy range (erg/cm^3).
 */
inline double planck_energy_density_group_integral(double const E_low, double const E_high, double const T) {
    using boost::math::pow;
    assert(E_low < E_high);

    if(T < 1e50*std::numeric_limits<double>::min()){
        return 0.0;
    }

    double const kT = units::k_boltz * T;
    double const a = E_low / kT;
    double const b = E_high / kT;
    return units::arad * pow<4>(T) * planck_integral(a, b);
}

} // namespace planck_integral

#endif //PLANK_INTEGRAL_HPP