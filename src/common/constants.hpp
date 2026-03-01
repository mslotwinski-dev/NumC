#pragma once

namespace numc {

/// @defgroup constants Mathematical Constants
/// @ingroup common
/// @brief Collection of fundamental mathematical constants (π, e, φ, γ, etc.) and their related values (e.g., logarithms, square roots) for use
/// across the library.
/// @{

/// @brief Liczba Pi (π) - stosunek obwodu koła do jego średnicy.
constexpr double PI = 3.141592653589793238462643383279502884;

/// @brief Dwa Pi (2π), znane również jako Tau (τ) - pełny kąt w radianach.
constexpr double TWO_PI = 6.283185307179586476925286766559005768;

/// @brief Pi podzielone przez 2 (π/2) - kąt prosty w radianach.
constexpr double PI_2 = 1.570796326794896619231321691639751442;

/// @brief Pi podzielone przez 4 (π/4) - kąt 45 stopni w radianach.
constexpr double PI_4 = 0.785398163397448309615660845819875721;

/// @brief Odwrotność Pi (1/π).
constexpr double INV_PI = 0.318309886183790671537767526745028724;

/// @brief Odwrotność dwóch Pi (1/2π).
constexpr double INV_TWO_PI = 0.159154943091895335768883763372514362;

// ============================================================================
// PODSTAWA LOGARYTMU NATURALNEGO I LOGARYTMY
// ============================================================================

/// @brief Liczba Eulera (e) - podstawa logarytmu naturalnego.
constexpr double E = 2.718281828459045235360287471352662497;

/// @brief Logarytm naturalny z 2 (ln 2).
constexpr double LN2 = 0.693147180559945309417232121458176568;

/// @brief Logarytm naturalny z 10 (ln 10).
constexpr double LN10 = 2.302585092994045684017991454684364207;

/// @brief Logarytm przy podstawie 2 z liczby Eulera (log2 e).
constexpr double LOG2E = 1.442695040888963407359924681001892137;

/// @brief Logarytm przy podstawie 10 z liczby Eulera (log10 e).
constexpr double LOG10E = 0.434294481903251827651128918916605082;

// ============================================================================
// PIERWIASTKI
// ============================================================================

/// @brief Pierwiastek kwadratowy z 2 (√2) - przekątna kwadratu o boku 1.
constexpr double SQRT2 = 1.414213562373095048801688724209698078;

/// @brief Odwrotność pierwiastka kwadratowego z 2 (1/√2), co jest równe √2/2.
constexpr double INV_SQRT2 = 0.707106781186547524400844362104849039;

/// @brief Pierwiastek kwadratowy z 3 (√3) - przekątna sześcianu o boku 1.
constexpr double SQRT3 = 1.732050807568877293527446341505872366;

// ============================================================================
// INNE WAŻNE STAŁE MATEMATYCZNE
// ============================================================================

/// @brief Złoty podział (φ, Phi) - (1 + √5)/2.
constexpr double PHI = 1.618033988749894848204586834365638117;

/// @brief Stała Eulera-Mascheroniego (γ, gamma).
constexpr double EGAMMA = 0.577215664901532860606512090082402431;

/// @brief Stała Catalana (G) - spotykana w kombinatoryce i całkach eliptycznych.
constexpr double CATALAN = 0.915965594177219015054603514932384110;

/// @} // Koniec grupy common

}  // namespace numc