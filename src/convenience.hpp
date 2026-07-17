#pragma once

/// @file convenience.hpp
/// @brief Sugar syntax and top-level re-exports for the most pleasant NumC experience.
/// Include via <numc> — no need to pull this in manually.

namespace numc {

/// @brief The independent variable x — shorthand for func<>::x().
/// @code
/// auto x = numc::x;
/// func f = sin(x) * exp(-0.1 * x);
/// @endcode

// --- Root finding (no extra using-namespace needed) ---
using optimalization::bisection;
using optimalization::regula_falsi;
using optimalization::fixed_point;
using optimalization::newton;
using optimalization::halley;
using optimalization::secant;
using optimalization::inverse_quadratic;

// --- Calculus ---
using analysis::derivative;
using analysis::integral;
using analysis::limit;
using analysis::tabulate;
using analysis::taylor_coeffs;
using analysis::taylor_eval;
using analysis::arc_length;
using analysis::simpson_1_3;
using analysis::adaptive_simpson;

// --- Statistics essentials ---
using statistics::mean;
using statistics::median;
using statistics::std_dev;
using statistics::variance;
using statistics::linear_regression;

}  // namespace numc
