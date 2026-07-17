#pragma once

#include "inc.hpp"

#include "statistics/histogram.hpp"
#include "statistics/sampling.hpp"
#include "statistics/random.hpp"
#include "statistics/error_propagation.hpp"
#include "statistics/param.hpp"
#include "statistics/distribution.hpp"
#include "statistics/fitting.hpp"

// Utility

#include "utility/log.hpp"
#include "utility/data_io.hpp"

// Common

#include "common/autodiff.hpp"
#include "common/complex.hpp"
#include "common/constants.hpp"
#include "common/function.hpp"
#include "common/polynomial.hpp"
#include "common/sparse.hpp"
#include "common/tensor.hpp"
#include "common/vector.hpp"

// Linear Algebra

#include "linear_algebra/solvers.hpp"
#include "linear_algebra/eigen.hpp"
#include "linear_algebra/matrix_ops.hpp"

// Analysis

#include "analysis/functions/elementary.hpp"
#include "analysis/functions/orthogonal.hpp"
#include "analysis/functions/signal.hpp"
#include "analysis/functions/special.hpp"
#include "analysis/complex/elementary.hpp"
#include "analysis/interpolation.hpp"
#include "analysis/differentiation.hpp"
#include "analysis/integration.hpp"
#include "analysis/fft.hpp"
#include "analysis/calculus.hpp"

// Multivariable

#include "analysis/multivariable/mfunc.hpp"
#include "analysis/multivariable/optimization.hpp"

// Optimalization

#include "optimalization/solvers/roots.hpp"
#include "optimalization/solvers/minimization.hpp"
#include "optimalization/solvers/least_squares.hpp"

// Differential Equations

#include "differential_equations/ode.hpp"
#include "differential_equations/bvp.hpp"
#include "differential_equations/pde.hpp"

// Statistics

#include "statistics/distributions/dist.hpp"
#include "statistics/visualization/plot.hpp"
#include "statistics/visualization/plt.hpp"
#include "statistics/stats.hpp"
#include "statistics/random.hpp"

// Machine Learning

#include "machine_learning/clustering.hpp"
#include "machine_learning/decomposition.hpp"
#include "machine_learning/neural_network.hpp"
#include "machine_learning/decision_tree.hpp"
#include "machine_learning/knn.hpp"

// Graphs

#include "graphs/pathfinding.hpp"
#include "graphs/mst.hpp"

// Number Theory

#include "number_theory/modular.hpp"
#include "number_theory/primes.hpp"

// Geometry

#include "geometry/spatial.hpp"
#include "geometry/curves.hpp"

// Convenience re-exports and sugar syntax (x(), newton(), mean(), …)
#include "convenience.hpp"

