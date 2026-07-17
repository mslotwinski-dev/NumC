#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../statistics/random.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace machine_learning {

/// @addtogroup machine_learning
/// @{

/// @brief Activation functions for neural networks.
enum class Activation {
  Sigmoid,
  ReLU,
  Tanh,
  Linear
};

/// @brief A fully connected (Dense) layer in a Multi-Layer Perceptron.
template <typename T = double>
class DenseLayer {
 private:
  size_t in_features;
  size_t out_features;
  tensor<T> weights; // shape (in_features, out_features)
  vector<T> biases;  // shape (out_features)
  Activation activation;

  // Caches for backpropagation
  vector<T> last_input;
  vector<T> last_z; // Before activation
  vector<T> last_a; // After activation

  // Adam parameters
  tensor<T> m_weights, v_weights;
  vector<T> m_biases, v_biases;
  int t_step;

 public:
  DenseLayer(size_t in_size, size_t out_size, Activation act = Activation::Sigmoid, unsigned seed = 42)
      : in_features(in_size), out_features(out_size), activation(act), t_step(0) {
    // Initialize weights (Xavier/Glorot uniform initialization)
    T limit = std::sqrt(T(6.0) / T(in_size + out_size));
    weights = numc::random::uniform_tensor<T>(in_size, out_size, -limit, limit, seed);
    biases = vector<T>(out_size); // initialized to 0

    m_weights = tensor<T>({in_size, out_size}); // zeros
    v_weights = tensor<T>({in_size, out_size}); // zeros
    m_biases = vector<T>(out_size);
    v_biases = vector<T>(out_size);
  }

  /// @brief Forward pass for a single sample.
  vector<T> forward(const vector<T>& input) {
    last_input = input;
    last_z = vector<T>(out_features);
    for (size_t j = 0; j < out_features; ++j) {
      T sum = biases[j];
      for (size_t i = 0; i < in_features; ++i) {
        sum += input[i] * weights(i, j);
      }
      last_z[j] = sum;
    }

    last_a = vector<T>(out_features);
    for (size_t j = 0; j < out_features; ++j) {
      switch (activation) {
        case Activation::Sigmoid:
          last_a[j] = T(1.0) / (T(1.0) + std::exp(-last_z[j]));
          break;
        case Activation::ReLU:
          last_a[j] = std::max(T(0.0), last_z[j]);
          break;
        case Activation::Tanh:
          last_a[j] = std::tanh(last_z[j]);
          break;
        case Activation::Linear:
        default:
          last_a[j] = last_z[j];
          break;
      }
    }
    return last_a;
  }

  /// @brief Backward pass for a single sample with Adam optimization.
  vector<T> backward(const vector<T>& dL_da, T learning_rate, bool use_adam = true) {
    vector<T> dL_dz(out_features);
    for (size_t j = 0; j < out_features; ++j) {
      T da_dz;
      switch (activation) {
        case Activation::Sigmoid: {
          T s = last_a[j];
          da_dz = s * (T(1.0) - s);
          break;
        }
        case Activation::ReLU:
          da_dz = last_z[j] > T(0.0) ? T(1.0) : T(0.0);
          break;
        case Activation::Tanh: {
          T t = last_a[j];
          da_dz = T(1.0) - t * t;
          break;
        }
        case Activation::Linear:
        default:
          da_dz = T(1.0);
          break;
      }
      dL_dz[j] = dL_da[j] * da_dz;
    }

    // Gradients for weights and biases
    tensor<T> dL_dw({in_features, out_features});
    for (size_t i = 0; i < in_features; ++i) {
      for (size_t j = 0; j < out_features; ++j) {
        dL_dw(i, j) = last_input[i] * dL_dz[j];
      }
    }
    vector<T> dL_db = dL_dz;

    // Gradient to pass back to the previous layer
    vector<T> dL_dinput(in_features);
    for (size_t i = 0; i < in_features; ++i) {
      T sum = T(0.0);
      for (size_t j = 0; j < out_features; ++j) {
        sum += dL_dz[j] * weights(i, j);
      }
      dL_dinput[i] = sum;
    }

    if (use_adam) {
      t_step++;
      T beta1 = T(0.9), beta2 = T(0.999), epsilon = T(1e-8);
      T lr_t = learning_rate * std::sqrt(T(1.0) - std::pow(beta2, t_step)) / (T(1.0) - std::pow(beta1, t_step));

      for (size_t i = 0; i < in_features; ++i) {
        for (size_t j = 0; j < out_features; ++j) {
          m_weights(i, j) = beta1 * m_weights(i, j) + (T(1.0) - beta1) * dL_dw(i, j);
          v_weights(i, j) = beta2 * v_weights(i, j) + (T(1.0) - beta2) * dL_dw(i, j) * dL_dw(i, j);
          weights(i, j) -= lr_t * m_weights(i, j) / (std::sqrt(v_weights(i, j)) + epsilon);
        }
      }
      for (size_t j = 0; j < out_features; ++j) {
        m_biases[j] = beta1 * m_biases[j] + (T(1.0) - beta1) * dL_db[j];
        v_biases[j] = beta2 * v_biases[j] + (T(1.0) - beta2) * dL_db[j] * dL_db[j];
        biases[j] -= lr_t * m_biases[j] / (std::sqrt(v_biases[j]) + epsilon);
      }
    } else {
      // Apply Gradient Descent update
      for (size_t i = 0; i < in_features; ++i) {
        for (size_t j = 0; j < out_features; ++j) {
          weights(i, j) -= learning_rate * dL_dw(i, j);
        }
      }
      for (size_t j = 0; j < out_features; ++j) {
        biases[j] -= learning_rate * dL_db[j];
      }
    }

    return dL_dinput;
  }
};

/// @brief Multi-Layer Perceptron (Basic Neural Network).
template <typename T = double>
class MLP {
 private:
  std::vector<DenseLayer<T>> layers;

 public:
  MLP() = default;

  /// @brief Add a Dense layer to the network.
  void add_layer(size_t in_size, size_t out_size, Activation act = Activation::Sigmoid, unsigned seed = 42) {
    layers.emplace_back(in_size, out_size, act, seed);
  }

  /// @brief Forward pass through the entire network.
  vector<T> predict(const vector<T>& input) {
    vector<T> out = input;
    for (auto& layer : layers) {
      out = layer.forward(out);
    }
    return out;
  }

  /// @brief Train the network using Mean Squared Error (MSE) loss and Stochastic Gradient Descent (SGD) or Adam.
  void train(const std::vector<vector<T>>& X, const std::vector<vector<T>>& Y, int epochs, T learning_rate, bool use_adam = true) {
    if (X.size() != Y.size() || X.empty()) {
      Log::Warn("MLP::train: Invalid data sizes.");
      return;
    }
    size_t n_samples = X.size();
    for (int epoch = 0; epoch < epochs; ++epoch) {
      T total_loss = T(0.0);
      for (size_t i = 0; i < n_samples; ++i) {
        // Forward
        vector<T> pred = predict(X[i]);
        
        // Compute loss (MSE) and initial gradient
        vector<T> dL_da(pred.size());
        for (size_t j = 0; j < pred.size(); ++j) {
          T diff = pred[j] - Y[i][j];
          total_loss += diff * diff;
          dL_da[j] = T(2.0) * diff; // derivative of (pred - y)^2 is 2*(pred - y)
        }

        // Backward
        for (int l = static_cast<int>(layers.size()) - 1; l >= 0; --l) {
          dL_da = layers[l].backward(dL_da, learning_rate, use_adam);
        }
      }
    }
  }
};

/// @}

}  // namespace machine_learning
}  // namespace numc
