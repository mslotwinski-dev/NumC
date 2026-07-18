<img src="public/logo.png" alt="NumC logo" width="320">

# NumC

**Numeryka w C++ tak piękna, jak w MATLAB-ie — ale szybsza, header-only i bez kompromisów.**

NumC to nowoczesna biblioteka obliczeń naukowych. Definiujesz funkcje jak w notatniku, rysujesz wykresy jednym wywołaniem, a reszta (całki, pochodne, układy równań, ML, FFT…) jest na wyciągnięcie ręki.

```cpp
#include <numc>
using namespace numc;

int main() {
  auto x = funcv::x();

  func f1 = sin(x) * (x / 2.0);
  func f2 = exp(-0.1 * x) * cos(x) * 5.0;

  vector data_x = {-8.0, -5.5, -2.1, 0.0, 3.14, 6.28, 9.0};
  vector data_y = {3.5, -2.1, -1.0, 5.2, -4.8, 5.1, -4.0};

  plot dashboard("Analiza drgań");
  dashboard.set_xlabel("Czas [ms]").set_ylabel("Amplituda [V]");
  dashboard.add(f1, -10.0, 10.0, "Model teoretyczny");
  dashboard.add(f2, -10.0, 10.0, "Model z tłumieniem");
  dashboard.add(data_x, data_y, "Pomiary z czujnika");
  dashboard.show();

  plt::show({f1, f2}, -10.0, 10.0, "Szybki wykres", "Czas [ms]", "Amplituda [V]",
            data_x, data_y, "Pomiary z czujnika");

  vector v1 = {1.0, 2.0, 3.0};
  vector v2 = {4.0, 5.0, 6.0};
  vector v3 = v1 + v2;    // dodawanie
  vector v4 = v1 * 2.0;   // skalar × wektor
  double dot = v1 * v2;   // iloczyn skalarny
  vector cross = v1 ^ v2; // iloczyn wektorowy (3D)
}
```

---

## Szybki start

```bash
git clone https://github.com/mslotwinski-dev/NumC.git
cd NumC
cmake -S . -B build
cmake --build build
./build/examples/example_showcase    # Linux / macOS
build\examples\Debug\example_showcase.exe   # Windows
```

Wystarczy `#include <numc>` — cała biblioteka to nagłówki.

---

## Co potrafi NumC?

| Obszar | Przykład użycia |
|--------|-----------------|
| **Wyrażenia** | `func f = sin(x) * exp(-x^2);` |
| **Wektory** | `v1 + v2`, `v1 * v2` (dot), `v1 ^ v2` (cross) |
| **Wykresy** | `plot`, `plt::show`, scatter, histogramy |
| **Rachunek** | `derivative(f)`, `f.integral(a,b)`, `limit(f, x0)` |
| **Korzenie** | `newton(f, x0)`, `bisection(f, a, b)`, `secant(…)` |
| **Statystyka** | `mean(v)`, `linear_regression(x, y)` |
| **Algebra liniowa** | tensory, SVD, QR, układy równań |
| **Równania różniczkowe** | ODE, BVP, PDE (Laplace) |
| **FFT & sygnały** | `fft(v)`, `convolve(a, b)` |
| **ML** | k-means, PCA, MLP, drzewa decyzyjne, k-NN |
| **Auto-diff** | `auto_diff([](auto x){ return sin(x)*exp(x); }, 1.0)` |

Pełny test wszystkich modułów: `examples/numerical_methods.cpp`.

---

## Filozofia API

NumC stawia na **przyjemność użytkowania**:

- **Domyślny `double`** — piszesz `vector`, `func`, `plot`, nie `vector<double>`.
- **Operator `^`** — potęgowanie funkcji i iloczyn wektorowy, jak w notatniku.
- **Łańcuchowe wywołania** — `plot.set_xlabel("t").set_ylabel("A").add(f, -10, 10).show()`.
- **Jedna linia na wykres** — `plt::show({f, g}, -10, 10, "Tytuł", "x", "y")`.
- **Bez zbędnych namespace'ów** — `newton`, `mean`, `derivative` dostępne od razu po `using namespace numc`.

---

## Przykłady w repozytorium

| Plik | Opis |
|------|------|
| `examples/showcase.cpp` | Flagowy demo — drgania, wektory, statystyka, wykresy |
| `examples/plots.cpp` | Wizualizacja funkcji i scatter |
| `examples/functions.cpp` | Budowanie wyrażeń funkcyjnych |
| `examples/vectors.cpp` | Arytmetyka wektorowa |
| `examples/roots.cpp` | Metody poszukiwania zer |
| `examples/numerical_methods.cpp` | Kompletny test wszystkich modułów |

---

## Licencja

MIT — używaj swobodnie w projektach naukowych i komercyjnych.

---

<p align="center"><i>NumC — bo obliczenia numeryczne mogą być piękne.</i></p>
