#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::number_theory;

// Funkcja pomocnicza: generowanie kluczy RSA (uproszczone z ustaloną wielkością)
// W praktyce liczby p i q powinny być losowane
void generuj_klucze_rsa(uint64_t p, uint64_t q, uint64_t& klucz_publiczny_e, uint64_t& modul_n, uint64_t& klucz_prywatny_d) {
  modul_n = p * q;
  uint64_t phi = (p - 1) * (q - 1); // Funkcja Eulera

  // Wybieramy e względnie pierwsze z phi. Tradycyjnie 65537, ale weźmy np. 17 dla małych liczb
  klucz_publiczny_e = 17;
  while (gcd(klucz_publiczny_e, phi) != 1) {
    klucz_publiczny_e += 2;
  }

  // d to odwrotność modularna e mod phi
  klucz_prywatny_d = mod_inverse(klucz_publiczny_e, phi);
}

int main() {
  std::cout << "=== Projekt Number Theory: Szyfrowanie RSA (Kryptografia) ===\n" << std::endl;

  std::cout << "1. Generowanie duzych liczb pierwszych (używamy Miller-Rabin test):" << std::endl;
  
  // Zakładamy, że wygenerowaliśmy te liczby korzystając np z losowania + is_prime_miller_rabin
  uint64_t p = 61; // Małe dla czytelności testu
  uint64_t q = 53;
  
  std::cout << "   p = " << p << " (czy pierwsza? " << (is_prime_miller_rabin(p) ? "TAK" : "NIE") << ")" << std::endl;
  std::cout << "   q = " << q << " (czy pierwsza? " << (is_prime_miller_rabin(q) ? "TAK" : "NIE") << ")" << std::endl;

  uint64_t e, n, d;
  generuj_klucze_rsa(p, q, e, n, d);

  std::cout << "\n2. Wygenerowano klucze RSA:" << std::endl;
  std::cout << "   Klucz publiczny (e, n) = (" << e << ", " << n << ")" << std::endl;
  std::cout << "   Klucz prywatny d = " << d << std::endl;

  std::cout << "\n3. Proces Szyfrowania (Modular Exponentiation):" << std::endl;
  uint64_t oryginalna_wiadomosc = 65; // 'A' w kodzie ASCII
  
  // Szyfrowanie: C = M^e (mod n)
  uint64_t szyfrogram = mod_pow(oryginalna_wiadomosc, e, n);
  
  std::cout << "   Wiadomosc (M) = " << oryginalna_wiadomosc << std::endl;
  std::cout << "   Szyfrogram (C) = " << szyfrogram << std::endl;

  std::cout << "\n4. Proces Deszyfrowania:" << std::endl;
  
  // Deszyfrowanie: M = C^d (mod n)
  uint64_t odszyfrowana_wiadomosc = mod_pow(szyfrogram, d, n);
  
  std::cout << "   Odszyfrowana Wiadomosc = " << odszyfrowana_wiadomosc << std::endl;
  if (oryginalna_wiadomosc == odszyfrowana_wiadomosc) {
    std::cout << "   [SUKCES] Wiadomosc zdeszyfrowano prawidlowo!" << std::endl;
  }

  return 0;
}
