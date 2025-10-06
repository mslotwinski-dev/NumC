#include "log.hpp"

std::string Log::Color(Colors i) {
  switch (i) {
    case RED: return "\033[31m";                 // Czerwony
    case GREEN: return "\033[32m";               // Zielony
    case YELLOW: return "\033[33m";              // Żółty
    case BLUE: return "\033[34m";                // Niebieski
    case MAGENTA: return "\033[35m";             // Fioletowy
    case CYAN: return "\033[36m";                // Cyjan
    case ORANGE: return "\033[38;2;255;165;0m";  // Pomarańczowy
    case LIGHT: return "\033[93m";               // Jasny zolty
    case NONE: return "\033[37m";                // Biały
    default: return "";                          // Brak koloru (domyślny)
  }
}

void Log::Info(std::string message) { std::cout << Color(CYAN) << "[Info] " << message << Color(NONE) << std::endl; }

void Log::Warn(std::string message) { std::cout << Color(YELLOW) << "[Warning] " << message << Color(NONE) << std::endl; }

void Log::Error(std::string message) { std::cout << Color(RED) << "[Error] " << message << Color(NONE) << std::endl; }