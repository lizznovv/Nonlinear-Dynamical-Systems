#pragma once
#include "System.h"

State euler(const State& state, std::vector<double>& coefs, double h);
State e_cromer(const State& state, std::vector<double>& coefs, double h);
State midpoint(const State& state, std::vector<double>& coefs, double h);
State vscd(const State& state, std::vector<double>& coefs, double h);
State rk4(const State& state, std::vector<double>& coefs, double h);