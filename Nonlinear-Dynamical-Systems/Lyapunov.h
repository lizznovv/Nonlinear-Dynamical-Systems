#pragma once

double LLE_calculation(int count, State start_state, double CT, double NT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
);

std::vector<LLEPoint> build_LLE_diagram(double start, double end, int resolution,
	int coefIndex, int count, double CT, double NT,
	std::vector<double>& coefs, State start_state,
	State(*method)(const State&, std::vector<double>&, double)
);

void save_lle_csv(const std::string& filename, const std::vector<LLEPoint>& data);
