#pragma once

double LLE_calculation(double h, State start_state, double CT, double NT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
);

std::vector<LLEPoint> build_LLE_diagram(double start, double end, int resolution,
    int coefIndex, double CT, double NT, std::vector<double>& coefs,
    State start_state, State(*method)(const State&, std::vector<double>&, double),
    bool varyStep = false
);


std::vector<LLEMapPoint> ac_build_LLE_diagram(double a_start, double a_end,
    double c_start, double c_end, int resolution,
    int a_index, int c_index, double CT, double NT, std::vector<double>& coefs,
    State start_state, State(*method)(const State&, std::vector<double>&, double)
);

std::vector<LLEMapPoint> ah_build_LLE_diagram(double a_start, double a_end,
    double h_start, double h_end, int resolution, int a_index, double CT, double NT,
    std::vector<double>& coefs, State start_state,
    State(*method)(const State&, std::vector<double>&, double)
);

void save_lle_csv(const std::string& filename, const std::vector<LLEPoint>& data);

void save_lle_map_csv(const std::string& filename,
    const std::vector<LLEMapPoint>& data);