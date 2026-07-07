#pragma once


double phase_volume(double h, State start_state, double CT, double TT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
);

std::vector<VolumePoint> build_volume_diagram(double start, double end, int resolution,
	int coefIndex, double CT, double NT, std::vector<double>& coefs,
	State start_state, State(*method)(const State&, std::vector<double>&, double),
	bool varyStep
);

void save_volume_csv(const std::string& filename,
	const std::vector<VolumePoint>& data
);

std::vector<BifPoint> build_contunue_diagram(double h, double windowTime,
	State start_state, double CT, double TT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
);