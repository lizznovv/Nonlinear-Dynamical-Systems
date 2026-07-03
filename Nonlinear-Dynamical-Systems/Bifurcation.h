#pragma once
#include "Methods.h"


void run_method(double h, double CT,
    State(*method)(const State& state, std::vector<double>& coefs, double h),
    std::vector<State>& trajectory, std::vector<double>& coefs
);

std::vector<double> findPeaks(std::vector<State>& trajectory);

std::vector<BifPoint> build_bifurcation(double start, double end, int resolution, int coefIndex, std::vector<double>& coefs,
    State(*method)(const State& state, std::vector<double>& coefs, double h)
);

std::vector<BifPoint> h_build_bifurcation(double start, double end, int resolution, std::vector<double>& coefs,
    State(*method)(const State&, std::vector<double>&, double)
);

std::vector<MapPoint> ac_build_bifurcation(double a_start, double a_end,
    double c_start, double c_end, int resolution,
    int a_index, int c_index, std::vector<double>& coefs,
    State(*method)(const State&, std::vector<double>&, double)
);

std::vector<MapPoint> ah_build_bifurcation(double a_start, double a_end,
    double h_start, double h_end, int resolution,
    int a_index, std::vector<double>& coefs,
    State(*method)(const State&, std::vector<double>&, double)
);

void save_csv(const std::string& filename, const std::vector<BifPoint>& data);

void save_2d_csv(const std::string& filename,
    const std::vector<MapPoint>& data,
    const std::string& xName,
    const std::string& yName
);
