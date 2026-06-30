// Dynamic_Systems_Bifurcation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "System.h"
#include "Methods.h"

double CT = 300, h = 0.01, h10 = 0.001, k = 10;
State state = { -0.7, -0.7, 1.0 };


void run_method(double h, double CT, 
    State (*method)(const State& state, std::vector<double>& coefs, double h),
    std::vector<State>& trajectory, std::vector<double>& coefs
) {
    int N = (CT / h);

    for (int i = 0; i < N; i++) {

        int length = trajectory.size();
        State state = method(trajectory[length - 1], coefs, h);
        trajectory.push_back(state);
    }
}

std::vector<double> findPeaks(std::vector<State>& trajectory) {

    std::vector<double> peaks;

    for (int i = 1; i < (trajectory.size() - 1); i++) {
        
        if (trajectory[i].x > trajectory[i - 1].x && trajectory[i].x > trajectory[i + 1].x)
            peaks.push_back(trajectory[i].x);
    }

    return peaks;
}

std::vector<BifPoint> build_bifurcation(double start, double end, int resolution, int coefIndex, std::vector<double>& coefs,
    State(*method)(const State& state, std::vector<double>& coefs, double h)
) {

    std::vector<BifPoint> diagram;
    double step = (end - start) / resolution;

    #pragma omp parallel
    {
        std::vector<BifPoint> localDiagram;

        #pragma omp for
        for (int i = 0; i < resolution; i++)
        {
            double param = start + step * i;

            std::vector<double> localCoefs = coefs;
            localCoefs[coefIndex] = param;

            std::vector<State> trajectory = { state };

            run_method(h, CT, method, trajectory, localCoefs);

            int startIndex = trajectory.size() * 0.8;

            std::vector<State> cutTrajectory(trajectory.begin() + startIndex, trajectory.end());

            auto peaks = findPeaks(cutTrajectory);

            for (double peak : peaks)
                localDiagram.push_back({ param, peak });
        }

        #pragma omp critical
        {
            diagram.insert(
                diagram.end(),
                localDiagram.begin(),
                localDiagram.end()
            );
        }
    }


    return diagram;
}

void save_csv(const std::string& filename, const std::vector<BifPoint>& data)
{
    std::ofstream file(filename);
    file << "parameter, peak\n";

    for (const auto& p : data)
    {
        file << p.param
            << ","
            << p.peak
            << "\n";
    }
}
int main()
{
    State state = { -0.7, -0.7, 1.0 };
    std::vector<double> coefs = { 0.2, 0.2, 5.7 };

    auto bif_a_euler = build_bifurcation(0.2, 0.35, 1000, 0, coefs, euler);
    save_csv("euler_a.csv", bif_a_euler);

    auto bif_b_euler = build_bifurcation(0.2, 0.35, 1000, 1, coefs, euler);
    save_csv("euler_b.csv", bif_b_euler);

    auto bif_c_euler = build_bifurcation(2.0, 10.0, 1000, 2, coefs, euler);
    save_csv("euler_c.csv", bif_c_euler);

    return 0;
}

