
#include <iostream>
#include <string>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "Bifurcation.h"
#include "System.h"


std::vector<double> findPeaks(std::vector<State>& trajectory) {

    std::vector<double> peaks;

    for (int i = 1; i < (trajectory.size() - 1); i++) {

        if (trajectory[i].x > trajectory[i - 1].x && trajectory[i].x > trajectory[i + 1].x)
            peaks.push_back(trajectory[i].x);
    }

    return peaks;
}

std::vector<BifPoint> build_bifurcation(double start, double end, int resolution, 
    int coefIndex, std::vector<double>& coefs,
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

            // переходный процесс
            run_method(h, TT, method, trajectory, localCoefs);

            // удалить все кроме последней точки
            trajectory = { trajectory.back() };

            run_method(h, CT, method, trajectory, localCoefs);

            auto peaks = findPeaks(trajectory);

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

std::vector<BifPoint> h_build_bifurcation(double start, double end, int resolution, std::vector<double>& coefs,
    State(*method)(const State&, std::vector<double>&, double)
) {
    std::vector<BifPoint> diagram;
    double k = pow(end / start, 1.0 / resolution);

#pragma omp parallel
    {
        std::vector<BifPoint> localDiagram;

#pragma omp for
        for (int i = 0; i <= resolution; i++)
        {
            double local_h = start * pow(k, i);

            std::vector<double> localCoefs = coefs;
            std::vector<State> trajectory = { state };

            // переходный процесс
            run_method(local_h, TT, method, trajectory, localCoefs);

            // удалить все кроме последней точки
            trajectory = { trajectory.back() };

            run_method(local_h, CT, method, trajectory, localCoefs);

            auto peaks = findPeaks(trajectory);

            for (double peak : peaks)
                localDiagram.push_back({ local_h, peak });
        }

#pragma omp critical
        {
            diagram.insert(
                diagram.end(),
                localDiagram.begin(),
                localDiagram.end());
        }
    }

    return diagram;
}


std::vector<MapPoint> ac_build_bifurcation(double a_start, double a_end,
    double c_start, double c_end, int resolution,
    int a_index, int c_index, std::vector<double>& coefs,
    State(*method)(const State&, std::vector<double>&, double)
) {

    std::vector<MapPoint> diagram;
    double a_step = (a_end - a_start) / resolution;
    double c_step = (c_end - c_start) / resolution;


#pragma omp parallel
    {
        std::vector<MapPoint> localDiagram;

#pragma omp for collapse(2)
        for (int i = 0; i < resolution; i++)
        {
            if (i % 10 == 0) std::cout << "\nIndex number: " << i;

            for (int j = 0; j < resolution; j++)
            {

                double a_param = a_start + a_step * i;
                double c_param = c_start + c_step * j;

                std::vector<double> localCoefs = coefs;
                localCoefs[a_index] = a_param;
                localCoefs[c_index] = c_param;

                std::vector<State> trajectory = { state };

                // переходный процесс
                run_method(h, TT, method, trajectory, localCoefs);

                // удалить все кроме последней точки
                trajectory = { trajectory.back() };

                run_method(h, CT, method, trajectory, localCoefs);

                auto peaks = findPeaks(trajectory);

                for (double peak : peaks)
                    localDiagram.push_back({ a_param, c_param, peak });
            }
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

std::vector<MapPoint> ah_build_bifurcation(double a_start, double a_end,
    double h_start, double h_end, int resolution,
    int a_index, std::vector<double>& coefs,
    State(*method)(const State&, std::vector<double>&, double)
) {
    std::vector<MapPoint> diagram;
    double a_step = (a_end - a_start) / resolution;
    double h_step = (h_end - h_start) / resolution;

#pragma omp parallel
    {
        std::vector<MapPoint> localDiagram;
#pragma omp for collapse(2)
        for (int i = 0; i < resolution; i++) {
            if (i % 10 == 0) std::cout << "\nIndex number: " << i;

            for (int j = 0; j < resolution; j++) {
                double a_param = a_start + a_step * i;
                double local_h = h_start + h_step * j;

                std::vector<double> localCoefs = coefs;
                localCoefs[a_index] = a_param;

                std::vector<State> trajectory = { state };

                run_method(local_h, TT, method, trajectory, localCoefs);
                trajectory = { trajectory.back() };

                run_method(local_h, CT, method, trajectory, localCoefs);

                auto peaks = findPeaks(trajectory);

                for (double peak : peaks) {
                    localDiagram.push_back({ a_param, local_h, peak });
                }
            }
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

    for (const auto& p : data) {
        file << p.param << "," << p.peak << "\n";
    }
}

void save_2d_csv(const std::string& filename,
    const std::vector<MapPoint>& data,
    const std::string& xName,
    const std::string& yName)
{
    std::ofstream file(filename);

    file << xName << "," << yName << ",peak\n";

    for (const auto& p : data)
    {
        file << p.p1 << ","
            << p.p2 << ","
            << p.peak << "\n";
    }
}
