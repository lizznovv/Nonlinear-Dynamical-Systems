// Dynamic_Systems_Bifurcation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "System.h"
#include "Methods.h"

double CT = 300, TT = 300, h = 0.01, h10 = 0.001, k = 10;
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
){
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
            if(i % 10 == 0) std::cout << "\nIndex number: " << i;

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

int main()
{
    State state = { -0.7, -0.7, 1.0 };
    std::vector<double> coefs = { 0.2, 0.2, 5.7, 0.5 };
    /*std::vector<std::string> names = {
        "euler_a.csv", "euler_b.csv", "euler_c.csv",
        "e_cromer_a.csv", "e_cromer_b.csv", "e_cromer_c.csv",
        "midpoint_a.csv", "midpoint_b.csv", "midpoint_c.csv",
        "vscd_a.csv", "vscd_b.csv", "vscd_c.csv", "vscd_s.csv",
        "rk4_a.csv", "rk4_b.csv", "rk4_c.csv",
    };
    */

    /*
    
    std::cout << "\nMethod: Euler";
    std::cout << "\na parameter processing...";
    auto bif_a_euler = build_bifurcation(0.0, 0.36, 1000, 0, coefs, euler);
    save_csv("euler_a.csv", bif_a_euler);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_euler = build_bifurcation(0.1, 2.0, 1000, 1, coefs, euler);
    save_csv("euler_b.csv", bif_b_euler);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_euler = build_bifurcation(2.0, 10.0, 1000, 2, coefs, euler);
    save_csv("euler_c.csv", bif_c_euler);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Euler-Cromer";
    std::cout << "\na parameter processing...";
    auto bif_a_ecromer = build_bifurcation(0.0, 0.36, 1000, 0, coefs, e_cromer);
    save_csv("e_cromer_a.csv", bif_a_ecromer);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_ecromer = build_bifurcation(0.1, 2.0, 1000, 1, coefs, e_cromer);
    save_csv("e_cromer_b.csv", bif_b_ecromer);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_ecromer = build_bifurcation(2.0, 10.0, 1000, 2, coefs, e_cromer);
    save_csv("e_cromer_c.csv", bif_c_ecromer);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Midpoint";
    std::cout << "\na parameter processing...";
    auto bif_a_midpoint = build_bifurcation(0.0, 0.36, 1000, 0, coefs, midpoint);
    save_csv("midpoint_a.csv", bif_a_midpoint);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_midpoint = build_bifurcation(0.1, 2.0, 1000, 1, coefs, midpoint);
    save_csv("midpoint_b.csv", bif_b_midpoint);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_midpoint = build_bifurcation(2.0, 10.0, 1000, 2, coefs, midpoint);
    save_csv("midpoint_c.csv", bif_c_midpoint);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: VSCD";
    std::cout << "\na parameter processing...";
    auto bif_a_vscd = build_bifurcation(0.0, 0.36, 1000, 0, coefs, vscd);
    save_csv("vscd_a.csv", bif_a_vscd);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_vscd = build_bifurcation(0.1, 2.0, 1000, 1, coefs, vscd);
    save_csv("vscd_b.csv", bif_b_vscd);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_vscd = build_bifurcation(2.0, 10.0, 1000, 2, coefs, vscd);
    save_csv("vscd_c.csv", bif_c_vscd);
    std::cout << "\ndone";

    std::cout << "\ns parameter processing...";
    auto bif_s_vscd = build_bifurcation(0.0, 1.0, 1000, 3, coefs, vscd);
    save_csv("vscd_s.csv", bif_s_vscd);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: RK4";
    std::cout << "\na parameter processing...";
    auto bif_a_rk4 = build_bifurcation(0.0, 0.36, 1000, 0, coefs, rk4);
    save_csv("rk4_a.csv", bif_a_rk4);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_rk4 = build_bifurcation(0.1, 2.0, 1000, 1, coefs, rk4);
    save_csv("rk4_b.csv", bif_b_rk4);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_rk4 = build_bifurcation(2.0, 10.0, 1000, 2, coefs, rk4);
    save_csv("rk4_c.csv", bif_c_rk4);
    std::cout << "\ndone\n";
    
    */


/*
    std::cout << "\nMethod: Euler";
    std::cout << "\nh exploring...";
    auto bif_h_euler = h_build_bifurcation(1e-3, 1e-1, 300, coefs, euler);
    save_csv("euler_h.csv", bif_h_euler);
    std::cout << "\ndone";

    std::cout << "\nMethod: Euler-Cromer";
    std::cout << "\nh exploring...";
    auto bif_h_ecromer = h_build_bifurcation(1e-4, 5e-2, 300, coefs, e_cromer);
    save_csv("e_cromer_h.csv", bif_h_ecromer);
    std::cout << "\ndone";

    std::cout << "\nMethod: Midpoint";
    std::cout << "\nh exploring...";
    auto bif_h_midpoint = h_build_bifurcation(1e-4, 5e-2, 300, coefs, midpoint);
    save_csv("midpoint_h.csv", bif_h_midpoint);
    std::cout << "\ndone";

    std::cout << "\nMethod: VSCD";
    std::cout << "\nh exploring...";
    auto bif_h_vscd = h_build_bifurcation(1e-4, 5e-2, 300, coefs, vscd);
    save_csv("vscd_h.csv", bif_h_vscd);
    std::cout << "\ndone";

    std::cout << "\nMethod: RK4";
    std::cout << "\nh exploring...";
    auto bif_h_rk4 = h_build_bifurcation(1e-4, 5e-2, 300, coefs, rk4);
    save_csv("rk4_h.csv", bif_h_rk4);
    std::cout << "\ndone";

    */


/*
    std::cout << "\nMethod: Euler";
    auto euler_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, euler);
    save_2d_csv("ac_euler.csv", euler_bif_ac, "a", "c");
    std::cout << "\ndone";

    std::cout << "\nMethod: Euler-Cromer";
    auto ecromer_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, e_cromer);
    save_2d_csv("ac_ecromer.csv", ecromer_bif_ac, "a", "c");
    std::cout << "\ndone";
    
    std::cout << "\nMethod: Midpoint";
    auto midpoint_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, midpoint);
    save_2d_csv("ac_midpoint.csv", midpoint_bif_ac, "a", "c");
    std::cout << "\ndone";
    
    std::cout << "\nMethod: VSCD";
    auto vscd_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, vscd);
    save_2d_csv("ac_vscd.csv", vscd_bif_ac, "a", "c");
    std::cout << "\ndone";

    std::cout << "\nMethod: RK4";
    auto rk4_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, rk4);
    save_2d_csv("ac_rk4.csv", rk4_bif_ac, "a", "c");
    std::cout << "\ndone";
    */


    std::cout << "\nMethod: Euler (a vs h)";
    auto euler_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, euler);
    save_2d_csv("ah_euler.csv", euler_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: Euler-Cromer (a vs h)";
    auto ecromer_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, e_cromer);
    save_2d_csv("ah_ecromer.csv", ecromer_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: Midpoint (a vs h)";
    auto midpoint_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, midpoint);
    save_2d_csv("ah_midpoint.csv", midpoint_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: VSCD (a vs h)";
    auto vscd_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, vscd);
    save_2d_csv("ah_vscd.csv", vscd_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: RK4 (a vs h)";
    auto rk4_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, rk4);
    save_2d_csv("ah_rk4.csv", rk4_bif_ah, "a", "h");
    std::cout << " done\n";

    return 0;
}

