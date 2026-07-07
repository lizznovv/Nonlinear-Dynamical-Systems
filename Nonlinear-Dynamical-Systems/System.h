#pragma once
#include <vector>

extern double CT;
extern double TT;
extern double h;
extern double h10;
extern double k;


struct State
{
    double x;
    double y;
    double z;
};

extern State state;

struct LLEPoint
{
    double param;
    double lambda;
};

struct LLEMapPoint
{
    double param1;
    double param2;
    double lambda;
};

struct BifPoint
{
    double param;
    double peak;
};

struct MapPoint
{
    double p1;
    double p2;
    double peak;
};


State Rossler_attractor(
    double x,
    double y,
    double z,
    std::vector<double>& coefs
); 
void run_method(double h, double CT,
    State(*method)(const State& state, std::vector<double>& coefs, double h),
    std::vector<State>& trajectory, std::vector<double>& coefs
);