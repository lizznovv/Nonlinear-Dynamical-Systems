#pragma once

struct State
{
    double x;
    double y;
    double z;
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