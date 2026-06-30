#include <iostream>
#include <vector>
#include "System.h"


State Rossler_attractor(double x, double y, double z, std::vector<double>& coefs) {

    return {
        -y - z,
        x + coefs[0] * y,
        coefs[1] + z * (x - coefs[2])
    };

}