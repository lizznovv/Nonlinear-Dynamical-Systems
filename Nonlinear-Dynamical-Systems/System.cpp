#include "System.h"

double CT = 300, TT = 300, h = 0.01, h10 = 0.001, k = 10;
State state = { -0.7, -0.7, 1.0 };


State Rossler_attractor(double x, double y, double z, std::vector<double>& coefs) {

    return {
        -y - z,
        x + coefs[0] * y,
        coefs[1] + z * (x - coefs[2])
    };

}


void run_method(double h, double CT,
    State(*method)(const State& state, std::vector<double>& coefs, double h),
    std::vector<State>& trajectory, std::vector<double>& coefs
) {
    int N = (CT / h);

    for (int i = 0; i < N; i++) {

        int length = trajectory.size();
        State state = method(trajectory[length - 1], coefs, h);
        trajectory.push_back(state);
    }
}