#include <iostream>
#include <vector>
#include "Methods.h"


State euler(const State& state, std::vector<double>& coefs, double h) {

	State derivatives = Rossler_attractor(state.x, state.y, state.z, coefs);

	return{
		state.x + h * derivatives.x,
		state.y + h * derivatives.y,
		state.z + h * derivatives.z
	};
}


State e_cromer(const State& state, std::vector<double>& coefs, double h) {

	State derivatives = Rossler_attractor(state.x, state.y, state.z, coefs);
	double x = state.x + h * derivatives.x;

	derivatives = Rossler_attractor(x, state.y, state.z, coefs);
	double y = state.y + h * derivatives.y;

	derivatives = Rossler_attractor(x, y, state.z, coefs);
	double z = state.z + h * derivatives.z;

	return{ x, y, z };
}


State midpoint(const State& state, std::vector<double>& coefs, double h) {

	State derivatives = Rossler_attractor(state.x, state.y, state.z, coefs);

	double half_x = state.x + (h / 2) * derivatives.x;
	double half_y = state.y + (h / 2) * derivatives.y;
	double half_z = state.z + (h / 2) * derivatives.z;

	derivatives = Rossler_attractor(half_x, half_y, half_z, coefs);

	double x = state.x + h * derivatives.x;
	double y = state.y + h * derivatives.y;
	double z = state.z + h * derivatives.z;

	return{ x, y, z };
}


State vscd(const State& state, std::vector<double>& coefs, double h) {
	double s = coefs[3];
	double h1 = s * h, h2 = (1 - s) * h;

	State derivatives = Rossler_attractor(state.x, state.y, state.z, coefs);
	double half_x = state.x + h1 * derivatives.x;

	derivatives = Rossler_attractor(half_x, state.y, state.z, coefs);
	double half_y = state.y + h1 * derivatives.y;

	derivatives = Rossler_attractor(half_x, half_y, state.z, coefs);
	double half_z = state.z + h1 * derivatives.z;

	double z = (half_z + h2 * coefs[1]) / (1 - h2 * half_x + h2 * coefs[2]);
	double y = (half_y + h2 * half_x) / (1 - h2 * coefs[0]);
	double x = half_x - h2 * y - h2 * z;

	return{ x, y, z };

}


State rk4(const State& state, std::vector<double>& coefs, double h) {

	State derivatives = Rossler_attractor(state.x, state.y, state.z, coefs);
	double k1X = derivatives.x;
	double k1Y = derivatives.y;
	double k1Z = derivatives.z;

	derivatives = Rossler_attractor(state.x + k1X * (h / 2), state.y + k1Y * (h / 2), state.z + k1Z * (h / 2), coefs);
	double k2X = derivatives.x;
	double k2Y = derivatives.y;
	double k2Z = derivatives.z;

	derivatives = Rossler_attractor(state.x + k2X * (h / 2), state.y + k2Y * (h / 2), state.z + k2Z * (h / 2), coefs);
	double k3X = derivatives.x;
	double k3Y = derivatives.y;
	double k3Z = derivatives.z;

	derivatives = Rossler_attractor(state.x + k3X * h, state.y + k3Y * h, state.z + k3Z * h, coefs);
	double k4X = derivatives.x;
	double k4Y = derivatives.y;
	double k4Z = derivatives.z;

	double x = state.x + (h / 6) * (k1X + 2 * k2X + 2 * k3X + k4X);
	double y = state.y + (h / 6) * (k1Y + 2 * k2Y + 2 * k3Y + k4Y);
	double z = state.z + (h / 6) * (k1Z + 2 * k2Z + 2 * k3Z + k4Z);

	return{ x, y, z };
}
