
#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "System.h"
#include "Lyapunov.h"

double LLE_calculation(int count, State start_state, double CT, double NT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
) {
	State new_state;
	double eps = 1e-8;
	double dt = CT / count;

	new_state.x = start_state.x + eps * (double)rand() / RAND_MAX;
	new_state.y = start_state.y + eps * (double)rand() / RAND_MAX;
	new_state.z = start_state.z + eps * (double)rand() / RAND_MAX;

	std::vector<double> localCoefs = coefs;

	std::vector<State> trajectory1 = { start_state };

	run_method(dt, NT, method, trajectory1, localCoefs);
	// удалить все кроме последней точки
	trajectory1 = { trajectory1.back() };

	std::vector<State> trajectory2 = { new_state };

	run_method(dt, NT, method, trajectory2, localCoefs);
	// удалить все кроме последней точки
	trajectory2 = { trajectory2.back() };

	double sum = 0;

	for (int i = 0; i < count; i++) {

		int length1 = trajectory1.size();
		State state1 = method(trajectory1[length1 - 1], localCoefs, dt);

		int length2 = trajectory2.size();
		State state2 = method(trajectory2[length2 - 1], localCoefs, dt);

		if (std::isnan(state1.x) || std::isnan(state1.y) || std::isnan(state1.z))
		{
			std::cout << "state1 is NaN\n";
			break;
		}

		if (std::isnan(state2.x) || std::isnan(state2.y) || std::isnan(state2.z))
		{
			std::cout << "state2 is NaN\n";
			break;
		}

		double distance = std::sqrt(pow(state1.x - state2.x, 2) + pow(state1.y - state2.y, 2) + pow(state1.z - state2.z, 2));
		
		if (distance < 1e-15)
			distance = 1e-15;

		sum += std::log(distance / eps);

		state2.x = state1.x + eps * ((state2.x - state1.x) / distance);
		state2.y = state1.y + eps * ((state2.y - state1.y) / distance);
		state2.z = state1.z + eps * ((state2.z - state1.z) / distance);

		trajectory1.push_back(state1);
		trajectory2.push_back(state2);

	}

	double lambda = sum / CT;

	return lambda;
}

std::vector<LLEPoint> build_LLE_diagram(double start, double end, int resolution,
	int coefIndex, int count, double CT, double NT,
	std::vector<double>& coefs, State start_state,
	State(*method)(const State&, std::vector<double>&, double)
) {
	std::vector<LLEPoint> diagram;
	double step = (end - start) / resolution;

#pragma omp parallel
	{

		std::vector<LLEPoint> localDiagram;

#pragma omp for

		for (int i = 0; i < resolution; i++)
		{
			if (i % 10 == 0) std::cout << "\nIndex number: " << i;

			double param = start + i * step;

			auto localCoefs = coefs;
			localCoefs[coefIndex] = param;

			double lambda = LLE_calculation(count, start_state, CT, NT, method, localCoefs);

			localDiagram.push_back({ param, lambda });
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

void save_lle_csv(const std::string& filename, const std::vector<LLEPoint>& data)
{
	std::ofstream file(filename);

	file << "parameter,lambda\n";
	std::cout << "Opened file...";

	for (const auto& p : data)
	{
		file << p.param
			<< ","
			<< p.lambda
			<< "\n";
	}
}