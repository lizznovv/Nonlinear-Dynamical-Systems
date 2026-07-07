
#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "System.h"
#include "Lyapunov.h"

double LLE_calculation(double h, State start_state, double CT, double NT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
) {
	State new_state;
	double eps = 1e-8;
	int count = static_cast<int>(CT / h);
	double dt = h;

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
            std::cout << "NaN! h = " << h
                << " a = " << coefs[0]
                << " b = " << coefs[1]
                << " c = " << coefs[2]
                << std::endl;
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
    int coefIndex, double CT, double NT, std::vector<double>& coefs,
    State start_state, State(*method)(const State&, std::vector<double>&, double),
    bool varyStep
)
{
    std::vector<LLEPoint> diagram;

    double linearStep = (end - start) / resolution;
    double k = 1.0;

    if (varyStep)
        k = pow(end / start, 1.0 / resolution);

#pragma omp parallel
    {
        std::vector<LLEPoint> localDiagram;

#pragma omp for
        for (int i = 0; i < resolution; i++)
        {
            if (i % 10 == 0)
                std::cout << "\nIndex number: " << i;

            auto localCoefs = coefs;

            double parameter;
            double lambda;
			double h;

            if (varyStep)
            {
                parameter = start * pow(k, i);
				h = parameter;
                lambda = LLE_calculation(h, start_state, CT, NT, method, localCoefs);
            }
            else
            {
                parameter = start + i * linearStep;
				localCoefs[coefIndex] = parameter;
				double h = 0.01;
                lambda = LLE_calculation(h, start_state, CT, NT, method, localCoefs);
            }

            localDiagram.push_back({ parameter, lambda });
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


std::vector<LLEMapPoint> ac_build_LLE_diagram(double a_start, double a_end,
    double c_start, double c_end, int resolution,
    int a_index, int c_index, double CT, double NT, std::vector<double>& coefs,
    State start_state, State(*method)(const State&, std::vector<double>&, double)
)
{
    std::vector<LLEMapPoint> diagram;

    double a_step = (a_end - a_start) / resolution;
    double c_step = (c_end - c_start) / resolution;

#pragma omp parallel
    {
        std::vector<LLEMapPoint> localDiagram;

#pragma omp for
        for (int i = 0; i < resolution; i++)
        {
            if (i % 10 == 0) std::cout << "\nIndex i number: " << i;

            double param1, param2;
            double lambda;

            param1 = a_start + i * a_step;
                
            for (int j = 0; j < resolution; j++) {

                if (j % 10 == 0) std::cout << "\nj = " << j;

                auto localCoefs = coefs;
                localCoefs[a_index] = param1;
                
                param2 = c_start + j * c_step;
                localCoefs[c_index] = param2;

                double h = 0.01;
                lambda = LLE_calculation(h, start_state, CT, NT, method, localCoefs);

                localDiagram.push_back({ param1, param2, lambda });
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

std::vector<LLEMapPoint> ah_build_LLE_diagram(double a_start, double a_end,
    double h_start, double h_end, int resolution, int a_index, double CT, double NT,
    std::vector<double>& coefs, State start_state,
    State(*method)(const State&, std::vector<double>&, double)
)
{
    std::vector<LLEMapPoint> diagram;

    double a_step = (a_end - a_start) / resolution;
    double k = pow(h_end / h_start, 1.0 / resolution);

#pragma omp parallel
    {
        std::vector<LLEMapPoint> localDiagram;

#pragma omp for
        for (int i = 0; i < resolution; i++)
        {
            double a = a_start + i * a_step;

            if (i % 10 == 0) std::cout << "\nIndex i number: " << i;

            for (int j = 0; j < resolution; j++)
            {
                if (j % 10 == 0) std::cout << "\nj = " << j;

                auto localCoefs = coefs;
                localCoefs[a_index] = a;
                double h = h_start * pow(k, j);

                double lambda = LLE_calculation(h, start_state, CT, NT, method, localCoefs);
                localDiagram.push_back({ a, h, lambda });
            }
        }

#pragma omp critical
        diagram.insert(diagram.end(),
            localDiagram.begin(),
            localDiagram.end());
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

void save_lle_map_csv(const std::string& filename,
    const std::vector<LLEMapPoint>& data)
{
    std::ofstream file(filename);

    file << "a,h,lambda\n";
    std::cout << "Opened file...";

    for (const auto& p : data)
    {
        file << p.param1
            << ","
            << p.param2
            << ","
            << p.lambda
            << "\n";
    }
}