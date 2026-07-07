#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "System.h"
#include "PhaseVolume.h"
#include "Bifurcation.h"

double phase_volume(double h, State start_state, double CT, double TT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs
) {
	std::vector<double> localCoefs = coefs;
	std::vector<State> trajectory = { start_state };

	run_method(h, TT, method, trajectory, localCoefs);
	trajectory = { trajectory.back() };

	double xmin = trajectory.back().x;
	double xmax = trajectory.back().x;

	double ymin = trajectory.back().y;
	double ymax = trajectory.back().y;

	double zmin = trajectory.back().z;
	double zmax = trajectory.back().z;


	run_method(h, CT, method, trajectory, localCoefs);

	for (const auto& state : trajectory){

		xmin = (state.x < xmin) ? state.x : xmin;
		xmax = (state.x > xmax) ? state.x : xmax;

		ymin = (state.y < ymin) ? state.y : ymin;
		ymax = (state.y > ymax) ? state.y : ymax;

		zmin = (state.z < zmin) ? state.z : zmin;
		zmax = (state.z > zmax) ? state.z : zmax;
	}

	double V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);

	return V;
}

std::vector<VolumePoint> build_volume_diagram(double start, double end, int resolution,
	int coefIndex, double CT, double NT, std::vector<double>& coefs,
	State start_state, State(*method)(const State&, std::vector<double>&, double),
	bool varyStep)
{
	std::vector<VolumePoint> diagram;

	double linearStep = (end - start) / resolution;
	double k = 1.0;

	if (varyStep)
		k = pow(end / start, 1.0 / resolution);

#pragma omp parallel
	{
		std::vector<VolumePoint> localDiagram;


#pragma omp for
		for (int i = 0; i < resolution; i++)
		{
			if (i % 10 == 0)
				std::cout << "\nIndex number: " << i;

			auto localCoefs = coefs;

			double parameter;
			double volume;
			double h;

			if (varyStep)
			{
				parameter = start * pow(k, i);
				h = parameter;
				volume = phase_volume(h, start_state, CT, NT, method, localCoefs);
			}
			else
			{
				parameter = start + i * linearStep;
				localCoefs[coefIndex] = parameter;
				h = 0.01;
				volume = phase_volume(h, start_state, CT, NT, method, localCoefs);
			}

			localDiagram.push_back({ parameter, volume });
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

void save_volume_csv(const std::string& filename,
	const std::vector<VolumePoint>& data)
{
	std::ofstream file(filename);

	file << "parameter,volume\n";
	std::cout << "Opened file...";

	for (const auto& p : data)
	{
		file << p.param
			<< ","
			<< p.volume
			<< "\n";
	}

	file.close();
}


std::vector<BifPoint> build_contunue_diagram(double h, double windowTime,
	State start_state, double CT, double TT,
	State(*method)(const State& state, std::vector<double>& coefs, double),
	std::vector<double>& coefs) 
{

	std::vector<BifPoint> diagram;
	
	std::vector<double> localCoefs = coefs;
	std::vector<State> trajectory = { start_state };

	// переходный процесс
	run_method(h, TT, method, trajectory, localCoefs);

	State current = trajectory.back();
	trajectory = { current };

	int windows = static_cast<int>(CT / windowTime);

	double time = 0;

	for (int i = 0; i < windows; i++)
	{
		run_method(h, windowTime, method, trajectory, localCoefs);

		current = trajectory.back();
		time += windowTime;

		auto peaks = findPeaks(trajectory);

		for (double peak : peaks)
			diagram.push_back({ time, peak });

		trajectory = { current };

	}


	return diagram;
}