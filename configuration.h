#pragma once
#include <string>
# include <omp.h>
#include <vector>
#include <string>

// Global parameters of simulations -> define iterations and steps
struct SimulationParameters
{
	// simulation parameters
	double totalTime;//
	double timeStep;//
};


// Classes for objects that store simulation configs: LoggerParameters, ModelParameters, InitialConditions
struct LoggerParameters
{
	enum class FilenameTemplate { PrefixName };
	FilenameTemplate filenametemplate;
	std::string filepath;
	std::string name;
};
struct ModelParameters
{
	double T ;
	double kT;
	double MAPsDiffusion ;
	double MAPstiffness ;
	double MAPforcesOn;
	double KINESINstiffness ;
	double KINESINforcesOn;
	double vUnloaded ;
	double omega ;
	double Fstall ;
	double kinesinPparam;
	double forceVelocityOn;
	double deltaPeriod ;
	double EtaMT;
	double dMT ;
	double hMT ;
	double gammaMT;
};

struct SystemState
{
	double MTposition;

	template <typename F>
	static void iterateFields(F&& f) {
		f(&SystemState::MTposition, "MTposition");
	}
};

struct InitialConditions
{
	SystemState initialState;
	double MAPdistance;
	double KINESINdistance;
	double surfaceLength;
	double MTlength;
	
	
};
// Composition of parameters
struct Configuration
{
	LoggerParameters loggerParameters;
	ModelParameters modelParameters;
	InitialConditions initialConditions;
	SystemState currentState;
};
