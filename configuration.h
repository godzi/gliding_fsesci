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
	int saveFrequency;//save evry saveFrequency step
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
	double thermalNoiseOn;
	double freeSpringLength;
	double MAPsDiffusion ;
	double MAPstiffness ;
	double MAPforcesOn;
	double MAPfsmPar;
	double MAPunbindBorderFreq;
	double KinesinunbindBorderFreq;
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
	double SummKINESINForces;
	double SummMAPForces;
	double currentTime;
	double MTpositionStep;
	double BoundedKinesins;
	double BoundedMAPs;
	double Monitorkinesin;
	double MonitorMAP;

	template <typename F>
	static void iterateFields(F&& f) {
		f(&SystemState::MTposition, "MTposition");
		f(&SystemState::SummKINESINForces, "SummKINESINForces");
		f(&SystemState::SummMAPForces, "SummMAPForces");
		f(&SystemState::currentTime, "currentTime");
		f(&SystemState::MTpositionStep, "MTpositionStep");
		f(&SystemState::BoundedKinesins, "BoundedKinesins");
		f(&SystemState::BoundedMAPs, "BoundedMAPs");
		f(&SystemState::Monitorkinesin, "Monitorkinesin");
		f(&SystemState::MonitorMAP, "MonitorMAP");
	}
};

struct InitialConditions
{
	SystemState initialState;
	double MAPdistance;
	double KINESINdistance;
	double surfaceLength;
	double MTlength;
	double surfaceKINESINstartPoint;
	double surfaceMAPstartPoint;
	
};
// Composition of parameters
struct Configuration
{
	LoggerParameters loggerParameters;
	ModelParameters modelParameters;
	InitialConditions initialConditions;
	SystemState currentState;
};
