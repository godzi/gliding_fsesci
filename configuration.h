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
	double kinesinForceUnbindingA;
	double kinesinForceUnbindingFd;
	double numberMAPsInOneSite;
	double numberKinesinsInOneSite;
	double useKinesinBindingKon;
	double useMAPBindingKon;
	double KinesinKon;
	double MAPKon;
	double MAPKoff;
	double MAPfsmParforKoff;
	double kinesinDparam;
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
	

	template <typename Callback>
	static void iterateFields(Callback&& callback) {
		callback(&SystemState::MTposition, "MTposition");
		callback(&SystemState::SummKINESINForces, "SummKINESINForces");
		callback(&SystemState::SummMAPForces, "SummMAPForces");
		callback(&SystemState::currentTime, "currentTime");
		callback(&SystemState::MTpositionStep, "MTpositionStep");
		callback(&SystemState::BoundedKinesins, "BoundedKinesins");
		callback(&SystemState::BoundedMAPs, "BoundedMAPs");
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
	double Monitorkinesin;
	double MonitorMAP;
	double watchMAPs;
	double watchKinesins;
	double fulllogMAPsKinesins;
	double kinesinUnbinding;
	double MAPUnbinding;
	double stopSimIfNoMAPsKinesins;
	double oldInitialMount;
	double numberKinesins;
	double numberMAPs;
};
// Composition of parameters
struct Configuration
{
	LoggerParameters loggerParameters;
	ModelParameters modelParameters;
	InitialConditions initialConditions;
	SystemState currentState;
};
