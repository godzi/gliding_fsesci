#include "configuration_loader.h"
#include "library.h"
#include <math.h>

using json = nlohmann::json;
const double E = std::exp(1.0);
const double kBoltz = 1.38064852e-5;// (*pN um *)



/// Assign Simulation Parameters
SimulationParameters assign_simulation_parameters_from_json(SimulationParameters simp, json jsonobjsimp) {
	if (!(jsonobjsimp["totalTime"].empty())) {
		simp.totalTime = stod(jsonobjsimp["totalTime"].get<std::string>());
	}
	if (!(jsonobjsimp["timeStep"].empty())) {
		simp.timeStep = stod(jsonobjsimp["timeStep"].get<std::string>());
	}
	if (!(jsonobjsimp["saveFrequency"].empty())) {
		simp.saveFrequency = std::stoi(jsonobjsimp["saveFrequency"].get<std::string>());
	}
	if (!(jsonobjsimp["useInitialSetup"].empty())) {
		simp.useInitialSetup = std::stoi(jsonobjsimp["useInitialSetup"].get<std::string>());
	}
	if (!(jsonobjsimp["initialSetupfile"].empty())) {
		simp.initialSetupfile = jsonobjsimp["initialSetupfile"].get<std::string>();
	}
	if (!(jsonobjsimp["stopSimIfNoMAPsKinesinsTmin"].empty())) {
		simp.stopSimIfNoMAPsKinesinsTmin = stod(jsonobjsimp["stopSimIfNoMAPsKinesinsTmin"].get<std::string>());
	}
	if (!(jsonobjsimp["extGaussianseed"].empty())) {
		simp.extGaussianseed = static_cast<unsigned>(std::stoul(jsonobjsimp["extGaussianseed"].get<std::string>()));
	}
	if (!(jsonobjsimp["extUniformseed"].empty())) {
		simp.extUniformseed = static_cast<unsigned>(std::stoul(jsonobjsimp["extUniformseed"].get<std::string>()));
	}
	if (!(jsonobjsimp["extMountseed"].empty())) {
		simp.extMountseed = static_cast<unsigned>(std::stoul(jsonobjsimp["extMountseed"].get<std::string>()));
	}
	if (!(jsonobjsimp["TimeMinLog"].empty())) {
		simp.TimeMinLog = stod(jsonobjsimp["TimeMinLog"].get<std::string>());
	}
	
	return simp;
}
// helper to process json value into double
double jsonValuetodouble(json obj)
{
	return stod(obj.get<std::string>());
}

InitialSetup assign_initial_setup_from_json(InitialSetup iSetup, json jsoniSetup)
{
	if (!(jsoniSetup["MTposition"].empty())) {
		iSetup.MTposition = stod(jsoniSetup["MTposition"].get<std::string>());
	}
	if (!(jsoniSetup["MTpositionStep"].empty())) {
		iSetup.MTpositionStep = stod(jsoniSetup["MTpositionStep"].get<std::string>());
	}
	if (!(jsoniSetup["boundKinesins"].empty())) {
		for (json::iterator jsonit = jsoniSetup["boundKinesins"].begin(); jsonit != jsoniSetup["boundKinesins"].end(); ++jsonit) {
			json jsonobj = *jsonit;
			iSetup.InitialBoundKINESINs.emplace_back(jsonValuetodouble(jsonobj["mountCoordinate"]), jsonValuetodouble(jsonobj["MTsite"]), jsonValuetodouble(jsonobj["springLength"]), jsonValuetodouble(jsonobj["id"]));
		}		
	}
	if (!(jsoniSetup["boundMaps"].empty())) {
		for (json::iterator jsonit = jsoniSetup["boundMaps"].begin(); jsonit != jsoniSetup["boundMaps"].end(); ++jsonit) {
			json jsonobj = *jsonit;
			iSetup.InitialBoundMAPs.emplace_back(jsonValuetodouble(jsonobj["mountCoordinate"]), jsonValuetodouble(jsonobj["MTsite"]), jsonValuetodouble(jsonobj["springLength"]), jsonValuetodouble(jsonobj["id"]));
		}
	}
	if (!(jsoniSetup["unboundKinesins"].empty())) {
		for (json::iterator jsonit = jsoniSetup["unboundKinesins"].begin(); jsonit != jsoniSetup["unboundKinesins"].end(); ++jsonit) {
			json jsonobj = *jsonit;
			iSetup.InitialUnboundKINESINs.emplace_back(jsonValuetodouble(jsonobj["mountCoordinate"]), jsonValuetodouble(jsonobj["id"]));
		}
	}
	if (!(jsoniSetup["unboundMaps"].empty())) {
		for (json::iterator jsonit = jsoniSetup["unboundMaps"].begin(); jsonit != jsoniSetup["unboundMaps"].end(); ++jsonit) {
			json jsonobj = *jsonit;
			iSetup.InitialUnboundMAPs.emplace_back(jsonValuetodouble(jsonobj["mountCoordinate"]), jsonValuetodouble(jsonobj["id"]));
		}
	}
	return iSetup;
}

// Assign configuration from json object
Configuration assign_config_from_json(Configuration conf, json jsonobj) {
	//// Assign Logger Parameters from json
	if (!(jsonobj["Name"].empty())) {
		conf.loggerParameters.name = jsonobj["Name"].get<std::string>();
	}
	if (!(jsonobj["LoggerParameters"]["FilePath"].empty())) {
		conf.loggerParameters.filepath= jsonobj["LoggerParameters"]["FilePath"].get<std::string>();
	}
	
	//// Assign Model Parameters from json
	if (!(jsonobj["ModelParameters"]["T"].empty())) {
		conf.modelParameters.T= stod(jsonobj["ModelParameters"]["T"].get<std::string>());
		conf.modelParameters.kT = kBoltz*conf.modelParameters.T;
	}
	if (!(jsonobj["ModelParameters"]["thermalNoiseOn"].empty())) {
		conf.modelParameters.thermalNoiseOn = stod(jsonobj["ModelParameters"]["thermalNoiseOn"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["freeSpringLength"].empty())) {
		conf.modelParameters.freeSpringLength = stod(jsonobj["ModelParameters"]["freeSpringLength"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["MAPsDiffusion"].empty())) {
		conf.modelParameters.MAPsDiffusion = stod(jsonobj["ModelParameters"]["MAPsDiffusion"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPstiffness"].empty())) {
		conf.modelParameters.MAPstiffness = stod(jsonobj["ModelParameters"]["MAPstiffness"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPforcesOn"].empty())) {
		conf.modelParameters.MAPforcesOn = stod(jsonobj["ModelParameters"]["MAPforcesOn"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["KINESINstiffness"].empty())) {
		conf.modelParameters.KINESINstiffness = stod(jsonobj["ModelParameters"]["KINESINstiffness"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["KINESINforcesOn"].empty())) {
		conf.modelParameters.KINESINforcesOn = stod(jsonobj["ModelParameters"]["KINESINforcesOn"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["vUnloaded"].empty())) {
		conf.modelParameters.vUnloaded = stod(jsonobj["ModelParameters"]["vUnloaded"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["forceVelocityOn"].empty())) {
		conf.modelParameters.forceVelocityOn = stod(jsonobj["ModelParameters"]["forceVelocityOn"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["omega"].empty())) {
		conf.modelParameters.omega = stod(jsonobj["ModelParameters"]["omega"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["Fstall"].empty())) {
		conf.modelParameters.Fstall = stod(jsonobj["ModelParameters"]["Fstall"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kinesinPparam"].empty())) {
		conf.modelParameters.kinesinPparam = stod(jsonobj["ModelParameters"]["kinesinPparam"].get<std::string>());
	}	
	if (!(jsonobj["ModelParameters"]["deltaPeriod"].empty())) {
		conf.modelParameters.deltaPeriod = stod(jsonobj["ModelParameters"]["deltaPeriod"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["EtaMT"].empty())) {
		conf.modelParameters.EtaMT =stod(jsonobj["ModelParameters"]["EtaMT"].get<std::string>()) ;
	}
	if (!(jsonobj["ModelParameters"]["dMT"].empty())) {
		conf.modelParameters.dMT = stod(jsonobj["ModelParameters"]["dMT"].get<std::string>()) ;
	}
	if (!(jsonobj["ModelParameters"]["hMT"].empty())) {
		conf.modelParameters.hMT = stod(jsonobj["ModelParameters"]["hMT"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPfsmPar"].empty())) {
		conf.modelParameters.MAPfsmPar = stod(jsonobj["ModelParameters"]["MAPfsmPar"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPunbindBorderFreq"].empty())) {
		conf.modelParameters.MAPunbindBorderFreq = stod(jsonobj["ModelParameters"]["MAPunbindBorderFreq"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["KinesinunbindBorderFreq"].empty())) {
		conf.modelParameters.KinesinunbindBorderFreq = stod(jsonobj["ModelParameters"]["KinesinunbindBorderFreq"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kinesinForceUnbindingFd"].empty())) {
		conf.modelParameters.kinesinForceUnbindingFd = stod(jsonobj["ModelParameters"]["kinesinForceUnbindingFd"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kinesinForceUnbindingA"].empty())) {
		conf.modelParameters.kinesinForceUnbindingA = stod(jsonobj["ModelParameters"]["kinesinForceUnbindingA"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["numberMAPsInOneSite"].empty())) {
		conf.modelParameters.numberMAPsInOneSite = stod(jsonobj["ModelParameters"]["numberMAPsInOneSite"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["numberKinesinsInOneSite"].empty())) {
		conf.modelParameters.numberKinesinsInOneSite = stod(jsonobj["ModelParameters"]["numberKinesinsInOneSite"].get<std::string>());
	}

	

	if (!(jsonobj["ModelParameters"]["useKinesinBindingKon"].empty())) {
		conf.modelParameters.useKinesinBindingKon = stod(jsonobj["ModelParameters"]["useKinesinBindingKon"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["useMAPBindingKon"].empty())) {
		conf.modelParameters.useMAPBindingKon = stod(jsonobj["ModelParameters"]["useMAPBindingKon"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["KinesinKon"].empty())) {
		conf.modelParameters.KinesinKon = stod(jsonobj["ModelParameters"]["KinesinKon"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPKon"].empty())) {
		conf.modelParameters.MAPKon = stod(jsonobj["ModelParameters"]["MAPKon"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPKoff"].empty())) {
		conf.modelParameters.MAPKoff = stod(jsonobj["ModelParameters"]["MAPKoff"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["MAPfsmParforKoff"].empty())) {
		conf.modelParameters.MAPfsmParforKoff = stod(jsonobj["ModelParameters"]["MAPfsmParforKoff"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["kinesinDparam"].empty())) {
		conf.modelParameters.kinesinDparam = stod(jsonobj["ModelParameters"]["kinesinDparam"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kinesinOneDet"].empty())) {
		conf.modelParameters.kinesinOneDet = stod(jsonobj["ModelParameters"]["kinesinOneDet"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kinesinOneExpfsp"].empty())) {
		conf.modelParameters.kinesinOneExpfsp = stod(jsonobj["ModelParameters"]["kinesinOneExpfsp"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["kinesinOneLinfsp"].empty())) {
		conf.modelParameters.kinesinOneLinfsp = stod(jsonobj["ModelParameters"]["kinesinOneLinfsp"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["constantForceonMT"].empty())) {
		conf.modelParameters.constantForceonMT = stod(jsonobj["ModelParameters"]["constantForceonMT"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["MTgrowthrate"].empty())) {
		conf.modelParameters.MTgrowthrate = stod(jsonobj["ModelParameters"]["MTgrowthrate"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MTshrinkrate"].empty())) {
		conf.modelParameters.MTshrinkrate = stod(jsonobj["ModelParameters"]["MTshrinkrate"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["dynamicMT"].empty())) {
		conf.modelParameters.dynamicMT = stod(jsonobj["ModelParameters"]["dynamicMT"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["useMAPunbindingThreshold"].empty())) {
		conf.modelParameters.useMAPunbindingThreshold = stod(jsonobj["ModelParameters"]["useMAPunbindingThreshold"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPThresholdRight"].empty())) {
		conf.modelParameters.MAPThresholdRight = stod(jsonobj["ModelParameters"]["MAPThresholdRight"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPThresholdLeft"].empty())) {
		conf.modelParameters.MAPThresholdLeft = stod(jsonobj["ModelParameters"]["MAPThresholdLeft"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPunbindingThresholdisExponential"].empty())) {
		conf.modelParameters.MAPunbindingThresholdisExponential = stod(jsonobj["ModelParameters"]["MAPunbindingThresholdisExponential"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPAssymDiffusion"].empty())) {
		conf.modelParameters.MAPAssymDiffusion = stod(jsonobj["ModelParameters"]["MAPAssymDiffusion"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["MAPDiffSmRight"].empty())) {
		conf.modelParameters.MAPDiffSmRight = stod(jsonobj["ModelParameters"]["MAPDiffSmRight"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPDiffSmLeft"].empty())) {
		conf.modelParameters.MAPDiffSmLeft = stod(jsonobj["ModelParameters"]["MAPDiffSmLeft"].get<std::string>());
	}
	
	if (!(jsonobj["ModelParameters"]["extensionCriticalLength"].empty())) {
		conf.modelParameters.extensionCriticalLength = stod(jsonobj["ModelParameters"]["extensionCriticalLength"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["criticalStiffness"].empty())) {
		conf.modelParameters.criticalStiffness = stod(jsonobj["ModelParameters"]["criticalStiffness"].get<std::string>());
	}
	
	//// Assign Initial Conditions from json
	if (!(jsonobj["InitialConditions"]["MAPdistance"].empty())) {
		conf.initialConditions.MAPdistance = stod(jsonobj["InitialConditions"]["MAPdistance"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["KINESINdistance"].empty())) {
		conf.initialConditions.KINESINdistance = stod(jsonobj["InitialConditions"]["KINESINdistance"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["surfaceLength"].empty())) {
		conf.initialConditions.surfaceLength = stod(jsonobj["InitialConditions"]["surfaceLength"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["MTlength"].empty())) {
		conf.initialConditions.MTlength = stod(jsonobj["InitialConditions"]["MTlength"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["MTposition"].empty())) {
		conf.initialConditions.initialState.MTposition = stod(jsonobj["InitialConditions"]["MTposition"].get<std::string>());
	}	
	if (!(jsonobj["InitialConditions"]["surfaceMAPstartPoint"].empty())) {
		conf.initialConditions.surfaceMAPstartPoint = stod(jsonobj["InitialConditions"]["surfaceMAPstartPoint"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["surfaceKINESINstartPoint"].empty())) {
		conf.initialConditions.surfaceKINESINstartPoint= stod(jsonobj["InitialConditions"]["surfaceKINESINstartPoint"].get<std::string>());
	}

	if (!(jsonobj["InitialConditions"]["Monitorkinesin"].empty())) {
		conf.initialConditions.Monitorkinesin = stod(jsonobj["InitialConditions"]["Monitorkinesin"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["MonitorMAP"].empty())) {
		conf.initialConditions.MonitorMAP = stod(jsonobj["InitialConditions"]["MonitorMAP"].get<std::string>());
	}
	
	if (!(jsonobj["InitialConditions"]["watchMAPs"].empty())) {
		conf.initialConditions.watchMAPs = stod(jsonobj["InitialConditions"]["watchMAPs"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["watchMAPs"].empty())) {
		conf.initialConditions.watchMAPs = stod(jsonobj["InitialConditions"]["watchMAPs"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["watchKinesins"].empty())) {
		conf.initialConditions.watchKinesins = stod(jsonobj["InitialConditions"]["watchKinesins"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["fulllogMAPsKinesins"].empty())) {
		conf.initialConditions.fulllogMAPsKinesins = stod(jsonobj["InitialConditions"]["fulllogMAPsKinesins"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["kinesinUnbinding"].empty())) {
		conf.initialConditions.kinesinUnbinding = stod(jsonobj["InitialConditions"]["kinesinUnbinding"].get<std::string>());
	}
	
	if (!(jsonobj["InitialConditions"]["MAPUnbinding"].empty())) {
		conf.initialConditions.MAPUnbinding = stod(jsonobj["InitialConditions"]["MAPUnbinding"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["stopSimIfNoMAPsKinesins"].empty())) {
		conf.initialConditions.stopSimIfNoMAPsKinesins = stod(jsonobj["InitialConditions"]["stopSimIfNoMAPsKinesins"].get<std::string>());
	}
	
	if (!(jsonobj["InitialConditions"]["oldInitialMount"].empty())) {
		conf.initialConditions.oldInitialMount = stod(jsonobj["InitialConditions"]["oldInitialMount"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["numberKinesins"].empty())) {
		conf.initialConditions.numberKinesins = stod(jsonobj["InitialConditions"]["numberKinesins"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["numberMAPs"].empty())) {
		conf.initialConditions.numberMAPs = stod(jsonobj["InitialConditions"]["numberMAPs"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["useKinesinOneparams"].empty())) {
		conf.initialConditions.useKinesinOneparams = stod(jsonobj["InitialConditions"]["useKinesinOneparams"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["watchMAPsCircularBuffer"].empty())) {
		conf.initialConditions.watchMAPsCircularBuffer = stod(jsonobj["InitialConditions"]["watchMAPsCircularBuffer"].get<std::string>());
	}
	if (!(jsonobj["InitialConditions"]["watchKinesinsCircularBuffer"].empty())) {
		conf.initialConditions.watchKinesinsCircularBuffer = stod(jsonobj["InitialConditions"]["watchKinesinsCircularBuffer"].get<std::string>());
	}

	//// Assign Dynamic Coordinates from json initial conditions
	if (!(jsonobj["InitialConditions"]["MTposition"].empty())) {
		conf.currentState.MTposition = stod(jsonobj["InitialConditions"]["MTposition"].get<std::string>());
	}
	
	////
	//conf.modelParameters.gammaMT = 2 * 3.14159*conf.modelParameters.EtaMT*conf.initialConditions.MTlength / (log(4 * conf.modelParameters.hMT / conf.modelParameters.dMT));
	conf.modelParameters.gammaMT = 2.0 * 3.14159*conf.modelParameters.EtaMT*4.0 / (log(4 * conf.modelParameters.hMT / conf.modelParameters.dMT));
	///
	return conf;
}


////Simulation params loader
SimulationParameters load_simulationparams(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json jsonsimp = fulljson["SimulationParameters"];
	SimulationParameters simp = {};
	return assign_simulation_parameters_from_json(simp, jsonsimp);
	 
}

//// Configuration creator
std::vector <Configuration> load_configuration(std::string paramInputFilename) {
	json fulljson = parse_json_string(readfile(paramInputFilename));
	json defaultjson = fulljson["Configuration"];
	Configuration default;
	default = assign_config_from_json(default, defaultjson);



	std::vector <Configuration> configurationsVector;
	Configuration iterate;

	// iterate the array of configurations
	for (json::iterator jsonit = defaultjson["Configurations"].begin(); jsonit != defaultjson["Configurations"].end(); ++jsonit) {
		iterate = default;
		iterate=assign_config_from_json(iterate, *jsonit);
		configurationsVector.push_back(iterate);
	}

	return configurationsVector;

}

//// Initial Setup
InitialSetup load_setup(std::string setupFilename,int useInitialSetup) {
	if (useInitialSetup == 1)
	{
		json fulljson = parse_json_string(readfile(setupFilename));
		//json jsonsimp = fulljson["SimulationParameters"];
		InitialSetup iSetup;
		return assign_initial_setup_from_json(iSetup, fulljson);
		//return iSetup;
	}
	else
	{
		InitialSetup iSetup;
		return iSetup;
	}
}