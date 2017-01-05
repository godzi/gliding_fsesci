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
	return simp;
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
	if (!(jsonobj["ModelParameters"]["MAPsDiffusion"].empty())) {
		conf.modelParameters.MAPsDiffusion = stod(jsonobj["ModelParameters"]["MAPsDiffusion"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["MAPstiffness"].empty())) {
		conf.modelParameters.MAPstiffness = stod(jsonobj["ModelParameters"]["MAPstiffness"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["KINESINstiffness"].empty())) {
		conf.modelParameters.KINESINstiffness = stod(jsonobj["ModelParameters"]["KINESINstiffness"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["vUnloaded"].empty())) {
		conf.modelParameters.vUnloaded = stod(jsonobj["ModelParameters"]["vUnloaded"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["omega"].empty())) {
		conf.modelParameters.omega = stod(jsonobj["ModelParameters"]["omega"].get<std::string>());
	}
	if (!(jsonobj["ModelParameters"]["Fstall"].empty())) {
		conf.modelParameters.Fstall = stod(jsonobj["ModelParameters"]["Fstall"].get<std::string>());
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
	//// Assign Dynamic Coordinates from json initial conditions
	if (!(jsonobj["InitialConditions"]["MTposition"].empty())) {
		conf.currentState.MTposition = stod(jsonobj["InitialConditions"]["MTposition"].get<std::string>());
	}
	
	////
	conf.modelParameters.gammaMT = 2 * 3.14159*conf.modelParameters.EtaMT*conf.initialConditions.MTlength / (log(4 * conf.modelParameters.hMT / conf.modelParameters.dMT));
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