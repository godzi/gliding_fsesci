#include <stdexcept>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string>
#include <iostream>
# include <cstdlib>
# include <iomanip>
# include <omp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>
#include <memory>
#include <immintrin.h>

#include <mkl_vsl.h>

#include "json.hpp"

#include "library.h"
#include "configuration.h"
#include "configuration_loader.h"
#include "mkl_gaussian_parallel_generator.h"
#include "mkl_flat_parallel_generator.h"
#include "binary_file_logger.h"
#include "dat_file_logger.h"

#include <fstream>

static constexpr unsigned nThreads = 5;

// Initialize global constants
//std::string inputparamfile = "C:\\Users\\Tolya\\Documents\\Visual Studio 2015\\Projects\\gliding_fsesci\\Release\\newconfig.json";
const double E = std::exp(1.0);
const double kBoltz= 1.38064852e-5;// (*pN um *)


/// ToDo try to use ofstream rawwrite




///
double mod(double a, double N)
{
	return a - N*floor(a / N); //return in range [0, N)
}

class KINESINvLoaded
{
public:
	double vUnloaded =0.0;
	double Fstall =0.0;
	double omega =0.0;
	double kinesinPparam = 0.0;
	double forceVelocityOn = 1.0;
	double calc(double Force) const
	{
		//
		//return vUnloaded*(1 - (pow((Force / Fstall) , omega)));
		return  forceVelocityOn*(vUnloaded / (kinesinPparam + ((1 - kinesinPparam)*exp(Force*2.79 / 4.11))))+(1.0- forceVelocityOn)*vUnloaded;
		
	}
};

class boundKINESIN
{
public:
	boundKINESIN(double coordinate, int MTsite, double springLength) :_mountCoordinate{ coordinate }, _MTsite{ MTsite }, _springLength{ springLength }
	{	};
	double _mountCoordinate;
	int _MTsite;
	double _springLength;
};
class unboundKINESIN
{
public:
	unboundKINESIN(double coordinate) {
		_mountCoordinate = coordinate;
	}
	double _mountCoordinate;
};
class boundMAP
{
public:
	boundMAP(double coordinate, int MTsite, double springLength) :_mountCoordinate{ coordinate }, _MTsite{ MTsite }, _springLength{ springLength }
	{	};
	double _mountCoordinate;
	int _MTsite;
	double _springLength;
};
class unboundMAP
{
public:
	unboundMAP(double coordinate) {
		_mountCoordinate = coordinate;
	}
	double _mountCoordinate;
};

class Task
{
public:
	Task(const SimulationParameters& simulationParameters, const Configuration& configuration):
		_sP( simulationParameters ),
		_mP( configuration.modelParameters ),
		_initC( configuration.initialConditions ),
		_state( configuration.initialConditions.initialState ),
		_kinesinStepsLog{ configuration.loggerParameters, "kinesinstepslog" },
		_MAPStepsLog{ configuration.loggerParameters, "mapstepslog" }
	{
		const auto& loggerParameters = configuration.loggerParameters;
		auto callback = [this, &loggerParameters](double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		};
		SystemState::iterateFields(callback);
		
		

	}
	int countTotalSteps() {
		return (int)ceil(_sP.totalTime / _sP.timeStep);
	}
	void saveStepingMap(double time, double proteinMountCoordinate,int MTsiteNum ) {
		if ((fabs(proteinMountCoordinate - _state.MonitorMAP)) <= _initC.MAPdistance / 2) {
			_MAPStepsLog.save(std::to_string(time) + "	" + std::to_string(proteinMountCoordinate) + "	" + std::to_string(MTsiteNum) + "\n");
		}
	}
	void saveStepingKinesin(double time, double proteinMountCoordinate, int MTsiteNum) {
		if ((fabs(proteinMountCoordinate - _state.Monitorkinesin)) <= _initC.KINESINdistance / 2) {
			_kinesinStepsLog.save(std::to_string(time) + "	" + std::to_string(proteinMountCoordinate) + "	" + std::to_string(MTsiteNum) + "\n");
		}
	}

	size_t initializeState() {
		_state.MTpositionStep = 0.0;
	//	_state.MTposition = 0.0;
		_NumberofMTsites = (int)floor(_initC.MTlength / _mP.deltaPeriod);
		
		
		_kinesinvLoaded.vUnloaded = _mP.vUnloaded;
		_kinesinvLoaded.Fstall = _mP.Fstall;
		_kinesinvLoaded.omega = _mP.omega;
		_kinesinvLoaded.kinesinPparam = _mP.kinesinPparam;
		_kinesinvLoaded.forceVelocityOn = _mP.forceVelocityOn;
		/// Initial binding
		for (double place = 0.0+ _initC.surfaceKINESINstartPoint; place < _initC.surfaceLength ; place = place + _initC.KINESINdistance) {
			_unboundKinesins.emplace_back(place);
		}
		for (double place = 0.0+ _initC.surfaceMAPstartPoint; place < _initC.surfaceLength ; place = place + _initC.MAPdistance) {
			_unboundMaps.emplace_back(place);
		}
	//	std::cout << "_unboundKinesins.size= " << _unboundKinesins.size() << std::endl;
	//	std::cout << "_unboundMaps.size= " << _unboundMaps.size() << std::endl;
		////////////////////////
		size_t neededFlatBufferSize = _unboundKinesins.size() + 2 * _unboundMaps.size();
		
		////
		/// test for binding
		for (int i = 1; i <= _NumberofMTsites;i++) {
			
			_SurfaceDistanceofiSite = _state.MTposition + _mP.deltaPeriod*(static_cast<double>(i) - 1.0);

			//std::cout << "_unboundMaps.size() = " << _unboundMaps.size() << std::endl;
			//std::cout << "iter->_mountCoordinate = " << _unboundMaps.begin()->_mountCoordinate << std::endl;
			for (auto iter = _unboundMaps.begin(); iter != _unboundMaps.end(); ) {
				// std::cout << "MAPS iter-> = " << iter- _unboundMaps.begin() << std::endl;
			//	std::cout << i << " "<< iter - _unboundMaps.begin() << " " << iter->_mountCoordinate << " " << _SurfaceDistanceofiSite << " " << ((fabs(iter->_mountCoordinate - _SurfaceDistanceofiSite)) <= (_mP.deltaPeriod / 2.0))  << std::endl;
				if ((fabs(iter->_mountCoordinate - _SurfaceDistanceofiSite)) <= 1.01*(_mP.deltaPeriod / 2.0)) {
					_boundMaps.emplace_back(iter->_mountCoordinate, i, iter->_mountCoordinate - _SurfaceDistanceofiSite);
					saveStepingMap(0.0, iter->_mountCoordinate,i);
					iter = _unboundMaps.erase(iter);
					
				}
				else
				{
					++iter;
				}					
			}
			for (auto iter = _unboundKinesins.begin(); iter != _unboundKinesins.end(); ) {
				if ((fabs(iter->_mountCoordinate - _SurfaceDistanceofiSite)) <= 1.01* (_mP.deltaPeriod / 2.0)) {
					_boundKinesins.emplace_back(iter->_mountCoordinate, i, iter->_mountCoordinate - _SurfaceDistanceofiSite);
					saveStepingKinesin(0.0, iter->_mountCoordinate, i);
					iter = _unboundKinesins.erase(iter);
					
				}
				else
				{
					++iter;
				}
			}
			
		}
	//	std::cout << "_boundKinesins.size() = " << _boundKinesins.size() << std::endl;
	//	std::cout << "_boundMaps.size() = " << _boundMaps.size() << std::endl;
		return neededFlatBufferSize;
	}
	double springForce(double stiffness, double length) {
		return stiffness*length;
		/*
		if (fabs(length) <= _mP.freeSpringLength) {
			return 0.0;
		}
		else {
			if (length >= 0.0) {
				return stiffness*(length - _mP.freeSpringLength);
			}
			if (length < 0.0) {
				return stiffness*(length + _mP.freeSpringLength);
			}
		}
		*/
		
	}
	//double taskStartTime,
	// rndNumbers must contain 3 * nSteps random numbers
	void advanceState(unsigned nSteps,  const double* rndNormalNumbers, const double* rndFlatNumbers, MklGaussianParallelGenerator* gaussGenerator, MklFlatParallelGenerator* flatGenerator) {

	//	std::cout << "MTgamma" << _mP.gammaMT << std::endl;
	//	std::cout << "kT" << _mP.kT << std::endl;

		int counterNormalRandomNumber = 0;
		int counterFlatRandomNumber = 0;
		
		// configurate force object
		


		//

		auto takeNormalRandomNumber = [rndNormalNumbers,&counterNormalRandomNumber,gaussGenerator,  this] () -> double {
			
			if (counterNormalRandomNumber == _optimalGaussianBuffer) {
				gaussGenerator->generateNumbers();
				//std::cout << ", NormalRandomCounter = " << counterNormalRandomNumber << std::endl;
				counterNormalRandomNumber = 0;

			}
			return rndNormalNumbers[counterNormalRandomNumber++];
		};
		auto takeFlatRandomNumber = [rndFlatNumbers, &counterFlatRandomNumber, flatGenerator, this]() -> double {
			//std::cout << counterFlatRandomNumber << " " << rndFlatNumbers[counterFlatRandomNumber] << std::endl;
			if (counterFlatRandomNumber == _optimalFlatBuffer) {
				flatGenerator->generateNumbers();
			//	std::cout << ", FlatRandomCounter = " << counterFlatRandomNumber << std::endl;
				counterFlatRandomNumber = 0;
			}
			return rndFlatNumbers[counterFlatRandomNumber++];
		};
		


		// Sites on MT to be checked for new bindings every iteration
		int sitestocheck[6] = { 1, 2, 3, _NumberofMTsites-2,_NumberofMTsites-1,_NumberofMTsites };
		///
		
		///
		/// Code below for single iteration
		double timeCompute = 0.0;
		_state.currentTime = 0.0;
		for (unsigned taskIteration = 0; taskIteration < nSteps; taskIteration++) {
			_state.currentTime += _sP.timeStep;
			_state.SummKINESINForces = 0.0;
			_state.SummMAPForces = 0.0;
			_SummForces = 0.0;
			
			
			if (taskIteration % 100000 == 0) {
				double procent = 100 * round(100000* (double)taskIteration / (double)nSteps) / 100000;
				std::cout << procent<< "%" << std::endl;
				
			//	std::cout << __rdtsc()-timeCompute << std::endl;
			//	timeCompute = __rdtsc();
				//std::cout << nst << std::endl;
			}
			
			
			//std::cout << "taskIteration = " << taskIteration << ", FlatRandomCounter = " << counterFlatRandomNumber << ", FlatRandom = " << rndFlatNumbers[counterFlatRandomNumber] << ", BoundMAPs = " << _boundMaps.size() <<std::endl;
			//takeFlatRandomNumber();

					/// Test for binding
				for (int it = 0; it < 6; it++) {
					int i = sitestocheck[it];
					_SurfaceDistanceofiSite = _state.MTposition + _mP.deltaPeriod*((double)i - 1.0);

					for (auto iter = _unboundMaps.begin(); iter != _unboundMaps.end(); ) {
						if ((fabs(iter->_mountCoordinate - _SurfaceDistanceofiSite)) <= 1.01*(_mP.deltaPeriod / 2.0)) {
							_boundMaps.emplace_back(iter->_mountCoordinate, i, iter->_mountCoordinate - _SurfaceDistanceofiSite);
							saveStepingMap(_state.currentTime, iter->_mountCoordinate, i);
							iter = _unboundMaps.erase(iter);
							//
							
						}
						else
						{
							++iter;
						}
					}
					for (auto iter = _unboundKinesins.begin(); iter != _unboundKinesins.end(); ) {
						if ((fabs(iter->_mountCoordinate - _SurfaceDistanceofiSite)) <= 1.01*(_mP.deltaPeriod / 2.0)) {
							_boundKinesins.emplace_back(iter->_mountCoordinate, i, iter->_mountCoordinate - _SurfaceDistanceofiSite);
							saveStepingKinesin(_state.currentTime, iter->_mountCoordinate, i);
							iter = _unboundKinesins.erase(iter);
							//
							
						}
						else
						{
							++iter;
						}
					}

				}
				/// Test for stepping for MT bound MAPs AND UPDATE SPRING EXTENSION LENGTHS
		
				for (auto iter = _boundMaps.begin(); iter != _boundMaps.end(); ) {
					
					// update spring extensions based on previous microtubule step
					iter->_springLength = iter->_springLength +  _state.MTpositionStep;
					//std::cout << "iter->_springLength*_mP.MAPstiffness " << iter->_springLength*_mP.MAPstiffness << std::endl;
					///
					_kprob = _mP.MAPsDiffusion / (_mP.deltaPeriod*_mP.deltaPeriod);
					_kpow = (-springForce(_mP.MAPstiffness, iter->_springLength))*_mP.MAPfsmPar / ( kBoltz*_mP.T);
					_kPlus = _kprob*exp(_kpow);
					_kMinus =_kprob*exp(-_kpow);
					
					if (takeFlatRandomNumber() > exp(-(_kMinus + _kPlus)*_sP.timeStep)) {
						//std::cout << "_state.MTpositionStep= " << _state.MTpositionStep << std::endl;
						
					// Make step
						if (takeFlatRandomNumber() > (_kPlus / (_kMinus + _kPlus))) {
						//	Make step to the left
							if (((iter->_MTsite) - 1) >= 1)
							{
							//test that we still have the site there
								iter->_MTsite = iter->_MTsite - 1;
								iter->_springLength = iter->_springLength - _mP.deltaPeriod;
								//
								saveStepingMap(_state.currentTime, iter->_mountCoordinate, iter->_MTsite);
								
								// Count new summ forces
								 _state.SummMAPForces =  _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
								// std::cout << " ((_kMinus + _kPlus)*_sP.timeStep)= " << ((_kMinus + _kPlus)*_sP.timeStep) << std::endl;
								++iter;
							}
							else
							{
								if(takeFlatRandomNumber()<=_mP.MAPunbindBorderFreq)
								{
									_unboundMaps.emplace_back(iter->_mountCoordinate);
									iter = _boundMaps.erase(iter);
								}
								else
								{
									_state.SummMAPForces = _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
									++iter;
								}

							}				
						}
						else
						{
						// Make step to the right
							if (((iter->_MTsite) + 1) <= _NumberofMTsites)
							{
								//test that we still have the site there
								iter->_MTsite = iter->_MTsite + 1;
								iter->_springLength = iter->_springLength + _mP.deltaPeriod;
								//
								saveStepingMap(_state.currentTime, iter->_mountCoordinate, iter->_MTsite);

								// Count new summ forces
								 _state.SummMAPForces =  _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
								++iter;
							}
							else
							{
								if (takeFlatRandomNumber() <= _mP.MAPunbindBorderFreq)
								{
									_unboundMaps.emplace_back(iter->_mountCoordinate);
									iter = _boundMaps.erase(iter);
								}
								else
								{
									_state.SummMAPForces = _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
									++iter;
								}
							}
						}
					
					}
					else
					{
						_state.SummMAPForces = _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
						++iter;
					}

				}
				////
				/// Test for stepping for MT bound MT kinesins AND UPDATE SPRING EXTENSION LENGTHS
				for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); ) {
					// update spring extensions based on previous microtubule step
					
					iter->_springLength = iter->_springLength +  _state.MTpositionStep;
					
					//std::cout << "iter->_springLength*_mP.KINESINstiffness " << iter->_springLength*_mP.KINESINstiffness << std::endl;
					//
					double pstep =exp(-_kinesinvLoaded.calc(fabs(springForce(_mP.KINESINstiffness, iter->_springLength)))*(_sP.timeStep / (2 * _mP.deltaPeriod)));
					if (takeFlatRandomNumber() > pstep)
					{
						//make step (always 2 site periods to the plus end of MT therefore to higher site number)
						if (iter->_MTsite + 2 <= _NumberofMTsites )
						{
							iter->_MTsite = iter->_MTsite + 2;
							iter->_springLength = iter->_springLength + 2* _mP.deltaPeriod;
							//
							
							
							saveStepingKinesin(_state.currentTime, iter->_mountCoordinate, iter->_MTsite);
							// Count new summ forces
							 _state.SummKINESINForces =  _state.SummKINESINForces + springForce(_mP.KINESINstiffness, iter->_springLength);
							++iter;
						}
						else
						{
							if (takeFlatRandomNumber() <= _mP.KinesinunbindBorderFreq)
							{
								_unboundKinesins.emplace_back(iter->_mountCoordinate);
								iter = _boundKinesins.erase(iter);
							}
							else
							{
								_state.SummKINESINForces = _state.SummKINESINForces + springForce(_mP.KINESINstiffness, iter->_springLength);
								++iter;
							}
						}
					}
					else {
						//
						 _state.SummKINESINForces =  _state.SummKINESINForces + springForce(_mP.KINESINstiffness, iter->_springLength);
						++iter;
					}
				}

		//////
			_SummForces = -(_mP.KINESINforcesOn* _state.SummKINESINForces+ _mP.MAPforcesOn* _state.SummMAPForces);
			
			////
			 _state.MTpositionStep= (_sP.timeStep / _mP.gammaMT)*_SummForces +	_mP.thermalNoiseOn*sqrt(2.0 * _mP.kT*_sP.timeStep / _mP.gammaMT) *	takeNormalRandomNumber();
			_state.MTposition = _state.MTposition +  _state.MTpositionStep;
			///
			_state.BoundedKinesins=_boundKinesins.size();
			_state.BoundedMAPs = _boundMaps.size();
			
			if (taskIteration %_sP.saveFrequency == 0) {
				writeStateTolog();
				
			}
		}
	}

	void writeStateTolog() const {
		for (const auto& logger : _loggers) {
			logger->save(&_state);
		}
	}

private:
	const SimulationParameters _sP;
	const ModelParameters _mP;
	const InitialConditions _initC;
	SystemState _state;
	std::vector<std::unique_ptr<BinaryFileLogger>> _loggers;
	
	//
	KINESINvLoaded _kinesinvLoaded;
	////
	int _NumberofMTsites;
	/*std::vector <double> _KINESINMountsites;
	std::vector <double> _MAPMountsites;
	std::vector <double> _MAPsitesonMT;
	std::vector <double> _KINESINsitesonMT;
	*/
	std::vector <unboundKINESIN> _unboundKinesins;
	std::vector <boundKINESIN> _boundKinesins;
	std::vector <unboundMAP> _unboundMaps;
	std::vector <boundMAP> _boundMaps;
	double _SurfaceDistanceofiSite;
	double _SummForces;
	double _kprob=0.0;
	double _kpow=0.0;
	double _kPlus=0.0;
	double _kMinus=0.0;
	DatFileLogger _kinesinStepsLog;
	DatFileLogger _MAPStepsLog;
public:
	int _testFlatBufferSizeFreq = 10;//iteration beetween tests
	int _testGaussBufferSizeFreq = 10000;
	size_t _optimalGaussianBuffer = 600000;
	size_t _optimalFlatBuffer = 600000;
};



int main(int argc, char *argv[])
{
	if (cmdOptionExists(argv, argv + argc, "-h"))
	{
		// Do stuff
		std::cout << "Sorry users, no help donations today." << std::endl;
	}
	char * param_input_filename = getCmdOption(argv, argv + argc, "-paramsfile");
	//char * output_filename = getCmdOption(argv, argv + argc, "-resultfile");
	
	std::string inputparamfile;
	inputparamfile.append(param_input_filename);
	//inputparamfile = "C:\\Users\\Tolya\\Documents\\Visual Studio 2015\\Projects\\gliding_fsesci\\x64\\Release\\NewConfig.json";
	
	// Create and load simulation parameters and configuration, values are taken from json file
	const auto simulationParameters = load_simulationparams(inputparamfile);
	const auto configurations = load_configuration(inputparamfile);

	std::vector<std::unique_ptr<Task>> tasks;
	for (const auto& configuration : configurations) {
		auto task = std::make_unique<Task>(simulationParameters, configuration);
		tasks.push_back(std::move(task));
	}
	
	// Check if number of configurations correspond to predefined threads number
	/*if (configurations.size() != nThreads) {
		throw std::runtime_error{ "Please check the number of configurations in json corresponds to the number of threads implemented in simulator" };
	}
	*/
	/*

	///////////////////////////////
	//////////////////////////////
	///////////// Iterations
	///////////////////////////
	/*
	const int frequency = (int)round(sP.microsteps);
	const int totalsimsteps = (int)round((sP.microsteps)*(sP.nTotal));
	const int nsteps = (int)round(3*(totalsimsteps / frequency) / intbufferSize);
		*/


		/*
		int frequency = sP.microsteps;
		int totalsimsteps = sP.microsteps*(sP.nTotal);
		int nsteps = 3 * (totalsimsteps / frequency) / intbufferSize;
		int dd = 2;

		frequency 100'000
		totalsimsteps 100'000'000'000
		buffer 1'000'000
		randms per simstep 3

		saved step range is 10'000 -> nTotal
		microsteps is macrostep*(buffersize/3)
		*/

	//size_t optimalFlatBufferSize = 900000;
	size_t flatbuffersizeperStep = tasks[0]->initializeState();
	int totalSteps = tasks[0]->countTotalSteps();
	int macrostepsNum = 1;
	//if (flatbuffersizeperStep*totalSteps > optimalFlatBufferSize) {
	//	macrostepsNum = (int)ceil(flatbuffersizeperStep*totalSteps / optimalFlatBufferSize);
	//}
	//int taskStepsNum = (int)floor(optimalFlatBufferSize/ flatbuffersizeperStep);
	int taskStepsNum = totalSteps / macrostepsNum;
	
	MklGaussianParallelGenerator gaussGenerator(0.0, 1.0, tasks[0]->_optimalGaussianBuffer, 1);
	MklFlatParallelGenerator flatGenerator;
	flatGenerator.initialize(0.0, 1.0, tasks[0]->_optimalFlatBuffer, 4);
	//std::cout << macrostepsNum << std::endl;
	//for (int savedstep = 0; savedstep < (100'000); savedstep++) {

		//for (int macrostep = 0; macrostep < macrostepsNum; macrostep++) {
			
			
			gaussGenerator.generateNumbers();
			const auto gaussBuffData = gaussGenerator.getNumbersBuffer();
			flatGenerator.generateNumbers();
			const auto flatBuffData = flatGenerator.getNumbersBuffer();

			tasks[0]->advanceState(taskStepsNum, gaussBuffData, flatBuffData,&gaussGenerator,&flatGenerator);
			
			/*
			#pragma omp parallel num_threads(nThreads) shared(buffData, tasks)
			{
				const auto begin = __rdtsc();
				tasks[omp_get_thread_num()]->advanceState(900'000 / 3, buffData);
				const auto end = __rdtsc();
				const auto cyclesPerStep = static_cast<double>(end - begin) / static_cast<double>(900'000 / 3);
				std::cout << "cyclesPerStep = " << cyclesPerStep << std::endl;
			} // end of openmp section

			for (const auto& task : tasks) {
				task->writeStateTolog();
			}
			*/
		
			/*if (macrostep % 10 == 0) {
				double procent = 100*round(1000*macrostep / macrostepsNum)/1000;
				std::cout << procent << "%" << std::endl;
				std::cout << __rdtsc() << std::endl;
				//std::cout << nst << std::endl;
			}*/
		//}

		
	//}

			std::cout << " end" << std::endl;
			std::cout << " end" << std::endl;
}

//std::chrono
//_rdtscp  in #immintrin.h

//		(*/),(+-)