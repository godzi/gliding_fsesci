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

#include <vector>
#include <deque>
#include <numeric> //std::accumulate

#include <algorithm>
#include <math.h>
#include <iomanip>
#include <map>
#include <random>
#include <chrono>
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
#include "simulation_agents.h"

#include <fstream>

static constexpr unsigned nThreads = 5;

// Initialize global constants
//std::string inputparamfile = "C:\\Users\\Tolya\\Documents\\Visual Studio 2015\\Projects\\gliding_fsesci\\Release\\newconfig.json";
const double E = std::exp(1.0);
const double kBoltz= 1.38064852e-5;// (* pN um *)


/// ToDo try to use ofstream rawwrite


double get_uniform_rnd_forMounts(double start, double end, std::default_random_engine& generatorMount)
{
	std::uniform_real_distribution<double> distribution(start, end);
	return distribution(generatorMount);

}



///
double mod(double a, double N)
{
	return a - N*floor(a / N); //return in range [0, N)
}
#pragma pack(push,1)
struct SingleProteinSteplog
{
	double time;
	double mountCoordinate;
	int MTsiteNum;
	double MTposition;
	int proteinstate;
};
#pragma pack(pop)
static_assert(sizeof(SingleProteinSteplog) ==32, "wrong SingleProteinSteplog packing");

class Task
{
public:
	Task(const SimulationParameters& simulationParameters, const Configuration& configuration, const InitialSetup& initSetup):
		_sP( simulationParameters ),
		_mP( configuration.modelParameters ),
		_initC( configuration.initialConditions ),
		_iSetup(initSetup),
		_state( configuration.initialConditions.initialState ),
		_loggerparams(configuration.loggerParameters),
		_kinesinStepsLog{ configuration.loggerParameters, "kinesinstepslog" },
		_MAPStepsLog{ configuration.loggerParameters, "mapstepslog" },
		//benchmark
		//_performanceLog{ _loggerparams, "performancelog" },
		//benchmark
		_kinesinPosLog{ configuration.loggerParameters, "kinposlog" },
		_MAPPosLog{ configuration.loggerParameters, "mapposlog" }
	{
		const auto& loggerParameters = configuration.loggerParameters;
		auto callback = [this, &loggerParameters](double(SystemState::* field), std::string fieldName) {
			auto logger = std::make_unique<BinaryFileLogger>(loggerParameters, field, fieldName);// creates object of class BinaryFileLogger but returns to logger variable the unique pointer to it. // for creation of object implicitly with arguments like this also remind yourself the vector.emplace(args) .
			this->_loggers.push_back(std::move(logger));// unique pointer can't be copied, only moved like this
		};
		SystemState::iterateFields(callback);
		
		

	}
	int countTotalSteps() {
		std::cout << "total steps" << (int)ceil(_sP.totalTime / _sP.timeStep) << std::endl;
		return (int)ceil(_sP.totalTime / _sP.timeStep);

	}


	void saveMapPosition(std::string positions) {
		_MAPPosLog.save(positions);
	}
	void saveKinesinPosition(std::string positions) {
		_kinesinPosLog.save(positions);
	}


	void saveStepingMap(double time, double proteinMountCoordinate,int MTsiteNum,int mstate ) {
	//	std::cout << proteinMountCoordinate << " " << _initC.MonitorMAP << std::endl;
		if ((_initC.watchMAPs == 1.0)&&(_initC.watchMAPsCircularBuffer==1.0)) {
			SingleProteinSteplog thisStepLog;
			thisStepLog.time= time;
			thisStepLog.mountCoordinate= proteinMountCoordinate;
			thisStepLog.MTsiteNum= MTsiteNum;
			thisStepLog.MTposition=_state.MTposition;
			thisStepLog.proteinstate = mstate;

			_MapBuffer.push_back(thisStepLog);
			if (_MapBuffer.size() > _saveStepingMapBufferSize)
			{
				_MapBuffer.pop_front();
			}
						
		}
	}
	void saveStepingKinesin(double time, double proteinMountCoordinate, int MTsiteNum, int kstate) {
		
			if ((_initC.watchKinesins == 1.0) && (_initC.watchKinesinsCircularBuffer == 1.0)) {
				SingleProteinSteplog thisStepLog;
				thisStepLog.time = time;
				thisStepLog.mountCoordinate = proteinMountCoordinate;
				thisStepLog.MTsiteNum = MTsiteNum;
				thisStepLog.MTposition = _state.MTposition;
				thisStepLog.proteinstate = kstate;

				_KinesinBuffer.push_back(thisStepLog);
				if (_KinesinBuffer.size() > _saveStepingMapBufferSize)
				{
					_KinesinBuffer.pop_front();
				}

			}		
	}
	
	




	size_t initializeState() {
		//Random for surface Mounts of MAPs and kinesins
		if (_sP.extMountseed==0)
		{
			_seedMount = std::chrono::system_clock::now().time_since_epoch().count();
		}
		else
		{
			_seedMount = _sP.extMountseed;
		}
		std::cout << "mount seed " << _seedMount << std::endl;
		std::default_random_engine generatorMount(_seedMount);
		///

		_state.currentTime = 0.0;
		_state.MTpositionStep = 0.0;
		_state.MTcurLength = _initC.MTlength;
		_intialMTposition = _state.MTposition;
	//	_state.MTposition = 0.0;
		_NumberofMTsites = (int)floor(_initC.MTlength / _mP.deltaPeriod);
		
		
		// loggers for initial binding

		
		DatFileLogger* kinesinCoordsLog;
		DatFileLogger* MAPCoordsLog;
			
		kinesinCoordsLog=new DatFileLogger(_loggerparams, "kinecoordslog"); // log of initial mounting
		MAPCoordsLog = new DatFileLogger(_loggerparams, "mapcoordslog"); //log of initial mounting
		
		
		/////
		/////
		if (_sP.useInitialSetup == 0)
		{
			/// Initial binding

			
			double mountcoord;
			//Kinesin binding

			for (double ind = 0.0; ind < _initC.numberKinesins; ind++)
			{
				mountcoord = get_uniform_rnd_forMounts(_initC.surfaceKINESINstartPoint, _initC.surfaceLength, generatorMount);
				_unboundKinesins.emplace_back(mountcoord, (int)ind);
				kinesinCoordsLog->save(std::to_string((int)ind) + "	" + std::to_string(mountcoord) + "	");
			}
			//MAPs binding
			for (double ind = 0.0; ind < _initC.numberMAPs; ind++)
			{
				mountcoord = get_uniform_rnd_forMounts(_initC.surfaceMAPstartPoint, _initC.surfaceLength, generatorMount);
				_unboundMaps.emplace_back(mountcoord, (int)ind);
				MAPCoordsLog->save(std::to_string((int)ind) + "	" + std::to_string(mountcoord) + "	");
			}
			

		}
		if (_sP.useInitialSetup == 1)
		{
			_state.MTposition = _iSetup.MTposition;
			_state.MTpositionStep = _iSetup.MTpositionStep;
			 _unboundKinesins= _iSetup.InitialUnboundKINESINs;
			 _boundKinesins = _iSetup.InitialBoundKINESINs;
			 _unboundMaps = _iSetup.InitialUnboundMAPs;
			 _boundMaps = _iSetup.InitialBoundMAPs;
		}



		delete MAPCoordsLog;
		delete kinesinCoordsLog;

		
	//	std::cout << "_unboundKinesins.size= " << _unboundKinesins.size() << std::endl;
	//	std::cout << "_unboundMaps.size= " << _unboundMaps.size() << std::endl;
		////////////////////////
		size_t neededFlatBufferSize = _unboundKinesins.size() + 2 * _unboundMaps.size();
		
		////
		
		
		

	//	std::cout << "_boundKinesins.size() = " << _boundKinesins.size() << std::endl;
	//	std::cout << "_boundMaps.size() = " << _boundMaps.size() << std::endl;
		return neededFlatBufferSize;
	}
	double springForce(double stiffness, double length) {
	//	return stiffness*length;
				
		if (fabs(length) <= _mP.extensionCriticalLength) {
			return stiffness*length;
		}
		else {
			if (length >= 0.0) {
				return (stiffness*length+_mP.criticalStiffness*(length - _mP.extensionCriticalLength));
			}
			if (length < 0.0) {
				return (stiffness*length + _mP.criticalStiffness*(length + _mP.extensionCriticalLength));
			}
		}
		
		
	}
	//double taskStartTime,
	// rndNumbers must contain 3 * nSteps random numbers
	std::string advanceState(unsigned nSteps,  const double* rndNormalNumbers, const double* rndFlatNumbers, MklGaussianParallelGenerator* gaussGenerator, MklFlatParallelGenerator* flatGenerator) {
		_totalsteps = nSteps;

	//	std::cout << "MTgamma" << _mP.gammaMT << std::endl;
	//	std::cout << "kT" << _mP.kT << std::endl;

		int counterNormalRandomNumber = 0;
		int counterFlatRandomNumber = 0;
		
		// configurate force object
		


		//

		auto takeNormalRandomNumber = [rndNormalNumbers,&counterNormalRandomNumber,gaussGenerator,  this] () -> double {
			
			if (counterNormalRandomNumber == _optimalGaussianBuffer-1) {
				gaussGenerator->generateNumbers();
				//std::cout << ", NormalRandomCounter = " << counterNormalRandomNumber << std::endl;
				counterNormalRandomNumber = 0;

			}
			return rndNormalNumbers[counterNormalRandomNumber++];
		};
		auto takeFlatRandomNumber = [rndFlatNumbers, &counterFlatRandomNumber, flatGenerator, this]() -> double {
			//std::cout << counterFlatRandomNumber << " " << rndFlatNumbers[counterFlatRandomNumber] << std::endl;
			if (counterFlatRandomNumber == _optimalFlatBuffer-1) {
				flatGenerator->generateNumbers();
			//	std::cout << ", FlatRandomCounter = " << counterFlatRandomNumber << std::endl;
				counterFlatRandomNumber = 0;
			}
			return rndFlatNumbers[counterFlatRandomNumber++];
		};
		


		
		///
		/// Code below for single iteration
		double timeCompute = 0.0;
		_state.currentTime = 0.0;
		_everboundedMAPsKinesins = 0.0;
		//dump initial state to json
		dumpFullStatetoJson("_initial_state_dump");

		for (unsigned taskIteration = 0; taskIteration < nSteps; taskIteration++) {

			// benchmark
			//auto begin = std::chrono::high_resolution_clock::now();
			// benchmark

			_currentstep = taskIteration;
			_state.currentTime += _sP.timeStep;
			_state.SummKINESINForces = 0.0;
			_state.SummMAPForces = 0.0;
			_SummForces = 0.0;


			if (taskIteration % 100000 == 0) {
				double procent = 100 * round(100000 * (double)taskIteration / (double)nSteps) / 100000;
				std::cout << procent << "%" << std::endl;

				//	std::cout << __rdtsc()-timeCompute << std::endl;
				//	timeCompute = __rdtsc();
					//std::cout << nst << std::endl;
			}


			//std::cout << "taskIteration = " << taskIteration << ", FlatRandomCounter = " << counterFlatRandomNumber << ", FlatRandom = " << rndFlatNumbers[counterFlatRandomNumber] << ", BoundMAPs = " << _boundMaps.size() <<std::endl;
			//takeFlatRandomNumber();

			/// Test for MT growth and shrinkage
			if (_mP.dynamicMT == 1.0)
			{
				if ((takeFlatRandomNumber() > exp(-(_mP.MTgrowthrate)*_sP.timeStep)))
				{
					_state.MTcurLength = _state.MTcurLength + 0.004;
					_NumberofMTsites = _NumberofMTsites + 1;
				}
				if ((takeFlatRandomNumber() > exp(-(_mP.MTshrinkrate)*_sP.timeStep)) && (_NumberofMTsites > 2))
				{
					_state.MTcurLength = _state.MTcurLength - 0.004;
					_NumberofMTsites = _NumberofMTsites - 1;
				}
			}


			/// test for binding
			
			for (auto iter = _unboundMaps.begin(); iter != _unboundMaps.end(); ) {
				_MTsystemCoordinate = iter->_mountCoordinate - _state.MTposition;
				if ((takeFlatRandomNumber() > exp(-(_mP.MAPKon)*_sP.timeStep)))
				{
					if ((_MTsystemCoordinate < -0.5*_mP.deltaPeriod) || (_MTsystemCoordinate > (double(_NumberofMTsites) - 0.5)*_mP.deltaPeriod))
					{
						++iter;
					}
					else
					{
						_fractpart = modf(((_MTsystemCoordinate + 0.5*_mP.deltaPeriod) / _mP.deltaPeriod), &_intpart);
						_sitenum = _intpart + 1.0;

						if ((_sitenum > 1.0) && (_sitenum < double(_NumberofMTsites)) && ((0.5 - fabs((_MTsystemCoordinate - (_mP.deltaPeriod*(_sitenum - 1.0))) / _mP.deltaPeriod)) < 0.01))
						{

							if (takeFlatRandomNumber() <= 0.5)
							{
								_sitenum = _sitenum + 1.0;
							}
						}
						_SurfaceDistanceofiSite = _state.MTposition + _mP.deltaPeriod*(_sitenum - 1.0);
						_boundMaps.emplace_back(iter->_mountCoordinate, static_cast<int>(_sitenum), _SurfaceDistanceofiSite - iter->_mountCoordinate, iter->_id);
						saveStepingMap(_state.currentTime, iter->_mountCoordinate, static_cast<int>(_sitenum), 1);
						_everboundedMAPsKinesins = _everboundedMAPsKinesins + 1.0;
						if (_everboundedMAPsKinesins==1.0)
						{
							std::cout << "first protein bound at " << std::to_string(_state.currentTime) << " time" << std::endl;
						}
						iter = _unboundMaps.erase(iter);

						//
					}
				}
				else
				{
					++iter;
				}
			}
		
			// Test for kinesin binding
			
				for (auto iter = _unboundKinesins.begin(); iter != _unboundKinesins.end(); ) {
					_MTsystemCoordinate = iter->_mountCoordinate - _state.MTposition;


					if ((takeFlatRandomNumber() > exp(-(_mP.KinesinKon)*_sP.timeStep)))
					{

						if ((_MTsystemCoordinate < -0.5*_mP.deltaPeriod) || (_MTsystemCoordinate >(double(_NumberofMTsites) + 1)*_mP.deltaPeriod))
						{
							++iter;
						}
						else
						{
							_fractpart = modf(((_MTsystemCoordinate + 0.5*_mP.deltaPeriod) / _mP.deltaPeriod), &_intpart);
							_sitenum = _intpart + 1.0;

							if ((_sitenum > 1.0) && (_sitenum < double(_NumberofMTsites)) && ((0.5 - fabs((_MTsystemCoordinate - (_mP.deltaPeriod*(_sitenum - 1.0))) / _mP.deltaPeriod)) < 0.01))
							{

								if (takeFlatRandomNumber() <= 0.5)
								{
									_sitenum = _sitenum + 1.0;
								}
							}


							//correct if binding is away MT
							if (_sitenum > double(_NumberofMTsites))
							{
								_sitenum = _NumberofMTsites;
							}

							_SurfaceDistanceofiSite = _state.MTposition + _mP.deltaPeriod*(_sitenum - 1.0);
							_boundKinesins.emplace_back(iter->_mountCoordinate, static_cast<int>(_sitenum), _SurfaceDistanceofiSite - iter->_mountCoordinate, iter->_id);
							saveStepingKinesin(_state.currentTime, iter->_mountCoordinate, static_cast<int>(_sitenum), 1);

							_everboundedMAPsKinesins = _everboundedMAPsKinesins + 1.0;
							if (_everboundedMAPsKinesins == 1.0)
							{
								std::cout << "first protein bound at " << std::to_string(_state.currentTime) << " time" << std::endl;
							}
							//	std::cout << _everboundedMAPsKinesins << std::endl;
							iter = _unboundKinesins.erase(iter);

							//
						}
					}
					else
					{
						++iter;
					}
				}
			
			// update spring extensions based on previous microtubule step
					//Update spring liength for MAPs (not to do it in stepping or unbinding)
					for (auto iter = _boundMaps.begin(); iter != _boundMaps.end(); )
					{
						iter->_springLength = iter->_springLength + _state.MTpositionStep;
						++iter;
					}

					//Update spring liength for kinesins (not to do it in stepping or unbinding)
					for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); )
					{
						iter->_springLength = iter->_springLength + _state.MTpositionStep;
						++iter;
					}

					// Test for MAP unbinding
					if (_initC.MAPUnbinding == 1.0)
					{


						for (auto iter = _boundMaps.begin(); iter != _boundMaps.end(); ) {



							_MAPsunbindProbability = _mP.MAPKoff*exp(fabs(springForce(_mP.MAPstiffness, (iter->_springLength)))*_mP.MAPfsmParforKoff / (kBoltz*_mP.T));

							if (takeFlatRandomNumber() > exp((-(_MAPsunbindProbability)*_sP.timeStep)))
							{
								saveStepingMap(_state.currentTime, iter->_mountCoordinate, iter->_MTsite, -1);
								std::cout << "this MAP was unbound " << iter->_mountCoordinate << " at time " << _state.currentTime << " s, under force " << springForce(_mP.MAPstiffness, iter->_springLength) << "pN" << std::endl;
								_unboundMaps.emplace_back(iter->_mountCoordinate, iter->_id);
								iter = _boundMaps.erase(iter);



							}
							else
							{
								++iter;
							}
						}


					}


					// Test for kinesin unbinding
					if (_initC.kinesinUnbinding==1.0)
					{
						for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); ) {
							
						
							if (_initC.useKinesinOneparams==1.0)
							{
								//use kinesin-1 params
								if ((iter->_springLength) <= 0)
								{
								//assisting force
									_KinesinunbindProbability =( _mP.kinesinOneDet + _mP.kinesinOneLinfsp*(-springForce(_mP.KINESINstiffness, iter->_springLength)));
									
								}
								else
								{
									//hindering force
									_KinesinunbindProbability =  _mP.kinesinOneDet*exp(springForce(_mP.KINESINstiffness, iter->_springLength) *_mP.kinesinOneExpfsp / _mP.kT);
								}
								
							}
							else
							{
								//use kinesin CENP-E params
								
								_KinesinunbindProbability = 1 / (_mP.kinesinForceUnbindingA*exp(-fabs(springForce(_mP.KINESINstiffness, iter->_springLength) / _mP.kinesinForceUnbindingFd)));
							}

							if (takeFlatRandomNumber() > exp((-(_KinesinunbindProbability)*_sP.timeStep)))
							{
								saveStepingKinesin(_state.currentTime, iter->_mountCoordinate,iter->_MTsite,-1);
								std::cout << "this kinesin was unbound " << iter->_mountCoordinate << " at time " << _state.currentTime << " s, under force " << springForce(_mP.KINESINstiffness, iter->_springLength) << "pN" << std::endl;
								_unboundKinesins.emplace_back(iter->_mountCoordinate, iter->_id);
								iter = _boundKinesins.erase(iter);
								
								

							}
							else
							{
								++iter;
							}
						}
					}
					
////////////////////// Kinesin unbinding from dynamic MT
					if (_mP.dynamicMT == 1.0)
					{
						for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); ) {
							if ((iter->_MTsite) > _NumberofMTsites)
							{
								saveStepingKinesin(_state.currentTime, iter->_mountCoordinate, iter->_MTsite, -1);
								std::cout << "this kinesin was unbound " << iter->_mountCoordinate << " at time " << _state.currentTime << " s, under force " << springForce(_mP.KINESINstiffness, iter->_springLength) << "pN. From srinked MT" << std::endl;
								_unboundKinesins.emplace_back(iter->_mountCoordinate, iter->_id);
								iter = _boundKinesins.erase(iter);
							}
							else
							{
								++iter;
							}
						}
					}



				/// Test for stepping for MT bound MAPs AND UPDATE SPRING EXTENSION LENGTHS
						//// save position
						if ((taskIteration %_sP.saveFrequency == 0) && (_state.currentTime >= _sP.TimeMinLog)) {
						//	_positionstosave = std::to_string(_state.MTposition) +"	"+std::to_string(_state.currentTime)+"	";
							_positionstosave = "";
						}
						////
				for (auto iter = _boundMaps.begin(); iter != _boundMaps.end(); ) {
					
					
					
					

					///// saving position log  
					if ((taskIteration %_sP.saveFrequency == 0)&&(_state.currentTime>=_sP.TimeMinLog)) {
						_positionstosave = _positionstosave + std::to_string((iter->_id)) + "	"+std::to_string((iter->_MTsite)) + "	";
					}

					/////
					

					//std::cout << "iter->_springLength*_mP.MAPstiffness " << iter->_springLength*_mP.MAPstiffness << std::endl;
					///
					_kprob = _mP.MAPsDiffusion / (_mP.deltaPeriod*_mP.deltaPeriod);
					
					//_kpow = (-springForce(_mP.MAPstiffness, iter->_springLength))*_mP.MAPfsmPar / ( kBoltz*_mP.T);

					if (_mP.MAPAssymDiffusion ==1.0)
					{
						if (iter->_springLength >= 0.0)
						{
							//MAP is pulled towards minus end							
							_kMinus = _kprob*exp((fabs(springForce(_mP.MAPstiffness, iter->_springLength)))*_mP.MAPDiffSmLeft / (kBoltz*_mP.T));
							//_kPlus = _kprob*exp((-fabs(springForce(_mP.MAPstiffness, iter->_springLength)))*(_mP.deltaPeriod - _mP.MAPDiffSmLeft) / (kBoltz*_mP.T));
							_kPlus = _kprob*exp((-fabs(springForce(_mP.MAPstiffness, iter->_springLength)))*_mP.MAPDiffSmLeft / (kBoltz*_mP.T));

						}
						else
						{
							//MAP is pulled towards plus end

							_kPlus = _kprob*exp((fabs(springForce(_mP.MAPstiffness, iter->_springLength)))*_mP.MAPDiffSmRight / (kBoltz*_mP.T));
							//_kMinus = _kprob*exp((-fabs(springForce(_mP.MAPstiffness, iter->_springLength)))*(_mP.deltaPeriod - _mP.MAPDiffSmRight) / (kBoltz*_mP.T));
							_kMinus = _kprob*exp((-fabs(springForce(_mP.MAPstiffness, iter->_springLength)))*_mP.MAPDiffSmRight / (kBoltz*_mP.T));

						}
					}


					
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
								saveStepingMap(_state.currentTime, iter->_mountCoordinate, iter->_MTsite,2);
								
								// Count new summ forces
								 _state.SummMAPForces =  _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
								// std::cout << " ((_kMinus + _kPlus)*_sP.timeStep)= " << ((_kMinus + _kPlus)*_sP.timeStep) << std::endl;
								++iter;
							}
							else
							{
								if(takeFlatRandomNumber()<=_mP.MAPunbindBorderFreq)
								{
									_unboundMaps.emplace_back(iter->_mountCoordinate, iter->_id);
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
								saveStepingMap(_state.currentTime, iter->_mountCoordinate, iter->_MTsite,2);

								// Count new summ forces
								 _state.SummMAPForces =  _state.SummMAPForces + springForce(_mP.MAPstiffness, iter->_springLength);
								++iter;
							}
							else
							{
								if (takeFlatRandomNumber() <= _mP.MAPunbindBorderFreq)
								{
									_unboundMaps.emplace_back(iter->_mountCoordinate, iter->_id);
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
				///// saving position log  saveMapPosition
				if ((taskIteration %_sP.saveFrequency == 0) && (_state.currentTime >= _sP.TimeMinLog)) {
					_positionstosave = _positionstosave + "\n";
					saveMapPosition(_positionstosave);
					_positionstosave = "";
					//_positionstosave = std::to_string(_state.MTposition) + "	"+std::to_string(_state.currentTime) + "	";
				}

				/////



				////
				/// Test for stepping for MT bound MT kinesins AND UPDATE SPRING EXTENSION LENGTHS
				for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); ) {
					

					///// saving position log  
					if ((taskIteration %_sP.saveFrequency == 0) && (_state.currentTime >= _sP.TimeMinLog)) {
						_positionstosave = _positionstosave + std::to_string((iter->_id)) + "	" + std::to_string((iter->_MTsite)) + "	";
					}

					/////
					
					
					
					//std::cout << "iter->_springLength*_mP.KINESINstiffness " << iter->_springLength*_mP.KINESINstiffness << std::endl;
					//
					double KinesinVelocity = (_mP.vUnloaded / (_mP.kinesinPparam + ((1 - _mP.kinesinPparam)*exp(springForce(_mP.KINESINstiffness, iter->_springLength)*_mP.kinesinDparam / (kBoltz*_mP.T)))));
					double pstep =exp(-KinesinVelocity*(_sP.timeStep / (2 * _mP.deltaPeriod)));
					//double pstep = 0.0;


					if (takeFlatRandomNumber() > pstep)
					{
						//make step (always 2 site periods to the plus end of MT therefore to higher site number)
						if (iter->_MTsite + 2 <= _NumberofMTsites )
						{
							iter->_MTsite = iter->_MTsite + 2;
							iter->_springLength = iter->_springLength + 2* _mP.deltaPeriod;
							//
							
							
							saveStepingKinesin(_state.currentTime, iter->_mountCoordinate, iter->_MTsite,2);
							// Count new summ forces
							 _state.SummKINESINForces =  _state.SummKINESINForces + (springForce(_mP.KINESINstiffness, iter->_springLength));
							++iter;
						}
						else
						{
							if (takeFlatRandomNumber() <= _mP.KinesinunbindBorderFreq)
							{
								_unboundKinesins.emplace_back(iter->_mountCoordinate, iter->_id);
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
				/// save Kinesin positions log
				if ((taskIteration %_sP.saveFrequency == 0) && (_state.currentTime >= _sP.TimeMinLog)) {
					_positionstosave = _positionstosave + "\n";
					saveKinesinPosition(_positionstosave);
					_positionstosave = "";
				}
				////
		//////
				// stop simulations if forces too high
				if ((fabs(_state.SummMAPForces)>1000000.0))
				{
					// basic coordinates log
					writeStateTolog();
					return "Calculations aborted at simulation time " + std::to_string(_state.currentTime) + " due to total MAP forces exceed 1000000.0" ;
					
				}
				if ((fabs(_state.SummKINESINForces)>1000000.0))
				{
					// basic coordinates log
					writeStateTolog();
					return "Calculations aborted at simulation time " + std::to_string(_state.currentTime) + " due to total Kinesin forces exceed 1000000.0";
					
				}
				////
				// stop simulation if no MAPs and Kinesins bound, but only if parameter stopSimIfNoMAPsKinesins set to 1.0
				if ((_initC.stopSimIfNoMAPsKinesins == 1.0) && (_everboundedMAPsKinesins > 0.0))
				{
					if ((_boundKinesins.size()==0)&&(_boundMaps.size()==0)&&(_state.currentTime>=_sP.stopSimIfNoMAPsKinesinsTmin))
					{
						
				
						terminate();
						return "Calculations aborted at simulation time " + std::to_string(_state.currentTime) + " due to no MAPs and Kinesins left on the MT and stopSimIfNoMAPsKinesins is on";

					}
				}

				//////////////////

			_SummForces = -_state.SummMAPForces -  _state.SummKINESINForces + _mP.constantForceonMT;
			
			////
			 
	
				_state.MTpositionStep = (_sP.timeStep / _mP.gammaMT)*_SummForces + _mP.thermalNoiseOn*sqrt(2.0 * _mP.kT*_sP.timeStep / _mP.gammaMT) *	takeNormalRandomNumber();
		
			
			 ////
			 
			 _state.MTposition = _state.MTposition +  _state.MTpositionStep;
			///
			_state.BoundedKinesins=_boundKinesins.size();
			_state.BoundedMAPs = _boundMaps.size();
			
			// benchmark
			//auto end = std::chrono::high_resolution_clock::now();


			//_timebenchmark = _timebenchmark + (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
			// benchmark


			if ((taskIteration %_sP.saveFrequency == 0) && (_state.currentTime >= _sP.TimeMinLog)) {
			//	std::cout << takeFlatRandomNumber() << std::endl;
			//	std::cout << takeNormalRandomNumber() << std::endl;
				//benchmark log
				//_performanceLog.save(std::to_string(_timebenchmark / _sP.saveFrequency) + "\n");
				//_timebenchmark = 0.0;
				//benchmark log
				
				// basic coordinates log
				writeStateTolog();
				// extnded stepping log if selected
				if(_initC.fulllogMAPsKinesins==1.0)
				{
					for (auto iter = _boundMaps.begin(); iter != _boundMaps.end(); )
					{
						saveStepingMap(_state.currentTime, iter->_mountCoordinate, iter->_MTsite,0);
						++iter;
					}
					for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); )
					{
						saveStepingKinesin(_state.currentTime, iter->_mountCoordinate, iter->_MTsite,0);
						++iter;
					}
				}
				///
				
				
				
			}

			
			
			
			
		}
		terminate();
		///////
		return "No errors encountered";
	}
	void saveSteppingDatFile()
	{
		_saveStepingMapBuffer = std::accumulate(begin(_MapBuffer), end(_MapBuffer), std::string{}, [](std::string result, const SingleProteinSteplog& element) {
			return
				std::move(result)
				+ std::to_string(element.time) + "	"
				+ std::to_string(element.mountCoordinate) + "	"
				+ std::to_string(element.MTsiteNum) + " "
				+ std::to_string(element.MTposition) + " "
				+ std::to_string(element.proteinstate) + "\n";
		});
		///
		_saveStepingKinesinBuffer = std::accumulate(begin(_KinesinBuffer), end(_KinesinBuffer), std::string{}, [](std::string result, const SingleProteinSteplog& element) {
			return
				std::move(result)
				+ std::to_string(element.time) + "	"
				+ std::to_string(element.mountCoordinate) + "	"
				+ std::to_string(element.MTsiteNum) + " "
				+ std::to_string(element.MTposition) + " "
				+ std::to_string(element.proteinstate) + "\n";
		});
		///
		writeFinalSteppingLogBuffer();
	}
	void saveSteppingBinaryFile()
	{
		// Binary MAPs stepping log
		std::ofstream MAPstepsfile{ _loggerparams.filepath + _loggerparams.name + "MAPsteps.binary" , std::ios::binary };
		const std::vector<SingleProteinSteplog> continiousMapStepsData{ begin(_MapBuffer), end(_MapBuffer) };
		//std::cout << "data vector size = " << continiousStepsData.size() << std::endl;
		MAPstepsfile.write(reinterpret_cast<const char*>(continiousMapStepsData.data()), continiousMapStepsData.size() * sizeof(SingleProteinSteplog));
		MAPstepsfile.close();
		// Binary Kinesins stepping log
		std::ofstream Kinesinstepsfile{ _loggerparams.filepath + _loggerparams.name + "Kinesinsteps.binary" , std::ios::binary };
		const std::vector<SingleProteinSteplog> continiousKinesinStepsData{ begin(_KinesinBuffer), end(_KinesinBuffer) };
		//std::cout << "data vector size = " << continiousStepsData.size() << std::endl;
		Kinesinstepsfile.write(reinterpret_cast<const char*>(continiousKinesinStepsData.data()), continiousKinesinStepsData.size() * sizeof(SingleProteinSteplog));
		Kinesinstepsfile.close();
	}
	void terminate()
	{
		//saveSteppingDatFile();
		saveSteppingBinaryFile();

		// basic coordinates log
		writeStateTolog();
		//log final state to json
		dumpFullStatetoJson("_final_state_dump");
			
	}
	void dumpFullStatetoJson(std::string jsonstatefilename)
	{
		/////// Json Dump of final state
		_stateDump["MTposition"] = (std::to_string(_state.MTposition));
		_stateDump["MTpositionStep"] = (std::to_string(_state.MTpositionStep));
		//_stateDump["boundMaps"].push_back(nullptr);
		//_stateDump["unboundMaps"].push_back(nullptr);
		//_stateDump["boundKinesins"].push_back(nullptr);
		//_stateDump["unboundKinesins"].push_back(nullptr);
		for (auto iter = _boundMaps.begin(); iter != _boundMaps.end(); )
		{
			json f;
			f["mountCoordinate"] = (std::to_string(iter->_mountCoordinate));
			f["MTsite"] = (std::to_string(iter->_MTsite));
			f["springLength"] = (std::to_string(iter->_springLength));
			f["id"] = (std::to_string(iter->_id));
			_stateDump["boundMaps"].push_back(f);
			///////
			++iter;
		}
		for (auto iter = _boundKinesins.begin(); iter != _boundKinesins.end(); )
		{
			json f;
			f["mountCoordinate"] = (std::to_string(iter->_mountCoordinate));
			f["MTsite"] = (std::to_string(iter->_MTsite));
			f["springLength"] = (std::to_string(iter->_springLength));
			f["id"] = (std::to_string(iter->_id));
			_stateDump["boundKinesins"].push_back(f);
			/////
			++iter;
		}
		for (auto iter = _unboundMaps.begin(); iter != _unboundMaps.end(); )
		{
			json f;
			f["mountCoordinate"] = (std::to_string(iter->_mountCoordinate));
			f["id"] = (std::to_string(iter->_id));
			_stateDump["unboundMaps"].push_back(f);
			///////
			++iter;
		}
		for (auto iter = _unboundKinesins.begin(); iter != _unboundKinesins.end(); )
		{
			json f;
			f["mountCoordinate"] = (std::to_string(iter->_mountCoordinate));
			f["id"] = (std::to_string(iter->_id));
			_stateDump["unboundKinesins"].push_back(f);
			/////
			++iter;
		}
		std::ofstream o(_loggerparams.filepath + _loggerparams.name + jsonstatefilename + ".json");
		o << std::setw(4) << _stateDump << std::endl;
		_stateDump.clear();
		o.close();
		/// json dump end
	}

	void writeStateTolog() const {
		for (const auto& logger : _loggers) {
			logger->save(&_state);
		}
	}
	void writeFinalSteppingLogBuffer()
	{
		_MAPStepsLog.save(_saveStepingMapBuffer);
		_kinesinStepsLog.save(_saveStepingKinesinBuffer);
	}

private:
	double _intialMTposition;
	unsigned _seedMount;
	json _stateDump;
	// benchmark
	//DatFileLogger _performanceLog;
	//double _timebenchmark=0.0;
	// benchmark
	const SimulationParameters _sP;
	const ModelParameters _mP;
	const InitialConditions _initC;
	const LoggerParameters _loggerparams;
	const InitialSetup _iSetup;
	SystemState _state;
	std::vector<std::unique_ptr<BinaryFileLogger>> _loggers;
	
	int _NumberofMTsites;
	
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
	
	double _MTsystemCoordinate;
	double  _fractpart, _intpart;
	double _KinesinunbindProbability;
	double _MAPsunbindProbability;
	double _sitenum;
	double _everboundedMAPsKinesins;
//	double _currentSpringLength;
	
	DatFileLogger _kinesinStepsLog; //log of stepping
	DatFileLogger _MAPStepsLog;//log of stepping

	

	DatFileLogger _kinesinPosLog;
	DatFileLogger _MAPPosLog;

	

	std::string _positionstosave="";

	int _saveStepingMapCounter = 0;
	int _saveStepingKinesinCounter = 0;
	std::string _saveStepingMapBuffer="";
	std::string _saveStepingKinesinBuffer="";
	int _saveStepingMapBufferSize=4000;
	int _saveStepingKinesinBufferSize=1000;

	unsigned _currentstep = 0;
	unsigned _totalsteps = 0;
public:
	size_t _optimalGaussianBuffer = 100000;//50000;
	size_t _optimalFlatBuffer = 500000;//200000;

	std::deque <SingleProteinSteplog> _MapBuffer;
	std::deque <SingleProteinSteplog> _KinesinBuffer;
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
	const auto initialSetup = load_setup(simulationParameters.initialSetupfile, simulationParameters.useInitialSetup);
	//const auto initialSetup= load_setup(inputparamfile);

	std::vector<std::unique_ptr<Task>> tasks;
	for (const auto& configuration : configurations) {
		auto task = std::make_unique<Task>(simulationParameters, configuration, initialSetup);
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

	MklGaussianParallelGenerator gaussGenerator(0.0, 1.0, tasks[0]->_optimalGaussianBuffer, 1, simulationParameters.extGaussianseed);
	MklFlatParallelGenerator flatGenerator;
	flatGenerator.initialize(0.0, 1.0, tasks[0]->_optimalFlatBuffer, 4, simulationParameters.extUniformseed);
	//std::cout << macrostepsNum << std::endl;
	//for (int savedstep = 0; savedstep < (100'000); savedstep++) {

		//for (int macrostep = 0; macrostep < macrostepsNum; macrostep++) {
			
			
			gaussGenerator.generateNumbers();
			const auto gaussBuffData = gaussGenerator.getNumbersBuffer();
			flatGenerator.generateNumbers();
			const auto flatBuffData = flatGenerator.getNumbersBuffer();

			std::string taskvalidation= tasks[0]->advanceState(taskStepsNum, gaussBuffData, flatBuffData,&gaussGenerator,&flatGenerator);
		
			std::cout << taskvalidation << std::endl;
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

			//std::cout << " end" << std::endl;
			std::cout << " end" << std::endl;
}

//std::chrono
//_rdtscp  in #immintrin.h

//		(*/),(+-)