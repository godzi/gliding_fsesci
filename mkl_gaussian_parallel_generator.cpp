#include "mkl_gaussian_parallel_generator.h"
#include <stdlib.h>
#include <stdexcept>
#include <omp.h>
#include <ctime>
#include <thread>
#include <chrono>
#include <random>
#include <iostream>

 unsigned build_random_gauss_seed() {
	unsigned seedMount = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generatorMount(seedMount);
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,99999999);
	unsigned result= static_cast<unsigned>(distribution(generatorMount));;
//	std::cout << result << std::endl;
	return result;

}


MklGaussianParallelGenerator::MklGaussianParallelGenerator(double mean, double stDeviation, std::size_t bufferSize, unsigned threadNum, unsigned extseed)
	:_mean{ mean }, _stDeviation{ stDeviation }, _bufferSize{ bufferSize }, _threadNum{ threadNum }, _extseed{extseed}
{
	///////////////// If reproducibility from launch to launch is required seed is const, eslse seed must be random
	//MKL_UINT seed = __rdtsc();
	//MKL_UINT unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
	//MKL_UINT seed = static_cast<unsigned>(std::time(0))*static_cast<unsigned>(std::hash<std::thread::id>()(std::this_thread::get_id()));
	MKL_UINT seed;
	if (extseed==0)
	{
		seed = build_random_gauss_seed();
	}
	else
	{
		seed = extseed;
	}
	std::cout << "gauss seed " <<  seed << std::endl;
	/////////////////
	for (unsigned i = 0; i < threadNum; i++) {
		_streamWrappers.emplace_back(VSL_BRNG_MT2203 + i, seed);
	}
	_nPerThread = _bufferSize / _threadNum;
	if (_bufferSize != _nPerThread * _threadNum) {
		throw std::logic_error{ "buffsize must be multiple of number of threads" };
	}
	_buffer.resize(_bufferSize);
}



void MklGaussianParallelGenerator::generateNumbers()
{
	//// Strange that it works without shared! check this out.
#pragma omp parallel num_threads(_threadNum) default(none)
	{
		const std::size_t threadId=omp_get_thread_num();
		const std::size_t begin = _nPerThread * threadId;
		const auto generationResult = vdRngGaussian(
			VSL_RNG_METHOD_GAUSSIAN_ICDF, _streamWrappers.at(threadId).get(),
			static_cast<MKL_INT>(_nPerThread), _buffer.data() + begin,
			_mean, _stDeviation
		);
		if (generationResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "can't generate numbers" };
		}
	}

	
}

const double * MklGaussianParallelGenerator::getNumbersBuffer() const
{
	return _buffer.data();
}

std::size_t MklGaussianParallelGenerator::getNumbersBufferSize() const
{
	return _bufferSize;
}
