#include "mkl_flat_parallel_generator.h"
#include <stdlib.h>
#include <stdexcept>
#include <omp.h>
#include <ctime>
#include <thread>
#include <iostream>
#include <random>

unsigned build_random_flat_seed() {
	unsigned seedMount = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generatorMount(seedMount);
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, 99999999);
	unsigned result = static_cast<unsigned>(distribution(generatorMount));;
	//	std::cout << result << std::endl;
	return result;

}

MklFlatParallelGenerator::MklFlatParallelGenerator() {}
	

void MklFlatParallelGenerator::initialize(double leftBound, double rightBound, std::size_t bufferSize, unsigned threadNum, unsigned extseed)
{
	_leftBound = leftBound;
	_rightBound = rightBound;
	_bufferSize = bufferSize;
	_threadNum = threadNum;
	_extseed = extseed;
	///////////////// If reproducibility from launch to launch is required seed is const, eslse seed must be random
	//MKL_UINT seed = __rdtsc();
	MKL_UINT seed;
	if (extseed == 0) {
		seed = build_random_flat_seed();
	}
	else
	{
		seed = extseed;
	}
	//MKL_UINT seed = static_cast<unsigned>(std::time(0))*static_cast<unsigned>(std::hash<std::thread::id>()(std::this_thread::get_id()));
	std::cout << "flat seed "<< seed << std::endl;
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

void MklFlatParallelGenerator::generateNumbers()
{
	//// Strange that it works without shared! check this out.
#pragma omp parallel num_threads(_threadNum) default(none)
	{
		const std::size_t threadId = omp_get_thread_num();
		const std::size_t begin = _nPerThread * threadId;
		const auto generationResult = vdRngUniform(
			VSL_RNG_METHOD_UNIFORM_STD, _streamWrappers.at(threadId).get(),
			static_cast<MKL_INT>(_nPerThread), _buffer.data() + begin,
			_leftBound, _rightBound
		);
		if (generationResult != VSL_STATUS_OK) {
			throw std::runtime_error{ "can't generate numbers" };
		}
	}


}

const double * MklFlatParallelGenerator::getNumbersBuffer() const
{
	return _buffer.data();
}

std::size_t MklFlatParallelGenerator::getNumbersBufferSize() const
{
	return _bufferSize;
}
