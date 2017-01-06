#include "mkl_flat_parallel_generator.h"
#include <stdlib.h>
#include <stdexcept>
#include <omp.h>



MklFlatParallelGenerator::MklFlatParallelGenerator() {}
	

void MklFlatParallelGenerator::initialize(double leftBound, double rightBound, std::size_t bufferSize, unsigned threadNum) 
{
	_leftBound = leftBound;
	_rightBound = rightBound;
	_bufferSize = bufferSize;
	_threadNum = threadNum;
	///////////////// If reproducibility from launch to launch is required seed is const, eslse seed must be random
	MKL_UINT seed = __rdtsc();
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
