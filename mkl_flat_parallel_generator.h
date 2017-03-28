#pragma once
#include "i_generator.h"
#include "vsl_stream_wrapper.h"
#include <vector>

unsigned build_random_flat_seed();
class MklFlatParallelGenerator :
	public IGenerator
{
public:
	MklFlatParallelGenerator();
	void initialize(double leftBound, double rightBound, std::size_t bufferSize, unsigned threadNum, unsigned extseed);
	virtual void generateNumbers() override;
	virtual const double* getNumbersBuffer() const override;
	virtual std::size_t getNumbersBufferSize() const override;
private:
	double _leftBound;
	double _rightBound;
	std::size_t _bufferSize;
	unsigned _threadNum;
	std::vector <VSLStreamWrapper> _streamWrappers;
	std::vector<double> _buffer;
	std::size_t _nPerThread;
	unsigned _extseed;
};

