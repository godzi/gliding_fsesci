#pragma once
#include <fstream>
#include "configuration.h"
class BinaryFileLogger
{
public:
	BinaryFileLogger(LoggerParameters loggerParams, double(SystemState::* loggedField), std::string coordinateName);
	~BinaryFileLogger();
	void save(const SystemState* systemState);

private:
	void flush();

	static constexpr std::size_t _buffsize = 4096 / sizeof(double); //previous buffer 4096
	std::ofstream _file;
	double(SystemState::* _loggedField);
	std::vector <double> _buffer;
};