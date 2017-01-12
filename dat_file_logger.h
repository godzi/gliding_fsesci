#pragma once
#include <fstream>
#include "configuration.h"
class DatFileLogger
{
public:
	DatFileLogger(LoggerParameters loggerParams, std::string coordinateName);
	
	void save(const std::string& dataToLog);

private:
	std::ofstream _file;
	
};