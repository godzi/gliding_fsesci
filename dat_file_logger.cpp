#pragma once
#include <iostream>
#include <fstream>
#include "configuration.h"
#include "dat_file_logger.h"


DatFileLogger::DatFileLogger(LoggerParameters loggerParams, std::string coordinateName) :
		_file{ loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".dat" }
	{
		if (!_file) {
			throw std::runtime_error{ "the file was not created" };
		}
	
	}
DatFileLogger::~DatFileLogger()
{
	_file.close();
}
void DatFileLogger::save(const std::string& dataToLog) {
	_file << dataToLog;
	}


	
