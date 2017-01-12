#pragma once
#include <iostream>
#include <fstream>
#include "configuration.h"
#include "binary_file_logger.h"


BinaryFileLogger::BinaryFileLogger(LoggerParameters loggerParams, double(SystemState::* loggedField), std::string coordinateName) :
		_file{ loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".binary", std::ios::binary },

		//_file{ loggerParams.filepath + loggerParams.name + "_results_" + coordinateName + ".dat", std::ofstream::out },
		_loggedField{ loggedField }
	{
		if (!_file) {
			throw std::runtime_error{ "the file was not created" };
		}
		_buffer.reserve(_buffsize);
	}
BinaryFileLogger::~BinaryFileLogger() {
		flush();
	}
void BinaryFileLogger::save(const SystemState* systemState) {
		_buffer.push_back(systemState->*_loggedField);
		if (_buffer.size() == _buffsize) {
			flush();
		}
	}

void BinaryFileLogger::flush() {
		if (_buffer.empty()) {
			return;
		}
		_file.write(reinterpret_cast<const char*>(_buffer.data()), _buffer.size() * sizeof(double));
		//_file.write(reinterpret_cast<const char*>(_buffer.data()), _buffer.size() * sizeof(double));


		if (!_file.good()) {
			throw std::runtime_error{ "not all data was written to file" };
		};
		_buffer.clear();
	}

	
