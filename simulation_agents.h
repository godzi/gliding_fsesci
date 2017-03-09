#pragma once

class boundKINESIN
{
public:
	boundKINESIN(double coordinate, int MTsite, double springLength, int id);
	double _mountCoordinate;
	int _MTsite;
	double _springLength;
	int _id;
};
class unboundKINESIN
{
public:
	unboundKINESIN(double coordinate, int id);
	double _mountCoordinate;
	int _id;
};
class boundMAP
{
public:
	boundMAP(double coordinate, int MTsite, double springLength, int id);
	double _mountCoordinate;
	int _MTsite;
	double _springLength;
	int _id;
};
class unboundMAP
{
public:
	unboundMAP(double coordinate, int id);
	double _mountCoordinate;
	int _id;
};


