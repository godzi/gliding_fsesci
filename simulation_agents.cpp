#include "simulation_agents.h"


boundKINESIN::boundKINESIN(double coordinate, int MTsite, double springLength, int id) :_mountCoordinate{ coordinate }, _MTsite{ MTsite }, _springLength{ springLength }, _id{ id }
{	

}

unboundKINESIN::unboundKINESIN(double coordinate, int id) 
{
		_mountCoordinate = coordinate;
		_id = id;
}

boundMAP::boundMAP(double coordinate, int MTsite, double springLength, int id) :_mountCoordinate{ coordinate }, _MTsite{ MTsite }, _springLength{ springLength }, _id{ id }
{	

}

unboundMAP::unboundMAP(double coordinate, int id)
{
		_mountCoordinate = coordinate;
		_id = id;
}
	