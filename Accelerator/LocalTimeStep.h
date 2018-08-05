/*
 * LocalTimeStep.h
 *
 *  Created on: Apr 18, 2018
 *      Author: arash
 */

#ifndef LOCALTIMESTEP_H_
#define LOCALTIMESTEP_H_

class Mesh;
class SteadyReader;

class LocalTimeStep {
public:
	LocalTimeStep(Mesh*,SteadyReader*,double*);
	void CalculateLTS();
};

#endif /* LOCALTIMESTEP_H_ */
