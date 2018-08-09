
#ifndef INPUTOUTPUT_SOLUTIONCALCULATOR_H_
#define INPUTOUTPUT_SOLUTIONCALCULATOR_H_

class Mesh;
class SteadyReader;
class HarmonicSolver;


struct Frequency{
	double* w1;double* w2;double* w3;double* w4;double* rho;double* u;double* v;double* P;double* M;
};


struct TimeSolution{
	double time;
	double* w1;double* w2;double* w3;double* w4;double* rho;double* u;double* v;double* P;double* M;
};


class SolutionCalculator {
public:
	SolutionCalculator(Mesh*,SteadyReader*,HarmonicSolver*);
	double PhiCalculator(double,double);
	Frequency* GetAmplitude();
	Frequency* GetPhi();
	TimeSolution* GetHarmonicSinus();
	TimeSolution* GetHarmonicCosinus();
	TimeSolution* GetUnsteadySinus();
	TimeSolution* GetUnsteadyCosinus();
	int GetInstantNumber();
	double GetDt();

private:
	int instantTotal=0;
	double dt=0;
	Frequency amplitude;
	Frequency phi;
	TimeSolution* harmonicSinus;
	TimeSolution* harmonicCosinus;
	TimeSolution* unsteadySinus;
	TimeSolution* unsteadyCosinus;

};

#endif
