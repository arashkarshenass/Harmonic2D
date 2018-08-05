/*===================================
 	 	 Steady Solution Class Header
====================================*/

#ifndef INPUTOUTPUT_STEADYREADER_H_
#define INPUTOUTPUT_STEADYREADER_H_

class Mesh;

class SteadyReader {
public:
	SteadyReader(char pubname[],Mesh*);
	double** GetPrimitiveDomain() const;
	double** GetUnDom() const;
	double** GetUtDom() const;

private:
	double** qsDom=nullptr;
	double** UnMatDom=nullptr;
	double** UtMatDom=nullptr;
};

#endif
