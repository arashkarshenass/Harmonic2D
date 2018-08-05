/*=================================================
                Wall Boundary Class header
=================================================*/
#ifndef WALL_H_
#define WALL_H_
class Mesh;
class SteadyJacob;
class SteadyReader;
class Wall
{
public:
	void Initializer(double*,double*,double*,double*);
	void CalculateGhost();

private:
	double*  wrFaceDomain=nullptr;
	double*  wiFaceDomain=nullptr;
	double*  wrFaceGhost=nullptr;
	double*  wiFaceGhost=nullptr;
};

#endif
