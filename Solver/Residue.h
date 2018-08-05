/*=================================================
                Residue Class header
=================================================*/

#ifndef RESIDUE_H_
#define RESIDUE_H_

class Solver;
class Utility;

class Residue
{
public:
	Residue(double**,double**,double*,double*, Utility*,int);
	void FirstIteration();
	void CalculateResidue(int);
	double findMax();

private:
	int celltot=0;
	double** dwR=nullptr;
	double** dwIm=nullptr;
	double* resR=nullptr;
	double* resIm=nullptr;
	double* dwRFirst=nullptr;
	double* dwImFirst=nullptr;
	Utility* up=nullptr;
};

#endif
