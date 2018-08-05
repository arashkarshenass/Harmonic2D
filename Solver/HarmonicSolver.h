/*========================================
 	   Harmonic Solver Class Header
=========================================*/

#ifndef SOLVER_HARMONICSOLVER_H_
#define SOLVER_HARMONICSOLVER_H_

class Mesh;
class SteadyReader;
class SteadyJacobian;
class Utility;
class OneCell;
class Residue;


class HarmonicSolver {
public:
	HarmonicSolver(Mesh*,SteadyReader*,SteadyJacobian*,Utility*);
	void Itterator();
	double** GetConservativeReal()const;
	double** GetConservativeImag()const;

private:
	int celltot=0;								//total cell number
	double* vol=nullptr;					//volume of cells
	double* dt=nullptr;						//time vector from local time step
	double** wr=nullptr;					//conservative variables (real)
	double** wi=nullptr;					//conservative variables (imaginary)
	double** wrOld=nullptr;				//conservative variables from previous iteration (real)
	double** wiOld=nullptr;				//conservative variables from previous iteration (imaginary)
	double** dwr=nullptr;					//cells total flux (real)
	double** dwi=nullptr;					//cells total flux (imaginary)
	double* cellFluxReal=nullptr;				//total real fluxes of a cell
	double* cellFluxImag=nullptr;			//total imag fluxes of a cell
	double* resR=nullptr;					//real residues
	double* resIm=nullptr;					//imaginary resudies
	double radiF=0;								//radial frequency
	OneCell* ocp=nullptr;					//pointer to one cell object
	Residue* rp=nullptr;						//pointer to residue object
	Utility* up=nullptr;						//ponter to utility object
	int itr=0;



};

#endif
