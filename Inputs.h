/*=================================================
                     inputs
=================================================*/

#ifndef INPUTS_H_
#define INPUTS_H_

//constants
#define cp 1006				//specific heat of air in constant pressure (j/kg.K)
#define cv 717.1			//specific heat of air in constant volume (j/kg.K)
#define gamma 1.4			//ratio of specific heat (-)
#define R 287.058			//gas constant of air (j/kg.K)
#define pi 3.1415926535		//value of pi

//inlet boundary conditions
#define tStagIn 55.555      	//stagnation inlet temperature (K)
#define pStagIn 6894.76   	//stagnation inlet pressure (Pa)
#define thetaIn 0			//flow angle at inlet (degree). It is the angle between flow direction and x axis. + for _\ and - for """/
							//below x-axis (+ for _\  and - for -/

//outlet boundary condition
#define stdPrsOut 6136.334		//static oulet pressure
#define dP 0.01					//harmonic outlet pressure amplitude ratio to steady outlet pressure
#define redF 1					//reduced frequency of harmonic BC
#define harmonicType 2		//type of harmonic BC (1->sinus	2->cosinus)

//controllers
#define cfl 0.01				//cfl number
#define resRef -7			//residue magnitude for stopping (10^resRef)
#define nfres 100			//frequency of calculation and printing residual
#define scaleRes 1			//use scaling of residues?		1->Yes		any oter value->No
#define prntmsh 0			//print mesh results? 0->No		1->YES
#define prntstd 0			//print steady solution? 0->No	1->Yes
#define prntSln 1			//print perturbed solution?		1->Yes		any oter value->No
#define periodDevision 5	//divide one period to how many time instances for solution in time

//initializers
#define wr0 0
#define wr1 0
#define wr2 0
#define wr3 0
#define wi0 0
#define wi1 0
#define wi2 0
#define wi3 0

#endif
