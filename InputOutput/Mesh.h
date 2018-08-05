/*=================================================
                     mesh class header
=================================================*/

#ifndef MESH_H_
#define MESH_H_

class Utility;

class Mesh
{
public:
	Mesh(char pub_name[]);
	void Geometry(Utility*);
	void Connectivity();
	void GhostCellStructure();
	int GetCelltot() const;
	int GetNodetot() const;
	int GetEdgeNumber() const;
	double* GetMinDis() const;
	int** GetCellmap() const;
	int** GetConnectivity() const;
	double* GetArea() const;
	double** GetDl() const;
	double** GetNx() const;
	double** GetNy() const;
	double* GetX() const;
	double* GetY() const;
	int** GetMapWall() const;
	int** GetMapInlet() const;
	int** GetMapOutlet() const;
	int GetCellTotalWall() const;
	int GetCellTotalInlet() const;
	int GetCellTotalOutlet() const;
	int GetNeighborFaceNum(int,int,int) const;


private:
	char fileName[100];
	int celltot=0;
	int nodetot=0;
	int edgeNumber=0;
	int **cellmap=nullptr;
	int **neigh=nullptr;
	double* x=nullptr;
	double* y=nullptr;
	int bouTypes=0;
	int* bouSizes=nullptr;
	int*** bouPoints=nullptr;
	double** dl=nullptr;
	double** normx=nullptr;
	double** normy=nullptr;
	double *area=nullptr;
	double *xc=nullptr;
	double *yc=nullptr;
	double* minDis=nullptr;
	int** mapWall=nullptr;
	int** mapInlet=nullptr;
	int** mapOutlet=nullptr;
	int cellTotWall=0;
	int cellTotInlet=0;
	int cellTotOutlet=0;
};

#endif
