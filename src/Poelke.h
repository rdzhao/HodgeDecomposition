#pragma once

#include "My_Triangulation.h"
#include "RKIntegrator.h"

using namespace std;

class Poelke
{
public:
	Poelke() {};

	My_Triangulation& Mesh();
	std::vector<Mesh_3_Point_3>& Vertices();
	std::vector<int>& Indices();
	std::vector<Mesh_3_Vector_3>& Normals();
	std::vector<Mesh_3_Vector_3>& Colors();
	std::vector<Mesh_3_Vector_3>& VfVertices();
	std::vector<Mesh_3_Vector_3>& VfFaces();
	std::vector<Mesh_3_Vector_3>& VfNormals();
	std::vector<Mesh_3_Vector_3>& VfColors();

	void init(double dsty, double sr, int stps, std::string sfx);
	void buildMeshFromSurface(std::string fn, double size);
	void prepare();
	void decompose();
	void visualize();
	void integrate();

public: 
	// analytical functions
	Mesh_3_Vector_3 normalGradient(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 tangentialCurl(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 normalHarmonic(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 tangentialHarmonic(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 centralHarmonic(Mesh_3_Vector_3 p);

public:
	void createInput();

	void solveNormalGradient();
	void solveTangentalCurl();
	void solveNormalHarmonic();
	void solveTangentialHarmonic();
	void solveCentralHarmonic();

	void solveTHBasis();
	void solveNGBasis();

	void computeArrows();
	void computeCrossSection();

	void evaluateError();
	void evaluateL2Norm();

public:
	void createCurlOperator_B();
	void createGradOperator_B();
	void createCurlOperator_NB();
	void createGradOperator_NB();

	void createCNbGFOperator();
	void createGF0GFOperator();

	void createMassMatrix();

public:
	void rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci);

public:
	EigenVector matlabLinearSolveQR(EigenSpMat A, EigenVector b);
	EigenVector matlabLinearSolveCholesky(EigenSpMat A, EigenVector b);
	void matlabSVDS(std::vector<EigenVector>& svec, std::vector<double>& sval, EigenSpMat A, int k);
	void matlabEIGS(std::vector<EigenVector>& evec, std::vector<double>& eval, EigenSpMat A, EigenSpMat B, int k);

	void concatenateVerticallySPMatrix(EigenSpMat& M, EigenSpMat A, EigenSpMat B);
	void concatenateHorizontallySPMatrix(EigenSpMat& M, EigenSpMat A, EigenSpMat B);
	void concatenateVector(EigenVector& V, EigenVector A, EigenVector B);

	void toEigenVector(EigenVector& v, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void fromEigenVector(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, EigenVector v);

	

public:
	double innerProduct(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf1, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf2);
	
	double computeFieldError(int fieldType);
	double computeFieldInnerProduct(int fieldType);
	double cellErrorQuadrature(Mesh_3_Cell_iterator ci, Mesh_3_Vector_3 v, int fieldType);
	double cellInnerProductQuadrature(Mesh_3_Cell_iterator ci, int fieldType);

	Mesh_3_Vector_3 averageNormalGradient(Mesh_3_Facet_iterator fi);
	Mesh_3_Vector_3 averageTangentialCurl(Mesh_3_Facet_iterator fi);
	Mesh_3_Vector_3 averageCentralHarmonic(Mesh_3_Facet_iterator fi);

public: //test functions
	void testOrthogonalityOfBasis();
	void checkCurlFieldBoundary();

private:
	std::string m_suffix;
	double m_density_ratio;
	double m_step_ratio;
	double m_steps;

private:
	My_Triangulation m_mesh_3;

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_input;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_original_input;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_normalGradient;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_tangentialCurl;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_normalHarmonic;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_tangentialHarmonic;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> m_centralHarmonic;

	EigenSpMat m_curlOp_B;
	EigenSpMat m_gradOp_B;
	EigenSpMat m_curlOp_NB;
	EigenSpMat m_gradOp_NB;
	EigenSpMat m_CNbGF;
	EigenSpMat m_GF0GF;
	EigenSpMat m_Mass;

	EigenSpMat m_thOp;
	EigenSpMat m_nhOp;

private:
	Mesh_3_Point_3 center;
	std::vector<Mesh_3_Point_3> vertices;
	std::vector<Mesh_3_Vector_3> colors;
	std::map<Mesh_3_Vertex_iterator, int> vertexMap;
	std::vector<int> indices;
	std::vector<Mesh_3_Vector_3> normals;
	std::vector<bool> realBoundary;

	std::vector<Mesh_3_Vector_3> vfVertices;
	std::vector<Mesh_3_Vector_3> vfFaces;
	std::vector<Mesh_3_Vector_3> vfNormals;
	std::vector<Mesh_3_Vector_3> vfColors; // face based
	std::vector<Mesh_3_Vector_3> vfvColors; // vertex based
	int vfvNum; // current total vertices in arrows.

	std::map<Mesh_3_Vertex_iterator, int> vertexBoundaryIdx; // components
	int numOfBoundaryComponents;

};