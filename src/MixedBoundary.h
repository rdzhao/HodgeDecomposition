#pragma once

#include "My_Triangulation.h"
#include "RKIntegrator.h"

using namespace std;

class MixedBoundary
{
public:
	MixedBoundary() {};
	 
	My_Triangulation& Mesh();
	std::vector<Mesh_3_Point_3>& Vertices();
	std::vector<int>& Indices();
	std::vector<Mesh_3_Vector_3>& Normals();
	std::vector<Mesh_3_Vector_3>& Colors();
	std::vector<Mesh_3_Vector_3>& VfVertices();
	std::vector<Mesh_3_Vector_3>& VfFaces();
	std::vector<Mesh_3_Vector_3>& VfNormals();
	std::vector<Mesh_3_Vector_3>& VfColors();

	// major functions
	void init(double dsty, double sr, int stps, std::string sfx);
	void buildMeshFromSurface(std::string fn, double size);
	void readSimulationData(std::string fn_tet, std::string fn_field);
	void prepare();
	void decompose();
	void visualize();
	void integrate();

private:
	void readField(std::string fn_field);

	void assignMixedElements();
	void assignMixedForwardBackward();
	void createOperators();

	void setOmega();
	void setOmegaFromPWCF();
	void computeHarmonic();
	void computeGradient();
	void computeCurl();

	void convertForms();
	void computeCrossSection();
	void computeArrows();

	void writeMesh();
	void writeMeshPLY();

private:
	void buildExteriorDerivative0Form();
	void buildExteriorDerivative1Form();
	void buildExteriorDerivative2Form();

	void buildHodgeStar0Form();
	void buildHodgeStar1Form();
	void buildHodgeStar2Form();
	void buildHodgeStar3Form();

	void buildLaplacian0();
	void buildLaplacian1();
	void buildLaplacian2();

private:
	EigenVector matlabLinearSolveCholesky(EigenSpMat A, EigenVector b);
	void matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
		std::vector<double> row, std::vector<double> col, std::vector<double> val,
		std::vector<double> brow, std::vector<double> bcol, std::vector<double> bval,
		double hm); // TODO: only evec is computed, fix eval if needed.

	void integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf,
		std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf,
		double dsty, double lth, int stps, int md, std::string sfx);

private:
	void mixedToFull1Form(EigenVector& form);
	void fulltoMixed1Form(EigenVector& form);

	void mixedToFull0Form(EigenVector& form);

	void convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);

	void computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci);

	void selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors);
	void correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors, EigenSpMat hs, std::vector<EigenVector> ef);

	void sortedIndices(std::vector<int>& sIndex, std::vector<double> len);
	Mesh_3_Vector_3 assignColorSF(double ratio);
	Mesh_3_Vector_3 assignColorVF(double ratio);
	void assignColorGradient(EigenVector ev);
	void assignColorVectorField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode);

private: //debug
	void assignVertexColorMixBoundary();
	void computeNorms();

private:
	std::string m_suffix;
	double m_density_ratio;
	double m_step_ratio;
	double m_steps;

private:
	My_Triangulation m_mesh_3;

	std::map<Mesh_3_Vertex_iterator, bool> vertexMixed;
	std::map<Mesh_3_Edge_iterator, bool, Mesh_3_Edge_Iterator_Comparator> edgeMixed;
	std::map<Mesh_3_Facet_iterator, bool, Mesh_3_Facet_Iterator_Comparator> facetMixed;
	
	std::map<int, int> vertexForwardMixed;
	std::map<int, int> edgeForwardMixed;
	std::map<int, int> facetForwardMixed;

	std::map<int, int> vertexBackwardMixed;
	std::map<int, int> edgeBackwardMixed;
	std::map<int, int> facetBackwardMixed;

	int numKV;
	int numKE;
	int numKF;
	int numMV;
	int numME;
	int numMF;

private:
	EigenSpMat m_ed0f;
	EigenSpMat m_ed1f;
	EigenSpMat m_ed2f;

	EigenSpMat m_hs0f;
	EigenSpMat m_hs1f;
	EigenSpMat m_hs2f;
	EigenSpMat m_hs3f;

	EigenSpMat m_laplacian0;
	EigenSpMat m_laplacian1;
	EigenSpMat m_laplacian2;

	EigenSpMat m_divEM;
	EigenSpMat m_curlEM;

private:
	EigenVector omega;
	vector<EigenVector> harmonicKnotBasis;
	EigenVector harmonic;
	EigenVector gradient;
	EigenVector curl;

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> omegaVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> gradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlVF;

	std::map<Mesh_3_Vertex_iterator, double> scalarPotential;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vectorPotential;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> omegaCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> gradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> curlCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> colorField;

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

private:
	std::map<Mesh_3_Vertex_iterator, double> boundaryVertexPotential;

private:
	std::map<Mesh_3_Vertex_iterator, int> vertexGroups; // for mix boundary condition

	std::map<int, int> curlFieldIdx;
	std::map<int, int> divFieldIdx;
};
