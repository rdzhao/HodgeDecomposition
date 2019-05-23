#pragma once

#include "My_Triangulation.h"
#include "RKIntegrator.h"

using namespace std;

class LaplacianAnalyzer
{
public:
	LaplacianAnalyzer() {};

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
	void buildDECOperators();
	void buildLaplacian();
	void setEigenfield(int mode, int idx, bool isCurl);
	void decompose();
	void visualize();
	void integrate();

	void process_Poelke();

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
	void buildLaplacian3();

	void computeBoundaryComponent();
	void setForm();
	void setFormRandom();
	void computeHarmonic();
	void computeCoexact();
	void computeExact();

	void convertForms();
	void computeCrossSection();
	void computeArrows();

	void setEigenfield_L0n(int idx, bool isCurl);
	void setEigenfield_L0t(int idx, bool isCurl);
	void setEigenfield_L1n(int idx, bool isCurl);
	void setEigenfield_L1t(int idx, bool isCurl);
	void setEigenfield_L2n(int idx, bool isCurl);
	void setEigenfield_L2t(int idx, bool isCurl);
	void setEigenfield_L3n(int idx, bool isCurl);
	void setEigenfield_L3t(int idx, bool isCurl);

private:
	void computeHarmonicKnot();
	void computeHarmonicKnotTest();
	void computeHarmonicGradient();
	void computeFluxlessKnot();
	void computeGroundedGradient();
	void computeCurlyGradient();

	void computeCurlFreeGradient();
	void computeRemainingDivFreeFieldPotential();
	void testSPDHS();
	void testHodgeStar1();

private:	
	void assignBoundaryVerticesGroup();
	void assignVertexColorMixBoundary();
	void computeHarmonicMixBoundary();
	void computeExactMixBoundary();
	void computeCoexactMixBoundary();
	void computeCurlyGradientMixBoundary();
	void assembleLaplacian1MixBoundary();
	void reassembleOmega();

private:
	void assignBoundaryVertexPotential(int idx);

	void matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
		std::vector<double> row, std::vector<double> col, std::vector<double> val,
		std::vector<double> brow, std::vector<double> bcol, std::vector<double> bval,
		double hm); // TODO: only evec is computed, fix eval if needed.

	void integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf,
		std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf,
		double dsty, double lth, int stps, int md, std::string sfx);

	void prepareEigenFeilds(std::vector<EigenVector>& eigenFields);
	void computeCurlFields2Form(std::vector<EigenVector>& eigenFields);
	void convert2FormFieldTo1FormField(std::vector<EigenVector>& eigenfields1Form, std::vector<EigenVector> eigenfields2Form);

	void groupFields_L1n(std::vector<EigenVector> eigenFields);
	void groupFields_L1t(std::vector<EigenVector> eigenFields);
	void groupFields_L2n(std::vector<EigenVector> eigenFields);
	void groupFields_L2t(std::vector<EigenVector> eigenFields);

	void printEigenvalues(std::vector<double> eigenvalues, string sfx);

private:
	// generate fields for omega
	void setConstantField(EigenVector& field, Mesh_3_Vector_3 d);
	void setPointChargeElectricField(EigenVector& field, Mesh_3_Point_3 rp, bool ori = true);
	void setCurrentMagneticField(EigenVector& field, Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, bool ori = true);
	void smooth(EigenVector& form);

	void boundaryToNoBoundary1Form(EigenVector& form);
	void noBoundaryToBoundary1Form(EigenVector& form);

	void boundaryToNoBoundary2Form(EigenVector& form);
	void noBoundaryToBoundary2Form(EigenVector& form);
	void computeScalarPotential(EigenVector form);

	void convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);

	void computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci);

	void selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors);
	void correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors);

	void assignColorStrength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode);
	void assignColorGradient(EigenVector ev);
	Mesh_3_Vector_3 assignColor(double ratio);
	Mesh_3_Vector_3 assignColorSF(double ratio);
	double medianLength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void sortedIndices(std::vector<int>& sIndex, std::vector<double> len);

	void writeMesh();
	void writeMeshPLY();

private: // functions for Peokle method
	void solveCurl_Poelke();
	void solveGradient_Poelke();

	void makeCurlOperator_Poelke(EigenSpMat& curlOp);
	void makeGradOperator_Poelke(EigenSpMat& gradOp);

	void concatenateSPMatrix(EigenSpMat& M, EigenSpMat A, EigenSpMat B);
	void concatenateVector(EigenVector& V, EigenVector A, EigenVector B);

private:
	void writeTestMesh();

private:
	std::string m_suffix;
	double m_density_ratio;
	double m_step_ratio;
	double m_steps;

private:
	My_Triangulation m_mesh_3;

	EigenSpMat m_laplacian0_NB;
	EigenSpMat m_laplacian0_B;
	EigenSpMat m_laplacian1_B;
	EigenSpMat m_laplacian1_NB;
	EigenSpMat m_laplacian2_B;
	EigenSpMat m_laplacian2_NB;
	EigenSpMat m_laplacian3_B;
	EigenSpMat m_laplacian3_NB;

	EigenSpMat m_divEM_L1n;
	EigenSpMat m_curlEM_L1n;
	EigenSpMat m_divEM_L1t;
	EigenSpMat m_curlEM_L1t;
	EigenSpMat m_divEM_L2n;
	EigenSpMat m_curlEM_L2n;
	EigenSpMat m_divEM_L2t;
	EigenSpMat m_curlEM_L2t;

private:
	EigenVector omega; // input form with boundary elements;
	vector<EigenVector> harmonicKnotBasis; // caution: For computation simplicity, these basis are without boundary elements
	vector<double> harmonicKnotCoeffs;
	vector<EigenVector> harmonicGradientPotentialBasis; // not orthogonal
	vector<EigenVector> harmonicGradientBasis; // orthogonal
	vector<double> harmonicGradientCoeffs;
	EigenVector harmonicKnot; // with boundary
	EigenVector harmonicKnotTest; // with boundary
	EigenVector harmonicGradient; // with boundary
	EigenVector fluxlessKnot; // with boundary
	EigenVector groundedGradient; //
	EigenVector curlyGradient;


	EigenVector curlyGradientVP;

	EigenVector curlFreeGradient; // 
	EigenVector remainingDivFreeField; // 

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> omegaVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicKnotVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> fluxlessKnotVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> groundedGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlyGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlFreeGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> targetVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlyGradientVPVF;


	std::map<Mesh_3_Vertex_iterator, double> scalarPotential;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vectorPotential;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> omegaCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicKnotCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> fluxlessKnotCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> groundedGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> curlyGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> scalarFieldCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> colorField;

	EigenVector gamma_1;
	EigenVector gamma_2;
	EigenVector gamma_3;
	EigenVector gamma_4;

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> gamma_1_VF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> gamma_2_VF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> gamma_3_VF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> gamma_4_VF;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> gamma_1_CF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> gamma_2_CF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> gamma_3_CF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> gamma_4_CF;


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

	int curlSize;
	int divSize;
	std::map<int, int> curlFieldIdx;
	std::map<int, int> divFieldIdx;
};
