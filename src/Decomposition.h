#pragma once

#include "My_Triangulation.h"
#include "RKIntegrator.h"
#include "JHTDB.h"

using namespace std;

class Decomposition
{
public:
	Decomposition() {};

	// accessor
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
	void decompose();
	void visualize();
	void integrate();
	void write();

public: //debug functions
	void testSurface();
	void testSurfaceAndArrowsPLY();

private:
	void buildLaplacian1();
	void buildLaplacian2();
	void buildLaplacian3();

	void setForm();
	void setFormRandom();
	void sampleFromJHTDB();

	void convertForms();
	void computeCrossSection(int mode);
	void computeArrows();

	void writeMesh();
	void writeMeshPLY(std::string sfx);

private:
	void computeHarmonicKnot();
	void computeHarmonicGradient();
	void computeFluxlessKnot();
	void computeGroundedGradient();
	void computeCurlyGradient();

	void computeCurlFreeGradient();

	void computeCurlyGradientSP();
	void computeCurlyGradientVP();
	void computeHarmonicKnotVP();
	
	void computeNewtonianDivField();

private:
	void computeNG();
	void computeTC();
	void computeNH();
	void computeTH();
	void computeCH();

private:
	EigenVector matlabLinearSolveCholesky(EigenSpMat A, EigenVector b);
	void matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
		std::vector<double> row, std::vector<double> col, std::vector<double> val,
		std::vector<double> brow, std::vector<double> bcol, std::vector<double> bval,
		double hm); // TODO: only evec is computed, fix eval if needed.

	void assignBoundaryFacetPotential(int idx);

	void integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf,
			std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf,
			double dsty, double lth, int stps, int md, std::string sfx);

	void groupFields(std::vector<EigenVector> eigenFields);

private:
	// generate fields for omega
	void setPointChargeElectricField(EigenVector& field, Mesh_3_Point_3 rp, bool ori = true);
	void setCurrentMagneticField(Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, EigenVector& field);
	void setConstantField(Mesh_3_Vector_3 d, EigenVector& field);
	void smooth(EigenVector& form);

	void boundaryToNoBoundary1Form(EigenVector& form);
	void noBoundaryToBoundary1Form(EigenVector& form);

	void boundaryToNoBoundary2Form(EigenVector& form);
	void noBoundaryToBoundary2Form(EigenVector& form);

	void convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void convert3Form(EigenVector form, std::map<Mesh_3_Vertex_iterator, double>& sf);
	
	void computeScalarPotential(EigenVector form);
	void assignHarmonicGradientCF();

	void computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci);

	void selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors);
	void correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors, EigenSpMat hs, std::vector<EigenVector> eigenfields);

	void assignColorStrength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode);
	void assignColorStrengthVP(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode);
	void assignColorStrengthSP(std::map<Mesh_3_Vertex_iterator, double> sf, int mode);
	
	Mesh_3_Vector_3 assignColor(double ratio);
	Mesh_3_Vector_3 assignColorVP(double ratio);
	Mesh_3_Vector_3 assignColorSP(double ratio);
	void sortedIndices(std::vector<int>& sIndex, std::vector<double> len);

	void testOrthogonality();
	void evaluateAccuracy();

	void correctBoundaryPotential(std::map<Mesh_3_Vertex_iterator, double>& sf);

private:
	void createDivTan(EigenVector& form);

private:
	void sampleInput();
	void computeError();
	void computeL2Norm();

	// debugging
	void testEnergy();
	void testNonDiagonalHodgeStar();
	void testProjectionMatrix();

private:
	// analytical functions
	Mesh_3_Vector_3 normalGradient(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 tangentialCurl(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 normalHarmonic(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 tangentialHarmonic(Mesh_3_Vector_3 p);
	Mesh_3_Vector_3 centralHarmonic(Mesh_3_Vector_3 p);

	double fluxNormalGradient(Mesh_3_Facet_iterator fi);
	double fluxTangentialCurl(Mesh_3_Facet_iterator fi);
	double fluxNormalHarmonic(Mesh_3_Facet_iterator fi);
	double fluxTangentialHarmonic(Mesh_3_Facet_iterator fi);
	double fluxCentralHarmonic(Mesh_3_Facet_iterator fi);

	EigenVector sampleField(int fieldType);
	double sampleFacetFlux(Mesh_3_Facet_iterator fi,int fieldType, int reso);
	void generateBaryCoords(std::vector<Mesh_3_Vector_3>& coords, int reso);

private:
	std::string m_suffix;
	double m_density_ratio;
	double m_step_ratio;
	double m_steps;

private:
	My_Triangulation m_mesh_3;

	EigenSpMat m_laplacian1_NB;
	EigenSpMat m_laplacian1_B;
	EigenSpMat m_laplacian2_NB;
	EigenSpMat m_laplacian2_B;
	EigenSpMat m_laplacian3_NB;
	EigenSpMat m_laplacian3_B;

	EigenSpMat m_divEM;
	EigenSpMat m_curlEM;

private:
	EigenVector noHarmonicComponent;

	EigenVector omega; // input form with boundary elements;
	EigenVector omegaOriginal;
	vector<EigenVector> harmonicKnotBasis; // caution: For computation simplicity, these basis are without boundary elements
	vector<double> harmonicKnotCoeffs;
	vector<EigenVector> harmonicGradientPotentialBasis; // not orthogonal
	vector<EigenVector> harmonicGradientBasis; // orthogonal
	vector<double> harmonicGradientCoeffs;
	EigenVector harmonicKnot; // with boundary
	EigenVector harmonicGradient; // with boundary
	EigenVector fluxlessKnot; // with boundary
	EigenVector groundedGradient;
	EigenVector curlyGradient;

	EigenVector harmonicGradientSP;
	EigenVector groundedGradientSP;
	EigenVector curlyGradientSP;
	EigenVector harmonicKnotVP;
	EigenVector fluxlessKnotVP;
	EigenVector curlyGradientVP;

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> omegaVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicKnotVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> fluxlessKnotVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> groundedGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlyGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlFreeGradientVF; /// !!!!!!!

	std::map<Mesh_3_Vertex_iterator, double> harmonicGradientSPSF;
	std::map<Mesh_3_Vertex_iterator, double> groundedGradientSPSF;
	std::map<Mesh_3_Vertex_iterator, double> curlyGradientSPSF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicKnotVPVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> fluxlessKnotVPVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlyGradientVPVF;

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> targetVF;

	std::map<Mesh_3_Vertex_iterator, double> scalarPotential;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vectorPotential;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> omegaCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicKnotCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> fluxlessKnotCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> groundedGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> curlyGradientCF;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicGradientSPCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> groundedGradientSPCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> curlyGradientSPCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicKnotVPCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> fluxlessKnotVPCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> curlyGradientVPCF;

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

private:
	std::map<Mesh_3_Facet_iterator, double, Mesh_3_Facet_Iterator_Comparator> boundaryFacetPotential;

	std::map<int, int> curlFieldIdx;
	std::map<int, int> divFieldIdx;

	double m_accu;
	double m_cond;

};