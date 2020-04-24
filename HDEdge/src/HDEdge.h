#pragma once

#include "My_Triangulation.h"
#include "RKIntegrator.h"

using namespace std;

class HDEdge
{
public:
	HDEdge() {};

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
	void initParameters(std::string sfx, double dr, double sr, int stps);
	void buildMeshFromSurface(std::string fn, double size);
	void buildLaplacian();
	void decompose();
	void visualize();
	void write();
	void integrate();

private:
	void buildLaplacian0();
	void buildLaplacian1();
	void buildLaplacian2();

	void computeBoundaryComponent();
	void setFormRandom();

	void convertForms();
	void computeCrossSection();
	void computeArrows();

	void writeVFasVTK(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, std::string sfx);
	
private:
	void computeHarmonicKnot();
	void computeHarmonicGradient();
	void computeFluxlessKnot();
	void computeGroundedGradient();
	void computeCurlyGradient();
	void computeCurlFreeGradient();

	void integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf, 
		double dsty, double len, int stps, std::string sfx);

private:
	void matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
		std::vector<double> row, std::vector<double> col, std::vector<double> val,
		std::vector<double> hs, double hm); // TODO: only evec is computed, fix eval if needed.
	
	void groupFields(std::vector<EigenVector> eigenFields);

private:
	// generate fields for omega
	void setPointChargeElectricField(EigenVector& field, Mesh_3_Point_3 rp, bool ori = true);
	void setCurrentMagneticField(EigenVector& field, Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, bool ori = true);
	void normalize(EigenVector& form);

	void boundaryToNoBoundary1Form(EigenVector& form);
	void noBoundaryToBoundary1Form(EigenVector& form);

	void boundaryToNoBoundary2Form(EigenVector& form);
	void noBoundaryToBoundary2Form(EigenVector& form);
	void computeScalarPotential(EigenVector form);

	void convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);
	void convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);

	void rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf);

	void selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors);
	void correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors);

	void assignColorStrength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode);
	void assignColorGradient();
	Mesh_3_Vector_3 assignColor(double ratio);
	double medianLength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void sortedIndices(std::vector<int>& sIndex, std::vector<double> len);

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

	EigenSpMat m_curlEnergyMatrix;
	EigenSpMat m_divEnergyMatrix;

private:
	EigenVector omega; // input form with boundary elements;
	vector<EigenVector> harmonicKnotBasis; // caution: For computation simplicity, these basis are without boundary elements
	vector<double> harmonicKnotCoeffs;
	vector<EigenVector> harmonicGradientPotentialBasis; // not orthogonal
	vector<EigenVector> harmonicGradientBasis; // orthogonal
	vector<double> harmonicGradientCoeffs;
	EigenVector harmonicKnot; // with boundary
	EigenVector harmonicGradient; // with boundary
	EigenVector fluxlessKnot; // with boundary
	EigenVector groundedGradient; //
	EigenVector curlyGradient;
	EigenVector curlFreeGradient; // 

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> omegaVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicKnotVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> harmonicGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> fluxlessKnotVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> groundedGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlyGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curlFreeGradientVF;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> targetVF;

	std::map<Mesh_3_Vertex_iterator, double> scalarPotential;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vectorPotential;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> omegaCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicKnotCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> harmonicGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> fluxlessKnotCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> groundedGradientCF;
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> curlyGradientCF;
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

	std::map<int, int> curlFieldIdx;
	std::map<int, int> divFieldIdx;
};
