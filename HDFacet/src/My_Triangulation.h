#pragma once

#ifndef _MY_TRIANGULATION_
#define _MY_TRIANGULATION_

#define BOOST_PARAMETER_MAX_ARITY  12

#define BOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/make_mesh_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

#include <stack>
#include <bitset>
#include <unordered_map>
#include <functional>
#include <utility>
#include <algorithm>
#include <random>
#include <cmath>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

const double EPS = 0.000000000001;

// Domain 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//typedef CGAL::Polyhedron_3<K> Polyhedron;
//typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Polyhedron_mesh_domain;

typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron_with_feature;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain_with_feature;

// Data Type 
typedef K::Vector_3 Mesh_3_Vector_3;
typedef K::Point_3 Mesh_3_Point_3;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
//typedef CGAL::Mesh_triangulation_3<Polyhedron_mesh_domain>::type Tr;
typedef CGAL::Mesh_triangulation_3<Mesh_domain_with_feature, CGAL::Default, Concurrency_tag>::type Tr;
//typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
	Tr, Mesh_domain_with_feature::Corner_index,
	Mesh_domain_with_feature::Curve_segment_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// Iterator 
typedef Tr::Facet Facet;
typedef Tr::Edge Edge;
typedef Tr::Vertex_iterator Mesh_3_Vertex_iterator;
typedef Tr::Edge_iterator Mesh_3_Edge_iterator;
typedef Tr::Facet_iterator Mesh_3_Facet_iterator;
typedef Tr::Cell_iterator Mesh_3_Cell_iterator;
typedef Tr::Cell_circulator Mesh_3_Cell_circulator;
typedef Tr::Facet_circulator Mesh_3_Facet_circulator;

// Non duplicate treatment
typedef std::tuple<Mesh_3_Cell_iterator, int> My_facet;
typedef std::tuple<Mesh_3_Cell_iterator, int, int> My_edge;

// Eigen 
typedef Eigen::SparseMatrix<double> EigenSpMat;
typedef Eigen::MatrixXd EigenMatrix;
typedef Eigen::VectorXd EigenVector;
typedef Eigen::Matrix3d EigenMatrix3d;
typedef Eigen::Matrix4d EigenMatrix4d;
typedef Eigen::Vector3d EigenVector3d;
typedef Eigen::Vector4d EigenVector4d;

// To avoid verbose function and named parameters call
using namespace std;
using namespace CGAL::parameters;

class Mesh_3_Facet_Iterator_Comparator {
public:
	bool operator() (const Mesh_3_Facet_iterator& iter1, const Mesh_3_Facet_iterator& iter2) const {
		return *iter1 < *iter2;
	}
};

class Mesh_3_Edge_Iterator_Comparator {
public:
	bool operator() (const Mesh_3_Edge_iterator& iter1, const Mesh_3_Edge_iterator& iter2) const {
		return *iter1 < *iter2;
	}
};

class My_Triangulation {
public:
	My_Triangulation() {};
	~My_Triangulation() {};

	// accessors
	C3t3& C3T3();

	EigenSpMat& ED0F_B();
	EigenSpMat& ED1F_B();
	EigenSpMat& ED2F_B();
	EigenSpMat& ED0F_NB();
	EigenSpMat& ED1F_NB();
	EigenSpMat& ED2F_NB();
	
	EigenSpMat& HS0F_B();
	EigenSpMat& HS1F_B();
	EigenSpMat& HS2F_B();
	EigenSpMat& HS3F();
	EigenSpMat& HS0F_NB();
	EigenSpMat& HS1F_NB();
	EigenSpMat& HS2F_NB();

	int NumV();
	int NumE();
	int NumF();
	int NumC();
	int NumIV();
	int NumBV();
	int NumIE();
	int NumBE();
	int NumIF();
	int NumBF();

	int VertexIdx(Mesh_3_Vertex_iterator vi);
	int EdgeIdx(Mesh_3_Edge_iterator ei);
	int FacetIdx(Mesh_3_Facet_iterator fi);
	int CellIdx(Mesh_3_Cell_iterator ci);

	Mesh_3_Vertex_iterator IdxToVertex(int idx);
	Mesh_3_Edge_iterator IdxToEdge(int idx);
	Mesh_3_Facet_iterator IdxToFacet(int idx);
	Mesh_3_Cell_iterator IdxToCell(int idx);

	Mesh_3_Facet_iterator MyCellFacetMap(My_facet mf);
	Mesh_3_Edge_iterator MyCellEdgeMap(My_edge me);

	bool VertexValid(Mesh_3_Vertex_iterator vi);
	bool EdgeValid(Mesh_3_Edge_iterator ei);
	bool FacetValid(Mesh_3_Facet_iterator fi);

	bool VertexOnBoundary(Mesh_3_Vertex_iterator vi);
	bool EdgeOnBoundary(Mesh_3_Edge_iterator ei);
	bool FacetOnBoundary(Mesh_3_Facet_iterator fi);

	Mesh_3_Point_3 EdgeCC(Mesh_3_Edge_iterator ei);
	Mesh_3_Point_3 FacetCC(Mesh_3_Facet_iterator fi);
	Mesh_3_Point_3 CellCC(Mesh_3_Cell_iterator ci);

	double VertexDual(Mesh_3_Vertex_iterator vi);
	double EdgePrimal(Mesh_3_Edge_iterator ei);
	double EdgeDual(Mesh_3_Edge_iterator ei);
	double FacetPrimal(Mesh_3_Facet_iterator fi);
	double FacetDual(Mesh_3_Facet_iterator fi);
	double CellPrimal(Mesh_3_Cell_iterator ci);

	int VertexBoundaryIdx(Mesh_3_Vertex_iterator vi);
	int NumBoundary();
	int NumHMLG1();
	int NumHMLG2();

	int VertexForward(int idx);
	int EdgeForward(int idx);
	int FacetForward(int idx);

	int VertexBackward(int idx);
	int EdgeBackward(int idx);
	int FacetBackward(int idx);

	double getRadius();
	Mesh_3_Point_3 getCenter();
	double getAveLen();
	Mesh_3_Point_3 BoundingBoxMin();
	Mesh_3_Point_3 BoundingBoxMax();

	void preprocessing();

private:
	void indexElements();
	void basicInit();
	void computeValidAndBoundary();
	void computeBoundaryComponent();
	void computeHomology();
	void buildDECOperators();

private:
	void computeElementCircumcenter();
	void computeElementVolume();

	void buildMyCellFacetEdgeMap();
	
	void computeEdgeCircumcenter();
	void computeFacetCircumcenter();
	void computeCellCircumcenter();

	void computeVertexVolumes();
	void computeEdgeVolumes();
	void computeFacetVolumes();
	void computeCellVolumes();

	void createForwardBackwardMap();

private:
	void buildHS0F_B();
	void buildHS1F_B();
	void buildHS2F_B();
	void buildHS3F();
	void buildHS0F_NB();
	void buildHS1F_NB();
	void buildHS2F_NB();

	void buildED0F_B();
	void buildED1F_B();
	void buildED2F_B();
	void buildED0F_NB();
	void buildED1F_NB();
	void buildED2F_NB();

private:
	// criteria
	double m_facet_angle;
	double m_facet_size;
	double m_facet_distance;
	double m_cell_radius_edge_ratio;
	double m_cell_size;

	// model
	std::vector<float> voxels;
	C3t3 m_c3t3; // tet mesh

	// mesh info
	Mesh_3_Point_3 center;
	Mesh_3_Point_3 bbMin;
	Mesh_3_Point_3 bbMax;
	double m_radius; // 
	double m_aveLen; //average edge length

	int nhmlg_1; // number of homolgy 1 generators
	int nhmlg_2; // number of relative homology 1 generator

	std::map<Mesh_3_Vertex_iterator, int> vertexIdx;
	std::map<Mesh_3_Edge_iterator, int, Mesh_3_Edge_Iterator_Comparator> edgeIdx;
	std::map<Mesh_3_Facet_iterator, int, Mesh_3_Facet_Iterator_Comparator> facetIdx;
	std::map<Mesh_3_Cell_iterator, int> cellIdx;

	std::map<int, Mesh_3_Vertex_iterator> idxToVertex;
	std::map<int, Mesh_3_Edge_iterator> idxToEdge;
	std::map<int, Mesh_3_Facet_iterator> idxToFacet;
	std::map<int, Mesh_3_Cell_iterator> idxToCell;

	std::map<Mesh_3_Vertex_iterator, bool> vertexValid;
	std::map<Mesh_3_Edge_iterator, bool, Mesh_3_Edge_Iterator_Comparator> edgeValid;
	std::map<Mesh_3_Facet_iterator, bool, Mesh_3_Facet_Iterator_Comparator> facetValid;

	std::map<Mesh_3_Vertex_iterator, bool> vertexOnBoundary;
	std::map<Mesh_3_Edge_iterator, bool, Mesh_3_Edge_Iterator_Comparator> edgeOnBoundary;
	std::map<Mesh_3_Facet_iterator, bool, Mesh_3_Facet_Iterator_Comparator> facetOnBoundary;
	std::map<Mesh_3_Cell_iterator, bool> cellOnBoundary;

	std::map<My_facet, Mesh_3_Facet_iterator> myCellFacetMap;
	std::map<My_edge, Mesh_3_Edge_iterator> myCellEdgeMap;

	// circumcenter
	std::map<Mesh_3_Edge_iterator, Mesh_3_Point_3, Mesh_3_Edge_Iterator_Comparator> edgeCC;
	std::map<Mesh_3_Facet_iterator, Mesh_3_Point_3, Mesh_3_Facet_Iterator_Comparator> facetCC;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Point_3> cellCC;

	// element volume
	std::map<Mesh_3_Vertex_iterator, double> vertexDual;
	std::map<Mesh_3_Edge_iterator, double, Mesh_3_Edge_Iterator_Comparator> edgePrimal;
	std::map<Mesh_3_Edge_iterator, double, Mesh_3_Edge_Iterator_Comparator> edgeDual;
	std::map<Mesh_3_Facet_iterator, double, Mesh_3_Facet_Iterator_Comparator> facetPrimal;
	std::map<Mesh_3_Facet_iterator, double, Mesh_3_Facet_Iterator_Comparator> facetDual;
	std::map<Mesh_3_Cell_iterator, double> cellPrimal;

	// boundary 
	std::map<Mesh_3_Vertex_iterator, int> vertexBoundaryIdx; // components
	int numBoundary;

	// discrete exterior derivatibe
	EigenSpMat m_ED0F_B;
	EigenSpMat m_ED1F_B;
	EigenSpMat m_ED2F_B;
	EigenSpMat m_ED0F_NB;
	EigenSpMat m_ED1F_NB;
	EigenSpMat m_ED2F_NB;

	// discrete hodge star
	EigenSpMat m_HS0F_B;
	EigenSpMat m_HS1F_B;
	EigenSpMat m_HS2F_B;
	EigenSpMat m_HS3F;
	EigenSpMat m_HS0F_NB;
	EigenSpMat m_HS1F_NB;
	EigenSpMat m_HS2F_NB;

private:
	std::map<int, int> vertexForward;
	std::map<int, int> edgeForward;
	std::map<int, int> facetForward;

	std::map<int, int> vertexBackward;
	std::map<int, int> edgeBackward;
	std::map<int, int> facetBackward;

	int numV;
	int numE;
	int numF;
	int numC;

	int numIV; // interior vertex
	int numIE; // interior edge
	int numIF; // interior face

	int numBV; // boundary vertex
	int numBE; // boundary edge
	int numBF; // boundary facet
};

#endif
