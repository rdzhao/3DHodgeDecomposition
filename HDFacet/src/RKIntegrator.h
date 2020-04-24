#pragma once

#include "My_Triangulation.h"

class PolyVertex
{
public:
	// barencentric coordinate to represent vertex
	Mesh_3_Cell_iterator ci;
	EigenVector3d bary;
	int lineTag;
	int nc; // num cell passed.
	Mesh_3_Vector_3 col;
};

class Poly
{
public:
	std::list<PolyVertex> vertices;
};

class TagNC
{
public:
	TagNC(int a, int b) : tag(a), nc(b) {}

	int tag;
	int nc;

};

typedef std::list<PolyVertex>::iterator PVIter;

class RKIntegrator
{
public:
	RKIntegrator(My_Triangulation& mesh):
		_mesh_3(mesh){};

	void init(double d, double l, int s, int m);
	void initColorField(std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> field);
	void convert(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	void integrate();
	//void clean();
	void write(std::string suffix = "default");

private:
	void _calculateVertexBasedVectorField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf);
	
	void _initLinesInCell();

	void _integrate();
	void _integrateFromCell(Mesh_3_Cell_iterator ci); // given a starting point

	void _clean();
	void _cleanShortLines();

	void _assignPolyVertexColor();

	void _writeStreamlines(std::string suffix);

private:
	bool _integrateAStep(PolyVertex rkv, double h, bool forward, PolyVertex& nrkv); // one step for a single line

	EigenVector3d _vertexGlobalToBarycentric(Mesh_3_Cell_iterator ci, Mesh_3_Point_3 p);
	Mesh_3_Point_3 _vertexBarycentricToGlobal(Mesh_3_Cell_iterator ci, EigenVector3d bary);
	EigenVector3d _vectorGlobalToBarycentric(Mesh_3_Cell_iterator ci, Mesh_3_Vector_3 v);
	Mesh_3_Vector_3 _vectorBarycentricToGlobal(Mesh_3_Cell_iterator ci, EigenVector3d bary);

	Mesh_3_Vector_3 _vectorInterpolation(Mesh_3_Cell_iterator ci, EigenVector3d bary); //given barycentric coordinate in cell
	EigenVector3d _baseTransformBarycentricCoordinate(Mesh_3_Cell_iterator ci, EigenVector3d inBary, int i);
	bool _computeIntersection(Mesh_3_Cell_iterator ci, EigenVector3d pBary, EigenVector3d vBary, double& th, int& i);
	bool _computeLocation(PolyVertex rkv, Mesh_3_Vector_3 v, PolyVertex& nrkv); // return the polyvertex after a translation v

	bool _checkNearBoundary(PolyVertex rkv);
	bool _checkLinePassed(PolyVertex rkv);

	double _length(Poly pl);

	double _distanceToFacet(Mesh_3_Point_3 p, Mesh_3_Facet_iterator fi);


private:
	My_Triangulation& _mesh_3;

	int _mode; // 0: long, 1: mid, 2: short
	int _spacing;
	double _step_length_ratio;
	int _step_size; // how many steps we do integration

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> colorField;

	std::vector<int> _line_density_in_cell;
	std::map<Mesh_3_Cell_iterator, std::vector<TagNC>> _lines_in_cell;
	std::map<Mesh_3_Cell_iterator, std::set<int>> _num_lines;

	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> _vertex_to_field;
	std::vector<Poly> _streamlines;
	std::vector<bool> _sl_valid; // whether streamlines are valid
};
