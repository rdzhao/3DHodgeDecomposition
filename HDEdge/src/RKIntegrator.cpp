#include "RKIntegrator.h"

void RKIntegrator::init(double d, double lr, int s)
{
	// d is streamline density ranging in (0, 1) for determining _spacing.
	// total number of cells times this parameter
	// 1/d almost means how many integral lines
	_spacing = static_cast<int>(ceil(_mesh_3.NumC()*d));
	// lr is length ratio for detemrining step length compared with average mesh edge length
	_step_length_ratio = lr;
	_step_size = s;
}

void RKIntegrator::initColorField(std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3> field)
{
	colorField = field;
}

void RKIntegrator::convert(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	_calculateVertexBasedVectorField(vf);
}

void RKIntegrator::integrate()
{
	_initLinesInCell();
	_integrate();
	_clean();
}

void RKIntegrator::write(std::string suffix)
{
	_writeStreamlines(suffix);
}

void RKIntegrator::_calculateVertexBasedVectorField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	for (Mesh_3_Vertex_iterator vi = _mesh_3.C3T3().triangulation().vertices_begin();
		vi != _mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (_mesh_3.C3T3().triangulation().is_infinite(vi))
			continue;

		_vertex_to_field[vi] = Mesh_3_Vector_3(0, 0, 0);
		double vol = 0;
		int count = 0;

		std::set<Mesh_3_Cell_iterator> visitedCells;
		std::stack<Mesh_3_Cell_iterator> NotVisitedCells;

		NotVisitedCells.push(vi->cell());

		while (!NotVisitedCells.empty()) {
			Mesh_3_Cell_iterator ci = NotVisitedCells.top();
			visitedCells.insert(ci);
			NotVisitedCells.pop();

			if (ci->subdomain_index() == 0) {
				int idx;
				if (!ci->has_vertex(vi, idx)) {
					std::cout << "Error: Cell is not retrived correctly from a given vertex" << std::endl;
				}

				// push adjacent unvisited cells
				for (int i = 0; i < 4; ++i) {
					if (i != idx && visitedCells.find(ci->neighbor(i)) == visitedCells.end()) {
						NotVisitedCells.push(ci->neighbor(i));
					}
				}
			}
			else {
				int idx;
				if (!ci->has_vertex(vi, idx)) {
					std::cout << "Error: Cell is not retrived correctly from a given vertex" << std::endl;
				}

				_vertex_to_field[vi] += vf[ci] * _mesh_3.CellPrimal(ci);
				vol += _mesh_3.CellPrimal(ci);
				++count;

				// push adjacent unvisited cells
				for (int i = 0; i < 4; ++i) {
					if (i != idx && visitedCells.find(ci->neighbor(i)) == visitedCells.end()) {
						NotVisitedCells.push(ci->neighbor(i));
					}
				}
			}
		}

		_vertex_to_field[vi] /= vol;
	
		//cout<< _vertex_to_field[vi].x()<<" " << _vertex_to_field[vi].y() << " " << _vertex_to_field[vi].z()<<endl;
	}
}

void RKIntegrator::_initLinesInCell()
{
	// fill density
	for (Mesh_3_Cell_iterator ci = _mesh_3.C3T3().triangulation().cells_begin();
		ci != _mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		std::vector<TagNC> ts;
		_lines_in_cell[ci] = ts;

		std::set<int> ss;
		_num_lines[ci] = ss;
	}
}

void RKIntegrator::_integrate()
{
	for (Mesh_3_Cell_iterator ci = _mesh_3.C3T3().triangulation().cells_begin();
		ci != _mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		if (_mesh_3.CellIdx(ci) % _spacing != 0)
			continue;

		if (!_lines_in_cell[ci].empty())
			continue;
		//cout << "cell: " << _mesh_3.CellIdx(ci) << endl;
		_integrateFromCell(ci);
	}

	_assignPolyVertexColor();
}

void RKIntegrator::_integrateFromCell(Mesh_3_Cell_iterator ci)
{
	Poly pl;

	PolyVertex rkv, nrkv;
	rkv.ci = ci;
	rkv.bary = EigenVector3d(0.25, 0.25, 0.25);
	rkv.lineTag = _streamlines.size();
	rkv.nc = 0;
	pl.vertices.push_back(rkv);
	
	//cout << pl.vertices.size() << endl;
	_lines_in_cell[rkv.ci].push_back(TagNC(rkv.lineTag, rkv.nc));
	_num_lines[rkv.ci].insert(rkv.lineTag);

	// forward
	bool terminate = false;
	bool terminateNearBoundary = false;
	bool terminateOverlap = false;
	int numc = 0;
	int i = 0;

	do {
		bool success = false;
		success = _integrateAStep(rkv, _step_length_ratio, true, nrkv);
		if (success) {
			if (rkv.ci != nrkv.ci)
				++numc;

			rkv = nrkv;
			rkv.lineTag = _streamlines.size();
			rkv.nc = numc;
			pl.vertices.push_back(rkv);

			terminateNearBoundary = _checkNearBoundary(rkv);
			terminateOverlap = _checkLinePassed(rkv);
			if (terminateNearBoundary || terminateOverlap) {
				//cout << terminateNearBoundary << " " << terminateOverlap << endl;			
				terminate = true;
			}

			_lines_in_cell[rkv.ci].push_back(TagNC(rkv.lineTag, rkv.nc));
			_num_lines[rkv.ci].insert(rkv.lineTag);
			++i;
		}
		else {
			terminate = true;
		}

	} while (!terminate && i < _step_size);

	// backward
	if (pl.vertices.front().ci != pl.vertices.back().ci) {
		rkv = pl.vertices.front();
		terminate = false;
		numc = 0;
		i = 0;
		do {
			bool success = false;
			success = _integrateAStep(rkv, _step_length_ratio, false, nrkv);

			if (success) {
				if (rkv.ci != nrkv.ci)
					--numc;

				rkv = nrkv;
				rkv.lineTag = _streamlines.size();
				rkv.nc = numc;
				pl.vertices.push_front(rkv);

				terminateNearBoundary = _checkNearBoundary(rkv);
				terminateOverlap = _checkLinePassed(rkv);
				if (terminateNearBoundary || terminateOverlap) {
					terminate = true;
				}

				_lines_in_cell[rkv.ci].push_back(TagNC(rkv.lineTag, rkv.nc));
				_num_lines[rkv.ci].insert(rkv.lineTag);
				++i;
			}
			else {
				terminate = true;
			}
		} while (!terminate && i < _step_size);
	}

	_streamlines.push_back(pl);
}

void RKIntegrator::_clean()
{
	_cleanShortLines();
}

void RKIntegrator::_cleanShortLines()
{
	std::vector<Poly>::iterator iter = _streamlines.begin();

	do {
		if (_length(*iter) < 2 * _mesh_3.getAveLen()) {
			_streamlines.erase(iter);
		}
		else {
			++iter;
		}
	} while (iter != _streamlines.end());
}

void RKIntegrator::_assignPolyVertexColor()
{
	for (int i = 0; i < _streamlines.size(); ++i) {
		for (auto iter = _streamlines[i].vertices.begin(); iter != _streamlines[i].vertices.end(); ++iter) {
			Mesh_3_Cell_iterator ci = iter->ci;
			EigenVector3d bary = iter->bary;

			std::vector<Mesh_3_Vector_3> colvecs;
			for (int i = 0; i < 4; ++i) {
				colvecs.push_back(colorField[ci->vertex(i)]);
			}

			Mesh_3_Vector_3 v = bary.x()*colvecs[0]
				+ bary.y()*colvecs[1]
				+ bary.z()*colvecs[2]
				+ (1 - bary.x() - bary.y() - bary.z())*colvecs[3];

			iter->col = v;
		}
	}
}

void RKIntegrator::_writeStreamlines(std::string suffix)
{
	std::string fn = "RKStreamlines_" + suffix + ".txt";

	std::ofstream out(fn.c_str());

	int idx = 0;
	for (int i = 0; i < _streamlines.size(); ++i) {

		//cout << _streamlines[i].vertices.size() << endl;

		for (auto iter = _streamlines[i].vertices.begin(); iter != _streamlines[i].vertices.end(); ++iter) {

			Mesh_3_Point_3 p = _vertexBarycentricToGlobal(iter->ci, iter->bary);

			out << idx << " "
				<< (p.x() /*- _mesh_3.getCenter().x()*/) / 1/*_mesh_3.getRadius()*/ << " "
				<< (p.y() /*- _mesh_3.getCenter().y()*/) / 1/*_mesh_3.getRadius()*/ << " "
				<< (p.z() /*- _mesh_3.getCenter().z()*/) / 1/*_mesh_3.getRadius()*/ << " "
				<< floor(iter->col.x() * 255) << " "
				<< floor(iter->col.y() * 255) << " "
				<< floor(iter->col.z() * 255) << std::endl;
		}

		++idx;
	}
}

bool RKIntegrator::_integrateAStep(PolyVertex rkv, double h, bool forward, PolyVertex& nrkv)
{
	bool success = false;

	if (forward) {
		bool valid;

		Mesh_3_Vector_3 v1 = _vectorInterpolation(rkv.ci, rkv.bary);

		if (sqrt(v1*v1) < 0.001*_mesh_3.getAveLen()) {
			success = false;
		}
		else {
			Mesh_3_Vector_3 k1 = h * v1;

			PolyVertex nrkv_2;
			valid = _computeLocation(rkv, k1 / 2, nrkv_2);

			//if (nrkv_2.ci->subdomain_index() == 0)
				//valid = false;

			if (valid) {
				Mesh_3_Vector_3 v2 = _vectorInterpolation(nrkv_2.ci, nrkv_2.bary);

				Mesh_3_Vector_3 k2 = h * v2;

				valid = _computeLocation(rkv, k2, nrkv);

				//if (nrkv.ci->subdomain_index() == 0)
					//valid = false;

				if (valid)
					success = true;
			}
		}
	}
	else {
		bool valid;

		Mesh_3_Vector_3 v1 = _vectorInterpolation(rkv.ci, rkv.bary);

		if (sqrt(v1*v1) < 0.001*_mesh_3.getAveLen()) {
			success = false;
		}
		else {
			Mesh_3_Vector_3 k1 = -h * v1;

			PolyVertex nrkv_2;
			valid = _computeLocation(rkv, k1 / 2, nrkv_2);

			//if (nrkv_2.ci->subdomain_index() == 0)
				//valid = false;

			if (valid) {
				Mesh_3_Vector_3 v2 = _vectorInterpolation(nrkv_2.ci, nrkv_2.bary);

				Mesh_3_Vector_3 k2 = -h * v2;

				valid = _computeLocation(rkv, k2, nrkv);

				//if (nrkv.ci->subdomain_index() == 0)
					//valid = false;

				if (valid)
					success = true;
			}
		}
	}

	return success;
}

EigenVector3d RKIntegrator::_vertexGlobalToBarycentric(Mesh_3_Cell_iterator ci, Mesh_3_Point_3 p)
{
	std::vector<Mesh_3_Point_3> points;
	for (int i = 0; i < 4; ++i) {
		points.push_back(ci->vertex(i)->point().point());
	}

	EigenMatrix3d A;
	A(0, 0) = points[0].x() - points[3].x(); A(0, 1) = points[1].x() - points[3].x(); A(0, 2) = points[2].x() - points[3].x();
	A(1, 0) = points[0].y() - points[3].y(); A(1, 1) = points[1].y() - points[3].y(); A(1, 2) = points[2].y() - points[3].y();
	A(2, 0) = points[0].z() - points[3].z(); A(2, 1) = points[1].z() - points[3].z(); A(2, 2) = points[2].z() - points[3].z();

	EigenVector3d b;
	b(0) = p.x() - points[3].x();
	b(1) = p.y() - points[3].y();
	b(2) = p.z() - points[3].z();

	EigenVector3d bary = A.inverse()*b;

	return bary;
}

Mesh_3_Point_3 RKIntegrator::_vertexBarycentricToGlobal(Mesh_3_Cell_iterator ci, EigenVector3d bary)
{
	std::vector<Mesh_3_Vector_3> vecs;
	for (int i = 0; i < 4; ++i) {
		vecs.push_back(ci->vertex(i)->point().point() - CGAL::ORIGIN);
	}

	Mesh_3_Vector_3 v = bary.x()*vecs[0] + bary.y()*vecs[1] + bary.z()*vecs[2] + (1 - bary.x() - bary.y() - bary.z())*vecs[3];

	Mesh_3_Point_3 p = CGAL::ORIGIN + v;

	return p;
}

EigenVector3d RKIntegrator::_vectorGlobalToBarycentric(Mesh_3_Cell_iterator ci, Mesh_3_Vector_3 v)
{
	std::vector<Mesh_3_Point_3> points;
	for (int i = 0; i < 4; ++i) {
		points.push_back(ci->vertex(i)->point().point());
	}

	Mesh_3_Point_3 p = points[0] + v;

	EigenMatrix3d A;
	A(0, 0) = points[0].x() - points[3].x(); A(0, 1) = points[1].x() - points[3].x(); A(0, 2) = points[2].x() - points[3].x();
	A(1, 0) = points[0].y() - points[3].y(); A(1, 1) = points[1].y() - points[3].y(); A(1, 2) = points[2].y() - points[3].y();
	A(2, 0) = points[0].z() - points[3].z(); A(2, 1) = points[1].z() - points[3].z(); A(2, 2) = points[2].z() - points[3].z();

	EigenVector3d b;
	b(0) = p.x() - points[3].x();
	b(1) = p.y() - points[3].y();
	b(2) = p.z() - points[3].z();

	EigenVector3d bary = A.inverse()*b;
	bary.x() = bary.x() - 1;

	return bary;
}

Mesh_3_Vector_3 RKIntegrator::_vectorBarycentricToGlobal(Mesh_3_Cell_iterator ci, EigenVector3d bary)
{
	bary.x() += 1;

	std::vector<Mesh_3_Vector_3> vecs;
	for (int i = 0; i < 4; ++i) {
		vecs.push_back(ci->vertex(i)->point().point() - CGAL::ORIGIN);
	}

	Mesh_3_Vector_3 v = bary.x()*vecs[0] + bary.y()*vecs[1] + bary.z()*vecs[2] + (1 - bary.x() - bary.y() - bary.z())*vecs[3];
	v -= vecs[0];

	return v;
}

Mesh_3_Vector_3 RKIntegrator::_vectorInterpolation(Mesh_3_Cell_iterator ci, EigenVector3d bary)
{
	Mesh_3_Vector_3 v(0, 0, 0);

	v += bary.x()*_vertex_to_field[ci->vertex(0)];
	v += bary.y()*_vertex_to_field[ci->vertex(1)];
	v += bary.z()*_vertex_to_field[ci->vertex(2)];
	v += (1 - bary.x() - bary.y() - bary.z())*_vertex_to_field[ci->vertex(3)];

	return v;
}

EigenVector3d RKIntegrator::_baseTransformBarycentricCoordinate(Mesh_3_Cell_iterator ci, EigenVector3d inBary, int i)
{
	EigenVector3d outBary(0, 0, 0);
	Mesh_3_Cell_iterator ci_neighbor = ci->neighbor(i);

	for (int k = 0; k < 4; ++k) {
		if (k != i) {
			double val;
			if (k == 0) {
				val = inBary.x();
			}
			else if (k == 1) {
				val = inBary.y();
			}
			else if (k == 2) {
				val = inBary.z();
			}
			else {
				val = 1 - inBary.x() - inBary.y() - inBary.z();
			}

			Mesh_3_Vertex_iterator vi = ci->vertex(k);

			int nk;
			ci_neighbor->has_vertex(vi, nk);

			if (nk == 0) {
				outBary.x() = val;
			}
			else if (nk == 1) {
				outBary.y() = val;
			}
			else if (nk == 2) {
				outBary.z() = val;
			}
		}
	}

	return outBary;
}

bool RKIntegrator::_checkNearBoundary(PolyVertex rkv)
{
	// check near boundary

	bool terminate = false;

	//if (rkv.ci->subdomain_index() == 0)
		//terminate = true;

	bool atBoundary = false;
	Mesh_3_Facet_iterator bfi; // boundary facet
	Mesh_3_Cell_iterator ci = rkv.ci;
	for (int i = 0; i < 4; ++i) {
		My_facet mf(ci, i);
		Mesh_3_Facet_iterator fi = _mesh_3.MyCellFacetMap(mf);

		if (_mesh_3.FacetOnBoundary(fi)) {
			atBoundary = true;
			//terminate = true;
			bfi = fi;
		}
	}

	if (atBoundary) {

		Mesh_3_Point_3 p = _vertexBarycentricToGlobal(rkv.ci, rkv.bary);
		double d = _distanceToFacet(p, bfi);

		if (d < _mesh_3.getAveLen()*0.5) {
			terminate = true;
		}
	}


	return terminate;
}

bool RKIntegrator::_checkLinePassed(PolyVertex rkv)
{
	bool terminate = false;

	Mesh_3_Cell_iterator ci = rkv.ci;

	for (int i = 0; i < _lines_in_cell[ci].size(); ++i) {
		if(_lines_in_cell[ci][i].tag == rkv.lineTag)
			if(_lines_in_cell[ci][i].nc != rkv.nc)
		{
				terminate = true;
				break;
			}
	}

	if (_num_lines[ci].size() > 4)
		terminate = true;

	return terminate;
}

bool RKIntegrator::_computeIntersection(Mesh_3_Cell_iterator ci, EigenVector3d pBary, EigenVector3d vBary, double& th, int& i)
{
	bool found = false;

	double ch = -1;

	if (fabs(vBary.x()) > 0 && pBary.x() != 0) {
		double h = -pBary.x() / vBary.x();
		double s = pBary.y() + h * vBary.y();
		double t = pBary.z() + h * vBary.z();
		double w = 1 - s - t;

		//std::cout << "In 1: " << s << " " << t << " " << w << " " << h << std::endl;
		if (s > 0 && t > 0 && w > 0 && h > 0) {
			if (h > ch) {
				ch = h;
				th = h;
				i = 0;
				found = true;
			}
		}
	}

	if (fabs(vBary.y()) > 0 && pBary.y() != 0) {
		double h = -pBary.y() / vBary.y();
		double r = pBary.x() + h * vBary.x();
		double t = pBary.z() + h * vBary.z();
		double w = 1 - r - t;

		//std::cout << "In 2: " << r << " " << t << " " << w << " " << h << std::endl;
		if (r > 0 && t > 0 && w > 0 && h > 0) {
			if (h > ch) {
				ch = h;
				th = h;
				i = 1;
				found = true;
			}
		}
	}

	if (fabs(vBary.z()) > 0 && pBary.z() != 0) {
		double h = -pBary.z() / vBary.z();
		double r = pBary.x() + h * vBary.x();
		double s = pBary.y() + h * vBary.y();
		double w = 1 - r - s;

		//std::cout << "In 3: " << r << " " << s << " " << w << " " << h << std::endl;
		if (r > 0 && s > 0 && w > 0 && h > 0) {
			if (h > ch) {
				ch = h;
				th = h;
				i = 2;
				found = true;
			}
		}
	}

	if (fabs(-vBary.x() - vBary.y() - vBary.z()) > 0 && 1 - pBary.x() - pBary.y() - pBary.z() != 0) {
		double h = -(1 - pBary.x() - pBary.y() - pBary.z()) / (-vBary.x() - vBary.y() - vBary.z());
		double r = pBary.x() + h * vBary.x();
		double s = pBary.y() + h * vBary.y();
		double t = pBary.z() + h * vBary.z();

		//std::cout << "In 4: " << r << " " << s << " " << t << " " << h << std::endl;
		if (r > 0 && s > 0 && t > 0 && h > 0) {
			if (h > ch) {
				ch = h;
				th = h;
				i = 3;
				found = true;
			}
		}
	}

	if (!found) {
		std::cout << "Error: Ray triangle intersection not found in tet.." << std::endl;
		cout << "Cell Info: " << ci->subdomain_index() << endl;
		int a;
		std::cin >> a;
	}

	return found;
}

bool RKIntegrator::_computeLocation(PolyVertex rkv, Mesh_3_Vector_3 v, PolyVertex& nrkv)
{
	bool terminate = false;
	PolyVertex pv = rkv;
	EigenVector3d vBary = _vectorGlobalToBarycentric(rkv.ci, v);

	double th;
	int i;
	//std::cout << "Searching ..." << std::endl;
	bool found = false;
	bool valid = false;
	do {
		found = _computeIntersection(rkv.ci, rkv.bary, vBary, th, i);
		
		if (found) {
			if (fabs(th) > 1) {
				if (rkv.ci->subdomain_index() == 0) {
					terminate = true;
					valid = false;
				}
				else {
					nrkv.ci = rkv.ci;
					nrkv.bary = rkv.bary + vBary;
					terminate = true;
					valid = true;
				}
			}
			else {
				EigenVector3d inBary = rkv.bary + th * vBary;
				EigenVector3d outBary = _baseTransformBarycentricCoordinate(rkv.ci, inBary, i);

				vBary -= th * vBary;
				v = _vectorBarycentricToGlobal(rkv.ci, vBary); // update v

				rkv.ci = rkv.ci->neighbor(i);
				rkv.bary = outBary;
				vBary = _vectorGlobalToBarycentric(rkv.ci, v); // update vBary
				
				if (rkv.ci->subdomain_index() == 0) {
					terminate = true;
					valid = false;
				}
				else
					terminate = false;
			}
		}
		else {
			terminate = true;
			valid = false;
		}
	} while (!terminate);

	//return nrkv;

	return valid;
}

double RKIntegrator::_length(Poly pl)
{
	double l = 0;
	auto end = --pl.vertices.end();
	for (auto iter = pl.vertices.begin(); iter != end; ++iter) {
		Mesh_3_Point_3 p1 = _vertexBarycentricToGlobal(iter->ci, iter->bary);
		auto it = iter;
		++it;
		Mesh_3_Point_3 p2 = _vertexBarycentricToGlobal(it->ci, it->bary);

		Mesh_3_Vector_3 v = (p1 - p2);
		double dl = sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z());

		l += dl;
	}

	return l;
}

double RKIntegrator::_distanceToFacet(Mesh_3_Point_3 p, Mesh_3_Facet_iterator fi)
{
	std::vector<Mesh_3_Point_3> vp;
	for (int i = 0; i < 4; ++i) {
		if (i != fi->second) {
			vp.push_back(fi->first->vertex(i)->point().point());
		}
	}

	Mesh_3_Vector_3 v1 = vp[1] - vp[0];
	Mesh_3_Vector_3 v2 = vp[2] - vp[0];

	Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
	n /= sqrt(n*n);

	Mesh_3_Vector_3 v = vp[0] - p;
	double d = v * n;

	return fabs(d);
}

Mesh_3_Cell_iterator RKIntegrator::_selectHighestDensityCell()
{
	Mesh_3_Cell_iterator ci = _mesh_3.IdxToCell(0);

	for (int i = 1; i < _mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator tci = _mesh_3.IdxToCell(i);

		if (_lines_in_cell[tci].size() > _lines_in_cell[ci].size())
			ci = tci;
	}

	return ci;
}
