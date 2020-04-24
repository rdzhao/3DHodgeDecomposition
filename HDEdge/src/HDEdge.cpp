#include "HDEdge.h"

My_Triangulation& HDEdge::Mesh()
{
	return m_mesh_3;
}

std::vector<Mesh_3_Point_3>& HDEdge::Vertices()
{
	return vertices;
}

std::vector<int>& HDEdge::Indices()
{
	return indices;
}

std::vector<Mesh_3_Vector_3>& HDEdge::Normals()
{
	return normals;
}

std::vector<Mesh_3_Vector_3>& HDEdge::Colors()
{
	return colors;
}

std::vector<Mesh_3_Vector_3>& HDEdge::VfVertices()
{
	return vfVertices;
}

std::vector<Mesh_3_Vector_3>& HDEdge::VfFaces()
{
	return vfFaces;
}

std::vector<Mesh_3_Vector_3>& HDEdge::VfNormals()
{
	return vfNormals;
}

std::vector<Mesh_3_Vector_3>& HDEdge::VfColors()
{
	return vfColors;
}

void HDEdge::initParameters(std::string sfx, double dr, double sr, int stps)
{
	m_suffix = sfx;
	m_density_ratio = dr;
	m_step_ratio = sr;
	m_steps = stps;
}

void HDEdge::buildMeshFromSurface(std::string fn, double size)
{
	//Polyhedron polyhedron;
	//std::ifstream input(fn.c_str());
	//input >> polyhedron;
	//input.close();
	//Polyhedron_mesh_domain domain(polyhedron);
	
	Polyhedron_with_feature polyhedron;
	std::ifstream input(fn.c_str());
	input >> polyhedron;
	input.close();
	Mesh_domain_with_feature domain(polyhedron);
	domain.detect_features();

	Mesh_criteria criteria(facet_angle = 30,
		facet_size = size,
		facet_distance = 0.2*size,
		cell_radius_edge_ratio = 2,
		cell_size = size);
	std::cout << "Criteria ... " << std::endl;

	//m_mesh_3.C3T3() = CGAL::makem_mesh_3<C3t3>(domain, criteria, CGAL::parameters::odt());
	//CGAL::perturbm_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, features(domain), odt());
	CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	std::cout << "Triangulation ... " << std::endl;

	m_mesh_3.preprocessing();
}

void HDEdge::buildLaplacian()
{
	buildLaplacian0();
	buildLaplacian1();
	buildLaplacian2();
}

void HDEdge::decompose()
{
	setFormRandom();

	cout << "Harmonic Knot ..." << endl;
	computeHarmonicKnot();
	cout << "Fluxless Knot ..." << endl;
	computeFluxlessKnot();
	cout << "Grounded Gradient ..." << endl;
	computeGroundedGradient();
	cout<<"Harmonic Gradient ..."<<endl;
	computeHarmonicGradient();
	cout<<"Curly Gradient ..."<<endl;
	computeCurlyGradient();
}

void HDEdge::visualize()
{
	convertForms();
	computeArrows();
	computeCrossSection();
}

void HDEdge::write()
{
	writeVFasVTK(omegaVF, "o");
	writeVFasVTK(harmonicKnotVF, "hk");
	writeVFasVTK(harmonicGradientVF, "hg");
	writeVFasVTK(fluxlessKnotVF, "fk");
	writeVFasVTK(groundedGradientVF, "gg");
	writeVFasVTK(curlyGradientVF, "cg");
}	

void HDEdge::integrate()
{
	cout<<"Omega Field ..."<<endl;
	integrateField(omegaVF, omegaCF, m_density_ratio, m_step_ratio, m_steps, m_suffix+"_o");

	if(m_mesh_3.NumHMLG1() != 0){
		cout<<"Harmonic Knot ..."<<endl;
		integrateField(harmonicKnotVF, harmonicKnotCF, m_density_ratio, m_step_ratio, m_steps, m_suffix+"_hk");
	}

	if(m_mesh_3.NumHMLG2() != 0){
		cout<<"Harmonic Gradient ..."<<endl;
		integrateField(harmonicGradientVF, harmonicGradientCF, m_density_ratio, m_step_ratio, m_steps, m_suffix+"_hg");
	}

	cout<<"Fluxless Knot ..."<<endl;
	integrateField(fluxlessKnotVF, fluxlessKnotCF, m_density_ratio, m_step_ratio, m_steps, m_suffix+"_fk");
	
	cout<<"Grounded Gradient ..."<<endl;
	integrateField(groundedGradientVF, groundedGradientCF, m_density_ratio, m_step_ratio, m_steps, m_suffix+"_gg");

	cout<<"Curly Gradient ..."<<endl;
	integrateField(curlyGradientVF, curlyGradientCF, m_density_ratio, m_step_ratio, m_steps, m_suffix+"_cg");
}

void HDEdge::buildLaplacian0()
{
	m_laplacian0_NB = m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB()*m_mesh_3.ED0F_NB();
	m_laplacian0_B = m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B()*m_mesh_3.ED0F_B();
	
}

void HDEdge::buildLaplacian1()
{	
	m_laplacian1_NB = m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB()
		+ m_mesh_3.ED0F_NB()*m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB();

	m_laplacian1_B = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()
		+ m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();


	m_curlEnergyMatrix = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B();
	m_divEnergyMatrix = m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();
	
}

void HDEdge::buildLaplacian2()
{
	m_laplacian2_B = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_B()
		+ m_mesh_3.ED1F_B()*m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B();
}

void HDEdge::computeBoundaryComponent()
{
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (m_mesh_3.C3T3().triangulation().is_infinite(vi))
			continue;

		vertexBoundaryIdx[vi] = -1;
	}

	std::map<Mesh_3_Facet_iterator, bool, Mesh_3_Facet_Iterator_Comparator> facetVisited;
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (m_mesh_3.FacetOnBoundary(fi)) {
			facetVisited[fi] = false;
		}
	}

	std::cout << "Classifying boundary" << std::endl;

	int seedIdx = 0;
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi) || !m_mesh_3.FacetOnBoundary(fi))
			continue;

		// find a seed
		// assume there is only one 3d mesh domain			
		bool visited = false;
		for (int i = 0; i < 4; ++i) {
			if (fi->second != i) {
				if (vertexBoundaryIdx[fi->first->vertex(i)] != -1)
					visited = true;
			}
		}

		if (!visited) {
			std::stack<Mesh_3_Facet_iterator> fstack;
			fstack.push(fi);
			facetVisited[fi] = true;

			while (!fstack.empty()) {
				Mesh_3_Facet_iterator fii = fstack.top();
				fstack.pop();

				// mark boundary vertices
				for (int i = 0; i < 4; ++i) {
					if (i != fii->second && vertexBoundaryIdx[fii->first->vertex(i)] == -1) {
						vertexBoundaryIdx[fii->first->vertex(i)] = seedIdx;
					}
				}

				// push adjacent boundary facet
				for (int i = 0; i < 4; ++i) {
					if (i != fii->second) {
						int idx1, idx2;
						idx1 = i;
						if ((i + 1) % 4 == fii->second)
							idx2 = (i + 2) % 4;
						else
							idx2 = (i + 1) % 4;

						if (idx1 > idx2)
							std::swap(idx1, idx2);

						My_edge me(fii->first, idx1, idx2);
						Edge e = *m_mesh_3.MyCellEdgeMap(me);

						// circulate e for adj faces
						Mesh_3_Facet_circulator fc = m_mesh_3.C3T3().triangulation().incident_facets(e);
						Mesh_3_Facet_circulator end = fc;
						do {
							if (m_mesh_3.C3T3().subdomain_index(fc->first)
								!= m_mesh_3.C3T3().subdomain_index(fc->first->neighbor(fc->second))
								&& *fc != *fii) {
								My_facet mf(fc->first, fc->second);

								if (!facetVisited[m_mesh_3.MyCellFacetMap(mf)]) {
									fstack.push(m_mesh_3.MyCellFacetMap(mf));
									facetVisited[m_mesh_3.MyCellFacetMap(mf)] = true;
								}
							}
							++fc;
						} while (fc != end);
					}
				}
			}

			++seedIdx;
		}
	}

	numOfBoundaryComponents = seedIdx;
	std::cout << "End " << seedIdx << std::endl;
}

void HDEdge::setFormRandom()
{
	omega.resize(m_mesh_3.NumE(), m_mesh_3.NumE());
	omega.setZero();

	EigenVector form;
	
	cout <<"Electric Field ..."<<endl;
	Mesh_3_Point_3 pcp(0.5, 0.7, 0.5);
	setPointChargeElectricField(form, pcp, true);
	omega += form;
	
	pcp = Mesh_3_Point_3(0.5, 0.5, 0.5);
	setPointChargeElectricField(form, pcp, false);
	omega += form;
	
	cout<<"Magnetic Field ..."<<endl;
	Mesh_3_Point_3 cp(0.2, 0.25, 0.7);
	Mesh_3_Vector_3 cv(0, 1, 0);
	setCurrentMagneticField(form, cp, cv, true);
	omega += form;

	cp = Mesh_3_Point_3(0.8, 0.5, 0.5);
    cv = Mesh_3_Vector_3(0, 0.5, 1);
    setCurrentMagneticField(form, cp, cv, false);
    omega += form;

	cp = Mesh_3_Point_3(0.6, 0.7, 0.2);
	cv = Mesh_3_Vector_3(1, 1, 0);
	setCurrentMagneticField(form, cp, cv, false);
	omega += form;
	
	cout<<"Eigen Field ..."<<endl;
	EigenSpMat s1l1 = m_mesh_3.HS1F_B()*m_laplacian1_B;

	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s1l1.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> hs;
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		hs.push_back(m_mesh_3.HS1F_B().coeff(i, i));
	}

	int num = 10;
	std::vector<double> emptyEvals;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEvals, eigenFields, row, col, val, hs, num);

	for(int i=0; i<eigenFields.size(); ++i){
		normalize(eigenFields[i]);
	}	
	
	groupFields(eigenFields);

	omega += eigenFields[curlFieldIdx[0]];
	omega += eigenFields[curlFieldIdx[2]];
}

void HDEdge::convertForms()
{
	convert1Form(omega, omegaVF);
	convert1Form(harmonicKnot, harmonicKnotVF);
	convert1Form(fluxlessKnot, fluxlessKnotVF);
	convert1Form(groundedGradient, groundedGradientVF);
	convert1Form(harmonicGradient, harmonicGradientVF);
	convert1Form(curlyGradient, curlyGradientVF);
}

void HDEdge::computeArrows()
{
	//assignColorGradient();
	assignColorStrength(omegaVF, 0);
	assignColorStrength(harmonicKnotVF, 1);
	assignColorStrength(harmonicGradientVF, 2);
	assignColorStrength(fluxlessKnotVF, 3);
	assignColorStrength(groundedGradientVF, 4);
	assignColorStrength(curlyGradientVF, 5);
	

	rescaleArrowsWithVF(omegaVF);
	rescaleArrowsWithVF(harmonicKnotVF);
	rescaleArrowsWithVF(harmonicGradientVF);
	rescaleArrowsWithVF(fluxlessKnotVF);
	rescaleArrowsWithVF(groundedGradientVF);
	rescaleArrowsWithVF(curlyGradientVF);	
}

void HDEdge::computeCrossSection()
{
	double cutThd = 10;

	// test for cross section of xy plane at center. If z>center.z, eliminate;
	int ii = 0;
	int faceCount = 0;
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {

		if (fi->first->subdomain_index() != fi->first->neighbor(fi->second)->subdomain_index()) { // boundary face

			Mesh_3_Cell_iterator validCell;

			bool side = true;
			for (int i = 0; i < 4; ++i) {
				if (i != fi->second) {
					if (fi->first->vertex(i)->point().point().z() < cutThd)
						side = true;
					else {
						side = false;
						break;
					}
				}
			}

			if (fi->first->subdomain_index() == 0) {
				if (m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point().z() > cutThd)
					side = false;

				validCell = fi->first->neighbor(fi->second);
			}
			else
			{
				if (fi->first->vertex(fi->second)->point().point().z() > cutThd)
					side = false;

				validCell = fi->first;
			}

			if (side) {
				std::vector<int> tmpIdx;

				for (int i = 0; i < 4; ++i) {
					if (i != fi->second) {
						if (vertexMap.find(fi->first->vertex(i)) == vertexMap.end()) {
							vertices.push_back(fi->first->vertex(i)->point().point());
							vertexMap[fi->first->vertex(i)] = ii;
							ii++;

							//colors.push_back(Mesh_3_Vector_3(0.6, 0.6, 0.6));
							colors.push_back(colorField[fi->first->vertex(i)]);

							if (!m_mesh_3.VertexOnBoundary(fi->first->vertex(i)))
								std::cout << "Error" << std::endl;
						}
						tmpIdx.push_back(vertexMap[fi->first->vertex(i)]);
					}
				}

				if ((fi->second % 2 == 1) == (fi->first->subdomain_index() == 0)) {
					indices.push_back(tmpIdx[0]);
					indices.push_back(tmpIdx[1]);
					indices.push_back(tmpIdx[2]);
				}
				else
				{
					indices.push_back(tmpIdx[1]);
					indices.push_back(tmpIdx[0]);
					indices.push_back(tmpIdx[2]);
				}
				realBoundary.push_back(true);

				// normal
				Mesh_3_Vector_3 v1, v2;
				v1 = vertices[indices[indices.size() - 2]] - vertices[indices[indices.size() - 1]];
				v2 = vertices[indices[indices.size() - 3]] - vertices[indices[indices.size() - 1]];

				Mesh_3_Vector_3 n = CGAL::cross_product(v2, v1);
				n /= sqrt(n*n);

				normals.push_back(n);
			}
		}
		else if (fi->first->subdomain_index() != 0) // non boundary face
		{
			int numTargetSide = 0;
			bool infinite = false;

			for (int i = 0; i < 4; ++i) {
				if (i != fi->second) {
					if (fi->first->subdomain_index() == 0)
						infinite = true;

					if (fi->first->vertex(i)->point().point().z() < cutThd)
						numTargetSide++;
				}
			}

			if (infinite)
				continue;

			bool onBoundary = false;
			Mesh_3_Cell_iterator boundaryCell;
			Mesh_3_Cell_iterator validCell; // for extracting vector field.
			if (numTargetSide == 3) {
				if (fi->first->vertex(fi->second)->point().point().z() > cutThd) {
					boundaryCell = fi->first;
					validCell = fi->first->neighbor(fi->second);
					onBoundary = !onBoundary;
				}

				if (m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point().z() > cutThd) {
					boundaryCell = fi->first->neighbor(fi->second);
					validCell = fi->first;
					onBoundary = !onBoundary;
				}
			}

			if (onBoundary) {
				std::vector<int> tmpIdx;
				for (int i = 0; i < 4; ++i) {
					if (i != fi->second) {
						if (vertexMap.find(fi->first->vertex(i)) == vertexMap.end()) {
							vertices.push_back(fi->first->vertex(i)->point().point());
							vertexMap[fi->first->vertex(i)] = ii;
							ii++;

							//colors.push_back(Mesh_3_Vector_3(0.6, 0.6, 0.6));
							colors.push_back(colorField[fi->first->vertex(i)]);
						}
						tmpIdx.push_back(vertexMap[fi->first->vertex(i)]);
					}
				}

				if ((fi->second % 2 == 1) == (fi->first != boundaryCell)) {
					indices.push_back(tmpIdx[1]);
					indices.push_back(tmpIdx[0]);
					indices.push_back(tmpIdx[2]);
				}
				else
				{
					indices.push_back(tmpIdx[0]);
					indices.push_back(tmpIdx[1]);
					indices.push_back(tmpIdx[2]);
				}
				realBoundary.push_back(false);

				// normal
				Mesh_3_Vector_3 v1, v2;
				v1 = vertices[indices[indices.size() - 2]] - vertices[indices[indices.size() - 1]];
				v2 = vertices[indices[indices.size() - 3]] - vertices[indices[indices.size() - 1]];

				Mesh_3_Vector_3 n = CGAL::cross_product(v2, v1);
				n /= sqrt(n*n);

				normals.push_back(n);
			}
		}
	}
}

void HDEdge::computeHarmonicKnot()
{
	EigenSpMat s1l1 = m_mesh_3.HS1F_B() * m_laplacian1_B;

	double hm1 = m_mesh_3.NumHMLG1();
	cout << "HM1: " << hm1 << endl;

	if (hm1 == 0) {
		harmonicKnot.setZero(m_mesh_3.NumE());
	}
	else {
		// assemble matrix to be eigen decomposed
		std::vector<double> row, col;
		std::vector<double> val;
		for (int i = 0; i < s1l1.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
				row.push_back(static_cast<double>(iter.row()) + 1);
				col.push_back(static_cast<double>(iter.col()) + 1);
				val.push_back(static_cast<double>(iter.value()));
			}
		}
		std::vector<double> hs;
		for (int i = 0; i < m_mesh_3.NumE(); ++i) {
			hs.push_back(m_mesh_3.HS1F_B().coeff(i, i));
		}

		std::vector<double> emptyEval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, harmonicKnotBasis, row, col, val, hs, hm1);

		// we have harmonic basis
		harmonicKnot.setZero(m_mesh_3.NumE());
		for (int i = 0; i < harmonicKnotBasis.size(); ++i) {
			double c = (omega.transpose()*m_mesh_3.HS1F_B()*harmonicKnotBasis[i]).coeff(0,0)
				/ (harmonicKnotBasis[i].transpose() * m_mesh_3.HS1F_B()*harmonicKnotBasis[i]).coeff(0, 0);
			harmonicKnotCoeffs.push_back(c);
			harmonicKnot += c * harmonicKnotBasis[i];
		}
	}
}

void HDEdge::computeHarmonicGradient()
{
	EigenSpMat s1l1nb = m_mesh_3.HS1F_NB() * m_laplacian1_NB;

	double hm2 = m_mesh_3.NumHMLG2();
	cout << "HM2: " << hm2 << endl;

	if (hm2 == 0) {
		harmonicGradient.setZero(m_mesh_3.NumE());
	}
	else {
		// assemble matrix to be eigen decomposed
		std::vector<double> row, col;
		std::vector<double> val;
		for (int i = 0; i < s1l1nb.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(s1l1nb, i); iter; ++iter) {
				row.push_back(static_cast<double>(iter.row()) + 1);
				col.push_back(static_cast<double>(iter.col()) + 1);
				val.push_back(static_cast<double>(iter.value()));
			}
		}
		std::vector<double> hs;
		for (int i = 0; i < m_mesh_3.NumIE(); ++i) {
			hs.push_back(m_mesh_3.HS1F_NB().coeff(i, i));
		}

		std::vector<double> emptyEval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, harmonicGradientBasis, row, col, val, hs, hm2);

		// we have harmonic basis
		harmonicGradient.setZero(m_mesh_3.NumIE());
		EigenVector omegaNB = omega;
		boundaryToNoBoundary1Form(omegaNB);
		for (int i = 0; i < harmonicGradientBasis.size(); ++i) {
			double c = (omegaNB.transpose()*m_mesh_3.HS1F_NB()*harmonicGradientBasis[i]).coeff(0,0)
				/ (harmonicGradientBasis[i].transpose() * m_mesh_3.HS1F_NB()*harmonicGradientBasis[i]).coeff(0,0);
			harmonicGradientCoeffs.push_back(c);
			harmonicGradient += c * harmonicGradientBasis[i];
		}
		noBoundaryToBoundary1Form(harmonicGradient);
	}
}

void HDEdge::computeFluxlessKnot()
{
	EigenSpMat s2l2 = m_mesh_3.HS2F_B() * m_laplacian2_B;

	double hm2 = m_mesh_3.NumHMLG2();

	if (hm2 == 0) {
		EigenVector b = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*omega;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(s2l2);

		EigenVector x = solver.solve(b); // potential 2-form for exact component

		fluxlessKnot = m_mesh_3.HS1F_B().cwiseInverse() *m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B() * x;
	}
	else {
		std::vector<double> row, col;
		std::vector<double> val;
		for (int i = 0; i < s2l2.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(s2l2, i); iter; ++iter) {
				row.push_back(static_cast<double>(iter.row()) + 1);
				col.push_back(static_cast<double>(iter.col()) + 1);
				val.push_back(static_cast<double>(iter.value()));
			}
		}
		std::vector<double> hs;
		for (int i = 0; i < m_mesh_3.NumF(); ++i) {
			hs.push_back(m_mesh_3.HS2F_B().coeff(i, i));
		}
		
		std::vector<double> emptyEval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, eigenFields, row, col, val, hs, hm2);

		EigenVector b = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*omega;

		std::set<int> anchors;
		selectAnchors(eigenFields, anchors);
		correctRankDeficiency(s2l2, b, anchors);

		// solve
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(s2l2);
		EigenVector x = solver.solve(b); // potential 2-form for exact component

		fluxlessKnot = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x;
	}
}

void HDEdge::computeGroundedGradient()
{
	EigenSpMat s0l0nb = m_mesh_3.HS0F_NB()*m_laplacian0_NB;
	
	EigenVector omegaNB = omega;
	boundaryToNoBoundary1Form(omegaNB);

	EigenVector b = m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB()*omegaNB;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(s0l0nb);

	EigenVector x = solver.solve(b);

	groundedGradient = m_mesh_3.ED0F_NB()*x;
	noBoundaryToBoundary1Form(groundedGradient);

	for (int j = 0; j < m_mesh_3.NumV(); ++j) {
		if (m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(j)))
			scalarPotential[m_mesh_3.IdxToVertex(j)] = 0;
		else
			scalarPotential[m_mesh_3.IdxToVertex(j)] = x(m_mesh_3.VertexForward(j));
	}
}

void HDEdge::computeCurlyGradient()
{
	curlyGradient = omega - harmonicKnot - harmonicGradient - fluxlessKnot - groundedGradient;
}

void HDEdge::computeCurlFreeGradient()
{
	EigenSpMat s0l0 = m_mesh_3.HS0F_B()*m_laplacian0_B;
	EigenVector b = m_mesh_3.ED0F_B().transpose()* m_mesh_3.HS1F_B()*omega;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(s0l0);
	EigenVector x = solver.solve(b);

	curlFreeGradient = m_mesh_3.ED0F_B()*x;
}

void HDEdge::integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf, 
	double dsty, double lth, int stps, std::string sfx)
{
	RKIntegrator integrator(m_mesh_3);

	
	integrator.init(dsty, lth, stps);
	integrator.initColorField(cf);
	integrator.convert(vf);
	integrator.integrate();
	integrator.write(sfx);
}

void HDEdge::matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
	std::vector<double> row, std::vector<double> col, std::vector<double> val,
	std::vector<double> hs, double hm)
{
	std::vector<double> vhm;
	vhm.push_back(hm);

	using namespace matlab::engine;

	std::unique_ptr<MATLABEngine> mptr = startMATLAB();
	matlab::data::ArrayFactory factory;

	matlab::data::TypedArray<double> mrow = factory.createArray<std::vector<double>::iterator>({ row.size(), 1 }, row.begin(), row.end());
	matlab::data::TypedArray<double> mcol = factory.createArray<std::vector<double>::iterator>({ col.size(), 1 }, col.begin(), col.end());
	matlab::data::TypedArray<double> mval = factory.createArray<std::vector<double>::iterator>({ val.size(), 1 }, val.begin(), val.end());
	matlab::data::TypedArray<double> mhs = factory.createArray<std::vector<double>::iterator>({ hs.size(), 1 }, hs.begin(), hs.end());
	matlab::data::TypedArray<double> mhm = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, vhm.begin(), vhm.end());

	mptr->setVariable(u"row", std::move(mrow));
	mptr->setVariable(u"col", std::move(mcol));
	mptr->setVariable(u"val", std::move(mval));
	mptr->setVariable(u"hs", std::move(mhs));
	mptr->setVariable(u"hm", std::move(mhm));

	mptr->eval(u"A=sparse(row, col, val);");
	mptr->eval(u"n=length(hs);");
	mptr->eval(u"B = spdiags(hs(:),0,n,n);");
	mptr->eval(u"[V,D]=eigs(A, B, hm, 'smallestabs');");

	matlab::data::TypedArray<double> dd = mptr->getVariable(u"D");
	matlab::data::TypedArray<double> vv = mptr->getVariable(u"V");

	EigenVector ev;
	ev.resize(vv.getDimensions()[0]);

	evec.clear();
	for (int j = 0; j < vv.getDimensions()[1]; ++j) {
		for (int i = 0; i < vv.getDimensions()[0]; ++i) {
			ev(i) = vv[i][j];
		}
		evec.push_back(ev);
	}
}

void HDEdge::groupFields(std::vector<EigenVector> eigenFields)
{
	int curlIdx = 0;
	int divIdx = 0;

	for(int i = 0; i < eigenFields.size(); ++i){
		if(eigenFields[i].dot(m_curlEnergyMatrix*eigenFields[i])
			> eigenFields[i].dot(m_divEnergyMatrix*eigenFields[i])){
			curlFieldIdx[divIdx] = i;
			++divIdx;
		}
		else{
			divFieldIdx[curlIdx] = i;
			++curlIdx;
		}
	}
}

void HDEdge::setPointChargeElectricField(EigenVector& field, Mesh_3_Point_3 rp, bool ori)
{
	// rp: relative position of bounding box
	Mesh_3_Point_3 bbMin, bbMax;
	bbMin = m_mesh_3.BoundingBoxMin();
	bbMax = m_mesh_3.BoundingBoxMax();

	Mesh_3_Point_3 p((1 - rp.x())*bbMin.x() + rp.x()*bbMax.x(), (1 - rp.y())*bbMin.y() + rp.y()*bbMax.y(), (1 - rp.z())*bbMin.z() + rp.z()*bbMax.z());

	field.setZero(m_mesh_3.NumE());
	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;
			
		Mesh_3_Point_3 ec = m_mesh_3.EdgeCC(ei);
		Mesh_3_Vector_3 d = ec - p;
		d /= sqrt(d * d);

		Mesh_3_Vertex_iterator v1, v2;
		v1 = ei->first->vertex(ei->second);
		v2 = ei->first->vertex(ei->third);

		Mesh_3_Vector_3 l;
		if (m_mesh_3.VertexIdx(v1) > m_mesh_3.VertexIdx(v2)) {
			l = v1->point().point() - v2->point().point();
		}
		else {
			l = v2->point().point() - v1->point().point();
		}


		if (ori)
			field(m_mesh_3.EdgeIdx(ei)) = d * l;
		else
			field(m_mesh_3.EdgeIdx(ei)) = -d * l;
	}

	normalize(field);
}

void HDEdge::setCurrentMagneticField(EigenVector& field, Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, bool ori)
{
	// rp: relative position of bounding box
	Mesh_3_Point_3 bbMin, bbMax;
	bbMin = m_mesh_3.BoundingBoxMin();
	bbMax = m_mesh_3.BoundingBoxMax();

	Mesh_3_Point_3 p((1 - rp.x())*bbMin.x() + rp.x()*bbMax.x(), (1 - rp.y())*bbMin.y() + rp.y()*bbMax.y(), (1 - rp.z())*bbMin.z() + rp.z()*bbMax.z());
	
	field.setZero(m_mesh_3.NumE());
	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;

		Mesh_3_Point_3 ec = m_mesh_3.EdgeCC(ei); // facet circum center
		Mesh_3_Vector_3 vpfc = ec - p; // vector p to fc
		Mesh_3_Vector_3 vpfcp = (vpfc * d)*d; // vector p to fc protected on d 
		Mesh_3_Vector_3 vpd = vpfc - vpfcp; // vector p to direction d
		
		Mesh_3_Vector_3 fD = CGAL::cross_product(d, vpd); // field direction
		fD /= sqrt(fD*fD); // normalize

		Mesh_3_Vertex_iterator v1, v2;
		v1 = ei->first->vertex(ei->second);
		v2 = ei->first->vertex(ei->third);

		Mesh_3_Vector_3 l;
		if (m_mesh_3.VertexIdx(v1) > m_mesh_3.VertexIdx(v2)) {
			l = v1->point().point() - v2->point().point();
		}
		else {
			l = v2->point().point() - v1->point().point();
		}

		if (ori)
			field(m_mesh_3.EdgeIdx(ei)) = fD * l;
		else
			field(m_mesh_3.EdgeIdx(ei)) = -fD * l;
	}

	normalize(field);
}

void HDEdge::normalize(EigenVector& form)
{
	form = form / sqrt(form.dot(form));
}

void HDEdge::boundaryToNoBoundary1Form(EigenVector& form)
{
	EigenVector formNB;
	formNB.resize(m_mesh_3.NumIE());

	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;

		if (!m_mesh_3.EdgeOnBoundary(ei))
			formNB(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei))) = form(m_mesh_3.EdgeIdx(ei));
	}

	form = formNB;
}

void HDEdge::noBoundaryToBoundary1Form(EigenVector& form)
{
	EigenVector formB;
	formB.resize(m_mesh_3.NumE());

	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;

		if (m_mesh_3.EdgeOnBoundary(ei))
			formB(m_mesh_3.EdgeIdx(ei)) = 0;
		else
			formB(m_mesh_3.EdgeIdx(ei)) = form(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)));
	}

	form = formB;
}

void HDEdge::boundaryToNoBoundary2Form(EigenVector& form)
{
	EigenVector formNB;
	formNB.resize(m_mesh_3.NumIF());

	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;

		if (!m_mesh_3.FacetOnBoundary(fi))
			formNB(m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi))) = form(m_mesh_3.FacetIdx(fi));
	}

	form = formNB;
}

void HDEdge::noBoundaryToBoundary2Form(EigenVector& form)
{
	EigenVector formB;
	formB.resize(m_mesh_3.NumF());

	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;

		if (m_mesh_3.FacetOnBoundary(fi))
			formB(m_mesh_3.FacetIdx(fi)) = 0;
		else
			formB(m_mesh_3.FacetIdx(fi)) = form(m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)));
	}

	form = formB;
}

void HDEdge::computeScalarPotential(EigenVector form)
{
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		scalarPotential[vi] = form(m_mesh_3.VertexIdx(vi));
	}
}

void HDEdge::convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
{
	// compute gradient of phi_i in cell j, corresponding a facet.
	std::map<My_facet, Mesh_3_Vector_3> grad;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		for (int i = 0; i < 4; ++i) {
			My_facet mf(ci, i);

			// magnitude
			double h = 3 * m_mesh_3.CellPrimal(ci) / m_mesh_3.FacetPrimal(m_mesh_3.MyCellFacetMap(mf));

			// direction normal
			std::vector<Mesh_3_Point_3> pv;
			for (int j = 0; j < 4; ++j) {
				if (j != i) {
					pv.push_back(ci->vertex(j)->point().point());
				}
			}

			Mesh_3_Vector_3 v1 = pv[1] - pv[0];
			Mesh_3_Vector_3 v2 = pv[2] - pv[0];

			Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
			n /= sqrt(n*n);

			Mesh_3_Vector_3 vc = pv[0] - ci->vertex(i)->point().point();
			if (vc*n < 0)
				n = -n;

			grad[mf] = n / h;
		}

		// iterate through edges in a cell to reconstruct vector field
		Mesh_3_Vector_3 of(0, 0, 0); // original vf on cell
		for (int i = 0; i < 4; ++i) {
			for (int j = i + 1; j < 4; ++j) {
				My_edge me(ci, i, j);
				My_facet mf1(ci, i);
				My_facet mf2(ci, j);

				if (m_mesh_3.VertexIdx(ci->vertex(i)) > m_mesh_3.VertexIdx(ci->vertex(j))) {
					of += form[m_mesh_3.EdgeIdx(m_mesh_3.MyCellEdgeMap(me))] * (1.0 / 3) * (-grad[mf1] + grad[mf2]);
				}
				else {
					of += form[m_mesh_3.EdgeIdx(m_mesh_3.MyCellEdgeMap(me))] * (1.0 / 3) * (grad[mf1] - grad[mf2]);
				}
			}
		}
		vf[ci] = of;
	}
}

void HDEdge::convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
{
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		std::map<int, Mesh_3_Vector_3> grad;
		for (int i = 0; i < 4; ++i) {
			My_facet mf(ci, i);

			// magnitude
			double h = 3 * m_mesh_3.CellPrimal(ci) / m_mesh_3.FacetPrimal(m_mesh_3.MyCellFacetMap(mf));

			// direction normal
			std::vector<Mesh_3_Point_3> pv;
			for (int j = 0; j < 4; ++j) {
				if (j != i) {
					pv.push_back(ci->vertex(j)->point().point());
				}
			}

			Mesh_3_Vector_3 v1 = pv[1] - pv[0];
			Mesh_3_Vector_3 v2 = pv[2] - pv[0];

			Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
			n /= sqrt(n*n);

			Mesh_3_Vector_3 vc = pv[0] - ci->vertex(i)->point().point();
			if (vc*n < 0)
				n = -n;

			grad[m_mesh_3.VertexIdx(ci->vertex(i))] = n / h;
			//std::cout<< grad[vertexIdx[ci->vertex(i)]] <<std::endl;
		}

		// evaluate piece wise constant vector in the cell
		Mesh_3_Vector_3 rv(0, 0, 0); // recovered vector
		for (int i = 0; i < 4; ++i) { //iterate through 4 faces;
			std::vector<int> vIdx;
			for (int j = 0; j < 4; ++j) {
				if (j != i) {
					vIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(j)));
				}
			}

			std::sort(vIdx.begin(), vIdx.end());

			My_facet mf(ci, i);

			rv += form[m_mesh_3.FacetIdx(m_mesh_3.MyCellFacetMap(mf))] * 2 * 0.25
				* (CGAL::cross_product(grad[vIdx[0]], grad[vIdx[1]])
					+ CGAL::cross_product(grad[vIdx[1]], grad[vIdx[2]])
					+ CGAL::cross_product(grad[vIdx[2]], grad[vIdx[0]]));
		}

		vf[ci] = rv;
	}
}

void HDEdge::writeVFasVTK(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, std::string sfx)
{
	std::string fn="vector_field_"+sfx+".vtk";
    std::ofstream out(fn.c_str());

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"3D triangulation with vector field"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_mesh_3.NumV()<<" double"<<std::endl;
    for(Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
        vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi){
            if(!m_mesh_3.VertexValid(vi)){
                continue;
            }

            Mesh_3_Point_3 p = vi->point().point();
            out<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
        }

    out<<"CELLS "<<m_mesh_3.NumC()<<" "<<5*m_mesh_3.NumC()<<std::endl;
    for(Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
        ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci){
            if(ci->subdomain_index() == 0){
                continue;
            }

            std::vector<int> ids;
            for(int i=0; i<4; ++i){
                ids.push_back(m_mesh_3.VertexIdx(ci->vertex(i)));
            }

            out<<"4 "<<ids[0]<<" "<<ids[1]<<" "<<ids[2]<<" "<<ids[3]<<std::endl;
        }

    out<<"CELL_TYPES "<<m_mesh_3.NumC()<<std::endl;
    for(int i=0; i<m_mesh_3.NumC(); ++i){
        out<<"10"<<std::endl;
    }

    out<<"CELL_DATA "<<m_mesh_3.NumC()<<std::endl;
    out<<"VECTORS U double"<<std::endl;
    for(Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
        ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci){
            if(ci->subdomain_index() == 0){
                continue;
            }

            Mesh_3_Vector_3 v = vf[ci];
            if(sqrt(v*v)<EPS)
                v = Mesh_3_Vector_3(0,0,0);
            out<<v.x()<<" "<<v.y()<<" "<<v.z()<<std::endl;
        }

    out.close();
}

void HDEdge::rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
{
	// set the max length to 0.7*aveLen 
	// get vector ave length 
	double aveVecLen = 0;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Vector_3 v = vf[ci];
		double len = sqrt(v*v);
		aveVecLen += len;
	}
	aveVecLen /= m_mesh_3.NumC();

	double pivotLen;
	std::vector<double> lenVec;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Vector_3 v = vf[ci];
		double len = sqrt(v*v);
		lenVec.push_back(len);
	}
	std::sort(lenVec.begin(), lenVec.end());
	pivotLen = lenVec[static_cast<int>(floor(lenVec.size()*0.6))];
	double aveLen = m_mesh_3.getAveLen();
	//if (pivotLen < aveLen)
	//	pivotLen = aveLen;

	
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Vector_3 v = vf[ci];
		double len = sqrt(v*v);
		v /= len;

		if (len > pivotLen)
			len = pivotLen;
		else if(len <0.1*pivotLen)
			len = 0.1*pivotLen;
			
		len = (aveLen)*(len / (pivotLen));

		vf[ci] = v * len;
	}
}

void HDEdge::selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors)
{
	cout << eigenfields.size() << endl;

	std::random_device rd;
	std::mt19937_64 mt(rd());
	std::uniform_int_distribution<int> distribution(0, m_mesh_3.NumIF());

	int condNumber = 0;

	int size = static_cast<int>(eigenfields.size());
	int iteration = 10;
	for (int i = 0; i < iteration; ++i) {
		//cout<<"---------------------"<<endl;
		std::set<int> rows;

		do {
			int d = distribution(mt);
			rows.insert(d);
			//cout <<"###: "<< rows.size() << endl;
		} while (rows.size() < size);
		//cout << "---------------------" << endl;
		EigenMatrix mat(size, size);
		int j = 0;
		for (auto iter = rows.begin(); iter != rows.end(); ++iter, ++j) {
			for (int k = 0; k < size; ++k) {
				//cout<<"@@@"<<endl;
				mat(j, k) = eigenfields[k](*iter);
				//cout << "@@@" << endl;
			}
		}
		//cout << "---------------------" << endl;
		//cout << mat(0, 0) << endl;

		if (mat.determinant() == 0) {
			--i;
			continue;
		}
		else {
			Eigen::JacobiSVD<EigenMatrix> svd(mat);
			double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

			if (cond > condNumber) {
				anchors = rows;
			}
		}
	}

	for (auto iter = anchors.begin(); iter != anchors.end(); ++iter) {
		cout << *iter << endl;
	}
}

void HDEdge::correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors)
{
	for (int i = 0; i < m.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m, i); iter; ++iter) {
			if (anchors.find(iter.row()) != anchors.end()
				|| anchors.find(iter.col()) != anchors.end()) {
				if (iter.row() == iter.col()) {
					iter.valueRef() = 1;
				}
				else {
					iter.valueRef() = 0;
				}
			}
		}
	}

	for (auto iter = anchors.begin(); iter != anchors.end(); ++iter) {
		b(*iter) = 0;
	}
}

void HDEdge::assignColorStrength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode)
{
	// calculate vertex based strength
	std::map<Mesh_3_Vertex_iterator, double> vertexStrength;
	std::vector<double> len;
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		std::vector<Mesh_3_Cell_iterator> cells;
		m_mesh_3.C3T3().triangulation().incident_cells(vi, std::back_inserter(cells));

		double p = 0;
		double vol = 0;
		for (int i = 0; i < cells.size(); ++i) {
			if (cells[i]->subdomain_index() == 0)
				continue;

			//p += potential[mm_mesh_3.CellIdx(cells[i])];
			Mesh_3_Vector_3 v = vf[cells[i]];
			p += sqrt(v*v) * m_mesh_3.CellPrimal(cells[i]);
			vol += m_mesh_3.CellPrimal(cells[i]);
		}
		p /= vol;

		vertexStrength[vi] = p;
		len.push_back(p);
	}

	std::vector<int> sIndices;
	sortedIndices(sIndices, len);

	for (int i = 0; i < sIndices.size(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(sIndices[i]);
		double ratio = i * 1.0 / sIndices.size();

		Mesh_3_Vector_3 col = assignColor(ratio);

		if (mode == 0) {
			omegaCF[vi] = col;
		}
		else if (mode == 1) {
			harmonicKnotCF[vi] = col;
		}
		else if (mode == 2) {
			harmonicGradientCF[vi] = col;
		}
		else if (mode == 3) {
			fluxlessKnotCF[vi] = col;
		}
		else if (mode == 4) {
			groundedGradientCF[vi] = col;
		}
		else if (mode == 5) {
			curlyGradientCF[vi] = col;
		}
	}

	if (mode == 0) {
		colorField = omegaCF;
	}
	else if (mode == 1) {
		colorField = harmonicKnotCF;
	}
	else if (mode == 2) {
		colorField = harmonicGradientCF;
	}
	else if (mode == 3) {
		colorField = fluxlessKnotCF;
	}
	else if (mode == 4) {
		colorField = groundedGradientCF;
	}
	else if (mode == 5) {
		colorField = curlyGradientCF;
	}


}

void HDEdge::assignColorGradient()
{
	// assign color 
	double pmin = scalarPotential[m_mesh_3.IdxToVertex(0)];
	double pmax = scalarPotential[m_mesh_3.IdxToVertex(0)];
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (scalarPotential[vi] > pmax)
			pmax = scalarPotential[vi];
		if (scalarPotential[vi] < pmin)
			pmin = scalarPotential[vi];
	}

	double diff = pmax - pmin;
	Mesh_3_Vector_3 cmax = Mesh_3_Vector_3(0.8, 0.4, 0.8);
	Mesh_3_Vector_3 cmin = Mesh_3_Vector_3(0.8, 0.8, 0.8);
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		double val = scalarPotential[vi];
		double d = val - (pmin + 0.0*diff);

		double ratio = d / (diff*1.0);
		if (ratio < 0)
			ratio = 0;
		if (ratio > 1)
			ratio = 1;

		groundedGradientCF[vi] = ratio * cmax + (1 - ratio)*cmin;
	}

	colorField = groundedGradientCF;
}

Mesh_3_Vector_3 HDEdge::assignColor(double ratio)
{
	// 0.0 - 0.25: blue to indigo
	// 0.25 - 0.5: indigo to green
	// 0.5 - 0.75: green to yellow
	// 0.75 - 1.0: yellow to red

	Mesh_3_Vector_3 color;

	Mesh_3_Vector_3 red(0.9, 0.7, 0.7);
	Mesh_3_Vector_3 yellow(0.9, 0.9, 0.7);
	Mesh_3_Vector_3 green(0.7, 0.9, 0.7);
	Mesh_3_Vector_3 indigo(0.7, 0.9, 0.9);
	Mesh_3_Vector_3 blue(0.7, 0.7, 0.9);



	if (ratio < 0.25) {
		double r = ratio / 0.25;
		color = (1 - r)*blue + r * indigo;
	}
	else if (ratio < 0.5) {
		double r = (ratio - 0.25) / 0.25;
		color = (1 - r)*indigo + r * green;
	}
	else if (ratio < 0.75) {
		double r = (ratio - 0.5) / 0.25;
		color = (1 - r)*green + r * yellow;
	}
	else {
		double r = (ratio - 0.75) / 0.25;
		color = (1 - r)*yellow + r * red;
	}

	return color;
}

double HDEdge::medianLength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	std::vector<double> vl;
	for (auto iter = vf.begin(); iter != vf.end(); ++iter) {
		vl.push_back(sqrt(iter->second*iter->second));
	}
	std::nth_element(vl.begin(), vl.begin() + vl.size() / 2, vl.end());

	return vl[vl.size() / 2];
}

void HDEdge::sortedIndices(std::vector<int>& sIndex, std::vector<double> len)
{
	sIndex.resize(len.size(), 0);
	std::iota(sIndex.begin(), sIndex.end(), 0);

	std::sort(sIndex.begin(), sIndex.end(), [&len](int i1, int i2) {return len[i1] < len[i2]; });
}
