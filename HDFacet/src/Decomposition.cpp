#include "Decomposition.h"


My_Triangulation& Decomposition::Mesh()
{
	return m_mesh_3;
}

std::vector<Mesh_3_Point_3>& Decomposition::Vertices()
{
	return vertices;
}

std::vector<int>& Decomposition::Indices()
{
	return indices;
}

std::vector<Mesh_3_Vector_3>& Decomposition::Normals()
{
	return normals;
}

std::vector<Mesh_3_Vector_3>& Decomposition::Colors()
{
	return colors;
}

std::vector<Mesh_3_Vector_3>& Decomposition::VfVertices()
{
	return vfVertices;
}

std::vector<Mesh_3_Vector_3>& Decomposition::VfFaces()
{
	return vfFaces;
}

std::vector<Mesh_3_Vector_3>& Decomposition::VfNormals()
{
	return vfNormals;
}

std::vector<Mesh_3_Vector_3>& Decomposition::VfColors()
{
	return vfColors;
}

void Decomposition::initParameters(std::string sfx, double dr, double sr, int stps)
{
	m_suffix = sfx;
	m_density_ratio = dr;
	m_step_ratio = sr;
	m_steps = stps;
}

void Decomposition::buildMeshFromSurface(std::string fn, double size)
{
	/*
	Polyhedron polyhedron;
	std::ifstream input(fn.c_str());
	input >> polyhedron;
	input.close();
	Polyhedron_mesh_domain domain(polyhedron);

	Mesh_criteria criteria(facet_angle = 30,
		facet_size = size,
		facet_distance = 0.1*size,
		cell_radius_edge_ratio = 2,
		cell_size = size);
	std::cout << "Criteria ... " << std::endl;

	m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::odt());
	CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	std::cout << "Triangulation ... " << std::endl;
	*/
	
	Polyhedron_with_feature polyhedron;
	std::ifstream input(fn.c_str());
	input >> polyhedron;
	input.close();
	Mesh_domain_with_feature domain(polyhedron);
	//domain.detect_features();

	Mesh_criteria criteria(edge_size = size,
		facet_angle = 30,
		facet_size = size,
		facet_distance = 0.1*size,
		cell_radius_edge_ratio = 2,
		cell_size = size);
	std::cout << "Criteria ... " << std::endl;

	m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, features(domain), odt());
	CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	std::cout << "Triangulation ... " << std::endl;
	
	m_mesh_3.preprocessing();
	cout << "Euler: " << m_mesh_3.NumV() - m_mesh_3.NumE() + m_mesh_3.NumF() - m_mesh_3.NumC() << endl;
}

void Decomposition::buildLaplacian()
{
	buildLaplacian1();
	buildLaplacian2(); 
	buildLaplacian3();
}

void Decomposition::decompose()
{
	setFormRandom();

	cout << "Harmonic Knot ..." << endl;
	computeHarmonicKnot();
	cout << "Harmonic Gradient ..." << endl;
	computeHarmonicGradient();
	cout << "Fluxless Knot ..." << endl;
	computeFluxlessKnot();
	cout << "Grounded Gradient ..."<< endl;
	computeGroundedGradient();
	cout << "Curly Gradient ..." << endl;
	computeCurlyGradient();
}

void Decomposition::visualize()
{
	convertForms();
	computeArrows();
	cout<<"Vector field visualization done ..."<<endl;
}

void Decomposition::integrate()
{
	cout<<"Omega Field ..."<<endl;
	integrateField(omegaVF, omegaCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix+"_o");

	if(m_mesh_3.NumHMLG1() != 0){
		cout<<"Harmonic Knot ..."<<endl;
		integrateField(harmonicKnotVF, harmonicKnotCF, m_density_ratio, m_step_ratio, m_steps, 0, m_suffix+"_hk");
	}

	if(m_mesh_3.NumHMLG2() != 0){
		cout<<"Harmonic Gradient ..."<<endl;
		integrateField(harmonicGradientVF, harmonicGradientCF, m_density_ratio/3, m_step_ratio, m_steps, 2, m_suffix+"_hg");
	}

	cout<<"Fluxless Knot ..."<<endl;
	integrateField(fluxlessKnotVF, fluxlessKnotCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix+"_fk");
	
	cout<<"Grounded Gradient ..."<<endl;
	integrateField(groundedGradientVF, groundedGradientCF, m_density_ratio/2, m_step_ratio, m_steps, 2, m_suffix+"_gg");

	cout<<"Curly Gradient ..."<<endl;
	integrateField(curlyGradientVF, curlyGradientCF, m_density_ratio, m_step_ratio, m_steps, 1 ,m_suffix+"_cg");
}

void Decomposition::write()
{
	writeVFasVTK(omegaVF, "o");

	if(m_mesh_3.NumHMLG1() != 0)
		writeVFasVTK(harmonicKnotVF, "hk");
	
	if(m_mesh_3.NumHMLG2() != 0)
		writeVFasVTK(harmonicGradientVF, "hg");
	
	writeVFasVTK(fluxlessKnotVF, "fk");
	writeVFasVTK(groundedGradientVF, "gg");
	writeVFasVTK(curlyGradientVF, "cg");
}

void Decomposition::buildLaplacian1()
{

	m_laplacian1_NB = m_mesh_3.HS1F_NB().cwiseInverse() *m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB()
		+ m_mesh_3.ED0F_NB() *  m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB(); 
	m_laplacian1_B = m_mesh_3.HS1F_B().cwiseInverse() *m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()
                + m_mesh_3.ED0F_B() *  m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();

	m_divEM1F = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B();
	m_curlEM1F = m_mesh_3.HS1F_B()*m_mesh_3.ED0F_B() *  m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();

}

void Decomposition::buildLaplacian2()
{
	m_laplacian2_NB = m_mesh_3.HS2F_NB().cwiseInverse() * m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_NB()
		+ m_mesh_3.ED1F_NB() *m_mesh_3.HS1F_NB().cwiseInverse() *m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB();

	m_divEnergyMatrix = m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_NB();
	m_curlEnergyMatrix = m_mesh_3.ED1F_NB()*m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB();
}

void Decomposition::buildLaplacian3()
{
	m_laplacian3_NB = m_mesh_3.ED2F_NB() * m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F();
	m_laplacian3_B = m_mesh_3.ED2F_B() * m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F();
}


void Decomposition::setFormRandom()
{
	omega.resize(m_mesh_3.NumF());
	omega.setZero();

	EigenVector form;
	
	cout <<"Electric Field ..."<< endl;
	Mesh_3_Point_3 pcp(0.3, 0.3, 0.5);
	setPointChargeElectricField(form, pcp);
	omega += 0.8*form;
	pcp = Mesh_3_Point_3(0.6, 0.6, 0.6);
	setPointChargeElectricField(form, pcp, false);
	omega += 0.6*form;

	cout<<"Eigen Field ..."<<endl;
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB()*m_laplacian2_NB;

	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s2l2nb.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s2l2nb, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> hs;
	for (int i = 0; i < m_mesh_3.NumIF(); ++i) {
		hs.push_back(m_mesh_3.HS2F_NB().coeff(i, i));
	}

	int num = 60;
	std::vector<double> emptyEvals;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEvals, eigenFields, row, col, val, hs, num);

	for(int i=0; i<eigenFields.size(); ++i){
		normalize(eigenFields[i]);
	}	
	
	groupFields(eigenFields);
	
	EigenVector cpnt;
	cpnt.setZero(m_mesh_3.NumIF());
	
	cpnt += 0.3*eigenFields[curlFieldIdx[0]];
        cpnt += 0.6*eigenFields[curlFieldIdx[1]];
	cpnt += 0.7*eigenFields[curlFieldIdx[2]];
        cpnt += 0.4*eigenFields[curlFieldIdx[3]];
		

	cpnt += 0.2*eigenFields[curlFieldIdx[10]];
	cpnt += 1.2*eigenFields[curlFieldIdx[12]];
	cpnt += 0.5*eigenFields[curlFieldIdx[14]];
	cpnt += 0.2*eigenFields[curlFieldIdx[16]];
	cpnt += 0.2*eigenFields[curlFieldIdx[18]];	

	noBoundaryToBoundary2Form(cpnt);

	omega += cpnt;
}


void Decomposition::computeHarmonicKnot()
{
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB() * m_laplacian2_NB;

	double hm1 = m_mesh_3.NumHMLG1();
	//double hm2 = m_mesh_3.NumHMLG2();
	cout<<"HM1: "<<hm1<<endl;

	if (hm1 == 0) {
		harmonicKnot.resize(m_mesh_3.NumF());
		harmonicKnot.setZero();
	}
	else {
		// assemble matrix to be eigen decomposed
		std::vector<double> row, col;
		std::vector<double> val;
		for (int i = 0; i < s2l2nb.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(s2l2nb, i); iter; ++iter) {
				row.push_back(static_cast<double>(iter.row()) + 1);
				col.push_back(static_cast<double>(iter.col()) + 1);
				val.push_back(static_cast<double>(iter.value()));
			}
		}
		std::vector<double> hs;
		for (int i = 0; i < m_mesh_3.NumIF(); ++i) {
			hs.push_back(m_mesh_3.HS2F_NB().coeff(i, i));
		}

		std::vector<double> emptyEval;
		matlabEIGS(emptyEval, harmonicKnotBasis, row, col, val, hs, hm1);

		//turn original form to no boundary
		EigenVector omegaNB = omega;
		boundaryToNoBoundary2Form(omegaNB);

		// we have harmonic basis
		EigenVector harmonicKnotNB;
		harmonicKnotNB.resize(m_mesh_3.NumIF());
		harmonicKnotNB.setZero();
		for (int i = 0; i < harmonicKnotBasis.size(); ++i) {
			double c = (omegaNB.transpose()*m_mesh_3.HS2F_NB()*harmonicKnotBasis[i])[0]
				/ (harmonicKnotBasis[i].transpose() * m_mesh_3.HS2F_NB()*harmonicKnotBasis[i])[0];
			harmonicKnotCoeffs.push_back(c);
			harmonicKnotNB += c * harmonicKnotBasis[i];
		}
		
		
		harmonicKnot = harmonicKnotNB;
		noBoundaryToBoundary2Form(harmonicKnot);
		
		// compute vector potential
		// need to fix beta_1 rank deficiency
		EigenSpMat s1l1 = m_mesh_3.HS1F_B()*m_laplacian1_B;
		EigenVector b = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*harmonicKnot;
			
		std::vector<double> row_vp, col_vp;
		std::vector<double> val_vp;
		for (int i = 0; i < s1l1.outerSize(); ++i) {
                        for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
                                row_vp.push_back(static_cast<double>(iter.row()) + 1);
                                col_vp.push_back(static_cast<double>(iter.col()) + 1);
                                val_vp.push_back(static_cast<double>(iter.value()));
                        }
                }
                std::vector<double> hs_vp;
                for (int i = 0; i < m_mesh_3.NumE(); ++i) {
                        hs_vp.push_back(m_mesh_3.HS1F_B().coeff(i, i));
                }

		std::vector<double> emptyEval_vp;
                std::vector<EigenVector> eigenFields_vp;
                matlabEIGS(emptyEval_vp, eigenFields_vp, row_vp, col_vp, val_vp, hs_vp, hm1);
		
		std::set<int> anchors_vp;
                selectAnchors(eigenFields_vp, anchors_vp);
                correctRankDeficiency(s1l1, b, anchors_vp);	
	
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_vp;
		solver_vp.compute(s1l1);
		EigenVector x = solver_vp.solve(b);

		harmonicKnotVP = x;
		
		// project
		cout<<"Potential tangential harmonic coeffs: "<<endl;
		for(int i=0; i<eigenFields_vp.size(); ++i){
			double c = harmonicKnotVP.dot(m_mesh_3.HS1F_B()*eigenFields_vp[i])
				/ eigenFields_vp[i].dot(m_mesh_3.HS1F_B()*eigenFields_vp[i]);
			harmonicKnotVP -= c*eigenFields_vp[i];
			//cout<<c<<endl;
		}

		// accuracy
		//double divEnergy = harmonicKnotVP.dot(m_divEM1F*harmonicKnotVP);
		//cout<<"Harmonic Knot Potential VP divergence Energy: "<<divEnergy<<endl;

	}
}

void Decomposition::computeHarmonicGradient()
{
	EigenSpMat A = m_laplacian3_B;
	EigenVector b;

	std::vector<EigenVector> tmpBasis;
	for (int i = 1; i < m_mesh_3.NumBoundary(); ++i) {
		assignBoundaryFacetPotential(i);

		b.setZero(m_mesh_3.NumC());
		for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
			ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
			if (ci->subdomain_index() == 0)
				continue;

			for (int j = 0; j < 4; ++j) {
				My_facet mf(ci, j);
				Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);

				if (m_mesh_3.FacetOnBoundary(fi)) {
					b(m_mesh_3.CellIdx(ci)) += boundaryFacetPotential[fi] * (1 / m_mesh_3.HS2F_B().coeff(m_mesh_3.FacetIdx(fi), m_mesh_3.FacetIdx(fi)));
				}
			}
		}

		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		EigenVector x = solver.solve(b); // potential 3 form for coexact component
		harmonicGradientPotentialBasis.push_back(x);

		EigenVector hmnc = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*x;
		for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
			fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
			if (!m_mesh_3.FacetValid(fi))
				continue;

			if (m_mesh_3.FacetOnBoundary(fi)) {
				int even = 1;
				if (fi->second % 2 == 1)
					even = -1;

				std::vector<int> lIdx;
				for (int j = 0; j < 4; ++j) {
					if (j != fi->second)
						lIdx.push_back(j);
				}

				std::vector<int> gIdx;
				gIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(lIdx[0])));
				gIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(lIdx[1])));
				gIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(lIdx[2])));

				for (int k = 0; k < gIdx.size(); ++k) {
					for (int j = k + 1; j < gIdx.size(); ++j) {
						if (gIdx[k] > gIdx[j]) {
							std::swap(gIdx[k], gIdx[j]);
							even = -even;
						}
					}
				}

				if (fi->first->subdomain_index() == 0)
					even = -even;

				hmnc(m_mesh_3.FacetIdx(fi)) -= even*boundaryFacetPotential[fi] * (1 / m_mesh_3.HS2F_B().coeff(m_mesh_3.FacetIdx(fi), m_mesh_3.FacetIdx(fi)));
			}
		}
		tmpBasis.push_back(hmnc);
	}

	// Gram schmidt
	for (int i = 0; i < tmpBasis.size(); ++i) {
		EigenVector orthoBasis = tmpBasis[i];

		for (int j = 0; j < harmonicGradientBasis.size(); ++j) {
			// compute coefficient
			double numerator, denominator;
			EigenVector currentBasis, normalizedBasis;
			currentBasis = tmpBasis[i];
			normalizedBasis = harmonicGradientBasis[j];

			numerator = currentBasis.transpose()*m_mesh_3.HS2F_B()*normalizedBasis;
			denominator = normalizedBasis.transpose()*m_mesh_3.HS2F_B()*normalizedBasis;

			// subtract
			orthoBasis -= (numerator / denominator)*normalizedBasis;
		}
		harmonicGradientBasis.push_back(orthoBasis);
	}
	
	harmonicGradient.setZero(m_mesh_3.NumF());
	for (int i = 0; i < harmonicGradientBasis.size(); ++i) {
		double numerator = (omega.transpose()*m_mesh_3.HS2F_B()*harmonicGradientBasis[i]);
		double denominator = (harmonicGradientBasis[i].transpose() * m_mesh_3.HS2F_B()*harmonicGradientBasis[i]);
		double c = numerator / denominator;
		cout << numerator << " " << denominator << endl;
		harmonicGradientCoeffs.push_back(c);
		
		harmonicGradient += c * harmonicGradientBasis[i];
	}
	
	if(m_mesh_3.NumHMLG2() != 0){
		// compute scalar potential
		EigenVector tmpGG_NB = harmonicGradient;
        	boundaryToNoBoundary2Form(tmpGG_NB);
        	EigenSpMat s3l3nb_sp = m_laplacian3_NB;
        	EigenVector b_sp = m_mesh_3.ED2F_NB()*tmpGG_NB;

		
        	for (int i = 0; i < s3l3nb_sp.outerSize(); ++i) {
                	for (EigenSpMat::InnerIterator iter(s3l3nb_sp, i); iter; ++iter) {
                        	if (iter.row() == 0 && iter.col() == 0) {
                                	iter.valueRef() = 1;
                        	}
                        	else if (iter.row() == 0 || iter.col() == 0) {
                                	iter.valueRef() = 0;
                        	}
                	}
        	}
		b_sp(0)=0;	
		
		
        	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_sp;
        	solver_sp.compute(s3l3nb_sp);
        	EigenVector x_vp = solver_sp.solve(b_sp);
	
		harmonicGradientSP = m_mesh_3.HS3F()*x_vp;
	}
}

void Decomposition::computeFluxlessKnot()
{
	// by solving potentials
	EigenSpMat s1l1nb = m_mesh_3.HS1F_NB() * m_laplacian1_NB;

	double hm2 = m_mesh_3.NumHMLG2();

	if (hm2 == 0) {
		EigenVector omegaNB = omega;
		boundaryToNoBoundary2Form(omegaNB);

		EigenVector b = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*omegaNB;

		// solve
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(s1l1nb);

		EigenVector x = solver.solve(b); // potential 1-form for exact component
		fluxlessKnotVP = x;
		noBoundaryToBoundary1Form(fluxlessKnotVP);
		//convert1Form(vp, vectorPotential);

		fluxlessKnot = m_mesh_3.ED1F_NB() * x;
		noBoundaryToBoundary2Form(fluxlessKnot);
		//cout<< fluxlessKnot <<endl;
	}
	else {
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
		matlabEIGS(emptyEval, eigenFields, row, col, val, hs, hm2);

		EigenVector omegaNB = omega;
		boundaryToNoBoundary2Form(omegaNB);

		EigenVector b = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*omegaNB;

		std::set<int> anchors;
		selectAnchors(eigenFields, anchors);
		correctRankDeficiency(s1l1nb, b, anchors);

		// solve
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(s1l1nb);

		EigenVector x = solver.solve(b); // potential 1-form for exact component
		fluxlessKnotVP = x;
		noBoundaryToBoundary1Form(fluxlessKnotVP);

		fluxlessKnot = m_mesh_3.ED1F_NB() * x;
		noBoundaryToBoundary2Form(fluxlessKnot);
	}
}

void Decomposition::computeGroundedGradient()
{
	EigenSpMat A = m_laplacian3_B;
	EigenVector b = m_mesh_3.ED2F_B()*omega;

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	EigenVector x = solver.solve(b); // potential 3 form

	groundedGradientSP = m_mesh_3.HS3F()*x;

	groundedGradient = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*x;
}

void Decomposition::computeCurlyGradient()
{
	curlyGradient = omega - harmonicKnot - harmonicGradient - fluxlessKnot - groundedGradient;
}

void Decomposition::testOrthogonality()
{
	cout<<harmonicKnot.dot(m_mesh_3.HS2F_B()*harmonicGradient)<<endl;
	cout<<harmonicKnot.dot(m_mesh_3.HS2F_B()*fluxlessKnot)<<endl;
	cout<<harmonicKnot.dot(m_mesh_3.HS2F_B()*groundedGradient)<<endl;
	cout<<harmonicKnot.dot(m_mesh_3.HS2F_B()*curlyGradient)<<endl;
	
	cout<<harmonicGradient.dot(m_mesh_3.HS2F_B()*fluxlessKnot)<<endl;
	cout<<harmonicGradient.dot(m_mesh_3.HS2F_B()*groundedGradient)<<endl;
	cout<<harmonicGradient.dot(m_mesh_3.HS2F_B()*curlyGradient)<<endl;
	
	cout<<fluxlessKnot.dot(m_mesh_3.HS2F_B()*groundedGradient)<<endl;
	cout<<fluxlessKnot.dot(m_mesh_3.HS2F_B()*curlyGradient)<<endl;

	cout<<groundedGradient.dot(m_mesh_3.HS2F_B()*curlyGradient)<<endl;

}

void Decomposition::integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf, 
	double dsty, double lth, int stps, int md, std::string sfx)
{
	RKIntegrator integrator(m_mesh_3);

	
	integrator.init(dsty, lth, stps, md);
	integrator.initColorField(cf);
	integrator.convert(vf);
	integrator.integrate();
	integrator.write(sfx);
}

void Decomposition::matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
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

void Decomposition::assignBoundaryFacetPotential(int idx)
{
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;

		if (m_mesh_3.FacetOnBoundary(fi)) {
			int pieceIdx;
			for (int i = 0; i < 4; ++i) {
				if (i != fi->second) {
					pieceIdx = m_mesh_3.VertexBoundaryIdx(fi->first->vertex(i));
					break;
				}
			}

			if (pieceIdx == idx) {
				boundaryFacetPotential[fi] = 1;
			}
			else{
				boundaryFacetPotential[fi] = 0;
			}
		}
	}
}

void Decomposition::groupFields(std::vector<EigenVector> eigenFields)
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

	cout<<"Curl Fields: "<<curlIdx<<endl;
	cout<<"Div Fields: "<<divIdx<<endl;
}

void Decomposition::convertForms()
{

	convert2Form(omega, omegaVF);
	if(m_mesh_3.NumHMLG1()!=0)
		convert2Form(harmonicKnot, harmonicKnotVF);
	if(m_mesh_3.NumHMLG2()!=0)
		convert2Form(harmonicGradient, harmonicGradientVF);
	convert2Form(fluxlessKnot, fluxlessKnotVF);
	convert2Form(groundedGradient, groundedGradientVF);
	convert2Form(curlyGradient, curlyGradientVF);
	cout<<"Field form conversion done ..."<<endl;
}

void Decomposition::computeArrows()
{
	assignColorStrength(omegaVF, 0);
	if(m_mesh_3.NumHMLG1()!=0)
		assignColorStrength(harmonicKnotVF, 1);
	if(m_mesh_3.NumHMLG2()!=0)
		assignColorStrength(harmonicGradientVF, 2);
	assignColorStrength(fluxlessKnotVF, 3);
	assignColorStrength(groundedGradientVF, 4);
	assignColorStrength(curlyGradientVF, 5);
	
	rescaleArrowsWithVF(omegaVF);
	if(m_mesh_3.NumHMLG1()!=0)
		rescaleArrowsWithVF(harmonicKnotVF);
	if(m_mesh_3.NumHMLG2()!=0)
		rescaleArrowsWithVF(harmonicGradientVF);
	rescaleArrowsWithVF(fluxlessKnotVF);
	rescaleArrowsWithVF(groundedGradientVF);
	rescaleArrowsWithVF(curlyGradientVF);
}

void Decomposition::computeCrossSection(int mode)
{
	vertices.clear();
	colors.clear();
	vertexMap.clear();
	indices.clear();
	normals.clear();
	realBoundary.clear();

	double cutThd = -10;

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
					if (fi->first->vertex(i)->point().point().y() > cutThd)
						side = true;
					else {
						side = false;
						break;
					}
				}
			}

			if (fi->first->subdomain_index() == 0) {
				if (m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point().y() < cutThd)
					side = false;

				validCell = fi->first->neighbor(fi->second);
			}
			else
			{
				if (fi->first->vertex(fi->second)->point().point().y() < cutThd)
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
							if(mode == 0)
								colors.push_back(harmonicGradientSPCF[fi->first->vertex(i)]);
							else if(mode == 1) 
								colors.push_back(groundedGradientSPCF[fi->first->vertex(i)]);
							else 
								colors.push_back(curlyGradientSPCF[fi->first->vertex(i)]);							

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

					if (fi->first->vertex(i)->point().point().y() > cutThd)
						numTargetSide++;
				}
			}

			if (infinite)
				continue;

			bool onBoundary = false;
			Mesh_3_Cell_iterator boundaryCell;
			Mesh_3_Cell_iterator validCell; // for extracting vector field.
			if (numTargetSide == 3) {
				if (fi->first->vertex(fi->second)->point().point().y() < cutThd) {
					boundaryCell = fi->first;
					validCell = fi->first->neighbor(fi->second);
					onBoundary = !onBoundary;
				}

				if (m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point().y() < cutThd) {
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

							if(mode == 0)
                                                                colors.push_back(harmonicGradientSPCF[fi->first->vertex(i)]);
                                                        else if(mode == 1)
                                                                colors.push_back(groundedGradientSPCF[fi->first->vertex(i)]);
                                                        else
                                                                colors.push_back(curlyGradientSPCF[fi->first->vertex(i)]);  
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

void Decomposition::writeVFasVTK(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, std::string sfx){
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

void Decomposition::rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

	
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Vector_3 v = vf[ci];
		double len = sqrt(v*v);
		v /= len;

		if (len > pivotLen)
			len = pivotLen;
		else if(len < 0.5*pivotLen)
			len = 0.5*pivotLen;
			
		len = (aveLen)*(len / (pivotLen));

		vf[ci] = v * len;
	}
}

void Decomposition::setPointChargeElectricField(EigenVector& field, Mesh_3_Point_3 rp, bool ori)
{
	// rp: relative position of bounding box
	Mesh_3_Point_3 bbMin, bbMax;
	bbMin = m_mesh_3.BoundingBoxMin();
	bbMax = m_mesh_3.BoundingBoxMax();

	Mesh_3_Point_3 p((1 - rp.x())*bbMin.x() + rp.x()*bbMax.x(), (1 - rp.y())*bbMin.y() + rp.y()*bbMax.y(), (1 - rp.z())*bbMin.z() + rp.z()*bbMax.z());

	field.setZero(m_mesh_3.NumF());
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;
		
		Mesh_3_Point_3 fc = m_mesh_3.FacetCC(fi); // cell circum center
		Mesh_3_Vector_3 d = fc - p;;
		double sq = d*d;
		d /= sqrt(sq);

		std::vector<int> vIdx;
		for (int i = 0; i < 4; ++i) {
			if (i != fi->second) {
				vIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i))); // get global vertex idx
			}
		}
		std::sort(vIdx.begin(), vIdx.end()); // set consistent orientation, increasing order of vertex indices

		Mesh_3_Vector_3 v1, v2;
		v1 = m_mesh_3.IdxToVertex(vIdx[1])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();
		v2 = m_mesh_3.IdxToVertex(vIdx[2])->point().point() - m_mesh_3.IdxToVertex(vIdx[1])->point().point();
		Mesh_3_Vector_3 av = CGAL::cross_product(v1, v2) / 2; // area normal vector

		if(ori)
			field(m_mesh_3.FacetIdx(fi)) = av * d;
		else 
			field(m_mesh_3.FacetIdx(fi)) = -av * d;

	}

	normalize(field);
}

void Decomposition::setCurrentMagneticField(EigenVector& field, Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, bool ori)
{
	// rp: relative position of bounding box
	Mesh_3_Point_3 bbMin, bbMax;
	bbMin = m_mesh_3.BoundingBoxMin();
	bbMax = m_mesh_3.BoundingBoxMax();

	Mesh_3_Point_3 p((1 - rp.x())*bbMin.x() + rp.x()*bbMax.x(), (1 - rp.y())*bbMin.y() + rp.y()*bbMax.y(), (1 - rp.z())*bbMin.z() + rp.z()*bbMax.z());

	field.setZero(m_mesh_3.NumF());
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;
		
		Mesh_3_Point_3 fc = m_mesh_3.FacetCC(fi); // facet circum center
		Mesh_3_Vector_3 vpfc = fc - p; // vector p to fc
		Mesh_3_Vector_3 vpfcp = (vpfc * d)*d; // vector p to fc protected on d 
		Mesh_3_Vector_3 vpd = vpfc - vpfcp; // vector p to direction d

		Mesh_3_Vector_3 fD = CGAL::cross_product(d, vpd); // field direction
		fD /= sqrt(fD*fD); // normalize
		
		std::vector<int> vIdx;
		for (int i = 0; i < 4; ++i) {
			if (i != fi->second) {
				vIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i))); // get global vertex idx
			}
		}
		std::sort(vIdx.begin(), vIdx.end()); // set consistent orientation, increasing order of vertex indices

		Mesh_3_Vector_3 v1, v2;
		v1 = m_mesh_3.IdxToVertex(vIdx[1])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();
		v2 = m_mesh_3.IdxToVertex(vIdx[2])->point().point() - m_mesh_3.IdxToVertex(vIdx[1])->point().point();
		Mesh_3_Vector_3 av = CGAL::cross_product(v1, v2) / 2; // area normal vector
	
		if(ori)
			field(m_mesh_3.FacetIdx(fi)) = av * fD;
		else
			field(m_mesh_3.FacetIdx(fi)) = -av * fD;
	}

	normalize(field);
}

void Decomposition::normalize(EigenVector& form)
{
	form = form / sqrt(form.dot(form));
}

void Decomposition::boundaryToNoBoundary1Form(EigenVector& form)
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

void Decomposition::noBoundaryToBoundary1Form(EigenVector& form)
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

void Decomposition::boundaryToNoBoundary2Form(EigenVector& form)
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

void Decomposition::noBoundaryToBoundary2Form(EigenVector& form)
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

void Decomposition::convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

void Decomposition::convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

void Decomposition::convert3Form(EigenVector form, std::map<Mesh_3_Vertex_iterator, double>& sf)
{
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

			p += form[m_mesh_3.CellIdx(cells[i])] * m_mesh_3.CellPrimal(cells[i]);
			vol += m_mesh_3.CellPrimal(cells[i]);
		}

		p /= vol;

		sf[vi] = p;
	}
}

void Decomposition::computeScalarPotential(EigenVector form)
{
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

			p += form[m_mesh_3.CellIdx(cells[i])] * m_mesh_3.CellPrimal(cells[i]);
			vol += m_mesh_3.CellPrimal(cells[i]);
		}

		p /= vol;

		scalarPotential[vi] = p;
	}
}

void Decomposition::selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors)
{
	cout << eigenfields.size() << endl;

	std::random_device rd;
	std::mt19937_64 mt(rd());
	std::uniform_int_distribution<int> distribution(0, m_mesh_3.NumIE());

	double condNumber = 0;

	int size = static_cast<int>(eigenfields.size());
	int iteration = 10;
	int flag = 0;
	for (int i = 0; i < iteration; ++i) {
		std::set<int> rows;

		do {
			int d = distribution(mt);
			rows.insert(d);
		} while (rows.size() < size);
		EigenMatrix mat(size, size);
		int j = 0;
		for (auto iter = rows.begin(); iter != rows.end(); ++iter, ++j) {
			for (int k = 0; k < size; ++k) {
				mat(j, k) = eigenfields[k](*iter);
			}
		}

		if (mat.determinant() < 10e-16) {
			continue;
		}
		else {
			Eigen::JacobiSVD<EigenMatrix> svd(mat);
			double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
			
			if(flag == 0){
				condNumber = fabs(cond);
				anchors = rows;
			}
			
			cout<<"Condition Number: "<<cond<<endl;
			if (fabs(cond) < fabs(condNumber)) {
				condNumber = fabs(cond);
				anchors = rows;
			}
		}
	}

	//for (auto iter = anchors.begin(); iter != anchors.end(); ++iter) {
	//	cout<<*iter<<endl;
	//}
}

void Decomposition::correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors)
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

void Decomposition::assignColorStrength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode)
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

void Decomposition::assignColorStrengthVP(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vp, int mode){	
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

                        Mesh_3_Vector_3 v = vp[cells[i]];
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

                Mesh_3_Vector_3 col = assignColorVP(ratio);

                if (mode == 0) { 
                        harmonicKnotVPCF[vi] = col;
                }
                else if (mode == 1) {
                        fluxlessKnotVPCF[vi] = col;
                }
                else if (mode == 2) {
                        curlyGradientVPCF[vi] = col;
                }
        }		

	if (mode == 0) {
                colorField = harmonicKnotVPCF;
        }
        else if (mode == 1) {
                colorField = fluxlessKnotVPCF;
        }
        else if (mode == 2) {
                colorField =  curlyGradientVPCF;
        }

}

void Decomposition::assignColorStrengthSP(std::map<Mesh_3_Vertex_iterator, double> sf, int mode){
	double minv = sf[m_mesh_3.C3T3().triangulation().vertices_begin()];
	double maxv = sf[m_mesh_3.C3T3().triangulation().vertices_begin()];
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;
		
		if (sf[vi] < minv)
			minv = sf[vi];
		if (sf[vi] > maxv)
			maxv = sf[vi];
	}

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		double diff = maxv - minv;

		double d = sf[vi] - (minv + 0*diff);

		double ratio = d / (diff);

		if (ratio < 0)
			ratio = 0;
		if (ratio > 1)
			ratio = 1;

		Mesh_3_Vector_3 col = assignColorSP(ratio);

		if (mode == 0) {
			harmonicGradientSPCF[vi] = col;
		}
		else if (mode == 1) {
			groundedGradientSPCF[vi] = col;
		}
		else if (mode == 2) {
			curlyGradientSPCF[vi] = col;
		}
	}

	if (mode == 0) {
		colorField = harmonicGradientSPCF;
	}
	else if (mode == 1) {
		colorField = groundedGradientSPCF;
	}
	else if (mode == 2) {
		colorField = curlyGradientSPCF;
	}
}

void Decomposition::assignHarmonicGradientCF()
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

		harmonicGradientCF[vi] = ratio * cmax + (1 - ratio)*cmin;
	}
	colorField = harmonicGradientCF;
}

Mesh_3_Vector_3 Decomposition::assignColor(double ratio)
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

Mesh_3_Vector_3 Decomposition::assignColorVP(double ratio)
{
	// 0.0 - 0.5 light pink to white
	// 0.5 - 1.0 white to light blue
	
	Mesh_3_Vector_3 color;

	Mesh_3_Vector_3 pink(0.9, 0.5, 0.7);
	Mesh_3_Vector_3 white(0.9, 0.9, 0.9);
	Mesh_3_Vector_3 blue(0.5, 0.7, 0.9);

	if(ratio < 0.5){
		double r = ratio / 0.5;
		color = (1-r)*blue+r*white;
	}	
	else{
		double r = (ratio - 0.5) / 0.5;
		color = (1-r)*white+r*pink;
	}

	return color;
}

Mesh_3_Vector_3 Decomposition::assignColorSP(double ratio)
{
	// 0.0 - 0.5 purple to white
	// 0.5 - 1.0 white to orange

	Mesh_3_Vector_3 color;

        Mesh_3_Vector_3 orange(0.9, 0.5, 0.3);
        Mesh_3_Vector_3 white(0.9, 0.9, 0.9);
        Mesh_3_Vector_3 purple(0.9, 0.3, 0.5);

        if(ratio < 0.5){
                double r = ratio / 0.5;
                color = (1-r)*purple+r*white;
        }
        else{
                double r = (ratio - 0.5) / 0.5;
                color = (1-r)*white+r*orange;
        }

        return color;
}

void Decomposition::sortedIndices(std::vector<int>& sIndex, std::vector<double> len)
{
	sIndex.resize(len.size(), 0);
	std::iota(sIndex.begin(), sIndex.end(), 0);

	std::sort(sIndex.begin(), sIndex.end(), [&len](int i1, int i2) {return len[i1] < len[i2]; });
}

