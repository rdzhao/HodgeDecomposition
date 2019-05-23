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

void Decomposition::init(double dsty, double sr, int stps, std::string sfx)
{
	m_density_ratio = dsty;
	m_step_ratio = sr;
	m_steps = stps;
	m_suffix = sfx;
}

void Decomposition::buildMeshFromSurface(std::string fn, double size)
{
	//Polyhedron polyhedron;
	//std::ifstream input(fn.c_str());
	//input >> polyhedron;
	//input.close();
	//Polyhedron_mesh_domain domain(polyhedron);

	//Mesh_criteria criteria(facet_angle = 25,
	//	facet_size = size,
	//	facet_distance = 0.1*size,
	//	cell_radius_edge_ratio = 2,
	//	cell_size = size);
	//std::cout << "Criteria ... " << std::endl;

	//m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::odt());
	//CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	//std::cout << "Triangulation ... " << std::endl;

	Polyhedron_with_feature polyhedron;
	std::ifstream input(fn.c_str());
	input >> polyhedron;
	input.close();
	Mesh_domain_with_feature domain(polyhedron);
	//domain.detect_features();

	Mesh_criteria criteria(/*edge_size = 1.2*size,*/
		facet_angle = 30,
		facet_size = size,
		facet_distance = 0.1*size,
		cell_radius_edge_ratio = 2,
		cell_size = size);
	std::cout << "Criteria ... " << std::endl;

	m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, features(domain), odt());
	//CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	std::cout << "Triangulation ... " << std::endl;

	m_mesh_3.preprocessing();

	cout << "Euler: " << m_mesh_3.NumV() - m_mesh_3.NumE() + m_mesh_3.NumF() - m_mesh_3.NumC() << endl;
}

void Decomposition::buildDECOperators()
{
	//buildExteriorDerivative0Form();
	//buildExteriorDerivative1Form();
	//buildExteriorDerivative2Form();

	//buildHodgeStar0Form();
	//buildHodgeStar1Form();
	//buildHodgeStar2Form();
	//buildHodgeStar3Form();
}

void Decomposition::buildLaplacian()
{
	cout << "@@@@@" << endl;
	buildLaplacian1();
	cout << "@@@@@" << endl;
	buildLaplacian2(); 
	cout << "@@@@@" << endl;
	buildLaplacian3();
}

void Decomposition::decompose()
{
	// ----------------------
	setFormRandom();
	//sampleInput();
	//sampleFromJHTDB();

	//testEnergy();
	//testNonDiagonalHodgeStar();
	
	
	//cout << "Compute TC ..." << endl;
	//computeTC();
	//cout << "Compute NH ..." << endl;
	//computeNH();
	//cout << "Compute TH ..." << endl;
	//computeTH();
	//cout << "Compute CH ..." << endl;
	//computeCH();
	//cout << "Compute NG ..." << endl;
	//computeNG();
	
	//testProjectionMatrix();
	

	cout<<"Harmonic Knot ..."<<endl;
	computeHarmonicKnot();
	cout << "Harmonic Gradient ..." << endl;
	computeHarmonicGradient();
	cout << "Fluxless Knot ..." << endl;
	computeFluxlessKnot();
	cout << "Grounded Gradient ..." << endl;
	computeGroundedGradient();
	cout << "Curly Gradient ..." << endl;
	computeCurlyGradient();

	//evaluateAccuracy();
	//testOrthogonality();
	//computeError();
	//computeL2Norm();
	
}

void Decomposition::visualize()
{
	//assignHarmonicGradientCF();

	//harmonicKnot += harmonicGradient;
	//curlyGradient += harmonicGradient;

	convertForms();
	computeArrows();
	cout << "Vector field visualization done ..." << endl;

	if (m_mesh_3.NumHMLG2() != 0) {
		computeCrossSection(0);
		writeMeshPLY("hg");
	}
	computeCrossSection(1);
	writeMeshPLY("gg");
	computeCrossSection(2);
	writeMeshPLY("cg");
	cout << "Scalar field visualization done ..." << endl;

	
	//writeMeshPLY();

	//testSurface();
	//testSurfaceAndArrowsPLY();
}

void Decomposition::integrate()
{
	cout << "Omega Field ..." << endl;
	integrateField(omegaVF, omegaCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_o");

	if (m_mesh_3.NumHMLG1() != 0) {
		cout << "Harmonic Knot ..." << endl;
		integrateField(harmonicKnotVF, harmonicKnotCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_hk");
	}

	if (m_mesh_3.NumHMLG2() != 0) {
		cout << "Harmonic Gradient ..." << endl;
		integrateField(harmonicGradientVF, harmonicGradientCF, m_density_ratio / 3, m_step_ratio, m_steps, 2, m_suffix + "_hg");
	}

	cout << "Fluxless Knot ..." << endl;
	integrateField(fluxlessKnotVF, fluxlessKnotCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_fk");

	cout << "Grounded Gradient ..." << endl;
	integrateField(groundedGradientVF, groundedGradientCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_gg");

	cout << "Curly Gradient ..." << endl;
	integrateField(curlyGradientVF, curlyGradientCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_cg");

	//// vector potential
	//if (m_mesh_3.NumHMLG1() != 0) {
	//	cout << "Harmonic Knot VP ..." << endl;
	//	integrateField(harmonicKnotVPVF, harmonicKnotVPCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_hkvp");
	//}

	//cout << "Fluxless Knot VP ..." << endl;
	//integrateField(fluxlessKnotVPVF, fluxlessKnotVPCF, m_density_ratio / 3, m_step_ratio, m_steps, 2, m_suffix + "_fkvp");

	//cout << "Curly Gradient VP ..." << endl;
	//integrateField(curlyGradientVPVF, curlyGradientVPCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_cgvp");

}

void Decomposition::write()
{
	//writeMesh();
	//writeMeshPLY();
}


void Decomposition::buildLaplacian1()
{
	//cout << m_mesh_3.HS1F_NB().innerSize() << " " << m_mesh_3.HS1F_NB().outerSize() << endl;

	m_laplacian1_NB = m_mesh_3.HS1F_NB().cwiseInverse() *m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB()
		+ m_mesh_3.ED0F_NB() *  m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB(); 

	m_laplacian1_B = m_mesh_3.HS1F_B().cwiseInverse() *m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()
		+ m_mesh_3.ED0F_B() *  m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();
}

void Decomposition::buildLaplacian2()
{
	m_laplacian2_NB = m_mesh_3.HS2F_NB().cwiseInverse() * m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_NB()
					+ m_mesh_3.ED1F_NB() *m_mesh_3.HS1F_NB().cwiseInverse() *m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB();
	m_laplacian2_B = m_mesh_3.HS2F_B().cwiseInverse() * m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_B()
		+ m_mesh_3.ED1F_B() *m_mesh_3.HS1F_B().cwiseInverse() *m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B();


	m_divEM = m_mesh_3.HS2F_NB().cwiseInverse() * m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_NB();
	m_curlEM = m_mesh_3.ED1F_NB() *m_mesh_3.HS1F_NB().cwiseInverse() *m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB();
}

void Decomposition::buildLaplacian3()
{
	m_laplacian3_B = m_mesh_3.ED2F_B() * m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F();
	m_laplacian3_NB = m_mesh_3.ED2F_NB() * m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F();
}

void Decomposition::setForm()
{
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB() * m_laplacian2_NB;

	Engine* eng;
	eng = engOpen("\0");

	std::string str;

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

	// transfer data
	mxArray *rowArr, *colArr, *valArr, *hsArr;
	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	hsArr = mxCreateDoubleMatrix(hs.size(), 1, mxREAL);
	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	memcpy(mxGetPr(hsArr), &hs[0], hs.size() * sizeof(double));
	engPutVariable(eng, "row", rowArr);
	engPutVariable(eng, "col", colArr);
	engPutVariable(eng, "val", valArr);
	engPutVariable(eng, "hs", hsArr);

	str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());

	str = "n=length(hs)";        engEvalString(eng, str.c_str());
	str = "B = spdiags(hs(:),0,n,n);";        engEvalString(eng, str.c_str());
	str = "[V,D]=eigs(A, B, 100, 'smallestabs')";        engEvalString(eng, str.c_str());

	mxArray* vv;
	vv = engGetVariable(eng, "V");
	double* v;
	v = reinterpret_cast<double*>(mxGetData(vv));
	const size_t* dims;
	dims = mxGetDimensions(vv);

	EigenVector ev; // no boundary
	ev.resize(dims[0]);

	std::vector<EigenVector> eigenFields;
	eigenFields.clear();
	for (int j = 0; j < dims[1]; ++j) {
		for (int i = 0; i < dims[0]; ++i) {
			ev(i) = v[i + j * dims[0]];
		}

		eigenFields.push_back(ev);
	}

	// randomly choose eigen basis for assembling omega.
	omega.setZero(m_mesh_3.NumIF());
	omega += eigenFields[3];

	noBoundaryToBoundary2Form(omega);
}

void Decomposition::setFormRandom()
{
	omega.resize(m_mesh_3.NumF());
	omega.setZero();

	EigenVector form;
	
	Mesh_3_Point_3 pcp(0.2, 0.65, 0.35);
	setPointChargeElectricField(form, pcp);
	omega += 0.5*form;
	pcp = Mesh_3_Point_3(0.8, 0.35, 0.65);
	setPointChargeElectricField(form, pcp, false);
	omega += 0.4*form;
	cout << "Electric Done ..." << endl;

	// normal example fields
	EigenSpMat s2l2 = m_mesh_3.HS2F_B()*m_laplacian2_B;
	std::vector<double> rowb, colb;
	std::vector<double> valb;
	for (int i = 0; i < s2l2.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s2l2, i); iter; ++iter) {
			rowb.push_back(static_cast<double>(iter.row()) + 1);
			colb.push_back(static_cast<double>(iter.col()) + 1);
			valb.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> hrowb, hcolb;
	std::vector<double> hvalb;
	for (int i = 0; i < m_mesh_3.HS2F_B().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_B(), i); iter; ++iter) {
			hrowb.push_back(static_cast<double>(iter.row()) + 1);
			hcolb.push_back(static_cast<double>(iter.col()) + 1);
			hvalb.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEvalsb;
	std::vector<EigenVector> eigenFieldsb;
	matlabEIGS(emptyEvalsb, eigenFieldsb, rowb, colb, valb, hrowb, hcolb, hvalb, 10);

	omega += 0.1*eigenFieldsb[0];
	omega += eigenFieldsb[2];
	//omega += eigenFieldsb[3];

	// tangential example fields
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB() * m_laplacian2_NB;
	std::vector<double> rownb, colnb;
	std::vector<double> valnb;
	for (int i = 0; i < s2l2nb.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s2l2nb, i); iter; ++iter) {
			rownb.push_back(static_cast<double>(iter.row()) + 1);
			colnb.push_back(static_cast<double>(iter.col()) + 1);
			valnb.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> hrownb, hcolnb;
	std::vector<double> hvalnb;
	for (int i = 0; i < m_mesh_3.HS2F_NB().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_NB(), i); iter; ++iter) {
			hrownb.push_back(static_cast<double>(iter.row()) + 1);
			hcolnb.push_back(static_cast<double>(iter.col()) + 1);
			hvalnb.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEvalsnb;
	std::vector<EigenVector> eigenFieldsnb;
	matlabEIGS(emptyEvalsnb, eigenFieldsnb, rownb, colnb, valnb, hrownb, hcolnb, hvalnb, 40);

	groupFields(eigenFieldsnb);

	EigenVector cpnt;
	cpnt.resize(m_mesh_3.NumIF());
	cpnt.setZero();

	cpnt += 0.1*eigenFieldsnb[0];

	cpnt += 0.5*eigenFieldsnb[curlFieldIdx[3]];
	cpnt += 0.8*eigenFieldsnb[curlFieldIdx[9]];
	//cpnt += 1.2*eigenFieldsnb[curlFieldIdx[8]];
	//cpnt += 0.5*eigenFieldsnb[curlFieldIdx[9]];

	noBoundaryToBoundary2Form(cpnt);

	omega += cpnt;

	//EigenVector divtan;
	//createDivTan(divtan);

	//divtan /= sqrt(divtan.dot(m_mesh_3.HS2F_B()*divtan));
	//divtan *= 0.8;
	//omega += divtan;

	//omega /= sqrt(omega.dot(omega));
}

void Decomposition::computeNG()
{
	EigenSpMat idF;
	idF.resize(m_mesh_3.NumF(), m_mesh_3.NumF());
	idF.setIdentity();

	EigenSpMat A = m_mesh_3.CHS3F()*m_mesh_3.ED2F_B()*idF*m_mesh_3.ED2F_B().transpose()*m_mesh_3.CHS3F();
	EigenVector b = m_mesh_3.CHS3F()*m_mesh_3.ED2F_B()*idF*m_mesh_3.CHS2F_B()*(omega);

	EigenVector x = matlabLinearSolveCholesky(A, b);

	EigenSpMat AA = m_mesh_3.CHS2F_B();
	EigenVector bb = m_mesh_3.ED2F_B().transpose()*m_mesh_3.CHS3F()*x;

	EigenVector xx = matlabLinearSolveCholesky(AA, bb);

	//EigenVector xx = m_mesh_3.ED2F_B().transpose()*m_mesh_3.CHS3F()*x;

	groundedGradient = xx;
}

void Decomposition::computeTC()
{
	EigenSpMat idV;
	idV.resize(m_mesh_3.NumV(), m_mesh_3.NumV());
	idV.setIdentity();
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (m_mesh_3.VertexOnBoundary(vi)) {
			idV.coeffRef(i, i) = 0;;
		}
	}
	

	EigenSpMat A = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*m_mesh_3.ED1F_B() +
		m_mesh_3.CHS1F_B()*m_mesh_3.ED0F_B()*idV*m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B();
	EigenVector b = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*omega;

	for (int i = 0; i < A.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
			int row = iter.row();
			int col = iter.col();

			Mesh_3_Edge_iterator eir = m_mesh_3.IdxToEdge(row);
			Mesh_3_Edge_iterator eic = m_mesh_3.IdxToEdge(col);

			if (m_mesh_3.EdgeOnBoundary(eir) || m_mesh_3.EdgeOnBoundary(eir)) {
				if (row == col)
					iter.valueRef() = 1;
				else
					iter.valueRef() = 0;
			}
		}
	}
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (m_mesh_3.EdgeOnBoundary(ei))
			b(i) = 0;
	}


	EigenVector x = matlabLinearSolveCholesky(A, b);

	fluxlessKnot = m_mesh_3.ED1F_B()*x;
	//noBoundaryToBoundary2Form(fluxlessKnot);
	omega -= fluxlessKnot;
}

void Decomposition::computeNH()
{
	int hm2 = m_mesh_3.NumHMLG2();

	if (hm2 == 0) {
		harmonicGradient.setZero(m_mesh_3.NumF());
	}
	else {
		EigenSpMat A = m_mesh_3.ED2F_B().transpose()*m_mesh_3.CHS3F()*m_mesh_3.ED2F_B()
			+ m_mesh_3.CHS2F_B()*m_mesh_3.ED1F_B() *m_mesh_3.HS1F_B().cwiseInverse() *m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B();
		EigenSpMat B = m_mesh_3.CHS2F_B();

		std::vector<double> row, col;
		std::vector<double> val;
		for (int i = 0; i < A.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
				row.push_back(static_cast<double>(iter.row()) + 1);
				col.push_back(static_cast<double>(iter.col()) + 1);
				val.push_back(static_cast<double>(iter.value()));
			}
		}

		std::vector<double> hrow, hcol;
		std::vector<double> hval;
		for (int i = 0; i < B.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(B, i); iter; ++iter) {
				hrow.push_back(static_cast<double>(iter.row()) + 1);
				hcol.push_back(static_cast<double>(iter.col()) + 1);
				hval.push_back(static_cast<double>(iter.value()));
			}
		}

		std::vector<double> eval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(eval, eigenFields, row, col, val, hrow, hcol, hval, hm2);

		EigenMatrix H;
		H.resize(m_mesh_3.NumF(), hm2);

		for (int i = 0; i < hm2; ++i) {
			EigenVector ev = eigenFields[i];
			H.col(i) = ev;
		}

		harmonicGradient = H * H.transpose()*(m_mesh_3.CHS2F_B()*omega);
		omega -= harmonicGradient;
	}
}

void Decomposition::computeTH()
{
	int hm1 = m_mesh_3.NumHMLG1();

	if (hm1 == 0) {
		harmonicKnot.setZero(m_mesh_3.NumF());
	}
	else {
		//EigenSpMat hs1fb_inv = m_mesh_3.HS1F_B().cwiseInverse();
		//for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		//	Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		//	if (m_mesh_3.EdgeOnBoundary(ei)) {
		//		hs1fb_inv.coeffRef(i, i) = 0;
		//	}
		//}
		EigenSpMat A = m_mesh_3.ED2F_NB().transpose()*m_mesh_3.CHS3F()*m_mesh_3.ED2F_NB()
			+ m_mesh_3.CHS2F_NB()*m_mesh_3.ED1F_NB() *m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.CHS2F_NB();
		EigenSpMat B = m_mesh_3.CHS2F_NB();
		
		std::vector<double> row, col;
		std::vector<double> val;
		for (int i = 0; i < A.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
				row.push_back(static_cast<double>(iter.row()) + 1);
				col.push_back(static_cast<double>(iter.col()) + 1);
				val.push_back(static_cast<double>(iter.value()));
			}
		}
		
		std::vector<double> hrow, hcol;
		std::vector<double> hval;
		for (int i = 0; i < B.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(B, i); iter; ++iter) {
				hrow.push_back(static_cast<double>(iter.row()) + 1);
				hcol.push_back(static_cast<double>(iter.col()) + 1);
				hval.push_back(static_cast<double>(iter.value()));
			}
		}
		
		std::vector<double> eval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(eval, eigenFields, row, col, val, hrow, hcol, hval, hm1);
		
		EigenMatrix H;
		H.resize(m_mesh_3.NumF(), hm1);

		for (int i = 0; i < hm1; ++i) {
			EigenVector ev = eigenFields[i];
			noBoundaryToBoundary2Form(ev);
			H.col(i) = ev;
		}

		harmonicKnot = H * H.transpose()*(m_mesh_3.CHS2F_B()*omega);
		omega -= harmonicKnot;
	}

	//int asd;
	//cin >> asd;
}

void Decomposition::computeCH()
{
	EigenSpMat idV;
	idV.resize(m_mesh_3.NumV(), m_mesh_3.NumV());
	idV.setIdentity();

	EigenSpMat A = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*m_mesh_3.ED1F_B() +
		m_mesh_3.CHS1F_B()*m_mesh_3.ED0F_B()*idV*m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B();
	EigenVector b = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*omega;

	EigenVector x = matlabLinearSolveCholesky(A, b);

	curlyGradient = m_mesh_3.ED1F_B()*x;

	omega -= curlyGradient;
}

void Decomposition::sampleFromJHTDB()
{	
	EigenMatrix points;
	points.resize(3, m_mesh_3.NumF());

	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);
		Mesh_3_Point_3 p = m_mesh_3.FacetCC(fi);

		points.coeffRef(0, i) = 2 * p.x()+1;
		points.coeffRef(1, i) = 2 * p.y();
		points.coeffRef(2, i) = 2 * p.z()+1;
	}

	JHTDB jh;
	EigenMatrix velocities;
	velocities = jh.getVelocities(points);

	// compute omega
	omega.resize(m_mesh_3.NumF());
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);

		Mesh_3_Point_3 fc = m_mesh_3.FacetCC(fi); // cell circum center
		Mesh_3_Vector_3 d(velocities.coeff(0, i), velocities.coeff(1, i), velocities.coeff(2, i));

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

		omega(i) = av * d;
	}
}

void Decomposition::computeHarmonicKnot()
{
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB() * m_laplacian2_NB;

	double hm1 = m_mesh_3.NumHMLG1();
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
		//std::vector<double> hs;
		//for (int i = 0; i < m_mesh_3.NumIF(); ++i) {
		//	hs.push_back(m_mesh_3.HS2F_NB().coeff(i, i));
		//}

		//std::vector<double> emptyEval;
		//matlabEIGS(emptyEval, harmonicKnotBasis, row, col, val, hs, hm1);
		std::vector<double> hrow, hcol;
		std::vector<double> hval;
		for (int i = 0; i < m_mesh_3.HS2F_NB().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_NB(), i); iter; ++iter) {
				hrow.push_back(static_cast<double>(iter.row()) + 1);
				hcol.push_back(static_cast<double>(iter.col()) + 1);
				hval.push_back(static_cast<double>(iter.value()));
			}
		}


		std::vector<double> emptyEval;
		matlabEIGS(emptyEval, harmonicKnotBasis, row, col, val, hrow, hcol, hval, hm1);


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
			cout << c << endl;
		}

		harmonicKnot = harmonicKnotNB;
		noBoundaryToBoundary2Form(harmonicKnot);


		cout << "Harmonic Knot potential ..." << endl;
		// compute vector potential for harmonic knot
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
		//std::vector<double> hs_vp;
		//for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		//	hs_vp.push_back(m_mesh_3.HS1F_B().coeff(i, i));
		//}

		//std::vector<double> emptyEval_vp;
		//std::vector<EigenVector> eigenFields_vp;
		//matlabEIGS(emptyEval_vp, eigenFields_vp, row_vp, col_vp, val_vp, hs_vp, hm1);
		std::vector<double> hrow_vp, hcol_vp;
		std::vector<double> hval_vp;
		for (int i = 0; i < m_mesh_3.HS1F_B().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_B(), i); iter; ++iter) {
				hrow_vp.push_back(static_cast<double>(iter.row()) + 1);
				hcol_vp.push_back(static_cast<double>(iter.col()) + 1);
				hval_vp.push_back(static_cast<double>(iter.value()));
			}
		}


		std::vector<double> emptyEval_vp;
		std::vector<EigenVector> eigenFields_vp;
		matlabEIGS(emptyEval_vp, eigenFields_vp, row_vp, col_vp, val_vp, hrow_vp, hcol_vp, hval_vp, hm1);


		std::set<int> anchors_vp;
		selectAnchors(eigenFields_vp, anchors_vp);
		correctRankDeficiency(s1l1, b, anchors_vp, m_mesh_3.HS1F_B(), eigenFields_vp);

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_vp;
		solver_vp.compute(s1l1);
		EigenVector x = solver_vp.solve(b);

		//EigenVector x = matlabLinearSolveCholesky(s1l1, b);

		harmonicKnotVP = x;

		// project
		cout << "Potential tangential harmonic coeffs: " << endl;
		for (int i = 0; i < eigenFields_vp.size(); ++i) {
			double c = harmonicKnotVP.dot(m_mesh_3.HS1F_B()*eigenFields_vp[i])
				/ eigenFields_vp[i].dot(m_mesh_3.HS1F_B()*eigenFields_vp[i]);
			harmonicKnotVP -= c * eigenFields_vp[i];
			cout << c << endl;
		}
	}
}

void Decomposition::computeHarmonicGradient()
{
	EigenSpMat A = m_laplacian3_B;
	EigenVector b;

	std::vector<EigenVector> tmpBasis;
	for (int i = 0; i < m_mesh_3.NumBoundary() - 1; ++i) {
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

				//cout << hmnc(m_mesh_3.FacetIdx(fi)) << " " << boundaryFacetPotential[fi] << " "<< fi->first->subdomain_index() <<" "
				//	<< even * boundaryFacetPotential[fi] * (1 / m_mesh_3.HS2F_B().coeff(m_mesh_3.FacetIdx(fi), m_mesh_3.FacetIdx(fi))) << endl;

				hmnc(m_mesh_3.FacetIdx(fi)) -= even*boundaryFacetPotential[fi] * (1 / m_mesh_3.HS2F_B().coeff(m_mesh_3.FacetIdx(fi), m_mesh_3.FacetIdx(fi)));
			}
		}
		tmpBasis.push_back(hmnc);
	}

	//// !!!!!!!!!!!!!! comment after
	//tmpBasis.clear();
	//EigenSpMat s2l2 = m_mesh_3.HS2F_B()*m_laplacian2_B;

	//std::vector<double> row, col;
	//std::vector<double> val;
	//for (int i = 0; i < s2l2.outerSize(); ++i) {
	//	for (EigenSpMat::InnerIterator iter(s2l2, i); iter; ++iter) {
	//		row.push_back(static_cast<double>(iter.row()) + 1);
	//		col.push_back(static_cast<double>(iter.col()) + 1);
	//		val.push_back(static_cast<double>(iter.value()));
	//	}
	//}
	//std::vector<double> hs;
	//for (int i = 0; i < m_mesh_3.NumF(); ++i) {
	//	hs.push_back(m_mesh_3.HS2F_B().coeff(i, i));
	//}

	//std::vector<double> emptyEval;
	//std::vector<EigenVector> eigenFields;
	//matlabEIGS(emptyEval, eigenFields, row, col, val, hs, 1);

	//tmpBasis.push_back(eigenFields[0]);

	//for (int i = 0; i < m_mesh_3.NumF(); ++i) {
	//	if (omega(i)*tmpBasis[0](i) > 0)
	//	//if(m_mesh_3.FacetOnBoundary(m_mesh_3.IdxToFacet(i)))
	//		cout << omega(i) << " " << tmpBasis[0](i) << " " << m_mesh_3.FacetOnBoundary(m_mesh_3.IdxToFacet(i)) << endl;
	//}

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

	//cout<<"Orthogonal? :"<<endl;
	//cout << tmpBasis[0].transpose() * m_mesh_3.HS2F_B()*tmpBasis[1] << endl;
	//cout<< harmonicGradientBasis[0].transpose() * m_mesh_3.HS2F_B()*harmonicGradientBasis[1] <<endl;

	// compute 
	//cout << omega.size() << endl;
	//cout << m_mesh_3.HS2F_B().outerSize() << " " << m_mesh_3.HS2F_B().innerSize() << endl;
	
	cout<<"HG coeffs: "<<endl;
	harmonicGradient.setZero(m_mesh_3.NumF());
	for (int i = 0; i < harmonicGradientBasis.size(); ++i) {
		double numerator = (omega.transpose()*m_mesh_3.HS2F_B()*harmonicGradientBasis[i]);
		double denominator = (harmonicGradientBasis[i].transpose() * m_mesh_3.HS2F_B()*harmonicGradientBasis[i]);
		double c = numerator / denominator;
		//cout << numerator << " " << denominator << endl;
		harmonicGradientCoeffs.push_back(c);
		
		cout << c << endl;
		harmonicGradient += c * harmonicGradientBasis[i];
	}

	if (m_mesh_3.NumHMLG2() != 0) {
		// compute scalar potential
		EigenVector tmpGG_NB = harmonicGradient;
		boundaryToNoBoundary2Form(tmpGG_NB);
		EigenSpMat s3l3nb_sp =  m_laplacian3_NB;
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
		b_sp(0) = 0;

		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_sp;
		solver_sp.compute(s3l3nb_sp);
		EigenVector x_vp = solver_sp.solve(b_sp);

		//EigenVector x_vp = matlabLinearSolveCholesky(s3l3nb_sp, b_sp);

		//computeScalarPotential(m_mesh_3.HS3F()*x_vp);
		harmonicGradientSP = m_mesh_3.HS3F()*x_vp;
		//harmonicGradient.setZero();
		//harmonicGradient = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*x_vp;
		//harmonicGradientSP = harmonicGradientPotentialBasis[0];

	}

}

void Decomposition::computeFluxlessKnot()
{
	// by solving potentials
	EigenSpMat s1l1nb = m_mesh_3.HS1F_NB() * m_laplacian1_NB;

	double hm2 = m_mesh_3.NumHMLG2();

	if (hm2 == 0) {
		cout<<"###########"<<endl;
		EigenVector omegaNB = omega;
		boundaryToNoBoundary2Form(omegaNB);

		EigenVector b = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*omegaNB;

		// solve
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(s1l1nb);
		EigenVector x = solver.solve(b); // potential 1-form for exact component
		
		//EigenVector x = matlabLinearSolveCholesky(s1l1nb, b);

		fluxlessKnotVP = x;
		noBoundaryToBoundary1Form(fluxlessKnotVP);

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
		//std::vector<double> hs;
		//for (int i = 0; i < m_mesh_3.NumIE(); ++i) {
		//	hs.push_back(m_mesh_3.HS1F_NB().coeff(i, i));
		//}

		//std::vector<double> emptyEval;
		//std::vector<EigenVector> eigenFields;
		//matlabEIGS(emptyEval, eigenFields, row, col, val, hs, hm2);

		std::vector<double> hrow, hcol;
		std::vector<double> hval;
		for (int i = 0; i < m_mesh_3.HS1F_NB().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_NB(), i); iter; ++iter) {
				hrow.push_back(static_cast<double>(iter.row()) + 1);
				hcol.push_back(static_cast<double>(iter.col()) + 1);
				hval.push_back(static_cast<double>(iter.value()));
			}
		}


		std::vector<double> emptyEval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, eigenFields, row, col, val, hrow, hcol, hval, hm2);

		EigenVector omegaNB = omega;
		boundaryToNoBoundary2Form(omegaNB);

		EigenVector b = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*omegaNB;

		std::set<int> anchors;
		EigenSpMat s1l1nb_fixed = s1l1nb;
		EigenVector b_fixed = b;
		
		selectAnchors(eigenFields, anchors);
		correctRankDeficiency(s1l1nb_fixed, b_fixed, anchors, m_mesh_3.HS1F_NB(), eigenFields);

		// solve
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(s1l1nb_fixed);
		EigenVector x = solver.solve(b_fixed); // potential 1-form for exact component
		
		//EigenVector x = matlabLinearSolveCholesky(s1l1nb_fixed, b_fixed);
		
		fluxlessKnotVP = x;
		noBoundaryToBoundary1Form(fluxlessKnotVP);

		 // check accuracy ....
		int r = eigenFields[0].size();
		int c = anchors.size();
		EigenSpMat H(r, c);
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				H.coeffRef(i, j) = eigenFields[j](i);
			}
		}

		EigenVector b_proj = b;
		b_proj -= m_mesh_3.HS1F_NB()*H*H.transpose()*b;

		//check accuracy
		EigenVector DM = s1l1nb * x - b_proj;
		double err = DM.dot(m_mesh_3.HS1F_NB().cwiseInverse()*DM) / (b_proj.dot(m_mesh_3.HS1F_NB().cwiseInverse()*b_proj));
		cout << "Accuracy of Solution: " << err << endl;
		m_accu = err;
		

		fluxlessKnot = m_mesh_3.ED1F_NB() * x;
		noBoundaryToBoundary2Form(fluxlessKnot);
	}
}

void Decomposition::computeGroundedGradient()
{
	EigenSpMat A =  m_laplacian3_B;
	EigenVector b = m_mesh_3.ED2F_B()*omega;

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	EigenVector x = solver.solve(b); // potential 3 form

	//EigenVector x = matlabLinearSolveCholesky(A, b);

	groundedGradientSP = m_mesh_3.HS3F()*x;

	groundedGradient = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*x;

	//computeScalarPotential(m_mesh_3.HS3F()*x);
}

void Decomposition::computeCurlyGradient()
{
	curlyGradient = omega - harmonicKnot - harmonicGradient - fluxlessKnot - groundedGradient;

	int hm1 = m_mesh_3.NumHMLG1();

	// vector potential
	if (hm1 == 0) {
		EigenVector tmpCG = curlyGradient;

		EigenSpMat s1l1 = m_mesh_3.HS1F_B()*m_laplacian1_B;
		EigenVector b = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*tmpCG;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_vp;
		solver_vp.compute(s1l1);
		EigenVector x = solver_vp.solve(b);

		//EigenVector x = matlabLinearSolveCholesky(s1l1, b);

		curlyGradientVP = x;
	}
	else {
		EigenVector tmpCG = curlyGradient;

		EigenSpMat s1l1 = m_mesh_3.HS1F_B()*m_laplacian1_B;
		EigenVector b = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*tmpCG;

		std::vector<double> row_vp, col_vp;
		std::vector<double> val_vp;
		for (int i = 0; i < s1l1.outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
				row_vp.push_back(static_cast<double>(iter.row()) + 1);
				col_vp.push_back(static_cast<double>(iter.col()) + 1);
				val_vp.push_back(static_cast<double>(iter.value()));
			}
		}
		//std::vector<double> hs_vp;
		//for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		//	hs_vp.push_back(m_mesh_3.HS1F_B().coeff(i, i));
		//}

		//std::vector<double> emptyEval_vp;
		//std::vector<EigenVector> eigenFields_vp;
		//matlabEIGS(emptyEval_vp, eigenFields_vp, row_vp, col_vp, val_vp, hs_vp, hm1);

		std::vector<double> hrow, hcol;
		std::vector<double> hval;
		for (int i = 0; i < m_mesh_3.HS1F_B().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_B(), i); iter; ++iter) {
				hrow.push_back(static_cast<double>(iter.row()) + 1);
				hcol.push_back(static_cast<double>(iter.col()) + 1);
				hval.push_back(static_cast<double>(iter.value()));
			}
		}


		std::vector<double> emptyEval_vp;
		std::vector<EigenVector> eigenFields_vp;
		matlabEIGS(emptyEval_vp, eigenFields_vp, row_vp, col_vp, val_vp, hrow, hcol, hval, hm1);

		std::set<int> anchors_vp;
		selectAnchors(eigenFields_vp, anchors_vp);
		correctRankDeficiency(s1l1, b, anchors_vp, m_mesh_3.HS1F_B(), eigenFields_vp);

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_vp;
		solver_vp.compute(s1l1);
		curlyGradientVP = solver_vp.solve(b);

		//curlyGradientVP = matlabLinearSolveCholesky(s1l1, b);

		cout << "Potential curly gradient projection coeffs: " << endl;
		for (int i = 0; i < eigenFields_vp.size(); ++i) {
			double c = curlyGradientVP.dot(m_mesh_3.HS1F_B()*eigenFields_vp[i])
				/ eigenFields_vp[i].dot(m_mesh_3.HS1F_B()*eigenFields_vp[i]);
			curlyGradientVP -= c * eigenFields_vp[i];
			cout << c << endl;
		}
	}

	/**************************************/
	// scalar potential
	EigenVector tmpCG = omega - harmonicKnot - harmonicGradient - fluxlessKnot - groundedGradient;
	EigenVector tmpCG_NB = tmpCG;
	boundaryToNoBoundary2Form(tmpCG_NB);
	EigenSpMat s3l3 = m_laplacian3_NB;
	EigenVector b =  m_mesh_3.ED2F_NB()*tmpCG_NB;

	
	for (int i = 0; i < s3l3.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s3l3, i); iter; ++iter) {
			if (iter.row() == 0 && iter.col() == 0) {
				iter.valueRef() = 1;
			}
			else if (iter.row() == 0 || iter.col() == 0) {
				iter.valueRef() = 0;
			}
		}
	}
	b(0) = 0;

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(s3l3);
	EigenVector x = solver.solve(b);

	//EigenVector x = matlabLinearSolveCholesky(s3l3, b);

	//EigenVector cg = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*x; // !!!!!! alternative

	EigenVector cg = m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*x;
	noBoundaryToBoundary2Form(cg);
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;

		if (m_mesh_3.FacetOnBoundary(fi)) {
			cg(m_mesh_3.FacetIdx(fi)) = tmpCG(m_mesh_3.FacetIdx(fi));
		}
	}

	EigenVector diff = cg - curlyGradient;
	cout << "WTF: " << diff.dot(m_mesh_3.ED2F_B()*diff)/ curlyGradient.dot(m_mesh_3.ED2F_B()*curlyGradient) << endl;

	//computeScalarPotential(m_mesh_3.HS3F()*x);
	curlyGradientSP = m_mesh_3.HS3F()*x;
}

void  Decomposition::computeCurlyGradientSP()
{

}

void  Decomposition::computeCurlyGradientVP()
{

}

void  Decomposition::computeHarmonicKnotVP()
{

}

void  Decomposition::computeNewtonianDivField()
{
	EigenSpMat A = m_mesh_3.HS3F() * m_laplacian3_NB;

	EigenVector noHarmonicComponentNB = omega;
	boundaryToNoBoundary2Form(noHarmonicComponentNB);

	EigenVector b = m_mesh_3.HS3F() * m_mesh_3.ED2F_NB() * noHarmonicComponentNB;

	for (int i = 0; i < A.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
			if (iter.row() == 0 && iter.col() == 0) {
				iter.valueRef() = 1;
			}
			else if (iter.row() == 0 || iter.col() == 0) {
				iter.valueRef() = 0;
			}
		}
	}

	b(0) = 0; // set fix potential root.

	std::cout << "A and b are set." << std::endl;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	EigenVector x = solver.solve(b); // potential 3 form for coexact component

	groundedGradientSP = m_mesh_3.HS3F()*x;

	groundedGradient = m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*x;
	noBoundaryToBoundary2Form(groundedGradient);

}

EigenVector Decomposition::matlabLinearSolveCholesky(EigenSpMat A, EigenVector b)
{
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < A.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> bv;
	for (int i = 0; i < b.size(); ++i) {
		bv.push_back(b(i));
	}

	Engine* eng;
	eng = engOpen("\0");

	std::string str;

	// transfer data
	mxArray *rowArr, *colArr, *valArr, *bvArr;
	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	bvArr = mxCreateDoubleMatrix(bv.size(), 1, mxREAL);

	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	memcpy(mxGetPr(bvArr), &bv[0], bv.size() * sizeof(double));

	engPutVariable(eng, "row", rowArr);
	engPutVariable(eng, "col", colArr);
	engPutVariable(eng, "val", valArr);
	engPutVariable(eng, "b", bvArr);

	const int buffersize = 256;
	char buffer[buffersize + 1];
	buffer[buffersize] = '\0';
	engOutputBuffer(eng, buffer, buffersize);
	
	str = "A=sparse(row, col, val);";       engEvalString(eng, str.c_str());
	str = "tic;";							engEvalString(eng, str.c_str());
	str = "x=A\\b;";						engEvalString(eng, str.c_str());
	
	//str = "R=chol(A);";						engEvalString(eng, str.c_str());
	//str = "x=R\\(R'\\b);";					engEvalString(eng, str.c_str());
	str = "elapsed=toc";					engEvalString(eng, str.c_str());
	printf("%s", buffer);

	mxArray* xx;
	xx = engGetVariable(eng, "x");
	double* x;
	x = reinterpret_cast<double*>(mxGetData(xx));
	const size_t* dims;
	dims = mxGetDimensions(xx);

	EigenVector solution;
	solution.resize(dims[0]);

	cout << dims[0] << " " << dims[1] << endl;

	for (int i = 0; i < dims[0]; ++i) {
		solution(i) = x[i];
	}

	return solution;
}

void Decomposition::matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
	std::vector<double> row, std::vector<double> col, std::vector<double> val,
	std::vector<double> brow, std::vector<double> bcol, std::vector<double> bval,
	double hm)
{
	Engine* eng;
	eng = engOpen("\0");

	std::string str;

	// transfer data
	mxArray *rowArr, *colArr, *valArr, *browArr, *bcolArr, *bvalArr, *hmArr;
	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	browArr = mxCreateDoubleMatrix(brow.size(), 1, mxREAL);
	bcolArr = mxCreateDoubleMatrix(bcol.size(), 1, mxREAL);
	bvalArr = mxCreateDoubleMatrix(bval.size(), 1, mxREAL);
	hmArr = mxCreateDoubleMatrix(1, 1, mxREAL);
	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	memcpy(mxGetPr(browArr), &brow[0], brow.size() * sizeof(double));
	memcpy(mxGetPr(bcolArr), &bcol[0], bcol.size() * sizeof(double));
	memcpy(mxGetPr(bvalArr), &bval[0], bval.size() * sizeof(double));
	memcpy(mxGetPr(hmArr), &hm, sizeof(double));
	engPutVariable(eng, "row", rowArr);
	engPutVariable(eng, "col", colArr);
	engPutVariable(eng, "val", valArr);
	engPutVariable(eng, "brow", browArr);
	engPutVariable(eng, "bcol", bcolArr);
	engPutVariable(eng, "bval", bvalArr);
	engPutVariable(eng, "hm", hmArr);

	const int buffersize = 256;
	char buffer[buffersize + 1];
	buffer[buffersize] = '\0';
	engOutputBuffer(eng, buffer, buffersize);

	str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());
	str = "B=sparse(brow, bcol, bval)";        engEvalString(eng, str.c_str());
	str = "tic";						engEvalString(eng, str.c_str());
	str = "[V,D]=eigs(A, B, hm, 'smallestabs')";        engEvalString(eng, str.c_str());
	str = "elapsed=toc";					engEvalString(eng, str.c_str());
	printf("%s", buffer);

	// read out eigenvector matrix D from matlab
	mxArray* vv;
	vv = engGetVariable(eng, "V");
	double* v;
	v = reinterpret_cast<double*>(mxGetData(vv));
	const size_t* dims;
	dims = mxGetDimensions(vv);

	EigenVector ev; // no boundary
	ev.resize(dims[0]);

	evec.clear();
	for (int j = 0; j < dims[1]; ++j) {
		for (int i = 0; i < dims[0]; ++i) {
			ev(i) = v[i + j * dims[0]];
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

void Decomposition::integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf,
	std::map<Mesh_3_Vertex_iterator, Mesh_3_Vector_3>& cf,
	double dsty, double lth, int stps, int md, std::string sfx)
{
	RKIntegrator integrator(m_mesh_3);

	integrator.init(dsty, lth, stps, md);
	integrator.initColorField(cf);
	integrator.convert(vf);
	integrator.integrate();
	integrator.write(sfx);

}

void Decomposition::convertForms()
{
	convert2Form(omega, omegaVF);
	if (m_mesh_3.NumHMLG1() != 0)
		convert2Form(harmonicKnot, harmonicKnotVF);
	if (m_mesh_3.NumHMLG2() != 0)
		convert2Form(harmonicGradient, harmonicGradientVF);
	convert2Form(fluxlessKnot, fluxlessKnotVF);
	convert2Form(groundedGradient, groundedGradientVF);
	convert2Form(curlyGradient, curlyGradientVF);
	cout << "Field form conversion done ..." << endl;

	//if (m_mesh_3.NumHMLG1() != 0)
	//	convert1Form(harmonicKnotVP, harmonicKnotVPVF);
	//convert1Form(fluxlessKnotVP, fluxlessKnotVPVF);
	//convert1Form(curlyGradientVP, curlyGradientVPVF);
	//cout << "Vector potential conversion done ..." << endl;

	//if (m_mesh_3.NumHMLG2() != 0) {
	//	convert3Form(harmonicGradientSP, harmonicGradientSPSF);
	//	correctBoundaryPotential(harmonicGradientSPSF);
	//}
	//convert3Form(groundedGradientSP, groundedGradientSPSF);
	//correctBoundaryPotential(groundedGradientSPSF);
	//convert3Form(curlyGradientSP, curlyGradientSPSF);
	//cout << "Scalar potential conversion done ..." << endl;
}

void Decomposition::computeArrows()
{
	assignColorStrength(omegaVF, 0);
	if (m_mesh_3.NumHMLG1() != 0)
		assignColorStrength(harmonicKnotVF, 1);
	if (m_mesh_3.NumHMLG2() != 0)
		assignColorStrength(harmonicGradientVF, 2);
	assignColorStrength(fluxlessKnotVF, 3);
	assignColorStrength(groundedGradientVF, 4);
	assignColorStrength(curlyGradientVF, 5);

	//if (m_mesh_3.NumHMLG1() != 0)
	//	assignColorStrengthVP(harmonicKnotVPVF, 0);
	//assignColorStrengthVP(fluxlessKnotVPVF, 1);
	//assignColorStrengthVP(curlyGradientVPVF, 2);

	//if (m_mesh_3.NumHMLG2() != 0)
	//	assignColorStrengthSP(harmonicGradientSPSF, 0);
	//assignColorStrengthSP(groundedGradientSPSF, 1);
	//assignColorStrengthSP(curlyGradientSPSF, 2);

	rescaleArrowsWithVF(omegaVF);
	if (m_mesh_3.NumHMLG1() != 0)
		rescaleArrowsWithVF(harmonicKnotVF);
	if (m_mesh_3.NumHMLG2() != 0)
		rescaleArrowsWithVF(harmonicGradientVF);
	rescaleArrowsWithVF(fluxlessKnotVF);
	rescaleArrowsWithVF(groundedGradientVF);
	rescaleArrowsWithVF(curlyGradientVF);

	//if (m_mesh_3.NumHMLG1() != 0)
	//	rescaleArrowsWithVF(harmonicKnotVPVF);
	//rescaleArrowsWithVF(fluxlessKnotVPVF);
	//rescaleArrowsWithVF(curlyGradientVPVF);

	computeArrowsWithVF(groundedGradientVF);
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

							colors.push_back(Mesh_3_Vector_3(0.6, 0.6, 0.6));
							//colors.push_back(colorField[fi->first->vertex(i)]);
							/*if (mode == 0)
								colors.push_back(harmonicGradientSPCF[fi->first->vertex(i)]);
							else if (mode == 1)
								colors.push_back(groundedGradientSPCF[fi->first->vertex(i)]);
							else
								colors.push_back(curlyGradientSPCF[fi->first->vertex(i)]);*/

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

							colors.push_back(Mesh_3_Vector_3(0.6, 0.6, 0.6));
							//colors.push_back(colorField[fi->first->vertex(i)]);
							/*if (mode == 0)
								colors.push_back(harmonicGradientSPCF[fi->first->vertex(i)]);
							else if (mode == 1)
								colors.push_back(groundedGradientSPCF[fi->first->vertex(i)]);
							else
								colors.push_back(curlyGradientSPCF[fi->first->vertex(i)]);*/

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

void Decomposition::writeMesh()
{
	std::ofstream out("mesh.obj");

	for (int i = 0; i < vertices.size(); ++i) {
		out << "v "
			<< (vertices[i].x() /*- m_mesh_3.getCenter().x()*/) / 1/*m_mesh_3.getRadius()*/ << " "
			<< (vertices[i].y() /*- m_mesh_3.getCenter().y()*/) / 1/*m_mesh_3.getRadius()*/ << " "
			<< (vertices[i].z() /*- m_mesh_3.getCenter().z()*/) / 1/*m_mesh_3.getRadius()*/ << endl;
	}

	for (int i = 0; i < indices.size(); i = i + 3) {
		out << "f "
			<< indices[i] + 1 << " "
			<< indices[i + 1] + 1 << " "
			<< indices[i + 2] + 1 << endl;
	}
}

void Decomposition::writeMeshPLY(std::string sfx)
{
	std::ofstream out("mesh_" + sfx + ".ply");

	out << "ply" << std::endl;
	out << "format ascii 1.0" << std::endl;
	out << "element vertex " << vertices.size() << std::endl;
	out << "property float32 x" << std::endl;
	out << "property float32 y" << std::endl;
	out << "property float32 z" << std::endl;
	out << "property uchar red" << std::endl;
	out << "property uchar green" << std::endl;
	out << "property uchar blue" << std::endl;
	out << "element face " << indices.size() / 3 << std::endl;
	out << "property list uchar int vertex_index" << std::endl;
	out << "end_header" << std::endl;

	for (int i = 0; i < vertices.size(); ++i) {
		out << vertices[i].x() << " "
			<< vertices[i].y() << " "
			<< vertices[i].z() << " "
			<< floor(colors[i].x() * 255) << " "
			<< floor(colors[i].y() * 255) << " "
			<< floor(colors[i].z() * 255) << " "
			<< std::endl;
	}
	for (int i = 0; i < indices.size(); i = i + 3) {
		out << "3 "
			<< indices[i] << " "
			<< indices[i + 1] << " "
			<< indices[i + 2] << std::endl;
	}
}

void Decomposition::computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	vfvNum = 0;

	int i = 0;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() != 0) {
			if (i % 4 == 0) {

				bool boundary = false;
				for (int k = 0; k < 4; ++k) {
					My_facet mf(ci, k);
					Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);
					if (m_mesh_3.FacetOnBoundary(fi))
						boundary = true;
				}

				if (!boundary)
					assembleArrow(m_mesh_3.CellCC(ci) - CGAL::ORIGIN, (m_mesh_3.CellCC(ci) - CGAL::ORIGIN) + vf[ci], 0.04*m_mesh_3.getAveLen(), ci);

			}
			++i;
		}
	}
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
	pivotLen = lenVec[static_cast<int>(floor(lenVec.size()*0.97))];


	double aveLen = m_mesh_3.getAveLen();
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Vector_3 v = vf[ci];
		double len = sqrt(v*v);
		v /= len;
		//len *= 10;
		//if (len > 0.7*m_aveLen) 
			//len = 0.7*m_aveLen; 
		if (len > pivotLen)
			len = pivotLen;
		if (len < 0.3*pivotLen)
			len = 0.3*pivotLen;

		len = (0.5*aveLen)*(len / (pivotLen));

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
		d /= sqrt(d*d);

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

		if (ori)
			field(m_mesh_3.FacetIdx(fi)) = av * d;
		else
			field(m_mesh_3.FacetIdx(fi)) = -av * d;
	}
}

void Decomposition::setCurrentMagneticField(Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, EigenVector& field)
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
	
		field(m_mesh_3.FacetIdx(fi)) = av * fD;
	}
}

void Decomposition::setConstantField(Mesh_3_Vector_3 d, EigenVector& field)
{
	field.setZero(m_mesh_3.NumF());
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (m_mesh_3.FacetValid(fi))
			continue;
		
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

		field(m_mesh_3.FacetIdx(fi)) = av * d;
	}
}

void Decomposition::smooth(EigenVector& form)
{
	EigenVector formNB = form;
	boundaryToNoBoundary2Form(formNB);

	EigenSpMat s2l2 = m_mesh_3.HS2F_NB() * m_laplacian2_NB;

	// identity
	//EigenSpMat id;
	//id.resize(m_mesh_3.NumIF(), m_mesh_3.NumIF());
	//id.setIdentity();

	double delta = 0.01; // smoothing time step

	EigenVector current = formNB; // form in current time step;
	EigenSpMat op = m_mesh_3.HS2F_NB() + delta * s2l2; // this is an symmetric matrix
	EigenVector next; // form after 1 step of smoothing
	Eigen::SimplicialLDLT<EigenSpMat> solver;
	solver.compute(op);
	int numIter = 5;
	int i = 0;
	while(i<numIter){
		current = m_mesh_3.HS2F_NB() * current;
		next = solver.solve(current);
		current = next;

		++i;
	}

	form = current;
	noBoundaryToBoundary2Form(form);
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

void Decomposition::assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci)
{
	Mesh_3_Vector_3 d = head - tail;
	if (sqrt(d*d) > 0.0000000001)
	{
		vec3 h(head.x(), head.y(), head.z());
		vec3 t(tail.x(), tail.y(), tail.z());
		Arrow ar(radius, t, h);

		std::vector<float> vertices = ar.getVertices();
		std::vector<int> faces = ar.getFaces();

		typedef CGAL::Simple_cartesian<double> Kernel;
		typedef Enriched_polyhedron<Kernel, Enriched_items> Mesh;
		typedef Mesh::HalfedgeDS HalfedgeDS;
		typedef Mesh::Vertex_iterator Vertex_iterator;
		typedef Mesh::Facet_iterator Facet_iterator;

		Mesh arw;
		ObjBuilder<HalfedgeDS> builder(vertices, faces);
		arw.delegate(builder);
		arw.basic_init();

		for (Vertex_iterator vi = arw.vertices_begin(); vi != arw.vertices_end(); ++vi) {
			vfVertices.push_back(Mesh_3_Vector_3(vi->point().x(), vi->point().y(), vi->point().z()));
			vfvColors.push_back(Mesh_3_Vector_3(0.7, 0.4, 0.1));
		}
		for (Facet_iterator fi = arw.facets_begin(); fi != arw.facets_end(); ++fi) {
			vfFaces.push_back(Mesh_3_Vector_3(fi->halfedge()->vertex()->idx() + vfvNum,
				fi->halfedge()->next()->vertex()->idx() + vfvNum,
				fi->halfedge()->next()->next()->vertex()->idx() + vfvNum));
			vfNormals.push_back(Mesh_3_Vector_3(fi->normal().x(), fi->normal().y(), fi->normal().z()));
			vfColors.push_back(Mesh_3_Vector_3(0.7, 0.4, 0.1));
		}

		vfvNum += static_cast<int>(arw.size_of_vertices());
	}
}

void Decomposition::selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors)
{
	cout << eigenfields.size() << endl;

	std::random_device rd;
	std::mt19937_64 mt(rd());
	std::uniform_int_distribution<int> distribution(0, m_mesh_3.NumIE());

	double condNumber;

	int size = static_cast<int>(eigenfields.size());
	int iteration = 1;
	int kk = 0;
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

		if (fabs(mat.determinant()) < 10e-15)
			continue;

		{
			Eigen::JacobiSVD<EigenMatrix> svd(mat);
			double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

			cout << "Cond: " << fabs(cond) << endl;
			cout << "Det: " << fabs(mat.determinant()) << endl;

			if (kk == 0){
				condNumber = fabs(cond);
				anchors = rows;

				m_cond = condNumber;
			}

			if (fabs(cond) > condNumber) {
				condNumber = fabs(cond);
				anchors = rows;

				m_cond = condNumber;
			}

			++kk;
		}
	}

	for (auto iter = anchors.begin(); iter != anchors.end(); ++iter) {
		cout<<*iter<<endl;
	}
}

void Decomposition::correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors, EigenSpMat hs, std::vector<EigenVector> eigenfields)
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

	int row = eigenfields[0].size();
	int col = anchors.size();
	EigenSpMat H(row, col);
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			H.coeffRef(i, j) = eigenfields[j](i);
		}
	}

	b -= hs * H*H.transpose()*b;

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

void Decomposition::assignColorStrengthVP(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vp, int mode) {
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
		colorField = curlyGradientVPCF;
	}
}

void Decomposition::assignColorStrengthSP(std::map<Mesh_3_Vertex_iterator, double> sf, int mode) {
	//std::map<Mesh_3_Vertex_iterator, double> vertexStrength;
	//std::vector<double> len;
	//for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
	//	vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
	//	if (!m_mesh_3.VertexValid(vi))
	//		continue;

	//	std::vector<Mesh_3_Cell_iterator> cells;
	//	m_mesh_3.C3T3().triangulation().incident_cells(vi, std::back_inserter(cells));

	//	double p = 0;
	//	double vol = 0;
	//	for (int i = 0; i < cells.size(); ++i) {
	//		if (cells[i]->subdomain_index() == 0)
	//			continue;

	//		p += form(m_mesh_3.CellIdx(cells[i]))* m_mesh_3.CellPrimal(cells[i]);
	//		vol += m_mesh_3.CellPrimal(cells[i]);
	//	}
	//	p /= vol;

	//	vertexStrength[vi] = p;
	//	len.push_back(p);
	//}

	//std::vector<int> sIndices;
	//sortedIndices(sIndices, len);

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

		double d = sf[vi] - (minv + 0.2*diff);

		double ratio = d / (diff*0.6);

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


/************************************************************\
 *															* 
 *	Debug functions.										*
 *															*
\************************************************************/

void Decomposition::testSurface()
{
	Mesh_3_Point_3 center = m_mesh_3.getCenter();
	double radius = m_mesh_3.getRadius();

	std::cout << "output" << std::endl;
	std::ofstream out("CrossSection.obj");

	for (int i = 0; i < vertices.size(); ++i) {
		out << "v "
			<< (vertices[i].x() - center.x()) / radius << " "
			<< (vertices[i].y() - center.y()) / radius << " "
			<< (vertices[i].z() - center.z()) / radius << std::endl;
	}
	for (int i = 0; i < indices.size(); i = i + 3) {
		out << "f " << indices[i] + 1 << " " << indices[i + 1] + 1 << " " << indices[i + 2] + 1 << std::endl;
	}
}

void Decomposition::testSurfaceAndArrowsPLY()
{
	Mesh_3_Point_3 center = m_mesh_3.getCenter();
	double radius = m_mesh_3.getRadius();

	Mesh_3_Vector_3 translation = -(center - CGAL::ORIGIN);
	double scale = 1 / radius;

	// model
	std::ofstream out("CrossSection.ply");
	// formating ply files
	out << "ply" << std::endl;
	out << "format ascii 1.0" << std::endl;
	out << "element vertex " << vertices.size() << std::endl;
	out << "property float32 x" << std::endl;
	out << "property float32 y" << std::endl;
	out << "property float32 z" << std::endl;
	out << "property uchar red" << std::endl;
	out << "property uchar green" << std::endl;
	out << "property uchar blue" << std::endl;
	out << "element face " << indices.size() / 3 << std::endl;
	out << "property list uchar int vertex_index" << std::endl;
	out << "end_header" << std::endl;
	// data
	for (int i = 0; i < vertices.size(); ++i) {
		out << (vertices[i].x() + translation.x())*scale << " "
			<< (vertices[i].y() + translation.y())*scale << " "
			<< (vertices[i].z() + translation.z())*scale << " "
			<< floor(colors[i].x() * 255) << " "
			<< floor(colors[i].y() * 255) << " "
			<< floor(colors[i].z() * 255) << " "
			<< std::endl;
	}
	for (int i = 0; i < indices.size(); i = i + 3) {
		out << "3 "
			<< indices[i] << " "
			<< indices[i + 1] << " "
			<< indices[i + 2] << std::endl;
	}

	//field
	std::ofstream out_field("Field.ply");
	// formating ply files
	out_field << "ply" << std::endl;
	out_field << "format ascii 1.0" << std::endl;
	out_field << "element vertex " << vfVertices.size() << std::endl;
	out_field << "property float32 x" << std::endl;
	out_field << "property float32 y" << std::endl;
	out_field << "property float32 z" << std::endl;
	out_field << "property uchar red" << std::endl;
	out_field << "property uchar green" << std::endl;
	out_field << "property uchar blue" << std::endl;
	out_field << "element face " << vfFaces.size() << std::endl;
	out_field << "property list uchar int vertex_index" << std::endl;
	out_field << "end_header" << std::endl;
	// data
	for (int i = 0; i < vfVertices.size(); ++i) {
		out_field << (vfVertices[i].x() + translation.x())*scale << " "
			<< (vfVertices[i].y() + translation.y())*scale << " "
			<< (vfVertices[i].z() + translation.z())*scale << " "
			<< floor(vfvColors[i].x() * 255) << " "
			<< floor(vfvColors[i].y() * 255) << " "
			<< floor(vfvColors[i].z() * 255) << " "
			<< std::endl;
	}
	for (int i = 0; i < vfFaces.size(); ++i) {
		out_field << "3 "
			<< vfFaces[i].x() << " "
			<< vfFaces[i].y() << " "
			<< vfFaces[i].z() << " "
			<< std::endl;
	}
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

	//color = (1 - ratio)*blue + ratio * pink;

	if (ratio < 0.5) {
		double r = ratio / 0.5;
		color = (1 - r)*blue + r * white;
	}
	else {
		double r = (ratio - 0.5) / 0.5;
		color = (1 - r)*white + r * pink;
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

	if (ratio < 0.5) {
		double r = ratio / 0.5;
		color = (1 - r)*purple + r * white;
	}
	else {
		double r = (ratio - 0.5) / 0.5;
		color = (1 - r)*white + r * orange;
	}

	return color;
}

void Decomposition::sortedIndices(std::vector<int>& sIndex, std::vector<double> len)
{
	sIndex.resize(len.size(), 0);
	std::iota(sIndex.begin(), sIndex.end(), 0);

	std::sort(sIndex.begin(), sIndex.end(), [&len](int i1, int i2) {return len[i1] < len[i2]; });
}

void Decomposition::evaluateAccuracy()
{
	ofstream out("stat.txt", std::ios::app);

	out << m_accu << " " << m_cond << endl;
}

void Decomposition::testOrthogonality()
{
	cout<<"Orthogonality: "<<endl;
	cout << harmonicKnot.dot(m_mesh_3.HS2F_B()*harmonicGradient) << endl;
	cout << harmonicKnot.dot(m_mesh_3.HS2F_B()*fluxlessKnot) << endl;
	cout << harmonicKnot.dot(m_mesh_3.HS2F_B()*groundedGradient) << endl;
	cout << harmonicKnot.dot(m_mesh_3.HS2F_B()*curlyGradient) << endl;

	cout << harmonicGradient.dot(m_mesh_3.HS2F_B()*fluxlessKnot) << endl;
	cout << harmonicGradient.dot(m_mesh_3.HS2F_B()*groundedGradient) << endl;
	cout << harmonicGradient.dot(m_mesh_3.HS2F_B()*curlyGradient) << endl;

	cout << fluxlessKnot.dot(m_mesh_3.HS2F_B()*groundedGradient) << endl;
	cout << fluxlessKnot.dot(m_mesh_3.HS2F_B()*curlyGradient) << endl;

	cout << groundedGradient.dot(m_mesh_3.HS2F_B()*curlyGradient) << endl;

	cout<<"Norm: "<<endl;
	cout << harmonicKnot.dot(m_mesh_3.HS2F_B()*harmonicKnot) << endl;
	cout << harmonicGradient.dot(m_mesh_3.HS2F_B()*harmonicGradient) << endl;
	cout << fluxlessKnot.dot(m_mesh_3.HS2F_B()*fluxlessKnot) << endl;
	cout << groundedGradient.dot(m_mesh_3.HS2F_B()*groundedGradient) << endl;
	cout << curlyGradient.dot(m_mesh_3.HS2F_B()*curlyGradient) << endl;
}


void Decomposition::groupFields(std::vector<EigenVector> eigenFields)
{
	int curlIdx = 0;
	int divIdx = 0;

	for (int i = 0; i < eigenFields.size(); ++i) {
		cout<<eigenFields[i].dot(m_curlEM*eigenFields[i])<<" "<<eigenFields[i].dot(m_divEM*eigenFields[i])<<endl;
		if (eigenFields[i].dot(m_curlEM*eigenFields[i]) > eigenFields[i].dot(m_divEM*eigenFields[i])) {
			curlFieldIdx[divIdx] = i;
			++divIdx;
		}
		else {
			divFieldIdx[curlIdx] = i;
			++curlIdx;
		}
	}

	cout << "Curl Fields: " << curlIdx << endl;
	cout << "Div Fields: " << divIdx << endl;
}

void Decomposition::correctBoundaryPotential(std::map<Mesh_3_Vertex_iterator, double>& sf)
{
	std::vector<double> bvt; // boundary vertex potential group
	std::vector<int> count;
	bvt.resize(m_mesh_3.NumBoundary(), 0);
	count.resize(m_mesh_3.NumBoundary(), 0);

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		bvt[m_mesh_3.VertexBoundaryIdx(vi)] += sf[vi];
		count[m_mesh_3.VertexBoundaryIdx(vi)]++;
	}

	for (int i = 0; i < bvt.size(); ++i) {
		bvt[i] /= count[i];
	}

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		int bidx = m_mesh_3.VertexBoundaryIdx(vi);

		sf[vi] = bvt[bidx];
	}
}

void Decomposition::createDivTan(EigenVector& form)
{
	// assign a scalar function --------------------------------------------
	Mesh_3_Point_3 s1(0, 0.05, -0.3);
	Mesh_3_Point_3 s2(0, 0.05, 0.3);

	Mesh_3_Vector_3 v1 = m_mesh_3.CellCC(m_mesh_3.C3T3().triangulation().cells_begin()) - s1;
	Mesh_3_Vector_3 v2 = m_mesh_3.CellCC(m_mesh_3.C3T3().triangulation().cells_begin()) - s2;

	double d1 = sqrt(v1*v1);
	double d2 = sqrt(v2*v2);

	int vidx1 = 0;
	int vidx2 = 0;

	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;
	
		Mesh_3_Vector_3 vv1 = m_mesh_3.CellCC(ci) - s1;
		Mesh_3_Vector_3 vv2 = m_mesh_3.CellCC(ci) - s2;

		if (sqrt(vv1*vv1) < d1) {
			d1 = sqrt(vv1*vv1);
			vidx1 = m_mesh_3.CellIdx(ci);
		}

		if (sqrt(vv2*vv2) < d2) {
			d2 = sqrt(vv2*vv2);
			vidx2 = m_mesh_3.CellIdx(ci);
		}
	}

	// solve for grounded gradient --------------------------------------------
	EigenVector b;
	b.resize(m_mesh_3.NumC());
	b.setZero();
	b(vidx1) = 1;
	b(vidx2) = -1;

	EigenSpMat A_gg = m_laplacian3_B;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_gg;
	solver_gg.compute(A_gg);
	EigenVector x_gg = solver_gg.solve(b); // potential 3 form

	EigenVector ggsf = m_mesh_3.HS3F()*x_gg;

	EigenVector gg = m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*x_gg;

	// solve for eta ------------------------------------------------------
	EigenVector cg;
	//cg.resize(m_mesh_3.NumF());
	//cg.setZero();

	EigenVector tmp_gg = gg;
	boundaryToNoBoundary2Form(tmp_gg);
	EigenSpMat s3l3_eta = m_laplacian3_NB;
	EigenVector b_eta;
	//EigenVector b_eta = m_mesh_3.ED2F_NB()*tmp_gg;

	//for (int i = 0; i < s3l3_eta.outerSize(); ++i) {
	//	for (EigenSpMat::InnerIterator iter(s3l3_eta, i); iter; ++iter) {
	//		if (iter.row() == 0 && iter.col() == 0) {
	//			iter.valueRef() = 1;
	//		}
	//		else if (iter.row() == 0 || iter.col() == 0) {
	//			iter.valueRef() = 0;
	//		}
	//	}
	//}
	//b_eta(0) = 0;

	b_eta.setZero(m_mesh_3.NumC());
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		for (int j = 0; j < 4; ++j) {
			My_facet mf(ci, j);
			Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);

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

			if (m_mesh_3.FacetOnBoundary(fi)) {
				b_eta(m_mesh_3.CellIdx(ci)) += even * gg(m_mesh_3.FacetIdx(fi));
			}
		}
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_eta;
	solver_eta.compute(s3l3_eta);
	EigenVector x_eta = solver_eta.solve(b_eta);

	cg = m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*x_eta;
	noBoundaryToBoundary2Form(cg);
	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;

		if (m_mesh_3.FacetOnBoundary(fi)) {
			cg(m_mesh_3.FacetIdx(fi)) = gg(m_mesh_3.FacetIdx(fi));
		}
	}

	// get final tangential divergence field

	EigenVector divtan = gg + cg;
	//groundedGradient = divtan;

	form = divtan;
}

Mesh_3_Vector_3 Decomposition::normalGradient(Mesh_3_Vector_3 p)
{
	return p*1;
}

Mesh_3_Vector_3 Decomposition::tangentialCurl(Mesh_3_Vector_3 p)
{
	return Mesh_3_Vector_3(p.y(), -p.x(), 0) * 1;
}

Mesh_3_Vector_3 Decomposition::normalHarmonic(Mesh_3_Vector_3 p)
{
	double de = pow((p*p), -3.0 / 2.0);

	return de * p;
}

Mesh_3_Vector_3 Decomposition::tangentialHarmonic(Mesh_3_Vector_3 p)
{
	double de = 1 / (p.x()*p.x() + p.y()*p.y());

	return de * Mesh_3_Vector_3(p.y(), -p.x(), 0);
}

Mesh_3_Vector_3 Decomposition::centralHarmonic(Mesh_3_Vector_3 p)
{
	return Mesh_3_Vector_3(1, 1, 1) *0.5;
}

double Decomposition::fluxNormalGradient(Mesh_3_Facet_iterator fi)
{
	// get normal vector
	Mesh_3_Vector_3 n = m_mesh_3.FacetNormal(fi);

	std::vector<int> vIdx;
	for (int i = 0; i < 4; ++i) {
		if (i != fi->second) {
			vIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i)));
		}
	}
	std::sort(vIdx.begin(), vIdx.end());
	
	Mesh_3_Vector_3 p = m_mesh_3.IdxToVertex(vIdx[0])->point().point() - CGAL::ORIGIN;
	Mesh_3_Vector_3 u = m_mesh_3.IdxToVertex(vIdx[1])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();
	Mesh_3_Vector_3 v = m_mesh_3.IdxToVertex(vIdx[2])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();

	double integral = (1.0 / 6.0)*(n.x()*(3 * p.x() + u.x() + v.x()) + n.y()*(3 * p.y() + u.y() + v.y()) + n.z()*(3 * p.z() + u.z() + v.z()));
	
	return integral;
}

double Decomposition::fluxTangentialCurl(Mesh_3_Facet_iterator fi)
{
	// get normal vector
	Mesh_3_Vector_3 n = m_mesh_3.FacetNormal(fi);

	std::vector<int> vIdx;
	for (int i = 0; i < 4; ++i) {
		if (i != fi->second) {
			vIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i)));
		}
	}
	std::sort(vIdx.begin(), vIdx.end());

	Mesh_3_Vector_3 p = m_mesh_3.IdxToVertex(vIdx[0])->point().point() - CGAL::ORIGIN;
	Mesh_3_Vector_3 u = m_mesh_3.IdxToVertex(vIdx[1])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();
	Mesh_3_Vector_3 v = m_mesh_3.IdxToVertex(vIdx[2])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();

	double integral = (1.0 / 6.0)*(-n.y()*(3*p.x()+u.x()+v.x())+ n.x()*(3 * p.y() + u.y() + v.y()));

	return integral;
}

double Decomposition::fluxNormalHarmonic(Mesh_3_Facet_iterator fi)
{
	return 0;
}

double Decomposition::fluxTangentialHarmonic(Mesh_3_Facet_iterator fi)
{
	return 0;
}

double Decomposition::fluxCentralHarmonic(Mesh_3_Facet_iterator fi)
{
	Mesh_3_Vector_3 n = m_mesh_3.FacetNormal(fi);

	return (n.x() + n.y() + n.z()) / 2;
}

void Decomposition::sampleInput()
{
	omega.resize(m_mesh_3.NumF());
	//omega = sampleField(0);
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);
		omega(i) = 0;
		//omega(i) += sampleFacetFlux(fi, 0, 3);
		//omega(i) += sampleFacetFlux(fi, 1, 3);
		//omega(i) += sampleFacetFlux(fi, 2, 3);
		omega(i) += sampleFacetFlux(fi, 3, 3);
		//omega(i) += sampleFacetFlux(fi, 4, 3);
		//omega(i) += sampleFacetFlux(fi, 5, 3);
	}

	omegaOriginal = omega;
}

void Decomposition::computeError()
{
	EigenVector ggana = sampleField(1);
	EigenVector fkana = sampleField(2);
	EigenVector hgana = sampleField(3);
	EigenVector hkana = sampleField(4);
	

	EigenVector ggdiff = ggana - groundedGradient;
	EigenVector fkdiff = fkana - fluxlessKnot;
	EigenVector hgdiff = hgana - harmonicGradient;
	EigenVector hkdiff = hkana - harmonicKnot;
	

	cout << "*** Numerical Errors ... ***" << endl;
	cout << "GG Error: " << sqrt(ggdiff.dot(m_mesh_3.CHS2F_B()*ggdiff)) / sqrt(ggana.dot(m_mesh_3.CHS2F_B()*ggana)) << endl;
	cout << "FK Error: " << sqrt(fkdiff.dot(m_mesh_3.CHS2F_B()*fkdiff)) / sqrt(fkana.dot(m_mesh_3.CHS2F_B()*fkana)) << endl;
	cout << "HG Error: " << sqrt(hgdiff.dot(m_mesh_3.HS2F_B()*hgdiff)) / sqrt(hgana.dot(m_mesh_3.HS2F_B()*hgana)) << endl;
	cout << "HK Error: " << sqrt(hkdiff.dot(m_mesh_3.HS2F_B()*hkdiff)) / sqrt(hkana.dot(m_mesh_3.HS2F_B()*hkana)) << endl;
}

void Decomposition::computeL2Norm()
{
	double inl2 = omegaOriginal.transpose()*m_mesh_3.CHS2F_B()*omegaOriginal;
	double ggl2 = groundedGradient.transpose()*m_mesh_3.CHS2F_B()*groundedGradient;
	double fkl2 = fluxlessKnot.transpose()*m_mesh_3.CHS2F_B()*fluxlessKnot;
	double hgl2 = harmonicGradient.transpose()*m_mesh_3.CHS2F_B()*harmonicGradient;
	double hkl2 = harmonicKnot.transpose()*m_mesh_3.CHS2F_B()*harmonicKnot;
	double cgl2 = curlyGradient.transpose()*m_mesh_3.CHS2F_B()*curlyGradient;

	double diff = inl2 - ggl2 - fkl2 - hgl2 - hkl2 - cgl2;

	cout << "Input L2-Norm: " << sqrt(inl2) << endl;
	cout << "Grounded Gradient L2-Norm: " << sqrt(ggl2) << endl;
	cout << "Fluxless Knot L2-Norm: " << sqrt(fkl2) << endl;
	cout << "Harmonic Gradient L2-Norm: " << sqrt(hgl2) << endl;
	cout << "Harmonic Knot L2-Norm: " << sqrt(hkl2) << endl;
	cout << "Curly Gradient L2-Norm: " << sqrt(cgl2) << endl;
	cout << "Summation Difference: " << diff << endl;

}

void Decomposition::testEnergy()
{
	EigenMatrix divEM = m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_B();
	EigenMatrix curlEM = m_mesh_3.CHS2F_NB()* m_mesh_3.ED1F_NB() *m_mesh_3.HS1F_NB().cwiseInverse() *m_mesh_3.ED1F_NB().transpose()*m_mesh_3.CHS2F_NB();

	EigenVector omegaNB = omega;
	boundaryToNoBoundary2Form(omegaNB);

	//EigenVector wtf = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.CHS2F_NB()*omegaNB;
	//for (int i = 0; i < wtf.size(); ++i) {
	//	if (wtf(i) > 10e-10) {
	//		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(m_mesh_3.EdgeBackward(i));
	//		Mesh_3_Vertex_iterator vi1 = ei->first->vertex(ei->second);
	//		Mesh_3_Vertex_iterator vi2 = ei->first->vertex(ei->third);
	//		cout << m_mesh_3.VertexOnBoundary(vi1) << " " << m_mesh_3.VertexOnBoundary(vi2) << endl;
	//	}
	//}
	//for (int i = 0; i < m_mesh_3.NumE(); ++i) {
	//	cout << wtf(i) << " " << m_mesh_3.EdgeOnBoundary(m_mesh_3.IdxToEdge(i)) << endl;
	//}

	double de = omega.dot(divEM*omega);
	double ce = omegaNB.dot(curlEM*omegaNB);

	cout << "Input Div Energy: " << de << endl;
	cout << "Input Curl Energy: " << ce << endl;
}

void Decomposition::testNonDiagonalHodgeStar()
{
	// -----------------------------------------------------
	double vol = 0;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		vol += m_mesh_3.CellPrimal(m_mesh_3.IdxToCell(i));
	}
	cout << "Cell Volume: " << vol << endl;

	// -----------------------------------------------------
	EigenVector form0;
	form0.resize(m_mesh_3.NumV());
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		form0(i) = 1;
	}
	double in0 = form0.dot(m_mesh_3.CHS0F_B()*form0);
	cout << "0-form L2: " << in0 << endl;

	// -----------------------------------------------------
	EigenVector form1;
	form1.resize(m_mesh_3.NumE());
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		Mesh_3_Vertex_iterator v1 = ei->first->vertex(ei->second);
		Mesh_3_Vertex_iterator v2 = ei->first->vertex(ei->third);

		int vidx1 = m_mesh_3.VertexIdx(v1);
		int vidx2 = m_mesh_3.VertexIdx(v2);
	
		Mesh_3_Vector_3 d = v1->point().point() - v2->point().point();

		if (vidx1 > vidx2)
			d = -d;

		Mesh_3_Vector_3 f(1, 0, 0);
		form1(i) = d * f;
	}
	double in1 = form1.dot(m_mesh_3.CHS1F_B()*form1);
	cout << "1-form L2: " << in1 << endl;

	double in2 = omega.dot(m_mesh_3.CHS2F_B()*omega);
	cout << "2-form L2: " << in2 << endl;

	EigenVector form3;
	form3.resize(m_mesh_3.NumC());
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		form3(i) = 1;
	}
	double in3 = form3.dot(m_mesh_3.CHS3F()*form3);
	cout << "3-form L2: " << in3 << endl;
}

void Decomposition::testProjectionMatrix()
{
	EigenSpMat idV;
	idV.resize(m_mesh_3.NumV(), m_mesh_3.NumV());
	idV.setIdentity();

	EigenSpMat idIV;
	idIV.resize(m_mesh_3.NumIV(), m_mesh_3.NumIV());
	idIV.setIdentity();


	//EigenSpMat A = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*m_mesh_3.ED1F_B() +
	//	m_mesh_3.CHS1F_B()*m_mesh_3.ED0F_B()*idV*m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B();
	//EigenVector b = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*omega;

	//for (int i = 0; i < A.outerSize(); ++i) {
	//	for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
	//		int row = iter.row();
	//		int col = iter.col();

	//		Mesh_3_Edge_iterator eir = m_mesh_3.IdxToEdge(row);
	//		Mesh_3_Edge_iterator eic = m_mesh_3.IdxToEdge(col);

	//		if (m_mesh_3.EdgeOnBoundary(eir) || m_mesh_3.EdgeOnBoundary(eir)) {
	//			if (row == col)
	//				iter.valueRef() = 1;
	//			else
	//				iter.valueRef() = 0;
	//		}
	//	}
	//}
	//for (int i = 0; i < m_mesh_3.NumE(); ++i) {
	//	Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

	//	if (m_mesh_3.EdgeOnBoundary(ei))
	//		b(i) = 0;
	//}


	//EigenVector x = matlabLinearSolveCholesky(A, b);

	//fluxlessKnot = m_mesh_3.ED1F_B()*x;
	////noBoundaryToBoundary2Form(fluxlessKnot);
	//omega -= fluxlessKnot;

	// --------------------------------------------------------------------

	// construct projection matrix
	std::vector<Eigen::Triplet<double>> triplet_P0;
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (m_mesh_3.VertexOnBoundary(vi))
			continue;

		int forwardIdx = m_mesh_3.VertexForward(i);
		triplet_P0.push_back(Eigen::Triplet<double>(forwardIdx, i, 1));
	}
	EigenSpMat P0;
	P0.resize(m_mesh_3.NumIV(), m_mesh_3.NumV());
	P0.setFromTriplets(triplet_P0.begin(), triplet_P0.end());

	std::vector<Eigen::Triplet<double>> triplet_P1;
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (m_mesh_3.EdgeOnBoundary(ei))
			continue;

		int forwardIdx = m_mesh_3.EdgeForward(i);
		triplet_P1.push_back(Eigen::Triplet<double>(forwardIdx, i, 1));
	}
	EigenSpMat P1;
	P1.resize(m_mesh_3.NumIE(), m_mesh_3.NumE());
	P1.setFromTriplets(triplet_P1.begin(), triplet_P1.end());

	std::vector<Eigen::Triplet<double>> triplet_P2;
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);

		if (m_mesh_3.FacetOnBoundary(fi))
			continue;

		int forwardIdx = m_mesh_3.FacetForward(i);
		triplet_P2.push_back(Eigen::Triplet<double>(forwardIdx, i, 1));
	}
	EigenSpMat P2;
	P2.resize(m_mesh_3.NumIF(), m_mesh_3.NumF());
	P2.setFromTriplets(triplet_P2.begin(), triplet_P2.end());

	// construct S_{2,int} from Projection
	//EigenSpMat S1 = P2 *m_mesh_3.CHS2F_B()*P2.transpose();
	//EigenSpMat S2 = m_mesh_3.CHS2F_NB();

	//EigenSpMat DM = S2 - S1;
	
	EigenSpMat L_1n_1 = 
		(P1*m_mesh_3.ED1F_B().transpose()*P2.transpose())
		*(P2*m_mesh_3.HS2F_B()*P2.transpose())
		*(P2*m_mesh_3.ED1F_B()*P1.transpose())
		+
		(P1*m_mesh_3.HS1F_B()*P1.transpose())
		*(P1*m_mesh_3.ED0F_B()*P0.transpose())
		*(P0*idV*P0.transpose())
		*(P0*m_mesh_3.ED0F_B().transpose()*P1.transpose())
		*(P1*m_mesh_3.HS1F_B()*P1.transpose());

	EigenSpMat L_1n_2 = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB() +
		m_mesh_3.HS1F_NB()*m_mesh_3.ED0F_NB()*idIV*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB();

	cout << "Size: " << L_1n_1.rows() << " " << L_1n_1.cols() << endl;
	cout << "Size: " << L_1n_2.rows() << " " << L_1n_2.cols() << endl;

	EigenSpMat DM = L_1n_1 - L_1n_2;

	cout << "Frobenius Norm: " << DM.norm() << endl;
}

EigenVector Decomposition::sampleField(int fieldType)
{
	EigenVector form;
	form.resize(m_mesh_3.NumF());

	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);
		form(i) = sampleFacetFlux(fi, fieldType, 3);
		/*form(i) = 0;
		
		if (fieldType == 0) {
			form(i) += fluxNormalGradient(fi);
			form(i) += fluxTangentialCurl(fi);
			form(i) += fluxCentralHarmonic(fi);
		}
		else if (fieldType == 1) {
			form(i) += fluxNormalGradient(fi);
		}
		else if (fieldType == 2) {
			form(i) += fluxTangentialCurl(fi);
		}
		else if (fieldType == 3) {
			form(i) += fluxCentralHarmonic(fi);
		}*/
	}

	return form;
}

double Decomposition::sampleFacetFlux(Mesh_3_Facet_iterator fi, int fieldType, int reso)
{
	// compute normal;
	std::vector<int> vidx;
	for (int i = 0; i < 4; ++i) {
		if (i != fi->second)
			vidx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i)));
	}

	std::sort(vidx.begin(), vidx.end());

	Mesh_3_Vector_3 p1 = m_mesh_3.IdxToVertex(vidx[0])->point().point() - CGAL::ORIGIN;
	Mesh_3_Vector_3 p2 = m_mesh_3.IdxToVertex(vidx[1])->point().point() - CGAL::ORIGIN;
	Mesh_3_Vector_3 p3 = m_mesh_3.IdxToVertex(vidx[2])->point().point() - CGAL::ORIGIN;

	Mesh_3_Vector_3 v1 = m_mesh_3.IdxToVertex(vidx[1])->point().point() - m_mesh_3.IdxToVertex(vidx[0])->point().point();
	Mesh_3_Vector_3 v2 = m_mesh_3.IdxToVertex(vidx[2])->point().point() - m_mesh_3.IdxToVertex(vidx[0])->point().point();

	Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
	n /= sqrt(n*n);

	//sampling
	double flux = 0;
	int num = reso * 3;

	std::vector<Mesh_3_Vector_3> coords;
	generateBaryCoords(coords, reso);

	//cout<<"@@@@@@@@@@@@@@@@@@@"<<endl;
	for (int i = 0; i < coords.size(); ++i) {	
		Mesh_3_Vector_3 bary = coords[i];
		//cout << "Bary Coord: " << bary.x() << " " << bary.y() << " " << bary.z() << endl;

		Mesh_3_Vector_3 p = p1 * bary.x() + p2 * bary.y() + p3 * bary.z();

		if (fieldType == 0)
			flux += n * (normalGradient(p) + tangentialCurl(p))*m_mesh_3.FacetPrimal(fi)*(1.0 / (reso*reso));
		else if (fieldType == 1)
			flux += n * normalGradient(p)*m_mesh_3.FacetPrimal(fi)*(1.0 / (reso*reso));
		else if (fieldType == 2)
			flux += n * tangentialCurl(p)*m_mesh_3.FacetPrimal(fi)*(1.0 / (reso*reso));
		else if (fieldType == 3)
			flux += n * normalHarmonic(p)*m_mesh_3.FacetPrimal(fi)*(1.0 / (reso*reso));
		else if (fieldType == 4)
			flux += n * tangentialHarmonic(p)*m_mesh_3.FacetPrimal(fi)*(1.0 / (reso*reso));
		else if (fieldType == 5)
			flux += n * centralHarmonic(p)*m_mesh_3.FacetPrimal(fi)*(1.0 / (reso*reso));
	}

	return flux;
}

void Decomposition::generateBaryCoords(std::vector<Mesh_3_Vector_3>& coords, int reso)
{
	std::vector<Mesh_3_Vector_3> cBary;
	cBary.push_back(Mesh_3_Vector_3((3 * reso - 2)*1.0 / (3 * reso), 1.0 / (3 * reso), 1.0 / (3 * reso)));
	coords.insert(coords.end(), cBary.begin(), cBary.end());

	double unit = 1.0 / (3 * reso);
	for (int i = 0; i < 2 * (reso - 1); ++i) {
		std::vector<Mesh_3_Vector_3> tBary;
		if (i % 2 == 0) {
			for (int j = 0; j < cBary.size(); ++j) {
				Mesh_3_Vector_3 b = cBary[j];
				Mesh_3_Vector_3 bb(b.x() - 2 * unit, b.y() + unit, b.z() + unit);
				tBary.push_back(bb);
			}
		}
		else {
			for (int j = 0; j < cBary.size(); ++j) {
				Mesh_3_Vector_3 b = cBary[j];
				Mesh_3_Vector_3 bb(b.x() - unit, b.y() + 2 * unit, b.z() - unit);
				tBary.push_back(bb);
			}
			Mesh_3_Vector_3 lb = cBary[cBary.size() - 1];
			Mesh_3_Vector_3 lbb(lb.x() - unit, lb.y() - unit, lb.z() + 2 * unit);
			tBary.push_back(lbb);
		}

		cBary.clear();
		cBary = tBary;
		coords.insert(coords.end(), cBary.begin(), cBary.end());
	}
}