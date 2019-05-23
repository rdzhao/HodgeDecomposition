#include "LaplacianAnalyzer.h"

My_Triangulation& LaplacianAnalyzer::Mesh()
{
	return m_mesh_3;
}

std::vector<Mesh_3_Point_3>& LaplacianAnalyzer::Vertices()
{
	return vertices;
}

std::vector<int>& LaplacianAnalyzer::Indices()
{
	return indices;
}

std::vector<Mesh_3_Vector_3>& LaplacianAnalyzer::Normals()
{
	return normals;
}

std::vector<Mesh_3_Vector_3>& LaplacianAnalyzer::Colors()
{
	return colors;
}

std::vector<Mesh_3_Vector_3>& LaplacianAnalyzer::VfVertices()
{
	return vfVertices;
}

std::vector<Mesh_3_Vector_3>& LaplacianAnalyzer::VfFaces()
{
	return vfFaces;
}

std::vector<Mesh_3_Vector_3>& LaplacianAnalyzer::VfNormals()
{
	return vfNormals;
}

std::vector<Mesh_3_Vector_3>& LaplacianAnalyzer::VfColors()
{
	return vfColors;
}

void LaplacianAnalyzer::init(double dsty, double sr, int stps, std::string sfx)
{
	m_density_ratio = dsty;
	m_step_ratio = sr;
	m_steps = stps;
	m_suffix = sfx;
}

void LaplacianAnalyzer::buildMeshFromSurface(std::string fn, double size)
{
	//Polyhedron polyhedron;
	//std::ifstream input(fn.c_str());
	//input >> polyhedron;
	//input.close();
	//Polyhedron_mesh_domain domain(polyhedron);
	//
	//Mesh_criteria criteria(
	//	facet_angle = 30,
	//	facet_size = size,
	//	facet_distance = size,
	//	cell_radius_edge_ratio = 2,
	//	cell_size = size);
	//std::cout << "Criteria ... " << std::endl;

	//m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::odt());
	//CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);

	// -----------------------------------

	//Polyhedron_with_feature polyhedron;
	//std::ifstream input(fn.c_str());
	//input >> polyhedron;
	//input.close();
	//Mesh_domain_with_feature domain(polyhedron);
	////domain.detect_features();

	//Mesh_criteria criteria(/*edge_size = size,*/
	//	facet_angle = 20,
	//	facet_size = size,
	//	facet_distance = 0.2*size,
	//	cell_radius_edge_ratio = 3,
	//	cell_size = size,
	//	facet_topology = CGAL::MANIFOLD);
	//std::cout << "Criteria ... " << std::endl;

	//m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, features(domain), no_lloyd(), no_odt(), no_perturb(), no_exude());
	////CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	//std::cout << "Triangulation ... " << std::endl;

	// -------------------------------------

	std::ifstream in("test1.tet", std::ios::in);
	in >> m_mesh_3.C3T3();
	cout<<"Reading Done ..."<<endl;
	
	m_mesh_3.preprocessing();

	//std::ofstream out("output.tet", std::ios::out);
	//out << m_mesh_3.C3T3();

	writeTestMesh();
}

void LaplacianAnalyzer::buildDECOperators()
{
	//buildExteriorDerivative0Form();
	//buildExteriorDerivative1Form();
	//buildExteriorDerivative2Form();

	//buildHodgeStar0Form();
	//buildHodgeStar1Form();
	//buildHodgeStar2Form();
	//buildHodgeStar3Form();
}

void LaplacianAnalyzer::buildLaplacian()
{
	buildLaplacian0();
	buildLaplacian1();
	buildLaplacian2();
	buildLaplacian3();
}

void LaplacianAnalyzer::decompose()
{
	setFormRandom();
	//sampleInput();

	//computeCurlFreeGradient();
	cout << "Harmonic Knot ..." << endl;
	computeHarmonicKnot();
	//computeHarmonicKnotTest();
	cout << "Fluxless Knot ..." << endl;
	computeFluxlessKnot();
	cout << "Grounded Gradient ..." << endl;
	computeGroundedGradient();
	cout<<"Harmonic Gradient ..."<<endl;
	computeHarmonicGradient();
	cout << "Curly Gradient ..." << endl;
	computeCurlyGradient();

	double inner = fluxlessKnot.transpose() * (m_mesh_3.HS1F_B() *groundedGradient);
	double fklen = fluxlessKnot.transpose() * (m_mesh_3.HS1F_B() *fluxlessKnot);
	double gglen = groundedGradient.transpose() * (m_mesh_3.HS1F_B() *groundedGradient);

	cout << "Inner: " << inner / (fklen*gglen) << endl;

	//computeCurlFreeGradient();
	//computeRemainingDivFreeFieldPotential();

	//// for mix boundary 
	//assignBoundaryVerticesGroup();
	//assignVertexColorMixBoundary();
	//computeHarmonicMixBoundary();
	//computeExactMixBoundary();
	//computeCoexactMixBoundary();

	//// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
}

void LaplacianAnalyzer::setEigenfield(int mode, int idx, bool isCurl)
{
	if (mode == 0) { // L_0n
		setEigenfield_L0n(idx, isCurl);
	}

	if (mode == 1) { // L_0t
		setEigenfield_L0t(idx, isCurl);
	}

	if (mode == 2) { // L_1n
		setEigenfield_L1n(idx, isCurl);
	}

	if (mode == 3) { // L_1t
		setEigenfield_L1t(idx, isCurl);
	}

	if (mode == 4) { // L_2n
		setEigenfield_L2n(idx, isCurl);
	}

	if (mode == 5) { // L_2t
		setEigenfield_L2t(idx, isCurl);
	}

	if (mode == 6) { // L_3n
		setEigenfield_L3n(idx, isCurl);
	}

	if (mode == 7) { // L_3t
		setEigenfield_L3t(idx, isCurl);
	}

	if (mode == -1) {
		EigenVector form;
		Mesh_3_Vector_3 d(1, 0, 0);
		setConstantField(form, d);

		omega.resize(m_mesh_3.NumE());
		omega.setZero();
		omega += form;

		double inner = omega.dot(m_mesh_3.HS1F_B()*omega);
		cout << "L2 Inner product: " << inner << endl;

		EigenVector form0;
		form0.resize(m_mesh_3.NumV());
		for (int i = 0; i < m_mesh_3.NumV(); ++i) {
			form0(i) = 1;
		}
		double inner0 = form0.dot(m_mesh_3.HS0F_B()*form0);
		cout << "L2 Inner product: " << inner0 << endl;

		EigenVector form3;
		form3.resize(m_mesh_3.NumC());
		for (int i = 0; i < m_mesh_3.NumC(); ++i) {
			form3(i) = m_mesh_3.CellPrimal(m_mesh_3.IdxToCell(i));
		}

		double inner3 = form3.dot(m_mesh_3.HS3F()*form3);
		cout << "L2 Inner product: " << inner3 << endl;
	}
}

void LaplacianAnalyzer::visualize()
{
	//convertForms();
	//computeArrows();
	computeCrossSection();

	//writeMesh();
	writeMeshPLY();
}

void LaplacianAnalyzer::integrate()
{
	//integrateField(omegaVF, omegaCF, m_density_ratio, m_step_ratio, m_steps, 2, m_suffix + "_o");

	//integrateField(fluxlessKnotVF, fluxlessKnotCF, m_density_ratio, m_step_ratio, m_steps, 2, m_suffix + "_fk");

	//integrateField(groundedGradientVF, groundedGradientCF, m_density_ratio, m_step_ratio, m_steps, 2, m_suffix + "_gg");

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		omegaCF[vi] = Mesh_3_Vector_3(0.4, 0.4, 0.4);
	}

	integrateField(gamma_1_VF, gamma_1_CF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_g1");
	integrateField(gamma_2_VF, gamma_2_CF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_g2");
	integrateField(gamma_3_VF, gamma_3_CF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_g3");
	integrateField(gamma_4_VF, gamma_4_CF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_g4");
}

void LaplacianAnalyzer::process_Poelke()
{
	// create example 1 form
	EigenVector field;
	field.setZero(m_mesh_3.NumF());
	setEigenfield_L2n(0, true);
	field += gamma_1;
	cout<<"###############"<<endl;
	setEigenfield_L2t(0, false);
	cout << "###############" << endl;
	field += gamma_1;
	gamma_1 = field;

	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf;
	convert2Form(gamma_1, vf);
	EigenVector vfcol;
	vfcol.resize(3 * m_mesh_3.NumC());
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);
		vfcol(3 * i) = vf[ci].x();
		vfcol(3 * i + 1) = vf[ci].y();
		vfcol(3 * i + 2) = vf[ci].z();
	}


	//// test ED
	//boundaryToNoBoundary1Form(gamma_1);
	//EigenVector dg1 = m_mesh_3.ED1F_NB()*gamma_1;
	//gamma_1 = dg1;
	//noBoundaryToBoundary2Form(gamma_1);
	//// ---------------------

	//// compute curl: gamma_1 is the form
	//EigenSpMat curlOp;
	//cout<<"Creating Curl Operator ..."<<endl;
	//makeCurlOperator_Poelke(curlOp);
	//cout << "Creating Curl Operator Done ..." << endl;
	//boundaryToNoBoundary1Form(gamma_1);
	//EigenVector fv = curlOp * gamma_1;

	//for (int i = 0; i < m_mesh_3.NumC(); ++i) {
	//	Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

	//	gamma_1_VF[ci] = Mesh_3_Vector_3(fv(3 * i), fv(3 * i + 1), fv(3 * i + 2));
	//}
	//// -------------------------

	// compute curl part
	EigenSpMat curlOp;
	cout<<"Creating Curl Operator ..."<<endl;
	makeCurlOperator_Poelke(curlOp);
	cout << "Creating Curl Operator Done ..." << endl;
	EigenSpMat A = curlOp.transpose()*curlOp;
	EigenSpMat B = m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB();
	EigenVector a = curlOp.transpose()*vfcol;
	EigenVector b;
	b.setZero(m_mesh_3.NumIV());

	//cout << "$$$$$$$$$$$$" << endl;
	EigenSpMat M;
	concatenateSPMatrix(M, A, B);
	//cout<<"$$$$$$$$$$$$ "<<M.innerSize()<<" "<<M.outerSize()<<endl;
	EigenVector v;
	concatenateVector(v,a,b);
	//cout << "$$$$$$$$$$$$ " << v.size() << endl;
	//cout<<v<<endl;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::StorageIndex>> solver;
	solver.compute(M);
	EigenVector x = solver.solve(v);
	//EigenVector d1 = A * x - a;
	//cout << "Diff: " << d1.transpose()*d1 << endl;
	
	auto diff = std::chrono::high_resolution_clock::now() - start;
	auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(diff);
	std::cout << "Time : " << 1.0*(t1.count()) / 1000000000 << std::endl;
	EigenVector cf = curlOp * x;
	//cout << x.size() << endl;
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> curl_VF;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		curl_VF[ci] = Mesh_3_Vector_3(cf(3 * i), cf(3 * i + 1), cf(3 * i + 2));
	}
	// ---------------------------------------

	// compute grad part
	EigenSpMat gradOp;
	makeGradOperator_Poelke(gradOp);
	EigenSpMat AA = gradOp.transpose()*gradOp;
	EigenVector bb = gradOp.transpose()*vfcol;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_grad;
	solver_grad.compute(AA);
	if (solver_grad.info() != Eigen::Success) {
		cout<<"Solver Fail ..."<<endl;
	}
	EigenVector xx = solver_grad.solve(bb);
	EigenVector d2 = AA * xx - bb;
	cout << "Diff: " << d2.transpose()*d2 << endl;
	EigenVector gf = gradOp * xx;

	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		gamma_1_VF[ci] = Mesh_3_Vector_3(gf(3 * i), gf(3 * i + 1), gf(3 * i + 2));
	}
	// -----------------------------------------

	// evaluate inner product.
	double inner = 0;
	double inner1 = 0;
	double inner2 = 0;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);
		
		//Mesh_3_Vector_3 dv = vf[ci] - gamma_1_VF[ci] - curl_VF[ci];

		inner += (gamma_1_VF[ci] * curl_VF[ci])*m_mesh_3.CellPrimal(ci);
		inner1 += (gamma_1_VF[ci]* gamma_1_VF[ci])*m_mesh_3.CellPrimal(ci);
		inner2 += (curl_VF[ci] * curl_VF[ci])*m_mesh_3.CellPrimal(ci);
		//inner += (dv*dv)*m_mesh_3.CellPrimal(ci);
	}
	cout << "Inner Product: " << inner / sqrt(inner1*inner2) << endl;

	EigenVector evf;
	evf.resize(m_mesh_3.NumIF());
	for (int i = 0; i < m_mesh_3.NumIF(); ++i) {
		evf(i) = 0;
	}

	EigenVector eve;
	eve.resize(m_mesh_3.NumIE());
	for (int i = 0; i < m_mesh_3.NumIE(); ++i) {
		eve(i) = 0;
	}


	Mesh_3_Facet_iterator fi;
	Mesh_3_Edge_iterator ei;
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		ei = m_mesh_3.IdxToEdge(i);

		if (m_mesh_3.EdgeOnBoundary(ei))
			continue;

		int idx;
		for (int i = 0; i < 4; ++i) {
			if (i != ei->second && i != ei->third) {
				idx = i;
			}
		}

		My_facet mf(ei->first, idx);
		fi = m_mesh_3.MyCellFacetMap(mf);

		break;
	}

	evf(m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi))) = 1;
	eve(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei))) = 1;

	EigenVector xf = gradOp * evf;
	EigenVector xe = curlOp * eve;

	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		curl_VF[ci] = Mesh_3_Vector_3(xf(3 * i), xf(3 * i + 1), xf(3 * i + 2));
	}

	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		gamma_1_VF[ci] = Mesh_3_Vector_3(xe(3 * i), xe(3 * i + 1), xe(3 * i + 2));
	}

	// evaluate inner product.
	double inn = 0;
	double inn1 = 0;
	double inn2 = 0;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);
		
		inn += (gamma_1_VF[ci] * curl_VF[ci])*m_mesh_3.CellPrimal(ci);
		inn1 += (gamma_1_VF[ci] * gamma_1_VF[ci])*m_mesh_3.CellPrimal(ci);
		inn2 += (curl_VF[ci] * curl_VF[ci])*m_mesh_3.CellPrimal(ci);
	}
	cout<<"Inner Product: "<< inn / sqrt(inn1*inn2) <<endl;

	//// my method for curl
	//EigenSpMat s1l1nb = m_mesh_3.HS1F_NB() * m_laplacian1_NB;

	//cout << "###########" << endl;
	//EigenVector omegaNB = gamma_1;
	//boundaryToNoBoundary2Form(omegaNB);

	//EigenVector b = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*omegaNB;

	//// solve
	//auto start = std::chrono::high_resolution_clock::now();
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//solver.compute(s1l1nb);
	//EigenVector x = solver.solve(b); // potential 1-form for exact component
	//auto diff = std::chrono::high_resolution_clock::now() - start;
	//auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(diff);
	//std::cout << "Time : " << 1.0*(t1.count()) / 1000000000 << std::endl;

	//gamma_1 = m_mesh_3.ED1F_NB() * x;
	//noBoundaryToBoundary2Form(gamma_1);

	//// -----------------------------------------

	//// compute div part
	//EigenSpMat gradOp;
	//cout<<"Creating Curl Operator ..."<<endl;
	//makeGradOperator_Poelke(gradOp);
	//cout << "Creating Curl Operator Done ..." << endl;
	//EigenSpMat A = gradOp.transpose()*gradOp;
	//EigenVector b = gradOp.transpose()*vfcol;

	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//solver.compute(A);
	//if (solver.info() != Eigen::Success) {
	//	cout<<"Solver Fail ..."<<endl;
	//}

	//EigenVector x = solver.solve(b);
	//EigenVector cf = gradOp * x;
	//for (int i = 0; i < m_mesh_3.NumC(); ++i) {
	//	Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

	//	gamma_1_VF[ci] = Mesh_3_Vector_3(cf(3 * i), cf(3 * i + 1), cf(3 * i + 2));
	//}


}

void LaplacianAnalyzer::buildExteriorDerivative0Form()
{
	//// orientation: increasing global indices
	//std::vector<Eigen::Triplet<double>> triplet;

	//int ii = 0;
	//for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
	//	ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
	//	if (!m_mesh_3.EdgeValid(ei))
	//		continue;

	//	if (!m_mesh_3.EdgeOnBoundary(ei)) {
	//		if (m_mesh_3.VertexIdx(ei->first->vertex(ei->second)) > m_mesh_3.VertexIdx(ei->first->vertex(ei->third))) {
	//			if (!m_mesh_3.VertexOnBoundary(ei->first->vertex(ei->second)))
	//				triplet.push_back(Eigen::Triplet<double>(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)),
	//					m_mesh_3.VertexForward(m_mesh_3.VertexIdx(ei->first->vertex(ei->second))), 1));
	//			if (!m_mesh_3.VertexOnBoundary(ei->first->vertex(ei->third)))
	//				triplet.push_back(Eigen::Triplet<double>(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)),
	//					m_mesh_3.VertexForward(m_mesh_3.VertexIdx(ei->first->vertex(ei->third))), -1));
	//		}
	//		else {
	//			if (!m_mesh_3.VertexOnBoundary(ei->first->vertex(ei->second)))
	//				triplet.push_back(Eigen::Triplet<double>(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)),
	//					m_mesh_3.VertexForward(m_mesh_3.VertexIdx(ei->first->vertex(ei->second))), -1));
	//			if (!m_mesh_3.VertexOnBoundary(ei->first->vertex(ei->third)))
	//				triplet.push_back(Eigen::Triplet<double>(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)),
	//					m_mesh_3.VertexForward(m_mesh_3.VertexIdx(ei->first->vertex(ei->third))), 1));
	//		}
	//	}
	//	++ii;
	//}

	//r_ed0f.resize(m_mesh_3.NumIE(), m_mesh_3.NumIV());
	//r_ed0f.setFromTriplets(triplet.begin(), triplet.end());
}

void LaplacianAnalyzer::buildExteriorDerivative1Form()
{
	//// orientation: increasing global indices 
	//std::vector<Eigen::Triplet<double>> triplet;

	//int ii = 0;
	//for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
	//	fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
	//	if (!m_mesh_3.FacetValid(fi))
	//		continue;

	//	if (m_mesh_3.FacetOnBoundary(fi)) // !!
	//		continue;

	//	std::vector<int> vGlobalIdx;
	//	std::vector<int> vLocalIdx;
	//	for (int i = 0; i < 4; ++i) {
	//		if (i == fi->second)
	//			continue;

	//		vGlobalIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i)));
	//		vLocalIdx.push_back(i);
	//	}

	//	// make increasing order to be consistent with orientation
	//	for (int i = 0; i < vGlobalIdx.size(); ++i) {
	//		for (int j = i + 1; j < vGlobalIdx.size(); ++j) {
	//			if (vGlobalIdx[i] > vGlobalIdx[j]) { // swap
	//				int tgi, tli;

	//				tgi = vGlobalIdx[i];
	//				vGlobalIdx[i] = vGlobalIdx[j];
	//				vGlobalIdx[j] = tgi;

	//				tli = vLocalIdx[i];
	//				vLocalIdx[i] = vLocalIdx[j];
	//				vLocalIdx[j] = tli;
	//			}
	//		}
	//	}

	//	// integrate based on oritation	
	//	for (int i = 0; i < vGlobalIdx.size(); ++i) {
	//		int localIdx1, localIdx2;
	//		localIdx1 = vLocalIdx[i];
	//		localIdx2 = vLocalIdx[(i + 1) % vLocalIdx.size()];

	//		if (localIdx1 > localIdx2)
	//			std::swap(localIdx1, localIdx2);

	//		int val = 1;
	//		if (i == vGlobalIdx.size() - 1)
	//			val = -1;

	//		My_edge me(fi->first, localIdx1, localIdx2);
	//		Mesh_3_Edge_iterator ei = m_mesh_3.MyCellEdgeMap(me);

	//		if (!m_mesh_3.EdgeOnBoundary(ei))
	//			triplet.push_back(Eigen::Triplet<double>(m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), val));

	//	}

	//	++ii;
	//}

	//r_ed1f.resize(m_mesh_3.NumIF(), m_mesh_3.NumIE()); // !!
	//r_ed1f.setFromTriplets(triplet.begin(), triplet.end());
}

void LaplacianAnalyzer::buildExteriorDerivative2Form()
{
	//// orientation inherited from local index 0123
	//std::vector<Eigen::Triplet<double>> triplet;

	//int ii = 0;
	//for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
	//	ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
	//	if (ci->subdomain_index() == 0)
	//		continue;

	//	for (int i = 0; i < 4; ++i) {
	//		My_facet mf(ci, i);

	//		int even = 1;
	//		if (i % 2 == 1)
	//			even = -1;

	//		std::vector<int> lIdx;
	//		for (int j = 0; j < 4; ++j) {
	//			if (j != i)
	//				lIdx.push_back(j);
	//		}

	//		std::vector<int> gIdx;
	//		gIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(lIdx[0])));
	//		gIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(lIdx[1])));
	//		gIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(lIdx[2])));

	//		//int even = 1;

	//		for (int i = 0; i < gIdx.size(); ++i) {
	//			for (int j = i + 1; j < gIdx.size(); ++j) {
	//				if (gIdx[i] > gIdx[j]) {
	//					std::swap(gIdx[i], gIdx[j]);
	//					even = -even;
	//				}
	//			}
	//		}

	//		Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);

	//		if (!m_mesh_3.FacetOnBoundary(fi)) {
	//			if (fi->first == ci) {
	//				triplet.push_back(Eigen::Triplet<double>(ii, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(m_mesh_3.MyCellFacetMap(mf))), even));
	//			}
	//			else {
	//				triplet.push_back(Eigen::Triplet<double>(ii, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(m_mesh_3.MyCellFacetMap(mf))), even));
	//			}
	//		}
	//	}

	//	++ii;
	//}

	//r_ed2f.resize(m_mesh_3.NumC(), m_mesh_3.NumIF()); // !!
	//r_ed2f.setFromTriplets(triplet.begin(), triplet.end());
}

void LaplacianAnalyzer::buildHodgeStar0Form()
{
	/*std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (m_mesh_3.VertexOnBoundary(vi))
			continue;

		triplet.push_back(Eigen::Triplet<double>(m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi)),
			m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi)),
			m_mesh_3.VertexDual(vi)));
	}
	r_hs0f.resize(m_mesh_3.NumIV(), m_mesh_3.NumIV());
	r_hs0f.setFromTriplets(triplet.begin(), triplet.end());*/
}

void LaplacianAnalyzer::buildHodgeStar1Form()
{
	/*std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;

		if (!m_mesh_3.EdgeOnBoundary(ei))
			triplet.push_back(Eigen::Triplet<double>(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)),
				m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)),
				m_mesh_3.EdgeDual(ei) / m_mesh_3.EdgePrimal(ei)));

	}
	r_hs1f.resize(m_mesh_3.NumIE(), m_mesh_3.NumIE());
	r_hs1f.setFromTriplets(triplet.begin(), triplet.end());*/
}

void LaplacianAnalyzer::buildHodgeStar2Form()
{
	//std::vector<Eigen::Triplet<double>> triplet;

	//for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
	//	fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
	//	if (!m_mesh_3.FacetValid(fi))
	//		continue;

	//	if (!m_mesh_3.FacetOnBoundary(fi)) // !!
	//		triplet.push_back(Eigen::Triplet<double>(m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)),
	//			m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)),
	//			m_mesh_3.FacetDual(fi) / m_mesh_3.FacetPrimal(fi)));
	//}

	//r_hs2f.resize(m_mesh_3.NumIF(), m_mesh_3.NumIF()); // !!
	//r_hs2f.setFromTriplets(triplet.begin(), triplet.end());
}

void LaplacianAnalyzer::buildHodgeStar3Form()
{
	/*std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		triplet.push_back(Eigen::Triplet<double>(m_mesh_3.CellIdx(ci), m_mesh_3.CellIdx(ci), 1 / m_mesh_3.CellPrimal(ci)));
	}

	r_hs3f.resize(m_mesh_3.NumC(), m_mesh_3.NumC());
	r_hs3f.setFromTriplets(triplet.begin(), triplet.end());*/
}

void LaplacianAnalyzer::buildLaplacian0()
{
	m_laplacian0_NB = m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB()*m_mesh_3.ED0F_NB();
	m_laplacian0_B = m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B()*m_mesh_3.ED0F_B();
	
}

void LaplacianAnalyzer::buildLaplacian1()
{	
	m_laplacian1_NB = m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB()
		+ m_mesh_3.ED0F_NB()*m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB();

	m_laplacian1_B = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()
		+ m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();

	m_divEM_L1n = m_mesh_3.HS1F_NB()*m_mesh_3.ED0F_NB()*m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB();
	m_curlEM_L1n = m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB();

	m_divEM_L1t = m_mesh_3.HS1F_B()*m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();
	m_curlEM_L1t = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B();
}

void LaplacianAnalyzer::buildLaplacian2()
{
	//_laplacian_2 = r_hs2f.cwiseInverse() * r_ed2f.transpose()*r_hs3f*r_ed2f
		//+ r_ed1f * r_hs1f.cwiseInverse() *r_ed1f.transpose()*r_hs2f;

	m_laplacian2_B = /*100 **/ m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_B()
		+ m_mesh_3.ED1F_B()*m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B();
	m_laplacian2_NB = m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_NB()
		+ m_mesh_3.ED1F_NB()*m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB();

	m_divEM_L2n = m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_NB();
	m_curlEM_L2n = m_mesh_3.HS2F_NB()*m_mesh_3.ED1F_NB()*m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB();

	m_divEM_L2t = /*100**/m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_B();
	m_curlEM_L2t = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B();
}

void LaplacianAnalyzer::buildLaplacian3()
{
	m_laplacian3_B = m_mesh_3.ED2F_B()*m_mesh_3.HS2F_B().cwiseInverse()*m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F();
	m_laplacian3_NB = m_mesh_3.ED2F_NB()*m_mesh_3.HS2F_NB().cwiseInverse()*m_mesh_3.ED2F_NB().transpose()*m_mesh_3.HS3F();
}

void LaplacianAnalyzer::computeBoundaryComponent()
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

void LaplacianAnalyzer::setForm()
{
	//EigenSpMat s1l1nb = r_hs1f * _laplacian_1;

	//Engine* eng;
	//eng = engOpen("\0");

	//std::string str;

	//// assemble matrix to be eigen decomposed
	//std::vector<double> row, col;
	//std::vector<double> val;
	//for (int i = 0; i < s1l1nb.outerSize(); ++i) {
	//	for (EigenSpMat::InnerIterator iter(s1l1nb, i); iter; ++iter) {
	//		row.push_back(static_cast<double>(iter.row()) + 1);
	//		col.push_back(static_cast<double>(iter.col()) + 1);
	//		val.push_back(static_cast<double>(iter.value()));
	//	}
	//}
	//std::vector<double> hs;
	//for (int i = 0; i < m_mesh_3.NumIE(); ++i) {
	//	hs.push_back(r_hs1f.coeff(i, i));
	//}

	//// transfer data
	//mxArray *rowArr, *colArr, *valArr, *hsArr;
	//rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	//colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	//valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	//hsArr = mxCreateDoubleMatrix(hs.size(), 1, mxREAL);
	//memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	//memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	//memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	//memcpy(mxGetPr(hsArr), &hs[0], hs.size() * sizeof(double));
	//engPutVariable(eng, "row", rowArr);
	//engPutVariable(eng, "col", colArr);
	//engPutVariable(eng, "val", valArr);
	//engPutVariable(eng, "hs", hsArr);

	//str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());

	//str = "n=length(hs)";        engEvalString(eng, str.c_str());
	//str = "B = spdiags(hs(:),0,n,n);";        engEvalString(eng, str.c_str());
	//str = "[V,D]=eigs(A, B, 50, 'smallestabs')";        engEvalString(eng, str.c_str());

	//// read out eigenvector matrix D from matlab
	//mxArray* vv;
	//vv = engGetVariable(eng, "V");
	//double* v;
	//v = reinterpret_cast<double*>(mxGetData(vv));
	//const size_t* dims;
	//dims = mxGetDimensions(vv);

	//EigenVector ev; // no boundary
	//ev.resize(dims[0]);

	//harmonicBasis.clear();
	//for (int j = 0; j < dims[1]; ++j) {
	//	for (int i = 0; i < dims[0]; ++i) {
	//		ev(i) = v[i + j * dims[0]];
	//	}

	//	harmonicBasis.push_back(ev);
	//}

	//// randomly choose eigen basis for assembling omega.
	//omega.setZero(m_mesh_3.NumIE());
	//omega += harmonicBasis[0];
	//omega += 0.5*harmonicBasis[2];
	////omega += 0.8*harmonicBasis[4];
	////omega += 0.2*harmonicBasis[6];
	//omega += 1.5*harmonicBasis[8];
	//omega += 0.7*harmonicBasis[10];

	////omega += 0.5*harmonicBasis[11];
	////omega += 0.7*harmonicBasis[12];
	//omega += 0.8*harmonicBasis[20];
	//omega += 0.4*harmonicBasis[22];
	//omega += 1*harmonicBasis[24];
	//omega += 0.8*harmonicBasis[26];

	//noBoundaryToBoundary1Form(omega);
}

void LaplacianAnalyzer::setFormRandom()
{
	omega.resize(m_mesh_3.NumE());
	omega.setZero();

	EigenVector form;

	Mesh_3_Point_3 pcp(0.2, 0.7, 0.6);
	setPointChargeElectricField(form, pcp, false);
	omega += 1.8*form;

	//pcp= Mesh_3_Point_3(0.4, 0.5, 0.5);
	//setPointChargeElectricField(form, pcp, true);
	//omega += 1.2*form;

	//pcp = Mesh_3_Point_3(0.6, 0.4, 0.5);
	//setPointChargeElectricField(form, pcp, true);
	//omega += 0.7*form;

	pcp = Mesh_3_Point_3(0.6, 0.3, 0.6);
	setPointChargeElectricField(form, pcp, true);
	omega += 1.3*form;

	cout << "Electric Done ..." << endl;

	//Mesh_3_Point_3 cp(0.2, 0.6, 0.6);
	//Mesh_3_Vector_3 cv(0, 0.2, 1);
	//setCurrentMagneticField(form, cp, cv, true);
	//omega += 0.5*form;

	//cp = Mesh_3_Point_3(0.7, 0.2, 0.5);
	//cv = Mesh_3_Vector_3(0.3, 0.7, 1);
	//setCurrentMagneticField(form, cp, cv, false);
	//omega += 0.7*form;

	//cp = Mesh_3_Point_3(2, 0, 0);
	//cv = Mesh_3_Vector_3(0, 0, 1);
	//setCurrentMagneticField(form, cp, cv, false);
	//omega += form;

	cout << "Magnetic Done ..." << endl;

	std::vector<EigenVector> eigenFields;
	prepareEigenFeilds(eigenFields);

	EigenVector cpnt;
	cpnt.resize(m_mesh_3.NumE());
	cpnt.setZero();
	//cpnt += 0.4*eigenFields[curlFieldIdx[8]];
	cpnt += 0.8*eigenFields[curlFieldIdx[12]];
	cpnt += 0.5*eigenFields[curlFieldIdx[14]];
	cpnt += 0.6*eigenFields[curlFieldIdx[16]];


	//cpnt += 0.8*eigenFields[divFieldIdx[3]];

	omega += cpnt;

	cout<<omega.size()<<endl;

	//smooth(omega);
	//cout << "Smooth Done" << endl;
}

void LaplacianAnalyzer::computeHarmonic()
{
	//EigenSpMat s1l1nb = r_hs1f * _laplacian_1;

	//double hm2 = m_mesh_3.NumHMLG2();

	//if (hm2 == 0) {
	//	harmonicComponent.resize(m_mesh_3.NumE(), m_mesh_3.NumE());
	//	harmonicComponent.setZero();

	//	noHarmonicComponent = omega - harmonicComponent;
	//}
	//else {
	//	Engine* eng;
	//	eng = engOpen("\0");

	//	std::string str;

	//	// assemble matrix to be eigen decomposed
	//	std::vector<double> row, col;
	//	std::vector<double> val;
	//	for (int i = 0; i < s1l1nb.outerSize(); ++i) {
	//		for (EigenSpMat::InnerIterator iter(s1l1nb, i); iter; ++iter) {
	//			row.push_back(static_cast<double>(iter.row()) + 1);
	//			col.push_back(static_cast<double>(iter.col()) + 1);
	//			val.push_back(static_cast<double>(iter.value()));
	//		}
	//	}
	//	std::vector<double> hs;
	//	for (int i = 0; i < m_mesh_3.NumIE(); ++i) {
	//		hs.push_back(r_hs1f.coeff(i, i));
	//	}

	//	// transfer data
	//	mxArray *rowArr, *colArr, *valArr, *hsArr, *hmArr;
	//	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	//	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	//	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	//	hsArr = mxCreateDoubleMatrix(hs.size(), 1, mxREAL);
	//	hmArr = mxCreateDoubleMatrix(1, 1, mxREAL);
	//	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	//	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	//	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	//	memcpy(mxGetPr(hsArr), &hs[0], hs.size() * sizeof(double));
	//	memcpy(mxGetPr(hmArr), &hm2, sizeof(double));
	//	engPutVariable(eng, "row", rowArr);
	//	engPutVariable(eng, "col", colArr);
	//	engPutVariable(eng, "val", valArr);
	//	engPutVariable(eng, "hs", hsArr);
	//	engPutVariable(eng, "hm2", hmArr);

	//	str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());
	//	str = "n=length(hs)";        engEvalString(eng, str.c_str());
	//	str = "B = spdiags(hs(:),0,n,n);";        engEvalString(eng, str.c_str());
	//	str = "[V,D]=eigs(A, B, hm2, 'smallestabs')";        engEvalString(eng, str.c_str());

	//	// read out eigenvector matrix D from matlab
	//	mxArray* vv;
	//	vv = engGetVariable(eng, "V");
	//	double* v;
	//	v = reinterpret_cast<double*>(mxGetData(vv));
	//	const size_t* dims;
	//	dims = mxGetDimensions(vv);

	//	EigenVector ev; // no boundary
	//	ev.resize(dims[0]);

	//	harmonicBasis.clear();
	//	for (int j = 0; j < dims[1]; ++j) {
	//		for (int i = 0; i < dims[0]; ++i) {
	//			ev(i) = v[i + j * dims[0]];
	//		}

	//		harmonicBasis.push_back(ev);
	//	}

	//	//turn original form to no boundary
	//	EigenVector omegaNB = omega;
	//	boundaryToNoBoundary1Form(omegaNB);

	//	// we have harmonic basis for tangential form with no boundary.
	//	EigenVector harmonicComponentNB;
	//	harmonicComponentNB.resize(m_mesh_3.NumIE());
	//	harmonicComponentNB.setZero();
	//	for (int i = 0; i < harmonicBasis.size(); ++i) {
	//		double c = (omegaNB.transpose()*r_hs1f*harmonicBasis[i])[0] / (harmonicBasis[i].transpose() * r_hs1f*harmonicBasis[i])[0];
	//		harmonicCoeffs.push_back(c);

	//		harmonicComponentNB += c * harmonicBasis[i];
	//	}

	//	harmonicComponent = harmonicComponentNB;
	//	noBoundaryToBoundary1Form(harmonicComponent);

	//	noHarmonicComponent = omega - harmonicComponent;
	//}
}

void LaplacianAnalyzer::computeCoexact()
{
	//EigenSpMat s2l2 = r_hs2f * _laplacian_2;

	//double hm1 = m_mesh_3.NumHMLG1();
	//if (hm1 == 0) {
	//	EigenVector noHarmonicComponentNB = noHarmonicComponent;
	//	boundaryToNoBoundary1Form(noHarmonicComponentNB);

	//	EigenVector b = r_hs2f * r_ed1f*noHarmonicComponentNB;

	//	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//	solver.compute(s2l2);

	//	EigenVector x = solver.solve(b); // potential 2-form for exact component
	//	EigenVector vp = x;
	//	noBoundaryToBoundary2Form(vp);
	//	convert2Form(vp, vectorPotential);

	//	coexactComponent = r_hs1f.cwiseInverse() *r_ed1f.transpose()*r_hs2f * x;
	//	noBoundaryToBoundary1Form(coexactComponent);
	//}
	//else {
	//	Engine* eng;
	//	eng = engOpen("\0");
	//	std::string str;

	//	std::vector<double> row, col;
	//	std::vector<double> val;
	//	for (int i = 0; i < s2l2.outerSize(); ++i) {
	//		for (EigenSpMat::InnerIterator iter(s2l2, i); iter; ++iter) {
	//			row.push_back(static_cast<double>(iter.row()) + 1);
	//			col.push_back(static_cast<double>(iter.col()) + 1);
	//			val.push_back(static_cast<double>(iter.value()));
	//		}
	//	}
	//	std::vector<double> hs;
	//	for (int i = 0; i < m_mesh_3.NumIF(); ++i) {
	//		hs.push_back(r_hs2f.coeff(i, i));
	//	}

	//	// transfer data
	//	mxArray *rowArr, *colArr, *valArr, *hsArr, *hmArr;
	//	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	//	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	//	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	//	hsArr = mxCreateDoubleMatrix(hs.size(), 1, mxREAL);
	//	hmArr = mxCreateDoubleMatrix(1, 1, mxREAL);
	//	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	//	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	//	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	//	memcpy(mxGetPr(hsArr), &hs[0], hs.size() * sizeof(double));
	//	memcpy(mxGetPr(hmArr), &hm1, sizeof(double));
	//	engPutVariable(eng, "row", rowArr);
	//	engPutVariable(eng, "col", colArr);
	//	engPutVariable(eng, "val", valArr);
	//	engPutVariable(eng, "hs", hsArr);
	//	engPutVariable(eng, "hm1", hmArr);

	//	str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());
	//	str = "n=length(hs)";        engEvalString(eng, str.c_str());
	//	str = "B = spdiags(hs(:),0,n,n);";        engEvalString(eng, str.c_str());
	//	str = "[V,D]=eigs(A, B, hm1, 'smallestabs')";        engEvalString(eng, str.c_str());

	//	mxArray* vv;
	//	vv = engGetVariable(eng, "V");
	//	double* v;
	//	v = reinterpret_cast<double*>(mxGetData(vv));
	//	const size_t* dims;
	//	dims = mxGetDimensions(vv);

	//	EigenVector ev; // no boundary
	//	ev.resize(dims[0]);

	//	std::vector<EigenVector> eigenfields;
	//	for (int j = 0; j < dims[1]; ++j) {
	//		for (int i = 0; i < dims[0]; ++i) {
	//			ev(i) = v[i + j * dims[0]];
	//		}
	//		eigenfields.push_back(ev);
	//	}

	//	EigenVector noHarmonicComponentNB = noHarmonicComponent;
	//	boundaryToNoBoundary1Form(noHarmonicComponentNB);

	//	EigenVector b = r_hs2f * r_ed1f*noHarmonicComponentNB;

	//	std::set<int> anchors;
	//	selectAnchors(eigenfields, anchors);
	//	correctRankDeficiency(s2l2, b, anchors);

	//	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//	solver.compute(s2l2);

	//	EigenVector x = solver.solve(b); // potential 1-form for exact component
	//	EigenVector vp = x;
	//	noBoundaryToBoundary2Form(vp);
	//	convert2Form(vp, vectorPotential);

	//	coexactComponent = r_hs1f.cwiseInverse() *r_ed1f.transpose()*r_hs2f * x;
	//	noBoundaryToBoundary1Form(coexactComponent);
	//}
}

void LaplacianAnalyzer::computeExact()
{
	//EigenSpMat A = r_hs0f * _laplacian_0;

	//EigenVector noHarmonicComponentNB = noHarmonicComponent;
	//boundaryToNoBoundary1Form(noHarmonicComponentNB);

	//EigenVector b = r_ed0f.transpose() * r_hs1f * noHarmonicComponentNB;

	////// potential has 1 degree of freedom. Fix one point at potential 0.
	////for (int i = 0; i < A.outerSize(); ++i) {
	////	for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
	////		if (iter.row() == 0 && iter.col() == 0) {
	////			iter.valueRef() = 1;
	////		}
	////		else if (iter.row() == 0 || iter.col() == 0) {
	////			iter.valueRef() = 0;
	////		}
	////	}
	////}

	////b(0) = 0; // set fix potential root.

	//std::cout << "A and b are set." << std::endl;

	//// solve
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//solver.compute(A);

	//EigenVector x = solver.solve(b); 
	//computeScalarPotential(x);

	//exactComponent = r_ed0f*x;
	//noBoundaryToBoundary1Form(exactComponent);
}

void LaplacianAnalyzer::convertForms()
{
	//convert1Form(omega, omegaVF);
	//convert1Form(harmonicKnot, harmonicKnotVF);
	//convert1Form(fluxlessKnot, fluxlessKnotVF);
	//convert1Form(groundedGradient, groundedGradientVF);
	//convert1Form(harmonicGradient, harmonicGradientVF);
	//convert1Form(curlyGradient, curlyGradientVF);
	//convert1Form(curlyGradientVP, curlyGradientVPVF);

	//convert1Form(curlFreeGradient, groundedGradientVF);
	//convert1Form(remainingDivFreeField, fluxlessKnotVF);

	//cout << "@@@@@@@@@@@@" << endl;
	//targetVF = fluxlessKnotVF;

	convert1Form(gamma_1, gamma_1_VF);
	convert1Form(gamma_2, gamma_2_VF);
	convert1Form(gamma_3, gamma_3_VF);
	convert1Form(gamma_4, gamma_4_VF);
}

void LaplacianAnalyzer::computeArrows()
{
	//assignColorGradient();
	//assignColorStrength(targetVF, 0);

	//assignColorStrength(omegaVF, 0);
	//assignColorStrength(fluxlessKnotVF, 3);
	//assignColorStrength(groundedGradientVF, 4);


	//rescaleArrowsWithVF(omegaVF);
	//rescaleArrowsWithVF(fluxlessKnotVF);
	//rescaleArrowsWithVF(groundedGradientVF);
	
	//computeArrowsWithVF(omegaVF);

	assignColorStrength(gamma_1_VF, 1);
	assignColorStrength(gamma_1_VF, 2);
	assignColorStrength(gamma_1_VF, 3);
	assignColorStrength(gamma_1_VF, 4);

	rescaleArrowsWithVF(gamma_1_VF);
	rescaleArrowsWithVF(gamma_2_VF);
	rescaleArrowsWithVF(gamma_3_VF);
	rescaleArrowsWithVF(gamma_4_VF);

	computeArrowsWithVF(gamma_2_VF);
}

void LaplacianAnalyzer::setEigenfield_L0n(int idx, bool isCurl)
{
	EigenSpMat s0l0nb = m_mesh_3.HS0F_NB()*m_laplacian0_NB;

	// assemble matrix to be eigen decomposed
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s0l0nb.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s0l0nb, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS0F_NB().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS0F_NB(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 101);

	printEigenvalues(emptyEval, m_suffix + "_L0n");

	//cout << "EigenValues:" << endl;
	//for (int i = 0; i < emptyEval.size(); ++i) {
		//cout << emptyEval[i] << endl;
	//}

	EigenVector evn = eigenFields[0];
	EigenVector ev;
	ev.setZero(m_mesh_3.NumV());
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		if (m_mesh_3.VertexOnBoundary(vi))
			ev(i) = 0;
		else
			ev(i) = evn(m_mesh_3.VertexForward(i));
	}

	assignColorGradient(ev);
}

void LaplacianAnalyzer::setEigenfield_L0t(int idx, bool isCurl)
{
	EigenSpMat s0l0 = m_mesh_3.HS0F_B()*m_laplacian0_B;

	// assemble matrix to be eigen decomposed
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s0l0.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s0l0, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS0F_B().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS0F_B(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 101);

	printEigenvalues(emptyEval, m_suffix + "_L0t");

	assignColorGradient(eigenFields[1]);
}

void LaplacianAnalyzer::setEigenfield_L1n(int idx, bool isCurl)
{
	EigenSpMat s1l1nb = m_mesh_3.HS1F_NB()*m_laplacian1_NB;

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

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS1F_NB().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_NB(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 30);

	groupFields_L1n(eigenFields);

	std::vector<double> divEV;
	std::vector<double> curlEV;
	for (int i = 0; i < curlSize; ++i) {
		curlEV.push_back(emptyEval[curlFieldIdx[i]]);
	}
	for (int i = 0; i < divSize; ++i) {
		divEV.push_back(emptyEval[divFieldIdx[i]]);
	}
	printEigenvalues(divEV, m_suffix + "_L1n_div");
	printEigenvalues(curlEV, m_suffix + "_L1n_curl");

	if (!isCurl) {
		//omega = eigenFields[divFieldIdx[idx]];
		
		gamma_1 = eigenFields[divFieldIdx[0]];
		gamma_2 = eigenFields[divFieldIdx[4]];
		gamma_3 = eigenFields[divFieldIdx[4]];
		gamma_4 = eigenFields[divFieldIdx[6]];
	}
	else {
		//omega = eigenFields[curlFieldIdx[idx]];

		gamma_1 = eigenFields[curlFieldIdx[0]];
		gamma_2 = eigenFields[curlFieldIdx[2]];
		gamma_3 = eigenFields[curlFieldIdx[4]];
		gamma_4 = eigenFields[curlFieldIdx[6]];
	}

	//noBoundaryToBoundary1Form(omega);
	noBoundaryToBoundary1Form(gamma_1);
	noBoundaryToBoundary1Form(gamma_2);
	noBoundaryToBoundary1Form(gamma_3);
	noBoundaryToBoundary1Form(gamma_4);
}

void LaplacianAnalyzer::setEigenfield_L1t(int idx, bool isCurl)
{
	EigenSpMat s1l1 = m_mesh_3.HS1F_B()*m_laplacian1_B;

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

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS1F_B().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_B(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 50);

	groupFields_L1t(eigenFields);

	std::vector<double> divEV;
	std::vector<double> curlEV;
	for (int i = 0; i < curlSize; ++i) {
		curlEV.push_back(emptyEval[curlFieldIdx[i]]);
	}
	for (int i = 0; i < divSize; ++i) {
		divEV.push_back(emptyEval[divFieldIdx[i]]);
	}
	printEigenvalues(divEV, m_suffix + "_L1t_div");
	printEigenvalues(curlEV, m_suffix + "_L1t_curl");

	if (!isCurl) {
		//omega = eigenFields[divFieldIdx[idx]];
		gamma_1 = eigenFields[divFieldIdx[0]];
		gamma_2 = eigenFields[divFieldIdx[3]];
		gamma_3 = eigenFields[divFieldIdx[6]];
		gamma_4 = eigenFields[divFieldIdx[9]];
	}
	else {
		gamma_1 = eigenFields[curlFieldIdx[0]];
		gamma_2 = eigenFields[curlFieldIdx[3]];
		gamma_3 = eigenFields[curlFieldIdx[6]];
		gamma_4 = eigenFields[curlFieldIdx[9]];

		/*gamma_1 = eigenFields[0];
		gamma_2 = eigenFields[1];
		gamma_3 = eigenFields[2];
		gamma_4 = eigenFields[3];*/
	}
}

void LaplacianAnalyzer::setEigenfield_L2n(int idx, bool isCurl)
{
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB()*m_laplacian2_NB;

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

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS2F_NB().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_NB(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 50);

	groupFields_L2n(eigenFields);

	std::vector<double> divEV;
	std::vector<double> curlEV;
	for (int i = 0; i < curlSize; ++i) {
		curlEV.push_back(emptyEval[curlFieldIdx[i]]);
	}
	for (int i = 0; i < divSize; ++i) {
		divEV.push_back(emptyEval[divFieldIdx[i]]);
	}
	printEigenvalues(divEV, m_suffix + "_L2n_div");
	printEigenvalues(curlEV, m_suffix + "_L2n_curl");

	if (!isCurl) {
		gamma_1 = eigenFields[divFieldIdx[0]];
	}
	else {
		gamma_1 = eigenFields[curlFieldIdx[0]];
	}

	noBoundaryToBoundary2Form(gamma_1);
}

void LaplacianAnalyzer::setEigenfield_L2t(int idx, bool isCurl)
{
	EigenSpMat s2l2 = m_mesh_3.HS2F_B()*m_laplacian2_B;

	// assemble matrix to be eigen decomposed
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s2l2.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s2l2, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS2F_B().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_B(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 50);

	groupFields_L2t(eigenFields);

	std::vector<double> divEV;
	std::vector<double> curlEV;
	for (int i = 0; i < curlSize; ++i) {
		curlEV.push_back(emptyEval[curlFieldIdx[i]]);
	}
	for (int i = 0; i < divSize; ++i) {
		divEV.push_back(emptyEval[divFieldIdx[i]]);
	}
	printEigenvalues(divEV, m_suffix + "_L2t_div");
	printEigenvalues(curlEV, m_suffix + "_L2t_curl");

	if (!isCurl) {
		gamma_1 = eigenFields[divFieldIdx[0]];
	}
	else {
		gamma_1 = eigenFields[curlFieldIdx[0]];
	}
}

void LaplacianAnalyzer::setEigenfield_L3n(int idx, bool isCurl)
{
	EigenSpMat s3l3nb = m_mesh_3.HS3F()*m_laplacian3_NB;

	// assemble matrix to be eigen decomposed
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s3l3nb.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s3l3nb, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS3F().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS3F(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 50);

	printEigenvalues(emptyEval, m_suffix + "_L3n");
}

void LaplacianAnalyzer::setEigenfield_L3t(int idx, bool isCurl)
{
	EigenSpMat s3l3 = m_mesh_3.HS3F()*m_laplacian3_B;

	// assemble matrix to be eigen decomposed
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s3l3.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s3l3, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS3F().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS3F(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, 101);

	printEigenvalues(emptyEval, m_suffix + "_L3t");
}

void LaplacianAnalyzer::computeCrossSection()
{
	double diff = m_mesh_3.BoundingBoxMax().y() - m_mesh_3.BoundingBoxMin().y();
	double pl = 0.5*diff + m_mesh_3.BoundingBoxMin().y();
	double cutThd = 0; // 2

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
					if (fi->first->vertex(i)->point().point().x() < cutThd)
						side = true;
					else {
						side = false;
						break;
					}
				}
			}

			if (fi->first->subdomain_index() == 0) {
				if (m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point().x() > cutThd)
					side = false;

				validCell = fi->first->neighbor(fi->second);
			}
			else
			{
				if (fi->first->vertex(fi->second)->point().point().x() > cutThd)
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
							//colors.push_back(scalarFieldCF[fi->first->vertex(i)]);

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

					if (fi->first->vertex(i)->point().point().x() < cutThd)
						numTargetSide++;
				}
			}

			if (infinite)
				continue;

			bool onBoundary = false;
			Mesh_3_Cell_iterator boundaryCell;
			Mesh_3_Cell_iterator validCell; // for extracting vector field.
			if (numTargetSide == 3) {
				if (fi->first->vertex(fi->second)->point().point().x() > cutThd) {
					boundaryCell = fi->first;
					validCell = fi->first->neighbor(fi->second);
					onBoundary = !onBoundary;
				}

				if (m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point().x() > cutThd) {
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
							//colors.push_back(scalarFieldCF[fi->first->vertex(i)]);
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

void LaplacianAnalyzer::computeHarmonicKnot()
{
	//EigenSpMat s1l1 = m_mesh_3.HS1F_B() * m_laplacian1_B;

	EigenSpMat eyeV;
	eyeV.resize(m_mesh_3.NumV(), m_mesh_3.NumV());
	eyeV.setIdentity();
	EigenSpMat s1l1 = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()
		+ m_mesh_3.HS1F_B() * m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();;

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
		//std::vector<double> hs;
		/*for (int i = 0; i < m_mesh_3.NumE(); ++i) {
			hs.push_back(m_mesh_3.HS1F_B().coeff(i, i));
		}*/
		std::vector<double> brow, bcol;
		std::vector<double> bval;
		for (int i = 0; i < m_mesh_3.HS1F_B().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_B(), i); iter; ++iter) {
				brow.push_back(static_cast<double>(iter.row()) + 1);
				bcol.push_back(static_cast<double>(iter.col()) + 1);
				bval.push_back(static_cast<double>(iter.value()));
			}
		}

		std::vector<double> emptyEval;
		//std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, harmonicKnotBasis, row, col, val, brow, bcol, bval, hm1);

		//harmonicKnot = harmonicKnotBasis[0];

		// we have harmonic basis
		harmonicKnot.setZero(m_mesh_3.NumE());
		for (int i = 0; i < harmonicKnotBasis.size(); ++i) {
			double c = (omega.transpose()*m_mesh_3.HS1F_B()*harmonicKnotBasis[i]).coeff(0,0)
				/ (harmonicKnotBasis[i].transpose() * m_mesh_3.HS1F_B()*harmonicKnotBasis[i]).coeff(0, 0);
			harmonicKnotCoeffs.push_back(c);
			harmonicKnot += c * harmonicKnotBasis[i];
			cout << c << endl;
		}




		//EigenVector aa = m_mesh_3.ED1F_B()*harmonicKnotBasis[0];
		//EigenVector bb = m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B()*harmonicKnotBasis[0];

		//cout << aa.dot(aa) << " " << bb.dot(bb) << endl;
	}
}

void LaplacianAnalyzer::computeHarmonicKnotTest()
{
	EigenSpMat eyeV;
	eyeV.resize(m_mesh_3.NumV(), m_mesh_3.NumV());
	eyeV.setIdentity();
	EigenSpMat s1l1 = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()
		+ m_mesh_3.CHS1F_B() * m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B();;

	double hm1 = m_mesh_3.NumHMLG1();
	cout << "HM1: " << hm1 << endl;

	if (hm1 == 0) {
		harmonicKnotTest.setZero(m_mesh_3.NumE());
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
		//std::vector<double> hs;
		/*for (int i = 0; i < m_mesh_3.NumE(); ++i) {
			hs.push_back(m_mesh_3.HS1F_B().coeff(i, i));
		}*/
		std::vector<double> brow, bcol;
		std::vector<double> bval;
		for (int i = 0; i < m_mesh_3.CHS1F_B().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.CHS1F_B(), i); iter; ++iter) {
				brow.push_back(static_cast<double>(iter.row()) + 1);
				bcol.push_back(static_cast<double>(iter.col()) + 1);
				bval.push_back(static_cast<double>(iter.value()));
			}
		}

		std::vector<double> emptyEval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, hm1);

		//harmonicKnotTest = eigenFields[0];

		// we have harmonic basis
		harmonicKnotTest.setZero(m_mesh_3.NumE());
		for (int i = 0; i < eigenFields.size(); ++i) {
			double c = (omega.transpose()*m_mesh_3.CHS1F_B()*eigenFields[i]).coeff(0,0)
				/ (eigenFields[i].transpose() * m_mesh_3.CHS1F_B()*eigenFields[i]).coeff(0, 0);
			//harmonicKnotCoeffs.push_back(c);
			harmonicKnotTest += c * eigenFields[i];
			cout << c << endl;
		}

		//EigenVector aa = m_mesh_3.ED1F_B()*eigenFields[0];
		//EigenVector bb = m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B()*eigenFields[0];

		//cout << aa.dot(aa) << " " << bb.dot(bb) << endl;
		//cout << harmonicKnot.dot(m_mesh_3.HS1F_B()*harmonicKnot) << endl;
		//cout << harmonicKnotTest.dot(m_mesh_3.CHS1F_B()*harmonicKnotTest) << endl;
		EigenVector cc = harmonicKnotTest - harmonicKnot;
		cout << cc.dot(m_mesh_3.HS1F_B()*cc) << endl;
	}
}

void LaplacianAnalyzer::computeHarmonicGradient()
{
	EigenSpMat A = m_laplacian0_NB;
	EigenVector b;
	
	std::vector<EigenVector> tmpBasis;
	for (int i = 0; i < m_mesh_3.NumBoundary() - 1; ++i) {
		assignBoundaryVertexPotential(i);

		b.setZero(m_mesh_3.NumIV());
		for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
			ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
			if (!m_mesh_3.EdgeValid(ei))
				continue;
			if (m_mesh_3.EdgeOnBoundary(ei))
				continue;

			Mesh_3_Vertex_iterator vi1 = ei->first->vertex(ei->second);
			Mesh_3_Vertex_iterator vi2 = ei->first->vertex(ei->third);

			if (m_mesh_3.VertexOnBoundary(vi1)) {
				b(m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi2))) += boundaryVertexPotential[vi1]
					* (1 / m_mesh_3.HS0F_NB().coeff(m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi2)), m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi2))))
					*m_mesh_3.HS1F_NB().coeff(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)));
			}
			else if (m_mesh_3.VertexOnBoundary(vi2)) {
				b(m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi1))) += boundaryVertexPotential[vi2]
					* (1 / m_mesh_3.HS0F_NB().coeff(m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi1)), m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi1))))
					* m_mesh_3.HS1F_NB().coeff(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)));
			}
		}
		
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		EigenVector x = solver.solve(b);
		
		EigenVector p;
		p.setZero(m_mesh_3.NumV());
		for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
			vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
			if (!m_mesh_3.VertexValid(vi))
				continue;

			if (m_mesh_3.VertexOnBoundary(vi)) {
				p(m_mesh_3.VertexIdx(vi)) = boundaryVertexPotential[vi];
			}
			else {
				p(m_mesh_3.VertexIdx(vi)) = x(m_mesh_3.VertexForward(m_mesh_3.VertexIdx(vi)));
			}
			
		}
		
		harmonicGradientPotentialBasis.push_back(p);
		//cout<< p <<endl;
		//for (int j = 0; j < p.size(); ++j) {
		//	scalarPotential[m_mesh_3.IdxToVertex(j)] = p(j);
		//}

		//harmonicGradient = m_mesh_3.ED0F_B()*p;

		EigenVector hmnc = m_mesh_3.ED0F_B()*p;
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

	harmonicGradient.setZero(m_mesh_3.NumE());
	for (int i = 0; i < harmonicGradientBasis.size(); ++i) {
		double numerator = (omega.transpose()*m_mesh_3.HS1F_B()*harmonicGradientBasis[i]);
		double denominator = (harmonicGradientBasis[i].transpose() * m_mesh_3.HS1F_B()*harmonicGradientBasis[i]);
		double c = numerator / denominator;
		cout << numerator << " " << denominator << endl;
		harmonicGradientCoeffs.push_back(c);

		cout << c << endl;
		harmonicGradient += c * harmonicGradientBasis[i];
	}


	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//solver.compute(s0l0nb);
	//EigenVector x = solver.solve(b);

	//groundedGradient = m_mesh_3.ED0F_NB()*x;
	//noBoundaryToBoundary1Form(groundedGradient);
}

void LaplacianAnalyzer::computeFluxlessKnot()
{
	EigenSpMat s2l2 = m_mesh_3.HS2F_B() * m_laplacian2_B;

	double hm2 = m_mesh_3.NumHMLG2();

	if (hm2 == 0) {
		EigenVector b = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*omega;

		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
		
		auto start = std::chrono::high_resolution_clock::now();
		solver.compute(s2l2);
		EigenVector x = solver.solve(b); // potential 2-form for exact component
		auto diff = std::chrono::high_resolution_clock::now() - start;
		auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(diff);
		std::cout << "Time 1: " << t1.count() << std::endl;
										 
										 //EigenVector vp = x;
		//noBoundaryToBoundary2Form(vp);
		//convert2Form(vp, vectorPotential);

		fluxlessKnot = m_mesh_3.HS1F_B().cwiseInverse() *m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B() * x;
		//noBoundaryToBoundary1Form(coexactComponent);
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
		//std::vector<double> hs;
		//for (int i = 0; i < m_mesh_3.NumF(); ++i) {
			//hs.push_back(m_mesh_3.HS2F_B().coeff(i, i));
		//}
		std::vector<double> brow, bcol;
		std::vector<double> bval;
		for (int i = 0; i < m_mesh_3.HS2F_B().outerSize(); ++i) {
			for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_B(), i); iter; ++iter) {
				brow.push_back(static_cast<double>(iter.row()) + 1);
				bcol.push_back(static_cast<double>(iter.col()) + 1);
				bval.push_back(static_cast<double>(iter.value()));
			}
		}
		
		std::vector<double> emptyEval;
		std::vector<EigenVector> eigenFields;
		matlabEIGS(emptyEval, eigenFields, row, col, val, brow, bcol, bval, hm2);

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

void LaplacianAnalyzer::computeGroundedGradient()
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

void LaplacianAnalyzer::computeCurlyGradient()
{
	//curlyGradient = omega - harmonicKnot - harmonicGradient - fluxlessKnot - groundedGradient;

	EigenVector tmpCG = omega - harmonicKnot - harmonicGradient - fluxlessKnot - groundedGradient;
	EigenVector tmpCG_NB = tmpCG;
	//EigenVector crucial = m_mesh_3.ED1F_B()*tmpCG;

	EigenSpMat s2l2nb = m_laplacian2_NB;
	boundaryToNoBoundary1Form(tmpCG_NB);
	EigenVector b = m_mesh_3.ED1F_NB()*tmpCG_NB;

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(s2l2nb);
	EigenVector x = solver.solve(b);
	curlyGradientVP = x;

	//noBoundaryToBoundary2Form(curlyGradientVP);
	//curlyGradient = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*curlyGradientVP;

	curlyGradient = m_mesh_3.HS1F_NB().cwiseInverse()*m_mesh_3.ED1F_NB().transpose()*m_mesh_3.HS2F_NB()*curlyGradientVP;
	noBoundaryToBoundary1Form(curlyGradient);
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		if (m_mesh_3.EdgeOnBoundary(m_mesh_3.IdxToEdge(i)))
			curlyGradient(i) = tmpCG(i);
	}
}

void LaplacianAnalyzer::computeCurlFreeGradient()
{
	//EigenSpMat s0l0 = m_mesh_3.HS0F_B()*m_laplacian0_B;
	//EigenVector b = m_mesh_3.ED0F_B().transpose()* m_mesh_3.HS1F_B()*omega;

	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	//solver.compute(s0l0);
	//EigenVector x = solver.solve(b);

	//curlFreeGradient = m_mesh_3.ED0F_B()*x;

	//cout << m_mesh_3.CHS1F_B() << endl;

	EigenSpMat A = m_mesh_3.ED0F_B().transpose()*m_mesh_3.CHS1F_B()*m_mesh_3.ED0F_B();
	EigenVector b = m_mesh_3.ED0F_B().transpose()* m_mesh_3.CHS1F_B()*omega;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;

	solver.compute(A);
	EigenVector x = solver.solve(b);


	curlFreeGradient = m_mesh_3.ED0F_B()*x;

	//cout<< curlFreeGradient <<endl;
	cout<<"Curl free gradient done ..."<<endl;
}

void LaplacianAnalyzer::computeRemainingDivFreeFieldPotential()
{
	EigenSpMat eyeE;
	eyeE.resize(m_mesh_3.NumE(), m_mesh_3.NumE());
	eyeE.setIdentity();

	EigenSpMat A = m_mesh_3.ED2F_B().transpose()*m_mesh_3.HS3F()*m_mesh_3.ED2F_B()
		+ m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B();

	EigenVector b1 = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.CHS1F_B()*(omega - curlFreeGradient);
	EigenVector b2 = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*(omega - curlFreeGradient);
	EigenVector b3 = m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B()*m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.CHS1F_B()*omega;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	EigenVector x1 = solver.solve(b1);
	EigenVector x2 = solver.solve(b2);
	EigenVector x3 = solver.solve(b3);

	EigenSpMat AA = m_mesh_3.CHS1F_B();
	EigenVector bb1 = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x1;
	EigenVector bb2 = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x2;
	EigenVector bb3 = m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x3;


	auto start = std::chrono::high_resolution_clock::now();
	solver.compute(AA);
	EigenVector delta_beta_1 = solver.solve(bb1);
	auto diff = std::chrono::high_resolution_clock::now() - start;
	auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(diff);
	std::cout << "Time 2: " << t1.count() << std::endl;

	EigenVector delta_beta_2 = solver.solve(bb2);
	EigenVector delta_beta_3 = solver.solve(bb3);
	EigenVector delta_beta_4 = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x1;
	EigenVector delta_beta_5 = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x2;
	EigenVector delta_beta_6 = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*x3;
	remainingDivFreeField = delta_beta_1;


	double jj = sqrt(curlFreeGradient.dot(m_mesh_3.CHS1F_B()*curlFreeGradient));
	double l1 = sqrt(delta_beta_1.dot(m_mesh_3.CHS1F_B()*delta_beta_1));
	double l2 = sqrt(delta_beta_2.dot(m_mesh_3.CHS1F_B()*delta_beta_2));
	double l3 = sqrt(delta_beta_3.dot(m_mesh_3.CHS1F_B()*delta_beta_3));
	double l4 = sqrt(delta_beta_4.dot(m_mesh_3.CHS1F_B()*delta_beta_4));
	double l5 = sqrt(delta_beta_5.dot(m_mesh_3.CHS1F_B()*delta_beta_5));
	double l6 = sqrt(delta_beta_6.dot(m_mesh_3.CHS1F_B()*delta_beta_6));

	cout << "Checking orthogonality: " << endl;
	cout << "with beta 1: " << curlFreeGradient.dot(m_mesh_3.CHS1F_B()*delta_beta_1) / (jj*l1) << endl;
	cout << "with beta 2: " << curlFreeGradient.dot(m_mesh_3.CHS1F_B()*delta_beta_2) / (jj*l2) << endl;
	cout << "with beta 3: " << curlFreeGradient.dot(m_mesh_3.CHS1F_B()*delta_beta_3) / (jj*l3) << endl;
	cout << "with beta 4: " << curlFreeGradient.dot(m_mesh_3.CHS1F_B()*delta_beta_4) / (jj*l4) << endl;
	cout << "with beta 5: " << curlFreeGradient.dot(m_mesh_3.CHS1F_B()*delta_beta_5) / (jj*l5) << endl;
	cout << "with beta 6: " << curlFreeGradient.dot(m_mesh_3.CHS1F_B()*delta_beta_6) / (jj*l6) << endl;

	EigenVector vv1 = omega - curlFreeGradient - delta_beta_1;
	EigenVector vv2 = omega - curlFreeGradient - delta_beta_2;
	EigenVector vv3 = omega - curlFreeGradient - delta_beta_3;
	EigenVector vv4 = omega - curlFreeGradient - delta_beta_4;
	EigenVector vv5 = omega - curlFreeGradient - delta_beta_5;
	EigenVector vv6 = omega - curlFreeGradient - delta_beta_6;

	//cout<<"J :"<< curlFreeGradient.dot(m_mesh_3.CHS1F_B()*curlFreeGradient)<<endl;
	//cout << "Delta_Beta 1: " << delta_beta_1.dot(m_mesh_3.CHS1F_B()*delta_beta_1) << endl;
	//cout << "Delta_Beta 2: " << delta_beta_2.dot(m_mesh_3.CHS1F_B()*delta_beta_2) << endl;
	//cout << "Delta_Beta 3: " << delta_beta_3.dot(m_mesh_3.CHS1F_B()*delta_beta_3) << endl;

	//cout << "Checking absolute error: " << endl;
	//cout << "with beta 1: " << vv1.dot(m_mesh_3.CHS1F_B()*vv1)  << endl;
	//cout << "with beta 2: " << vv2.dot(m_mesh_3.CHS1F_B()*vv2)  << endl;
	//cout << "with beta 3: " << vv3.dot(m_mesh_3.CHS1F_B()*vv3)  << endl;

	cout<<"Checking relative error: "<<endl;
	cout << "with beta 1: " << vv1.dot(m_mesh_3.CHS1F_B()*vv1) / omega.dot(m_mesh_3.CHS1F_B()*omega) << endl;
	cout << "with beta 2: " << vv2.dot(m_mesh_3.CHS1F_B()*vv2) / omega.dot(m_mesh_3.CHS1F_B()*omega) << endl;
	cout << "with beta 3: " << vv3.dot(m_mesh_3.CHS1F_B()*vv3) / omega.dot(m_mesh_3.CHS1F_B()*omega) << endl;
	cout << "with beta 4: " << vv4.dot(m_mesh_3.CHS1F_B()*vv4) / omega.dot(m_mesh_3.CHS1F_B()*omega) << endl;
	cout << "with beta 5: " << vv5.dot(m_mesh_3.CHS1F_B()*vv5) / omega.dot(m_mesh_3.CHS1F_B()*omega) << endl;
	cout << "with beta 6: " << vv6.dot(m_mesh_3.CHS1F_B()*vv6) / omega.dot(m_mesh_3.CHS1F_B()*omega) << endl;
	//remainingDivFreeField = m_mesh_3.ED1F_B().transpose()*m_mesh_3.CHS2F_B()*x;

	curlFreeGradient -= 0.8*curlyGradient;
}

void LaplacianAnalyzer::testSPDHS()
{
	cout << "Start Testing ..." << endl;

	EigenSpMat A = m_mesh_3.CHS1F_B();

	for (int i = 0; i < A.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
			double r, c, v;
			r = iter.row();
			c = iter.col();
			v = iter.value();

			if (v != A.coeff(c, r)) {
				cout << v << " " << A.coeff(c, r) << endl;
			}
		}
	}

	Eigen::SimplicialLLT<EigenSpMat> llt(A);
	if (llt.info() == Eigen::NumericalIssue)
	{
		throw std::runtime_error("Possibly non semi-positive definitie matrix!");
	}

	cout<< A <<endl;
}

void LaplacianAnalyzer::testHodgeStar1()
{

}

void LaplacianAnalyzer::assignBoundaryVertexPotential(int idx)
{
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (m_mesh_3.VertexOnBoundary(vi)) {
			if (m_mesh_3.VertexBoundaryIdx(vi) == idx) {
				boundaryVertexPotential[vi] = 1;
			}
			else {
				boundaryVertexPotential[vi] = 0;
			}
		}
	}
}

void LaplacianAnalyzer::matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
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

	str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());
	//str = "n=length(hs)";        engEvalString(eng, str.c_str());
	//str = "B = spdiags(hs(:),0,n,n);";        engEvalString(eng, str.c_str());
	str = "B=sparse(brow, bcol, bval)";        engEvalString(eng, str.c_str());
	str = "[V,D]=eigs(A, B, hm, 'smallestabs')";        engEvalString(eng, str.c_str());

	// read out eigenvector matrix V from matlab
	mxArray* vv;
	vv = engGetVariable(eng, "V");
	double* v;
	v = reinterpret_cast<double*>(mxGetData(vv));
	const size_t* dims;
	dims = mxGetDimensions(vv);

	// read out eigenvector matrix V from matlab
	mxArray* dd;
	dd = engGetVariable(eng, "D");
	double* d;
	d = reinterpret_cast<double*>(mxGetData(dd));

	EigenVector ev; // no boundary
	ev.resize(dims[0]);

	evec.clear();
	for (int j = 0; j < dims[1]; ++j) {
		for (int i = 0; i < dims[0]; ++i) {
			ev(i) = v[i + j * dims[0]];
		}
		evec.push_back(ev);
	}

	eval.clear();
	for (int i = 0; i < dims[1]; ++i) {
		eval.push_back(d[i + i * dims[1]]);
	}
}

void LaplacianAnalyzer::integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf,
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

void LaplacianAnalyzer::prepareEigenFeilds(std::vector<EigenVector>& eigenFields)
{
	std::vector<EigenVector> eigenfields2Form;
	computeCurlFields2Form(eigenfields2Form);
	groupFields_L2n(eigenfields2Form);
	cout<< "EigenFields Num: "<<eigenfields2Form.size()<<endl;
	convert2FormFieldTo1FormField(eigenFields, eigenfields2Form);
}

void LaplacianAnalyzer::computeCurlFields2Form(std::vector<EigenVector>& eigenFields)
{
	EigenSpMat s2l2nb = m_mesh_3.HS2F_NB() * m_laplacian2_NB;

	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s2l2nb.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s2l2nb, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS2F_NB().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS2F_NB(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> eigenVals;
	matlabEIGS(eigenVals, eigenFields, row, col, val, brow, bcol, bval, 50);
}

void LaplacianAnalyzer::groupFields_L1n(std::vector<EigenVector> eigenFields)
{
	int curlIdx = 0;
	int divIdx = 0;

	int start = m_mesh_3.NumHMLG2();

	for (int i = start; i < eigenFields.size(); ++i) {
		cout << eigenFields[i].dot(m_curlEM_L1n*eigenFields[i]) << " " << eigenFields[i].dot(m_divEM_L1n*eigenFields[i]) << endl;
		if (eigenFields[i].dot(m_curlEM_L1n*eigenFields[i]) > eigenFields[i].dot(m_divEM_L1n*eigenFields[i])) {
			curlFieldIdx[curlIdx] = i;
			++curlIdx;
		}
		else {
			divFieldIdx[divIdx] = i;
			++divIdx;
		}
	}

	curlSize = curlIdx;
	divSize = divIdx;

	cout << "Curl Fields: " << curlIdx << endl;
	cout << "Div Fields: " << divIdx << endl;
}

void LaplacianAnalyzer::groupFields_L1t(std::vector<EigenVector> eigenFields)
{
	int curlIdx = 0;
	int divIdx = 0;

	int start = m_mesh_3.NumHMLG1();

	for (int i = start; i < eigenFields.size(); ++i) {
		cout << eigenFields[i].dot(m_curlEM_L1t*eigenFields[i]) << " " << eigenFields[i].dot(m_divEM_L1t*eigenFields[i]) << endl;
		if (eigenFields[i].dot(m_curlEM_L1t*eigenFields[i]) > eigenFields[i].dot(m_divEM_L1t*eigenFields[i])) {
			curlFieldIdx[curlIdx] = i;
			++curlIdx;
		}
		else {
			divFieldIdx[divIdx] = i;
			++divIdx;
		}
	}

	curlSize = curlIdx;
	divSize = divIdx;

	cout << "Curl Fields: " << curlIdx << endl;
	cout << "Div Fields: " << divIdx << endl;
}

void LaplacianAnalyzer::groupFields_L2n(std::vector<EigenVector> eigenFields)
{
	int curlIdx = 0;
	int divIdx = 0;

	for (int i = 0; i < eigenFields.size(); ++i) {
		cout << eigenFields[i].dot(m_curlEM_L2n*eigenFields[i]) << " " << eigenFields[i].dot(m_divEM_L2n*eigenFields[i]) << endl;
		if (eigenFields[i].dot(m_curlEM_L2n*eigenFields[i]) > eigenFields[i].dot(m_divEM_L2n*eigenFields[i])) {
			curlFieldIdx[curlIdx] = i;
			++curlIdx;
		}
		else {
			divFieldIdx[divIdx] = i;
			++divIdx;
		}
	}

	curlSize = curlIdx;
	divSize = divIdx;

	cout << "Curl Fields: " << curlIdx << endl;
	cout << "Div Fields: " << divIdx << endl;
}

void LaplacianAnalyzer::groupFields_L2t(std::vector<EigenVector> eigenFields)
{
	int curlIdx = 0;
	int divIdx = 0;

	for (int i = 0; i < eigenFields.size(); ++i) {
		cout << eigenFields[i].dot(m_curlEM_L2t*eigenFields[i]) << " " << eigenFields[i].dot(m_divEM_L2t*eigenFields[i]) << endl;
		if (eigenFields[i].dot(m_curlEM_L2t*eigenFields[i]) > eigenFields[i].dot(m_divEM_L2t*eigenFields[i])) {
			curlFieldIdx[curlIdx] = i;
			++curlIdx;
		}
		else {
			divFieldIdx[divIdx] = i;
			++divIdx;
		}
	}

	curlSize = curlIdx;
	divSize = divIdx;

	cout << "Curl Fields: " << curlIdx << endl;
	cout << "Div Fields: " << divIdx << endl;
}

void LaplacianAnalyzer::printEigenvalues(std::vector<double> eigenvalues, string sfx)
{
	std::string fn = "Eigenvalues_" + sfx + ".txt";

	std::ofstream file(fn.c_str());

	for (int i = 0; i < eigenvalues.size(); ++i) {
		file << eigenvalues[i] << endl;
	}

	file.close();
}

void LaplacianAnalyzer::convert2FormFieldTo1FormField(std::vector<EigenVector>& eigenfields1Form, std::vector<EigenVector> eigenfields2Form)
{
	std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf;
	for (int i = 0; i < eigenfields2Form.size(); ++i) {
		noBoundaryToBoundary2Form(eigenfields2Form[i]);
		//cout << eigenfields2Form[i].size() << endl;
		convert2Form(eigenfields2Form[i], vf);
	
		EigenVector field;
		field.resize(m_mesh_3.NumE());
		field.setZero();
		for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
			ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
			if (!m_mesh_3.EdgeValid(ei))
				continue;

			Mesh_3_Cell_circulator cir = m_mesh_3.C3T3().triangulation().incident_cells(*ei);
			Mesh_3_Cell_circulator end = cir;

			Mesh_3_Vector_3 v(0, 0, 0);
			int num = 0;
			do {
				if (cir->subdomain_index() != 0) {
					v += vf[cir];
					++num;
				}
				++cir;
			} while (cir != end);
			
			v /= num;

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

			field(m_mesh_3.EdgeIdx(ei)) = l * v;
		}

		eigenfields1Form.push_back(field);
		//cout<< field <<endl;
	}
}

void LaplacianAnalyzer::setConstantField(EigenVector& field, Mesh_3_Vector_3 d)
{
	field.setZero(m_mesh_3.NumE());
	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;

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

		field(m_mesh_3.EdgeIdx(ei)) = d * l;
	}

	//field.setZero(m_mesh_3.NumF());
	//for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
	//	fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
	//	if (!m_mesh_3.FacetValid(fi))
	//		continue;

	//	std::vector<int> vIdx;
	//	for (int i = 0; i < 4; ++i) {
	//		if (i != fi->second) {
	//			vIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i))); // get global vertex idx
	//		}
	//	}
	//	std::sort(vIdx.begin(), vIdx.end()); // set consistent orientation, increasing order of vertex indices

	//	Mesh_3_Vector_3 v1, v2;
	//	v1 = m_mesh_3.IdxToVertex(vIdx[1])->point().point() - m_mesh_3.IdxToVertex(vIdx[0])->point().point();
	//	v2 = m_mesh_3.IdxToVertex(vIdx[2])->point().point() - m_mesh_3.IdxToVertex(vIdx[1])->point().point();
	//	Mesh_3_Vector_3 av = CGAL::cross_product(v1, v2) / 2; // area normal vector

	//	field(m_mesh_3.FacetIdx(fi)) = av * d;
	//}
}

void LaplacianAnalyzer::setPointChargeElectricField(EigenVector& field, Mesh_3_Point_3 rp, bool ori)
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

		//if (sqrt(d*d) > 0.3*m_mesh_3.getAveLen())
			//d = 0.3*m_mesh_3.getAveLen()*(d / sqrt(d*d));

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
}

void LaplacianAnalyzer::setCurrentMagneticField(EigenVector& field, Mesh_3_Point_3 rp, Mesh_3_Vector_3 d, bool ori)
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
}

void LaplacianAnalyzer::smooth(EigenVector& form)
{
	//EigenVector formNB = form;
	//boundaryToNoBoundary1Form(formNB);

	//EigenSpMat s1l1 = m_mesh_3.HS1F_NB() * m_laplacian_1;

	//double delta = 0.01; // smoothing time step
	//EigenVector current = formNB; // form in current time step;
	//EigenSpMat op = m_mesh_3.HS1F_NB() + delta * s1l1; // this is an symmetric matrix
	//EigenVector next; // form after 1 step of smoothing
	//Eigen::SimplicialLDLT<EigenSpMat> solver;
	//solver.compute(op);
	//int numIter = 1;
	//int i = 0;
	//while (i < numIter) {
	//	current = m_mesh_3.HS1F_NB() * current;
	//	next = solver.solve(current);
	//	current = next;

	//	++i;
	//}

	//form = current;
	//noBoundaryToBoundary1Form(form);
}

void LaplacianAnalyzer::boundaryToNoBoundary1Form(EigenVector& form)
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

void LaplacianAnalyzer::noBoundaryToBoundary1Form(EigenVector& form)
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

void LaplacianAnalyzer::boundaryToNoBoundary2Form(EigenVector& form)
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

void LaplacianAnalyzer::noBoundaryToBoundary2Form(EigenVector& form)
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

void LaplacianAnalyzer::computeScalarPotential(EigenVector form)
{
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		scalarPotential[vi] = form(m_mesh_3.VertexIdx(vi));
	}
}

void LaplacianAnalyzer::convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

void LaplacianAnalyzer::convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

void LaplacianAnalyzer::computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	vfvNum = 0;

	int i = 0;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() != 0) {
			if (i % 10 == 0) {

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

void LaplacianAnalyzer::rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
{
	// set the max length to 0.5*aveLen 
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
	//pivotLen = lenVec[static_cast<int>(lenVec.size() - 1)];
	double aveLen = m_mesh_3.getAveLen();
	if (pivotLen < aveLen)
		pivotLen = aveLen;

	
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Vector_3 v = vf[ci];
		double len = sqrt(v*v);
		v /= len;

		if (len > pivotLen)
			len = pivotLen;
		if (len < 0.4*pivotLen)
			len = 0.4*pivotLen;

		len = (0.6*aveLen)*(len / (pivotLen));

		vf[ci] = v * len;
	}
}

void LaplacianAnalyzer::assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci)
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

void LaplacianAnalyzer::selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors)
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

void LaplacianAnalyzer::correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors)
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

void LaplacianAnalyzer::assignColorStrength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode)
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

		if (mode == 1) {
			gamma_1_CF[vi] = col;
		}
		else if (mode == 2) {
			gamma_2_CF[vi] = col;
		}
		else if (mode == 3) {
			gamma_3_CF[vi] = col;
		}
		else if (mode == 4) {
			gamma_4_CF[vi] = col;
		}
	}

	if (mode == 1) {
		colorField = gamma_1_CF;
	}
	else if (mode == 2) {
		colorField = gamma_2_CF;
	}
	else if (mode == 3) {
		colorField = gamma_3_CF;
	}
	else if (mode == 4) {
		colorField = gamma_4_CF;
	}


}

void LaplacianAnalyzer::assignColorGradient(EigenVector ev)
{
	//std::vector<double> val;
	//for (int i = 0; i < ev.size(); ++i) {
	//	val.push_back(ev(i));
	//}

	//std::vector<int> sIndices;
	//sortedIndices(sIndices, val);

	////defien color range
	//Mesh_3_Vector_3 colorMin(0.2f, 0.4f, 0.8f);
	//Mesh_3_Vector_3 colorMax(0.8f, 0.4f, 0.2f);

	//for (int i = 0; i < sIndices.size(); ++i) {
	//	Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(sIndices[i]);
	//	double ratio = i * 1.0 / sIndices.size();

	//	scalarFieldCF[vi] = assignColorSF(ratio);
	//}

	// assign color 
	double pmin = ev(0);
	double pmax = ev(0);
	for (int i = 1; i < m_mesh_3.NumV(); ++i) {
		if (ev(i) > pmax)
			pmax = ev(i);
		if (ev(i) < pmin)
			pmin = ev(i);
	}

	double diff = pmax - pmin;
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		double val = ev(i);
		double d = val - (pmin + 0.0*diff);

		double ratio = d / (diff*1.0);
		if (ratio < 0)
			ratio = 0;
		if (ratio > 1)
			ratio = 1;

		scalarFieldCF[vi] = assignColorSF(ratio);
	}

	colorField = groundedGradientCF;
}

Mesh_3_Vector_3 LaplacianAnalyzer::assignColor(double ratio)
{
	Mesh_3_Vector_3 color;

	//Mesh_3_Vector_3 red(0.9, 0.7, 0.7);
	//Mesh_3_Vector_3 yellow(0.9, 0.9, 0.7);
	//Mesh_3_Vector_3 green(0.7, 0.9, 0.7);
	//Mesh_3_Vector_3 indigo(0.7, 0.9, 0.9);
	//Mesh_3_Vector_3 blue(0.7, 0.7, 0.9);

	Mesh_3_Vector_3 red(0.8, 0.3, 0.6);
	Mesh_3_Vector_3 white(0.7, 0.7, 0.7);
	Mesh_3_Vector_3 blue(0.2, 0.6, 0.8);

	if (ratio < 0.6) {
		double r = ratio / 0.6;
		color = (1 - r)*blue + r * white;
	}
	else if (ratio < 0.7) {
		color = white;
	}
	else{
		double r = (ratio - 0.7) / 0.3;
		color = (1 - r)*white + r * red;
	}

	return color;
}

Mesh_3_Vector_3 LaplacianAnalyzer::assignColorSF(double ratio)
{
	Mesh_3_Vector_3 color;

	Mesh_3_Vector_3 blue(0.2, 0.7, 0.9);
	Mesh_3_Vector_3 white(0.9, 0.9, 0.9);
	Mesh_3_Vector_3 red(0.9, 0.5, 0.2);

	if (ratio < 0.55) {
		double r = ratio / 0.55;
		color = (1 - r)*blue + r * white;
	}
	else if (ratio < 0.7) {
		color = white;
	}
	else {
		double r = (ratio - 0.7) / 0.3;
		color = (1 - r)*white + r * red;
	}

	return color;
}

double LaplacianAnalyzer::medianLength(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	std::vector<double> vl;
	for (auto iter = vf.begin(); iter != vf.end(); ++iter) {
		vl.push_back(sqrt(iter->second*iter->second));
	}
	std::nth_element(vl.begin(), vl.begin() + vl.size() / 2, vl.end());

	return vl[vl.size() / 2];
}

void LaplacianAnalyzer::sortedIndices(std::vector<int>& sIndex, std::vector<double> len)
{
	sIndex.resize(len.size(), 0);
	std::iota(sIndex.begin(), sIndex.end(), 0);

	std::sort(sIndex.begin(), sIndex.end(), [&len](int i1, int i2) {return len[i1] < len[i2]; });
}

void LaplacianAnalyzer::assignBoundaryVerticesGroup()
{
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (vi->point().point().y() > 0.35 - 0.001)
			vertexGroups[vi] = 1;
		else if(vi->point().point().y() < -0.35 + 0.001)
			vertexGroups[vi] = 2;
		else
			vertexGroups[vi] = 0;
	}
}

void LaplacianAnalyzer::assignVertexColorMixBoundary()
{
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;
		
		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] == 1) {
			colorField[vi] = Mesh_3_Vector_3(0.4, 0.4, 0.8);
		}
		else if (vertexGroups[vi] == 2) {
			colorField[vi] = Mesh_3_Vector_3(0.8, 0.4, 0.4);
		}
		else{
			colorField[vi] = Mesh_3_Vector_3(0.6, 0.6, 0.6);
		}
	}
}

void LaplacianAnalyzer::computeHarmonicMixBoundary()
{
	std::map<Mesh_3_Vertex_iterator, double> boundaryPotential;
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] == 1) {
			boundaryPotential[vi] = 0;
		}
		else if (vertexGroups[vi] == 2) {
			boundaryPotential[vi] = 1;
		}

	}

	EigenSpMat A = m_laplacian0_B;
	EigenVector b;
	b.setZero(m_mesh_3.NumV());

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] != 0) {
			EigenVector help;
			help.setZero(m_mesh_3.NumV());
			help(m_mesh_3.VertexIdx(vi)) = 1;

			b -= boundaryPotential[vi] * (A * help);
		}
	}

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] != 0) {
			b(m_mesh_3.VertexIdx(vi)) = boundaryPotential[vi];
		}
	}

	for (int i = 0; i < A.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
			if (iter.row() == iter.col()) {
				if (!m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(iter.row())))
					continue;

				if (vertexGroups[m_mesh_3.IdxToVertex(iter.row())] != 0) {
					iter.valueRef() = 1;
				}
			}
			else {
				if (!m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(iter.row()))
					&& !m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(iter.col())))
					continue;

				if (vertexGroups[m_mesh_3.IdxToVertex(iter.row())] != 0 
					|| vertexGroups[m_mesh_3.IdxToVertex(iter.col())] != 0) {
					iter.valueRef() = 0;
				}
			}
		}
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	EigenVector x = solver.solve(b);

	harmonicGradient = m_mesh_3.ED0F_B()*x;

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		scalarPotential[vi] = x(m_mesh_3.VertexIdx(vi));
	}

	//assignColorGradient();
}

void LaplacianAnalyzer::computeExactMixBoundary()
{
	std::map<Mesh_3_Vertex_iterator, double> boundaryPotential;
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] == 1) {
			boundaryPotential[vi] = 0;
		}
		else if (vertexGroups[vi] == 2) {
			boundaryPotential[vi] = 0;
		}

	}

	EigenSpMat A = m_laplacian0_B;
	EigenVector b;
	b = m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B()*omega;

	for (int i = 0; i < b.size(); ++i) {
	
		if (m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(i))) {
			if (vertexGroups[m_mesh_3.IdxToVertex(i)] == 0)
				b(i) = 0;
		}
	}

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] != 0) {
			EigenVector help;
			help.setZero(m_mesh_3.NumV());
			help(m_mesh_3.VertexIdx(vi)) = 1;

			b -= boundaryPotential[vi] * (A * help);
		}
	}

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] != 0) {
			b(m_mesh_3.VertexIdx(vi)) = boundaryPotential[vi];
		}
	}

	for (int i = 0; i < A.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(A, i); iter; ++iter) {
			if (iter.row() == iter.col()) {
				if (!m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(iter.row())))
					continue;

				if (vertexGroups[m_mesh_3.IdxToVertex(iter.row())] != 0) {
					iter.valueRef() = 1;
				}
			}
			else {
				if (!m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(iter.row()))
					&& !m_mesh_3.VertexOnBoundary(m_mesh_3.IdxToVertex(iter.col())))
					continue;

				if (vertexGroups[m_mesh_3.IdxToVertex(iter.row())] != 0
					|| vertexGroups[m_mesh_3.IdxToVertex(iter.col())] != 0) {
					iter.valueRef() = 0;
				}
			}
		}
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	EigenVector x = solver.solve(b);

	groundedGradient = m_mesh_3.ED0F_B()*x;

	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		scalarPotential[vi] = x(m_mesh_3.VertexIdx(vi));
	}
}

void LaplacianAnalyzer::computeCoexactMixBoundary()
{
	EigenSpMat s1l1 = m_mesh_3.HS1F_B() * m_laplacian1_B;

	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s1l1.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> brow, bcol;
	std::vector<double> bval;
	for (int i = 0; i < m_mesh_3.HS1F_B().outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_mesh_3.HS1F_B(), i); iter; ++iter) {
			brow.push_back(static_cast<double>(iter.row()) + 1);
			bcol.push_back(static_cast<double>(iter.col()) + 1);
			bval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> eigenVals;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(eigenVals, eigenFields, row, col, val, brow, bcol, bval, 10);

	fluxlessKnot = eigenFields[3];


}

void LaplacianAnalyzer::assembleLaplacian1MixBoundary()
{
	EigenSpMat ed0f_m = m_mesh_3.ED0F_B();
	EigenSpMat helper_ed0f;
	helper_ed0f.resize(m_mesh_3.NumV(), m_mesh_3.NumV());
	for (Mesh_3_Vertex_iterator vi = m_mesh_3.C3T3().triangulation().vertices_begin();
		vi != m_mesh_3.C3T3().triangulation().vertices_end(); ++vi) {
		if (!m_mesh_3.VertexValid(vi))
			continue;

		if (!m_mesh_3.VertexOnBoundary(vi))
			continue;

		if (vertexGroups[vi] != 0){
			helper_ed0f.coeffRef(m_mesh_3.VertexIdx(vi), m_mesh_3.VertexIdx(vi)) = 1;
		}
	}
}

void LaplacianAnalyzer::writeMesh()
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

void LaplacianAnalyzer::writeMeshPLY()
{
	std::ofstream out("mesh.ply");

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << vertices.size() << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "element face " << indices.size() / 3 << endl;
	out << "property list uchar int vertex_index" << endl;
	out << "end_header" << endl;

	for (int i = 0; i < vertices.size(); ++i) {
		out << vertices[i].x() << " "
			<< vertices[i].y() << " "
			<< vertices[i].z() << " "
			<< static_cast<int>(round(colors[i].x() * 255)) << " "
			<< static_cast<int>(round(colors[i].y() * 255)) << " "
			<< static_cast<int>(round(colors[i].z() * 255)) << endl;
	}

	for (int i = 0; i < indices.size(); i += 3) {
		out << "3 "
			<< indices[i] << " "
			<< indices[i + 1] << " "
			<< indices[i + 2] << endl;
	}

	out.close();
}

// functions for Poelke method
void LaplacianAnalyzer::solveCurl_Poelke()
{

}

void LaplacianAnalyzer::solveGradient_Poelke()
{

}

void LaplacianAnalyzer::makeCurlOperator_Poelke(EigenSpMat& curlOp)
{
	// Curl operator is a span of basis
	// It should be an 3|C| by |IE| matrix mapping 1-form to PWC vector field in each tet.
	// Each column should be a basis of curl space.

	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;
		if (m_mesh_3.EdgeOnBoundary(ei))
			continue;
		//cout << "Edge: " << m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)) << endl;
		// iterate through all the adjacent tets.
		// for each tet, assign 2 form value for 4 facets.
		// compute vector at center of tet

		Mesh_3_Cell_circulator ccir = m_mesh_3.C3T3().triangulation().incident_cells(*ei);
		Mesh_3_Cell_circulator end = ccir;

		do {
			std::map<int, Mesh_3_Vector_3> grad;
			for (int i = 0; i < 4; ++i) {
				My_facet mf(ccir, i);

				// magnitude
				double h = 3 * m_mesh_3.CellPrimal(ccir) / m_mesh_3.FacetPrimal(m_mesh_3.MyCellFacetMap(mf));

				// direction normal
				std::vector<Mesh_3_Point_3> pv;
				for (int j = 0; j < 4; ++j) {
					if (j != i) {
						pv.push_back(ccir->vertex(j)->point().point());
					}
				}

				Mesh_3_Vector_3 v1 = pv[1] - pv[0];
				Mesh_3_Vector_3 v2 = pv[2] - pv[0];

				Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
				n /= sqrt(n*n);

				Mesh_3_Vector_3 vc = pv[0] - ccir->vertex(i)->point().point();
				if (vc*n < 0)
					n = -n;

				grad[m_mesh_3.VertexIdx(ccir->vertex(i))] = n / h;
				//std::cout<< grad[vertexIdx[ci->vertex(i)]] <<std::endl;
			}

			// get 4 facet 2-form value after curl (d_1)
			std::vector<int> facetVal;
			for (int i = 0; i < 4; ++i) {
				My_facet mf(ccir, i);
				Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);
				double val = 0;

				EigenVector colv;
				colv.setZero(m_mesh_3.NumIE());
				colv(m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei))) = 1;
				EigenMatrix rowv = m_mesh_3.ED1F_NB().row(m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)));
				val = (rowv *colv)(0,0);

				facetVal.push_back(val);
				//cout<<"Facet Val: "<<val<<endl;
			}

			// get vector for this tet
			Mesh_3_Vector_3 rv(0, 0, 0); // recovered vector
			Mesh_3_Vector_3 rvt(0, 0, 0);
			std::vector<double> bary{ 0.1,0.2,0.3,0.4 };
			for (int i = 0; i < 4; ++i) { //iterate through 4 faces;
				std::vector<int> vIdx;
				std::vector<double> br;
				for (int j = 0; j < 4; ++j) {
					if (j != i) {
						vIdx.push_back(m_mesh_3.VertexIdx(ccir->vertex(j)));
						br.push_back(bary[j]);
					}
				}

				//std::sort(vIdx.begin(), vIdx.end());
				for (int k = 0; k < vIdx.size(); ++k) {
					for (int l = k + 1; l < vIdx.size(); ++l) {
						if (vIdx[k] > vIdx[l]) {
							std::swap(vIdx[k], vIdx[l]);
							std::swap(br[k], br[l]);
						}
					}
				}
				My_facet mf(ccir, i);

				rv += facetVal[i] * 2 * 0.25
					* (CGAL::cross_product(grad[vIdx[0]], grad[vIdx[1]])
						+ CGAL::cross_product(grad[vIdx[1]], grad[vIdx[2]])
						+ CGAL::cross_product(grad[vIdx[2]], grad[vIdx[0]]));
				
				/*rvt += facetVal[i] * 2
					* (CGAL::cross_product(grad[vIdx[0]], grad[vIdx[1]])*br[2]
						+ CGAL::cross_product(grad[vIdx[1]], grad[vIdx[2]])*br[0]
						+ CGAL::cross_product(grad[vIdx[2]], grad[vIdx[0]])*br[1]);*/
			}

			//Mesh_3_Vector_3 diff = rv - rvt;
			//cout<<rv<<endl;
			//cout << rvt << endl;
			//cout << "Difference :" << diff * diff << endl;

			//std::vector<int> opidx;
			//for (int i = 0; i < 4; ++i) {
			//	if (i != ei->second && i != ei->third)
			//		opidx.push_back(i);
			//}
			//Mesh_3_Point_3 p1 = ei->first->vertex(opidx[0])->point().point();
			//Mesh_3_Point_3 p2 = ei->first->vertex(opidx[0])->point().point();

			//cout << CGAL::cross_product((p1 - p2), rv) << endl;

			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir), m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), rv.x()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir) + 1, m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), rv.y()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir) + 2, m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), rv.z()));

			++ccir;
		} while (ccir != end);
	}

	curlOp.resize(3 * m_mesh_3.NumC(), m_mesh_3.NumIE());
	curlOp.setFromTriplets(triplet.begin(), triplet.end());
}
void LaplacianAnalyzer::makeGradOperator_Poelke(EigenSpMat& gradOp)
{
	// Div operator is a span of basis
	// It should be an 3|C| by |IF| matrix mapping 1-form to PWC vector field in each tet.
	// Each column should be a basis of curl space.

	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;
		if (m_mesh_3.FacetOnBoundary(fi))
			continue;

		Mesh_3_Cell_iterator ci1 = fi->first;
		Mesh_3_Cell_iterator ci2 = fi->first->neighbor(fi->second);

		double mag1 = m_mesh_3.FacetPrimal(fi) / m_mesh_3.CellPrimal(ci1);
		double mag2 = m_mesh_3.FacetPrimal(fi) / m_mesh_3.CellPrimal(ci2);
	
		// direction normal
		std::vector<Mesh_3_Point_3> pv;
		for (int i = 0; i < 4; ++i) {
			if (i != fi->second) {
				pv.push_back(ci1->vertex(i)->point().point());
			}
		}

		Mesh_3_Vector_3 v1 = pv[1] - pv[0];
		Mesh_3_Vector_3 v2 = pv[2] - pv[0];

		Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
		n /= sqrt(n*n);

		Mesh_3_Vector_3 vc = pv[0] - ci1->vertex(fi->second)->point().point();
		if (vc*n > 0)
			n = -n;

		Mesh_3_Vector_3 rv1 = mag1 * n;
		Mesh_3_Vector_3 rv2 = -mag2 * n;

		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1), m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv1.x()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1) + 1, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv1.y()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1) + 2, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv1.z()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2), m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv2.x()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2) + 1, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv2.y()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2) + 2, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv2.z()));
	}

	gradOp.resize(3 * m_mesh_3.NumC(), m_mesh_3.NumIF());
	gradOp.setFromTriplets(triplet.begin(), triplet.end());
}

void LaplacianAnalyzer::concatenateSPMatrix(EigenSpMat& M, EigenSpMat A, EigenSpMat B)
{
	M.resize(A.rows() + B.rows(), A.cols());
	M.reserve(A.nonZeros() + B.nonZeros());
	for (Eigen::Index c = 0; c < A.cols(); ++c)
	{
		M.startVec(c); // Important: Must be called once for each column before inserting!
		for (EigenSpMat::InnerIterator itL(A, c); itL; ++itL)
			M.insertBack(itL.row(), c) = itL.value();
		for (EigenSpMat::InnerIterator itC(B, c); itC; ++itC)
			M.insertBack(itC.row() + A.rows(), c) = itC.value();
	}
	M.finalize();
}

void LaplacianAnalyzer::concatenateVector(EigenVector& V, EigenVector A, EigenVector B)
{
	V.resize(A.size() + B.size());

	for (int i = 0; i < A.size(); ++i) {
		V(i) = A(i);
	}

	for (int i = 0; i < B.size(); ++i) {
		V(i + A.size()) = B(i);
	}
}

void LaplacianAnalyzer::writeTestMesh()
{
	ofstream out("compare.tet");

	out << "3" << endl;
	out << m_mesh_3.NumV() << endl;

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		int k;
		if (m_mesh_3.VertexOnBoundary(vi))
			k = 2;
		else 
			k = 3;

		out << vi->point().point().x() << " "
			<< vi->point().point().y() << " "
			<< vi->point().point().z() << " "
			<< "0 " << k << " 1" << endl;
	}

	std::map<Mesh_3_Cell_iterator, int> canonicalIdx;
	int kk = 0;
	out << m_mesh_3.C3T3().triangulation().number_of_cells() << endl;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		canonicalIdx[ci] = kk;
		
		for (int i = 0; i < 4; ++i) {
			if (m_mesh_3.C3T3().triangulation().is_infinite(ci->vertex(i)))
				out << "0" << " ";
			else
				out << m_mesh_3.VertexIdx(ci->vertex(i)) + 1 << " ";
		}
		out << endl;

		++kk;
	}
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		for (int i = 0; i < 4; ++i) {
			out<<canonicalIdx[ci->neighbor(i)]<<" ";
		}
		out << endl;
	}
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		out << ci->subdomain_index() << " ";
		for (int i = 0; i < 4; ++i) {
			My_facet mf(ci, i);
			Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);
			int status;
			if (m_mesh_3.FacetOnBoundary(fi))
				status = 1;
			else
				status = 0;
			out << status << " ";
		}
		out << endl;
	}

}