#include "MixedBoundary.h"

My_Triangulation& MixedBoundary::Mesh()
{
	return m_mesh_3;
}

std::vector<Mesh_3_Point_3>& MixedBoundary::Vertices()
{
	return vertices;
}

std::vector<int>& MixedBoundary::Indices()
{
	return indices;
}

std::vector<Mesh_3_Vector_3>& MixedBoundary::Normals()
{
	return normals;
}

std::vector<Mesh_3_Vector_3>& MixedBoundary::Colors()
{
	return colors;
}

std::vector<Mesh_3_Vector_3>& MixedBoundary::VfVertices()
{
	return vfVertices;
}

std::vector<Mesh_3_Vector_3>& MixedBoundary::VfFaces()
{
	return vfFaces;
}

std::vector<Mesh_3_Vector_3>& MixedBoundary::VfNormals()
{
	return vfNormals;
}

std::vector<Mesh_3_Vector_3>& MixedBoundary::VfColors()
{
	return vfColors;
}

void MixedBoundary::init(double dsty, double sr, int stps, std::string sfx)
{
	m_density_ratio = dsty;
	m_step_ratio = sr;
	m_steps = stps;
	m_suffix = sfx;
}

void MixedBoundary::buildMeshFromSurface(std::string fn, double size)
{
	Polyhedron_with_feature polyhedron;
	std::ifstream input(fn.c_str());
	input >> polyhedron;
	input.close();
	Mesh_domain_with_feature domain(polyhedron);
	domain.detect_features();

	Mesh_criteria criteria(edge_size = size,
		facet_angle = 30,
		facet_size = size,
		facet_distance = 0.1*size,
		cell_radius_edge_ratio = 3,
		cell_size = size);
	std::cout << "Criteria ... " << std::endl;

	m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, features(domain), odt());
	//CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	std::cout << "Triangulation ... " << std::endl;

	m_mesh_3.preprocessing();
}

void MixedBoundary::readSimulationData(std::string fn_tet, std::string fn_field)
{
	std::ifstream in(fn_tet.c_str());
	in >> m_mesh_3.C3T3();
	cout << "Reading Mesh Done ..." << endl;

	m_mesh_3.preprocessing();
	cout << "Preprocessing Done ..." << endl;

	readField(fn_field);
	cout << "Reading Field Done ..." << endl;
}

void MixedBoundary::prepare()
{
	assignMixedElements();
	cout<<"Mixed Elements Done ..."<<endl;
	assignMixedForwardBackward();
	cout << "Mixed Forward Backward Done ..." << endl;
	createOperators();
	cout << "Mixed Operators Done ..." << endl;
}

void MixedBoundary::decompose()
{
	
	setOmega();
	//setOmegaFromPWCF();
	cout << "Omega Set ..." << endl;

	computeHarmonic();
	cout<<"Harmonic Done ..."<<endl;
	computeGradient();
	cout << "Gradient Done ..." << endl;
	computeCurl();
	cout << "Curl Done ..." << endl;

	assignVertexColorMixBoundary();
	//computeNorms();
}

void MixedBoundary::visualize()
{
	convertForms();
	computeArrows();
	computeCrossSection();

	//writeMesh();
	writeMeshPLY();
}

void MixedBoundary::integrate()
{
	integrateField(omegaVF, omegaCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_o");
	integrateField(harmonicVF, harmonicCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_h");
	integrateField(curlVF, curlCF, m_density_ratio, m_step_ratio, m_steps, 1, m_suffix + "_c");
}

void MixedBoundary::readField(std::string fn_field)
{
	ifstream in(fn_field.c_str());

	std::string line;

	std::getline(in, line);
	std::stringstream ss(line);

	int nc;
	ss >> nc;

	for (int i = 0; i < nc; ++i) {
		std::getline(in, line);
		std::stringstream sss(line);

		double v1, v2, v3;
		sss >> v1 >> v2 >> v3;
		//cout << v1 << " " << v2 << " " << v3 << endl;
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);
		omegaVF[ci] = Mesh_3_Vector_3(v1, v2, v3);
	}
}

void MixedBoundary::assignMixedElements()
{
	// assign mixed vertices in X direction
	double cmin, cmax;
	cmin = -0.3683;
	cmax = 0.6317;
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (vi->point().point().y() < cmin + 0.0001 || vi->point().point().y() > cmax - 0.0001)
			vertexMixed[vi] = true;
		else
			vertexMixed[vi] = false;
	}

	// compute mixed other elements
	// edges
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		Mesh_3_Vertex_iterator vi1 = ei->first->vertex(ei->second);
		Mesh_3_Vertex_iterator vi2 = ei->first->vertex(ei->third);

		if (vertexMixed[vi1] && vertexMixed[vi2])
			edgeMixed[ei] = true;
		else
			edgeMixed[ei] = false;
	}

	//facets
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);

		std::vector<Mesh_3_Vertex_iterator> viv;
		for (int k = 0; k < 4; ++k) {
			if (fi->second != k) {
				viv.push_back(fi->first->vertex(k));
			}
		}

		if (vertexMixed[viv[0]] && vertexMixed[viv[1]] && vertexMixed[viv[2]])
			facetMixed[fi] = true;
		else
			facetMixed[fi] = false;
	}
}

void MixedBoundary::assignMixedForwardBackward()
{
	// K:kept M:mixed

	// vertex
	numKV = 0;
	numMV = 0;
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (!vertexMixed[vi]) {
			vertexForwardMixed[i] = numKV;
			vertexBackwardMixed[numKV] = i;
			++numKV;
		}
		else {
			++numMV;
		}
	}

	// edge
	numKE = 0;
	numME = 0;
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (!edgeMixed[ei]) {
			edgeForwardMixed[i] = numKE;
			vertexBackwardMixed[numKE] = i;
			++numKE;
		}
		else {
			++numME;
		}
	}

	// facet
	numKF = 0;
	numMF = 0;
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);

		if (!facetMixed[fi]) {
			facetForwardMixed[i] = numKF;
			facetBackwardMixed[numKF] = i;
			++numKF;
		}
		else {
			++numMF;
		}
	}
}

void MixedBoundary::createOperators()
{
	buildExteriorDerivative0Form();
	buildExteriorDerivative1Form();
	buildExteriorDerivative2Form();

	buildHodgeStar0Form();
	buildHodgeStar1Form();
	buildHodgeStar2Form();
	buildHodgeStar3Form();

	buildLaplacian0();
	buildLaplacian1();
	buildLaplacian2();
}

void MixedBoundary::setOmega()
{
	EigenSpMat s1l1 = m_hs1f * m_laplacian1;
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s1l1.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> hsrow, hscol;
	std::vector<double> hsval;
	for (int i = 0; i < m_hs1f.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_hs1f, i); iter; ++iter) {
			hsrow.push_back(static_cast<double>(iter.row()) + 1);
			hscol.push_back(static_cast<double>(iter.col()) + 1);
			hsval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, hsrow, hscol, hsval, 20);

	EigenVector form;
	form.setZero(m_mesh_3.NumE());

	form += 0.5*eigenFields[0];
	form += eigenFields[1];
	//form += eigenFields[3];
	//form += eigenFields[5];
	//form += eigenFields[7];
	
	mixedToFull1Form(form);

	omega = form;
}

void MixedBoundary::setOmegaFromPWCF()
{
	omega.resize(m_mesh_3.NumE());

	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (edgeMixed[ei]) {
			omega(i) = 0;
		}
		else {
			Mesh_3_Cell_circulator cc = m_mesh_3.C3T3().triangulation().incident_cells(*ei);
			Mesh_3_Cell_circulator end = cc;

			Mesh_3_Vector_3 v(0, 0, 0);
			double w = 0;
			do {
				if (cc->subdomain_index() != 0) {
					v += omegaVF[cc] * m_mesh_3.CellPrimal(cc);
					w += m_mesh_3.CellPrimal(cc);
				}

				++cc;
			} while (cc != end);
			v /= w;

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

			omega(i) = l * v;
		}
	}
}

void MixedBoundary::computeHarmonic()
{
	int hm = 2; // harmonic dimension
	
	EigenSpMat s1l1 = m_hs1f * m_laplacian1;
	std::vector<double> row, col;
	std::vector<double> val;
	for (int i = 0; i < s1l1.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(s1l1, i); iter; ++iter) {
			row.push_back(static_cast<double>(iter.row()) + 1);
			col.push_back(static_cast<double>(iter.col()) + 1);
			val.push_back(static_cast<double>(iter.value()));
		}
	}
	std::vector<double> hsrow, hscol;
	std::vector<double> hsval;
	for (int i = 0; i < m_hs1f.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_hs1f, i); iter; ++iter) {
			hsrow.push_back(static_cast<double>(iter.row()) + 1);
			hscol.push_back(static_cast<double>(iter.col()) + 1);
			hsval.push_back(static_cast<double>(iter.value()));
		}
	}

	std::vector<double> emptyEval;
	std::vector<EigenVector> eigenFields;
	matlabEIGS(emptyEval, eigenFields, row, col, val, hsrow, hscol, hsval, hm);

	EigenMatrix H;
	H.resize(m_mesh_3.NumE(), hm);

	for (int i = 0; i < hm; ++i) {
		EigenVector ev = eigenFields[i];
		mixedToFull1Form(ev);

		H.col(i) = ev;
	}

	harmonic = H * H.transpose()*(m_mesh_3.HS1F_B()*omega);
}

void MixedBoundary::computeGradient()
{
	EigenSpMat s0l0 = m_hs0f * m_laplacian0;

	// solving
	EigenVector omegaMixed = omega;
	fulltoMixed1Form(omegaMixed);

	EigenVector b = m_ed0f.transpose()*m_hs1f * omegaMixed;

	EigenVector x = matlabLinearSolveCholesky(s0l0, b);

	gradient = m_ed0f * x;
	mixedToFull1Form(gradient);

	mixedToFull0Form(x);
	assignColorGradient(x);
}

void MixedBoundary::computeCurl()
{
	EigenVector omegaMixed = omega;
	fulltoMixed1Form(omegaMixed);

	EigenSpMat s2l2 = m_hs2f * m_laplacian2;
	EigenVector b = m_hs2f * m_ed1f*omegaMixed;

	EigenVector x = matlabLinearSolveCholesky(s2l2, b);

	curl = m_hs1f.cwiseInverse()*m_ed1f.transpose()*m_hs2f*x;
	mixedToFull1Form(curl);

	gradient = omega - harmonic - curl;
}

void MixedBoundary::buildExteriorDerivative0Form()
{
	// orientation: increasing global indices
	std::vector<Eigen::Triplet<double>> triplet;

	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);
		if (edgeMixed[ei])
			continue;
		
		if (m_mesh_3.VertexIdx(ei->first->vertex(ei->second)) > m_mesh_3.VertexIdx(ei->first->vertex(ei->third))) {
			if (!vertexMixed[ei->first->vertex(ei->second)])
				triplet.push_back(Eigen::Triplet<double>(edgeForwardMixed[i],
					vertexForwardMixed[m_mesh_3.VertexIdx(ei->first->vertex(ei->second))], 1));
			if (!vertexMixed[ei->first->vertex(ei->third)])
				triplet.push_back(Eigen::Triplet<double>(edgeForwardMixed[i],
					vertexForwardMixed[m_mesh_3.VertexIdx(ei->first->vertex(ei->third))], -1));
		}
		else {
			if (!vertexMixed[ei->first->vertex(ei->second)])
				triplet.push_back(Eigen::Triplet<double>(edgeForwardMixed[i],
					vertexForwardMixed[m_mesh_3.VertexIdx(ei->first->vertex(ei->second))], -1));
			if (!vertexMixed[ei->first->vertex(ei->third)])
				triplet.push_back(Eigen::Triplet<double>(edgeForwardMixed[i],
					vertexForwardMixed[m_mesh_3.VertexIdx(ei->first->vertex(ei->third))], 1));
		}
		
	}

	m_ed0f.resize(numKE, numKV);
	m_ed0f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildExteriorDerivative1Form()
{
	// orientation: increasing global indices 
	std::vector<Eigen::Triplet<double>> triplet;

	for (int k = 0; k < m_mesh_3.NumF(); ++k) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(k);
		if (facetMixed[fi])
			continue;

		std::vector<int> vGlobalIdx;
		std::vector<int> vLocalIdx;
		for (int i = 0; i < 4; ++i) {
			if (i == fi->second)
				continue;

			vGlobalIdx.push_back(m_mesh_3.VertexIdx(fi->first->vertex(i)));
			vLocalIdx.push_back(i);
		}

		// make increasing order to be consistent with orientation
		for (int i = 0; i < vGlobalIdx.size(); ++i) {
			for (int j = i + 1; j < vGlobalIdx.size(); ++j) {
				if (vGlobalIdx[i] > vGlobalIdx[j]) { // swap
					int tgi, tli;

					tgi = vGlobalIdx[i];
					vGlobalIdx[i] = vGlobalIdx[j];
					vGlobalIdx[j] = tgi;

					tli = vLocalIdx[i];
					vLocalIdx[i] = vLocalIdx[j];
					vLocalIdx[j] = tli;
				}
			}
		}
		
		// integrate based on oritation	
		for (int i = 0; i < vGlobalIdx.size(); ++i) {
			int localIdx1, localIdx2;
			localIdx1 = vLocalIdx[i];
			localIdx2 = vLocalIdx[(i + 1) % vLocalIdx.size()];

			if (localIdx1 > localIdx2)
				std::swap(localIdx1, localIdx2);

			int val = 1;
			if (i == vGlobalIdx.size() - 1)
				val = -1;

			My_edge me(fi->first, localIdx1, localIdx2);
			Mesh_3_Edge_iterator ei = m_mesh_3.MyCellEdgeMap(me);

			if (!edgeMixed[ei])
				triplet.push_back(Eigen::Triplet<double>(facetForwardMixed[k], edgeForwardMixed[m_mesh_3.EdgeIdx(ei)], val));
		}
	}

	m_ed1f.resize(numKF, numKE);
	m_ed1f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildExteriorDerivative2Form()
{
	// orientation inherited from local index 0123
	std::vector<Eigen::Triplet<double>> triplet;
	
	for (int k = 0; k < m_mesh_3.NumC(); ++k) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(k);

		for (int i = 0; i < 4; ++i) {
			My_facet mf(ci, i);

			int even = 1;
			if (i % 2 == 1)
				even = -1;

			std::vector<int> lIdx;
			for (int j = 0; j < 4; ++j) {
				if (j != i)
					lIdx.push_back(j);
			}

			std::vector<int> gIdx;
			gIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(lIdx[0])));
			gIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(lIdx[1])));
			gIdx.push_back(m_mesh_3.VertexIdx(ci->vertex(lIdx[2])));

			for (int l = 0; l < gIdx.size(); ++l) {
				for (int j = l + 1; j < gIdx.size(); ++j) {
					if (gIdx[l] > gIdx[j]) {
						std::swap(gIdx[l], gIdx[j]);
						even = -even;
					}
				}
			}

			Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);

			if (!facetMixed[fi]) {
				if (fi->first == ci) {
					triplet.push_back(Eigen::Triplet<double>(k, facetForwardMixed[m_mesh_3.FacetIdx(m_mesh_3.MyCellFacetMap(mf))], even));
				}
				else {
					triplet.push_back(Eigen::Triplet<double>(k, facetForwardMixed[m_mesh_3.FacetIdx(m_mesh_3.MyCellFacetMap(mf))], even));
				}
			}
		}
	}

	m_ed2f.resize(m_mesh_3.NumC(), numKF);
	m_ed2f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildHodgeStar0Form()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (vertexMixed[vi])
			continue;

		triplet.push_back(Eigen::Triplet<double>(vertexForwardMixed[i],
			vertexForwardMixed[i],
			m_mesh_3.VertexDual(vi)));
	}
	m_hs0f.resize(numKV, numKV);
	m_hs0f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildHodgeStar1Form()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (edgeMixed[ei])
			continue;

		triplet.push_back(Eigen::Triplet<double>(edgeForwardMixed[i],
			edgeForwardMixed[i],
			m_mesh_3.EdgeDual(ei) / m_mesh_3.EdgePrimal(ei)));

	}

	m_hs1f.resize(numKE, numKE);
	m_hs1f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildHodgeStar2Form()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);

		if (facetMixed[fi])
			continue;

		triplet.push_back(Eigen::Triplet<double>(facetForwardMixed[i],
			facetForwardMixed[i],
			m_mesh_3.FacetDual(fi) / m_mesh_3.FacetPrimal(fi)));
	}

	m_hs2f.resize(numKF, numKF); // !!
	m_hs2f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildHodgeStar3Form()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		triplet.push_back(Eigen::Triplet<double>(i, i, 1 / m_mesh_3.CellPrimal(ci)));
	}

	m_hs3f.resize(m_mesh_3.NumC(), m_mesh_3.NumC());
	m_hs3f.setFromTriplets(triplet.begin(), triplet.end());
}

void MixedBoundary::buildLaplacian0()
{
	m_laplacian0 = m_hs0f.cwiseInverse()*m_ed0f.transpose()*m_hs1f*m_ed0f;
}

void MixedBoundary::buildLaplacian1()
{
	m_laplacian1 = m_hs1f.cwiseInverse() *m_ed1f.transpose()*m_hs2f*m_ed1f
		+ m_ed0f * m_hs0f.cwiseInverse()*m_ed0f.transpose()*m_hs1f;

	//m_divEM = m_mesh_3.ED0F_B()*m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();
	//m_curlEM = m_mesh_3.HS1F_B().cwiseInverse()*m_mesh_3.ED1F_B().transpose()*m_mesh_3.HS2F_B()*m_mesh_3.ED1F_B();
}

void MixedBoundary::buildLaplacian2()
{
	m_laplacian2 = m_hs2f.cwiseInverse() * m_ed2f.transpose()*m_hs3f*m_ed2f
		+ m_ed1f * m_hs1f.cwiseInverse() *m_ed1f.transpose()*m_hs2f;
}


void MixedBoundary::convertForms()
{
	convert1Form(omega, omegaVF);
	convert1Form(harmonic, harmonicVF);
	convert1Form(gradient, gradientVF);
	convert1Form(curl, curlVF);
}

void MixedBoundary::computeArrows()
{
	assignColorVectorField(omegaVF, 0);
	assignColorVectorField(harmonicVF, 1);
	assignColorVectorField(curlVF, 2);
	//assignColorVectorField(gradientVF, 2);

	rescaleArrowsWithVF(omegaVF);
	rescaleArrowsWithVF(harmonicVF);
	rescaleArrowsWithVF(curlVF);
	//rescaleArrowsWithVF(gradientVF);

	computeArrowsWithVF(harmonicVF);
}

void MixedBoundary::writeMesh()
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

void MixedBoundary::writeMeshPLY()
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

void MixedBoundary::computeCrossSection()
{
	double cutThd = 1;

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

EigenVector MixedBoundary::matlabLinearSolveCholesky(EigenSpMat A, EigenVector b)
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
	//str = "x=pcg(A,b);";						engEvalString(eng, str.c_str());
	str = "elapsed=toc";					engEvalString(eng, str.c_str());
	//str = "R=chol(A);";						engEvalString(eng, str.c_str());
	//str = "x=R\\(R'\\b);";					engEvalString(eng, str.c_str());
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

void MixedBoundary::matlabEIGS(std::vector<double>& eval, std::vector<EigenVector>& evec,
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

void MixedBoundary::integrateField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf,
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

void MixedBoundary::mixedToFull1Form(EigenVector& form)
{
	EigenVector formFull;
	formFull.resize(m_mesh_3.NumE());

	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (edgeMixed[ei])
			formFull[i] = 0;
		else
			formFull[i] = form(edgeForwardMixed[i]);
	}

	form = formFull;
}

void  MixedBoundary::mixedToFull0Form(EigenVector& form)
{
	EigenVector formFull;
	formFull.resize(m_mesh_3.NumV());

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (vertexMixed[vi])
			formFull[i] = 0;
		else
			formFull[i] = form(vertexForwardMixed[i]);
	}

	form = formFull;
}

void MixedBoundary::fulltoMixed1Form(EigenVector& form)
{
	EigenVector formMixed;
	formMixed.resize(numKE);

	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (edgeMixed[ei])
			continue;

		formMixed(edgeForwardMixed[i]) = form(i);
	}

	form = formMixed;
}

void MixedBoundary::convert1Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

void MixedBoundary::convert2Form(EigenVector form, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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

void MixedBoundary::computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	vfvNum = 0;

	int i = 0;
	for (Mesh_3_Cell_iterator ci = m_mesh_3.C3T3().triangulation().cells_begin();
		ci != m_mesh_3.C3T3().triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() != 0) {
			if (i % 5 == 0) {

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

void MixedBoundary::rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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
	pivotLen = lenVec[static_cast<int>(floor(lenVec.size()*0.7))];
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
		if (len < 0.1*pivotLen)
			len = 0.1*pivotLen;

		len = (0.5*aveLen)*(len / (pivotLen));

		vf[ci] = v * len;

		//cout<< vf[ci] <<endl;
	}
}

void MixedBoundary::assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci)
{
	Mesh_3_Vector_3 d = head - tail;
	//cout<< sqrt(d*d) <<endl;
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

void MixedBoundary::selectAnchors(std::vector<EigenVector> eigenfields, std::set<int>& anchors)
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

void MixedBoundary::correctRankDeficiency(EigenSpMat& m, EigenVector& b, std::set<int> anchors, EigenSpMat hs, std::vector<EigenVector> ef)
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

	int row = ef[0].size();
	int col = anchors.size();
	EigenSpMat H(row, col);
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			H.coeffRef(i, j) = ef[j](i);
		}
	}

	b -= hs * H*H.transpose()*b;

	for (auto iter = anchors.begin(); iter != anchors.end(); ++iter) {
		b(*iter) = 0;
	}
}

void MixedBoundary::sortedIndices(std::vector<int>& sIndex, std::vector<double> len)
{
	sIndex.resize(len.size(), 0);
	std::iota(sIndex.begin(), sIndex.end(), 0);

	std::sort(sIndex.begin(), sIndex.end(), [&len](int i1, int i2) {return len[i1] < len[i2]; });
}

Mesh_3_Vector_3  MixedBoundary::assignColorSF(double ratio)
{
	Mesh_3_Vector_3 color;

	//Mesh_3_Vector_3 red(0.9, 0.7, 0.7);
	//Mesh_3_Vector_3 yellow(0.9, 0.9, 0.7);
	//Mesh_3_Vector_3 green(0.7, 0.9, 0.7);
	//Mesh_3_Vector_3 indigo(0.7, 0.9, 0.9);
	//Mesh_3_Vector_3 blue(0.7, 0.7, 0.9);

	Mesh_3_Vector_3 red(0.8, 0.3, 0.6);
	Mesh_3_Vector_3 white(0.9, 0.9, 0.9);
	Mesh_3_Vector_3 blue(0.2, 0.6, 0.8);

	if (ratio < 0.4) {
		double r = ratio / 0.4;
		color = (1 - r)*blue + r * white;
	}
	else if (ratio < 0.6) {
		color = white;
	}
	else {
		double r = (ratio - 0.6) / 0.4;
		color = (1 - r)*white + r * red;
	}

	return color;
}

Mesh_3_Vector_3 MixedBoundary::assignColorVF(double ratio)
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

void MixedBoundary::assignColorGradient(EigenVector ev)
{
	std::vector<double> val;
	for (int i = 0; i < ev.size(); ++i) {
		val.push_back(ev(i));
	}
	std::sort(val.begin(), val.end());

	// assign color 
	//double pmin = ev(0);
	//double pmax = ev(0);
	//for (int i = 1; i < m_mesh_3.NumV(); ++i) {
	//	if (ev(i) > pmax)
	//		pmax = ev(i);
	//	if (ev(i) < pmin)
	//		pmin = ev(i);
	//}
	double pmin = val[floor(val.size()*0.2)];
	double pmax = val[floor(val.size()*0.8)];


	double diff = pmax - pmin;
	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		double val = ev(i);
		double d = val - (pmin);

		double ratio = d / (diff);
		if (ratio < 0)
			ratio = 0;
		if (ratio > 1)
			ratio = 1;

		colorField[vi] = assignColorVF(ratio);
	}

	//colorField = groundedGradientCF;
}

void MixedBoundary::assignColorVectorField(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf, int mode)
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

		Mesh_3_Vector_3 col = assignColorVF(ratio);

		if (mode==0) {
			omegaCF[vi] = col;
		}
		else if (mode == 1) {
			harmonicCF[vi] = col;
		}
		else if (mode == 2) {
			curlCF[vi] = col;
		}
		else if (mode == 3) {
			gradientCF[vi] = col;
		}
	}

	if (mode == 0) {
		colorField = omegaCF;
	}
	else if (mode == 1) {
		colorField = harmonicCF;
	}
	else if (mode == 2) {
		colorField = curlCF;
	}
	else if (mode == 3) {
		colorField = gradientCF;
	}
}

void MixedBoundary::assignVertexColorMixBoundary()
{
	/*for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);

		if (vertexMixed[vi])
			colorField[vi] = Mesh_3_Vector_3(0.8, 0.2, 0.2);
		else
			colorField[vi] = Mesh_3_Vector_3(0.6, 0.6, 0.6);
	}*/

	for (int i = 0; i < m_mesh_3.NumV(); ++i) {
		Mesh_3_Vertex_iterator vi = m_mesh_3.IdxToVertex(i);
		colorField[vi] = Mesh_3_Vector_3(0.6, 0.6, 0.6);
	}

	//for (int i = 0; i < m_mesh_3.NumE(); ++i) {
	//	Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

	//	if (edgeMixed[ei]) {
	//		Mesh_3_Vertex_iterator vi1 = ei->first->vertex(ei->second);
	//		Mesh_3_Vertex_iterator vi2 = ei->first->vertex(ei->third);

	//		colorField[vi1] = Mesh_3_Vector_3(0.8, 0.2, 0.2);
	//		colorField[vi2] = Mesh_3_Vector_3(0.8, 0.2, 0.2);
	//	}

	//}

	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);
	
		if (facetMixed[fi]) {
			for (int k = 0; k < 4; ++k) {
				if (k != fi->second)
					colorField[fi->first->vertex(k)] = Mesh_3_Vector_3(0.8, 0.2, 0.2);
			}
		}
	}
}

void MixedBoundary::computeNorms()
{
	double hn = sqrt((harmonic.transpose()*m_mesh_3.HS1F_B()*harmonic)(0, 0));
	double gn = sqrt((gradient.transpose()*m_mesh_3.HS1F_B()*gradient)(0, 0));
	double cn = sqrt((curl.transpose()*m_mesh_3.HS1F_B()*curl)(0, 0));

	cout << "Hamonic Norm: "<< hn << endl;
	cout << "Gradient Norm: "<< gn << endl;
	cout << "Curl Norm: " << cn << endl;
}