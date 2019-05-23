#include "Poelke.h"

My_Triangulation& Poelke::Mesh()
{
	return m_mesh_3;
}

std::vector<Mesh_3_Point_3>& Poelke::Vertices()
{
	return vertices;
}

std::vector<int>& Poelke::Indices()
{
	return indices;
}

std::vector<Mesh_3_Vector_3>& Poelke::Normals()
{
	return normals;
}

std::vector<Mesh_3_Vector_3>& Poelke::Colors()
{
	return colors;
}

std::vector<Mesh_3_Vector_3>& Poelke::VfVertices()
{
	return vfVertices;
}

std::vector<Mesh_3_Vector_3>& Poelke::VfFaces()
{
	return vfFaces;
}

std::vector<Mesh_3_Vector_3>& Poelke::VfNormals()
{
	return vfNormals;
}

std::vector<Mesh_3_Vector_3>& Poelke::VfColors()
{
	return vfColors;
}


void Poelke::init(double dsty, double sr, int stps, std::string sfx)
{
	m_density_ratio = dsty;
	m_step_ratio = sr;
	m_steps = stps;
	m_suffix = sfx;
}

void Poelke::buildMeshFromSurface(std::string fn, double size)
{
	Polyhedron_with_feature polyhedron;
	std::ifstream input(fn.c_str());
	input >> polyhedron;
	input.close();
	Mesh_domain_with_feature domain(polyhedron);
	//domain.detect_features();

	Mesh_criteria criteria(facet_angle = 30,
		facet_size = size,
		facet_distance = 0.1*size,
		cell_radius_edge_ratio = 2,
		cell_size = size);
	std::cout << "Criteria ... " << std::endl;

	m_mesh_3.C3T3() = CGAL::make_mesh_3<C3t3>(domain, criteria, features(domain), odt());
	//CGAL::perturb_mesh_3(m_mesh_3.C3T3(), domain, time_limit = 10);
	std::cout << "Triangulation ... " << std::endl;

	m_mesh_3.preprocessing();
}

void Poelke::prepare()
{
	cout<<"Creating Operators ..."<<endl;
	createCurlOperator_B();
	createGradOperator_B();
	createCurlOperator_NB();
	createGradOperator_NB();
	createMassMatrix();
	cout<<"Operators Done ..."<<endl;

	testOrthogonalityOfBasis();
	
}

void Poelke::decompose()
{
	createInput();

	cout<<"Compute TC ..."<<endl;
	solveTangentalCurl();
	cout << "Compute NG ..." << endl;
	solveNormalGradient();
	cout << "Compute TH ..." << endl;
	solveTangentialHarmonic();
	cout << "Compute CH ..." << endl;
	solveCentralHarmonic();
	cout << "Compute NH ..." << endl;
	solveNormalHarmonic();

	//solveTHBasis();
	solveNGBasis();

	//checkCurlFieldBoundary();

	evaluateError();
	evaluateL2Norm();
}

void Poelke::visualize()
{
	computeArrows();
	computeCrossSection();
}

void Poelke::integrate()
{

}

void Poelke::createInput()
{
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		Mesh_3_Vector_3 v(0, 0, 0);
		Mesh_3_Vector_3 p = m_mesh_3.CellCC(ci) - CGAL::ORIGIN;

		//v += normalGradient(p);
		//v += tangentialCurl(p);
		v += normalHarmonic(p);
		//v += tangentialHarmonic(p);
		//v += centralHarmonic(p);

		m_input[ci] = v;
	}

	m_original_input = m_input;
}

void Poelke::solveNormalGradient()
{
	int a;

	EigenVector f;
	toEigenVector(f, m_input);

	EigenSpMat A = m_gradOp_NB.transpose()*m_Mass*m_gradOp_NB;
	EigenVector b = m_gradOp_NB.transpose()*m_Mass*f;
	
	EigenVector x = matlabLinearSolveCholesky(A, b);
	EigenVector ngf = m_gradOp_NB * x;

	f -= ngf;
	fromEigenVector(m_normalGradient, ngf);
	fromEigenVector(m_input, f);
}

void Poelke::solveTangentalCurl()
{
	EigenVector f;
	toEigenVector(f, m_input);

	EigenSpMat A = m_curlOp_NB.transpose()*m_Mass*m_curlOp_NB;
	EigenSpMat B = m_mesh_3.HS0F_NB().cwiseInverse()*m_mesh_3.ED0F_NB().transpose()*m_mesh_3.HS1F_NB();
	EigenVector a = m_curlOp_NB.transpose()*m_Mass*f;
	EigenVector b;
	b.setZero(m_mesh_3.NumIV());

	EigenSpMat M;
	concatenateVerticallySPMatrix(M, A, B);
	EigenVector v;
	concatenateVector(v, a, b);

	EigenVector x = matlabLinearSolveQR(M, v);
	EigenVector tcf = m_curlOp_NB * x;

	//EigenSpMat M;
	////concatenateVerticallySPMatrix(M, m_curlOp_NB.transpose()*m_Mass, m_gradOp_B.transpose()*m_Mass);
	//concatenateHorizontallySPMatrix(M, m_curlOp_NB, m_gradOp_B);
	//EigenVector a = m_curlOp_NB.transpose()*m_Mass*f;
	//EigenVector b;
	//b.setZero(3 * m_mesh_3.NumC() - m_mesh_3.NumIE());

	//cout << M.rows() << " " << M.cols() << endl;

	//EigenVector v;
	//concatenateVector(v, a, b);
	//cout << v.size() << endl;

	//EigenVector z = matlabLinearSolveCholesky(M*M.transpose(), v);
	//EigenVector x = M.transpose()*z;

	////int asd;
	////cin >> asd;

	//EigenVector tcf;
	//tcf.resize(3 * m_mesh_3.NumC());
	//for (int i = 0; i < 3 * m_mesh_3.NumC(); ++i) {
	//	tcf(i) = x(i);
	//}

	f -= tcf;
	fromEigenVector(m_tangentialCurl, tcf);
	fromEigenVector(m_input, f);
}

void Poelke::solveNormalHarmonic()
{
	m_normalHarmonic = m_input;
}

void Poelke::solveTangentialHarmonic()
{
	EigenVector f;
	toEigenVector(f, m_input);

	EigenSpMat A = m_gradOp_B.transpose()*m_Mass*m_gradOp_B;
	EigenVector b = m_gradOp_B.transpose()*m_Mass*f;

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
	b(0) = 0;

	EigenVector x = matlabLinearSolveCholesky(A, b);
	EigenVector rf = m_gradOp_B * x;

	f -= rf;
	fromEigenVector(m_tangentialHarmonic, f);
	fromEigenVector(m_input, rf);

	EigenVector ff;
	toEigenVector(ff, m_original_input);

	/*for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Vector_3 a(f(3 * i), f(3 * i + 1), f(3 * i + 2));
		Mesh_3_Vector_3 b(ff(3 * i), ff(3 * i + 1), ff(3 * i + 2));

		cout << sqrt(a*a) / sqrt(b*b) << endl;
	}*/
}

void Poelke::solveCentralHarmonic()
{
	EigenVector f;
	toEigenVector(f, m_input);

	EigenSpMat A = m_curlOp_B.transpose()*m_Mass*m_curlOp_B;
	EigenSpMat B = m_mesh_3.HS0F_B().cwiseInverse()*m_mesh_3.ED0F_B().transpose()*m_mesh_3.HS1F_B();
	EigenVector a = m_curlOp_B.transpose()*m_Mass*f;
	EigenVector b;
	b.setZero(m_mesh_3.NumV());

	EigenSpMat M;
	concatenateVerticallySPMatrix(M, A, B);
	EigenVector v;
	concatenateVector(v, a, b);

	EigenVector x = matlabLinearSolveQR(M, v);
	EigenVector chf = m_curlOp_B * x;

	f -= chf;
	fromEigenVector(m_centralHarmonic, chf);
	fromEigenVector(m_input, f);
}

void Poelke::solveTHBasis()
{
	EigenSpMat L;
	concatenateVerticallySPMatrix(L, m_curlOp_NB.transpose()*m_Mass, m_gradOp_B.transpose()*m_Mass);

	std::vector<EigenVector> svec;
	std::vector<double> sval;
	matlabSVDS(svec, sval, L, 1);
	//cout << svec[0] << endl;
	//cout << svec.size() << endl;
	//fromEigenVector(m_tangentialHarmonic, svec[0]);

	EigenVector f;
	toEigenVector(f, m_original_input);

	EigenVector H = svec[0];

	/*for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Vector_3 a(f(3 * i), f(3 * i + 1), f(3 * i + 2));
		Mesh_3_Vector_3 b(H(3 * i), H(3 * i + 1), H(3 * i + 2));
	
		cout << sqrt(a*a) / sqrt(b*b) << endl;
	}*/

	//EigenVector PP = m_curlOp_NB.transpose()*m_Mass*H1;
	//EigenVector QQ = m_gradOp_B.transpose()*m_Mass*H1;
	//EigenVector RR = L * H1;

	//cout << PP.dot(PP) << endl;
	//cout << QQ.dot(QQ) << endl;
	//cout << RR.dot(RR) << endl;
	//cout << H0.transpose()*m_Mass*H1 << endl;

	EigenMatrix A = H.transpose()*m_Mass*H;
	EigenVector b = H.transpose()*m_Mass*f;

	double x = b(0) / A.coeff(0, 0);

	EigenVector h = x * H;
	fromEigenVector(m_tangentialHarmonic, h);
}

void Poelke::solveNGBasis()
{
	EigenSpMat L;
	concatenateVerticallySPMatrix(L, m_curlOp_B.transpose()*m_Mass, m_gradOp_NB.transpose()*m_Mass);

	std::vector<EigenVector> svec;
	std::vector<double> sval;
	matlabSVDS(svec, sval, L, 1);
	//cout << svec.size() << endl;

	//fromEigenVector(m_normalHarmonic, svec[0]);

	EigenVector f;
	toEigenVector(f, m_original_input);
	EigenVector H = svec[0];

	/*for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Vector_3 a(f(3 * i), f(3 * i + 1), f(3 * i + 2));
		Mesh_3_Vector_3 b(H(3 * i), H(3 * i + 1), H(3 * i + 2));

		cout << sqrt(a*a) / sqrt(b*b) << endl;
	}*/

	EigenVector PP = m_curlOp_B.transpose()*m_Mass*H;
	EigenVector QQ = m_gradOp_NB.transpose()*m_Mass*H;
	EigenVector RR = L * H;

	cout << PP.dot(PP) << endl;
	cout << QQ.dot(QQ) << endl;
	cout << RR.dot(RR) << endl;

	EigenMatrix A = H.transpose()*m_Mass*H;
	EigenVector b = H.transpose()*m_Mass*f;

	double x = b(0) / A.coeff(0, 0);

	EigenVector h = x * H;
	fromEigenVector(m_normalHarmonic, h);
}

void Poelke::computeArrows()
{
	rescaleArrowsWithVF(m_original_input);
	computeArrowsWithVF(m_original_input);
}

void Poelke::computeCrossSection()
{
	double diff = m_mesh_3.BoundingBoxMax().y() - m_mesh_3.BoundingBoxMin().y();
	double pl = 0.5*diff + m_mesh_3.BoundingBoxMin().y();
	double cutThd = -0.05;

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

void Poelke::evaluateError()
{
	//double ttErr = computeFieldError(0);
	//double ttip = computeFieldInnerProduct(0);
	//double ngErr = computeFieldError(1);
	//double ngip = computeFieldInnerProduct(1);
	//double tcErr = computeFieldError(2);
	//double tcip = computeFieldInnerProduct(2);

	double nhErr = computeFieldError(3);
	double nhip = computeFieldInnerProduct(3);
	//double thErr = computeFieldError(4);
	//double thip = computeFieldInnerProduct(4);
	//double chErr = computeFieldError(5);
	//double chip = computeFieldInnerProduct(5);

	cout<<"*** Numerical Errors ... ***"<<endl;
	//cout << "TT Err: " << ttErr / ttip << endl;
	//cout << "NG Err: " << ngErr / ngip << endl;
	//cout << "TC Err: " << tcErr / tcip << endl;
	cout << "NH Err: " << nhErr / nhip << endl;
	//cout << "TH Err: " << thErr / thip << endl;
	//cout << "CC Err: " << chErr / chip << endl;
}

void Poelke::evaluateL2Norm()
{
	EigenVector in;
	toEigenVector(in, m_original_input);
	EigenVector ng;
	toEigenVector(ng, m_normalGradient);
	EigenVector tc;
	toEigenVector(tc, m_tangentialCurl);
	EigenVector nh;
	toEigenVector(nh, m_normalHarmonic);
	EigenVector th;
	toEigenVector(th, m_tangentialHarmonic);
	EigenVector ch;
	toEigenVector(ch, m_centralHarmonic);


	double inl2 = (in.transpose()*m_Mass*in).coeff(0, 0);
	double ngl2 = (ng.transpose()*m_Mass*ng).coeff(0, 0);
	double tcl2 = (tc.transpose()*m_Mass*tc).coeff(0, 0);
	double nhl2 = (nh.transpose()*m_Mass*nh).coeff(0, 0);
	double thl2 = (th.transpose()*m_Mass*th).coeff(0, 0);
	double chl2 = (ch.transpose()*m_Mass*ch).coeff(0, 0);

	double diff = inl2 - ngl2 - tcl2 - nhl2 - thl2 - chl2;

	cout << "Input L2-Norm: " << sqrt(inl2) << endl;
	cout << "Normal Gradient L2-Norm: " << sqrt(ngl2) << endl;
	cout << "Tagential Curl L2-Norm: " << sqrt(tcl2) << endl;
	cout << "Normal Harmonic L2-Norm: " << sqrt(nhl2) << endl;
	cout << "Tangential Harmonic L2-Norm: " << sqrt(thl2) << endl;
	cout << "Central Harmonic L2-Norm: " << sqrt(chl2) << endl;
	cout << "Summation Difference: " << diff << endl;
}

void Poelke::createCurlOperator_B()
{
	// Curl operator is a span of basis
	// It should be an 3|C| by |E| matrix mapping 1-form to PWC vector field in each tet.
	// Each column should be a basis of curl space.

	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_mesh_3.C3T3().triangulation().edges_begin();
		ei != m_mesh_3.C3T3().triangulation().edges_end(); ++ei) {
		if (!m_mesh_3.EdgeValid(ei))
			continue;

		Mesh_3_Cell_circulator ccir = m_mesh_3.C3T3().triangulation().incident_cells(*ei);
		Mesh_3_Cell_circulator end = ccir;
		
		do {
			if (ccir->subdomain_index() == 0) {
				++ccir;
			}
			else {
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
					if (vc*n > 0)
						n = -n;

					//if (i % 2 == 0)
						//n = -n;

					grad[m_mesh_3.VertexIdx(ccir->vertex(i))] = n / h;
				}
				
				// get 4 facet 2-form value after curl (d_1)
				std::vector<int> facetVal;
				for (int i = 0; i < 4; ++i) {
					My_facet mf(ccir, i);
					Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);
					double val = 0;

					EigenVector colv;
					colv.setZero(m_mesh_3.NumE());
					colv(m_mesh_3.EdgeIdx(ei)) = 1;
					EigenMatrix rowv = m_mesh_3.ED1F_B().row(m_mesh_3.FacetIdx(fi));
					val = (rowv *colv)(0, 0);

					facetVal.push_back(val);
				}
				
				// get vector for this tet
				Mesh_3_Vector_3 rv(0, 0, 0); // recovered vector
				for (int i = 0; i < 4; ++i) { //iterate through 4 faces;
					std::vector<int> vIdx;
					for (int j = 0; j < 4; ++j) {
						if (j != i) {
							vIdx.push_back(m_mesh_3.VertexIdx(ccir->vertex(j)));
						}
					}

					//std::sort(vIdx.begin(), vIdx.end());
					for (int k = 0; k < vIdx.size(); ++k) {
						for (int l = k + 1; l < vIdx.size(); ++l) {
							if (vIdx[k] > vIdx[l]) {
								std::swap(vIdx[k], vIdx[l]);
							}
						}
					}
					My_facet mf(ccir, i);

					rv += facetVal[i] * 2 * 0.25
						* (CGAL::cross_product(grad[vIdx[0]], grad[vIdx[1]])
							+ CGAL::cross_product(grad[vIdx[1]], grad[vIdx[2]])
							+ CGAL::cross_product(grad[vIdx[2]], grad[vIdx[0]]));

				}

				triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir), m_mesh_3.EdgeIdx(ei), rv.x()));
				triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir) + 1, m_mesh_3.EdgeIdx(ei), rv.y()));
				triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir) + 2, m_mesh_3.EdgeIdx(ei), rv.z()));

				++ccir;
			}
		} while (ccir != end);
	}

	m_curlOp_B.resize(3 * m_mesh_3.NumC(), m_mesh_3.NumE());
	m_curlOp_B.setFromTriplets(triplet.begin(), triplet.end());
}

void Poelke::createGradOperator_B()
{
	// Div operator is a span of basis
	// It should be an 3|C| by |F| matrix mapping 1-form to PWC vector field in each tet.
	// Each column should be a basis of curl space.

	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Facet_iterator fi = m_mesh_3.C3T3().triangulation().facets_begin();
		fi != m_mesh_3.C3T3().triangulation().facets_end(); ++fi) {
		if (!m_mesh_3.FacetValid(fi))
			continue;

		/*if (m_mesh_3.FacetIdx(fi) == m_mesh_3.NumF() - 1)
			continue;*/

		Mesh_3_Cell_iterator ci1 = fi->first;
		Mesh_3_Cell_iterator ci2 = fi->first->neighbor(fi->second);

		double mag1;
		double mag2;
		if (ci1->subdomain_index() != 0)
			mag1 = m_mesh_3.FacetPrimal(fi) / m_mesh_3.CellPrimal(ci1);
		if (ci2->subdomain_index() != 0)
			mag2 = m_mesh_3.FacetPrimal(fi) / m_mesh_3.CellPrimal(ci2);

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

		//Mesh_3_Vector_3 vc = pv[0] - ci1->vertex(fi->second)->point().point();
		//if (vc*n > 0)
			//n = -n;
		
		if (fi->second % 2 == 0)
			n = -n;

		if (ci1->subdomain_index() != 0) {
			Mesh_3_Vector_3 rv1 = mag1 * n;
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1), m_mesh_3.FacetIdx(fi), rv1.x()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1) + 1, m_mesh_3.FacetIdx(fi), rv1.y()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1) + 2, m_mesh_3.FacetIdx(fi), rv1.z()));
		}
		if (ci2->subdomain_index() != 0) {
			Mesh_3_Vector_3 rv2 = -mag2 * n;
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2), m_mesh_3.FacetIdx(fi), rv2.x()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2) + 1, m_mesh_3.FacetIdx(fi), rv2.y()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2) + 2, m_mesh_3.FacetIdx(fi), rv2.z()));
		}
	}

	m_gradOp_B.resize(3 * m_mesh_3.NumC(), m_mesh_3.NumF());
	m_gradOp_B.setFromTriplets(triplet.begin(), triplet.end());
}

void Poelke::createCurlOperator_NB()
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
				if (vc*n > 0)
					n = -n;

				//if (i % 2 == 0)
					//n = -n;

				grad[m_mesh_3.VertexIdx(ccir->vertex(i))] = n / h;
			}

			// get 4 facet 2-form value after curl (d_1)
			std::vector<int> facetVal;
			for (int i = 0; i < 4; ++i) {
				My_facet mf(ccir, i);
				Mesh_3_Facet_iterator fi = m_mesh_3.MyCellFacetMap(mf);
				double val = 0;

				EigenVector colv;
				colv.setZero(m_mesh_3.NumE());
				colv(m_mesh_3.EdgeIdx(ei)) = 1;
				EigenMatrix rowv = m_mesh_3.ED1F_B().row(m_mesh_3.FacetIdx(fi));
				val = (rowv *colv)(0, 0);

				facetVal.push_back(val);
			}

			// get vector for this tet
			Mesh_3_Vector_3 rv(0, 0, 0); // recovered vector
			for (int i = 0; i < 4; ++i) { //iterate through 4 faces;
				std::vector<int> vIdx;
				for (int j = 0; j < 4; ++j) {
					if (j != i) {
						vIdx.push_back(m_mesh_3.VertexIdx(ccir->vertex(j)));
					}
				}

				//std::sort(vIdx.begin(), vIdx.end());
				for (int k = 0; k < vIdx.size(); ++k) {
					for (int l = k + 1; l < vIdx.size(); ++l) {
						if (vIdx[k] > vIdx[l]) {
							std::swap(vIdx[k], vIdx[l]);
						}
					}
				}
				My_facet mf(ccir, i);

				rv += facetVal[i] * 2 * 0.25
					* (CGAL::cross_product(grad[vIdx[0]], grad[vIdx[1]])
						+ CGAL::cross_product(grad[vIdx[1]], grad[vIdx[2]])
						+ CGAL::cross_product(grad[vIdx[2]], grad[vIdx[0]]));

			}

			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir), m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), rv.x()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir) + 1, m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), rv.y()));
			triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ccir) + 2, m_mesh_3.EdgeForward(m_mesh_3.EdgeIdx(ei)), rv.z()));

			++ccir;
		} while (ccir != end);
	}

	m_curlOp_NB.resize(3 * m_mesh_3.NumC(), m_mesh_3.NumIE());
	m_curlOp_NB.setFromTriplets(triplet.begin(), triplet.end());
}

void Poelke::createGradOperator_NB()
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

		/*Mesh_3_Vector_3 vc = pv[0] - ci1->vertex(fi->second)->point().point();
		if (vc*n < 0)
			n = -n;*/

		if (fi->second % 2 == 0)
			n = -n;

		Mesh_3_Vector_3 rv1 = mag1 * n;
		Mesh_3_Vector_3 rv2 = -mag2 * n;

		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1), m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv1.x()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1) + 1, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv1.y()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci1) + 2, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv1.z()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2), m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv2.x()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2) + 1, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv2.y()));
		triplet.push_back(Eigen::Triplet<double>(3 * m_mesh_3.CellIdx(ci2) + 2, m_mesh_3.FacetForward(m_mesh_3.FacetIdx(fi)), rv2.z()));
	
		////test
		//Mesh_3_Point_3 cc = m_mesh_3.FacetCC(fi);
		//Mesh_3_Point_3 p1 = cc - 3 * (1 / mag1)*n;
		//Mesh_3_Point_3 p2 = cc + 3 * (1 / mag2)*n;
		//Mesh_3_Vector_3 tv1 = p1 - ci1->vertex(fi->second)->point().point();
		//Mesh_3_Vector_3 tv2 = p2 - m_mesh_3.C3T3().triangulation().mirror_vertex(fi->first, fi->second)->point().point();
		//cout << mag1 << " " << mag2 << endl;
		//cout << "Ortho? : " << tv1*rv1 << endl;
		//cout << "Ortho? : " << tv2*rv2 << endl;
	}

	m_gradOp_NB.resize(3 * m_mesh_3.NumC(), m_mesh_3.NumIF());
	m_gradOp_NB.setFromTriplets(triplet.begin(), triplet.end());

	//cout<< m_gradOp_NB <<endl;
}

void Poelke::createCNbGFOperator()
{
	std::vector<Eigen::Triplet<double>> triplet;

	int ebidx = 0;
	for (int i = 0; i < m_mesh_3.NumE(); ++i) {
		Mesh_3_Edge_iterator ei = m_mesh_3.IdxToEdge(i);

		if (m_mesh_3.EdgeOnBoundary(ei)) {
			


			ebidx++;
		}
	}
}

void Poelke::createGF0GFOperator()
{

}

void Poelke::createMassMatrix()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		triplet.push_back(Eigen::Triplet<double>(3 * i, 3 * i, m_mesh_3.CellPrimal(ci)));
		triplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, m_mesh_3.CellPrimal(ci)));
		triplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, m_mesh_3.CellPrimal(ci)));		
	}

	m_Mass.resize(3 * m_mesh_3.NumC(), 3 * m_mesh_3.NumC());
	m_Mass.setFromTriplets(triplet.begin(), triplet.end());
}

void Poelke::rescaleArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf)
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
	pivotLen = lenVec[static_cast<int>(floor(lenVec.size()*0.99))];
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
		if (len < 0.6*pivotLen)
			len = 0.6*pivotLen;

		len = (0.6*aveLen)*(len / (pivotLen));

		vf[ci] = v * len;
	}
}

void Poelke::computeArrowsWithVF(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
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

void Poelke::assembleArrow(Mesh_3_Vector_3 tail, Mesh_3_Vector_3 head, double radius, Mesh_3_Cell_iterator ci)
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

Mesh_3_Vector_3 Poelke::normalGradient(Mesh_3_Vector_3 p)
{
	return p;

	/*double a = p1 * (p2 - p1);
	double b = (p2 - p1)*(p2 - p1);

	return a + b / 2;*/
}

Mesh_3_Vector_3 Poelke::tangentialCurl(Mesh_3_Vector_3 p)
{
	return Mesh_3_Vector_3(p.y(), -p.x(), 0);

	/*Mesh_3_Vector_3 d = p2 - p1;
	double a = d.x()*p1.y() - d.y()*p1.x();
	double b = d.x()*d.y() - d.y()*d.x();

	return a + b / 2;*/
}

Mesh_3_Vector_3 Poelke::normalHarmonic(Mesh_3_Vector_3 p)
{
	double de = 1 / pow((p*p), 3.0 / 2.0);
	return de * p;

	/*Mesh_3_Vector_3 d = p2 - p1;
	double c1 = p1 * p1;
	double c2 = p1 * d;
	double c3 = d * d;

	double f0 = 2 * (c2*c2 - 2 * c1*c3) / ((c2*c2 - 4 * c1*c3)*sqrt(c1));
	double f1 = 2 * (c2*c3 + c2 * c2 - 2 * c1*c3) / ((c2*c2 - 4 * c1*c3)*sqrt(c1 + c2 + c3));

	return f1 - f0;*/
}

Mesh_3_Vector_3 Poelke::tangentialHarmonic(Mesh_3_Vector_3 p)
{
	double de = 1 / (p.x()*p.x() + p.y()*p.y());
	return de * Mesh_3_Vector_3(p.y(), -p.x(), 0);

	/*Mesh_3_Vector_3 d = p2 - p1;
	double c1 = p1.x()*p1.x() + p2.x()*p2.x();
	double c2 = p1.x()*d.x() + p1.y()*d.y();
	double c3 = d.x() * d.x() + d.x() * d.x();
	double c4 = d.x()*p1.y() - d.y()*p1.x();
	double c5 = d.x()*d.y() - d.y()*d.x();
	
	double f0 = (c5 * log(c1) - (2 * (c2*c5 - 2 * c3*c4)*atan(c2) / sqrt(4 * c1*c3 - c2 * c2)) / sqrt(4 * c1*c3 - c2 * c2)) / (2 * c3);
	double f1 = (c5 * log(c1 + c2 + c3) - (2 * (c2*c5 - 2 * c3*c4)*atan(2 * c3 + c2) / sqrt(4 * c1*c3 - c2 * c2)) / sqrt(4 * c1*c3 - c2 * c2)) / (2 * c3);

	return f1 - f0;*/
}

Mesh_3_Vector_3 Poelke::centralHarmonic(Mesh_3_Vector_3 p)
{
	return Mesh_3_Vector_3(1, 1, 1)*0.5;

	//Mesh_3_Vector_3 cv(1, 1, 1);
	//double a = (p2 - p1)*cv;

	//return a;
}

EigenVector Poelke::matlabLinearSolveQR(EigenSpMat A, EigenVector b)
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
	str = "tic;";					engEvalString(eng, str.c_str());
	str = "[C,R]=qr(A,b);";					engEvalString(eng, str.c_str());
	str = "x=R\\C;";							engEvalString(eng, str.c_str());
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

	for (int i = 0; i < dims[0]; ++i) {
		solution(i) = x[i];
	}

	return solution;
}

EigenVector Poelke::matlabLinearSolveCholesky(EigenSpMat A, EigenVector b)
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
	printf("%s", buffer);
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

	cout<< dims[0] <<" "<<dims[1]<<endl;

	for (int i = 0; i < dims[0]; ++i) {
		solution(i) = x[i];
	}

	return solution;
}

void Poelke::matlabSVDS(std::vector<EigenVector>& svec, std::vector<double>& sval, EigenSpMat A, int k)
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

	Engine* eng;
	eng = engOpen("\0");

	std::string str;

	// transfer data
	mxArray *rowArr, *colArr, *valArr;
	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);

	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));

	engPutVariable(eng, "row", rowArr);
	engPutVariable(eng, "col", colArr);
	engPutVariable(eng, "val", valArr);

	const int buffersize = 256;
	char buffer[buffersize + 1];
	buffer[buffersize] = '\0';
	engOutputBuffer(eng, buffer, buffersize);
	str = "A=sparse(row, col, val);";       engEvalString(eng, str.c_str());
	//str = "[r,c] = size(A);";       engEvalString(eng, str.c_str());
	//str = "A=A(1:r-1,:);";       engEvalString(eng, str.c_str());
	str = "tic;";							engEvalString(eng, str.c_str());
	str = "[U,S,V]=svds(A, 1, 'smallest');";       engEvalString(eng, str.c_str());
	str = "elapsed=toc";					engEvalString(eng, str.c_str());
	printf("%s", buffer);

	mxArray* vv;
	vv = engGetVariable(eng, "V");
	double* v;
	v = reinterpret_cast<double*>(mxGetData(vv));
	const size_t* dims;
	dims = mxGetDimensions(vv);

	// read out eigenvector matrix V from matlab
	mxArray* dd;
	dd = engGetVariable(eng, "S");
	double* d;
	d = reinterpret_cast<double*>(mxGetData(dd));

	EigenVector ev; 
	ev.resize(dims[0]);

	svec.clear();
	for (int j = 0; j < dims[1]; ++j) {
		for (int i = 0; i < dims[0]; ++i) {
			ev(i) = v[i + j * dims[0]];
		}
		svec.push_back(ev);
	}

	sval.clear();
	for (int i = 0; i < dims[1]; ++i) {
		sval.push_back(d[i + i * dims[1]]);
	}
}

void Poelke::matlabEIGS(std::vector<EigenVector>& evec, std::vector<double>& eval, EigenSpMat A, EigenSpMat B, int k)
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

	std::vector<double> hrow, hcol;
	std::vector<double> hval;
	for (int i = 0; i < B.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(B, i); iter; ++iter) {
			hrow.push_back(static_cast<double>(iter.row()) + 1);
			hcol.push_back(static_cast<double>(iter.col()) + 1);
			hval.push_back(static_cast<double>(iter.value()));
		}
	}

	Engine* eng;
	eng = engOpen("\0");

	std::string str;

	// transfer data
	mxArray *rowArr, *colArr, *valArr, *hrowArr, *hcolArr, *hvalArr;
	rowArr = mxCreateDoubleMatrix(row.size(), 1, mxREAL);
	colArr = mxCreateDoubleMatrix(col.size(), 1, mxREAL);
	valArr = mxCreateDoubleMatrix(val.size(), 1, mxREAL);
	hrowArr = mxCreateDoubleMatrix(hrow.size(), 1, mxREAL);
	hcolArr = mxCreateDoubleMatrix(hcol.size(), 1, mxREAL);
	hvalArr = mxCreateDoubleMatrix(hval.size(), 1, mxREAL);
	memcpy(mxGetPr(rowArr), &row[0], row.size() * sizeof(double));
	memcpy(mxGetPr(colArr), &col[0], col.size() * sizeof(double));
	memcpy(mxGetPr(valArr), &val[0], val.size() * sizeof(double));
	memcpy(mxGetPr(hrowArr), &hrow[0], hrow.size() * sizeof(double));
	memcpy(mxGetPr(hcolArr), &hcol[0], hcol.size() * sizeof(double));
	memcpy(mxGetPr(hvalArr), &hval[0], hval.size() * sizeof(double));
	engPutVariable(eng, "row", rowArr);
	engPutVariable(eng, "col", colArr);
	engPutVariable(eng, "val", valArr);
	engPutVariable(eng, "hrow", hrowArr);
	engPutVariable(eng, "hcol", hcolArr);
	engPutVariable(eng, "hval", hvalArr);
	

	str = "A=sparse(row, col, val)";        engEvalString(eng, str.c_str());
	str = "B=sparse(brow, bcol, bval)";        engEvalString(eng, str.c_str());
	str = "[V,D]=eigs(A, B, 3, 'smallestabs')";        engEvalString(eng, str.c_str());

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

void Poelke::concatenateVerticallySPMatrix(EigenSpMat& M, EigenSpMat A, EigenSpMat B)
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


void Poelke::concatenateHorizontallySPMatrix(EigenSpMat& M, EigenSpMat A, EigenSpMat B)
{
	M.resize(A.rows(), A.cols() + B.cols());
	M.reserve(A.nonZeros() + B.nonZeros());

	for (Eigen::Index c = 0; c < A.cols(); ++c){
		M.startVec(c); // Important: Must be called once for each column before inserting!
		for (EigenSpMat::InnerIterator itL(A, c); itL; ++itL)
			M.insertBack(itL.row(), c) = itL.value();
	}
	for (Eigen::Index c = 0; c < B.cols(); ++c){
		M.startVec(c + A.cols());
		for (EigenSpMat::InnerIterator itC(B, c); itC; ++itC)
			M.insertBack(itC.row(), c + A.cols()) = itC.value();
	}
	M.finalize();
}

void Poelke::concatenateVector(EigenVector& V, EigenVector A, EigenVector B)
{
	V.resize(A.size() + B.size());

	for (int i = 0; i < A.size(); ++i) {
		V(i) = A(i);
	}

	for (int i = 0; i < B.size(); ++i) {
		V(i + A.size()) = B(i);
	}
}

void Poelke::toEigenVector(EigenVector& v, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf)
{
	v.resize(3 * m_mesh_3.NumC());

	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		v(3 * i) = vf[ci].x();
		v(3 * i + 1) = vf[ci].y();
		v(3 * i + 2) = vf[ci].z();
	}
}

void Poelke::fromEigenVector(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3>& vf, EigenVector v)
{
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		vf[ci] = Mesh_3_Vector_3(v(3 * i), v(3 * i + 1), v(3 * i + 2));
	}
}

double Poelke::innerProduct(std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf1, std::map<Mesh_3_Cell_iterator, Mesh_3_Vector_3> vf2)
{
	double inner = 0;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);
	
		inner += vf1[ci] * vf2[ci] * m_mesh_3.CellPrimal(ci);
	}

	return inner;
}

double Poelke::cellErrorQuadrature(Mesh_3_Cell_iterator ci, Mesh_3_Vector_3 v, int fieldType)
{
	// quadrature rule: 3/4 1/12 1/12 1/12
	double val = 0;
	
	//for (int i = 0; i < 4; ++i) {
		Mesh_3_Vector_3 p(0, 0, 0);
		for (int j = 0; j < 4; ++j) {
			//if (j == i) {
				p += (1.0 / 4)*(ci->vertex(j)->point().point() - CGAL::ORIGIN);
			//}
			//else {
			//	p += (1.0 / 8)*(ci->vertex(j)->point().point() - CGAL::ORIGIN);
			//}
		}

		Mesh_3_Vector_3 anav;

		if (fieldType == 0) { // normalGradient
			anav = normalGradient(p);
		}
		else if (fieldType == 1) {
			anav = normalGradient(p);
		}
		else if (fieldType == 2) {
			anav = tangentialCurl(p);
		}
		else if (fieldType == 3) { 
			anav = normalHarmonic(p);
		}
		else if (fieldType == 4) { 
			anav = tangentialHarmonic(p);
		}
		else if (fieldType == 5) { 
			anav = centralHarmonic(p);
		}

		Mesh_3_Vector_3 diff = anav - v;
		//cout << "##########" << endl;
		//cout << anav << endl;
		//cout << v << endl;
		//cout << diff << endl;
		val += diff*diff*m_mesh_3.CellPrimal(ci);
	//}
	
	return val;
}

double Poelke::cellInnerProductQuadrature(Mesh_3_Cell_iterator ci, int fieldType)
{
	// quadrature rule: 3/4 1/12 1/12 1/12
	double val = 0;

	//for (int i = 0; i < 4; ++i) {
		Mesh_3_Vector_3 p(0, 0, 0);
		for (int j = 0; j < 4; ++j) {
			//if (j == i) {
				p += (1.0 / 4)*(ci->vertex(j)->point().point() - CGAL::ORIGIN);
			//}
			//else {
				//p += (1.0 / 8)*(ci->vertex(j)->point().point() - CGAL::ORIGIN);
			//}
		}

		Mesh_3_Vector_3 anav;

		if (fieldType == 0) { // normalGradient
			anav = normalGradient(p);
		}
		else if (fieldType == 1) { // normalGradient
			anav = normalGradient(p);
		}
		else if (fieldType == 2) { // tangentialCurl
			anav = tangentialCurl(p);
		}
		else if (fieldType == 3) { // tangentialCurl
			anav = normalHarmonic(p);
		}
		else if (fieldType == 4) { // tangentialCurl
			anav = tangentialHarmonic(p);
		}
		else if (fieldType == 5) { // tangentialCurl
			anav = centralHarmonic(p);
		}

		val += anav*anav*m_mesh_3.CellPrimal(ci);
	//}

	return val;
}

Mesh_3_Vector_3 Poelke::averageNormalGradient(Mesh_3_Facet_iterator fi)
{
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

	Mesh_3_Vector_3 f(3 * p.x() + u.x() + v.x(), 3 * p.y() + u.y() + v.y(), 3 * p.z() + u.z() + v.z());

	return f / 6;
}

Mesh_3_Vector_3 Poelke::averageTangentialCurl(Mesh_3_Facet_iterator fi)
{
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

	Mesh_3_Vector_3 f(3 * p.y() + u.y() + v.y(), -(3 * p.x() + u.x() + v.x()), 0);

	return f / 6;
}

Mesh_3_Vector_3 Poelke::averageCentralHarmonic(Mesh_3_Facet_iterator fi)
{
	return Mesh_3_Vector_3(1, 1, 1);
}

double Poelke::computeFieldError(int fieldType)
{
	double totalError = 0;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);

		if (fieldType == 0) {
			double cer = cellErrorQuadrature(ci, m_original_input[ci], fieldType);
			totalError += cer;
		}
		else if (fieldType == 1) { // normalGradient
			double cer = cellErrorQuadrature(ci, m_normalGradient[ci], fieldType);
			totalError += cer;
		}
		else if (fieldType == 2) { // tangentialCurl
			double cer = cellErrorQuadrature(ci, m_tangentialCurl[ci], fieldType);
			totalError += cer;
		}
		else if (fieldType == 3) { // normalHarmonic
			double cer = cellErrorQuadrature(ci, m_normalHarmonic[ci], fieldType);
			totalError += cer;
		}
		else if (fieldType == 4) { // tangentialHarmonic
			double cer = cellErrorQuadrature(ci, m_tangentialHarmonic[ci], fieldType);
			totalError += cer;
		}
		else if (fieldType == 5) { // centralHarmonic
			double cer = cellErrorQuadrature(ci, m_centralHarmonic[ci], fieldType);
			totalError += cer;
		}
	}

	return sqrt(totalError);
}

double Poelke::computeFieldInnerProduct(int fieldType)
{
	double total = 0;
	for (int i = 0; i < m_mesh_3.NumC(); ++i) {
		Mesh_3_Cell_iterator ci = m_mesh_3.IdxToCell(i);
		double inc = cellInnerProductQuadrature(ci, fieldType);
		total += inc;
	}

	return sqrt(total);
}

void Poelke::testOrthogonalityOfBasis()
{
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);

		if (m_mesh_3.FacetOnBoundary(fi))
			continue;

		int cIdx1 = m_mesh_3.CellIdx(fi->first);
		int cIdx2 = m_mesh_3.CellIdx(fi->first->neighbor(fi->second));

		std::vector<int> ids;
		for (int j = 0; j < 4; ++j) {
			if (j != fi->second)
				ids.push_back(j);
		}

		for (int j = 0; j < ids.size(); ++j) {
			int idx1 = ids[j];
			int idx2 = ids[(j + 1) % ids.size()];

			if (idx1 > idx2)
				std::swap(idx1, idx2);

			My_edge me(fi->first, idx1, idx2);
			Mesh_3_Edge_iterator ei = m_mesh_3.MyCellEdgeMap(me);
			int eIdx = m_mesh_3.EdgeIdx(ei);

			if (m_mesh_3.EdgeOnBoundary(ei))
				continue;

			Mesh_3_Vector_3 c1(m_curlOp_NB.coeff(3 * cIdx1, m_mesh_3.EdgeForward(eIdx)),
				m_curlOp_NB.coeff(3 * cIdx1 + 1, m_mesh_3.EdgeForward(eIdx)),
				m_curlOp_NB.coeff(3 * cIdx1 + 2, m_mesh_3.EdgeForward(eIdx)));
			Mesh_3_Vector_3 c2(m_curlOp_NB.coeff(3 * cIdx2, m_mesh_3.EdgeForward(eIdx)),
				m_curlOp_NB.coeff(3 * cIdx2 + 1, m_mesh_3.EdgeForward(eIdx)),
				m_curlOp_NB.coeff(3 * cIdx2 + 2, m_mesh_3.EdgeForward(eIdx)));

			Mesh_3_Vector_3 g1(m_gradOp_NB.coeff(3 * cIdx1, m_mesh_3.FacetForward(i)),
				m_gradOp_NB.coeff(3 * cIdx1+1, m_mesh_3.FacetForward(i)),
				m_gradOp_NB.coeff(3 * cIdx1+2, m_mesh_3.FacetForward(i)));
			Mesh_3_Vector_3 g2(m_gradOp_NB.coeff(3 * cIdx2, m_mesh_3.FacetForward(i)),
				m_gradOp_NB.coeff(3 * cIdx2 + 1, m_mesh_3.FacetForward(i)),
				m_gradOp_NB.coeff(3 * cIdx2 + 2, m_mesh_3.FacetForward(i)));
		
			double inner = c1 * g1*m_mesh_3.CellPrimal(fi->first) + c2 * g2*m_mesh_3.CellPrimal(fi->first->neighbor(fi->second));
			if (fabs(inner) > 10e-10) {
				cout << "EIP: " << inner << endl;
				cout << c1 << endl;
				cout << c2 << endl;
				cout << g1 << endl;
				cout << g2 << endl;
			}

		}


	}
}

void Poelke::checkCurlFieldBoundary()
{
	for (int i = 0; i < m_mesh_3.NumF(); ++i) {
		Mesh_3_Facet_iterator fi = m_mesh_3.IdxToFacet(i);
	
		if (!m_mesh_3.FacetOnBoundary(fi))
			continue;
	
		Mesh_3_Cell_iterator ci = fi->first;
		if (ci->subdomain_index() == 0)
			ci = ci->neighbor(fi->second);

		Mesh_3_Vector_3 v = m_tangentialCurl[ci];
		Mesh_3_Vector_3 n = m_mesh_3.FacetNormal(fi);

		if ((n * v) / (sqrt(v*v)*sqrt(n*n)) > 10e-14)
			cout << "Orthogonal?: " << (n * v) / (sqrt(v*v)*sqrt(n*n)) << endl;
		
	}
}