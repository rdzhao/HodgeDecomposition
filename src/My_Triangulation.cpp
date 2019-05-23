#include "My_Triangulation.h"

C3t3& My_Triangulation::C3T3()
{
	return m_c3t3;
}

EigenSpMat& My_Triangulation::ED0F_B()
{
	return m_ED0F_B;
}

EigenSpMat& My_Triangulation::ED1F_B()
{
	return m_ED1F_B;
}

EigenSpMat& My_Triangulation::ED2F_B()
{
	return m_ED2F_B;
}

EigenSpMat& My_Triangulation::ED0F_NB()
{
	return m_ED0F_NB;
}

EigenSpMat& My_Triangulation::ED1F_NB()
{
	return m_ED1F_NB;
}

EigenSpMat& My_Triangulation::ED2F_NB()
{
	return m_ED2F_NB;
}

EigenSpMat& My_Triangulation::HS0F_B()
{
	return m_HS0F_B;
}

EigenSpMat& My_Triangulation::HS1F_B()
{
	return m_HS1F_B;
}

EigenSpMat& My_Triangulation::HS2F_B()
{
	return m_HS2F_B;
}

EigenSpMat& My_Triangulation::HS3F()
{
	return m_HS3F;
}

EigenSpMat& My_Triangulation::HS0F_NB()
{
	return m_HS0F_NB;
}

EigenSpMat& My_Triangulation::HS1F_NB()
{
	return m_HS1F_NB;
}

EigenSpMat& My_Triangulation::HS2F_NB()
{
	return m_HS2F_NB;
}

EigenSpMat& My_Triangulation::CHS0F_B()
{
	return m_CHS0F_B;
}

EigenSpMat& My_Triangulation::CHS1F_B()
{
	return m_CHS1F_B;
}

EigenSpMat& My_Triangulation::CHS2F_B()
{
	return m_CHS2F_B;
}

EigenSpMat& My_Triangulation::CHS3F()
{
	return m_CHS3F;
}

EigenSpMat& My_Triangulation::CHS0F_NB()
{
	return m_CHS0F_NB;
}

EigenSpMat& My_Triangulation::CHS1F_NB()
{
	return m_CHS1F_NB;
}

EigenSpMat& My_Triangulation::CHS2F_NB()
{
	return m_CHS2F_NB;
}

int My_Triangulation::NumV()
{
	return numV;
}

int My_Triangulation::NumE()
{
	return numE;
}

int My_Triangulation::NumF()
{
	return numF;
}

int My_Triangulation::NumC()
{
	return numC;
}

int My_Triangulation::NumIV()
{
	return numIV;
}

int My_Triangulation::NumBV()
{
	return numBV;
}

int My_Triangulation::NumIE()
{
	return numIE;
}

int My_Triangulation::NumBE()
{
	return numBE;
}

int My_Triangulation::NumIF()
{
	return numIF;
}

int My_Triangulation::NumBF()
{
	return numBF;
}

int My_Triangulation::VertexIdx(Mesh_3_Vertex_iterator vi)
{
	return vertexIdx[vi];
}

int My_Triangulation::EdgeIdx(Mesh_3_Edge_iterator ei)
{
	return edgeIdx[ei];
}

int My_Triangulation::FacetIdx(Mesh_3_Facet_iterator fi)
{
	return facetIdx[fi];
}

int My_Triangulation::CellIdx(Mesh_3_Cell_iterator ci)
{
	return cellIdx[ci];
}

Mesh_3_Vertex_iterator My_Triangulation::IdxToVertex(int idx)
{
	return idxToVertex[idx];
}

Mesh_3_Edge_iterator My_Triangulation::IdxToEdge(int idx)
{
	return idxToEdge[idx];
}

Mesh_3_Facet_iterator My_Triangulation::IdxToFacet(int idx)
{
	return idxToFacet[idx];
}

Mesh_3_Cell_iterator My_Triangulation::IdxToCell(int idx)
{
	return idxToCell[idx];
}

Mesh_3_Facet_iterator My_Triangulation::MyCellFacetMap(My_facet mf)
{
	return myCellFacetMap[mf];
}

Mesh_3_Edge_iterator My_Triangulation::MyCellEdgeMap(My_edge me)
{
	return myCellEdgeMap[me];
}

bool My_Triangulation::VertexValid(Mesh_3_Vertex_iterator vi)
{
	return vertexValid[vi];
}

bool My_Triangulation::EdgeValid(Mesh_3_Edge_iterator ei)
{
	return edgeValid[ei];
}

bool My_Triangulation::FacetValid(Mesh_3_Facet_iterator fi)
{
	return facetValid[fi];
}

bool My_Triangulation::VertexOnBoundary(Mesh_3_Vertex_iterator vi)
{
	return vertexOnBoundary[vi];
}

bool My_Triangulation::EdgeOnBoundary(Mesh_3_Edge_iterator ei)
{
	return edgeOnBoundary[ei];
}

bool My_Triangulation::FacetOnBoundary(Mesh_3_Facet_iterator fi)
{
	return facetOnBoundary[fi];
}

Mesh_3_Point_3 My_Triangulation::EdgeCC(Mesh_3_Edge_iterator ei)
{
	return edgeCC[ei];
}

Mesh_3_Point_3 My_Triangulation::FacetCC(Mesh_3_Facet_iterator fi)
{
	return facetCC[fi];
}

Mesh_3_Point_3 My_Triangulation::CellCC(Mesh_3_Cell_iterator ci)
{
	return cellCC[ci];
}

double My_Triangulation::VertexDual(Mesh_3_Vertex_iterator vi)
{
	return vertexDual[vi];
}

double My_Triangulation::EdgePrimal(Mesh_3_Edge_iterator ei)
{
	return edgePrimal[ei];
}

double My_Triangulation::EdgeDual(Mesh_3_Edge_iterator ei)
{
	return edgeDual[ei];
}

double My_Triangulation::FacetPrimal(Mesh_3_Facet_iterator fi)
{
	return facetPrimal[fi];
}

double My_Triangulation::FacetDual(Mesh_3_Facet_iterator fi)
{
	return facetDual[fi];
}

double My_Triangulation::CellPrimal(Mesh_3_Cell_iterator ci)
{
	return cellPrimal[ci];
}

Mesh_3_Vector_3 My_Triangulation::FacetNormal(Mesh_3_Facet_iterator fi)
{
	return facetNormal[fi];
}

int My_Triangulation::VertexBoundaryIdx(Mesh_3_Vertex_iterator vi)
{
	return vertexBoundaryIdx[vi];
}

int My_Triangulation::NumBoundary()
{
	return numBoundary;
}

int My_Triangulation::NumHMLG1()
{
	return nhmlg_1;
}

int My_Triangulation::NumHMLG2()
{
	return nhmlg_2;
}

int My_Triangulation::VertexForward(int idx)
{
	return vertexForward[idx];
}

int My_Triangulation::EdgeForward(int idx)
{
	return edgeForward[idx];
}

int My_Triangulation::FacetForward(int idx)
{
	return facetForward[idx];
}

int My_Triangulation::VertexBackward(int idx)
{
	return vertexBackward[idx];
}

int My_Triangulation::EdgeBackward(int idx)
{
	return edgeBackward[idx];
}

int My_Triangulation::FacetBackward(int idx)
{
	return facetBackward[idx];
}

double My_Triangulation::getRadius()
{
	return m_radius;
}

Mesh_3_Point_3 My_Triangulation::getCenter()
{
	return center;
}

double My_Triangulation::getAveLen()
{
	return m_aveLen;
}

Mesh_3_Point_3 My_Triangulation::BoundingBoxMin()
{
	return bbMin;
}

Mesh_3_Point_3 My_Triangulation::BoundingBoxMax()
{
	return bbMax;
}

void My_Triangulation::processRawData()
{
	//// get grid dimention and width
	//double xmin, ymin, zmin, xmax, ymax, zmax;
	//unsigned int rx, ry, rz;
	//double gs;
	//std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
	//std::ifstream dataFile("proteinwhole.vox");
	//dataFile >> xmin >> ymin >> zmin;
	//dataFile >> xmax >> ymax >> zmax;
	//dataFile >> rx >> ry >> rz;
	//gs = (xmax - xmin) / rx;

	//voxels.reserve(rx*ry*rz);

	//std::string str;
	//while (std::getline(dataFile, str)) {
	//	float val;
	//	std::stringstream ss(str);
	//	ss >> val;
	//	
	//	float t = val;
	//	//double t = (val + 1) / 2;
	//	//if (t < 0)
	//	//	t = 0;
	//	//if (t > 1)
	//	//	t = 1;

	//	voxels.push_back(t);
	//}
	//
	//_image* image = ::_initImage();
	//image->vectMode = VM_SCALAR;
	//image->xdim = rx;
	//image->ydim = ry;
	//image->zdim = rz;
	//image->vdim = 1;
	//image->vx = gs;
	//image->vy = gs;
	//image->vz = gs;
	//image->endianness = ::_getEndianness();
	//image->wdim = sizeof(float);
	////std::cout << "Size of float: " << sizeof(float) << std::endl;
	//image->wordKind = WK_FLOAT;
	//image->sign = SGN_SIGNED;
	//image->data = &voxels[0];

	////ofstream test("test.txt");
	////int* ptr = (int*)image->data;
	////for (int i = 0; i < voxels.size(); ++i)
	////	test << ptr[i] << endl;

	//m_image = CGAL::Image_3(image);
}

void My_Triangulation::setCriteria(double fa, double fs, double fd, double cr, double cs)
{
	m_facet_angle = fa;
	m_facet_size = fs;
	m_facet_distance = fd;
	m_cell_radius_edge_ratio = cr;
	m_cell_size = cs;

	//m_facet_angle = 30;
	//m_facet_size = 2;
	//m_facet_distance = 2;
	//m_cell_radius_edge_ratio = 2;
	//m_cell_size = 2;
}

void My_Triangulation::buildTetMesh()
{
	//Polyhedron polyhedron;
	//std::ifstream input("models\\sphere.off");
	//input >> polyhedron;
	//input.close();
	//Polyhedron_mesh_domain domain(polyhedron);

	//Mesh_criteria criteria(facet_angle = m_facet_angle,
	//	facet_size = m_facet_size,
	//	facet_distance = m_facet_distance,
	//	cell_radius_edge_ratio = m_cell_radius_edge_ratio,
	//	cell_size = m_cell_size);
	//std::cout<<"Criteria ... "<< std::endl;

	//m_c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::odt());
	//CGAL::perturb_mesh_3(m_c3t3, domain, time_limit = 10);
	//std::cout << "Triangulation ... " << std::endl;


	//Polyhedron_with_feature polyhedron;
	//std::ifstream input("models\\sphere.off");
	//input >> polyhedron;
	//input.close();
	//Mesh_domain_with_feature domain(polyhedron);
	//domain.detect_features();

	//Mesh_criteria criteria(facet_angle = m_facet_angle,
	//	facet_size = m_facet_size,
	//	facet_distance = m_facet_distance,
	//	cell_radius_edge_ratio = m_cell_radius_edge_ratio,
	//	cell_size = m_cell_size);
	//std::cout << "Criteria ... " << std::endl;

	////r_c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::odt());
	////CGAL::perturb_mesh_3(r_c3t3, domain, time_limit = 10);
	//m_c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
	//std::cout << "Triangulation ... " << std::endl;
}

void My_Triangulation::preprocessing()
{
	computeValidAndBoundary();
	indexElements();
	basicInit();
	
	buildMyCellFacetEdgeMap();
	computeBoundaryComponent();
	computeHomology();
	createForwardBackwardMap();

	computeElementCircumcenter();
	computeElementVolume();

	computeFacetNormal();

	buildDECOperators();
}

void My_Triangulation::buildDECOperators()
{
	buildHS0F_B();
	buildHS1F_B();
	buildHS2F_B();
	buildHS3F();
	buildHS0F_NB();
	buildHS1F_NB();
	buildHS2F_NB();

	buildED0F_B();
	buildED1F_B();
	buildED2F_B();
	buildED0F_NB();
	buildED1F_NB();
	buildED2F_NB();

	buildCHS0F_B();
	buildCHS1F_B();
	buildCHS2F_B();
	buildCHS3F();
	buildCHS0F_NB();
	buildCHS1F_NB();
	buildCHS2F_NB();
}

void My_Triangulation::computeElementCircumcenter()
{
	std::cout << "Circumcenter ..." << std::endl;
	computeEdgeCircumcenter();
	computeFacetCircumcenter();
	computeCellCircumcenter();
}

void My_Triangulation::computeElementVolume()
{
	std::cout<<"Element Volume ..."<<std::endl;
	computeEdgeVolumes();
	computeFacetVolumes();
	computeCellVolumes();
	computeVertexVolumes(); // dependent on cell volume
}

void My_Triangulation::clean()
{

}

void My_Triangulation::indexElements()
{
	int ii;

	ii = 0;
	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!vertexValid[vi])
			continue;

		vertexIdx[vi] = ii;
		idxToVertex[ii] = vi;
		++ii;
	}
	numV = ii;

	ii = 0;
	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;
		
		edgeIdx[ei] = ii;
		idxToEdge[ii] = ei;
		++ii;
	}
	numE = ii;

	ii = 0;
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		facetIdx[fi] = ii;
		idxToFacet[ii] = fi;
		++ii;
	}
	numF = ii;

	ii = 0;
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;
		
		if (ci->subdomain_index() != 1)
			cout << ci->subdomain_index() << endl;

		cellIdx[ci] = ii;
		idxToCell[ii] = ci;
		++ii;
		//cout << ci->subdomain_index() << endl;
	}
	numC = ii;

	cout << "NumC: " << numC << endl;

}

void My_Triangulation::basicInit()
{
	// center and radius
	Mesh_3_Point_3 p(0, 0, 0);
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!vertexValid[vi])
			continue;

		p += (vi->point().point() - CGAL::ORIGIN) / static_cast<int>(m_c3t3.triangulation().number_of_vertices());
	
		if (vertexIdx[vi] == 0) {
			xmin = vi->point().point().x();
			ymin = vi->point().point().y();
			zmin = vi->point().point().z();
			xmax = vi->point().point().x();
			ymax = vi->point().point().y();
			zmax = vi->point().point().z();
		}
		else {
			if (xmin > vi->point().point().x())
				xmin = vi->point().point().x();
			if (ymin > vi->point().point().y())
				ymin = vi->point().point().y();
			if (zmin > vi->point().point().z())
				zmin = vi->point().point().z();
			if (xmax < vi->point().point().x())
				xmax = vi->point().point().x();
			if (ymax < vi->point().point().y())
				ymax = vi->point().point().y();
			if (zmax < vi->point().point().z())
				zmax = vi->point().point().z();
		}
	}
	center = p;
	m_radius = (xmax + ymax + zmax - xmin - ymin - zmin) / 6;
	bbMin = Mesh_3_Point_3(xmin, ymin, zmin);
	bbMax = Mesh_3_Point_3(xmax, ymax, zmax);

	// average edge length
	double tl = 0;
	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;
		
		Mesh_3_Vector_3 v = ei->first->vertex(ei->second)->point().point() - ei->first->vertex(ei->third)->point().point();

		tl += sqrt(v*v);
	}
	m_aveLen = tl / numE;

}

void My_Triangulation::buildMyCellFacetEdgeMap() {
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		My_facet mf1(fi->first, fi->second);
		myCellFacetMap[mf1] = fi;
		
		Facet opf = m_c3t3.triangulation().mirror_facet(*fi);
		My_facet mf2(opf.first, opf.second);
		myCellFacetMap[mf2] = fi;

		//std::cout << cellIdx[fi->first] << " " << cellIdx[opf.first] << " " << fi->second << " " << opf.second << std::endl;
	}

	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		Mesh_3_Cell_circulator cc = m_c3t3.triangulation().incident_cells(*ei);
		Mesh_3_Cell_circulator end = cc;
		Mesh_3_Vertex_iterator vi1, vi2;
		vi1 = ei->first->vertex(ei->second);
		vi2 = ei->first->vertex(ei->third);

		do {
			int i, j;

			cc->has_vertex(vi1, i);
			cc->has_vertex(vi2, j);

			if (i > j)
				std::swap(i, j);

			My_edge me(cc, i, j);
			myCellEdgeMap[me] = ei;

			++cc;
		} while (cc != end);
	}
}

void My_Triangulation::computeEdgeCircumcenter()
{
	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		Mesh_3_Point_3 p1 = ei->first->vertex(ei->second)->point().point();
		Mesh_3_Point_3 p2 = ei->first->vertex(ei->third)->point().point();

		Mesh_3_Point_3 p = CGAL::circumcenter(p1, p2);
		edgeCC[ei] = p;
	}
}

void My_Triangulation::computeFacetCircumcenter()
{
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		std::vector<Mesh_3_Point_3> vec;
		for (int i = 0; i < 4; ++i) {
			if (fi->second != i) {
				vec.push_back(fi->first->vertex(i)->point().point());
			}
		}

		Mesh_3_Point_3 p = CGAL::circumcenter(vec[0], vec[1], vec[2]);
		facetCC[fi] = p;
	}
}

void My_Triangulation::computeCellCircumcenter()
{
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		std::vector<Mesh_3_Point_3> vec;
		for (int i = 0; i < 4; ++i) {
			vec.push_back(ci->vertex(i)->point().point());
		}

		Mesh_3_Point_3 p = CGAL::circumcenter(vec[0], vec[1], vec[2], vec[3]);
		cellCC[ci] = p;

		//// test whether circumcenter in tet.
		//double a, a1, a2, a3, a4;
		//a = CGAL::volume(vec[0], vec[1], vec[2], vec[3]);
		//a1 = CGAL::volume(p, vec[1], vec[2], vec[3]);
		//a2 = CGAL::volume(vec[0], p, vec[2], vec[3]);
		//a3 = CGAL::volume(vec[0], vec[1], p, vec[3]);
		//a4 = CGAL::volume(vec[0], vec[1], vec[2], p);
		////std::cout << a << " " << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
		//if (a1 < -0.0001 || a4 < -0.0001 || a4 < -0.0001 || a4 < -0.0001){
		//	std::cout << a << " " << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
		//}
	}
}

void My_Triangulation::computeVertexVolumes()
{
	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!vertexValid[vi])
			continue;

		vertexDual[vi] = 0;

		/*
		std::vector<Mesh_3_Cell_iterator> cells;
		m_c3t3.triangulation().incident_cells(vi, std::back_inserter(cells));

		for (int i = 0; i < cells.size(); ++i) {
			if (cells[i]->subdomain_index() == 0)
				continue;

			// compute volume
			vertexDual[vi] += cellPrimal[cells[i]] / 4;
		}
		*/

		std::vector<Edge> edges;
		m_c3t3.triangulation().incident_edges(vi, std::back_inserter(edges));

		for (int i = 0; i < edges.size(); ++i) {
			Mesh_3_Cell_iterator ci = edges[i].first;
			int idx1 = edges[i].second;
			int idx2 = edges[i].third;
			if (idx1 > idx2)
				std::swap(idx1, idx2);

			My_edge me(ci, idx1, idx2);
			Mesh_3_Edge_iterator ei = myCellEdgeMap[me];

			vertexDual[vi] += edgeDual[ei] * edgePrimal[ei] / 6;
		}
	}
}

void My_Triangulation::computeEdgeVolumes()
{
	//primal
	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		Mesh_3_Vector_3 v = ei->first->vertex(ei->second)->point().point() - ei->first->vertex(ei->third)->point().point();
		edgePrimal[ei] = sqrt(v*v);
		//std::cout<< edgePrimal[ei] <<std::endl;
	}

	//dual
	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		std::vector<Mesh_3_Point_3> pVec;
		Mesh_3_Cell_circulator cc = m_c3t3.triangulation().incident_cells(*ei);
		Mesh_3_Cell_circulator end = cc;

		edgeDual[ei] = 0;
		do {
			if (cc->subdomain_index() != 0) {
				int idx1, idx2;
				bool b1= cc->has_vertex(ei->first->vertex(ei->second), idx1);
				bool b2 = cc->has_vertex(ei->first->vertex(ei->third), idx2);
				//std::cout << b1 << " " << b2 << std::endl;

				std::vector<Mesh_3_Point_3> fcv;
				for (int i = 0; i < 4; ++i) {
					if (i != idx1 && i != idx2) {
						My_facet mf(cc, i);
						fcv.push_back(facetCC[myCellFacetMap[mf]]);
					}
				}

				Mesh_3_Point_3 cCenter = cellCC[cc];

				Mesh_3_Vector_3 vfc1, vfc2, vcc;
				vfc1 = fcv[0] - edgeCC[ei];
				vfc2 = fcv[1] - edgeCC[ei];
				vcc = cCenter - edgeCC[ei];

				Mesh_3_Vector_3 av1 = CGAL::cross_product(vfc1, vcc);
				Mesh_3_Vector_3 av2 = CGAL::cross_product(vfc2, vcc);


				edgeDual[ei] += sqrt(av1*av1) / 2;
				edgeDual[ei] += sqrt(av2*av2) / 2;
			}
			++cc;
		} while (cc != end);

		//std::cout << edgeDual[ei] << std::endl;
	}
}

void My_Triangulation::computeFacetVolumes()
{
	// primal
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		std::vector<Mesh_3_Point_3> pVec;
		for (int i = 0; i < 4; ++i) {
			if (i != fi->second) {
				pVec.push_back(fi->first->vertex(i)->point().point());
			}
		}

		std::vector<Mesh_3_Vector_3> vVec;
		vVec.push_back(pVec[1] - pVec[0]);
		vVec.push_back(pVec[2] - pVec[0]);

		Mesh_3_Vector_3 av = CGAL::cross_product(vVec[0], vVec[1]);

		facetPrimal[fi] = sqrt(av*av) / 2;
	}

	// dual
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		facetDual[fi] = 0;
		Mesh_3_Point_3 fc = facetCC[fi];
		
		Mesh_3_Cell_iterator ci1 = fi->first;
		Mesh_3_Cell_iterator ci2 = fi->first->neighbor(fi->second); /*m_c3t3.triangulation().mirror_facet(*fi).first*/;
		if (ci1->subdomain_index() != 0) {
			Mesh_3_Vector_3 v1 = cellCC[ci1] - fc;
			facetDual[fi] += sqrt(v1*v1);
		}
		if (ci2->subdomain_index() != 0) {
			Mesh_3_Vector_3 v2 = cellCC[ci2] - fc;
			facetDual[fi] += sqrt(v2*v2);
		}

		//Mesh_3_Point_3 p1, p2;
		//p1 = cellCC[fi->first];
		//p2 = cellCC[fi->first->neighbor(fi->second)];

		//facetDual[fi] = sqrt((p1 - p2)*(p1 - p2));
	}
}

void My_Triangulation::computeCellVolumes()
{
	for (Mesh_3_Cell_iterator ci = m_c3t3.cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		Mesh_3_Point_3 p0, p1, p2, p3;
		p0 = ci->vertex(0)->point().point();
		p1 = ci->vertex(1)->point().point();
		p2 = ci->vertex(2)->point().point();
		p3 = ci->vertex(3)->point().point();

		Mesh_3_Vector_3 v0, v1, v2;
		v0 = p0 - p3;
		v1 = p1 - p3;
		v2 = p2 - p3;

		Mesh_3_Vector_3 av = CGAL::cross_product(v0, v1);

		cellPrimal[ci] = fabs(av*v2) / 6;
	}
}

void My_Triangulation::computeFacetNormal()
{
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		std::vector<int> pIdx;
		for (int i = 0; i < 4; ++i) {
			if (i != fi->second) {
				pIdx.push_back(vertexIdx[fi->first->vertex(i)]);
			}
		}

		std::sort(pIdx.begin(), pIdx.end());

		Mesh_3_Vector_3 v1 = idxToVertex[pIdx[1]]->point().point() - idxToVertex[pIdx[0]]->point().point();
		Mesh_3_Vector_3 v2 = idxToVertex[pIdx[2]]->point().point() - idxToVertex[pIdx[0]]->point().point();

		Mesh_3_Vector_3 n = CGAL::cross_product(v1, v2);
		n /= sqrt(n*n);

		facetNormal[fi] = n;
	}
}

void My_Triangulation::createForwardBackwardMap()
{
	numIV = 0;
	numBV = 0;
	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!m_c3t3.triangulation().is_infinite(vi)){
			if (vertexValid[vi]) {
				if (!vertexOnBoundary[vi]) {
					vertexForward[vertexIdx[vi]] = numIV;
					vertexBackward[numIV] = vertexIdx[vi];
					++numIV;
				}
				else
					++numBV;
			}
		}
	}

	numIE = 0;
	numBE = 0;
	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (edgeValid[ei]) {
			if (!edgeOnBoundary[ei]) {
				edgeForward[edgeIdx[ei]] = numIE;
				edgeBackward[numIE] = edgeIdx[ei];
				++numIE;
			}
			else
				++numBE;
		}
	}

	numIF = 0;
	numBF = 0;
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (facetValid[fi]) {
			if (!facetOnBoundary[fi]) {
				facetForward[facetIdx[fi]] = numIF;
				facetBackward[numIF] = facetIdx[fi];
				++numIF;
			}
			else
				++numBF;
		}
	}

	std::cout << "Vertex: " << numIV << " " << numBV << std::endl;
	std::cout << "Edge: " << numIE << " " << numBE << std::endl;
	std::cout << "Facet: " << numIF << " " << numBF << std::endl;
}

void My_Triangulation::computeValidAndBoundary()
{
	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (m_c3t3.triangulation().is_infinite(vi)) {
			vertexValid[vi] = false;
			vertexOnBoundary[vi] = false;
		}
		else {
			std::vector<Mesh_3_Cell_iterator> cells;
			m_c3t3.triangulation().incident_cells(vi, std::back_inserter(cells));

			bool valid = false;
			bool hasOutside = false;
			for (int i = 0; i < cells.size(); ++i) {
				if (cells[i]->subdomain_index() != 0)
					valid = true;
				else
					hasOutside = true;
			}
			vertexValid[vi] = valid;
			vertexOnBoundary[vi] = (valid && hasOutside);
		}
	}

	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (fi->first->subdomain_index() == 0 && fi->first->neighbor(fi->second)->subdomain_index() == 0) {
			facetValid[fi] = false;
			facetOnBoundary[fi] = false;
		}
		else {
			facetValid[fi] = true;
			if (fi->first->subdomain_index() != fi->first->neighbor(fi->second)->subdomain_index()) {
				facetOnBoundary[fi] = true;
			
				//for (int i = 0; i < 4; ++i) {
					//if (i != fi->second) {
						//vertexOnBoundary[fi->first->vertex(i)] = true;
					//}
				//}
			}
			else {
				facetOnBoundary[fi] = false;
			}
		}
	}

	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		//Mesh_3_Vertex_iterator v1, v2;
		//v1 = ei->first->vertex(ei->second);
		//v2 = ei->first->vertex(ei->third);

		//if (!vertexValid[v1] || !vertexValid[v2]) {
		//	edgeValid[ei] = false;
		//	edgeOnBoundary[ei] = false;
		//}
		//else
		//{
		//	edgeValid[ei] = true;
		//	if (vertexOnBoundary[v1] && vertexOnBoundary[v2])
		//		edgeOnBoundary[ei] = true;
		//	else
		//		edgeOnBoundary[ei] = false;
		//}

		Mesh_3_Cell_circulator cc = m_c3t3.triangulation().incident_cells(*ei);
		Mesh_3_Cell_circulator end = cc;

		bool valid = false;
		bool hasOutside = false;
		do {

			if (cc->subdomain_index() != 0)
				valid = true;
			else
				hasOutside = true;
			
			++cc;
		} while (cc!=end);

		edgeValid[ei] = valid;
		edgeOnBoundary[ei] = (valid && hasOutside);
	}

	//for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
	//	ei != m_c3t3.triangulation().edges_end(); ++ei) {
	//	if (!edgeValid[ei])
	//		continue;

	//	//std::cout<< edgeOnBoundary[ei] <<" "<< vertexOnBoundary[ei->first->vertex(ei->second)] <<" "<< vertexOnBoundary[ei->first->vertex(ei->third)] <<std::endl;
	//	//if (edgeOnBoundary[ei]) {
	//		//if (!vertexOnBoundary[ei->first->vertex(ei->second)] || !vertexOnBoundary[ei->first->vertex(ei->third)])
	//			//std::cout<<"WTF!"<<std::endl;
	//	//}
	//}

}

void My_Triangulation::computeBoundaryComponent()
{
	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!vertexValid[vi])
			continue;
		
		vertexBoundaryIdx[vi] = -1;
	}

	std::map<Mesh_3_Facet_iterator, bool, Mesh_3_Facet_Iterator_Comparator> facetVisited;
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (facetOnBoundary[fi]) {
			facetVisited[fi] = false;
		}
	}

	std::cout<<"Classifying boundary"<< std::endl;

	int seedIdx = 0;
	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi] || !facetOnBoundary[fi])
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
						Edge e = *myCellEdgeMap[me];
						
						// circulate e for adj faces
						Mesh_3_Facet_circulator fc = m_c3t3.triangulation().incident_facets(e);
						Mesh_3_Facet_circulator end = fc;
						do {
							if (m_c3t3.subdomain_index(fc->first)
								!= m_c3t3.subdomain_index(fc->first->neighbor(fc->second))
								&& *fc != *fii) {
								My_facet mf(fc->first, fc->second);
								
								if (!facetVisited[myCellFacetMap[mf]]) {
									fstack.push(myCellFacetMap[mf]);
									facetVisited[myCellFacetMap[mf]] = true;
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

	numBoundary = seedIdx;
	std::cout << "End " << seedIdx << std::endl;
}

void My_Triangulation::computeHomology()
{
	nhmlg_2 = numBoundary - 1;
	nhmlg_1 = numBoundary - (numV - numE + numF - numC);

	//nhmlg_2 = 4;
	//nhmlg_1 = 0;

	cout <<"Homology: "<< nhmlg_1 << " " << nhmlg_2 << endl;
}

void My_Triangulation::buildHS0F_B()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!vertexValid[vi])
			continue;

		triplet.push_back(Eigen::Triplet<double>(vertexIdx[vi],
			vertexIdx[vi],
			vertexDual[vi]));
	}
	m_HS0F_B.resize(numV, numV);
	m_HS0F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildHS1F_B()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		triplet.push_back(Eigen::Triplet<double>(edgeIdx[ei],
			edgeIdx[ei],
			edgeDual[ei] / edgePrimal[ei]));

	}
	m_HS1F_B.resize(numE, numE);
	m_HS1F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildHS2F_B()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		triplet.push_back(Eigen::Triplet<double>(facetIdx[fi],
			facetIdx[fi],
			facetDual[fi] / facetPrimal[fi]));
	}

	m_HS2F_B.resize(numF, numF); // !!
	m_HS2F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildHS3F()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		triplet.push_back(Eigen::Triplet<double>(cellIdx[ci],
			cellIdx[ci],
			1 / cellPrimal[ci]));
	}

	m_HS3F.resize(numC, numC);
	m_HS3F.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildHS0F_NB()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
		vi != m_c3t3.triangulation().vertices_end(); ++vi) {
		if (!vertexValid[vi])
			continue;
		if (vertexOnBoundary[vi])
			continue;

		triplet.push_back(Eigen::Triplet<double>(vertexForward[vertexIdx[vi]],
			vertexForward[vertexIdx[vi]],
			vertexDual[vi]));
	}
	m_HS0F_NB.resize(numIV, numIV);
	m_HS0F_NB.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildHS1F_NB()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		if (!edgeOnBoundary[ei])
			triplet.push_back(Eigen::Triplet<double>(edgeForward[edgeIdx[ei]],
				edgeForward[edgeIdx[ei]],
				edgeDual[ei] / edgePrimal[ei]));

	}
	m_HS1F_NB.resize(numIE, numIE);
	m_HS1F_NB.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildHS2F_NB()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		if (!facetOnBoundary[fi]) // !!
			triplet.push_back(Eigen::Triplet<double>(facetForward[facetIdx[fi]],
				facetForward[facetIdx[fi]],
				facetDual[fi] / facetPrimal[fi]));
	}

	m_HS2F_NB.resize(numIF, numIF); // !!
	m_HS2F_NB.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildCHS0F_B()
{
	std::map<std::pair<int, int>, double> tripmap;
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;
	
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				int idx1 = vertexIdx[ci->vertex(i)];
				int idx2 = vertexIdx[ci->vertex(j)];

				if (idx1 > idx2)
					std::swap(idx1, idx2);
				std::pair<int, int> p(idx1, idx2);
				
				tripmap[p] = 0;
			}
		}
	}
	
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;
		
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				int idx1 = vertexIdx[ci->vertex(i)];
				int idx2 = vertexIdx[ci->vertex(j)];

				if (idx1 > idx2)
					std::swap(idx1, idx2);
				std::pair<int, int> p(idx1, idx2);

				if (idx1 == idx2) {
					tripmap[p] += (1.0 / 10)*cellPrimal[ci];
				}
				else {
					tripmap[p] += (1.0 / 20)*cellPrimal[ci];
				}
			}
		}
	}

	std::vector<Eigen::Triplet<double>> triplet;
	for (auto iter = tripmap.begin(); iter != tripmap.end(); ++iter) {
		if (iter->first.first == iter->first.second) {
			triplet.push_back(Eigen::Triplet<double>(iter->first.first, iter->first.second, iter->second));
		}
		else{
			triplet.push_back(Eigen::Triplet<double>(iter->first.first, iter->first.second, iter->second / 2));
			triplet.push_back(Eigen::Triplet<double>(iter->first.second, iter->first.first, iter->second / 2));
		}
	}

	m_CHS0F_B.resize(numV, numV);
	m_CHS0F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildCHS1F_B()
{
	std::map<std::pair<int, int>, double> tripmap;
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;
	
		// build a simple map to index 6 edges
		std::map<int, My_edge> localmap;
		localmap[0] = My_edge(ci, 0, 1);
		localmap[1] = My_edge(ci, 0, 2);
		localmap[2] = My_edge(ci, 0, 3);
		localmap[3] = My_edge(ci, 1, 2);
		localmap[4] = My_edge(ci, 1, 3);
		localmap[5] = My_edge(ci, 2, 3);
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				My_edge me1 = localmap[i];
				My_edge me2 = localmap[j];

				Mesh_3_Edge_iterator ei1 = myCellEdgeMap[me1];
				Mesh_3_Edge_iterator ei2 = myCellEdgeMap[me2];

				int eidx1 = edgeIdx[ei1];
				int eidx2 = edgeIdx[ei2];
			
				if (eidx1 > eidx2)
					std::swap(eidx1, eidx2);
				std::pair<int, int> p(eidx1, eidx2);

				tripmap[p] = 0;
			}
		}
	}
	
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		// height
		std::map<My_facet, Mesh_3_Vector_3> grad;
		for (int i = 0; i < 4; ++i) {
			My_facet mf(ci, i);

			// magnitude
			double h = 3 * CellPrimal(ci) / FacetPrimal(MyCellFacetMap(mf));

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
			//cout<< n / h <<endl;
		}

		// build a simple map to index 6 edges
		std::map<int, My_edge> localmap;
		localmap[0] = My_edge(ci, 0, 1);
		localmap[1] = My_edge(ci, 0, 2);
		localmap[2] = My_edge(ci, 0, 3);
		localmap[3] = My_edge(ci, 1, 2);
		localmap[4] = My_edge(ci, 1, 3);
		localmap[5] = My_edge(ci, 2, 3);
		//cout<<"WTF ... "<<endl;
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				My_edge me0 = localmap[i];
				My_edge me1 = localmap[j];

				int l0 = std::get<1>(me0);
				int l1 = std::get<2>(me0);
				int l2 = std::get<1>(me1);
				int l3 = std::get<2>(me1);
				//cout << l0 << " " << l1 << " " << l2 << " " << l3 << endl;
				// get correct orientation
				if (VertexIdx(ci->vertex(l0)) > VertexIdx(ci->vertex(l1)))
					std::swap(l0, l1);
				if (VertexIdx(ci->vertex(l2)) > VertexIdx(ci->vertex(l3)))
					std::swap(l2, l3);

				My_facet mf0 = My_facet(ci, l0);
				My_facet mf1 = My_facet(ci, l1);
				My_facet mf2 = My_facet(ci, l2);
				My_facet mf3 = My_facet(ci, l3);

				double val = 0;

				if (l0 == l2) {
					val += (1.0 / 10)*(grad[mf1] * grad[mf3]);
				}
				else {
					val += (1.0 / 20)*(grad[mf1] * grad[mf3]);
				}

				if (l0 == l3) {
					val -= (1.0 / 10)*(grad[mf1] * grad[mf2]);
				}
				else{
					val -= (1.0 / 20)*(grad[mf1] * grad[mf2]);
				}

				if (l1 == l2) {
					val -= (1.0 / 10)*(grad[mf0] * grad[mf3]);
				}
				else {
					val -= (1.0 / 20)*(grad[mf0] * grad[mf3]);
				}

				if (l1 == l3) {
					val += (1.0 / 10)*(grad[mf0] * grad[mf2]);
				}
				else {
					val += (1.0 / 20)*(grad[mf0] * grad[mf2]);
				}

				// i-th edge on j-th edge ----- matrix (i, j)
				Mesh_3_Edge_iterator ei1 = myCellEdgeMap[me0];
				Mesh_3_Edge_iterator ei2 = myCellEdgeMap[me1];

				int eidx1 = edgeIdx[ei1];
				int eidx2 = edgeIdx[ei2];

				if (eidx1 > eidx2)
					std::swap(eidx1, eidx2);
				std::pair<int, int> p(eidx1, eidx2);
				tripmap[p] += val * CellPrimal(ci);
				
				//Mesh_3_Edge_iterator ei0 = MyCellEdgeMap(me0);
				//Mesh_3_Edge_iterator ei1 = MyCellEdgeMap(me1);
				//cout << EdgeIdx(ei0) << " " << EdgeIdx(ei1) << " " << val << endl;
				//m_CHS1F_B.coeffRef(EdgeIdx(ei0), EdgeIdx(ei1)) += val * CellPrimal(ci);
			}
		}
	}

	std::vector<Eigen::Triplet<double>> triplet;
	for (auto iter = tripmap.begin(); iter != tripmap.end(); ++iter) {
		if (iter->first.first == iter->first.second) {
			triplet.push_back(Eigen::Triplet<double>(iter->first.first, iter->first.second, iter->second));
		}
		else {
			triplet.push_back(Eigen::Triplet<double>(iter->first.first, iter->first.second, iter->second / 2));
			triplet.push_back(Eigen::Triplet<double>(iter->first.second, iter->first.first, iter->second / 2));
		}
	}

	m_CHS1F_B.resize(numE, numE);
	m_CHS1F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildCHS2F_B()
{
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

		std::map<int, Mesh_3_Vector_3> grad;
		for (int i = 0; i < 4; ++i) {
			My_facet mf(ci, i);

			// magnitude
			double h = 3 * CellPrimal(ci) / FacetPrimal(MyCellFacetMap(mf));

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

			grad[vertexIdx[ci->vertex(i)]] = n / h;
			//cout<< n / h <<endl;
		}

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				My_facet mfi(ci, i);
				My_facet mfj(ci, j);

				Mesh_3_Facet_iterator fi = myCellFacetMap[mfi];
				Mesh_3_Facet_iterator fj = myCellFacetMap[mfj];
			
				double in = innerProductInCellHS2F(ci, fi, fj, grad);
				in *= CellPrimal(ci);

				if (i == j) {
					triplet.push_back(Eigen::Triplet<double>(facetIdx[fi], facetIdx[fj], in));
				}
				else{
					triplet.push_back(Eigen::Triplet<double>(facetIdx[fi], facetIdx[fj], in / 2));
					triplet.push_back(Eigen::Triplet<double>(facetIdx[fj], facetIdx[fi], in / 2));
				}
			}
		}
	}
	
	m_CHS2F_B.resize(numF, numF);
	m_CHS2F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildCHS3F()
{
	std::vector<Eigen::Triplet<double>> triplet;
	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;
		
		triplet.push_back(Eigen::Triplet<double>(cellIdx[ci], cellIdx[ci], cellPrimal[ci]));
	}

	m_CHS3F.resize(numC, numC);
	m_CHS3F.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildCHS0F_NB()
{
	std::vector<Eigen::Triplet<double>> triplet;
	for (int i = 0; i < m_CHS0F_B.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_CHS0F_B, i); iter; ++iter) {
			double row = iter.row();
			double col = iter.col();
			double val = iter.value();
		
			Mesh_3_Vertex_iterator vir = idxToVertex[row];
			Mesh_3_Vertex_iterator vic = idxToVertex[col];

			if (!vertexOnBoundary[vir] && !vertexOnBoundary[vic]) {
				double nrow = vertexForward[row];
				double ncol = vertexForward[col];

				triplet.push_back(Eigen::Triplet<double>(nrow, ncol, val));
			}
		}
	}
	m_CHS0F_NB.resize(numIV, numIV);
	m_CHS0F_NB.setFromTriplets(triplet.begin(), triplet.end());


	//m_CHS0F_NB.resize(numIV, numIV);

	//for (Mesh_3_Vertex_iterator vi = m_c3t3.triangulation().vertices_begin();
	//	vi != m_c3t3.triangulation().vertices_end(); ++vi) {
	//	if (!vertexValid[vi])
	//		continue;
	//	if (vertexOnBoundary[vi])
	//		continue;

	//	std::vector<Mesh_3_Vertex_iterator> verts;
	//	m_c3t3.triangulation().incident_vertices(vi, std::back_inserter(verts));

	//	for (int i = 0; i < verts.size(); ++i) {
	//		if (!vertexValid[verts[i]])
	//			continue;
	//		if (vertexOnBoundary[verts[i]])
	//			continue;

	//		m_CHS0F_B.coeffRef(vertexForward[vertexIdx[vi]], vertexForward[vertexIdx[verts[i]]]) = 1;
	//	}

	//	m_CHS0F_B.coeffRef(vertexIdx[vi], vertexIdx[vi]) = verts.size() + 1;
	//}
}

void My_Triangulation::buildCHS1F_NB()
{
	std::vector<Eigen::Triplet<double>> triplet;
	for (int i = 0; i < m_CHS1F_B.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_CHS1F_B, i); iter; ++iter) {
			double row = iter.row();
			double col = iter.col();
			double val = iter.value();

			Mesh_3_Edge_iterator eir = idxToEdge[row];
			Mesh_3_Edge_iterator eic = idxToEdge[col];

			if (!edgeOnBoundary[eir] && !edgeOnBoundary[eic]) {
				double nrow = edgeForward[row];
				double ncol = edgeForward[col];

				triplet.push_back(Eigen::Triplet<double>(nrow, ncol, val));
			}
		}
	}
	m_CHS1F_NB.resize(numIE, numIE);
	m_CHS1F_NB.setFromTriplets(triplet.begin(), triplet.end());

	//m_CHS1F_NB.resize(numIE, numIE);

	//for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
	//	ei != m_c3t3.triangulation().edges_end(); ++ei) {
	//	if (!edgeValid[ei])
	//		continue;
	//	if (edgeOnBoundary[ei])
	//		continue;

	//	int idx = edgeForward[edgeIdx[ei]];
	//	m_CHS1F_NB.coeffRef(idx, idx) = 0;
	//}

	//for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
	//	fi != m_c3t3.triangulation().facets_end(); ++fi) {
	//	if (!facetValid[fi])
	//		continue;
	//	if (facetOnBoundary[fi])
	//		continue;

	//	std::vector<int> localIdx;
	//	for (int i = 0; i < 4; ++i) {
	//		if (i == fi->second)
	//			continue;

	//		localIdx.push_back(i);
	//	}

	//	My_edge me1(fi->first, localIdx[0], localIdx[1]);
	//	My_edge me2(fi->first, localIdx[1], localIdx[2]);
	//	My_edge me3(fi->first, localIdx[2], localIdx[0]);
	//	Mesh_3_Edge_iterator ei1 = myCellEdgeMap[me1];
	//	Mesh_3_Edge_iterator ei2 = myCellEdgeMap[me2];
	//	Mesh_3_Edge_iterator ei3 = myCellEdgeMap[me3];

	//	if (!edgeOnBoundary[ei1] && !edgeOnBoundary[ei2]) {
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei1]], edgeForward[edgeIdx[ei2]]) = 0;
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei2]], edgeForward[edgeIdx[ei1]]) = 0;
	//	}
	//	if (!edgeOnBoundary[ei2] && !edgeOnBoundary[ei3]) {
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei2]], edgeForward[edgeIdx[ei3]]) = 0;
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei3]], edgeForward[edgeIdx[ei2]]) = 0;
	//	}
	//	if (!edgeOnBoundary[ei3] && !edgeOnBoundary[ei1]) {
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei3]], edgeForward[edgeIdx[ei1]]) = 0;
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei1]], edgeForward[edgeIdx[ei3]]) = 0;
	//	}
	//	if(!edgeOnBoundary[ei1])
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei1]], edgeForward[edgeIdx[ei1]]) += 10;
	//	if (!edgeOnBoundary[ei2])
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei2]], edgeForward[edgeIdx[ei2]]) += 10;
	//	if (!edgeOnBoundary[ei3])
	//		m_CHS1F_NB.coeffRef(edgeForward[edgeIdx[ei2]], edgeForward[edgeIdx[ei3]]) += 10;
	//}
}

void My_Triangulation::buildCHS2F_NB()
{
	std::vector<Eigen::Triplet<double>> triplet;
	for (int i = 0; i < m_CHS2F_B.outerSize(); ++i) {
		for (EigenSpMat::InnerIterator iter(m_CHS2F_B, i); iter; ++iter) {
			int row = iter.row();
			int col = iter.col();
			double val = iter.value();

			Mesh_3_Facet_iterator fir = idxToFacet[row];
			Mesh_3_Facet_iterator fic = idxToFacet[col];

			if (!facetOnBoundary[fir] && !facetOnBoundary[fic]) {
				int nrow = facetForward[row];
				int ncol = facetForward[col];

				triplet.push_back(Eigen::Triplet<double>(nrow, ncol, val));
			}
		}
	}
	m_CHS2F_NB.resize(numIF, numIF);
	m_CHS2F_NB.setFromTriplets(triplet.begin(), triplet.end());

	//m_CHS2F_NB.resize(numIF, numIF);

	//for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
	//	fi != m_c3t3.triangulation().facets_end(); ++fi) {
	//	if (!facetValid[fi])
	//		continue;
	//	if (facetOnBoundary[fi])
	//		continue;

	//	m_CHS2F_NB.coeffRef(facetForward[facetIdx[fi]], facetForward[facetIdx[fi]]) = 0;
	//}

	//for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
	//	ci != m_c3t3.triangulation().cells_end(); ++ci) {
	//	if (ci->subdomain_index() == 0)
	//		continue;

	//	for (int i = 0; i < 4; ++i) {
	//		for (int j = i + 1; j < 4; ++j) {
	//			My_facet mf1(ci, i);
	//			My_facet mf2(ci, j);

	//			Mesh_3_Facet_iterator fi1 = myCellFacetMap[mf1];
	//			Mesh_3_Facet_iterator fi2 = myCellFacetMap[mf2];

	//			if (facetOnBoundary[fi1] || facetOnBoundary[fi2])
	//				continue;

	//			m_CHS2F_NB.coeffRef(facetForward[facetIdx[fi1]], facetForward[facetIdx[fi2]]) = 1;
	//			m_CHS2F_NB.coeffRef(facetForward[facetIdx[fi2]], facetForward[facetIdx[fi1]]) = 1;

	//			m_CHS2F_NB.coeffRef(facetForward[facetIdx[fi1]], facetForward[facetIdx[fi1]]) += 1;
	//			m_CHS2F_NB.coeffRef(facetForward[facetIdx[fi2]], facetForward[facetIdx[fi2]]) += 1;
	//		}
	//	}
	//}
}

void My_Triangulation::buildED0F_B()
{
	// orientation: increasing global indices
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		if (vertexIdx[ei->first->vertex(ei->second)] > vertexIdx[ei->first->vertex(ei->third)]) {
			triplet.push_back(Eigen::Triplet<double>(edgeIdx[ei],
				vertexIdx[ei->first->vertex(ei->second)], 1));
			triplet.push_back(Eigen::Triplet<double>(edgeIdx[ei],
				vertexIdx[ei->first->vertex(ei->third)], -1));
			}
		else {
			triplet.push_back(Eigen::Triplet<double>(edgeIdx[ei],
				vertexIdx[ei->first->vertex(ei->second)], -1));
			triplet.push_back(Eigen::Triplet<double>(edgeIdx[ei],
				vertexIdx[ei->first->vertex(ei->third)], 1));
		}
	}

	m_ED0F_B.resize(numE, numV);
	m_ED0F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildED1F_B()
{
	// orientation: increasing global indices 
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		std::vector<int> vGlobalIdx;
		std::vector<int> vLocalIdx;
		for (int i = 0; i < 4; ++i) {
			if (i == fi->second)
				continue;

			vGlobalIdx.push_back(vertexIdx[fi->first->vertex(i)]);
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
			Mesh_3_Edge_iterator ei = myCellEdgeMap[me];

			triplet.push_back(Eigen::Triplet<double>(facetIdx[fi],
				edgeIdx[ei], val));

		}
	}

	m_ED1F_B.resize(numF, numE); // !!
	m_ED1F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildED2F_B()
{
	// orientation inherited from local index 0123
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

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
			gIdx.push_back(vertexIdx[ci->vertex(lIdx[0])]);
			gIdx.push_back(vertexIdx[ci->vertex(lIdx[1])]);
			gIdx.push_back(vertexIdx[ci->vertex(lIdx[2])]);


			for (int k = 0; k < gIdx.size(); ++k) {
				for (int j = k + 1; j < gIdx.size(); ++j) {
					if (gIdx[k] > gIdx[j]) {
						std::swap(gIdx[k], gIdx[j]);
						even = -even;
					}
				}
			}

			Mesh_3_Facet_iterator fi = myCellFacetMap[mf];

			triplet.push_back(Eigen::Triplet<double>(cellIdx[ci], facetIdx[fi], even));
		}
	}

	m_ED2F_B.resize(numC, numF); // !!
	m_ED2F_B.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildED0F_NB()
{
	// orientation: increasing global indices
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Edge_iterator ei = m_c3t3.triangulation().edges_begin();
		ei != m_c3t3.triangulation().edges_end(); ++ei) {
		if (!edgeValid[ei])
			continue;

		if (!edgeOnBoundary[ei]) {
			if (vertexIdx[ei->first->vertex(ei->second)] > vertexIdx[ei->first->vertex(ei->third)]) {
				if (!vertexOnBoundary[ei->first->vertex(ei->second)])
					triplet.push_back(Eigen::Triplet<double>(edgeForward[edgeIdx[ei]],
						vertexForward[vertexIdx[ei->first->vertex(ei->second)]], 1));
				if (!vertexOnBoundary[ei->first->vertex(ei->third)])
					triplet.push_back(Eigen::Triplet<double>(edgeForward[edgeIdx[ei]],
						vertexForward[vertexIdx[ei->first->vertex(ei->third)]], -1));
			}
			else {
				if (!vertexOnBoundary[ei->first->vertex(ei->second)])
					triplet.push_back(Eigen::Triplet<double>(edgeForward[edgeIdx[ei]],
						vertexForward[vertexIdx[ei->first->vertex(ei->second)]], -1));
				if (!vertexOnBoundary[ei->first->vertex(ei->third)])
					triplet.push_back(Eigen::Triplet<double>(edgeForward[edgeIdx[ei]],
						vertexForward[vertexIdx[ei->first->vertex(ei->third)]], 1));
			}
		}
	}

	m_ED0F_NB.resize(numIE, numIV);
	m_ED0F_NB.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildED1F_NB()
{
	// orientation: increasing global indices 
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Facet_iterator fi = m_c3t3.triangulation().facets_begin();
		fi != m_c3t3.triangulation().facets_end(); ++fi) {
		if (!facetValid[fi])
			continue;

		if (facetOnBoundary[fi]) // !!
			continue;

		std::vector<int> vGlobalIdx;
		std::vector<int> vLocalIdx;
		for (int i = 0; i < 4; ++i) {
			if (i == fi->second)
				continue;

			vGlobalIdx.push_back(vertexIdx[fi->first->vertex(i)]);
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
			Mesh_3_Edge_iterator ei = myCellEdgeMap[me];

			if (!edgeOnBoundary[ei])
				triplet.push_back(Eigen::Triplet<double>(facetForward[facetIdx[fi]], 
					edgeForward[edgeIdx[ei]], val));

		}
	}

	m_ED1F_NB.resize(numIF, numIE); // !!
	m_ED1F_NB.setFromTriplets(triplet.begin(), triplet.end());
}

void My_Triangulation::buildED2F_NB()
{
	// orientation inherited from local index 0123
	std::vector<Eigen::Triplet<double>> triplet;

	for (Mesh_3_Cell_iterator ci = m_c3t3.triangulation().cells_begin();
		ci != m_c3t3.triangulation().cells_end(); ++ci) {
		if (ci->subdomain_index() == 0)
			continue;

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
			gIdx.push_back(vertexIdx[ci->vertex(lIdx[0])]);
			gIdx.push_back(vertexIdx[ci->vertex(lIdx[1])]);
			gIdx.push_back(vertexIdx[ci->vertex(lIdx[2])]);

			//int even = 1;

			for (int i = 0; i < gIdx.size(); ++i) {
				for (int j = i + 1; j < gIdx.size(); ++j) {
					if (gIdx[i] > gIdx[j]) {
						std::swap(gIdx[i], gIdx[j]);
						even = -even;
					}
				}
			}

			Mesh_3_Facet_iterator fi = myCellFacetMap[mf];

			if (!facetOnBoundary[fi])
				triplet.push_back(Eigen::Triplet<double>(cellIdx[ci], facetForward[facetIdx[fi]], even));
		}
	}

	m_ED2F_NB.resize(numC, numIF); // !!
	m_ED2F_NB.setFromTriplets(triplet.begin(), triplet.end());
}

double My_Triangulation::innerProductInCellHS2F(Mesh_3_Cell_iterator ci, Mesh_3_Facet_iterator f1, Mesh_3_Facet_iterator f2,
	std::map<int, Mesh_3_Vector_3> grad)
{
	// get orientation of facets
	std::vector<int> idxsf1;
	for (int i = 0; i < 4; ++i) {
		if (i != f1->second)
			idxsf1.push_back(vertexIdx[f1->first->vertex(i)]);
	}
	std::vector<int> idxsf2;
	for (int i = 0; i < 4; ++i) {
		if (i != f2->second)
			idxsf2.push_back(vertexIdx[f2->first->vertex(i)]);
	}

	std::sort(idxsf1.begin(), idxsf1.end());
	std::sort(idxsf2.begin(), idxsf2.end());


	double val = 0;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			std::vector<int> vec;
			vec.push_back(idxsf1[i]);
			vec.push_back(idxsf1[(i + 1) % 3]);
			vec.push_back(idxsf1[(i + 2) % 3]);
			vec.push_back(idxsf2[j]);
			vec.push_back(idxsf2[(j + 1) % 3]);
			vec.push_back(idxsf2[(j + 2) % 3]);

			val += helperHS2F(grad, vec);
		}
	}

	return val;
}

double My_Triangulation::helperHS2F(std::map<int, Mesh_3_Vector_3> grad, std::vector<int> vec)
{
	double scaling;
	if (vec[0] == vec[3]) {
		scaling = 1.0 / 10;
	}
	else {
		scaling = 1.0 / 20;
	}

	double val = (grad[vec[1]] * grad[vec[4]])*(grad[vec[2]] * grad[vec[5]]) -
		(grad[vec[1]] * grad[vec[5]])*(grad[vec[2]] * grad[vec[4]]);
	val *= 4;

	return scaling * val;
}