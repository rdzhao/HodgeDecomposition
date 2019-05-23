#include "MyOpenGLWidget.h"

void OGLWidget::initializeGL()
{
	mode = 0;

	// general initialize
	leftPressed = false;
	faceSelection = false;
	edgeSelection = false;
	vertSelection = false;

	// opengl initialize
	initializeOpenGLFunctions();

	//connect(this, SIGNAL(frameSwapped()), this, SLOT(update()));

	printContextInformation();

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 1.0f);
	
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void OGLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT);

	for (int i = 0; i < rModules.size(); ++i)
	{
		if (rModules[i]->visible())
		{
			rModules[i]->setCamera(camera);
			
			if(mode == 0)
				rModules[i]->setLightDistance(10 * decomtan.Mesh().getRadius());
			if(mode == 1)
				rModules[i]->setLightDistance(10 * decomnor.Mesh().getRadius());
			if (mode == 2)
				rModules[i]->setLightDistance(10 * lapana.Mesh().getRadius());
			if (mode == 3)
				rModules[i]->setLightDistance(10 * poelke.Mesh().getRadius());
			if (mode == 4)
				rModules[i]->setLightDistance(10 * mxb.Mesh().getRadius());

			rModules[i]->render();
		}
	}
}

void OGLWidget::resizeGL(int width, int height)
{
	camera.setWinWidth(width);
	camera.setWinHeight(height);
	//cout << width << " " << height << endl;
}

void OGLWidget::mousePressEvent(QMouseEvent* event)
{
	if (faceSelection && event->button() == Qt::LeftButton)
	{
		cout << "At Face Selection..."<< endl;

		selectFace(event->x(), event->y());
	}
	else if (edgeSelection && event->button() == Qt::LeftButton)
	{
		cout << "At Edge Selection..." << endl;
	
		selectEdge(event->x(), event->y());
	}
	else if (vertSelection && event->button() == Qt::LeftButton)
	{
		cout << "At Vertex Selection..." << endl;

		selectVertex(event->x(), event->y());
	}
	else if (event->button() == Qt::LeftButton)
	{
		//cout << event->x()<<" "<< event->y() << endl;
		leftPressed = true;

		camera.setPWX(event->x());
		camera.setPWY(event->y());
		camera.updatePUnitCoord();
		camera.setPRotationIdentity();
	}
	else if (event->button() == Qt::RightButton)
	{
		rightPressed = true;

		camera.setPWX(event->x());
		camera.setPWY(event->y());
		camera.updatePUnitCoord();
		camera.setPTranslationIdentity();
	}
}

void OGLWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
	setCamera();
	update();
}

void OGLWidget::mouseReleaseEvent(QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton)
	{
		leftPressed = false;
	}
	else if (event->button() == Qt::RightButton)
	{
		rightPressed = false;
	}
}

void OGLWidget::mouseMoveEvent(QMouseEvent* event)
{
	if (leftPressed)
	{
		//camera.moveUnitCoordToPre();

		camera.setWX(event->x());
		camera.setWY(event->y());
		camera.updateUnitCoord();

		camera.arcballRotate();
		update();
	}
	else if (rightPressed)
	{
		camera.setWX(event->x());
		camera.setWY(event->y());
		camera.updateUnitCoord();

		double ratio = mesh->radius();

		camera.move(ratio);
		update();
	}
}

void OGLWidget::wheelEvent(QWheelEvent* event)
{
	if (event->delta() != 0)
	{
		camera.setscroll(event->delta());
		camera.zoom();
		update();
	}
}

void OGLWidget::clear()
{
	cout<<"@@@@@@@@@@@ "<< rModules.size() <<endl;
	for (int i = 0; i < rModules.size(); ++i) {
		rModules[i]->clearData();
	}
	cout << "@@@@@@@@@@@" << endl;
	rModules.clear();
	cout << "@@@@@@@@@@@" << endl;
	//mesh->clear();
	//delete kdTree;
	cout << "@@@@@@@@@@@" << endl;
	selectedVerts.clear();
}

void OGLWidget::setMesh(Mesh* m)
{
	mesh = m;
}

void OGLWidget::setBall(Mesh* b)
{
	ball = b;
}

void OGLWidget::setKDTree()
{
	auto start = std::chrono::system_clock::now();

	int num_verts, num_tris;
	glm::vec3 *verts, *tris;
	// import mesh data
	num_verts = static_cast<int>(mesh->size_of_vertices());
	verts = new glm::vec3[num_verts];
	for (Vertex_iterator vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		verts[vi->idx()] = 
		glm::vec3(vi->point().x() - mesh->xcenter(),
			vi->point().y() - mesh->ycenter(),
			vi->point().z() - mesh->zcenter());
	num_tris = static_cast<int>(mesh->size_of_facets());
	tris = new glm::vec3[num_tris];
	for (Facet_iterator fi = mesh->facets_begin(); fi != mesh->facets_end(); ++fi)
		tris[fi->idx()] =
		glm::vec3(fi->halfedge()->vertex()->idx(),
			fi->halfedge()->next()->vertex()->idx(),
			fi->halfedge()->next()->next()->vertex()->idx());

	kdTree = new KDTreeCPU(num_tris, tris, num_verts, verts);
	
	cout << "Num Verts: " << num_verts << endl;
	cout << "Num Faces: " << num_tris << endl;
	cout << "KD tree construction complete ..." << endl;

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "Construction Time: " << elapsed_seconds.count() << endl;

}

void OGLWidget::setRenderContexts()
{
	meshModule = new MeshModule();
	rModules.push_back(meshModule);
	
	wireFrameModule = new WireFrameModule();
	rModules.push_back(wireFrameModule);

	//vertexHLModule = new MeshModule();
	//rModules.push_back(vertexHLModule);

	vectorFieldModule = new VectorFieldModule();
	rModules.push_back(vectorFieldModule);

	setMeshModule(meshModule);
	setWireFrameModule(wireFrameModule);
	setVectorFieldModule(vectorFieldModule);

}

void OGLWidget::setMeshModule(RenderModule* rm)
{
	std::vector<float> v;
	std::vector<int> idx;
	std::vector<float> n;
	std::vector<float> c;

	//int mode = 1; // 0: tan 1: nor 3:poelke

	if(mode ==0)
	{ // decomtan

		for (int i = 0; i < decomtan.Indices().size() / 3; ++i) {
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].z() - decomtan.Mesh().getCenter().z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i]].x());
			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i]].y());
			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i]].z());

			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i + 1]].x());
			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i + 1]].y());
			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i + 1]].z());

			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i + 2]].x());
			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i + 2]].y());
			c.push_back(decomtan.Colors()[decomtan.Indices()[3 * i + 2]].z());
		}
	}

	if(mode ==1)
	{ // decomnor

		for (int i = 0; i < decomnor.Indices().size() / 3; ++i) {
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].z() - decomnor.Mesh().getCenter().z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i]].x());
			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i]].y());
			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i]].z());

			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i + 1]].x());
			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i + 1]].y());
			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i + 1]].z());

			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i + 2]].x());
			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i + 2]].y());
			c.push_back(decomnor.Colors()[decomnor.Indices()[3 * i + 2]].z());
		}
	}

	if (mode == 2)
	{ // analyzer

		for (int i = 0; i < lapana.Indices().size() / 3; ++i) {
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].z() - lapana.Mesh().getCenter().z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			c.push_back(lapana.Colors()[lapana.Indices()[3 * i]].x());
			c.push_back(lapana.Colors()[lapana.Indices()[3 * i]].y());
			c.push_back(lapana.Colors()[lapana.Indices()[3 * i]].z());

			c.push_back(lapana.Colors()[lapana.Indices()[3 * i + 1]].x());
			c.push_back(lapana.Colors()[lapana.Indices()[3 * i + 1]].y());
			c.push_back(lapana.Colors()[lapana.Indices()[3 * i + 1]].z());

			c.push_back(lapana.Colors()[lapana.Indices()[3 * i + 2]].x());
			c.push_back(lapana.Colors()[lapana.Indices()[3 * i + 2]].y());
			c.push_back(lapana.Colors()[lapana.Indices()[3 * i + 2]].z());
		}
	}

	if (mode == 3)
	{ // poelke

		for (int i = 0; i < poelke.Indices().size() / 3; ++i) {
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].z() - poelke.Mesh().getCenter().z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			c.push_back(poelke.Colors()[poelke.Indices()[3 * i]].x());
			c.push_back(poelke.Colors()[poelke.Indices()[3 * i]].y());
			c.push_back(poelke.Colors()[poelke.Indices()[3 * i]].z());

			c.push_back(poelke.Colors()[poelke.Indices()[3 * i + 1]].x());
			c.push_back(poelke.Colors()[poelke.Indices()[3 * i + 1]].y());
			c.push_back(poelke.Colors()[poelke.Indices()[3 * i + 1]].z());

			c.push_back(poelke.Colors()[poelke.Indices()[3 * i + 2]].x());
			c.push_back(poelke.Colors()[poelke.Indices()[3 * i + 2]].y());
			c.push_back(poelke.Colors()[poelke.Indices()[3 * i + 2]].z());
		}
	}

	if (mode == 4)
	{ // mixedBoundary

		for (int i = 0; i < mxb.Indices().size() / 3; ++i) {
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].z() - mxb.Mesh().getCenter().z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			c.push_back(mxb.Colors()[mxb.Indices()[3 * i]].x());
			c.push_back(mxb.Colors()[mxb.Indices()[3 * i]].y());
			c.push_back(mxb.Colors()[mxb.Indices()[3 * i]].z());

			c.push_back(mxb.Colors()[mxb.Indices()[3 * i + 1]].x());
			c.push_back(mxb.Colors()[mxb.Indices()[3 * i + 1]].y());
			c.push_back(mxb.Colors()[mxb.Indices()[3 * i + 1]].z());

			c.push_back(mxb.Colors()[mxb.Indices()[3 * i + 2]].x());
			c.push_back(mxb.Colors()[mxb.Indices()[3 * i + 2]].y());
			c.push_back(mxb.Colors()[mxb.Indices()[3 * i + 2]].z());
		}
	}

	rm->setData(v, idx, n, c);
}

void OGLWidget::setWireFrameModule(RenderModule* rm)
{
	std::vector<float> v;
	std::vector<int> idx;
	std::vector<float> n;
	std::vector<float> c;


	//int mode = 0; // 0: tan 1: nor

	if (mode == 0)
	{ // decomtan

		for (int i = 0; i < decomtan.Indices().size() / 3; ++i) {
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 1]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i + 2]].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.Vertices()[decomtan.Indices()[3 * i]].z() - decomtan.Mesh().getCenter().z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			n.push_back(decomtan.Normals()[i].x());
			n.push_back(decomtan.Normals()[i].y());
			n.push_back(decomtan.Normals()[i].z());

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
		}
	}

	if(mode == 1)
	{ // decomnor

		for (int i = 0; i < decomnor.Indices().size() / 3; ++i) {
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 1]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i + 2]].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.Vertices()[decomnor.Indices()[3 * i]].z() - decomnor.Mesh().getCenter().z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			n.push_back(decomnor.Normals()[i].x());
			n.push_back(decomnor.Normals()[i].y());
			n.push_back(decomnor.Normals()[i].z());

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
		}
	}

	if (mode == 2)
	{ // lapana

		for (int i = 0; i < lapana.Indices().size() / 3; ++i) {
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 1]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i + 2]].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.Vertices()[lapana.Indices()[3 * i]].z() - lapana.Mesh().getCenter().z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			n.push_back(lapana.Normals()[i].x());
			n.push_back(lapana.Normals()[i].y());
			n.push_back(lapana.Normals()[i].z());

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
		}
	}

	if (mode == 3)
	{ // poelke

		for (int i = 0; i < poelke.Indices().size() / 3; ++i) {
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 1]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i + 2]].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.Vertices()[poelke.Indices()[3 * i]].z() - poelke.Mesh().getCenter().z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			n.push_back(poelke.Normals()[i].x());
			n.push_back(poelke.Normals()[i].y());
			n.push_back(poelke.Normals()[i].z());

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
		}
	}

	if (mode == 4)
	{ // mixedBoundary

		for (int i = 0; i < mxb.Indices().size() / 3; ++i) {
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 1]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i + 2]].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.Vertices()[mxb.Indices()[3 * i]].z() - mxb.Mesh().getCenter().z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			n.push_back(mxb.Normals()[i].x());
			n.push_back(mxb.Normals()[i].y());
			n.push_back(mxb.Normals()[i].z());

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);

			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
			c.push_back(0); c.push_back(0); c.push_back(0);
		}
	}

	rm->setData(v, idx, n, c);
}

void OGLWidget::setVectorFieldModule(RenderModule* rm)
{
	std::vector<float> v;
	std::vector<int> idx;
	std::vector<float> n;
	std::vector<float> c;


	// insert arrow mesh for each vector

	//int mode = 0; // 0: tan 1: nor

	if (mode == 0)
	{ // decomtan

		for (int i = 0; i < decomtan.VfFaces().size(); ++i) {
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].x()].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].x()].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].x()].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].y()].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].y()].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].y()].z() - decomtan.Mesh().getCenter().z());

			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].z()].x() - decomtan.Mesh().getCenter().x());
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].z()].y() - decomtan.Mesh().getCenter().y());
			v.push_back(decomtan.VfVertices()[decomtan.VfFaces()[i].z()].z() - decomtan.Mesh().getCenter().z());

			n.push_back(decomtan.VfNormals()[i].x());
			n.push_back(decomtan.VfNormals()[i].y());
			n.push_back(decomtan.VfNormals()[i].z());

			n.push_back(decomtan.VfNormals()[i].x());
			n.push_back(decomtan.VfNormals()[i].y());
			n.push_back(decomtan.VfNormals()[i].z());

			n.push_back(decomtan.VfNormals()[i].x());
			n.push_back(decomtan.VfNormals()[i].y());
			n.push_back(decomtan.VfNormals()[i].z());

			c.push_back(decomtan.VfColors()[i].x());
			c.push_back(decomtan.VfColors()[i].y());
			c.push_back(decomtan.VfColors()[i].z());

			c.push_back(decomtan.VfColors()[i].x());
			c.push_back(decomtan.VfColors()[i].y());
			c.push_back(decomtan.VfColors()[i].z());

			c.push_back(decomtan.VfColors()[i].x());
			c.push_back(decomtan.VfColors()[i].y());
			c.push_back(decomtan.VfColors()[i].z());
		}
	}

	if(mode == 1)
	{ // decomnor
		for (int i = 0; i < decomnor.VfFaces().size(); ++i) {
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].x()].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].x()].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].x()].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].y()].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].y()].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].y()].z() - decomnor.Mesh().getCenter().z());

			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].z()].x() - decomnor.Mesh().getCenter().x());
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].z()].y() - decomnor.Mesh().getCenter().y());
			v.push_back(decomnor.VfVertices()[decomnor.VfFaces()[i].z()].z() - decomnor.Mesh().getCenter().z());

			n.push_back(decomnor.VfNormals()[i].x());
			n.push_back(decomnor.VfNormals()[i].y());
			n.push_back(decomnor.VfNormals()[i].z());

			n.push_back(decomnor.VfNormals()[i].x());
			n.push_back(decomnor.VfNormals()[i].y());
			n.push_back(decomnor.VfNormals()[i].z());

			n.push_back(decomnor.VfNormals()[i].x());
			n.push_back(decomnor.VfNormals()[i].y());
			n.push_back(decomnor.VfNormals()[i].z());

			c.push_back(decomnor.VfColors()[i].x());
			c.push_back(decomnor.VfColors()[i].y());
			c.push_back(decomnor.VfColors()[i].z());

			c.push_back(decomnor.VfColors()[i].x());
			c.push_back(decomnor.VfColors()[i].y());
			c.push_back(decomnor.VfColors()[i].z());

			c.push_back(decomnor.VfColors()[i].x());
			c.push_back(decomnor.VfColors()[i].y());
			c.push_back(decomnor.VfColors()[i].z());
		}
	}

	if (mode == 2)
	{ // lapana
		for (int i = 0; i < lapana.VfFaces().size(); ++i) {
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].x()].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].x()].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].x()].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].y()].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].y()].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].y()].z() - lapana.Mesh().getCenter().z());

			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].z()].x() - lapana.Mesh().getCenter().x());
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].z()].y() - lapana.Mesh().getCenter().y());
			v.push_back(lapana.VfVertices()[lapana.VfFaces()[i].z()].z() - lapana.Mesh().getCenter().z());

			n.push_back(lapana.VfNormals()[i].x());
			n.push_back(lapana.VfNormals()[i].y());
			n.push_back(lapana.VfNormals()[i].z());

			n.push_back(lapana.VfNormals()[i].x());
			n.push_back(lapana.VfNormals()[i].y());
			n.push_back(lapana.VfNormals()[i].z());

			n.push_back(lapana.VfNormals()[i].x());
			n.push_back(lapana.VfNormals()[i].y());
			n.push_back(lapana.VfNormals()[i].z());

			c.push_back(lapana.VfColors()[i].x());
			c.push_back(lapana.VfColors()[i].y());
			c.push_back(lapana.VfColors()[i].z());

			c.push_back(lapana.VfColors()[i].x());
			c.push_back(lapana.VfColors()[i].y());
			c.push_back(lapana.VfColors()[i].z());

			c.push_back(lapana.VfColors()[i].x());
			c.push_back(lapana.VfColors()[i].y());
			c.push_back(lapana.VfColors()[i].z());
		}
	}

	if (mode == 3)
	{ // poelke
		for (int i = 0; i < poelke.VfFaces().size(); ++i) {
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].x()].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].x()].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].x()].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].y()].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].y()].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].y()].z() - poelke.Mesh().getCenter().z());

			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].z()].x() - poelke.Mesh().getCenter().x());
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].z()].y() - poelke.Mesh().getCenter().y());
			v.push_back(poelke.VfVertices()[poelke.VfFaces()[i].z()].z() - poelke.Mesh().getCenter().z());

			n.push_back(poelke.VfNormals()[i].x());
			n.push_back(poelke.VfNormals()[i].y());
			n.push_back(poelke.VfNormals()[i].z());

			n.push_back(poelke.VfNormals()[i].x());
			n.push_back(poelke.VfNormals()[i].y());
			n.push_back(poelke.VfNormals()[i].z());

			n.push_back(poelke.VfNormals()[i].x());
			n.push_back(poelke.VfNormals()[i].y());
			n.push_back(poelke.VfNormals()[i].z());

			c.push_back(poelke.VfColors()[i].x());
			c.push_back(poelke.VfColors()[i].y());
			c.push_back(poelke.VfColors()[i].z());

			c.push_back(poelke.VfColors()[i].x());
			c.push_back(poelke.VfColors()[i].y());
			c.push_back(poelke.VfColors()[i].z());

			c.push_back(poelke.VfColors()[i].x());
			c.push_back(poelke.VfColors()[i].y());
			c.push_back(poelke.VfColors()[i].z());
		}
	}

	if (mode == 4)
	{ // mixedBoundary
		for (int i = 0; i < mxb.VfFaces().size(); ++i) {
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].x()].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].x()].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].x()].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].y()].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].y()].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].y()].z() - mxb.Mesh().getCenter().z());

			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].z()].x() - mxb.Mesh().getCenter().x());
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].z()].y() - mxb.Mesh().getCenter().y());
			v.push_back(mxb.VfVertices()[mxb.VfFaces()[i].z()].z() - mxb.Mesh().getCenter().z());

			n.push_back(mxb.VfNormals()[i].x());
			n.push_back(mxb.VfNormals()[i].y());
			n.push_back(mxb.VfNormals()[i].z());

			n.push_back(mxb.VfNormals()[i].x());
			n.push_back(mxb.VfNormals()[i].y());
			n.push_back(mxb.VfNormals()[i].z());

			n.push_back(mxb.VfNormals()[i].x());
			n.push_back(mxb.VfNormals()[i].y());
			n.push_back(mxb.VfNormals()[i].z());

			c.push_back(mxb.VfColors()[i].x());
			c.push_back(mxb.VfColors()[i].y());
			c.push_back(mxb.VfColors()[i].z());

			c.push_back(mxb.VfColors()[i].x());
			c.push_back(mxb.VfColors()[i].y());
			c.push_back(mxb.VfColors()[i].z());

			c.push_back(mxb.VfColors()[i].x());
			c.push_back(mxb.VfColors()[i].y());
			c.push_back(mxb.VfColors()[i].z());
		}
	}

	rm->setData(v, idx, n, c);
}

void OGLWidget::addVertexHL(RenderModule* rm, float x, float y, float z)
{
	std::vector<float> v;
	std::vector<int> idx;
	std::vector<float> n;
	std::vector<float> c;

	double ave=0;
	for (Halfedge_iterator hei = mesh->halfedges_begin(); hei != mesh->halfedges_end(); ++hei)
		ave += hei->length();
	ave /= 2 * mesh->size_of_halfedges();
	float ratio = ave / ball->radius()*0.2;
	
	for (Facet_iterator fi = ball->facets_begin(); fi != ball->facets_end(); ++fi)
	{
		Halfedge_around_facet_circulator he = fi->facet_begin();
		Halfedge_around_facet_circulator end = he;

		do {
			v.push_back((he->vertex()->point().x() - ball->xcenter())*ratio + x - mesh->xcenter());
			v.push_back((he->vertex()->point().y() - ball->ycenter())*ratio + y - mesh->ycenter());
			v.push_back((he->vertex()->point().z() - ball->zcenter())*ratio + z - mesh->zcenter());

			n.push_back(he->vertex()->normal().x());
			n.push_back(he->vertex()->normal().y());
			n.push_back(he->vertex()->normal().z());

			c.push_back(0.8f);
			c.push_back(0.6f);
			c.push_back(0.2f);

			//idx.push_back(he->vertex()->idx());
			++he;
		} while (he != end);
	}

	rm->appData(v, idx, n, c);
}

void OGLWidget::delVertexHL(RenderModule* rm, int pos)
{
	int posBegin = pos*static_cast<int>(ball->size_of_facets()) * 9;
	rm->delData(posBegin, static_cast<int>(ball->size_of_facets()) * 9);
}

void OGLWidget::setCamera()
{
	if(mode == 0)
		camera.setView(QVector3D(0, 0, 6* decomtan.Mesh().getRadius()), QVector3D(0, 0, 0), QVector3D(0, 1, 0));
	if (mode == 1)
		camera.setView(QVector3D(0, 0, 6 * decomnor.Mesh().getRadius()), QVector3D(0, 0, 0), QVector3D(0, 1, 0));
	if (mode == 2)
		camera.setView(QVector3D(0, 0, 6 * lapana.Mesh().getRadius()), QVector3D(0, 0, 0), QVector3D(0, 1, 0));
	if (mode == 3)
		camera.setView(QVector3D(0, 0, 6 * poelke.Mesh().getRadius()), QVector3D(0, 0, 0), QVector3D(0, 1, 0));
	if (mode == 4)
		camera.setView(QVector3D(0, 0, 6 * mxb.Mesh().getRadius()), QVector3D(0, 0, 0), QVector3D(0, 1, 0));


	camera.setProject(45.0, 4.0 / 3.0, 0.01, 10000);
	
	camera.init();
}

void OGLWidget::printContextInformation()
{
	QString glType;
	QString glVersion;
	QString glProfile;

	// Get Version Information
	glType = (context()->isOpenGLES()) ? "OpenGL ES" : "OpenGL";
	glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

	// Get Profile Information
#define CASE(c) case QSurfaceFormat::c: glProfile = #c; break
	switch (format().profile())
	{
		CASE(NoProfile);
		CASE(CoreProfile);
		CASE(CompatibilityProfile);
	}
#undef CASE

	// qPrintable() will print our QString w/o quotes around it.
	qDebug() << qPrintable(glType) << qPrintable(glVersion) << "(" << qPrintable(glProfile) << ")";
}

void OGLWidget::setFaceVisible(bool v)
{
	meshModule->setVisible(v);
	update();
}

void OGLWidget::setEdgeVisible(bool v)
{
	wireFrameModule->setVisible(v);
	update();
}

void OGLWidget::setVectorFieldVisible(bool v)
{
	vectorFieldModule->setVisible(v);
	update();
}

void OGLWidget::setFaceSelection(bool b)
{
	faceSelection = b;
}

void OGLWidget::setEdgeSelection(bool b)
{
	edgeSelection = b;
}

void OGLWidget::setVertSelection(bool b)
{
	vertSelection = b;
}

void OGLWidget::selectFace(int wx, int wy)
{
	QVector3D nearP, farP, d;
	camera.getFarNearPointWorld(wx, wy, nearP, farP);
	d = (farP - nearP).normalized();

	glm::vec3 rayO, rayD, hitP, normal;
	float t;
	int idx;
	rayO = glm::vec3(nearP.x(), nearP.y(), nearP.z());
	rayD = glm::vec3(d.x(), d.y(), d.z());
	bool intersected = kdTree->intersectNew(rayO, rayD, t, hitP, normal, idx);

	if (intersected)
	{
		mesh->facet(idx)->selected() = !mesh->facet(idx)->selected();
		meshModule->highlightFace(idx, mesh->facet(idx)->selected());
		update();
	}
}

void OGLWidget::selectEdge(int wx, int wy)
{
	QVector3D nearP, farP, d;
	camera.getFarNearPointWorld(wx, wy, nearP, farP);
	d = (farP - nearP).normalized();

	glm::vec3 rayO, rayD, hitP, normal;
	float t;
	int idx;
	rayO = glm::vec3(nearP.x(), nearP.y(), nearP.z());
	rayD = glm::vec3(d.x(), d.y(), d.z());
	bool intersected = kdTree->intersectNew(rayO, rayD, t, hitP, normal, idx);

	if (intersected)
	{
		// choose the closest edge to the hit point
		double distance = std::numeric_limits<double>::max();
		int eIdx;

		Point_3 hp(hitP.x + mesh->xcenter(), hitP.y + mesh->ycenter(), hitP.z + mesh->zcenter());
		Facet_iterator fi = mesh->facet(idx);
		Halfedge_iterator hei = fi->halfedge();
		do {
			Vector_3 ev, hv;
			ev = hei->opposite()->vertex()->point() - hei->vertex()->point();
			hv = hp - hei->vertex()->point();
			glm::vec3 gev, ghv;
			gev = glm::vec3(ev.x(), ev.y(), ev.z());
			ghv = glm::vec3(hv.x(), hv.y(), hv.z());
			double dis = glm::length((ghv - (glm::dot(ghv, glm::normalize(gev))*glm::normalize(gev))));

			if (dis < distance)
			{
				distance = dis;
				eIdx = hei->idx();
			}

			hei = hei->next();
		} while (hei != fi->halfedge());

		mesh->halfedge(eIdx)->selected() = !mesh->halfedge(eIdx)->selected();
		mesh->halfedge(eIdx)->opposite()->selected() = mesh->halfedge(eIdx)->selected();
		wireFrameModule->highlightEdge(eIdx, mesh->halfedge(eIdx)->selected());
		update();
	}
}

void OGLWidget::selectVertex(int wx, int wy)
{
	QVector3D nearP, farP, d;
	camera.getFarNearPointWorld(wx, wy, nearP, farP);
	d = (farP - nearP).normalized();

	glm::vec3 rayO, rayD, hitP, normal;
	float t;
	int idx;
	rayO = glm::vec3(nearP.x(), nearP.y(), nearP.z());
	rayD = glm::vec3(d.x(), d.y(), d.z());
	bool intersected = kdTree->intersectNew(rayO, rayD, t, hitP, normal, idx);

	if (intersected){
		double distance = std::numeric_limits<double>::max();
		int vIdx;

		glm::vec3 hp(hitP.x + mesh->xcenter(), hitP.y + mesh->ycenter(), hitP.z + mesh->zcenter());
		Facet_iterator fi = mesh->facet(idx);
		Halfedge_iterator hei = fi->halfedge();
		do {
			glm::vec3 gp(hei->vertex()->point().x(), hei->vertex()->point().y(), hei->vertex()->point().z());
			double dis = glm::length(hp - gp);

			if (dis < distance){
				distance = dis;
				vIdx = hei->vertex()->idx();
			}

			hei = hei->next();
		} while (hei != fi->halfedge());
		//cout << "Vertex Idx: " << vIdx << endl;

		if (!mesh->vertex(vIdx)->selected()) {
			selectedVerts.push_back(vIdx);
			mesh->vertex(vIdx)->selected() = !mesh->vertex(vIdx)->selected();
			mesh->vertex(vIdx)->pos() = static_cast<int>(selectedVerts.size()) - 1;
			addVertexHL(vertexHLModule, mesh->vertex(vIdx)->point().x(), mesh->vertex(vIdx)->point().y(), mesh->vertex(vIdx)->point().z());
		}
		else{
			//update selected pos
			selectedVerts.erase(selectedVerts.begin() + mesh->vertex(vIdx)->pos());
			//update pos
			for (int i = 0; i < selectedVerts.size(); ++i)
				mesh->vertex(selectedVerts[i])->pos() = i;
			mesh->vertex(vIdx)->selected() = !mesh->vertex(vIdx)->selected();
			delVertexHL(vertexHLModule, mesh->vertex(vIdx)->pos());
		}

		update();
	}
}

void OGLWidget::computeTriangulation_0()
{
	//mesh_3.processRawData();
	//mesh_3.counter = 0;
	//mesh_3.setCriteria(30, 2, 2, 2, 2);
	//mesh_3.buildTetMesh();
	//mesh_3.preprocessing();
	//mesh_3.buildDECOperators_0();

	//mesh_3.decompose();  // orthogonal form
}

void OGLWidget::computeTriangulation_2()
{
	//mesh_3.counter = 0;
	//mesh_3.setCriteria(30, 0.2, 0.2, 2, 0.2);
	//mesh_3.buildTetMesh();
	//mesh_3.preprocessing();
	//mesh_3.buildDECOperators_2();

	//mesh_3.decomposeTangential(); // tangential form

}

void OGLWidget::computeLaplacian0()
{
	//lapana_0.buildMeshFromSurface("proteins\\1ajj.off", 2);
	//lapana_0.buildDECOperators();
	//cout<<"Laplacian ..."<<endl;
	//lapana_0.buildLaplacian0();
	//cout<<"Laplacian done ..."<<endl;
	////lapana_0.analyze();
	//lapana_0.analyzeSpectra();
	//lapana_0.visualize();
}

void OGLWidget::computeLaplacian1()
{
	//lapana.buildMeshFromSurface("models\\bimba.off", 0.06);
	////lapana_1.buildMeshFromSurface("models\\moai.off", 0.06);

	//lapana_1.buildDECOperators();
	//lapana_1.buildLaplacian1();
	//lapana_1.analyze();
	//lapana_1.group();
	//lapana_1.visualize(5, false);
	//lapana_1.integrate();
	//lapana_1.write();

	lapana.init(0.004, 0.2, 40, "sphere");
	lapana.buildMeshFromSurface("models\\sphere.off", 0.2);
	//lapana.buildMeshFromSurface("emd\\emd_1776.off", 0.04);
	
	cout << "Triangulation Done ..." << endl;
	//lapana.buildLaplacian();
	cout << "Laplacian Done ..." << endl;
	//lapana.decompose();
	//lapana.setEigenfield(0, 0, true);
	//lapana.setEigenfield(1, 0, false);
	//lapana.setEigenfield(2, 0, false);
	//lapana.setEigenfield(3, 2, true);
	//lapana.setEigenfield(4, 10, false);
	//lapana.setEigenfield(5, 2, false);
	//lapana.setEigenfield(6, 0, false);
	//lapana.setEigenfield(7, 0, false);
	//lapana.process_mxb();
	cout << "Decomposition Done ..." << endl;
	lapana.visualize();
	//lapana.integrate();
	cout << "Post Processing Done ..." << endl;
}

void OGLWidget::computeLaplacian2()
{
	////lapana_2.buildMeshFromSurface("models\\kitten100K.off", 0.04);
	//lapana_2.buildMeshFromSurface("models\\cheese.off", 0.02);
	//
	//lapana_2.buildDECOperators();
	//lapana_2.buildLaplacian2();
	//lapana_2.analyze();
	//lapana_2.group();
	////lapana_2.batchIntegrate(); // problem, need debug ...
	//lapana_2.visualize(0, false);
	//lapana_2.integrate();
	//lapana_2.write();


}

void OGLWidget::computeDecompositionTangential()
{
	decomtan.init(0.004, 0.2, 50, "emd");
	decomtan.buildMeshFromSurface("emd\\emd_1590.off", 0.04);
	cout<<"Triangulation Done ..."<<endl;
	decomtan.buildLaplacian();
	cout<<"Laplacian Done ..."<<endl;
	decomtan.decompose();
	cout<<"Decomposition Done ..."<<endl;
	decomtan.visualize();
	decomtan.integrate();
}

void OGLWidget::computeDecompositionNormal()
{
	decomnor.init(0.004, 0.2, 500, "bunny");
	decomnor.buildMeshFromSurface("models\\torus.off", 0.03);
	cout << "Triangulation Done ..." << endl;
	decomnor.buildLaplacian();
	cout << "Laplacian Done ..." << endl;
	decomnor.decompose();
	cout << "Decomposition Done ..." << endl;
	decomnor.visualize();
	decomnor.integrate();
	cout<<"Post Processing Done ..."<<endl;
}

void OGLWidget::computePoelke()
{
	poelke.init(0.01, 0.2, 500, "sphere");
	poelke.buildMeshFromSurface("compare\\2sphere_plot.off", 0.27);
	cout << "Triangulation Done ..." << endl;
	poelke.prepare();
	cout << "Preparation Done ..." << endl;
	poelke.decompose();
	cout << "Decomposition Done ..." << endl;
	poelke.visualize();
}

void OGLWidget::computeMixedBoundary()
{
	mxb.init(0.01, 0.2, 500, "joint");
	mxb.buildMeshFromSurface("models\\joint.off", 0.05);
	//mxb.readSimulationData("simulation\\channel_0.tet", "simulation\\channel_0.field");
	cout << "Triangulation Done ..." << endl;
	mxb.prepare();
	cout << "Preparation Done ..." << endl;
	mxb.decompose();
	cout << "Decomposition Done ..." << endl;
	mxb.visualize();
	cout << "Visualization Done ..." << endl;
	mxb.integrate();
	cout << "Integration Done ..." << endl;
}

void OGLWidget::analyzeFlexibility()
{
	//flex.init("proteins\\1AKG_CA_A2.pdb", "proteins\\mesh_1AKG_CA_A2.off", 100, 2.0, 0);
	//flex.compute();
	//flex.write();
}

void OGLWidget::visualizeTriangulation()
{	
	//mesh_3.outputBoundaryMesh();
	//mesh_3.normalizeInputMesh();
	//mesh_3.computeCrossSection();
	

	//std::cout << "Compute RK streamlines ..." << std::endl;
	//mesh_3.computeRKStreamlines();
	//std::cout << "Done ..." << std::endl;

	//std::cout << "Compute VMD ...... " << std::endl;
	//mesh_3.writeVMDNMD();
	//std::cout << "Done ......" << std::endl;

	//std::cout << "Compute Streamlines ...... " << std::endl;
	//mesh_3.computeStreamlines();
	//std::cout << "Done ......" << std::endl;


	
	//mesh_3.testInfinite();
	//mesh_3.testBoundaryPieces();
	//mesh_3.testVertex();
	//mesh_3.testOrthogonality();
}