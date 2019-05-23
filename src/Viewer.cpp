#include "Viewer.h"

Viewer::Viewer(QWidget *parent)
	: QMainWindow(parent)
{
	createCanvas();
	createMenu();
	createStatus();
	createToolBar();
}

void Viewer::createMenu()
{
	menuBar()->setMinimumHeight(25);

	fileMenu = menuBar()->addMenu("File");
	fileMenu->setFixedWidth(250);
	{
		newAct = new QAction("New", this);
		newAct->setShortcuts(QKeySequence::New);
		fileMenu->addAction(newAct);

		openAct = new QAction("Open", this);
		openAct->setShortcuts(QKeySequence::Open);
		fileMenu->addAction(openAct);
		connect(openAct, SIGNAL(triggered()), this, SLOT(load()));

		saveAct = new QAction("Save", this);
		saveAct->setShortcuts(QKeySequence::Save);
		fileMenu->addAction(saveAct);

		exitAct = new QAction("Exit", this);
		exitAct->setShortcuts(QKeySequence::Quit);
		fileMenu->addAction(exitAct);
		connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
	}

	sceneMenu = menuBar()->addMenu("Scene");
	sceneMenu->setFixedWidth(250);
	{
		showOriginalFieldAct = new QAction("Show Original Field", this);
		showOriginalFieldAct->setShortcut(QKeySequence(tr("Ctrl+F")));
		sceneMenu->addAction(showOriginalFieldAct);
		connect(showOriginalFieldAct, SIGNAL(triggered()), this, SLOT(showVectorFieldSlot()));

		ClearAct = new QAction("Clear Scene", this);
		ClearAct->setShortcut(QKeySequence(tr("Ctrl+P")));
		sceneMenu->addAction(ClearAct);
		connect(ClearAct, SIGNAL(triggered()), this, SLOT(clearSceneSlot()));
	}

	opMenu = menuBar()->addMenu("Operation");
	opMenu->setFixedWidth(250);
	{
		//decomposeAct = new QAction("Vector Field Decomposition", this);
		//decomposeAct->setShortcut(tr("Crtl+D"));
		//opMenu->addAction(decomposeAct);
		//connect(decomposeAct, SIGNAL(triggered()), this, SLOT(decomposeSlot()));

		analyzeLaplacian0Act = new QAction("Analyze Laplacian 0", this);
		analyzeLaplacian0Act->setShortcut(tr("Crtl+Q"));
		opMenu->addAction(analyzeLaplacian0Act);
		connect(analyzeLaplacian0Act, SIGNAL(triggered()), this, SLOT(analyzeLaplacian0Slot()));

		analyzeLaplacian1Act = new QAction("Analyze Laplacian 1", this);
		analyzeLaplacian1Act->setShortcut(tr("Crtl+W"));
		opMenu->addAction(analyzeLaplacian1Act);
		connect(analyzeLaplacian1Act, SIGNAL(triggered()), this, SLOT(analyzeLaplacian1Slot()));

		analyzeLaplacian2Act = new QAction("Analyze Laplacian 2", this);
		analyzeLaplacian2Act->setShortcut(tr("Crtl+E"));
		opMenu->addAction(analyzeLaplacian2Act);
		connect(analyzeLaplacian2Act, SIGNAL(triggered()), this, SLOT(analyzeLaplacian2Slot()));

		decomposeNormalAct = new QAction("Decompose Normal", this);
		decomposeNormalAct->setShortcut(tr("Crtl+N"));
		opMenu->addAction(decomposeNormalAct);
		connect(decomposeNormalAct, SIGNAL(triggered()), this, SLOT(decomposeNormalSlot()));

		decomposeTangentialAct = new QAction("Decompose Tangential", this);
		decomposeTangentialAct->setShortcut(tr("Crtl+T"));
		opMenu->addAction(decomposeTangentialAct);
		connect(decomposeTangentialAct, SIGNAL(triggered()), this, SLOT(decomposeTangentialSlot()));

		flexibilityAct = new QAction("Analyze Flexibility", this);
		flexibilityAct->setShortcut(tr("Crtl+F"));
		opMenu->addAction(flexibilityAct);
		connect(flexibilityAct, SIGNAL(triggered()), this, SLOT(flexibilitySlot()));

		poelkeAct = new QAction("Poelke!", this);
		poelkeAct->setShortcut(tr("Crtl+P"));
		opMenu->addAction(poelkeAct);
		connect(poelkeAct, SIGNAL(triggered()), this, SLOT(poelkeSlot()));

		mixedBoundaryAct = new QAction("Mixed Boundary", this);
		mixedBoundaryAct->setShortcut(tr("Crtl+M"));
		opMenu->addAction(mixedBoundaryAct);
		connect(mixedBoundaryAct, SIGNAL(triggered()), this, SLOT(mixedBoundarySlot()));


	}
}

void Viewer::createCanvas()
{
	canvas = new OGLWidget();
	setCentralWidget(canvas);

	QSurfaceFormat format;
	format.setRenderableType(QSurfaceFormat::OpenGL);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setVersion(3, 3);
	format.setSamples(4);

	canvas->setFormat(format);
	canvas->show();

	string fn = "../../resources/ball.obj";
	ObjReader<Kernel, Enriched_items> reader(fn);
	reader.read();

	ball = new Enriched_polyhedron<Kernel, Enriched_items>;
	ObjBuilder<HalfedgeDS> builder(reader.vertices(), reader.facets());
	ball->delegate(builder);
	ball->basic_init();
	canvas->setBall(ball);
}

void Viewer::createStatus()
{
	statusBar()->show();
}

void Viewer::createToolBar()
{
	toolBar = addToolBar("Tools");
	toolBar->setMinimumHeight(50);
	toolBar->setMinimumWidth(50);
	toolBar->setIconSize(QSize(45, 45));

	tbImportAct = new QAction(this);
	tbImportAct->setIcon(QIcon("../../resources/import.png"));
	tbImportAct->setText("Import Mesh");
	toolBar->addAction(tbImportAct);
	connect(tbImportAct, SIGNAL(triggered()), this, SLOT(load()));

	tbExportAct = new QAction(this);
	tbExportAct->setIcon(QIcon("../../resources/export.png"));
	tbExportAct->setText("Export Mesh");
	toolBar->addAction(tbExportAct);

	toolBar->addSeparator();

	showFaceAct = new QAction(this);
	showFaceAct->setIcon(QIcon("../../resources/show_face.png"));
	showFaceAct->setText("Show Face");
	showFaceAct->setCheckable(true);
	showFaceAct->setChecked(true);
	toolBar->addAction(showFaceAct);
	connect(showFaceAct, SIGNAL(triggered()), this, SLOT(showFaceSlot()));

	showEdgeAct = new QAction(this);
	showEdgeAct->setIcon(QIcon("../../resources/show_edge.png"));
	showEdgeAct->setText("Show Wireframe");
	showEdgeAct->setCheckable(true);
	showEdgeAct->setChecked(false);
	toolBar->addAction(showEdgeAct);
	connect(showEdgeAct, SIGNAL(triggered()), this, SLOT(showEdgeSlot()));

	toolBar->addSeparator();

	tbSelectFaceAct = new QAction();
	tbSelectFaceAct->setIcon(QIcon("../../resources/select_face.png"));
	tbSelectFaceAct->setText("Select Face");
	tbSelectFaceAct->setCheckable(true);
	tbSelectFaceAct->setChecked(false);
	toolBar->addAction(tbSelectFaceAct);
	connect(tbSelectFaceAct, SIGNAL(triggered()), this, SLOT(faceSelectionSlot()));
	
	tbSelectEdgeAct = new QAction();
	tbSelectEdgeAct->setIcon(QIcon("../../resources/select_edge.png"));
	tbSelectEdgeAct->setText("Select Edge");
	tbSelectEdgeAct->setCheckable(true);
	tbSelectEdgeAct->setChecked(false);
	toolBar->addAction(tbSelectEdgeAct);
	connect(tbSelectEdgeAct, SIGNAL(triggered()), this, SLOT(edgeSelectionSlot()));

	tbSelectVertAct = new QAction();
	tbSelectVertAct->setIcon(QIcon("../../resources/select_vertex.png"));
	tbSelectVertAct->setText("Select Vertex");
	tbSelectVertAct->setCheckable(true);
	tbSelectVertAct->setChecked(false);
	toolBar->addAction(tbSelectVertAct);
	connect(tbSelectVertAct, SIGNAL(triggered()), this, SLOT(vertSelectionSlot()));
}

// Slots
void Viewer::load()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open"), "",
		tr("WaveFront Obj (*.obj);;All Files (*)"));

	if (fileName.isEmpty())
		return;
	else 
	{
		string fn = fileName.toLocal8Bit().constData();
		ObjReader<Kernel, Enriched_items> reader(fn);
		reader.read();

		mesh = new Enriched_polyhedron<Kernel, Enriched_items>;
		ObjBuilder<HalfedgeDS> builder(reader.vertices(), reader.facets());
		mesh->delegate(builder);
		mesh->basic_init();

		canvas->setMesh(mesh);
		canvas->setKDTree();
		canvas->setRenderContexts();
		canvas->setCamera();
		canvas->update();
	}

}

void Viewer::showFaceSlot()
{
	canvas->setFaceVisible(showFaceAct->isChecked());
}

void Viewer::showEdgeSlot()
{
	canvas->setEdgeVisible(showEdgeAct->isChecked());
}

void Viewer::faceSelectionSlot()
{
	if (tbSelectFaceAct->isChecked())
	{
		tbSelectEdgeAct->setChecked(false);
		tbSelectVertAct->setChecked(false);

		canvas->setFaceSelection(true);
		canvas->setEdgeSelection(false);
		canvas->setVertSelection(false);
	}
	else
	{
		canvas->setFaceSelection(false);
	}
}

void Viewer::edgeSelectionSlot()
{
	if (tbSelectEdgeAct->isChecked())
	{
		tbSelectFaceAct->setChecked(false);
		tbSelectVertAct->setChecked(false);

		canvas->setFaceSelection(false);
		canvas->setEdgeSelection(true);
		canvas->setVertSelection(false);
	}
	else
	{
		canvas->setEdgeSelection(false);
	}
}

void Viewer::vertSelectionSlot()
{
	if (tbSelectVertAct->isChecked())
	{
		tbSelectFaceAct->setChecked(false);
		tbSelectEdgeAct->setChecked(false);

		canvas->setFaceSelection(false);
		canvas->setEdgeSelection(false);
		canvas->setVertSelection(true);
	}
	else
	{
		canvas->setVertSelection(false);
	}
}

void Viewer::showVectorFieldSlot()
{
	canvas->clear();

	std::cout<<"Clear. "<<std::endl;
	//canvas->visualizeTriangulation();
	//std::cout<<"Visial Data Done."<<std::endl;
	
	canvas->setRenderContexts();
	cout << "Rendering Context Complete" << endl;
	canvas->setCamera();
	canvas->update();

	canvas->setVectorFieldVisible(true); // comment this if do not display vector field.
}

void Viewer::clearSceneSlot()
{
	canvas->clear();
	canvas->update();
}

void Viewer::decomposeSlot()
{
	canvas->computeTriangulation_2();
	cout << "Triangulation Complete ..." << endl;
}

void Viewer::analyzeLaplacian0Slot()
{
	canvas->computeLaplacian0();
}

void Viewer::analyzeLaplacian1Slot()
{
	canvas->computeLaplacian1();
}

void Viewer::analyzeLaplacian2Slot()
{
	canvas->computeLaplacian2();
}

void Viewer::decomposeNormalSlot()
{
	canvas->computeDecompositionNormal();
}

void Viewer::decomposeTangentialSlot()
{
	canvas->computeDecompositionTangential();
}

void Viewer::poelkeSlot()
{
	canvas->computePoelke();
}

void Viewer::mixedBoundarySlot()
{
	canvas->computeMixedBoundary();
}

void Viewer::flexibilitySlot()
{
	canvas->analyzeFlexibility();
}