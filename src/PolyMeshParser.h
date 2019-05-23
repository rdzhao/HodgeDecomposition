#include <stack>
#include <bitset>
#include <map>
#include <utility>
#include <functional>
#include <utility>
#include <algorithm>
#include <random>
#include <cmath>
#include <chrono>

typedef std::vector<double> Point; // coordinates of points
typedef std::vector<int> Face; // vertex indices list of face
typedef std::vector<int> Cell; // vertex indices list of cell
typedef std::pair<int, int> CellPoint; // cell vertex pair for face

class PolyMeshParser 
{
public:
	PolyMeshParser() {};

	void readPolyMesh();
	void parse();
	void write();

private:
	void readPoints();
	void readFaces();
	void readOwners();
	void readNeighbours();

private:
	void initFacesAndCells(); // for boundary and subdomain index
	void buildCellFaceMap();
	void buildCellNeighbourMap();
	void extractCellPointsList();

private:
	void writeCGAL();

private:
	std::vector<Point> m_points;
	std::vector<Face> m_faces;
	std::vector<int> m_owners;
	std::vector<int> m_neighbours;

	std::map<int, bool> m_faceOnBoundary;
	std::map<int, bool> m_cellSubdomainIndex;
	std::map<CellPoint, int> m_cellFaceMap;
	std::map<CellPoint, int> m_cellNeighbourMap;

	std::vector<Point> m_points_CGAL;
	std::vector<Cell> m_cells_CGAL;
};