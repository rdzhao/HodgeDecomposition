#include <stack>
#include <bitset>
#include <unordered_map>
#include <functional>
#include <utility>
#include <algorithm>
#include <random>
#include <cmath>
#include <chrono>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include "engine.h"

// Eigen 
typedef Eigen::SparseMatrix<double> EigenSpMat;
typedef Eigen::MatrixXd EigenMatrix;
typedef Eigen::VectorXd EigenVector;
typedef Eigen::Matrix3d EigenMatrix3d;
typedef Eigen::Matrix4d EigenMatrix4d;
typedef Eigen::Vector3d EigenVector3d;
typedef Eigen::Vector4d EigenVector4d;

class JHTDB
{
public:
	JHTDB() {};

	EigenMatrix getVelocities(EigenMatrix pts);
};