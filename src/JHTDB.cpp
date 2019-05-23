#include "JHTDB.h"

EigenMatrix JHTDB::getVelocities(EigenMatrix pts)
{
	Engine* eng;
	eng = engOpen("\0");
	std::string str;

	mxArray* ptsArray;
	ptsArray = mxCreateDoubleMatrix(pts.rows(), pts.cols(), mxREAL);
	memcpy(mxGetPr(ptsArray), pts.data(), pts.rows() *pts.cols() * sizeof(double));
	engPutVariable(eng, "points", ptsArray);

	str = "addpath('C:\\Users\\zhaor\\Desktop\\Project\\turbmat');";        engEvalString(eng, str.c_str());
	str = "velocities = VelocityInquiry(points)";		engEvalString(eng, str.c_str());

	mxArray* vv;
	vv = engGetVariable(eng, "velocities");
	double* v;
	v = reinterpret_cast<double*>(mxGetData(vv));
	const size_t* dims;
	dims = mxGetDimensions(vv);

	EigenMatrix velocities;
	velocities.resize(dims[0], dims[1]);


	std::cout << "Dimension: " << dims[0] << " " << dims[1] << std::endl;

	for (int j = 0; j < dims[1]; ++j) {
		for (int i = 0; i < dims[0]; ++i) {
			velocities.coeffRef(i, j) = v[i + j * dims[0]];
		}
	}

	return velocities;
}