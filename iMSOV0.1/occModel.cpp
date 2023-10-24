#include "occModel.h"

void OccSurface::SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum, 
	const vector<double>& uKnots, const vector<double>& vKnots)
{
	m_uDegree = uDegree;
	m_vDegree = vDegree;
	m_uNum = uNum;
	m_vNum = vNum;
	m_uKnots = uKnots;
	m_vKnots = vKnots;
}
