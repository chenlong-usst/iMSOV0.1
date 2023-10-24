#pragma once

#include "vector"

using namespace std;
class  point {
public:
	double x, y, z, w;
};


//NURBS曲线
struct OccCurve 
{
	int m_Degree;//曲线次数
	vector<double> m_Knots;//曲线节点矢量
	vector<point> m_CtrlPts;//控制点

};

//NURBS曲面
struct OccSurface 
{
	int m_uDegree;//v方向次数
	int m_vDegree;//v方向次数
	int m_uNum;//v方向控制点数量
	int m_vNum;//v方向控制点数量
	vector<double> m_uKnots;//v方向节点矢量
	vector<double> m_vKnots;//v方向节点矢量
	vector<point> m_CtrlPts;//控制点，按v-u方向存储

	void SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum,
		const vector<double>& uKnots, const vector<double>& vKnots);
};

//NURBS体
struct OccVol 
{
	int m_uDegree;//v方向次数
	int m_vDegree;//v方向次数
	int m_wDegree;//w方向次数
	int m_uNum;//v方向控制点数量
	int m_vNum;//v方向控制点数量
	int m_wNum;//w方向控制点数量
	vector<double> m_uKnots;//v方向节点矢量
	vector<double> m_vKnots;//v方向节点矢量
	vector<double> m_wKnots;//w方向节点矢量
	vector<point> m_CtrlPts;//控制点，按v-u-w方向存储

};

