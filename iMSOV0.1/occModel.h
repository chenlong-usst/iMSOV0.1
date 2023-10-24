#pragma once

#include "vector"

using namespace std;
class  point {
public:
	double x, y, z, w;
};


//NURBS����
struct OccCurve 
{
	int m_Degree;//���ߴ���
	vector<double> m_Knots;//���߽ڵ�ʸ��
	vector<point> m_CtrlPts;//���Ƶ�

};

//NURBS����
struct OccSurface 
{
	int m_uDegree;//v�������
	int m_vDegree;//v�������
	int m_uNum;//v������Ƶ�����
	int m_vNum;//v������Ƶ�����
	vector<double> m_uKnots;//v����ڵ�ʸ��
	vector<double> m_vKnots;//v����ڵ�ʸ��
	vector<point> m_CtrlPts;//���Ƶ㣬��v-u����洢

	void SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum,
		const vector<double>& uKnots, const vector<double>& vKnots);
};

//NURBS��
struct OccVol 
{
	int m_uDegree;//v�������
	int m_vDegree;//v�������
	int m_wDegree;//w�������
	int m_uNum;//v������Ƶ�����
	int m_vNum;//v������Ƶ�����
	int m_wNum;//w������Ƶ�����
	vector<double> m_uKnots;//v����ڵ�ʸ��
	vector<double> m_vKnots;//v����ڵ�ʸ��
	vector<double> m_wKnots;//w����ڵ�ʸ��
	vector<point> m_CtrlPts;//���Ƶ㣬��v-u-w����洢

};

