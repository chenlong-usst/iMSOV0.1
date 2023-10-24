//ȫ�ֺ���
#pragma once

#include <vector>
#include <string>

#include "Eigen/Dense"

#include "varray.h"
#include "pointXd.h"
#include "XVec.h"

#ifdef IMSOBASE_EXPORT
#define GLOBAL_API _declspec(dllexport)
#else
#define GLOBAL_API _declspec(dllimport)
#endif

using Eigen::Vector4d;
using Eigen::Matrix4d;
using std::string;
using std::vector;
using namespace base;

class GLOBAL_API Global {
public:

	/////////////////////////////////////////////���ı�    2018///////////////////////////////
	/*ά�Ƚ�һ
	hightDimVarray����άvarray
	lowDimVarray����άvarray
	clear���Ƿ������άvarrayԭ������
	*/
	//template<class _T>
	static void ReduceVarrayDim(const varray<varray<point4d>>& hightDimVarray, varray<point4d>& lowDimVarray, const bool clear);

	//��ά����ת��
	template<class _T>
	void Transpose(varray<varray<_T>>& Varr);

	/*����ʽϵ��
	Bin(x,y)=x!/(y!*(x-y)!)
	*/
	static double Bin(const int x, const int y);

	//vectorת��Ϊvarray
	//template<typename _T>
	//static varray<_T> transVectorToVarray(vector<_T> a);

	static varray<double> transVectorToVarray(vector<double> a);
	static varray<point4d> transVectorToVarray(vector<point4d> a);
	static varray<point3d> transVectorToVarray(vector<point3d> a);

	//varrayת��Ϊvector
	//template<typename _T>
	//vector<_T> transVarrayToVector(varray<_T> a);
	static vector<double> transVarrayToVector(varray<double> a);
	static vector<point4d> transVarrayToVector(varray<point4d> a);
	static vector<point3d> transVarrayToVector(varray<point3d> a);

	//��ʾ��淶����[a��b]
	void Normalization(varray<point3d>& pts, const double a, const double b);

	//��ʾ��淶����[a��b]
	void Normalization(varray<varray<point3d>>& pts, const double a, const double b);

	//����㼯����
	static void GetCenter(varray<point3d>& pts, point3d & center);

	//�㼯����ƽ��
	//mode=1 vecΪƽ���յ�,��ʱ��Ҫ����㼯����center
	//mode=0 vecΪƽ������
	template<class T>
	static void TransPts(varray<T> &pts, const T &vec, bool mode, const T &center/* = T()*/);

	//���ʽת��
	//point4dתVector4d
	static Vector4d P4dToV4d(const point4d P4d);

	static Vec4 Point4dToVec4(const point4d p4d);
	static varray<Vec4>Point4dToVec4(const varray<point4d> VP4d);

	//���ʽת��
	//Vector4dתpoint4d
	static point4d V4dToP4d(const Vector4d V4d);

	//Vector4d��λ��
	void Vector4dToUnit(Vector4d &v);

	//Vector4d���
	void AcrossB(const Vector4d & a, const Vector4d & b, Vector4d & c);

	//����н�
	double angleAB(const Vector4d & va, const Vector4d & vb);

	//�������߳���,LΪ���������
	double CurveLength(const varray<point3d>& L);

	//�ж�a�Ƿ���S��
	//template<class numtype>
	//bool IsInSet(const numtype a, const vector<numtype>& S);
	static bool IsInSet(const double a, const vector<double>& S);

	//��������ϵobj��target�Ĺ��ɾ���
	//Coordinates��X�ᣬY�ᣬZ��, ����ԭ��
	static void CalCoordinatesTransMat(const varray<point3d>& objCoordinates, const varray<point3d>& targetCoordinates, Matrix4d& mat);

	//����ʸ���任����
	static void CalTransMat(const Vector4d & targetPts, const Vector4d & targetVec, const Vector4d & objPts, const Vector4d & objVec, Matrix4d & mat);

	static void CalTransMat(const point4d & targetPts, const point4d & targetVec, const point4d & objPts, const point4d & objVec, Matrix4d & mat);

	//����任
	static void TransByMat(varray<Vector4d>& varr, const Matrix4d & mat);

	static void TransByMat(varray<point4d>& varr, const Matrix4d& mat);

	static void TransByMat(varray<Vec4>& varr, const Matrix4d & mat);

	static void TransByMat(vector<point4d>& varr, const Matrix4d& mat);

	//ȡ���ź���
	int Sgn(double a);

	/*����������㼰����
	beginIdx:����Ϊ���Ԫ�صĵ�ǰ�±�
	dir���Ƿ񰴵�ǰ˳�򴢴�
	*/
	template<class _T>
	void SetBeginDir(varray<_T>& varr, const size_t beginIdx, const bool dir/* = true*/);

	//���Լ��
	static int GCD(int a, int b);

	//����
	//template<class _T>
	static void KnotUnify(const varray<double>& KnotsA, const varray<double>& KnotsB, varray<double>& NewKnots);

	//�
	//diffKnots=KnotsL-KnotsS
	//KnotsL��С��С��KnotsS
	//template<class _T>
	static void KnotsDiff(const varray<double>& KnotsL, const varray<double>& KnotsS, varray<double>& diffKnots);

	//��pts������v1ͶӰ��ƽ��(v,p)
	//v1Ĭ��Ϊ0��������ʱΪ��ͶӰ���൱��v1 = -v
	//ƽ�з���0,���򷵻�1
	static int Project2Plane(const point3d& v, const point3d& p, const point3d& pts, point3d& res, const point3d& v1/* = point3d()*/);

	//�㼯Cpts������v1ͶӰ��ƽ��(v,p)
	//v1Ĭ��Ϊ0��������ʱΪ��ͶӰ���൱��v1 = -v
	static int Project2Plane(const point3d & v, const point3d & p, const varray<point4d>& Cpts, varray<point4d>& ProPts, const point3d & v1/* = point3d()*/);




};
