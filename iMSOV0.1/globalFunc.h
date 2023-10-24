//全局函数
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

	/////////////////////////////////////////////何文彬    2018///////////////////////////////
	/*维度降一
	hightDimVarray：高维varray
	lowDimVarray：低维varray
	clear：是否清除低维varray原有内容
	*/
	//template<class _T>
	static void ReduceVarrayDim(const varray<varray<point4d>>& hightDimVarray, varray<point4d>& lowDimVarray, const bool clear);

	//二维数组转置
	template<class _T>
	void Transpose(varray<varray<_T>>& Varr);

	/*二项式系数
	Bin(x,y)=x!/(y!*(x-y)!)
	*/
	static double Bin(const int x, const int y);

	//vector转换为varray
	//template<typename _T>
	//static varray<_T> transVectorToVarray(vector<_T> a);

	static varray<double> transVectorToVarray(vector<double> a);
	static varray<point4d> transVectorToVarray(vector<point4d> a);
	static varray<point3d> transVectorToVarray(vector<point3d> a);

	//varray转换为vector
	//template<typename _T>
	//vector<_T> transVarrayToVector(varray<_T> a);
	static vector<double> transVarrayToVector(varray<double> a);
	static vector<point4d> transVarrayToVector(varray<point4d> a);
	static vector<point3d> transVarrayToVector(varray<point3d> a);

	//显示点规范化至[a，b]
	void Normalization(varray<point3d>& pts, const double a, const double b);

	//显示点规范化至[a，b]
	void Normalization(varray<varray<point3d>>& pts, const double a, const double b);

	//计算点集中心
	static void GetCenter(varray<point3d>& pts, point3d & center);

	//点集中心平移
	//mode=1 vec为平移终点,此时需要传入点集中心center
	//mode=0 vec为平移向量
	template<class T>
	static void TransPts(varray<T> &pts, const T &vec, bool mode, const T &center/* = T()*/);

	//点格式转换
	//point4d转Vector4d
	static Vector4d P4dToV4d(const point4d P4d);

	static Vec4 Point4dToVec4(const point4d p4d);
	static varray<Vec4>Point4dToVec4(const varray<point4d> VP4d);

	//点格式转换
	//Vector4d转point4d
	static point4d V4dToP4d(const Vector4d V4d);

	//Vector4d单位化
	void Vector4dToUnit(Vector4d &v);

	//Vector4d叉乘
	void AcrossB(const Vector4d & a, const Vector4d & b, Vector4d & c);

	//计算夹角
	double angleAB(const Vector4d & va, const Vector4d & vb);

	//计算曲线长度,L为曲线坐标点
	double CurveLength(const varray<point3d>& L);

	//判断a是否在S内
	//template<class numtype>
	//bool IsInSet(const numtype a, const vector<numtype>& S);
	static bool IsInSet(const double a, const vector<double>& S);

	//计算坐标系obj到target的过渡矩阵
	//Coordinates：X轴，Y轴，Z轴, 坐标原点
	static void CalCoordinatesTransMat(const varray<point3d>& objCoordinates, const varray<point3d>& targetCoordinates, Matrix4d& mat);

	//计算矢量变换矩阵
	static void CalTransMat(const Vector4d & targetPts, const Vector4d & targetVec, const Vector4d & objPts, const Vector4d & objVec, Matrix4d & mat);

	static void CalTransMat(const point4d & targetPts, const point4d & targetVec, const point4d & objPts, const point4d & objVec, Matrix4d & mat);

	//坐标变换
	static void TransByMat(varray<Vector4d>& varr, const Matrix4d & mat);

	static void TransByMat(varray<point4d>& varr, const Matrix4d& mat);

	static void TransByMat(varray<Vec4>& varr, const Matrix4d & mat);

	static void TransByMat(vector<point4d>& varr, const Matrix4d& mat);

	//取符号函数
	int Sgn(double a);

	/*更改容器起点及方向
	beginIdx:设置为起点元素的当前下标
	dir：是否按当前顺序储存
	*/
	template<class _T>
	void SetBeginDir(varray<_T>& varr, const size_t beginIdx, const bool dir/* = true*/);

	//最大公约数
	static int GCD(int a, int b);

	//并集
	//template<class _T>
	static void KnotUnify(const varray<double>& KnotsA, const varray<double>& KnotsB, varray<double>& NewKnots);

	//差集
	//diffKnots=KnotsL-KnotsS
	//KnotsL大小不小于KnotsS
	//template<class _T>
	static void KnotsDiff(const varray<double>& KnotsL, const varray<double>& KnotsS, varray<double>& diffKnots);

	//点pts沿向量v1投影至平面(v,p)
	//v1默认为0向量，此时为正投影，相当于v1 = -v
	//平行返回0,否则返回1
	static int Project2Plane(const point3d& v, const point3d& p, const point3d& pts, point3d& res, const point3d& v1/* = point3d()*/);

	//点集Cpts沿向量v1投影至平面(v,p)
	//v1默认为0向量，此时为正投影，相当于v1 = -v
	static int Project2Plane(const point3d & v, const point3d & p, const varray<point4d>& Cpts, varray<point4d>& ProPts, const point3d & v1/* = point3d()*/);




};
