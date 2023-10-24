#pragma once


#include <iostream>
#include <vector>
#include <thread>
#include <future>
#include <vector>
#include <list>
#include <stack>
#include <string>
#include <math.h>
#include <algorithm>
#include <memory>

//公用全局头文件
#include "pointXd.h"
#include "varray.h"
#include "XVec.h"
#include "globalFunc.h"

#include "threadStruct.h"

#ifndef minDisBetweenSurfce
#define minDisBetweenSurfce 0.01
#endif

using namespace std;
using namespace base;

class  VolumeVertex//体模型中的一个控制点所有信息
{
public:
	Vec4  m_pt;      //坐标点
	Vec4  m_oript;   //原始控制点，用于备份形体变性前的控制点。
	Vec4  m_matval;  //材料参数
	Vec4  m_displacement;  //在施加载荷后形成的位移。  分别表示ux，uy,uz
	Vec4  m_stress;    //改点外力载荷。     分别表示sigma(x),sigma(y),sigma(z).
	float m_fvonMissStressMat, m_fDisplaceMat;    //表示根据最大最小值得到的应力颜色值和位移颜色值，范围在[0,1].需要显示位移或者应力时，只需将相应的值赋予颜色变量即可。
public:
	VolumeVertex();
	VolumeVertex(Vec3& vt);
	VolumeVertex(Vec3& vt, Vec3& matt);
	VolumeVertex(Vec4& vt);
	VolumeVertex(Vec4& vt, Vec4& matt);
	float getVonmissStress() { return sqrt((m_stress.x - m_stress.y)*(m_stress.x - m_stress.y) / 2 + (m_stress.y - m_stress.z)*(m_stress.y - m_stress.z) / 2 + (m_stress.z - m_stress.x)*(m_stress.z - m_stress.x) / 2); };
	~VolumeVertex();
};


class NurbsLine
{
public:
	NurbsLine();
	NurbsLine(const NurbsLine &);
	NurbsLine& operator=(const NurbsLine&);
	//控制点个数，阶数，节点矢量，控制点
	NurbsLine(int, int, vector<double>, vector<point4d>);
	virtual ~NurbsLine();

	int GetUDegree()const;
	shared_ptr<vector<double>> GetUKonts()const;
	int GetUNum()const;
	shared_ptr<vector<point4d>> GetControlPointer()const;

	bool SetControlPoint(const vector<point4d>);
	bool SetControlPoint(shared_ptr<vector<point4d>>);
	bool SetUDegree(const int);
	bool SetUKonts(const vector<double>);
	bool SetUKonts(shared_ptr<vector<double>>);
	bool SetUNum(const int);
	//检查参数是否正确
	//控制点个数+1+阶数=节点矢量个数
	bool isParameterCorrect(int, int, int);

	//两点绘制直线
	void CreatLineWithTwoPoints(const point4d &p1, const point4d &p2, int degree = 2);

	/*计算节点下标
	x：节点值
	degree：次数
	CtrlPtsNum：控制点数量
	knots：节点矢量*/
	static int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots);
	static int FindSpan(const double x, const int degree, const int CtrlPtsNum, const varray<double>& knots);

	/*根据参数值，节点下标，计算基函数
	u：参数值
	k：节点下标
	degree：次数
	knots：节点矢量
	N：返回的基函数*/
	static void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,vector<double>& N);
	static void BasisFuns(const double u, const int k, const int degree, const varray<double>& knots, varray<double>& N);

	/*根据参数值，节点下标，计算一次基函数*/
	float OneBasisFun(int p, int m, const varray<double>& U, int i, double u);

	/*u处所有基函数
	u：参数值
	k：节点下标
	degree：次数
	knots：节点矢量
	ndu：返回的所有基函数*/
	void AllBasisFuns(const double u, const int k, const int degree,
		const vector<double>& knots, vector<vector<double>>& ndu);

	/*基函数n阶导
	u：参数值
	k：节点下标
	degree：次数
	n：导矢阶数
	knots：节点矢量
	basisDu：基函数n阶导*/
	static void DerBasisFuns(const double u, const int k, const int degree, const int n, const vector<double>& knots,
		vector<vector<double>>& basisDu);
	static void DerBasisFuns(const double u, const int k, const int degree, const int n, const varray<double>& knots,
		varray<varray<double>>& basisDu);


	/*计算u节点对应的曲线坐标点
	u：节点参数*/
	point3d GetLinePoint(const double u);

	/*获取线上的点
	Unum：曲线点数量
	linePts：曲线点	*/
	void CalLinePoint(const int Unum, vector<point3d>& linePts);

	/*曲线矢值函数A(u)所有p阶导，矢值函数A(u)即NURBS曲线分子式
	u：参数
	n：导矢阶数
	Der：A(u)所有p阶导*/
	void PtsDerivsAu(const double u, const int n, vector<point4d>& Der);

	/*曲线上u处所有n阶导
	u：参数
	n：导矢阶数
	Der：所有n阶导 */
	void PtsDerivs(const double u, const int n, vector<point4d>& Der);

	/*曲线所有n阶导
	step：参数域等分步进量
	Der：曲线点对应的n阶导*/
	void CurveDerivs(const int n, const double step, vector<vector<point3d>>& Der);

	/*节点插入
	u：需要插入的节点
	r：插入次数*/
	void KnotInsert(const double u, const int r);

	/*节点去除
	u：需要删除的节点
	r：删除次数*/
	int KnotRemove(const double u, const int r);

	/*节点细化
	u:需要插入的节点*/
	virtual void KnotsRefine(const vector<double>& u);

	/* 曲线升阶
	degree：升阶后次数*/
	virtual void DegreeElevate(const int degree);

	//在节点u处对曲线进行分割
	void Segmentation(const double u, vector<NurbsLine>& lines);

	/*NURBS曲线分解为BEZIER
	Qw：输出的BEZIER曲线段控制点*/
	void Decompose(varray<varray<point4d>>& Qw);

	//反转曲线方向
	void CruveReverse();	
public:
	//u方向次数
	int _u_Degree;
	//u方向节点矢量
	shared_ptr<vector<double>> _u_Knots;
	//u方向控制点个数
	int _u_Num;
	//控制点
	shared_ptr<vector<point4d>> _ControlPts;
};

class NurbsSurface
{
public:
	NurbsSurface();
	NurbsSurface(int, int, int, int, vector<double>, vector<double>, vector<point4d>);
	NurbsSurface(const NurbsSurface&);
	NurbsSurface& operator=(const NurbsSurface&);
	virtual ~NurbsSurface();

	int GetUDegree()const;
	shared_ptr<vector<double>> GetUKonts()const;
	int GetUNum()const;
	shared_ptr<vector<point4d>> GetControlPointer()const;

	int GetVDegree()const;
	shared_ptr<vector<double>> GetVKonts()const;
	int GetVNum()const;
	point4d getControlPoint(int u, int v);

	bool SetControlPoint(const vector<point4d>);
	bool SetControlPoint(shared_ptr<vector<point4d>>);
	bool SetUDegree(const int);
	bool SetUKonts(const vector<double>);
	bool SetUKonts(shared_ptr<vector<double>>);
	bool SetUNum(const int);

	bool SetVDegree(const int);
	bool SetVKonts(const vector<double>);
	bool SetVKonts(shared_ptr<vector<double>>);
	bool SetVNum(const int);

	void SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum,
		const vector<double>& uKnots, const vector<double>& vKnots);

	/*根据参数值，节点下标，计算一次基函数*/
	float OneBasisFun(int p, int m, const varray<double>& U, int i, double u);

	/*u处所有基函数
	u：参数值
	k：节点下标
	degree：次数
	knots：节点矢量
	ndu：返回的所有基函数*/
	void AllBasisFuns(const double u, const int k, const int degree, const vector<double>& knots, vector<vector<double>>& ndu);

	/*计算节点下标
	x：节点值
	degree：次数
	CtrlPtsNum：控制点数量
	knots：节点矢量*/
	int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const;

	/*根据参数值，节点下标，计算基函数
	u：参数值
	k：节点下标
	degree：次数
	knots：节点矢量
	N：返回的基函数*/
	void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
		vector<double>& N)const;

	//判断两个面是否相同
	bool isTwoSurfaceSame(NurbsSurface);

	//返回两个面之间的哈斯多夫距离
	float GetDistanceBetweenTwoSurface(NurbsSurface);

	//计算（u，v）对应的曲面上的点
	point3d GetSurFacePoint(const double u, const double v)const;

	//计算四边形面片显示数据
	//num:该方向采样点数量
	//quads:面片数据
	//lines:等参线数据
	threadParam CalQuads(const int Unum, const int Vnum)const;

	/*Coons插值
	EndgCtrlPts:边界控制点（e0,e1,e2,e3顺序）*/
	void CoonsInterpolate(const varray<varray<point4d>>& EdgeCtrlPts);

	/*Coons插值
	EdgeLines:边界曲线（e0,e1,e2,e3顺序）*/
	void CoonsInterpolate(const varray<NurbsLine>& EdgeLines);

	//曲面升阶
	//Udegree,Vdegree:升阶后次数
	void DegreeElevate(const int Udegree, const int Vdegree);

	//曲面节点插入
	//Uknot,Vknot:需要插入的节点
	void KnotsRefine(const varray<double>& Uknot, const varray<double>& Vknot);

	//根据控制点二维序号计算一维序号
	int CtrlPtsIdx(const int uIdx, const int vIdx);

	//控制点排序转换为U-V
	void OrderCtrlPts();

	//控制点排序转换为U-V
	void OrderCtrlPts(NurbsSurface& sf);

	//提取边界线,排序已经是U-V情况下
	void GetEdgeLines(varray<NurbsLine>& EdgeLines);

public:
	//u方向次数
	int _u_Degree;
	//u方向节点矢量
	shared_ptr<vector<double>> _u_Knots;
	//u方向控制点个数
	int _u_Num;
	//控制点
	shared_ptr<vector<point4d>> _ControlPts;
	//v方向次数
	int _v_Degree;
	//v方向节点矢量
	shared_ptr<vector<double>> _v_Knots;
	//v方向控制点个数
	int _v_Num;
};

class NurbsVol
{
public:
	NurbsVol();
	~NurbsVol();
	NurbsVol(const NurbsVol &);
	NurbsVol& operator=(const NurbsVol&);

	int GetUDegree()const;
	shared_ptr<vector<double>> GetUKonts()const;
	int GetUNum()const;
	shared_ptr<vector<point4d>> GetControlPointer()const;

	int GetVDegree()const;
	shared_ptr<vector<double>> GetVKonts()const;
	int GetVNum()const;
	point4d getControlPoint(int u, int v);
		
	shared_ptr<vector<double>> GetWKonts()const;
	int GetWNum()const;
	point4d getControlPoint(int u, int v, int w);

	bool SetControlPoint(const vector<point4d>);
	bool SetControlPoint(shared_ptr<vector<point4d>>);
	bool SetUDegree(const int);
	bool SetUKonts(const vector<double>);
	bool SetUKonts(shared_ptr<vector<double>>);
	bool SetUNum(const int);

	bool SetVDegree(const int);
	bool SetVKonts(const vector<double>);
	bool SetVKonts(shared_ptr<vector<double>>);
	bool SetVNum(const int);

	int GetWDegree()const;
	bool SetWDegree(const int);
	bool SetWKonts(const vector<double>);
	bool SetWKonts(shared_ptr<vector<double>>);
	bool SetWNum(const int);
	void SetVol(const int uDegree, const int vDegree, const int wDegree, const int uNum, const int vNum, const int wNum,
		const vector<double>& uKnots, const vector<double>& vKnots, const vector<double>& wKnots);

	/*计算节点下标
	x：节点值
	degree：次数
	CtrlPtsNum：控制点数量
	knots：节点矢量*/
	int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const;

	/*根据参数值，节点下标，计算基函数
	u：参数值
	k：节点下标
	degree：次数
	knots：节点矢量
	N：返回的基函数*/
	void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
		vector<double>& N)const;

	/*根据参数值，节点下标，计算一次基函数*/
	float OneBasisFun(int p, int m, const varray<double>& U, int i, double u);

	/*u处所有基函数
	u：参数值
	k：节点下标
	degree：次数
	knots：节点矢量
	ndu：返回的所有基函数*/
	void AllBasisFuns(const double u, const int k, const int degree,
		const vector<double>& knots, vector<vector<double>>& ndu);

	//(u,v,w)对应的体上的点,
	point3d GetVolPoint(const double u, const double v, const double w)const;


	//计算四边形面片显示数据
	//num:该方向采样点数量
	//quads:面片数据
	//lines:等参线数据
	threadParamVOL CalQuads(const int Unum, const int Vnum, const int Wnum,
		varray<varray<varray<point3d>>>& quads, varray<varray<varray<point3d>>>& lines)const;

	threadParamVOL CalQuads2(const int Unum, const int Vnum, const int Wnum)const;


	//根据控制点三维序号计算一维序号
	int CtrlPtsIdx(const int uIdx, const int vIdx, const int wIdx);

	//控制点排序转换为U-V-W
	void OrderCtrlPts();

	//控制点排序转换为U-V-W
	void OrderCtrlPts(NurbsVol& vol);

	//节点插入
	//Uknot,Vknot,Wknot:需要插入的节点
	void KnotsRefine(varray<double> Uknot, varray<double> Vknot, varray<double> Wknot);

	//升阶
	//Udegree,Vdegree,Wdegree:升阶后阶数
	void DegreeElevate(const int Udegree, const int Vdegree, const int Wdegree);

	/*扫描生成Nurbs体模型，截面垂直于路径
	pathT：扫描路径
	nurbsSF：起始截面
	K：截面实例数量,一般取路径控制点数量减1
	*/
	void CreateSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface& nurbsSF, const int K);

	/*平移扫描生成Nurbs体模型，截面不垂直于路径
	pathT：扫描路径
	nurbsSF：起始截面*/
	void CreateTransSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface& nurbsSF);

	//放样
	//path:路径
	//surfaces:放样截面
	void LoftingNurbsVol(const NurbsLine& path, const varray<NurbsSurface>& surfaces);

	//放样
	//path:路径
	//surfaces0:路径起点放样截面
	//surfaces1:路径终点放样截面
	void LoftingNurbsVol(const NurbsLine& path, const NurbsSurface& surfaces0, const NurbsSurface& surfaces1);

	//测试网格质量
	void  TestMeshQuanity(int segmentNum);

	//提取单个边界面
		//dir:1->6表示方向
	NurbsSurface GetSingleSurfaces(int dir) const;

private:
	//计算等参面
	//uvw:0=u,1=v,2=w
	//t:参数
	//num:以u-v-w顺序
	//L:等参面四边形面片集
	void CalIsoSurface(const int uvw, const double t, const int num1, const int num2,
		varray<varray<point3d>>& quads, varray<varray<point3d>>& lines)const;

	/*扫描时的截面实例位置
	pathT：扫描路径
	K：截面实例数量
	pos：截面实例位置
	NewKnots：新节点矢量*/
	void InsLocation(const NurbsLine & pathT, int K, varray<double> & pos);

	/*计算扫描路径上的局部坐标系
	pathT：扫描路径
	pos：截面实例位置
	TranMat：pos处的局部坐标系*/
	void LocalCoordinates(const NurbsLine & pathT, const varray<double> & pos,
		varray<varray<point4d>> & TranMat);

	/*沿扫描路径对截面进行变换
	nurbsSF：截面控制点
	TranMat：变换矩阵
	OriCoordinates:截面的局部坐标系
	allNurbsSF：得到的所有截面实例控制点
	*/
	void MatrixTran(const varray<varray<point4d>> & nurbsSF, const varray<varray<point4d>>& TranMat,
		varray<varray<varray<point4d>>>& allNurbsSF);

	//扫描生成Nurbs体
	void SweepSurface(const varray<varray<varray<point4d>>>& allNurbsSF,
		const varray<varray<varray<double>>>& SFw, const varray<double>& pos);

	//取最高次
	void MaxDegree(const varray<NurbsSurface>& surfaces, int& uDegree, int& vDegree);

	//节点矢量并集
	void KnotsUnify(const varray<NurbsSurface>& surfaces, varray<double>& NewUKnots, varray<double>& NewVKnots);

	//获得单个点的雅可比值
	float GetJacobianValue(float u, float v, float w);

	float GetFirstPartialDirivatives(float u, float v, float w, Vec3& pu, Vec3& pv, Vec3& pw);

public:
	//u方向次数
	int _u_Degree;
	//u方向节点矢量
	shared_ptr<vector<double>> _u_Knots;
	//u方向控制点个数
	int _u_Num;
	//控制点
	shared_ptr<vector<point4d>> _ControlPts;
	//v方向次数
	int _v_Degree;
	//v方向节点矢量
	shared_ptr<vector<double>> _v_Knots;
	//v方向控制点个数
	int _v_Num;
	//w方向次数
	int _w_Degree;
	//w方向节点矢量
	shared_ptr<vector<double>> _w_Knots;
	//w方向控制点个数
	int _w_Num;

	//雅可比控制值
	int m_minorJocbianNum;
	float m_isovalAdd;
	float m_minJocbian;
	float m_maxJocbian;


	

	//20230715 直接复制于SplineVolume.h    ，还需要修改
	//int   m_uNum, m_uDegree, m_uRenderNum;
	//int   m_vNum, m_vDegree, m_vRenderNum;
	//int   m_wNum, m_wDegree, m_wRenderNum;
	int m_uRenderNum, m_vRenderNum, m_wRenderNum;
	bool  m_bHeterogeneous;//？？？
	//varray<Vec4>			 m_CtrlPts;
	varray<VolumeVertex> m_vAllCtrlPts;  //控制点
	varray<VolumeVertex> m_vVolumePts;   //用于显示的节点。
	varray<int>  m_CtrlPtsIDinKKKMatrix;   //每个控制点在总纲矩阵中的位置。 如果是单片，默认就是0~m_uNum*m_vNum*m_wNum-1;如果是多片就需要重新编号。
	bool  m_bHasInterpolated;//？？？
	//varray<double> m_uKnots;
	//varray<double> m_vKnots;
	//varray<double> m_wKnots;
	varray<NurbsSurface> m_6boundarySurface;  //排序前后都是对应m_6PatchIdx，对应于拟合的surface,存储的是旋转后的曲面
	bool  m_bHasIGASolution;
	//float m_maxJocbian, m_minJocbian;
	//int   m_minorJocbianNum;
	//double m_isovalAdd;
	//以下用于IGA3D
	varray<point4d>m_CPdu;
	varray<double>m_dis;
	varray<double>m_stress;
	varray<varray<varray<Vec3>>>m_show3Dcpts;//用于3dIGA显示的控制点
	varray<varray<varray<double>>>m_show3DDisp;//用于3dIGA显示的位移
	varray<varray<varray<Vec3>>>m_show3DEdge;//用于3dIGA显示边界线
	varray<varray<varray<Vec3>>>m_show3DcptsOfStress;//用于3dIGA显示的控制点
	varray<varray<varray<double>>>m_show3DStrainOrStress;//用于3dIGA显示的应力
};