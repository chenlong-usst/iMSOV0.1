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

//����ȫ��ͷ�ļ�
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

class  VolumeVertex//��ģ���е�һ�����Ƶ�������Ϣ
{
public:
	Vec4  m_pt;      //�����
	Vec4  m_oript;   //ԭʼ���Ƶ㣬���ڱ����������ǰ�Ŀ��Ƶ㡣
	Vec4  m_matval;  //���ϲ���
	Vec4  m_displacement;  //��ʩ���غɺ��γɵ�λ�ơ�  �ֱ��ʾux��uy,uz
	Vec4  m_stress;    //�ĵ������غɡ�     �ֱ��ʾsigma(x),sigma(y),sigma(z).
	float m_fvonMissStressMat, m_fDisplaceMat;    //��ʾ���������Сֵ�õ���Ӧ����ɫֵ��λ����ɫֵ����Χ��[0,1].��Ҫ��ʾλ�ƻ���Ӧ��ʱ��ֻ�轫��Ӧ��ֵ������ɫ�������ɡ�
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
	//���Ƶ�������������ڵ�ʸ�������Ƶ�
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
	//�������Ƿ���ȷ
	//���Ƶ����+1+����=�ڵ�ʸ������
	bool isParameterCorrect(int, int, int);

	//�������ֱ��
	void CreatLineWithTwoPoints(const point4d &p1, const point4d &p2, int degree = 2);

	/*����ڵ��±�
	x���ڵ�ֵ
	degree������
	CtrlPtsNum�����Ƶ�����
	knots���ڵ�ʸ��*/
	static int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots);
	static int FindSpan(const double x, const int degree, const int CtrlPtsNum, const varray<double>& knots);

	/*���ݲ���ֵ���ڵ��±꣬���������
	u������ֵ
	k���ڵ��±�
	degree������
	knots���ڵ�ʸ��
	N�����صĻ�����*/
	static void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,vector<double>& N);
	static void BasisFuns(const double u, const int k, const int degree, const varray<double>& knots, varray<double>& N);

	/*���ݲ���ֵ���ڵ��±꣬����һ�λ�����*/
	float OneBasisFun(int p, int m, const varray<double>& U, int i, double u);

	/*u�����л�����
	u������ֵ
	k���ڵ��±�
	degree������
	knots���ڵ�ʸ��
	ndu�����ص����л�����*/
	void AllBasisFuns(const double u, const int k, const int degree,
		const vector<double>& knots, vector<vector<double>>& ndu);

	/*������n�׵�
	u������ֵ
	k���ڵ��±�
	degree������
	n����ʸ����
	knots���ڵ�ʸ��
	basisDu��������n�׵�*/
	static void DerBasisFuns(const double u, const int k, const int degree, const int n, const vector<double>& knots,
		vector<vector<double>>& basisDu);
	static void DerBasisFuns(const double u, const int k, const int degree, const int n, const varray<double>& knots,
		varray<varray<double>>& basisDu);


	/*����u�ڵ��Ӧ�����������
	u���ڵ����*/
	point3d GetLinePoint(const double u);

	/*��ȡ���ϵĵ�
	Unum�����ߵ�����
	linePts�����ߵ�	*/
	void CalLinePoint(const int Unum, vector<point3d>& linePts);

	/*����ʸֵ����A(u)����p�׵���ʸֵ����A(u)��NURBS���߷���ʽ
	u������
	n����ʸ����
	Der��A(u)����p�׵�*/
	void PtsDerivsAu(const double u, const int n, vector<point4d>& Der);

	/*������u������n�׵�
	u������
	n����ʸ����
	Der������n�׵� */
	void PtsDerivs(const double u, const int n, vector<point4d>& Der);

	/*��������n�׵�
	step��������ȷֲ�����
	Der�����ߵ��Ӧ��n�׵�*/
	void CurveDerivs(const int n, const double step, vector<vector<point3d>>& Der);

	/*�ڵ����
	u����Ҫ����Ľڵ�
	r���������*/
	void KnotInsert(const double u, const int r);

	/*�ڵ�ȥ��
	u����Ҫɾ���Ľڵ�
	r��ɾ������*/
	int KnotRemove(const double u, const int r);

	/*�ڵ�ϸ��
	u:��Ҫ����Ľڵ�*/
	virtual void KnotsRefine(const vector<double>& u);

	/* ��������
	degree�����׺����*/
	virtual void DegreeElevate(const int degree);

	//�ڽڵ�u�������߽��зָ�
	void Segmentation(const double u, vector<NurbsLine>& lines);

	/*NURBS���߷ֽ�ΪBEZIER
	Qw�������BEZIER���߶ο��Ƶ�*/
	void Decompose(varray<varray<point4d>>& Qw);

	//��ת���߷���
	void CruveReverse();	
public:
	//u�������
	int _u_Degree;
	//u����ڵ�ʸ��
	shared_ptr<vector<double>> _u_Knots;
	//u������Ƶ����
	int _u_Num;
	//���Ƶ�
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

	/*���ݲ���ֵ���ڵ��±꣬����һ�λ�����*/
	float OneBasisFun(int p, int m, const varray<double>& U, int i, double u);

	/*u�����л�����
	u������ֵ
	k���ڵ��±�
	degree������
	knots���ڵ�ʸ��
	ndu�����ص����л�����*/
	void AllBasisFuns(const double u, const int k, const int degree, const vector<double>& knots, vector<vector<double>>& ndu);

	/*����ڵ��±�
	x���ڵ�ֵ
	degree������
	CtrlPtsNum�����Ƶ�����
	knots���ڵ�ʸ��*/
	int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const;

	/*���ݲ���ֵ���ڵ��±꣬���������
	u������ֵ
	k���ڵ��±�
	degree������
	knots���ڵ�ʸ��
	N�����صĻ�����*/
	void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
		vector<double>& N)const;

	//�ж��������Ƿ���ͬ
	bool isTwoSurfaceSame(NurbsSurface);

	//����������֮��Ĺ�˹������
	float GetDistanceBetweenTwoSurface(NurbsSurface);

	//���㣨u��v����Ӧ�������ϵĵ�
	point3d GetSurFacePoint(const double u, const double v)const;

	//�����ı�����Ƭ��ʾ����
	//num:�÷������������
	//quads:��Ƭ����
	//lines:�Ȳ�������
	threadParam CalQuads(const int Unum, const int Vnum)const;

	/*Coons��ֵ
	EndgCtrlPts:�߽���Ƶ㣨e0,e1,e2,e3˳��*/
	void CoonsInterpolate(const varray<varray<point4d>>& EdgeCtrlPts);

	/*Coons��ֵ
	EdgeLines:�߽����ߣ�e0,e1,e2,e3˳��*/
	void CoonsInterpolate(const varray<NurbsLine>& EdgeLines);

	//��������
	//Udegree,Vdegree:���׺����
	void DegreeElevate(const int Udegree, const int Vdegree);

	//����ڵ����
	//Uknot,Vknot:��Ҫ����Ľڵ�
	void KnotsRefine(const varray<double>& Uknot, const varray<double>& Vknot);

	//���ݿ��Ƶ��ά��ż���һά���
	int CtrlPtsIdx(const int uIdx, const int vIdx);

	//���Ƶ�����ת��ΪU-V
	void OrderCtrlPts();

	//���Ƶ�����ת��ΪU-V
	void OrderCtrlPts(NurbsSurface& sf);

	//��ȡ�߽���,�����Ѿ���U-V�����
	void GetEdgeLines(varray<NurbsLine>& EdgeLines);

public:
	//u�������
	int _u_Degree;
	//u����ڵ�ʸ��
	shared_ptr<vector<double>> _u_Knots;
	//u������Ƶ����
	int _u_Num;
	//���Ƶ�
	shared_ptr<vector<point4d>> _ControlPts;
	//v�������
	int _v_Degree;
	//v����ڵ�ʸ��
	shared_ptr<vector<double>> _v_Knots;
	//v������Ƶ����
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

	/*����ڵ��±�
	x���ڵ�ֵ
	degree������
	CtrlPtsNum�����Ƶ�����
	knots���ڵ�ʸ��*/
	int FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const;

	/*���ݲ���ֵ���ڵ��±꣬���������
	u������ֵ
	k���ڵ��±�
	degree������
	knots���ڵ�ʸ��
	N�����صĻ�����*/
	void BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
		vector<double>& N)const;

	/*���ݲ���ֵ���ڵ��±꣬����һ�λ�����*/
	float OneBasisFun(int p, int m, const varray<double>& U, int i, double u);

	/*u�����л�����
	u������ֵ
	k���ڵ��±�
	degree������
	knots���ڵ�ʸ��
	ndu�����ص����л�����*/
	void AllBasisFuns(const double u, const int k, const int degree,
		const vector<double>& knots, vector<vector<double>>& ndu);

	//(u,v,w)��Ӧ�����ϵĵ�,
	point3d GetVolPoint(const double u, const double v, const double w)const;


	//�����ı�����Ƭ��ʾ����
	//num:�÷������������
	//quads:��Ƭ����
	//lines:�Ȳ�������
	threadParamVOL CalQuads(const int Unum, const int Vnum, const int Wnum,
		varray<varray<varray<point3d>>>& quads, varray<varray<varray<point3d>>>& lines)const;

	threadParamVOL CalQuads2(const int Unum, const int Vnum, const int Wnum)const;


	//���ݿ��Ƶ���ά��ż���һά���
	int CtrlPtsIdx(const int uIdx, const int vIdx, const int wIdx);

	//���Ƶ�����ת��ΪU-V-W
	void OrderCtrlPts();

	//���Ƶ�����ת��ΪU-V-W
	void OrderCtrlPts(NurbsVol& vol);

	//�ڵ����
	//Uknot,Vknot,Wknot:��Ҫ����Ľڵ�
	void KnotsRefine(varray<double> Uknot, varray<double> Vknot, varray<double> Wknot);

	//����
	//Udegree,Vdegree,Wdegree:���׺����
	void DegreeElevate(const int Udegree, const int Vdegree, const int Wdegree);

	/*ɨ������Nurbs��ģ�ͣ����洹ֱ��·��
	pathT��ɨ��·��
	nurbsSF����ʼ����
	K������ʵ������,һ��ȡ·�����Ƶ�������1
	*/
	void CreateSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface& nurbsSF, const int K);

	/*ƽ��ɨ������Nurbs��ģ�ͣ����治��ֱ��·��
	pathT��ɨ��·��
	nurbsSF����ʼ����*/
	void CreateTransSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface& nurbsSF);

	//����
	//path:·��
	//surfaces:��������
	void LoftingNurbsVol(const NurbsLine& path, const varray<NurbsSurface>& surfaces);

	//����
	//path:·��
	//surfaces0:·������������
	//surfaces1:·���յ��������
	void LoftingNurbsVol(const NurbsLine& path, const NurbsSurface& surfaces0, const NurbsSurface& surfaces1);

	//������������
	void  TestMeshQuanity(int segmentNum);

	//��ȡ�����߽���
		//dir:1->6��ʾ����
	NurbsSurface GetSingleSurfaces(int dir) const;

private:
	//����Ȳ���
	//uvw:0=u,1=v,2=w
	//t:����
	//num:��u-v-w˳��
	//L:�Ȳ����ı�����Ƭ��
	void CalIsoSurface(const int uvw, const double t, const int num1, const int num2,
		varray<varray<point3d>>& quads, varray<varray<point3d>>& lines)const;

	/*ɨ��ʱ�Ľ���ʵ��λ��
	pathT��ɨ��·��
	K������ʵ������
	pos������ʵ��λ��
	NewKnots���½ڵ�ʸ��*/
	void InsLocation(const NurbsLine & pathT, int K, varray<double> & pos);

	/*����ɨ��·���ϵľֲ�����ϵ
	pathT��ɨ��·��
	pos������ʵ��λ��
	TranMat��pos���ľֲ�����ϵ*/
	void LocalCoordinates(const NurbsLine & pathT, const varray<double> & pos,
		varray<varray<point4d>> & TranMat);

	/*��ɨ��·���Խ�����б任
	nurbsSF��������Ƶ�
	TranMat���任����
	OriCoordinates:����ľֲ�����ϵ
	allNurbsSF���õ������н���ʵ�����Ƶ�
	*/
	void MatrixTran(const varray<varray<point4d>> & nurbsSF, const varray<varray<point4d>>& TranMat,
		varray<varray<varray<point4d>>>& allNurbsSF);

	//ɨ������Nurbs��
	void SweepSurface(const varray<varray<varray<point4d>>>& allNurbsSF,
		const varray<varray<varray<double>>>& SFw, const varray<double>& pos);

	//ȡ��ߴ�
	void MaxDegree(const varray<NurbsSurface>& surfaces, int& uDegree, int& vDegree);

	//�ڵ�ʸ������
	void KnotsUnify(const varray<NurbsSurface>& surfaces, varray<double>& NewUKnots, varray<double>& NewVKnots);

	//��õ�������ſɱ�ֵ
	float GetJacobianValue(float u, float v, float w);

	float GetFirstPartialDirivatives(float u, float v, float w, Vec3& pu, Vec3& pv, Vec3& pw);

public:
	//u�������
	int _u_Degree;
	//u����ڵ�ʸ��
	shared_ptr<vector<double>> _u_Knots;
	//u������Ƶ����
	int _u_Num;
	//���Ƶ�
	shared_ptr<vector<point4d>> _ControlPts;
	//v�������
	int _v_Degree;
	//v����ڵ�ʸ��
	shared_ptr<vector<double>> _v_Knots;
	//v������Ƶ����
	int _v_Num;
	//w�������
	int _w_Degree;
	//w����ڵ�ʸ��
	shared_ptr<vector<double>> _w_Knots;
	//w������Ƶ����
	int _w_Num;

	//�ſɱȿ���ֵ
	int m_minorJocbianNum;
	float m_isovalAdd;
	float m_minJocbian;
	float m_maxJocbian;


	

	//20230715 ֱ�Ӹ�����SplineVolume.h    ������Ҫ�޸�
	//int   m_uNum, m_uDegree, m_uRenderNum;
	//int   m_vNum, m_vDegree, m_vRenderNum;
	//int   m_wNum, m_wDegree, m_wRenderNum;
	int m_uRenderNum, m_vRenderNum, m_wRenderNum;
	bool  m_bHeterogeneous;//������
	//varray<Vec4>			 m_CtrlPts;
	varray<VolumeVertex> m_vAllCtrlPts;  //���Ƶ�
	varray<VolumeVertex> m_vVolumePts;   //������ʾ�Ľڵ㡣
	varray<int>  m_CtrlPtsIDinKKKMatrix;   //ÿ�����Ƶ����ܸپ����е�λ�á� ����ǵ�Ƭ��Ĭ�Ͼ���0~m_uNum*m_vNum*m_wNum-1;����Ƕ�Ƭ����Ҫ���±�š�
	bool  m_bHasInterpolated;//������
	//varray<double> m_uKnots;
	//varray<double> m_vKnots;
	//varray<double> m_wKnots;
	varray<NurbsSurface> m_6boundarySurface;  //����ǰ���Ƕ�Ӧm_6PatchIdx����Ӧ����ϵ�surface,�洢������ת�������
	bool  m_bHasIGASolution;
	//float m_maxJocbian, m_minJocbian;
	//int   m_minorJocbianNum;
	//double m_isovalAdd;
	//��������IGA3D
	varray<point4d>m_CPdu;
	varray<double>m_dis;
	varray<double>m_stress;
	varray<varray<varray<Vec3>>>m_show3Dcpts;//����3dIGA��ʾ�Ŀ��Ƶ�
	varray<varray<varray<double>>>m_show3DDisp;//����3dIGA��ʾ��λ��
	varray<varray<varray<Vec3>>>m_show3DEdge;//����3dIGA��ʾ�߽���
	varray<varray<varray<Vec3>>>m_show3DcptsOfStress;//����3dIGA��ʾ�Ŀ��Ƶ�
	varray<varray<varray<double>>>m_show3DStrainOrStress;//����3dIGA��ʾ��Ӧ��
};