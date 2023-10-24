#pragma once
//#include <Geom2d_Line.hxx>
//#include <TopoDS_Edge.hxx>
#include "iMSO_NURBS.h"
#include <Geom2d_Point.hxx>
#include <Geom2d_Circle.hxx>
#include <BRepLib_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>



//������
class FeaturePoint
{
	
};

//������
class FeatureCurve
{
public:
	//������ת��anglesΪ�Ƕȣ�mode123�ֱ�Ϊ��xyz����
	/*void Rolate(NurbsLine &NL, double angles, int mode);*/
public:
	NurbsLine nl;
	//BRepLib_MakeEdge ol;
};


//���㹹���߶�

class FeatureCurve_Line
{
public:
	FeatureCurve_Line(double x1, double y1, double x2, double y2);

public:
	FeatureCurve FC;
	TopoDS_Edge ol1;
};


//�߽Ǿ��Σ����������Դ���Լ����
class FeatureCurve_Rectangle
{
public:
	enum InitType
	{
		Edges,
		Wire
	};
	FeatureCurve_Rectangle(double x1, double y1, double x2, double y2, InitType type = Edges);

public:
	varray<FeatureCurve> FC;
	//������
	varray<TopoDS_Edge> aEdges;
	//�߿�
	TopoDS_Wire aWire;
private:
	void InitEdges(double x1, double y1, double x2, double y2);
	void InitWire(double x1, double y1, double x2, double y2);

};

//class FeatureCurve_Circle_arc
//{
//public:
//	FeatureCurve_Circle_arc();
//private:
//	FeatureCurve FC;
//};

class FeatureCurve_Circle
{
public:
	FeatureCurve_Circle();

public:
	varray<FeatureCurve> FC;
};


