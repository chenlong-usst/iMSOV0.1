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



//特征点
class FeaturePoint
{
	
};

//特征线
class FeatureCurve
{
public:
	//曲线旋转，angles为角度，mode123分别为绕xyz方向
	/*void Rolate(NurbsLine &NL, double angles, int mode);*/
public:
	NurbsLine nl;
	//BRepLib_MakeEdge ol;
};


//两点构造线段

class FeatureCurve_Line
{
public:
	FeatureCurve_Line(double x1, double y1, double x2, double y2);

public:
	FeatureCurve FC;
	TopoDS_Edge ol1;
};


//边角矩形（包含矩形自带的约束）
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
	//四条边
	varray<TopoDS_Edge> aEdges;
	//线框
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


