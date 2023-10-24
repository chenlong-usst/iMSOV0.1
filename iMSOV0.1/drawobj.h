#pragma once
#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

#include <gp_Lin2d.hxx>

#include <GCE2d_MakeSegment.hxx>

#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <BRepLib.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

#include <AIS_Shape.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <vector>
#include <V3d_View.hxx>


class DrawObj
{
public:
	DrawObj();
	//绘制一个矩形
	//返回值为AIS_Shape
	Handle(AIS_Shape) drawRectangle();
	//绘制一个圆
	//返回值为AIS_Shape
	Handle(AIS_Shape) drawCircle();
	//绘制一个曲线
	Handle(AIS_Shape) drawCurve();
	Handle(AIS_Shape) prevShape;
private:
	Standard_Integer x0;    //鼠标初始点
	Standard_Integer y0;    //
	Standard_Integer x1;    //鼠标移动点
	Standard_Integer y1;    //

	TColStd_Array1OfReal knots;     //节点矢量
	TColStd_Array1OfInteger mulits; //重复度
	TColStd_Array1OfReal wights;    //权值
	Standard_Integer degree;        //次数

public:
	void setInitXY(Standard_Integer x, Standard_Integer y);
	void setMoveXY(Standard_Integer x, Standard_Integer y);
	void setDegree(Standard_Integer degree);
	void setKnots(TColStd_Array1OfReal knots, TColStd_Array1OfInteger mulits);

public:
	//绘制模式
	enum DrawMode {
		DrawMode_Rectangle, //绘制矩形
		DrawMode_Circle     //绘制圆
	};


};