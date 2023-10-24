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
	//����һ������
	//����ֵΪAIS_Shape
	Handle(AIS_Shape) drawRectangle();
	//����һ��Բ
	//����ֵΪAIS_Shape
	Handle(AIS_Shape) drawCircle();
	//����һ������
	Handle(AIS_Shape) drawCurve();
	Handle(AIS_Shape) prevShape;
private:
	Standard_Integer x0;    //����ʼ��
	Standard_Integer y0;    //
	Standard_Integer x1;    //����ƶ���
	Standard_Integer y1;    //

	TColStd_Array1OfReal knots;     //�ڵ�ʸ��
	TColStd_Array1OfInteger mulits; //�ظ���
	TColStd_Array1OfReal wights;    //Ȩֵ
	Standard_Integer degree;        //����

public:
	void setInitXY(Standard_Integer x, Standard_Integer y);
	void setMoveXY(Standard_Integer x, Standard_Integer y);
	void setDegree(Standard_Integer degree);
	void setKnots(TColStd_Array1OfReal knots, TColStd_Array1OfInteger mulits);

public:
	//����ģʽ
	enum DrawMode {
		DrawMode_Rectangle, //���ƾ���
		DrawMode_Circle     //����Բ
	};


};