#include "drawobj.h"
#include <cmath>
#include <QDebug>
#include <gp_Trsf.hxx>
//#include "Conv/RWGeometry.h"

DrawObj::DrawObj() :x0(0), y0(0), x1(0), y1(0), degree(2)
{
	TColStd_Array1OfReal knot(0, 1);
	knot.SetValue(0, 0);
	knot.SetValue(1, 1);
	TColStd_Array1OfInteger mulit(0, 1);
	for (int i = 0; i <= 1; i++) {
		mulit.SetValue(i, 3);
	}
	TColStd_Array1OfReal wights(0, 2);
	for (int i = 0; i <= 2; i++) {
		wights.SetValue(i, 1);
	}

	this->wights.Move(wights);
	this->knots.Move(knot);
	this->mulits.Move(mulit);

}


Handle(AIS_Shape) DrawObj::drawRectangle()
{
	TColgp_Array1OfPnt cpts(0, 2);
	Standard_Real a0 = static_cast<Standard_Real>(x0);
	Standard_Real a1 = static_cast<Standard_Real>(x1);
	Standard_Real b0 = static_cast<Standard_Real>(y0);
	Standard_Real b1 = static_cast<Standard_Real>(y1);
	gp_Pnt p0 = { a0,b0,0 };
	gp_Pnt p1 = { (a0 + a1) / 2, b0, 0 };
	gp_Pnt p2 = { a1, b0, 0 };
	cpts.SetValue(0, p0);
	cpts.SetValue(1, p1);
	cpts.SetValue(2, p2);
	Handle(Geom_BSplineCurve) line0 = new Geom_BSplineCurve(cpts, knots, mulits, degree);
	p0 = { a0,b0,0 };
	p1 = { a0, (b0 + b1) / 2, 0 };
	p2 = { a0, b1, 0 };
	cpts.SetValue(0, p0);
	cpts.SetValue(1, p1);
	cpts.SetValue(2, p2);

	Handle(Geom_BSplineCurve) line1 = new Geom_BSplineCurve(cpts, knots, mulits, degree);

	p0 = { a0,b1,0 };
	p1 = { (a0 + a1) / 2, b1, 0 };
	p2 = { a1, b1, 0 };
	cpts.SetValue(0, p0);
	cpts.SetValue(1, p1);
	cpts.SetValue(2, p2);
	Handle(Geom_BSplineCurve) line2 = new Geom_BSplineCurve(cpts, knots, mulits, degree);

	p0 = { a1,b0,0 };
	p1 = { a1, (b0 + b1) / 2, 0 };
	p2 = { a1, b1, 0 };
	cpts.SetValue(0, p0);
	cpts.SetValue(1, p1);
	cpts.SetValue(2, p2);
	Handle(Geom_BSplineCurve) line3 = new Geom_BSplineCurve(cpts, knots, mulits, degree);

	TopoDS_Edge l1 = BRepBuilderAPI_MakeEdge(line0);
	TopoDS_Edge l2 = BRepBuilderAPI_MakeEdge(line1);
	TopoDS_Edge l3 = BRepBuilderAPI_MakeEdge(line2);
	TopoDS_Edge l4 = BRepBuilderAPI_MakeEdge(line3);

	TopoDS_Wire r = BRepBuilderAPI_MakeWire(l1, l2, l3, l4);
	Handle(AIS_Shape) rectangle = new AIS_Shape(r);


	//point3d p00 = { a0,b0,0 };
	//point3d p01 = { a1,b0,0 };
	//point3d p02 = { a0,b1,0 };
	//point3d p03 = { a1,b1,0 };

	//NurbsLine line01;
	//line01.CreatLineWithTwoPoints(p00, p01);

	//NurbsLine line02;
	//line02.CreatLineWithTwoPoints(p01, p03);

	//NurbsLine line03;
	//line03.CreatLineWithTwoPoints(p00, p02);

	//NurbsLine line04;
	//line04.CreatLineWithTwoPoints(p02, p03);

	//varray<NurbsLine> v1;
	//v1.push_back(line01);
	//v1.push_back(line02);
	//v1.push_back(line03);
	//v1.push_back(line04);

	//RWGeometric rwg;
	//rwg.WriteNurbsLine("C://oldiMSO//Models_IGA//QtGuiApplication1/data//Sketch.txt", v1);

	return rectangle;

}


Handle(AIS_Shape) DrawObj::drawCurve()
{
	TColgp_Array1OfPnt cpts(0, 2);
	Standard_Real a0 = static_cast<Standard_Real>(x0);
	Standard_Real a1 = static_cast<Standard_Real>(x1);
	Standard_Real b0 = static_cast<Standard_Real>(y0);
	Standard_Real b1 = static_cast<Standard_Real>(y1);

	gp_Pnt p0 = { a0,b0,0 };
	gp_Pnt p1 = { (a0 + a1) / 2, (b0 + b1) / 2, 0 };
	gp_Pnt p2 = { a1, b1, 0 };
	cpts.SetValue(0, p0);
	cpts.SetValue(1, p1);
	cpts.SetValue(2, p2);
	Handle(Geom_BSplineCurve) line0 = new Geom_BSplineCurve(cpts, knots, mulits, degree);
	p0 = { a0,b0,0 };
	p1 = { a0, (b0 + b1) / 2, 0 };
	p2 = { a0, b1, 0 };
	cpts.SetValue(0, p0);
	cpts.SetValue(1, p1);
	cpts.SetValue(2, p2);


	TopoDS_Edge l1 = BRepBuilderAPI_MakeEdge(line0);
	
	TopoDS_Wire r = BRepBuilderAPI_MakeWire(l1);

	Handle(AIS_Shape) line = new AIS_Shape(r);
	   
	return line;

}


Handle(AIS_Shape) DrawObj::drawCircle()
{

	TColgp_Array1OfPnt cpts(0, 2);
	Standard_Real a0 = static_cast<Standard_Real>(x0);
	Standard_Real a1 = static_cast<Standard_Real>(x1);
	Standard_Real b0 = static_cast<Standard_Real>(y0);
	Standard_Real b1 = static_cast<Standard_Real>(y1);
	gp_Pnt c = { a0,b0,0 };
	gp_Pnt p0 = { a1,b1,0 };
	Standard_Real r = p0.Distance(c);
	Standard_Real alph = M_PI / 4;
	gp_Pnt p1 = { a0 - sqrt(2)*r / 2,b0 - sqrt(2)*r / 2,0 };
	gp_Pnt p2 = { a0, b0 - sqrt(2)*r, 0 };
	gp_Pnt p3 = { a0 + sqrt(2)*r / 2,b0 - sqrt(2)*r / 2,0 };
	wights.SetValue(1, cos(alph));
	cpts.SetValue(0, p1);
	cpts.SetValue(1, p2);
	cpts.SetValue(2, p3);
	Handle(Geom_BSplineCurve) line0 = new Geom_BSplineCurve(cpts, wights, knots, mulits, degree);
	gp_Ax1 xc;
	xc.SetLocation(c);
	TopoDS_Edge l1 = BRepBuilderAPI_MakeEdge(line0);
	p1 = { a0 + sqrt(2)*r / 2,b0 - sqrt(2)*r / 2,0 };
	p2 = { a0 + sqrt(2)*r, b0, 0 };
	p3 = { a0 + sqrt(2)*r / 2,b0 + sqrt(2)*r / 2,0 };
	cpts.SetValue(0, p1);
	cpts.SetValue(1, p2);
	cpts.SetValue(2, p3);
	Handle(Geom_BSplineCurve) line1 = new Geom_BSplineCurve(cpts, wights, knots, mulits, degree);
	TopoDS_Edge l2 = BRepBuilderAPI_MakeEdge(line1);
	p1 = { a0 + sqrt(2)*r / 2,b0 + sqrt(2)*r / 2,0 };
	p2 = { a0, b0 + sqrt(2)*r, 0 };
	p3 = { a0 - sqrt(2)*r / 2,b0 + sqrt(2)*r / 2,0 };
	cpts.SetValue(0, p1);
	cpts.SetValue(1, p2);
	cpts.SetValue(2, p3);
	Handle(Geom_BSplineCurve) line2 = new Geom_BSplineCurve(cpts, wights, knots, mulits, degree);
	TopoDS_Edge l3 = BRepBuilderAPI_MakeEdge(line2);

	p1 = { a0 - sqrt(2)*r / 2,b0 + sqrt(2)*r / 2,0 };
	p2 = { a0 - sqrt(2)*r, b0, 0 };
	p3 = { a0 - sqrt(2)*r / 2,b0 - sqrt(2)*r / 2,0 };
	cpts.SetValue(0, p1);
	cpts.SetValue(1, p2);
	cpts.SetValue(2, p3);
	Handle(Geom_BSplineCurve) line3 = new Geom_BSplineCurve(cpts, wights, knots, mulits, degree);
	TopoDS_Edge l4 = BRepBuilderAPI_MakeEdge(line3);

	TopoDS_Wire cir = BRepBuilderAPI_MakeWire(l1, l2, l3, l4);
	Handle(AIS_Shape) circle = new AIS_Shape(cir);

	return circle;
}


void DrawObj::setInitXY(Standard_Integer x, Standard_Integer y)
{
	x0 = x;
	y0 = y;
}

void DrawObj::setMoveXY(Standard_Integer x, Standard_Integer y)
{
	x1 = x;
	y1 = y;
}

void DrawObj::setDegree(Standard_Integer degree)
{
	this->degree = degree;
}

void DrawObj::setKnots(TColStd_Array1OfReal knots, TColStd_Array1OfInteger mulits)
{
	this->knots = knots;
	this->mulits = mulits;
}
