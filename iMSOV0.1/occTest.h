#pragma once
#include <iostream>
using namespace std;

#include <gp_Cylinder.hxx>
#include <gp_Torus.hxx>

#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include <AIS_Shape.hxx>
#include <AIS_ColoredShape.hxx>
#include <AIS_TextLabel.hxx>

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>

#include <TopExp.hxx>
#include <BRep_Tool.hxx>

#include "FeatureCurve.h"

class occTest
{
public:

	occTest();
	Handle(AIS_Shape)myDraw(int drawType);
	Handle(AIS_Shape)FeatureLineDraw(double x1, double y1, double x2, double y2);
	varray<Handle(AIS_Shape)>FeatureRecDraw(double x1, double y1, double x2, double y2, FeatureCurve_Rectangle::InitType type = FeatureCurve_Rectangle::InitType::Edges);


	Handle(AIS_TextLabel)myText;

	enum drawType {
		myDraw_Vertex,
		myDraw_Edge,
		myDraw_Face,
		myDraw_Wire,
		myDraw_Shell,
		myDraw_Solid
	};

private:
	Handle(AIS_Shape) drawVertex();
	Handle(AIS_Shape) drawEdge();
	Handle(AIS_Shape) drawFace();
	Handle(AIS_Shape) drawWire();
	Handle(AIS_Shape) drawShell();
	Handle(AIS_Shape) drawSolid();
};

