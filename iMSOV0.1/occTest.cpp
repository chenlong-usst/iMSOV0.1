#include "occTest.h"

occTest::occTest()
{
}

Handle(AIS_Shape) occTest::myDraw(int drawType)
{
	switch (drawType)
	{
	case 0:
		return drawVertex();
		break;
	case 1:
		return drawEdge();
		break;
	case 2:
		return drawFace();
		break;
	case 3:
		return drawWire(); 
		break;
	case 4:
		return drawShell();
		break;
	case 5:
		return drawSolid();
		break;
	default:
		break;
	}

	return Handle(AIS_Shape)();
}

Handle(AIS_Shape) occTest::FeatureLineDraw(double x1, double y1, double x2, double y2)
{
	FeatureCurve_Line fl(x1, y1, x2, y2);
	Handle(AIS_Shape) myShape = new AIS_Shape(fl.ol1);
	cout << "FetureCurve_Line was created in red" << endl;
	return myShape;
}



varray<Handle(AIS_Shape)> occTest::FeatureRecDraw(double x1, double y1, double x2, double y2, FeatureCurve_Rectangle::InitType type)
{
	varray<Handle(AIS_Shape)> v;
	FeatureCurve_Rectangle fr(x1, y1, x2, y2, type);
	if (type == FeatureCurve_Rectangle::InitType::Edges)
	{
		for (int i = 0; i < fr.aEdges.size(); i++)
		{
			Handle(AIS_Shape) myShape = new AIS_Shape(fr.aEdges[i]);
			v.push_back(myShape);
		}
	}
	else if (type == FeatureCurve_Rectangle::InitType::Wire) 
	{
		Handle(AIS_Shape) myShape = new AIS_Shape(fr.aWire);
		v.push_back(myShape);
	}

	cout << "FetureCurve_Rectangle was created in red" << endl;


	return v;
}


Handle(AIS_Shape) occTest::drawVertex()
{
	// Make a vertex from a 3D point.
	gp_Pnt aPnt(0.0, 0.0, 10.0);
	TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(aPnt);
	cout << "TopoDS_Vertex was created at [ "
		<< aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z()
		<< " ]" << endl;
	Handle(AIS_Shape) aAisVertex = new AIS_Shape(aVertex);
	Handle(AIS_TextLabel) anAisLabel = new AIS_TextLabel();
	Standard_SStream aSS;
	aSS << "TopoDS_Vertex [" << aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z() << "]" << std::endl;
	anAisLabel->SetText(aSS.str().c_str());
	anAisLabel->SetPosition(aPnt);
	myText = anAisLabel;
	return aAisVertex;
}

Handle(AIS_Shape) occTest::drawEdge()
{
	// Make an edge from a 3D curve (BSpline).
	// Define points.
	gp_Pnt aPole1(0.0, 0.0, 10.0);
	gp_Pnt aPole2(5.0, 5.0, 5.0);
	gp_Pnt aPole3(10.0, 10.0, 15.0);
	gp_Pnt aPole4(15.0, 5.0, 20.0);
	// Add points to the curve poles array.
	TColgp_Array1OfPnt aPoles(1, 4);
	aPoles.SetValue(1, aPole1);
	aPoles.SetValue(2, aPole2);
	aPoles.SetValue(3, aPole3);
	aPoles.SetValue(4, aPole4);
	// Make a BSpline curve from the points array
	Handle(Geom_BSplineCurve) aBSplineCurve = GeomAPI_PointsToBSpline(aPoles).Curve();
	// Make an edge between two point on the BSpline curve.
	gp_Pnt aPntOnCurve1, aPntOnCurve2;
	aBSplineCurve->D0(0.75 * aBSplineCurve->FirstParameter()
		+ 0.25 * aBSplineCurve->LastParameter(),
		aPntOnCurve1);
	aBSplineCurve->D0(0.25 * aBSplineCurve->FirstParameter()
		+ 0.75 * aBSplineCurve->LastParameter(),
		aPntOnCurve2);
	TopoDS_Edge anEdgeBSpline = BRepBuilderAPI_MakeEdge(aBSplineCurve, aPntOnCurve1, aPntOnCurve2);
	cout << "TopoDS_Edge on the BSpline curve" << endl
		<< "between [ "
		<< aPntOnCurve1.X() << ", " << aPntOnCurve1.Y() << ", " << aPntOnCurve1.Z()
		<< " ] and [ "
		<< aPntOnCurve2.X() << ", " << aPntOnCurve2.Y() << ", " << aPntOnCurve2.Z()
		<< " ]" << std::endl
		<< "was created in green" << endl;
	Handle(AIS_ColoredShape) anAisEdgeBSpline = new AIS_ColoredShape(anEdgeBSpline);
	anAisEdgeBSpline->SetColor(Quantity_Color(Quantity_NOC_GREEN));

	Handle(AIS_TextLabel) anAisEdgeBSplineLabel = new AIS_TextLabel();
	anAisEdgeBSplineLabel->SetText("BSpline edge");
	anAisEdgeBSplineLabel->SetPosition(aPole3);
	anAisEdgeBSplineLabel->SetColor(Quantity_Color(Quantity_NOC_GREEN));
	myText = anAisEdgeBSplineLabel;

	return anAisEdgeBSpline;
}

Handle(AIS_Shape) occTest::drawFace()
{
	// Make a face from a BSpline surface.
  // Define a 4x4 grid of points for BSpline surface.
	TColgp_Array2OfPnt aPoints(1, 4, 1, 4);
	for (Standard_Integer i = 1; i <= 4; ++i)
	{
		gp_Pnt aPnt;
		aPnt.SetX(5.0 * i);
		for (Standard_Integer j = 1; j <= 4; ++j)
		{
			aPnt.SetY(5.0 * j);
			if (1 < i && i < 4 && 1 < j && j < 4)
			{
				aPnt.SetZ(15.0);
			}
			else
			{
				aPnt.SetZ(10.0);
			}
			aPoints.SetValue(i, j, aPnt);
		}
	}
	// Make a BSpline surface from the points array.
	Handle(Geom_BSplineSurface) aBSplineSurf = GeomAPI_PointsToBSplineSurface(aPoints).Surface();
	Standard_Real aU1, aU2, aV1, aV2;
	aBSplineSurf->Bounds(aU1, aU2, aV1, aV2);
	TopoDS_Face aFaceBSpline = BRepBuilderAPI_MakeFace(aBSplineSurf, aU1, aU2, aV1, aV2, Precision::Confusion());
	cout << "TopoDS_Face on the BSpline surface was created in green" << endl << endl;

	Handle(AIS_ColoredShape) anAisFaceBSpline = new AIS_ColoredShape(aFaceBSpline);
	anAisFaceBSpline->SetColor(Quantity_Color(Quantity_NOC_GREEN));
	Handle(AIS_TextLabel) anAisFaceBSplineLabel = new AIS_TextLabel();
	anAisFaceBSplineLabel->SetText("BSpline face");
	anAisFaceBSplineLabel->SetPosition(aPoints(4, 4));
	anAisFaceBSplineLabel->SetColor(Quantity_Color(Quantity_NOC_GREEN));
	myText = anAisFaceBSplineLabel;

	return anAisFaceBSpline;
}

Handle(AIS_Shape) occTest::drawWire()
{
	// Make a wire from edges created on a set of points.
  // Add points to the curve poles array.
	TColgp_Array1OfPnt aPoints(1, 4);
	aPoints.SetValue(1, gp_Pnt(0.0, 0.0, 0.0));
	aPoints.SetValue(2, gp_Pnt(20.0, 0.0, 0.0));
	aPoints.SetValue(3, gp_Pnt(20.0, 10.0, 0.0));
	aPoints.SetValue(4, gp_Pnt(0.0, 10.0, 0.0));
	// A wire maker contains an empty wire.
	BRepBuilderAPI_MakeWire aMakeWire;
	for (Standard_Integer i = 1; i <= 4; ++i)
	{
		Standard_Integer i1 = i;
		Standard_Integer i2 = i1 < 4 ? i1 + 1 : 1;
		const gp_Pnt& aPnt1 = aPoints.Value(i1);
		const gp_Pnt& aPnt2 = aPoints.Value(i2);
		TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(aPnt1, aPnt2);
		// Add an edge to the wire under construction.
		// The edge must be connectible to the wire under construction, and,
		// unless it is the first edge of the wire, must satisfy the following
		// condition: one of its vertices must be geometrically coincident
		// with one of the vertices of the wire (provided that the highest
		// tolerance factor is assigned to the two vertices).
		// It could also be the same vertex.
		// Warning
		// If the edge is not connectible to the wire under construction it is not added.
		// The function IsDone will return false and the function
		// Wire will raise an error, until a new connectible edge is added.
		aMakeWire.Add(anEdge);
		Standard_ASSERT_VOID(aMakeWire.IsDone(), "Added edge isn't connectible!");
	}
	// Retrieve a constructed wire.
	TopoDS_Wire aWire = aMakeWire.Wire();
	cout << "TopoDS_Wire was created. Vertices :" << endl;
	// Retrieve wire vertices. 4 vertices are expected, because of
	// edges connecting during wire constructing.
	TopTools_IndexedMapOfShape aVertices;
	TopExp::MapShapes(aWire, TopAbs_VERTEX, aVertices);
	for (TopTools_IndexedMapOfShape::Iterator anIt(aVertices); anIt.More(); anIt.Next())
	{
		TopoDS_Vertex aVertex = TopoDS::Vertex(anIt.Value());
		gp_Pnt aPnt = BRep_Tool::Pnt(aVertex);
		cout << "[ " << aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z() << " ]" << endl;
		Handle(AIS_Shape) anAisVertex = new AIS_Shape(aVertex);
		//myObject3d.Append(anAisVertex);
	}

	Handle(AIS_Shape) anAisWire = new AIS_Shape(aWire);
	//myObject3d.Append(anAisWire);

	return anAisWire;
}

Handle(AIS_Shape) occTest::drawShell()
{
	// Make a shell from a cylinder with R = 5 and directed along Z axis
	gp_Cylinder aCyl(gp::XOY(), 5.0);
	Handle(Geom_Surface) aCylSurf = new Geom_CylindricalSurface(aCyl);
	TopoDS_Shell aCylShell = BRepBuilderAPI_MakeShell(aCylSurf, 0.0, 2.0 * M_PI, -10.0, +10.0);
	cout << "TopoDS_Shell on the cylinder R = " << aCyl.Radius() << endl
		<< "with axis [ "
		<< aCyl.Position().Direction().X() << ", "
		<< aCyl.Position().Direction().Y() << ", "
		<< aCyl.Position().Direction().Z() << " ]" << std::endl
		<< "limited in length [-10 ... +10] was created" << std::endl;

	Handle(AIS_Shape) anAisShell = new AIS_Shape(aCylShell);
	return anAisShell;
}

Handle(AIS_Shape) occTest::drawSolid()
{
	// Make a torus from a shell.
	gp_Torus aTorus(gp::XOY(), 20.0, 7.5);

	Handle(Geom_Surface) aTorusSurf = new Geom_ToroidalSurface(aTorus);
	TopoDS_Shell aTorusShell = BRepBuilderAPI_MakeShell(aTorusSurf, 0.0, 2.0 * M_PI, 0.0, 2.0 * M_PI);
	// Make a solid on the torus shell.
	TopoDS_Solid aTorusSolid = BRepBuilderAPI_MakeSolid(aTorusShell);
	cout << "TopoDS_Solid on the torus with" << endl
		<< "R major = " << aTorus.MajorRadius() << endl
		<< "R minor = " << aTorus.MinorRadius() << endl
		<< "was created" << endl;

	Handle(AIS_Shape) anAisSolid = new AIS_Shape(aTorusSolid);
	return anAisSolid;
}
