#pragma once
#include <V3d_Viewer.hxx>
#include <AIS_Shape.hxx>
#include <AIS_InteractiveContext.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>



class SelectObj
{

public:

	SelectObj();
	~SelectObj();

public:
	//Activate参数：0为整体选择，1为选择点，2为选择线，3为选择线框，4为选择面
	enum SelectMode {
		SelectMode_Whole,
		SelectMode_Point,
		SelectMode_Line,
		SelectMode_Boundary,
		SelectMode_Face
	};
private:
	TopoDS_Face selectedFace;
	Standard_Integer m_currentmode;
	Handle(AIS_Shape) m_selectedShape;


public:
	bool SelectModel(Handle(AIS_InteractiveContext)& myContext,
		Handle(V3d_View)& myView,
		Standard_Integer x,
		Standard_Integer y,
		Standard_Integer mode);
	Handle(AIS_Shape) getSelectedShape();

};
