#include "selectobj.h"
#include <utility>
#include <Prs3d_LineAspect.hxx>

SelectObj::SelectObj()
{

}

SelectObj::~SelectObj()
{

}

Handle(AIS_Shape) SelectObj::getSelectedShape()
{
	return m_selectedShape;
}

bool SelectObj::SelectModel(Handle(AIS_InteractiveContext)& myContext,
	Handle(V3d_View)& myView,
	Standard_Integer x,
	Standard_Integer y,
	Standard_Integer mode)
{
	myContext->Activate(mode);
	myContext->MoveTo(x, y, myView, true);
	auto r = myContext->Select(true);
	//若选中了一个对象，返回true
	if (r == AIS_SOP_OneSelected) {
		//m_selectedShape->SetShape(myContext->SelectedShape());
		return true;
	}

	return false;

}



