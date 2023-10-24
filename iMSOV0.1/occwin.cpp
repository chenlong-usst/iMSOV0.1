#include "occwin.h"
#include "ReadModelDialog.h"

#include <TColgp_Array1OfPnt2d.hxx>
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
#include <QPalette>
#include <fstream>
#include <STEPCAFControl_Writer.hxx>

OccWin::OccWin(QWidget *parent)
	: QDialog(parent)
{
	//ui.setupUi(this);
	occView = new OccView(this);
	menuBar = new QMenuBar(this);
	layout = new QVBoxLayout(this);
	myAction = new QAction(this);
	toolBar = new QToolBar(this);
	this->setGeometry(50, 50, 1600, 900);
	menuBar->setGeometry(0, 0, 800, 40);
	occView->setGeometry(0, 0, 1600, 900);
	toolBar->setGeometry(0, 40, 800, 40);
	toolBar->setStyleSheet("background-color:rgb(255,255,255);");
	toolBar->hide();
	drawObj = new DrawObj();
	myTest = new occTest();

	//加载菜单栏
	InitMenu();

	connect(menuBar, SIGNAL(triggered(QAction*)), this, SLOT(trigerMenu(QAction*)));
	connect(toolBar, SIGNAL(actionTriggered(QAction*)), this, SLOT(trigerToolBar(QAction*)));

}

OccWin::~OccWin()
{
	delete occView;
	delete menuBar;
	delete layout;
	delete myAction;
	delete selectObj;
	delete toolBar;
	delete drawObj;
	delete myTest;
}


void OccWin::Convert_Spline_Occtype(OccCurve & spline, Handle(Geom_BSplineCurve)& occline)
{

	TColStd_Array1OfReal occknot;

	TColStd_Array1OfInteger occmulit;

	Convert_knot_occtype(spline.m_Knots, occknot, occmulit);

	TColStd_Array1OfReal occwights(0, spline.m_CtrlPts.size() - 1);
	TColgp_Array1OfPnt cpts(0, spline.m_CtrlPts.size() - 1);
	for (int i = 0; i < spline.m_CtrlPts.size(); i++) {
		gp_Pnt pt = { spline.m_CtrlPts[i].x,spline.m_CtrlPts[i].y,spline.m_CtrlPts[i].z };
		cpts.SetValue(i, pt);
		occwights.SetValue(i, spline.m_CtrlPts[i].w);
	}

	occline = new Geom_BSplineCurve(cpts, occwights, occknot, occmulit, spline.m_Degree);

}

void OccWin::Convert_SplineSurface_Occtype(OccSurface & splinesurface, Handle(Geom_BSplineSurface)& occsurface)
{
	TColStd_Array1OfReal occuknot;
	TColStd_Array1OfInteger occumulit;
	TColStd_Array1OfReal occvknot;
	TColStd_Array1OfInteger occvmulit;
	Convert_knot_occtype(splinesurface.m_uKnots, occuknot, occumulit);
	Convert_knot_occtype(splinesurface.m_vKnots, occvknot, occvmulit);
	TColStd_Array2OfReal occwights(0, splinesurface.m_uNum - 1, 0, splinesurface.m_vNum - 1);
	TColgp_Array2OfPnt cpts(0, splinesurface.m_uNum - 1, 0, splinesurface.m_vNum - 1);
	for (int i = 0; i < splinesurface.m_vNum; i++) {
		for (int j = 0; j < splinesurface.m_uNum; j++) {
			int idx = i * splinesurface.m_vNum + j;
			gp_Pnt pt = { splinesurface.m_CtrlPts[idx].x,splinesurface.m_CtrlPts[idx].y,splinesurface.m_CtrlPts[idx].z };
			cpts.SetValue(j, i, pt);
			occwights.SetValue(j, i, splinesurface.m_CtrlPts[idx].w);
		}
	}

	occsurface = new Geom_BSplineSurface(cpts, occwights, occuknot, occvknot,
		occumulit, occvmulit, splinesurface.m_uDegree, splinesurface.m_vDegree, false, false);
}

void OccWin::Convert_SplineVolume_Occtype(OccVol & splinevolume, vector<Handle(Geom_BSplineSurface)>& occsurfaces)
{
	auto surfaces = GetSurfaces(splinevolume);
	for (int i = 0; i < surfaces.size(); i++) {
		Handle(Geom_BSplineSurface) sf;
		Convert_SplineSurface_Occtype(surfaces[i], sf);
		occsurfaces.push_back(sf);
	}
}


void OccWin::Convert_knot_occtype(vector<double>& knot, TColStd_Array1OfReal& occknot, TColStd_Array1OfInteger& occmulits)
{
	vector<double> newknot;
	newknot.push_back(knot[0]);
	vector<int> mulits;
	int mulit = 1;
	for (int i = 0; i < knot.size(); i++) {
		if (i == 0) continue;
		if (i > 0 && abs(knot[i] - knot[i - 1]) > 0.001) {
			mulits.push_back(mulit);
			newknot.push_back(knot[i]);
			mulit = 1;
		}
		else {
			mulit++;
		}
	}
	mulits.push_back(mulit);

	TColStd_Array1OfReal occknot1(0, newknot.size() - 1);
	for (int i = 0; i < newknot.size(); i++) {
		occknot1.SetValue(i, newknot[i]);
	}
	TColStd_Array1OfInteger occmulit1(0, mulits.size() - 1);
	for (int i = 0; i < mulits.size(); i++) {
		occmulit1.SetValue(i, mulits[i]);
	}
	occknot = std::move(occknot1);
	occmulits = std::move(occmulit1);

}

vector<OccSurface> OccWin::GetSurfaces(const OccVol& NVS)			//将所有的面提取出来
{
	OccSurface U0;		//0
	OccSurface U1;		//1
	OccSurface V0;		//2
	OccSurface V1;		//3
	OccSurface W0;		//4
	OccSurface W1;		//5
	vector<vector<OccSurface>> Res;
	vector<OccSurface> Sum;
	int m = NVS.m_vNum;
	int n = NVS.m_uNum;
	int p = NVS.m_wNum;
	OccVol Vtemp = NVS;
	int idx = 0;
	for (int x = 0; x < m*n; x++)
	{
		W0.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
	}
	W0.SetSurface(Vtemp.m_uDegree, Vtemp.m_vDegree, Vtemp.m_uNum, Vtemp.m_vNum, Vtemp.m_uKnots, Vtemp.m_vKnots);
	for (int x = m * n*(p - 1); x < m*n*p; x++)
	{
		W1.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
	}
	W1.SetSurface(Vtemp.m_uDegree, Vtemp.m_vDegree, Vtemp.m_uNum, Vtemp.m_vNum, Vtemp.m_uKnots, Vtemp.m_vKnots);

	for (int y = 0; y <= (p - 1)*m*n; y += m * n) {
		for (int x = y; x <= y + m * n - m; x += m) {
			V0.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
		}
	}
	V0.SetSurface(Vtemp.m_wDegree, Vtemp.m_uDegree, Vtemp.m_wNum, Vtemp.m_uNum, Vtemp.m_wKnots, Vtemp.m_uKnots);
	for (int y = m - 1; y <= m * n*(p - 1) + m - 1; y += m * n)
		for (int x = y; x <= y + m * n - m; x += m) {
			V1.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
		}
	V1.SetSurface(Vtemp.m_wDegree, Vtemp.m_uDegree, Vtemp.m_wNum, Vtemp.m_uNum, Vtemp.m_wKnots, Vtemp.m_uKnots);

	for (int y = 0; y < p; y++) {
		for (int x = y * m*n; x < y*m*n + m; x++)
		{
			U0.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
		}
	}
	U0.SetSurface(Vtemp.m_wDegree, Vtemp.m_vDegree, Vtemp.m_wNum, Vtemp.m_vNum, Vtemp.m_wKnots, Vtemp.m_vKnots);

	for (int y = 0; y < p; y++) {
		for (int x = y * m*n + m * n - m; x < y*m*n + m * n; x++)
		{
			U1.m_CtrlPts.push_back(NVS.m_CtrlPts[x]);
		}
	}
	U1.SetSurface(Vtemp.m_wDegree, Vtemp.m_vDegree, Vtemp.m_wNum, Vtemp.m_vNum, Vtemp.m_wKnots, Vtemp.m_vKnots);


	Sum.push_back(V1);
	V1.m_CtrlPts.clear();
	Sum.push_back(V0);
	V0.m_CtrlPts.clear();

	Sum.push_back(U0);
	U0.m_CtrlPts.clear();

	Sum.push_back(U1);
	U1.m_CtrlPts.clear();

	Sum.push_back(W0);
	W0.m_CtrlPts.clear();

	Sum.push_back(W1);
	W1.m_CtrlPts.clear();

	return Sum;
}
//void OccWin::Convert_Spline_Occtype(NurbsLine & spline, Handle(Geom_BSplineCurve)& occline)
//{
//
//	TColStd_Array1OfReal occknot;
//
//	TColStd_Array1OfInteger occmulit;
//	vector<double> v = *spline._u_Knots;
//
//	varray<double>knot;
//	for (auto&v0 : v) {
//		knot.push_back(v0);
//	}
//
//	Convert_knot_occtype(knot, occknot, occmulit);
//
//	TColStd_Array1OfReal occwights(0, (*(spline._ControlPts)).size() - 1);
//	TColgp_Array1OfPnt cpts(0, (*(spline._ControlPts)).size() - 1);
//	for (int i = 0; i < (*(spline._ControlPts)).size(); i++) {
//
//		gp_Pnt pt = { (*(spline._ControlPts))[i].x,(*(spline._ControlPts))[i].y,(*(spline._ControlPts))[i].z };
//		cpts.SetValue(i, pt);
//		occwights.SetValue(i, (*(spline._ControlPts))[i].w);
//	}
//
//	occline = new Geom_BSplineCurve(cpts, occwights, occknot, occmulit, spline._u_Degree);
//
//}
//
//void OccWin::Convert_SplineSurface_Occtype(NurbsSurface & splinesurface, Handle(Geom_BSplineSurface)& occsurface)
//{
//	TColStd_Array1OfReal occuknot;
//	TColStd_Array1OfInteger occumulit;
//	TColStd_Array1OfReal occvknot;
//	TColStd_Array1OfInteger occvmulit;
//	auto uknot = *splinesurface._u_Knots;
//	auto vknot = *splinesurface._v_Knots;
//	varray<double>uknots, vknots;
//	for (auto&v : uknot) {
//		uknots.push_back(v);
//	}
//	for (auto&v : vknot) {
//		vknots.push_back(v);
//	}
//	Convert_knot_occtype(uknots, occuknot, occumulit);
//	Convert_knot_occtype(vknots, occvknot, occvmulit);
//	TColStd_Array2OfReal occwights(0, splinesurface._u_Num - 1, 0, splinesurface._v_Num - 1);
//	TColgp_Array2OfPnt cpts(0, splinesurface._u_Num - 1, 0, splinesurface._v_Num - 1);
//	for (int i = 0; i < splinesurface._v_Num; i++) {
//		for (int j = 0; j < splinesurface._u_Num; j++) {
//			int idx = i * splinesurface._v_Num + j;
//			gp_Pnt pt = { (*(splinesurface._ControlPts))[idx].x,(*(splinesurface._ControlPts))[idx].y,(*(splinesurface._ControlPts))[idx].z };
//			cpts.SetValue(j, i, pt);
//			occwights.SetValue(j, i, (*(splinesurface._ControlPts))[idx].w);
//		}
//	}
//
//	occsurface = new Geom_BSplineSurface(cpts, occwights, occuknot, occvknot,
//		occumulit, occvmulit, splinesurface._u_Num, splinesurface._v_Num, false, false);
//}
//
//void OccWin::Convert_SplineVolume_Occtype(NurbsVol & splinevolume, varray<Handle(Geom_BSplineSurface)>& occsurfaces)
//{
//	auto surfaces = GetSurfaces(splinevolume);
//	for (int i = 0; i < surfaces.size(); i++) {
//		Handle(Geom_BSplineSurface) sf;
//		Convert_SplineSurface_Occtype(surfaces[i], sf);
//		occsurfaces.push_back(sf);
//	}
//}
//
//void OccWin::Convert_Occline_Spline(NurbsLine & spline, Handle(Geom_BSplineCurve)& occline)
//{
//
//}
//
//void OccWin::Convert_OccSurface_SplineSurface(NurbsSurface & splinesurface, Handle(Geom_BSplineSurface)& occsurface)
//{
//}
//
//void OccWin::Convert_knot_occtype(varray<double>& knot, TColStd_Array1OfReal& occknot, TColStd_Array1OfInteger& occmulits)
//{
//	varray<double> newknot;
//	newknot.push_back(knot[0]);
//	varray<int> mulits;
//	int mulit = 1;
//	for (int i = 0; i < knot.size(); i++) {
//		if (i == 0) continue;
//		if (i > 0 && abs(knot[i] - knot[i - 1]) > 0.001) {
//			mulits.push_back(mulit);
//			newknot.push_back(knot[i]);
//			mulit = 1;
//		}
//		else {
//			mulit++;
//		}
//	}
//	mulits.push_back(mulit);
//
//	TColStd_Array1OfReal occknot1(0, newknot.size() - 1);
//	for (int i = 0; i < newknot.size(); i++) {
//		occknot1.SetValue(i, newknot[i]);
//	}
//	TColStd_Array1OfInteger occmulit1(0, mulits.size() - 1);
//	for (int i = 0; i < mulits.size(); i++) {
//		occmulit1.SetValue(i, mulits[i]);
//	}
//	occknot = std::move(occknot1);
//	occmulits = std::move(occmulit1);
//
//}

void OccWin::SelectModel(Standard_Integer x, Standard_Integer y)
{
	Handle(AIS_InteractiveContext) myContext = occView->GetInteractiveContext();
	Handle(V3d_View) myView = occView->GetView();
	bool flag = selectObj->SelectModel(myContext, myView, x, y, occView->getSelectMode());
	if (flag == true)
		occView->modelSelected();
}

void OccWin::InitMenu()
{
	QMenu* menuCreate = menuBar->addMenu(tr("&Model"));
	menuCreate->addAction(tr("&Read Model"));
	//menuCreate->addAction(tr("&save"));

	menuCreate = menuBar->addMenu(tr("&View"));
	menuCreate->addAction(tr("&Reset"));
	menuCreate->addAction(tr("&Fitall"));
	menuCreate->addAction(tr("&Clear"));

	menuCreate = menuBar->addMenu(tr("&InterActive"));
	menuCreate->addAction(tr("&Select"));
	menuCreate->addAction(tr("&Draw"));

	menuCreate = menuBar->addMenu(tr("&Test"));
	menuCreate = menuCreate->addMenu(tr("Topology"));
	menuCreate->addAction(tr("&Vertex"));
	menuCreate->addAction(tr("&Edge"));
	menuCreate->addAction(tr("&Face"));
	menuCreate->addAction(tr("&Wire"));
	menuCreate->addAction(tr("&Shell"));
	menuCreate->addAction(tr("&Solid"));

	menuCreate = menuBar->addMenu(tr("&FeatureTest"));
	menuCreate->addAction(tr("&FeatrueLineTest"));
	menuCreate->addAction(tr("&FeatrueRecTest_Edges"));
	menuCreate->addAction(tr("&FeatrueRecTest_Wire"));

	toolBar->addAction(tr("&Rectangle"));
	toolBar->addAction(tr("&Circle"));
	toolBar->addAction(tr("&OCCTest"));
	toolBar->addAction(tr("&Exit"));
	
}

void OccWin::trigerMenu(QAction* act)
{
	if (act->text() == "&Read Model") {
		OccReadModel();
		
		occView->GetView()->FitAll();
		occView->setViewMode(occView->ViewMode_View);
	}
	//else if (act->text() == "&Save") {
	//	actionsave();
	//}
	else if (act->text() == "&Reset") {
		occView->setViewMode(occView->ViewMode_View);
		viewReset();
	}
	else if (act->text() == "&Clear") {
		occView->setViewMode(occView->ViewMode_View);
		viewClear();
	}
	else if (act->text() == "&Fitall") {
		occView->setViewMode(occView->ViewMode_View);
		occView->GetView()->FitAll();
	}
	else if (act->text() == "&Select") {
		occView->setViewMode(occView->ViewMode_Select);
		occView->callSelectDialog();
		connect(occView, SIGNAL(myMouseClickXY(Standard_Integer, Standard_Integer)),
			this, SLOT(SelectModel(Standard_Integer, Standard_Integer)));
	}
	else if (act->text() == "&Draw") {
		occView->setViewMode(occView->ViewMode_Draw);
		toolBar->show();
		occView->GetView()->SetProj(V3d_Zpos);
		occView->GetView()->SetBackgroundColor(m_DrawBackground_Color);
	}
	else if (act->text() == "&Vertex") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->myDraw(myTest->myDraw_Vertex);
		DrawShape(myShape);
		DrawShape(myTest->myText);

	}
	else if (act->text() == "&Edge") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->myDraw(myTest->myDraw_Edge);
		DrawShape(myShape);
		DrawShape(myTest->myText);
	}
	else if (act->text() == "&Face") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->myDraw(myTest->myDraw_Face);
		DrawShape(myShape);
		DrawShape(myTest->myText);
	}
	else if (act->text() == "&Wire") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->myDraw(myTest->myDraw_Wire);
		DrawShape(myShape);
		DrawShape(myTest->myText);
	}
	else if (act->text() == "&Shell") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->myDraw(myTest->myDraw_Shell);
		DrawShape(myShape);
		DrawShape(myTest->myText);
	}
	else if (act->text() == "&Solid") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->myDraw(myTest->myDraw_Solid);
		DrawShape(myShape);
		DrawShape(myTest->myText);
	}
	else if (act->text() == "&FeatrueLineTest") {
		occView->setViewMode(occView->ViewMode_View);
		Handle(AIS_Shape) myShape = myTest->FeatureLineDraw(1, 1, 2, 2);
		DrawShape(myShape);
		DrawShape(myTest->myText);
	}
	else if (act->text() == "&FeatrueRecTest_Edges") {
		
		varray<Handle(AIS_Shape)> myShape = myTest->FeatureRecDraw(-10, 5, 20, 10);
		DrawShape(myShape);
	}
	
	else if (act->text() == "&FeatrueRecTest_Wire") {

		varray<Handle(AIS_Shape)> myShape = myTest->FeatureRecDraw(-30, -15, 0, -10, FeatureCurve_Rectangle::InitType::Wire);
		DrawShape(myShape);
	}


}

void OccWin::viewReset()
{
	if (!occView->GetView().IsNull())
		occView->GetView()->Reset();
}

void OccWin::viewClear()
{
	occView->GetInteractiveContext()->EraseAll(Standard_True);
}

void OccWin::OccReadModel()
{
	readModelDialog readFile(this);
	readFile.exec();
	if (readFile.isRead == false)
		return;
	string path = readFile.path;
	int modeltype = readFile.modelType;
	vector<OccCurve> lines;
	vector<OccSurface> surfaces;
	vector<OccVol> volumes;

	switch (modeltype) {
		//0,1是点，2是贝塞尔曲线，3是NURBS曲线，4是NURBS曲面，5是NURBS体
	case 0:

		break;
	case 1:

		break;
	case 2:

		break;
	case 3:
		ReadOccCurve(path, lines);
		for (int i = 0; i < lines.size(); i++) {
			Handle(Geom_BSplineCurve) occline;
			Convert_Spline_Occtype(lines[i], occline);
			occlines.push_back(occline);
		}
		for (int i = 0; i < occlines.size(); i++) {
			TopoDS_Edge L = BRepBuilderAPI_MakeEdge(occlines[i]);
			Handle(AIS_Shape) r = new AIS_Shape(L);
			occView->GetInteractiveContext()->Display(r, true);
		}
		break;
	case 4:
		ReadOccSurface(path, surfaces);
		for (int i = 0; i < surfaces.size(); i++) {
			Handle(Geom_BSplineSurface) occsurface;
			Convert_SplineSurface_Occtype(surfaces[i], occsurface);
			occsurfaces.push_back(occsurface);
		}
		for (int i = 0; i < occsurfaces.size(); i++) {
			TopoDS_Face NURBS = BRepBuilderAPI_MakeFace(occsurfaces[i], 0.3);
			Handle(AIS_Shape) NurbsShape = new AIS_Shape(NURBS);
			NurbsShape->SetColor(Quantity_NOC_RED);
			occView->GetInteractiveContext()->Display(NurbsShape, Standard_True);
		}
		break;
	case 5:
		ReadOccVol(path, volumes);
		for (int i = 0; i < volumes.size(); i++) {
			vector<Handle(Geom_BSplineSurface)> occsurface;
			Convert_SplineVolume_Occtype(volumes[i], occsurface);
			for (int j = 0; j < occsurface.size(); j++) {
				occsurfaces.push_back(occsurface[j]);
			}
		}
		for (int i = 0; i < occsurfaces.size(); i++) {
			TopoDS_Face NURBS = BRepBuilderAPI_MakeFace(occsurfaces[i], 0.0001);
			Handle(AIS_Shape) NurbsShape = new AIS_Shape(NURBS);
			NurbsShape->SetColor(Quantity_NOC_RED);
			occView->GetInteractiveContext()->Display(NurbsShape, Standard_True);
		}
		break;
	default:

		break;
	}

}

//void OccWin::actionsave()//保存为step
//{
//	STEPCAFControl_Writer aWriter = STEPCAFControl_Writer();
//	//Interface_Static::SetCVal("write.step.schema", "AP242DIS"));
//
//	const char* test;
//	IFSelect_ReturnStatus aStat = aWriter.Write("test.stp");
//	//return aStat;
//}

//varray<NurbsSurface> OccWin::GetSurfaces(const NurbsVol& NVS)			//将所有的面提取出来
//{
//	NurbsSurface U0;		//0
//	NurbsSurface U1;		//1
//	NurbsSurface V0;		//2
//	NurbsSurface V1;		//3
//	NurbsSurface W0;		//4
//	NurbsSurface W1;		//5
//	varray<varray<NurbsSurface>> Res;
//	varray<NurbsSurface> Sum;
//	int m = NVS._v_Num;
//	int n = NVS._u_Num;
//	int p = NVS._w_Num;
//	NurbsVol Vtemp = NVS;
//	int idx = 0;
//	for (int x = 0; x < m*n; x++)
//	{
//		(*(W0._ControlPts)).push_back((*(NVS._ControlPts))[x]);
//	}
//	W0.SetSurface(Vtemp._u_Degree, Vtemp._v_Degree, Vtemp._u_Num, Vtemp._v_Num, *Vtemp._u_Knots, *Vtemp._v_Knots);
//	for (int x = m * n*(p - 1); x < m*n*p; x++)
//	{
//		(*(W1._ControlPts)).push_back((*(NVS._ControlPts))[x]);
//	}
//	W1.SetSurface(Vtemp._u_Degree, Vtemp._v_Degree, Vtemp._u_Num, Vtemp._v_Num, *Vtemp._u_Knots, *Vtemp._v_Knots);
//
//	for (int y = 0; y <= (p - 1)*m*n; y += m * n) {
//		for (int x = y; x <= y + m * n - m; x += m) {
//			(*(V0._ControlPts)).push_back((*(NVS._ControlPts))[x]);
//		}
//	}
//	V0.SetSurface(Vtemp._w_Degree, Vtemp._u_Degree, Vtemp._w_Num, Vtemp._u_Num, *Vtemp._w_Knots, *Vtemp._u_Knots);
//	for (int y = m - 1; y <= m * n*(p - 1) + m - 1; y += m * n)
//		for (int x = y; x <= y + m * n - m; x += m) {
//			(*(V1._ControlPts)).push_back((*(NVS._ControlPts))[x]);
//		}
//	V1.SetSurface(Vtemp._w_Degree, Vtemp._u_Degree, Vtemp._w_Num, Vtemp._u_Num, *Vtemp._w_Knots, *Vtemp._u_Knots);
//
//	for (int y = 0; y < p; y++) {
//		for (int x = y * m*n; x < y*m*n + m; x++)
//		{
//			(*(U0._ControlPts)).push_back((*(NVS._ControlPts))[x]);
//		}
//	}
//	U0.SetSurface(Vtemp._w_Degree, Vtemp._v_Degree, Vtemp._w_Num, Vtemp._v_Num, *Vtemp._w_Knots, *Vtemp._v_Knots);
//
//	for (int y = 0; y < p; y++) {
//		for (int x = y * m*n + m * n - m; x < y*m*n + m * n; x++)
//		{
//			(*(U1._ControlPts)).push_back((*(NVS._ControlPts))[x]);
//		}
//	}
//	U1.SetSurface(Vtemp._w_Degree, Vtemp._v_Degree, Vtemp._w_Num, Vtemp._v_Num, *Vtemp._w_Knots, *Vtemp._v_Knots);
//
//
//	Sum.push_back(V1);
//	(*(V1._ControlPts)).clear();
//	Sum.push_back(V0);
//	(*(V0._ControlPts)).clear();
//
//	Sum.push_back(U0);
//	(*(U0._ControlPts)).clear();
//
//	Sum.push_back(U1);
//	(*(U1._ControlPts)).clear();
//
//	Sum.push_back(W0);
//	(*(W0._ControlPts)).clear();
//
//	Sum.push_back(W1);
//	(*(W1._ControlPts)).clear();
//
//	return Sum;
//}

void OccWin::draw2D(Standard_Integer mode)
{
	connect(occView, static_cast<void(OccView::*)(Standard_Integer, Standard_Integer)>(&OccView::myMouseClickXY),
		[=](Standard_Integer x, Standard_Integer y) {
		drawObj->setInitXY(x, y);
	});

	connect(occView, static_cast<void(OccView::*)(Standard_Integer, Standard_Integer)>(&OccView::myMouseMoveXY),
		[=](Standard_Integer x, Standard_Integer y) {
		drawObj->setMoveXY(x, y);
		Handle(AIS_Shape) myShape;
		switch (mode) {
		case 0:
			myShape = drawObj->drawRectangle();
			break;
		case 1:
			myShape = drawObj->drawCircle();
			break;
		}
		if (!myShape.IsNull())
			if (!drawObj->prevShape.IsNull())
				occView->GetInteractiveContext()->Remove(drawObj->prevShape, Standard_True);
		myShape->SetColor(m_Draw_Color);
		occView->GetInteractiveContext()->Display(myShape, Standard_True);
		drawObj->prevShape = myShape;
		occView->GetInteractiveContext()->UpdateCurrentViewer();
	});

	connect(occView, &OccView::drawComplete, [=] {
		drawObj->prevShape.Nullify();
	});
}

void OccWin::exitPrevDraw()
{
	disconnect(occView, &OccView::myMouseClickXY, 0, 0);
	disconnect(occView, &OccView::myMouseMoveXY, 0, 0);
	disconnect(occView, &OccView::drawComplete, 0, 0);
}

void OccWin::DrawShape(Handle(AIS_InteractiveObject) shape)
{
	occView->GetInteractiveContext()->Display(shape, Standard_True);
}

void OccWin::DrawShape(varray<Handle(AIS_Shape)> shapes)
{
	for (int i = 0; i < shapes.size(); i++)
	{
		occView->GetInteractiveContext()->Display(shapes[i], Standard_True);
	}
}



void OccWin::trigerToolBar(QAction* action)
{
	if (action->text() == "&Rectangle") {
		actiondrawRectangle();
	}
	else if (action->text() == "&Circle") {
		actiondrawCircle();
	}
	else if (action->text() == "&OCCTest") {
		actionOCCTest();
	}
	else if (action->text() == "&Exit") {
		actionexit();
	}
	
}

void OccWin::actiondrawRectangle()
{
	exitPrevDraw();
	draw2D(drawObj->DrawMode_Rectangle);
}

void OccWin::actiondrawCircle()
{
	exitPrevDraw();
	draw2D(drawObj->DrawMode_Circle);
}


void OccWin::actionOCCTest()
{
	

	// Make a vertex from a 3D point.
	gp_Pnt aPnt(0.0, 0.0, 10.0);
	TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(aPnt);
	std::cout << "TopoDS_Vertex was created at [ "
		<< aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z()
		<< " ]" << std::endl;

	Handle(AIS_Shape) aAisVertex = new AIS_Shape(aVertex);
	Handle(AIS_TextLabel) anAisLabel = new AIS_TextLabel();
	Standard_SStream aSS;
	aSS << "TopoDS_Vertex [" << aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z() << "]" << std::endl;
	anAisLabel->SetText(aSS.str().c_str());
	anAisLabel->SetPosition(aPnt);
	occView->GetInteractiveContext()->Display(aAisVertex, Standard_True);
	occView->GetInteractiveContext()->Display(anAisLabel, Standard_True);

}

void OccWin::actionexit()
{
	exitPrevDraw();
	toolBar->hide();
	occView->GetView()->SetBackgroundColor(m_View_Color);
	occView->setViewMode(occView->ViewMode_View);
}


//读取NURBS曲线
//返回曲线数量
int OccWin::ReadOccCurve(const string& path, vector<OccCurve>& lines)
{
	lines.clear();
	ifstream inf(path, ios::binary);
	if (!inf)return 0;
	string temp = "";
	int flag = 0;
	OccCurve line;
	while (!inf.eof())
	{
		inf >> temp;
		//判断标识
		if (temp == "<idx>") { inf >> temp; continue; }
		if (temp == "<line/>")flag = 1;
		else if (temp == "<degree>")flag = 2;
		else if (temp == "<knots/>") { flag = 3; continue; }
		else if (temp == "</knots>")continue;
		else if (temp == "<ctrlPts/>") { flag = 4; continue; }
		else if (temp == "</ctrlPts>")continue;
		else if (temp == "</line>")flag = 5;
		else if (temp == "<allend>")break;
		//根据标识对temp操作
		switch (flag)
		{
		case 1://<line/>
			line.m_Knots.clear();
			line.m_CtrlPts.clear();
			break;
		case 2://<degree>
			inf >> line.m_Degree;
			break;
		case 3://<knots/>
		{
			double k = atof(temp.c_str());
			line.m_Knots.push_back(k);
			break;
		}
		case 4://<ctrlPts/>
		{
			point pts;
			pts.x = atof(temp.c_str());
			inf >> pts.y;
			inf >> pts.z;
			inf >> pts.w;
			line.m_CtrlPts.push_back(pts);
			break;
		}
		case 5://</line>
		{
			lines.push_back(line);
			break;
		}
		default:
			break;
		}
	}
	inf.close();
	return lines.size();
}

//读取NURBS曲面
//返回曲面数量
int OccWin::ReadOccSurface(const string& path, vector<OccSurface>& surfaces)
{
	surfaces.clear();
	ifstream inf(path, ios::binary);
	if (!inf)return 0;
	string temp = "";
	int flag = 0;
	OccSurface surface;
	while (!inf.eof())
	{
		inf >> temp;
		//判断标识
		if (temp == "<idx>") { inf >> temp; continue; }
		if (temp == "<surface/>")flag = 1;
		else if (temp == "<udegree>")flag = 2;
		else if (temp == "<vdegree>")flag = 3;
		else if (temp == "<uNum>")flag = 4;
		else if (temp == "<vNum>")flag = 5;
		else if (temp == "<uknots/>") { flag = 6; continue; }
		else if (temp == "</uknots>")continue;
		else if (temp == "<vknots/>") { flag = 7; continue; }
		else if (temp == "</vknots>")continue;
		else if (temp == "<ctrlPts/>") { flag = 8; continue; }
		else if (temp == "</ctrlPts>")continue;
		else if (temp == "</surface>")flag = 9;
		else if (temp == "<allend>")break;
		//根据标识对temp操作
		switch (flag)
		{
		case 1://<surface/>
			surface.m_uKnots.clear();
			surface.m_vKnots.clear();
			surface.m_CtrlPts.clear();
			break;
		case 2://<udegree>
			inf >> surface.m_uDegree;
			break;
		case 3://<vdegree>
			inf >> surface.m_vDegree;
			break;
		case 4://<uNum>
			inf >> surface.m_uNum;
			break;
		case 5://<vNum>
			inf >> surface.m_vNum;
			break;
		case 6://<uknots/>
		{
			double k = atof(temp.c_str());
			surface.m_uKnots.push_back(k);
			break;
		}
		case 7://<vknots/>
		{
			double k = atof(temp.c_str());
			surface.m_vKnots.push_back(k);
			break;
		}
		case 8://<ctrlPts/>
		{
			point pts;
			pts.x = atof(temp.c_str());
			inf >> pts.y;
			inf >> pts.z;
			inf >> pts.w;
			surface.m_CtrlPts.push_back(pts);
			break;
		}
		case 9://</surface>
		{
			surfaces.push_back(surface);
			break;
		}
		default:
			break;
		}
	}
	inf.close();
	return surfaces.size();
}

//读取NURBS体
//返回体模型数量
int OccWin::ReadOccVol(const string& path, vector<OccVol>& vols)
{
	vols.clear();
	ifstream inf(path, ios::binary);
	if (!inf)return 0;
	string temp = "";
	int flag = 0;
	OccVol vol;
	while (!inf.eof())
	{
		inf >> temp;
		//判断标识
		if (temp == "<idx>") { inf >> temp; continue; }
		if (temp == "<vol/>")flag = 1;
		else if (temp == "<udegree>")flag = 2;
		else if (temp == "<vdegree>")flag = 3;
		else if (temp == "<wdegree>")flag = 4;
		else if (temp == "<uNum>")flag = 5;
		else if (temp == "<vNum>")flag = 6;
		else if (temp == "<wNum>")flag = 7;
		else if (temp == "<uknots/>") { flag = 8; continue; }
		else if (temp == "</uknots>")continue;
		else if (temp == "<vknots/>") { flag = 9; continue; }
		else if (temp == "</vknots>")continue;
		else if (temp == "<wknots/>") { flag = 10; continue; }
		else if (temp == "</wknots>")continue;
		else if (temp == "<ctrlPts/>") { flag = 11; continue; }
		else if (temp == "</ctrlPts>")continue;
		else if (temp == "</vol>")flag = 12;
		else if (temp == "<allend>")break;
		//根据标识对temp操作
		switch (flag)
		{
		case 1://<surface/>
			vol.m_uKnots.clear();
			vol.m_vKnots.clear();
			vol.m_wKnots.clear();
			vol.m_CtrlPts.clear();
			break;
		case 2://<udegree>
			inf >> vol.m_uDegree;
			break;
		case 3://<vdegree>
			inf >> vol.m_vDegree;
			break;
		case 4://<wdegree>
			inf >> vol.m_wDegree;
			break;
		case 5://<uNum>
			inf >> vol.m_uNum;
			break;
		case 6://<vNum>
			inf >> vol.m_vNum;
			break;
		case 7://<wNum>
			inf >> vol.m_wNum;
			break;
		case 8://<uknots/>
		{
			double k = atof(temp.c_str());
			vol.m_uKnots.push_back(k);
			break;
		}
		case 9://<vknots/>
		{
			double k = atof(temp.c_str());
			vol.m_vKnots.push_back(k);
			break;
		}
		case 10://<wknots/>
		{
			double k = atof(temp.c_str());
			vol.m_wKnots.push_back(k);
			break;
		}
		case 11://<ctrlPts/>
		{
			point pts;
			pts.x = atof(temp.c_str());
			inf >> pts.y;
			inf >> pts.z;
			inf >> pts.w;
			vol.m_CtrlPts.push_back(pts);
			break;
		}
		case 12://</vol>
		{
			vols.push_back(vol);
			break;
		}
		default:
			break;
		}
	}
	inf.close();
	return vols.size();
}