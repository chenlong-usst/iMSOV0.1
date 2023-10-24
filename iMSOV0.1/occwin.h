#pragma once


#include <QtWidgets/qdialog.h>
#include <QtWidgets/qtoolbar.h>
#include <QtWidgets/qmenubar.h>
#include <QtWidgets/qboxlayout.h>
#include <QtWidgets/qdialog.h>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <AIS_Shape.hxx>
#include <AIS_TextLabel.hxx>

#include "ui_occwin.h"
#include "occview.h"
#include "occModel.h"
#include "selectobj.h"
#include "occTest.h"
#include "varray.h"
#define m_View_Color Quantity_NOC_GRAY60
//画图时的背景颜色
#define m_DrawBackground_Color Quantity_NOC_BLACK
//草图颜色
#define m_Draw_Color Quantity_NOC_WHITE
//被选中时线的颜色
#define m_Select_LineColor Quantity_NOC_GREEN


class OccWin : public QDialog
{
	Q_OBJECT

public:
	OccWin(QWidget *parent = Q_NULLPTR);
	~OccWin();
public:

	//将本程序的NURBS数据转换成Occ需要的NURBS数据，以便进行显示
	void Convert_Spline_Occtype(OccCurve& spline, Handle(Geom_BSplineCurve)& occline);
	void Convert_SplineSurface_Occtype(OccSurface& splinesurface, Handle(Geom_BSplineSurface)& occsurface);
	void Convert_SplineVolume_Occtype(OccVol& splinevolume, vector<Handle(Geom_BSplineSurface)>& occsurfaces);
	void Convert_knot_occtype(vector<double>& knot, TColStd_Array1OfReal& occknot, TColStd_Array1OfInteger& occmulits);
private slots:
	void trigerMenu(QAction*);
	void trigerToolBar(QAction*);
	void SelectModel(Standard_Integer x, Standard_Integer y);

private:

	void InitMenu();
	void OccReadModel();
	//void OccShowModel();
	void viewReset();
	void viewClear();
	vector<OccSurface>GetSurfaces(const OccVol& NVS);
	void actiondrawRectangle();
	void actiondrawCircle();

	void actionOCCTest();

	void actionexit();
	/*void actionsave();*/
	//读取NURBS曲线
//返回曲线数量
	int ReadOccCurve(const string& path, vector<OccCurve>& lines);

	//读取NURBS曲面
	//返回曲面数量
	int ReadOccSurface(const string& path, vector<OccSurface>& surfaces);

	//读取NURBS体
	//返回体模型数量
	int ReadOccVol(const string& path, vector<OccVol>& vols);
	//绘制二维图形
	void draw2D(Standard_Integer mode);
	//退出上一次的绘制
	void exitPrevDraw();

	void DrawShape(Handle(AIS_InteractiveObject)shape);
	void DrawShape(varray<Handle(AIS_Shape)>shapes);
private:
	vector<Handle(Geom_BSplineCurve)> occlines;
	vector<Handle(Geom_BSplineSurface)> occsurfaces;
	OccView* occView;
	QMenuBar* menuBar;
	QVBoxLayout* layout;
	QAction* myAction;
	SelectObj* selectObj;
	DrawObj* drawObj;
	QToolBar* toolBar;
	occTest* myTest;
private:
	Ui::OccWin ui;
};
