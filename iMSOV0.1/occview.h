#pragma once

#include <QtWidgets/qwidget.h>
#include <QMouseEvent>
#include <AIS_InteractiveContext.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <V3d_Viewer.hxx>
#include <V3d_View.hxx>
#include "selectdialog.h"
#include "drawobj.h"
#include <Prs3d_Drawer.hxx>
#include <Prs3d_LineAspect.hxx>

namespace Ui {
	class OccView;
}

class OccView : public QWidget
{
	Q_OBJECT

public:
	explicit OccView(QWidget *parent = nullptr);
	~OccView();
private:
	Ui::OccView *ui;
	Handle(AIS_InteractiveContext) m_context;
	Handle(V3d_Viewer) m_viewer;
	Handle(V3d_View) m_view;
	Handle(Graphic3d_GraphicDriver) m_graphic_driver;
public:
	QPaintEngine* paintEngine() const;
	Handle(AIS_InteractiveContext) GetInteractiveContext();
	Handle(V3d_View) GetView();

	gp_Pnt ConvertClickToPoint(Standard_Real theX, Standard_Real theY, Handle(V3d_View) theView);

protected:

	void paintEvent(QPaintEvent* event);
	void resizeEvent(QResizeEvent* event);
	void mousePressEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void mouseMoveEvent(QMouseEvent *event);

protected:
	enum CurrentAction
	{
		CurrentAction_Nothing,
		CurrentAction_Panning,  //平移，鼠标左键
		CurrentAction_Zooming,  //缩放，鼠标滚轮
		CurrentAction_Rotation  //旋转，按住鼠标中间并拖动
	};
public:
	//当前视野处于哪个模式
	enum ViewMode
	{
		ViewMode_View,   //浏览模式
		ViewMode_Draw,   //画草图模式
		ViewMode_Select  //选择模式
	};
private:
	Standard_Integer m_x;
	Standard_Integer m_y;
	Standard_Integer m_select_mode;
	CurrentAction m_current_mode;
	ViewMode m_view_mode;
	Standard_Boolean m_isSelected;

signals:
	void myMouseClickXY(Standard_Integer x, Standard_Integer y);
	void myMouseMoveXY(Standard_Integer x, Standard_Integer y);
	void drawComplete();


public:
	void setViewMode(ViewMode viewmode);
	void callSelectDialog();
	void setProj();
	void modelSelected();
	Standard_Integer getSelectMode();
	Standard_Integer getViewMode();
private:
	SelectDialog *selectDialog;

	void selectedModelMove();

private slots:
	void setSelectMode(Standard_Integer x);




};

// OCCVIEW_H
