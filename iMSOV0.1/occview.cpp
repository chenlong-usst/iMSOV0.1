#include "occview.h"
#include <ui_occview.h>
#include <WNT_Window.hxx>
#include <ProjLib.hxx>
#include <ElSLib.hxx>
#include <QDebug>
OccView::OccView(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::OccView), m_isSelected(Standard_False)
{
	ui->setupUi(this);
	selectDialog = new SelectDialog(this);
	if (m_context.IsNull()) {
		Handle(Aspect_DisplayConnection) m_display_connection = new Aspect_DisplayConnection();
		if (m_graphic_driver.IsNull())
		{
			m_graphic_driver = new OpenGl_GraphicDriver(m_display_connection);
		}
		WId window_handle = static_cast<WId>(winId());
		Handle(WNT_Window) wind = new WNT_Window(reinterpret_cast<Aspect_Handle>(window_handle));

		m_viewer = new V3d_Viewer(m_graphic_driver);
		m_view = m_viewer->CreateView();
		m_view->SetWindow(wind);
		if (!wind->IsMapped()) wind->Map();
		m_context = new AIS_InteractiveContext(m_viewer);
		m_viewer->SetDefaultLights();
		m_viewer->SetLightOn();
		m_view->SetBackgroundColor(Quantity_NOC_GRAY60);
		m_view->MustBeResized();
		m_view->TriedronDisplay(Aspect_TOTP_LEFT_LOWER, Quantity_NOC_GOLD, 0.08, V3d_ZBUFFER);
		m_context->SetDisplayMode(AIS_Shaded, Standard_True);
		//m_context->SetSelectionStyle();

		// 设置选择模型的风格
		Handle(Prs3d_Drawer) t_select_style = m_context->SelectionStyle();  // 获取选择风格
		t_select_style->SetMethod(Aspect_TOHM_COLOR);  // 颜色显示方式
		t_select_style->SetColor(Quantity_NOC_LIGHTSEAGREEN);   // 设置选择后颜色
		t_select_style->SetDisplayMode(1); // 整体高亮
		t_select_style->SetTransparency(0.4f);
	}
	setAttribute(Qt::WA_PaintOnScreen);
	setAttribute(Qt::WA_NoSystemBackground);
	setBackgroundRole(QPalette::NoRole);
	setFocusPolicy(Qt::StrongFocus);

	setMouseTracking(true);

	connect(selectDialog, SIGNAL(sendSelectMode(Standard_Integer)), this, SLOT(setSelectMode(Standard_Integer)));
	connect(selectDialog, &SelectDialog::cancelSelectMode, [=]() {
		selectDialog->close();
	});

}

OccView::~OccView()
{
	delete selectDialog;
	delete ui;
}

QPaintEngine* OccView::paintEngine() const
{
	return 0;
}

void OccView::setViewMode(ViewMode viewmode)
{
	m_view_mode = viewmode;
}

void OccView::callSelectDialog()
{
	selectDialog->show();
}

Standard_Integer OccView::getSelectMode()
{
	return m_select_mode;
}

Standard_Integer OccView::getViewMode()
{
	return m_view_mode;
}

void OccView::selectedModelMove()
{
	qDebug() << "success";
}

void OccView::setSelectMode(Standard_Integer x)
{
	m_select_mode = x;
}


void OccView::setProj()
{
	m_view->SetFront();
}

void OccView::modelSelected()
{
	m_isSelected = Standard_True;
}

Handle(AIS_InteractiveContext) OccView::GetInteractiveContext()
{
	return m_context;
}

Handle(V3d_View) OccView::GetView()
{
	return m_view;
}

void OccView::paintEvent(QPaintEvent *event)
{
	m_view->Redraw();
}

void OccView::resizeEvent(QResizeEvent *event)
{
	if (!m_view.IsNull()) {
		m_view->MustBeResized();
	}
}

void OccView::mousePressEvent(QMouseEvent *event)
{
	switch (m_view_mode) {
	case ViewMode_View:
		//左键拖动
		if (event->buttons() & Qt::LeftButton)
		{
			m_current_mode = CurrentAction_Panning;
			m_x = event->pos().x();
			m_y = event->pos().y();
		}
		else if (event->buttons() & Qt::MidButton)
		{
			//按住中键
			m_current_mode = CurrentAction_Rotation;
			m_view->StartRotation(event->pos().x(), event->pos().y());
		}
		break;
	case ViewMode_Draw:
		//TODO
		if (event->buttons() & Qt::LeftButton) {
			gp_Pnt p = ConvertClickToPoint(event->pos().x(), event->pos().y(), m_view);
			myMouseClickXY(static_cast<int>(p.X()), static_cast<int>(p.Y()));
		}
		break;
	case ViewMode_Select:
		if (event->buttons() & Qt::LeftButton)
			myMouseClickXY(event->pos().x(), event->pos().y());
		break;
	default:
		break;
	}


}

void OccView::mouseReleaseEvent(QMouseEvent *event)
{
	m_current_mode = CurrentAction_Nothing;
	if (m_view_mode == ViewMode_Draw)
		emit drawComplete();
	else if (m_view_mode == ViewMode_Select)
		m_isSelected = Standard_False;
}

void OccView::wheelEvent(QWheelEvent *event)
{
	if (m_view_mode == ViewMode_View) {
		m_view->StartZoomAtPoint(event->pos().x(), event->pos().y());
		m_view->ZoomAtPoint(0, 0, event->angleDelta().y(), 0);
	}
}

void OccView::mouseMoveEvent(QMouseEvent *event)
{
	if (m_view_mode == ViewMode_View)
		switch (m_current_mode) {
		case CurrentAction_Panning:
			m_view->Pan(event->pos().x() - m_x, m_y - event->pos().y());
			m_x = event->pos().x();
			m_y = event->pos().y();
			break;
		case CurrentAction_Rotation:
			m_view->Rotation(event->pos().x(), event->pos().y());
			break;
		default:
			break;
		}
	else if (m_view_mode == ViewMode_Draw) {
		if (event->buttons() & Qt::LeftButton) {
			gp_Pnt p = ConvertClickToPoint(event->pos().x(), event->pos().y(), m_view);
			myMouseMoveXY(static_cast<int>(p.X()), static_cast<int>(p.Y()));
		}
		m_x = event->pos().x();
		m_y = event->pos().y();
	}
	else if (m_view_mode == ViewMode_Select) {
		if ((event->buttons() & Qt::LeftButton) && m_isSelected)
			selectedModelMove();
	}
}

//二维画图点转换，窗口坐标转换到世界坐标
gp_Pnt OccView::ConvertClickToPoint(Standard_Real theX, Standard_Real theY, Handle(V3d_View) theView) {
	Standard_Real XEye, YEye, ZEye, XAt, YAt, ZAt;
	theView->Eye(XEye, YEye, ZEye);
	theView->At(XAt, YAt, ZAt);
	gp_Pnt EyePoint(XEye, YEye, ZEye);
	gp_Pnt AtPoint(XAt, YAt, ZAt);

	gp_Vec EyeVector(EyePoint, AtPoint);
	gp_Dir EyeDir(EyeVector);

	gp_Pln PlaneOfTheView = gp_Pln(AtPoint, EyeDir);
	Standard_Real X, Y, Z;
	theView->Convert(static_cast<int>(theX), static_cast<int>(theY), X, Y, Z);
	gp_Pnt ConvertedPoint(X, Y, Z);
	gp_Pnt2d ConvertedPointOnPlane = ProjLib::Project(PlaneOfTheView, ConvertedPoint);

	gp_Pnt ResultPoint = ElSLib::Value(ConvertedPointOnPlane.X(),
		ConvertedPointOnPlane.Y(),
		PlaneOfTheView);
	return ResultPoint;
}

