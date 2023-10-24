#include "occview.h"
#include "ui_occview.h"
#include <WNT_Window.hxx>

OccView::OccView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::OccView)
{
    ui->setupUi(this);
    if(m_context.IsNull()){
        Handle(Aspect_DisplayConnection) m_display_connection = new Aspect_DisplayConnection();
        if(m_graphic_driver.IsNull())
        {
            m_graphic_driver = new OpenGl_GraphicDriver(m_display_connection);

        }
        WId window_handle = (WId)winId();
        Handle(WNT_Window) wind = new WNT_Window((Aspect_Handle)window_handle);

        m_viewer = new V3d_Viewer(m_graphic_driver);
        m_view = m_viewer->CreateView();
        m_view->SetWindow(wind);
        if(!wind->IsMapped()) wind->Map();
        m_context = new AIS_InteractiveContext(m_viewer);
        m_viewer->SetDefaultLights();
        m_viewer->SetLightOn();
        m_view->SetBackgroundColor(Quantity_NOC_GRAY60);
        m_view->MustBeResized();
        m_view->TriedronDisplay(Aspect_TOTP_LEFT_LOWER, Quantity_NOC_GOLD, 0.08, V3d_ZBUFFER);
        m_context->SetDisplayMode(AIS_Shaded, Standard_True);
    }
    setAttribute(Qt::WA_PaintOnScreen);
    setAttribute(Qt::WA_NoSystemBackground);
    setBackgroundRole(QPalette::NoRole);
    setFocusPolicy(Qt::StrongFocus);

    setMouseTracking(true);
}

OccView::~OccView()
{
    delete ui;
}

QPaintEngine* OccView::paintEngine() const
{
    return 0;
}

void OccView::paintEvent(QPaintEvent *event)
{
    m_view->Redraw();
}

void OccView::resizeEvent(QResizeEvent *event)
{
    if(!m_view.IsNull()){
        m_view->MustBeResized();
    }
}
