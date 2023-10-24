#ifndef OCCVIEW_H
#define OCCVIEW_H

#include <QWidget>
#include <AIS_InteractiveContext.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <V3d_Viewer.hxx>
#include <V3d_View.hxx>
#include <Quantity_NameOfColor.hxx>

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


protected:

    void paintEvent(QPaintEvent* event);
    void resizeEvent(QResizeEvent* event);

};

#endif // OCCVIEW_H
