#include "OCCMainWindow.h"

OCCMainWindow::OCCMainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	setWindowState(Qt::WindowMaximized);
	//QMainWindow *mainWindow = new QMainWindow();    // 创建一个QWidget对象
	QDialog *myDialog = new OccWin(this);
	this->setCentralWidget(myDialog);
	myDialog->exec();
	//QWidget *myWidget = new OccWin(this);    // 将myWidget设置为mainWindow的中心窗口部件
	//this->setCentralWidget(myWidget);    // 显示主窗口

	//this->show();
}

OCCMainWindow::~OCCMainWindow()
{}
