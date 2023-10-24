#include "OCCMainWindow.h"

OCCMainWindow::OCCMainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	setWindowState(Qt::WindowMaximized);
	//QMainWindow *mainWindow = new QMainWindow();    // ����һ��QWidget����
	QDialog *myDialog = new OccWin(this);
	this->setCentralWidget(myDialog);
	myDialog->exec();
	//QWidget *myWidget = new OccWin(this);    // ��myWidget����ΪmainWindow�����Ĵ��ڲ���
	//this->setCentralWidget(myWidget);    // ��ʾ������

	//this->show();
}

OCCMainWindow::~OCCMainWindow()
{}
