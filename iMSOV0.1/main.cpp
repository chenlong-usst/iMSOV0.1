#include <iostream>
//#include <QApplication>

#include "OCCMainWindow.h"
int main(int argc, char *argv[])
{
	QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
	QApplication a(argc, argv);

	//OCCMainWindow o;
	//o.show();
	QDialog *myDialog = new OccWin();
	//this->setCentralWidget(myDialog);
	myDialog->exec();


	//system("pause");
	return a.exec();

}