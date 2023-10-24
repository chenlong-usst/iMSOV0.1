#pragma once

#include <QMainWindow>
#include "ui_OCCMainWindow.h"
#include "occwin.h"

class OCCMainWindow : public QMainWindow
{
	Q_OBJECT

public:
	OCCMainWindow(QWidget *parent = nullptr);
	~OCCMainWindow();

private:
	Ui::OCCMainWindowClass ui;
};
