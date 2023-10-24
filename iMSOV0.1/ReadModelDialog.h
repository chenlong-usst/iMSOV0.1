#pragma once

#include <QtWidgets/qdialog.h>
#include "ui_readmodeldialog.h"
#include <string>
class readModelDialog : public QDialog
{
	Q_OBJECT

public:
	readModelDialog(QWidget *parent = Q_NULLPTR);
	~readModelDialog();
	bool isRead = false;
	int modelType;
	std::string path;
private:
	Ui::readModelDialog ui;
};
