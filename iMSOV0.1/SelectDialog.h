#pragma once

#include <QtWidgets/qdialog.h>
#include "ui_selectdialog.h"
#include "drawobj.h"

class SelectDialog : public QDialog
{
	Q_OBJECT

public:
	explicit SelectDialog(QWidget *parent = nullptr);
	~SelectDialog();
signals:
	void sendSelectMode(Standard_Integer mode);
	void cancelSelectMode();
private:
	Ui::SelectDialog *ui;
	Standard_Integer curmode;
};
