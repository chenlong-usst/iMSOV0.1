#include <QtWidgets/qpushbutton.h>

#include "selectdialog.h"
#include <ui_selectdialog.h>
#include "selectobj.h"

SelectDialog::SelectDialog(QWidget *parent) :
	QDialog(parent),
	ui(new Ui::SelectDialog)
{

	ui->setupUi(this);
	connect(ui->select_face, &QRadioButton::clicked, [=]() {
		curmode = SelectObj::SelectMode_Face;
	});

	connect(ui->select_line, &QRadioButton::clicked, [=]() {
		curmode = SelectObj::SelectMode_Line;
	});


	connect(ui->select_lines, &QRadioButton::clicked, [=]() {
		curmode = SelectObj::SelectMode_Boundary;
	});

	connect(ui->select_point, &QRadioButton::clicked, [=]() {
		curmode = SelectObj::SelectMode_Point;
	});

	connect(ui->select_whole, &QRadioButton::clicked, [=]() {
		curmode = SelectObj::SelectMode_Whole;
	});

	connect(ui->button_ok, &QPushButton::clicked, [=]() {
		emit sendSelectMode(curmode);
		emit cancelSelectMode();
	});

	connect(ui->button_cancel, &QPushButton::clicked, [=]() {
		emit cancelSelectMode();
	});



}

SelectDialog::~SelectDialog()
{
	delete ui;
}
