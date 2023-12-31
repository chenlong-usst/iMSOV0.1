#include "ReadModelDialog.h"
#include <QtWidgets/qfiledialog.h>
#include <QtWidgets/qmessagebox.h>

readModelDialog::readModelDialog(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);

	connect(ui.cancle, &QPushButton::clicked, this, &QWidget::close);

	connect(ui.selectFile, &QPushButton::clicked, [=]() {
		QString filePath = QFileDialog::getOpenFileName(this, QString::fromLocal8Bit("选择文件"), "", "*.txt");
		ui.filePath->setText(filePath);
	});

	connect(ui.conform, &QPushButton::clicked, [=]() {
		if (ui.filePath->text() == NULL)
		{
			QMessageBox box(QMessageBox::Warning, QString::fromLocal8Bit("警告"), QString::fromLocal8Bit("路径不能为空"));
			box.setStandardButtons(QMessageBox::Ok);
			box.setButtonText(QMessageBox::Ok, QString::fromLocal8Bit("确 定"));
			box.exec();
		}
		else {
			modelType = ui.modelType->currentIndex();
			path = ui.filePath->text().toStdString();
			isRead = true;
			this->close();
		}
	});
}

readModelDialog::~readModelDialog()
{
	
}
