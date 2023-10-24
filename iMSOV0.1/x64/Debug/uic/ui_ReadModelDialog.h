/********************************************************************************
** Form generated from reading UI file 'ReadModelDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_READMODELDIALOG_H
#define UI_READMODELDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_readModelDialog
{
public:
    QPushButton *conform;
    QPushButton *cancle;
    QWidget *widget_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QComboBox *modelType;
    QLabel *label;
    QLineEdit *filePath;
    QPushButton *selectFile;

    void setupUi(QDialog *readModelDialog)
    {
        if (readModelDialog->objectName().isEmpty())
            readModelDialog->setObjectName(QStringLiteral("readModelDialog"));
        readModelDialog->resize(419, 300);
        conform = new QPushButton(readModelDialog);
        conform->setObjectName(QStringLiteral("conform"));
        conform->setGeometry(QRect(120, 210, 93, 28));
        cancle = new QPushButton(readModelDialog);
        cancle->setObjectName(QStringLiteral("cancle"));
        cancle->setGeometry(QRect(220, 210, 93, 28));
        widget_3 = new QWidget(readModelDialog);
        widget_3->setObjectName(QStringLiteral("widget_3"));
        widget_3->setGeometry(QRect(10, 120, 378, 54));
        horizontalLayout_2 = new QHBoxLayout(widget_3);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_2 = new QLabel(widget_3);
        label_2->setObjectName(QStringLiteral("label_2"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy);

        horizontalLayout_2->addWidget(label_2);

        modelType = new QComboBox(widget_3);
        modelType->setObjectName(QStringLiteral("modelType"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(modelType->sizePolicy().hasHeightForWidth());
        modelType->setSizePolicy(sizePolicy1);

        horizontalLayout_2->addWidget(modelType);

        label = new QLabel(readModelDialog);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(27, 58, 60, 33));
        filePath = new QLineEdit(readModelDialog);
        filePath->setObjectName(QStringLiteral("filePath"));
        filePath->setGeometry(QRect(94, 64, 189, 21));
        selectFile = new QPushButton(readModelDialog);
        selectFile->setObjectName(QStringLiteral("selectFile"));
        selectFile->setGeometry(QRect(290, 60, 93, 28));

        retranslateUi(readModelDialog);

        QMetaObject::connectSlotsByName(readModelDialog);
    } // setupUi

    void retranslateUi(QDialog *readModelDialog)
    {
        readModelDialog->setWindowTitle(QApplication::translate("readModelDialog", "readModelDialog", Q_NULLPTR));
        conform->setText(QApplication::translate("readModelDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        cancle->setText(QApplication::translate("readModelDialog", "\345\217\226\346\266\210", Q_NULLPTR));
        label_2->setText(QApplication::translate("readModelDialog", "\346\250\241\345\236\213\345\261\236\346\200\247", Q_NULLPTR));
        modelType->clear();
        modelType->insertItems(0, QStringList()
         << QApplication::translate("readModelDialog", "point3d", Q_NULLPTR)
         << QApplication::translate("readModelDialog", "point4d", Q_NULLPTR)
         << QApplication::translate("readModelDialog", "bezierLine", Q_NULLPTR)
         << QApplication::translate("readModelDialog", "nurbusLine", Q_NULLPTR)
         << QApplication::translate("readModelDialog", "nurbusSurface", Q_NULLPTR)
         << QApplication::translate("readModelDialog", "nurbusVol", Q_NULLPTR)
        );
        label->setText(QApplication::translate("readModelDialog", "\346\226\207\344\273\266\350\267\257\345\276\204", Q_NULLPTR));
        selectFile->setText(QApplication::translate("readModelDialog", "\351\200\211\346\213\251\346\226\207\344\273\266", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class readModelDialog: public Ui_readModelDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_READMODELDIALOG_H
