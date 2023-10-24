/********************************************************************************
** Form generated from reading UI file 'SelectDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SELECTDIALOG_H
#define UI_SELECTDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>

QT_BEGIN_NAMESPACE

class Ui_SelectDialog
{
public:
    QRadioButton *select_whole;
    QRadioButton *select_point;
    QRadioButton *select_face;
    QRadioButton *select_lines;
    QPushButton *button_cancel;
    QPushButton *button_ok;
    QRadioButton *select_line;

    void setupUi(QDialog *SelectDialog)
    {
        if (SelectDialog->objectName().isEmpty())
            SelectDialog->setObjectName(QStringLiteral("SelectDialog"));
        SelectDialog->resize(270, 300);
        select_whole = new QRadioButton(SelectDialog);
        select_whole->setObjectName(QStringLiteral("select_whole"));
        select_whole->setGeometry(QRect(50, 20, 151, 41));
        select_point = new QRadioButton(SelectDialog);
        select_point->setObjectName(QStringLiteral("select_point"));
        select_point->setGeometry(QRect(50, 60, 151, 41));
        select_face = new QRadioButton(SelectDialog);
        select_face->setObjectName(QStringLiteral("select_face"));
        select_face->setGeometry(QRect(50, 180, 151, 41));
        select_lines = new QRadioButton(SelectDialog);
        select_lines->setObjectName(QStringLiteral("select_lines"));
        select_lines->setGeometry(QRect(50, 140, 151, 41));
        button_cancel = new QPushButton(SelectDialog);
        button_cancel->setObjectName(QStringLiteral("button_cancel"));
        button_cancel->setGeometry(QRect(150, 230, 93, 28));
        button_ok = new QPushButton(SelectDialog);
        button_ok->setObjectName(QStringLiteral("button_ok"));
        button_ok->setGeometry(QRect(40, 230, 93, 28));
        select_line = new QRadioButton(SelectDialog);
        select_line->setObjectName(QStringLiteral("select_line"));
        select_line->setGeometry(QRect(50, 100, 151, 41));

        retranslateUi(SelectDialog);

        QMetaObject::connectSlotsByName(SelectDialog);
    } // setupUi

    void retranslateUi(QDialog *SelectDialog)
    {
        SelectDialog->setWindowTitle(QApplication::translate("SelectDialog", "SelectDialog", Q_NULLPTR));
        select_whole->setText(QApplication::translate("SelectDialog", "\351\200\211\346\213\251\346\225\264\344\275\223", Q_NULLPTR));
        select_point->setText(QApplication::translate("SelectDialog", "\351\200\211\346\213\251\347\202\271", Q_NULLPTR));
        select_face->setText(QApplication::translate("SelectDialog", "\351\200\211\346\213\251\351\235\242", Q_NULLPTR));
        select_lines->setText(QApplication::translate("SelectDialog", "\351\200\211\346\213\251\350\276\271\346\241\206", Q_NULLPTR));
        button_cancel->setText(QApplication::translate("SelectDialog", "\345\217\226\346\266\210", Q_NULLPTR));
        button_ok->setText(QApplication::translate("SelectDialog", "\347\241\256\345\256\232", Q_NULLPTR));
        select_line->setText(QApplication::translate("SelectDialog", "\351\200\211\346\213\251\347\272\277", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class SelectDialog: public Ui_SelectDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SELECTDIALOG_H
