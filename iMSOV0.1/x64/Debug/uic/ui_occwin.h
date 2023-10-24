/********************************************************************************
** Form generated from reading UI file 'occwin.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OCCWIN_H
#define UI_OCCWIN_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_OccWin
{
public:

    void setupUi(QDialog *OccWin)
    {
        if (OccWin->objectName().isEmpty())
            OccWin->setObjectName(QStringLiteral("OccWin"));
        OccWin->resize(1644, 1111);

        retranslateUi(OccWin);

        QMetaObject::connectSlotsByName(OccWin);
    } // setupUi

    void retranslateUi(QDialog *OccWin)
    {
        OccWin->setWindowTitle(QApplication::translate("OccWin", "OccWin", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class OccWin: public Ui_OccWin {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OCCWIN_H
