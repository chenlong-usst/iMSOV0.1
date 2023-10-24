/********************************************************************************
** Form generated from reading UI file 'occview.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OCCVIEW_H
#define UI_OCCVIEW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_OccView
{
public:

    void setupUi(QWidget *OccView)
    {
        if (OccView->objectName().isEmpty())
            OccView->setObjectName(QStringLiteral("OccView"));
        OccView->resize(1760, 1196);

        retranslateUi(OccView);

        QMetaObject::connectSlotsByName(OccView);
    } // setupUi

    void retranslateUi(QWidget *OccView)
    {
        OccView->setWindowTitle(QApplication::translate("OccView", "OccView", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class OccView: public Ui_OccView {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OCCVIEW_H
