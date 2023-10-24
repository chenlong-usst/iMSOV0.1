/********************************************************************************
** Form generated from reading UI file 'OCCMainWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.9.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OCCMAINWINDOW_H
#define UI_OCCMAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_OCCMainWindowClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *OCCMainWindowClass)
    {
        if (OCCMainWindowClass->objectName().isEmpty())
            OCCMainWindowClass->setObjectName(QStringLiteral("OCCMainWindowClass"));
        OCCMainWindowClass->resize(600, 400);
        menuBar = new QMenuBar(OCCMainWindowClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        OCCMainWindowClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(OCCMainWindowClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        OCCMainWindowClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(OCCMainWindowClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        OCCMainWindowClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(OCCMainWindowClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        OCCMainWindowClass->setStatusBar(statusBar);

        retranslateUi(OCCMainWindowClass);

        QMetaObject::connectSlotsByName(OCCMainWindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *OCCMainWindowClass)
    {
        OCCMainWindowClass->setWindowTitle(QApplication::translate("OCCMainWindowClass", "OCCMainWindow", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class OCCMainWindowClass: public Ui_OCCMainWindowClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OCCMAINWINDOW_H
