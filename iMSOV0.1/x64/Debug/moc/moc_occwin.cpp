/****************************************************************************
** Meta object code from reading C++ file 'occwin.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../occwin.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'occwin.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_OccWin_t {
    QByteArrayData data[9];
    char stringdata0[75];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_OccWin_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_OccWin_t qt_meta_stringdata_OccWin = {
    {
QT_MOC_LITERAL(0, 0, 6), // "OccWin"
QT_MOC_LITERAL(1, 7, 10), // "trigerMenu"
QT_MOC_LITERAL(2, 18, 0), // ""
QT_MOC_LITERAL(3, 19, 8), // "QAction*"
QT_MOC_LITERAL(4, 28, 13), // "trigerToolBar"
QT_MOC_LITERAL(5, 42, 11), // "SelectModel"
QT_MOC_LITERAL(6, 54, 16), // "Standard_Integer"
QT_MOC_LITERAL(7, 71, 1), // "x"
QT_MOC_LITERAL(8, 73, 1) // "y"

    },
    "OccWin\0trigerMenu\0\0QAction*\0trigerToolBar\0"
    "SelectModel\0Standard_Integer\0x\0y"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_OccWin[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x08 /* Private */,
       4,    1,   32,    2, 0x08 /* Private */,
       5,    2,   35,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, 0x80000000 | 6, 0x80000000 | 6,    7,    8,

       0        // eod
};

void OccWin::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        OccWin *_t = static_cast<OccWin *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->trigerMenu((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 1: _t->trigerToolBar((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 2: _t->SelectModel((*reinterpret_cast< Standard_Integer(*)>(_a[1])),(*reinterpret_cast< Standard_Integer(*)>(_a[2]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QAction* >(); break;
            }
            break;
        case 1:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QAction* >(); break;
            }
            break;
        }
    }
}

const QMetaObject OccWin::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_OccWin.data,
      qt_meta_data_OccWin,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *OccWin::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *OccWin::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_OccWin.stringdata0))
        return static_cast<void*>(const_cast< OccWin*>(this));
    return QDialog::qt_metacast(_clname);
}

int OccWin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
