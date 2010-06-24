/****************************************************************************
** Meta object code from reading C++ file 'EdgeRatioThresholdDlg.h'
**
** Created: Thu Jun 24 12:50:26 2010
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "EdgeRatioThresholdDlg.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'EdgeRatioThresholdDlg.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_EdgeRatioThresholdDlg[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // slots: signature, parameters, type, tag, flags
      25,   23,   22,   22, 0x08,
      51,   23,   22,   22, 0x08,
      78,   23,   22,   22, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_EdgeRatioThresholdDlg[] = {
    "EdgeRatioThresholdDlg\0\0t\0"
    "setSmallThreshold(double)\0"
    "setMeidanThreshold(double)\0"
    "setLargeThreshold(double)\0"
};

const QMetaObject EdgeRatioThresholdDlg::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_EdgeRatioThresholdDlg,
      qt_meta_data_EdgeRatioThresholdDlg, 0 }
};

const QMetaObject *EdgeRatioThresholdDlg::metaObject() const
{
    return &staticMetaObject;
}

void *EdgeRatioThresholdDlg::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_EdgeRatioThresholdDlg))
        return static_cast<void*>(const_cast< EdgeRatioThresholdDlg*>(this));
    return QDialog::qt_metacast(_clname);
}

int EdgeRatioThresholdDlg::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: setSmallThreshold((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: setMeidanThreshold((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 2: setLargeThreshold((*reinterpret_cast< double(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 3;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
