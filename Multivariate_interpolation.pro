QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += datavisualization
TARGET = Multivariate-interpolation
TEMPLATE = app


SOURCES += main.cpp\
        interpolationgraph.cpp \
    smoothingspline.cpp \
    spline.cpp

HEADERS  += interpolationgraph.h \
    smoothingspline.h \
    spline.h



unix|win32: LIBS += -L$$PWD/../GnuWin32/lib/ -lgsl

INCLUDEPATH += $$PWD/../GnuWin32/include
DEPENDPATH += $$PWD/../GnuWin32/include

unix|win32: LIBS += -L$$PWD/../GnuWin32/lib/ -lgslcblas

INCLUDEPATH += $$PWD/../GnuWin32/include
DEPENDPATH += $$PWD/../GnuWin32/include
