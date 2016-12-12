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



LIBS += -LC:\Qt\Projects\Multivariate_interpolation\GnuWin32\lib/ -lgsl

INCLUDEPATH += C:\Qt\Projects\Multivariate_interpolation\GnuWin32\include
DEPENDPATH += C:\Qt\Projects\Multivariate_interpolation\GnuWin32\include

LIBS += -LC:\Qt\Projects\Multivariate_interpolation\GnuWin32\lib/ -lgslcblas

INCLUDEPATH += C:\Qt\Projects\Multivariate_interpolation\GnuWin32\include
DEPENDPATH += C:\Qt\Projects\Multivariate_interpolation\GnuWin32\include
