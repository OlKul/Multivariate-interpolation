#include "interpolationgraph.h"
#include "alglib/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "alglib/linalg.h"
#include "alglib/interpolation.h"

using namespace alglib;

#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/Q3DTheme>
#include <QtGui/QImage>
#include <QtCore/qmath.h>
#include <QtDebug>
#include <QTextStream>
#include <QFile>
using namespace QtDataVisualization;

const int sampleCountX = 50;
const int sampleCountZ = 50;
const double sampleMin = -8.0f;
const double sampleMax = 8.0f;

InterpolationGraph::InterpolationGraph(Q3DSurface *surface)
    : m_graph(surface)
{
    m_graph->setAxisX(new QValue3DAxis);
    m_graph->setAxisY(new QValue3DAxis);
    m_graph->setAxisZ(new QValue3DAxis);

    //! [0]
    m_sqrtSinProxy = new QSurfaceDataProxy();
    m_sqrtSinSeries = new QSurface3DSeries(m_sqrtSinProxy);
    //! [0]
    fillSqrtSinProxy();

}

InterpolationGraph::~InterpolationGraph()
{
    delete m_graph;
}

//! Строим массив значений:
void InterpolationGraph::fillSqrtSinProxy()
{
    qDebug() << "1234   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    QVector<double>
       xQt, // координаты x
       yQt, // координаты y
       zQt; // координаты z
    QFile InputFile("Input.txt");
    int i;
    if (!InputFile.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        QTextStream inXY(&InputFile);
        while (!inXY.atEnd()) {
            QString line = inXY.readLine();
            QStringList list = line.split("\t");
            QString strX = list.at(0);
            QString strY = list.at(1);
            QString strZ = list.at(2);
            xQt.append(strX.toDouble());
            yQt.append(strY.toDouble());
            zQt.append(strZ.toDouble());
            qDebug() << "1234   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            i++;
        }
        alglib::ae_int_t  xsizeAlglib=xQt.size();
        int xsize=xQt.size();
        alglib::ae_int_t  ysizeAlglib=yQt.size();
        int ysize=yQt.size();
        //alglib::ae_int_t  zsizeAlglib=zQt.size();
        int zsize=zQt.size();
        alglib::real_1d_array xAlglib,yAlglib,zAlglib;

        double xmin=xQt[0],xmax=xQt[0],ymin=yQt[0],ymax=yQt[0];
        double x[xsize-1], y[ysize-1], z[zsize-1];
    for (i=0; i<xsize; i++)
    {
        if (xQt[i]<xmin) xmin=xQt[i];
        if (xQt[i]>xmax) xmax=xQt[i];
        x[i]=xQt[i];
    }
    xAlglib.setcontent(xsize,x);

    for (i=0; i<ysize; i++)
    {
        if (yQt[i]<ymin) ymin=yQt[i];
        if (yQt[i]>ymax) ymax=yQt[i];
        y[i]=yQt[i];
    }

    yAlglib.setcontent(ysize,y);
    for (i=0; i<xsize; i++)
    {
//      if (z[i]<zmin) zmin=z[i];
//      if (z[i]>zmax) zmax=z[i];
        z[i]=zQt[i];
    }
    zAlglib.setcontent(zsize,z);
    spline2dinterpolant s;

    spline2dbuildbicubicv(xAlglib, xsizeAlglib, yAlglib, ysizeAlglib, zAlglib, 1, s);
    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve((ymax - ymin) / 0.01);
    double x2,y2,z2;
    for ( y2=ymin ; y2 <= ymax  ; y2=y2+0.01) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow((xmax - xmin) / 0.01);
        // Keep values within range bounds, since just adding step can cause minor drift due
        // to the rounding errors.

        int index = 0;
        for (x2=xmin ; i <= xmax  ; x2=x2+0.01) {
            z2 = spline2dcalc(s, x2, y2);
            (*newRow)[index++].setPosition(QVector3D(x2, y2, z2));
        }
        *dataArray << newRow;
    }

    m_sqrtSinProxy->resetArray(dataArray);
}
//!

void InterpolationGraph::enableSqrtSinModel(bool enable)
{
    if (enable) {
        //! Сначала мы устанавливаем декоративные проблемы, любят, включают сетку
        //!  для поверхности и выбирают режим закрашивания плоскостями. Следующие строки
        //! определяют формат метки оси и диапазоны значений. Автоматическое вращение метки
        //!  собирается улучшить удобочитаемость метки в низких ракурсах.
        //!  Наконец мы удостоверяемся, что корректный ряд добавлен к графику:
        m_sqrtSinSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_sqrtSinSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(sampleMin, sampleMax);
        m_graph->axisY()->setRange(0.0f, 2.0f);
        m_graph->axisZ()->setRange(sampleMin, sampleMax);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);

        m_graph->addSeries(m_sqrtSinSeries);
        //!

        //! The example has four slider controls for adjusting
        //!  the min and max values for X and Z axis. When selecting
        //! the proxy these sliders are adjusted so that one step on the slider moves the range by one segment step:
        // Reset range sliders for Sqrt&Sin
        m_rangeMinX = sampleMin;
        m_rangeMinZ = sampleMin;
        m_stepX = (sampleMax - sampleMin) / double(sampleCountX - 1);
        m_stepZ = (sampleMax - sampleMin) / double(sampleCountZ - 1);
        m_axisMinSliderX->setMaximum(sampleCountX - 2);
        m_axisMinSliderX->setValue(0);
        m_axisMaxSliderX->setMaximum(sampleCountX - 1);
        m_axisMaxSliderX->setValue(sampleCountX - 1);
        m_axisMinSliderZ->setMaximum(sampleCountZ - 2);
        m_axisMinSliderZ->setValue(0);
        m_axisMaxSliderZ->setMaximum(sampleCountZ - 1);
        m_axisMaxSliderZ->setValue(sampleCountZ - 1);
        //!
    }
}



void InterpolationGraph::adjustXMin(int min)
{
    double minX = m_stepX * double(min) + m_rangeMinX;

    int max = m_axisMaxSliderX->value();
    if (min >= max) {
        max = min + 1;
        m_axisMaxSliderX->setValue(max);
    }
    double maxX = m_stepX * max + m_rangeMinX;

    setAxisXRange(minX, maxX);
}

void InterpolationGraph::adjustXMax(int max)
{
    double maxX = m_stepX * double(max) + m_rangeMinX;

    int min = m_axisMinSliderX->value();
    if (max <= min) {
        min = max - 1;
        m_axisMinSliderX->setValue(min);
    }
    double minX = m_stepX * min + m_rangeMinX;

    setAxisXRange(minX, maxX);
}

void InterpolationGraph::adjustZMin(int min)
{
    double minZ = m_stepZ * double(min) + m_rangeMinZ;

    int max = m_axisMaxSliderZ->value();
    if (min >= max) {
        max = min + 1;
        m_axisMaxSliderZ->setValue(max);
    }
    double maxZ = m_stepZ * max + m_rangeMinZ;

    setAxisZRange(minZ, maxZ);
}

void InterpolationGraph::adjustZMax(int max)
{
    double maxX = m_stepZ * double(max) + m_rangeMinZ;

    int min = m_axisMinSliderZ->value();
    if (max <= min) {
        min = max - 1;
        m_axisMinSliderZ->setValue(min);
    }
    double minX = m_stepZ * min + m_rangeMinZ;

    setAxisZRange(minX, maxX);
}

//! Установка диапазонов:
void InterpolationGraph::setAxisXRange(double min, double max)
{
    m_graph->axisX()->setRange(min, max);
}

void InterpolationGraph::setAxisZRange(double min, double max)
{
    m_graph->axisZ()->setRange(min, max);
}
//!

//! Выбор темы:
void InterpolationGraph::changeTheme(int theme)
{
    m_graph->activeTheme()->setType(Q3DTheme::Theme(theme));
}
//!

void InterpolationGraph::setBlackToYellowGradient()
{
    //! Один из градиентов:
    QLinearGradient gr;
    gr.setColorAt(0.0, Qt::black);
    gr.setColorAt(0.33, Qt::blue);
    gr.setColorAt(0.67, Qt::red);
    gr.setColorAt(1.0, Qt::yellow);

    m_graph->seriesList().at(0)->setBaseGradient(gr);
    m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
    //!
}

void InterpolationGraph::setGreenToRedGradient()
{
    QLinearGradient gr;
    gr.setColorAt(0.0, Qt::darkGreen);
    gr.setColorAt(0.5, Qt::yellow);
    gr.setColorAt(0.8, Qt::red);
    gr.setColorAt(1.0, Qt::darkRed);

    m_graph->seriesList().at(0)->setBaseGradient(gr);
    m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
}
