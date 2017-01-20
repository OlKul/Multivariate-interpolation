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

int sampleCountX=5;
int sampleCountZ=5;
double xmin=5;
double xmax=5;
double zmin=5;
double zmax=5;
double ymin=5;
double ymax=5;

InterpolationGraph::InterpolationGraph(Q3DSurface *surface)
    : m_graph(surface)
{
    m_graph->setAxisX(new QValue3DAxis);
    m_graph->setAxisY(new QValue3DAxis);
    m_graph->setAxisZ(new QValue3DAxis);


    m_sqrtSinProxy = new QSurfaceDataProxy();
    m_sqrtSinSeries = new QSurface3DSeries(m_sqrtSinProxy);

    fillSqrtSinProxy();

}

InterpolationGraph::~InterpolationGraph()
{
    delete m_graph;
}

/*!
 * /brief Функций вычисления, построения сплайна
 * Функция считывает точки из файла, заполняет ими массивы, строит сплайн функции,
 * после чего строит их.
 */
void InterpolationGraph::fillSqrtSinProxy()
{

    qDebug() << "1234   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    QVector <double> xQt, yQt,zQt;
    QFile InputFile("C:\\Qt\\Projects\\Multivariate_interpolation\\Input.txt");
    int i;
    if (!InputFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "Can not open Input.txt\n";
        return;
    }
    qDebug() << "#1\n";
    QTextStream inXY(&InputFile);
    while (!inXY.atEnd()) {
        qDebug() << "#3\n";
        QString line = inXY.readLine();
        QStringList list = line.split(" ");
        qDebug() << "#3.1" << list.size() << " " << list << "\n";
        QString strX = list.at(0);
        QString strY = list.at(1);
        QString strZ = list.at(2);
        qDebug() << "#3.2\n";
        //xQt.append(strX.toDouble());
        //yQt.append(strY.toDouble());
        //zQt.append(strZ.toDouble());
        qDebug() << "#3.3\n";
    }
    qDebug() << "#2\n";
    xQt = {1, 2, 3};
    yQt = {1, 2, 3};
    zQt = {3.2, 2.5, 5.1, 4.4, 4.7, 3.6, 6.5, 5.8, 2.9};
    qDebug() << "#2.1\n";
    alglib::ae_int_t  xsizeAlglib=xQt.size();
    qDebug() << "#2.2\n";
    int xsize=xQt.size();
    alglib::ae_int_t  ysizeAlglib=yQt.size();
    int ysize=yQt.size();
    //alglib::ae_int_t  zsizeAlglib=zQt.size();
    int zsize=zQt.size();
    qDebug() << "#2.3\n";
    alglib::real_1d_array xAlglib,yAlglib,zAlglib;

    xmin=xQt[0],xmax=xQt[0],ymin=yQt[0],ymax=yQt[0], zmin=zQt[0],zmax=zQt[0];
    double x[xsize-1], y[xsize-1], z[zsize-1];
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
    for (i=0; i<zsize; i++)
    {
        if (zQt[i]<zmin) zmin=zQt[i];
        if (zQt[i]>zmax) zmax=zQt[i];
        z[i]=zQt[i];
    }
    zAlglib.setcontent(zsize,z);
    float stepX = 0.01;
    float stepZ = 0.01;

    spline2dinterpolant s;
    qDebug() << "#3\n";

        spline2dbuildbicubicv(xAlglib, xsizeAlglib, yAlglib, ysizeAlglib, zAlglib, 1, s);

    qDebug() << "#3.1\n";
    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(sampleCountZ);
    double x2,y2,z2;
    qDebug() << "#3.2\n";
    sampleCountX=(xmax - xmin) / 0.01;
    sampleCountZ=(ymax - ymin) / 0.01;
    for ( i=0 ; i < sampleCountZ  ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(sampleCountX);
        // Keep values within range bounds, since just adding step can cause minor drift due
        // to the rounding errors.
        y2 = ymin+i*0.01;
        float z3 = qMin(ymax, (i * stepZ + ymin));
        int index = 0;
        qDebug() << "#3.3\n";
        for (int j=0 ; j < sampleCountX  ; j++) {
            x2 = xmin+j*0.01;
            qDebug() << "#3.4\n";
            //z2 = spline2dcalc(s, x2, y2);
            z2= x2*x2+y2*y2;
            float x3 = qMin(xmax, (j * stepX + xmin));
            //float R = qSqrt(z3 * z3 + x3 * x3) + 0.01f;
            float y3 = spline2dcalc(s, x3, z3);
            //qDebug() << "#3.5 " << ((xmax - xmin) / 0.1) << " " << (*newRow).size() << " " << index << "\n";
            (*newRow)[index++].setPosition(QVector3D(x3, y3, z3));
            qDebug() << "#3.6\n";
        }
        *dataArray << newRow;
    }
    qDebug() << "#3.7\n" << xmax << ymax << zmax;
    m_sqrtSinProxy->resetArray(dataArray);

}

/// Улучшение интерфеса
void InterpolationGraph::enableSqrtSinModel(bool enable)
{
    if (enable) {

        m_sqrtSinSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_sqrtSinSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f X");
        m_graph->axisZ()->setLabelFormat("%.2f Y");
        m_graph->axisY()->setLabelFormat("%.2f Z");
        m_graph->axisX()->setRange(xmin, xmax);
        m_graph->axisY()->setRange(zmin, zmax);
        m_graph->axisZ()->setRange(ymin, ymax);
        m_graph->axisX()->setTitle("X");
        m_graph->axisY()->setTitle("Y");
        m_graph->axisZ()->setTitle("Z");
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);

        m_graph->addSeries(m_sqrtSinSeries);
        //!

        //! The example has four slider controls for adjusting
        //!  the min and max values for X and Z axis. When selecting
        //! the proxy these sliders are adjusted so that one step on the slider moves the range by one segment step:
        // Reset range sliders for Sqrt&Sin
        m_rangeMinX = xmin;
        m_rangeMinZ = ymin;
        m_stepX = (xmax - xmin) / double(sampleCountX - 1);
        m_stepZ = (ymax - ymin) / double(sampleCountZ - 1);
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

/// Установка диапазонов:
void InterpolationGraph::setAxisXRange(double min, double max)
{
    m_graph->axisX()->setRange(min, max);
}

void InterpolationGraph::setAxisZRange(double min, double max)
{
    m_graph->axisZ()->setRange(min, max);
}
//!

/// Выбор темы
void InterpolationGraph::changeTheme(int theme)
{
    m_graph->activeTheme()->setType(Q3DTheme::Theme(theme));
}

/// Первый градиент
void InterpolationGraph::setBlackToYellowGradient()
{

    QLinearGradient gr;
    gr.setColorAt(0.0, Qt::black);
    gr.setColorAt(0.33, Qt::blue);
    gr.setColorAt(0.67, Qt::red);
    gr.setColorAt(1.0, Qt::yellow);

    m_graph->seriesList().at(0)->setBaseGradient(gr);
    m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

}
/// Второй градиент
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
