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
#include <QFileDialog>
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


    m_BicubicSplineProxy = new QSurfaceDataProxy();
    m_BicubicSplineSeries = new QSurface3DSeries(m_BicubicSplineProxy);

    fillBicubicSplineProxy();

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
void InterpolationGraph::fillBicubicSplineProxy()
{


    QVector <double> xQt, yQt,zQt;
    QString Input_X,Input_Y,Input_Z;
    Input_X = QFileDialog::getOpenFileName(0,
        tr("Open Input data of X"), "", tr("Text Files (*.txt)"));
    Input_Y = QFileDialog::getOpenFileName(0,
        tr("Open Input data of Y"), "", tr("Text Files (*.txt)"));
    Input_Z = QFileDialog::getOpenFileName(0,
        tr("Open Input data of Z"), "", tr("Text Files (*.txt)"));
    QFile InputFile_X(Input_X),
          InputFile_Y(Input_Y),
          InputFile_Z(Input_Z);
    int i;
    if (!InputFile_X.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "Can not open Input_X.txt\n";
        return;
    }
    if (!InputFile_Y.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "Can not open Input_Y.txt\n";
        return;
    }
    if (!InputFile_Z.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "Can not open Input_Z.txt\n";
        return;
    }
    qDebug() << "#1\n";
    QTextStream inX(&InputFile_X);
    QTextStream inY(&InputFile_Y);
    QTextStream inZ(&InputFile_Z);
    while (!inX.atEnd()) {
        qDebug() << "#3\n";
        QString line = inX.readLine();
        qDebug() << "#3.1" << line << "\n";
        qDebug() << "#3.2\n";
        xQt.append(line.toDouble());
        qDebug() << "#3.3\n";
    }
    while (!inY.atEnd()) {
        qDebug() << "#3\n";
        QString line = inY.readLine();
        qDebug() << "#3.1" << line << "\n";
        qDebug() << "#3.2\n";
        yQt.append(line.toDouble());
        qDebug() << "#3.3\n";
    }
    while (!inZ.atEnd()) {
        qDebug() << "#3\n";
        QString line = inZ.readLine();
        qDebug() << "#3.1" << line << "\n";
        qDebug() << "#3.2\n";
        zQt.append(line.toDouble());
        qDebug() << "#3.3\n";
    }

    qDebug() << "#2\n";

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

    qDebug() << "#3.2\n";
    sampleCountX=(xmax - xmin) / 0.01;
    sampleCountZ=(ymax - ymin) / 0.01;
    for ( i=0 ; i < sampleCountZ  ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(sampleCountX);


        float z3 = qMin(ymax, (i * stepZ + ymin));
        int index = 0;
        qDebug() << "#3.3\n";
        for (int j=0 ; j < sampleCountX  ; j++) {

            qDebug() << "#3.4\n";

            float x3 = qMin(xmax, (j * stepX + xmin));

            float y3 = spline2dcalc(s, x3, z3);
            //qDebug() << "#3.5 " << ((xmax - xmin) / 0.1) << " " << (*newRow).size() << " " << index << "\n";
            (*newRow)[index++].setPosition(QVector3D(x3, y3, z3));
            qDebug() << "#3.6\n";
        }
        *dataArray << newRow;
    }
    qDebug() << "#3.7\n" << xmax << ymax << zmax;
    m_BicubicSplineProxy->resetArray(dataArray);

}

/// Улучшение интерфеса
void InterpolationGraph::enableBicubicSplineModel(bool enable)
{
    if (enable) {

        m_BicubicSplineSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_BicubicSplineSeries->setFlatShadingEnabled(true);

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

        m_graph->addSeries(m_BicubicSplineSeries);
        //!

        //! The example has four slider controls for adjusting
        //!  the min and max values for X and Z axis. When selecting
        //! the proxy these sliders are adjusted so that one step on the slider moves the range by one segment step:

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
