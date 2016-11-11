#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
double g(double m0, double m1, double f0, double f1, double u, double x2, double x3, double w);
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

double g(double m0, double m1, double f0, double f1, double u, double x2, double x3, double w)
{
    return (m0*((x3-w)*(x3-w)*(x3-w)))/(6*u)+m1*((w-x2)*(w-x2)*(w-x2))/(6*u)+(f0-m0*u*u/6)*(x3-w)/u+(f1-m1*u*u/6)*(w-x2)/u;
}

void MainWindow::on_pushButton_clicked()
{
        int n=6;
        int l=7;
        double x[l], y[l], p[l];
        QVector<double> xv(l-l),yv(l-l);
        double t=0.01;
        x[0]=0;
        x[5]=6;
        int N = (x[n]-x[0])/t+2;
        QVector<double> X_Out(N-N), Y_Out(N-N);
        QFile Input_X ("Input_X.txt");
        QFile Input_Y ("Input_Y.txt");
        QFile Input_P ("Input_P.txt");
        int i=0;
        if (!Input_X.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        QTextStream inX(&Input_X);
        while (!inX.atEnd()) {
            QString line = inX.readLine();
            x[i] = line.toDouble();
            //qDebug() << x[i];
            i++;
        }


        i=0;
        if (!Input_Y.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        QTextStream inY(&Input_Y);
        while (!inY.atEnd()) {
            QString line = inY.readLine();
            y[i] = line.toDouble();
            //qDebug() << y[i];
            i++;
        }

        i=0;
        if (!Input_P.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        QTextStream inP(&Input_P);
        while (!inP.atEnd()) {
            QString line = inP.readLine();
            p[i] = line.toDouble();
            //qDebug() << p[i];
            i++;
        }
        //Spline(x, y, p);
        int k = l;
        double h[k-1];

        for (i=0; i<=k-2; i++)
        {
            h[i]=x[i+1]-x[i];
        }
        gsl_vector  *X = gsl_vector_alloc (k);

        for (i=0; i<=k-1;i++)
            {
                gsl_vector_set (X, i, x[i]);
                xv.push_back(x[i]);
            }
        //gsl_vector_fprintf (stdout, X, "%g");

        /*!
         * \brief Задаем Y
         */
        gsl_vector  *Y = gsl_vector_alloc (k);
        for (i=0; i<=k-1;i++)
            {
                gsl_vector_set (Y, i, y[i]);
                yv.push_back(y[i]);
            }
        //gsl_vector_fprintf (stdout, Y, "%g");
        /*!
         * \brief Составляем матрицу А
         */
        gsl_matrix *A = gsl_matrix_alloc (k-2,k-2);
        gsl_matrix_set_zero (A);
        int j;
        for (i = 0; i <= k-3; i++)
        {
            for (j = 0; j <= k-3; j++)
            {
                if (j==i-1) gsl_matrix_set (A, i, j, h[i]/6);

                if (i==j) gsl_matrix_set (A, i, j, (h[i]+h[i+1])/3);
                if (j==i+1) gsl_matrix_set (A, i, j, h[j]/6);

            }
        }
        //gsl_matrix_fprintf (stdout, A, "%g");
        /*!
         * \brief Составляем матрицу H
         */
        gsl_matrix *H = gsl_matrix_alloc (n-1,n+1);
        gsl_matrix_set_zero (H);
        for (i = 0; i <= n-2; i++)
        {
            for (j = 0; j <= n; j++)
            {
                if (j==i) gsl_matrix_set (H, i, j, 1/h[i]);
                if (j==i+1) gsl_matrix_set (H, i, j, -(1/h[i]+1/h[j]));
                if (j==i+2) gsl_matrix_set (H, i, j, 1/h[i+1]);
            }
        }
        //gsl_matrix_fprintf (stdout, H, "%g");
        /*!
         * \brief Составляем матрицу P
         */
        gsl_matrix *P = gsl_matrix_alloc (l,l);
        gsl_matrix_set_zero (P);
        for (i = 0; i <=l-1; i++)
        {
            for (j = 0; j <= l-1; j++)
            {

                if (i==j) gsl_matrix_set (P, i, j, 1/p[i]);

            }

        }
        //gsl_matrix_fprintf (stdout, P, "%g");

        /*!
         * \brief Вычисляем матрицу H*(H1)
         */
        gsl_matrix *H1 = gsl_matrix_alloc (n+1,n-1);
        gsl_matrix_transpose_memcpy (H1,H);

        //gsl_matrix_fprintf (stdout, H1, "%g");
        /*!
         * \brief Решаем пятидиагоняльную систему (À+HPH1)m=HY и нахождим вектор M
         */
        gsl_matrix *HP = gsl_matrix_alloc (n-1,n+1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, H, P, 0.0, HP);
        //gsl_matrix_fprintf (stdout, HP, "%g");
        gsl_matrix *HPH1 = gsl_matrix_alloc (n-1,n-1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, HP, H1, 0.0, HPH1);
        //gsl_matrix_fprintf (stdout, HPH1, "%g");
        gsl_matrix_add (A, HPH1);
        //gsl_matrix_fprintf (stdout, A, "%g");
        gsl_vector *HY = gsl_vector_alloc (n-1);
        gsl_blas_dgemv (CblasNoTrans, 1.0, H, Y, 0.0, HY);
        //gsl_vector_fprintf (stdout, HY, "%g");
        gsl_vector *m = gsl_vector_alloc (n-1);
        int s;
        gsl_permutation *per = gsl_permutation_alloc (n-1);
        gsl_linalg_LU_decomp (A, per, &s);
        gsl_linalg_LU_solve (A, per, HY, m);;
        //gsl_vector_fprintf (stdout, m, "%g");
        /*!
         * brief Теперьнеобходимо найти вектор сеточных значений- YS
         * Для этого решаем систему YS=Y-PH1m
         */
        gsl_matrix *PH1 = gsl_matrix_alloc (n+1,n-1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, P, H1, 0.0, PH1);
        gsl_vector *PH1m = gsl_vector_alloc (n+1);
        gsl_blas_dgemv (CblasNoTrans, 1.0, PH1, m, 0.0, PH1m);
        gsl_vector *YS = gsl_vector_alloc (n+1);
        gsl_vector_sub (Y, PH1m);
        YS = Y;
        //gsl_vector_fprintf (stdout, YS, "%g");

        gsl_vector *M = gsl_vector_alloc (n+1);
        gsl_vector_set (M, 0,  0);
        gsl_vector_set (M, n,0);
        for (i=1; i<=k-2; i++)
        {
            gsl_vector_set (M, i, gsl_vector_get (m, i-1));
        }



        /*!
         * \brief Далее вычисляем значения функции g(x)
         *  В точках на интервале [x1,xn] на расстоянии t друг от друга
         */



        double x2, x3, f0,f1,m0,m1;
        i=0;
        double w;
        int q;
        FILE *Output_X = fopen("Output_X.txt", "w");
        FILE *Output_Y = fopen("Output_Y.txt", "w");
        char str1[150];
        char str2[150];
        for (q=0; q<=k-2; q++)
        {
            for (w = gsl_vector_get(X,q); w <= gsl_vector_get(X,q+1); w=w+t)
            {
                m0= gsl_vector_get (M,q);
                m1=gsl_vector_get(M,q+1);
                f0=gsl_vector_get (YS,q);
                f1=gsl_vector_get(YS,q+1);
                x2=gsl_vector_get(X,q);
                x3=gsl_vector_get(X,q+1);
                X_Out.push_back(w);
                Y_Out.push_back(g( m0, m1,f0, f1, h[q], x2, x3,w));
                sprintf(str1, "%f", X_Out[i]);
                fputs(str1, Output_X);
                fputs("\n", Output_X);
                sprintf(str2, "%f", Y_Out[i]);
                fputs(str2, Output_Y);
                fputs("\n", Output_Y);
                i++;
             }
        }


        fclose(Output_X);
        fclose(Output_Y);



        ui->widget->clearGraphs();
        ui->widget->addGraph();
        ui->widget->graph(0)->setPen(QPen(Qt::blue));
        //ui->widget->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20)));
        ui->widget->addGraph();
        ui->widget->graph(1)->setPen(QPen(Qt::red));
        ui->widget->xAxis2->setVisible(true);
        ui->widget->xAxis2->setTickLabels(false);
        ui->widget->yAxis2->setVisible(true);
        ui->widget->yAxis2->setTickLabels(false);
        connect(ui->widget->xAxis, SIGNAL(rangeChanged(QCPRange)), ui->widget->xAxis2, SLOT(setRange(QCPRange)));
        connect(ui->widget->yAxis, SIGNAL(rangeChanged(QCPRange)), ui->widget->yAxis2, SLOT(setRange(QCPRange)));
        ui->widget->graph(0)->setData(X_Out, Y_Out);
        ui->widget->graph(1)->setData(xv, yv);
        ui->widget->graph(1)->setLineStyle(QCPGraph::lsNone);//убираем линии
        ui->widget->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
        ui->widget->graph(0)->rescaleAxes();
        ui->widget->graph(1)->rescaleAxes(true);

        ui->widget->xAxis->setLabel("X");
        ui->widget->yAxis->setLabel("Y");

        ui->widget->replot();

}


void MainWindow::on_pushButton_2_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this,tr("Choose a filename to save under"),QString(),tr("PNG(*.png);;JPG(*.jpg);;PDF(*.pdf);;BMP(*.bmp);;All Files(*)"));
       if(!fileName.isEmpty())
       {
           if(fileName.endsWith(".png", Qt::CaseInsensitive))
           {
               ui->widget->savePng(fileName);
           }
           else if(fileName.endsWith(".jpg", Qt::CaseInsensitive))
           {
               ui->widget->saveJpg(fileName);
           }
           else if(fileName.endsWith(".pdf", Qt::CaseInsensitive))
           {
               ui->widget->savePdf(fileName);
           }
           else if(fileName.endsWith(".bmp", Qt::CaseInsensitive))
           {
               ui->widget->saveBmp(fileName);
           }
           else
           {
               fileName += ".png";
               ui->widget->savePng(fileName);
           }
       }
}



void MainWindow::on_pushButton_3_clicked()
{

        qApp->quit();

}

void MainWindow::on_pushButton_4_clicked()
{
    ui->widget->clearGraphs();
}

void MainWindow::on_pushButton_5_clicked()
{
        double a = -1;
        double b =  1;
        double h = 0.01;

        int N=(b-a)/h + 2;
        QVector<double> x(N), y(N);

        int i=0;
        for (double X=a; X<=b; X+=h)
        {
            x[i] = X;
            y[i] = sin(X);
            i++;
        }

        ui->widget->clearGraphs();
        ui->widget->addGraph();
        ui->widget->graph(0)->setData(x, y);
        ui->widget->xAxis->setLabel("x");
        ui->widget->yAxis->setLabel("y");
        ui->widget->xAxis->setRange(a, b);
        double minY = y[0], maxY = y[0];
        for (int i=1; i<N; i++)
        {
            if (y[i]<minY) minY = y[i];
            if (y[i]>maxY) maxY = y[i];
        }
        ui->widget->yAxis->setRange(minY, maxY);
        ui->widget->replot();
}

void MainWindow::on_pushButton_6_clicked()
{


}
