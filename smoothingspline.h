#ifndef SMOOTHINGSPLINE_H
#define SMOOTHINGSPLINE_H


#include <QFile>
#include <QVector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

class SmoothingSpline // имя класса
{
private: // спецификатор доступа private
    QFile Input;
    QFile Output;
    int l; // количество точек
    QVector<double>
    x, // координаты x
    y, // координаты y
    z, // координаты z
    X_Out,
    Y_Out,
    Z_Out,
    m, // коэффициент одномерный случай
    h; // интервалы
    double p; // погрешность точек
public: // спецификатор доступа public
    double g(double , double , double , double , double , double , double , double );
    void setSP2D(QFile, double); // чтение координат
    void setSP3D(QFile, double);
    void setP(double); // установка даты в формате дд.мм.гг
    void getX(); // отобразить значения новых X
    void getY(); // отобразить значения новых Y
    void getZ(); // отобразить значения новых Z
    void buildSP2D();

};
#endif // SMOOTHINGSPLINE_H
