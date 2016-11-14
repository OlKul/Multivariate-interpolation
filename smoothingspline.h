#ifndef SMOOTHINGSPLINE_H
#define SMOOTHINGSPLINE_H

#include <QFile>


class SmoothingSpline // имя класса
{
private: // спецификатор доступа private
    QFile Input;
    QFile Output;
    int l; // количество точек
    double
    *x = new double [l],// координаты x
    *y = new double [l],// координаты y
    *z = new double [l],// координаты z
           p; // погрешность точек
public: // спецификатор доступа public
    double g(double , double , double , double , double , double , double , double );
    void setSP2D(QFile, double); // чтение координат
    void setSP3D(QFile, double);
    void setP(double); // установка даты в формате дд.мм.гг
    void getX(); // отобразить значения новых X
    void getY(); // отобразить значения новых Y
    void getZ(); // отобразить значения новых Z

};
#endif // SMOOTHINGSPLINE_H
