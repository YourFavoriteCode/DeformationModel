#ifndef __DISTRIB_H 
#define __DISTRIB_H

namespace model
{
	/*
	Различные распределения случайной величины
	Для двухпараметрических - первый параметр - 
	мат.ожидание, второй - дисперсия
	*/
	double UniformDistrib(double, double);		//Равномерное распределение
	double NormalDistrib(double m, double d);	//Нормальное распределение
	double LogNormalDistrib(double m, double d);//Логнормальное распределение
	double ExpDistrib(double l);				//Показательное распределение
}

#endif __DISTRIB_H