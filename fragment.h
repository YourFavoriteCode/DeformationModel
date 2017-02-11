#ifndef __FRAGMENT_H 
#define __FRAGMENT_H

#include "MathCore.h"

namespace model
{

	class Fragment
	{
	public:
		Tensor d;						//Тензор деформации скорости
		Tensor w;						//Тензор вихря
		Tensor d_in;					//Тензор неупругой составляющей деформации
		Tensor o;						//Ориентационный тензор
		Tensor om;						//Тензор спина решётки
		Tensor sgm;						//Тензор напряжений
		Tensor dsgm;					//Тензор скоростей напряжений
		Tensor e;						//Тензор деформаций
		Tensor4 p;						//Тензор упругих свойств
		
		int SS_count;					//Кол-во систем скольжения
		SlipSystem *SS;					//Системы скольжения

		int material;					//Тип материала
		double size;					//Линейный размер фрагмента
		double volume;					//Объём фрагмента с учётом срезов
		double stress;					//Интенсивность напряжений
		double strain;					//Интенсивность деформаций
		
		double norm;//СОТРИ ЭТО ДЕРЬМО
		double dmc;	//И это наверное

		Fragment *surrounds;			//Ссылки на граничащие фрагменты
		Vector *normals;				//Вектора нормали к граничащим фрагментам
		Vector moment;				
		int *contact;					

		int position;					//Порядковый номер данного элемента в поликристалле
		bool crystall_center;			//Является ли фрагмент центром кристаллизации
		/*
		contact - массив с информацией о том, как соприкасаются фрагменты
		возможные значения:
		|-1	| Контакт ещё не задан						|
		| 0	| Нет контакта								|
		| 1	| Контакт на полной площади фасетки (100 %)	|
		| 2	| Контакт на ребре (10 %)					|
		| 3	| Контакт на вершине (5 %)					|
		*/

		bool isRotate;					//Вращается ли решётка на текущем шаге
		double rot_speed;				//Скорость вращения решётки
		double sum_angle;				//Накопленный угол поворота решётки
		double rot_energy;				//Энергия ротаций фрагмента
		
		double rot_Mc;					//Начальный критический момент
		double rot_H;					//Коэффициент при необратимой составляющей поворота
		double rot_A;					//Коэффициент при обратимой составляющей поворота
		double rot_L;					//Коэффициент при скорости поверхностных моментов

		void setMaterialParams(int);	//Задание материальных параметров
		void Orientate(double, double,
			double, double);			//Задание начальной ориентации
		void NDScalc();					//Вычисление НДС фрагмента
		
		double DisorientMeasure(int);	//Возвращает меру разориентации с выбранным соседом

		Fragment();
		~Fragment();


	private:

	};


}

#endif