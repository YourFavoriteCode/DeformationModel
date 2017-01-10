#ifndef __POLYCRYST_H 
#define __POLYCRYST_H

#include "fragment.h"

namespace model
{
	class Polycrystall
	{
	public:
		Fragment ***C;				//Массив элементов поликристалла
		int fragm_count;			//Кол-во фрагментов на ребре
		int total_fragm_count;		//Общее кол-во фрагментов

		Tensor D;					//Тензор деформации скорости
		Tensor W;					//Тензор вихря
		Tensor D_in;				//Тензор неупругой части деформации скорости
		Tensor E;					//Тензор деформаций
		Tensor dSgm;				//Тензор скоростей напряжений
		Tensor Sgm;					//Тензор напряжений
		Tensor4 P;					//Усреднённый тензор упругих констант

		int cycle;
		int CURR_STEP;				//Текущий шаг интегрирования
		int PLOT_STEP;				//Шаг сохранения графиков
		int POLUS_STEP;				//Шаг сохранения ПФ
		int PROC_STEP;				//Шаг отображения прогресса
		int DEBUG_STEP;				//Шаг записи отладочных данных
		int proc_period;			//Период обновления процента выполнения

		double Strain;				//Интенсивность деформаций
		double Stress;				//Интенсивность напряжений

		int file_count;				//Кол-во отладочных файлов
		std::ofstream *dbgstream;	//Массив файлов для отладки
		std::ofstream *Datastream;	//Массив файлов с кривыми НДС
		std::ofstream *TestStream;	//Массив временных (тестовых) файлов

		Polycrystall();
		~Polycrystall();

		double tension_component;	//Вытягивающая компонента в одноосном нагружении
		double final_stress;		//Значение, до которого разгружать
		double lam;					//Коэффициент в разгрузке
		double addition_strain;		//Добавочный коэффициент для продолжения циклического нагружения

		void Init(int);				//Выделение памяти под зёрна
		void MakeStruct();			//Распределение нормалей и фасеток всех фрагментов
		void setParams();			//Распределение параметров фрагментов
		void Deformate();			//Деформирование поликристалла
		
		void Fragmentate();			//Процедура фрагментации
		void SavePoleFig();			//Сохранение ПФ
		void SaveDbgInfo();			//Сохранение отладочных данных
		void OpenFiles();			//Открытие всех файлов для записи
		void CloseFiles();			//Закрытие всех файлов для записи

		int getNumArr(int, int, int);	//Возвращает уникальный номер фрагмента в общей структуре
		void getMultyNumArr(int, int&, int&, int&);	//Возвращает позицию фрагмента в трёхмерном массиве 
	private:
		void Load(bool unload);		//Нагружение поликристалла
	};
}

#endif __POLYCRYST_H