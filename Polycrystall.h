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

		Tensor D;					//Тензор деформации скорости макроуровня
		Tensor W;					//Тензор вихря макроуровня
		Tensor D_in;				//Тензор неупругой части деформации макроуровня
		Tensor E;					//Тензор макродеформаций
		Tensor dSgm;				//Тензор скоростей макронапряжений
		Tensor Sgm;					//Тензор макронапряжений
		Tensor4 P;					//Усреднённый тензор упругих констант

		int cycle;
		int CURR_STEP;				//Текущий шаг интегрирования
		int PLOT_STEP;					//Шаг сохранения графиков
		int POLUS_STEP;					//Шаг сохранения ПФ
		int PROC_STEP;				//Шаг отображения прогресса
		int DEBUG_STEP;				//Шаг записи отладочных данных
		int proc_period;			//Период обновления процента выполнения

		double Strain;				//Интенсивность макродеформаций
		double Stress;				//Интенсивность макронапряжений

		int file_count;
		std::ofstream *dbgstream;
		std::ofstream *TestStream;
		std::ofstream *Datastream;
	
		Polycrystall();
		~Polycrystall();

		double tension_component;	//Вытягивающая компонента в одноосье
		double final_stress;		//Значение, до которого разгружать
		double lam;				//Коэффициент в разгрузке
		double addition_strain;	//Добавочный множитель для продолжения циклики

		void Init(int);				//Выделение памяти под зёрна
		void MakeStruct();			//Распределение нормалей и фасеток
		void setParams();
		void Deformate();
		
		void Fragmentate();
		void SavePoleFig();
		void SaveDbgInfo();
		void OpenFiles();
		void CloseFiles();

		int getNumArr(int, int, int);
		void getMultyNumArr(int, int&, int&, int&);
	private:
		void Load(bool unload);				//Нагружение
	};
}

#endif __POLYCRYST_H