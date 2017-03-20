#ifndef __MATHCORE_H 
#define __MATHCORE_H

namespace model
{
	/*
	*Описание классов "тензор", "вектор", "система скольжения"
	*и методов работы с ними
	*/

	/**********************************************************************
	***********				Математические константы			***********
	**********************************************************************/
	const double SQRT3 = 1.732050807568;		//Корень из 3-х
	const double SQRT2 = 1.414213562373;		//Корень из 2-х
	const double SQRT2_3 = 0.816496580927;		//Корень из 2/3
	const double SQRT3_2 = 1.224744871391;		//Корень из 3/2

	const double PI = 3.141592653589;			//Число Пи
	const double PI_2 = 1.570796326794;			//Пи/2
	const double PIx2 = 6.283185307179;			//2Пи

	const double EPS = 1e-10;					//Малая величина

	const int DIM = 3;	//Размерность пространства


	/***********************************************************************
	****************                  Вектор               *****************
	***********************************************************************/

	class Vector
	{
	public:
		double C[DIM];					//Компоненты вектора

		double getNorm();				//Возвращает норму вектора
		int Normalize();				//Нормализует вектор
		void setZero();					//Обнуляет компоненты вектора
		void set(double,
			double, double);			//Задаёт значения компонент
		double ScalMult(Vector v);		//Скалярное произведение векторов
		Vector VectMult(Vector v);		//Векторное произведение векторов
			
		Vector operator + (Vector);		//Оператор сложения
		Vector operator - (Vector);		//Оператор вычитания
		Vector operator * (double);
		void operator += (Vector);		//Оператор прибавления вектора
		void operator -= (Vector);		//Оператор вычитания вектора
		void operator *= (double);		//Оператор умножения вектора на число
		int operator /= (double);		//Оператор деления вектора на число

		Vector();
		~Vector();

	private:

	};


	/***********************************************************************
	****************                  Тензор               *****************
	***********************************************************************/

	class Tensor
	{
	public:
		double C[DIM][DIM];				//Компоненты тензора

		double getDet();				//Возвращает определитель матрицы компонент
		void set(double, double, double,
			double, double, double,
			double, double, double);	//Задание всех компонент тензора
		void setZero();					//Зануляет компоненты тензора
		void setUnit();					//Делает матрицу компонент тензора единичной
		void Transp();					//Транспонирует матрицу компонент тензора
		double doubleScalMult(Tensor);	//Свёртка (двойное скалярное произведение тензоров)
		Tensor getSymmetryPart();		//Возвращает симметричную часть тензора
		Tensor getAntiSymmetryPart();	//Возвращает антисимметричную часть тензора
		Vector getRow(int);				//Вектор из компонент заданной строки
		Vector getCol(int);				//Вектор из компонент заданного стобца
		
		Tensor operator + (Tensor);		//Оператор сложения тензоров
		Tensor operator - (Tensor);		//Оператор вычитания тензоров
		Tensor operator * (Tensor);		//Оператор умножения тензоров
		void operator += (Tensor);		//Оператор прибавления тензора
		void operator -= (Tensor);		//Оператор убавления тензора
		void operator *= (Tensor);		//Оператор домножения на тензор
		void operator *= (double);		//Оператор умножения тензора на число
		int operator /= (double);		//Оператор деления тензора на число
		
		Tensor();
		~Tensor();

	private:

	};

	/***********************************************************************
	****************         Система скольжения            *****************
	***********************************************************************/

	class SlipSystem
	{
	public:
		Vector n;						//Вектор нормали
		Vector b;						//Вектор Бюргерса
		double t;						//Действующее касательное напряжение
		double tc;						//Критическое касательное напряжение
		double dgm;						//Скорость сдвига
		double gmm;						//Накопленный сдвиг

		void Initialize(float, float, float,
			float, float, float);	//Инициализация значений
		void Initialize(float, float, float, float,
			float, float, float, float, float);
		
		SlipSystem();
		~SlipSystem();
	};

	/***********************************************************************
	**************         Тензор четвёртого ранга            **************
	***********************************************************************/

	class Tensor4
	{
	public:
		double C[DIM][DIM][DIM][DIM];	//Компоненты тензора

		void setZero();					//Обнуление компонент тензора
		void Symmetrize();				//Симметризация компонент тензора

		Tensor4 ToLSK(Tensor O);		//Перевод компонент тензора в ЛСК

		void operator += (Tensor4);		//Оператор прибавления тензора
		void operator -= (Tensor4);		//Оператор отнимания тензора
		void operator *= (double);		//Оператор умножения тензора на число
		int operator /= (double);		//Оператор деления тензора на число

		Tensor4();
		~Tensor4();
	private:
	};


	int LeviCivit(int i, int j, int k);	//Псевдо-тензор Леви-Чивита


	Tensor VectMult(Vector, Tensor);	//Векторное произведение вектора на тензор
	Tensor VectMult(Tensor, Vector);	//Векторное произведение тензора на вектор
	Vector ScalMult(Vector, Tensor);	//Скалярное произведение вектора на тензор
	Vector ScalMult(Tensor, Vector);	//Скалярное произведение тензора на вектор
}

#endif __MATHCORE_H