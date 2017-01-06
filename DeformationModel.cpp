#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Windows.h>
#include <omp.h>

#include "MathCore.h"
#include "params.h"
#include "fragment.h"
#include "rotations.h"
#include "hardening.h"
#include "functions.h"
#include "tension.h"
#include "distributions.h"
#include "Fragmentation.h"

using namespace model;

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc == 1) return 1;	//Программа закроется, если вызвана без аргументов
	
	char* param_file = new char[256];
	wcstombs(param_file, argv[1], 256);//Получили имя файла с параметрами
	ReadParams(param_file);		//Считали параметры из файла
	
	/****************************************************************
	*********	  Создание несуществующих директорий		*********
	****************************************************************/
	if (!isDirectoryExists(L"Plot"))
	{
		CreateDirectory(L"Plot", NULL);//Для графиков
	}
	if (!isDirectoryExists(L"Polus"))
	{
		CreateDirectory(L"Polus", NULL);//Для ПФ
	}
	if (!isDirectoryExists(L"DBG"))
	{
		CreateDirectory(L"DBG", NULL);//Для отладочных данных
	}
	const int total_fragm_count = (int)pow(fragm_count, 3);	//Общее кол-во фрагментов

	/*****************************************************************
	********     Интерфейс ввода/вывода параметров модели     ********
	*****************************************************************/
	std::cout << " Build on " << __DATE__ << " " << __TIME__ << std::endl;
	std::cout << " Parameters file: " << param_file << std::endl;
	delete param_file;//Больше не нужен
	std::cout << " __________________________________________" << std::endl;
	std::cout << " Fragments count: " << total_fragm_count << std::endl;
	std::cout << " Max. strain: " << strain_max << std::endl;
	std::cout << " Integration step: " << dt << std::endl;
	if (ROTATIONS_HARDENING)
	{
		std::cout << " Using rotation hardening model" << std::endl;
	}
	if (ROTATIONS_TAYLOR)
	{
		std::cout << " Using Taylor's rotation model" << std::endl;
	}
	if (ROTATIONS_TRUSOV)
	{
		std::cout << " Using Trusov's rotation model" << std::endl;
		std::cout << "       A: " << ROT_A << std::endl;
		std::cout << "       H: " << ROT_H << std::endl;
		std::cout << "       L: " << ROT_L << std::endl;
		std::cout << "       MC: " << ROT_MC << std::endl;
	}
	if (HARDENING_BASE)
	{
		std::cout << " Using basic hardening" << std::endl;
		std::cout << "       A: " << HARD_BASE_A << std::endl;
		std::cout << "       Delta: " << HARD_BASE_DELTA << std::endl;
		std::cout << "       Psi: " << HARD_BASE_PSI << std::endl;
	}
	if (HARDENING_BOUND)
	{
		std::cout << " Using boundary hardening" << std::endl;
		std::cout << "       K: " << HARD_BOUND_K << std::endl;
	}
	if (debug_period > 0)
	{
		std::cout << " ++++++++++++++DEBUG MODE++++++++++++++" << std::endl;
		std::cout << "       Period: " << debug_period << std::endl;
		std::cout << "       Start: " << DEBUG_START << std::endl;
		std::cout << "       Stop: " << DEBUG_STOP << std::endl;
	}
	if (fix_orient == 1)
	{
		std::cout << " Saving current orientations and normals" << std::endl;
	}
	if (fix_orient == 2)
	{
		std::cout << " Reading saved orientations and normals" << std::endl;
	}

	Fragment ***PC = new Fragment**[fragm_count];		//Поликристалл	
	for (int i = 0; i < fragm_count; i++)
	{
		PC[i] = new Fragment*[fragm_count];
		for (int j = 0; j < fragm_count; j++)
		{
			PC[i][j] = new Fragment[fragm_count];
		}
	}
	unsigned long t1, t2;			//Отсечки времени

	Tensor macro_D;					//Тензор деформации скорости макроуровня
	Tensor macro_W;					//Тензор вихря макроуровня
	Tensor macro_D_in;				//Тензор неупругой части деформации макроуровня
	Tensor macro_E;					//Тензор макродеформаций
	Tensor macro_dSgm;				//Тензор скоростей макронапряжений
	Tensor macro_Sgm;				//Тензор макронапряжений
	Tensor4 macro_P;				//Усреднённый тензор упругих констант
	
	std::cout << " Initializing all fragments... ";
	t1 = clock();

	macro_D = gradV.getSymmetryPart();
	macro_W = gradV.getAntiSymmetryPart();

	std::ofstream TestStream[6];	//Отладочные файловые потоки
	TestStream[0].open("Test0.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[1].open("Test1.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[2].open("Test2.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[3].open("Test3.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[4].open("Test4.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[5].open("Test5.txt", std::ios_base::out | std::ios_base::trunc);
	
	std::srand(time(NULL));

	switch (SurroundsGrade)			//Степень учёта соседних элементов
	{
	case 0:
	{
		surround_count = 6;			//Обычный уровень
		break;
	}
	case 1:
	{
		surround_count = 18;		//Повышенный уровень
		break;
	}
	case 2:
	{
		surround_count = 26;		//Самый высокий уровень
		break;
	}
	}
	
	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{
				//Задание материала 
				int another_material;//Примесная фаза
				another_material = (material == 1) ? 0 : 1;

				int a = (int)(((double)rand() / RAND_MAX) * 100);//На всё воля божья
				if (a <= material_purity)
				{
					PC[q1][q2][q3].setMaterialParams(material);
				}
				else
				{
					PC[q1][q2][q3].setMaterialParams(another_material);
				}
				
				PC[q1][q2][q3].rot_Mc = ROT_MC;//Раздача начальных критических моментов
				PC[q1][q2][q3].rot_A = ROT_A;
				PC[q1][q2][q3].rot_H = ROT_H;
				PC[q1][q2][q3].rot_L = ROT_L;
				PC[q1][q2][q3].position = get1DPos(q1, q2, q3);//Получение порядкового номера
				
				if (RAND_ORIENT)//Случайный равномерный закон
				{

					double a = ((double)rand() / RAND_MAX) * (PI);
					double g = ((double)rand() / RAND_MAX) * (PI);
					double y1 = ((double)rand() / RAND_MAX);
					double y2 = ((double)rand() / RAND_MAX);
					PC[q1][q2][q3].Orientate(a, g, y1, y2);
				}
				else//КСК=ЛСК
				{
					PC[q1][q2][q3].o.setUnit();
				}
				//Задание размеров фрагментов
				switch (fragm_size_law)
				{
				case 0://Равномерное
				{
					PC[q1][q2][q3].size = UniformDistrib(fragm_size_m,fragm_size_dsp);
					break;
				}
				case 1://Нормальное
				{
					PC[q1][q2][q3].size = NormalDistrib(fragm_size_m, fragm_size_dsp);
					break;
				}
				case 2://Логнормальное
				{
					PC[q1][q2][q3].size = LogNormalDistrib(fragm_size_m, fragm_size_dsp);
					break;
				}
				case 3://Показательное
				{
					PC[q1][q2][q3].size = ExpDistrib(fragm_size_m);//Только один параметр
					break;
				}
				}
				PC[q1][q2][q3].volume = pow(PC[q1][q2][q3].size, 3);
				PC[q1][q2][q3].surrounds = new Fragment[surround_count];
				PC[q1][q2][q3].normals = new Vector[surround_count];
				PC[q1][q2][q3].moments = new Vector[surround_count];
				PC[q1][q2][q3].contact = new int[surround_count];

				for (int h = 0; h < surround_count; h++)
				{
					PC[q1][q2][q3].contact[h] = -1;		//Изначально контакт не задан
				}
			}
		}
	}
	
	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{

				for (int h = 0; h < surround_count; h++)
				{
					//Если контакт уже был задан - пропускаем
					if (PC[q1][q2][q3].contact[h] != -1) continue;
					//Определяем, граничат ли фрагменты
					//Первые 6, т.е. боковые грани, граничат всегда
					double a = h < 6 ? 1 : ((double)rand() / RAND_MAX);//На всё воля божья
					if (a < 0.5)
					{
						//Контакта нет - тоже пропускаем
						PC[q1][q2][q3].contact[h] = 0;
						continue;
					}

					int qq1 = q1, qq2 = q2, qq3 = q3, y;
					//qq1, qq2, qq3 - координаты зерна соседа
					//y - номер нормали в соседнем зерне в направлении данного зерна
					double fi = ((double)rand() / RAND_MAX) * (M_PI / 12);
					double cf = cos(fi);
					double sf = sin(fi);
					//TODO: предвычислить наиболее распространенные слагаемые для удобства чтения
					switch (h)
					{
					case 0://Вверх
					{
						PC[q1][q2][q3].normals[h].set(-sf, sf / cf, 1 / cf);
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 5;
						break;
					}
					case 1://От нас
					{
						PC[q1][q2][q3].normals[h].set(-1 / cf, sf, sf / cf);
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						y = 3;
						break;
					}
					case 2://Вправо
					{
						PC[q1][q2][q3].normals[h].set(sf / cf, 1 / cf, -sf);
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						y = 4;
						break;
					}
					case 3://На нас
					{
						PC[q1][q2][q3].normals[h].set(1 / cf, -sf, sf / cf);
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						y = 1;
						break;
					}
					case 4://Влево
					{
						PC[q1][q2][q3].normals[h].set(sf, -1 / cf, -sf / cf);
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 2;
						break;
					}
					case 5://Вниз
					{
						PC[q1][q2][q3].normals[h].set(sf / cf, sf, -1 / cf);
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 0;
						break;
					}
					//Далее идут уже необязательные соседи
					/**************           Рёбра куба          ***************************/
					case 6://Лево от нас
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cf*cos(PI_2), -cos(fi + PI_2) * cf*cos(PI_2), cf*cf*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 9;
						break;
					}
					case 7://Лево на нас
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cf*cos(PI_2), -cos(fi + PI_2) * cf*cos(PI_2), cf*cf*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 8;
						break;
					}
					case 8://право от нас
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cf*cos(PI_2), cos(fi + PI_2) * cf*cos(PI_2), cf*cf*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						y = 7;
						break;
					}
					case 9://Право на нас
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cf*cos(PI_2), cos(fi + PI_2) * cf*cos(PI_2), cf*cf*cos(PI_2));
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						y = 6;
						break;
					}
					case 10://Верх лево
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sf, -cos(fi + PI_2) * cf, sin(fi + PI_2)*cf);
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 15;
						break;
					}
					case 11://Верх право
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sf, cos(fi + PI_2) * cf, sin(fi + PI_2)*cf);
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						y = 14;
						break;
					}
					case 12://Верх на нас
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cf, cos(fi + PI_2) * sf, sin(fi + PI_2)*cf);
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						y = 17;
						break;
					}
					case 13://Верх от нас
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cf, cos(fi + PI_2) * sf, sin(fi + PI_2)*cf);
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						y = 16;
						break;
					}
					case 14://Низ лево
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sf, -cos(fi + PI_2) * cf, -sin(fi + PI_2)*cf);
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 11;
						break;
					}
					case 15://Низ право
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sf, cos(fi + PI_2) * cf, -sin(fi + PI_2)*cf);
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 10;
						break;
					}
					case 16://Низ на нас
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cf, -cos(fi + PI_2) * sf, -sin(fi + PI_2)*cf);
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 13;
						break;
					}
					case 17://Низ от нас
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cf, -cos(fi + PI_2) * sf, -sin(fi + PI_2)*cf);
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 12;
						break;
					}
					/**************      Вершины     *****************/
					case 18://верх лево от нас
					{
						PC[q1][q2][q3].normals[h].set(-cf * cf*cos(PI_2)*cos(PI_2), -cf*cf * cos(PI_2)*cos(PI_2), cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 25;
						break;
					}
					case 19://верх лево на нас
					{
						PC[q1][q2][q3].normals[h].set(cf * cf*cos(PI_2)*cos(PI_2), -cf*cf * cos(PI_2)*cos(PI_2), cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 24;
						break;
					}
					case 20://верх право от нас
					{
						PC[q1][q2][q3].normals[h].set(-cf * cf*cos(PI_2)*cos(PI_2), cf*cf * cos(PI_2)*cos(PI_2), cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 23;
						break;
					}
					case 21://верх право на нас
					{
						PC[q1][q2][q3].normals[h].set(cf * cf*cos(PI_2)*cos(PI_2), cf*cf * cos(PI_2)*cos(PI_2), cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 22;
						break;
					}
					case 22://низ лево от нас
					{
						PC[q1][q2][q3].normals[h].set(-cf * cf*cos(PI_2)*cos(PI_2), -cf*cf * cos(PI_2)*cos(PI_2), -cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 21;
						break;
					}
					case 23://низ лево на нас
					{
						PC[q1][q2][q3].normals[h].set(cf * cf*cos(PI_2)*cos(PI_2), -cf*cf * cos(PI_2)*cos(PI_2), -cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 20;
						break;
					}
					case 24://низ право от нас
					{
						PC[q1][q2][q3].normals[h].set(-cf * cf*cos(PI_2)*cos(PI_2), cf*cf * cos(PI_2)*cos(PI_2), -cf*cf*cos(PI_2)*cos(PI_2));

						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 19;
						break;
					}
					case 25://низ право на нас
					{
						PC[q1][q2][q3].normals[h].set(cf * cf*cos(PI_2)*cos(PI_2), cf*cf * cos(PI_2)*cos(PI_2), -cf*cf*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 18;
						break;
					}
					}

					PC[q1][q2][q3].surrounds[h] = PC[qq1][qq2][qq3];//Здравствуй, сосед!
					PC[qq1][qq2][qq3].surrounds[y] = PC[q1][q2][q3];//Приятно познакомиться!
					PC[q1][q2][q3].normals[h].Normalize();

					for (int i = 0; i < DIM; i++)
					{
						PC[qq1][qq2][qq3].normals[y].C[i] = -PC[q1][q2][q3].normals[h].C[i];//Поделись нормалью
					}

					if (h < 6) PC[q1][q2][q3].contact[h] = 1;//100 % соприкосновения
					else if (h < 14) PC[q1][q2][q3].contact[h] = 3;//5 % соприкосновения
					else PC[q1][q2][q3].contact[h] = 2;//10 % соприкосновения
				}
				if (surround_count > 6)	//Уменьшение объёма
				{
					double a = PC[q1][q2][q3].size * 0.1;				//Длина срезанной части вдоль ребра
					double vol_edge = a*a*PC[q1][q2][q3].size / 2.0;	//Объём, срезанный рёбрами
					double vol_vertex = a*a*a / SQRT3;					//Объём, срезанный вершинами
					int cut_edge = 0;		//Кол-во срезанных рёбер
					int cut_vertex = 0;		//Кол-во срезанных вершин
					for (int h = 6; h < surround_count; h++)
					{
						if (PC[q1][q2][q3].contact[h] != 0)
						{
							if (h < 14) cut_vertex++;
							else cut_edge++;
						}
					}
					PC[q1][q2][q3].volume -= (cut_edge*vol_edge + cut_vertex*vol_vertex);//Вычитание
				}

			}
		}
	}
	
	if (fix_orient == 2)	//Считывание записанных ориентаций
	{
		std::ifstream StreamO("DBG\\o.txt", std::ios_base::in);
		std::ifstream StreamNorm("DBG\\Norm.txt", std::ios_base::in);
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					for (int i = 0; i < DIM; i++)	//Считываем значения ориентационных тензоров
					{
						for (int j = 0; j < DIM; j++)
						{
							StreamO >> PC[q1][q2][q3].o.C[i][j];
						}
					}
					for (int h = 0; h < surround_count; h++)//Считываем значения нормалей
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm >> PC[q1][q2][q3].normals[h].C[i];
						}
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (fix_orient == 1)//Запоминание начальных ориентаций
	{
		std::ofstream StreamO("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
		std::ofstream StreamNorm("DBG\\Norm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					WriteDebugInfo(StreamO, PC[q1][q2][q3].o.C);//Записываем значения тензоров ориентации
					for (int h = 0; h < surround_count; h++)//Записываем значения нормалей
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm << PC[q1][q2][q3].normals[h].C[i] << " ";
						}
						StreamNorm << std::endl;
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (read_init_stress)	//Считывание остаточных напряжений
	{
			std::ifstream StreamSgm("DBG\\sgm.txt", std::ios_base::in);
			
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						for (int i = 0; i < DIM; i++)	//Считываем значения тензоров
						{
							for (int j = 0; j < DIM; j++)
							{
								StreamSgm >> PC[q1][q2][q3].sgm.C[i][j];
							}
						}
						
					}
				}
			}
			StreamSgm.close();
			
	}
	t2 = clock();
	std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;

	int CURR_STEP = 0;				//Текущий шаг интегрирования
	int PLOT_STEP;					//Шаг сохранения графиков
	int POLUS_STEP;					//Шаг сохранения ПФ
	int PROC_STEP = 0;				//Шаг отображения прогресса
	int DEBUG_STEP = 0;				//Шаг записи отладочных данных
	int proc_period = 400;			//Период обновления процента выполнения
	double macro_stress = 0;		//Интенсивность макронапряжений
	double macro_strain = 0;		//Интенсивность макродеформаций

	/*Работа с упругой разгрузкой*/
	double tension_component = 0;	//Вытягивающая компонента в одноосье
	double final_stress = 1e3;		//Значение, до которого разгружать
	double lam = 2.5;				//Коэффициент в разгрузке
	double addition_strain = 1e-4;	//Добавочный множитель для продолжения циклики
	if (REAL_UNIAX)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				if (macro_D.C[i][j] != 0)	//Определение растягивающей компоненты
				{
					tension_component = macro_D.C[i][j];
					break;
				}
			}
		}
	}
	

	/*******************************************************
	**********       Работа с файлами данных       *********
	*******************************************************/

	TruncPoleFiles();				//Очистка всех файлов полюсных фигур
	TruncSSTFiles();
	std::ofstream StrainStream("Plot\\X.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream StressStream("Plot\\Y.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream StrainStreamAll("Plot\\Xall.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream StressStreamAll("Plot\\Yall.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream ActiveSysStream("Plot\\ActiveSS.dat", std::ios_base::out | std::ios_base::trunc | std::ios::binary);
	const int file_count = 16;
	std::ofstream dbgstream[file_count];
	if (debug_period > 0)				//Открытие файлов для отладочных данных
	{
		dbgstream[0].open("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[1].open("DBG\\e.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[2].open("DBG\\d.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[4].open("DBG\\om.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[5].open("DBG\\dsgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[6].open("DBG\\din.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[7].open("DBG\\w.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[8].open("DBG\\dgamma.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[9].open("DBG\\t.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[10].open("DBG\\Macro_D.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[11].open("DBG\\Macro_Din.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[12].open("DBG\\Macro_Sgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[13].open("DBG\\Macro_dSgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[14].open("DBG\\Macro_E.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[15].open("DBG\\VOL_M.txt", std::ios_base::out | std::ios_base::trunc);
	}
	
	//Сохранение начальных полюсных фигур и ССТ
	if (polus_period > 0)
	{
		std::cout << " Saving pole figures... ";
		t1 = clock();
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					GetPoleFig(&PC[q1][q2][q3]);
					if (SST_SAVING) GetSST(&PC[q1][q2][q3]);
				}
			}
		}
		t2 = clock();
		std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;
	}
	
	{
		std::ofstream dbg;
		dbg.open("scal.txt", std::ios_base::out | std::ios_base::trunc);
		for (int i = 0; i < PC[0][0][0].SS_count; i++)
		{
			double sc = PC[0][0][0].SS[i].n.ScalMult(PC[0][0][0].SS[i].b);
			dbg << sc << " ";
		}
		dbg.close();
	}

	/*******************************************************
	**********      Цикл по этапам нагружения      *********
	*******************************************************/
	t1 = clock();
	
	for (int cycle = 0; cycle < cycle_count; cycle++)
	{
		PLOT_STEP = 0;
		POLUS_STEP = 0;
		if (cycle_count > 1) std::cout << std::endl << " Cycle # " << cycle+1;
		std::cout << std::endl << "        0.00%";

		/***************************************
		********    Цикл по времени    *********
		***************************************/
		while (macro_strain < strain_max)
		{
			//Раскомментировать строчки внутри куска
			if (REAL_UNIAX)	//Одноосное растяжение
			{
				//Осреднение
				macro_P.setZero();
				macro_D_in.setZero();
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_D_in += PC[q1][q2][q3].d_in;
							macro_P += PC[q1][q2][q3].p.ToLSK(PC[q1][q2][q3].o);
						}
					}
				}

				macro_D_in /= (total_fragm_count);
				macro_P /= (total_fragm_count);

				//Симметризация
				macro_P.Symmetrize();

				macro_D = TensionStrainCalc(macro_P, macro_D_in, macro_D.C[0][0]);

				Tensor b = macro_D;
				b *= dt;
				macro_E += b;
				Tensor buf = macro_E;
				macro_strain = macro_E.doubleScalMult(buf);
				macro_strain = SQRT2_3*sqrt(macro_strain);

				macro_dSgm = TensionStressCalc(macro_P, macro_D_in, macro_D);
				macro_dSgm *= dt;				//Приращение напряжений на шаге
				macro_Sgm += macro_dSgm;
				buf = macro_Sgm;
				macro_stress = macro_Sgm.doubleScalMult(buf);
				macro_stress = SQRT3_2*sqrt(macro_stress);
			}
			else
			{
				macro_stress = 0;		//Вычисление интенсивностей осреднением
				macro_strain = 0;
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_strain += PC[q1][q2][q3].strain;
							macro_stress += PC[q1][q2][q3].stress;
						}
					}
				}
				macro_strain /= total_fragm_count;
				macro_stress /= total_fragm_count;
			}
			omp_set_num_threads(thread_count);
			#pragma omp parallel for
			//Часть, которую можно паралелить
			//Здесь необходимо гарантировать защиту данных каждого фрагмента
			//от перезаписи другими фрагментами
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{

						/**************************************************
						************       Переходим в КСК       **********
						**************************************************/

						Tensor O = PC[q1][q2][q3].o;
						Tensor OT = O;
						OT.Transp();
						PC[q1][q2][q3].d = O*macro_D*OT;//Гипотеза Фойгта
						PC[q1][q2][q3].w = O*macro_W*OT;//Расширенная

						PC[q1][q2][q3].sgm = O*PC[q1][q2][q3].sgm*OT;
						PC[q1][q2][q3].d_in = O*PC[q1][q2][q3].d_in*OT;


						/***************************************************
						***********       Пересчитываем НДС      ***********
						***************************************************/

						PC[q1][q2][q3].NDScalc();

						if (HARDENING_BASE)			//Базовое упрочнение
						{
							Base_hardening(&PC[q1][q2][q3]);
						}
						if (ROTATIONS_TAYLOR)		//Ротации по Тейлору
						{
							Taylor_rotations(&PC[q1][q2][q3]);
						}
						if (ROTATIONS_TRUSOV && ROTATIONS_HARDENING)	//Ротационное упрочнение
						{
							Rotation_hardening(&PC[q1][q2][q3]);
						}

						/**************************************************
						************       Переходим в ЛСК       **********
						**************************************************/

						PC[q1][q2][q3].sgm = OT*PC[q1][q2][q3].sgm*O;
						PC[q1][q2][q3].d_in = OT*PC[q1][q2][q3].d_in*O;

					}
				}
			}
			/*Механизмы, учитывающие соседние элементы*/
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						if (HARDENING_BOUND)	//Зернограничное упрочнение
						{
							Boundary_hardening(&PC[q1][q2][q3]);
						}

						if (ROTATIONS_TRUSOV)		//Ротации по Трусову
						{
							Trusov_rotations(&PC[q1][q2][q3]);
						}
						
					}
				}
			}


			if (!REAL_UNIAX)		//Этот блок нужен исключительно для работы с энергией!
			{
				macro_Sgm.setZero();
				macro_D_in.setZero();
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_Sgm += PC[q1][q2][q3].sgm;
							macro_D_in += PC[q1][q2][q3].d_in;
						}
					}
				}
				macro_Sgm /= total_fragm_count;
				macro_D_in /= total_fragm_count;
			}


			/************************************************************
			***********	        Прогресс выполнения 	      ***********
			************************************************************/
						
			double progress = (macro_strain / strain_max) * 100;
			
			if (!(cycle_count == 1 || cycle == 0))	//Многоцикловые нагружения
			{
				progress /= 2;
				if (macro_E.C[0][0] > 0)			//Этот код позволяет корректно
				{									//отображать прогресс выполнения,
					if (macro_Sgm.C[0][0] > 0)		//когда на графике петли циклические
					{
						progress += 50.0;
					}
					else
					{
						progress = 50.0 - progress;
					}
				}
				else
				{
					if (macro_Sgm.C[0][0] > 0)
					{
						progress = 50.0 - progress;
					}
					else
					{
						progress += 50.0;
					}
				}
			}

			if (PROC_STEP == proc_period)
			{
				PROC_STEP = 0;
				std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
				printf("        %2.2f", progress);
				std::cout << "%";
			}

			/************************************************************
			***********	    Запись данных для графиков НДС    ***********
			************************************************************/
			
			if (progress - PLOT_STEP > plot_period && plot_period > 0)
			{
				if (!REAL_UNIAX)
				{
					StrainStream.write((char *)&macro_strain, sizeof macro_strain);
					StressStream.write((char *)&macro_stress, sizeof macro_stress);
				}
				else
				{
					StrainStream.write((char *)&macro_E.C[0][0], sizeof macro_E.C[0][0]);
					StressStream.write((char *)&macro_Sgm.C[0][0], sizeof macro_Sgm.C[0][0]);
				}
				double ActiveSysCount = 0;			//Среднее кол-во активных систем скольжения на шаге
				double RotEnergy = 0;				//Энергия ротаций на шаге
				double RotSpeed = 0;				//Средняя скорость вращения на шаге
				int RotCount = 0;					//Кол-во вращающихся фрагментов
				double norma = 0;
				double Mc = 0;
				double dmc = 0;
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							if (!REAL_UNIAX)
							{
								StrainStreamAll.write((char *)&PC[q1][q2][q3].strain, sizeof PC[q1][q2][q3].strain);
								StressStreamAll.write((char *)&PC[q1][q2][q3].stress, sizeof PC[q1][q2][q3].stress);
							}
							else
							{
								StrainStreamAll.write((char *)&PC[q1][q2][q3].e.C[0][0], sizeof PC[q1][q2][q3].e.C[0][0]);
								StressStreamAll.write((char *)&PC[q1][q2][q3].sgm.C[0][0], sizeof PC[q1][q2][q3].sgm.C[0][0]);
							}
							
							for (int i = 0; i < PC[q1][q2][q3].SS_count; i++)
							{
								if (PC[q1][q2][q3].SS[i].dgm > EPS) ActiveSysCount++;//Подсчёт активных СС
							}
							
							if (PC[q1][q2][q3].isRotate) RotCount++;		//Подсчёт вращающихся решёток
							RotEnergy += PC[q1][q2][q3].rot_energy;				//Суммирование энергий вращения
							RotSpeed += PC[q1][q2][q3].rot_speed;				//Суммирование скоростей вращения
							norma += PC[q1][q2][q3].norm;
							Mc += PC[q1][q2][q3].rot_Mc;
							dmc += PC[q1][q2][q3].dmc;
							
						}
					}
				}
				//dbgstream[15] << std::endl;
				norma /= total_fragm_count;
				Mc /= total_fragm_count;
				dmc /= total_fragm_count;
				ActiveSysCount /= total_fragm_count;
				ActiveSysStream.write((char *)&ActiveSysCount, sizeof ActiveSysCount);//Запись кол-ва активных СС
				if (RotCount != 0)
				{
					RotSpeed /= RotCount;
				}
				else RotSpeed = 0;
				
				/*******************************************************
				********* 		   Работа с энергией          **********
				*******************************************************/

				//Полная энергия деформирования - сумма элементарных энергий на каждом шаге
				//Элементарная энергия - свёртка напряжений с приращением деформации
				//Энергия ротаций - момент*приращение угла

				Tensor dE = macro_D;
				dE *= dt;			//Приращение деформации на шаге

				double StepEnergy = macro_Sgm.doubleScalMult(dE);	//Полная энергия на шаге

				TestStream[0] << RotCount << std::endl;	
				TestStream[1] << RotSpeed << std::endl;
				TestStream[2] << RotEnergy << std::endl;
				TestStream[3] << StepEnergy << std::endl;
				TestStream[4] << Mc << std::endl;
				TestStream[5] << norma << std::endl;

				PLOT_STEP =  progress;
			}
			
			/************************************************************
			***********	      Сохранение полюсных фигур	      ***********
			************************************************************/
			if (progress - POLUS_STEP > polus_period && polus_period > 0)
			{
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							GetPoleFig(&PC[q1][q2][q3]);
							if (SST_SAVING) GetSST(&PC[q1][q2][q3]);
						}
					}
				}
				POLUS_STEP = progress;
			}

			/************************************************************
			***********	       Запись пошаговых данных	      ***********
			************************************************************/
			if (CURR_STEP >= DEBUG_START && CURR_STEP <= DEBUG_STOP && DEBUG_STEP == debug_period)
			{
				DEBUG_STEP = 0;

				for (int i = 0; i < file_count; i++)
				{
					dbgstream[i] << "********      STEP " << CURR_STEP << "      ********" << std::endl << std::endl;
				}

				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							WriteDebugInfo(dbgstream[0], PC[q1][q2][q3].o.C);
							WriteDebugInfo(dbgstream[1], PC[q1][q2][q3].e.C);
							WriteDebugInfo(dbgstream[2], PC[q1][q2][q3].d.C);
							WriteDebugInfo(dbgstream[3], PC[q1][q2][q3].sgm.C);
							WriteDebugInfo(dbgstream[4], PC[q1][q2][q3].om.C);
							WriteDebugInfo(dbgstream[5], PC[q1][q2][q3].dsgm.C);
							WriteDebugInfo(dbgstream[6], PC[q1][q2][q3].d_in.C);
							WriteDebugInfo(dbgstream[7], PC[q1][q2][q3].w.C);
							for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
							{
								dbgstream[8] << PC[q1][q2][q3].SS[f].dgm << " ";
							}
							dbgstream[8] << std::endl << std::endl;
							for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
							{
								dbgstream[9] << PC[q1][q2][q3].SS[f].t << " ";
							}
							dbgstream[9] << std::endl << std::endl;
							dbgstream[15] << PC[q1][q2][q3].dmoment.C[0] << " " << PC[q1][q2][q3].dmoment.C[1] << " " << PC[q1][q2][q3].dmoment.C[2]<<std::endl;
						}
					}
				}
				WriteDebugInfo(dbgstream[10], macro_D.C);
				WriteDebugInfo(dbgstream[11], macro_D_in.C);
				WriteDebugInfo(dbgstream[12], macro_Sgm.C);
				WriteDebugInfo(dbgstream[13], macro_dSgm.C);
				WriteDebugInfo(dbgstream[14], macro_E.C);
			}
			CURR_STEP++;
			PROC_STEP++;
			if (CURR_STEP >= DEBUG_START && CURR_STEP <= DEBUG_STOP) DEBUG_STEP++;
		}

		if (UNLOADING)	//Упругая разгрузка представительного объёма
		{
			std::cout << std::endl << " Unloading # " << cycle;
			std::cout << std::endl << "        0.00%";
			while (abs(macro_stress) > final_stress)
			{
				//Осреднение
				macro_P.setZero();
				macro_D_in.setZero();
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_D_in += PC[q1][q2][q3].d_in;
							macro_P += PC[q1][q2][q3].p.ToLSK(PC[q1][q2][q3].o);
						}
					}
				}
				macro_D_in /= total_fragm_count;
				macro_P /= total_fragm_count;

				//Симметризация
				macro_P.Symmetrize();
				macro_D = UnloadingStrainCalc(macro_P, macro_D_in, macro_Sgm, lam);
				Tensor b = macro_D;
				b *= dt;
				macro_E += b;
				Tensor buf = macro_E;
				macro_strain = macro_E.doubleScalMult(buf);
				macro_strain = SQRT2_3*sqrt(macro_strain);
				
				macro_dSgm = TensionStressCalc(macro_P, macro_D_in, macro_D);
				macro_dSgm *= dt;
				macro_Sgm += macro_dSgm;
				buf = macro_Sgm;
				macro_stress = macro_Sgm.doubleScalMult(buf);
				macro_stress = SQRT3_2*sqrt(macro_stress);

				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{

							/**************************************************
							************       Переходим в КСК       **********
							**************************************************/

							Tensor O = PC[q1][q2][q3].o;	//Гипотеза Фойгта
							Tensor OT = O;
							OT.Transp();
							PC[q1][q2][q3].d = O*macro_D*OT;
							PC[q1][q2][q3].w = O*macro_W*OT;
							PC[q1][q2][q3].sgm = O*PC[q1][q2][q3].sgm*OT;
							PC[q1][q2][q3].d_in = O*PC[q1][q2][q3].d_in*OT;


							/***************************************************
							***********       Пересчитываем НДС      ***********
							***************************************************/

							PC[q1][q2][q3].NDScalc();

							if (HARDENING_BASE)			//Базовое упрочнение
							{
								Base_hardening(&PC[q1][q2][q3]);
							}
							if (HARDENING_BOUND)
							{
								Boundary_hardening(&PC[q1][q2][q3]);
							}
							if (ROTATIONS_TAYLOR)		//Ротации по Тейлору
							{
								Taylor_rotations(&PC[q1][q2][q3]);
							}
							if (ROTATIONS_TRUSOV)		//Ротации по Трусову
							{
								Trusov_rotations(&PC[q1][q2][q3]);
							}
							if (ROTATIONS_TRUSOV && ROTATIONS_HARDENING)	//Ротационное упрочнение
							{
								Rotation_hardening(&PC[q1][q2][q3]);
							}

							/**************************************************
							************       Переходим в ЛСК       **********
							**************************************************/
							PC[q1][q2][q3].sgm = OT*PC[q1][q2][q3].sgm*O;
							PC[q1][q2][q3].d_in = OT*PC[q1][q2][q3].d_in*O;

						}
					}
				}

				/************************************************************
				***********	        Прогресс выполнения 	      ***********
				************************************************************/

				double progress = (final_stress / macro_stress) * 100;

				if (PROC_STEP == proc_period/4)
				{
					PROC_STEP = 0;
					std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
					printf("        %2.2f", progress);
					std::cout << "%";
				}

				/************************************************************
				***********	    Запись данных для графиков НДС    ***********
				************************************************************/
				if (plot_period > 0)
				{
					
					StrainStream.write((char *)&macro_E.C[0][0], sizeof macro_E.C[0][0]);
					StressStream.write((char *)&macro_Sgm.C[0][0], sizeof macro_Sgm.C[0][0]);
					double ActiveSysCount = 0;		//Среднее кол-во активных систем скольжения
					for (int q1 = 0; q1 < fragm_count; q1++)
					{
						for (int q2 = 0; q2 < fragm_count; q2++)
						{
							for (int q3 = 0; q3 < fragm_count; q3++)
							{
								StrainStreamAll.write((char *)&PC[q1][q2][q3].e.C[0][0], sizeof PC[q1][q2][q3].e.C[0][0]);
								StressStreamAll.write((char *)&PC[q1][q2][q3].sgm.C[0][0], sizeof PC[q1][q2][q3].sgm.C[0][0]);
								for (int i = 0; i < PC[q1][q2][q3].SS_count; i++)
								{
									if (PC[q1][q2][q3].SS[i].dgm != 0) ActiveSysCount++;//Подсчёт активных СС
								}
							
							}
						}
					}
					ActiveSysCount /= total_fragm_count;
					ActiveSysStream.write((char *)&ActiveSysCount, sizeof ActiveSysCount);
					PLOT_STEP = progress;
				}

				/************************************************************
				***********	      Сохранение полюсных фигур	      ***********
				************************************************************/
				if (progress - POLUS_STEP > polus_period && polus_period > 0)
				{
					for (int q1 = 0; q1 < fragm_count; q1++)
					{
						for (int q2 = 0; q2 < fragm_count; q2++)
						{
							for (int q3 = 0; q3 < fragm_count; q3++)
							{
								GetPoleFig(&PC[q1][q2][q3]);
								if (SST_SAVING) GetSST(&PC[q1][q2][q3]);
							}
						}
					}
					POLUS_STEP = progress;
				}
				/************************************************************
				***********	       Запись пошаговых данных	      ***********
				************************************************************/
				if (CURR_STEP >= DEBUG_START && CURR_STEP <= DEBUG_STOP && DEBUG_STEP == debug_period)
				{
					DEBUG_STEP = 0;

					for (int i = 0; i < file_count; i++)
					{
						dbgstream[i] << "********      STEP " << CURR_STEP << "      ********" << std::endl << std::endl;
					}

					for (int q1 = 0; q1 < fragm_count; q1++)
					{
						for (int q2 = 0; q2 < fragm_count; q2++)
						{
							for (int q3 = 0; q3 < fragm_count; q3++)
							{
								WriteDebugInfo(dbgstream[0], PC[q1][q2][q3].o.C);
								WriteDebugInfo(dbgstream[1], PC[q1][q2][q3].e.C);
								WriteDebugInfo(dbgstream[2], PC[q1][q2][q3].d.C);
								WriteDebugInfo(dbgstream[3], PC[q1][q2][q3].sgm.C);
								WriteDebugInfo(dbgstream[4], PC[q1][q2][q3].om.C);
								WriteDebugInfo(dbgstream[5], PC[q1][q2][q3].dsgm.C);
								WriteDebugInfo(dbgstream[6], PC[q1][q2][q3].d_in.C);
								WriteDebugInfo(dbgstream[7], PC[q1][q2][q3].w.C);
								for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
								{
									dbgstream[8] << PC[q1][q2][q3].SS[f].dgm << " ";
								}
								dbgstream[8] << std::endl << std::endl;
								for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
								{
									dbgstream[9] << PC[q1][q2][q3].SS[f].t << " ";
								}
								dbgstream[9] << std::endl << std::endl;
							}
						}
					}
					WriteDebugInfo(dbgstream[10], macro_D.C);
					WriteDebugInfo(dbgstream[11], macro_D_in.C);
					WriteDebugInfo(dbgstream[12], macro_Sgm.C);
					WriteDebugInfo(dbgstream[13], macro_dSgm.C);
					WriteDebugInfo(dbgstream[14], macro_E.C);
				}
				
				CURR_STEP++;
				PROC_STEP++;
				if (CURR_STEP >= DEBUG_START && CURR_STEP <= DEBUG_STOP) DEBUG_STEP++;
			}
		}

		if (cycle_count > 1 && REAL_UNIAX)
		{
			macro_D.C[0][0] = pow(-1, cycle + 1) * tension_component;	//Меняем знак растягивающей компоненты
			strain_max += strain_max * addition_strain;					//Повышаем предел интенсивности
		}

		/****************************************************************************************************
		*****************					    Фрагментация и дробление					*****************
		*****************************************************************************************************/

		if (FRAGMENTATION)
		{
			float** sm_matrix = new float*[total_fragm_count];			//Матрица смежности
			for (int i = 0; i < total_fragm_count; i++)
			{
				sm_matrix[i] = new float[total_fragm_count];
			}

			for (int i = 0; i < total_fragm_count; i++)
			{
				for (int j = 0; j < total_fragm_count; j++)
				{
					sm_matrix[i][j] = -1;
				}
			}

			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						int pos1 = get1DPos(q1, q2, q3);				//Позиция первого элемента
						
						for (int h = 0; h < surround_count; h++)
						{
							int pos2 = PC[q1][q2][q3].surrounds[h].position;//Позиция второго элемента
							//if (pos2 < pos1) continue;			//Раз матрица диагональная - нижнюю половину не нужно отдельно считать
							/*if (PC[q1][q2][q3].contact[h] == 0)
							{
								sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = -1;//Если фрагменты не контактируют
								continue;
							}
							else if (sm_matrix[pos1][pos2]==0)//Если ещё не прошли
							{
								sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = (float)rand() / RAND_MAX;
							}*/
							if (pos1<pos2) sm_matrix[pos1][pos2] = PC[q1][q2][q3].DisorientMeasure(h);
							//sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = PC[q1][q2][q3].contact[h];
						
							//if (pos1<pos2)sm_matrix[pos1][pos2] = PC[q1][q2][q3].surrounds[h].position;
						}


					}
				}
			}
			std::ofstream MatrStream("DBG\\SmMatrix.txt", std::ios_base::out | std::ios_base::trunc);
			for (int i = 0; i < total_fragm_count; i++)
			{
				for (int j = 0; j < total_fragm_count; j++)
				{
					MatrStream << sm_matrix[i][j] << " ";
				}
				MatrStream << std::endl;
			}
			MatrStream.close();

			for (int i = 0; i < total_fragm_count; i++)
			{
				delete [] sm_matrix[i];
			}
			delete[]sm_matrix;
		}


	}

	if (read_init_stress)	//Сохранение остаточных напряжений
	{ 
		dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
				
					WriteDebugInfo(dbgstream[3], PC[q1][q2][q3].sgm.C);
				
				}
			}
		}
		dbgstream[3].close();
	}
	//Сохранение конечных полюсных фигур
/*	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{
				GetPoleFig(PC[q1][q2][q3]);
				if (SST_SAVING) GetSST(PC[q1][q2][q3]);
			}
		}
	}
	*/
	t2 = clock();//Финальная отсечка времени

	/***********************************************************
	*********      Закрытие всех файловых потоков      *********
	***********************************************************/

	StrainStream.close();
	StressStream.close();
	StrainStreamAll.close();
	StressStreamAll.close();
	ActiveSysStream.close();
	for (int i = 0; i < 6; i++)
	{
		TestStream[i].close();
	}


	if (debug_period > 0)
	{
		for (int i = 0; i < file_count; i++)
		{
			dbgstream[i].close();
		}
	}
	/***********************************************************
	*********       Информация о времени и шагах		********
	*********	(автоматически отправляется в панель	********
	*********	управления и программа закрывается,		********
	*********		если она была запущена из неё		********
	***********************************************************/

	if (isnan(macro_strain)) std::cout << std::endl << " Calculation ERROR!" << std::endl;
	else std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << " Done    " << std::endl;
	std::cout << " __________________________________________________" << std::endl;
	std::cout << " Processing time: " << (t2 - t1) / 1000.0 << " sec" << std::endl;
	std::cout << " Number of steps: " << CURR_STEP << std::endl;
	if (isnan(macro_strain))//Если не зафиксированы ошибки - закрытие
	{
		std::cout << " __________________________________________________" << std::endl;
		std::cout << " Press any key or STOP button to exit...";
		std::system("title Done");
		std::cin.get();
	}

	/************************************************************
	***********      Освобождение занятой памяти      ***********
	************************************************************/

	for (int i = 0; i < fragm_count; i++)
	{
		for (int j = 0; j < fragm_count; j++)
		{
			delete[] PC[i][j];
		}
		delete[] PC[i];
	}
	delete[] PC;


	/************************************************************
	*******      Передача данных в панель управления      *******
	************************************************************/
	HWND hwnd;
	hwnd = ::FindWindow(NULL, L"Панель управления моделью");
	if (hwnd != NULL)
	{
		::SendMessage(hwnd, WM_USER + 1, CURR_STEP, (t2 - t1));
		//Аргумент 1 - количество шагов
		//Аргумент 2 - затраченное время (мс)
	}
	return 0;
}

