#include "stdafx.h"

#include <fstream>
#include <iostream>
#include <ctime>
#include <Windows.h>

#include "params.h"
#include "functions.h"
#include "Polycrystall.h"

using namespace model;

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc == 1) return 1;	//Программа закроется, если вызвана без аргументов
	
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
	
	/*****************************************************************
	********     Интерфейс ввода/вывода параметров модели     ********
	*****************************************************************/
	std::cout << " Build on " << __DATE__ << " " << __TIME__ << std::endl;
	char* param_file = new char[256];
	wcstombs(param_file, argv[1], 256);//Получили имя файла с параметрами
	std::cout << " Parameters file: " << param_file << std::endl;
	if (prms::ReadParams(param_file) == 1) std::cout << " Error in file!" << std::endl;		//Считали параметры из файла
	delete param_file;//Больше не нужен
	std::cout << " ==========================================" << std::endl;
	const int total_fragm_count = (int)pow(prms::fragm_count, 3);	//Общее кол-во фрагментов
	std::cout << " Fragments count: " << total_fragm_count << std::endl;
	std::cout << " Max. strain: " << prms::strain_max << std::endl;
	std::cout << " Integration step: " << prms::dt << std::endl;
	if (prms::ROTATIONS_HARDENING)
	{
		std::cout << " Using rotation hardening model" << std::endl;
	}
	if (prms::ROTATIONS_TAYLOR)
	{
		std::cout << " Using Taylor's rotation model" << std::endl;
	}
	if (prms::ROTATIONS_TRUSOV)
	{
		std::cout << " Using Trusov's rotation model" << std::endl;
		std::cout << "       A: " << prms::ROT_A << std::endl;
		std::cout << "       H: " << prms::ROT_H << std::endl;
		std::cout << "       L: " << prms::ROT_L << std::endl;
		std::cout << "       MC: " << prms::ROT_MC << std::endl;
	}
	if (prms::HARDENING_BASE)
	{
		std::cout << " Using basic hardening" << std::endl;
		std::cout << "       A: " << prms::HARD_BASE_A << std::endl;
		std::cout << "       Delta: " << prms::HARD_BASE_DELTA << std::endl;
		std::cout << "       Psi: " << prms::HARD_BASE_PSI << std::endl;
	}
	if (prms::HARDENING_BOUND)
	{
		std::cout << " Using boundary hardening" << std::endl;
		std::cout << "       K: " << prms::HARD_BOUND_K << std::endl;
	}
	if (prms::debug_period > 0)
	{
		std::cout << " ============= DEBUG MODE =============" << std::endl;
		std::cout << "       Period: " << prms::debug_period << std::endl;
		std::cout << "       Start: " << prms::DEBUG_START << std::endl;
		std::cout << "       Stop: " << prms::DEBUG_STOP << std::endl;
	}
	if (prms::fix_orient == 1)
	{
		std::cout << " Saving current orientations and normals" << std::endl;
	}
	if (prms::fix_orient == 2)
	{
		std::cout << " Reading saved orientations and normals" << std::endl;
	}
	
	Polycrystall PC;				//Создание и инициализация поликристалла
	PC.Init(prms::fragm_count);

	unsigned long t1, t2;			//Отсечки времени

	std::cout << " Initializing all fragments... ";
	t1 = clock();

	PC.D = prms::gradV.getSymmetryPart();
	PC.W = prms::gradV.getAntiSymmetryPart();

	std::srand(time(NULL));

	switch (prms::SurroundsGrade)			//Степень учёта соседних элементов
	{
	case 0:
	{
		prms::surround_count = 6;			//Обычный уровень
		break;
	}
	case 1:
	{
		prms::surround_count = 18;		//Повышенный уровень
		break;
	}
	case 2:
	{
		prms::surround_count = 26;		//Самый высокий уровень
		break;
	}
	}

	PC.setParams();					//Заполнение всех параметров поликристалла
	PC.MakeStruct();				//Формирование зёренной структуры
	PC.MakeGrains();


	if (prms::fix_orient == 2)	//Считывание записанных ориентаций
	{
		std::ifstream StreamO("DBG\\o.txt", std::ios_base::in);
		std::ifstream StreamNorm("DBG\\Norm.txt", std::ios_base::in);
		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{
					for (int i = 0; i < DIM; i++)	//Считываем значения ориентационных тензоров
					{
						for (int j = 0; j < DIM; j++)
						{
							StreamO >> PC.C[q1][q2][q3].o.C[i][j];
						}
					}
					for (int h = 0; h < prms::surround_count; h++)//Считываем значения нормалей
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm >> PC.C[q1][q2][q3].normals[h].C[i];
						}
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (prms::fix_orient == 1)//Запоминание начальных ориентаций
	{
		std::ofstream StreamO("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
		std::ofstream StreamNorm("DBG\\Norm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{
					WriteDebugInfo(StreamO, PC.C[q1][q2][q3].o.C);//Записываем значения тензоров ориентации
					for (int h = 0; h < prms::surround_count; h++)//Записываем значения нормалей
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm << PC.C[q1][q2][q3].normals[h].C[i] << " ";
						}
						StreamNorm << std::endl;
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (prms::read_init_stress)	//Считывание остаточных напряжений
	{
		std::ifstream StreamSgm("DBG\\sgm.txt", std::ios_base::in);

		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{
					for (int i = 0; i < DIM; i++)	//Считываем значения тензоров
					{
						for (int j = 0; j < DIM; j++)
						{
							StreamSgm >> PC.C[q1][q2][q3].sgm.C[i][j];
						}
					}

				}
			}
		}
		StreamSgm.close();

	}
	t2 = clock();
	std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;

	PC.OpenFiles();				//Открытие и очистка файлов для вывода
	//Сохранение начальных полюсных фигур и ССТ
	if (prms::polus_period > 0)
	{
		std::cout << " Saving pole figures... ";
		t1 = clock();
		PC.SavePoleFig();
		t2 = clock();
		std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;
	}

	/*{
		std::ofstream dbg;
		dbg.open("scal.txt", std::ios_base::out | std::ios_base::trunc);
		for (int i = 0; i < PC.C[0][0][0].SS_count; i++)
		{
			double sc = PC.C[0][0][0].SS[i].n.ScalMult(PC.C[0][0][0].SS[i].b);
			dbg << sc << " ";
		}
		dbg.close();
	}*/

	t1 = clock();		//Начальная отсечка времени
	PC.Deformate();		//Деформирование
	t2 = clock();		//Финальная отсечка времени

	if (prms::read_init_stress)	//Сохранение остаточных напряжений
	{
		PC.dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{

					WriteDebugInfo(PC.dbgstream[3], PC.C[q1][q2][q3].sgm.C);

				}
			}
		}
		PC.dbgstream[3].close();
	}

	/***********************************************************
	*********       Информация о времени и шагах		********
	*********	(автоматически отправляется в панель	********
	*********	управления и программа закрывается,		********
	*********		если она была запущена из неё		********
	***********************************************************/

	if (!isnormal(PC.Strain)) std::cout << std::endl << " Calculation ERROR!" << std::endl;
	else std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << " Done    " << std::endl;
	std::cout << " ==================================================" << std::endl;
	std::cout << " Processing time: " << (t2 - t1) / 1000.0 << " sec" << std::endl;
	std::cout << " Number of steps: " << PC.CURR_STEP << std::endl;
	if (!isnormal(PC.Strain))//Если не зафиксированы ошибки - закрытие
	{
		std::cout << " ==================================================" << std::endl;
		std::cout << " Press any key or STOP button to exit...";
		std::system("title Done");
		std::cin.get();
	}
	
	PC.CloseFiles();				//Сохранение всех полученных данных

	/************************************************************
	*******      Передача данных в панель управления      *******
	************************************************************/
	HWND hwnd;
	hwnd = ::FindWindow(NULL, L"Панель управления моделью");
	if (hwnd != NULL)
	{
		::SendMessage(hwnd, WM_USER + 1, PC.CURR_STEP, (t2 - t1));
		//Аргумент 1 - количество шагов
		//Аргумент 2 - затраченное время (мс)
	}

	return 0;
}

