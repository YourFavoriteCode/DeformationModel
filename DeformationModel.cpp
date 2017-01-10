﻿#include "stdafx.h"
#define _USE_MATH_DEFINES

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

	Polycrystall PC;
	PC.Init(fragm_count);

	unsigned long t1, t2;			//Отсечки времени

	std::cout << " Initializing all fragments... ";
	t1 = clock();

	PC.D = gradV.getSymmetryPart();
	PC.W = gradV.getAntiSymmetryPart();

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

	PC.setParams();
	PC.MakeStruct();



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
							StreamO >> PC.C[q1][q2][q3].o.C[i][j];
						}
					}
					for (int h = 0; h < surround_count; h++)//Считываем значения нормалей
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
					WriteDebugInfo(StreamO, PC.C[q1][q2][q3].o.C);//Записываем значения тензоров ориентации
					for (int h = 0; h < surround_count; h++)//Записываем значения нормалей
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

	PC.OpenFiles();
	//Сохранение начальных полюсных фигур и ССТ
	if (polus_period > 0)
	{
		std::cout << " Saving pole figures... ";
		t1 = clock();
		PC.SavePoleFig();
		t2 = clock();
		std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;
	}

	{
		std::ofstream dbg;
		dbg.open("scal.txt", std::ios_base::out | std::ios_base::trunc);
		for (int i = 0; i < PC.C[0][0][0].SS_count; i++)
		{
			double sc = PC.C[0][0][0].SS[i].n.ScalMult(PC.C[0][0][0].SS[i].b);
			dbg << sc << " ";
		}
		dbg.close();
	}

	t1 = clock();
	PC.Deformate();
	t2 = clock();//Финальная отсечка времени

	if (read_init_stress)	//Сохранение остаточных напряжений
	{
		PC.dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
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

	if (isnan(PC.Strain)) std::cout << std::endl << " Calculation ERROR!" << std::endl;
	else std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << " Done    " << std::endl;
	std::cout << " __________________________________________________" << std::endl;
	std::cout << " Processing time: " << (t2 - t1) / 1000.0 << " sec" << std::endl;
	std::cout << " Number of steps: " << PC.CURR_STEP << std::endl;
	if (isnan(PC.Strain))//Если не зафиксированы ошибки - закрытие
	{
		std::cout << " __________________________________________________" << std::endl;
		std::cout << " Press any key or STOP button to exit...";
		std::system("title Done");
		std::cin.get();
	}

	PC.CloseFiles();

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

