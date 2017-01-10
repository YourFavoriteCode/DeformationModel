#include "stdafx.h"
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <fstream>

#include "Polycrystall.h"
#include "params.h"
#include "MathCore.h"
#include "Fragmentation.h"
#include "distributions.h"
#include "tension.h"
#include "rotations.h"
#include "hardening.h"
#include "functions.h"

namespace model
{
	Polycrystall::Polycrystall()
	{
		Strain = 0;
		Stress = 0;

		cycle = 0;
		CURR_STEP = 0;
		PROC_STEP = 0;
		POLUS_STEP = 0;
		PLOT_STEP = 0;
		DEBUG_STEP = 0;
		proc_period = 400;
		file_count = 16;

		tension_component = 0;
		final_stress = 1e3;
		lam = 2.5;
		addition_strain = 1e-5;
	}

	Polycrystall::~Polycrystall()
	{
		for (int i = 0; i < fragm_count; i++)
		{
			for (int j = 0; j < fragm_count; j++)
			{
				delete[] C[i][j];
			}
			delete[] C[i];
		}
		delete[] C;
	}

	void Polycrystall::OpenFiles()
	{
		TruncPoleFiles();				//Очистка всех файлов полюсных фигур
		TruncSSTFiles();

		Datastream = new std::ofstream[5];//Открытие файлов для записи кривых НДС
		Datastream[0].open("Plot\\X.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		Datastream[1].open("Plot\\Y.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		Datastream[2].open("Plot\\Xall.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		Datastream[3].open("Plot\\Yall.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		Datastream[4].open("Plot\\ActiveSS.dat", std::ios_base::out | std::ios_base::trunc | std::ios::binary);

		dbgstream = new std::ofstream[file_count];
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

		TestStream = new std::ofstream[6];
		TestStream[0].open("Test0.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[1].open("Test1.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[2].open("Test2.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[3].open("Test3.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[4].open("Test4.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[5].open("Test5.txt", std::ios_base::out | std::ios_base::trunc);
	}

	void Polycrystall::CloseFiles()
	{
		for (int i = 0; i < 5; i++)
		{
			Datastream[i].close();
		}
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
	}

	void Polycrystall::Init(int count)
	{
		fragm_count = count;
		total_fragm_count = (int)pow(count, 3);

		C = new Fragment**[count];		//Выделение памяти под массив
		for (int i = 0; i < count; i++)
		{
			C[i] = new Fragment*[count];
			for (int j = 0; j < count; j++)
			{
				C[i][j] = new Fragment[count];
			}
		}
	}

	void Polycrystall::setParams()
	{
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
						C[q1][q2][q3].setMaterialParams(material);
					}
					else
					{
						C[q1][q2][q3].setMaterialParams(another_material);
					}

					C[q1][q2][q3].rot_Mc = ROT_MC;	//Раздача начальных критических моментов
					C[q1][q2][q3].rot_A = ROT_A;	//и параметров модели ротаций
					C[q1][q2][q3].rot_H = ROT_H;
					C[q1][q2][q3].rot_L = ROT_L;
					C[q1][q2][q3].position = get1DPos(q1, q2, q3);//Получение порядкового номера фрагмента

					if (RAND_ORIENT)//Получение ориентационного тензора (случайный равномерный закон распределения)
					{
						double a = ((double)rand() / RAND_MAX) * (PI);
						double g = ((double)rand() / RAND_MAX) * (PI);
						double y1 = ((double)rand() / RAND_MAX);
						double y2 = ((double)rand() / RAND_MAX);
						C[q1][q2][q3].Orientate(a, g, y1, y2);
					}
					else//Получение ориентационного тензора (КСК=ЛСК)
					{
						C[q1][q2][q3].o.setUnit();
					}
					//Задание размеров фрагментов
					switch (fragm_size_law)
					{
					case 0://Равномерное
					{
						C[q1][q2][q3].size = UniformDistrib(fragm_size_m, fragm_size_dsp);
						break;
					}
					case 1://Нормальное
					{
						C[q1][q2][q3].size = NormalDistrib(fragm_size_m, fragm_size_dsp);
						break;
					}
					case 2://Логнормальное
					{
						C[q1][q2][q3].size = LogNormalDistrib(fragm_size_m, fragm_size_dsp);
						break;
					}
					case 3://Показательное
					{
						C[q1][q2][q3].size = ExpDistrib(fragm_size_m);//Только один параметр
						break;
					}
					}
					C[q1][q2][q3].volume = pow(C[q1][q2][q3].size, 3);	//Объём фрагмента
					
					//Выделение памяти под массивы, необходимые для работы с окружением
					C[q1][q2][q3].surrounds = new Fragment[surround_count];
					C[q1][q2][q3].normals = new Vector[surround_count];
					C[q1][q2][q3].moments = new Vector[surround_count];
					C[q1][q2][q3].contact = new int[surround_count];

					for (int h = 0; h < surround_count; h++)
					{
						C[q1][q2][q3].contact[h] = -1;		//Изначально контакт не задан
					}
				}
			}
		}
	}

	void Polycrystall::MakeStruct()
	{
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{

					for (int h = 0; h < surround_count; h++)
					{
						//Если контакт уже был задан - пропускаем
						if (C[q1][q2][q3].contact[h] != -1) continue;
						//Определяем, граничат ли фрагменты
						//Первые 6, т.е. боковые грани, граничат всегда
						double a = h < 6 ? 1 : ((double)rand() / RAND_MAX);//На всё воля божья
						if (a < 0.5)
						{
							//Контакта нет - тоже пропускаем
							C[q1][q2][q3].contact[h] = 0;
							continue;
						}

						int qq1 = q1, qq2 = q2, qq3 = q3, y;
						//qq1, qq2, qq3 - координаты зерна соседа
						//y - номер нормали в соседнем зерне в направлении данного зерна
						double fi = ((double)rand() / RAND_MAX) * (PI / 12);
						//TODO: предвычислить наиболее распространенные слагаемые для удобства чтения
						switch (h)
						{
						case 0://Вверх
						{
							C[q1][q2][q3].normals[h].set(-sin(fi), sin(fi) / cos(fi), 1 / cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 5;
							break;
						}
						case 1://От нас
						{
							C[q1][q2][q3].normals[h].set(-1 / cos(fi), sin(fi), sin(fi) / cos(fi));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							y = 3;
							break;
						}
						case 2://Вправо
						{
							C[q1][q2][q3].normals[h].set(sin(fi) / cos(fi), 1 / cos(fi), -sin(fi));
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							y = 4;
							break;
						}
						case 3://На нас
						{
							C[q1][q2][q3].normals[h].set(1 / cos(fi), -sin(fi), sin(fi) / cos(fi));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							y = 1;
							break;
						}
						case 4://Влево
						{
							C[q1][q2][q3].normals[h].set(sin(fi), -1 / cos(fi), -sin(fi) / cos(fi));
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 2;
							break;
						}
						case 5://Вниз
						{
							C[q1][q2][q3].normals[h].set(sin(fi) / cos(fi), sin(fi), -1 / cos(fi));
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 0;
							break;
						}
						//Далее идут уже необязательные соседи
						/**************           Рёбра куба          ***************************/
						case 6://Лево от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 9;
							break;
						}
						case 7://Лево на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 8;
							break;
						}
						case 8://право от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							y = 7;
							break;
						}
						case 9://Право на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							y = 6;
							break;
						}
						case 10://Верх лево
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 15;
							break;
						}
						case 11://Верх право
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							y = 14;
							break;
						}
						case 12://Верх на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							y = 17;
							break;
						}
						case 13://Верх от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							y = 16;
							break;
						}
						case 14://Низ лево
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 11;
							break;
						}
						case 15://Низ право
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 10;
							break;
						}
						case 16://Низ на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 13;
							break;
						}
						case 17://Низ от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 12;
							break;
						}
						/**************      Вершины     *****************/
						case 18://верх лево от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 25;
							break;
						}
						case 19://верх лево на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 24;
							break;
						}
						case 20://верх право от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 23;
							break;
						}
						case 21://верх право на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 22;
							break;
						}
						case 22://низ лево от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 21;
							break;
						}
						case 23://низ лево на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 20;
							break;
						}
						case 24://низ право от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));

							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 19;
							break;
						}
						case 25://низ право на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 18;
							break;
						}
						}

						C[q1][q2][q3].surrounds[h] = C[qq1][qq2][qq3];//Здравствуй, сосед!
						C[qq1][qq2][qq3].surrounds[y] = C[q1][q2][q3];//Приятно познакомиться!
						C[q1][q2][q3].normals[h].Normalize();

						for (int i = 0; i < DIM; i++)
						{
							C[qq1][qq2][qq3].normals[y].C[i] = -C[q1][q2][q3].normals[h].C[i];//Поделись нормалью
						}

						if (h < 6) C[q1][q2][q3].contact[h] = 1;		//Контакт на грани октаэдра
						else if (h < 14) C[q1][q2][q3].contact[h] = 3;	//Контакт на вершине октаэдра
						else C[q1][q2][q3].contact[h] = 2;				//Контакт на ребре
					}
					if (surround_count > 6)	//Уменьшение объёма из-за отсечений
					{
						double a = C[q1][q2][q3].size * 0.1;			//Длина срезанной части вдоль ребра
						double vol_edge = a*a*C[q1][q2][q3].size / 2.0;	//Объём, срезанный рёбрами
						double vol_vertex = a*a*a / SQRT3;				//Объём, срезанный вершинами
						int cut_edge = 0;		//Кол-во срезанных рёбер
						int cut_vertex = 0;		//Кол-во срезанных вершин
						for (int h = 6; h < surround_count; h++)
						{
							if (C[q1][q2][q3].contact[h] != 0)
							{
								if (h < 14) cut_vertex++;
								else cut_edge++;
							}
						}
						C[q1][q2][q3].volume -= (cut_edge*vol_edge + cut_vertex*vol_vertex);//Вычитание
					}

				}
			}
		}
	}

	void Polycrystall::SavePoleFig()
	{
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					GetPoleFig(&C[q1][q2][q3]);
					if (SST_SAVING) GetSST(&C[q1][q2][q3]);
				}
			}
		}
	}

	void Polycrystall::SaveDbgInfo()
	{
		for (int i = 0; i < file_count; i++)//Визуальное разделение шагов 
		{
			dbgstream[i] << "#########################      STEP " << CURR_STEP << "      #########################" << std::endl << std::endl;
		}

		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					WriteDebugInfo(dbgstream[0], C[q1][q2][q3].o.C);
					WriteDebugInfo(dbgstream[1], C[q1][q2][q3].e.C);
					WriteDebugInfo(dbgstream[2], C[q1][q2][q3].d.C);
					WriteDebugInfo(dbgstream[3], C[q1][q2][q3].sgm.C);
					WriteDebugInfo(dbgstream[4], C[q1][q2][q3].om.C);
					WriteDebugInfo(dbgstream[5], C[q1][q2][q3].dsgm.C);
					WriteDebugInfo(dbgstream[6], C[q1][q2][q3].d_in.C);
					WriteDebugInfo(dbgstream[7], C[q1][q2][q3].w.C);
					for (int f = 0; f < C[q1][q2][q3].SS_count; f++)
					{
						dbgstream[8] << C[q1][q2][q3].SS[f].dgm << " ";
					}
					dbgstream[8] << std::endl << std::endl;
					for (int f = 0; f < C[q1][q2][q3].SS_count; f++)
					{
						dbgstream[9] << C[q1][q2][q3].SS[f].t << " ";
					}
					dbgstream[9] << std::endl << std::endl;
					dbgstream[15] << C[q1][q2][q3].moment.C[0] << " " << C[q1][q2][q3].moment.C[1] << " " << C[q1][q2][q3].moment.C[2] << std::endl;
				}
			}
		}
		WriteDebugInfo(dbgstream[10], D.C);
		WriteDebugInfo(dbgstream[11], D_in.C);
		WriteDebugInfo(dbgstream[12], Sgm.C);
		WriteDebugInfo(dbgstream[13], dSgm.C);
		WriteDebugInfo(dbgstream[14], E.C);
	}

	void Polycrystall::Load(bool unload)
	{
		/*Параметр unload включает разгрузку представительного объёма*/
		if (REAL_UNIAX || unload)	//Одноосное растяжение
		{
			//Осреднение
			P.setZero();
			D_in.setZero();
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						D_in += C[q1][q2][q3].d_in;
						P += C[q1][q2][q3].p.ToLSK(C[q1][q2][q3].o);
					}
				}
			}

			D_in /= (total_fragm_count);
			P /= (total_fragm_count);

			//Симметризация
			P.Symmetrize();

			D = !unload ? TensionStrainCalc(P, D_in, D.C[0][0]) : UnloadingStrainCalc(P, D_in, Sgm, lam);

			Tensor b = D;	//Вычисление интенсивности деформаций
			b *= dt;
			E += b;
			Tensor buf = E;
			Strain = E.doubleScalMult(buf);
			Strain = SQRT2_3*sqrt(Strain);

			dSgm = TensionStressCalc(P, D_in, D);//Вычисление интенсивности напряжений
			dSgm *= dt;				//Приращение напряжений на шаге
			Sgm += dSgm;
			buf = Sgm;
			Stress = Sgm.doubleScalMult(buf);
			Stress = SQRT3_2*sqrt(Stress);
		}
		else
		{
			Stress = 0;		//Вычисление интенсивностей осреднением
			Strain = 0;
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						Strain += C[q1][q2][q3].strain;
						Stress += C[q1][q2][q3].stress;
					}
				}
			}
			Strain /= total_fragm_count;
			Stress /= total_fragm_count;
		}

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

					Tensor O = C[q1][q2][q3].o;
					Tensor OT = O;
					OT.Transp();
					C[q1][q2][q3].d = O*D*OT;//Гипотеза Фойгта
					C[q1][q2][q3].w = O*W*OT;//Расширенная

					C[q1][q2][q3].sgm = O*C[q1][q2][q3].sgm*OT;
					C[q1][q2][q3].d_in = O*C[q1][q2][q3].d_in*OT;


					/***************************************************
					***********       Пересчитываем НДС      ***********
					***************************************************/

					C[q1][q2][q3].NDScalc();

					if (HARDENING_BASE)			//Базовое упрочнение
					{
						Base_hardening(&C[q1][q2][q3]);
					}
				
					if (ROTATIONS_TAYLOR)		//Ротации по Тейлору
					{
						Taylor_rotations(&C[q1][q2][q3]);
					}
				
					if (ROTATIONS_TRUSOV && ROTATIONS_HARDENING)	//Ротационное упрочнение
					{
						Rotation_hardening(&C[q1][q2][q3]);
					}

					if (HARDENING_BOUND)	//Зернограничное упрочнение
					{
						Boundary_hardening(&C[q1][q2][q3]);
					}

					if (ROTATIONS_TRUSOV)		//Ротации по Трусову
					{
						Trusov_rotations(&C[q1][q2][q3]);
					}
					/**************************************************
					************       Переходим в ЛСК       **********
					**************************************************/

					C[q1][q2][q3].sgm = OT*C[q1][q2][q3].sgm*O;
					C[q1][q2][q3].d_in = OT*C[q1][q2][q3].d_in*O;
				}
			}
		}


		if (!REAL_UNIAX && !unload)		//Этот блок нужен исключительно для работы с энергией!
		{
			Sgm.setZero();
			D_in.setZero();
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						Sgm += C[q1][q2][q3].sgm;
						D_in += C[q1][q2][q3].d_in;
					}
				}
			}
			Sgm /= total_fragm_count;
			D_in /= total_fragm_count;
		}


		/************************************************************
		***********	        Прогресс выполнения 	      ***********
		************************************************************/

		double progress;
		if (!unload)
		{
			progress = Strain / strain_max * 100.0;

			if (!(cycle_count == 1 || cycle == 0))	//Для многоцикловых нагружений
			{
				progress /= 2.0;
				if (E.C[0][0] > 0)					//Этот код позволяет корректно
				{									//отображать прогресс выполнения,
					if (Sgm.C[0][0] > 0)			//когда на графике петли циклические
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
					if (Sgm.C[0][0] > 0)
					{
						progress = 50.0 - progress;
					}
					else
					{
						progress += 50.0;
					}
				}
			}
		}
		else
		{
			progress = final_stress / Stress * 100.0;	//Индикация при разгрузке
		}

		int period = unload ? proc_period/10 : proc_period;

		if (PROC_STEP == period)
		{
			PROC_STEP = 0;
			std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
			printf("        %2.2f", progress);
			std::cout << "%";
		}
		
		/************************************************************
		***********	    Запись данных для графиков НДС    ***********
		************************************************************/

		if ((progress - PLOT_STEP > plot_period || unload) && plot_period > 0)
		{
			if (!REAL_UNIAX && !unload)
			{
				Datastream[0].write((char *)&Strain, sizeof Strain);
				Datastream[1].write((char *)&Stress, sizeof Stress);
			}
			else
			{
				Datastream[0].write((char *)&E.C[0][0], sizeof E.C[0][0]);
				Datastream[1].write((char *)&Sgm.C[0][0], sizeof Sgm.C[0][0]);
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
						if (!REAL_UNIAX && !unload)
						{
							Datastream[2].write((char *)&C[q1][q2][q3].strain, sizeof C[q1][q2][q3].strain);
							Datastream[3].write((char *)&C[q1][q2][q3].stress, sizeof C[q1][q2][q3].stress);
						}
						else
						{
							Datastream[2].write((char *)&C[q1][q2][q3].e.C[0][0], sizeof C[q1][q2][q3].e.C[0][0]);
							Datastream[3].write((char *)&C[q1][q2][q3].sgm.C[0][0], sizeof C[q1][q2][q3].sgm.C[0][0]);
						}

						for (int i = 0; i < C[q1][q2][q3].SS_count; i++)
						{
							if (C[q1][q2][q3].SS[i].dgm > EPS) ActiveSysCount++;//Подсчёт активных СС
						}

						if (C[q1][q2][q3].isRotate) RotCount++;		//Подсчёт вращающихся решёток
						RotEnergy += C[q1][q2][q3].rot_energy;		//Суммирование энергий вращения
						RotSpeed += C[q1][q2][q3].rot_speed;		//Суммирование скоростей вращения
						norma += C[q1][q2][q3].norm;
						Mc += C[q1][q2][q3].rot_Mc;
						dmc += C[q1][q2][q3].dmc;

					}
				}
			}
			norma /= total_fragm_count;
			Mc /= total_fragm_count;
			dmc /= total_fragm_count;
			ActiveSysCount /= total_fragm_count;
			Datastream[4].write((char *)&ActiveSysCount, sizeof ActiveSysCount);//Запись кол-ва активных СС
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

			Tensor dE = D;
			dE *= dt;			//Приращение деформации на шаге

			double StepEnergy = Sgm.doubleScalMult(dE);	//Полная энергия на шаге

			TestStream[0] << RotCount << std::endl;
			TestStream[1] << RotSpeed << std::endl;
			TestStream[2] << RotEnergy << std::endl;
			TestStream[3] << StepEnergy << std::endl;
			TestStream[4] << Mc << std::endl;
			TestStream[5] << norma << std::endl;

			PLOT_STEP = progress;
		}

		/************************************************************
		***********	      Сохранение полюсных фигур	      ***********
		************************************************************/
		if (progress - POLUS_STEP > polus_period && polus_period > 0)
		{
			SavePoleFig();
			POLUS_STEP = progress;
		}

		/************************************************************
		***********	       Запись пошаговых данных	      ***********
		************************************************************/
		if (CURR_STEP >= DEBUG_START && CURR_STEP <= DEBUG_STOP && DEBUG_STEP == debug_period)
		{
			DEBUG_STEP = 0;
			SaveDbgInfo();
		}
		CURR_STEP++;
		PROC_STEP++;
		if (CURR_STEP >= DEBUG_START && CURR_STEP <= DEBUG_STOP) DEBUG_STEP++;

	}

	void Polycrystall::Fragmentate()
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
						int pos2 = C[q1][q2][q3].surrounds[h].position;//Позиция второго элемента
						//if (pos2 < pos1) continue;			//Раз матрица диагональная - нижнюю половину не нужно отдельно считать
						/*if (C[q1][q2][q3].contact[h] == 0)
						{
						sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = -1;//Если фрагменты не контактируют
						continue;
						}
						else if (sm_matrix[pos1][pos2]==0)//Если ещё не прошли
						{
						sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = (float)rand() / RAND_MAX;
						}*/
						if (pos1 < pos2) sm_matrix[pos1][pos2] = C[q1][q2][q3].DisorientMeasure(h);
						//sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = C[q1][q2][q3].contact[h];

						//if (pos1<pos2)sm_matrix[pos1][pos2] = C[q1][q2][q3].surrounds[h].position;
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
			delete[] sm_matrix[i];
		}
		delete[]sm_matrix;
	}
	void Polycrystall::Deformate()
	{
		omp_set_num_threads(thread_count);

		if (REAL_UNIAX)
		{
			tension_component = D.C[0][0];
		}
		for (cycle = 0; cycle < cycle_count; cycle++)
		{
			PLOT_STEP = 0;
			POLUS_STEP = 0;
			PROC_STEP = 0;
			if (cycle_count > 1)
			{
				std::cout << std::endl << " Cycle # " << cycle + 1;
			}
			std::cout << std::endl << "        0.00%";
			while (Strain < strain_max)
			{
				Load(false);
			}
		
			if (UNLOADING)
			{
				std::cout << std::endl << " Unloading # " << cycle + 1;
				std::cout << std::endl << "        0.00%";
				while (fabs(Stress) > final_stress)
				{
					 Load(true);
				}
			}
			if (cycle_count > 1)
			{
				D.C[0][0] = pow(-1, cycle + 1) * tension_component;	//Меняем знак растягивающей компоненты
				strain_max += strain_max * addition_strain;					//Повышаем предел интенсивности
			}
			if (FRAGMENTATION) Fragmentate();
		
		}
	}
}