// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include "params.h"
#include <fstream>
#include <vector>

#include "Polycrystall.h"
#include "Fragmentation.h"

namespace model
{
	int get1DPos(int q1, int q2, int q3)
	{
		//По трём координатам в объёме поликристалла возвращает уникальный номер фрагмента
		int res = q1*prms::fragm_count*prms::fragm_count + q2*prms::fragm_count + q3;
		return res;
	}

	void get3DPos(int pos, int* q1, int* q2, int* q3)
	{
		//Восстанавливает пространственные координаты поликристалла по уникальному номеру
		int C2d = prms::fragm_count*prms::fragm_count;
		int C3d = C2d*prms::fragm_count;

		int qq1 = pos / C2d;
		int qq2 = (pos - qq1*C2d) / prms::fragm_count;
		int qq3 = (pos - qq1*C2d) % prms::fragm_count;

		*q1 = qq1;
		*q2 = qq2;
		*q3 = qq3;
	}

	int*** mass;				//Массив принадлежности фрагментов определенным зернам
	std::vector<Grain> G;		//Массив зерен

	void Save_Sort_Size()
	{
		std::vector<int> G_Size, G_Count;
		for (int i = 0; i < G.size(); i++)//Сортировка
		{
			for (int j = 0; j < G.size()-1; j++)
			{
				if (G[j + 1].size < G[j].size)//По возрастанию
				{
					int buf = G[j + 1].size;
					G[j + 1].size = G[j].size;
					G[j].size = buf;
				}
			}
		}
		//Нужно повторяющиеся просуммировать
		int l = 0;
		G_Size.push_back(G[0].size);
		G_Count.push_back(1);
		for (int i = 0; i < G.size()-1; i++)//Сортировка
		{
			if (G[i].size != G[i + 1].size)//Нет повтора
			{
				++l;
				G_Size.push_back(G[i + 1].size);
				G_Count.push_back(1);
			}
			else//Повтор
			{
				G_Count[l]++;
			}
		}
		//Сохранение в файлы
		FILE* G_File=fopen("Plot\\GrainsSize.txt","w");
		for (int i = 0; i < G_Size.size(); i++)//Сортировка
		{
			fprintf(G_File, "%d %d ", G_Size[i], G_Count[i]);
		}
		fclose(G_File);
	}

	void Polycrystall::Fragmentate()
	{
		/*
		Поиск новых большеугловых границ ведется
		исключительно внутри ранее образованных зерен
		*/
		for (int i = 0; i < G.size(); i++)
		{
			if (G[i].size == 1) continue;		//Не ищем внутри единичных фрагментов

		}

	}

	void Polycrystall::MakeGrains()
	{
		int gsz = prms::Grain_size;					//Желаемый размер зерна (во фрагментах на ребере)
		int cnt = pow(int(fragm_count / gsz), 3);	//Желаемое кол-во зёрен
		mass = new int**[fragm_count];
		for (int i = 0; i < fragm_count; i++)
		{
			mass[i] = new int*[fragm_count];
			for (int j = 0; j < fragm_count; j++)
			{
				mass[i][j] = new int[fragm_count];
			}
		}
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					mass[q1][q2][q3] = -1;//Вначале все фрагменты никому не принадлежат
				}
			}
		}

		int* cryst_center = new int[cnt];		//Массив центров кристаллизации зерен
		for (int i = 0; i < cnt; i++)
		{
			bool good = false;
			int rnd;
			while (!good)
			{
				int rnd1 = rand() % fragm_count;
				int rnd2 = rand() % fragm_count;
				int rnd3 = rand() % fragm_count;
				rnd = get1DPos(rnd1, rnd2, rnd3);	//Случайным образом выбирается центр
				good = true;
				for (int j = 0; j < i; j++)
				{
					if (rnd == cryst_center[j])//Исключение возможных повторений
					{
						good = false;
						break;
					}
				}
			}
			cryst_center[i] = rnd;
		}

		for (int i = 0; i < cnt; i++)//Сортировка для оптимального обхода массива
		{
			for (int j = 0; j < cnt - 1; j++)
			{
				if (cryst_center[j] > cryst_center[j + 1])
				{
					int buf = cryst_center[j + 1];
					cryst_center[j + 1] = cryst_center[j];
					cryst_center[j] = buf;
				}
			}
		}

		for (int i = 0; i < cnt; i++)		//Вывод отладочных данных
		{
			printf("\n Center %d: %d\n", i, cryst_center[i]);
			int q1, q2, q3;
			get3DPos(cryst_center[i], &q1, &q2, &q3);
			C[q1][q2][q3].crystall_center = true;
		}

		int f = 0;	//Сколько зерен уже было обнаружено
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{

				//	if (mass[q1][q2][q3] != -1) continue;	//Фрагмент уже принадлежит какому-то зерну

					for (int i = f; i < cnt; i++)
					{
						if (get1DPos(q1, q2, q3) == cryst_center[i])//Попали в центр кристаллизации
						{
							int dsp = 2 * (gsz - 2) + 1;//Диапазон вариации размеров (минимальный размер зерна равен 2)
							int grain_size = rand() % dsp + 2;//Вариация размеров зёрен
							Grain new_gr;
							new_gr.center = cryst_center[i];
							new_gr.size = grain_size;
							new_gr.num = i;
							G.push_back(new_gr);

							double a = ((double)rand() / RAND_MAX) * (PI);//Общая ориентация для данного зерна
							double g = ((double)rand() / RAND_MAX) * (PI);
							double y1 = ((double)rand() / RAND_MAX);
							double y2 = ((double)rand() / RAND_MAX);

							for (int qq1 = q1; qq1 < q1 + grain_size; qq1++)
							{
								for (int qq2 = q2; qq2 < q2 + grain_size; qq2++)
								{
									for (int qq3 = q3; qq3 < q3 + grain_size; qq3++)
									{
									//	if (mass[qq1][qq2][qq3] != -1) continue;
										int i1 = qq1 > fragm_count - 1 ? qq1 - fragm_count : qq1;
										int i2 = qq2 > fragm_count - 1 ? qq2 - fragm_count : qq2;
										int i3 = qq3 > fragm_count - 1 ? qq3 - fragm_count : qq3;
										double delta1 = ((double)rand() / RAND_MAX) * (PI / 60);
										double delta2 = ((double)rand() / RAND_MAX) * (PI / 60);
										C[i1][i2][i3].Orientate(a - delta1, g - delta2, y1, y2);//Малые разориентации
										mass[i1][i2][i3] = i;//Данный фрагмент теперь принадлежит зерну i
									}
								}
							}
							f++;
							break;	//Раз нашли, то нет смысла проверять другие
						}

					}


				}
			}
		}

		for (int q1 = 0; q1 < fragm_count; q1++)//Оставшиеся ни при делах просто становятся отдельными зёрнами
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					if (mass[q1][q2][q3] == -1)
					{
						printf(" Single fragm [%d_%d_%d]\n", q1, q2, q3);
						double a = ((double)rand() / RAND_MAX) * (PI);
						double g = ((double)rand() / RAND_MAX) * (PI);
						double y1 = ((double)rand() / RAND_MAX);
						double y2 = ((double)rand() / RAND_MAX);
						C[q1][q2][q3].Orientate(a, g, y1, y2);
						mass[q1][q2][q3] = f++;

						Grain new_gr;
						new_gr.center = get1DPos(q1, q2, q3);
						new_gr.size = 1;
						new_gr.num = mass[q1][q2][q3];
						G.push_back(new_gr);
					}
				}
			}
		}

		delete[] cryst_center;

		/* (int i = 0; i < fragm_count; i++)
		{
		for (int j = 0; j < fragm_count; j++)
		{
		delete[] mass[i][j];
		}
		delete[] mass[i];
		}
		delete[] mass;*/
		Save_Sort_Size();
	}

	void Polycrystall::Illustrate()
	{
		float** sm_matrix = new float*[total_fragm_count];	//Весовая матрица
		for (int i = 0; i < total_fragm_count; i++)
		{
			sm_matrix[i] = new float[total_fragm_count];
		}

		for (int i = 0; i < total_fragm_count; i++)
		{
			for (int j = 0; j < total_fragm_count; j++)
			{
				sm_matrix[i][j] = -1;	//"Зануление"
			}
		}

		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					int pos1 = get1DPos(q1, q2, q3);	//Позиция первого элемента
					if (C[q1][q2][q3].crystall_center)
					{
						sm_matrix[pos1][pos1] = -2;		//Обозначение центра кристаллизации зерна
					}
					for (int h = 0; h < prms::surround_count; h++)
					{
						if (C[q1][q2][q3].contact[h] == 0) continue;//Если фрагменты не контактируют
						int pos2 = C[q1][q2][q3].surrounds[h].position;//Позиция второго элемента
						if (pos1 < pos2) sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = C[q1][q2][q3].DisorientMeasure(h);
						//sm_matrix[pos1][pos2] = sm_matrix[pos2][pos1] = C[q1][q2][q3].contact[h];

					}


				}
			}
		}
		//Сохранение весовой матрицы в файл
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
		//Очистка памяти


		for (int i = 0; i < total_fragm_count; i++)
		{
			delete[] sm_matrix[i];
		}
		delete[]sm_matrix;
	}
}