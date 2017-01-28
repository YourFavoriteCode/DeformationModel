#include "stdafx.h"
#include "params.h"
#include <fstream>

#include "Polycrystall.h"

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

	void Polycrystall::Fragmentate()
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
						if (pos1 < pos2) sm_matrix[pos1][pos2] = C[q1][q2][q3].DisorientMeasure(h);
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

	void Polycrystall::MakeGrains()
	{
		int gsz = prms::Grain_size;	//Желаемый размер зерна (во фрагментах на ребере)
		int cnt = pow(int(fragm_count / gsz), 3);	//Желаемое кол-во зёрен
		int*** mass = new int**[fragm_count];		//Массив принадлежности фрагментов каким-либо зёрнам
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

		
		int* cryst_center = new int[cnt];
		for (int i = 0; i < cnt; i++)
		{
			bool good = false;
			int rnd, rnd1, rnd2, rnd3;
			while (!good)
			{
				//rnd = rand() % total_fragm_count;
				rnd1 = rand() % fragm_count;
				rnd2 = rand() % fragm_count;
				rnd3 = rand() % fragm_count;
				rnd = get1DPos(rnd1, rnd2, rnd3);
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
		//cryst_center[0] = 1;
		//cryst_center[1] = 26;
		//cryst_center[2] = 98;
		//cryst_center[3] = 164;
		//cryst_center[4] = 201;
		for (int i = 0; i < cnt; i++)//Сортировка для оптимального поиска
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
		for (int i = 0; i < cnt; i++)
		{
			printf(" Center %d: %d\n", i, cryst_center[i]);
			int q1, q2, q3;
			get3DPos(cryst_center[i], &q1, &q2, &q3);
			C[q1][q2][q3].crystall_center = true;
		}

		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{

					for (int i = 0; i < cnt; i++)
					{
						if (get1DPos(q1, q2, q3) == cryst_center[i])//попали в центр кристаллизации
						{
							int dsp = 2*(gsz - 2) + 1;//Диапазон вариации размеров (минимальный размер зерна равен 2)
							int grain_size = rand() % dsp + 2;//Вариация размеров зёрен
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
										int i1 = qq1 > fragm_count - 1 ? qq1 - fragm_count : qq1;
										int i2 = qq2 > fragm_count - 1 ? qq2 - fragm_count : qq2;
										int i3 = qq3 > fragm_count - 1 ? qq3 - fragm_count : qq3;
										double delta1 = ((double)rand() / RAND_MAX) * (PI / 60);
										double delta2 = ((double)rand() / RAND_MAX) * (PI / 60);
										C[i1][i2][i3].Orientate(a - delta1, g - delta2, y1-delta1, y2-delta2);//Малые разориентации
										mass[i1][i2][i3] = i;//Данный фрагмент теперь принадлежит зерну i
									}
								}
							}
							break;
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
						printf(" fuck [%d_%d_%d]\n", q1, q2, q3);
						double a = ((double)rand() / RAND_MAX) * (PI);
						double g = ((double)rand() / RAND_MAX) * (PI);
						double y1 = ((double)rand() / RAND_MAX);
						double y2 = ((double)rand() / RAND_MAX);
						C[q1][q2][q3].Orientate(a, g, y1, y2);
						mass[q1][q2][q3] = cnt++;
					}
				}
			}
		}

		delete[] cryst_center;

		for (int i = 0; i < fragm_count; i++)
		{
			for (int j = 0; j < fragm_count; j++)
			{
				delete[] mass[i][j];
			}
			delete[] mass[i];
		}
		delete[] mass;
	}


}