// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include "params.h"
#include "rotations.h"
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
			int i1, i2, i3;
			get3DPos(G[i].center, &i1, &i2, &i3);//Получаем координаты угла зерна
			for (int q1 = i1; q1 < i1 + G[i].size; q1++)
			{
				for (int q2 = i2; q2 < i2 + G[i].size; q2++)
				{
					for (int q3 = i3; q3 < i3 + G[i].size; q3++)
					{
						
					}
				}
			}
		}

	}

	void Polycrystall::GrainRotate()
	{
		//Бежим по всем внешним фасеткам зерна и вычисляем несовместность
		//Здесь только большие фасетки учитываем
		for (int i = 0; i < G.size(); i++)
		{
			if (G[i].size == 1) continue; //Единичные зерна сами по себе вращаются
			int i1, i2, i3;
			get3DPos(G[i].center, &i1, &i2, &i3);//Получаем координаты угла зерна
			Tensor GrLp;	//Тензор скачка деформаций
			Tensor d_in1[6], d_in2[6];
			Vector normals[6];
			//6 площадок
			//0 - вверх, 1 - от нас, 2 - вправо, 3 - на нас, 4 - влево, 5 - вниз
			Vector dM;
			for (int qq1 = i1; qq1 < i1 + G[i].size; qq1++)
			{
				for (int qq2 = i2; qq2 < i2 + G[i].size; qq2++)
				{
					int q1 = qq1 > fragm_count - 1 ? qq1 - fragm_count : qq1;
					int q2 = qq2 > fragm_count - 1 ? qq2 - fragm_count : qq2;
					

					normals[0] += C[q1][q2][i3].normals[5];
					normals[1] += C[q1][q2][i3 + G[i].size].normals[0];
					/*for (int i = 0; i < DIM; i++)
					{
						for (int j = 0; j < DIM; j++)
						{
							//Вниз
							for (int k = 0; k < C[q1][q2][i3].SS_count; k++)
							{
							
								if (C[q1][q2][i3].SS[k].b.ScalMult(C[q1][q2][i3].normals[5]) < 0) continue; //Скольжение от границы - вклад не вносится
								d_in1[0].C[i][j] += C[q1][q2][i3].SS[k].dgm * (C[q1][q2][i3].SS[k].n.C[i] * C[q1][q2][i3].SS[k].b.C[j]);
								
							}
							for (int k = 0; k < C[q1][q2][i3].surrounds[5].SS_count; k++)
							{
								d_in2[0].C[i][j] += C[q1][q2][i3].surrounds[5].SS[k].dgm * (C[q1][q2][i3].surrounds[5].SS[k].n.C[i] * C[q1][q2][i3].surrounds[5].SS[k].b.C[j]);
							}
							//Вверх
							int k3 = i3 + G[i].size;
							for (int k = 0; k < C[q1][q2][k3].SS_count; k++)
							{

								if (C[q1][q2][k3].SS[k].b.ScalMult(C[q1][q2][k3].normals[0]) < 0) continue; //Скольжение от границы - вклад не вносится
								d_in1[1].C[i][j] += C[q1][q2][k3].SS[k].dgm * (C[q1][q2][k3].SS[k].n.C[i] * C[q1][q2][k3].SS[k].b.C[j]);

							}
							for (int k = 0; k < C[q1][q2][k3].surrounds[0].SS_count; k++)
							{
								d_in2[1].C[i][j] += C[q1][q2][k3].surrounds[0].SS[k].dgm * (C[q1][q2][k3].surrounds[0].SS[k].n.C[i] * C[q1][q2][k3].surrounds[0].SS[k].b.C[j]);
							}
							
						}
					}*/
					
				}
			}

		/*	for (int q1 = i1; q1 < i1 + G[i].size; q1++)
			{
				for (int q3 = i3; q3 < i3 + G[i].size; q3++)
				{
					normals[2] += C[q1][i2][q3].normals[4];
					normals[3] += C[q1][i2 + G[i].size][q3].normals[2];
					for (int i = 0; i < DIM; i++)
					{
						for (int j = 0; j < DIM; j++)
						{
							//Влево
							for (int k = 0; k < C[q1][i2][q3].SS_count; k++)
							{

								if (C[q1][i2][q3].SS[k].b.ScalMult(C[q1][i2][q3].normals[4]) < 0) continue; //Скольжение от границы - вклад не вносится
								d_in1[2].C[i][j] += C[q1][i2][q3].SS[k].dgm * (C[q1][i2][q3].SS[k].n.C[i] * C[q1][i2][q3].SS[k].b.C[j]);

							}
							for (int k = 0; k < C[q1][i2][q3].surrounds[4].SS_count; k++)
							{
								d_in2[2].C[i][j] += C[q1][i2][q3].surrounds[4].SS[k].dgm * (C[q1][i2][q3].surrounds[4].SS[k].n.C[i] * C[q1][i2][q3].surrounds[4].SS[k].b.C[j]);
							}
							//Вправо
							int k2 = i2 + G[i].size;
							for (int k = 0; k < C[q1][k2][i3].SS_count; k++)
							{

								if (C[q1][k2][i3].SS[k].b.ScalMult(C[q1][k2][i3].normals[2]) < 0) continue; //Скольжение от границы - вклад не вносится
								d_in1[3].C[i][j] += C[q1][k2][i3].SS[k].dgm * (C[q1][k2][i3].SS[k].n.C[i] * C[q1][k2][i3].SS[k].b.C[j]);

							}
							for (int k = 0; k < C[q1][k2][i3].surrounds[2].SS_count; k++)
							{
								d_in2[3].C[i][j] += C[q1][k2][i3].surrounds[2].SS[k].dgm * (C[q1][k2][i3].surrounds[2].SS[k].n.C[i] * C[q1][k2][i3].surrounds[2].SS[k].b.C[j]);
							}

						}
					}

				}
			}
			for (int q2 = i1; q2 < i2 + G[i].size; q2++)
			{
				for (int q3 = i3; q3 < i3 + G[i].size; q3++)
				{
					normals[4] += C[i1][q2][q3].normals[1];
					normals[5] += C[i1 + G[i].size][q2][q3].normals[3];
					for (int i = 0; i < DIM; i++)
					{
						for (int j = 0; j < DIM; j++)
						{
							//От нас
							for (int k = 0; k < C[i1][q2][q3].SS_count; k++)
							{

								if (C[i1][q2][q3].SS[k].b.ScalMult(C[i1][q2][q3].normals[1]) < 0) continue; //Скольжение от границы - вклад не вносится
								d_in1[4].C[i][j] += C[i1][q2][q3].SS[k].dgm * (C[i1][q2][q3].SS[k].n.C[i] * C[i1][q2][q3].SS[k].b.C[j]);

							}
							for (int k = 0; k < C[i1][q2][q3].surrounds[1].SS_count; k++)
							{
								d_in2[4].C[i][j] += C[i1][q2][q3].surrounds[1].SS[k].dgm * (C[i1][q2][q3].surrounds[1].SS[k].n.C[i] * C[i1][q2][q3].surrounds[1].SS[k].b.C[j]);
							}
							//На нас
							int k1 = i1 + G[i].size;
							for (int k = 0; k < C[k1][q2][q3].SS_count; k++)
							{

								if (C[k1][q2][q3].SS[k].b.ScalMult(C[k1][q2][q3].normals[3]) < 0) continue; //Скольжение от границы - вклад не вносится
								d_in1[5].C[i][j] += C[k1][q2][q3].SS[k].dgm * (C[k1][q2][q3].SS[k].n.C[i] * C[k1][q2][q3].SS[k].b.C[j]);

							}
							for (int k = 0; k < C[k1][q2][q3].surrounds[3].SS_count; k++)
							{
								d_in2[5].C[i][j] += C[k1][q2][q3].surrounds[3].SS[k].dgm * (C[k1][q2][q3].surrounds[3].SS[k].n.C[i] * C[k1][q2][q3].surrounds[3].SS[k].b.C[j]);
							}

						}
					}

				}
			}*/

		/*	double S = prms::fragm_size_m*prms::fragm_size_m*G[i].size*G[i].size;//Площадь крупной фасетки зерна
			int volume = C[i1][i2][i3].volume*pow(G[i].size, 3);//Объем зерна
			for (int p = 0; p < 6; p++)
			{
				normals[p].Normalize();
				Tensor Lp = d_in1[p] - d_in2[p];
				Lp.Transp();
				Tensor buf = VectMult(normals[p], Lp);
				Vector dm = ScalMult(buf, normals[p]);//Поверхностный вектор-момент
				dm *= prms::ROT_L;
				dM += dm*S;
			}
			
			dM /= volume;
			double dMnorm = dM.getNorm();*/
			Vector M = /*f->moment + */dM*prms::dt;//Пока что не накапливается
			//f->moment = M;
			/*double norm = M.getNorm();
			
			if (norm > prms::ROT_MC || norm == -1)
			{
				norm = prms::ROT_MC;
			}

			double pr = M.ScalMult(dM);

			double dFi = 0;		//Скорость вращения
			if (norm == prms::ROT_MC && pr >= 0)	//Пластические и упругие развороты
			{
				dFi = prms::ROT_A * dMnorm + prms::ROT_H * norm;
			}
			else
			{
				dFi = prms::ROT_A * dMnorm;		//Только упругие развороты
			}

			Vector e = M;
			Tensor OM;
			if (dFi > EPS*1e4)
			{
				e.Normalize();
				dFi *= prms::dt;
			
				for (int i = 0; i < DIM; i++)
				{
					for (int j = 0; j < DIM; j++)
					{
						for (int k = 0; k < DIM; k++)
						{
							OM.C[i][j] -= LeviCivit(i, j, k) * e.C[k] * dFi;
						}
					}
				}

			}*/
		
			/*for (int q1 = i1; q1 < i1 + G[i].size; q1++)
			{
				for (int q2 = i2; q2 < i2 + G[i].size; q2++)
				{
					for (int q3 = i3; q3 < i3 + G[i].size; q3++)
					{
						C[q1][q2][q3].om += OM;
						Rotate(&C[q1][q2][q3], dFi, e);
					}
				}
			}*/

		}
	}

	void Polycrystall::MakeGrains()
	{
		int gsz = prms::Grain_size;					//Желаемый размер зерна (во фрагментах на ребере)
		int cnt = 2*pow(int(fragm_count / gsz), 3);	//Желаемое кол-во зёрен
		
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

		for (int i = 0; i < cnt; i++)
		{

			bool good = false;
			int rnd, rnd1, rnd2, rnd3;
			while (!good)
			{
				rnd1 = rand() % fragm_count;
				rnd2 = rand() % fragm_count;
				rnd3 = rand() % fragm_count;
				rnd = get1DPos(rnd1, rnd2, rnd3);	//Случайным образом выбирается центр
				good = true;
				for (int j = 0; j < G.size(); j++)
				{
					if (mass[rnd1][rnd2][rnd3] != -1)//Исключение наложений
					{
						good = false;
						break;
					}
				}

			}
			Grain new_gr;
			new_gr.center = rnd;					//Назначили центр зерна
			new_gr.size = rand() % (2 * (gsz - 2) + 1) + 2;	//Назначили размер зерна
			new_gr.num = i;
			G.push_back(new_gr);
			double a = ((double)rand() / RAND_MAX) * (PI);//Общая ориентация для данного зерна
			double g = ((double)rand() / RAND_MAX) * (PI);
			double y1 = ((double)rand() / RAND_MAX);
			double y2 = ((double)rand() / RAND_MAX);
			for (int q1 = rnd1; q1 < rnd1 + G[i].size; q1++)
			{
				for (int q2 = rnd2; q2 < rnd2 + G[i].size; q2++)
				{
					for (int q3 = rnd3; q3 < rnd3 + G[i].size; q3++)
					{
						int i1 = q1 > fragm_count - 1 ? q1 - fragm_count : q1;
						int i2 = q2 > fragm_count - 1 ? q2 - fragm_count : q2;
						int i3 = q3 > fragm_count - 1 ? q3 - fragm_count : q3;
						double delta1 = ((double)rand() / RAND_MAX) * (PI / 60);
						double delta2 = ((double)rand() / RAND_MAX) * (PI / 60);
						C[i1][i2][i3].Orientate(a + delta1, g + delta2, y1, y2);
						mass[i1][i2][i3] = i;
					}
				}
			}

		}

		for (int q1 = 0; q1 < fragm_count; q1++)//Обработка фрагментов, не вошедших в состав зерен
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					if (mass[q1][q2][q3] == -1)
					{
						double a = ((double)rand() / RAND_MAX) * (PI);
						double g = ((double)rand() / RAND_MAX) * (PI);
						double y1 = ((double)rand() / RAND_MAX);
						double y2 = ((double)rand() / RAND_MAX);
						
						Grain new_gr;
						new_gr.center = get1DPos(q1, q2, q3);
						//Попытка образовать зерно размером 2
						if (q1 < fragm_count - 1 && mass[q1 + 1][q2][q3] == -1 &&
							q2 < fragm_count - 1 && mass[q1][q2 + 1][q3] == -1 &&
							q3 < fragm_count - 1 && mass[q1][q2][q3 + 1] == -1)
						{
							for (int qq1 = q1; qq1 < q1 + 1; qq1++)
							{
								for (int qq2 = q2; qq2 < q2 + 1; qq2++)
								{
									for (int qq3 = q3; qq3 < q3 + 1; qq3++)
									{
										mass[qq1][qq2][qq3] = cnt++;
										double delta1 = ((double)rand() / RAND_MAX) * (PI / 60);
										double delta2 = ((double)rand() / RAND_MAX) * (PI / 60);
										C[qq1][qq2][qq3].Orientate(a+delta1, g+delta2, y1, y2);
									}
								}
							}
							new_gr.size = 2;
							new_gr.num = mass[q1][q2][q3];
						}
						else//В противном случае образуются единичные зерна
						{
							C[q1][q2][q3].Orientate(a, g, y1, y2);
							mass[q1][q2][q3] = cnt++;
							new_gr.size = 1;
							new_gr.num = mass[q1][q2][q3];
						}
											
						G.push_back(new_gr);
					}
				}
			}
		}

		Save_Sort_Size();//Вывод в файл данных о зернах
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