#include "stdafx.h"
#include <cmath>

#include "fragment.h"
#include "params.h"
#include "materials.h"

namespace model
{

	Fragment::Fragment()
	{
		//Обнуление всех переменных при инициализации
		//кроме тех, которые имеют свой конструктор
		stress = 0;
		strain = 0;
		SS_count = 0;
		size = 0;
		rot_speed = 0;
		isRotate = false;
		sum_angle = 0;
		mc = 0;
		rot_energy = 0;
	}

	Fragment::~Fragment()
	{
	
	}

	int inline H(double a)		//Функция Хэвисайда
	{
		return (a >= 0 ? 1 : 0);
	}

	void Fragment::Orientate(double a, double g,
		double y1, double y2)
	{
		/*
		* Функция задаёт ориентацию решётки на основе двух
		* углов и двух случайных чисел, генерирую равномерное 
		* распределение косинуса третьего угла
		*/
		double cb = y1 > 0.5 ? y2 : -y2;
		double sb = sqrt(1.0 - cb * cb);

		double sa = sin(a);
		double sg = sin(g);
		double ca = cos(a);
		double cg = cos(g);

		o.C[0][0] = ca * cg - sa * cb * sg;
		o.C[0][1] = -ca * sg - sa * cb * cg;
		o.C[0][2] = sa * sb;

		o.C[1][0] = sa * cg + ca * cb * sg;
		o.C[1][1] = -sa * sg + ca * cb * cg;
		o.C[1][2] = -ca * sb;

		o.C[2][0] = sb * sg;
		o.C[2][1] = sb * cg;
		o.C[2][2] = cb;

	}

	void Fragment::setMaterialParams(int material_type)
	{
		/*
		* Задание всех материальных параметров фрагмента
		* в зависимости от выбранного материала
		*/
		this->material = material_type;
		double P1, P2, P3;
		switch (material_type)
		{
		case 0://Сталь 45
		{
			SS_count = SS_COUNT_BCC;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < 24; i++)
			{
				SS[i].tc = STEEL_TC1;
			}
			for (int i = 24; i < 48; i++)
			{
				SS[i].tc = STEEL_TC2;
			}
			for (int i = 48; i < 96; i++)
			{
				SS[i].tc = STEEL_TC3;
			}
			P1 = STEEL_P1;
			P2 = STEEL_P2;
			P3 = STEEL_P3;

			
			//-----------[110]
			SS[0].Initialize(0, 1, 1, 1, -1, 1);
			SS[1].Initialize(0, 1, 1, 1, 1, -1);
			SS[2].Initialize(0, -1, 1, 1, 1, 1);
			SS[3].Initialize(0, -1, 1, -1, 1, 1);

			SS[4].Initialize(-1, 1, 0, 1, 1, 1);
			SS[5].Initialize(-1, 1, 0, 1, 1, -1);
			SS[6].Initialize(1, 1, 0, -1, 1, 1);
			SS[7].Initialize(1, 1, 0, 1, -1, 1);

			SS[8].Initialize(1, 0, 1, 1, 1, -1);
			SS[9].Initialize(1, 0, 1, -1, 1, 1);
			SS[10].Initialize(-1, 0, 1, 1, 1, 1);
			SS[11].Initialize(-1, 0, 1, 1, -1, 1);

			SS[12].Initialize(0, 1, 1, -1, 1, -1);
			SS[13].Initialize(0, 1, 1, -1, -1, 1);
			SS[14].Initialize(0, -1, 1, -1, -1, -1);
			SS[15].Initialize(0, -1, 1, 1, -1, -1);

			SS[16].Initialize(-1, 1, 0, -1, -1, -1);
			SS[17].Initialize(-1, 1, 0, -1, -1, 1);
			SS[18].Initialize(1, 1, 0, 1, -1, -1);
			SS[19].Initialize(1, 1, 0, -1, 1, -1);

			SS[20].Initialize(1, 0, 1, -1, -1, 1);
			SS[21].Initialize(1, 0, 1, 1, -1, -1);
			SS[22].Initialize(-1, 0, 1, -1, -1, -1);
			SS[23].Initialize(-1, 0, 1, -1, 1, -1);
			//-----------[112]
			SS[24].Initialize(2, 1, 1, 1, -1, -1);
			SS[25].Initialize(2, -1, 1, 1, 1, -1);
			SS[26].Initialize(2, -1, -1, 1, 1, 1);
			SS[27].Initialize(2, 1, -1, 1, -1, 1);

			SS[28].Initialize(1, 2, 1, -1, 1, -1);
			SS[29].Initialize(-1, 2, 1, 1, 1, -1);
			SS[30].Initialize(-1, 2, -1, 1, 1, 1);
			SS[31].Initialize(1, 2, -1, -1, 1, 1);

			SS[32].Initialize(1, 1, 2, -1, -1, 1);
			SS[33].Initialize(-1, 1, 2, 1, -1, 1);
			SS[34].Initialize(-1, -1, 2, 1, 1, 1);
			SS[35].Initialize(1, -1, 2, -1, 1, 1);

			SS[36].Initialize(2, 1, 1, -1, 1, 1);
			SS[37].Initialize(2, -1, 1, -1, -1, 1);
			SS[38].Initialize(2, -1, -1, -1, -1, -1);
			SS[39].Initialize(2, 1, -1, -1, 1, -1);

			SS[40].Initialize(1, 2, 1, 1, -1, 1);
			SS[41].Initialize(-1, 2, 1, -1, -1, 1);
			SS[42].Initialize(-1, 2, -1, -1, -1, -1);
			SS[43].Initialize(1, 2, -1, 1, -1, -1);

			SS[44].Initialize(1, 1, 2, 1, 1, -1);
			SS[45].Initialize(-1, 1, 2, -1, 1, -1);
			SS[46].Initialize(-1, -1, 2, -1, -1, -1);
			SS[47].Initialize(1, -1, 2, 1, -1, -1);
			//-----------[123]
			SS[48].Initialize(1, 2, 3, -1, -1, 1);
			SS[49].Initialize(-1, 2, 3, 1, -1, 1);
			SS[50].Initialize(-1, -2, 3, 1, 1, 1);
			SS[51].Initialize(1, -2, 3, -1, 1, 1);

			SS[52].Initialize(3, 1, 2, 1, -1, -1);
			SS[53].Initialize(3, -1, 2, 1, 1, -1);
			SS[54].Initialize(3, -1, -2, 1, 1, 1);
			SS[55].Initialize(3, 1, -2, 1, -1, 1);

			SS[56].Initialize(2, 3, 1, -1, 1, -1);
			SS[57].Initialize(2, 3, -1, -1, 1, 1);
			SS[58].Initialize(-2, 3, -1, 1, 1, 1);
			SS[59].Initialize(-2, 3, 1, 1, 1, -1);

			SS[60].Initialize(2, 1, 3, -1, -1, 1);
			SS[61].Initialize(-2, 1, 3, 1, -1, 1);
			SS[62].Initialize(-2, -1, 3, 1, 1, 1);
			SS[63].Initialize(2, -1, 3, -1, 1, 1);

			SS[64].Initialize(1, 3, 2, -1, 1, -1);
			SS[65].Initialize(-1, 3, 2, 1, 1, -1);
			SS[66].Initialize(-1, 3, -2, 1, 1, 1);
			SS[67].Initialize(1, 3, -2, -1, 1, 1);

			SS[68].Initialize(3, 2, 1, 1, -1, -1);
			SS[69].Initialize(3, -2, 1, 1, 1, -1);
			SS[70].Initialize(3, -2, -1, 1, 1, 1);
			SS[71].Initialize(3, 2, -1, 1, -1, 1);

			SS[72].Initialize(1, 2, 3, 1, 1, -1);
			SS[73].Initialize(-1, 2, 3, -1, 1, -1);
			SS[74].Initialize(-1, -2, 3, -1, -1, -1);
			SS[75].Initialize(1, -2, 3, 1, -1, -1);

			SS[76].Initialize(3, 1, 2, -1, 1, 1);
			SS[77].Initialize(3, -1, 2, -1, -1, 1);
			SS[78].Initialize(3, -1, -2, -1, -1, -1);
			SS[79].Initialize(3, 1, -2, -1, 1, -1);

			SS[80].Initialize(2, 3, 1, 1, -1, 1);
			SS[81].Initialize(2, 3, -1, 1, -1, -1);
			SS[82].Initialize(-2, 3, -1, -1, -1, -1);
			SS[83].Initialize(-2, 3, 1, -1, -1, 1);

			SS[84].Initialize(2, 1, 3, 1, 1, -1);
			SS[85].Initialize(-2, 1, 3, -1, 1, -1);
			SS[86].Initialize(-2, -1, 3, -1, -1, -1);
			SS[87].Initialize(2, -1, 3, 1, -1, -1);

			SS[88].Initialize(1, 3, 2, 1, -1, 1);
			SS[89].Initialize(-1, 3, 2, -1, -1, 1);
			SS[90].Initialize(-1, 3, -2, -1, -1, -1);
			SS[91].Initialize(1, 3, -2, 1, -1, -1);

			SS[92].Initialize(3, 2, 1, -1, 1, 1);
			SS[93].Initialize(3, -2, 1, -1, -1, 1);
			SS[94].Initialize(3, -2, -1, -1, -1, -1);
			SS[95].Initialize(3, 2, -1, -1, 1, -1);

			break;
		}
		case 1://Медь
		{
			SS_count = SS_COUNT_FCC;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < SS_count; i++)
			{
				SS[i].tc = CUPR_TC;
			}
			P1 = CUPR_P1;
			P2 = CUPR_P2;
			P3 = CUPR_P3;

			

			SS[0].Initialize(1, -1, 1, 0, 1, 1);
			SS[1].Initialize(1, 1, -1, 0, 1, 1);
			SS[2].Initialize(1, 1, 1, 0, -1, 1);
			SS[3].Initialize(-1, 1, 1, 0, -1, 1);

			SS[4].Initialize(1, 1, 1, -1, 1, 0);
			SS[5].Initialize(1, 1, -1, -1, 1, 0);
			SS[6].Initialize(-1, 1, 1, 1, 1, 0);
			SS[7].Initialize(1, -1, 1, 1, 1, 0);

			SS[8].Initialize(1, 1, -1, 1, 0, 1);
			SS[9].Initialize(-1, 1, 1, 1, 0, 1);
			SS[10].Initialize(1, 1, 1, -1, 0, 1);
			SS[11].Initialize(1, -1, 1, -1, 0, 1);

			SS[12].Initialize(1, -1, 1, 0, -1, -1);
			SS[13].Initialize(1, 1, -1, 0, -1, -1);
			SS[14].Initialize(1, 1, 1, 0, 1, -1);
			SS[15].Initialize(-1, 1, 1, 0, 1, -1);

			SS[16].Initialize(1, 1, 1, 1, -1, 0);
			SS[17].Initialize(1, 1, -1, 1, -1, 0);
			SS[18].Initialize(-1, 1, 1, -1, -1, 0);
			SS[19].Initialize(1, -1, 1, -1, -1, 0);

			SS[20].Initialize(1, 1, -1, -1, 0, -1);
			SS[21].Initialize(-1, 1, 1, -1, 0, -1);
			SS[22].Initialize(1, 1, 1, 1, 0, -1);
			SS[23].Initialize(1, -1, 1, 1, 0, -1);

			break;
		}
		case 2://Титан
		{
			SS_count = SS_COUNT_HCP;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < 6; i++)
			{
				SS[i].tc = TITAN_TC1;
			}
			for (int i = 6; i < 12; i++)
			{
				SS[i].tc = TITAN_TC2;
			}
			for (int i = 12; i < 36; i++)
			{
				SS[i].tc = TITAN_TC3;
			}
				
			const double c = TITAN_C;
			//Базисное скольжение
			SS[0].Initialize(0, 0, 0, 1, 1, 1, -2, 0, c);// {0001} <11-20> 
			SS[1].Initialize(0, 0, 0, 1, -1, -1, 2, 0, c);
			SS[2].Initialize(0, 0, 0, 1, -2, 1, 1, 0, c);
			SS[3].Initialize(0, 0, 0, 1, 2, -1, -1, 0, c);
			SS[4].Initialize(0, 0, 0, 1, 1, -2, 1, 0, c);
			SS[5].Initialize(0, 0, 0, 1, -1, 2, -1, 0, c);
			//Призматическое скольжение
			SS[6].Initialize(1, -1, 0, 0, 1, 1, -2, 0, c);// {10-10} <11-20>
			SS[7].Initialize(1, -1, 0, 0, -1, -1, 2, 0, c);
			SS[8].Initialize(0, -1, 1, 0, -2, 1, 1, 0, c);
			SS[9].Initialize(0, -1, 1, 0, 2, -1, -1, 0, c);
			SS[10].Initialize(-1, 0, 1, 0, 1, -2, 1, 0, c);
			SS[11].Initialize(-1, 0, 1, 0, -1, 2, -1, 0, c);
			//Пирамидальное скольжение <a>
			SS[12].Initialize(1, -1, 0, 1, 1, 1, -2, 0, c);// {10-11} <11-20>
			SS[13].Initialize(1, -1, 0, 1, -1, -1, 2, 0, c);
			SS[14].Initialize(0, -1, 1, 1, -2, 1, 1, 0, c);
			SS[15].Initialize(0, -1, 1, 1, 2, -1, -1, 0, c);
			SS[16].Initialize(-1, 0, 1, 1, 1, -2, 1, 0, c);
			SS[17].Initialize(-1, 0, 1, 1, -1, 2, -1, 0, c);
			//Пирамидальное скольжение <a + c>
			SS[18].Initialize(1, -1, 0, 0, 1, 1, -2, 3, c);// {10-10} <11-23>
			SS[19].Initialize(1, -1, 0, 0, -1, -1, 2, 3, c);
			SS[20].Initialize(0, -1, 1, 0, -2, 1, 1, 3, c);
			SS[21].Initialize(0, -1, 1, 0, 2, -1, -1, 3, c);
			SS[22].Initialize(-1, 0, 1, 0, 1, -2, 1, 3, c);
			SS[23].Initialize(-1, 0, 1, 0, -1, 2, -1, 3, c);

			SS[24].Initialize(1, -2, 1, 1, 1, 1, -2, 3, c);// {11-21} <11-23>
			SS[25].Initialize(1, -2, 1, 1, -1, -1, 2, 3, c);
			SS[26].Initialize(1, 1, -2, 1, -2, 1, 1, 3, c);
			SS[27].Initialize(1, 1, -2, 1, 2, -1, -1, 3, c);
			SS[28].Initialize(-2, 1, 1, 1, 1, -2, 1, 3, c);
			SS[29].Initialize(-2, 1, 1, 1, -1, 2, -1, 3, c);

			SS[24].Initialize(1, -2, 1, 2, 1, 1, -2, 3, c);// {11-22} <11-23>
			SS[25].Initialize(1, -2, 1, 2, -1, -1, 2, 3, c);
			SS[26].Initialize(1, 1, -2, 2, -2, 1, 1, 3, c);
			SS[27].Initialize(1, 1, -2, 2, 2, -1, -1, 3, c);
			SS[28].Initialize(-2, 1, 1, 2, 1, -2, 1, 3, c);
			SS[29].Initialize(-2, 1, 1, 2, -1, 2, -1, 3, c);
			break;
		}
		
		}
		if (material == 0 || material == 1)
		{
			//Подобная симметрия допустима только для кубических решёток
			p.C[0][0][0][0] = p.C[1][1][1][1] = p.C[2][2][2][2] = P1;

			p.C[0][0][1][1] = p.C[1][1][0][0] = p.C[2][2][1][1] =
				p.C[1][1][2][2] = p.C[2][2][0][0] = p.C[0][0][2][2] = P2;

			p.C[0][1][1][0] = p.C[1][2][1][2] = p.C[0][1][0][1] =
				p.C[2][0][2][0] = p.C[0][2][0][2] = p.C[2][1][2][1] =
				p.C[1][0][1][0] = p.C[1][0][0][1] = p.C[2][1][1][2] =
				p.C[0][2][2][0] = p.C[1][2][2][1] = p.C[2][0][0][2] = P3;
		}
		else if (material == 2)//ГПУ, альфа-титан
		{
			p.C[2][2][2][2] = TITAN_P1;

			p.C[0][0][0][0] = p.C[1][1][1][1] = TITAN_P2;
			
			p.C[0][0][1][1] = p.C[1][1][0][0] = TITAN_P4;

			p.C[0][0][2][2] = p.C[2][2][0][0] = p.C[1][1][2][2] = p.C[2][2][1][1] = TITAN_P3;

			p.C[1][2][1][2] = p.C[2][1][1][2] = p.C[1][2][2][1] = p.C[2][1][2][1] =
				p.C[2][0][2][0] = p.C[2][0][0][2] = p.C[0][2][0][2] = p.C[0][2][2][0] = TITAN_P5;

			p.C[0][1][1][0] = p.C[1][0][1][0] = p.C[1][0][0][1] = p.C[0][1][0][1] = TITAN_P6;
			
		}
	}

	void Fragment::NDScalc()
	{
		/*
		* Функция вычисляет все характеристики НДС фрагмента
		*/
		for (int k = 0; k < SS_count; k++)
		{
			/************************************************************
			***********              Закон Шмида             ************
			************************************************************/
			SS[k].t = 0;
			if (!SYMMETRY)
			{
				for (int i = 0; i < DIM; i++)
				{
					for (int j = 0; j < DIM; j++)
					{
						SS[k].t += sgm.C[i][j] * SS[k].n.C[i] * SS[k].b.C[j];
					}
				}
			}
			else
			{
				
				for (int i = 0; i < DIM; i++)
				{
					for (int j = 0; j < DIM; j++)
					{
						SS[k].t += sgm.C[i][j] * (SS[k].n.C[i] * SS[k].b.C[j] + SS[k].n.C[j] * SS[k].b.C[i]);
					}
				}
				SS[k].t /= 2.0;
			}

			/************************************************************
			***********      Соотношение Хатчинсона          ************
			************************************************************/
			SS[k].dgm = dgm0 * pow(fabs(SS[k].t / SS[k].tc), m) * H(SS[k].t - SS[k].tc);

			SS[k].gmm += SS[k].dgm * dt; //Сдвиг по СС
		}
		/************************************************************
		**********    Вычисление неупругих деформаций     ***********
		************************************************************/
		d_in.setZero();
		if (!SYMMETRY)
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < SS_count; k++)
					{
						d_in.C[i][j] += SS[k].dgm * SS[k].n.C[i] * SS[k].b.C[j];
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < SS_count; k++)

					{
						d_in.C[i][j] += SS[k].dgm * (SS[k].n.C[i] * SS[k].b.C[j] + SS[k].n.C[j] * SS[k].b.C[i]);
					}
				}
			}
			d_in /= 2.0;
		}
		

		/************************************************************
		**********               Закон Гука               ***********
		************************************************************/
		Tensor buf1 = om*sgm;	//Для коротационных производных
		Tensor buf2 = sgm*om;

		dsgm.setZero();
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int l = 0; l < DIM; l++)
				{
					for (int n = 0; n < DIM; n++)
					{
						dsgm.C[i][j] += p.C[i][j][l][n] * (d.C[n][l] - d_in.C[n][l]);
					}
				}
			}
		}
		dsgm += buf1;
		dsgm -= buf2;

		/************************************************************
		**********  Интенсивности деформаций и напряжений  **********
		************************************************************/

		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				e.C[i][j] += d.C[i][j] * dt;
				sgm.C[i][j] += dsgm.C[i][j] * dt;
			}
		}
		//strain = 0;
		//stress = 0;
		strain = e.doubleScalMult(e);
		stress = sgm.doubleScalMult(sgm);
		/*for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				strain += e.C[i][j] * e.C[j][i];
				stress += sgm.C[i][j] * sgm.C[j][i];
			}
		}*/
		strain = SQRT2_3 * sqrt(strain);
		stress = SQRT3_2 * sqrt(stress);
	}

	void Fragment::Rotate(double dFi, const Vector a)
	{
		/*
		* Поворот решётки фрагмента вокруг
		* Заданной оси на заданный угол
		*/
		double CosFi = cos(dFi);
		double SinFiX = sin(dFi)*a.C[0];
		double SinFiY = sin(dFi)*a.C[1];
		double SinFiZ = sin(dFi)*a.C[2];
		double COS = (1.0 - CosFi);
		
		Tensor dO;			//Тензор поворота на шаге

		dO.C[0][0] = CosFi + COS * a.C[0] * a.C[0];
		dO.C[0][1] = COS * a.C[0]*a.C[1] - SinFiZ;
		dO.C[0][2] = COS * a.C[0] * a.C[2] + SinFiY;

		dO.C[1][0] = COS * a.C[0] * a.C[1] + SinFiZ;
		dO.C[1][1] = CosFi + COS * a.C[1] * a.C[1];
		dO.C[1][2] = COS * a.C[1] * a.C[2] - SinFiX;

		dO.C[2][0] = COS * a.C[0] * a.C[2] - SinFiY;
		dO.C[2][1] = COS * a.C[1] * a.C[2] + SinFiX;
		dO.C[2][2] = CosFi + COS * a.C[2] * a.C[2];

		Tensor buf = dO*o;
		o = buf;

	}

}