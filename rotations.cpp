// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <fstream>

#include "rotations.h"
#include "params.h"
#include "fragment.h"
#include "distributions.h"


namespace model
{
	void Rotate(Fragment* f, double dFi, const Vector a)
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
		dO.C[0][1] = COS * a.C[0] * a.C[1] - SinFiZ;
		dO.C[0][2] = COS * a.C[0] * a.C[2] + SinFiY;

		dO.C[1][0] = COS * a.C[0] * a.C[1] + SinFiZ;
		dO.C[1][1] = CosFi + COS * a.C[1] * a.C[1];
		dO.C[1][2] = COS * a.C[1] * a.C[2] - SinFiX;

		dO.C[2][0] = COS * a.C[0] * a.C[2] - SinFiY;
		dO.C[2][1] = COS * a.C[1] * a.C[2] + SinFiX;
		dO.C[2][2] = CosFi + COS * a.C[2] * a.C[2];

		Tensor buf = dO*f->o;
		f->o = buf;

	}
	
	void Taylor_rotations(Fragment *f)
	{
		f->om.setZero();
		//f->om = f->w - f->d_in.getAntiSymmetryPart();
		for (int i = 0; i < DIM; i++)//Спин решётки
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < f->SS_count; k++)
				{
					f->om.C[i][j] -= f->SS[k].dgm * (f->SS[k].n.C[i] * f->SS[k].b.C[j] - f->SS[k].b.C[i] * f->SS[k].n.C[j]);
				}
			}
		}
		f->om /= 2.0;
		f->om += f->w;
		
		Vector e;				//Ось вращения решётки
		e.set(f->om.C[1][2], f->om.C[2][0], f->om.C[0][1]);
		
		
		double dFi = e.getNorm();
		f->isRotate = dFi > EPS;
		if (f->isRotate)
		{
			f->rot_speed = dFi;
			dFi *= prms::dt;			//Получаем угол поворота
			f->sum_angle += dFi;
			e.Normalize();
			Rotate(f, dFi, e);
	
		}
		else
		{
			f->rot_speed = 0;
		
		}
	}

	void Rotation_hardening(Fragment *f)
	{
		if (f->sum_angle > 100*EPS)
		{
			double vol = pow(f->size, 3);//Объём элемента
			double dmc = prms::ROT_HARD_K1 / vol * exp( - prms::ROT_HARD_K2 * f->sum_angle);//Скорость приращения
			f->dmc = dmc;
			f->rot_Mc += dmc*prms::dt;//Приращение критического момента
		}
	}

	void Trusov_rotations(Fragment *f)
	{
		Vector dM;						//Производная вектор-момента
		double S = f->size*f->size;		//Площадь фасетки (полная)
		for (int h = 0; h < prms::surround_count; h++)//Пробегаем по всем соседям фрагмента
		{
			if (f->contact[h] == 0) continue;//Если нет контакта - пропускаем
			Tensor d_in1, d_in2;
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < f->SS_count; k++)
					{
						if (f->SS[k].b.ScalMult(f->normals[h]) < 0) continue; //Скольжение от границы - вклад не вносится
						d_in1.C[i][j] += f->SS[k].dgm * (f->SS[k].n.C[i] * f->SS[k].b.C[j]);
					}
					for (int k = 0; k < f->surrounds[h].SS_count; k++)
					{
						//if (f->SS[k].b.ScalMult(f->normals[h]) < 0) continue; //Скольжение от границы - вклад не вносится
						d_in2.C[i][j] += f->surrounds[h].SS[k].dgm * (f->surrounds[h].SS[k].n.C[i] * f->surrounds[h].SS[k].b.C[j]);
					}
				}
			}
			Tensor Lp = d_in1 - d_in2;
			Lp.Transp();
			Tensor buf = VectMult(f->normals[h], Lp);
			Vector dm = ScalMult(buf, f->normals[h]);//Поверхностный вектор-момент
			dm *= f->rot_L;
			Vector b1 = ScalMult(f->om, dm);//(коротационная производная)
			Vector b2 = ScalMult(dm, f->om);
			dm = dm + b1 - b2;
			double c;		//Определяет площадь контакта (в долях от полной площади стороны куба)
			if (h < 6) c = UniformDistrib(0.93, 0.07);
			else if (h < 14) c =  UniformDistrib(0.1, 0.05);
			else c = UniformDistrib(0.01, 0.005);
			dM += dm*S*c;
		}
		dM /= f->volume;
		double dMnorm = dM.getNorm();
		Vector M = f->moment + dM*prms::dt;
		f->moment = M;
		double norm = M.getNorm();
		
		if (norm > f->rot_Mc || norm == -1)
		{
			norm = f->rot_Mc;
		}
		f->norm = norm;
	
		double pr = M.ScalMult(dM);
	
		double dFi = 0;		//Скорость вращения
		if (norm == f->rot_Mc && pr >= 0)	//Пластические и упругие развороты
		{
			dFi = f->rot_A * dMnorm + f->rot_H * norm;
		}
		else
		{
			dFi = f->rot_A * dMnorm;		//Только упругие развороты
		}
		
		f->isRotate = (dFi > EPS*1e4);
		if (f->isRotate)
		{
			Vector e = M;					//Ось вращения решётки сонаправлена с вектором момента
			e.Normalize();
			f->rot_speed = dFi;
			dFi *= prms::dt;
			f->sum_angle += dFi;		//Накопленный угол вращения увеличивается
			Rotate(f, dFi, e);			//Вращение решётки
			f->rot_energy = norm*dFi;	//Считаем энергию ротаций

			f->om.setZero();			//Спин решётки
			/*for (int i = 0; i < DIM; i++)//Спин решётки (Часть от Тейлора)
			{
				for (int j = 0; j < DIM; j++)
				{
					
					for (int k = 0; k < f->SS_count; k++)
					{
						f->om.C[i][j] -= f->SS[k].dgm * (f->SS[k].n.C[i] * f->SS[k].b.C[j] - f->SS[k].b.C[i] * f->SS[k].n.C[j]);
					}
					f->om.C[i][j] /= 2.0;
				
				}
			}
			f->om += f->w;
			*/
			for (int i = 0; i < DIM; i++)    //Спин решётки (Вторая часть)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < DIM; k++)
					{
						f->om.C[i][j] -= LeviCivit(i, j, k) * e.C[k] * f->rot_speed;
					}
				}
			}
			
		}
		else
		{
			f->om.setZero();
			f->rot_speed = 0;		//Решётка не вращается
			f->rot_energy = 0;		//Энергия вращения равна нулю
		}
	}

	void inline SavePoints(Tensor O, const char *file, const int i, const int j, const int k)
	{
		Vector e;
		e.set(i, j, k);
		e.Normalize();
		Vector e1 = ScalMult(e, O);//Перевод в ЛСК
		e1.Normalize();
		std::ofstream of;
		of.open(file, std::ios::out | std::ios_base::app | std::ios::binary);
		//Запись в виде x,y,z
		of.write((char *)&e1.C[0], sizeof(double));
		of.write((char *)&e1.C[1], sizeof(double));
		of.write((char *)&e1.C[2], sizeof(double));
		of.close();
	}

	void GetPoleFig(Fragment *f)
	{
		/*---Семейство направлений [001]---*/
		SavePoints(f->o, "Polus\\S001.dat", 0, 0, 1);
		SavePoints(f->o, "Polus\\S010.dat", 0, 1, 0);
		SavePoints(f->o, "Polus\\S100.dat", 1, 0, 0);

		/*---Семейство направлений [011]---*/
		SavePoints(f->o, "Polus\\S011.dat", 0, 1, 1);
		SavePoints(f->o, "Polus\\S110.dat", 1, 1, 0);
		SavePoints(f->o, "Polus\\S101.dat", 1, 0, 1);

		/*---Семейство направлений [111]---*/
		SavePoints(f->o, "Polus\\S1-11.dat", 1, -1, 1);
		SavePoints(f->o, "Polus\\S-111.dat", -1, 1, 1);
		SavePoints(f->o, "Polus\\S111.dat", 1, 1, 1);
	}

	void inline GetSSTPoint(Tensor O, double dFi, const char *FileName, const int i, const int j, const int k)
	{
		Vector a;
		a.set(i, j, k);
		a.Normalize();

		double CosFi = cos(dFi);
		double SinFiX = sin(dFi) * a.C[0];
		double SinFiY = sin(dFi) * a.C[1];
		double SinFiZ = sin(dFi) * a.C[2];
		double COS = (1.0 - CosFi);

		Tensor dO;			//Тензор поворота на шаге

		dO.C[0][0] = CosFi + COS * a.C[0] * a.C[0];
		dO.C[0][1] = COS * a.C[0] * a.C[1] - SinFiZ;
		dO.C[0][2] = COS * a.C[0] * a.C[2] + SinFiY;

		dO.C[1][0] = COS * a.C[0] * a.C[1] + SinFiZ;
		dO.C[1][1] = CosFi + COS * a.C[1] * a.C[1];
		dO.C[1][2] = COS * a.C[1] * a.C[2] - SinFiX;

		dO.C[2][0] = COS * a.C[0] * a.C[2] - SinFiY;
		dO.C[2][1] = COS * a.C[1] * a.C[2] + SinFiX;
		dO.C[2][2] = CosFi + COS * a.C[2] * a.C[2];
		
		Tensor R = dO * O;
		Vector v = ScalMult(a, R);
		v.Normalize();
	
		std::ofstream Of;
		Of.open(FileName, std::ios::out | std::ios_base::app | std::ios::binary);
		//Запись координат в формате [xyz] подряд для всех зёрен 
		Of.write((char *)&v.C[0], sizeof(double));
		Of.write((char *)&v.C[1], sizeof(double));
		Of.write((char *)&v.C[2], sizeof(double));
		Of.close();
	}

	void GetSST(Fragment *f)
	{
		for (double fi = 0; fi < 2 * PI; fi += PI_2)//Ось четвёртого порядка
		{
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 1, 0, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", -1, 0, 0);
		}
		for (double fi = 0; fi < 2 * PI; fi += PI)//Ось второго порядка
		{
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, -1, -1);
		}
		for (double fi = 0; fi < 2 * PI; fi += 2.0 * PI / 3.0)//Ось третьего порядка
		{
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, -1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, -1, -1);
		}

	}
}