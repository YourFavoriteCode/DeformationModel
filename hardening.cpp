﻿// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <cmath>

#include "params.h"
#include "materials.h"
#include "hardening.h"
namespace model
{
	void Base_hardening(Fragment *f)
	{
		double nsd = 0;	//Накопленный сдвиг
		for (int k = 0; k < f->SS_count; k++)
		{
			nsd += f->SS[k].gmm;		//Суммируем по всем СС
		}
		for (int k = 0; k < f->SS_count; k++)
		{
			double osn = 0;
			for (int j = 0; j < f->SS_count; j++)
			{
				double a = (j == k ? prms::HARD_BASE_A : 1.25 * prms::HARD_BASE_A);
				if ((f->SS[j].dgm > EPS) && (fabs(nsd) > EPS))
				{
					osn += a*pow(f->SS[j].gmm, prms::HARD_BASE_PSI)*pow((f->SS[j].dgm / prms::dgm0), prms::HARD_BASE_DELTA) / nsd;
				}
			}
			double napr = 0;
			switch (f->material)
			{
			case 0://Сталь 45 (ОЦК)
			{
				if (k < 24) napr = STEEL_TC1;
				else if (k < 48) napr = STEEL_TC2;
				else napr = STEEL_TC3;
				break;
			}
			case 1://Медь (ГЦК)
			{
				napr = CUPR_TC;
				break;
			}
			case 2://Титан (ГПУ)
			{
				if (k < 6) napr = TITAN_TC1;
				else if (k < 12) napr = TITAN_TC2;
				else napr = TITAN_TC3;
				break;
			}

			}
			f->SS[k].tc += napr*osn*prms::dt;
		}
	}
	
	void Boundary_hardening(Fragment *f)
	{
		for (int k = 0; k < f->SS_count; k++)	//Цикл по СС текущего фрагмента
		{
			Vector b1 = ScalMult(f->o, f->SS[k].b);//Перевели вектор b текущей СС данного зерна в ЛСК
			double zgu = 0;
			for (int h = 0; h < prms::surround_count; h++)	//Цикл по фасеткам			
			{
				if (f->contact[h] == 0) continue;//Если нет контакта - пропускаем
				if (f->SS[k].b.ScalMult(f->normals[h]) < 0) continue; //Скольжение от границы - пропускаем
				double zguk = prms::HARD_BOUND_K * f->SS[k].dgm * f->SS[k].gmm / f->size;
				double min = 1.0;//Минимум
				min = f->DisorientMeasure(h);
				/*for (int p = 0; p < f->surrounds[h].SS_count; p++)	//Цикл по системам соседнего зерна
				{
					Vector b2 = ScalMult(f->surrounds[h].o, f->surrounds[h].SS[p].b);//Перевели вектор b p-ой СС соседнего зерна в ЛСК
					Vector diff = b1 - b2;
					diff.Normalize();
					double M = fabs(diff.ScalMult(f->normals[h]));

					if (M < min)
					{
						min = M;
					}
				}*/
				zgu += zguk*min;
			}

			f->SS[k].tc += zgu;
		}
	}

}