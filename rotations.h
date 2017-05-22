#ifndef __ROTATIONS_H 
#define __ROTATIONS_H

#include "fragment.h"

/*
*Ротационные механизмы
*/

namespace model
{
	
	void Taylor_rotations(Fragment*);		//Модель стеснённого поворота по Тейлору
	void Trusov_rotations(Fragment*);		//Модель, связанная с несовместностью сдвигов
	void Rotation_hardening(Fragment*);		//Модель ротационного упрочнения
	void Rotate(Fragment* f, double dFi, const Vector a);	//Поворот вокруг заданной оси на заданный угол
	/*********************************************************
	********	    Получение полюсных фигур и ССТ	   *******
	*********************************************************/

	void SavePoints(Tensor,	char*,
		int, int, int);			//Запись проекций ПФ в файл
	void GetPoleFig(Fragment*);				//Сохранение полюсной фигуры
	
	void SaveSSTPoints(Tensor&, float,		//Запись проекций ССТ в файл
		char*, int, int, int);
	void GetSST(Fragment*);					//Сохранение ССТ
}

#endif __ROTATIONS_H