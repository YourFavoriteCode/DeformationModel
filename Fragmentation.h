#ifndef __FRAGMENTATION_H 
#define __FRAGMENTATION_H

#include "MathCore.h"
#include <vector>

/*
Работа с фрагментацией
*/
namespace model
{
	int get1DPos(int, int, int);
	void get3DPos(int, int*, int*, int*);

	class Grain			//Зерно (совокупность фрагментов)
	{
	public:
		int num;		//Порядковый номер
		int center;		//Центр кристаллизации
		int size;		//Характерный размер
	};


}

#endif