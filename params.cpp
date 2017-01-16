#include "stdafx.h"
#include <iostream>

#include "params.h"
#include "TinyXML\tinyxml2.h"

namespace model
{
	/*****************************************************************
	*****		Задание значений параметров по умолчанию		******
	*****************************************************************/
	bool SYMMETRY = true;
	bool REAL_UNIAX = false;
	bool UNLOADING = false;
	bool RAND_ORIENT = true;
	int fix_orient = 0;
	double dt = 5e-4;
	int material = 0;
	double strain_max = 1e-1;
	Tensor gradV;
	
	int surround_count = 6;
	int material_purity = 100;
	bool read_init_stress = false;
	
	int fragm_size_law = 0;
	double fragm_size_m = 5e-5;
	double fragm_size_dsp = 0;

	int fragm_count = 4;
	int cycle_count = 1;
	int thread_count = 1;

	double plot_period = 2;
	double polus_period = 25;
	int debug_period = 0;
	int DEBUG_START = 0;
	int DEBUG_STOP = INT_MAX;

	double dgm0 = 1e-5;
	double m = 100;

	bool ROTATIONS_TAYLOR = false;
	bool ROTATIONS_TRUSOV = false;
	bool ROTATIONS_HARDENING = false;

	double ROT_A = 3e-9;
	double ROT_H = 1e-7;
	double ROT_L = 10;
	double ROT_MC = 3e4;

	bool HARDENING_BASE = false;
	bool HARDENING_BOUND = false;

	double HARD_BOUND_K = 0;
	double HARD_BASE_DELTA = 0;
	double HARD_BASE_PSI = 0;
	double HARD_BASE_A = 0;

	int SurroundsGrade = 1;
	bool SST_SAVING = false;

	bool FRAGMENTATION = false;

	int getValue(tinyxml2::XMLElement *rootnode, const char* name, int* var)//Целочисленные
	{
		if (rootnode->FirstChildElement(name) != NULL)//Проверка, есть ли в файле такой элемент
		{
			const char* title;
			title = rootnode->FirstChildElement(name)->GetText();
			*var = atoi(title);
		}
		else//Если нет, то предупреждаем и выводим значение по умолчанию
		{
			printf(" Parameter <%s> not found. Default value: %d.\n", name, *var);
		}
		return 1;
	}
	
	int getValue(tinyxml2::XMLElement *rootnode, const char* name, double* var)//Вещественные
	{
		if (rootnode->FirstChildElement(name) != NULL)//Проверка, есть ли в файле такой элемент
		{
			const char* title;
			title = rootnode->FirstChildElement(name)->GetText();
			*var = atof(title);
		}
		else//Если нет, то предупреждаем и выводим значение по умолчанию
		{
			printf(" Parameter <%s> not found. Default value: %e.\n", name, *var);
		}
		return 1;
	}

	int getValue(tinyxml2::XMLElement *rootnode, const char* name, bool* var)//Логические
	{
		if (rootnode->FirstChildElement(name) != NULL)//Проверка, есть ли в файле такой элемент
		{
			const char* title;
			title = rootnode->FirstChildElement(name)->GetText();
			*var = (bool)atoi(title);
		}
		else//Если нет, то предупреждаем и выводим значение по умолчанию
		{
			printf(" Parameter <%s> not found. Default value: %d.\n", name,*var);
		}
		return 1;
	}

	int ReadParams(const char * filename)
	{
		tinyxml2::XMLDocument doc;
		doc.LoadFile(filename);//Загрузка файла

		tinyxml2::XMLElement *rootnode = doc.FirstChildElement("Parameters");
		if (rootnode == NULL) return 1;
		getValue(rootnode, "GrainCount", &fragm_count);
		getValue(rootnode, "Material", &material);
		getValue(rootnode, "RandomOrientations", &RAND_ORIENT);
		getValue(rootnode, "M", &m);
		getValue(rootnode, "dGamma0", &dgm0);
		getValue(rootnode, "DeformationLimit", &strain_max);
		getValue(rootnode, "IntegrationStep", &dt);
		getValue(rootnode, "RealUniaxial", &REAL_UNIAX);
		getValue(rootnode, "Symmetrisation", &SYMMETRY);
		/*Упрочнение*/
		getValue(rootnode, "HardBase", &HARDENING_BASE);
		getValue(rootnode, "HardBound", &HARDENING_BOUND);
		getValue(rootnode, "HardBaseDelta", &HARD_BASE_DELTA);
		getValue(rootnode, "HardBasePsi", &HARD_BASE_PSI);
		getValue(rootnode, "HardBaseA", &HARD_BASE_A);
		getValue(rootnode, "HardBoundK", &HARD_BOUND_K);
		/*Ротации*/
		int rt;
		getValue(rootnode, "Rotate", &rt);
		if (rt == 0) ROTATIONS_TAYLOR = true;
		else if (rt == 1) ROTATIONS_TRUSOV = true;
		getValue(rootnode, "RotateA", &ROT_A);
		getValue(rootnode, "RotateH", &ROT_H);
		getValue(rootnode, "RotateLambda", &ROT_L);
		getValue(rootnode, "RotateMc", &ROT_MC);

		getValue(rootnode, "CycleCount", &cycle_count);
		getValue(rootnode, "Unloading", &UNLOADING);
		getValue(rootnode, "SurroundsDegree", &SurroundsGrade);
		getValue(rootnode, "RotationHardening", &ROTATIONS_HARDENING);
		getValue(rootnode, "MaterialPurity", &material_purity);
		getValue(rootnode, "PlotPeriod", &plot_period);
		getValue(rootnode, "PolusPeriod", &polus_period);
		getValue(rootnode, "DebugPeriod", &debug_period);
		getValue(rootnode, "ThreadCount", &thread_count);
		getValue(rootnode, "FixedOrientations", &fix_orient);
		
		/*Градиент скорости*/
		gradV.setZero();
		getValue(rootnode, "gradV00", &gradV.C[0][0]);
		getValue(rootnode, "gradV01", &gradV.C[0][1]);
		getValue(rootnode, "gradV02", &gradV.C[0][2]);
		getValue(rootnode, "gradV10", &gradV.C[1][0]);
		getValue(rootnode, "gradV11", &gradV.C[1][1]);
		getValue(rootnode, "gradV12", &gradV.C[1][2]);
		getValue(rootnode, "gradV20", &gradV.C[2][0]);
		getValue(rootnode, "gradV21", &gradV.C[2][1]);
		getValue(rootnode, "gradV22", &gradV.C[2][2]);
	
		getValue(rootnode, "FragmSizeLaw", &fragm_size_law);
		getValue(rootnode, "FragmSizeM", &fragm_size_m);
		getValue(rootnode, "FragmSizeDsp", &fragm_size_dsp);
		getValue(rootnode, "StartWritingDbgInfo", &DEBUG_START);
		getValue(rootnode, "StopWritingDbgInfo", &DEBUG_STOP);
		
		/*Экспериментальные параметры*/
		getValue(rootnode, "ReadInitStress", &read_init_stress);
		getValue(rootnode, "SaveSST", &SST_SAVING);
		getValue(rootnode, "Fragmentation", &FRAGMENTATION);

		return 0;
	}
}