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
	double dt = 1e-4;
	int material = 0;
	double strain_max = 5e-2;
	Tensor gradV;
	int surround_count = 6;
	int material_purity = 100;
	bool read_init_stress;
	
	int fragm_size_law = 0;
	double fragm_size_m = 5e-5;
	double fragm_size_dsp = 0;

	int fragm_count = 3;
	int cycle_count = 1;
	int thread_count = 1;

	double plot_period = 2;
	double polus_period = 20;
	int debug_period = 1000;
	int DEBUG_START = 0;
	int DEBUG_STOP = INT_MAX;

	double dgm0 = 1e-5;
	double m = 100;

	bool ROTATIONS_TAYLOR = false;
	bool ROTATIONS_TRUSOV = false;
	bool ROTATIONS_HARDENING = false;

	double ROT_A = 1e-8;
	double ROT_H = 2e-8;
	double ROT_L = 3e6;
	double ROT_MC = 2e5;

	bool HARDENING_BASE = false;
	bool HARDENING_BOUND = false;

	double HARD_BOUND_K = 0;
	double HARD_BASE_DELTA = 0;
	double HARD_BASE_PSI = 0;
	double HARD_BASE_A = 0;

	int SurroundsGrade = 1;
	bool SST_SAVING = false;

	bool FRAGMENTATION = false;

	int ReadParams(const char * filename)
	{
		tinyxml2::XMLDocument doc;
		doc.LoadFile(filename);//Загрузка файла

		tinyxml2::XMLElement *rootnode = doc.FirstChildElement("Parameters");
		const char* title;
		
		title = rootnode->FirstChildElement("GrainCount")->GetText();
		fragm_count = atoi(title);

		title = rootnode->FirstChildElement("Material")->GetText();
		material = atoi(title);

		title = rootnode->FirstChildElement("Rotate")->GetText();
		if (atoi(title) == 0) ROTATIONS_TAYLOR = true;
		else if (atoi(title) == 1) ROTATIONS_TRUSOV = true;

		title = rootnode->FirstChildElement("RandomOrientations")->GetText();
		RAND_ORIENT = (bool)atoi(title);

		title = rootnode->FirstChildElement("M")->GetText();
		m = atoi(title);

		title = rootnode->FirstChildElement("dGamma0")->GetText();
		dgm0 = atof(title);
		
		title = rootnode->FirstChildElement("DeformationLimit")->GetText();
		strain_max = atof(title);

		title = rootnode->FirstChildElement("IntegrationStep")->GetText();
		dt = atof(title);
		
		title = rootnode->FirstChildElement("RealUniaxial")->GetText();
		REAL_UNIAX = (bool) atoi(title);

		title = rootnode->FirstChildElement("Symmetrisation")->GetText();
		SYMMETRY = (bool) atoi(title);

		title = rootnode->FirstChildElement("HardBase")->GetText();
		HARDENING_BASE = (bool) atoi(title);

		title = rootnode->FirstChildElement("HardBound")->GetText();
		HARDENING_BOUND = (bool) atoi(title);

		title = rootnode->FirstChildElement("HardBaseDelta")->GetText();
		HARD_BASE_DELTA = atof(title);

		title = rootnode->FirstChildElement("HardBasePsi")->GetText();
		HARD_BASE_PSI = atof(title);

		title = rootnode->FirstChildElement("HardBaseA")->GetText();
		HARD_BASE_A = atof(title);

		title = rootnode->FirstChildElement("HardBoundK")->GetText();
		HARD_BOUND_K = atof(title);

		title = rootnode->FirstChildElement("RotateA")->GetText();
		ROT_A = atof(title);

		title = rootnode->FirstChildElement("RotateH")->GetText();
		ROT_H = atof(title);

		title = rootnode->FirstChildElement("RotateLambda")->GetText();
		ROT_L = atof(title);

		title = rootnode->FirstChildElement("RotateMc")->GetText();
		ROT_MC = atof(title);

		title = rootnode->FirstChildElement("CycleCount")->GetText();
		cycle_count = atoi(title);

		title = rootnode->FirstChildElement("Unloading")->GetText();
		UNLOADING = (bool) atoi(title);

		title = rootnode->FirstChildElement("SurroundsDegree")->GetText();
		SurroundsGrade = atoi(title);

		title = rootnode->FirstChildElement("RotationHardening")->GetText();
		ROTATIONS_HARDENING = (bool) atoi(title);

		title = rootnode->FirstChildElement("MaterialPurity")->GetText();
		material_purity = atoi(title);

		title = rootnode->FirstChildElement("PlotPeriod")->GetText();
		plot_period = atof(title);

		title = rootnode->FirstChildElement("PolusPeriod")->GetText();
		polus_period = atof(title);

		title = rootnode->FirstChildElement("DebugPeriod")->GetText();
		debug_period = atoi(title);

		title = rootnode->FirstChildElement("ThreadCount")->GetText();
		thread_count = atoi(title);

		title = rootnode->FirstChildElement("FixedOrientations")->GetText();
		fix_orient = atoi(title);
		
		gradV.setZero();
		
		title = rootnode->FirstChildElement("gradV00")->GetText();
		gradV.C[0][0] = atof(title);

		title = rootnode->FirstChildElement("gradV01")->GetText();
		gradV.C[0][1] = atof(title);

		title = rootnode->FirstChildElement("gradV02")->GetText();
		gradV.C[0][2] = atof(title);

		title = rootnode->FirstChildElement("gradV10")->GetText();
		gradV.C[1][0] = atof(title);

		title = rootnode->FirstChildElement("gradV11")->GetText();
		gradV.C[1][1] = atof(title);

		title = rootnode->FirstChildElement("gradV12")->GetText();
		gradV.C[1][2] = atof(title);
		
		title = rootnode->FirstChildElement("gradV20")->GetText();
		gradV.C[2][0] = atof(title);

		title = rootnode->FirstChildElement("gradV21")->GetText();
		gradV.C[2][1] = atof(title);

		title = rootnode->FirstChildElement("gradV22")->GetText();
		gradV.C[2][2] = atof(title);
		
		title = rootnode->FirstChildElement("FragmSizeLaw")->GetText();
		fragm_size_law = atoi(title);
		title = rootnode->FirstChildElement("FragmSizeM")->GetText();
		fragm_size_m = atof(title);
		title = rootnode->FirstChildElement("FragmSizeDsp")->GetText();
		fragm_size_dsp = atof(title);

		title = rootnode->FirstChildElement("StartWritingDbgInfo")->GetText();
		DEBUG_START = atoi(title);
		title = rootnode->FirstChildElement("StopWritingDbgInfo")->GetText();
		DEBUG_STOP = atoi(title);
		
		title = rootnode->FirstChildElement("ReadInitStress")->GetText();//TODO: Проверка на NULL
		read_init_stress = atoi(title);

		title = rootnode->FirstChildElement("SaveSST")->GetText();
		SST_SAVING = atoi(title);
		
		title = rootnode->FirstChildElement("Fragmentation")->GetText();
		FRAGMENTATION = atoi(title);

		return 0;
	}
}