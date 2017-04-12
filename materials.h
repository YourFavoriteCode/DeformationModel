#ifndef __MATERIALS_H 
#define __MATERIALS_H
/*
*В этом файле определены все материальные параметры-константы
*/
namespace model
{
	const int SS_COUNT_BCC = 96;		//Кол-во СС в ОЦК
	const int SS_COUNT_FCC = 24;		//Кол-во СС в ГЦК
	const int SS_COUNT_HCP = 36;		//Кол-во СС в ГПУ

	/********************************************
	***********			ОЦК		    *************
	********************************************/

	//Сталь 45
	/******       Упругие константы        *****/
	const double STEEL_P1 = 2.2e11;
	const double STEEL_P2 = 1.66e11;
	const double STEEL_P3 = 8.7e10;

	/***** Начальные критические напряжения *****/
	const double STEEL_TC1 = 0.1e9;
	const double STEEL_TC2 = 0.47e9;
	const double STEEL_TC3 = 2.47e9;

	/********************************************
	***********			ГЦК		    *************
	********************************************/

	//Чистая медь
	/******       Упругие константы        *****/
	const double CUPR_P1 = 1.684e11;
	const double CUPR_P2 = 1.214e11;
	const double CUPR_P3 = 7.54e10;

	/***** Начальные критические напряжения *****/
	const double CUPR_TC = 1.75e7;

	//Алюминий
	const double ALLUM_P1 = 1.08e11;
	const double ALLUM_P2 = 6.1e10;
	const double ALLUM_P3 = 2.8e10;

	/********************************************
	***********			ГПУ		    *************
	********************************************/
	//Альфа-титан
	const double TITAN_TC1 = 15e7;
	const double TITAN_TC2 = 3e7;
	const double TITAN_TC3 = 12e7;

	const double TITAN_P1 = 1.807e11;
	const double TITAN_P2 = 1.624e11;
	const double TITAN_P3 = 6.9e10;
	const double TITAN_P4 = 9.2e10;
	const double TITAN_P5 = 4.67e10;
	const double TITAN_P6 = 3.5e10;

	const double TITAN_C = 1.587;

}

#endif __MATERIALS_H