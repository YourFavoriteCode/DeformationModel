#ifndef __POLYCRYST_H 
#define __POLYCRYST_H

#include "fragment.h"

namespace model
{
	class Polycrystall
	{
	public:
		Fragment ***C;				//������ ��������� �������������
		int fragm_count;			//���-�� ���������� �� �����
		int total_fragm_count;		//����� ���-�� ����������

		Tensor D;					//������ ���������� �������� �����������
		Tensor W;					//������ ����� �����������
		Tensor D_in;				//������ ��������� ����� ���������� �����������
		Tensor E;					//������ ���������������
		Tensor dSgm;				//������ ��������� ���������������
		Tensor Sgm;					//������ ���������������
		Tensor4 P;					//����������� ������ ������� ��������

		int cycle;
		int CURR_STEP;				//������� ��� ��������������
		int PLOT_STEP;					//��� ���������� ��������
		int POLUS_STEP;					//��� ���������� ��
		int PROC_STEP;				//��� ����������� ���������
		int DEBUG_STEP;				//��� ������ ���������� ������
		int proc_period;			//������ ���������� �������� ����������

		double Strain;				//������������� ���������������
		double Stress;				//������������� ���������������

		int file_count;
		std::ofstream *dbgstream;
		std::ofstream *TestStream;
		std::ofstream *Datastream;
	
		Polycrystall();
		~Polycrystall();

		double tension_component;	//������������ ���������� � ��������
		double final_stress;		//��������, �� �������� ����������
		double lam;				//����������� � ���������
		double addition_strain;	//���������� ��������� ��� ����������� �������

		void Init(int);				//��������� ������ ��� ����
		void MakeStruct();			//������������� �������� � �������
		void setParams();
		void Deformate();
		
		void Fragmentate();
		void SavePoleFig();
		void SaveDbgInfo();
		void OpenFiles();
		void CloseFiles();

		int getNumArr(int, int, int);
		void getMultyNumArr(int, int&, int&, int&);
	private:
		void Load(bool unload);				//����������
	};
}

#endif __POLYCRYST_H