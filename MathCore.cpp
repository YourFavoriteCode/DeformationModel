// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <cmath>

#include "MathCore.h"
#include "Eigen/Eigen"

namespace model
{
	/***********************************************************************
	*****************      Псевдо-тензор Леви-Чивита     ******************
	***********************************************************************/
	int LeviCivit(const int i, const int j, const int k)
	{
		if (i == 0 && j == 1 && k == 2 || i == 2 && j == 0 && k == 1 || i == 1 && j == 2 && k == 0) return 1;
		else if (i == 2 && j == 1 && k == 0 || i == 0 && j == 2 && k == 1 || i == 1 && j == 0 && k == 2) return -1;
		else return 0;
	}

	/***********************************************************************
	****************      Методы для работы с вектором     *****************
	***********************************************************************/
	Vector::Vector()
	{
		C[0] = 0;
		C[1] = 0;
		C[2] = 0;
	}

	void Vector::set(const double x, const double y, const double z)
	{
		C[0] = x;
		C[1] = y;
		C[2] = z;
	}

	double Vector::getNorm()
	{
		double sum = 0;
		for (int i = 0; i < DIM; i++)
		{
			sum += C[i] * C[i];
		}
		double res;	
		if (isnormal(sum)) res = sqrt(sum); 
		else res = -1;//В случае ошибки вернёт -1
		return res;
	}

	int Vector::Normalize()
	{
		double norm = getNorm();
		if (norm > 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				C[i] /= norm;
			}
			return 0;
		}
		else
		{
			return -1;		//Код ошибки деления на 0 и QNaN
		}
		
	}

	void Vector::setZero()
	{
		for (int i = 0; i < DIM; i++)
		{
			C[i] = 0;
		}
	}

	double Vector::ScalMult(const Vector v)
	{
		double res = 0;
		for (int i = 0; i < DIM; i++)
		{
			res += C[i] * v.C[i];
		}
		return res;
	}
	
	Vector Vector::VectMult(const Vector v)
	{
		Vector res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					res.C[i] += LeviCivit(i, j, k) * C[j] * v.C[k];
				}
			}
		}
		return res;
	}

	Vector Vector::operator+(const Vector v)
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.C[i] = C[i] + v.C[i];
		}
		return buf;
	}

	Vector Vector::operator-(const Vector v)
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.C[i] = C[i] - v.C[i];
		}
		return buf;
	}
	
	Vector Vector::operator-()
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.C[i] = -C[i];
		}
		return buf;
	}

	Vector Vector::operator*(const double r)
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.C[i] = C[i] * r;
		}
		return buf;
	}

	void Vector::operator += (const Vector v)
	{
		for (int i = 0; i < DIM; i++)
		{
			C[i] += v.C[i];
		}
	}
	void Vector::operator -= (const Vector v)
	{
		for (int i = 0; i < DIM; i++)
		{
			C[i] -= v.C[i];
		}
	}
	void Vector::operator *= (const double r)
	{
		for (int i = 0; i < DIM; i++)
		{
			C[i] *= r;
		}
	}

	int Vector::operator /= (const double r)
	{
		if (r != 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				C[i] /= r;
			}
		}
		else
		{
			return -1;
		}
		return 0;
	}

	Vector::~Vector()
	{
	}

	/***********************************************************************
	****************      Методы для работы с тензором     *****************
	***********************************************************************/

	void Tensor::Transp()
	{
		double buf;
		for (int i = DIM-1; i >= 0; --i)
		{
			for (int j = i - 1; j >= 0; --j)
			{
				buf = C[i][j];
				C[i][j] = C[j][i];
				C[j][i] = buf;
			}
		}
	}

	double Tensor::getDet()
	{
		return (C[0][0] * C[1][1] * C[2][2] + C[0][1] * C[1][2] * C[2][0] +
			C[1][0] * C[2][1] * C[0][2] - C[2][0] * C[1][1] * C[0][2] -
			C[0][0] * C[1][2] * C[2][1] - C[2][2] * C[0][1] * C[1][0]);
	}

	void Tensor::set(double C00, double C01, double C02,
		double C10, double C11, double C12,
		double C20, double C21, double C22)
	{
		C[0][0] = C00;
		C[0][1] = C01;
		C[0][2] = C02;

		C[1][0] = C10;
		C[1][1] = C11;
		C[1][2] = C12;

		C[2][0] = C20;
		C[2][1] = C21;
		C[2][2] = C22;
	}

	void Tensor::setZero()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				C[i][j] = 0;
			}
		}
	}

	Vector Tensor::getRow(const int n)
	{
		Vector res;
		res.set(C[n][0], C[n][1], C[n][2]);
		return res;
	}

	Vector Tensor::getCol(const int n)
	{
		Vector res;
		res.set(C[0][n], C[1][n], C[2][n]);
		return res;
	}

	void Tensor::setUnit()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				C[i][j] = (i == j) ? 1.0 : 0;
			}
		}
	}

	double Tensor::doubleScalMult(Tensor t)
	{
		double res = 0;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res += C[i][j] * t.C[j][i];
			}
		}
		return res;
	}

	Tensor Tensor::getSymmetryPart()
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = C[i][j] + C[j][i];
			}
		}
		res /= 2.0;
		return res;
	}

	Tensor Tensor::getAntiSymmetryPart()
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = C[i][j] - C[j][i];
			}
		}
		res /= 2.0;
		return res;
	}
	
	Tensor::Tensor()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				C[i][j] = 0;
			}
		}
	}

	Tensor Tensor::operator + (const Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = C[i][j] + t.C[i][j];
			}
		}
		return res;
	}

	void Tensor::operator += (const Tensor t)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				C[i][j] += t.C[i][j];
			}
		}
	}

	Tensor Tensor::operator - (const Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = C[i][j] - t.C[i][j];
			}
		}
		return res;
	}

	Tensor Tensor::operator - ()
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = -C[i][j];
			}
		}
		return res;
	}

	void Tensor::operator -= (const Tensor t)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				C[i][j] -= t.C[i][j];
			}
		}
	}

	Tensor operator *(const double r, Tensor &t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = t.C[i][j]*r;
			}
		}
		return res;
	}

	Tensor operator *(Tensor &t, const double r)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = t.C[i][j] * r;
			}
		}
		return res;
	}

	Tensor Tensor::operator *(const Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					res.C[i][j] += C[i][k] * t.C[k][j];
				}
			}
		}
		return res;
	}

	void Tensor::operator *=(const Tensor t)
	{
		Tensor buf;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					buf.C[i][j] += C[i][k] * t.C[k][j];
				}
				C[i][j] = buf.C[i][j];
			}
		}
	}

	void Tensor::operator *= (const double r)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				C[i][j] *= r;
			}
		}
	}

	int Tensor::operator /= (const double r)
	{
		if (r != 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					C[i][j] /= r;
				}
			}
		}
		else
		{
			return -1;	//Код ошибки деления на 0
		}
		return 0;
	}

	Tensor::~Tensor()
	{
	}

	/***********************************************************************
	***********      Методы для работы с системой скольжения     ***********
	***********************************************************************/

	void SlipSystem::Initialize(float nx, float ny,
		float nz, float bx, float by, float bz)
	{
		/*
		Инициализация трёхиндексовых СС
		*/
		n.set(nx, ny, nz);
		b.set(bx, by, bz);
		
		n.Normalize();
		b.Normalize();
	}

	void SlipSystem::Initialize(float nx1, float nx2, float nx3, float nx4,
		float bx1, float bx2, float bx3, float bx4, float c)
	{
		/*
		Инициализация четырёхиндексовых СС
		*/
		n.set(nx1, nx2, nx4);

		double u = (2.0 * bx1 + bx2) / 3.0;
		double v = (bx1 + 2.0 * bx2) / 3.0;
		double w = bx4 / 3.0;
		b.set(u, v, w);

		n.Normalize();
		b.Normalize();
	}

	SlipSystem::SlipSystem()
	{
		t = 0;
		tc = 0;
		tbs = 0;
		dgm = 0;
		gmm = 0;
	}

	SlipSystem::~SlipSystem()
	{

	}

	/***********************************************************************
	********      Методы для работы с тензором четвёртого ранга    *********
	***********************************************************************/
	
	Tensor4::Tensor4()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] = 0;
					}
				}
			}
		}

	}

	Tensor4::~Tensor4()
	{

	}
	

	void Tensor4::setZero()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] = 0;
					}
				}
			}
		}
	}

	void Tensor4::Symmetrize()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						double buf = (C[i][j][k][l] + C[j][i][k][l] + C[i][j][l][k] + C[j][i][l][k]) / 4.0;
						C[i][j][k][l] = C[j][i][k][l] = C[i][j][l][k] = C[j][i][l][k] = buf;

					}
				}
			}
		}
	}
	
	Tensor4 Tensor4::ToLSK(Tensor O)
	{
		Tensor4 e;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						for (int i1 = 0; i1 < DIM; i1++)
						{
							for (int j1 = 0; j1 < DIM; j1++)
							{
								for (int k1 = 0; k1 < DIM; k1++)
								{
									for (int l1 = 0; l1 < DIM; l1++)
									{
										e.C[i][j][k][l] += O.C[i][i1] * O.C[j][j1] * O.C[k][k1] * O.C[l][l1] * C[i1][j1][k1][l1];
									}
								}
							}
						}
					}
				}
			}
		}
		return e;
	}

	void Tensor4::operator += (Tensor4 e)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] += e.C[i][j][k][l];
					}
				}
			}
		}
	}

	void Tensor4::operator -= (Tensor4 e)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] -= e.C[i][j][k][l];
					}
				}
			}
		}
	}

	void Tensor4::operator *= (double r)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] *= r;
					}
				}
			}
		}
	}

	int Tensor4::operator /= (double r)
	{
		if (r != 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < DIM; k++)
					{
						for (int l = 0; l < DIM; l++)
						{
							C[i][j][k][l] /= r;
						}
					}
				}
			}
		}
		else
		{
			return -1;
		}
		return 0;
	}

	/***********************************************************************
	***********   Прочие функции, не вошедшие в состав классов   ***********
	***********************************************************************/
	
	Tensor VectMult(Vector v, Tensor t)
	{
		Tensor res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getCol(k);
			Vector r = v.VectMult(buf);
			for (int i = 0; i < DIM; i++)
			{
				res.C[i][k] = r.C[i];
			}
		}
		return res;
	}
	
	Tensor VectMult(Tensor t, const Vector v)
	{
		Tensor res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getRow(k);
			Vector r = buf.VectMult(v);
			for (int i = 0; i < DIM; i++)
			{
				res.C[k][i] = r.C[i];
			}
		}
		return res;
	}

	Vector ScalMult(Vector v, Tensor t)
	{
		Vector res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getCol(k);
			res.C[k] = v.ScalMult(buf);
		}
		return res;
	}

	Vector ScalMult(Tensor t, const Vector v)
	{
		Vector res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getRow(k);
			res.C[k] = buf.ScalMult(v);
		}
		return res;
	}

	Tensor Transp(Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.C[i][j] = t.C[j][i];
			}
		}
		return res;
	}

	double Tensor::getL(int n)
	{
		Eigen::MatrixXf m(DIM, DIM);
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				m(i, j) = C[i][j];
			}
		}
		Eigen::EigenSolver<Eigen::MatrixXf> solver(m);
		return solver.eigenvalues().col(0)[n].real();
	}
}