#include "stdafx.h"
#include "params.h"

namespace model
{
	int get1DPos(int q1, int q2, int q3)
	{
		int res = q1*fragm_count*fragm_count + q2*fragm_count + q3;
		return res;
	}

	void get3DPos(int pos, int* q1, int* q2, int* q3)
	{
		int C2d = fragm_count*fragm_count;
		int C3d = C2d*fragm_count;
		
		int qq1 = pos / C2d;
		int qq2 = (pos - qq1*C2d) / fragm_count;
		int qq3 = (pos - qq1*C2d) % fragm_count;

		*q1 = qq1;
		*q2 = qq2;
		*q3 = qq3;
	}
}