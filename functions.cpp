// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <fstream>
#include "windows.h"

#include "functions.h"

namespace model
{
	bool isDirectoryExists(LPCWSTR filename)
	{
		DWORD dwFileAttributes = GetFileAttributes(filename);
		if (dwFileAttributes == 0xFFFFFFFF)
			return false;
		return dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
	}

	void WriteDebugInfo(std::ofstream& Stream, double Matrix[3][3])
	{
		for (int i = 0; i < 3; i++)
		{
			Stream << Matrix[i][0] << " " << Matrix[i][1] << " " << Matrix[i][2] << std::endl;
		}
		Stream << std::endl;
	}
	void TruncSSTFiles()
	{
		std::ofstream of1;
		of1.open("Polus\\SST001.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of1.close();
		of1.open("Polus\\SST011.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of1.close();
		of1.open("Polus\\SST111.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of1.close();
	}

	void TruncPoleFiles()
	{
		std::ofstream of;
		of.open("Polus\\S001.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S010.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S100.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();

		of.open("Polus\\S011.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S110.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S101.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S01-1.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S1-10.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S10-1.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();

		of.open("Polus\\S1-11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S-111.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S11-1.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
		of.open("Polus\\S111.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		of.close();
	}

}