#pragma once
#define MAX_UNK_PARA_NUM 12
#define uint unsigned int 

#include <vector>

#ifdef _WIN32
#include "io.h"
#include "windows.h"
#else
#include "sys/io.h"
#include "unistd.h"
#include "string.h"
#include "pthread.h"
#include <sys/stat.h>
#include <stdarg.h>

#define _access access
#define printf_s printf
#define sprintf_s sprintf
#define sscanf_s sscanf
#define fprintf_s fprintf
#define LPVOID void*
#define LPTHREAD_START_ROUTINE LPVOID
#define DWORD int
#define HANDLE pthread_t
#define INFINITE NULL
#define WAIT_OBJECT_0 0
#define FALSE 0
#define _fseeki64 fseeko64
#define _isnan isnan
//#define __int64 long long
typedef void* (*start_routine)(void*);
typedef long long __int64;
struct SYSTEMTIME
{
	unsigned int wYear;
	unsigned int wMonth;
	unsigned int wDay;
	unsigned int wHour;
	unsigned int wMinute;
	unsigned int wSecond;
};
struct MEMORYSTATUSEX
{
	int dwLength;
	long long ullAvailPhys;
};
struct SYSTEM_INFO
{
	int dwNumberOfProcessors;
};
inline int GetSystemInfo(SYSTEM_INFO* info) {
	info->dwNumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
	return 1;
}
inline int GlobalMemoryStatusEx(MEMORYSTATUSEX* sts) {
	long long pages = sysconf(_SC_AVPHYS_PAGES);
	long long page_size = sysconf(_SC_PAGE_SIZE);
	sts->ullAvailPhys = pages * page_size;

	return 1;
}
inline int fopen_s(FILE** fp, char* strFile, const char* mode)
{
	FILE* fpNew = NULL;
	fpNew = fopen(strFile, mode);
	*fp = fpNew;
	if (fpNew == NULL) {
		return -1;
	}
	else return 0;
}
inline int CreateDirectory(char* strPath, int mode)
{
	return mkdir(strPath, S_IRWXU) + 1;
}
inline double max(double a, double b) {
	if (a >= b) return a;
	else return b;
}
inline pthread_t CreateThread(int thread_attrib, int size, LPTHREAD_START_ROUTINE lpStartAddress, LPVOID lpParameter, DWORD dwCreationFlags, LPVOID lpThreadid)
{
	pthread_t t = pthread_t(lpThreadid);
	start_routine thread_fun = (start_routine)lpStartAddress;
	int test = pthread_create(&t, 0, thread_fun, lpParameter);
	return t;
}
inline int WaitForSingleObject(HANDLE handle, DWORD time)
{
	int test = pthread_join(handle, NULL);
	return test;
}
inline int CloseHandle(HANDLE handle) {
	return 1;
}

inline int GetLocalTime(SYSTEMTIME* st)
{
	time_t nSeconds;
	struct tm* pMt;

	time(&nSeconds);
	pMt = localtime(&nSeconds);
	st->wYear = pMt->tm_year + 1900;
	st->wMonth = pMt->tm_mon + 1;
	st->wDay = pMt->tm_mday;
	st->wHour = pMt->tm_hour;
	st->wMinute = pMt->tm_min;
	st->wSecond = pMt->tm_sec;

	return 1;
}

inline int CopyFile(char* strSrc, char* strDst, int mode)
{
	char buff[512] = { 0 };
	sprintf(buff, "cp %s %s", strSrc, strDst);
	system(buff);
	return 1;
}
inline double GetTickCount()
{
	double t;
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	t = ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
	return t;
}
#endif

typedef struct SUBMATRIXROWINDEX
{
	__int64 sNum;
	std::vector<int> veID;

	SUBMATRIXROWINDEX()
	{
		sNum = 0;
	}
}S_RINDEX;
typedef struct WEIGHTEDGE
{
	std::vector<int> veID;
	std::vector<double> veWei;
}WEDGE;
typedef struct SUBMATRIX
{
	int sR;
	int sC;
	int nR;
	int nC;

	double* pData;
	SUBMATRIX()
	{
		sR = sC = nR = nC = 0;
		pData = NULL;
	}
}SMA;
typedef struct ALIGNEDSUBMATRIX
{
	int sR;
	int sC;
	int nR;
	int nC;

	double pData[MAX_UNK_PARA_NUM * MAX_UNK_PARA_NUM];
	ALIGNEDSUBMATRIX()
	{
		sR = sC = nR = nC = 0;
		memset(pData, 0, sizeof(double) * MAX_UNK_PARA_NUM * MAX_UNK_PARA_NUM);
	}
}ALIGNED_SMA;
typedef struct OBJECTPOINT3D
{
	std::vector<int> veID;
}P3D;
struct THREAD_BMSV_CPU
{
	int nThreadID;
	double coef;
	unsigned int* pID;
	void* pBssm;
	SMA* pSma;
	int nHeight;
	double* pCst;
};
class CBlockSymmetrySparseMatrix
{
public:
	CBlockSymmetrySparseMatrix(void);
	~CBlockSymmetrySparseMatrix(void);
public:

	int FindIndex(int nID, std::vector<int>& veIndex, int& nNewID);
	int MultiMatrix(double* pA, int aR, int aC, double* pB, int bR, int bC, double* pAB);
	int TransMultiMatrix(double* pA, int aR, int aC, double* pB, int bR, int bC, double* pAB);
	int MultiTransMatrix(double* pA, int aR, int aC, double* pB, int bR, int bC, double* pAB);

	int ReadBinNormalEquation(char *strFile, std::vector<SMA>& veSMA, double* pCst);
	int ReadBinNormalEquation(FILE* fp, std::vector<SMA>& veSMA, double* pCst);

	int ResetSparseMatrix();
	int InitialSparseMatrix(uint nSmaNum, uint nGroupNum, uint* pBlockNum, uint* pBlockRow, S_RINDEX* psrIndex);
	

	int ComputeBlockIndex(S_RINDEX* psrIndex, int nPtNum, P3D* p3);
	int AddToSubMatrixIndex(S_RINDEX* psrIndex, int row, int col);
	int InitialSubMatrixIndex(S_RINDEX* psrIndex);

	int GetSparseBlockIDWithFullBlockID(int nBlockRowID, int nBlockColID);
	int GetEntryIDWithBlockID(int nBlockRowID, int nBlockColID, int& rowID, int& colID);

	int SparseMatrixToFull(uint nSmaNum, SMA* pS, uint nRow, uint nCol, double* pF);

	int UpdateBatchSparseMatrix(int sign, int nNum, SMA* pSub, SMA* pS);
	int UpdateBatchSparseMatrix(int sign, int nNum, SMA* pSub);
	int UpdateSparseMatrix(int sign, double* pSub, uint nSubR, uint nSubC, uint nBegR, uint nBegC, SMA* pS);
	int UpdateSparseMatrix(int sign, double* pSub, uint nSubR, uint nSubC, uint nBegR, uint nBegC);

	int SparseMatrixVectorProduct(uint nSmaNum, SMA* pSA, uint nRowA, uint nColA, double* pB, uint nRowB, double* pAB);
	int SparseUpperMatrixVectorProduct(double coef, uint nSmaNum, SMA* pSA, uint nRowA, uint nColA, double* pB, uint nRowB, double* pAB);
	int SparseUpperMatrixVectorProduct(double coef, int nStartID, int nNum, double* pB, uint nRowB, double* pAB);
	int SparseUpperMatrixVectorProduct_MultiThread(int nThreadNum, uint* pID, THREAD_BMSV_CPU* pB, double coef, int nSubMatNum, SMA* pSma, int nHeight, double* pCst, double* pRst);
	int SparseUpperMatrixVectorProduct_MultiThread(int nThreadNum, uint* pID, THREAD_BMSV_CPU* pB, double coef, int nHeight, double* pCst, double* pRst);
	int AddDampingCoef(double coef);
	SMA* GetBlockData(int blockID);
	SMA* GetBlockData(int rowID, int colID);
	SMA* GetData() { return m_pSma; };

private:

	uint m_nSize;
	uint m_nGroupNum;
	uint m_nBlockNum;
	
	uint* m_pBlockNum;
	uint* m_pBlockSize;

	uint m_nSmaNum;
	S_RINDEX* m_psrIndex;
	SMA* m_pSma;
};

