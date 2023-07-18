#include "CBlockSymmetrySparseMatrix.h"


double* pTmpRst = NULL;

void SparseUpperMatrixVectorProduct_Kernel(LPVOID lpParam)
{
	THREAD_BMSV_CPU* pB = (THREAD_BMSV_CPU*)lpParam;
	CBlockSymmetrySparseMatrix* pBssm = (CBlockSymmetrySparseMatrix*)pB->pBssm;
	double* pRst = pTmpRst + pB->nThreadID * pB->nHeight;
	pBssm->SparseUpperMatrixVectorProduct(pB->coef,  pB->pID[0], pB->pID[1], pB->pCst, pB->nHeight, pRst);
	return;
}
CBlockSymmetrySparseMatrix::CBlockSymmetrySparseMatrix()
{
	m_nSize = 0;
	m_nSmaNum = 0;
	m_nBlockNum = 0;
	m_nGroupNum = 1;

	m_pBlockNum = NULL;
	m_pBlockSize = NULL;

	m_psrIndex = NULL;
	m_pSma = NULL;
}
CBlockSymmetrySparseMatrix::~CBlockSymmetrySparseMatrix()
{
	if (m_pBlockNum) {
		delete[] m_pBlockNum;
		m_pBlockNum = NULL;
	}
	if (m_pBlockSize) {
		delete[] m_pBlockSize;
		m_pBlockSize = NULL;
	}

	if (m_pSma) {
		for (int i = 0; i < m_nSmaNum; i++) {
			delete[] m_pSma[i].pData;
			m_pSma[i].pData = NULL;
		}
		delete[] m_pSma;
		m_pSma = NULL;
	}
	if (m_psrIndex) {
		for (int i = 0; i < m_nBlockNum; i++) {
			std::vector<int>().swap(m_psrIndex[i].veID);
		}
		delete[] m_psrIndex;
		m_psrIndex = NULL;
	}
}
int CBlockSymmetrySparseMatrix::ReadBinNormalEquation(FILE* fp, std::vector<SMA>& veSMA, double* pCst)
{
	if (fp == NULL) {
		printf("Read binary sub normal equation file failed!%s\n");
		return 0;
	}
	__int64 num[2] = { 0 };
	fread(num, sizeof(__int64), 2, fp);
	__int64 i = 0; ALIGNED_SMA align_sma;
	veSMA.resize(num[0]);
	for (i = 0; i < num[0]; i++) {
		fread(&align_sma, sizeof(ALIGNED_SMA), 1, fp);
		veSMA[i].nR = align_sma.nR;
		veSMA[i].nC = align_sma.nC;
		veSMA[i].sR = align_sma.sR;
		veSMA[i].sC = align_sma.sC;
		if (align_sma.nR < 0 || align_sma.nR > MAX_UNK_PARA_NUM || align_sma.nC < 0 || align_sma.nC > MAX_UNK_PARA_NUM) {
			printf("Normal equation data unexpected!%s\n");
			continue;
		}
		veSMA[i].pData = new double[align_sma.nR * align_sma.nC];
		memcpy(veSMA[i].pData, align_sma.pData, sizeof(double) * align_sma.nR * align_sma.nC);
	}
	fread(pCst, sizeof(double), num[1], fp);
	return 1;
}
int CBlockSymmetrySparseMatrix::ReadBinNormalEquation(char* strFile, std::vector<SMA>& veSMA, double* pCst)
{
	FILE* fp = NULL;
	fopen_s(&fp, strFile, "rb");
	if (fp == NULL) {
		printf("Read binary sub normal equation file failed!%s\n", strFile);
		return 0;
	}
	ReadBinNormalEquation(fp, veSMA, pCst);

	fclose(fp);
	return 1;
}
int CBlockSymmetrySparseMatrix::ResetSparseMatrix()
{
	for (int i = 0; i < m_nSmaNum; i++) {
		if (m_pSma[i].pData) {
			memset(m_pSma[i].pData, 0, sizeof(double) * m_pSma[i].nR * m_pSma[i].nC);
		}

	}
	return 0;
}
int CBlockSymmetrySparseMatrix::InitialSparseMatrix(uint nSmaNum, uint nGroupNum, uint* pBlockNum, uint* pBlockSize, S_RINDEX* psrIndex)
{
	if (nGroupNum == 0 || pBlockNum == NULL || pBlockSize == NULL) {
		return -1;
	}
	m_nGroupNum = nGroupNum;
	m_pBlockNum = new uint[nGroupNum];
	m_pBlockSize = new uint[nGroupNum];
	memcpy(m_pBlockNum, pBlockNum, sizeof(uint) * nGroupNum);
	memcpy(m_pBlockSize, pBlockSize, sizeof(uint) * nGroupNum);

	int i = 0, nSize = 0, nBlock = 0;
	for (i = 0; i < nGroupNum; i++) {
		nSize += m_pBlockNum[i] * pBlockSize[i];
		nBlock += m_pBlockNum[i];
	}
	m_nSize = nSize;

	m_nBlockNum = nBlock;
	m_nSmaNum = nSmaNum;
	m_pSma = new SMA[m_nSmaNum];
	memset(m_pSma, 0, sizeof(SMA) * m_nSmaNum);
	m_psrIndex = new S_RINDEX[m_nBlockNum];
	memset(m_psrIndex, 0, sizeof(S_RINDEX) * m_nBlockNum);
	if (psrIndex) {
		InitialSubMatrixIndex(psrIndex);
	}

	return 0;
}
int CBlockSymmetrySparseMatrix::InitialSubMatrixIndex(S_RINDEX* psrIndex)
{
	int i = 0;
	if (m_psrIndex == NULL) {
		m_psrIndex = new S_RINDEX[m_nBlockNum];
		memset(m_psrIndex, 0, sizeof(S_RINDEX) * m_nBlockNum);
	}
	for (i = 0; i < m_nBlockNum; i++)
	{
		m_psrIndex[i].sNum = psrIndex[i].sNum;
		m_psrIndex[i].veID.assign(psrIndex[i].veID.begin(), psrIndex[i].veID.end());
	}
	return 0;
}
int CBlockSymmetrySparseMatrix::AddToSubMatrixIndex(S_RINDEX* psrIndex, int row, int col)
{
	if (col < row) {
		int t = col;
		col = row;
		row = t;
	}

	__int64 nStart = __int64(2 * m_nBlockNum + 1 - row) * row / 2;
	__int64 nID = nStart + col - row;
	__int64 nUpperSubMatrixNum = (__int64(m_nBlockNum) * m_nBlockNum + m_nBlockNum) / 2;


	if (nID >= nUpperSubMatrixNum) {
		printf("Error!Sub-normal-matrix ID out of range %d %d %d %d!\n", row, col, nID, nUpperSubMatrixNum);
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	int nNewCol = 0;
	int nFlag = FindIndex(col, psrIndex[row].veID, nNewCol);

	if (nFlag == 0) {
		if (nNewCol > psrIndex[row].veID.size() - 1) {
			psrIndex[row].veID.push_back(col);
		}
		else {
			psrIndex[row].veID.insert(psrIndex[row].veID.begin() + nNewCol, col);
		}
	}
	else if (nFlag == -1) {
		psrIndex[row].veID.push_back(col);
	}


	return 1;
}
int CBlockSymmetrySparseMatrix::ComputeBlockIndex(S_RINDEX* psrIndex, int nPtNum, P3D* p3)
{
	int i = 0, j = 0, k = 0, j1 = 0, j2 = 0, k1 = 0, k2 = 0, s1 = 0, s2 = 0;
	int rowID = 0, colID = 0, nImgNum = 0;
	uint* pSum = new uint[m_nGroupNum];
	memset(pSum, 0, sizeof(uint) * m_nGroupNum);
	for (i = 0; i < m_nGroupNum; i++) {
		pSum[i] = 0;
		for (j = 0; j < i; j++) {
			pSum[i] += m_pBlockNum[j];
		}
	}
	for (i = 0; i < nPtNum; i++) {
		nImgNum = p3[i].veID.size()/m_nGroupNum;

		for (j = 0; j < nImgNum*nImgNum; j++) {
			j1 = j / nImgNum;
			j2 = j % nImgNum;
			rowID = p3[i].veID[j1];
			colID = p3[i].veID[j2];
			for (k = 0; k < m_nGroupNum* m_nGroupNum; k++) {
				k1 = k / m_nGroupNum;
				k2 = k % m_nGroupNum;

				rowID = p3[i].veID[pSum[k1]+k1];
				colID = p3[i].veID[pSum[k2]+k2];
				AddToSubMatrixIndex(psrIndex, rowID, colID);
			}
		}
	}

	delete[] pSum; pSum = NULL;
	return 1;
}
int CBlockSymmetrySparseMatrix::FindIndex(int nID, std::vector<int>& veIndex, int& nNewID)
{
	int nNum = veIndex.size();
	if (nNum == 0) {
		return -1;
	}
	int nBegID = 0;
	int nEndID = nNum - 1;
	int nMidID = 0;
	if (nID < veIndex[nBegID]) {
		nNewID = nBegID;
		return 0;
	}
	else if (nID == veIndex[nBegID]) {
		nNewID = nBegID;
		return 1;
	}
	else if (nID > veIndex[nEndID]) {
		nNewID = nEndID + 1;
		return 0;
	}
	else if (nID == veIndex[nEndID]) {
		nNewID = nEndID;
		return 1;
	}

	int count = 0;
	while (nEndID - nBegID > 1)
	{

		if (nID == veIndex[nBegID]) {
			nNewID = nBegID;
			return 1;
		}
		else if (nID == veIndex[nEndID]) {
			nNewID = nEndID;
			return 1;
		}
		else {
			nMidID = (nBegID + nEndID) / 2;
			if (nID < veIndex[nMidID]) {
				nEndID = nMidID;
			}
			else if (nID == veIndex[nMidID]) {
				nNewID = nMidID;
				return 1;
			}
			else {
				nBegID = nMidID;
			}
		}
		count++;
	};
	nNewID = nEndID;
	return 0;
}
int CBlockSymmetrySparseMatrix::GetSparseBlockIDWithFullBlockID(int nBlockRowID, int nBlockColID)
{
	if (nBlockRowID < 0 || nBlockRowID > m_nBlockNum-1 || nBlockColID < 0 || nBlockColID > m_nBlockNum -1) {
		return -1;
	}
	int nNewID = 0;
	int nFlag = FindIndex(nBlockColID, m_psrIndex[nBlockRowID].veID, nNewID);
	if (nFlag == 1) {
		return m_psrIndex[nBlockRowID].sNum + nNewID;
	}
	else {
		return -1;
	}
}
int CBlockSymmetrySparseMatrix::GetEntryIDWithBlockID(int nBlockRowID, int nBlockColID, int& rowID, int& colID)
{
	if (nBlockRowID < 0 || nBlockRowID > m_nBlockNum - 1 || nBlockColID < 0 || nBlockColID > m_nBlockNum - 1) {
		return -1;
	}
	int i = 0, nBlockSum = 0, nBlockRowSum = 0, nBlockColSum = 0;
	int nEntryRowSum = 0, nEntryColSum = 0, nGroupRowID = -1, nGroupColID = -1;

	for (i = 0; i < m_nGroupNum; i++) {
		if (nGroupRowID == -1 ) {
			nBlockRowSum += m_pBlockNum[i];
			nEntryRowSum += m_pBlockNum[i] * m_pBlockSize[i];
			if (nBlockRowID < nBlockRowSum) {
				nGroupRowID = i;
				nBlockRowSum -= m_pBlockNum[i];
				nEntryRowSum -= m_pBlockNum[i] * m_pBlockSize[i];
			}
		}

		if (nGroupColID == -1 ) {
			nBlockColSum += m_pBlockNum[i];
			nEntryColSum += m_pBlockNum[i] * m_pBlockSize[i];
			if (nBlockColID < nBlockColSum) {
				nGroupColID = i;
				nBlockColSum -= m_pBlockNum[i];
				nEntryColSum -= m_pBlockNum[i] * m_pBlockSize[i];
			}
		}
	}

	rowID = nEntryRowSum + (nBlockRowID - nBlockRowSum) * m_pBlockSize[nGroupRowID];
	colID = nEntryColSum + (nBlockColID - nBlockColSum) * m_pBlockSize[nGroupColID];


	return 0;
}
int CBlockSymmetrySparseMatrix::SparseMatrixToFull(uint nSmaNum, SMA* pS, uint nRow, uint nCol, double* pF)
{
	int i = 0, j = 0, k = 0;
	int sR = 0, sC = 0;

	double* pData = new double[nRow * nCol];
	memset(pData, 0, sizeof(double) * nRow * nCol);
	for (i = 0; i < nSmaNum; i++) {
		SMA sma = pS[i];
		GetEntryIDWithBlockID(sma.sR, sma.sC, sR, sC);
		for (j = 0; j < sma.nR; j++) {
			for (k = 0; k < sma.nC; k++) {
				pData[(sR + j) * nCol + sC + k] = sma.pData[j * sma.nC + k];
			}
		}
	}

	for (i = 0; i < nRow; i++) {
		for (j = 0; j < nCol; j++) {
			if (i > j) {
				pData[i * nCol + j] = pData[j * nRow + i];
			}
		}
	}
	memcpy(pF, pData, sizeof(double) * nRow * nCol);

	delete[] pData;
	pData = NULL;
	return 1;
}
int CBlockSymmetrySparseMatrix::UpdateBatchSparseMatrix(int sign, int nNum, SMA* pSub, SMA* pS)
{
	double TotMem = 0;
	for (int i = 0; i < nNum; i++) {
		UpdateSparseMatrix(sign, pSub[i].pData, pSub[i].nR, pSub[i].nC, pSub[i].sR, pSub[i].sC, pS);

		if (pSub[i].pData) {
			TotMem += pSub[i].nR * pSub[i].nC * sizeof(double);
			pSub[i].nR = 0;
			pSub[i].nC = 0;
			delete[] pSub[i].pData;
			pSub[i].pData = NULL;
		}
	}
	return 1;
}
int CBlockSymmetrySparseMatrix::UpdateBatchSparseMatrix(int sign, int nNum, SMA* pSub)
{
	return UpdateBatchSparseMatrix(sign, nNum, pSub, m_pSma);
}
int CBlockSymmetrySparseMatrix::UpdateSparseMatrix(int sign, double* pSub, uint nSubR, uint nSubC, uint nBegR, uint nBegC, SMA* pS)
{
	uint i, j;
	if (pSub == NULL || pS == NULL) return -1;
	if (nBegR >= m_nBlockNum || nBegC >= m_nBlockNum) return -1;


	uint nID = GetSparseBlockIDWithFullBlockID(nBegR, nBegC);

	if (nID < 0 || nID >= m_nSmaNum) {
		//	WriteLog("Sub matrix ID out of range:%d %d %d!\n", nBegR, nBegC, nID);
		return -1;
	}

	if (pS[nID].pData == NULL) {

		pS[nID].sR = nBegR;
		pS[nID].sC = nBegC;
		pS[nID].nR = nSubR;
		pS[nID].nC = nSubC;
		pS[nID].pData = new double[nSubR * nSubC];
		memset(pS[nID].pData, 0, sizeof(double) * nSubR * nSubC);

	}

	for (i = 0; i < nSubR; i++) {
		for (j = 0; j < nSubC; j++) {
			pS[nID].pData[i * pS[nID].nC + j] += sign * pSub[i * nSubC + j];
		}
	}

	return 1;
}
int CBlockSymmetrySparseMatrix::UpdateSparseMatrix(int sign, double* pSub, uint nSubR, uint nSubC, uint nBegR, uint nBegC)
{
	return UpdateSparseMatrix(sign, pSub, nSubR, nSubC, nBegR, nBegC, m_pSma);
}
int CBlockSymmetrySparseMatrix::SparseUpperMatrixVectorProduct(double coef, uint nSmaNum, SMA* pSA, uint nRowA, uint nColA, double* pB, uint nRowB, double* pAB)
{
	if (pSA == NULL || pB == NULL || pAB == NULL) return -1;
	int i = 0, j = 0; int sR, sC, nSize = 0;
	memset(pAB, 0, sizeof(double) * nRowA);
	double pT[MAX_UNK_PARA_NUM * MAX_UNK_PARA_NUM] = {0};
	for (i = 0; i < nSmaNum; i++) {
		SMA sma = pSA[i];
		nSize = sma.nR * sma.nC;
		memcpy(pT, sma.pData, sizeof(double) * nSize);

		GetEntryIDWithBlockID(sma.sR, sma.sC, sR, sC);

		if (sma.sR == sma.sC) {
			for (j = 0; j < sma.nR; j++) {
				pT[j * sma.nC + j] *= (1 + coef);
			}
		}

		MultiMatrix(pT, sma.nR, sma.nC, pB + sC, sma.nC, 1, pT+nSize);
		for (j = 0; j < sma.nR; j++) { pAB[sR + j] += pT[nSize+j]; }

		if (sR != sC) {
			TransMultiMatrix(pT, sma.nR, sma.nC, pB + sR, sma.nR, 1, pT + nSize + sma.nR);
			for (j = 0; j < sma.nC; j++) { pAB[sC + j] += pT[nSize + sma.nR + j]; }
		}
	}
	return 0;
}
int CBlockSymmetrySparseMatrix::SparseUpperMatrixVectorProduct(double coef, int nStartID, int nNum, double* pB, uint nRowB, double* pAB)
{
	return SparseUpperMatrixVectorProduct(coef, nNum, m_pSma + nStartID, m_nBlockNum, m_nBlockNum, pB, nRowB, pAB);
}

int CBlockSymmetrySparseMatrix::SparseUpperMatrixVectorProduct_MultiThread(int nThreadNum, uint *pID, THREAD_BMSV_CPU * pB, double coef, int nSubMatNum, SMA* pSma, int nHeight, double* pCst, double* pRst)
{
	int i = 0, j = 0, k = 0;
//	nThreadNum = 1;
	memset(pRst, 0, sizeof(double) * nHeight);
	pTmpRst = new double[nThreadNum * nHeight];
	memset(pTmpRst, 0, sizeof(double) * nThreadNum * nHeight);

	double *pTmpCst = new double[nThreadNum * nHeight];
	memset(pTmpCst, 0, sizeof(double) * nThreadNum * nHeight);

	HANDLE* pThread = new HANDLE[nThreadNum];
	memset(pThread, 0, sizeof(HANDLE) * nThreadNum);

	for (i = 0; i < nThreadNum; i++) {
		DWORD lpThreadID = i;
		pB[i].nThreadID = i;
		pB[i].nHeight = nHeight;
		pB[i].coef = coef;
		pB[i].pCst = pCst;
		pB[i].pSma = pSma;
		pB[i].pID = pID + 2 * i;
		pB[i].pBssm = this;
		//	printf("Sub task %d: %d/%d %d/%d\n", i, pID[2*i+0], nSubMatNum, pID[2*i+1], nSubMatNum);
		pThread[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)SparseUpperMatrixVectorProduct_Kernel, pB + i, 0, &lpThreadID);

	}
	for (i = 0; i < nThreadNum; i++) {
		DWORD dwRet = WaitForSingleObject(pThread[i], INFINITE);
		if (dwRet == WAIT_OBJECT_0)
		{
			//	printf("Thread %d completed task\n", i);
		}
		CloseHandle(pThread[i]);

		for (j = 0; j < nHeight; j++) {
			pRst[j] += pTmpRst[i * nHeight + j];
		}
	}
	delete[] pThread; pThread = NULL;
	delete[] pTmpCst; pTmpCst = NULL;
	delete[] pTmpRst; pTmpRst = NULL;

	return 1;
}
int CBlockSymmetrySparseMatrix::SparseUpperMatrixVectorProduct_MultiThread(int nThreadNum, uint* pID, THREAD_BMSV_CPU* pB, double coef, int nHeight, double* pCst, double* pRst)
{
	return SparseUpperMatrixVectorProduct_MultiThread(nThreadNum, pID, pB, coef, m_nSmaNum, m_pSma, nHeight, pCst, pRst);
}
int CBlockSymmetrySparseMatrix::SparseMatrixVectorProduct(uint nSmaNum, SMA* pSA, uint nRowA, uint nColA, double* pB, uint nRowB, double* pAB)
{

	if (pSA == NULL || pB == NULL || pAB == NULL) return -1;
	int i = 0, j = 0; int sR, sC;
	memset(pAB, 0, sizeof(double) * nRowA);

	for (i = 0; i < nSmaNum; i++) {
		SMA sma = pSA[i];

		double* pT = new double[sma.nR];
		memset(pT, 0, sizeof(double) * sma.nR);
		GetEntryIDWithBlockID(sma.sR, sma.sC, sR, sC);

		MultiMatrix(sma.pData, sma.nR, sma.nC, pB + sC, sma.nC, 1, pT);
		for (j = 0; j < sma.nR; j++) { pAB[sR + j] += pT[j]; }

		delete[] pT; pT = NULL;
	}

	return 0;
}
int CBlockSymmetrySparseMatrix::AddDampingCoef(double coef)
{
	int i = 0, j = 0;

	for (i = 0; i < m_nBlockNum; i++) {

		int nID = GetSparseBlockIDWithFullBlockID(i, i);

		if (nID < 0) {
			printf("Error!EOP Sub Matrix %d %d ID invalid!\n", i, i);
			return 0;
		}
		int nR = m_pSma[nID].nR;
		int nC = m_pSma[nID].nC;

		if (nR != nC || nR == 0 || nC == 0) {
			printf("Error!EOP Sparse normal matrix size invalid:%d-%d-%d!\n", i, nR, nC);
			return 0;
		}
		for (j = 0; j < nR; j++) {

			m_pSma[nID].pData[j * nR + j] *= coef + 1;
		};

	}
		
	return 0;
}
SMA* CBlockSymmetrySparseMatrix::GetBlockData(int blockID)
{
	if (blockID < 0 || blockID > m_nSmaNum - 1 )
	{
		return NULL;
	}
	return m_pSma + blockID;
}
SMA * CBlockSymmetrySparseMatrix::GetBlockData(int rowID, int colID)
{
	if (rowID < 0 || rowID > m_nBlockNum - 1 || colID < 0 || colID > m_nBlockNum - 1)
	{
		return NULL;
	}
	int nID = GetSparseBlockIDWithFullBlockID(rowID, colID);
	return m_pSma+nID;
}
int CBlockSymmetrySparseMatrix::MultiMatrix(double* pA, int aR, int aC, double* pB, int bR, int bC, double* pAB)
{
	if (pA == 0 || pB == 0 || pAB == 0 || aC != bR) return -1;
	memset(pAB, 0, sizeof(double) * aR * bC);
	int i = 0; int j = 0; int k = 0;
	for (i = 0; i < aR; i++) {
		for (j = 0; j < bC; j++) {
			for (k = 0; k < aC; k++) {
				pAB[i * bC + j] += pA[i * aC + k] * pB[k * bC + j];
			}
		}
	}
	return 0;
}

int CBlockSymmetrySparseMatrix::TransMultiMatrix(double* pA, int aR, int aC, double* pB, int bR, int bC, double* pAB)
{
	if (pA == 0 || pB == 0 || pAB == 0 || aR != bR) return -1;
	memset(pAB, 0, sizeof(double) * aC * bC);
	int i = 0; int j = 0; int k = 0;
	for (i = 0; i < aC; i++) {
		for (j = 0; j < bC; j++) {
			for (k = 0; k < aR; k++) {
				pAB[i * bC + j] += pA[k * aC + i] * pB[k * bC + j];
			}
		}
	}
	return 0;
}

int CBlockSymmetrySparseMatrix::MultiTransMatrix(double* pA, int aR, int aC, double* pB, int bR, int bC, double* pAB)
{
	if (pA == 0 || pB == 0 || pAB == 0 || aC != bC) return -1;
	memset(pAB, 0, sizeof(double) * aR * bR);
	int i = 0; int j = 0; int k = 0;
	for (i = 0; i < aR; i++) {
		for (j = 0; j < bR; j++) {
			for (k = 0; k < aC; k++) {
				pAB[i * bR + j] += pA[i * aC + k] * pB[j * bC + k];
			}
		}
	}
	return 0;
}