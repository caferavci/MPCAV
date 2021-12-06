template <typename T>
T* AllocateDynamicVector(int nRows)
{
	T* dynamicVector;

	dynamicVector = new (std::nothrow) T[nRows]();

	if (dynamicVector == NULL)
	{
		exit(1);

	}
	return dynamicVector;
}

template <typename T>
void DeallocateDynamicVector(T* dVector, int nRows)
{
	if (!dVector)
		return;
	delete[] dVector;
}

template <typename T>
T** AllocateDynamicArray(int nRows, int nCols)
{
	T** dynamicArray;

	dynamicArray = new (std::nothrow) T * [nRows]();

	if (dynamicArray == NULL)
	{
		exit(1);
	}

	for (int i = 0; i < nRows; i++)
	{
		dynamicArray[i] = new (std::nothrow) T[nCols]();

		if (dynamicArray[i] == NULL)
		{
			exit(1);
		}
	}

	return dynamicArray;
}

template <typename T>
void DeallocateDynamicArray(T** dArray, int nRows, int nCols)
{
	if (!dArray)
		return;

	for (int x = 0; x < nRows; x++)
	{
		delete[] dArray[x];
	}

	delete[] dArray;

}

template <typename T>
T*** Allocate3DDynamicArray(int nX, int nY, int nZ)
{
	T*** dynamicArray;

	dynamicArray = new (std::nothrow) T * *[nX]();

	if (dynamicArray == NULL)
	{
		exit(1);
	}

	for (int x = 0; x < nX; x++)
	{
		dynamicArray[x] = new (std::nothrow) T * [nY]();

		if (dynamicArray[x] == NULL)
		{
			exit(1);
		}

		for (int y = 0; y < nY; y++)
		{
			dynamicArray[x][y] = new (std::nothrow) T[nZ]();
			if (dynamicArray[x][y] == NULL)
			{
				exit(1);
			}
		}
	}

	return dynamicArray;

}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
	if (!dArray)
		return;
	for (int x = 0; x < nX; x++)
	{
		for (int y = 0; y < nY; y++)
		{
			delete[] dArray[x][y];
		}

		delete[] dArray[x];
	}

	delete[] dArray;

}

template <typename T>
T**** Allocate4DDynamicArray(int nM, int nX, int nY, int nZ)
{
	T**** dynamicArray;

	dynamicArray = new (std::nothrow) T * **[nX]();

	if (dynamicArray == NULL)
	{
		exit(1);
	}
	for (int m = 0; m < nM; m++)
	{
		if (m % 100 == 0)
			//cout << "allocating 4D memory for " << m << endl;

			dynamicArray[m] = new (std::nothrow) T * *[nX]();

		if (dynamicArray[m] == NULL)
		{
			exit(1);
		}

		for (int x = 0; x < nX; x++)
		{
			dynamicArray[m][x] = new (std::nothrow) T * [nY]();

			if (dynamicArray[m][x] == NULL)
			{
				exit(1);
			}

			for (int y = 0; y < nY; y++)
			{
				dynamicArray[m][x][y] = new (std::nothrow) T[nZ]();
				if (dynamicArray[m][x][y] == NULL)
				{
					exit(1);
				}
			}
		}
	}
	return dynamicArray;

}

template <typename T>
void Deallocate4DDynamicArray(T**** dArray, int nM, int nX, int nY)
{
	if (!dArray)
		return;
	for (int m = 0; m < nM; m++)
	{
		for (int x = 0; x < nX; x++)
		{
			for (int y = 0; y < nY; y++)
			{
				delete[] dArray[m][x][y];
			}
			delete[] dArray[m][x];
		}
		delete[] dArray[m];
	}
	delete[] dArray;

}