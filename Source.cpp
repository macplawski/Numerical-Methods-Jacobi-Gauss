#include<iostream>
#include<cmath>
#include<ctime>

using namespace std;

const int N = 967;

void zeruj_wektor(double *wektor)
{
	for (int i = 0; i < N; i++)
	{
		wektor[i]=1;
	}
}

double oblicz_residuum(double macierzA[N][N], double *wektorX, double *wektorB)
{
	static double tmp[N];
	double res=0;

	for (int i = 0; i < N; i++)
	{
		double suma = 0;
		for (int j = 0; j < N; j++)
		{
			suma += macierzA[i][j] * wektorX[j];
		}

		suma -= wektorB[i];
		
		tmp[i] = suma;

	}

	for (int i = 0; i < N; i++)
	{
		res += tmp[i] * tmp[i];
	}



	return sqrt(res);

}

int jacobbi(double macierzA[N][N], double *wektorB, double *wektorX)
{
	double res;
	int iteracje=0;
	static double tmp[N];

	do 
	{
		iteracje++;

		for (int i = 0; i < N; i++)
		{
			tmp[i] = wektorX[i];
		}

		for (int i = 0; i < N; i++)
		{
			double suma = wektorB[i];

			for (int j = 0; j < i; j++)
			{
				suma -= macierzA[i][j] * tmp[j];
			}
			for (int j = i+1; j < N; j++)
			{
				suma -= macierzA[i][j] * tmp[j];
			} 
			suma /= macierzA[i][i];
			wektorX[i] = suma;
		}
	
		res = oblicz_residuum(macierzA, wektorX, wektorB);		
	
	} while (res >= 1e-9 && res < INFINITY); //10^-9
	
	if (res < 1e-9)
		return iteracje;
	else
		return INT_MAX;
}

int gauss(double macierzA[N][N], double *wektorB, double *wektorX)
{
	double res = 1;
	int iteracje = 0;

	while (res > 1e-9 && res<INFINITY) //10^-9
	{
		iteracje++;
		for (int i = 0; i < N; i++)
		{
			double suma = 0, suma2 = 0;
			for (int j = 0; j < N; j++)
			{
				if (j != i)
					suma += macierzA[i][j] * wektorX[j];
			}
			wektorX[i] = (wektorB[i] - suma) / macierzA[i][i];
		}

		res = oblicz_residuum(macierzA, wektorX, wektorB);

	}

	if (res < 1e-9)
		return iteracje;
	else
		return INT_MAX;
}

double faktoryzacjaLU(double macierzA[N][N], double *wektorB, double *wektorX)
{
	static double L[N][N], U[N][N];
	static double wektorD[N];
	double res;

	for (int i = 0; i < N; i++)
	{
		for (int h = 0; h < N; h++)
		{
			if (h == i)
				L[i][h] = 1;

			U[i][h] = macierzA[i][h];
		}
	}



	for (int i = 0; i < N-1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			L[j][i] = U[j][i] / U[i][i];

			for (int m = i; m < N; m++)
				U[j][m] -= L[j][i] * U[i][m];
			
		}
	}

	//podstawianie wprzod

	for (int i = 0; i < N; i++)
	{
		wektorD[i] = wektorB[i];
		for (int j = 0; j < i; j++)
		{
			wektorD[i] -= L[i][j] * wektorD[j];
		}

		wektorD[i] /= L[i][i];
	}

	// podstawianie wstecz

	for (int i = N - 1; i >= 0; i--)
	{
		wektorX[i] = wektorD[i];
		for (int j = i+1; j < N; j++)
		{
			wektorX[i] -= U[i][j] * wektorX[j];
		}

		wektorX[i] /= U[i][i];
	}

	res = oblicz_residuum(macierzA, wektorX, wektorB);

	return res;
}

void zapeln_macierz(double macierzA[N][N], double *wektorB, int a1, int a2, int a3)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			if (i == j)
			{
				macierzA[i][j] = a1;
			}
			else if (i - j >= -2 && i - j <= 2)
			{
				macierzA[i][j] = a2;
			}
			/*else if (i == (j + 1) || (i + 1) == j)
			{
				macierzA[i][j] = a2;
			}
			else if (i == (j + 2) || (i + 2) == j)
			{
				macierzA[i][j] = a3;
			}*/

		}

		wektorB[i] = sin((i + 1) * 5);
	}
}

int main()
{
	int a1 = 5+1, a2 = -1, a3 = -1;
	static double wektorB[N];
	static double macierzA[N][N];
	static double residuum[N];
	static double wektorX[N];
	clock_t start, stop;
	double czas_j, czas_g;
	int iteracje;
	double res;

	//zadanie A

	zapeln_macierz(macierzA, wektorB, a1, a2, a3);
	zeruj_wektor(wektorX);

	//zadanie B

	cout << "	Zadanie B\n\n";

	start = clock();
	iteracje = jacobbi(macierzA, wektorB, wektorX);
	stop = clock();
	czas_j = (double)(stop - start) / CLOCKS_PER_SEC;
	if (iteracje == INT_MAX)
		cout << "Metoda sie nie zbiega!\n";
	else
		cout << "Jacobbi potrzebuje " << iteracje << " iteracji i wykonuje sie w czasie " << czas_j <<endl;
	
	zeruj_wektor(wektorX);
	start = clock();
	iteracje = gauss(macierzA, wektorB, wektorX);
	stop = clock();
	czas_g = (double)(stop - start) / CLOCKS_PER_SEC;
	if (iteracje == INT_MAX)
		cout << "Metoda sie nie zbiega!\nBrak mozliwosci porownania dwoch metod.\n\n";
	else
	{
		cout << "Gauss potrzebuje " << iteracje << " iteracji i wykonuje sie w czasie " << czas_g << endl;
		cout << "\nMetoda Gaussa jest szybsza o " << czas_j - czas_g<<"s\n\n";
	}
	
	start = clock();
	iteracje = faktoryzacjaLU(macierzA, wektorB, wektorX);
	stop = clock();
	czas_j = (double)(stop - start) / CLOCKS_PER_SEC;
	if (iteracje == INT_MAX)
		cout << "Metoda sie nie zbiega!\n\n";
	else
		cout << "Faktoryzacja LU potrzebuje " << iteracje << " iteracji i wykonuje sie w czasie " << czas_j << endl<<endl;

	//zadanie C
	cout << "	Zadanie C\n\n";

	zeruj_wektor(wektorX);
	a1 = 3;
	zapeln_macierz(macierzA, wektorB, a1, a2, a3);
	iteracje = jacobbi(macierzA, wektorB, wektorX);
	if (iteracje == INT_MAX)
		cout << "Dla a1=3 metoda jacobiego sie nie zbiega!\n";
	else
		cout << "Dla a1=3 metoda jacobiego sie zbiega\n\n";

	zeruj_wektor(wektorX);
	iteracje = gauss(macierzA, wektorB, wektorX);
	if (iteracje == INT_MAX)
		cout << "Dla a1=3 metoda gaussa sie nie zbiega!\n\n";
	else
		cout << "Dla a1=3 metoda gaussa sie zbiega\n\n";
	
	//zadanie D

	cout << "	Zadanie D\n\n";

	zeruj_wektor(wektorX);
	res = faktoryzacjaLU(macierzA, wektorB, wektorX);
	cout << "Norma residuum dla faktoryzacji wynosi " << res << endl << endl;


	return 0;
}