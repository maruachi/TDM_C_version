#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct _cubeInfo {
	int 	N[3];
	int 	numAtom;
	double 	a[3], b[3], c[3];
	double 	infVol;
} cubeInfo, *pCubeInfo;

//parse cube file
pCubeInfo 		readCubeInfo(char* filein);
double*** 		readCubeWf(char* filein, pCubeInfo cInfo);

//fucntions to calculate the potential
double 			nomalization(double*** wf, pCubeInfo cInfo);
void 			findPositionFromIndex(int i, int j, int k, double* position, pCubeInfo cInfo);
double*			calStaticDipole(double*** wf1, pCubeInfo cInfo);

//tools for vector calculation
double 			calVolume(double* a, double* b, double* c);

//tools to process data
double*** 		dynAllocate3D(int size1, int size2, int size3);
double 			sciNotToDouble(char* str);
void 			lineSkip(int number, FILE* fCube);

//print functions to test this program
void printcInfo(pCubeInfo cInfo);
void printWf(double*** wf, pCubeInfo cInfo, char* text);

int main(int argc, char* argv[])
{
	char* fileWf1;
	pCubeInfo cInfo;
	double*** wf1;
	double to_Debye;

	if (argc != 2)
	{
		printf("static_dipole [charge density in cube file]\n");
		printf("this program calculate the dipole matrix integral(r rho(r))\n");
	}
	
	fileWf1 = argv[1];
	cInfo = readCubeInfo(fileWf1);

	wf1 = readCubeWf(fileWf1, cInfo);
	

	double* dipoles = calStaticDipole(wf1, cInfo);

	printf("--dipole matrix elements--\n");
	printf("X matrix : %.6lf (# of electrons*bohr)\n", dipoles[0]);
	printf("Y matrix : %.6lf (# of electrons*bohr)\n", dipoles[1]);
	printf("Z matrix : %.6lf (# of electrons*bohr)\n", dipoles[2]);

	to_Debye = 1.6 * 5.2918 / 3.33564;

	printf("--dipole matrix elements--\n");
	printf("X matrix : %.6lf Debye\n", dipoles[0]*to_Debye);
	printf("Y matrix : %.6lf Debye\n", dipoles[1]*to_Debye);
	printf("Z matrix : %.6lf Debye\n", dipoles[2]*to_Debye);

	return 0;
}

//parse cube file
pCubeInfo readCubeInfo(char* filein)
{
	pCubeInfo cInfo = (pCubeInfo)malloc(sizeof(cubeInfo));
	int numAtom = 0;
	int Nx = 0;
	int Ny = 0;
	int Nz = 0;
	FILE* fCube;
	
	if(!(fCube = fopen(filein, "r")))
	{
		printf("error2 : the file can not be opend.\n");
		return NULL;
	}
	
	lineSkip(2, fCube); // the ahead tow lines are just comments.

	fscanf(fCube, "%d", &numAtom);
	cInfo->numAtom = numAtom;
	
	lineSkip(1, fCube); // the next to the number of atom is the origin.
	
	fscanf(fCube, "%d",  &(cInfo->N[0]));
	fscanf(fCube, "%lf", &(cInfo->a[0]));
	fscanf(fCube, "%lf", &(cInfo->a[1]));
	fscanf(fCube, "%lf", &(cInfo->a[2]));
	fscanf(fCube, "%d",  &(cInfo->N[1]));
	fscanf(fCube, "%lf", &(cInfo->b[0]));
	fscanf(fCube, "%lf", &(cInfo->b[1]));
	fscanf(fCube, "%lf", &(cInfo->b[2]));
	fscanf(fCube, "%d",  &(cInfo->N[2]));
	fscanf(fCube, "%lf", &(cInfo->c[0]));
	fscanf(fCube, "%lf", &(cInfo->c[1]));
	fscanf(fCube, "%lf", &(cInfo->c[2]));
	cInfo->infVol = calVolume(cInfo->a, cInfo->b, cInfo->c);
	
	return cInfo;
}

double*** readCubeWf(char* filein, pCubeInfo cInfo)
{
	int Nx = cInfo->N[0];
	int Ny = cInfo->N[1];
	int Nz = cInfo->N[2];
	int numAtom = cInfo->numAtom;
	char temp[20];
	double*** wf = dynAllocate3D(Nx, Ny, Nz);
	FILE* fCube;
	
	if(!(fCube = fopen(filein, "r")))
	{
		printf("error2 : the file can not be opend.\n");
		return NULL;
	}
	
	lineSkip(2, fCube); // the ahead tow lines are just comments.
	lineSkip(1, fCube); // coordinate of the origin
	lineSkip(3, fCube); // unit vectors of infinitisimal volume
	lineSkip(numAtom, fCube); // charge and coordinates of atoms
	
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				fscanf(fCube, "%s", temp);
				if (sciNotToDouble(temp) >= 0)
				{
					wf[i][j][k] = sqrt(fabs(sciNotToDouble(temp)));
				}
				else
				{
					wf[i][j][k] = -sqrt(fabs(sciNotToDouble(temp)));
				}
			}
		}
	}
	
	return wf;
}

//fucntions to calculate the potential
void findPositionFromIndex(int i, int j, int k, double* position, pCubeInfo cInfo)
{
    for (int n = 0; n < 3; n++)
    {
        position[n] = cInfo->a[n] * i + cInfo->b[n] * j + cInfo->c[n] * k;
    }
    
    return;
}


double* calStaticDipole(double*** wf1, pCubeInfo cInfo)
{
	int Nx = cInfo->N[0];
	int Ny = cInfo->N[1];
	int Nz = cInfo->N[2];
	double infVol = cInfo->infVol;
	double position[3];
	double Xdipole = 0;
	double Ydipole = 0;
	double Zdipole = 0;
	double* dipoles = (double*) calloc(3, sizeof(double));
	double origin[3];
	
	findPositionFromIndex(Nx/3, Ny/3, Nz/3, origin, cInfo);
	

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				findPositionFromIndex(i, j, k, position, cInfo);
				for (int m = 0; m < 3; m++)
				{
					position[m] = position[m] - origin[m];
				}
				Xdipole += wf1[i][j][k] * position[0] * infVol;
				Ydipole += wf1[i][j][k] * position[1] * infVol;
				Zdipole += wf1[i][j][k] * position[2] * infVol;
			}
		}
	}


	dipoles[0] = Xdipole;
	dipoles[1] = Ydipole;
	dipoles[2] = Zdipole;

	return dipoles;
}

//tools for vector calculation
double calVolume(double* a, double* b, double* c)
{
	double volume = 0;
	double cross[3];

	cross[0] = a[1]*b[2] - a[2]*b[1];
	cross[1] = a[2]*b[0] - a[0]*b[2];
	cross[2] = a[0]*b[1] - a[1]*b[0];
	
	for (int i = 0; i < 3; i++)
	{
		volume += cross[i] * c[i];
	}
	volume = fabs(volume);

	return volume;
}

//tools to process data
double*** dynAllocate3D(int size1, int size2, int size3)
{
	double*** arr;

	arr = (double***)calloc(size1, sizeof(double**));
	for(int i = 0; i < size1; i++)
	{
		arr[i] = (double**)calloc(size2, sizeof(double*));
		for (int j = 0; j < size2; j++)
		{
			arr[i][j] = (double*)calloc(size3, sizeof(double));
		}
	}

	return arr;
}

double sciNotToDouble(char* str)
{
	int length = 0;
	int flag = 0;
	int det = 0;
	int sign = 0; // 0: plus, 1: minus
	char mantissa[20];
	char exponent[20];
	double mant, exp;
	double doubleVal = 0;
	
	length = strlen(str);
	for (int i = 0; i <= length; i++)
	{	
		if (flag == 0) 
		{
			if (str[i] != 'E')
			{
				mantissa[i] = str[i];
			}
			else
			{
				flag = 1;
				det = i;
				mantissa[i] = '\0';
			}
		}
		else
		{
			exponent[i-det-1] = str[i];
		}
	}

	doubleVal = atof(mantissa)*pow(10, atof(exponent));

	return doubleVal;
}

void lineSkip(int number, FILE* fCube)
{
	char temp[200];

	for (int i = 0; i < number; i++)
	{
		fgets(temp, 100, fCube);
	}

	return;
}

void printWf(double*** wf, pCubeInfo cInfo, char* text)
{
	int Nx = cInfo->N[0];
	int Ny = cInfo->N[1];
	int Nz = cInfo->N[2];
	FILE* p;
	
	p = fopen(text, "w");
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				fprintf(p, "%lf ", wf[i][j][k]);
				if(k % 6 == 5) fprintf(p, "\n");
			}
		}
	}

	fclose(p);
}

void printcInfo(pCubeInfo cInfo)
{
	FILE* p;
	
	p = fopen("cInfoTest.txt", "w");
	fprintf(p, "%d ", cInfo->numAtom);
	fprintf(p, "%lf ", cInfo->infVol);
	for (int i = 0; i < 3; i++)
	{
		fprintf(p, "%d ", cInfo->N[i]);
	}
	fprintf(p, "\n");
	
	for (int i = 0; i < 3; i++)
	{
		fprintf(p, "%lf ", cInfo->a[i]);
	}
	fprintf(p, "\n");
	
	for (int i = 0; i < 3; i++)
	{
		fprintf(p, "%lf ", cInfo->b[i]);
	}
	fprintf(p, "\n");
	
	for (int i = 0; i < 3; i++)
	{
		fprintf(p, "%lf ", cInfo->c[i]);
	}
	fprintf(p, "\n");

	fclose(p);
}
