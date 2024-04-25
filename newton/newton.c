#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include"tinyexpr.h"
#include "tinyexpr.c"

#include "from_pastebin_com.h"

// Задать:
// 1. формулы: строки
// 2. начальное приближение в виде: x1,x2,...,xN

// Получить: решение - вектор в виде x1,x2,...,xN

// Dimension for Compile Time example
#define DIM_FOR_COMPILE_TIME_EXAMPLE 2

bool bRuntime = true;

//////////////////////////////////////////////////////////////

// массив строк, задающих формулы функций т.е. типа "x1^2 + x2^2 - 1"
char** formulae = NULL;

char** createFormulaeArray(int count)
{
    if (count > 0)
    {
        formulae = (char**)malloc(count * sizeof(char*));
    
        char *pFormula = formulae[0];
    
        for (int i = 0; i < count; i++)
        {
            int formulaLength = 0; 

            // TODO retrieve formula length
            // formulaLength = ....;
            pFormula = (char*)malloc(formulaLength + 1);
             
            // TODO retrieve and fill the formula string as 0-terminated
            // ....(pFormula)...
            
            formulae[i] = pFormula;
            pFormula += strlen(pFormula);
        }
    }
    
    return formulae;
}

int getFormulae()
{
    int dim = 0;
    
    // TODO - retrieve number of functions i.e. the dimension
    // int dim = ...
    
    formulae = createFormulaeArray(dim);
    
    return dim;
}

//////////////////////////////////////////////////////////////

typedef double* tVectorValue;

// compile-time example of the functions
double getFuncValue_CompileTime(int fIndex, tVectorValue vectorX)
{
    double funcValue = 0;
    
    switch (fIndex)
    {
        case 0:
        {
            // Example from task pdf. Has a solution.
            // f1: x1^2 + x2^2 - 1
            funcValue = vectorX[0] * vectorX[0] + vectorX[1] * vectorX[1] - 1;
            
            // Example. Has a solution.
            // f1: x1^2 + (x2 - 1.25)^2 - 1
            // funcValue = vectorX[0] * vectorX[0] + (vectorX[1] - 1.25) * (vectorX[1] - 1.25) - 1;            
            
            // Example. Does NOT have a solution.
            // f1: x1^2 + (x2 - 1.26)^2 - 1
            // funcValue = vectorX[0] * vectorX[0] + (vectorX[1] - 1.26) * (vectorX[1] - 1.26) - 1;
            break;
        }
            
        case 1:
        {           
            // f2: x1^2 - x2
            funcValue = vectorX[0] * vectorX[0] - vectorX[1];
            break;
        }
        
        default: 
            break;
    }
    
    return funcValue;
}

/*double getFuncValue_CompileTime(int fIndex, tVectorValue vectorX)
{
    double funcValue = 0;

    switch (fIndex)
    {
        case 0:
        {
            // Example from task pdf. Has a solution.
            // f1: 3x1+2x2-x3-13
            funcValue = vectorX[0] *3 + vectorX[1] * 2 - vectorX[2]-13;

            // Example. Has a solution.
            // f1: x1^2 + (x2 - 1.25)^2 - 1
            // funcValue = vectorX[0] * vectorX[0] + (vectorX[1] - 1.25) * (vectorX[1] - 1.25) - 1;

            // Example. Does NOT have a solution.
            // f1: x1^2 + (x2 - 1.26)^2 - 1
            // funcValue = vectorX[0] * vectorX[0] + (vectorX[1] - 1.26) * (vectorX[1] - 1.26) - 1;
            break;
        }

        case 1:
        {
            // f2: 2*x1-x2+3*x3+2
            funcValue = vectorX[0] * 2 - vectorX[1]+3*vectorX[2]+2;
            break;
        }
        case 2:
        {
            // f3: x1+2*x2-x3-9
            funcValue=vectorX[0]+2*vectorX[1]-vectorX[2]-9;
            break;
        }
        default:
            break;
    }

    return funcValue;
}*/
/*
double getFuncValue_CompileTime(int fIndex, tVectorValue vectorX) {
    double funcValue = 0;
    // Example from task pdf. Has a solution.
    // f1: sin(x)-1=0;
    // f1: 2^x-4=0;
    //funcValue = sin(vectorX[0])-1;
    funcValue=pow(2,vectorX[0])-4;
    return funcValue;
}*/


// TODO - better to avoid using of this "magic number" at all, but use some another criteria.
#define MAX_ITERATIONS 1000

// количество уравнений(а также неизвестных)
int n = 0;

// массив x - начальное приближение
tVectorValue vectorXstart = NULL;

// массив x - решение
tVectorValue vectorXsolution = NULL;

void getVectorXstart()
{   
    // TBD - the start approximation can be input in different ways as well as just hardcoded
    for(int i = 0; i < n; i++)
    {
        //if (!bRuntime)
        {
            vectorXstart[i] = -0.5;
            te_variable vars[] ={"vectorXstart[i]",&vectorXstart[i]};
        }
    }
}

//TBD
double DX_FOR_DERIVATIVE = 1E-5;
double F_ZERO_APPROXIMATION = 1E-5;

void freeNonNull(void ** ptr)
{
    if (ptr != NULL) 
    { 
        free(*ptr);
        ptr = NULL;
    }        
}

void freeMemory()
{    
    freeNonNull((void**)&vectorXstart);
    freeNonNull((void**)&vectorXsolution);
    freeNonNull((void**)&formulae);
}

// fIndex - INPUT - индекс функции, для которой считаем
// vectorX - INPUT  - вектор аргументов X, т.е. в какой точке считаем
// OUTPUT - fi(x1,...,xN), where i = fIndex
double getFuncValue(int fIndex, tVectorValue vectorX)
{
    if (!bRuntime) 
    {
        return getFuncValue_CompileTime(fIndex, vectorX);
    }

    double funcValue = 0;
    
    char* formula = formulae[fIndex];
    
    // TODO - use some libs to get func value by formula string and args
    // funcValue = ...(formula, vectorX);
    
    return funcValue;
}

// vectorX - INPUT
// vectorF - OUTPUT
void getVectorF(tVectorValue vectorX, tVectorValue vectorF)
{
    for (int i = 0; i < n; i++)
    {
        vectorF[i] = getFuncValue(i, vectorX);
    }
}


// Required for Yakobi matrix - считаем частную производную 
// заданной функции в заданной точке по заданному аргументу
//
// fIndex - INPUT - индекс функции, для которой считаем
// vectorX - INPUT  - вектор аргументов X, т.е. в какой точке считаем
// xIndex - INPUT   - по какому аргументу 
double getPartialDerivativeByXindex(int fIndex, tVectorValue vectorX, int xIndex) 
{
    double xToRestore = vectorX[xIndex];

    double deltaX = DX_FOR_DERIVATIVE;
    double funcValueBeforeDeltaX = getFuncValue(fIndex, vectorX);
    
    vectorX[xIndex] += deltaX;
    double deltaF = getFuncValue(fIndex, vectorX) - funcValueBeforeDeltaX;

    // restore x in vectorX
    vectorX[xIndex] = xToRestore;
    
    return deltaF/deltaX;
}

// vectorX - INPUT  - вектор аргументов X, т.е. в какой точке считаем матрицу
// matrix - OUTPUT - matrix to be filled with values
void getMatrixYakobi(double** matrix, tVectorValue vectorX)
{
    for(int fIndex = 0; fIndex < n; fIndex++)    
    {
        for (int xIndex = 0; xIndex < n; xIndex++)
        {
            matrix[fIndex][xIndex] = getPartialDerivativeByXindex(fIndex, vectorX, xIndex);
        }
    }
}

// Транспонированная матрица. За основу взят код из https://pastebin.com/TccZvAEW, адаптирован.
// matrix - INPUT
// dim - INPUT
// matrixTransposed - OUTPUT
void getTransposedMatrix(double ** matrix, int dim, double ** matrixTransposed) 
{
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            matrixTransposed[j][i] = matrix[i][j];
        }
    }
}

// dim - INPUT - dimension
// matrix - INPUT
// matrixReversed - OUTPUT
void getMatrixReversed(double ** matrix, int dim)
{
    if (dim==1)
    {
        matrix[0][0]=((double)1)/matrix[0][0];
    }
    else
    {
    double determinant = get_determinant(matrix, dim);
    
    // матрица алгебраических дополнений
    double** matrixOfCofactors = get_matrix_of_cofactors(matrix, dim); 
    
    // матрица алгебраических дополнений транспонированная
    getTransposedMatrix(matrixOfCofactors, dim, matrix);
    
    free_memory(&matrixOfCofactors, dim);
    
    multiply_by_number(&matrix, dim, dim, ((double)1)/determinant);
}
}

// можно использовать сторонние функции
// matrix - INPUT
// vector - INPUT
// dim - INPUT
// vectorResult - OUTPUT
void multiplyMatrixByVector(double ** matrix, tVectorValue vector, tVectorValue vectorResult, int dim)
{
     for (int i = 0; i < dim; i++)
     {
         double x = 0;

         double * pRow = matrix[i];
         
         for (int j = 0; j < dim; j++)
         {
             x += pRow[j] * vector[j];
         }
         
         vectorResult[i] = x;
         pRow += dim;     
     }     
}

// v1, v2 - INPUT
// dim - INPUT
// as result: v1 has same coordinates as v2
void assignVector(tVectorValue v1, tVectorValue v2, int dim)
{
    memcpy((void*)v1, (void*)v2, dim * sizeof(double));
}

// v1, v2 - INPUT
// vResult - OUTPUT
// dim - INPUT
void substractVectors(tVectorValue v1, tVectorValue v2, tVectorValue vResult, int dim)
{
    for(int i = 0; i < dim; i++)
    {
        vResult[i] = v1[i] - v2[i];
    }
}

// оценка значений {f1, ..., fn} на предмет близости к 0
bool isVectorFnearZero(tVectorValue vectorF)
{
    for(int i = 0; i < n; i++)
    {
        printf("\n isVectorFnearZero vectorF[%d]: %f", i, vectorF[i]);
        
        if (fabs(vectorF[i]) >= F_ZERO_APPROXIMATION)         
        {
             return false;    
        }
    }
    
    return true;
}

// OUTPUT - true if a solution has been found
bool newton()
{    
    bool bRet = false;

    double** matrixYakobiReversed = NULL;
    
    allocate_memory(&matrixYakobiReversed, n, n);    

    // Xk = {x1,...xN} on step k
    tVectorValue vectorXk = (tVectorValue)malloc(n * sizeof(tVectorValue));
    
    // Diff between Xk and Xk+1    
    tVectorValue vectorDeltaXk = (tVectorValue)malloc(n * sizeof(tVectorValue));

    // value of vector Fk = {f1,...fN}(x1,...xN) on step k    
    tVectorValue vectorFk = (tVectorValue)malloc(n * sizeof(tVectorValue));    

    // k - actually may be not needed at all - depends on stopping criteria to be used in case of the system doesn't have a solution...
    int k = 1;

    assignVector(vectorXk, vectorXstart, n);

    do    
    { 
        getVectorF(vectorXk, vectorFk);
        
        // Check if we've found the solution
        if (isVectorFnearZero(vectorFk)) 
        {
            assignVector(vectorXsolution, vectorXk, n);
            bRet = true;
            break;
        }
        else
        {
            printf("\n\n newton iteration k = %d", k);
            
            // operations for formula (3) from pdf     
            double** matrix = matrixYakobiReversed;
            
            getMatrixYakobi(matrix, vectorXk);
            getMatrixReversed(matrix, n);    
            multiplyMatrixByVector(matrix, vectorFk, vectorDeltaXk, n);
            substractVectors(vectorXk, vectorDeltaXk, vectorXk, n);
            
            k++;
            printf("\n ");
        }
    }
    while(k <= MAX_ITERATIONS); //TODO. Condition here is just to stop iterations somewhen in this life... Actually some another criteria is required.
    
    free(vectorXk);
    free(vectorDeltaXk);
    free(vectorFk);
           
    free_memory(&matrixYakobiReversed, n);
    
    return bRet;
}

int main()
{
    bRuntime = false;
    
    if (!bRuntime)
    {
        n = DIM_FOR_COMPILE_TIME_EXAMPLE;
    }
    else
    {
        n = getFormulae();
    }
    
    vectorXsolution = (tVectorValue)malloc(n * sizeof(tVectorValue));
    vectorXstart = (tVectorValue)malloc(n * sizeof(tVectorValue));
    
    getVectorXstart(); 
    
    if (newton())
    {
        printf("\n\n Approximate solution is:");
        for (int i=0; i < n; i++)
        {
            printf("\n x%d  = %f", i, vectorXsolution[i]);
        }
    }
    else
    {
        printf("\n Does not have a solution!");
    }
    
    freeMemory();
}
