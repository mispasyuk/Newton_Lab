#include <stdlib.h>
 
typedef double Type;//Тип значений в матрице(лучше не менять)
 
void allocate_memory(Type *** matrix, int row, int col) {//Выделяет память для матрицы
    (*matrix) = (Type **)malloc(sizeof(Type *) * (row));
    for (int i = 0; i < row; i++) {
        (*matrix)[i] = (Type *)malloc(col * sizeof(Type));
    }
}

// Not in https://pastebin.com/TccZvAEW 
// Added.
void free_memory(Type*** pMatrix, int row)
{
    Type** matrix = (*pMatrix);
    
    for (int i = 0; i < row; i++)
    {
       free(matrix[i]);
    }
    
    free(matrix);
    (*pMatrix) = NULL;    
}
 
Type ** get_cofactor(Type ** matrix, int order, int row, int col) {//Алгебраическое дополнение
    Type ** ans;
    allocate_memory(&ans, order - 1, order - 1);
    for (int i = 0; i < order - 1; ++i) {
        for (int j = 0; j < order - 1; ++j) {
            int x, y;
            x = i + (i >= row);
            y = j + (j >= col);
            ans[i][j] = matrix[x][y];
        }
    }
    return ans;
}
 
Type get_determinant(Type ** matrix, int order) {//Определитель
    if (order == 1) {
        return matrix[0][0];
    }
    Type ans = 0;
    for (int i = 0; i < order; ++i) {
        int pm = ((i & 1) ? -1 : 1);
        ans += (Type)pm * matrix[0][i] * get_determinant(get_cofactor(matrix, order, 0, i), order - 1);
    }
    return ans;
}
 
Type ** get_matrix_of_cofactors(Type ** matrix, int order) {//Матрица алгебраических дополнений
    Type ** ans;
    allocate_memory(&ans, order, order);
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            int pm = (((i + j) & 1) ? -1 : 1);
            ans[i][j] = (Type)pm * get_determinant(get_cofactor(matrix, order, i, j), order - 1);
        }
    }
    return ans;
}
 
void multiply_by_number(Type *** matrix, int row, int col, Type number) {//Умножение матрицы на число
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            (*matrix)[i][j] *= number;
        }
    }
}