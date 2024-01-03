/*
@Author: Andy
@Date: 2024/1/1
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#pragma ones

typedef double elem_type;

typedef struct __Matrix{
    int rows;
    int columns;
    elem_type **data;
} Matrix;

Matrix *mat_init(int rows, int columns);
Matrix *mat_init_like(Matrix *M);
Matrix *mat_create( elem_type *data, int rows, int columns);
int mat_free(Matrix *M);
int mat_print(Matrix *M);

Matrix *mat_zeros(int rows, int columns);
Matrix *mat_zeros_like(Matrix *M);
Matrix *mat_ones(int rows, int columns);
Matrix *mat_ones_like(Matrix *M);
Matrix *mat_I(int order);
Matrix *mat_I_like(Matrix *M);
Matrix *mat_copy(Matrix *M);
Matrix *mat_row_stake(Matrix *M1, Matrix *M2);
Matrix *mat_column_stake(Matrix *M1, Matrix *M2);
Matrix *mat_concatenate(Matrix *M1, Matrix *M2, int axis);

Matrix *mat_get_row(Matrix *M, int row);
Matrix *mat_get_column(Matrix *M, int column);
Matrix *mat_get_elem(Matrix *M, int row, int column);
Matrix *mat_reshape(Matrix *M, int rows, int columns);
int mat_reshape_dst(Matrix *M, int rows, int columns);
Matrix *mat_abs(Matrix *M);
int mat_abs_dst(Matrix *M);
Matrix *mat_min(Matrix *M, int axis);
Matrix *mat_max(Matrix *M, int axis);
Matrix *mat_sum(Matrix *M, int axis);

Matrix *mat_add(Matrix *M1, Matrix *M2);
int mat_add_dst(Matrix *M_Dst, Matrix *M_Src);
Matrix *mat_neg(Matrix *M);
int mat_neg_dst(Matrix *M_Dst);
Matrix *mat_sub(Matrix *M1, Matrix *M2);
int mat_sub_dst(Matrix *M_Dst, Matrix *M_Src);
Matrix *mat_numul(Matrix *M, elem_type k);
int mat_numul_dst(Matrix *M_Dst, elem_type k);
Matrix *mat_mul(Matrix *M1, Matrix *M2);
Matrix *mat_pow(Matrix *M, int n);

Matrix *mat_T(Matrix *M);
Matrix *mat_tr(Matrix *M);
Matrix *mat_diag(Matrix *M);


Matrix *mat_init(int rows, int columns)
{

    Matrix *M=NULL;
    int k=0;
    M=(Matrix *)malloc(sizeof(Matrix));
    M->rows=rows;
    M->columns=columns;
    M->data=(elem_type**)malloc(sizeof(elem_type*)*(rows+1));
    for (k=0;k<rows;k++)
        M->data[k]=(elem_type*)malloc(sizeof(elem_type)*(columns+1));
    return M;
}

Matrix *mat_init_like(Matrix *M)
{

    return mat_init(M->rows,M->columns);
}

Matrix *mat_create(elem_type *data, int rows, int columns)
{

    Matrix *M=mat_init(rows, columns);
    int i,j;
    for(i=0;i<rows;i++)
        for(j=0;j<columns;j++)
            M->data[i][j]=data[i*columns+j];
    return M;
}

int mat_free(Matrix *M)
{
    if (M == NULL){
        printf("WARNING:Empty Matrix! from function %s",__func__);
        return -1;
    }
    int k=0;
    for(k=0;k<M->rows;k++)
        free(M->data[k]);
    free(M->data);
    free(M);
    printf("success free\n");
    return 0;
}

int mat_print(Matrix *M)
{
    int i=0,j=0;
    if(M->rows>1 || M->columns>1)
		printf("rows=%d, columns=%d\n",M->rows,M->columns);
    for(i=0;i<M->rows;i++){
		printf(i==0?"[[":"\n[");
        for(j=0;j<M->columns;j++)
            if(fabs(M->data[i][j]-(int)M->data[i][j])<1e-8)
                printf(" %4ld ",(int)M->data[i][j]);
            else
                printf(" %4.2f ",M->data[i][j]);
        printf("]");
    }
    printf("]\n");
    return 0;
}

Matrix *mat_zeros(int rows, int columns)
{

    Matrix *M=mat_init(rows,columns);
    int i,j;
    for(i=0;i<rows;i++)
        for(j=0;j<columns;j++)
            M->data[i][j]=0;
    return M;
}

Matrix *mat_zeros_like(Matrix *M)
{
    return mat_zeros(M->rows,M->columns);
}

Matrix *mat_ones(int rows, int columns)
{

    Matrix *M=mat_init(rows,columns);
    int i,j;
    for(i=0;i<rows;i++)
        for(j=0;j<columns;j++)
            M->data[i][j]=1;
    return M;
}

Matrix *mat_ones_like(Matrix *M)
{
    return mat_ones(M->rows,M->columns);
}

Matrix *mat_I(int order)
{
    
    Matrix *M=mat_zeros(order,order);
    int k;
    for(k=0;k<order;k++)
        M->data[k][k]=1;
    return M;
}

Matrix *mat_I_like(Matrix *M)
{

    return mat_I(M->rows);
}

Matrix *mat_copy(Matrix *M)
{

    Matrix *M_res=mat_init_like(M);
    int i,j;
    for(i=0;i<M->rows;i++)
        for(j=0;j<M->columns;j++)
            M_res->data[i][j]=M->data[i][j];
    return M_res;
}

Matrix *mat_row_stake(Matrix *M1, Matrix *M2)
{

    Matrix *M_res=mat_init(M1->rows+M2->rows,M2->columns);
    int i,j;
    for(i=0;i<M_res->rows;i++)
        for(j=0;j<M_res->columns;j++)
            M_res->data[i][j]=
            i<M1->rows?M1->data[i][j]:M2->data[i-M1->rows][j];
    return M_res;
}

Matrix *mat_column_stake(Matrix *M1, Matrix *M2)
{

    Matrix *M_res=mat_init(M1->rows,M1->columns+M2->columns);
    int i,j;
    for(i=0;i<M_res->rows;i++)
        for(j=0;j<M_res->columns;j++)
            M_res->data[i][j]=
            j<M1->columns?M1->data[i][j]:M2->data[i][j-M1->columns];
    return M_res;
            
}

Matrix *mat_concatenate(Matrix *M1, Matrix *M2, int axis)
{

    Matrix *M_res=NULL;
    if(axis==0)
        M_res=mat_row_stake(M1,M2);
    else if(axis==1)
        M_res=mat_column_stake(M1,M2);
    return M_res;
}

Matrix *mat_get_row(Matrix *M, int row)
{

    Matrix *M_res=mat_init(1,M->columns);
    int j;
    row+=row>=0?0:M->rows;
    for(j=0;j<M->columns;j++)
        M_res->data[0][j]=M->data[row][j];
    return M_res;
}

Matrix *mat_get_column(Matrix *M, int column)
{

    Matrix *M_res=mat_init(M->rows,1);
    int i;
    column+=column>=0?0:M->columns;
    for(i=0;i<M->rows;i++)
        M_res->data[i][0]=M->data[i][column];
    return M_res;
}

Matrix *mat_get_elem(Matrix *M, int row, int column)
{

    row+=row>=0?0:M->rows;
    column+=column>=0?0:M->columns;
    Matrix *M_res=mat_init(1,1);
    M_res->data[0][0]=M->data[row][column];
    return M_res;
}

Matrix *mat_reshape(Matrix *M, int rows, int columns)
{

    Matrix *M_res=mat_init(rows,columns);
    int i,j,k;
    for(i=0,k=0;i<rows;i++)
        for(j=0;j<columns;j++,k++)
            M_res->data[i][j]=M->data[k/M->columns][k%M->columns];
    return M_res;

}

int mat_reshape_dst(Matrix *M, int rows, int columns)
{

    Matrix *M_res=mat_reshape(M, rows, columns);
    mat_free(M);
    M=M_res;
    return 0;
}

Matrix *mat_abs(Matrix *M)
{

    Matrix *M_res=mat_copy(M);
    mat_abs_dst(M_res);
    return M_res;
}

int mat_abs_dst(Matrix *M)
{

    int i,j;
    for(i=0;i<M->rows;i++)
        for(j=0;j<M->columns;j++)
            M->data[i][j]=M->data[i][j]>=0?M->data[i][j]:-M->data[i][j];
    return 0;
}

Matrix *mat_min(Matrix *M, int axis)
{

    Matrix *M_res=NULL;
    int i,j;
    if(axis==0){
        M_res=mat_get_row(M,0);
        for(i=0;i<M->rows;i++)
            for(j=0;j<M->columns;j++)
                if(M_res->data[i][0]>M->data[i][j])
                    M_res->data[i][0]=M->data[i][j];
    }
    else if(axis==1){
        M_res=mat_get_column(M,0);
        for(i=0;i<M->rows;i++)
            for(j=0;j<M->columns;j++)
                if(M_res->data[0][j]>M->data[i][j])
                    M_res->data[0][j]=M->data[i][j];
    }
    return M_res;
}

Matrix *mat_max(Matrix *M, int axis)
{

    Matrix *M_res=NULL;
    int i,j;
    if(axis==0){
        M_res=mat_get_row(M,0);
        for(i=0;i<M->rows;i++)
            for(j=0;j<M->columns;j++)
                if(M_res->data[0][j]<M->data[i][j])
                    M_res->data[0][j]=M->data[i][j];
    }
    else if(axis==1){
        M_res=mat_get_column(M,0);
        for(i=0;i<M->rows;i++)
            for(j=0;j<M->columns;j++)
                if(M_res->data[i][0]<M->data[i][j])
                    M_res->data[i][0]=M->data[i][j];
    }
    return M_res;    
}

Matrix *mat_sum(Matrix *M, int axis)
{

    Matrix *M_res=NULL;
    int i,j;
    if(axis==0){
        M_res=mat_zeros(1,M->columns);
        for(i=0;i<M->rows;i++)
            for(j=0;j<M->columns;j++)
                M_res->data[0][j]+=M->data[i][j];
    }
    else if(axis==1){
        M_res=mat_zeros(M->rows,1);
        for(i=0;i<M->rows;i++)
            for(j=0;j<M->columns;j++)
                M_res->data[i][0]+=M->data[i][j];
    }
    return M_res;    
}

Matrix *mat_add(Matrix *M1, Matrix *M2)
{
    
    Matrix *M_res=mat_init_like(M1);
    int i,j;
    for(i=0;i<M1->rows;i++)
        for(j=0;j<M1->columns;j++)
            M_res->data[i][j]=M1->data[i][j]+M2->data[i][j];
    return M_res;
}

int mat_add_dst(Matrix *M_Dst, Matrix *M_Src)
{

    int i,j;
    for(i=0;i<M_Dst->rows;i++)
        for(j=0;j<M_Dst->columns;j++)
            M_Dst->data[i][j]+=M_Src->data[i][j];
    return 0;
}

Matrix *mat_neg(Matrix *M)
{

    Matrix *M_res=mat_init_like(M);
    int i,j;
    for(i=0;i<M->rows;i++)
        for(j=0;j<M->columns;j++)
            M_res->data[i][j]=-M->data[i][j];
    return M_res;
}

int mat_neg_dst(Matrix *M)
{

    int i,j;
    for(i=0;i<M->rows;i++)
        for(j=0;j<M->columns;j++)
            M->data[i][j]=-M->data[i][j];
    return 0;
}

Matrix *mat_sub(Matrix *M1, Matrix *M2)
{

    Matrix *M_res=mat_init_like(M1);
    int i,j;
    for(i=0;i<M1->rows;i++)
        for(j=0;j<M1->columns;j++)
            M_res->data[i][j]=M1->data[i][j]-M2->data[i][j];
    return M_res;
}

int mat_sub_dst(Matrix *M_Dst, Matrix *M_Src)
{

    int i,j;
    for(i=0;i<M_Dst->rows;i++)
        for(j=0;j<M_Dst->columns;j++)
            M_Dst->data[i][j]-=M_Src->data[i][j];
    return 0;
}

Matrix *mat_numul(Matrix *M, elem_type k)
{

    Matrix *M_res=mat_init_like(M);
    int i,j;
    for(i=0;i<M->rows;i++)
        for(j=0;j<M->columns;j++)
            M_res->data[i][j]=k*M->data[i][j];
    return M_res;
}

int mat_numul_dst(Matrix *M_Dst, elem_type k)
{

    int i,j;
    for(i=0;i<M_Dst->rows;i++)
        for(j=0;j<M_Dst->columns;j++)
            M_Dst->data[i][j]*=k;
    return 0;
}

Matrix *mat_mul(Matrix *M1, Matrix *M2)
{

	int i,j,k;
	Matrix *M_res=mat_zeros(M1->rows,M2->columns);
	for(i=0;i<M1->rows;i++)
		for(j=0;j<M2->columns;j++)
			for(k=0;k<M1->columns;k++)
				M_res->data[i][j]+=M1->data[i][k]*M2->data[k][j];
	return M_res;  
}

Matrix *mat_pow(Matrix *M, int n)
{

    int k=0;
	Matrix *M_tmp=mat_copy(M);
	Matrix *M_res=NULL;
	for(k=1;k<n;k++){
		M_res=mat_mul(M_tmp,M);
		mat_free(M_tmp);
		M_tmp=M_res;
	}
	return M_res;
}

Matrix *mat_T(Matrix *M)
{

    Matrix *M_res=mat_init(M->columns,M->rows);
    int i,j;
    for(i=0;i<M_res->rows;i++)
        for(j=0;j<M_res->columns;j++)
            M_res->data[i][j]=M->data[j][i];
    return M_res;
}

Matrix *mat_tr(Matrix *M)
{

    Matrix *M_res=mat_zeros(1,1);
    int k;
    for(k=0;k<M->rows;k++)
        M_res->data[0][0]+=M->data[k][k];
    return M_res;
}

Matrix *mat_diag(Matrix *M)
{


    Matrix *M_res=NULL;
    int k;
    if(M->columns==1){
        M_res=mat_zeros(M->rows,M->rows);
        for(k=0;k<M_res->rows;k++)
            M_res->data[k][k]=M->data[k][0];
    }
    else if(M->rows==1){
        M_res=mat_zeros(M->columns,M->columns);
        for(k=0;k<M_res->columns;k++)
            M_res->data[k][k]=M->data[0][k];
    }
    else if(M->rows==M->columns){
        M_res=mat_zeros(1,M->columns);
        for(k=0;k<M->columns;k++)
            M_res->data[0][k]=M->data[k][k];
    }
    return M_res;
}


