#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fillMatrixAndY(double **matrix, double *y, int r, int c, FILE *fp);

void fillIdentity(double **identity, int r);

void gaussJordan(double **mat, double **identity, int r);

void transposeM(double **matrix, double **transpose, int r, int c);

void printMatrix(double **temp, int r, int c);

void vectorProduct(double **squareXtrans, double *y, double *w, int r, int c);

void printVector(double *y, int r);

void squarify(double **matrix, double **squareMat, double **transpose, int r1, int c1, int r2, int c2);

void product(double **squareMat, double **transpose, double **squareXtrans, int r1, int c1, int r2, int c2);

void findPrice(double **data, double *w, double *prices, int r, int c);

int main(int argc, char *argv[])
{

    FILE *fp;
    int r, c;

    // this will most likely be changed to args[1] and args[2], instead of two txt file names
    // WILL THE TRAINING DATA EVER COME AFTER DATA DATA

    char *process = malloc(sizeof(char) * 6);
    fp = fopen(argv[1], "r");

    if (fp == NULL)
    {
        printf("Error opening file\n");
        return EXIT_FAILURE;
    }

    int train = 0;

    fscanf(fp, "%s", &*process);

    if (strcmp(process, "train") == 0)
    {
        train = 1;
    }

    fscanf(fp, "%d", &c);

    fscanf(fp, "%d", &r);

    c = c + 1; // this is because of the first column being all 1's

    // first layer of malloc-ing
    double **transpose = malloc(sizeof(double *) * c);
    double **matrix = malloc(sizeof(double *) * r);
    double **identity = malloc(sizeof(double *) * c);
    double **squareMat = malloc(sizeof(double *) * c);
    double **squareXtrans = malloc(sizeof(double *) * c);
    double *y = malloc(sizeof(double) * r);
    double *w = malloc(sizeof(double) * c);

    for (int i = 0; i < c; i++)
    {
        identity[i] = malloc(c * sizeof(double));
    }

    // based on TA's email
    for (int i = 0; i < c; i++)
    {
        squareMat[i] = malloc(c * sizeof(double));
        squareXtrans[i] = malloc(r * sizeof(double));
        transpose[i] = malloc(r * sizeof(double));
    }
    for (int i = 0; i < r; i++)
    {
        matrix[i] = malloc(c * sizeof(double));
        // identity is malloced in its own method
    }

    // delete all print statements later cause theres gonna be a lot
    if (train)
    {

        fillMatrixAndY(matrix, y, r, c, fp);
        transposeM(matrix, transpose, r, c);
        fillIdentity(identity, c);
        squarify(transpose, squareMat, matrix, c, r, r, c);

        gaussJordan(squareMat, identity, c);
        printf("\n Inverse \n");
        printMatrix(squareMat, c, c);
        product(squareMat, transpose, squareXtrans, c, c, c, r);
        vectorProduct(squareXtrans, y, w, c, r);
    }

    fclose(fp);

    FILE *fp2;

    fp2 = fopen(argv[2], "r");

    if (fp2 == NULL)
    {
        printf("Error opening file\n");
        return EXIT_FAILURE;
    }
    fscanf(fp2, "%s", &*process);

    int dataR = 0;
    int dataC = 0;

    fscanf(fp2, "%d", &dataC);
    fscanf(fp2, "%d", &dataR);
    dataC = dataC + 1;

    double **data = malloc(sizeof(double *) * dataR); // matrix of second file
    double *prices = malloc(sizeof(double) * dataR);

    for (int i = 0; i < dataR; i++)
    {
        {
            data[i] = malloc(dataC * sizeof(double));
            // identity is malloced in its own method
        }
    }

    double content2;
    for (int i = 0; i < dataR; i++)
    {
        data[i][0] = 1; // ASK IF THIS SHOULD BE W[0] BASED ON INSTRUCTIONS
    }
    for (int row = 0; row < dataR; row++)
    {
        for (int col = 1; col < dataC; col++)
        {
            fscanf(fp2, "%lf", &content2);
            data[row][col] = content2;
        }
    }

    findPrice(data, w, prices, dataR, dataC);

    printVector(prices, dataR);

    for (int i = 0; i < r; i++)
    {
        double *del = matrix[i];

        free(del);
    }
    for (int i = 0; i < c; i++)
    {
        free(squareMat[i]);
        free(transpose[i]);
        free(identity[i]);
        free(squareXtrans[i]);
    }
    for (int i = 0; i < dataR; i++)
    {
        free(data[i]);
    }

    fclose(fp2);

    free(squareMat);
    free(data);
    free(matrix);
    free(transpose);
    free(identity);
    free(squareXtrans);

    free(y);
    free(w);
    free(prices);

    free(process);
    return EXIT_SUCCESS;
}

// DO NOT CHANGE
void fillMatrixAndY(double **matrix, double *y, int r, int c, FILE *fp)
{

    double content;
    int i = 0;
    for (int row = 0; row < r; row++)
    {
        matrix[row][0] = 1;
    }
    for (int row = 0; row < r; row++)
    {
        for (int col = 1; col < c; col++)
        {
            fscanf(fp, "%lf", &content);
            matrix[row][col] = content;
        }
        fscanf(fp, "%lf", &content);
        y[i] = content;
        i++;
    }
}

// dont worry about this
void printVector(double *temp, int s)
{
    for (int i = 0; i < s; i++)
    {
        printf("%.0f\n", temp[i]);
        // printf("%lf ", temp[i]); does not round, one above does
    }
}

// DO NOT TOUCH
void fillIdentity(double **identity, int c)
{
    // if row = col add 1, else add 0
    for (int row = 0; row < c; row++)
    {
        for (int col = 0; col < c; col++)
        {
            if (row == col)
            {
                identity[row][col] = 1;
            }
            else
            {
                identity[row][col] = 0;
            }
        }
    }
}

// DO NOT TOUCH
void transposeM(double **matrix, double **transpose, int r, int c)
{
    for (int i = 0; i < c; i++)
    {
        for (int j = 0; j < r; j++)
        {
            transpose[i][j] = matrix[j][i];
        }
    }
}

// DO NOT TOUCH
void squarify(double **transpose, double **squareMat, double **matrix, int r1, int c1, int r2, int c2)
{
    for (int row = 0; row < r1; row++)
    {
        for (int col = 0; col < c2; col++)
        {
            double sum = 0;
            for (int k = 0; k < r2; k++)
            {
                sum += (transpose[row][k] * matrix[k][col]);
            }
            squareMat[row][col] = sum;
            sum = 0;
        }
    }
}

// this is the (xTx)^-1 * xT step
// this MIGHT be the problem child
void product(double **squareMat, double **transpose, double **squareXtrans, int r1, int c1, int r2, int c2)
{

    for (int row = 0; row < r1; row++)
    {
        for (int col = 0; col < c2; col++)
        {
            double sum = 0;
            for (int k = 0; k < r2; k++)
            {
                sum += (squareMat[row][k] * transpose[k][col]);
            }
            squareXtrans[row][col] = sum;
        }
    }
}

// DO NOT TOUCH
void gaussJordan(double **squareMat, double **identity, int r)
{

    double factor;
    for (int p = 0; p < r; p++)
    {
        factor = squareMat[p][p];
        for (int c = 0; c < r; c++)
        {
            squareMat[p][c] /= factor;
            identity[p][c] /= factor;
        }
        for (int i = p + 1; i < r; i++) // re add r-1 if this breaks
        {
            factor = squareMat[i][p];
            for (int c = 0; c < r; c++)
            {
                squareMat[i][c] -= squareMat[p][c] * factor;
                identity[i][c] -= identity[p][c] * factor;
            }
        }
    }

    for (int p = r - 1; p >= 0; p--)
    {
        for (int i = p - 1; i >= 0; i--)
        {
            factor = squareMat[i][p];
            for (int c = 0; c < r; c++)
            {
                squareMat[i][c] -= squareMat[p][c] * factor;
                identity[i][c] -= identity[p][c] * factor;
            }
        }
    }

    for (int row = 0; row < r; row++)
    {
        for (int col = 0; col < r; col++)
        {
            squareMat[row][col] = identity[row][col];
        }
    }
}

// DO NOT TOUCH
void printMatrix(double **temp, int r, int c)
{
    for (int row = 0; row < r; row++)
    {
        for (int col = 0; col < c; col++)
        {
            printf("%lf ", temp[row][col]);
        }
        printf("\n");
    }
}

// gets w vector
void vectorProduct(double **squareXtrans, double *y, double *w, int r, int c)
{
    int wi = 0;
    double sum = 0;
    for (int row = 0; row < r; row++)
    {
        for (int col = 0; col < c; col++)
        {
            sum += y[col] * squareXtrans[row][col];
        }
        w[wi] = sum;
        wi++;
        sum = 0;
    }
}

void findPrice(double **data, double *w, double *prices, int r, int c)
{
    int pi = 0;
    double sum = 0;
    for (int row = 0; row < r; row++)
    {
        for (int col = 0; col < c; col++)
        {
            sum += w[col] * data[row][col];
        }
        prices[pi] = sum;
        pi++;
        sum = 0;
    }
}