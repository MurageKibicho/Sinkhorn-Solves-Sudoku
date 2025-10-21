#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#define STR_INDEX(x, y, cols) ((x) * (cols) + (y))
//clear && gcc ParseCSV.c -lm -o m.o && ./m.o
typedef struct sudoku_grid_struct *Sudoku;
struct sudoku_grid_struct
{
	int **question;
	int **sinkhornTrialAnswer;
	int **answer;
	double **cellProbabilities;
	double **newProbabilities;
	double **constraintMatrices;
	double **constraintBalanced;
	int *updateCount;
	int rating;
	int gridSquareDimension;//9*9 classic sudoku grid
	int symbolCount;//9 symbols[1,2,3,4,5,6,7,8,9]
	int constraintMatrixCount;//27 possible constraints
	
	//Accuracy metrics
	double mse;
	double accuracy;
	int correctPredictions;
	int totalPredictions;
};

Sudoku CreateSudokuGrid()
{
	Sudoku sudoku = malloc(sizeof(struct sudoku_grid_struct));
	sudoku->rating = 0;
	sudoku->symbolCount = 9;
	sudoku->gridSquareDimension = 9;
	sudoku->constraintMatrixCount = 27;
	sudoku->question = malloc(sudoku->gridSquareDimension * sizeof(int *));
	sudoku->answer   = malloc(sudoku->gridSquareDimension * sizeof(int *));
	sudoku->sinkhornTrialAnswer   = malloc(sudoku->gridSquareDimension * sizeof(int *));
	
	//81 cellProbabilities
	sudoku->updateCount   = calloc(sudoku->gridSquareDimension * sudoku->gridSquareDimension, sizeof(int));
	sudoku->cellProbabilities   = malloc(sudoku->gridSquareDimension * sudoku->gridSquareDimension * sizeof(double *));
	sudoku->newProbabilities   = malloc(sudoku->gridSquareDimension * sudoku->gridSquareDimension * sizeof(double *));
	for(int i = 0; i < sudoku->gridSquareDimension; i++)
	{
		sudoku->question[i] = calloc(sudoku->gridSquareDimension, sizeof(int));
		sudoku->answer[i] = calloc(sudoku->gridSquareDimension, sizeof(int));
		sudoku->sinkhornTrialAnswer[i] = calloc(sudoku->gridSquareDimension, sizeof(int));
	}
	for(int i = 0; i < sudoku->gridSquareDimension * sudoku->gridSquareDimension; i++)
	{
		sudoku->cellProbabilities[i] = calloc(sudoku->symbolCount, sizeof(double));
		sudoku->newProbabilities[i] = calloc(sudoku->symbolCount, sizeof(double));
	}
	//27 constraint matrices
	sudoku->constraintMatrices = malloc(sudoku->constraintMatrixCount * sizeof(double*));
	sudoku->constraintBalanced = malloc(sudoku->constraintMatrixCount * sizeof(double*));
	for(int i = 0; i < sudoku->constraintMatrixCount; i++)
	{
		sudoku->constraintMatrices[i] = calloc(sudoku->gridSquareDimension * sudoku->gridSquareDimension, sizeof(double));
		sudoku->constraintBalanced[i] = calloc(sudoku->gridSquareDimension * sudoku->gridSquareDimension, sizeof(double));	
	}
	return sudoku;
}

void DestroySudoku(Sudoku sudoku)
{
	if(sudoku)
	{
		for(int i = 0; i < sudoku->constraintMatrixCount; i++)
		{
			free(sudoku->constraintBalanced[i]);
			free(sudoku->constraintMatrices[i]);
		}
		free(sudoku->constraintBalanced);
		free(sudoku->constraintMatrices);
		for(int i = 0; i < sudoku->gridSquareDimension * sudoku->gridSquareDimension; i++)
		{
			free(sudoku->newProbabilities[i]);
			free(sudoku->cellProbabilities[i]);
		}
		for(int i = 0; i < sudoku->gridSquareDimension; i++)
		{
			free(sudoku->question[i]);
			free(sudoku->sinkhornTrialAnswer[i]);
			free(sudoku->answer[i]);
		}
		free(sudoku->newProbabilities);
		free(sudoku->cellProbabilities);
		free(sudoku->question);
		free(sudoku->sinkhornTrialAnswer);
		free(sudoku->answer);
		free(sudoku->updateCount);
		free(sudoku);
	}
}

void PrintGrid(Sudoku sudoku, int printQuestion)
{
	if(printQuestion == 0)
	{
		printf("Question rating(%d)\n", sudoku->rating);
		for(int i = 0; i < sudoku->gridSquareDimension; i++)
		{
			for(int j = 0; j < sudoku->gridSquareDimension; j++)
			{
				printf("%d,", sudoku->question[i][j]);
			}
			printf("\n");
		}
	}
	else if(printQuestion == 1)
	{
		printf("Answer\n");
		for(int i = 0; i < sudoku->gridSquareDimension; i++)
		{
			for(int j = 0; j < sudoku->gridSquareDimension; j++)
			{
				printf("%d,", sudoku->answer[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	else if(printQuestion == 2)
	{
		printf("Trial\n");
		for(int i = 0; i < sudoku->gridSquareDimension; i++)
		{
			for(int j = 0; j < sudoku->gridSquareDimension; j++)
			{
				printf("%d,", sudoku->sinkhornTrialAnswer[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n");
}

void UpdateSudokuGrid(Sudoku sudoku, char *question, char *answer, char *rating)
{
	sudoku->rating = atoi(rating);
	for(int i = 0; i < sudoku->gridSquareDimension; i++)
	{
		for(int j = 0; j < sudoku->gridSquareDimension; j++)
		{	
			int index = STR_INDEX(i,j, sudoku->gridSquareDimension);
			
			if(question[index] != '.')
			{
				sudoku->question[i][j] = question[index] - '0';	
			}
 			if(answer[index] != '.')
			{
				sudoku->answer[i][j] = answer[index] - '0';	
			}
			if(sudoku->question[i][j] != sudoku->answer[i][j])
			{
				sudoku->question[i][j] = 0;
			}
		}
	}
}

void FindCellProbabilities(Sudoku sudoku)
{
	for(int i = 0; i < sudoku->gridSquareDimension; i++)
	{
		for(int j = 0; j < sudoku->gridSquareDimension; j++)
		{	
			int index = STR_INDEX(i,j, sudoku->gridSquareDimension);
			int symbol = sudoku->question[i][j];
			if(symbol != 0)
			{
				//Set all probabilities to zero and set the clue to 1
				for(int k = 0; k < sudoku->symbolCount; k++){sudoku->cellProbabilities[index][k] = 0;}
				sudoku->cellProbabilities[index][symbol-1] = 1.0;
			}
			else
			{
				//Set all probabilities to one and set the clues to 0
				for(int k = 0; k < sudoku->symbolCount; k++){sudoku->cellProbabilities[index][k] = 1;}
				//Loop vertically and see what clues given
				for(int k = 0; k < sudoku->gridSquareDimension; k++)
				{
					symbol = sudoku->question[k][j];
					if(symbol != 0){sudoku->cellProbabilities[index][symbol-1] = 0.0;}
				}
				//Loop horizontally and see what clues given
				for(int k = 0; k < sudoku->gridSquareDimension; k++)
				{
					symbol = sudoku->question[i][k];
					if(symbol != 0){sudoku->cellProbabilities[index][symbol-1] = 0.0;}
				}
				//Loop through 3*3 tile and see what clues given
				int boxSize = 3;  
				int boxStartRow = (i / boxSize) * boxSize;
				int boxStartCol = (j / boxSize) * boxSize;
				for(int row = 0; row < boxSize; row++)
				{
					for(int col = 0; col < boxSize; col++)
					{
						symbol = sudoku->question[boxStartRow + row][boxStartCol + col];
						if(symbol != 0){sudoku->cellProbabilities[index][symbol - 1] = 0.0;}
						
					}
				}
				//Count number of ones
				double oneCount = 0;
				for(int k = 0; k < sudoku->symbolCount; k++){if(sudoku->cellProbabilities[index][k] == 1){oneCount += 1;}}
				//Normalize by one count
				for(int k = 0; k < sudoku->symbolCount; k++){sudoku->cellProbabilities[index][k] /= oneCount;}
			}
			
		}
	}
}

void UpdateNewCellProbabilities(Sudoku sudoku)
{
	int dim = sudoku->gridSquareDimension;
	int symbolCount = sudoku->symbolCount;
	int boxSize = 3;

	//Reset update count
	for(int i = 0; i < dim * dim; i++){sudoku->updateCount[i] = 0;}
	//Reset new cell probabilities
	for(int i = 0; i < sudoku->gridSquareDimension * sudoku->gridSquareDimension; i++)
	{
		for(int j = 0; j < sudoku->symbolCount; j++)
		{
			sudoku->newProbabilities[i][j] = 0.0f;
		}
	}
	//Sum all probabilities from BALANCED constraints
	//Row constraints(0–8)
	for(int row = 0; row < dim; row++)
	{
		int constraintIndex = row;
		for(int col = 0; col < dim; col++)
		{
			int globalCellIndex = STR_INDEX(row, col, dim);
			int localCellIndex = col;
			for(int k = 0; k < symbolCount; k++)
			{
				sudoku->newProbabilities[globalCellIndex][k] += sudoku->constraintBalanced[constraintIndex][localCellIndex * symbolCount + k];
			}
			sudoku->updateCount[globalCellIndex]++;
		}
	}

	//Column constraints(9–17)
	for(int col = 0; col < dim; col++)
	{
		int constraintIndex = 9 + col;
		for(int row = 0; row < dim; row++)
		{
			int globalCellIndex = STR_INDEX(row, col, dim);
			int localCellIndex = row;
			for (int k = 0; k < symbolCount; k++)
			{
				sudoku->newProbabilities[globalCellIndex][k] += sudoku->constraintBalanced[constraintIndex][localCellIndex * symbolCount + k];
			}
			sudoku->updateCount[globalCellIndex]++;
		}
	}

	//Tile constraints (18–26)
	for(int box = 0; box < dim; box++)
	{
		int constraintIndex = 18 + box;
		int startRow = (box / boxSize) * boxSize;
		int startCol = (box % boxSize) * boxSize;

		int localCellIndex = 0;
		for(int r = 0; r < boxSize; r++)
		{
			for(int c = 0; c < boxSize; c++)
			{
				int row = startRow + r;
				int col = startCol + c;
				int globalCellIndex = STR_INDEX(row, col, dim);
				for (int k = 0; k < symbolCount; k++)
				{
					sudoku->newProbabilities[globalCellIndex][k] += sudoku->constraintBalanced[constraintIndex][localCellIndex * symbolCount + k];
				}
				localCellIndex++;
				sudoku->updateCount[globalCellIndex]++;
			}
		}
	}
	
	//Average the probabilities from all constraints and update
	for(int i = 0; i < sudoku->gridSquareDimension * sudoku->gridSquareDimension; i++)
	{
		int row = i / dim;
		int col = i % dim;
		if(sudoku->question[row][col] == 0 && sudoku->updateCount[i] > 0)
		{
			for(int j = 0; j < sudoku->symbolCount; j++)
			{
				sudoku->newProbabilities[i][j] /= (double) sudoku->updateCount[i];
				//Also update cellProbabilities
				sudoku->cellProbabilities[i][j] = sudoku->newProbabilities[i][j];
			}
		}	
	}
}
void PrintCellProbabilities(Sudoku sudoku)
{
	for(int i = 0; i < sudoku->gridSquareDimension; i++)
	{
		for(int j = 0; j < sudoku->gridSquareDimension; j++)
		{
			int index = STR_INDEX(i, j, sudoku->gridSquareDimension);
			printf("Cell (%d, %d): ", i, j);
			for(int k = 0; k < sudoku->symbolCount; k++)
			{
				printf("%.1f ", sudoku->cellProbabilities[index][k]);
			}
			printf("\n");
		}
	}
}
void FindConstraintMatrices(Sudoku sudoku)
{
	int dim = sudoku->gridSquareDimension;
	int symbolCount = sudoku->symbolCount;
	int boxSize = 3;

	//Row constraints(0–8)
	for(int row = 0; row < dim; row++)
	{
		int constraintIndex = row;
		for(int col = 0; col < dim; col++)
		{
			int globalCellIndex = STR_INDEX(row, col, dim);
			int localCellIndex = col;
			for(int k = 0; k < symbolCount; k++)
			{
				sudoku->constraintMatrices[constraintIndex][localCellIndex * symbolCount + k] = sudoku->cellProbabilities[globalCellIndex][k];
			}
		}
	}

	//Column constraints(9–17)
	for(int col = 0; col < dim; col++)
	{
		int constraintIndex = 9 + col;
		for(int row = 0; row < dim; row++)
		{
			int globalCellIndex = STR_INDEX(row, col, dim);
			int localCellIndex = row;
			for (int k = 0; k < symbolCount; k++)
			{
				sudoku->constraintMatrices[constraintIndex][localCellIndex * symbolCount + k] = sudoku->cellProbabilities[globalCellIndex][k];
			}
		}
	}

	//Tile constraints (18–26)
	for(int box = 0; box < dim; box++)
	{
		int constraintIndex = 18 + box;
		int startRow = (box / boxSize) * boxSize;
		int startCol = (box % boxSize) * boxSize;

		int localCellIndex = 0;
		for(int r = 0; r < boxSize; r++)
		{
			for(int c = 0; c < boxSize; c++)
			{
				int row = startRow + r;
				int col = startCol + c;
				int globalCellIndex = STR_INDEX(row, col, dim);
				for (int k = 0; k < symbolCount; k++)
				{
					sudoku->constraintMatrices[constraintIndex][localCellIndex * symbolCount + k] = sudoku->cellProbabilities[globalCellIndex][k];
				}
				localCellIndex++;
			}
		}
	}
}

void PrintSingleConstraintMatrix(double *constraintMatrix)
{
	for(int i = 0; i < 9; i++)
	{	
		for(int j = 0; j < 9; j++)
		{
			printf("%.1f,", constraintMatrix[STR_INDEX(i, j, 9)]);
		}
		printf("\n");
	}
}

double SinkhornKnoppAlgorithm(int rows, int cols, int sinkhornIterations, double sinkhornTolerance, double *matrix, double *matrixCopy, double *rowScaling, double *colScaling)
{
	//We assume you already checked for total support 
	assert(rows == cols);
	memcpy(matrixCopy, matrix, rows * cols * sizeof(double));
	
	//Initialize Scaling Vectors to 1
	for(int i = 0; i < rows; i++)
	{
		rowScaling[i] = 1.0;
		colScaling[i] = 1.0;
	}
	int n = rows;
	double maxError = 0.0;
	for(int iteration = 0; iteration < sinkhornIterations; iteration++)
	{
		maxError = 0.0;
		//Update Row Scaling
		for(int i = 0; i < n; i++)
		{
			double rowSum = 0.0;
			for(int j = 0; j < n; j++)
			{
				rowSum += matrix[i * n + j] * colScaling[j];
			}
			if(rowSum > 1e-12)
			{
				rowScaling[i] = 1.0 / rowSum;
			}
			maxError = fmax(maxError, fabs(rowSum - 1.0));
		}
		//Update column scaling
		for(int j = 0; j < n; j++)
		{
			double colSum = 0.0;
			for(int i = 0; i < n; i++)
			{
				colSum += rowScaling[i] * matrix[i * n + j];
			}
			if(colSum > 1e-12)
			{
				colScaling[j] = 1.0 / colSum;
			}
			maxError = fmax(maxError, fabs(colSum - 1.0));
		}
		
		if(maxError < sinkhornTolerance){break;}
	}	
	//Apply scaling
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			matrixCopy[i * n + j] = rowScaling[i] * matrix[i * n + j] * colScaling[j];
		}
	}
	return maxError;	
}

void ExtractSinkhornAnswer(Sudoku sudoku)
{
	sudoku->mse = 0.0f;
	sudoku->accuracy = 0.0f;
	sudoku->correctPredictions = 0;
	sudoku->totalPredictions = 0;
	int totalCells = sudoku->gridSquareDimension * sudoku->gridSquareDimension;
	int cellsCompared = 0;
	for(int i = 0; i < sudoku->gridSquareDimension; i++)
	{
		for(int j = 0; j < sudoku->gridSquareDimension; j++)
		{
			int index = STR_INDEX(i, j, sudoku->gridSquareDimension);
			//Keep original clue cell value
			if(sudoku->question[i][j] != 0){sudoku->sinkhornTrialAnswer[i][j] = sudoku->question[i][j];}
			else
			{
				//Find symbol with highest probability
				double maxProbabibility = -1.0;
				int bestSymbol = 0;
				for(int k = 0; k < sudoku->symbolCount; k++)
				{
					if(sudoku->cellProbabilities[index][k] > maxProbabibility)
					{
						maxProbabibility = sudoku->cellProbabilities[index][k];
						bestSymbol = k + 1;
					}
				}
				sudoku->sinkhornTrialAnswer[i][j] = bestSymbol;
				if(sudoku->sinkhornTrialAnswer[i][j] == sudoku->answer[i][j])
				{
					sudoku->correctPredictions += 1;
				}
				double diff = sudoku->sinkhornTrialAnswer[i][j] - sudoku->answer[i][j];
				sudoku->mse += diff * diff;
				sudoku->totalPredictions += 1;
			}	
		}
	}
	if(sudoku->totalPredictions > 0)
	{
		sudoku->accuracy = (double) sudoku->correctPredictions / sudoku->totalPredictions * 100.0f;
		sudoku->mse /= sudoku->totalPredictions;
	}
}

void PrintMetrics(int index, Sudoku sudoku)
{
	printf("Sudoku Index: %3d\n", index);
	printf("Correct predictions: %3d\n", sudoku->correctPredictions);
	printf("Total   predictions: %3d\n", sudoku->totalPredictions);
	printf("Prediction accuracy: %.3f\n", sudoku->accuracy);
	printf("Prediction MSE     : %.3f\n", sudoku->mse);
	printf("Rating             : %3d\n\n", sudoku->rating);
}


void ParseCSV()
{
	Sudoku sudoku = CreateSudokuGrid();
	FILE *fp = fopen("sudoku-data/test.csv", "r");
	assert(fp != NULL);
	char line[512];
	int lineCount = 0;
	double *rowScaling = calloc(9, sizeof(double));
	double *colScaling = calloc(9, sizeof(double));
	int sinkhornIterations = 400;
	double sinkhornTolerance = 1e-10;
	int maxSolverIterations = 50;
	while(fgets(line, sizeof(line), fp))
	{
		//Skip header
		if(lineCount == 0){lineCount++;continue;}
		line[strcspn(line, "\r\n")] = 0;

		//Tokenize
		char *source = strtok(line, ",");
		char *question = strtok(NULL, ",");
		char *answer = strtok(NULL, ",");  
		char *ratingString = strtok(NULL, ",");
		if(!source || !question || !ratingString) {fprintf(stderr, "Malformed line at %d\n", lineCount);continue;}

		UpdateSudokuGrid(sudoku, question, answer, ratingString);
		
		//Initialize probabilities from the original question grid
		FindCellProbabilities(sudoku);
		for(int epoch = 0; epoch < 100; epoch += 1)
		{
			FindConstraintMatrices(sudoku);
			//Balance all constraints
			for(int constraintIndex = 0; constraintIndex < 27; constraintIndex++)
			{
				double sinkhornError = SinkhornKnoppAlgorithm(9, 9, sinkhornIterations, sinkhornTolerance, sudoku->constraintMatrices[constraintIndex], sudoku->constraintBalanced[constraintIndex], rowScaling, colScaling);
			}
			//Update all cell probabilities
			UpdateNewCellProbabilities(sudoku);
			//Find Sinkhorn answer
			ExtractSinkhornAnswer(sudoku);
			
		}
		PrintMetrics(lineCount-1,sudoku);
		//PrintGrid(sudoku, 1);
		//PrintGrid(sudoku, 2);
		//PrintSingleConstraintMatrix(sudoku->constraintMatrices[0]);
		//PrintSingleConstraintMatrix(sudoku->constraintBalanced[0]);
		//PrintCellProbabilities(sudoku);
		//PrintGrid(sudoku, false);
		lineCount += 1;
		if(lineCount == 30){break;}
	}
	fclose(fp);
	DestroySudoku(sudoku);
	free(rowScaling);free(colScaling);
}

int main()
{	
	ParseCSV();
	return 0;
}

