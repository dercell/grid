#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include "_mt19937.c"

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"

mt19937_state state;
int size = 54;
float freq = 0.9;
int N = 500;
int R = 1;
float T;

void f_gen(int gr[size][size]){
	int i,j;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			if(mt19937_generate_uniform_float_(&state)<freq)
				gr[i][j]=1;
			else
				gr[i][j]=0;
			}
	}
}
void display(int gr[size][size]){
	int i,j;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			if(gr[i][j]==1)
				printf("%s%d",KGRN, gr[i][j]);
			else
				printf("%s%d",KRED, gr[i][j]);
			}
		printf("\n");
	}
}

void frequency(int gr[size][size], int counter, float T){
	int i,j,k=0;
	char str[15], t_str[10], path[10];
	FILE* file=fopen("data.txt","a");
	sprintf(t_str, "%f", T);
	sprintf(path,"%f", T);
	strcat(path, ".txt");
	
	FILE* file2=fopen(path, "a");
	for(i=2;i<size-2;i++){
		for(j=2;j<size-2;j++){
			if(gr[i][j]==1)
				k++;}}
	//printf("%f\n", k*1.0/((size-4)*(size-4)));
	sprintf(str, "%f", k*1.0/((size-4)*(size-4)));
	fputs(t_str, file2);
	fputs(" ", file2);
	fputs(str, file2);
	fputs("\n", file2);
	if(counter==499){
		fputs("\n", file2);
	}
	fputs(str, file);
	fputs("\n", file);
	fclose(file);
	fclose(file2);

}

void boundary(int gr[size][size]){
	int i;
	for(i=2;i<size-2;i++){
		gr[0][i]=gr[size-4][i];
		gr[1][i]=gr[size-3][i];
		gr[size-2][i]=gr[2][i];
		gr[size-1][i]=gr[3][i];

		gr[i][0]=gr[i][size-4];
		gr[i][1]=gr[i][size-3];
		gr[i][size-2]=gr[i][2];
		gr[i][size-1]=gr[i][3];
	}
	gr[0][0]=gr[size-4][size-4];
	gr[1][1]=gr[size-3][size-3];
	gr[0][1]=gr[size-4][size-3];
	gr[1][0]=gr[size-3][size-4];

	gr[size-1][size-1]=gr[3][3];
	gr[size-2][size-2]=gr[2][2];
	gr[size-2][size-1]=gr[2][3];
	gr[size-1][size-2]=gr[3][2];
	
	gr[0][size-2]=gr[size-4][2];
	gr[0][size-1]=gr[size-4][3];
	gr[1][size-2]=gr[size-3][2];
	gr[1][size-1]=gr[size-3][3];

	gr[size-2][0]=gr[2][size-4];
	gr[size-1][0]=gr[3][size-4];
	gr[size-2][1]=gr[2][size-3];
	gr[size-1][1]=gr[3][size-3];
}

float game(int gr[size][size], int i, int j, float T){
	int p,q; 
	float payoff=0;
	//printf("%d, %d\n", i, j);
	for(p=0;p<3;p++){
		for(q=0;q<3;q++){
			//printf("%d", gr[i+p-1][j+q-1]);
			if(gr[i][j]==1){
				if(gr[i+p-1][j+q-1]==gr[i][j])
					payoff += R;
				}			
			else {
				if(gr[i+p-1][j+q-1]!=gr[i][j])
					payoff += T;
				}
			}
		//printf("\n");
		
		}
	//printf("payoff %f\n", payoff); 
	return payoff;
	
}

void generations(int gr[size][size], int counter, float T)
{
	int new_grid[size][size];
	int i,j,p,q, i_max, j_max; 
	float payoff, payoff_max;
#pragma omp parrallel shared(gr, new_grid, size)
{
#pragma omp for private(i, payoff, j, j_max, i_max, payoff_max, p,q)
	for(i=2;i<size-2;i++){
		for(j=2;j<size-2;j++){
			payoff=0;
			payoff_max=0;
			i_max=-1;
			j_max=-1;
			for(p=0;p<3;p++){
				for(q=0;q<3;q++){
					payoff=game(gr, i+p-1, j+q-1, T);
					if(payoff > payoff_max){
						payoff_max=payoff;
						i_max=i+p-1;
						j_max=j+q-1;
					}
 				}
			}
			//printf("max %d , %f\n ", gr[i_max][j_max], payoff_max);
			new_grid[i][j]=gr[i_max][j_max];
		}
	}
}
//#pragma omp barrier
//#pragma omp master
//{
	boundary(new_grid);
//}
	//#pragma omp for private(i,j)
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			gr[i][j]=new_grid[i][j];
		}
	}
//}
	if(counter>=200){
		frequency(gr, counter, T);
	}
}

int main(int argc, char *argv[]){
	int grid[size][size], i,j, k=0;
	float fc[N];
	mt19937_init_sequence_(&state, time(NULL));
	unlink("data.txt");
	T=1.8;

		
		for(j=0;j<50;j++)
		{	
			f_gen(grid);
			boundary(grid);
			for(i=0;i<N;i++)
			{
				generations(grid, i,T);
				display(grid);
				getchar();
				//printf("\n");
				//sleep(1);
				system("clear");
			}
		}
	//T+=0.0001;
	
	
	return 0;
}
