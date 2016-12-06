#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "_mt19937.c"
//#include "graphics.h"

#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"
#define YEL  "\x1B[33m"
#define BLU   "\x1B[34m"
#define RST  "\x1b[37m"

mt19937_state state;

int size = 44;
float freq = 0.9;
int N = 500;
int R = 1;
float T;
int cluster_number=0;
int *LL=NULL;
int *value;
int *perc_arr=NULL;
int *val;
int perc_size=0;
int *per_arr=NULL;
int *v;
int per_num=0;

//добавление нового кластера
void add()
{ 	      
	value=(int*)realloc(LL,cluster_number*sizeof(int));
	LL=value;
	LL[cluster_number-1]=1;
}

void perc_add()
{	
	val=(int*)realloc(perc_arr, perc_size*sizeof(int));
	perc_arr=val;

}

void per_add(int per)
{
	v=(int*)realloc(per_arr, per_num*sizeof(int));
	per_arr=v;
	per_arr[per_num-1]=per;

}

void display(short int gr[size][size], int cluster1[size][size], int cluster0[size][size]){
	int i,j,k;
	for(i=2;i<size-2;i++){
		for(j=2;j<size-2;j++){
			if(gr[i][j]==1)
			{
				printf(BLU "%c" RST, 219);
			}
				
			else
				printf(RED "%c" RST, 254);
			
			//printf("%d", gr[i][j]);
			}
		printf("\t");
		for(j=2;j<size-2;j++)
		{	
			if(cluster1[i][j]==-1)
				printf(YEL "%c"RST, 254);
			if(cluster1[i][j]==0)
				//printf("%d", find_root(cluster1[i][j], LL));
				printf( RED "%c" RST, 254);
			if(cluster1[i][j]!=0 && cluster1[i][j]!=-1)
				printf(BLU "%c" RST, 254);
		}
		printf("\t");
		for(j=2;j<size-2;j++)
		{
			printf("%d", find_root(cluster0[i][j], LL));
		}
		printf("\n");
	}
	printf("\n\n");
}




int find_root(int num, int* mass)
{	
	
	int root=num;
	while(mass[root-1]<0) 
	{  
		root=-mass[root-1]; 
	}
	return root;
}

//остальные строки

void other_raws(short int grid[size][size], int cluster[size][size], int param)
{
	int i, j;
	
	for(i=2;i<size-2;i++)
	{	

		for(j=2;j<size-2;j++)
		{	
			if(grid[i][j]==param)
			{	
		
				int n, min=0, sz=8, root, neighbour;
				int site[8];
				//соседние клетки	
				site[0]=cluster[i][j-1];
				site[1]=cluster[i-1][j-1];
				site[2]=cluster[i-1][j];
				site[3]=cluster[i-1][j+1];
				site[4]=0;
				site[5]=0;
				site[6]=0;
				site[7]=0;
				//если последняя строка, то сравниваем с элементами первой
				if(i==size-3)
				{
					site[5]=cluster[2][j-1];
					site[6]=cluster[2][j];
					site[7]=cluster[2][j+1];
					sz=8;
				}
				//если последний столбец, то сравниваем с элементами вначале 
				if(j==size-3)
				{
					site[3]=cluster[i-1][2];
					site[4]=cluster[i][2];
					sz=5;
				}
				//если первый столбец, то сравниваем с последним элементом строки выше(фактически он по диагонали влево)
				if(j==2)
				{
					site[0]=0;
					site[1]=cluster[i-1][size-3];
					sz=3;
				}

				//поиск кластера с минимальным номером среди соседей
				for(n=0;n<sz;n++)
				{	
					if(site[n]!=0)
					{	
						neighbour=find_root(site[n], LL);
						if(min==0) { min=neighbour; }
						else { if(neighbour<min) { min=neighbour; } }
					}
				}
				//если кластер с мин номером найден, объединяем  кластеры
				if(min!=0)
				{	
					//текущая клетка получает лейбл минимального из соседей
					cluster[i][j]=min;
					//ищем номер корневого кластера

					//root=find_root(cluster[i][j]);

					//увеличивает его "массу" на единицу
					LL[min-1]++;
					//делаем номер корневого минимумом, чтобы сравнить с корневыми соседних кластеров
					//min=root;
					

					for(n=0;n<sz;n++)
					{	

						//если это ненулевой элемен и не элемент этого же кластера
						if(site[n]!=0 && site[n]!=min)
						{	

							//ищем номер его корневого кластера
							root=find_root(site[n], LL);
							//если он больше, то прибавляем его "массу" к первому
							//и записываем вместо его массы отрицательный номер кластера,
							//в который он слился
							if(root>min)
							{	

								LL[min-1]+=LL[root-1];
								LL[root-1]=-min;
							}
							//если он оказался меньше, то та же процедура,
							//только первый кластер сливается во второй
							else
							{	

								//условие, если вдруг это один и тот же кластер, на всякий случай
								if(root!=min)
								{	

								LL[root-1]+=LL[min-1];
								LL[min-1]=-root;
								}
							}
						}
					}
				}
				//если не найден - создаем новый кластер
				else { cluster_number++; cluster[i][j]=cluster_number; add(); }
	 
			}
			else cluster[i][j]=0;
		}
	}

}




void f_gen(short int gr[size][size]){
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

float frequency(short int gr[size][size]){
	int i,j,k=0;
	for(i=2;i<size-2;i++){
		for(j=2;j<size-2;j++){
			if(gr[i][j]==1)
				k++;}}
	return k*1.0/((size-4)*(size-4));

}

void boundary(short int gr[size][size]){
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

float game(short int gr[size][size], int i, int j, float T){
	int p,q; 
	float payoff=0;
	for(p=0;p<3;p++){
		for(q=0;q<3;q++){
			if(gr[i][j]==1){
				if(gr[i+p-1][j+q-1]==gr[i][j])
					payoff += R;
				}			
			else {
				if(gr[i+p-1][j+q-1]!=gr[i][j])
					payoff += T;
				}
			}
		}
	return payoff;
	
}

void generations(short int gr[size][size], int counter, float T){
	short int new_grid[size][size];
	int i,j,p,q, i_max, j_max; 
	float payoff, payoff_max;

#pragma omp parallel default(shared)
{
	#pragma omp for private(i,j, payoff, payoff_max, i_max, j_max, p,q)
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
			new_grid[i][j]=gr[i_max][j_max];
		}
	}
#pragma omp barrier
}
	boundary(new_grid);

	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			gr[i][j]=new_grid[i][j];
		}
	}
}


void route(int cluster1[size][size], int cluster0[size][size], int root1, int root2, int* cop, int* def)
{
	int root, i, j, nb[4];

	for(i=2;i<size-2;i++)
	{
		for(j=2;j<size-2;j++)
		{	
			if(cluster1[i][j]!=0)
			{	
				//printf("ifclus1\t");
				nb[0]=cluster0[i][j-1];
				nb[1]=cluster0[i][j+1];
				nb[2]=cluster0[i-1][j];
				nb[3]=cluster0[i+1][j];
				if(i==2)
					nb[2]=cluster0[size-3][j];
				if(i==size-3)
					nb[3]=cluster0[2][j];
				if(j==2)
					nb[0]=cluster0[i][size-3];
				if(j==size-3)
					nb[1]=cluster0[i][2];
				//printf("root\t");
				root=find_root(cluster1[i][j], cop);
				//printf("ifroot\n");
				if(root==root1)
				{	
					//printf("iside_if_root\t");
					if(find_root(nb[0], def)==root2 || find_root(nb[1], def)==root2 \
					|| find_root(nb[2], def)==root2 || find_root(nb[3], def)==root2)
					cluster1[i][j]=-1;
					//printf("after_per++\n");
				}
			}

		}
	}

}

int edge(short int grid[size][size], int cluster1[size][size], int cluster0[size][size], int i1, int j1, int len)
{
	
	int flag=0, j=0;
	cluster1[i1][j1] = 1;
	
	display(grid, cluster1, cluster0);
	getchar();	
	

	if(i1 == size-3)
		return len;
	
	if(j1 == size-3)
		j = 2;
	else
		j = j1 + 1;
	
	if((cluster1[i1][j] == -1)  && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1, j, len);
	}
	
	if((cluster1[i1+1][j] == -1) && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1+1, j, len);
	}
	
	if((cluster1[i1+1][j1] == -1) && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1+1, j1, len);
	}
	
	if(j1 == 2)
		j = size - 3;
	else
		j = j1 - 1;	

	if((cluster1[i1+1][j] == -1) && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1+1, j, len);
	}
	
	if((cluster1[i1][j] == -1) && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1, j, len);
	}
	
	if((cluster1[i1-1][j] == -1) && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1-1, j, len);
	}
	
	if((cluster1[i1-1][j1] == -1) && (flag == 0)){
		len++;
		flag=edge(grid, cluster1, cluster0, i1-1, j1, len);
	}
	
	if(j1 == size -3)
		j = 2;
	else 
		j = j1 + 1;

	if((cluster1[i1-1][j] == -1) && (flag == 0 )){
		len++;
		flag=edge(grid, cluster1, cluster0, i1-1, j, len);
	}
		   
	if(flag != 0 ){
		return flag;
	}
	else
		return 0;
}

int perimeter(short int grid[size][size], int cluster1[size][size], int cluster0[size][size], int root1, int root2, int* def, int* cop)
{
	
	int i1, j1,k=0, flag=0, flag2=0, flag3=0, k1, jstart=2, cluster_edge=0, i, j, ilast, jlast;
	
	
	route(cluster1, cluster0, root1, root2, cop, def);

	for(j=2;j<size-2;j++)
	{
		if(cluster1[2][j]==-1)
		{
			i1=2; j1=j;
			jstart=j;
			//printf("%d %d\n", i1, j1);
			cluster_edge=edge(grid, cluster1, cluster0, i1, j1, cluster_edge);
			break;
		}
		if(j==size-3)
			return 0;
	}
	
	
	return cluster_edge;
}

int perc(int cluster[size][size], int i, int j, int root)
{	
	int flg=0;
	
	if(i!=size-3)
	{
	if(cluster[i+1][j]!=0 && find_root(cluster[i+1][j], LL)==root && flg<0)
		flg=perc(cluster, i+1,j, root);
	if(cluster[i+1][j+1]!=0 && find_root(cluster[i+1][j+1], LL)==root && flg<0)
		flg=perc(cluster, i+1, j+1, root);
	if(cluster[i+1][j-1]!=0 && find_root(cluster[i+1][j-1],LL)==root && flg<0)
		flg=perc(cluster, i+1, j-1, root);
	if(cluster[i][j+1]!=0 && find_root(cluster[i][j+1],LL)==root && flg<0)
		flg=perc(cluster, i, j+1, root);
	if(cluster[i][j-1]!=0 && find_root(cluster[i][j-1],LL)==root && flg<0)
		flg=perc(cluster, i, j-1, root);
	if(cluster[i-1][j+1]!=0 && find_root(cluster[i-1][j+1],LL)==root && flg<0)
		flg=perc(cluster, i-1, j+1, root);
	if(cluster[i-1][j-1]!=0 && find_root(cluster[i-1][j-1],LL)==root && flg<0)
		flg=perc(cluster, i-1,j-1, root);
	if(cluster[i-1][j]!=0 && find_root(cluster[i-1][j],LL)==root && flg<0)
		flg=perc(cluster, i-1, j, root);
	
		return -1;
	}
	else
		return root;
}

int percolate(int cluster[size][size])
{	
	perc_size=0;
	int i, j, k, flag=0, root1, root2, per=0;
	
	for(j=2;j<size-2;j++)
	{	
		flag=0;
		if(cluster[2][j]!=0)
		{	
			root1=find_root(cluster[2][j], LL);
			for(i=0;i<perc_size;i++)
			{	
				if(root1==perc_arr[i])
				flag=1;
			}
			
			if(flag!=1)
			{	
				for(k=2;k<size-2;k++)
					{	
						if(cluster[size-3][k]!=0)
						root2=find_root(cluster[size-3][k], LL);
						if(root1==root2 && LL[root2-1]>size)
						{	
							return root1;
						}	
					}
				//root2=perc(cluster,2,j, root1);
			}	
		}
	//	if(root2>0)
			//return root2;
		
	}
	return -1;

}

double average(double* mass, int n, int flag)
{
	int i;
	double sum=0, av=0, up=0, var=0;
	for(i=0;i<n;i++)
	{	
		sum+=mass[i];
	}
	
	if(n!=0)
	av=sum/n;
	else
	av=0;
	
	for(i=0;i<n;i++)
	{
		up+=(mass[i]-av)*(mass[i]-av);
	}
	
	if(n!=0)
	var=sqrt(up/n);	
	else
	var=0;
	//if(flag)
	//printf("%f\t%f\t", av, var);
	return av;	
}


double int_average(int* mass, int n)
{
	int i;
	double sum=0, av=0, up=0, var=0, down=0;
	char aver[10], variance[10], t_str[10], path[10];
	sprintf(path, "%d", size);
	strcat(path, ".txt");
	FILE* fin=fopen(path, "a");
	sprintf(t_str, "%f", T);
	fputs(t_str, fin);
	fputs("\t", fin);
	
	//printf("sum\t");
	for(i=0;i<n;i++)
	{	
		if(mass[i]!=0)
		{
			sum+=mass[i];
			down++;
		}
	}
	//printf("av\t");
	if(n!=0 &&  down!=0)
		av=sum/down;
	else
		av=0;
	//printf("up\t");
	for(i=0;i<n;i++)
	{	
		if(mass[i]!=0)
			up+=(mass[i]-av)*(mass[i]-av);
	}
	//printf("var\n");
	if(n!=0 && down!=0)
		var=sqrt(up*1.0/down);	
	else
		var=0;
	printf("%f\t %f\n", av, var);
	sprintf(aver, "%f", av);
	sprintf(variance, "%f", var);

	fputs(aver, fin);
	fputs("\t", fin);
	fputs(variance, fin);
	fputs("\n", fin);
	//printf("files\n");
	fclose(fin);
	return av;	
}


int count(short int grid[size][size], int param)
{	
	int i,j,k, res=0, c, interval=0, max=0, coop, def;
	int cluster1[size][size], cluster0[size][size]; 
	
	for(i=0;i<size;i++)
	{ for(j=0;j<size;j++)
	  { if(i==0 || i==1 || i==size-1 || i==size-2 || j==0 || j==1 || j==size-1 || j==size-2)
		cluster0[i][j]=0;
		cluster1[i][j]=0;
	  }
	}
	//printf("coop hosh\t");
	cluster_number=0;
	
	other_raws(grid, cluster1, 1);
	//printf("coop perc\n");
	coop=percolate(cluster1);

	int cop[cluster_number], cop_size;
	cop_size=cluster_number;
	for(i=0;i<cluster_number;i++)
		cop[i]=LL[i];
	
	cluster_number=0;
	other_raws(grid, cluster0, 0);
	def=percolate(cluster0);
	
	int deff[cluster_number], deff_size;
	deff_size=cluster_number;
	for(i=0;i<cluster_number;i++)
		deff[i]=LL[i];
	
	if(coop>0 && def>0)
	{
		max=perimeter(grid, cluster1, cluster0, coop, def, deff, cop);
		
		//printf("%d\n", max);
	}
	
	return max;
	
}


int main(int argc, char *argv[]){
	short int grid[size][size], i,j, k=0, m=0, c=0;
	int c_clusters[15000], d_clusters[15000];
	double fc[N], freq[300], av_freq[50];
	char str[15], t_str[10], path[10], val[10];
	mt19937_init_sequence_(&state, time(NULL));
	//FILE *fin=fopen("edge20.txt", "a");

	T=1.799999;
		for(j=0;j<50;j++)
		{	
			f_gen(grid);
			boundary(grid);
			for(i=0;i<N;i++)
			{	
				//printf("%d\n", i);
				generations(grid, i,T);
				if(i>=200)
					{ 

					freq[i-200]=frequency(grid);

					c_clusters[m]=count(grid,1);

					m++;
					//getchar(); 	
					}
			}
		}
		int_average(c_clusters, 15000);		
		m=0;
	
//	fclose(fin);
	free(LL);
	free(perc_arr);
	return 0;
}
	
