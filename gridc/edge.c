#include <stdio.h>

int size=24;

int edge(int matrix[size][size], int i, int j)
{
	int i1, j1,k=0, cluster_edge=1, m,n,p,q, zeros, zeros_max, imax, jmax, iprev, jprev;
	int pt[8];	

	i1=i;
	j1=j;
	
	//printf("%d %d %d\n", i1, j1, matrix[i1][j1]);
	iprev=-1;
	jprev=-1;
	while(i1<size-3)
	{	
		//printf("%d %d\n", i1, j1);	
		zeros_max=0;
		imax=-1;
		jmax=-1;
		for(n=0;n<3;n++){
			for(m=0;m<3;m++)
			{	
				//printf("%d", matrix[i1+n-1][j1+m-1]);		
				if(matrix[i1+n-1][j1+m-1]==1)
				{	
					//printf("%d\n", matrix[i1+n-1][j1+m-1]);
				
					zeros=0;
					for(p=0;p<3;p++)
					{
						for(q=0;q<3;q++)
						{
							if(matrix[i1+n+p-2][j1+m+q-2]==0)
								zeros++;
						}				
					}
					
					//printf("i1+n-1=%d n=%d j1+m-1=%d, m=%d, zeros=%d, max=%d\n", i1+n-1, n, j1+m-1, m, zeros, zeros_max);	
					if(zeros > zeros_max && (i1+n-1!=iprev || j1+m-1!=jprev) && (n!=1 || m!=1))
					{
						
						//printf("%d %d=>%d\n",i1+n-1, j1+m-1, zeros);	
						zeros_max=zeros;
						imax=i1+n-1;
						jmax=j1+m-1;
					}
				}	
 			}
			//printf("\n");
		}
		iprev=i1;
		jprev=j1;
		i1=imax;
		j1=jmax;
		cluster_edge++;
		printf("%d %d, %d\n", i1,j1, matrix[i1][j1]);
		
	}

	return cluster_edge;	
}

int main()
{	
	int matrix[size][size], i, j, ed;

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			if(j<=i && j>=2 && i< size/2 && i!=1 && i!=0 && i!=size-2 && i!=size-1)
				matrix[i][j]=1;
			else matrix[i][j]=0;
			if(i==12 && j>=i && j>=2 && j <=18)
				matrix[i][j]=1;
			if(i>12 &&  i<22 && j>=2 && j<=18)
				matrix[i][j]=1;
		printf("%d", matrix[i][j]);
		}
		printf("\n");
	}
	getchar();	
	ed=edge(matrix, 2, 2);
	printf("%d\n", ed);
	return 0;
}
