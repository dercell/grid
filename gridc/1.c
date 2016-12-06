#include <stdio.h>


int size=10;

int main()
{	
	int i=0, j=0, elem;
	
	int grid[size][size];
	FILE* fin=fopen("matrix.txt","r");
	char str;	

	while(str != EOF)
	{
		str=fgetc(fin);
		if((str != '\r') && (str != '\n'))
		{	
			elem=(int)str-48;
			if(i<10)
			{
				if(j<10)
				{	
					grid[i][j]=elem;
					j++;
				}
				else
				{
					i++;
					j=0;
					grid[i][j]=elem;
					j++;
				}
			}
			

		}

	
	}

	 for(i=0;i<size;i++)
                        {
                                for(j=0;j<size;j++)
                                        printf("%d", grid[i][j]);
				printf("\n");
                        }
		
	
	fclose(fin);
	return 0;
}
