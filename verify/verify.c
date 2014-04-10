#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
	char infile1[256], infile2[256];
	FILE  *fp1, *fp2;
	double *buf, *dat, sub;
	int   i,k1,k2,ns,num=0;
	double diffsum=0.0,diffsum2=0.0,file1sum=0.0,file1sum2=0.0,file2sum=0.0,diffmax=0.0;
	double L1,L2,RMS;
	if(argc<2)
	{
		printf("please add 2 parameter: file1 file2\n");
		return 0;
	}
	strcpy(infile1,argv[1]);
	strcpy(infile2,argv[2]);

	fp1 = fopen(infile1,"r");
	if(fp1 == NULL)
	{
		printf("file %s is  not exist\n",infile1);
		exit(0);
	}

	fp2 = fopen(infile2,"r");
	if(fp2 == NULL)	
	{
		printf("file %s is  not exist\n",infile2);
		exit(0);
	}


	ns=1;
	buf=(double *)malloc(sizeof(double)*ns);
	dat=(double *)malloc(sizeof(double)*ns);

	fseek(fp1, 0 ,0);
	fseek(fp2, 0 ,0);
	for(;;)
	{
		k1=fread(buf,sizeof(double), ns, fp1);
		k2=fread(dat,sizeof(double), ns, fp2);
		if(k1!=ns) break;
		for(i=0;i<ns;i++)
		{
                        file1sum += fabs(buf[i]);
			file2sum += fabs(dat[i]);
		        file1sum2 += buf[i]*buf[i];
			sub=(buf[i]-dat[i]);
			if(fabs(sub)>diffmax) 
			{
				diffmax=fabs(sub);
			}
		  	diffsum += fabs(sub);
		  	diffsum2 += sub*sub;
		}
		num+=ns;
	}
	if(file1sum!=0)
	{
		L1=diffsum/file1sum;
		L2=diffsum2/file1sum2;
		RMS=sqrt(L2);
	}
	else
	{
		L1=0.0;
		L2=0.0;
		RMS=0.0;
	}
	printf(" L1 model = %.10lf\n",L1);
	printf(" L2 model = %.10lf\n",L2);
	printf(" RMS model = %.10lf\n",RMS);
	printf(" file1sum = %.10lf\n",file1sum);
	printf(" file2sum = %.10lf\n",file2sum);

	free(buf);
	fclose(fp1);
	fclose(fp2);

}

