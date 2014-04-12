#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#include <omp.h>

#define PIE 3.1415926

#define DEBUG_CPU_RUNNING

#define DEBUG_NO_PARALLEL

#define POSITION_INDEX_HOST_X(_z,_y,_x)        ((_z)*ny*nx + (_y)*nx + (_x))
#define POSITION_INDEX_HOST_Y(_z,_y,_x)        ((_x)*nz*ny + (_z)*ny + (_y))
#define POSITION_INDEX_HOST_Z(_z,_y,_x)        ((_y)*nx*nz + (_x)*nz + (_z))

#define POSITION_INDEX_X(_z,_y,_x)   ((_z)*(nMicMaxYLength+10)*(nMicMaxXLength+10) + (_y)*(nMicMaxXLength+10) + (_x))
#define POSITION_INDEX_Y(_z,_y,_x)   ((_x)*(nMicMaxZLength+10)*(nMicMaxYLength+10) + (_z)*(nMicMaxYLength+10) + (_y))
#define POSITION_INDEX_Z(_z,_y,_x)   ((_y)*(nMicMaxXLength+10)*(nMicMaxZLength+10) + (_x)*(nMicMaxZLength+10) + (_z))

#ifndef	DEBUG_CPU_RUNNING
#	define MIC_ALLOC       alloc_if(1) free_if(0)
#	define MIC_FREE        alloc_if(0) free_if(1)
#	define MIC_REUSE       alloc_if(0) free_if(0)
#	define MIC_SELF_MANAGE alloc_if(0) free_if(0)
#	define MIC_VAR         __attribute__((target(mic)))
#else
#	define MIC_VAR
#endif

#define ROUND_TO_SIZE(_length, _alignment)    \
			(((_length) + ((_alignment)-1)) & ~((_alignment) - 1))

#define ROUND_TO_SIZE_LESS(_length,_alignment)  \
			((_length) &~ ((_alignment)-1))

#define vpp2(_z) (((_z+ntop-5)<210)?(5290000):((_z+ntop-5)>=260?12250000:7840000))
#define vss2(_z) (((_z+ntop-5)<210)?(1517824):((_z+ntop-5)>=260?3644281:2277081))

typedef struct _MEMORY_BLOCKS{

 double * up, * up1, * up2, * vp, * vp1, * vp2, * wp, * wp1, \
	  * wp2, * us, * us1, * us2, * vs, * vs1, * vs2, * ws, * ws1, * ws2;
 double *u,*v,*w;
 double * to_write;
}MEMORY_BLOCKS,*PMEMORY_BLOCKS;


MIC_VAR int i,j,k,kk,kkk,l;
MIC_VAR int ishot,ncy_shot,ncx_shot;

MIC_VAR double *wave;
MIC_VAR double nshot,t0,tt,c0;
MIC_VAR double dtx,dtz;
MIC_VAR double px;

MIC_VAR	double vvp2_dtx_dtx;
MIC_VAR	double vvs2_dtz_dtz;
MIC_VAR	double vvs2_dtx_dtx;
MIC_VAR	double vvp2_dtz_dtz;
MIC_VAR	double vvp2_dtz_dtx;
MIC_VAR	double vvs2_dtz_dtx;

int nx,ny,nz,lt,nedge;
double frequency;
double velmax;
double dt;
int ncx_shot1,ncy_shot1,ncz_shot;
double unit;
int nxshot,nyshot,dxshot,dyshot;

char infile[80],outfile[80],logfile[80],tmp[80];
FILE  *fin, *fout, *flog;
struct timeval start,end;
double all_time;
double c[7][11];

int nSize;
int nSliceSize;

double * to_write;

int mic_used_size;
int mic_slice_size;

double *up_out;

void initailize(){
	strcpy(tmp,"date ");
	strncat(tmp, ">> ",3);
	strncat(tmp, logfile, strlen(logfile));
	flog = fopen(logfile,"w");
	fprintf(flog,"------------start time------------\n");
	fclose(flog);
	system(tmp);
	gettimeofday(&start,NULL);

	fin = fopen(infile,"r");
	if(fin == NULL)
	{
		printf("file %s is  not exist\n",infile);
		exit(0);
	}
	fscanf(fin,"nx=%d\n",&nx);
	fscanf(fin,"ny=%d\n",&ny);
	fscanf(fin,"nz=%d\n",&nz);
	fscanf(fin,"lt=%d\n",&lt);
	fscanf(fin,"nedge=%d\n",&nedge);
	fscanf(fin,"ncx_shot1=%d\n",&ncx_shot1);
	fscanf(fin,"ncy_shot1=%d\n",&ncy_shot1);
	fscanf(fin,"ncz_shot=%d\n",&ncz_shot);
	fscanf(fin,"nxshot=%d\n",&nxshot);
	fscanf(fin,"nyshot=%d\n",&nyshot);
	fscanf(fin,"frequency=%lf\n",&frequency);
	fscanf(fin,"velmax=%lf\n",&velmax);
	fscanf(fin,"dt=%lf\n",&dt);
	fscanf(fin,"unit=%lf\n",&unit);
	fscanf(fin,"dxshot=%d\n",&dxshot);
	fscanf(fin,"dyshot=%d\n",&dyshot);
	fclose(fin);

	printf("\n--------workload parameter--------\n");
	printf("nx=%d\n",nx);
	printf("ny=%d\n",ny);
	printf("nz=%d\n",nz);
	printf("lt=%d\n",lt);
	printf("nedge=%d\n",nedge);
	printf("ncx_shot1=%d\n",ncx_shot1);
	printf("ncy_shot1=%d\n",ncy_shot1);
	printf("ncz_shot=%d\n",ncz_shot);
	printf("nxshot=%d\n",nxshot);
	printf("nyshot=%d\n",nyshot);
	printf("frequency=%lf\n",frequency);
	printf("velmax=%lf\n",velmax);
	printf("dt=%lf\n",dt);
	printf("unit=%lf\n",unit);
	printf("dxshot=%d\n",dxshot);
	printf("dyshot=%d\n\n",dyshot);
	flog = fopen(logfile,"a");
	fprintf(flog,"\n--------workload parameter--------\n");
	fprintf(flog,"nx=%d\n",nx);
	fprintf(flog,"ny=%d\n",ny);
	fprintf(flog,"nz=%d\n",nz);
	fprintf(flog,"lt=%d\n",lt);
	fprintf(flog,"nedge=%d\n",nedge);
	fprintf(flog,"ncx_shot1=%d\n",ncx_shot1);
	fprintf(flog,"ncy_shot1=%d\n",ncy_shot1);
	fprintf(flog,"ncz_shot=%d\n",ncz_shot);
	fprintf(flog,"nxshot=%d\n",nxshot);
	fprintf(flog,"nyshot=%d\n",nyshot);
	fprintf(flog,"frequency=%lf\n",frequency);
	fprintf(flog,"velmax=%lf\n",velmax);
	fprintf(flog,"dt=%lf\n",dt);
	fprintf(flog,"unit=%lf\n",unit);
	fprintf(flog,"dxshot=%d\n",dxshot);
	fprintf(flog,"dyshot=%d\n\n",dyshot);
	fclose(flog);

	nSize = nz*ny*nx;
	nSliceSize = ny*nx;

	mic_used_size = pow(2.*((lt*dt*velmax)/unit+10.)+10.,3.);
	mic_slice_size = pow(2.*((lt*dt*velmax)/unit+10.)+10.,2.);

	printf("length: %lfmic_slice_size:%d mic_used_size:%d\n",2.*((lt*dt*velmax)/unit+10.)+10.,mic_slice_size, mic_used_size);

	up_out = (double *)malloc(mic_slice_size * sizeof(double));

	wave = (double*)malloc(sizeof(double)*lt);

	t0=1.0/frequency;
	for(l=0;l<lt;l++)
	{
		tt=l*dt;
		tt=tt-t0;
		double sp=PIE*frequency*tt;
		double fx=100000.*exp(-sp*sp)*(1.-2.*sp*sp);
		wave[l]=fx;
		// printf("wave l:%d:%lf",l,wave[l]);
	}

	c0=-2.927222164;
	c[0][0]=1.66666665;
	c[0][1]=-0.23809525;
	c[0][2]=0.03968254;
	c[0][3]=-0.004960318;
	c[0][4]=0.0003174603;
	c[1][0]=0.83333;
	c[1][1]=-0.2381;
	c[1][2]=0.0595;
	c[1][3]=-0.0099;
	c[1][4]=0.0008;

	for(i=0;i<5;i++)
		for(j=0;j<5;j++)
			c[2+i][j]=c[1][i]*c[1][j];

    for(j=0;j<7;++j){
        for(i=0;i<5;++i){
            c[j][10-i] = c[j][i];
        }
        c[j][5] = c0;
    }
}

void calc_shot(
	int ncx_shot,
	int ncy_shot,
	int lStart,int lEnd,
	PMEMORY_BLOCKS pMemBlocks
	)
{
	int i_begin,i_end,j_begin,j_end,k_begin,k_end;
	int nMicXLength,nMicYLength,nMicZLength;
	int nMicMaxXLength,nMicMaxYLength,nMicMaxZLength;

	double current_c;

	int n_mic_left,n_mic_right,n_mic_front,n_mic_back,n_mic_top,n_mic_bottom;

	double * to_write 	= pMemBlocks->to_write;
	double * up_out 	= ((double *)malloc(mic_slice_size * sizeof(double)));

	double * up  = pMemBlocks->up;
	double * up1 = pMemBlocks->up1;
	double * up2 = pMemBlocks->up2;
	double * vp  = pMemBlocks->vp;
	double * vp1 = pMemBlocks->vp1;
	double * vp2 = pMemBlocks->vp2;
	double * wp  = pMemBlocks->wp;
	double * wp1 = pMemBlocks->wp1;
	double * wp2 = pMemBlocks->wp2;
	double * us  = pMemBlocks->us;
	double * us1 = pMemBlocks->us1;
	double * us2 = pMemBlocks->us2;
	double * vs  = pMemBlocks->vs;
	double * vs1 = pMemBlocks->vs1;
	double * vs2 = pMemBlocks->vs2;
	double * ws  = pMemBlocks->ws;
	double * ws1 = pMemBlocks->ws1;
	double * ws2 = pMemBlocks->ws2;
	double * u   = pMemBlocks->u;
	double * v   = pMemBlocks->v;
	double * w   = pMemBlocks->w;

	double xmax=lEnd*dt*velmax;
	int nleft=ncx_shot-xmax/unit-10;
	int nright=ncx_shot+xmax/unit+10;
	int nfront=ncy_shot-xmax/unit-10;
	int nback=ncy_shot+xmax/unit+10;
	int ntop=ncz_shot-xmax/unit-10;
	int nbottom=ncz_shot+xmax/unit+10;

	ntop = ntop-1;
	nfront = nfront-1;
	nleft = nleft-1;

	if(nleft<5) nleft=5;
	if(nright>nx-5) nright=nx-5;
	if(nfront<5) nfront=5;
	if(nback>ny-5) nback=ny-5;
	if(ntop<5) ntop=5;
	if(nbottom>nz-5) nbottom=nz-5;

	nMicMaxXLength = nright  - nleft;
	nMicMaxYLength = nback   - nfront;
	nMicMaxZLength = nbottom - ntop;

	// printf("MAX_X: %d, MAX_Y: %d, MAX_Z: %d", nMicMaxXLength, nMicMaxYLength, nMicMaxZLength);

	memset(u  ,0,sizeof(double)*mic_used_size);
	memset(v  ,0,sizeof(double)*mic_used_size);
	memset(w  ,0,sizeof(double)*mic_used_size);
	memset(up ,0,sizeof(double)*mic_used_size);
	memset(up1,0,sizeof(double)*mic_used_size);
	memset(up2,0,sizeof(double)*mic_used_size);
	memset(vp ,0,sizeof(double)*mic_used_size);
	memset(vp1,0,sizeof(double)*mic_used_size);
	memset(vp2,0,sizeof(double)*mic_used_size);
	memset(wp ,0,sizeof(double)*mic_used_size);
	memset(wp1,0,sizeof(double)*mic_used_size);
	memset(wp2,0,sizeof(double)*mic_used_size);
	memset(us ,0,sizeof(double)*mic_used_size);
	memset(us1,0,sizeof(double)*mic_used_size);
	memset(us2,0,sizeof(double)*mic_used_size);
	memset(vs ,0,sizeof(double)*mic_used_size);
	memset(vs1,0,sizeof(double)*mic_used_size);
	memset(vs2,0,sizeof(double)*mic_used_size);
	memset(ws ,0,sizeof(double)*mic_used_size);
	memset(ws1,0,sizeof(double)*mic_used_size);
	memset(ws2,0,sizeof(double)*mic_used_size);

	for(l=lStart;l<=lEnd;l++) {
		printf("l %d started\n",l);
		xmax=l*dt*velmax;
		n_mic_left=ncx_shot-xmax/unit-10;
		n_mic_right=ncx_shot+xmax/unit+10;
		n_mic_front=ncy_shot-xmax/unit-10;
		n_mic_back=ncy_shot+xmax/unit+10;
		n_mic_top=ncz_shot-xmax/unit-10;
		n_mic_bottom=ncz_shot+xmax/unit+ 10;

		--n_mic_left;
		--n_mic_front;
		--n_mic_top;

		if(n_mic_left<5) n_mic_left=5;
		if(n_mic_right>nx-5) n_mic_right=nx-5;
		if(n_mic_front<5) n_mic_front=5;
		if(n_mic_back>ny-5) n_mic_back=ny-5;
		if(n_mic_top<5) n_mic_top=5;
		if(n_mic_bottom>nz-5) n_mic_bottom=nz-5;

		// printf("[L]Starting l xmax:%d n_mic_left:%d n_mic_right:%d n_mic_top:%d n_mic_bottom:%d n_mic_front:%d n_mic_back:%d....%d\n",
		// 	xmax,n_mic_left,n_mic_right,n_mic_top,n_mic_bottom,n_mic_front,n_mic_back,l);

		//
		//	此处n_mic_XXX 系列变量同Host上的实际值相等。
		//  此前申请控空间时已经考虑留出了5的边界
		//	故此处循环应该从5开始。
		//
		nMicXLength = n_mic_right  - n_mic_left   ;
		nMicYLength = n_mic_back   - n_mic_front  ;
		nMicZLength = n_mic_bottom - n_mic_top 	  ;

		i_begin = 5+n_mic_left-nleft;
		i_end	= n_mic_left-nleft+nMicXLength+5;

		k_begin	= 5+n_mic_top - ntop;
		k_end	= n_mic_top - ntop+nMicZLength+5;

		j_begin = 5+n_mic_front-nfront;
		j_end	= n_mic_front-nfront+nMicYLength+5;

        #ifndef DEBUG_NO_PARALLEL
        #pragma omp parallel for private(i,j,k)
        #endif
		for(k=k_begin;k<k_end;k++)
		{

		    double vvp2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
							tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;

		    double _tempuxz,_tempuxy,_tempvyz,_tempvxy,_tempwxz,_tempwyz;

			vvp2=vpp2(k);
			vvs2=vss2(k);

			vvs2_dtz_dtz = vvs2*dtz*dtz;
			vvp2_dtx_dtx = vvp2*dtx*dtx;
			vvs2_dtx_dtx = vvs2*dtx*dtx;
			vvp2_dtz_dtz = vvp2*dtz*dtz;
			vvp2_dtz_dtx = vvp2*dtz*dtx;
			vvs2_dtz_dtx = vvs2*dtz*dtx;

			for(j=j_begin;j<j_end;j++) {
				for(i=i_begin;i<i_end;i++) {

					// printf("i:%d j:%d k:%d\n",i-5+nleft,j-5+nfront,k-5+ntop);
					int nIndex = POSITION_INDEX_X(k,j,i);

					if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
					{
						px=1.;
					}
					else
					{
						px=0.;
					}

					tempux2=0.0f;
					tempuy2=0.0f;
					tempuz2=0.0f;
					tempvx2=0.0f;
					tempvy2=0.0f;
					tempvz2=0.0f;
					tempwx2=0.0f;
					tempwy2=0.0f;
                	tempwz2=0.0f;
					tempuxz=0.0f;
					tempuxy=0.0f;
					tempvyz=0.0f;
					tempvxy=0.0f;
					tempwxz=0.0f;
					tempwyz=0.0f;

					for(kk=1;kk<=5;kk++)
					{
						tempux2=tempux2+c[0][kk-1]*(u[nIndex+kk]+u[nIndex-kk]);

						tempvx2=tempvx2+c[0][kk-1]*(v[nIndex+kk]+v[nIndex-kk]);

						tempwx2=tempwx2+c[0][kk-1]*(w[nIndex+kk]+w[nIndex-kk]);

						tempuy2=tempuy2+c[0][kk-1]*(u[POSITION_INDEX_X(k,j+kk,i)]+u[POSITION_INDEX_X(k,j-kk,i)]);

						tempvy2=tempvy2+c[0][kk-1]*(v[POSITION_INDEX_X(k,j+kk,i)]+v[POSITION_INDEX_X(k,j-kk,i)]);

						tempwy2=tempwy2+c[0][kk-1]*(w[POSITION_INDEX_X(k,j+kk,i)]+w[POSITION_INDEX_X(k,j-kk,i)]);

						tempuz2=tempuz2+c[0][kk-1]*(u[POSITION_INDEX_X(k+kk,j,i)]+u[POSITION_INDEX_X(k-kk,j,i)]);

						tempvz2=tempvz2+c[0][kk-1]*(v[POSITION_INDEX_X(k+kk,j,i)]+v[POSITION_INDEX_X(k-kk,j,i)]);

						tempwz2=tempwz2+c[0][kk-1]*(w[POSITION_INDEX_X(k+kk,j,i)]+w[POSITION_INDEX_X(k-kk,j,i)]);

					} //for(kk=1;kk<=5;kk++) end

					tempux2=(tempux2+c0*u[nIndex])*vvp2_dtx_dtx;

					tempvx2=(tempvx2+c0*v[nIndex])*vvs2_dtx_dtx;

					tempwx2=(tempwx2+c0*w[nIndex])*vvs2_dtx_dtx;

					tempuy2=(tempuy2+c0*u[nIndex])*vvs2_dtx_dtx;

					tempvy2=(tempvy2+c0*v[nIndex])*vvp2_dtx_dtx;

					tempwy2=(tempwy2+c0*w[nIndex])*vvs2_dtx_dtx;

					tempuz2=(tempuz2+c0*u[nIndex])*vvs2_dtz_dtz;

					tempvz2=(tempvz2+c0*v[nIndex])*vvs2_dtz_dtz;

					tempwz2=(tempwz2+c0*w[nIndex])*vvp2_dtz_dtz;

					for(kk=1;kk<=5;kk++)
					{
						for(kkk=1;kkk<=5;kkk++)
						{
							current_c = c[1+kk][kkk-1];

							_tempuxy = u[POSITION_INDEX_X(k,j+kkk,i+kk)];
							_tempvxy = v[POSITION_INDEX_X(k,j+kkk,i+kk)];

							_tempuxy -= u[POSITION_INDEX_X(k,j-kkk,i+kk)];
							_tempvxy -= v[POSITION_INDEX_X(k,j-kkk,i+kk)];

							_tempuxy += u[POSITION_INDEX_X(k,j-kkk,i-kk)];
							_tempvxy += v[POSITION_INDEX_X(k,j-kkk,i-kk)];

							_tempuxy -= u[POSITION_INDEX_X(k,j+kkk,i-kk)];
							_tempvxy -= v[POSITION_INDEX_X(k,j+kkk,i-kk)];

							_tempvyz = v[POSITION_INDEX_X(k+kkk,j+kk,i)];
							_tempwyz = w[POSITION_INDEX_X(k+kkk,j+kk,i)];

							_tempvyz -= v[POSITION_INDEX_X(k-kkk,j+kk,i)];
							_tempwyz -= w[POSITION_INDEX_X(k-kkk,j+kk,i)];

							_tempvyz += v[POSITION_INDEX_X(k-kkk,j-kk,i)];
							_tempwyz += w[POSITION_INDEX_X(k-kkk,j-kk,i)];

							_tempvyz -= v[POSITION_INDEX_X(k+kkk,j-kk,i)];
							_tempwyz -= w[POSITION_INDEX_X(k+kkk,j-kk,i)];

							_tempuxz = u[POSITION_INDEX_X(k+kkk,j,i+kk)];
							_tempwxz = w[POSITION_INDEX_X(k+kkk,j,i+kk)];

							_tempuxz -= u[POSITION_INDEX_X(k-kkk,j,i+kk)];
							_tempwxz -= w[POSITION_INDEX_X(k-kkk,j,i+kk)];

							_tempuxz += u[POSITION_INDEX_X(k-kkk,j,i-kk)];
							_tempwxz += w[POSITION_INDEX_X(k-kkk,j,i-kk)];

							_tempuxz -= u[POSITION_INDEX_X(k+kkk,j,i-kk)];
							_tempwxz -= w[POSITION_INDEX_X(k+kkk,j,i-kk)];

							tempuxz = tempuxz + (current_c*_tempuxz);
							tempwxz = tempwxz + (current_c*_tempwxz);

							tempvyz = tempvyz + (current_c*_tempvyz);
							tempwyz = tempwyz + (current_c*_tempwyz);

							tempuxy = tempuxy + (current_c*_tempuxy);
							tempvxy = tempvxy + (current_c*_tempvxy);

						} // for(kkk=1;kkk<=5;kkk++) end
					} //for(kk=1;kk<=5;kk++) end
					up[nIndex] = tempux2 + tempvxy * vvp2_dtz_dtx										+ tempwxz * vvp2_dtz_dtx;
					vp[nIndex] = tempuxy * vvp2_dtz_dtx				+ tempvy2 + tempwyz*vvp2_dtz_dtx 	;
					wp[nIndex] = px * wave[l-1]						+ tempvyz * vvp2_dtz_dtx			+ tempwz2 + tempuxz * vvp2_dtz_dtx;

					us[nIndex] = - tempvxy * vvs2_dtz_dtx			+ tempuy2							+ tempuz2 - tempwxz * vvs2_dtz_dtx;;
					vs[nIndex] = tempvx2 - tempuxy * vvs2_dtz_dtx   - tempwyz * vvs2_dtz_dtx			+ tempvz2;
					ws[nIndex] = tempwx2							+ tempwy2 - tempvyz*vvs2_dtz_dtx	- tempuxz * vvs2_dtz_dtx;
					// //Debug
				// if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
					// 	printf("[X]%lf wave:%lf px:%lf should be %lf\n", wp[nIndex],wave[l-1],px,wave[l-1]*px);
				}
			}
		}

        #ifndef DEBUG_NO_PARALLEL
        #pragma omp parallel for private(i,j,k)
        #endif
		for(k=k_begin;k<k_end;k++)
			for(j=j_begin;j<j_end;j++)
				for(i=i_begin;i<i_end;i++)
				{
					int nIndex              = POSITION_INDEX_X(k,j,i);

					up[nIndex] += 2 * up1[nIndex] - up2[nIndex];
					vp[nIndex] += 2 * vp1[nIndex] - vp2[nIndex];
					wp[nIndex] += 2 * wp1[nIndex] - wp2[nIndex];
					us[nIndex] += 2 * us1[nIndex] - us2[nIndex];
					vs[nIndex] += 2 * vs1[nIndex] - vs2[nIndex];
					ws[nIndex] += 2 * ws1[nIndex] - ws2[nIndex];
					// if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
					// 	printf("[Final]%lf\n", wp[nIndex]);
					u[nIndex] = up[nIndex] + us[nIndex];
					v[nIndex] = vp[nIndex] + vs[nIndex];
					w[nIndex] = wp[nIndex] + ws[nIndex];

					// if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
					// 	printf("[Final]%lf\n", w_x[nIndex]);
				}//for(i=nleft;i<nright;i++) end

		// printf("Start waiting....%d\n",l);

		double *swap_temp;
		swap_temp = up2; up2 = up1; up1 = up; up = swap_temp;
		swap_temp = vp2; vp2 = vp1; vp1 = vp; vp = swap_temp;
		swap_temp = wp2; wp2 = wp1; wp1 = wp; wp = swap_temp;
		swap_temp = us2; us2 = us1; us1 = us; us = swap_temp;
		swap_temp = vs2; vs2 = vs1; vs1 = vs; vs = swap_temp;
		swap_temp = ws2; ws2 = ws1; ws1 = ws; ws = swap_temp;
		// printf("[L]Finished %d\n",l);
	}//for(l=1;l<=lt;l++) end

	for(j=5;j<nMicYLength+5;j++)
		for(i=5;i<nMicXLength+5;i++)
		{
			up_out[POSITION_INDEX_X(0,j,i)] = up1[POSITION_INDEX_X(169-ntop+5,j,i)];
		}

    for(j=nfront;j<nback;j++)
        for(i=nleft;i<nright;i++){
        	to_write[POSITION_INDEX_HOST_X(0,j,i)] =up_out[POSITION_INDEX_X(0,j+5-nfront,i+5-nleft)];
        }
    free(up_out);
}



int calc_slice_on_mic();


int main(int argc, char **argv) {

	MEMORY_BLOCKS memory_blocks;

	if(argc<4)
	{
		printf("please add 3 parameter: inpurfile, outfile, logfile\n");
		exit(0);
	}

	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);
	strcpy(logfile,argv[3]);

	initailize();

	memory_blocks.u   = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.v   = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.w   = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.up  = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.up1 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.up2 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.vp  = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.vp1 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.vp2 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.wp  = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.wp1 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.wp2 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.us  = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.us1 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.us2 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.vs  = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.vs1 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.vs2 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.ws  = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.ws1 = (double *)malloc(sizeof(double)*mic_used_size);
	memory_blocks.ws2 = (double *)malloc(sizeof(double)*mic_used_size);

	memory_blocks.to_write = (double*)calloc(nSliceSize,sizeof(double));

	nshot=nxshot*nyshot;

	dtx=dt/unit;
	dtz=dt/unit;

	fout=fopen(outfile,"wb");

	// shot is divided to cluster, MPI
	for(ishot=1;ishot<=nshot;ishot++)
	{
		printf("shot=%d\n",ishot);
		flog = fopen(logfile,"a");
//		fprintf(flog,"shot=%d\n",ishot);
//		fclose(flog);

		ncy_shot=ncy_shot1+(ishot/nxshot)*dyshot;
		ncx_shot=ncx_shot1+(ishot%nxshot)*dxshot;

		calc_shot(
			ncx_shot,
			ncy_shot,
			1,lt,
			&memory_blocks
			);


		fwrite(to_write,sizeof(double),nSliceSize,fout);

	}//for(ishot=1;ishot<=nshot;ishot++) end

	fclose(fout);

	gettimeofday(&end,NULL);
	all_time = (end.tv_sec-start.tv_sec)+(double)(end.tv_usec-start.tv_usec)/1000000.0;
	printf("run time:\t%lf s\n",all_time);
	flog = fopen(logfile,"a");
	fprintf(flog,"\nrun time:\t%lf s\n\n",all_time);
	fclose(flog);
	flog = fopen(logfile,"a");
	fprintf(flog,"------------end time------------\n");
	fclose(flog);
	system(tmp);
	return 1;
}
