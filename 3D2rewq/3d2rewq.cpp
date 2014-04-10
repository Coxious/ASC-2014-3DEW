#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#include <omp.h>

#define PIE 3.1415926

#define DEBUG_NO_PARALLEL

#define DEBUG_CPU_RUNNING

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

#define _DEBUG_LEVEL_1_

#ifdef _DEBUG_LEVEL_1_

#define _MAX_DIFF_

#endif

// #define MIC_MEMORY_SIZE             (8*1024*1024*1024)
// #define MEMORY_REQUIRED_EACH_POint  (27 * sizeof(double))

// typedef struct _BLOCK_ARRANGE{
//     int nx;
//     int ny;
//     int nz;

//     int nxBlocksNormal;
//     int nyBlocksNormal;
//     int nzBlocksNormal;

//     int nEachNormal;

//     int nxEachEdgeEnd;
//     int nyEachEdgeEnd;
//     int nzEachEdgeEnd;

// }BLOCK_ARRANGE,*PBLOCK_ARRANGE;


// void CalcBlockSize(PBLOCK_ARRANGE pBlockArrange){

//     int nMemoryOffloadEachLimit = MIC_MEMORY_SIZE/2;

//     int nMemoryOffloadEachPointLimit
//             = nMemoryOffloadEachLimit/MEMORY_REQUIRED_EACH_POint;

//     int nPointLimitEachDirection
//         = pow(nMemoryOffloadEachPointLimit,1./3);

//     pBlockArrange -> nEachNormal
//         = ROUND_TO_SIZE_LESS(nPointLimitEachDirection,16);

//     //
//     // Here we program the seperation arrangement staticly for now.
//     //
//     pBlockArrange->nxBlocksNormal= pBlockArrange->nx / pBlockArrange->nEachNormal;
//     pBlockArrange->nyBlocksNormal= pBlockArrange->ny / pBlockArrange->nEachNormal;
//     pBlockArrange->nzBlocksNormal= pBlockArrange->nz / pBlockArrange->nEachNormal;

//     pBlockArrange->nxEachEdgeEnd = pBlockArrange->nx - pBlockArrange->nxBlocksNormal * pBlockArrange->nEachNormal;
//     pBlockArrange->nxEachEdgeEnd = pBlockArrange->ny - pBlockArrange->nyBlocksNormal * pBlockArrange->nEachNormal;
//     pBlockArrange->nxEachEdgeEnd = pBlockArrange->nz - pBlockArrange->nzBlocksNormal * pBlockArrange->nEachNormal;


// }

typedef struct _MEMORY_BLOCKS{

 double * up, * up1, * up2, * vp, * vp1, * vp2, * wp, * wp1, \
	  * wp2, * us, * us1, * us2, * vs, * vs1, * vs2, * ws, * ws1, * ws2;
 double *u_x, *u_y, *u_z, *v_x, *v_y, *v_z, *w_x, *w_y, *w_z;

}MEMORY_BLOCKS,*PMEMORY_BLOCKS;

// MIC_VAR double * up, * up1, * up2, * vp, * vp1, * vp2, * wp, * wp1, \
// 	  * wp2, * us, * us1, * us2, * vs, * vs1, * vs2, * ws, * ws1, * ws2;
// MIC_VAR double *u_x, *u_y, *u_z, *v_x, *v_y, *v_z, *w_x, *w_y, *w_z;


MIC_VAR int i,j,k,kk,kkk,l;
MIC_VAR int ishot,ncy_shot,ncx_shot;

MIC_VAR    double *wave;
MIC_VAR    double nshot,t0,tt,c0;
MIC_VAR    double dtx,dtz;
MIC_VAR    double px;

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


	to_write = (double*)calloc(nSliceSize,sizeof(double));

	mic_used_size = pow(2.*((lt*dt*velmax)/unit+10.)+10.,3.);
	mic_slice_size = pow(2.*((lt*dt*velmax)/unit+10.)+10.,2.);

	printf("length: %lfmic_slice_size:%d mic_used_size:%d\n",2.*((lt*dt*velmax)/unit+10.)+10.,mic_slice_size, mic_used_size);

	*up_out = (double *)malloc(mic_slice_size * sizeof(double));

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

    double vvp2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
					tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;

    double _tempuxz,_tempuxy,_tempvyz,_tempvxy,_tempwxz,_tempwyz;

	int n_mic_left,n_mic_right,n_mic_front,n_mic_back,n_mic_top,n_mic_bottom;
	int nMicXLength,nMicYLength,nMicZLength;
	int nMicMaxXLength,nMicMaxYLength,nMicMaxZLength;

	double * up  = pMemBlocks->up,
	double * up1 = pMemBlocks->up1,
	double * up2 = pMemBlocks->up2,
	double * vp, = pMemBlocks->vp,,
	double * vp1 = pMemBlocks->vp1,
	double * vp2 = pMemBlocks->vp2,
	double * wp, = pMemBlocks->wp,,
	double * wp1 = pMemBlocks->wp1,
	double * wp2 = pMemBlocks->wp2,
	double * us, = pMemBlocks->us,,
	double * us1 = pMemBlocks->us1,
	double * us2 = pMemBlocks->us2,
	double * vs, = pMemBlocks->vs,,
	double * vs1 = pMemBlocks->vs1,
	double * vs2 = pMemBlocks->vs2,
	double * ws, = pMemBlocks->ws,,
	double * ws1 = pMemBlocks->ws1,
	double * ws2 = pMemBlocks->ws2,
	double * u_x = pMemBlocks->u_x,
	double * u_y = pMemBlocks->u_y,
	double * u_z = pMemBlocks->u_z,
	double * v_x = pMemBlocks->v_x,
	double * v_y = pMemBlocks->v_y,
	double * v_z = pMemBlocks->v_z,
	double * w_x = pMemBlocks->w_x,
	double * w_y = pMemBlocks->w_y,
	double * w_z = pMemBlocks->w_z;

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

#ifndef DEBUG_CPU_RUNNING

	#pragma offload target(mic:0) \
		nocopy(u_x,v_x,w_x,u_y,v_y,w_y,u_z,v_z,w_z,\
			up ,up1,up2,vp ,vp1,vp2,wp ,wp1,wp2,us ,us1,us2,vs,\
			vs1,vs2,ws ,ws1,ws2) \
		nocopy(wave) \
		nocopy(c) wait(c)
	{
		// printf("MIC started!\n");
#else
	{
		// printf("CPU started!\n");
#endif
        // omp_set_num_threads(200);

		memset(u_x,0,sizeof(double)*mic_used_size);
		memset(v_x,0,sizeof(double)*mic_used_size);
		memset(w_x,0,sizeof(double)*mic_used_size);
		memset(u_y,0,sizeof(double)*mic_used_size);
		memset(v_y,0,sizeof(double)*mic_used_size);
		memset(w_y,0,sizeof(double)*mic_used_size);
		memset(u_z,0,sizeof(double)*mic_used_size);
		memset(v_z,0,sizeof(double)*mic_used_size);
		memset(w_z,0,sizeof(double)*mic_used_size);
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

		for(l=lStart;l<=lEnd;l++)
		{
			printf("l %d started\n",l);
			xmax=l*dt*velmax;
			n_mic_left=ncx_shot-xmax/unit-10;
			n_mic_right=ncx_shot+xmax/unit+10;
			n_mic_front=ncy_shot-xmax/unit-10;
			n_mic_back=ncy_shot+xmax/unit+10;
			n_mic_top=ncz_shot-xmax/unit-10;
			n_mic_bottom=ncz_shot+xmax/unit+10;

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
			// printf("Len_X: %d, Len_Y: %d, Len_Z: %d\n", nMicXLength, nMicYLength, nMicZLength );
	 	// 	printf("Size to use now : %d\n",(nMicXLength + 10)*(nMicZLength + 10)*(nMicYLength + 10));

			// printf("nX: %d nY: %d nZ: %d \n",nright - nleft,nback - nfront,nbottom - ntop);

			// ZYX
            #ifndef DEBUG_NO_PARALLEL
            #pragma omp parallel for private(i,j,k)
            #endif
			for(k=5+n_mic_top - ntop;k<n_mic_top - ntop + nMicZLength+5;k++)
			{
				vvp2=vpp2(k);
				vvs2=vss2(k);

				vvs2_dtz_dtz = vvs2*dtz*dtz;
				vvp2_dtx_dtx = vvp2*dtx*dtx;
				vvs2_dtx_dtx = vvs2*dtx*dtx;
				vvp2_dtz_dtz = vvp2*dtz*dtz;
				vvp2_dtz_dtx = vvp2*dtz*dtx;
				vvs2_dtz_dtx = vvs2*dtz*dtx;

				for(j=5+n_mic_front-nfront;j<n_mic_front-nfront+nMicYLength+5;j++) {
					for(i=5+n_mic_left-nleft;i<n_mic_left-nleft+nMicXLength+5;i++) {

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

 //#ifdef DEBUG_CPU_RUNNING
						for(kk=1;kk<=5;kk++)
						{
							tempux2=tempux2+c[0][kk-1]*(u_x[POSITION_INDEX_X(k,j,i+kk)]+u_x[POSITION_INDEX_X(k,j,i-kk)]);

							tempvx2=tempvx2+c[0][kk-1]*(v_x[POSITION_INDEX_X(k,j,i+kk)]+v_x[POSITION_INDEX_X(k,j,i-kk)]);

							tempwx2=tempwx2+c[0][kk-1]*(w_x[POSITION_INDEX_X(k,j,i+kk)]+w_x[POSITION_INDEX_X(k,j,i-kk)]);
						} //for(kk=1;kk<=5;kk++) end
//#else

//#endif
						tempux2=(tempux2+c0*u_x[POSITION_INDEX_X(k,j,i)])*vvp2_dtx_dtx;

						tempvx2=(tempvx2+c0*v_x[POSITION_INDEX_X(k,j,i)])*vvs2_dtx_dtx;

						tempwx2=(tempwx2+c0*w_x[POSITION_INDEX_X(k,j,i)])*vvs2_dtx_dtx;

						for(kk=1;kk<=5;kk++)
						{
							for(kkk=1;kkk<=5;kkk++)
							{
								current_c = c[1+kk][kkk-1];

								_tempuxy = u_x[POSITION_INDEX_X(k,j+kkk,i+kk)];
								_tempvxy = v_x[POSITION_INDEX_X(k,j+kkk,i+kk)];

								_tempuxy -= u_x[POSITION_INDEX_X(k,j-kkk,i+kk)];
								_tempvxy -= v_x[POSITION_INDEX_X(k,j-kkk,i+kk)];

								_tempuxy += u_x[POSITION_INDEX_X(k,j-kkk,i-kk)];
								_tempvxy += v_x[POSITION_INDEX_X(k,j-kkk,i-kk)];

								_tempuxy -= u_x[POSITION_INDEX_X(k,j+kkk,i-kk)];
								_tempvxy -= v_x[POSITION_INDEX_X(k,j+kkk,i-kk)];

								tempuxy = tempuxy + (current_c*_tempuxy);
								tempvxy = tempvxy + (current_c*_tempvxy);

							} // for(kkk=1;kkk<=5;kkk++) end
						} //for(kk=1;kk<=5;kk++) end
						up[nIndex] = tempux2 + tempvxy * vvp2_dtz_dtx;
						vp[nIndex] = tempuxy * vvp2_dtz_dtx;
						us[nIndex] = - tempvxy * vvs2_dtz_dtx;
						vs[nIndex] = tempvx2 - tempuxy * vvs2_dtz_dtx;
						ws[nIndex] = tempwx2;
						wp[nIndex] = px * wave[l-1];

						// //Debug
					// if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
						// 	printf("[X]%lf wave:%lf px:%lf should be %lf\n", wp[nIndex],wave[l-1],px,wave[l-1]*px);

					}
				}
			}

			// X Z Y
            #ifndef DEBUG_NO_PARALLEL
            #pragma omp parallel for private(i,j,k)
            #endif
			for(i=5+n_mic_left-nleft;i<n_mic_left-nleft+nMicXLength+5;i++) {
				for(k=5+n_mic_top - ntop;k<nMicZLength+n_mic_top - ntop+5;k++) {
					vvp2=vpp2(k);
					vvs2=vss2(k);

					vvs2_dtz_dtz = vvs2*dtz*dtz;
					vvp2_dtx_dtx = vvp2*dtx*dtx;
					vvs2_dtx_dtx = vvs2*dtx*dtx;
					vvp2_dtz_dtz = vvp2*dtz*dtz;
					vvp2_dtz_dtx = vvp2*dtz*dtx;
					vvs2_dtz_dtx = vvs2*dtz*dtx;

					for(j=5+n_mic_front-nfront;j<n_mic_front-nfront+nMicYLength+5;j++) {

						int nIndex = POSITION_INDEX_Y(k,j,i);
						int nIndex_X = POSITION_INDEX_X(k,j,i);

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
							//if ( u_x[POSITION_INDEX_X(k,j+kk,i) != u_y[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
							tempuy2=tempuy2+c[0][kk-1]*(u_y[POSITION_INDEX_Y(k,j+kk,i)]+u_y[POSITION_INDEX_Y(k,j-kk,i)]);

							//if ( v_x[POSITION_INDEX_X(k,j+kk,i) != v_y[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
							tempvy2=tempvy2+c[0][kk-1]*(v_y[POSITION_INDEX_Y(k,j+kk,i)]+v_y[POSITION_INDEX_Y(k,j-kk,i)]);

							//if ( w_x[POSITION_INDEX_X(k,j+kk,i) != w_y[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
							tempwy2=tempwy2+c[0][kk-1]*(w_y[POSITION_INDEX_Y(k,j+kk,i)]+w_y[POSITION_INDEX_Y(k,j-kk,i)]);
						} //for(kk=1;kk<=5;kk++) end

						tempuy2=(tempuy2+c0*u_y[POSITION_INDEX_Y(k,j,i)])*vvs2_dtx_dtx;

						tempvy2=(tempvy2+c0*v_y[POSITION_INDEX_Y(k,j,i)])*vvp2_dtx_dtx;

						tempwy2=(tempwy2+c0*w_y[POSITION_INDEX_Y(k,j,i)])*vvs2_dtx_dtx;

						for(kk=1;kk<=5;kk++)
						{
							for(kkk=1;kkk<=5;kkk++)
							{
								current_c = c[1+kk][kkk-1];

								_tempvyz = v_y[POSITION_INDEX_Y(k+kkk,j+kk,i)];
								_tempwyz = w_y[POSITION_INDEX_Y(k+kkk,j+kk,i)];

								_tempvyz -= v_y[POSITION_INDEX_Y(k-kkk,j+kk,i)];
								_tempwyz -= w_y[POSITION_INDEX_Y(k-kkk,j+kk,i)];

								_tempvyz += v_y[POSITION_INDEX_Y(k-kkk,j-kk,i)];
								_tempwyz += w_y[POSITION_INDEX_Y(k-kkk,j-kk,i)];

								_tempvyz -= v_y[POSITION_INDEX_Y(k+kkk,j-kk,i)];
								_tempwyz -= w_y[POSITION_INDEX_Y(k+kkk,j-kk,i)];

								tempvyz = tempvyz + (current_c*_tempvyz);
								tempwyz = tempwyz + (current_c*_tempwyz);

							} // for(kkk=1;kkk<=5;kkk++) end
						} //for(kk=1;kk<=5;kk++) end

						up[nIndex_X] += 0;
						vp[nIndex_X] += tempvy2 + tempwyz*vvp2_dtz_dtx;
						wp[nIndex_X] += tempvyz * vvp2_dtz_dtx;
						us[nIndex_X] += tempuy2;
						vs[nIndex_X] += - tempwyz * vvs2_dtz_dtx;
						ws[nIndex_X] += tempwy2 - tempvyz*vvs2_dtz_dtx;

					}
				}
			}

			// YXZ
            #ifndef DEBUG_NO_PARALLEL
            #pragma omp parallel for private(i,j,k)
            #endif
			for(j=5+n_mic_front-nfront;j<n_mic_front-nfront+nMicYLength+5;j++) {
				for(i=5+n_mic_left-nleft;i<n_mic_left-nleft+nMicXLength+5;i++) {
					for(k=5+n_mic_top - ntop;k<n_mic_top - ntop+nMicZLength+5;k++)
					{
						vvp2=vpp2(k);
						vvs2=vss2(k);

						vvs2_dtz_dtz = vvs2*dtz*dtz;
						vvp2_dtx_dtx = vvp2*dtx*dtx;
						vvs2_dtx_dtx = vvs2*dtx*dtx;
						vvp2_dtz_dtz = vvp2*dtz*dtz;
						vvp2_dtz_dtx = vvp2*dtz*dtx;
						vvs2_dtz_dtx = vvs2*dtz*dtx;

						int nIndex = POSITION_INDEX_Z(k,j,i);
						int nIndex_X = POSITION_INDEX_X(k,j,i);

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
							//if ( u_x[POSITION_INDEX_X(k,j+kk,i) != u_z[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
							tempuz2=tempuz2+c[0][kk-1]*(u_z[POSITION_INDEX_Z(k+kk,j,i)]+u_z[POSITION_INDEX_Z(k-kk,j,i)]);

							// if ( v_x[POSITION_INDEX_X(k,j+kk,i) != v_z[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
							tempvz2=tempvz2+c[0][kk-1]*(v_z[POSITION_INDEX_Z(k+kk,j,i)]+v_z[POSITION_INDEX_Z(k-kk,j,i)]);

							// if ( w_x[POSITION_INDEX_X(k,j+kk,i) != w_z[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
							tempwz2=tempwz2+c[0][kk-1]*(w_z[POSITION_INDEX_Z(k+kk,j,i)]+w_z[POSITION_INDEX_Z(k-kk,j,i)]);
						} //for(kk=1;kk<=5;kk++) end

						tempuz2=(tempuz2+c0*u_z[POSITION_INDEX_Z(k,j,i)])*vvs2_dtz_dtz;

						tempvz2=(tempvz2+c0*v_z[POSITION_INDEX_Z(k,j,i)])*vvs2_dtz_dtz;

						tempwz2=(tempwz2+c0*w_z[POSITION_INDEX_Z(k,j,i)])*vvp2_dtz_dtz;

						for(kk=1;kk<=5;kk++)
						{
							for(kkk=1;kkk<=5;kkk++)
							{
								current_c = c[1+kk][kkk-1];

								_tempuxz = u_z[POSITION_INDEX_Z(k+kkk,j,i+kk)];
								_tempwxz = w_z[POSITION_INDEX_Z(k+kkk,j,i+kk)];

								_tempuxz -= u_z[POSITION_INDEX_Z(k-kkk,j,i+kk)];
								_tempwxz -= w_z[POSITION_INDEX_Z(k-kkk,j,i+kk)];

								_tempuxz += u_z[POSITION_INDEX_Z(k-kkk,j,i-kk)];
								_tempwxz += w_z[POSITION_INDEX_Z(k-kkk,j,i-kk)];

								_tempuxz -= u_z[POSITION_INDEX_Z(k+kkk,j,i-kk)];
								_tempwxz -= w_z[POSITION_INDEX_Z(k+kkk,j,i-kk)];

								tempuxz = tempuxz + (current_c*_tempuxz);
								tempwxz = tempwxz + (current_c*_tempwxz);

							} // for(kkk=1;kkk<=5;kkk++) end
						} //for(kk=1;kk<=5;kk++) end

						up[nIndex_X] += tempwxz * vvp2_dtz_dtx;
						vp[nIndex_X] += 0;
						wp[nIndex_X] += tempwz2 + tempuxz * vvp2_dtz_dtx;
						us[nIndex_X] += tempuz2 - tempwxz * vvs2_dtz_dtx;
						vs[nIndex_X] += tempvz2;
						ws[nIndex_X] += - tempuxz * vvs2_dtz_dtx;
						// DEBUG
						// if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
						// 	printf("[Z]%lf\n", wp[nIndex]);
					}
				}
			}

            #ifndef DEBUG_NO_PARALLEL
            #pragma omp parallel for private(i,j,k)
            #endif
			for(k=5+n_mic_top - ntop;k<n_mic_top - ntop+nMicZLength+5;k++)
				for(j=5+n_mic_front-nfront;j<n_mic_front-nfront+nMicYLength+5;j++)
					for(i=5+n_mic_left-nleft;i<n_mic_left-nleft+nMicXLength+5;i++)
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
						u_x[POSITION_INDEX_X(k,j,i)] = up[POSITION_INDEX_X(k,j,i)] + us[POSITION_INDEX_X(k,j,i)];
						v_x[POSITION_INDEX_X(k,j,i)] = vp[POSITION_INDEX_X(k,j,i)] + vs[POSITION_INDEX_X(k,j,i)];
						w_x[POSITION_INDEX_X(k,j,i)] = wp[POSITION_INDEX_X(k,j,i)] + ws[POSITION_INDEX_X(k,j,i)];

						u_z[POSITION_INDEX_Z(k,j,i)] = u_y[POSITION_INDEX_Y(k,j,i)] = u_x[POSITION_INDEX_X(k,j,i)];
						v_z[POSITION_INDEX_Z(k,j,i)] = v_y[POSITION_INDEX_Y(k,j,i)] = v_x[POSITION_INDEX_X(k,j,i)];
						w_z[POSITION_INDEX_Z(k,j,i)] = w_y[POSITION_INDEX_Y(k,j,i)] = w_x[POSITION_INDEX_X(k,j,i)];
						// if(i-5+nleft==ncx_shot-1&&j-5+nfront==ncy_shot-1&&k-5+ntop==ncz_shot-1)
						// 	printf("[Final]%lf\n", w_x[nIndex]);
					}//for(i=nleft;i<nright;i++) end

			// printf("Start waiting....%d\n",l);
#ifndef DEBUG_CPU_RUNNING
#			pragma offload_wait target(mic:0) wait(up_out)
			{}
#endif
			double *swap_temp;
			swap_temp = up2; up2 = up1; up1 = up; up = swap_temp;
			swap_temp = vp2; vp2 = vp1; vp1 = vp; vp = swap_temp;
			swap_temp = wp2; wp2 = wp1; wp1 = wp; wp = swap_temp;
			swap_temp = us2; us2 = us1; us1 = us; us = swap_temp;
			swap_temp = vs2; vs2 = vs1; vs1 = vs; vs = swap_temp;
			swap_temp = ws2; ws2 = ws1; ws1 = ws; ws = swap_temp;
			// printf("[L]Finished %d\n",l);
		}//for(l=1;l<=lt;l++) end
	}//MIC END
}

int main(int argc, char **argv) {

	if(argc<4)
	{
		printf("please add 3 parameter: inpurfile, outfile, logfile\n");
		exit(0);
	}

	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);
	strcpy(logfile,argv[3]);

	initailize();

	double current_c;

#ifndef DEBUG_CPU_RUNNING
		u_x = (double *)malloc(sizeof(double)*mic_used_size);
		v_x = (double *)malloc(sizeof(double)*mic_used_size);
		w_x = (double *)malloc(sizeof(double)*mic_used_size);
		u_y = (double *)malloc(sizeof(double)*mic_used_size);
		v_y = (double *)malloc(sizeof(double)*mic_used_size);
		w_y = (double *)malloc(sizeof(double)*mic_used_size);
		u_z = (double *)malloc(sizeof(double)*mic_used_size);
		v_z = (double *)malloc(sizeof(double)*mic_used_size);
		w_z = (double *)malloc(sizeof(double)*mic_used_size);
		up  = (double *)malloc(sizeof(double)*mic_used_size);
		up1 = (double *)malloc(sizeof(double)*mic_used_size);
		up2 = (double *)malloc(sizeof(double)*mic_used_size);
		vp  = (double *)malloc(sizeof(double)*mic_used_size);
		vp1 = (double *)malloc(sizeof(double)*mic_used_size);
		vp2 = (double *)malloc(sizeof(double)*mic_used_size);
		wp  = (double *)malloc(sizeof(double)*mic_used_size);
		wp1 = (double *)malloc(sizeof(double)*mic_used_size);
		wp2 = (double *)malloc(sizeof(double)*mic_used_size);
		us  = (double *)malloc(sizeof(double)*mic_used_size);
		us1 = (double *)malloc(sizeof(double)*mic_used_size);
		us2 = (double *)malloc(sizeof(double)*mic_used_size);
		vs  = (double *)malloc(sizeof(double)*mic_used_size);
		vs1 = (double *)malloc(sizeof(double)*mic_used_size);
		vs2 = (double *)malloc(sizeof(double)*mic_used_size);
		ws  = (double *)malloc(sizeof(double)*mic_used_size);
		ws1 = (double *)malloc(sizeof(double)*mic_used_size);
		ws2 = (double *)malloc(sizeof(double)*mic_used_size);
	#pragma offload_transfer target(mic:0) \
        nocopy(u_x:length(mic_used_size) MIC_ALLOC) \
        nocopy(v_x:length(mic_used_size) MIC_ALLOC) \
        nocopy(w_x:length(mic_used_size) MIC_ALLOC) \
        nocopy(u_y:length(mic_used_size) MIC_ALLOC) \
        nocopy(v_y:length(mic_used_size) MIC_ALLOC) \
        nocopy(w_y:length(mic_used_size) MIC_ALLOC) \
        nocopy(u_z:length(mic_used_size) MIC_ALLOC) \
        nocopy(v_z:length(mic_used_size) MIC_ALLOC) \
        nocopy(w_z:length(mic_used_size) MIC_ALLOC) \
        nocopy(up :length(mic_used_size) MIC_ALLOC) \
        nocopy(up1:length(mic_used_size) MIC_ALLOC) \
        nocopy(up2:length(mic_used_size) MIC_ALLOC) \
        nocopy(vp :length(mic_used_size) MIC_ALLOC) \
        nocopy(vp1:length(mic_used_size) MIC_ALLOC) \
        nocopy(vp2:length(mic_used_size) MIC_ALLOC) \
        nocopy(wp :length(mic_used_size) MIC_ALLOC) \
        nocopy(wp1:length(mic_used_size) MIC_ALLOC) \
        nocopy(wp2:length(mic_used_size) MIC_ALLOC) \
        nocopy(us :length(mic_used_size) MIC_ALLOC) \
        nocopy(us1:length(mic_used_size) MIC_ALLOC) \
        nocopy(us2:length(mic_used_size) MIC_ALLOC) \
        nocopy(vs :length(mic_used_size) MIC_ALLOC) \
        nocopy(vs1:length(mic_used_size) MIC_ALLOC) \
        nocopy(vs2:length(mic_used_size) MIC_ALLOC) \
        nocopy(ws :length(mic_used_size) MIC_ALLOC) \
        nocopy(ws1:length(mic_used_size) MIC_ALLOC) \
        nocopy(ws2:length(mic_used_size) MIC_ALLOC) \
        nocopy(up_out:length(mic_slice_size) MIC_ALLOC)\
        in(wave:length(lt) MIC_ALLOC) \
        in(c: MIC_ALLOC) signal(c)
	{
		// printf("Transfer finished!\n");
	}

#else
		u_x = (double *)malloc(sizeof(double)*mic_used_size);
		v_x = (double *)malloc(sizeof(double)*mic_used_size);
		w_x = (double *)malloc(sizeof(double)*mic_used_size);
		u_y = (double *)malloc(sizeof(double)*mic_used_size);
		v_y = (double *)malloc(sizeof(double)*mic_used_size);
		w_y = (double *)malloc(sizeof(double)*mic_used_size);
		u_z = (double *)malloc(sizeof(double)*mic_used_size);
		v_z = (double *)malloc(sizeof(double)*mic_used_size);
		w_z = (double *)malloc(sizeof(double)*mic_used_size);
		up  = (double *)malloc(sizeof(double)*mic_used_size);
		up1 = (double *)malloc(sizeof(double)*mic_used_size);
		up2 = (double *)malloc(sizeof(double)*mic_used_size);
		vp  = (double *)malloc(sizeof(double)*mic_used_size);
		vp1 = (double *)malloc(sizeof(double)*mic_used_size);
		vp2 = (double *)malloc(sizeof(double)*mic_used_size);
		wp  = (double *)malloc(sizeof(double)*mic_used_size);
		wp1 = (double *)malloc(sizeof(double)*mic_used_size);
		wp2 = (double *)malloc(sizeof(double)*mic_used_size);
		us  = (double *)malloc(sizeof(double)*mic_used_size);
		us1 = (double *)malloc(sizeof(double)*mic_used_size);
		us2 = (double *)malloc(sizeof(double)*mic_used_size);
		vs  = (double *)malloc(sizeof(double)*mic_used_size);
		vs1 = (double *)malloc(sizeof(double)*mic_used_size);
		vs2 = (double *)malloc(sizeof(double)*mic_used_size);
		ws  = (double *)malloc(sizeof(double)*mic_used_size);
		ws1 = (double *)malloc(sizeof(double)*mic_used_size);
		ws2 = (double *)malloc(sizeof(double)*mic_used_size);
	#endif

	nshot=nxshot*nyshot;

	dtx=dt/unit;
	dtz=dt/unit;

	fout=fopen(outfile,"wb");

	// shot is divided to cluster, MPI
	for(ishot=1;ishot<=nshot;ishot++)
	{
		printf("shot=%d\n",ishot);
		flog = fopen(logfile,"a");
		fprintf(flog,"shot=%d\n",ishot);
		fclose(flog);

		ncy_shot=ncy_shot1+(ishot/nxshot)*dyshot;
		ncx_shot=ncx_shot1+(ishot%nxshot)*dxshot;

		calc_shot(
			ncx_shot,
			ncy_shot,
			1,lt,
			u_x,u_y,u_z,
			v_x,v_y,v_z,
			w_x,w_y,w_z
			);

#ifndef DEBUG_CPU_RUNNING
		#pragma offload target(mic:0)\
			nocopy(up1:MIC_REUSE)\
			out(up_out :length(mic_slice_size) MIC_REUSE) signal(up_out)
#endif
		{

		// printf("To reindex the data!");
		for(j=5;j<nMicYLength+5;j++)
			for(i=5;i<nMicXLength+5;i++)
			{
				up_out[POSITION_INDEX_X(0,j,i)] = up1[POSITION_INDEX_X(169-ntop+5,j,i)];
//				up_out[i]=((double*)(up1+(169-ntop+5)*mic_slice_size))[i];
				// printf("up_out i:%d j:%d %lf\n",i,j,up_out[POSITION_INDEX_X(0,j,i)]);
			}
		}

#ifndef DEBUG_CPU_RUNNING
		#pragma offload_wait target(mic:0) wait(up_out)
#endif
		// printf("To write the data!");
        for(j=nfront;j<nback;j++)
            for(i=nleft;i<nright;i++){
            	to_write[POSITION_INDEX_HOST_X(0,j,i)] =up_out[POSITION_INDEX_X(0,j+5-nfront,i+5-nleft)];
//		        fprintf(fout,"i:%d j:%d :%lf\n",i,j,to_write[POSITION_INDEX_HOST_X(0,j,i)]);
            }
		fwrite(to_write,sizeof(double),nSliceSize,fout);

	}//for(ishot=1;ishot<=nshot;ishot++) end

	fclose(fout);

#ifdef DEBUG_CPU_RUNNING
	free(u_x);
	free(v_x);
	free(w_x);
	free(u_y);
	free(v_y);
	free(w_y);
	free(u_z);
	free(v_z);
	free(w_z);
	free(up );
	free(up1);
	free(up2);
	free(vp );
	free(vp1);
	free(vp2);
	free(wp );
	free(wp1);
	free(wp2);
	free(us );
	free(us1);
	free(us2);
	free(vs );
	free(vs1);
	free(vs2);
	free(ws );
	free(ws1);
	free(ws2);
#endif

	free(up_out);
	free(wave);

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
