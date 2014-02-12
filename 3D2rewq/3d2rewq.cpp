#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#define PIE 3.1415926

#define POSITION_INDEX_X(_z,_y,_x)        ((_z)*ny*nx + (_y)*nx + (_x))
#define POSITION_INDEX_Y(_z,_y,_x)        ((_x)*nz*ny + (_z)*ny + (_y))
#define POSITION_INDEX_Z(_z,_y,_x)        ((_y)*nx*nz + (_x)*nz + (_z))

#define vpp2(_z) ((_z<210)?(5290000):(_z>=260?12250000:7840000))
#define vss2(_z) ((_z<210)?(1517824):(_z>=260?3644281:2277081))

//#define _DEBUG_TEST
#define _DEBUG_LEVEL_2

#ifdef _DEBUG_TEST
#define _DEBUG_MAX_DIFF	1e-1
#endif

#ifdef _DEBUG_LEVEL_2
#define _DEBUG_MAX_DIFF	1e-1
#endif



float *u_x, *v_x, *w_x, *u_y, *v_y, *w_y, *u_z, *v_z, *w_z, *up , *up1, *up2, *vp , *vp1, *vp2,
		 *wp , *wp1, *wp2, *us , *us1, *us2, *vs , *vs1, *vs2, *ws , *ws1, *ws2, *swap_temp;

#ifdef _DEBUG_TEST
float *debug_up,*debug_vp,*debug_wp,*debug_us,*debug_vs,*debug_ws;
#endif

#ifdef _DEBUG_LEVEL_2
    float * debug_2_u  ;
    float * debug_2_v  ;
    float * debug_2_w  ;
    float * debug_2_up ;
    float * debug_2_up1;
    float * debug_2_up2;
    float * debug_2_vp ;
    float * debug_2_vp1;
    float * debug_2_vp2;
    float * debug_2_wp ;
    float * debug_2_wp1;
    float * debug_2_wp2;
    float * debug_2_us ;
    float * debug_2_us1;
    float * debug_2_us2;
    float * debug_2_vs ;
    float * debug_2_vs1;
    float * debug_2_vs2;
    float * debug_2_ws ;
    float * debug_2_ws1;
    float * debug_2_ws2;
#endif

int main(int argc, char **argv)
{
	int i,j,k,kk,kkk,l,mm=5;
	int nx,ny,nz,lt,nedge;
	int nleft,nright,nfront,nback,ntop,nbottom;
	float frequency;
	float velmax;
	float dt;
	int ncx_shot1,ncy_shot1,ncz_shot;
	int ishot,ncy_shot,ncx_shot;
	float unit;
	int nxshot,nyshot,dxshot,dyshot;
	char infile[80],outfile[80],logfile[80],tmp[80];
	FILE  *fin, *fout, *flog;
	struct timeval start,end;
	float all_time;
	int px, sx;
	float c[5][7];
	float *wave;
	float nshot,t0,tt,c0;
	float dtx,dtz;
	float xmax;
	float vvp2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
		  tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;

	float _tempuxz,_tempuxy,_tempvyz,_tempvxy,_tempwxz,_tempwyz;

	if(argc<4)
	{
		printf("please add 3 parameter: inpurfile, outfile, logfile\n");
		exit(0);
	}

	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);
	strcpy(logfile,argv[3]);

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
	fscanf(fin,"frequency=%f\n",&frequency);
	fscanf(fin,"velmax=%f\n",&velmax);
	fscanf(fin,"dt=%f\n",&dt);
	fscanf(fin,"unit=%f\n",&unit);
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
	printf("frequency=%f\n",frequency);
	printf("velmax=%f\n",velmax);
	printf("dt=%f\n",dt);
	printf("unit=%f\n",unit);
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
	fprintf(flog,"frequency=%f\n",frequency);
	fprintf(flog,"velmax=%f\n",velmax);
	fprintf(flog,"dt=%f\n",dt);
	fprintf(flog,"unit=%f\n",unit);
	fprintf(flog,"dxshot=%d\n",dxshot);
	fprintf(flog,"dyshot=%d\n\n",dyshot);
	// fclose(flog);

	const int nSize = nz*ny*nx;
	const int nSliceSize = ny*nx;

	u_x = (float *) calloc(nSize,sizeof(float));
	v_x = (float *) calloc(nSize,sizeof(float));
	w_x = (float *) calloc(nSize,sizeof(float));
	u_y = (float *) calloc(nSize,sizeof(float));
	v_y = (float *) calloc(nSize,sizeof(float));
	w_y = (float *) calloc(nSize,sizeof(float));
	u_z = (float *) calloc(nSize,sizeof(float));
	v_z = (float *) calloc(nSize,sizeof(float));
	w_z = (float *) calloc(nSize,sizeof(float));
	up  = (float *) calloc(nSize,sizeof(float));
	up1 = (float *) calloc(nSize,sizeof(float));
	up2 = (float *) calloc(nSize,sizeof(float));
	vp  = (float *) calloc(nSize,sizeof(float));
	vp1 = (float *) calloc(nSize,sizeof(float));
	vp2 = (float *) calloc(nSize,sizeof(float));
	wp  = (float *) calloc(nSize,sizeof(float));
	wp1 = (float *) calloc(nSize,sizeof(float));
	wp2 = (float *) calloc(nSize,sizeof(float));
	us  = (float *) calloc(nSize,sizeof(float));
	us1 = (float *) calloc(nSize,sizeof(float));
	us2 = (float *) calloc(nSize,sizeof(float));
	vs  = (float *) calloc(nSize,sizeof(float));
	vs1 = (float *) calloc(nSize,sizeof(float));
	vs2 = (float *) calloc(nSize,sizeof(float));
	ws  = (float *) calloc(nSize,sizeof(float));
	ws1 = (float *) calloc(nSize,sizeof(float));
	ws2 = (float *) calloc(nSize,sizeof(float));

	wave    = (float*) calloc(lt,sizeof(float));

#ifdef _DEBUG_TEST
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
	debug_up= (float *) calloc(nSize,sizeof(float));
	debug_vp= (float *) calloc(nSize,sizeof(float));
	debug_wp= (float *) calloc(nSize,sizeof(float));
	debug_us= (float *) calloc(nSize,sizeof(float));
	debug_vs= (float *) calloc(nSize,sizeof(float));
	debug_ws= (float *) calloc(nSize,sizeof(float));
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif

#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug level 2 purpose code start

    debug_2_u       = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_v       = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_w       = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_up      = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_up1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_up2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_vp      = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_vp1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_vp2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_wp      = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_wp1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_wp2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_us      = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_us1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_us2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_vs      = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_vs1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_vs2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_ws      = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_ws1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    debug_2_ws2     = (float*)malloc(sizeof(float)*nz*ny*nx);
	//	Debug level 2 purpose code end
	/////////////////////////////////////////////////////////////
#endif

	nshot=nxshot*nyshot;
	t0=1.0/frequency;
	for(l=0;l<lt;l++)
	{
		tt=l*dt;
		tt=tt-t0;
		float sp=PIE*frequency*tt;
		float fx=100000.*exp(-sp*sp)*(1.-2.*sp*sp);
		wave[l]=fx;
	}

	if(mm==5)
	{
		c0=-2.927222164;
		c[0][0]=1.66666665;
		c[1][0]=-0.23809525;
		c[2][0]=0.03968254;
		c[3][0]=-0.004960318;
		c[4][0]=0.0003174603;
	}

	c[0][1]=0.83333;
	c[1][1]=-0.2381;
	c[2][1]=0.0595;
	c[3][1]=-0.0099;
	c[4][1]=0.0008;

	for(i=0;i<5;i++)
		for(j=0;j<5;j++)
			c[j][2+i]=c[i][1]*c[j][1];

	dtx=dt/unit;
	dtz=dt/unit;

	float vvp2_dtx_dtx;
	float vvs2_dtz_dtz;
	float vvs2_dtx_dtx;
	float vvp2_dtz_dtz;
	float vvp2_dtz_dtx;
	float vvs2_dtz_dtx;

	float current_c;
	long long nIndex, nIndex_X, nIndex_Y, nIndex_Z;

	fout=fopen(outfile,"wb");
	// shot is divided to cluster, MPI
	for(ishot=1;ishot<=nshot;ishot++)
	{
		printf("shot=%d\n",ishot);
		// flog = fopen(logfile,"a");
		fprintf(flog,"shot=%d\n",ishot);
		// fclose(flog);
		ncy_shot=ncy_shot1+(ishot/nxshot)*dyshot;
		ncx_shot=ncx_shot1+(ishot%nxshot)*dxshot;
//		printf("x: %d, y: %d, z: %d", ncx_shot, ncy_shot, ncz_shot);

		memset(u_x ,0,nSize*sizeof(float));
		memset(v_x ,0,nSize*sizeof(float));
		memset(w_x ,0,nSize*sizeof(float));
		memset(u_y ,0,nSize*sizeof(float));
		memset(v_y ,0,nSize*sizeof(float));
		memset(w_y ,0,nSize*sizeof(float));
		memset(u_z ,0,nSize*sizeof(float));
		memset(v_z ,0,nSize*sizeof(float));
		memset(w_z ,0,nSize*sizeof(float));
		memset(up  ,0,nSize*sizeof(float));
		memset(up1 ,0,nSize*sizeof(float));
		memset(up2 ,0,nSize*sizeof(float));
		memset(vp  ,0,nSize*sizeof(float));
		memset(vp1 ,0,nSize*sizeof(float));
		memset(vp2 ,0,nSize*sizeof(float));
		memset(wp  ,0,nSize*sizeof(float));
		memset(wp1 ,0,nSize*sizeof(float));
		memset(wp2 ,0,nSize*sizeof(float));
		memset(us  ,0,nSize*sizeof(float));
		memset(us1 ,0,nSize*sizeof(float));
		memset(us2 ,0,nSize*sizeof(float));
		memset(vs  ,0,nSize*sizeof(float));
		memset(vs1 ,0,nSize*sizeof(float));
		memset(vs2 ,0,nSize*sizeof(float));
		memset(ws  ,0,nSize*sizeof(float));
		memset(ws1 ,0,nSize*sizeof(float));
		memset(ws2 ,0,nSize*sizeof(float));

		for(l=1;l<=lt;l++)
		{
			xmax=l*dt*velmax;
			nleft=ncx_shot-xmax/unit-10;
			nright=ncx_shot+xmax/unit+10;
			nfront=ncy_shot-xmax/unit-10;
			nback=ncy_shot+xmax/unit+10;
			ntop=ncz_shot-xmax/unit-10;
			nbottom=ncz_shot+xmax/unit+10;
			if(nleft<5) nleft=5;
			if(nright>nx-5) nright=nx-5;
			if(nfront<5) nfront=5;
			if(nback>ny-5) nback=ny-5;
			if(ntop<5) ntop=5;
			if(nbottom>nz-5) nbottom=nz-5;

			ntop = ntop-1;
			nfront = nfront-1;
			nleft = nleft-1;
			// cout << ntop-nbottom << ' ' << nback-nfront << ' ' << nright-nleft << endl;

#ifdef _DEBUG_TEST
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
			for(k=ntop;k<nbottom;k++)
				for(j=nfront;j<nback;j++)
					for(i=nleft;i<nright;i++)
					{
						if(i==ncx_shot-1&&j==ncy_shot-1&&k==ncz_shot-1)
						{
							px=1.;
							sx=0.;
						}
						else
						{
							px=0.;
							sx=0.;
						}

						vvp2=vpp2(k);
						vvs2=vss2(k);


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
						for(kk=1;kk<=mm;kk++)
						{
							tempux2=tempux2+c[kk-1][0]*(u_x[k*ny*nx+j*nx+(i+kk)]+u_x[k*ny*nx+j*nx+(i-kk)]);
							tempuy2=tempuy2+c[kk-1][0]*(u_x[k*ny*nx+(j+kk)*nx+i]+u_x[k*ny*nx+(j-kk)*nx+i]);
							tempuz2=tempuz2+c[kk-1][0]*(u_x[(k+kk)*ny*nx+j*nx+i]+u_x[(k-kk)*ny*nx+j*nx+i]);

							tempvx2=tempvx2+c[kk-1][0]*(v_x[k*ny*nx+j*nx+(i+kk)]+v_x[k*ny*nx+j*nx+(i-kk)]);
							tempvy2=tempvy2+c[kk-1][0]*(v_x[k*ny*nx+(j+kk)*nx+i]+v_x[k*ny*nx+(j-kk)*nx+i]);
							tempvz2=tempvz2+c[kk-1][0]*(v_x[(k+kk)*ny*nx+j*nx+i]+v_x[(k-kk)*ny*nx+j*nx+i]);

							tempwx2=tempwx2+c[kk-1][0]*(w_x[k*ny*nx+j*nx+(i+kk)]+w_x[k*ny*nx+j*nx+(i-kk)]);
							tempwy2=tempwy2+c[kk-1][0]*(w_x[k*ny*nx+(j+kk)*nx+i]+w_x[k*ny*nx+(j-kk)*nx+i]);
							tempwz2=tempwz2+c[kk-1][0]*(w_x[(k+kk)*ny*nx+j*nx+i]+w_x[(k-kk)*ny*nx+j*nx+i]);

						} //for(kk=1;kk<=mm;kk++) end

						tempux2=(tempux2+c0*u_x[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
						tempuy2=(tempuy2+c0*u_x[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempuz2=(tempuz2+c0*u_x[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;

						tempvx2=(tempvx2+c0*v_x[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempvy2=(tempvy2+c0*v_x[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
						tempvz2=(tempvz2+c0*v_x[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;

						tempwx2=(tempwx2+c0*w_x[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempwy2=(tempwy2+c0*w_x[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempwz2=(tempwz2+c0*w_x[k*ny*nx+j*nx+i])*vvp2*dtz*dtz;

						for(kk=1;kk<=mm;kk++)
						{
							for(kkk=1;kkk<=mm;kkk++)
							{
								tempuxz=tempuxz+c[kkk-1][1+kk]*(u_x[(k+kkk)*ny*nx+j*nx+(i+kk)] -u_x[(k-kkk)*ny*nx+j*nx+(i+kk)] +u_x[(k-kkk)*ny*nx+j*nx+(i-kk)] -u_x[(k+kkk)*ny*nx+j*nx+(i-kk)]);
								tempuxy=tempuxy+c[kkk-1][1+kk]*(u_x[k*ny*nx+(j+kkk)*nx+(i+kk)] -u_x[k*ny*nx+(j-kkk)*nx+(i+kk)] +u_x[k*ny*nx+(j-kkk)*nx+(i-kk)] -u_x[k*ny*nx+(j+kkk)*nx+(i-kk)]);

								tempvyz=tempvyz+c[kkk-1][1+kk]*(v_x[(k+kkk)*ny*nx+(j+kk)*nx+i] -v_x[(k-kkk)*ny*nx+(j+kk)*nx+i] +v_x[(k-kkk)*ny*nx+(j-kk)*nx+i] -v_x[(k+kkk)*ny*nx+(j-kk)*nx+i]);
								tempvxy=tempvxy+c[kkk-1][1+kk]*(v_x[k*ny*nx+(j+kkk)*nx+(i+kk)] -v_x[k*ny*nx+(j-kkk)*nx+(i+kk)] +v_x[k*ny*nx+(j-kkk)*nx+(i-kk)] -v_x[k*ny*nx+(j+kkk)*nx+(i-kk)]);

								tempwyz=tempwyz+c[kkk-1][1+kk]*(w_x[(k+kkk)*ny*nx+(j+kk)*nx+i] -w_x[(k-kkk)*ny*nx+(j+kk)*nx+i] +w_x[(k-kkk)*ny*nx+(j-kk)*nx+i] -w_x[(k+kkk)*ny*nx+(j-kk)*nx+i]);
								tempwxz=tempwxz+c[kkk-1][1+kk]*(w_x[(k+kkk)*ny*nx+j*nx+(i+kk)] -w_x[(k-kkk)*ny*nx+j*nx+(i+kk)] +w_x[(k-kkk)*ny*nx+j*nx+(i-kk)] -w_x[(k+kkk)*ny*nx+j*nx+(i-kk)]);
							} // for(kkk=1;kkk<=mm;kkk++) end
						} //for(kk=1;kk<=mm;kk++) end
						debug_up[k*ny*nx+j*nx+i]=2.*up1[k*ny*nx+j*nx+i]-up2[k*ny*nx+j*nx+i]
										  +tempux2+tempwxz*vvp2*dtz*dtx
										  +tempvxy*vvp2*dtz*dtx;
						debug_vp[k*ny*nx+j*nx+i]=2.*vp1[k*ny*nx+j*nx+i]-vp2[k*ny*nx+j*nx+i]
										  +tempvy2+tempuxy*vvp2*dtz*dtx
										  +tempwyz*vvp2*dtz*dtx;
						debug_wp[k*ny*nx+j*nx+i]=2.*wp1[k*ny*nx+j*nx+i]-wp2[k*ny*nx+j*nx+i]
										  +tempwz2+tempuxz*vvp2*dtz*dtx
										  +tempvyz*vvp2*dtz*dtx
										  +px*wave[l-1];
						debug_us[k*ny*nx+j*nx+i]=2.*us1[k*ny*nx+j*nx+i]-us2[k*ny*nx+j*nx+i]+tempuy2+tempuz2
										  -tempvxy*vvs2*dtz*dtx-tempwxz*vvs2*dtz*dtx;
						debug_vs[k*ny*nx+j*nx+i]=2.*vs1[k*ny*nx+j*nx+i]-vs2[k*ny*nx+j*nx+i]+tempvx2+tempvz2
										  -tempuxy*vvs2*dtz*dtx-tempwyz*vvs2*dtz*dtx;
						debug_ws[k*ny*nx+j*nx+i]=2.*ws1[k*ny*nx+j*nx+i]-ws2[k*ny*nx+j*nx+i]+tempwx2+tempwy2
										  -tempuxz*vvs2*dtz*dtx-tempvyz*vvs2*dtz*dtx;
					}//for(i=nleft;i<nright;i++) end
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif


#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug level 2 purpose code start

			for(k=ntop;k<nbottom;k++)
				for(j=nfront;j<nback;j++)
					for(i=nleft;i<nright;i++)
					{
						if(i==ncx_shot-1&&j==ncy_shot-1&&k==ncz_shot-1)
						{
							px=1.;
							sx=0.;
						}
						else
						{
							px=0.;
							sx=0.;
						}

						vvp2=vpp2(k);
						vvs2=vss2(k);


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
						for(kk=1;kk<=mm;kk++)
						{
							tempux2=tempux2+c[kk-1][0]*(debug_2_u[k*ny*nx+j*nx+(i+kk)]+debug_2_u[k*ny*nx+j*nx+(i-kk)]);
							tempuy2=tempuy2+c[kk-1][0]*(debug_2_u[k*ny*nx+(j+kk)*nx+i]+debug_2_u[k*ny*nx+(j-kk)*nx+i]);
							tempuz2=tempuz2+c[kk-1][0]*(debug_2_u[(k+kk)*ny*nx+j*nx+i]+debug_2_u[(k-kk)*ny*nx+j*nx+i]);

							tempvx2=tempvx2+c[kk-1][0]*(debug_2_v[k*ny*nx+j*nx+(i+kk)]+debug_2_v[k*ny*nx+j*nx+(i-kk)]);
							tempvy2=tempvy2+c[kk-1][0]*(debug_2_v[k*ny*nx+(j+kk)*nx+i]+debug_2_v[k*ny*nx+(j-kk)*nx+i]);
							tempvz2=tempvz2+c[kk-1][0]*(debug_2_v[(k+kk)*ny*nx+j*nx+i]+debug_2_v[(k-kk)*ny*nx+j*nx+i]);

							tempwx2=tempwx2+c[kk-1][0]*(debug_2_w[k*ny*nx+j*nx+(i+kk)]+debug_2_w[k*ny*nx+j*nx+(i-kk)]);
							tempwy2=tempwy2+c[kk-1][0]*(debug_2_w[k*ny*nx+(j+kk)*nx+i]+debug_2_w[k*ny*nx+(j-kk)*nx+i]);
							tempwz2=tempwz2+c[kk-1][0]*(debug_2_w[(k+kk)*ny*nx+j*nx+i]+debug_2_w[(k-kk)*ny*nx+j*nx+i]);

						} //for(kk=1;kk<=mm;kk++) end

						tempux2=(tempux2+c0*debug_2_u[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
						tempuy2=(tempuy2+c0*debug_2_u[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempuz2=(tempuz2+c0*debug_2_u[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;

						tempvx2=(tempvx2+c0*debug_2_v[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempvy2=(tempvy2+c0*debug_2_v[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
						tempvz2=(tempvz2+c0*debug_2_v[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;

						tempwx2=(tempwx2+c0*debug_2_w[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempwy2=(tempwy2+c0*debug_2_w[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						tempwz2=(tempwz2+c0*debug_2_w[k*ny*nx+j*nx+i])*vvp2*dtz*dtz;

						for(kk=1;kk<=mm;kk++)
						{
							for(kkk=1;kkk<=mm;kkk++)
							{
								tempuxz=tempuxz+c[kkk-1][1+kk]*(debug_2_u[(k+kkk)*ny*nx+j*nx+(i+kk)] -debug_2_u[(k-kkk)*ny*nx+j*nx+(i+kk)] +debug_2_u[(k-kkk)*ny*nx+j*nx+(i-kk)] -debug_2_u[(k+kkk)*ny*nx+j*nx+(i-kk)]);
								tempuxy=tempuxy+c[kkk-1][1+kk]*(debug_2_u[k*ny*nx+(j+kkk)*nx+(i+kk)] -debug_2_u[k*ny*nx+(j-kkk)*nx+(i+kk)] +debug_2_u[k*ny*nx+(j-kkk)*nx+(i-kk)] -debug_2_u[k*ny*nx+(j+kkk)*nx+(i-kk)]);

								tempvyz=tempvyz+c[kkk-1][1+kk]*(debug_2_v[(k+kkk)*ny*nx+(j+kk)*nx+i] -debug_2_v[(k-kkk)*ny*nx+(j+kk)*nx+i] +debug_2_v[(k-kkk)*ny*nx+(j-kk)*nx+i] -debug_2_v[(k+kkk)*ny*nx+(j-kk)*nx+i]);
								tempvxy=tempvxy+c[kkk-1][1+kk]*(debug_2_v[k*ny*nx+(j+kkk)*nx+(i+kk)] -debug_2_v[k*ny*nx+(j-kkk)*nx+(i+kk)] +debug_2_v[k*ny*nx+(j-kkk)*nx+(i-kk)] -debug_2_v[k*ny*nx+(j+kkk)*nx+(i-kk)]);

								tempwyz=tempwyz+c[kkk-1][1+kk]*(debug_2_w[(k+kkk)*ny*nx+(j+kk)*nx+i] -debug_2_w[(k-kkk)*ny*nx+(j+kk)*nx+i] +debug_2_w[(k-kkk)*ny*nx+(j-kk)*nx+i] -debug_2_w[(k+kkk)*ny*nx+(j-kk)*nx+i]);
								tempwxz=tempwxz+c[kkk-1][1+kk]*(debug_2_w[(k+kkk)*ny*nx+j*nx+(i+kk)] -debug_2_w[(k-kkk)*ny*nx+j*nx+(i+kk)] +debug_2_w[(k-kkk)*ny*nx+j*nx+(i-kk)] -debug_2_w[(k+kkk)*ny*nx+j*nx+(i-kk)]);
							} // for(kkk=1;kkk<=mm;kkk++) end
						} //for(kk=1;kk<=mm;kk++) end
						debug_2_up[k*ny*nx+j*nx+i]=2.*debug_2_up1[k*ny*nx+j*nx+i]-debug_2_up2[k*ny*nx+j*nx+i]
										  +tempux2+tempwxz*vvp2*dtz*dtx
										  +tempvxy*vvp2*dtz*dtx;
						debug_2_vp[k*ny*nx+j*nx+i]=2.*debug_2_vp1[k*ny*nx+j*nx+i]-debug_2_vp2[k*ny*nx+j*nx+i]
										  +tempvy2+tempuxy*vvp2*dtz*dtx
										  +tempwyz*vvp2*dtz*dtx;
						debug_2_wp[k*ny*nx+j*nx+i]=2.*debug_2_wp1[k*ny*nx+j*nx+i]-debug_2_wp2[k*ny*nx+j*nx+i]
										  +tempwz2+tempuxz*vvp2*dtz*dtx
										  +tempvyz*vvp2*dtz*dtx
										  +px*wave[l-1];
						debug_2_us[k*ny*nx+j*nx+i]=2.*debug_2_us1[k*ny*nx+j*nx+i]-debug_2_us2[k*ny*nx+j*nx+i]+tempuy2+tempuz2
										  -tempvxy*vvs2*dtz*dtx-tempwxz*vvs2*dtz*dtx;
						debug_2_vs[k*ny*nx+j*nx+i]=2.*debug_2_vs1[k*ny*nx+j*nx+i]-debug_2_vs2[k*ny*nx+j*nx+i]+tempvx2+tempvz2
										  -tempuxy*vvs2*dtz*dtx-tempwyz*vvs2*dtz*dtx;
						debug_2_ws[k*ny*nx+j*nx+i]=2.*debug_2_ws1[k*ny*nx+j*nx+i]-debug_2_ws2[k*ny*nx+j*nx+i]+tempwx2+tempwy2
										  -tempuxz*vvs2*dtz*dtx-tempvyz*vvs2*dtz*dtx;
					}//for(i=nleft;i<nright;i++) end

	//	Debug level 2 purpose code end
	/////////////////////////////////////////////////////////////
#endif
			// #pragma omp parallel for
			// Z Y X
			for ( k=ntop; k<nbottom; k++ ) {

				vvp2=vpp2(k);
				vvs2=vss2(k);

				vvp2_dtx_dtx = vvp2*dtx*dtx;
				vvs2_dtx_dtx = vvs2*dtx*dtx;
				vvs2_dtz_dtx = vvs2*dtz*dtx;
				vvp2_dtz_dtx = vvp2*dtz*dtx;
				for ( j=nfront; j<nback; j++ ) {
					nIndex = POSITION_INDEX_X(k,j,nleft);
					for ( i=nleft; i<nright; i++ ) {

						tempux2 = c0*u_x[nIndex];
						tempvx2 = c0*v_x[nIndex];
						tempwx2 = c0*w_x[nIndex];
						for(kk=1;kk<=mm;kk++)
						{
							tempux2=tempux2+c[kk-1][0]*(u_x[POSITION_INDEX_X(k,j,i+kk)]+u_x[POSITION_INDEX_X(k,j,i-kk)]);

							tempvx2=tempvx2+c[kk-1][0]*(v_x[POSITION_INDEX_X(k,j,i+kk)]+v_x[POSITION_INDEX_X(k,j,i-kk)]);

							tempwx2=tempwx2+c[kk-1][0]*(w_x[POSITION_INDEX_X(k,j,i+kk)]+w_x[POSITION_INDEX_X(k,j,i-kk)]);
						} //for(kk=1;kk<=mm;kk++) end

						tempux2 *= vvp2_dtx_dtx;
						tempvx2 *= vvs2_dtx_dtx;
						tempwx2 *= vvs2_dtx_dtx;
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						float debug_tempux2 = 0;
						float debug_tempvx2 = 0;
						float debug_tempwx2 = 0;
						for(kk=1;kk<=mm;kk++)
						{
							debug_tempux2=debug_tempux2+c[kk-1][0]*(debug_2_u[k*ny*nx+j*nx+(i+kk)]+debug_2_u[k*ny*nx+j*nx+(i-kk)]);

							debug_tempvx2=debug_tempvx2+c[kk-1][0]*(debug_2_v[k*ny*nx+j*nx+(i+kk)]+debug_2_v[k*ny*nx+j*nx+(i-kk)]);

							debug_tempwx2=debug_tempwx2+c[kk-1][0]*(debug_2_w[k*ny*nx+j*nx+(i+kk)]+debug_2_w[k*ny*nx+j*nx+(i-kk)]);

						} //for(kk=1;kk<=mm;kk++) end
						debug_tempux2=(debug_tempux2+c0*debug_2_u[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
						debug_tempvx2=(debug_tempvx2+c0*debug_2_v[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						debug_tempwx2=(debug_tempwx2+c0*debug_2_w[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						if ( fabs(debug_tempux2 - tempux2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Ux2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff:%f\n", l, i, j, k,debug_tempux2, tempux2 ,fabs(debug_tempux2-tempux2));
						if ( fabs(debug_tempvx2 - tempvx2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Vx2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff:%f\n", l, i, j, k,debug_tempvx2, tempvx2 ,fabs(debug_tempvx2-tempvx2));
						if ( fabs(debug_tempwx2 - tempwx2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Wx2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff:%f\n", l, i, j, k,debug_tempwx2, tempwx2 ,fabs(debug_tempwx2-tempwx2));
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif

						tempuxy = 0.0;
						tempvxy = 0.0;

						for(kkk=1;kkk<=mm;kkk++)
						{
							for(kk=1;kk<=mm;kk++)
							{
								current_c = c[kkk-1][1+kk];

								_tempuxy  = u_x[POSITION_INDEX_X(k,j+kkk,i+kk)];
								_tempvxy  = v_x[POSITION_INDEX_X(k,j+kkk,i+kk)];

								_tempuxy -= u_x[POSITION_INDEX_X(k,j-kkk,i+kk)];
								_tempvxy -= v_x[POSITION_INDEX_X(k,j-kkk,i+kk)];

								_tempuxy += u_x[POSITION_INDEX_X(k,j-kkk,i-kk)];
								_tempvxy += v_x[POSITION_INDEX_X(k,j-kkk,i-kk)];

								_tempuxy -= u_x[POSITION_INDEX_X(k,j+kkk,i-kk)];
								_tempvxy -= v_x[POSITION_INDEX_X(k,j+kkk,i-kk)];

								tempuxy += current_c*_tempuxy;
								tempvxy += current_c*_tempvxy;
							}
						} //for(kk=1;kk<=mm;kk++) end
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						float debug_tempuxy = 0;
						float debug_tempvxy = 0;
						for(kk=1;kk<=mm;kk++)
						{
							for(kkk=1;kkk<=mm;kkk++)
							{
								debug_tempuxy=debug_tempuxy+c[kkk-1][1+kk]*(debug_2_u[k*ny*nx+(j+kkk)*nx+(i+kk)] -debug_2_u[k*ny*nx+(j-kkk)*nx+(i+kk)] +debug_2_u[k*ny*nx+(j-kkk)*nx+(i-kk)] -debug_2_u[k*ny*nx+(j+kkk)*nx+(i-kk)]);
								debug_tempvxy=debug_tempvxy+c[kkk-1][1+kk]*(debug_2_v[k*ny*nx+(j+kkk)*nx+(i+kk)] -debug_2_v[k*ny*nx+(j-kkk)*nx+(i+kk)] +debug_2_v[k*ny*nx+(j-kkk)*nx+(i-kk)] -debug_2_v[k*ny*nx+(j+kkk)*nx+(i-kk)]);
							} // for(kkk=1;kkk<=mm;kkk++) end
						} //for(kk=1;kk<=mm;kk++) end

						if ( fabs(tempuxy-debug_tempuxy) > _DEBUG_MAX_DIFF ) printf("[Warining tempuxy] l: %d i: %d j: %d k: %d debugValue: %f, Value: %f diff:%f\n", l, i, j, k, debug_tempuxy, tempuxy,fabs(debug_tempuxy-tempuxy));
						if ( fabs(tempvxy-debug_tempvxy) > _DEBUG_MAX_DIFF ) printf("[Warining tempvxy] l: %d i: %d j: %d k: %d debugValue: %f, Value: %f diff:%f\n", l, i, j, k, debug_tempvxy, tempvxy,fabs(debug_tempvxy-tempvxy));
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif

						up[nIndex] = tempux2+tempvxy*vvp2_dtz_dtx;
						vp[nIndex] = tempuxy*vvp2_dtz_dtx;
						us[nIndex] = -tempvxy*vvs2_dtz_dtx;
						vs[nIndex] = tempvx2-tempuxy*vvs2_dtz_dtx;
						ws[nIndex] = tempwx2;
						nIndex ++;
					}
				}
			}
			// X Z Y
			for ( i=nleft; i<nright; i++ ) {
				for ( k=ntop; k<nbottom; k++ ) {
					vvp2=vpp2(k);
					vvs2=vss2(k);

					vvp2_dtx_dtx = vvp2*dtx*dtx;
					vvs2_dtx_dtx = vvs2*dtx*dtx;
					vvs2_dtz_dtx = vvs2*dtz*dtx;
					vvp2_dtz_dtx = vvp2*dtz*dtx;

					nIndex 		= POSITION_INDEX_Y(k,nfront,i);
					nIndex_X 	= POSITION_INDEX_X(k,nfront,i);
					for ( j=nfront; j<nback; j++ ) {
						tempuy2 = c0*u_y[nIndex];
						tempvy2 = c0*v_y[nIndex];
						tempwy2 = c0*w_y[nIndex];

						for(kk=1;kk<=mm;kk++)
						{
							tempuy2=tempuy2+c[kk-1][0]*(u_y[POSITION_INDEX_Y(k, j+kk, i)]+u_y[POSITION_INDEX_Y(k, j-kk, i)]);
							tempvy2=tempvy2+c[kk-1][0]*(v_y[POSITION_INDEX_Y(k, j+kk, i)]+v_y[POSITION_INDEX_Y(k, j-kk, i)]);
							tempwy2=tempwy2+c[kk-1][0]*(w_y[POSITION_INDEX_Y(k, j+kk, i)]+w_y[POSITION_INDEX_Y(k, j-kk, i)]);
						} //for(kk=1;kk<=mm;kk++) end

						tempuy2 *= vvs2_dtx_dtx;
						tempvy2 *= vvp2_dtx_dtx;
						tempwy2 *= vvs2_dtx_dtx;
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						float debug_tempuy2 = 0;
						float debug_tempvy2 = 0;
						float debug_tempwy2 = 0;
						for(kk=1;kk<=mm;kk++)
						{
							debug_tempuy2=debug_tempuy2+c[kk-1][0]*(debug_2_u[k*ny*nx+(j+kk)*nx+i]+debug_2_u[k*ny*nx+(j-kk)*nx+i]);

							debug_tempvy2=debug_tempvy2+c[kk-1][0]*(debug_2_v[k*ny*nx+(j+kk)*nx+i]+debug_2_v[k*ny*nx+(j-kk)*nx+i]);

							debug_tempwy2=debug_tempwy2+c[kk-1][0]*(debug_2_w[k*ny*nx+(j+kk)*nx+i]+debug_2_w[k*ny*nx+(j-kk)*nx+i]);

						} //for(kk=1;kk<=mm;kk++) end
						debug_tempuy2=(debug_tempuy2+c0*debug_2_u[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						debug_tempvy2=(debug_tempvy2+c0*debug_2_v[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
						debug_tempwy2=(debug_tempwy2+c0*debug_2_w[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
						if ( fabs(debug_tempuy2 - tempuy2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Uy2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff %f\n", l, i, j, k, debug_tempuy2, tempuy2,fabs(debug_tempuy2-tempuy2) );
						if ( fabs(debug_tempvy2 - tempvy2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Vy2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff %f\n", l, i, j, k, debug_tempvy2, tempvy2,fabs(debug_tempvy2-tempvy2) );
						if ( fabs(debug_tempwy2 - tempwy2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Wy2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff %f\n", l, i, j, k, debug_tempwy2, tempwy2,fabs(debug_tempwy2-tempwy2) );
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif

						tempvyz = 0.0;
						tempwyz = 0.0;

						for(kkk=1;kkk<=mm;kkk++)
						{
							for(kk=1;kk<=mm;kk++)
							{
								current_c = c[kkk-1][1+kk];

								_tempvyz  = v_y[POSITION_INDEX_Y(k+kkk,j+kk,i)];
								_tempwyz  = w_y[POSITION_INDEX_Y(k+kkk,j+kk,i)];

								_tempvyz -= v_y[POSITION_INDEX_Y(k-kkk,j+kk,i)];
								_tempwyz -= w_y[POSITION_INDEX_Y(k-kkk,j+kk,i)];

								_tempvyz += v_y[POSITION_INDEX_Y(k-kkk,j-kk,i)];
								_tempwyz += w_y[POSITION_INDEX_Y(k-kkk,j-kk,i)];

								_tempvyz -= v_y[POSITION_INDEX_Y(k+kkk,j-kk,i)];
								_tempwyz -= w_y[POSITION_INDEX_Y(k+kkk,j-kk,i)];

								tempvyz += current_c*_tempvyz;
								tempwyz += current_c*_tempwyz;
							}
						} //for(kk=1;kk<=mm;kk++) end
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						float debug_tempvyz = 0;
						float debug_tempwyz = 0;
						for(kk=1;kk<=mm;kk++)
						{
							for(kkk=1;kkk<=mm;kkk++)
							{
								debug_tempvyz=debug_tempvyz+c[kkk-1][1+kk]*(debug_2_v[(k+kkk)*ny*nx+(j+kk)*nx+i] -debug_2_v[(k-kkk)*ny*nx+(j+kk)*nx+i] +debug_2_v[(k-kkk)*ny*nx+(j-kk)*nx+i] -debug_2_v[(k+kkk)*ny*nx+(j-kk)*nx+i]);
								debug_tempwyz=debug_tempwyz+c[kkk-1][1+kk]*(debug_2_w[(k+kkk)*ny*nx+(j+kk)*nx+i] -debug_2_w[(k-kkk)*ny*nx+(j+kk)*nx+i] +debug_2_w[(k-kkk)*ny*nx+(j-kk)*nx+i] -debug_2_w[(k+kkk)*ny*nx+(j-kk)*nx+i]);
							} // for(kkk=1;kkk<=mm;kkk++) end
						} //for(kk=1;kk<=mm;kk++) end

						if ( fabs(tempvyz-debug_tempvyz) > _DEBUG_MAX_DIFF ) printf("[Warining tempvyz] l: %d i: %d j: %d k: %d debugValue: %f, Value: %f diff:%f\n", l, i, j, k, debug_tempvyz, tempvyz,fabs(debug_tempvyz- tempvyz));
						if ( fabs(tempwyz-debug_tempwyz) > _DEBUG_MAX_DIFF ) printf("[Warining tempwyz] l: %d i: %d j: %d k: %d debugValue: %f, Value: %f diff:%f\n", l, i, j, k, debug_tempwyz, tempwyz,fabs(debug_tempwyz- tempwyz));
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif

						vp[nIndex_X] += tempvy2 + tempwyz*vvp2_dtz_dtx;
						us[nIndex_X] += tempuy2;
						wp[nIndex_X]  = tempvyz*vvp2_dtz_dtx;
						vs[nIndex_X] -= tempwyz*vvs2_dtz_dtx;
						ws[nIndex_X] += tempwy2	- tempvyz*vvs2_dtz_dtx;

						nIndex++;
						nIndex_X += nx;
					}
				}
			}
			// Y X Z
			for ( j=nfront; j<nback; j++ ) {
				for ( i=nleft; i<nright; i++ ) {
					nIndex 		= POSITION_INDEX_Z(ntop, j, i);
					nIndex_X 	= POSITION_INDEX_X(ntop, j, i);
					for ( k=ntop; k<nbottom; k++ ) {
						vvp2=vpp2(k);
						vvs2=vss2(k);

						vvs2_dtz_dtz = vvs2*dtz*dtz;
						vvp2_dtz_dtz = vvp2*dtz*dtz;
						vvs2_dtz_dtx = vvs2*dtz*dtx;
						vvp2_dtz_dtx = vvp2*dtz*dtx;

						tempuz2 = c0*u_z[nIndex];
						tempvz2 = c0*v_z[nIndex];
						tempwz2 = c0*w_z[nIndex];
						for(kk=1;kk<=mm;kk++)
						{
							tempuz2=tempuz2+c[kk-1][0]*(u_z[POSITION_INDEX_Z(k+kk, j, i)]+u_z[POSITION_INDEX_Z(k-kk, j, i)]);

							tempvz2=tempvz2+c[kk-1][0]*(v_z[POSITION_INDEX_Z(k+kk, j, i)]+v_z[POSITION_INDEX_Z(k-kk, j, i)]);

							tempwz2=tempwz2+c[kk-1][0]*(w_z[POSITION_INDEX_Z(k+kk, j, i)]+w_z[POSITION_INDEX_Z(k-kk, j, i)]);
						}

						tempuz2 *= vvs2_dtz_dtz;
						tempvz2 *= vvs2_dtz_dtz;
						tempwz2 *= vvp2_dtz_dtz;
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						float debug_tempuz2 = 0;
						float debug_tempvz2 = 0;
						float debug_tempwz2 = 0;
						for(kk=1;kk<=mm;kk++)
						{
							debug_tempuz2=debug_tempuz2+c[kk-1][0]*(debug_2_u[(k+kk)*ny*nx+j*nx+i]+debug_2_u[(k-kk)*ny*nx+j*nx+i]);

							debug_tempvz2=debug_tempvz2+c[kk-1][0]*(debug_2_v[(k+kk)*ny*nx+j*nx+i]+debug_2_v[(k-kk)*ny*nx+j*nx+i]);

							debug_tempwz2=debug_tempwz2+c[kk-1][0]*(debug_2_w[(k+kk)*ny*nx+j*nx+i]+debug_2_w[(k-kk)*ny*nx+j*nx+i]);

						} //for(kk=1;kk<=mm;kk++) end
						debug_tempuz2=(debug_tempuz2+c0*debug_2_u[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;
						debug_tempvz2=(debug_tempvz2+c0*debug_2_v[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;
						debug_tempwz2=(debug_tempwz2+c0*debug_2_w[k*ny*nx+j*nx+i])*vvp2*dtz*dtz;
						if ( fabs(debug_tempuz2 - tempuz2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Uz2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff: %f \n", l, i, j, k, debug_tempuz2, tempuz2 ,fabs(debug_tempuz2-tempuz2));
						if ( fabs(debug_tempvz2 - tempvz2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Vz2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff: %f \n", l, i, j, k, debug_tempvz2, tempvz2 ,fabs(debug_tempvz2-tempvz2));
						if ( fabs(debug_tempwz2 - tempwz2) > _DEBUG_MAX_DIFF ) printf("[Temp Warining Wz2] l: %d i: %d j: %d k: %d debugvalue: %f value: %f diff: %f \n", l, i, j, k, debug_tempwz2, tempwz2 ,fabs(debug_tempwz2-tempwz2));
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif
						tempuxz = 0.0;
						tempwxz = 0.0;

						for(kk=1;kk<=mm;kk++)
						{
							for(kkk=1;kkk<=mm;kkk++)
							{
								current_c = c[kkk-1][1+kk];
								_tempuxz  = u_z[POSITION_INDEX_Z(k+kkk,j,i+kk)];
								_tempwxz  = w_z[POSITION_INDEX_Z(k+kkk,j,i+kk)];

								_tempuxz -= u_z[POSITION_INDEX_Z(k-kkk,j,i+kk)];
								_tempwxz -= w_z[POSITION_INDEX_Z(k-kkk,j,i+kk)];

								_tempuxz += u_z[POSITION_INDEX_Z(k-kkk,j,i-kk)];
								_tempwxz += w_z[POSITION_INDEX_Z(k-kkk,j,i-kk)];

								_tempuxz -= u_z[POSITION_INDEX_Z(k+kkk,j,i-kk)];
								_tempwxz -= w_z[POSITION_INDEX_Z(k+kkk,j,i-kk)];

								tempuxz += current_c*_tempuxz;
								tempwxz += current_c*_tempwxz;
							}
						} //for(kk=1;kk<=mm;kk++) end
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						float debug_tempuxz = 0;
						float debug_tempwxz = 0;
						for(kk=1;kk<=mm;kk++)
						{
							for(kkk=1;kkk<=mm;kkk++)
							{
								debug_tempuxz=debug_tempuxz+c[kkk-1][1+kk]*(debug_2_u[(k+kkk)*ny*nx+j*nx+(i+kk)] -debug_2_u[(k-kkk)*ny*nx+j*nx+(i+kk)] +debug_2_u[(k-kkk)*ny*nx+j*nx+(i-kk)] -debug_2_u[(k+kkk)*ny*nx+j*nx+(i-kk)]);
								debug_tempwxz=debug_tempwxz+c[kkk-1][1+kk]*(debug_2_w[(k+kkk)*ny*nx+j*nx+(i+kk)] -debug_2_w[(k-kkk)*ny*nx+j*nx+(i+kk)] +debug_2_w[(k-kkk)*ny*nx+j*nx+(i-kk)] -debug_2_w[(k+kkk)*ny*nx+j*nx+(i-kk)]);
							} // for(kkk=1;kkk<=mm;kkk++) end
						} //for(kk=1;kk<=mm;kk++) end

						if ( fabs(tempuxz-debug_tempuxz) > _DEBUG_MAX_DIFF ) printf("[Warining tempuxz] l: %d i: %d j: %d k: %d debugValue: %f, Value: %f diff:%f\n", l, i, j, k, debug_tempuxz, tempuxz,fabs(debug_tempuxz- tempuxz));
						if ( fabs(tempwxz-debug_tempwxz) > _DEBUG_MAX_DIFF ) printf("[Warining tempwxz] l: %d i: %d j: %d k: %d debugValue: %f, Value: %f diff:%f\n", l, i, j, k, debug_tempwxz, tempwxz,fabs(debug_tempwxz- tempwxz));
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif

						up[nIndex_X] += tempwxz*vvp2_dtz_dtx;
						wp[nIndex_X] += tempwz2+tempuxz*vvp2_dtz_dtx +(int)(i==ncx_shot-1&&j==ncy_shot-1&&k==ncz_shot-1)*wave[l-1];
						us[nIndex_X] += tempuz2 - tempwxz*vvs2_dtz_dtx;
						vs[nIndex_X] += tempvz2;
						ws[nIndex_X] -= tempuxz*vvs2_dtz_dtx;
#ifdef _DEBUG_TEST
	/////////////////////////////////////////////////////////////
	//	Debug purpose code start
						if(nIndex_X != POSITION_INDEX_X(k,j,i)){printf("[Warining] Index X : %d does not match %d.\n",nIndex_X,POSITION_INDEX_X(k,j,i));}
	//	Debug purpose code end
	/////////////////////////////////////////////////////////////
#endif
						nIndex++;
						nIndex_X += nSliceSize;
					}
				}
			}

			int nIndex_X = POSITION_INDEX_X(ntop, nfront, nleft);

			// #pragma omp parallel for
			nIndex = POSITION_INDEX_X(ntop, nfront, nleft);
			for(k=ntop;k<nbottom;k++) {
				for(j=nfront;j<nback;j++) {
					for(i=nleft;i<nright;i++) {

						u_x[nIndex] = (up[nIndex] = (2*up1[nIndex] - up2[nIndex] + up[nIndex])) + ((us[nIndex] = 2*us1[nIndex] - us2[nIndex] + us[nIndex]));
						v_x[nIndex] = (vp[nIndex] = (2*vp1[nIndex] - vp2[nIndex] + vp[nIndex])) + ((vs[nIndex] = 2*vs1[nIndex] - vs2[nIndex] + vs[nIndex]));
						w_x[nIndex] = (wp[nIndex] = (2*wp1[nIndex] - wp2[nIndex] + wp[nIndex])) + ((ws[nIndex] = 2*ws1[nIndex] - ws2[nIndex] + ws[nIndex]));
#ifdef _DEBUG_LEVEL_2
	/////////////////////////////////////////////////////////////
	//	Debug level 2 code start
						float diff;
                        debug_2_u[k*ny*nx+j*nx+i]=debug_2_up[k*ny*nx+j*nx+i]+debug_2_us[k*ny*nx+j*nx+i];
                        debug_2_v[k*ny*nx+j*nx+i]=debug_2_vp[k*ny*nx+j*nx+i]+debug_2_vs[k*ny*nx+j*nx+i];
                        debug_2_w[k*ny*nx+j*nx+i]=debug_2_wp[k*ny*nx+j*nx+i]+debug_2_ws[k*ny*nx+j*nx+i];

						if((diff = fabs(debug_2_u[k*ny*nx+j*nx+i] - u_x[nIndex])) > _DEBUG_MAX_DIFF) { printf("[Warining 2 U] l: %d nIndex:%d index: %d k: %d j: %d i: %d diff: %f\n", l, nIndex, k*ny*nx+j*nx+i, k, j, i, diff);} //u[k*ny*nx+j*nx+i]=
						if((diff = fabs(debug_2_v[k*ny*nx+j*nx+i] - v_x[nIndex])) > _DEBUG_MAX_DIFF) { printf("[Warining 2 V] l: %d nIndex:%d index: %d k: %d j: %d i: %d diff: %f\n", l, nIndex, k*ny*nx+j*nx+i, k, j, i, diff);} //v[k*ny*nx+j*nx+i]=
						if((diff = fabs(debug_2_w[k*ny*nx+j*nx+i] - w_x[nIndex])) > _DEBUG_MAX_DIFF) { printf("[Warining 2 W] l: %d nIndex:%d index: %d k: %d j: %d i: %d diff: %f\n", l, nIndex, k*ny*nx+j*nx+i, k, j, i, diff);} //w[k*ny*nx+j*nx+i]=


                        debug_2_up2[k*ny*nx+j*nx+i]=debug_2_up1[k*ny*nx+j*nx+i];
                        debug_2_up1[k*ny*nx+j*nx+i]=debug_2_up[k*ny*nx+j*nx+i];
                        debug_2_us2[k*ny*nx+j*nx+i]=debug_2_us1[k*ny*nx+j*nx+i];
                        debug_2_us1[k*ny*nx+j*nx+i]=debug_2_us[k*ny*nx+j*nx+i];
                        debug_2_vp2[k*ny*nx+j*nx+i]=debug_2_vp1[k*ny*nx+j*nx+i];
                        debug_2_vp1[k*ny*nx+j*nx+i]=debug_2_vp[k*ny*nx+j*nx+i];
                        debug_2_vs2[k*ny*nx+j*nx+i]=debug_2_vs1[k*ny*nx+j*nx+i];
                        debug_2_vs1[k*ny*nx+j*nx+i]=debug_2_vs[k*ny*nx+j*nx+i];
                        debug_2_wp2[k*ny*nx+j*nx+i]=debug_2_wp1[k*ny*nx+j*nx+i];
                        debug_2_wp1[k*ny*nx+j*nx+i]=debug_2_wp[k*ny*nx+j*nx+i];
                        debug_2_ws2[k*ny*nx+j*nx+i]=debug_2_ws1[k*ny*nx+j*nx+i];
                        debug_2_ws1[k*ny*nx+j*nx+i]=debug_2_ws[k*ny*nx+j*nx+i];
	//	Debug level 2 code end
	/////////////////////////////////////////////////////////////
#endif

// #ifdef _DEBUG_TEST
// 	/////////////////////////////////////////////////////////////
// 	//	Debug purpose code start
// 			float diff;
// 						if((diff = fabs(debug_up[k*ny*nx+j*nx+i]+debug_us[k*ny*nx+j*nx+i] - u_x[nIndex])) > _DEBUG_MAX_DIFF) { printf("[Warining U] l: %d nIndex:%d index: %d k: %d j: %d i: %d diff: %f\n", l, nIndex, k*ny*nx+j*nx+i, k, j, i, diff);} //u[k*ny*nx+j*nx+i]=
// 						if((diff = fabs(debug_vp[k*ny*nx+j*nx+i]+debug_vs[k*ny*nx+j*nx+i] - v_x[nIndex])) > _DEBUG_MAX_DIFF) { printf("[Warining V] l: %d nIndex:%d index: %d k: %d j: %d i: %d diff: %f\n", l, nIndex, k*ny*nx+j*nx+i, k, j, i, diff);} //v[k*ny*nx+j*nx+i]=
// 						if((diff = fabs(debug_wp[k*ny*nx+j*nx+i]+debug_ws[k*ny*nx+j*nx+i] - w_x[nIndex])) > _DEBUG_MAX_DIFF) { printf("[Warining W] l: %d nIndex:%d index: %d k: %d j: %d i: %d diff: %f\n", l, nIndex, k*ny*nx+j*nx+i, k, j, i, diff);} //w[k*ny*nx+j*nx+i]=
// 	//	Debug purpose code end
// 	/////////////////////////////////////////////////////////////
// #endif
						nIndex++;
					}//for(i=nleft;i<nright;i++) end
					nIndex += nx+nleft-nright;
				}
				nIndex+=nx*(ny + nfront - nback);
			}

//			#pragma omp parallel for
			for(k=ntop;k<nbottom;k++) {
				for(j=nfront;j<nback;j++) {
					for(i=nleft;i<nright;i++) {
						nIndex_X = POSITION_INDEX_X(k, j ,i);
						nIndex_Y = POSITION_INDEX_Y(k, j, i);
						nIndex_Z = POSITION_INDEX_Z(k, j ,i);
						u_y[nIndex_Y] = u_z[nIndex_Z] = u_x[nIndex_X];
						v_y[nIndex_Y] = v_z[nIndex_Z] = v_x[nIndex_X];
						w_y[nIndex_Y] = w_z[nIndex_Z] = w_x[nIndex_X];
					}
				}
			}

			swap_temp = up2; up2 = up1; up1 = up; up = swap_temp;
			swap_temp = vp2; vp2 = vp1; vp1 = vp; vp = swap_temp;
			swap_temp = wp2; wp2 = wp1; wp1 = wp; wp = swap_temp;

		}//for(l=1;l<=lt;l++) end

//#ifndef _DEBUG_TEST
		for(int i=0;i<nx;++i){
			for(int j=0;j<ny;++j)
			fprintf(flog,"%f\n",up[POSITION_INDEX_X(169,i,j)]);
		}
		fwrite(up1+POSITION_INDEX_X(169,0,0),sizeof(float),nSliceSize,fout);
//#endif

	}//for(ishot=1;ishot<=nshot;ishot++) end
	fclose(fout);

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
	free(wave);

#ifdef _DEBUG_TEST
	free(debug_up);
	free(debug_vp);
	free(debug_wp);
	free(debug_us);
	free(debug_vs);
	free(debug_ws);
#endif

#ifdef _DEBUG_LEVEL_2
	free(debug_2_u  );
    free(debug_2_v  );
    free(debug_2_w  );
    free(debug_2_up );
    free(debug_2_up1);
    free(debug_2_up2);
    free(debug_2_vp );
    free(debug_2_vp1);
    free(debug_2_vp2);
    free(debug_2_wp );
    free(debug_2_wp1);
    free(debug_2_wp2);
    free(debug_2_us );
    free(debug_2_us1);
    free(debug_2_us2);
    free(debug_2_vs );
    free(debug_2_vs1);
    free(debug_2_vs2);
    free(debug_2_ws );
    free(debug_2_ws1);
    free(debug_2_ws2);
#endif

	gettimeofday(&end,NULL);
	all_time = (end.tv_sec-start.tv_sec)+(float)(end.tv_usec-start.tv_usec)/1000000.0;
	printf("run time:\t%f s\n",all_time);
	// flog = fopen(logfile,"a");
	fprintf(flog,"\nrun time:\t%f s\n\n",all_time);
	fclose(flog);
	flog = fopen(logfile,"a");
	fprintf(flog,"------------end time------------\n");
	fclose(flog);

	system(tmp);
	return 1;
}
