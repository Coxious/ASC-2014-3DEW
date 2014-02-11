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

float *u_x, *v_x, *w_x, *u_y, *v_y, *w_y, *u_z, *v_z, *w_z, *up , *up1, *up2, *vp , *vp1, *vp2,
		 *wp , *wp1, *wp2, *us , *us1, *us2, *vs , *vs1, *vs2, *ws , *ws1, *ws2, *swap_temp;

float *debug_up,*debug_vp,*debug_wp,*debug_us,*debug_vs,*debug_ws;

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

	u_x = (float *) calloc(sizeof(float), nSize);
	v_x = (float *) calloc(sizeof(float), nSize);
	w_x = (float *) calloc(sizeof(float), nSize);
	u_y = (float *) calloc(sizeof(float), nSize);
	v_y = (float *) calloc(sizeof(float), nSize);
	w_y = (float *) calloc(sizeof(float), nSize);
	u_z = (float *) calloc(sizeof(float), nSize);
	v_z = (float *) calloc(sizeof(float), nSize);
	w_z = (float *) calloc(sizeof(float), nSize);
	up  = (float *) calloc(sizeof(float), nSize);
	up1 = (float *) calloc(sizeof(float), nSize);
	up2 = (float *) calloc(sizeof(float), nSize);
	vp  = (float *) calloc(sizeof(float), nSize);
	vp1 = (float *) calloc(sizeof(float), nSize);
	vp2 = (float *) calloc(sizeof(float), nSize);
	wp  = (float *) calloc(sizeof(float), nSize);
	wp1 = (float *) calloc(sizeof(float), nSize);
	wp2 = (float *) calloc(sizeof(float), nSize);
	us  = (float *) calloc(sizeof(float), nSize);
	us1 = (float *) calloc(sizeof(float), nSize);
	us2 = (float *) calloc(sizeof(float), nSize);
	vs  = (float *) calloc(sizeof(float), nSize);
	vs1 = (float *) calloc(sizeof(float), nSize);
	vs2 = (float *) calloc(sizeof(float), nSize);
	ws  = (float *) calloc(sizeof(float), nSize);
	ws1 = (float *) calloc(sizeof(float), nSize);
	ws2 = (float *) calloc(sizeof(float), nSize);

	wave    = (float*) calloc(sizeof(float), lt);


	/////////////////////////////////////////////////////////////
	//	Debug purpos code start

	debug_up= (float *) calloc(sizeof(float), nSize);
	debug_vp= (float *) calloc(sizeof(float), nSize);
	debug_wp= (float *) calloc(sizeof(float), nSize);
	debug_us= (float *) calloc(sizeof(float), nSize);
	debug_vs= (float *) calloc(sizeof(float), nSize);
	debug_ws= (float *) calloc(sizeof(float), nSize);

	//	Debug purpos code end
	/////////////////////////////////////////////////////////////

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


	/////////////////////////////////////////////////////////////
	//	Debug purpos code start
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
                                tempuxz=tempuxz+c[kkk-1][1+kk]*(u_x[(k+kkk)*ny*nx+j*nx+(i+kk)]
                                                   -u_x[(k-kkk)*ny*nx+j*nx+(i+kk)]
                                                   +u_x[(k-kkk)*ny*nx+j*nx+(i-kk)]
                                                   -u_x[(k+kkk)*ny*nx+j*nx+(i-kk)]);
                                tempuxy=tempuxy+c[kkk-1][1+kk]*(u_x[k*ny*nx+(j+kkk)*nx+(i+kk)]
                                                   -u_x[k*ny*nx+(j-kkk)*nx+(i+kk)]
                                                   +u_x[k*ny*nx+(j-kkk)*nx+(i-kk)]
                                                   -u_x[k*ny*nx+(j+kkk)*nx+(i-kk)]);

                                tempvyz=tempvyz+c[kkk-1][1+kk]*(v_x[(k+kkk)*ny*nx+(j+kk)*nx+i]
                                                   -v_x[(k-kkk)*ny*nx+(j+kk)*nx+i]
                                                   +v_x[(k-kkk)*ny*nx+(j-kk)*nx+i]
                                                   -v_x[(k+kkk)*ny*nx+(j-kk)*nx+i]);
                                tempvxy=tempvxy+c[kkk-1][1+kk]*(v_x[k*ny*nx+(j+kkk)*nx+(i+kk)]
                                                   -v_x[k*ny*nx+(j-kkk)*nx+(i+kk)]
                                                   +v_x[k*ny*nx+(j-kkk)*nx+(i-kk)]
                                                   -v_x[k*ny*nx+(j+kkk)*nx+(i-kk)]);

                                tempwyz=tempwyz+c[kkk-1][1+kk]*(w_x[(k+kkk)*ny*nx+(j+kk)*nx+i]
                                                   -w_x[(k-kkk)*ny*nx+(j+kk)*nx+i]
                                                   +w_x[(k-kkk)*ny*nx+(j-kk)*nx+i]
                                                   -w_x[(k+kkk)*ny*nx+(j-kk)*nx+i]);
                                tempwxz=tempwxz+c[kkk-1][1+kk]*(w_x[(k+kkk)*ny*nx+j*nx+(i+kk)]
                                                   -w_x[(k-kkk)*ny*nx+j*nx+(i+kk)]
                                                   +w_x[(k-kkk)*ny*nx+j*nx+(i-kk)]
                                                   -w_x[(k+kkk)*ny*nx+j*nx+(i-kk)]);
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
	//	Debug purpos code end
	/////////////////////////////////////////////////////////////



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
						nIndex++;
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

						up[nIndex] = tempux2+tempvxy*vvs2_dtz_dtx;
						vp[nIndex] = tempuxy;
						us[nIndex] = -tempvxy*vvs2_dtz_dtx;
						vs[nIndex] = tempvx2-tempuxy*vvs2_dtz_dtx;
						ws[nIndex] = tempwx2;
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
							tempuy2=tempuy2+c[kk-1][0]*(u_y[POSITION_INDEX_Y(i, k, j+kk)]+u_y[POSITION_INDEX_Y(i, k, j-kk)]);

							tempvy2=tempvy2+c[kk-1][0]*(v_y[POSITION_INDEX_Y(i, k, j+kk)]+v_y[POSITION_INDEX_Y(i, k, j-kk)]);

							tempwy2=tempwy2+c[kk-1][0]*(w_y[POSITION_INDEX_Y(i, k, j+kk)]+w_y[POSITION_INDEX_Y(i, k, j-kk)]);
						} //for(kk=1;kk<=mm;kk++) end

						tempuy2 *= vvs2_dtx_dtx;
						tempvy2 *= vvp2_dtx_dtx;
						tempwy2 *= vvs2_dtx_dtx;

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

						vp[nIndex_X] += tempvy2 + tempwyz*vvp2_dtz_dtx;
						us[nIndex_X] += tempuy2;
						wp[nIndex_X]  = +tempvyz*vvp2_dtz_dtx;
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
					nIndex 		= POSITION_INDEX_Z(ntop, i, j);
					nIndex_X 	= POSITION_INDEX_X(ntop, i, j);
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
							tempuz2=tempuz2+c[kk-1][0]*(u_z[POSITION_INDEX_Z(j, i, k+kk)]+u_z[POSITION_INDEX_Z(j, i, k-kk)]);

							tempvz2=tempvz2+c[kk-1][0]*(v_z[POSITION_INDEX_Z(j, i, k+kk)]+v_z[POSITION_INDEX_Z(j, i, k-kk)]);

							tempwz2=tempwz2+c[kk-1][0]*(w_z[POSITION_INDEX_Z(j, i, k+kk)]+w_z[POSITION_INDEX_Z(j, i, k-kk)]);
						}

						tempuz2 *= vvs2_dtz_dtz;
						tempvz2 *= vvs2_dtz_dtz;
						tempwz2 *= vvp2_dtz_dtz;

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

						up[nIndex_X] += tempwxz*vvp2_dtz_dtx;
						wp[nIndex_X] += tempwz2+tempuxz*vvp2_dtz_dtx +(i==ncx_shot-1&&j==ncy_shot-1&&k==ncz_shot-1)*wave[l-1];
						us[nIndex_X] += tempuz2 - tempwxz*vvs2_dtz_dtx;
						vs[nIndex_X] += tempvz2;
						ws[nIndex_X] -= tempuxz*vvs2_dtz_dtx;

						nIndex++;
						nIndex_X += nSliceSize;
					}
				}
			}

			int nIndex_X = POSITION_INDEX_X(ntop, nfront, nleft);

			// #pragma omp parallel for
			for(k=ntop;k<nbottom;k++) {
				for(j=nfront;j<nback;j++) {
					for(i=nleft;i<nright;i++) {
						u_x[nIndex] = 2*up1[nIndex] - up2[nIndex] + up[nIndex] + 2*us1[nIndex] - us2[nIndex] + us[nIndex];
						v_x[nIndex] = 2*vp1[nIndex] - vp2[nIndex] + vp[nIndex] + 2*vs1[nIndex] - vs2[nIndex] + vs[nIndex];
						w_x[nIndex] = 2*wp1[nIndex] - wp2[nIndex] + wp[nIndex] + 2*ws1[nIndex] - ws2[nIndex] + ws[nIndex];
						nIndex++;



	/////////////////////////////////////////////////////////////
	//	Debug purpos code start
	                    if(up[k*ny*nx+j*nx+i]+us[k*ny*nx+j*nx+i] - u_x[nIndex] > 1e-6) { printf("[Warining] lt: %d nIndex:%d index: %d k: %d j: %d i: %d", lt, nIndex, k*ny*nx+j*nx+i, k, j, i);} //u[k*ny*nx+j*nx+i]=
                        if(vp[k*ny*nx+j*nx+i]+vs[k*ny*nx+j*nx+i] - v_x[nIndex] > 1e-6) { printf("[Warining] lt: %d nIndex:%d index: %d k: %d j: %d i: %d", lt, nIndex, k*ny*nx+j*nx+i, k, j, i);} //v[k*ny*nx+j*nx+i]=
                        if(wp[k*ny*nx+j*nx+i]+ws[k*ny*nx+j*nx+i] - w_x[nIndex] > 1e-6) { printf("[Warining] lt: %d nIndex:%d index: %d k: %d j: %d i: %d", lt, nIndex, k*ny*nx+j*nx+i, k, j, i);} //w[k*ny*nx+j*nx+i]=

	//	Debug purpos code end
	/////////////////////////////////////////////////////////////
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

		// fwrite(up+POSITION_INDEX_X(169,0,0),sizeof(float),nSliceSize,fout);
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
