#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#define PIE 3.1415926

#define POSITION_INDEX(_z,_y,_x)          ((_z)*ny*nx + (_y)*nx + (_x))
#define POSITION_INDEX_X(_z,_y,_x)        ((_z)*ny*nx + (_y)*nx + (_x))
#define POSITION_INDEX_Y(_z,_y,_x)        ((_x)*nz*ny + (_z)*ny + (_y)) 
#define POSITION_INDEX_Z(_z,_y,_x)        ((_y)*nx*nz + (_x)*nz + (_z))

#define vpp2(_z) ((_z<210)?(5290000):(_z>=260?12250000:7840000))
#define vss2(_z) ((_z<210)?(1517824):(_z>=260?3644281:2277081))

#define _DEBUG_LEVEL_1_

#ifdef _DEBUG_LEVEL_1_

#define _MAX_DIFF_

#endif

float * u, * v, * w, * up, * up1, * up2, * vp, * vp1, * vp2, * wp, * wp1,
      * wp2, * us, * us1, * us2, * vs, * vs1, * vs2, * ws, * ws1, * ws2; 
float *u_x, *u_y, *u_z, *v_x, *v_y, *v_z, *w_x, *w_y, *w_z;

int main(int argc, char **argv) {
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

    float *up_out;
    float c[5][7];
    float *wave;
    float nshot,t0,tt,c0;
    float dtx,dtz;
    float xmax,px;
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
    fclose(flog);

    const int nSize = nz*ny*nx;
    const int nSliceSize = ny*nx;


    u   = ( float * )malloc(sizeof(float)*nSize);
    v   = ( float * )malloc(sizeof(float)*nSize);
    w   = ( float * )malloc(sizeof(float)*nSize);
    u_x = ( float * )malloc(sizeof(float)*nSize);
    v_x = ( float * )malloc(sizeof(float)*nSize);
    w_x = ( float * )malloc(sizeof(float)*nSize);
    u_y = ( float * )malloc(sizeof(float)*nSize);
    v_y = ( float * )malloc(sizeof(float)*nSize);
    w_y = ( float * )malloc(sizeof(float)*nSize);
    u_z = ( float * )malloc(sizeof(float)*nSize);
    v_z = ( float * )malloc(sizeof(float)*nSize);
    w_z = ( float * )malloc(sizeof(float)*nSize);
    up  = ( float * )malloc(sizeof(float)*nSize);
    up1 = ( float * )malloc(sizeof(float)*nSize); 
    up2 = ( float * )malloc(sizeof(float)*nSize); 
    vp  = ( float * )malloc(sizeof(float)*nSize);
    vp1 = ( float * )malloc(sizeof(float)*nSize);
    vp2 = ( float * )malloc(sizeof(float)*nSize); 
    wp  = ( float * )malloc(sizeof(float)*nSize); 
    wp1 = ( float * )malloc(sizeof(float)*nSize); 
    wp2 = ( float * )malloc(sizeof(float)*nSize); 
    us  = ( float * )malloc(sizeof(float)*nSize);
    us1 = ( float * )malloc(sizeof(float)*nSize); 
    us2 = ( float * )malloc(sizeof(float)*nSize); 
    vs  = ( float * )malloc(sizeof(float)*nSize); 
    vs1 = ( float * )malloc(sizeof(float)*nSize); 
    vs2 = ( float * )malloc(sizeof(float)*nSize); 
    ws  = ( float * )malloc(sizeof(float)*nSize); 
    ws1 = ( float * )malloc(sizeof(float)*nSize); 
    ws2 = ( float * )malloc(sizeof(float)*nSize); 

    wave    = (float*)malloc(sizeof(float)*lt);
    up_out  = (float*)malloc(sizeof(float)*nx*ny);

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

        // memset
        memset(u_x, 0, sizeof(float)*nSize); 
        memset(v_x, 0, sizeof(float)*nSize); 
        memset(w_x, 0, sizeof(float)*nSize); 
        memset(u_y, 0, sizeof(float)*nSize); 
        memset(v_y, 0, sizeof(float)*nSize); 
        memset(w_y, 0, sizeof(float)*nSize); 
        memset(u_z, 0, sizeof(float)*nSize); 
        memset(v_z, 0, sizeof(float)*nSize); 
        memset(w_z, 0, sizeof(float)*nSize); 
        memset(up , 0, sizeof(float)*nSize); 
        memset(up1, 0, sizeof(float)*nSize); 
        memset(up2, 0, sizeof(float)*nSize); 
        memset(vp , 0, sizeof(float)*nSize); 
        memset(vp1, 0, sizeof(float)*nSize); 
        memset(vp2, 0, sizeof(float)*nSize); 
        memset(wp , 0, sizeof(float)*nSize); 
        memset(wp1, 0, sizeof(float)*nSize); 
        memset(wp2, 0, sizeof(float)*nSize); 
        memset(us , 0, sizeof(float)*nSize); 
        memset(us1, 0, sizeof(float)*nSize); 
        memset(us2, 0, sizeof(float)*nSize); 
        memset(vs , 0, sizeof(float)*nSize); 
        memset(vs1, 0, sizeof(float)*nSize); 
        memset(vs2, 0, sizeof(float)*nSize); 
        memset(ws , 0, sizeof(float)*nSize); 
        memset(ws1, 0, sizeof(float)*nSize); 
        memset(ws2, 0, sizeof(float)*nSize); 

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

            // ZYX
            for(k=ntop;k<nbottom;k++)
            {
                vvp2=vpp2(k);
                vvs2=vss2(k);

                vvs2_dtz_dtz = vvs2*dtz*dtz;
                vvp2_dtx_dtx = vvp2*dtx*dtx;
                vvs2_dtx_dtx = vvs2*dtx*dtx;
                vvp2_dtz_dtz = vvp2*dtz*dtz;
                vvp2_dtz_dtx = vvp2*dtz*dtx;
                vvs2_dtz_dtx = vvs2*dtz*dtx;

                for(j=nfront;j<nback;j++) {
                    for(i=nleft;i<nright;i++) {

                        int nIndex = POSITION_INDEX_X(k,j,i);

                        if(i==ncx_shot-1&&j==ncy_shot-1&&k==ncz_shot-1)
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

                        for(kk=1;kk<=mm;kk++)
                        {
                            tempux2=tempux2+c[kk-1][0]*(u_x[POSITION_INDEX_X(k,j,i+kk)]+u_x[POSITION_INDEX_X(k,j,i-kk)]);

                            tempvx2=tempvx2+c[kk-1][0]*(v_x[POSITION_INDEX_X(k,j,i+kk)]+v_x[POSITION_INDEX_X(k,j,i-kk)]);

                            tempwx2=tempwx2+c[kk-1][0]*(w_x[POSITION_INDEX_X(k,j,i+kk)]+w_x[POSITION_INDEX_X(k,j,i-kk)]);
                        } //for(kk=1;kk<=mm;kk++) end

                        tempux2=(tempux2+c0*u_x[POSITION_INDEX_X(k,j,i)])*vvp2_dtx_dtx;

                        tempvx2=(tempvx2+c0*v_x[POSITION_INDEX_X(k,j,i)])*vvs2_dtx_dtx;

                        tempwx2=(tempwx2+c0*w_x[POSITION_INDEX_X(k,j,i)])*vvs2_dtx_dtx;

                        for(kk=1;kk<=mm;kk++)
                        {
                            for(kkk=1;kkk<=mm;kkk++)
                            {
                                current_c = c[kkk-1][1+kk];

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

                            } // for(kkk=1;kkk<=mm;kkk++) end
                        } //for(kk=1;kk<=mm;kk++) end
                        up[nIndex] = tempux2 + tempvxy * vvp2_dtz_dtx;
                        vp[nIndex] = tempuxy * vvp2_dtz_dtx;
                        wp[nIndex] = px * wave[l-1];
                        us[nIndex] = - tempvxy * vvs2_dtz_dtx;
                        vs[nIndex] = tempvx2 - tempuxy * vvs2_dtz_dtx;
                        ws[nIndex] = tempwx2;
                    }
                }
            }

            // X Z Y
            for(i=nleft;i<nright;i++) {
                for(k=ntop;k<nbottom;k++) {
                    vvp2=vpp2(k);
                    vvs2=vss2(k);

                    vvs2_dtz_dtz = vvs2*dtz*dtz;
                    vvp2_dtx_dtx = vvp2*dtx*dtx;
                    vvs2_dtx_dtx = vvs2*dtx*dtx;
                    vvp2_dtz_dtz = vvp2*dtz*dtz;
                    vvp2_dtz_dtx = vvp2*dtz*dtx;
                    vvs2_dtz_dtx = vvs2*dtz*dtx;

                    for(j=nfront;j<nback;j++) {

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

                        for(kk=1;kk<=mm;kk++)
                        {
                            if ( u_x[POSITION_INDEX_X(k,j+kk,i) != u_y[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
                            tempuy2=tempuy2+c[kk-1][0]*(u_y[POSITION_INDEX_Y(k,j+kk,i)]+u_y[POSITION_INDEX_Y(k,j-kk,i)]);

                            if ( v_x[POSITION_INDEX_X(k,j+kk,i) != v_y[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
                            tempvy2=tempvy2+c[kk-1][0]*(v_y[POSITION_INDEX_Y(k,j+kk,i)]+v_y[POSITION_INDEX_Y(k,j-kk,i)]);

                            if ( w_x[POSITION_INDEX_X(k,j+kk,i) != w_y[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
                            tempwy2=tempwy2+c[kk-1][0]*(w_y[POSITION_INDEX_Y(k,j+kk,i)]+w_y[POSITION_INDEX_Y(k,j-kk,i)]);
                        } //for(kk=1;kk<=mm;kk++) end

                        tempuy2=(tempuy2+c0*u_y[POSITION_INDEX_Y(k,j,i)])*vvs2_dtx_dtx;

                        tempvy2=(tempvy2+c0*v_y[POSITION_INDEX_Y(k,j,i)])*vvp2_dtx_dtx;

                        tempwy2=(tempwy2+c0*w_y[POSITION_INDEX_Y(k,j,i)])*vvs2_dtx_dtx;

                        for(kk=1;kk<=mm;kk++)
                        {
                            for(kkk=1;kkk<=mm;kkk++)
                            {
                                current_c = c[kkk-1][1+kk];

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

                            } // for(kkk=1;kkk<=mm;kkk++) end
                        } //for(kk=1;kk<=mm;kk++) end

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
            for(j=nfront;j<nback;j++) {
                for(i=nleft;i<nright;i++) {
                    for(k=ntop;k<nbottom;k++)
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

                        for(kk=1;kk<=mm;kk++)
                        {
                            if ( u_x[POSITION_INDEX_X(k,j+kk,i) != u_z[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
                            tempuz2=tempuz2+c[kk-1][0]*(u_z[POSITION_INDEX_Z(k+kk,j,i)]+u_z[POSITION_INDEX_Z(k-kk,j,i)]);

                            if ( v_x[POSITION_INDEX_X(k,j+kk,i) != v_z[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
                            tempvz2=tempvz2+c[kk-1][0]*(v_z[POSITION_INDEX_Z(k+kk,j,i)]+v_z[POSITION_INDEX_Z(k-kk,j,i)]);

                            if ( w_x[POSITION_INDEX_X(k,j+kk,i) != w_z[POSITION_INDEX_Y(k,j+kk,i)]] ) printf("[Warining] i: %d j: %d k: %d", i, j, k);
                            tempwz2=tempwz2+c[kk-1][0]*(w_z[POSITION_INDEX_Z(k+kk,j,i)]+w_z[POSITION_INDEX_Z(k-kk,j,i)]);
                        } //for(kk=1;kk<=mm;kk++) end

                        tempuz2=(tempuz2+c0*u_z[POSITION_INDEX_Z(k,j,i)])*vvs2_dtz_dtz;

                        tempvz2=(tempvz2+c0*v_z[POSITION_INDEX_Z(k,j,i)])*vvs2_dtz_dtz;

                        tempwz2=(tempwz2+c0*w_z[POSITION_INDEX_Z(k,j,i)])*vvp2_dtz_dtz;

                        for(kk=1;kk<=mm;kk++)
                        {
                            for(kkk=1;kkk<=mm;kkk++)
                            {
                                current_c = c[kkk-1][1+kk];

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

                            } // for(kkk=1;kkk<=mm;kkk++) end
                        } //for(kk=1;kk<=mm;kk++) end

                        up[nIndex_X] += tempwxz * vvp2_dtz_dtx;
                        vp[nIndex_X] += 0;
                        wp[nIndex_X] += tempwz2 + tempuxz * vvp2_dtz_dtx;
                        us[nIndex_X] += tempuz2 - tempwxz * vvs2_dtz_dtx; 
                        vs[nIndex_X] += tempvz2;
                        ws[nIndex_X] += - tempuxz * vvs2_dtz_dtx;
                    }
                }
            }

            // Original
            // for(k=ntop;k<nbottom;k++)
            // {
            //     vvp2=vpp2(k);
            //     vvs2=vss2(k);

            //     vvs2_dtz_dtz = vvs2*dtz*dtz;
            //     vvp2_dtx_dtx = vvp2*dtx*dtx;
            //     vvs2_dtx_dtx = vvs2*dtx*dtx;
            //     vvp2_dtz_dtz = vvp2*dtz*dtz;
            //     vvp2_dtz_dtx = vvp2*dtz*dtx;
            //     vvs2_dtz_dtx = vvs2*dtz*dtx;

            //     for(j=nfront;j<nback;j++)
            //         for(i=nleft;i<nright;i++)
            //         {
            //             int nIndex              = POSITION_INDEX(k,j,i);

            //             if(i==ncx_shot-1&&j==ncy_shot-1&&k==ncz_shot-1)
            //             {
            //                 px=1.;
            //             }
            //             else
            //             {
            //                 px=0.;
            //             }

            //             tempux2=0.0f;
            //             tempuy2=0.0f;
            //             tempuz2=0.0f;
            //             tempvx2=0.0f;
            //             tempvy2=0.0f;
            //             tempvz2=0.0f;
            //             tempwx2=0.0f;
            //             tempwy2=0.0f;
            //             tempwz2=0.0f;
            //             tempuxz=0.0f;
            //             tempuxy=0.0f;
            //             tempvyz=0.0f;
            //             tempvxy=0.0f;
            //             tempwxz=0.0f;
            //             tempwyz=0.0f;

            //             for(kk=1;kk<=mm;kk++)
            //             {
            //                 tempux2=tempux2+c[kk-1][0]*(u[POSITION_INDEX(k,j,i+kk)]+u[POSITION_INDEX(k,j,i-kk)]);
            //                 tempuy2=tempuy2+c[kk-1][0]*(u[POSITION_INDEX(k,j+kk,i)]+u[POSITION_INDEX(k,j-kk,i)]);
            //                 tempuz2=tempuz2+c[kk-1][0]*(u[POSITION_INDEX(k+kk,j,i)]+u[POSITION_INDEX(k-kk,j,i)]);

            //                 tempvx2=tempvx2+c[kk-1][0]*(v[POSITION_INDEX(k,j,i+kk)]+v[POSITION_INDEX(k,j,i-kk)]);
            //                 tempvy2=tempvy2+c[kk-1][0]*(v[POSITION_INDEX(k,j+kk,i)]+v[POSITION_INDEX(k,j-kk,i)]);
            //                 tempvz2=tempvz2+c[kk-1][0]*(v[POSITION_INDEX(k+kk,j,i)]+v[POSITION_INDEX(k-kk,j,i)]);

            //                 tempwx2=tempwx2+c[kk-1][0]*(w[POSITION_INDEX(k,j,i+kk)]+w[POSITION_INDEX(k,j,i-kk)]);
            //                 tempwy2=tempwy2+c[kk-1][0]*(w[POSITION_INDEX(k,j+kk,i)]+w[POSITION_INDEX(k,j-kk,i)]);
            //                 tempwz2=tempwz2+c[kk-1][0]*(w[POSITION_INDEX(k+kk,j,i)]+w[POSITION_INDEX(k-kk,j,i)]);
            //             } //for(kk=1;kk<=mm;kk++) end

            //             tempux2=(tempux2+c0*u[POSITION_INDEX(k,j,i)])*vvp2_dtx_dtx;
            //             tempuy2=(tempuy2+c0*u[POSITION_INDEX(k,j,i)])*vvs2_dtx_dtx;
            //             tempuz2=(tempuz2+c0*u[POSITION_INDEX(k,j,i)])*vvs2_dtz_dtz;

            //             tempvx2=(tempvx2+c0*v[POSITION_INDEX(k,j,i)])*vvs2_dtx_dtx;
            //             tempvy2=(tempvy2+c0*v[POSITION_INDEX(k,j,i)])*vvp2_dtx_dtx;
            //             tempvz2=(tempvz2+c0*v[POSITION_INDEX(k,j,i)])*vvs2_dtz_dtz;

            //             tempwx2=(tempwx2+c0*w[POSITION_INDEX(k,j,i)])*vvs2_dtx_dtx;
            //             tempwy2=(tempwy2+c0*w[POSITION_INDEX(k,j,i)])*vvs2_dtx_dtx;
            //             tempwz2=(tempwz2+c0*w[POSITION_INDEX(k,j,i)])*vvp2_dtz_dtz;

            //             for(kk=1;kk<=mm;kk++)
            //             {
            //                 for(kkk=1;kkk<=mm;kkk++)
            //                 {
            //                     current_c = c[kkk-1][1+kk];

            //                     _tempuxz = POSITION_VALUE(k+kkk,j,i+kk,u);
            //                     _tempwxz = POSITION_VALUE(k+kkk,j,i+kk,w);

            //                     _tempuxy = POSITION_VALUE(k,j+kkk,i+kk,u);
            //                     _tempvxy = POSITION_VALUE(k,j+kkk,i+kk,v);

            //                     _tempvyz = POSITION_VALUE(k+kkk,j+kk,i,v);
            //                     _tempwyz = POSITION_VALUE(k+kkk,j+kk,i,w);


            //                     _tempuxz -= POSITION_VALUE(k-kkk,j,i+kk,u);
            //                     _tempwxz -= POSITION_VALUE(k-kkk,j,i+kk,w);

            //                     _tempuxy -= POSITION_VALUE(k,j-kkk,i+kk,u);
            //                     _tempvxy -= POSITION_VALUE(k,j-kkk,i+kk,v);

            //                     _tempvyz -= POSITION_VALUE(k-kkk,j+kk,i,v);
            //                     _tempwyz -= POSITION_VALUE(k-kkk,j+kk,i,w);


            //                     _tempuxz += POSITION_VALUE(k-kkk,j,i-kk,u);
            //                     _tempwxz += POSITION_VALUE(k-kkk,j,i-kk,w);

            //                     _tempuxy += POSITION_VALUE(k,j-kkk,i-kk,u);
            //                     _tempvxy += POSITION_VALUE(k,j-kkk,i-kk,v);

            //                     _tempvyz += POSITION_VALUE(k-kkk,j-kk,i,v);
            //                     _tempwyz += POSITION_VALUE(k-kkk,j-kk,i,w);


            //                     _tempuxz -= POSITION_VALUE(k+kkk,j,i-kk,u);
            //                     _tempwxz -= POSITION_VALUE(k+kkk,j,i-kk,w);

            //                     _tempuxy -= POSITION_VALUE(k,j+kkk,i-kk,u);
            //                     _tempvxy -= POSITION_VALUE(k,j+kkk,i-kk,v);

            //                     _tempvyz -= POSITION_VALUE(k+kkk,j-kk,i,v);
            //                     _tempwyz -= POSITION_VALUE(k+kkk,j-kk,i,w);

            //                     tempuxz = tempuxz + (current_c*_tempuxz);
            //                     tempwxz = tempwxz + (current_c*_tempwxz);

            //                     tempuxy = tempuxy + (current_c*_tempuxy);
            //                     tempvxy = tempvxy + (current_c*_tempvxy);

            //                     tempvyz = tempvyz + (current_c*_tempvyz);
            //                     tempwyz = tempwyz + (current_c*_tempwyz);

            //                 } // for(kkk=1;kkk<=mm;kkk++) end
            //             } //for(kk=1;kk<=mm;kk++) end

            //             // DEBUG MODULE
            //             

            //             
            //             
            //             // DEBUG MODULE END

            //             POSITION_VALUE(k,j,i,up)=2.*POSITION_VALUE(k,j,i,up1)-POSITION_VALUE(k,j,i,up2) +tempux2+tempwxz*vvp2_dtz_dtx +tempvxy*vvp2_dtz_dtx;

            //             POSITION_VALUE(k,j,i,vp)=2.*POSITION_VALUE(k,j,i,vp1)-POSITION_VALUE(k,j,i,vp2) +tempvy2+tempuxy*vvp2_dtz_dtx +tempwyz*vvp2_dtz_dtx;

            //             POSITION_VALUE(k,j,i,wp)=2.*POSITION_VALUE(k,j,i,wp1)-POSITION_VALUE(k,j,i,wp2) +tempwz2+tempuxz*vvp2_dtz_dtx +tempvyz*vvp2_dtz_dtx +px*wave[l-1];

            //             POSITION_VALUE(k,j,i,us)=2.*POSITION_VALUE(k,j,i,us1)-POSITION_VALUE(k,j,i,us2) +tempuy2+tempuz2 -tempvxy*vvs2_dtz_dtx-tempwxz*vvs2_dtz_dtx;

            //             POSITION_VALUE(k,j,i,vs)=2.*POSITION_VALUE(k,j,i,vs1)-POSITION_VALUE(k,j,i,vs2) +tempvx2+tempvz2 -tempuxy*vvs2_dtz_dtx-tempwyz*vvs2_dtz_dtx;

            //             POSITION_VALUE(k,j,i,ws)=2.*POSITION_VALUE(k,j,i,ws1)-POSITION_VALUE(k,j,i,ws2) +tempwx2+tempwy2 -tempuxz*vvs2_dtz_dtx-tempvyz*vvs2_dtz_dtx;
            //         }//for(i=nleft;i<nright;i++) end
            // }//for(k)

            for(k=ntop;k<nbottom;k++)
                for(j=nfront;j<nback;j++)
                    for(i=nleft;i<nright;i++)
                    {
                        int nIndex              = POSITION_INDEX(k,j,i);

                        up[nIndex] += 2 * up1[nIndex] - up2[nIndex];
                        vp[nIndex] += 2 * vp1[nIndex] - vp2[nIndex];
                        wp[nIndex] += 2 * wp1[nIndex] - wp2[nIndex];
                        us[nIndex] += 2 * us1[nIndex] - us2[nIndex];
                        vs[nIndex] += 2 * vs1[nIndex] - vs2[nIndex];
                        ws[nIndex] += 2 * ws1[nIndex] - ws2[nIndex];

                        u_x[POSITION_INDEX(k,j,i)] = up[POSITION_INDEX(k,j,i)] + us[POSITION_INDEX(k,j,i)];
                        v_x[POSITION_INDEX(k,j,i)] = vp[POSITION_INDEX(k,j,i)] + vs[POSITION_INDEX(k,j,i)];
                        w_x[POSITION_INDEX(k,j,i)] = wp[POSITION_INDEX(k,j,i)] + ws[POSITION_INDEX(k,j,i)];

                        u_z[POSITION_INDEX_Z(k,j,i)] = u_y[POSITION_INDEX_Y(k,j,i)] = u_x[POSITION_INDEX(k,j,i)];
                        v_z[POSITION_INDEX_Z(k,j,i)] = v_y[POSITION_INDEX_Y(k,j,i)] = v_x[POSITION_INDEX(k,j,i)];
                        w_z[POSITION_INDEX_Z(k,j,i)] = w_y[POSITION_INDEX_Y(k,j,i)] = w_x[POSITION_INDEX(k,j,i)];

                    }//for(i=nleft;i<nright;i++) end

            float *swap_temp;
            swap_temp = up2; up2 = up1; up1 = up; up = swap_temp;
            swap_temp = vp2; vp2 = vp1; vp1 = vp; vp = swap_temp;
            swap_temp = wp2; wp2 = wp1; wp1 = wp; wp = swap_temp;
            swap_temp = us2; us2 = us1; us1 = us; us = swap_temp;
            swap_temp = vs2; vs2 = vs1; vs1 = vs; vs = swap_temp;
            swap_temp = ws2; ws2 = ws1; ws1 = ws; ws = swap_temp;
        }//for(l=1;l<=lt;l++) end

        fwrite(up1+169*nSliceSize,sizeof(float),nSliceSize,fout);
    }//for(ishot=1;ishot<=nshot;ishot++) end
    fclose(fout);

    free(u);
    free(v);
    free(w);
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
    free(up_out);

    gettimeofday(&end,NULL);
    all_time = (end.tv_sec-start.tv_sec)+(float)(end.tv_usec-start.tv_usec)/1000000.0;
    printf("run time:\t%f s\n",all_time);
    flog = fopen(logfile,"a");
    fprintf(flog,"\nrun time:\t%f s\n\n",all_time);
    fclose(flog);
    flog = fopen(logfile,"a");
    fprintf(flog,"------------end time------------\n");
    fclose(flog);
    system(tmp);
    return 1;
}

