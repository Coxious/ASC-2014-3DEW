#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#define PIE 3.1415926

#define POSITION_INDEX(_z,_y,_x)        ((_z)*ny*nx + (_y)*nx + (_x))
#define POSITION_VALUE(_z,_y,_x,_attr)  (( pPositionData + ((_z)*ny*nx + (_y)*nx + (_x)))->_attr)


#define vpp2(_z) ((_z<210)?(5290000):(_z>=260?12250000:7840000))
#define vss2(_z) ((_z<210)?(1517824):(_z>=260?3644281:2277081))
 
typedef struct _POSITION_DATA{
    float u  ;
    float v  ;
    float w  ;
    float up ;
    float up1; 
    float up2; 
    float vp ;
    float vp1;
    float vp2; 
    float wp ; 
    float wp1; 
    float wp2; 
    float us ;
    float us1; 
    float us2; 
    float vs ; 
    float vs1; 
    float vs2; 
    float ws ; 
    float ws1; 
    float ws2; 
}POSITION_DATA,*PPOSITION_DATA;

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

    float *vpp, /* *density,**/ *vss ,*up_out;
    float c[5][7];
    float *wave;
    float nshot,t0,tt,c0;
    float dtx,dtz,dtxz;
    float xmax,px,sx;
    float vvp2,drd1,drd2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
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

    PPOSITION_DATA pPositionData = (PPOSITION_DATA)malloc(sizeof(POSITION_DATA)*(long long)(nSize));

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
    dtxz=dtx*dtz;

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

        memset(pPositionData,0,nSize*sizeof(POSITION_DATA));

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

            #pragma omp parallel for
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

                for(j=nfront;j<nback;j++)
                    for(i=nleft;i<nright;i++)
                    {
                        int nIndex              = POSITION_INDEX(k,j,i);
                        PPOSITION_DATA  pCurPos = pPositionData + nIndex;


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
                            tempux2=tempux2+c[kk-1][0]*(POSITION_VALUE(k,j,i+kk,u)+POSITION_VALUE(k,j,i-kk,u));
                            tempuy2=tempuy2+c[kk-1][0]*(POSITION_VALUE(k,j+kk,i,u)+POSITION_VALUE(k,j-kk,i,u));
                            tempuz2=tempuz2+c[kk-1][0]*(POSITION_VALUE(k+kk,j,i,u)+POSITION_VALUE(k-kk,j,i,u));

                            tempvx2=tempvx2+c[kk-1][0]*(POSITION_VALUE(k,j,i+kk,v)+POSITION_VALUE(k,j,i-kk,v));
                            tempvy2=tempvy2+c[kk-1][0]*(POSITION_VALUE(k,j+kk,i,v)+POSITION_VALUE(k,j-kk,i,v));
                            tempvz2=tempvz2+c[kk-1][0]*(POSITION_VALUE(k+kk,j,i,v)+POSITION_VALUE(k-kk,j,i,v));

                            tempwx2=tempwx2+c[kk-1][0]*(POSITION_VALUE(k,j,i+kk,w)+POSITION_VALUE(k,j,i-kk,w));
                            tempwy2=tempwy2+c[kk-1][0]*(POSITION_VALUE(k,j+kk,i,w)+POSITION_VALUE(k,j-kk,i,w));
                            tempwz2=tempwz2+c[kk-1][0]*(POSITION_VALUE(k+kk,j,i,w)+POSITION_VALUE(k-kk,j,i,w));
                        } //for(kk=1;kk<=mm;kk++) end

                        tempux2=(tempux2+c0*pCurPos->u)*vvp2_dtx_dtx;
                        tempuy2=(tempuy2+c0*pCurPos->u)*vvs2_dtx_dtx;
                        tempuz2=(tempuz2+c0*pCurPos->u)*vvs2_dtz_dtz;

                        tempvx2=(tempvx2+c0*pCurPos->v)*vvs2_dtx_dtx;
                        tempvy2=(tempvy2+c0*pCurPos->v)*vvp2_dtx_dtx;
                        tempvz2=(tempvz2+c0*pCurPos->v)*vvs2_dtz_dtz;

                        tempwx2=(tempwx2+c0*pCurPos->w)*vvs2_dtx_dtx;
                        tempwy2=(tempwy2+c0*pCurPos->w)*vvs2_dtx_dtx;
                        tempwz2=(tempwz2+c0*pCurPos->w)*vvp2_dtz_dtz;

                        for(kk=1;kk<=mm;kk++)
                        {
                            for(kkk=1;kkk<=mm;kkk++)
                            {
                                current_c = c[kkk-1][1+kk];

                                _tempuxz = POSITION_VALUE(k+kkk,j,i+kk,u);
                                _tempwxz = POSITION_VALUE(k+kkk,j,i+kk,w);

                                _tempuxy = POSITION_VALUE(k,j+kkk,i+kk,u);
                                _tempvxy = POSITION_VALUE(k,j+kkk,i+kk,v);

                                _tempvyz = POSITION_VALUE(k+kkk,j+kk,i,v);
                                _tempwyz = POSITION_VALUE(k+kkk,j+kk,i,w);


                                _tempuxz -= POSITION_VALUE(k-kkk,j,i+kk,u);
                                _tempwxz -= POSITION_VALUE(k-kkk,j,i+kk,w);

                                _tempuxy -= POSITION_VALUE(k,j-kkk,i+kk,u);
                                _tempvxy -= POSITION_VALUE(k,j-kkk,i+kk,v);

                                _tempvyz -= POSITION_VALUE(k-kkk,j+kk,i,v);
                                _tempwyz -= POSITION_VALUE(k-kkk,j+kk,i,w);


                                _tempuxz += POSITION_VALUE(k-kkk,j,i-kk,u);
                                _tempwxz += POSITION_VALUE(k-kkk,j,i-kk,w);

                                _tempuxy += POSITION_VALUE(k,j-kkk,i-kk,u);
                                _tempvxy += POSITION_VALUE(k,j-kkk,i-kk,v);

                                _tempvyz += POSITION_VALUE(k-kkk,j-kk,i,v);
                                _tempwyz += POSITION_VALUE(k-kkk,j-kk,i,w);


                                _tempuxz -= POSITION_VALUE(k+kkk,j,i-kk,u);
                                _tempwxz -= POSITION_VALUE(k+kkk,j,i-kk,w);

                                _tempuxy -= POSITION_VALUE(k,j+kkk,i-kk,u);
                                _tempvxy -= POSITION_VALUE(k,j+kkk,i-kk,v);

                                _tempvyz -= POSITION_VALUE(k+kkk,j-kk,i,v);
                                _tempwyz -= POSITION_VALUE(k+kkk,j-kk,i,w);

                                tempuxz = tempuxz + (current_c*_tempuxz);
                                tempwxz = tempwxz + (current_c*_tempwxz);

                                tempuxy = tempuxy + (current_c*_tempuxy);
                                tempvxy = tempvxy + (current_c*_tempvxy);

                                tempvyz = tempvyz + (current_c*_tempvyz);
                                tempwyz = tempwyz + (current_c*_tempwyz);

                            } // for(kkk=1;kkk<=mm;kkk++) end
                        } //for(kk=1;kk<=mm;kk++) end
                        pCurPos->up=2.*pCurPos->up1-pCurPos->up2
                                          +tempux2+tempwxz*vvp2_dtz_dtx
                                          +tempvxy*vvp2_dtz_dtx;

                        pCurPos->vp=2.*pCurPos->vp1-pCurPos->vp2
                                          +tempvy2+tempuxy*vvp2_dtz_dtx
                                          +tempwyz*vvp2_dtz_dtx;
                        
                        pCurPos->wp=2.*pCurPos->wp1-pCurPos->wp2
                                          +tempwz2+tempuxz*vvp2_dtz_dtx
                                          +tempvyz*vvp2_dtz_dtx
                                          +px*wave[l-1];
                        
                        pCurPos->us=2.*pCurPos->us1-pCurPos->us2
                                          +tempuy2+tempuz2
                                          -tempvxy*vvs2_dtz_dtx-tempwxz*vvs2_dtz_dtx;
                        
                        pCurPos->vs=2.*pCurPos->vs1-pCurPos->vs2
                                          +tempvx2+tempvz2
                                          -tempuxy*vvs2_dtz_dtx-tempwyz*vvs2_dtz_dtx;
                        
                        pCurPos->ws=2.*pCurPos->ws1-pCurPos->ws2
                                          +tempwx2+tempwy2
                                          -tempuxz*vvs2_dtz_dtx-tempvyz*vvs2_dtz_dtx;
                    }//for(i=nleft;i<nright;i++) end
            }//for(k)
            

            #pragma omp parallel for 
            for(k=ntop;k<nbottom;k++)
                for(j=nfront;j<nback;j++)
                    for(i=nleft;i<nright;i++)
                    {
                        int nIndex              = POSITION_INDEX(k,j,i);
                        PPOSITION_DATA  pCurPos = pPositionData + nIndex;

                        pCurPos->u = pCurPos->up + pCurPos->us;
                        pCurPos->v = pCurPos->vp + pCurPos->vs;
                        pCurPos->w = pCurPos->wp + pCurPos->ws;

                        pCurPos->up2 = pCurPos->up1;
                        pCurPos->up1 = pCurPos->up;

                        pCurPos->us2 = pCurPos->us1;
                        pCurPos->us1 = pCurPos->us;

                        pCurPos->vp2 = pCurPos->vp1;
                        pCurPos->vp1 = pCurPos->vp;

                        pCurPos->vs2 = pCurPos->vs1;
                        pCurPos->vs1 = pCurPos->vs;

                        pCurPos->wp2 = pCurPos->wp1;
                        pCurPos->wp1 = pCurPos->wp;

                        pCurPos->ws2 = pCurPos->ws1;
                        pCurPos->ws1 = pCurPos->ws;

                    }//for(i=nleft;i<nright;i++) end
        }//for(l=1;l<=lt;l++) end
        
        PPOSITION_DATA  pCurPos;
        for(i=0,pCurPos=pPositionData+POSITION_INDEX(169,0,0);i<nSliceSize;++i){
            up_out[i] = pCurPos->up;
            ++pCurPos;
        }
        fwrite(up_out,sizeof(float),nSliceSize,fout);
    }//for(ishot=1;ishot<=nshot;ishot++) end
    fclose(fout);

    free(pPositionData);

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
