# Proposal of the 4th problem: Design and Optimization of 3D-EW

## 1. Introduction to the software
This program is a typical application of intensive calculation in physical called 3D-EW. It’s a kind of wave field extrapolation method to simulate the propagation of elastic wave. In this code ,P-wave and S-wave are simulated separately but interact to each other. So the calculation cannot be separated. 
The code dealing with the problem mainly has three parts:

1. Initialization
1. Data Processing
1. Print Results.

The main time-costing part is the second part. The Data processing section is a Multilevel model. It scan the whole area where the starting point affect and apply the effect to the whole area by three loops stands for x, y and z directions. Then the loop of l extend the area that the prior waves affect, then the loop of shot change the starting points several times. Each iteration of shots is individual which means it does not depend on the prior result.

## 2. Introduction to the platform

Thanks for Inspur providing machines we test our programs mainly on. It' Servers named Inspur NF5280M3 which has two Intel Xeon E5-2692 12 cores CPUs each Server, 16G * 8 DDR3 RAM and a Intel XEON PHI-7110P corprocessor. We used for of them as a cluster to implement MPI parallel calculating.

## 3. The test result of the original program

Beyond the pre-prepared two workloads, we do a little change in 'para1.in', reducing shots to only 1 remains, as 'para0.in'. The test results of orinal programs are as listed:

para0.in | para1.in | para2.in
--- | --- | ---
32.03532 s | 533.644775 s | 9753.830078 s

The original code is pretty coarse, but we still need a guideline to evaluate the effects of our optimization.

## 4. An adjustment of data structure

We designed two ways to adjust the data structure storing waves' information.

1. Chagne the arrays into block structure.
2. Sperate the u,v,w arrays into 3 direction -- x,y,z.

The first changes is fully in consideration of CPU cache miss. The caculations of area is one point by one point adjacently. So if we stored all information about one point together, it wil be extracted one time instead of 21 times read from memory -- RAM read cost a lot.

For more detailed description, we desined a data structure to store information in one 'struct', and then we get those points one by one, use properties instead of arrays.

```c
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
```

That's a intutive optimization and easy understood, it shrinks the time cost by memory allocation, memory addressing and address calculation. It's effective but not notable, the test result is behind.

The second changes is based on the same thought, the loop of kk and kkk used the continues memory block, and it's also the preparation of vectorization. Another adventrue of the second way is that we can use swap of pointers instead of swap value of whole array.

The implementation of this optimization is quite difficult. The calculation is happend on a cube which is stored as a line. The original code use three loops of different directions to traverse the space by generate the index of specific point of the line.

The calculation of a point needs a accumulation of value of points near it in different directions. The accumulation is calculate together in original code. But we seperate it this time, the accumulation happens in three different directions. This may have two benefits.

1. It used to need several multiplications to generate indexes, but the in order accumulation only need one self-incresment.
1. The read and write of RAM become sequential so cache-miss may be reduced.
1. Preparation of vectorization.

The second benefit is also a great optimization, but this chapter is only to compare the optimization of memory structure, so we didn't implement it. It will be recalled on the last chapter.

The most important shortage is that we need to maintain three copy of waves' values.

The test result of two different adjustment in 'para0.in' is as follow:

| Workload | Adjustment 1 | Adjustmen 2 |
| --- | --- | --- |
| **para0.in** | 27.296394 | 22.960428 |
| **para1.in** | 449.614471 | 366.522797 |

It seems that the second way is a little faster than the first. And in consideration of the convinience the second way provided for verctorization, our following test is based on the second way.

## 5. Mult-thread test result

Our cluster provide 2 cpus with 24 cores each pre node, so to make full use of the resource, we need to alter program to adjust the massive parallel processors enviroment. First we try to use OpenMP to parallel our calculation tasks. The only data dependency exist in loop for 'l'. So we add complier directives before the rest of loops. There are the test results of the three workloads.

To adjust the 

| para0.in | para1.in | para2.in |
| --- | --- | --- |
| 16.163277 | 19.320847 | 320.572784 |

The simple test shows that parallel running reducing the runing time notablely, and the accelerate times is nearly linear to cores of CPUs. Thanks to OpenMP, it's easy to implement our code to make full use of mult-thread resources. The key snapshots is as following:

First we need to decouple the use of space, so we use ```calloc``` beneath the loop for shots instead of ```malloc``` at the beginning and ```memset``` it to zero every iterations of shots.

And then, the temporary variables need to be declared in block of loop as local variables instead of global variables because they are all private to the parallel of loop.

Third, the output satement used to be like this:

```c
fwrite(up1+169*ny*nx,sizeof(float),ny*nx,fout);
```
The results are write to disk at the end of the loop of each shot, but as the result of parallel running, th sequential goes disorder. The we can't write the result at then end of the loop of each shot. We applied a space called ```up_out``` and copy the result from ```up``` to ```up_out``` begin with where is should start at the end of the loop of each shot. The statements look like this:

```c
for(int ishot = 1; ishot <= nshot; ishot++) {
    ... // content of loop of shots.
    memcpy(up_out+(ishot-1)*nSliceSize, up1+169*nSliceSize,nSliceSize*sizeof(float));
    ... // ebd of loop of shots.
}
fwrite(up_out,sizeof(float),nshot*nSliceSize,fout);
fclose(fout);
```

And the fatal change of code to adjust the mult-processor enviroment is directive instructions of OpenMP, we a statement above the loop of shot:

```c
#pragma omp parallel for private(u_x, v_x, w_x, u_y, v_y, w_y, u_z, v_z, w_z,  up, up1, up2, vp, vp1, vp2, wp, wp1, wp2, us, us1, us2, vs, vs1, vs2, ws, ws1, ws2, ishot, i, j, k, kk, kkk)
for(int ishot=1; ishot <= nshot; ishot++) {...}
```

Finished add `-openmp` switch of compilier, the program now will distribute the loop of shots to different processors.

**Then** we tried another standard and allowed way, **MPI**, to implement parallel calculation. MPI is not as simple as OpenMP. It needs to do some allocation of tasks, and it's hard to use shared-memory to communicating each process. So some specific function need to use to gather the answer.
 
Compared with many-to-one model in parallel transimission, we used serial transimission to gather results every processors generate. It's in consideration of that the only need of parallel transimission is to gather resutls which only happened at the end of the program. So parallel transimission doesn't have too much adnvantages and optimization, and the limit of cache may cause some overflow problems.

There are results of three test workload tested in one node:

| para0.in | para1.in | para2.in |
| --- | --- | --- |
| - | 19.245184 | 316.847473 |

First we need to initialize MPI by calling specific functions like this:

```c
char processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
MPI_Comm_size(MPI_COMM_WORLD, &mpi_sum);
printf("Hello, I'm %d, MPI starting!\n", mpi_id);
```

And the part of initialization of origin code like reading files will be done only in the main process which ```mpi_id``` is **0**. And pass the arguments and distribute tasks by sending messages. I used serial sending rather than boardcast because the data is short and only need to do only once.

The rest of code is similar to original code except the part of output. The result need to be passed to the main process. To avoid the message clogging and stack overflowing, each process except the main will waiting for a signal from the main process, then it will send the result back to the main process. The main process will send signal to each process in order once it finished its work. Aiming to reduce waiting time, the main process will be assigned less tasks.

The data gather snapshots is as following:

```c
if ( mpi_id == 0 ) {
	fout=fopen(outfile,"wb");
	float *buffer = (float *) malloc(sizeof(float)*nSliceSize*(buf[1]-buf[0]+3));
	float *ans = (float *) calloc(nshot*nSliceSize, sizeof(float));
	memcpy(ans+nSliceSize*(buf[0]-1), to_write, nSliceSize*(buf[1]-buf[0]+1)*sizeof(float));

	int iii = 0, flag = 1;
	for(iii = 1; iii < mpi_sum; iii++) {
		MPI_Send(&flag, 1, MPI_INT, iii, 2, MPI_COMM_WORLD);
		MPI_Recv(buffer, nSliceSize*(buf[1]-buf[0]+3), MPI_FLOAT, iii, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		int length, source = status.MPI_TAG;
		MPI_Get_count(&status, MPI_FLOAT, &length);
		memcpy(ans+nSliceSize*(source-1), buffer, length*sizeof(float));
	}
	
	fwrite(ans,sizeof(float),nshot*nSliceSize,fout);
	fclose(fout);

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

	MPI_Finalize();
	free(ans);
	free(buffer);
} else {
	int flag = 0;
	MPI_Recv(&flag, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
	if ( flag ) {
		MPI_Send(to_write, nSliceSize*(buf[1]-buf[0]+1), MPI_FLOAT, 0, buf[0], MPI_COMM_WORLD);
	}
	MPI_Finalize();
}
```

The code need to be complied by the command on the server;

```sh
# mpicxx 3d2rewq.cpp -o 3d2rewq -O3 -openmp
```

And run with command:

```sh
# mpirun -n 24 ./3d2rewq input output log
```

## 6. Implementation of cluster by MPI

Once we finished the implementation of MPI before, the implementation of cluster is such a easy thing -- all we need to do is change the running parameters. 

There are the test results of 'para1.in' and 'para2.in':

| para1.in | para2.in |
| --- | --- |
| 0 | 0 |

We don't provide the test resutl of 'para0.in' because there is only one shot need to be calculate in it which means no parallel calculation happend.

## 7. Acceleration in MIC platform



## 8. Vectorrization



## 9. Details of Further Optimizations

This chapter will show the details of optimizations, some of them may didn't provide notable promotability, but all of them together will have big effect on large workload.

