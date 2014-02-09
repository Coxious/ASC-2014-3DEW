#include <iostream>

int main() {

	long *a;
	long out_p[200];

	#pragma offload target(mic:0) nocopy(a:length(2000) alloc_if(1) free_if(0))
	#pragma omp parallel for
	for(int i=0; i< 2000; i++) {
		a[i] = 1;
	}

	for(int j=0; j<100; j++) {
	
		#pragma offload target(mic:0) out(out_p:length(200)) nocopy(a:length(2000) alloc_if(0) free_if(1))
		{
			#pragma omp parallel for
			for( int i=0; i < 2000; i++) 
				a[i]*=i;
		
			for(int i=0; i<2000; i+=10)
				out_p[i/10] = a[i];
		}
		for(int i=0; i<200; i++) cout << out_p[i] << ' ';
		cout << endl;
	
	}
	

	return 0;
}