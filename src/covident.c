/*
for i in $IMS/_distribution; do ./covident $i computed/$(basename $i); done
gcc -O3 -D_FILE_OFFSET_BITS=64 covident.c -o covident
**********/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define MAXLID	200
#define OUT_MATSIZE 101



int main(int argc, char ** av){
	
	if(argc != 3){ printf("USE:covident in output.MAT\n"); exit(-1); }
	
	
	FILE * f, * g;
	f = fopen(av[1], "rt");
	if(f == NULL){
		printf("Error opening input file\n");
		exit(-1);
	}
	g = fopen(av[2], "wt");
	if(g == NULL){
		printf("	Error opening output file %s\n", av[2]);
		exit(-1);
	}
	
	
	uint64_t matrix[OUT_MATSIZE][OUT_MATSIZE];	
	uint64_t i,j;
	for(i=0;i<OUT_MATSIZE;i++){
		for(j=0;j<OUT_MATSIZE;j++){
			matrix[i][j] = 0;
		}
	}
	uint64_t cov, ident, len;
	while(!feof(f)){
		
		if(3 != fscanf(f, "%"PRIu64" %"PRIu64" %"PRIu64, &cov, &ident, &len) && !feof(f)){ printf("Did not read 3\n"); exit(-1);}
		matrix[cov][ident]++;
		
		
	}
	
	
	for(i=0;i<OUT_MATSIZE;i++){
		for(j=0;j<OUT_MATSIZE;j++){
			fprintf(g, "%"PRIu64"\t", matrix[i][j]);
		}
		fprintf(g, "\n");
	}
	
	fclose(f);
	fclose(g);
	return 0;
}




