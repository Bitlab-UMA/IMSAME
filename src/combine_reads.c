#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes
#define MATVAL 101

// void terror(char *s) {
//     fprintf(stdout, "ERR**** %s ****\n", s);
//     exit(-1);
// }


int main(int argc, char ** av){
    
    if(argc != 2) terror("USE:: combine_reads <file>");

	
    char names[10][100] = {"14np", "47np", "BR0716C", "BR11044P", "CA01044C", "CA01067C", "CA11149P", "CA11154P", "MA07066C", "MA09032P"};
    uint64_t llen[10] = {1758451, 959020, 1616163, 884193, 934256, 1007436, 1670521, 584151, 671010, 516632};
    uint64_t nseqs = 10;

    FILE * results, * data;
    data = fopen(av[1], "rt");
    if(data == NULL) terror("Could not open input file");

    uint64_t * mat = (uint64_t *) calloc(MATVAL, sizeof(uint64_t));
    if(mat == NULL) terror("Could not allocate matrix array");
    long double * mat_e = (long double *) calloc(MATVAL, sizeof(long double));
    if(mat_e == NULL) terror("Could not allocate float table");

    char current_file[100];
    uint64_t i=5, idx=0;
    while(av[1][i] != '_' && av[1][i+1] != 'g'){
	    current_file[idx++] = av[1][i];
	    i++;
    }
    current_file[idx] = '\0';
    fprintf(stdout, "Current file is %s\n", current_file);

    char buffer[MAXLID];
    if ((results = fopen("accu.log", "r")) == NULL){
        results = fopen("accu.log", "wt");
    }else{
        // Load the matrix
       
        
        for(i=0;i<100;i++){
            if(0 == fgets(buffer, MAXLID, results)) terror("Missing number on load");
            
            //fprintf(stdout, "Have %s\n", buffer);
            buffer[strlen(buffer)-1] = '\0';
            //mat[i] = asciiToUint64(buffer);
	    mat_e[i] = (long double) atof(buffer);
            //fprintf(stdout, "%"PRIu64"\n", mat[i]);
            //getchar();
        }
        fclose(results);
        results = fopen("accu.log", "wt"); // Re open
    }

    // Read file 
    uint64_t read_id_1, read_id_2, coverage, identity, length, current, currmax, j, last = 100000, lastmax = 0xFFFFFFFFFFFFFFFF;
    while(!feof(data)){
        if(0 == fgets(buffer, MAXLID, data) && !feof(data)) terror("Missing values");
        // 2 77277 89 64 213
        sscanf(buffer, "%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64, &read_id_1, &read_id_2, &coverage, &identity, &length);
        //fprintf(stdout, "%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n", read_id_1, read_id_2, coverage, identity, length);
        currmax = MIN(identity, coverage);
        //fprintf(stdout, "%"PRIu64"\n", currmax);
        current = read_id_1;
        /*
        for(j=currmax; j > 1; j--){
            if(current != lasts[j]){
                mat[j]++;
                lasts[j] = current;
            }
        }
        */
	if(lastmax == 0xFFFFFFFFFFFFFFFF){
		last = current;
		lastmax = 0;
	}
	if(current == last){
		//mat[currmax] += 1;
		if(lastmax < currmax){
			lastmax = currmax;
		}
	}else{
		mat[lastmax]++;
		last = current;
		lastmax = 0;
	}

    }

    uint64_t divider = 0;
    for(i=0; i < nseqs; i++){

	    if(strcmp(current_file, names[i]) == 0){
		    divider = llen[i];
	    }
    }
    if(divider == 0) terror("Could not find file for divider"); else fprintf(stdout, "The divider is %"PRIu64"\n", divider);

    for(j=99; j>0; j--){
        mat[j] += mat[j+1];
    }
    mat[0] = mat[1];

    for(j=0; j<100; j++){
        fprintf(stdout, "%"PRIu64" --> %.5f\n", mat[j], (float)((100*(long double)mat[j])/(long double)divider));
	
        fprintf(results, "%.5f\n", (float)(mat_e[j] + ((100*(long double)mat[j])/(long double)divider))/2);
    }


    fclose(results);
    fclose(data);
    free(mat);
    return 0;
}
