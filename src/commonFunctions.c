#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"

void terror(char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

void generate_queue(Head * queue_head, uint64_t t_reads, uint64_t n_threads, uint64_t levels){
    uint64_t i, j, k;
    uint64_t reads_per_thread;
    uint64_t current_piece = t_reads;
    uint64_t from, to;
    
    for(i=0;i<levels;i++){

        current_piece = t_reads/(powl(2, i+1));
        printf("current_piece: %"PRIu64"\n", current_piece);
        
        /*
        reads_per_thread = (uint64_t) (floorl((long double) current_piece / (long double) n_threads));
        for(j=0;j<n_threads;j++){

            
            from = j * reads_per_thread;
            to = (j + 1) * reads_per_thread;

            if(j==n_threads-1) to = t_reads;

            //Add to queue

        }
        */
    }
    
}