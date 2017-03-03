#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>
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
    uint64_t pieces = t_reads/levels;
    uint64_t from, to, t_queues = 0, current_queue = 0;
    for(i=0;i<levels;i++) t_queues += ((i+1)*n_threads);
    Queue * queues = (Queue *) malloc(t_queues * sizeof(Queue));
    queue_head->head = &queues[0];
    
    for(i=0;i<levels;i++){

        reads_per_thread = (uint64_t) (floorl((long double) pieces / (long double) ((i+1)*n_threads)));
        

        for(j=0;j<(i+1)*n_threads;j++){
            from = j * reads_per_thread + (pieces*i);
            to = (j + 1) * reads_per_thread + (pieces*i);
            
            if(i==levels - 1 && j == (i+1)*n_threads - 1){
                to = t_reads;
                queues[current_queue].next = NULL;
            }else{
                queues[current_queue].next = &queues[current_queue+1];
            }

            queues[current_queue].r1 = from;
            queues[current_queue].r2 = to;
            current_queue++;
            //printf("current_piece: %"PRIu64"-%"PRIu64"\n", from, to);

        }

    }
    //printf("TREADS was %"PRIu64"\n", t_reads);    
}

Queue * get_task_from_queue(Head * queue_head, pthread_mutex_t * lock){
    pthread_mutex_lock(lock);
    
    printf("Mutex was locked\n");
    printf("Mutex was UN-locked\n");

    pthread_mutex_unlock(lock);

    return NULL;
}