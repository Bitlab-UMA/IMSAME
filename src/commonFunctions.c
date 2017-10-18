#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>
#include <emmintrin.h>
#include <math.h>
#include "structs.h"

#define ALIGNED_BYTES 16

void terror(char *s) {
    fprintf(stdout, "ERR**** %s ****\n", s);
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

void get_num_seqs_and_length(char * seq_buffer, uint64_t * n_seqs, uint64_t * t_len, LoadingDBArgs * ldbargs){
    char check[16] = ">>>>>>>>>>>>>>>>";
    uint64_t a_fourth = *t_len / 4;

    uint64_t i=0;
    
    __m128i * ptr1;
    __m128i * ptr2;
    __m128i res = _mm_setzero_si128();
    ptr2 = (__m128i *) check;

    ldbargs[0].read_from = 0;
    for(i=0; i<a_fourth - 16; i+=16){
        ptr1 = (__m128i *) &seq_buffer[i];
        res = _mm_cmpeq_epi8(*ptr1, *ptr2);
        // TODO parallelize this
        *n_seqs += (uint64_t) (((unsigned char *) &res)[0]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[1]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[2]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[3]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[4]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[5]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[6]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[7]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[8]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[9]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[10]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[11]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[12]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[13]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[14]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[15]) >> 7;
    }
    while(seq_buffer[i] != '>'){ ++i; }
    ldbargs[0].read_to = i;
    ldbargs[0].data_database->n_seqs = *n_seqs;
    ldbargs[1].read_from = i;
    while(i % 16 != 0){ ++i; if(seq_buffer[i] == '>') terror("Sequence of length < 16"); }
    if(i > ldbargs[1].read_from) ++(*n_seqs); // For the advance to align bytes

    //printf("seqs at 0: %"PRIu64" readfrom: %"PRIu64", readto: %"PRIu64", \n", ldbargs[0].data_database->n_seqs, ldbargs[0].read_from, ldbargs[0].read_to);

    
    for(i=i; i<2*(a_fourth) - 16; i+=16){
        ptr1 = (__m128i *) &seq_buffer[i];
        res = _mm_cmpeq_epi8(*ptr1, *ptr2);
        // TODO parallelize this
        *n_seqs += (uint64_t) (((unsigned char *) &res)[0]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[1]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[2]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[3]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[4]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[5]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[6]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[7]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[8]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[9]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[10]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[11]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[12]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[13]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[14]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[15]) >> 7;
    }
    while(seq_buffer[i] != '>') ++i;
    ldbargs[1].read_to = i;
    ldbargs[1].data_database->n_seqs = *n_seqs - ldbargs[0].data_database->n_seqs;
    ldbargs[2].read_from = i;
    while(i % 16 != 0){ ++i; if(seq_buffer[i] == '>') terror("Sequence of length < 16"); }
    if(i > ldbargs[2].read_from) ++(*n_seqs); // For the advance to align bytes

    //printf("seqs at 1: %"PRIu64" readfrom: %"PRIu64", readto: %"PRIu64", \n", ldbargs[1].data_database->n_seqs, ldbargs[1].read_from, ldbargs[1].read_to);

    
    for(i=i; i<3*(a_fourth) - 16; i+=16){
        ptr1 = (__m128i *) &seq_buffer[i];
        res = _mm_cmpeq_epi8(*ptr1, *ptr2);
        // TODO parallelize this
        *n_seqs += (uint64_t) (((unsigned char *) &res)[0]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[1]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[2]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[3]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[4]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[5]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[6]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[7]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[8]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[9]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[10]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[11]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[12]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[13]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[14]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[15]) >> 7;
    }
    while(seq_buffer[i] != '>') ++i;
    ldbargs[2].read_to = i;
    ldbargs[2].data_database->n_seqs = *n_seqs - (ldbargs[1].data_database->n_seqs + ldbargs[0].data_database->n_seqs);
    ldbargs[3].read_from = i;
    while(i % 16 != 0){ ++i; if(seq_buffer[i] == '>') terror("Sequence of length < 16"); }
    if(i > ldbargs[3].read_from) ++(*n_seqs); // For the advance to align bytes

    //printf("seqs at 2: %"PRIu64" readfrom: %"PRIu64", readto: %"PRIu64", \n", ldbargs[2].data_database->n_seqs, ldbargs[2].read_from, ldbargs[2].read_to);

    
    for(i=i; i<(*t_len) - 16; i+=16){
        ptr1 = (__m128i *) &seq_buffer[i];
        res = _mm_cmpeq_epi8(*ptr1, *ptr2);
        // TODO parallelize this
        *n_seqs += (uint64_t) (((unsigned char *) &res)[0]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[1]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[2]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[3]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[4]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[5]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[6]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[7]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[8]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[9]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[10]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[11]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[12]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[13]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[14]) >> 7;
        *n_seqs += (uint64_t) (((unsigned char *) &res)[15]) >> 7;
    }
    while(i < *t_len){
        if(seq_buffer[i] == '>') ++(*n_seqs);
        ++i;
    }
    ldbargs[3].read_to = *t_len;
    ldbargs[3].data_database->n_seqs = *n_seqs - (ldbargs[2].data_database->n_seqs + ldbargs[1].data_database->n_seqs + ldbargs[0].data_database->n_seqs);

    //printf("seqs at 3: %"PRIu64" readfrom: %"PRIu64", readto: %"PRIu64", \n", ldbargs[3].data_database->n_seqs, ldbargs[3].read_from, ldbargs[3].read_to);

}

Queue * generate_queue(Head * queue_head, uint64_t t_reads, uint64_t n_threads, uint64_t levels){
    uint64_t i, j;
    uint64_t reads_per_thread;
    if(levels > t_reads) levels = t_reads;
    uint64_t pieces = t_reads/levels;
    uint64_t from, to, t_queues = 0, current_queue = 0;
    for(i=0;i<levels;i++) t_queues += ((i+1)*n_threads);
    if(t_queues > t_reads) t_queues = t_reads;
    Queue * queues = (Queue *) malloc(t_queues * sizeof(Queue));
    if(queues == NULL) terror("Could not allocate queue tasks");
    queue_head->head = &queues[0];

    
    for(i=0;i<levels;i++){

        //reads_per_thread = (uint64_t) (floorl((long double) pieces / (long double) ((i+1)*n_threads)));
        reads_per_thread = (uint64_t) (ceill((long double) pieces / (long double) ((i+1)*n_threads)));
        

        for(j=0;j<(i+1)*n_threads;j++){
            from = j * reads_per_thread + (pieces*i);
            to = (j + 1) * reads_per_thread + (pieces*i);
            
            if(j == (i+1)*n_threads - 1) to = pieces*(i+1);
		
	    if(from >= t_reads){
		    queues[current_queue-1].next = NULL;
                    return &queues[0];
            }


            if(i == levels - 1 && j == (i+1)*n_threads - 1){
                //If its the last 
                queues[current_queue].next = NULL;
            }else{
                //Else add it to the queue
                queues[current_queue].next = &queues[current_queue+1];
            }
            

            queues[current_queue].r1 = from;
            queues[current_queue].r2 = (to <= t_reads) ? (to) : (t_reads);
            current_queue++;
            //printf("current_piece: %"PRIu64"-%"PRIu64" diff: %"PRIu64"\n", from, to, to - from);

        }

    }
    //printf("TREADS was %"PRIu64"\n", t_reads);    
    queues[current_queue-1].r2 = t_reads;
    return &queues[0];
}

void print_queue(Queue * q){
    fprintf(stdout, "Task: %"PRIu64"-%"PRIu64"\n", q->r1, q->r2);
}

Queue * get_task_from_queue(Head * queue_head, pthread_mutex_t * lock){
    pthread_mutex_lock(lock);

    Queue * ptr = queue_head->head;
    if(queue_head->head != NULL) queue_head->head = queue_head->head->next;
    //if(ptr != NULL){ printf("Taking "); /*print_queue(ptr);*/ }

    pthread_mutex_unlock(lock);


    return ptr;
}

uint64_t quick_pow4(uint64_t n){
    static uint64_t pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L, 
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L, 
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L, 
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L, 
    72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};
    return pow4[n];
}

uint64_t quick_pow4byLetter(uint64_t n, const char c){
    static uint64_t pow4_G[33]={2*1L, 2*4L, 2*16L, 2*64L, 2*256L, 2*1024L, 2*4096L, 2*16384L, 2*65536L, 
    (uint64_t)2*262144L, (uint64_t)2*1048576L,(uint64_t)2*4194304L, (uint64_t)2*16777216L, (uint64_t)2*67108864L, (uint64_t)2*268435456L, (uint64_t)2*1073741824L, (uint64_t)2*4294967296L, 
    (uint64_t)2*17179869184L, (uint64_t)2*68719476736L, (uint64_t)2*274877906944L, (uint64_t)2*1099511627776L, (uint64_t)2*4398046511104L, (uint64_t)2*17592186044416L, 
    (uint64_t)2*70368744177664L, (uint64_t)2*281474976710656L, (uint64_t)2*1125899906842624L, (uint64_t)2*4503599627370496L, (uint64_t)2*18014398509481984L, 
    (uint64_t)2*72057594037927936L, (uint64_t) 2*288230376151711744L, (uint64_t) 2*1152921504606846976L, (uint64_t) 2*4611686018427387904L};
    
    static uint64_t pow4_T[33]={3*1L, 3*4L, 3*16L, 3*64L, 3*256L, 3*1024L, 3*4096L, 3*16384L, 3*65536L, 
    (uint64_t)3*262144L, (uint64_t) 3*1048576L, (uint64_t)3*4194304L, (uint64_t)3*16777216L, (uint64_t)3*67108864L, (uint64_t)3*268435456L, (uint64_t)3*1073741824L, (uint64_t)3*4294967296L, 
    (uint64_t)3*17179869184L, (uint64_t)3*68719476736L, (uint64_t)3*274877906944L, (uint64_t)3*1099511627776L, (uint64_t)3*4398046511104L, (uint64_t)3*17592186044416L, 
    (uint64_t)3*70368744177664L, (uint64_t)3*281474976710656L, (uint64_t)3*1125899906842624L, (uint64_t)3*4503599627370496L, (uint64_t)3*18014398509481984L, 
    (uint64_t)3*72057594037927936L, (uint64_t) 3*288230376151711744L, (uint64_t) 3*1152921504606846976L, (uint64_t) 3*4611686018427387904L};
    
    if(c == 'A') return 0;
    if(c == 'C') return quick_pow4(n);
    if(c == 'G') return pow4_G[n];
    if(c == 'T') return pow4_T[n];
    return 0;
}

uint64_t hashOfWord(const unsigned char * word, uint32_t k){
    
    uint64_t value = 0, jIdx;
    for(jIdx=0;jIdx<k;jIdx++){
        value += quick_pow4byLetter(k-(jIdx+1), (char) word[jIdx]);
    }
    return value;
}


uint64_t asciiToUint64(const char *text){
	uint64_t number=0;

	for(;*text;text++){
		char digit=*text-'0';           
		number=(number*10)+digit;
	}
	return number;
}
