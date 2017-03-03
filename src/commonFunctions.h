#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H
#include "structs.h"
/**
 * Print the error message 's' and exit(-1)
 */
void terror(char *s);


/**
 * Function to read char by char buffered from a FILE
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

/*
    Generates a queue of tasks for threads
*/
void generate_queue(Head * queue_head, uint64_t t_reads, uint64_t n_threads, uint64_t levels);

#endif /* COMMON_FUNCTIONS_H */
