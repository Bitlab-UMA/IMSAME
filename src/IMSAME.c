/*********

File        IMSAME.c
Author      EPW <estebanpw@uma.es>
Description Computes an incremental alignment on reads versus reads using n threads

USAGE       Usage is described by calling ./IMSAME --help



**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes


void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t  * n_threads, long double * minevalue, long double * mincoverage, int * igap, int * egap, long double * minidentity, long double * window);

int VERBOSE_ACTIVE = 0;

int main(int argc, char ** av){
    

    clock_t begin, end;
    
    int error; //To tell if threads could not be launched
    uint64_t i,j;

    //query to read kmers from, database to find seeds
    FILE * query = NULL, * database = NULL, * out_database = NULL;
    long double minevalue = 1/powl(10, 20); //Default 1 * 10^-20
    
    long double mincoverage = 0.5, minidentity = 0.5, window = 0.15; //Default
    int igap = -5, egap = -2;

    uint64_t n_threads = 4;
    init_args(argc, av, &query, &database, &out_database, &n_threads, &minevalue, &mincoverage, &igap, &egap, &minidentity, &window);
    
    //uint64_t reads_per_thread;
    uint64_t sum_accepted_reads = 0;

    unsigned char char_converter[91];
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;

    begin = clock();

    fprintf(stdout, "[INFO] Init. quick table\n");

    pthread_t * threads = (pthread_t *) malloc(n_threads * sizeof(pthread_t));
    if(threads == NULL) terror("Could not create threads");


    HashTableArgs * hta = (HashTableArgs *) malloc(n_threads*sizeof(HashTableArgs));
    if(hta == NULL) terror("Could not allocate arguments for hash table");
    
    pthread_mutex_t lock; //The mutex to lock the queue
    if (pthread_mutex_init(&lock, NULL) != 0) terror("Could not init mutex");

    unsigned char ** my_x = (unsigned char **) malloc(n_threads * sizeof(unsigned char*));
    unsigned char ** my_y = (unsigned char **) malloc(n_threads * sizeof(unsigned char*));

    struct positioned_cell ** mc = (struct positioned_cell **) malloc(n_threads * sizeof(struct positioned_cell *));
    struct cell *** table = (struct cell ***) malloc(n_threads * sizeof(struct cell **));
    if(table == NULL) terror("Could not allocate NW table");
    char ** reconstruct_X = (char **) malloc(n_threads * sizeof(char *));
    char ** reconstruct_Y = (char **) malloc(n_threads * sizeof(char *));
    if(reconstruct_Y == NULL || reconstruct_X == NULL) terror("Could not allocate output alignment sequences");
    char ** writing_buffer_alignment = (char **) malloc(n_threads * sizeof(char*));
    for(i=0;i<n_threads;i++){
	
    	table[i] = (struct cell **) malloc(MAX_READ_SIZE * sizeof(struct cell *));
	    for(j=0;j<MAX_READ_SIZE;j++){
		    table[i][j] = (struct cell *) malloc(MAX_READ_SIZE*sizeof(struct cell));
		    if(table[i][j] == NULL) terror("Could not allocate memory for second loop of table");
	    }
    	mc[i] = (struct positioned_cell *) malloc(MAX_READ_SIZE * sizeof(struct positioned_cell));
    	my_x[i] = (unsigned char *) malloc(MAX_READ_SIZE * sizeof(unsigned char));
    	my_y[i] = (unsigned char *) malloc(MAX_READ_SIZE * sizeof(unsigned char));
    	reconstruct_X[i] = (char *) malloc(2*MAX_READ_SIZE * sizeof(char));
    	reconstruct_Y[i] = (char *) malloc(2*MAX_READ_SIZE * sizeof(char));
	    writing_buffer_alignment[i] = (char *) malloc(MAX_READ_SIZE*MAX_READ_SIZE*sizeof(char));
	    if(table[i] == NULL || mc[i] == NULL || my_x[i] == NULL || my_y[i] == NULL || reconstruct_X[i] == NULL || reconstruct_Y[i] == NULL || writing_buffer_alignment[i] == NULL) terror("Could not allocate buffer for alignment output");
    }



    end = clock();
    fprintf(stdout, "[INFO] Initialization took %e seconds \n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Variables to account for positions
    //Print info
    fprintf(stdout, "[INFO] Loading database\n");
    //Variables to read kmers
    char c = 'N'; //Char to read character
    //Current length of array and variables for the buffer
    uint64_t idx = 0, r = 0;
    
    //Vector to read in batches
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }
    //Vector to store database seq
    unsigned char * seq_vector_database = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_database == NULL) terror("Could not allocate memory for database vector");
    uint64_t n_realloc_database = 1;
    uint64_t * database_positions = (uint64_t *) malloc(INITSEQS*sizeof(uint64_t));
    if(database_positions == NULL) terror("Could not allocate database sequences positions");

    unsigned char curr_kmer[FIXED_K];
    curr_kmer[0] = '\0';
    uint64_t word_size = 0, pos_in_database = 0, n_seqs_database_realloc = 1;

    //To hold all information related to database
    SeqInfo data_database;
    SeqInfo data_query;
    data_database.sequences = seq_vector_database;
    data_database.start_pos = database_positions;
    data_database.total_len = 0;
    data_database.n_seqs = 0;
    
    //To force reading from the buffer
    idx = READBUF + 1;

    //Store positions of kmers
    uint64_t n_pools_used = 0;
    //Mempool_l * mp = (Mempool_l *) malloc(MAX_MEM_POOLS*sizeof(Mempool_l));
    //if(mp == NULL) terror("Could not allocate vectors for memory pools");
    Mempool_l mp[MAX_MEM_POOLS];
    init_mem_pool_llpos(&mp[n_pools_used]);
    llpos * aux, * pointer;

    unsigned char aux_kmer[FIXED_K+1];
    
    //Vector to store query seq
    unsigned char * seq_vector_query = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_query == NULL) terror("Could not allocate memory for query vector");
    uint64_t n_realloc_query = 1, pos_in_query = 0, n_seqs_query_realloc = 1;
    uint64_t * query_positions = (uint64_t *) malloc(INITSEQS*sizeof(uint64_t));
    if(query_positions == NULL) terror("Could not allocate query sequences positions");


    Container * ct = (Container *) calloc(1, sizeof(Container));
    if(ct == NULL) terror("Could not allocate container");
    
    /*
    uint64_t w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11;
    for(w0=0;w0<4;w0++){
        for(w1=0;w1<4;w1++){
            for(w2=0;w2<4;w2++){
                for(w3=0;w3<4;w3++){
                    for(w4=0;w4<4;w4++){
                        for(w5=0;w5<4;w5++){
                            for(w6=0;w6<4;w6++){
                                for(w7=0;w7<4;w7++){
                                    for(w8=0;w8<4;w8++){
                                        for(w9=0;w9<4;w9++){
                                            for(w10=0;w10<4;w10++){
                                                for(w11=0;w11<4;w11++){
                                                    ct->table[w0][w1][w2][w3][w4][w5][w6][w7][w8][w9][w10][w11] = NULL;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    */







    begin = clock();

    c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
    while((!feof(database) || (feof(database) && idx < r))){

        if(c == '>'){
            data_database.start_pos[data_database.n_seqs++] = pos_in_database;
            
            if(pos_in_database == READBUF*n_realloc_database){ 
                n_realloc_database++; data_database.sequences = (unsigned char *) realloc(data_database.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                if(data_database.sequences == NULL) terror("Could not reallocate temporary database");
            }

            if(data_database.n_seqs == INITSEQS*n_seqs_database_realloc){
                n_seqs_database_realloc++; data_database.start_pos =  (uint64_t *) realloc(data_database.start_pos, INITSEQS*n_seqs_database_realloc*sizeof(uint64_t));
            }


            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);  //Skip ID
                

            while(c != '>' && (!feof(database) || (feof(database) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < FIXED_K) word_size++;
                    data_database.sequences[pos_in_database++] = (unsigned char) c;
            
                    if(pos_in_database == READBUF*n_realloc_database){ 
                        n_realloc_database++; data_database.sequences = (unsigned char *) realloc(data_database.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                        if(data_database.sequences == NULL) terror("Could not reallocate temporary database");
                    }


                }else{ //It can be anything (including N, Y, X ...)

                    if(c != '\n'){
                        word_size = 0;
                        data_database.sequences[pos_in_database++] = (unsigned char) 'N'; //Convert to N
                        if(pos_in_database == READBUF*n_realloc_database){ 
                            n_realloc_database++; data_database.sequences = (unsigned char *) realloc(data_database.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                        if(data_database.sequences == NULL) terror("Could not reallocate temporary database");
                        }
                    } 
                }
                if(word_size == FIXED_K){
                    //write to hash table
                    
		
                    pointer = ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];


                    if(pointer == NULL){

                        pointer = getNewLocationllpos(mp, &n_pools_used);
                        

                        pointer->pos = pos_in_database;

                        pointer->s_id = data_database.n_seqs-1;

                        pointer->next = NULL;

                    

                    }else{

                        
                        aux = ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                        [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                        [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                        [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                        pointer = getNewLocationllpos(mp, &n_pools_used);

                        pointer->pos = pos_in_database;
                        pointer->s_id = data_database.n_seqs-1;
                        pointer->next = aux;

                        

                    }

                    ct->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                    [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                    [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                    [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]] = pointer;
		
		            memcpy(aux_kmer, &curr_kmer[1], FIXED_K-1);
                    memcpy(curr_kmer, aux_kmer, FIXED_K-1);
                    word_size--;
                }
            }
            word_size = 0;
            
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, database);    
        }
        
    }

    end = clock();

    data_database.total_len = pos_in_database;

    fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64". Hash table building took %e seconds\n", data_database.total_len, (double)(end-begin)/CLOCKS_PER_SEC);
    //close database
    fclose(database);


    
    begin = clock();

    

    //To force reading from the buffer
    idx = READBUF + 1;

    
    data_query.sequences = seq_vector_query;
    data_query.start_pos = query_positions;
    data_query.total_len = 0;
    data_query.n_seqs = 0;
    n_realloc_database = 1;
    
    
    //Print info
    fprintf(stdout, "[INFO] Loading query.\n");   

    
    c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
    while((!feof(query) || (feof(query) && idx < r))){

        if(c == '>'){
            data_query.start_pos[data_query.n_seqs++] = pos_in_query;
            word_size = 0;
            
            if(pos_in_query == READBUF*n_realloc_query){
                n_realloc_query++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_query*sizeof(unsigned char));
                if(data_query.sequences == NULL) terror("Could not reallocate temporary query");
            }

            if(data_query.n_seqs == INITSEQS*n_seqs_query_realloc){
                n_seqs_query_realloc++; data_query.start_pos =  (uint64_t *) realloc(data_query.start_pos, INITSEQS*n_seqs_query_realloc*sizeof(uint64_t));
            }


            while(c != '\n'){ c = buffered_fgetc(temp_seq_buffer, &idx, &r, query); } //Skip ID
                

            while(c != '>' && (!feof(query) || (feof(query) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    data_query.sequences[pos_in_query++] = (unsigned char) c;
                    curr_kmer[word_size++] = (unsigned char) c;
                    
                    if(word_size == FIXED_K){
                        memcpy(aux_kmer, curr_kmer, FIXED_K);
                        aux_kmer[FIXED_K] = '\0';
                        
                        memcpy(aux_kmer, &curr_kmer[1], FIXED_K-1);
                        memcpy(curr_kmer, aux_kmer, FIXED_K-1);
                        word_size--;
                    }
            
                    if(pos_in_query == READBUF*n_realloc_database){ 
                        n_realloc_database++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                        if(data_query.sequences == NULL) terror("Could not reallocate temporary query");
                    }
                }else{
                    if(c != '\n'){
                        word_size = 0;
                        data_query.sequences[pos_in_query++] = (unsigned char) 'N'; //Convert to N
                        if(pos_in_query == READBUF*n_realloc_database){ 
                            n_realloc_database++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_database*sizeof(unsigned char));
                            if(data_query.sequences == NULL) terror("Could not reallocate temporary query");
                        }
                    }
                }
            }
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, query);    
        }
        
    }

    
    
    end = clock();

    
    /*
    for(w0=0;w0<4;w0++){
        for(w1=0;w1<4;w1++){
            for(w2=0;w2<4;w2++){
                for(w3=0;w3<4;w3++){
                    for(w4=0;w4<4;w4++){
                        for(w5=0;w5<4;w5++){
                            for(w6=0;w6<4;w6++){
                                for(w7=0;w7<4;w7++){
                                    for(w8=0;w8<4;w8++){
                                        for(w9=0;w9<4;w9++){
                                            for(w10=0;w10<4;w10++){
                                                for(w11=0;w11<4;w11++){
                                                    fprintf(stdout, "%p\n", ct->table[w0][w1][w2][w3][w4][w5][w6][w7][w8][w9][w10][w11]);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    */
    data_query.total_len = pos_in_query;

    fprintf(stdout, "[INFO] Query loaded and of length %"PRIu64". Took %e seconds\n", data_query.total_len, (double)(end-begin)/CLOCKS_PER_SEC);

    begin = clock();
        
    Head queue_head;
    Queue * first_task = generate_queue(&queue_head, data_query.n_seqs, n_threads, 3);

    /*
    Queue * traverse = queue_head.head;
    while(traverse != NULL){
        printf("current_piece: %"PRIu64"-%"PRIu64"\n", traverse->r1, traverse->r2);
        traverse = traverse->next;
    }
    */


    //reads_per_thread = (uint64_t) (floorl((long double) data_query.n_seqs / (long double) n_threads));
    
    fprintf(stdout, "[INFO] Computing alignments.\n");

    

    /*
    uint64_t z;
    for(z=0; z<POOL_SIZE; z++){
        aux = mp[0].base + z;
        fprintf(stdout, "%p\n", aux);
        fflush(stdout);
    }
    */
    

    for(i=0;i<n_threads;i++){
        hta[i].database = &data_database;
        hta[i].query = &data_query;
        //hta[i].from = i * reads_per_thread;
        //hta[i].to = (i + 1) * reads_per_thread;
        hta[i].container = ct;
        hta[i].accepted_query_reads = 0;
        hta[i].min_e_value = minevalue;
        hta[i].min_coverage = mincoverage;
        hta[i].min_identity = minidentity;
        hta[i].out = out_database;
        hta[i].igap = igap;
        hta[i].egap = egap;
        hta[i].window = window;
        hta[i].mc = mc[i];
        hta[i].table = table[i];
        hta[i].reconstruct_X = reconstruct_X[i];
        hta[i].reconstruct_Y = reconstruct_Y[i];
        hta[i].writing_buffer_alignment = writing_buffer_alignment[i];
        hta[i].my_x = my_x[i];
        hta[i].my_y = my_y[i];
        hta[i].queue_head = &queue_head;
        hta[i].lock = &lock;

        //if(i==n_threads-1) hta[i].to = data_query.n_seqs;

        if( 0 != (error = pthread_create(&threads[i], NULL, computeAlignmentsByThread, (void *) (&hta[i])) )){
            fprintf(stdout, "Thread %"PRIu64" returned %d:", i, error); terror("Could not launch");
        }
    }
    
    //Wait for threads to finish
    for(i=0;i<n_threads;i++){
        pthread_join(threads[i], NULL);
    }

    
    for(i=0;i<n_threads;i++){
        sum_accepted_reads += hta[i].accepted_query_reads;
    }
    
    end = clock();
    fprintf(stdout, "[INFO] Alignments computed in %e seconds.\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fprintf(stdout, "[INFO] %"PRIu64" reads (%"PRIu64") from the query were found in the database (%"PRIu64") at a minimum e-value of %Le and minimum coverage of %d%%.\n", sum_accepted_reads, data_query.n_seqs, data_database.n_seqs, (long double)minevalue, (int) (100*mincoverage));
    fprintf(stdout, "[INFO] The Jaccard-index is: %Le\n", (long double)sum_accepted_reads/((data_database.n_seqs+data_query.n_seqs)-sum_accepted_reads));
    fprintf(stdout, "[INFO] Deallocating heap memory.\n");

    fclose(query);
    
    if(out_database != NULL) fclose(out_database);
    free(temp_seq_buffer);
    
    free(data_database.sequences);
    free(data_database.start_pos);
    free(data_query.sequences);
    free(data_query.start_pos);
    free(ct->table);
    //free(ct);
    free(threads);
    free(hta);


    for(i=0;i<n_threads;i++){

        for(j=0;j<MAX_READ_SIZE;j++){
		    free(table[i][j]);
        }
        free(table[i]);
        free(mc[i]);
        free(reconstruct_X[i]);
        free(reconstruct_Y[i]);
        free(my_x[i]);	
        free(my_y[i]);
        free(writing_buffer_alignment[i]);
    }
    free(table);
    free(mc);
    free(reconstruct_X);
    free(reconstruct_Y);
    free(my_y);
    free(my_x);
    free(writing_buffer_alignment);

    pthread_mutex_destroy(&lock);

    for(i=0;i<=n_pools_used;i++){
        free(mp[i].base);
    }
    //Deallocate queue (its allocated as an array)
    free(first_task);

    //free(mp);
    
    return 0;
}

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t  * n_threads, long double * minevalue, long double * mincoverage, int * igap, int * egap, long double * minidentity, long double * window){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--verbose") == 0) VERBOSE_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           IMSAME -query [query] -db [database]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -n_threads  [Integer:   0<n_threads] (default 4)\n");
            fprintf(stdout, "           -evalue     [Double:    0<=pval<1] (default: 1 * 10^-20)\n");
            fprintf(stdout, "           -coverage   [Double:    0<coverage<=1 (default: 0.5)\n");
            fprintf(stdout, "           -identity   [Double:    0<identity<=1 (default: 0.5)\n");
            fprintf(stdout, "           -igap       [Integer:   (default: 5)\n");
            fprintf(stdout, "           -egap       [Integer:   (default: 2)\n");
            fprintf(stdout, "           -window     [Double:    0<window<=0.5 (default: 0.15)");
            fprintf(stdout, "           -out        [File path]\n");
            fprintf(stdout, "           --verbose   Turns verbose on\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            exit(1);
        }
        if(strcmp(av[pNum], "-query") == 0){
            *query = fopen64(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
        }
        if(strcmp(av[pNum], "-db") == 0){
            *database = fopen64(av[pNum+1], "rt");
            if(database==NULL) terror("Could not open database file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_database = fopen64(av[pNum+1], "wt");
            if(out_database==NULL) terror("Could not open output database file");
        }
        if(strcmp(av[pNum], "-evalue") == 0){
            *minevalue = (long double) atof(av[pNum+1]);
            if(*minevalue < 0) terror("Min-e-value must be larger than zero");
        }
        if(strcmp(av[pNum], "-window") == 0){
            *window = (long double) atof(av[pNum+1]);
            if(*window <= 0 || *window > 0.5) terror("Window percentage size must lie between 0<window<=0.5");
        }
        if(strcmp(av[pNum], "-coverage") == 0){
            *mincoverage = (long double) atof(av[pNum+1]);
            if(*mincoverage <= 0) terror("Min-coverage must be larger than zero");
        }
        if(strcmp(av[pNum], "-identity") == 0){
            *minidentity = (long double) atof(av[pNum+1]);
            if(*minidentity <= 0) terror("Min-identity must be larger than zero");
        }
        if(strcmp(av[pNum], "-igap") == 0){
            *igap = - (atoi(av[pNum+1]));
        }
        if(strcmp(av[pNum], "-egap") == 0){
            *egap = - (atoi(av[pNum+1]));
        }
        if(strcmp(av[pNum], "-n_threads") == 0){
            *n_threads = (uint64_t) atoi(av[pNum+1]);
            if(*n_threads < 0) terror("Number of threads must be larger than zero");
        }
        pNum++;
    }
    
    if(*query==NULL || *database==NULL) terror("A query and database is required");
}

