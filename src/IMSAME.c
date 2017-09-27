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


uint64_t custom_kmer = 12; // Defined as external in structs.h

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t  * n_threads, long double * minevalue, long double * mincoverage, int * igap, int * egap, long double * minidentity, long double * window, unsigned char * full_comp, uint64_t * custom_kmer, unsigned char * hits_only, uint64_t * n_parts);

int VERBOSE_ACTIVE = 0;

int main(int argc, char ** av){


    


    clock_t begin, end;
    
    int error; //To tell if threads could not be launched
    uint64_t i,j;

    //query to read kmers from, database to find seeds
    FILE * query = NULL, * database = NULL, * out_database = NULL;
    long double minevalue = 1/powl(10, 10); //Default 1 * 10^-10
    
    long double mincoverage = 0.5, minidentity = 0.5, window = 0.15; //Default
    int igap = -5, egap = -2;
    unsigned char full_comp = FALSE;
    unsigned char hits_only = FALSE;

    uint64_t n_threads = 4;
    uint64_t n_parts = 3; // Default is 3

    init_args(argc, av, &query, &database, &out_database, &n_threads, &minevalue, &mincoverage, &igap, &egap, &minidentity, &window, &full_comp, &custom_kmer, &hits_only, &n_parts);
    
    //uint64_t reads_per_thread;
    uint64_t sum_accepted_reads = 0;

    begin = clock();

    fprintf(stdout, "[INFO] Init. quick table\n");

    pthread_t * threads = (pthread_t *) malloc(n_threads * sizeof(pthread_t));
    if(threads == NULL) terror("Could not create threads");

    pthread_t * loading_threads = (pthread_t *) malloc(FIXED_LOADING_THREADS * sizeof(pthread_t));
    if(loading_threads == NULL) terror("Could not create loading threads");


    HashTableArgs * hta = (HashTableArgs *) malloc(n_threads*sizeof(HashTableArgs));
    if(hta == NULL) terror("Could not allocate arguments for hash table");
    
    pthread_mutex_t lock; //The mutex to lock the queue
    if (pthread_mutex_init(&lock, NULL) != 0) terror("Could not init mutex");

    // To be used if only computing hits
    uint64_t ** hits_table = NULL;

    unsigned char ** my_x = (unsigned char **) malloc(n_threads * sizeof(unsigned char*));
    unsigned char ** my_y = (unsigned char **) malloc(n_threads * sizeof(unsigned char*));

    unsigned char ** marker_taggs = NULL; // To be used with full comparison

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
            table[i][j] = (struct cell *) malloc((1+MAX_WINDOW_SIZE)*sizeof(struct cell));
		    if(table[i][j] == NULL) terror("Could not allocate memory for second loop of table");
            // Delete this 
            /*
            uint64_t r;
            for(r=0;r<MAX_WINDOW_SIZE+1;r++){
                table[i][j][r].score = INT64_MIN;
                table[i][j][r].xfrom = 10000000000;
                table[i][j][r].yfrom = 10000000000;
            }
            */
	    }
    	mc[i] = (struct positioned_cell *) malloc(MAX_READ_SIZE * sizeof(struct positioned_cell));
    	my_x[i] = (unsigned char *) malloc(2*MAX_READ_SIZE * sizeof(unsigned char));
    	my_y[i] = (unsigned char *) malloc(2*MAX_READ_SIZE * sizeof(unsigned char));
    	reconstruct_X[i] = (char *) malloc(2*MAX_READ_SIZE * sizeof(char));
    	reconstruct_Y[i] = (char *) malloc(2*MAX_READ_SIZE * sizeof(char));
	    writing_buffer_alignment[i] = (char *) malloc(6*MAX_READ_SIZE*sizeof(char)); //6 times because of 2 times the length of the max of the read, and that happens 3 times (seqX,align,Y)
	    if(table[i] == NULL || mc[i] == NULL || my_x[i] == NULL || my_y[i] == NULL || reconstruct_X[i] == NULL || reconstruct_Y[i] == NULL || writing_buffer_alignment[i] == NULL) terror("Could not allocate buffer for alignment output");
    }



    end = clock();
    fprintf(stdout, "[INFO] Initialization took %e seconds.\n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Variables to account for positions
    //Print info
    fprintf(stdout, "[INFO] Loading database.\n");
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
    unsigned char ** seq_vector_database = (unsigned char **) malloc(FIXED_LOADING_THREADS*sizeof(unsigned char *));
    if(seq_vector_database == NULL) terror("Could not allocate memory for database vector");
    uint64_t ** database_positions = (uint64_t **) malloc(FIXED_LOADING_THREADS*sizeof(uint64_t));
    if(database_positions == NULL) terror("Could not allocate database sequences positions");

    //Mempool_l * mp = (Mempool_l *) malloc(MAX_MEM_POOLS*sizeof(Mempool_l));
    //if(mp == NULL) terror("Could not allocate vectors for memory pools");
    //Mempool_l mp[FIXED_LOADING_THREADS][MAX_MEM_POOLS];

    
    Mempool_l ** mp = (Mempool_l **) malloc(FIXED_LOADING_THREADS*sizeof(Mempool_l *));
    if(mp == NULL) terror("Could not allocate memory pools");
    for(i=0; i<FIXED_LOADING_THREADS; i++){
        mp[i] = (Mempool_l *) malloc(MAX_MEM_POOLS*sizeof(Mempool_l));
        if(mp[i] == NULL) terror("Could not allocate individual memory pools");
    }
    

    Container * ct_A = (Container *) calloc(1, sizeof(Container));
    if(ct_A == NULL) terror("Could not allocate container A");
    Container * ct_B = (Container *) calloc(1, sizeof(Container));
    if(ct_B == NULL)    terror("Could not allocate container B");
    Container * ct_C = (Container *) calloc(1, sizeof(Container));
    if(ct_C == NULL) terror("Could not allocate container C");
    Container * ct_D = (Container *) calloc(1, sizeof(Container));
    if(ct_D == NULL) terror("Could not allocate container D");

    SeqInfo data_database[FIXED_LOADING_THREADS];
    uint64_t full_db_n_seqs = 0;
    
    
    unsigned char curr_kmer[custom_kmer];
    unsigned char aux_kmer[custom_kmer+1];
    curr_kmer[0] = '\0';
    uint64_t word_size = 0;

    
    SeqInfo data_query;
    
    //To force reading from the buffer
    idx = READBUF + 1;
    
    //Vector to store query seq
    unsigned char * seq_vector_query = (unsigned char *) malloc(READBUF*sizeof(unsigned char));
    if(seq_vector_query == NULL) terror("Could not allocate memory for query vector");
    uint64_t n_realloc_query = 1, pos_in_query = 0, n_seqs_query_realloc = 1;
    uint64_t * query_positions = (uint64_t *) malloc(INITSEQS*sizeof(uint64_t));
    if(query_positions == NULL) terror("Could not allocate query sequences positions");
    
    
    
    // Read number of sequences and load into RAM
    begin = clock();
    fseek(database, 0L, SEEK_END);
    uint64_t db_temp_size = ftell(database);
    char * load_buffer = (char *) malloc(db_temp_size * sizeof(char));
    if(load_buffer == NULL) terror("Could not allocate intermediate buffer for threads sequence array");
    fseek(database, 0L, SEEK_SET);

    if(db_temp_size != fread(load_buffer, sizeof(char), db_temp_size, database)) terror("Could not read full sequence");

    LoadingDBArgs args_DB_load[FIXED_LOADING_THREADS];


    args_DB_load[0].data_database = &data_database[0];
    args_DB_load[1].data_database = &data_database[1];
    args_DB_load[2].data_database = &data_database[2];
    args_DB_load[3].data_database = &data_database[3];

    args_DB_load[0].ct = ct_A;
    args_DB_load[1].ct = ct_B;
    args_DB_load[2].ct = ct_C;
    args_DB_load[3].ct = ct_D;
    for(i=0; i<FIXED_LOADING_THREADS; i++){
        args_DB_load[i].read_to = 0;
        args_DB_load[i].read_from = 0;
    }

    //uint64_t a_fourth = db_temp_size / 4;
    

    //get_num_seqs_and_length(load_buffer, &full_db_n_seqs, &db_temp_size, args_DB_load);

    end = clock();
    fprintf(stdout, "[INFO] Loading into RAM took %e seconds.\n", (double)(end-begin)/CLOCKS_PER_SEC);
    begin = clock();
    
    /*
    char * temp_seq_buffer;
    SeqInfo * data_database;
    uint64_t t_len;
    uint64_t word_size;
    uint64_t read_from;
    uint64_t read_to;
    char thread_id;
    Mempool_l * mp;
    uint64_t n_pools_used;
    */

    // Launch threads to process database

    args_DB_load[0].thread_id = 'A';
    args_DB_load[1].thread_id = 'B';
    args_DB_load[2].thread_id = 'C';
    args_DB_load[3].thread_id = 'D';
    

    for(i=0; i<FIXED_LOADING_THREADS; i++){

        //seq_vector_database[i] = (unsigned char *) malloc((args_DB_load[i].read_to - args_DB_load[i].read_from)*sizeof(unsigned char));
        seq_vector_database[i] = (unsigned char *) malloc((READBUF)*sizeof(unsigned char));
        //database_positions[i] = (uint64_t *) malloc((1+data_database[i].n_seqs)*sizeof(uint64_t));
        database_positions[i] = (uint64_t *) malloc((INITSEQS)*sizeof(uint64_t));
        if(seq_vector_database[i] == NULL || database_positions[i] == NULL) terror("Could not allocate memory for individual database vectors");
        data_database[i].sequences = seq_vector_database[i];
        
        
        //To hold all information related to database
        args_DB_load[i].n_pools_used = 0;
        init_mem_pool_llpos(&mp[i][args_DB_load[i].n_pools_used]);

        data_database[i].start_pos = database_positions[i];
        
        args_DB_load[i].n_allocs = 1;
        args_DB_load[i].read_from = i * (db_temp_size / 4);
        args_DB_load[i].read_to = (i+1) * (db_temp_size / 4);
        args_DB_load[i].temp_seq_buffer = load_buffer;
        args_DB_load[i].t_len = db_temp_size;
        args_DB_load[i].word_size = custom_kmer;
        args_DB_load[i].mp = mp[i];
        args_DB_load[i].offloaded = 0;
        
        if( 0 != (error = pthread_create(&loading_threads[i], NULL, load_input, (void *) (&args_DB_load[i])) )){
            fprintf(stdout, "[@loading] Thread %"PRIu64" returned %d:", i, error); terror("Could not launch");
        }
    }
    
    //Wait for threads to finish
    for(i=0;i<FIXED_LOADING_THREADS;i++){
        pthread_join(loading_threads[i], NULL);
    }

    // Deallocate memory not needed anymore
    free(load_buffer);
    free(loading_threads);

    
    
    //fprintf(stdout, "[INFO] WARNING!!!!!!!!! USING NON OVERLAPPING MERS, WHICH IS NOT INCLUDED AS OPTION!!!! DISABLE\n");

    

    

    //data_database.total_len = pos_in_database;
    //printf("tables have %"PRIu64" %"PRIu64", %"PRIu64" %"PRIu64"\n", data_database[0].total_len, data_database[1].total_len, data_database[2].total_len, data_database[3].total_len);
    uint64_t full_db_len = data_database[0].total_len + data_database[1].total_len + data_database[2].total_len + data_database[3].total_len;
    full_db_n_seqs = args_DB_load[0].contained_reads + args_DB_load[1].contained_reads + args_DB_load[2].contained_reads + args_DB_load[3].contained_reads;

    end = clock();
    fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64" (%"PRIu64" sequences). Hash table building took %e seconds\n", full_db_len, full_db_n_seqs, (double)(end-begin)/CLOCKS_PER_SEC);
    //close database
    fclose(database);

    // If FULL comparison is performed, allocate memory for the marker taggs
    if(full_comp == TRUE){
        marker_taggs = (unsigned char **) malloc(n_threads * sizeof(unsigned char *));
        if(marker_taggs == NULL) terror("Could not allocate marker taggs");
        for(i=0;i<n_threads;i++){
            marker_taggs[i] = (unsigned char *) malloc(full_db_n_seqs);
            if(marker_taggs[i] == NULL) terror("Could not allocate second loop of marker taggs");
        }
    }
    if(hits_only == TRUE){
        hits_table = (uint64_t **) calloc(n_threads, sizeof(uint64_t *));
        if(hits_table == NULL) terror("Could not allocate hits table");
        uint64_t z;
        for(z=0;z<n_threads;z++){
            hits_table[z] = (uint64_t *) calloc(full_db_n_seqs, sizeof(uint64_t));
            if(hits_table[z] == NULL) terror("Could not allocate sub table of hits");
        }
    }

    
    begin = clock();

    

    //To force reading from the buffer
    idx = READBUF + 1;

    
    data_query.sequences = seq_vector_query;
    data_query.start_pos = query_positions;
    data_query.total_len = 0;
    data_query.n_seqs = 0;
    
    
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
                    
                    if(word_size == custom_kmer){
                        memcpy(aux_kmer, curr_kmer, custom_kmer);
                        aux_kmer[custom_kmer] = '\0';
                        
                        memcpy(aux_kmer, &curr_kmer[1], custom_kmer-1);
                        memcpy(curr_kmer, aux_kmer, custom_kmer-1);
                        word_size--;
                    }
            
                    if(pos_in_query == READBUF*n_realloc_query){ 
                        n_realloc_query++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_query*sizeof(unsigned char));
                        if(data_query.sequences == NULL) terror("Could not reallocate temporary query");
                    }
                }else{
                    if(c != '\n' && c != '\r' && c != '>'){
                        word_size = 0;
                        data_query.sequences[pos_in_query++] = (unsigned char) 'N'; //Convert to N
                        if(pos_in_query == READBUF*n_realloc_query){ 
                            n_realloc_query++; data_query.sequences = (unsigned char *) realloc(data_query.sequences, READBUF*n_realloc_query*sizeof(unsigned char));
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

    
    data_query.total_len = pos_in_query;

    fprintf(stdout, "[INFO] Query loaded and of length %"PRIu64". Took %e seconds\n", data_query.total_len, (double)(end-begin)/CLOCKS_PER_SEC);

    begin = clock();
        
    Head queue_head;
    Queue * first_task = generate_queue(&queue_head, data_query.n_seqs, n_threads, n_parts);

    /*
    Queue * traverse = queue_head.head;
    while(traverse != NULL){
        printf("current_piece: %"PRIu64"-%"PRIu64"\n", traverse->r1, traverse->r2);
        traverse = traverse->next;
    }
    getchar();
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

    // Make the full db
    uint64_t contained_reads[FIXED_LOADING_THREADS] = {0,0,0,0};
    uint64_t base_coordinates[FIXED_LOADING_THREADS] = {0,0,0,0};
    //contained_reads[0] = 0;
    //base_coordinates[0] = 0;
    for(i=1;i<FIXED_LOADING_THREADS;i++){
        contained_reads[i] = args_DB_load[i-1].contained_reads + contained_reads[i-1];
        base_coordinates[i] = args_DB_load[i-1].base_coordinates + base_coordinates[i-1];
    }

    /*
    for(i = 0; i < 4 ; i++){
        printf("total len: %"PRIu64"\n", data_database[i].total_len);
        
    }
    getchar();
    */
    /*
    for(i = 0; i < 4 ; i++){
        
        printf("c:%"PRIu64" - b:%"PRIu64"\n", contained_reads[i], base_coordinates[i]);
    }
    getchar();
    */
    



    SeqInfo final_db;
    final_db.sequences = (unsigned char *) malloc(full_db_len * sizeof(unsigned char));
    if(final_db.sequences == NULL) terror ("Could not allocate final database sequences");
    memcpy(&final_db.sequences[0], &data_database[0].sequences[0], data_database[0].total_len);
    memcpy(&final_db.sequences[data_database[0].total_len], &data_database[1].sequences[0], data_database[1].total_len);
    memcpy(&final_db.sequences[data_database[0].total_len+data_database[1].total_len], &data_database[2].sequences[0], data_database[2].total_len);
    memcpy(&final_db.sequences[data_database[0].total_len+data_database[1].total_len+data_database[2].total_len], &data_database[3].sequences[0], data_database[3].total_len);
    final_db.start_pos = (uint64_t *) malloc(full_db_n_seqs * sizeof(uint64_t));
    if(final_db.start_pos == NULL) terror("Could not allocate final db starting positions");
    // Copy them with offset
    i=0;
    while(i<data_database[0].n_seqs){ final_db.start_pos[i] = data_database[0].start_pos[i]; ++i; }
    j=0;
    //printf("switch\n");
    while(j<data_database[1].n_seqs){ final_db.start_pos[i] = data_database[1].start_pos[j] + data_database[0].total_len; ++i; ++j; }
    j=0;
    //printf("switch\n");
    while(j<data_database[2].n_seqs){ final_db.start_pos[i] = data_database[2].start_pos[j] + data_database[1].total_len + data_database[0].total_len; ++i; ++j; }
    j=0;
    //printf("switch\n");
    while(j<data_database[3].n_seqs){ final_db.start_pos[i] = data_database[3].start_pos[j] + data_database[2].total_len + data_database[1].total_len + data_database[0].total_len; ++i; ++j; }
    
    final_db.total_len = full_db_len;
    final_db.n_seqs = full_db_n_seqs;
    
    for(i=0;i<FIXED_LOADING_THREADS;i++){
        free(data_database[i].sequences);
        free(data_database[i].start_pos);
    }
    

    // Debug
    /*
    for(i=0; i<full_db_n_seqs-1; i++){
        printf("%"PRIu64" - %"PRIu64"\n", final_db.start_pos[i], final_db.start_pos[i+1]);
        getchar();
    }
    */
    
    uint64_t idx_tablespaces = 0; // Indexes for loading the data tables
    

    Container * ptr_table_redirect[4];
    ptr_table_redirect[0] = hta->container_a;
    ptr_table_redirect[1] = hta->container_b;
    ptr_table_redirect[2] = hta->container_c;
    ptr_table_redirect[3] = hta->container_d;

    while(idx_tablespaces <= args_DB_load[0].offloaded || idx_tablespaces <= args_DB_load[1].offloaded || idx_tablespaces <= args_DB_load[2].offloaded || idx_tablespaces <= args_DB_load[3].offloaded){
        for(i=0;i<n_threads;i++){
            hta[i].id = i;
            hta[i].database = &final_db;
            hta[i].query = &data_query;
            
            
            // The appropriate data table and container must be loaded
            char read_name[MAXPATH]; read_name[0] = '\0';
            for(j=0; j<FIXED_LOADING_THREADS; j++){
                sprintf(&read_name[0], "tablespace-%c-%"PRIu64, args_DB_load[j].thread_id, idx_tablespaces);
                FILE * tablespace_handler = fopen(read_name, "rb");
                if(tablespace_handler == NULL){ fprintf(stdout, "At %s\n", read_name); terror("Could not open tablespace"); } 
                

                fclose(tablespace_handler);

            }
            

            hta[i].container_a = ct_A;
            hta[i].container_b = ct_B;
            hta[i].container_c = ct_C;
            hta[i].container_d = ct_D;
            hta[i].contained_reads = contained_reads;
            hta[i].base_coordinates = base_coordinates;
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
            hta[i].full_comp = full_comp;
            hta[i].mp_pools = mp;
            if(full_comp){
                hta[i].markers = marker_taggs[i];
            }
            if(hits_only) hta[i].hits = hits_table[i]; else hta[i].hits = NULL;
            

            //if(i==n_threads-1) hta[i].to = data_query.n_seqs;

            if( 0 != (error = pthread_create(&threads[i], NULL, computeAlignmentsByThread, (void *) (&hta[i])) )){
                fprintf(stdout, "Thread %"PRIu64" returned %d:", i, error); terror("Could not launch");
            }


            ++idx_tablespaces;
        }
        
        //Wait for threads to finish
        for(i=0;i<n_threads;i++){
            pthread_join(threads[i], NULL);
        }

        
        for(i=0;i<n_threads;i++){
            sum_accepted_reads += hta[i].accepted_query_reads;
        }
    }
    
    // Accumulate hits just in case
    if(hits_only == TRUE){
        for(i=0;i<full_db_n_seqs;i++){
            for(j=1;j<n_threads;j++){
                hta[0].hits[i] += hta[j].hits[i];
            }
            if(out_database != NULL){
                
                if(i == final_db.n_seqs - 1){
                    fprintf(out_database, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", i, hta[0].hits[i], full_db_len - final_db.start_pos[i]);
                }else{
                    fprintf(out_database, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", i, hta[0].hits[i], final_db.start_pos[i+1] - final_db.start_pos[i]);
                }
            }
        }
    }
    
    
    end = clock();
    fprintf(stdout, "[INFO] Alignments computed in %e seconds.\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fprintf(stdout, "[INFO] %"PRIu64" reads (%"PRIu64") from the query were found in the database (%"PRIu64") at a minimum e-value of %Le and minimum coverage of %d%%.\n", sum_accepted_reads, data_query.n_seqs, final_db.n_seqs, (long double)minevalue, (int) (100*mincoverage));
    fprintf(stdout, "[INFO] The Jaccard-index is: %Le\n", (long double)sum_accepted_reads/((final_db.n_seqs+data_query.n_seqs)-sum_accepted_reads));
    fprintf(stdout, "[INFO] Deallocating heap memory.\n");

    fclose(query);
    
    if(out_database != NULL) fclose(out_database);
    free(temp_seq_buffer);
    

    if(hits_only == TRUE){
        // Write hits here
        for(i=0;i<n_threads;i++){
            free(hits_table[i]);
        }
        free(hits_table);
    }

    free(final_db.sequences);
    free(final_db.start_pos);
    free(data_query.sequences);
    free(data_query.start_pos);
    free(ct_A->table);
    free(ct_B->table);
    free(ct_C->table);
    free(ct_D->table);
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

        if(full_comp == TRUE){
            free(marker_taggs[i]);
        }
        
    }
    if(full_comp == TRUE) free(marker_taggs);

    
    free(table);
    free(mc);
    free(reconstruct_X);
    free(reconstruct_Y);
    free(my_y);
    free(my_x);
    free(writing_buffer_alignment);

    pthread_mutex_destroy(&lock);

    

    for(i=0; i<FIXED_LOADING_THREADS; i++){
        for(j=0;j<=args_DB_load[i].n_pools_used;j++){
            free(mp[i][j].base);
        }
        free(mp[i]);
    }
    free(mp);
    free(seq_vector_database);
    free(database_positions);
    
    //Deallocate queue (its allocated as an array)
    free(first_task);

    //free(mp);
    
    return 0;
}

void init_args(int argc, char ** av, FILE ** query, FILE ** database, FILE ** out_database, uint64_t  * n_threads, long double * minevalue, long double * mincoverage, int * igap, int * egap, long double * minidentity, long double * window, unsigned char * full_comp, uint64_t * custom_kmer, unsigned char * hits_only, uint64_t * n_parts){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--verbose") == 0) VERBOSE_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           IMSAME -query [query] -db [database]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -n_threads  [Integer:   0<n_threads] (default 4)\n");
            fprintf(stdout, "           -evalue     [Double:    0<=pval<1] (default: 1 * 10^-10)\n");
            fprintf(stdout, "           -coverage   [Double:    0<coverage<=1 (default: 0.5)\n");
            fprintf(stdout, "           -identity   [Double:    0<identity<=1 (default: 0.5)\n");
            fprintf(stdout, "           -igap       [Integer:   (default: 5)\n");
            fprintf(stdout, "           -egap       [Integer:   (default: 2)\n");
            fprintf(stdout, "           -window     [Double:    0<window<=0.5 (default: 0.15)\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 12)]\n");
	        fprintf(stdout, "		-n_parts    [Integer:	n>0 (default 3)]\n");
            fprintf(stdout, "           -out        [File path]\n");
            fprintf(stdout, "           --full      Does not stop at first match and reports all equalities\n");
            fprintf(stdout, "           --verbose   Turns verbose on\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "           --hits      Compute only the non-overlapping hits\n");
            exit(1);
        }
        if(strcmp(av[pNum], "--full") == 0) *full_comp = TRUE;
        if(strcmp(av[pNum], "--hits") == 0) *hits_only = TRUE;
	if(strcmp(av[pNum], "-n_parts") == 0) *n_parts = (uint64_t) atoi(av[pNum+1]);
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
        if(strcmp(av[pNum], "-kmer") == 0){
            *custom_kmer = (uint64_t) atoi(av[pNum+1]);
            if(*custom_kmer < 2) terror("K-mer size must be larger than 1");
        }
        pNum++;
    }
    
    if(*query==NULL || *database==NULL) terror("A query and database is required");
}

