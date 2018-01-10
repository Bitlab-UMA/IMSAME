#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))

int64_t compare_letters(unsigned char a, unsigned char b){
    if(a != (unsigned char) 'N' && a != (unsigned char) '>') return (a == b) ? POINT : -POINT;
    return -POINT;
}

llpos * getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used){

    if(mp[*n_pools_used].current == POOL_SIZE){
        *n_pools_used += 1;
        if(*n_pools_used == MAX_MEM_POOLS) terror("Reached max pools");
        init_mem_pool_llpos(&mp[*n_pools_used]);
        
    }

    llpos * new_pos = mp[*n_pools_used].base + mp[*n_pools_used].current;
    mp[*n_pools_used].current++;

    
    return new_pos;
}

void init_mem_pool_llpos(Mempool_l * mp){
    mp->base = (llpos *) calloc(POOL_SIZE, sizeof(llpos));
    if(mp->base == NULL) terror("Could not request memory pool");
    mp->current = 0;
}


void * load_input(void * a){

    LoadingDBArgs * ldbargs = (LoadingDBArgs *) a;
    
    // Requires
    /*
    char * temp_seq_buffer;
    SeqInfo * data_database; 1 per thread
    uint64_t t_len;
    uint64_t word_size;
    uint64_t read_from;
    uint64_t read_to;
    char thread_id;
    */

    uint64_t c_pos;
    
    unsigned char curr_kmer[custom_kmer];
    unsigned char aux_kmer[custom_kmer+1];
    curr_kmer[0] = '\0';
    uint64_t word_size = 0, pos_in_database = 0;
    unsigned char char_converter[91];
    uint64_t curr_seq = 0;
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;
    llpos * aux, * pointer;

    char c;

    /*
    if(ldbargs->thread_id == 'A'){
        printf("read to is: %"PRIu64"\n", ldbargs->read_to);
        printf("make sure: %c %c %c\n", ldbargs->temp_seq_buffer[ldbargs->read_to-1], ldbargs->temp_seq_buffer[ldbargs->read_to], ldbargs->temp_seq_buffer[ldbargs->read_to+1]);
        
        uint64_t z = ldbargs->read_to-1;
        while(ldbargs->temp_seq_buffer[z] != '>'){
            printf("%c", ldbargs->temp_seq_buffer[z]);
            z--;
        }
        getchar();
    }
    if(ldbargs->thread_id == 'C'){
        printf("HELLOOOOOOO im going from %"PRIu64"\n", ldbargs->read_from);
    }
    */

    c_pos = ldbargs->read_from;
    while(ldbargs->temp_seq_buffer[c_pos] != '>') ++c_pos;
    ldbargs->read_from = c_pos;
    c_pos = ldbargs->read_to;
    while(c_pos < ldbargs->t_len && ldbargs->temp_seq_buffer[c_pos] != '>') ++c_pos;
    ldbargs->read_to = c_pos;

    c_pos = ldbargs->read_from;
    c = ldbargs->temp_seq_buffer[c_pos];

    
    //printf("thread going from %"PRIu64" to %"PRIu64"\n", ldbargs->read_from, ldbargs->read_to);
    

    while(c_pos < ldbargs->read_to){
        

        if(c == '>'){
            
            //if(ldbargs->thread_id == 'G') printf("putting in %"PRIu64" @ %"PRIu64"\n", curr_seq, c_pos);
            ldbargs->data_database->start_pos[curr_seq] = pos_in_database; ++curr_seq;

            // REalloc sequences and sequence index
            if(pos_in_database == READBUF*ldbargs->n_allocs){
                ldbargs->n_allocs++; ldbargs->data_database->sequences = (unsigned char *) realloc(ldbargs->data_database->sequences, READBUF*ldbargs->n_allocs*sizeof(unsigned char));
                if(ldbargs->data_database->sequences == NULL) terror("Could not reallocate temporary database");
            }

            if(curr_seq == INITSEQS*ldbargs->n_allocs){
                ldbargs->n_allocs++; ldbargs->data_database->start_pos =  (uint64_t *) realloc(ldbargs->data_database->start_pos, INITSEQS*ldbargs->n_allocs*sizeof(uint64_t));
            }



            while(c != '\n'){ c = ldbargs->temp_seq_buffer[c_pos]; ++c_pos; }  //Skip ID

            while(c != '>' && c_pos < ldbargs->read_to){ //Until next id

                //if(ldbargs->thread_id == 'A') printf("!!!!!!%"PRIu64" from:%"PRIu64", to %"PRIu64"\n", c_pos, ldbargs->read_from, ldbargs->read_to);
                c = ldbargs->temp_seq_buffer[c_pos]; ++c_pos;
                c = toupper(c);
                if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
                    curr_kmer[word_size] = (unsigned char) c;
                    if(word_size < custom_kmer) ++word_size;

                    ldbargs->data_database->sequences[pos_in_database] = (unsigned char) c; ++pos_in_database;
                    // REalloc sequences and sequence index
                    if(pos_in_database == READBUF*ldbargs->n_allocs){
                        ldbargs->n_allocs++; ldbargs->data_database->sequences = (unsigned char *) realloc(ldbargs->data_database->sequences, READBUF*ldbargs->n_allocs*sizeof(unsigned char));
                        if(ldbargs->data_database->sequences == NULL) terror("Could not reallocate temporary database");
                    }

                    if(curr_seq == INITSEQS*ldbargs->n_allocs){
                        ldbargs->n_allocs++; ldbargs->data_database->start_pos =  (uint64_t *) realloc(ldbargs->data_database->start_pos, INITSEQS*ldbargs->n_allocs*sizeof(uint64_t));
                    }
                }else{ //It can be anything (including N, Y, X ...)

                    if(c != '\n' && c != '\r' && c != '>'){
                        word_size = 0;
                        ldbargs->data_database->sequences[pos_in_database] = (unsigned char) 'N'; ++pos_in_database; //Convert to N
                        // REalloc sequences and sequence index
                        if(pos_in_database == READBUF*ldbargs->n_allocs){
                            ldbargs->n_allocs++; ldbargs->data_database->sequences = (unsigned char *) realloc(ldbargs->data_database->sequences, READBUF*ldbargs->n_allocs*sizeof(unsigned char));
                            if(ldbargs->data_database->sequences == NULL) terror("Could not reallocate temporary database");
                        }

                        if(curr_seq == INITSEQS*ldbargs->n_allocs){
                            ldbargs->n_allocs++; ldbargs->data_database->start_pos =  (uint64_t *) realloc(ldbargs->data_database->start_pos, INITSEQS*ldbargs->n_allocs*sizeof(uint64_t));
                        }
                    } 
                }
                if(word_size == custom_kmer){
                    //write to hash table
                    
		
                    pointer = ldbargs->ct->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                    

                    if(pointer == NULL){

                        pointer = getNewLocationllpos(ldbargs->mp, &ldbargs->n_pools_used);
                        pointer->pos = pos_in_database;
                        pointer->extended_hash = hashOfWord(&curr_kmer[FIXED_K], custom_kmer - FIXED_K);
                        pointer->s_id = curr_seq-1;
                        pointer->next = NULL;

                    

                    }else{

                        
                        aux = ldbargs->ct->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                        pointer = getNewLocationllpos(ldbargs->mp, &ldbargs->n_pools_used);

                        pointer->pos = pos_in_database;
                        pointer->extended_hash = hashOfWord(&curr_kmer[FIXED_K], custom_kmer - FIXED_K);
                        pointer->s_id = curr_seq-1;
                        pointer->next = aux;

                    }

                    ldbargs->ct->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]] = pointer;
    

                    // CURRENTLY USING OVERLAPPING
                    
                    memcpy(aux_kmer, &curr_kmer[1], custom_kmer-1);
                    memcpy(curr_kmer, aux_kmer, custom_kmer-1);
                    word_size--;
                        
                    // For NON OVERLAPPING ENABLE THIS
                    //word_size = 0;
                }

            }
            word_size = 0;
        }else{
            c = ldbargs->temp_seq_buffer[c_pos]; ++c_pos;
        }

    }
    /*
    if(ldbargs->thread_id == 'T'){
        uint64_t j;
        for(j=0; j < curr_seq-1; j++){
            printf("%"PRIu64" - %"PRIu64"\n", ldbargs->data_database->start_pos[j], ldbargs->data_database->start_pos[j+1]);
        }
    }
    */
    /*
    if(ldbargs->thread_id == 'A'){
        printf("\nAT %"PRIu64", and pos_database = %"PRIu64"\n", ldbargs->data_database->start_pos[curr_seq-1], pos_in_database);        
        uint64_t z = pos_in_database;
        while(z > ldbargs->data_database->start_pos[curr_seq-1]){
            printf("%c", ldbargs->data_database->sequences[z]);
            z--;
        }

        getchar();
    }
    */

    ldbargs->data_database->start_pos[curr_seq] = pos_in_database;
    ldbargs->data_database->total_len = pos_in_database;
    ldbargs->contained_reads = curr_seq;
    ldbargs->data_database->n_seqs = curr_seq;
    ldbargs->base_coordinates = pos_in_database;
    return NULL;
}


void * computeAlignmentsByThread(void * a){

/*
typedef struct {
    SeqInfo * database; //Database sequence and lengths
    SeqInfo * query;    //Query sequence and lengths
    uint64_t from;      //Starting READ to compute alignments from
    uint64_t to;        //End READ to compute alignments from
    Container * container; //Container to hold the multidimensional array
    uint64_t accepted_query_reads; //Number of reads that have a fragment with evalue less than specified
    long double min_e_value;    //Minimum evalue to accept read
} HashTableArgs;
*/

    HashTableArgs * hta = (HashTableArgs *) a;
    Queue * my_current_task = NULL;

    unsigned char char_converter[91];
    char_converter[(unsigned char)'A'] = 0;
    char_converter[(unsigned char)'C'] = 1;
    char_converter[(unsigned char)'G'] = 2;
    char_converter[(unsigned char)'T'] = 3;
    Quickfrag qf;
    int64_t * cell_path_y = (int64_t *) malloc(MAX_READ_SIZE*sizeof(int64_t));
    if(cell_path_y == NULL) terror("Could not allocate cell paths");
    

    Point p0, p1, p2, p3; //Points for NW anchored
    p0.x = 0; p0.y = 0;

    Container * ptr_table_redirect[4];
    ptr_table_redirect[0] = hta->container_a;
    ptr_table_redirect[1] = hta->container_b;
    ptr_table_redirect[2] = hta->container_c;
    ptr_table_redirect[3] = hta->container_d;
    unsigned char current_table = 0;

    

    //To keep track of which reads are we reading
    uint64_t curr_read, curr_db_seq, xlen, ylen;
    uint64_t crrSeqL, pos_of_hit = 0xFFFFFFFFFFFFFFFF;

    //Reading from buffer
    char c;
    unsigned char curr_kmer[custom_kmer], b_aux[custom_kmer];
    llpos * aux;

    // For NW-alignment
    int NWaligned;
    uint64_t n_hits, alignments_tried;

    BasicAlignment ba; //The resulting alignment from the NW
    uint64_t curr_pos = 0; //Reading-head position
    uint64_t up_to = 0;

    int64_t last_diagonal = INT64_MIN; // Diagonal to skip repeated hits
    unsigned char already_aligned = FALSE; // To not count more times the same read
    

    

    //Get next operation in queue
    while(NULL != ( my_current_task = get_task_from_queue(hta->queue_head, hta->lock))){
        //Initialize all variables
        
        qf.x_start = qf.y_start = qf.t_len = 0;
        qf.e_value = LDBL_MAX;
        last_diagonal = INT64_MIN;

        //Starting from
        curr_read = my_current_task->r1;
        crrSeqL = 0; pos_of_hit = 0;

        curr_kmer[0] = '\0'; b_aux[0] = '\0';
        aux = NULL;

        NWaligned = 0;
        n_hits = 0;
        alignments_tried = 0;

        ba.identities = 0; ba.length = 0; ba.igaps = 0xFFFFFFFFFFFFFFFF; ba.egaps = 0xFFFFFFFFFFFFFFFF;
        memset(&hta->markers[0], 0, hta->database->n_seqs); // Reset used tags
        already_aligned = FALSE;

        //Set current header position at the position of the read start (the ">")
        curr_pos = hta->query->start_pos[curr_read]; //Skip the ">"
        c = (char) hta->query->sequences[curr_pos];

        //printf("Im doing from %"PRIu64" to %"PRIu64", nseqs=%"PRIu64"\n", my_current_task->r1, my_current_task->r2, hta->query->n_seqs);
        //getchar();

        while(curr_read < my_current_task->r2 && curr_pos < hta->query->total_len){

            
            
            if(curr_read != hta->query->n_seqs) up_to = hta->query->start_pos[curr_read+1]-1; else up_to = hta->query->total_len;
            //printf("Currrpos: %"PRIu64" up to: %"PRIu64" on read: %"PRIu64"\n", curr_pos, up_to, curr_read);

            if (curr_pos == up_to) { // Comment, empty or quality (+) line
                crrSeqL = 0; // Reset buffered sequence length
                #ifdef VERBOSE
                printf("Read: %"PRIu64" yielded (%d)\n", curr_read, NWaligned);
                #endif
                //if(NWaligned == 0){ printf("Read: %"PRIu64" yielded (%d)\n", curr_read, NWaligned);}
                //printf("Read: %"PRIu64" yielded (%d)\n", curr_read, NWaligned); if(NWaligned == 0) getchar();
                NWaligned = 0;
                //fprintf(stdout, "Seq %"PRIu64" has %"PRIu64" hits and tried to align %"PRIu64" times\n", curr_read, n_hits, alignments_tried);
                //fflush(stdout);
                n_hits = 0;
                already_aligned = FALSE;
                alignments_tried = 0;
                last_diagonal = INT64_MIN; // This is not perfect but if the diagonal reaches the value then we have an overflow anyway
                qf.x_start = 0;
                qf.t_len = 0;

                //if(hta->full_comp == TRUE) memset(&hta->markers[my_current_task->r1], 0, my_current_task->r2 - my_current_task->r1 + 1); // Reset used tags
                memset(&hta->markers[0], FALSE, hta->database->n_seqs); // Reset used tags
                
                curr_read++;
                //printf("On current read %"PRIu64"\n", curr_read);
                continue;
            }

            if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
                curr_kmer[crrSeqL] = (unsigned char) c;
                crrSeqL++;
            }else{
                crrSeqL = 0;
            }

            if (crrSeqL >= custom_kmer) { // Full well formed sequence
            
                
                //printf("comparing hit: %.11s\n", (char *)&curr_kmer[0]);
                //getchar();

                // Choose table
                /*
                if(curr_kmer[0] == (unsigned char) 'A'){
                    ptr_table_redirect = hta->container_A;
                }else if(curr_kmer[0] == (unsigned char) 'C'){
                    ptr_table_redirect = hta->container_C;
                }else if(curr_kmer[0] == (unsigned char) 'G'){
                    ptr_table_redirect = hta->container_G;
                }else{
                    ptr_table_redirect = hta->container_T;
                }
                */

                //fprintf(stdout, "%s\n", curr_kmer);
                //fflush(stdout);
                current_table = 0;
                aux = ptr_table_redirect[current_table]->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                if(aux == NULL){
                    ++current_table;
                    aux = ptr_table_redirect[current_table]->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
                }
                if(aux == NULL){
                    ++current_table;
                    aux = ptr_table_redirect[current_table]->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
                }
                if(aux == NULL){
                    ++current_table;
                    aux = ptr_table_redirect[current_table]->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
                }

                //While there are hits
                //fprintf(stdout, "%p\n", aux);
                //fflush(stdout);

                

                while(aux != NULL && aux->extended_hash == hashOfWord(&curr_kmer[FIXED_K], custom_kmer - FIXED_K) && ((hta->full_comp == FALSE && NWaligned == 0 && hta->markers[aux->s_id+ hta->contained_reads[current_table]] == 0) || (hta->full_comp && hta->markers[aux->s_id+ hta->contained_reads[current_table]] == 0))){

                    n_hits++;
                    //fprintf(stdout, "%p\n", aux);
                    //fflush(stdout);
                    // ADD OFFSET CUCOOOOOOOOOOOOOOOOOOOOOO!!!!!!!!!!!!!!!
                    //printf("my current table is %u\n", current_table);
                    //printf("check this woop: %"PRIu64" - %"PRIu64" - %"PRIu64" - %"PRIu64"\n", aux->s_id, hta->contained_reads[current_table], aux->pos, hta->base_coordinates[current_table]);
                    curr_db_seq = aux->s_id + hta->contained_reads[current_table];
                    pos_of_hit = aux->pos + hta->base_coordinates[current_table];
                    if(hta->hits != NULL){
                        hta->hits[curr_db_seq]++;
                        goto only_hits; // Count only hits and skip the rest
                    } 
                    
                    //fprintf(stdout, "Launching curr_read: %"PRIu64" @ %"PRIu64", vs curr_db_read: %"PRIu64" @ %"PRIu64": ", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                    /*
                    if(curr_read == 534) fprintf(stdout, "Launching %"PRIu64" @ %"PRIu64", vs %"PRIu64" @ %"PRIu64": ", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                    */
                    #ifdef VERBOSE 
                    fprintf(stdout, "Launching %"PRIu64" @ %"PRIu64", vs %"PRIu64" @ %"PRIu64": ", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                    #endif 
                    int64_t curr_diagonal = (int64_t)(curr_pos+1) - (int64_t) pos_of_hit;

                    
                    
                    if( (last_diagonal != curr_diagonal && !(qf.x_start <= (pos_of_hit + custom_kmer) && pos_of_hit <= (qf.x_start + qf.t_len)))){
                        
                        /*
                        if(curr_db_seq == hta->database->n_seqs-1){
                            xlen = hta->database->total_len - hta->database->start_pos[curr_db_seq];
                        }else{
                            xlen = hta->database->start_pos[curr_db_seq+1] - hta->database->start_pos[curr_db_seq];
                        }
                        if(curr_read == hta->query->n_seqs-1){
                            ylen = hta->query->total_len - hta->query->start_pos[curr_read];
                        }else{
                            ylen = hta->query->start_pos[curr_read+1] - hta->query->start_pos[curr_read];
                        }
                        */
                        
                        /*if(curr_read == 35){
                            fprintf(stdout, "Launching %"PRIu64" @ %"PRIu64", vs %"PRIu64" @ %"PRIu64": \n", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                            printf("what do you think will happen: %"PRIu64"\n", hta->database->start_pos[curr_db_seq]);
                            //fprintf(stdout, "lengths: x: %"PRIu64", y: %"PRIu64"\n", xlen, ylen);
                            getchar();
                        }*/
                        

                        //printf("accepted because: \n");
                        /*if(curr_read % 100 == 0){
                            fprintf(stdout, "Launching %"PRIu64" @ %"PRIu64", vs %"PRIu64" @ %"PRIu64": \n", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                        }
                        */
                        //printf("prev_diag: %"PRId64"- currdiag: %"PRId64"\n", last_diagonal, curr_diagonal);
                        //printf("\t covers to [%"PRIu64"+%"PRIu64"=%"PRIu64"] [%"PRIu64"] \n", qf.x_start, qf.t_len, qf.x_start+qf.t_len, pos_of_hit); getchar();
                        if(hta->markers[curr_db_seq] == FALSE){

                            //if(current_table > 0) getchar();
                            alignmentFromQuickHits(hta->database, hta->query, pos_of_hit, curr_pos+1, curr_read, curr_db_seq, &qf, hta->contained_reads[current_table], hta->base_coordinates[current_table]);
                            last_diagonal = curr_diagonal;
                        }else{
                            qf.e_value = 100000000;
                        }
                        
                        //if(curr_read == 35) printf(" evalue: %Le from %"PRIu64", %"PRIu64" with l: %"PRIu64"\n", qf.e_value, qf.x_start, qf.y_start, qf.t_len);
                    }else{
                        
                        /*if(curr_read == 35){
                            printf("rejected because: \n");
                            fprintf(stdout, "UNLaunching %"PRIu64" @ %"PRIu64", vs %"PRIu64" @ %"PRIu64": \n", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                            printf("prev_diag: %"PRId64"- currdiag: %"PRId64"\n", last_diagonal, curr_diagonal);
                            printf("\t covers to [%"PRIu64"+%"PRIu64"=%"PRIu64"] [%"PRIu64"] \n", qf.x_start, qf.t_len, qf.x_start+qf.t_len, pos_of_hit); getchar();
                        }*/
                        
                        qf.e_value = 100000000;
                    }

                    #ifdef VERBOSE 
                    printf(" evalue: %Le %"PRIu64"\n", qf.e_value, qf.t_len);
                    #endif
                    //getchar();


                    


                    //If e-value of current frag is good, then we compute a good gapped alignment
                    if(qf.e_value < hta->min_e_value /*&& xlen == 799 && ylen == 2497*/){
                        alignments_tried++;
                        ba.identities = ba.length = ba.igaps = ba.egaps = 0;
                        //Compute lengths of reads
                        if(curr_db_seq == hta->database->n_seqs-1){
                            xlen = hta->database->total_len - hta->database->start_pos[curr_db_seq];
                        }else{
                            xlen = hta->database->start_pos[curr_db_seq+1] - hta->database->start_pos[curr_db_seq];
                            //printf("!!!\n%"PRIu64", %"PRIu64" :: %"PRIu64"; its db->start_pos[curr_db_seq]  db->start_pos[curr_db_seq+1]  curr_db_seq\n", hta->database->start_pos[curr_db_seq], hta->database->start_pos[curr_db_seq+1], curr_db_seq);
                        }
                        if(curr_read == hta->query->n_seqs-1){
                            ylen = hta->query->total_len - hta->query->start_pos[curr_read];
                        }else{
                            ylen = hta->query->start_pos[curr_read+1] - hta->query->start_pos[curr_read];
                        }
                        //fprintf(stdout, "lengths: x: %"PRIu64", y: %"PRIu64"\n", xlen, ylen);
                        //Perform alignment plus backtracking
                        //void build_alignment(char * reconstruct_X, char * reconstruct_Y, uint64_t curr_db_seq, uint64_t curr_read, HashTableArgs * hta, char * my_x, char * my_y, struct cell ** table, struct cell * mc, char * writing_buffer_alignment, BasicAlignment * ba, uint64_t xlen, uint64_t ylen)
                        if(xlen > MAX_READ_SIZE || ylen > MAX_READ_SIZE){ printf("(%"PRIu64",%"PRIu64")\n", xlen, ylen); terror("Read size reached for gapped alignment."); }
                        //fprintf(stdout, "R0 %"PRIu64", %"PRIu64"\n", curr_db_seq, curr_read);
                        
                        
                        #ifdef VERBOSE 
                        fprintf(stdout, "qfxs %"PRIu64", dbs %"PRIu64", qfys %"PRIu64" qys %"PRIu64"\n", qf.x_start, hta->database->start_pos[curr_db_seq], qf.y_start, hta->query->start_pos[curr_read]);
                        #endif

                        //fprintf(stdout, "dbFragxs %"PRIu64", dbs %"PRIu64", rFragys %"PRIu64" rys %"PRIu64"\n", qf.x_start, hta->database->start_pos[curr_db_seq], qf.y_start, hta->query->start_pos[curr_read]);
                        /*
                        printf("at table: %u\n", current_table);
                        fprintf(stdout, "Launching curr_read: %"PRIu64" @ %"PRIu64", vs curr_db_read: %"PRIu64" @ %"PRIu64": ", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                        fprintf(stdout, "Launching NW %"PRIu64" @ %"PRIu64", vs %"PRIu64" @ %"PRIu64": \n", curr_read, curr_pos+1, curr_db_seq, pos_of_hit);
                        printf("have len: %"PRIu64", %"PRIu64"\n", xlen, ylen);
                        printf("Quickfrag (xs, ys): qf::%"PRIu64", %"PRIu64", tlen:%"PRIu64"\n", qf.x_start, qf.y_start, qf.t_len);
                        printf("Yea, but start and end of read in db is: %"PRIu64" - %"PRIu64"\n", hta->database->start_pos[curr_db_seq], hta->database->start_pos[curr_db_seq+1]);
                        */

                        p1.x = qf.x_start - hta->database->start_pos[curr_db_seq];
                        //p1.y = qf.y_start - hta->query->start_pos[curr_read];
                        p1.y = qf.y_start - (hta->query->start_pos[curr_read] -1);
                        p2.x = p1.x + qf.t_len;
                        p2.y = p1.y + qf.t_len;
                        p3.x = xlen;
                        p3.y = ylen;

                        #ifdef VERBOSE
                        fprintf(stdout, "p0 (%"PRIu64", %"PRIu64") p1 (%"PRIu64", %"PRIu64") p2 (%"PRIu64", %"PRIu64") p3 (%"PRIu64", %"PRIu64")\n", p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
                        #endif
                        //fprintf(stdout, "p0 (%"PRIu64", %"PRIu64") p1 (%"PRIu64", %"PRIu64") p2 (%"PRIu64", %"PRIu64") p3 (%"PRIu64", %"PRIu64")\n", p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
                        
                        calculate_y_cell_path(p0, p1, p2, p3, cell_path_y);

                        // REMOVE
                        /*
                        uint64_t r1,r2;
                        for(r1=0;r1<MAX_WINDOW_SIZE;r1++){
                            for(r2=0;r2<MAX_WINDOW_SIZE;r2++){
                                hta->table[r1][r2].score = INT64_MIN;
                            }
                        }
                        */
                        
                        build_alignment(hta->reconstruct_X, hta->reconstruct_Y, curr_db_seq, curr_read, hta, hta->my_x, hta->my_y, hta->table, hta->mc, hta->writing_buffer_alignment, &ba, xlen, ylen, cell_path_y, &hta->window);
                        
                        // Set the read to already aligned so that it does not repeat
                        hta->markers[curr_db_seq] = 1;
                        
                        #ifdef VERBOSE
                        printf("len 1 %"PRIu64", len 2 %"PRIu64"\n", ba.length, ylen);
                        printf("ident %"PRIu64"\n", ba.identities);
                        #endif

                        //If is good
                        if(((long double)(ba.length-(ba.igaps+ba.egaps))/ylen) >= hta->min_coverage && ((long double)ba.identities/(ba.length-(ba.igaps+ba.egaps))) >=  hta->min_identity){
                            if(already_aligned == FALSE){
                                hta->accepted_query_reads++;
                                already_aligned = TRUE;
                                //printf("accepted: %"PRIu64"\n", hta->accepted_query_reads);
                            }
                            
                            hta->markers[curr_db_seq] = 1;
                            if(hta->out != NULL){
                                //printf("Last was: (%"PRIu64", %"PRIu64")\n", curr_read, curr_db_seq);
                                fprintf(hta->out, "(%"PRIu64", %"PRIu64") : %d%% %d%% %"PRIu64"\n $$$$$$$ \n", curr_read, curr_db_seq, MIN(100,(int)(100*(ba.length-(ba.igaps+ba.egaps))/ylen)), MIN(100,(int)((long double)100*ba.identities/(ba.length-(ba.igaps+ba.egaps)))), ylen);
                                fprintf(hta->out, "%s", hta->writing_buffer_alignment);
                                //fprintf(stdout, "(%"PRIu64", %"PRIu64") : %d%% %d%% %"PRIu64"\n $$$$$$$ \n", curr_read, curr_db_seq, MIN(100,(int)(100*ba.identities/ba.length)), MIN(100,(int)(100*ba.length/ylen)), ylen);
                                //fprintf(stdout, "%s", hta->writing_buffer_alignment);
                            }
                            NWaligned = 1;
                        }/*else{
                            printf("what: ");
                            printf("len x %"PRIu64", len y %"PRIu64"\n", xlen, ylen);
                            printf("ident %"PRIu64" len %"PRIu64"\n", ba.identities, ba.length); getchar();
                        }*/
                    
                    }

                    //strncpy(get_from_db, &hta->database->sequences[qf.x_start], qf.t_len);
                    //strncpy(get_from_query, &hta->query->sequences[qf.y_start], qf.t_len);
                    //fprintf(hta->out, "%s\n%s\n%Le\t%d\n-------------------\n", get_from_db, get_from_query, qf.e_value, (int)(100*qf.coverage));
                    //fprintf(hta->out, "%"PRIu64", %"PRIu64", %"PRIu64"\n", qf.x_start, qf.y_start, qf.t_len);

                    //printf("Hit comes from %"PRIu64", %"PRIu64"\n", pos_of_hit, curr_pos);
                    only_hits:
                    aux = aux->next;
                    while(aux == NULL && current_table < FIXED_LOADING_THREADS-1){
                        ++current_table;
                        aux = ptr_table_redirect[current_table]->table[char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]][char_converter[curr_kmer[3]]]
                        [char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]][char_converter[curr_kmer[6]]]
                        [char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]][char_converter[curr_kmer[9]]]
                        [char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];
                    }
                    //fprintf(stdout, "%p\n", aux);
                    //fflush(stdout);
                }
                //printf("SWITCHED\n");

                if(NWaligned == 1 && hta->full_comp == FALSE){
                    if(curr_read < hta->query->n_seqs) curr_pos = hta->query->start_pos[curr_read+1]-2;
                }else{
                    memcpy(b_aux, curr_kmer, custom_kmer);
                    memcpy(curr_kmer, &b_aux[1], custom_kmer-1);
                    crrSeqL -= 1;
                }
            }
            //printf("current pos: %"PRIu64"\n", curr_pos);
            curr_pos++;
            if(curr_pos < hta->query->total_len) c = (char) hta->query->sequences[curr_pos];
            


        }





    }
    
    
    

    
    //fprintf(stdout, "Going from %"PRIu64" to %"PRIu64"\n", hta->from, hta->to);
    //fflush(stdout);

    free(cell_path_y);

    return NULL;

}

void build_alignment(char * reconstruct_X, char * reconstruct_Y, uint64_t curr_db_seq, uint64_t curr_read, HashTableArgs * hta, unsigned char * my_x, unsigned char * my_y, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, BasicAlignment * ba, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window){
 

    //Do some printing of alignments here
    uint64_t maximum_len, i, j, curr_pos_buffer, curr_window_size;

    maximum_len = 2*MAX(xlen,ylen);
    memcpy(my_x, &hta->database->sequences[hta->database->start_pos[curr_db_seq]], xlen);
    memcpy(my_y, &hta->query->sequences[hta->query->start_pos[curr_read]], ylen);

    struct best_cell bc = NW(my_x, 0, xlen, my_y, 0, ylen, (int64_t) hta->igap, (int64_t) hta->egap, table, mc, 0, cell_path_y, window, &curr_window_size);
    backtrackingNW(my_x, 0, xlen, my_y, 0, ylen, table, reconstruct_X, reconstruct_Y, &bc, &i, &j, ba, cell_path_y, curr_window_size);
    uint64_t offset = 0, before_i = 0, before_j = 0;
    i++; j++;
    
    #ifdef VERBOSE
    uint64_t z=0;
    for(z=0;z<maximum_len;z++) printf("%c", reconstruct_X[z]);
    printf("\n");
    for(z=0;z<maximum_len;z++) printf("%c", reconstruct_Y[z]);
    #endif
    
    curr_pos_buffer = 0;
    while(i <= maximum_len && j <= maximum_len){
        offset = 0;
        before_i = i;
        writing_buffer_alignment[curr_pos_buffer++] = 'D';
        writing_buffer_alignment[curr_pos_buffer++] = '\t';
        while(offset < ALIGN_LEN && i <= maximum_len){
            //fprintf(stdout, "%c", reconstruct_X[i]);
            writing_buffer_alignment[curr_pos_buffer++] = (char) reconstruct_X[i];
            i++;
            offset++;
        }
        //fprintf(out, "\n");
        
        writing_buffer_alignment[curr_pos_buffer++] = '\n';
        offset = 0;
        before_j = j;
        
        //fprintf(stdout, "\n");

        writing_buffer_alignment[curr_pos_buffer++] = 'Q';
        writing_buffer_alignment[curr_pos_buffer++] = '\t';

        while(offset < ALIGN_LEN && j <= maximum_len){
            //fprintf(stdout, "%c", reconstruct_Y[j]);
            writing_buffer_alignment[curr_pos_buffer++] = (char) reconstruct_Y[j];
            j++;
            offset++;
        }
        //fprintf(out, "\n");
        writing_buffer_alignment[curr_pos_buffer++] = '\n';
        writing_buffer_alignment[curr_pos_buffer++] = ' ';
        writing_buffer_alignment[curr_pos_buffer++] = '\t';
        while(before_i < i){
            if(reconstruct_X[before_i] != '-' && reconstruct_Y[before_j] != '-' && reconstruct_X[before_i] == reconstruct_Y[before_j]){
                //fprintf(out, "*");
                writing_buffer_alignment[curr_pos_buffer++] = '*';
                ba->identities++;
            }else{
                //fprintf(out, " ");
                writing_buffer_alignment[curr_pos_buffer++] = ' ';
            }
            before_j++;
            before_i++;
        }
        writing_buffer_alignment[curr_pos_buffer++] = '\n';

    }
    //fprintf(out, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    writing_buffer_alignment[curr_pos_buffer++] = '\n';
    writing_buffer_alignment[curr_pos_buffer++] = '\0';


}

void alignmentFromQuickHits(SeqInfo * database, SeqInfo * query, uint64_t pos_database, uint64_t pos_query, uint64_t curr_read, uint64_t curr_db_seq, Quickfrag * qf, uint64_t offset_db_reads, uint64_t offset_db_coordinates){

    int64_t read_x_start, read_x_end, read_y_start, read_y_end;

    if(curr_db_seq == database->n_seqs-1){
        read_x_start = database->start_pos[curr_db_seq];
	    read_x_end = database->total_len;
    }else{     
        read_x_start = database->start_pos[curr_db_seq];
	    read_x_end = database->start_pos[curr_db_seq+1] - 1;
    }

    //printf("read x start -> %"PRId64", end -> %"PRId64" btw: %"PRIu64"\n", read_x_start, read_x_end, database->n_seqs);

    if(curr_read == query->n_seqs-1){
        read_y_start = query->start_pos[curr_read];
        read_y_end = query->total_len;
    }else{
        read_y_start = query->start_pos[curr_read];
        read_y_end = query->start_pos[curr_read+1] - 1;
    }

    //printf("db_end %"PRId64" query_end %"PRId64"\n", read_x_end, read_y_end);
    //printf("pos database: %"PRIu64"\n", pos_database);
    int64_t curr_pos_db = (int64_t) pos_database;
    int64_t curr_pos_qy = (int64_t) pos_query;
    int64_t final_end_x = (int64_t) pos_database - 1, final_start_x = final_end_x - custom_kmer + 1, final_start_y = pos_query - custom_kmer;
    int64_t score_right = custom_kmer * POINT;
    int64_t score_left = score_right;
    int64_t high_left = score_left, high_right = score_right;
    qf->t_len = custom_kmer;
    uint64_t idents = custom_kmer;

    /*
    char le_hit[1000];
    memcpy(le_hit, &database->sequences[final_start_x], FIXED_K);
    
	fprintf(stdout, "HIT: %s\n", le_hit);
	fflush(stdout);
    */
    //printf("final start x: %"PRId64"\n", final_start_x);
    int keep_going = 1;

    //Forward search
    while(keep_going == 1){
        
        
        if(score_right > 0 && curr_pos_db < database->total_len && curr_pos_qy < query->total_len){
            if(curr_pos_db  > read_x_end ||  curr_pos_qy > read_y_end) break;
            //if(database->sequences[curr_pos_db] == query->sequences[curr_pos_qy]){ score_right+=POINT; idents++; }else{ score_right-=POINT;}
            if(compare_letters(database->sequences[curr_pos_db], query->sequences[curr_pos_qy]) == POINT){ score_right+=POINT; idents++; }else{ score_right-=POINT;}
            if(high_right <= score_right){
                final_end_x = curr_pos_db;
                high_right = score_right;
            }
            curr_pos_db++;
            curr_pos_qy++;
        }else{
            keep_going = 0;
        }
    }

    //printf("pos here %"PRIu64" curr_pos_db, curr_pos_query %"PRIu64"\n", curr_pos_db, curr_pos_qy);
    //printf("final start x: %"PRId64"\n", final_start_x);
    keep_going = 1;
    curr_pos_db = pos_database - custom_kmer - 1;
    curr_pos_qy = pos_query - custom_kmer - 1;

    score_left = high_right;

    //Backward search
    while(keep_going == 1){
        
        if(score_left > 0 && curr_pos_db >= 0 && curr_pos_qy >= 0){
            if(curr_pos_db < read_x_start || curr_pos_qy < read_y_start ) break;
            //if(database->sequences[curr_pos_db] == query->sequences[curr_pos_qy]){ score_left+=POINT; idents++; }else{ score_left-=POINT;}
            if(compare_letters(database->sequences[curr_pos_db], query->sequences[curr_pos_qy]) == POINT){ score_left+=POINT; idents++; }else{ score_left-=POINT;}
            if(high_left <= score_left){
                final_start_x = curr_pos_db;
                final_start_y = curr_pos_qy;
                //printf("got %"PRIu64" when min is %"PRIu64"\n", final_start_y, read_y_start);
                high_left = score_left;
            }
            curr_pos_db--;
            curr_pos_qy--;
        }else{
            keep_going = 0;
        }
    }

    qf->t_len = final_end_x - final_start_x;
    
    /*
    char s1[1000];
    char s2[1000];

    memcpy(s1, &database->sequences[final_start_x], qf->t_len);
    memcpy(s2, &query->sequences[final_start_y], qf->t_len);
    s1[qf->t_len] = '\0';
    s2[qf->t_len] = '\0';
    
	fprintf(stdout, "%s\n%s\n------\n", s1, s2);
	fflush(stdout);

    printf("the real hit was:\n");
    memcpy(s1, &database->sequences[pos_database-12+1], 12);
    memcpy(s2, &query->sequences[pos_query-12+1], 12);
    s1[12] = '\0';
    s2[12] = '\0';
    
	fprintf(stdout, "%s\n%s\n------\n", s1, s2);
	fflush(stdout);
    //getchar();
    printf("%"PRIu64"\n", idents);
    */
    
    long double rawscore = (idents*POINT) - (qf->t_len - idents)*(POINT);

    long double t_len;
    if(curr_read == query->n_seqs-1){
        t_len = (long double) query->total_len - query->start_pos[curr_read];
    }else{
        t_len = (long double) query->start_pos[curr_read+1] - query->start_pos[curr_read];
    }
    //printf("final start x: %"PRId64"\n", final_start_x);
    qf->x_start = final_start_x;
    qf->y_start = final_start_y;
    qf->e_value = (long double) QF_KARLIN*t_len*database->total_len*expl(-QF_LAMBDA * rawscore);
    qf->coverage = qf->t_len / t_len;

}

void calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3, int64_t * y_points){
    
    //Calculate lines between points
    uint64_t i;

    #ifdef VERBOSE
    printf("Built on\n");
    printf("(%"PRIu64", %"PRIu64")\n", p0.x, p0.y);
    printf("(%"PRIu64", %"PRIu64")\n", p1.x, p1.y);
    printf("(%"PRIu64", %"PRIu64")\n", p2.x, p2.y);
    printf("(%"PRIu64", %"PRIu64")\n", p3.x, p3.y);
    #endif

    

    if(p0.x > MAX_READ_SIZE){ fprintf(stdout, "LEN error %"PRIu64"\n", p0.x); terror("Reached max length in read for anchoring procedure (1)"); }
    if(p1.x > MAX_READ_SIZE){ fprintf(stdout, "LEN error %"PRIu64"\n", p1.x); terror("Reached max length in read for anchoring procedure (2)"); }
    if(p2.x > MAX_READ_SIZE){ fprintf(stdout, "LEN error %"PRIu64"\n", p2.x); terror("Reached max length in read for anchoring procedure (3)"); }
    if(p3.x > MAX_READ_SIZE){ fprintf(stdout, "LEN error %"PRIu64"\n", p3.x); terror("Reached max length in read for anchoring procedure (4)"); }

    long double deltax, deltay, deltaerr, error;
    uint64_t y;

    //P0 to P1
    deltax = p1.x - p0.x;
    deltay = p1.y - p0.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p0.y;

    for(i=p0.x;i<p1.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }

    //P1 to P2

    deltax = p2.x - p1.x;
    deltay = p2.y - p1.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p1.y;

    for(i=p1.x;i<p2.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }
    
    //P2 to P3

    deltax = p3.x - p2.x;
    deltay = p3.y - p2.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p2.y;

    for(i=p2.x;i<p3.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }

    /*
    if(p3.x == 799 && p3.y == 2497){
        for(i=0;i<p3.x;i++){
            printf("%"PRIu64": %"PRIu64"\n", i, y_points[i]);
            if(i % 50 == 0) getchar();
        }    
    }
    */

    #ifdef VERBOSE
    for(i=0;i<p3.x;i++){
        printf("%"PRIu64" -> ", y_points[i]);
        if(i % 50 == 0) getchar();
    }
    
    #endif
    
    

}

struct best_cell NW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show, int64_t * cell_path_y, long double * window, uint64_t * current_window_size){
    

    uint64_t i, j, j_prime;
    int64_t scoreDiagonal = INT64_MIN, scoreLeft = INT64_MIN, scoreRight = INT64_MIN, score = INT64_MIN, delta_dif_1_0, delta_dif_2_1, limit_left, limit_right, j_right_prime = 1, j_left_prime = 1, j_diag_prime = 1;

    struct best_cell bc;
    bc.c.score = INT64_MIN;
    bc.c.xpos = 0; bc.c.ypos = 0;
    
    //The window size will be a +-15% of the square root of the product of lengths
    int64_t window_size = MIN(MAX_WINDOW_SIZE/2, (uint64_t) (*window * sqrtl((long double) Xend * (long double) Yend)));
    //printf("xlen: %"PRIu64", ylen: %"PRIu64" w-size: %"PRId64"\n", Xend, Yend, window_size);    
    *current_window_size = (uint64_t) window_size;

    //The limits to the window
    limit_left = 0;
    limit_right = 2*window_size + 1;
    if(limit_right > MAX_WINDOW_SIZE) limit_right = MAX_WINDOW_SIZE;
    
    struct positioned_cell mf;
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");
    //Zero will always be
    table[0][0].score = compare_letters(X[0], Y[0]);
    mc[0].score = table[0][0].score;
    mc[0].xpos = 0;
    mc[0].ypos = 0;
    
    //if(Xend == 799 && Yend == 2497) printf("I am %p The count is real %.5s %.5s %p %p \n", &table[0][0], X, Y, X, Y);

    
    for(i=1;i<Yend;i++){
        //table[0][i].score = (X[0] == Y[i]) ? POINT : -POINT;
        if(i < MAX_WINDOW_SIZE) table[0][i].score = compare_letters(X[0], Y[i]) + iGap + (i-1)*eGap;
        //table[Xstart][i].xfrom = Xstart;
        //table[Xstart][i].yfrom = i;
        //Set every column max
        mc[i].score = compare_letters(X[0], Y[i]) + iGap + (i-1)*eGap;
        #ifdef VERBOSE
        printf("%02"PRId64" ", mc[i].score);
        #endif
        mc[i].xpos = 0;
        mc[i].ypos = i;

    }
    #ifdef VERBOSE
    printf("\n");
    #endif
    //Set row max
    mf.score = table[0][0].score;
    mf.xpos = 0;
    mf.ypos = 0;
    //Init j
    j = MAX(1,(cell_path_y[1] - window_size));

    //Go through full matrix
    for(i=1;i<Xend;i++){
        //Fill first rowcell
        if(cell_path_y[i-1]+window_size < cell_path_y[i]) return bc; //terror("Sequence proportions make window shift too large");
        //Conversion for the j-coordinate
        j_prime = 1;

        //table[i][0].score = (X[i] == Y[0]) ? POINT : -POINT;
        if(cell_path_y[i] - window_size <= 0){
            table[i][0].score = compare_letters(X[i], Y[0]) + iGap + (i-1)*eGap;
            mf.score = table[i][0].score;
        }else{
            mf.score = compare_letters(X[i], Y[0]) + iGap + (i-1)*eGap;
        }

        mf.xpos = i-1;
        mf.ypos = 0;

        delta_dif_1_0 = MAX(1, (cell_path_y[i] - window_size)) - MAX(1,(cell_path_y[i-1] - window_size)); //j-1
        if(i>1) delta_dif_2_1 = MAX(1, (cell_path_y[i-1] - window_size)) - MAX(1, (cell_path_y[i-2] - window_size)); //j-2

        #ifdef VERBOSE 
        printf("D1_0: %"PRId64" D2_1: %"PRId64"\n", delta_dif_1_0, delta_dif_2_1);
        #endif

        #ifdef VERBOSE
        printf("%02"PRId64" ", mf.score);
        #endif
        //printf("Check on i: (%"PRIu64") from - to (%"PRIu64", %"PRIu64")\n", i, 0L, Xend);
        /*
        if(1||i==262){
            printf("I will go from %"PRIu64" to %"PRIu64" and I am %"PRIu64", %"PRIu64"\n", (uint64_t) MAX(1,(cell_path_y[i] - (int64_t)window_size)), (uint64_t) MIN((int64_t)Yend,(cell_path_y[i] + (int64_t)window_size)), i, j);
            //printf("lengs: %"PRIu64", %"PRIu64"\n", Xend, Yend);
            //printf("cp[i]: %"PRId64", cp[i-1] %"PRId64"\n", cell_path_y[i], cell_path_y[i-1]);
            //printf("min(%"PRId64", %"PRId64" + %"PRId64")-------------------\n", Yend, cell_path_y[i] ,(int64_t)window_size);
            
        }
        */
        //getchar();

        //printf("@%"PRIu64"[%"PRId64"] -> (%"PRIu64", %"PRIu64") jp %"PRIu64", lright %"PRIu64"\n", i, cell_path_y[i], MAX(1,(cell_path_y[i] - window_size)), MIN(Yend,(cell_path_y[i] + window_size)), j_prime, limit_right);
        //printf("M:@%"PRIu64"-> %"PRIu64"\n", i, MIN(Yend,(cell_path_y[i] + window_size)));
        #ifdef VERBOSE
        int64_t r;
        for(r=0;r<MAX(0,(cell_path_y[i] - window_size)); r++){
            printf("  ");
        }
        #endif

        /*
        if(Xend == 799 && Yend == 2497 && i >= 145 && i <= 155){
                printf("them limits @i %"PRIu64"::: %"PRIu64", %"PRIu64"\n", i, MAX(1,(cell_path_y[i] - window_size)), MIN(Yend,(cell_path_y[i] + window_size)));
                getchar();
        }
        */
        

        for(j=MAX(1,(cell_path_y[i] - window_size));j<MIN(Yend,(cell_path_y[i] + window_size)) && j_prime < limit_right;j++){
            //if(i == 8302){ printf("Doing on : (%"PRIu64",%"PRIu64" and jprime=%"PRIu64"\n", i,j,j_prime); getchar(); }
            //Check if max in row has changed
            //if(j > MAX(1, cell_path_y[i-1] - window_size +1) && mf.score <= table[i][j-2].score){
            //if(j_prime == MAX_WINDOW_SIZE) break;
            //Calculate the real j position in the windowed table
            
            j_left_prime = ((int64_t)j_prime - (2 - delta_dif_1_0));
            //j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            if(i > 1){
                j_right_prime = ((int64_t)j_prime - (1 - (delta_dif_1_0 + delta_dif_2_1)));
            }

            if(j > MAX(1, cell_path_y[i-1] - window_size +1) && j < MIN(Yend,(cell_path_y[i-1] + window_size)) && j_left_prime < limit_right && table[i-1][j_left_prime].score >= mf.score){
                //mf.score = table[i-1][j-2].score;
                mf.score = table[i-1][j_left_prime].score;
                mf.xpos = i-1;
                mf.ypos = j-2;
                if(table[i-1][j_left_prime].score == INT64_MIN){ printf("A: mf.x\t%"PRIu64"\tmf.y\t%"PRIu64"\ts%"PRId64"\n", mf.xpos, mf.ypos, mf.score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
                
            }
            //printf("RowMax: %"PRId64"@(%"PRIu64", %"PRIu64")\t", mf.score, mf.xpos, mf.ypos);
            
            //score = (X[i] == Y[j]) ? POINT : -POINT;
            score = compare_letters(X[i], Y[j]);

            //Precondition: Upper row needs to reach up to diagonal
            //if((cell_path_y[i-1]+window_size) >= j-1){
            if(i > 1 && j >= 1 && j-1 >= MAX(1,(cell_path_y[i-2] - window_size)) && j-1 < MIN(Yend,(cell_path_y[i-2] + window_size)) && j_right_prime >= limit_left && j_right_prime < limit_right && table[i-2][j_right_prime].score >= mc[j-1].score ){
                //mc[j-1].score = table[i-2][j-(1+j_prime)].score;
                //Should be the j_prime we had at cell_path_y
                //MAX(1,(cell_path_y[i] - window_size));j<MIN(Yend,(cell_path_y[i] + window_size))
                
                mc[j-1].score = table[i-2][j_right_prime].score;
                mc[j-1].xpos = i-2;
                mc[j-1].ypos = j-1;

                if(table[i-2][j_right_prime].score == INT64_MIN){ printf("A: j-1\t%"PRIu64"\tmc.xpos\t%"PRIu64"\ts%"PRId64"\n", j-1, mc[j-1].xpos, mc[j-1].score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
    
            }
            

            if(j-1 >= MAX(0, (cell_path_y[i-1]-window_size)) && (cell_path_y[i-1]+window_size) >= j-1 && j_diag_prime >= limit_left && j_diag_prime < limit_right && j_diag_prime < cell_path_y[i-1]+window_size){
                //scoreDiagonal = table[i-1][j-1].score + score;
                //printf("prevdiag: %"PRId64"\n", table[i-1][j_diag_prime].score);
                scoreDiagonal = table[i-1][j_diag_prime].score + score;                
                if(table[i-1][j_diag_prime].score == INT64_MIN){ printf("A: i-1\t%"PRIu64"\tj_diag\t%"PRIu64"\ts%"PRId64"\n", i-1, j_diag_prime, table[i-1][j_diag_prime].score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
            }else{
                scoreDiagonal = INT64_MIN;
            }
            
            if(i>=1 && j>1){
                scoreLeft = mf.score + iGap + (j - (mf.ypos+2))*eGap + score;
                
                if(mf.score == INT64_MIN){ printf("A: mf.x\t%"PRIu64"\tmf.y\t%"PRIu64"\ts%"PRId64"\n", mf.xpos, mf.ypos, mf.score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
            }else{
                scoreLeft = INT64_MIN;
            }

            if(j>=1 && i>1){
                scoreRight = mc[j-1].score + iGap + (i - (mc[j-1].xpos+2))*eGap + score;
                //if(scoreRight == -12) printf("MC: %"PRId64", From: %"PRIu64", %"PRIu64"->", mc[j-1].score, mc[j-1].xpos, mc[j-1].ypos);
                
                if(mc[j-1].score == INT64_MIN){ printf("A: j-1\t%"PRIu64"\tmc.xpos\t%"PRIu64"\ts%"PRId64"\n", j-1, mc[j-1].xpos, mc[j-1].score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
            }else{
                scoreRight = INT64_MIN;
            }
            
            /*
            if(Xend == 799 && Yend == 2497 && i >= 152 && i == 153){
                printf("@%"PRIu64", %"PRIu64" -> scores: %"PRId64", %"PRId64", %"PRId64"\n", i, j, scoreDiagonal, scoreRight, scoreLeft);
                printf("in position @ jprime= %"PRIu64" cellpaths [i-2, i-1, i] are %"PRId64", %"PRId64", %"PRId64", window_size: %"PRId64", j_diag_prime: %"PRId64"\n", j_prime, cell_path_y[i-2], cell_path_y[i-1], cell_path_y[i], window_size, j_diag_prime);
                printf("Mfs from scoreLeft: mf.x\t%"PRIu64"\tmf.y\t%"PRIu64"\ts%"PRId64"\n", mf.xpos, mf.ypos, mf.score);
                getchar();
            }
            */

            //Choose maximum
            /*
            #ifdef VERBOSE
            printf("The game starts at %"PRId64"\n", MAX(0, cell_path_y[i] - window_size));
            printf("from %c %c and I get to %"PRIu64" while j=%"PRIu64"\n", X[i], Y[j], j_prime, j);
            printf("j_prime: %"PRId64"\n", j_prime);
            printf("j_diag_prime: %"PRId64" limits[%"PRId64", %"PRId64"]\n", j_diag_prime, limit_left, limit_right);
            printf("Score DIAG: %"PRId64"; LEFT: %"PRId64"; RIGHT: %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
            printf("currmf: %"PRId64" mc: %"PRId64"\n", mf.score, mc[j-1].score);
            #endif
            */
            

            //if(i >= MAX_READ_SIZE){ printf("i=%"PRIu64"\n", i); terror("i overflowed\n");}
            //if(j_prime >= MAX_WINDOW_SIZE){ printf("upper : %"PRId64"\n", MIN(Yend,(cell_path_y[i] + window_size-1))); printf("jp=%"PRIu64"\n", j_prime); terror("j overflowed\n"); }
          


            if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight){
                //Diagonal
                
                //fprintf(stdout, "The JPRIME: %"PRId64" actual pos: %"PRIu64"\n", j_prime, j); getchar();
                table[i][j_prime].score = scoreDiagonal;
                table[i][j_prime].xfrom = i-1;
                table[i][j_prime].yfrom = j-1;
                
                                
            }else if(scoreRight > scoreLeft){
                table[i][j_prime].score = scoreRight;
                table[i][j_prime].xfrom = mc[j-1].xpos;
                table[i][j_prime].yfrom = mc[j-1].ypos;
                
            }else{
                //printf("Scores %"PRId64", %"PRId64", %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
                table[i][j_prime].score = scoreLeft;
                table[i][j_prime].xfrom = mf.xpos;
                table[i][j_prime].yfrom = mf.ypos;
            }
            //printf("F: i\t%"PRIu64"\tj_prime\t%"PRIu64"\n", i, j_prime);
            //getchar();
            //if(i == 94){ printf("showing j %"PRIu64" jprime %"PRIu64" lleft %"PRIu64", llright %"PRIu64"\n", j, j_prime, limit_left, limit_right); getchar(); }
            //if(i == 94 && j == 374){ printf("stopped at 94, 374 s %"PRId64"\n", table[i][j_prime].score); getchar(); }
            
            
                       

            /*
            if(i == 264 && j == 176){
                    printf("@%"PRIu64", %"PRIu64"\n", i, j);
                    printf("my score is %"PRId64"\n", mc[j-1].score);
                    printf("in position @ jprime= %"PRIu64" cellpaths [i-1, i] are %"PRId64", %"PRId64"\n", j_prime, cell_path_y[i-1], cell_path_y[i]);
                    printf("Scores %"PRId64", %"PRId64", %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
                    printf("check j_right_prime == %"PRIu64"\n", j_right_prime);
                    getchar();
                    //exit(-1);
            }
            */
            
            //check if column max has changed
            //New condition: check if you filled i-2, j-1
            
            
            if(i == Xend-1 || j == Yend-1){

                if(i == Xend-1 && j != Yend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Yend - j)*eGap;
            	}else if(j == Yend-1 && i != Xend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Xend - i)*eGap;
            	}
                //Check for best cell
                if(table[i][j_prime].score >= bc.c.score){ 
                    
                    /*
                    if(i == 798 && j == 1052){ // yields 799, 2497
                        printf("in position @ jprime= %"PRIu64" cellpaths [i-1, i] are %"PRId64", %"PRId64"\n", j_prime, cell_path_y[i-1], cell_path_y[i]);
                        printf("Scores %"PRId64", %"PRId64", %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
                        printf("score comes from %"PRIu64", %"PRIu64",\n", mc[j-1].xpos, mc[j-1].ypos);
                        printf("IDlengths: %"PRIu64", %"PRIu64"\n", Xend, Yend);
                        
                        //exit(-1);
                    }
                    */
                    
                    bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime; 
                }
                //bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime;
            }
            
            
            #ifdef VERBOSE
            //printf("Put score: %"PRId64"\n\n", table[i][j_prime].score);
            //printf("(%"PRId64")%02"PRId64" ", j_diag_prime, table[i][j_prime].score); //printf("->(%"PRIu64", %"PRIu64")", i, j); printf("[%c %c]", X[i], Y[j]);
            //if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight) printf("*\t");
            //else if(scoreRight > scoreLeft) printf("{\t"); else printf("}\t");
            //getchar();

            #endif
            j_prime++;
        }
        #ifdef VERBOSE
        printf("\n");
        getchar();
        #endif
    }
        
    return bc;
}



void backtrackingNW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, struct cell ** table, char * rec_X, char * rec_Y, struct best_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, BasicAlignment * ba, int64_t * cell_path_y, uint64_t window_size){
    uint64_t curr_x, curr_y, prev_x, prev_y, head_x, head_y, limit_x, limit_y;
    int64_t k, j_prime, delta_diff = 0;

    //limit_x = 2*MAX_READ_SIZE-1;
    //limit_y = limit_x;
    limit_x = 2*MAX(Xend, Yend);
    limit_y = limit_x;
    //head_x = 2*MAX(Xend, Yend);
    //head_y = 2*MAX(Xend, Yend);
    head_x = limit_x;
    head_y = limit_y;
    curr_x = bc->c.xpos;
    curr_y = bc->c.ypos;
    #ifdef VERBOSE
    printf("Optimum : %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    #endif
    //printf("Optimum : %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    
    prev_x = curr_x;
    prev_y = curr_y;
    int show = 0;
   
    for(k=Xend-1; k>curr_x; k--) rec_X[head_x--] = '-';
    for(k=Yend-1; k>curr_y; k--) rec_Y[head_y--] = '-';

    j_prime = bc->j_prime;
    //printf("init prime: %"PRIu64"\n", j_prime);
    unsigned char first_track = 1;
    
    while(curr_x > 0 && curr_y > 0){

        
        if(first_track == 0){
            delta_diff = MAX(1, cell_path_y[prev_x] - (int64_t) window_size) - MAX(1, cell_path_y[curr_x] - (int64_t)window_size); //j-1
            j_prime = MAX(0, j_prime - (int64_t)(prev_y - curr_y) + (int64_t) delta_diff);
            

            if(/*(bc->c.xpos == 630 && bc->c.ypos == 541 )||*/ j_prime > MAX_WINDOW_SIZE){
                
                printf("from %"PRIu64", %"PRIu64"\nto   %"PRIu64", %"PRIu64"\n", prev_x, prev_y, curr_x, curr_y);
                printf("jp: %"PRIu64", py,cy %"PRIu64", %"PRIu64", delta: %"PRId64"\n", j_prime, prev_y, curr_y, (int64_t)delta_diff);
                printf("currx curry : %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
                printf("window size: %"PRIu64"\n", window_size);
                printf("cp[prev, curr] : %"PRId64", %"PRId64"\n", cell_path_y[prev_x], cell_path_y[curr_x]);
                printf("my cell path: %"PRId64"\n", cell_path_y[curr_x]);
                printf("Optimum : %"PRIu64", %"PRIu64"\n", bc->c.xpos, bc->c.ypos);
                getchar();
            }
            
            //j_prime = j_prime - (int64_t)(prev_y - curr_y) + (int64_t) delta_diff;

            prev_x = curr_x;
            prev_y = curr_y;

            /*
            if(bc->c.xpos == 798 && bc->c.ypos == 1052){
                printf("[%c %c]", X[prev_x], Y[prev_y]);
                printf("(%"PRIu64", %"PRIu64") ::: \n", curr_x, curr_y);

            }
            */

            #ifdef VERBOSE
            //printf("Jprime: %"PRId64" :DELTADIF:%"PRId64"\n", j_prime, delta_diff);
            printf("[%c %c]", X[prev_x], Y[prev_y]);
            printf("(%"PRIu64", %"PRIu64") ::: \n", curr_x, curr_y);
            //printf("(%"PRIu64", %"PRIu64") ::: \n", prev_x, prev_y);
            //printf("cellp Prev: %"PRId64" Post: %"PRId64"\n", cell_path_y[prev_x], cell_path_y[curr_x]);
            //printf("the difs? %"PRId64" the other: %"PRId64"\n", MAX(0, cell_path_y[prev_x] - (int64_t) window_size), MAX(0, cell_path_y[curr_x] - (int64_t)window_size));
            getchar();
            #endif

        }
        
        //if(table[prev_x][j_prime].xfrom > MAX_READ_SIZE || table[prev_x][j_prime].yfrom > MAX_WINDOW_SIZE) fprintf(stdout, "OH NOES !! %"PRIu64"\t%"PRId64"\t%"PRIu64"\t%"PRIu64" dangers: %"PRIu64", %"PRIu64"\n", prev_x, j_prime, Xend, Yend, table[prev_x][j_prime].xfrom, table[prev_x][j_prime].yfrom);

        /*
        if(table[prev_x][j_prime].xfrom > MAX_READ_SIZE || table[prev_x][j_prime].yfrom > MAX_WINDOW_SIZE){
            fprintf(stdout, "OH NOES !! %"PRIu64"\t%"PRId64"\t%"PRIu64"\t%"PRIu64" dangers: %"PRIu64", %"PRIu64"\n", prev_x, j_prime, Xend, Yend, table[prev_x][j_prime].xfrom, table[prev_x][j_prime].yfrom);
            uint64_t k;
            for(k=0;k<Xend;k++){
                fprintf(stdout, "%c", X[k]);
            }
            fprintf(stdout, "\n");
            for(k=0;k<Yend;k++){
                fprintf(stdout, "%c", Y[k]);
            }
            fprintf(stdout, "\n");
            show = 1;
        } 
        */
        if(j_prime >= MAX_WINDOW_SIZE) printf("j_prime:overflow %"PRIu64"\n", j_prime);

 
        curr_x = table[prev_x][j_prime].xfrom;
        curr_y = table[prev_x][j_prime].yfrom;
        first_track = 0;
        
        //printf("w: %"PRIu64"- %"PRIu64"\n", curr_x, curr_y);
        

        if((curr_x == (prev_x - 1)) && (curr_y == (prev_y -1))){
            //Diagonal case
            //printf("DIAG\n");
            if(head_x == 0 || head_y == 0) goto exit_point;
            rec_X[head_x--] = (char) X[prev_x];
            rec_Y[head_y--] = (char) Y[prev_y];
            ba->length++;
            
        }else if((prev_x - curr_x) > (prev_y - curr_y)){
            //Gap in X
            //printf("Gap X\n");
            if(head_x == 0 || head_y == 0) goto exit_point;
            if(bc->c.xpos != prev_x && bc->c.ypos != prev_y){
                rec_Y[head_y--] = Y[prev_y];
                rec_X[head_x--] = X[prev_x];
            }else{
                rec_Y[head_y--] = '-';
                rec_X[head_x--] = X[prev_x];
            }
            ba->length++;
            
            for(k=prev_x-1;k>curr_x;k--){
                if(head_x == 0 || head_y == 0) goto exit_point;
                #ifdef VERBOSE 
                if(head_x == 0 || head_y == 0){
                    printf("%"PRIu64" %"PRIu64" and prevs are %"PRIu64" %"PRIu64"\n", head_x, head_y, prev_x, prev_y);
                    printf("origin is %"PRIu64", %"PRIu64"\n", bc->c.xpos, bc->c.ypos);
                    uint64_t z;
                    for(z=head_x;z<limit_x;z++){
                        fprintf(stdout, "%c", (char) rec_X[z]);
                    }
                    printf("\n");
                    for(z=head_y;z<limit_y;z++){
                        fprintf(stdout, "%c", (char) rec_Y[z]);
                    }
                    getchar();
                }
                #endif
                rec_Y[head_y--] = '-';
                rec_X[head_x--] = (char) X[k];
                ba->length++;
                ba->egaps++;
            }
            ba->igaps += 1;
            ba->egaps--;
        }else{
            //Gap in Y
            //printf("GAP Y\n");
            //10, 0, 401, 281
            if(head_x == 0 || head_y == 0) goto exit_point;
            if(bc->c.xpos != prev_x && bc->c.ypos != prev_y){
                rec_Y[head_y--] = Y[prev_y];
                rec_X[head_x--] = X[prev_x];
            }else{
                rec_Y[head_y--] = Y[prev_y];
                rec_X[head_x--] = '-';
            }
            ba->length++;

            for(k=prev_y-1;k>curr_y;k--){
                if(head_x == 0 || head_y == 0) goto exit_point;
                #ifdef VERBOSE 
                if(head_x == 0 || head_y == 0){
                    printf("%"PRIu64" %"PRIu64" and prevs are %"PRIu64" %"PRIu64"\n", head_x, head_y, prev_x, prev_y);
                    printf("origin is %"PRIu64", %"PRIu64"\n", bc->c.xpos, bc->c.ypos);
                    uint64_t z;
                    for(z=head_x;z<limit_x;z++){
                        fprintf(stdout, "%c", (char) rec_X[z]);
                    }
                    printf("\n");
                    for(z=head_y;z<limit_y;z++){
                        fprintf(stdout, "%c", (char) rec_Y[z]);
                    }
                    getchar();
                }
                #endif
                rec_X[head_x--] = '-';
                rec_Y[head_y--] = (char) Y[k];
                ba->length++;
                ba->egaps++;
            }
            
            ba->igaps += 1;
            ba->egaps--;
        }
        
    }
    
    if(curr_x == 0 && curr_y == 0 && (curr_x == (prev_x - 1)) && (curr_y == (prev_y -1))){
        rec_X[head_x--] = (char) X[curr_x];
        rec_Y[head_y--] = (char) Y[curr_y];
        ba->length++;
    }
    
    exit_point:

    //printf("curr: %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    //printf("Heads: %"PRIu64", %"PRIu64"\n", head_x, head_y);
    if(show == 1)fprintf(stdout, "%"PRIu64", %"PRIu64"\n", head_x, head_y);
    uint64_t huecos_x = 0, huecos_y = 0;
    k=(int64_t)curr_x-1;
    while(k>=0){ if(head_x == 0) break; rec_X[head_x--] = '-'; huecos_x++;  k--; }
    k=(int64_t)curr_y-1;
    while(k>=0){ if(head_y == 0) break; rec_Y[head_y--] = '-'; huecos_y++; k--; }
    
    if(show == 1)fprintf(stdout, "%"PRIu64", %"PRIu64"\n", head_x, head_y);

    if(huecos_x >= huecos_y){
        while(huecos_x > 0) { if(head_y == 0) break; rec_Y[head_y--] = ' '; huecos_x--;}
    }else{
        while(huecos_y > 0) { if(head_x == 0) break; rec_X[head_x--] = ' '; huecos_y--;}
    }

    if(show == 1){
        fprintf(stdout, "%"PRIu64", %"PRIu64"\n", head_x, head_y);
        fprintf(stdout, "%"PRIu64", %"PRIu64"\n", 2*Xend, 2*Yend);
        uint64_t k;
        for(k=head_x;k<limit_x;k++){
            fprintf(stdout, "%c", (char) rec_X[k]);
        }
        printf("\n");
        for(k=head_y;k<limit_y;k++){
            fprintf(stdout, "%c", (char) rec_Y[k]);
        }
        printf("\n");
        getchar();
    }

    *ret_head_x = head_x;
    *ret_head_y = head_y;
    #ifdef VERBOSE
    printf("hx hy: %"PRIu64", %"PRIu64"\n", head_x, head_y);
    #endif
}
