#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))


inline int64_t compare_letters(unsigned char a, unsigned char b){
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

    

    //To keep track of which reads are we reading
    uint64_t curr_read, curr_db_seq, xlen, ylen;
    uint64_t crrSeqL, pos_of_hit;

    //Reading from buffer
    char c;
    unsigned char curr_kmer[FIXED_K], b_aux[FIXED_K];
    llpos * aux;

    // For NW-alignment
    int NWaligned;
    uint64_t n_hits, alignments_tried;

    BasicAlignment ba; //The resulting alignment from the NW
    uint64_t curr_pos; //Reading-head position
    uint64_t up_to;

    //Get next operation in queue
    while(NULL != ( my_current_task = get_task_from_queue(hta->queue_head, hta->lock))){
        //Initialize all variables
        qf.x_start = qf.y_start = qf.t_len = 0;
        qf.e_value = LDBL_MAX;

        //Starting from
        curr_read = my_current_task->r1;
        crrSeqL = 0; pos_of_hit = 0;

        curr_kmer[0] = '\0'; b_aux[0] = '\0';
        aux = NULL;

        NWaligned = 0;
        n_hits = 0;
        alignments_tried = 0;

        ba.identities = 0; ba.length = 0; ba.igaps = 0xFFFFFFFFFFFFFFFF; ba.egaps = 0xFFFFFFFFFFFFFFFF;


        //Set current header position at the position of the read start (the ">")
        curr_pos = hta->query->start_pos[curr_read]; //Skip the ">"
        c = (char) hta->query->sequences[curr_pos];


        while(curr_read < my_current_task->r2 && curr_pos < hta->query->total_len){

            if(curr_read < hta->query->n_seqs - 1) up_to = hta->query->start_pos[curr_read+1]-1; else up_to = hta->query->total_len;
            //printf("Currrpos: %"PRIu64" up to: %"PRIu64" on read: %"PRIu64"\n", curr_pos, up_to, curr_read);

            if (curr_pos == up_to) { // Comment, empty or quality (+) line
                crrSeqL = 0; // Reset buffered sequence length
                #ifdef VERBOSE
                printf("Read: %"PRIu64" yielded (%d)\n", curr_read, NWaligned);
                #endif
                NWaligned = 0;
                //fprintf(stdout, "Seq %"PRIu64" has %"PRIu64" hits and tried to align %"PRIu64" times\n", curr_read, n_hits, alignments_tried);
                //fflush(stdout);
                n_hits = 0;
                alignments_tried = 0;
                curr_read++;
                continue;
            }

            if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
                curr_kmer[crrSeqL] = (unsigned char) c;
                crrSeqL++;
            }else{
                crrSeqL = 0;
            }

            if (crrSeqL >= FIXED_K) { // Full well formed sequence

                //fprintf(stdout, "%s\n", curr_kmer);
                //fflush(stdout);
                aux = hta->container->table[char_converter[curr_kmer[0]]][char_converter[curr_kmer[1]]][char_converter[curr_kmer[2]]]
                        [char_converter[curr_kmer[3]]][char_converter[curr_kmer[4]]][char_converter[curr_kmer[5]]]
                        [char_converter[curr_kmer[6]]][char_converter[curr_kmer[7]]][char_converter[curr_kmer[8]]]
                        [char_converter[curr_kmer[9]]][char_converter[curr_kmer[10]]][char_converter[curr_kmer[11]]];

                //While there are hits
                //fprintf(stdout, "%p\n", aux);
                //fflush(stdout);
                while(aux != NULL && NWaligned == 0){
                    n_hits++;
                    //fprintf(stdout, "%p\n", aux);
                    //fflush(stdout);
                    curr_db_seq = aux->s_id;
                    pos_of_hit = aux->pos;
                    #ifdef VERBOSE 
                    fprintf(stdout, "Launching %"PRIu64" @ %"PRIu64", %"PRIu64"\n", curr_read, curr_pos+1, pos_of_hit);
                    #endif 
                    alignmentFromQuickHits(hta->database, hta->query, pos_of_hit, curr_pos+1, curr_read, curr_db_seq, &qf);

                    #ifdef VERBOSE 
                    printf("curr evalue: %Le %"PRIu64"\n", qf.e_value, qf.t_len);
                    #endif
                    //getchar();
                    

                    //If e-value of current frag is good, then we compute a good gapped alignment
                    if(qf.e_value < hta->min_e_value){
                        alignments_tried++;
                        ba.identities = ba.length = ba.igaps = ba.egaps = 0;
                        //Compute lengths of reads
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
                        //Perform alignment plus backtracking
                        //void build_alignment(char * reconstruct_X, char * reconstruct_Y, uint64_t curr_db_seq, uint64_t curr_read, HashTableArgs * hta, char * my_x, char * my_y, struct cell ** table, struct cell * mc, char * writing_buffer_alignment, BasicAlignment * ba, uint64_t xlen, uint64_t ylen)
                        if(xlen > MAX_READ_SIZE || ylen > MAX_READ_SIZE) terror("Read size reached for gapped alignment.");
                        //fprintf(stdout, "R0 %"PRIu64", %"PRIu64"\n", curr_db_seq, curr_read);
                        
                        #ifdef VERBOSE 
                        fprintf(stdout, "qfxs %"PRIu64", dbs %"PRIu64", qfys %"PRIu64" qys %"PRIu64"\n", qf.x_start, hta->database->start_pos[curr_db_seq], qf.y_start, hta->query->start_pos[curr_read]);
                        #endif

                        p1.x = qf.x_start - hta->database->start_pos[curr_db_seq];
                        //p1.y = qf.y_start - hta->query->start_pos[curr_read];
                        p1.y = qf.y_start - (hta->query->start_pos[curr_read] -1);
                        p2.x = p1.x + qf.t_len;
                        p2.y = p1.y + qf.t_len;
                        p3.x = xlen;
                        p3.y = ylen;
                        calculate_y_cell_path(p0, p1, p2, p3, cell_path_y);
                        build_alignment(hta->reconstruct_X, hta->reconstruct_Y, curr_db_seq, curr_read, hta, hta->my_x, hta->my_y, hta->table, hta->mc, hta->writing_buffer_alignment, &ba, xlen, ylen, cell_path_y, &hta->window);
                        #ifdef VERBOSE
                        printf("len 1 %"PRIu64", len 2 %"PRIu64"\n", ba.length, ylen);
                        printf("ident %"PRIu64"\n", ba.identities);
                        #endif

                        //If is good
                        if(((long double)ba.length/xlen) >= hta->min_coverage && ((long double)ba.identities/ba.length) >=  hta->min_identity){
                            hta->accepted_query_reads++;   
                            if(hta->out != NULL){
                                //printf("Last was: (%"PRIu64", %"PRIu64")\n", curr_read, curr_db_seq);
                                fprintf(hta->out, "(%"PRIu64", %"PRIu64") : %d%% %d%% %"PRIu64"\n $$$$$$$ \n", curr_read, curr_db_seq, MIN(100,(int)(100*ba.identities/ba.length)), MIN(100,(int)(100*ba.length/xlen)), xlen);
                                fprintf(hta->out, "%s", hta->writing_buffer_alignment);
                                //fprintf(stdout, "(%"PRIu64", %"PRIu64") : %d%% %d%% %"PRIu64"\n $$$$$$$ \n", curr_read, curr_db_seq, MIN(100,(int)(100*ba.identities/ba.length)), MIN(100,(int)(100*ba.length/ylen)), ylen);
                                //fprintf(stdout, "%s", hta->writing_buffer_alignment);
                            }
                            NWaligned = 1;
                        }
                    
                    }

                    //strncpy(get_from_db, &hta->database->sequences[qf.x_start], qf.t_len);
                    //strncpy(get_from_query, &hta->query->sequences[qf.y_start], qf.t_len);
                    //fprintf(hta->out, "%s\n%s\n%Le\t%d\n-------------------\n", get_from_db, get_from_query, qf.e_value, (int)(100*qf.coverage));
                    //fprintf(hta->out, "%"PRIu64", %"PRIu64", %"PRIu64"\n", qf.x_start, qf.y_start, qf.t_len);

                    //printf("Hit comes from %"PRIu64", %"PRIu64"\n", pos_of_hit, curr_pos);
                    aux = aux->next;
                    //fprintf(stdout, "%p\n", aux);
                    //fflush(stdout);
                }
                //printf("SWITCHED\n");

                if(NWaligned == 1){
                    if(curr_read < hta->query->n_seqs) curr_pos = hta->query->start_pos[curr_read+1]-2;
                }else{
                    memcpy(b_aux, curr_kmer, FIXED_K);
                    memcpy(curr_kmer, &b_aux[1], FIXED_K-1);
                    crrSeqL -= 1;
                }
            }
        
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

void alignmentFromQuickHits(SeqInfo * database, SeqInfo * query, uint64_t pos_database, uint64_t pos_query, uint64_t curr_read, uint64_t curr_db_seq, Quickfrag * qf){

    int64_t read_x_start, read_x_end, read_y_start, read_y_end;

    if(curr_db_seq == database->n_seqs-1){
        read_x_start = database->start_pos[curr_db_seq];
	    read_x_end = database->total_len;
    }else{     
        read_x_start = database->start_pos[curr_db_seq];
	    read_x_end = database->start_pos[curr_db_seq+1] - 1;
    }

    if(curr_read == query->n_seqs-1){
        read_y_start = query->start_pos[curr_read];
        read_y_end = query->total_len;
    }else{
        read_y_start = query->start_pos[curr_read];
        read_y_end = query->start_pos[curr_read+1] - 1;
    }



    int64_t curr_pos_db = (int64_t) pos_database;
    int64_t curr_pos_qy = (int64_t) pos_query;
    int64_t final_end_x = pos_database - 1, final_start_x = final_end_x - FIXED_K + 1, final_start_y = pos_query - FIXED_K;
    int64_t score_right = FIXED_K * POINT;
    int64_t score_left = score_right;
    int64_t high_left = score_left, high_right = score_right;
    qf->t_len = FIXED_K;
    uint64_t idents = FIXED_K;

    /*
    char le_hit[1000];
    memcpy(le_hit, &database->sequences[final_start_x], FIXED_K);
    
	fprintf(stdout, "HIT: %s\n", le_hit);
	fflush(stdout);
    */

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

    keep_going = 1;
    curr_pos_db = pos_database - FIXED_K - 1;
    curr_pos_qy = pos_query - FIXED_K - 1;

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

    
	fprintf(stdout, "%s\n%s\n------\n", s1, s2);
	fflush(stdout);
    getchar();
    */

    long double rawscore = (idents*POINT) - (qf->t_len - idents)*(POINT);

    long double t_len;
    if(curr_read == query->n_seqs-1){
        t_len = (long double) query->total_len - query->start_pos[curr_read];
    }else{
        t_len = (long double) query->start_pos[curr_read+1] - query->start_pos[curr_read];
    }

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

    #ifdef VERBOSE
    for(i=0;i<p3.x;i++){
        printf("%"PRIu64" -> ", y_points[i]);
        if(i % 50 == 0) getchar();
    }
    #endif
    

}

struct best_cell NW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show, int64_t * cell_path_y, long double * window, uint64_t * current_window_size){
    

    uint64_t i, j, j_prime;
    int64_t scoreDiagonal,scoreLeft,scoreRight,score, delta_dif_1_0, delta_dif_2_1, limit_left, limit_right, j_right_prime = 1, j_left_prime, j_diag_prime;

    struct best_cell bc;
    bc.c.score = INT64_MIN;
    
    //The window size will be a +-15% of the square root of the product of lengths
    int64_t window_size = (uint64_t) (*window * sqrtl((long double) Xend * (long double) Yend));
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

    for(i=1;i<Yend;i++){
        //table[0][i].score = (X[0] == Y[i]) ? POINT : -POINT;
        if(i < cell_path_y[0] + window_size) table[0][i].score = compare_letters(X[0], Y[i]) + iGap + (i-1)*eGap;
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
        //printf("I will go from %"PRIu64" to %"PRIu64"\n", (uint64_t) MAX(1,(cell_path_y[i] - window_size)), (uint64_t) MIN(Yend,(cell_path_y[i] + window_size)));
        //getchar();

        #ifdef VERBOSE
        int64_t r;
        for(r=0;r<MAX(0,(cell_path_y[i] - window_size)); r++){
            printf("  ");
        }
        #endif

        

        

        for(j=MAX(1,(cell_path_y[i] - window_size));j<MIN(Yend,(cell_path_y[i] + window_size));j++){
            //printf("Doing on : (%"PRIu64",%"PRIu64" and jprime=%"PRIu64"\n", i,j,j_prime);
            //Check if max in row has changed
            //if(j > MAX(1, cell_path_y[i-1] - window_size +1) && mf.score <= table[i][j-2].score){

            //Calculate the real j position in the windowed table
            
            j_left_prime = ((int64_t)j_prime - (2 - delta_dif_1_0));
            //j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            if(i > 1){
                j_right_prime = ((int64_t)j_prime - (1 - (delta_dif_1_0 + delta_dif_2_1)));
            }

            if(j > MAX(1, cell_path_y[i-1] - window_size +1) && j_left_prime >= limit_left && j_left_prime < limit_right && table[i-1][j_left_prime].score >= mf.score){
                //mf.score = table[i-1][j-2].score;
                mf.score = table[i-1][j_left_prime].score;
                mf.xpos = i-1;
                mf.ypos = j-2;
            }
            //printf("RowMax: %"PRId64"@(%"PRIu64", %"PRIu64")\t", mf.score, mf.xpos, mf.ypos);
            
            //score = (X[i] == Y[j]) ? POINT : -POINT;
            score = compare_letters(X[i], Y[j]);

            //Precondition: Upper row needs to reach up to diagonal
            //if((cell_path_y[i-1]+window_size) >= j-1){
            if(j-1 >= MAX(0, (cell_path_y[i-1]-window_size)) && (cell_path_y[i-1]+window_size) >= j-1 && j_diag_prime >= limit_left && j_diag_prime < limit_right){
                //scoreDiagonal = table[i-1][j-1].score + score;
                //printf("prevdiag: %"PRId64"\n", table[i-1][j_diag_prime].score);
                scoreDiagonal = table[i-1][j_diag_prime].score + score;
                //printf("j_diag: %"PRId64":", j_diag_prime);
            }else{
                scoreDiagonal = INT64_MIN;
            }
            
            if(i>=1 && j>1){
                scoreLeft = mf.score + iGap + (j - (mf.ypos+2))*eGap + score;
                
            }else{
                scoreLeft = INT64_MIN;
            }

            if(j>=1 && i>1){
                scoreRight = mc[j-1].score + iGap + (i - (mc[j-1].xpos+2))*eGap + score;
                //if(scoreRight == -12) printf("MC: %"PRId64", From: %"PRIu64", %"PRIu64"->", mc[j-1].score, mc[j-1].xpos, mc[j-1].ypos);
            }else{
                scoreRight = INT64_MIN;
            }
            
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
            
            if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight){
                //Diagonal
                table[i][j_prime].score = scoreDiagonal;
                table[i][j_prime].xfrom = i-1;
                table[i][j_prime].yfrom = j-1;
                
                                
            }else if(scoreRight > scoreLeft){
                table[i][j_prime].score = scoreRight;
                table[i][j_prime].xfrom = mc[j-1].xpos;
                table[i][j_prime].yfrom = mc[j-1].ypos;
                
            }else{
                table[i][j_prime].score = scoreLeft;
                table[i][j_prime].xfrom = mf.xpos;
                table[i][j_prime].yfrom = mf.ypos;
            }
        
        
            //check if column max has changed
            //New condition: check if you filled i-2, j-1
            
            if(i > 1 && j >= 1 && j_right_prime >= limit_left && j_right_prime < limit_right && table[i-2][j_right_prime].score >= mc[j-1].score){
                //mc[j-1].score = table[i-2][j-(1+j_prime)].score;
                //Should be the j_prime we had at cell_path_y
                
                mc[j-1].score = table[i-2][j_right_prime].score;
                mc[j-1].xpos = i-2;
                mc[j-1].ypos = j-1;
            }
            if(i == Xend-1 || j == Yend-1){

                if(i == Xend-1 && j != Yend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Yend - j)*eGap;
            	}else if(j == Yend-1 && i != Xend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Xend - i)*eGap;
            	}
                //Check for best cell
                if(table[i][j_prime].score >= bc.c.score){ bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime; }
                //bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime;
            }
            #ifdef VERBOSE
            //printf("Put score: %"PRId64"\n\n", table[i][j_prime].score);
            printf("(%"PRId64")%02"PRId64" ", j_diag_prime, table[i][j_prime].score); //printf("->(%"PRIu64", %"PRIu64")", i, j); printf("[%c %c]", X[i], Y[j]);
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
    prev_x = curr_x;
    prev_y = curr_y;
    int show = 0;
   
    for(k=Xend-1; k>curr_x; k--) rec_X[head_x--] = '-';
    for(k=Yend-1; k>curr_y; k--) rec_Y[head_y--] = '-';

    j_prime = bc->j_prime;
    unsigned char first_track = 1;
    
    while(curr_x > 0 && curr_y > 0){

        
        if(first_track == 0){
            delta_diff = MAX(1, cell_path_y[prev_x] - (int64_t) window_size) - MAX(1, cell_path_y[curr_x] - (int64_t)window_size); //j-1
            j_prime = MAX(0, j_prime - (int64_t)(prev_y - curr_y) + (int64_t) delta_diff);
            //j_prime = j_prime - (int64_t)(prev_y - curr_y) + (int64_t) delta_diff;

            prev_x = curr_x;
            prev_y = curr_y;

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

        curr_x = table[prev_x][j_prime].xfrom;
        curr_y = table[prev_x][j_prime].yfrom;
        first_track = 0;
        

        

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
