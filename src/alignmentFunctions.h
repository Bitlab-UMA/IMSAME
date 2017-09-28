#define QF_LAMBDA 0.275
#define QF_KARLIN 0.333



typedef struct {
    uint64_t id;        //The thread id
    SeqInfo * database; //Database sequence and lengths
    SeqInfo * query;    //Query sequence and lengths
    uint64_t from;      //Starting READ to compute alignments from
    uint64_t to;        //End READ to compute alignments from
    Container * container_a; //Container to hold the multidimensional array
    Container * container_b; //Container to hold the multidimensional array
    Container * container_c; //Container to hold the multidimensional array
    Container * container_d; //Container to hold the multidimensional array
    uint64_t * contained_reads; //Number of reads contained in each of the tablespace chunks
    uint64_t * base_coordinates; //Coordinates where each thread started to read their DB chunk
    unsigned char * accepted_query_reads; //reads that have a fragment with evalue less than specified
    long double min_e_value;    //Minimum evalue to accept read
    long double min_coverage;    //Minimum coverage percentage to accept read
    long double min_identity;    //Minimum identity percentage to accept read
    long double window;         //Percentage of window that will be explored (+-)
    FILE * out; //File to write alignments out
    int igap;
    int egap;
    uint64_t * hits;        // To work in hits mode only
    struct positioned_cell * mc;
    struct cell ** table;
    char * reconstruct_X;
    char * reconstruct_Y;
    char * writing_buffer_alignment;
    unsigned char * my_x;
    unsigned char * my_y;
    Head * queue_head;  //To tell where the queue starts after modifications
    pthread_mutex_t * lock;
    unsigned char full_comp; // Tells whether read reporting should stop at first match or keep reporting
    unsigned char * markers; // To tell which sequences were already used
    Mempool_l ** mp_pools; // All memory pools
} HashTableArgs;



/*
    Nucleotides matching function
*/
int64_t compare_letters(unsigned char a, unsigned char b);

/**
 * Initialize the memory pool to later retrieve individual memory addresses for llpos
 * 
 */
void init_mem_pool_llpos(Mempool_l * mp);

/**
 * Get a new memory address from the pool mp for a type llpos
 * 
 */
uint64_t getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used);

/*
    Load input database using 4 threads
*/
void * load_input(void * a);
/*
    Compute alignments by thread given a hash table argument
*/
void * computeAlignmentsByThread(void * a);


/*
    Performs NW and backtracking to recover alignment
*/
void build_alignment(char * reconstruct_X, char * reconstruct_Y, uint64_t curr_db_seq, uint64_t curr_read, HashTableArgs * hta, unsigned char * my_x, unsigned char * my_y, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, BasicAlignment * ba, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window);

/*
    Compute the alignment and evalue of a given hit
    The positions pos_database and pos_query refer to the last match in the hit
*/
void alignmentFromQuickHits(SeqInfo * database, SeqInfo * query, uint64_t pos_database, uint64_t pos_query, uint64_t curr_read, uint64_t curr_db_seq, Quickfrag * qf, uint64_t offset_db_reads, uint64_t offset_db_coordinates);

/*
    Computes the cell path for the y points given incremental x
    Only add +- window size to each to know which path to go through
*/
void calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3, int64_t * cell_path_y);
/*
    Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/
struct best_cell NW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show, int64_t * cell_path_y, long double * window, uint64_t * curr_window_size);

/*
    Computes the alignment given a NW table
*/
void backtrackingNW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, struct cell ** table, char * rec_X, char * rec_Y, struct best_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, BasicAlignment * ba, int64_t * cell_path_y, uint64_t window_size);

