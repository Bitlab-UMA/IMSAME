#define QF_LAMBDA 0.275
#define QF_KARLIN 0.333

typedef struct container{
    llpos * table[4][4][4][4][4][4][4][4][4][4][4][4];
} Container;



typedef struct {
    SeqInfo * database; //Database sequence and lengths
    SeqInfo * query;    //Query sequence and lengths
    uint64_t from;      //Starting READ to compute alignments from
    uint64_t to;        //End READ to compute alignments from
    Container * container; //Container to hold the multidimensional array
    uint64_t accepted_query_reads; //Number of reads that have a fragment with evalue less than specified
    long double min_e_value;    //Minimum evalue to accept read
    long double min_coverage;    //Minimum coverage percentage to accept read
    long double min_identity;    //Minimum identity percentage to accept read
    FILE * out; //File to write alignments out
    int igap;
    int egap;
    struct positioned_cell * mc;
    struct cell ** table;
    char * reconstruct_X;
    char * reconstruct_Y;
    char * writing_buffer_alignment;
    unsigned char * my_x;
    unsigned char * my_y;
} HashTableArgs;




/**
 * Initialize the memory pool to later retrieve individual memory addresses for llpos
 * 
 */
void init_mem_pool_llpos(Mempool_l * mp);

/**
 * Get a new memory address from the pool mp for a type llpos
 * 
 */
llpos * getNewLocationllpos(Mempool_l * mp, uint64_t * n_pools_used);


/*
    Compute alignments by thread given a hash table argument
*/
void * computeAlignmentsByThread(void * a);


/*
    Performs NW and backtracking to recover alignment
*/
void build_alignment(char * reconstruct_X, char * reconstruct_Y, uint64_t curr_db_seq, uint64_t curr_read, HashTableArgs * hta, unsigned char * my_x, unsigned char * my_y, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, BasicAlignment * ba, uint64_t xlen, uint64_t ylen, uint64_t * cell_path_y);

/*
    Compute the alignment and evalue of a given hit
    The positions pos_database and pos_query refer to the last match in the hit
*/
void alignmentFromQuickHits(SeqInfo * database, SeqInfo * query, uint64_t pos_database, uint64_t pos_query, uint64_t curr_read, uint64_t curr_db_seq, Quickfrag * qf);

/*
    Computes the cell path for the y points given incremental x
    Only add +- window size to each to know which path to go through
*/
uint64_t * calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3);
/*
    Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/
struct positioned_cell NW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show, uint64_t * cell_path_y);

/*
    Computes the alignment given a NW table
*/
void backtrackingNW(unsigned char * X, uint64_t Xstart, uint64_t Xend, unsigned char * Y, uint64_t Ystart, uint64_t Yend, struct cell ** table, char * rec_X, char * rec_Y, struct positioned_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, BasicAlignment * ba);

