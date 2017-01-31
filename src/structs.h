#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>

//Structs required for the dotplot workflow
#define MAXLID 200
//#define READBUF 2000000000 //2 GB
#define READBUF 50000000 //50MB
#define INITSEQS 3000 //Number of starting sequences (in database)
#define POINT 4

#define FIXED_K 12

#define MAXLID 200
#define ALIGN_LEN 60 //For NW alignment
#define MAX_READ_SIZE 3000 //Beware: A NW table is computed of cuadratic size
//#define POOL_SIZE 2500000000 // by 16 bytes it is 40 GB
#define POOL_SIZE 12500000 // 1 GB if 16 bytes
#define MAX_MEM_POOLS 256 


//Struct for linked list of positions
typedef struct linked_list_pos{
    uint64_t pos;
    uint64_t s_id;
    struct linked_list_pos * next;
} llpos;

//Struct for memory pool por lists
typedef struct mempool_l{
    llpos * base;
    uint64_t current;
} Mempool_l;


//Struct for a whole sequence(s) data
typedef struct seqinfo{
    unsigned char * sequences;
    uint64_t * start_pos;
    uint64_t total_len;
    uint64_t n_seqs;
} SeqInfo;

//Struct for the alignment of a quick hit
typedef struct quickfrag{
    uint64_t x_start;
    uint64_t y_start;
    uint64_t t_len;
    long double coverage;
    long double e_value;
} Quickfrag;

typedef struct point{
    uint64_t x;
    uint64_t y;
} Point;



struct cell{
    int64_t score;
    uint64_t xfrom;
    uint64_t yfrom;
};

struct positioned_cell{
    int64_t score;
    uint64_t xpos;
    uint64_t ypos;
};

typedef struct{
    uint64_t identities;
    uint64_t length;
    uint64_t igaps;
    uint64_t egaps;
} BasicAlignment;


#endif
