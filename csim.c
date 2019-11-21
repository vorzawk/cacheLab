/* csim.c - A cache simulator which outputs the hits, misses and evictions for a given sequence 
 * of memory references.
 * Required inputs : number of bits in the set index(s), associativity(E), number of bits in
 * the offset(b) and a trace file containing memory accesses */

#include "cachelab.h"
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>

typedef struct {
    long long tag;
    short valid;
    short lru_cntr;
} cache_line;

cache_line *cache_lookup(cache_line set[], int assoc, long long tag);
void update_lru_cntr(cache_line set[], int assoc, short lru_cntr_accessed);
void usage(char *argv[]);

int main(int argc, char *argv[]) {
    extern char* optarg;
    int sflag = 0, Eflag = 0, bflag = 0, tflag = 0, err_flag = 0;
    char *sname, *Ename, *trace_file, *bname;
    char c;

/* Use getopt to read the commandline arguments */
    while((c = getopt(argc, argv, "s:E:b:t:")) != -1) {
        switch(c) {
            case 's':
                sflag = 1;
                sname = optarg;
                break;
            case 'E':
                Eflag = 1;
                Ename = optarg;
                break;
            case 'b':
                bflag = 1;
                bname = optarg;
                break;
            case 't':
                tflag = 1;
                trace_file = optarg;
                break;
            case '?':
                err_flag = 1;
                break;
          }
    }

    if ((sflag == 0) || (Eflag == 0) || (bflag == 0) || (tflag == 0)) {
        fprintf(stderr, "required parameter missing, check usage\n");
        usage(argv);
        return -1;
    }
    
    if (err_flag) {
        fprintf(stderr, "unknown parameter encountered, check usage\n");
        usage(argv);
        return -2;
    }

    int s, b, num_sets, assoc;
    s = atoi(sname);
    b = atoi(bname);
    num_sets = (1 << s);
    assoc = atoi(Ename);
    const int INDEX_MASK = (1 << s) - 1;

    cache_line *cache[num_sets];

    for (int i = 0; i < num_sets; i++) {
        cache[i] = (cache_line *)malloc(assoc*sizeof(cache_line));
    }

    for (int i = 0; i < num_sets; i++) {
        for (int j = 0; j < assoc; j++) {
            cache[i][j].valid = 0;
            /* Initializing with j ensures that all the lru_cntr values are distinct as is the case with normal
             * operation */
            cache[i][j].lru_cntr = j; 
        }
    }

    FILE *tracefp;
    tracefp = fopen(trace_file, "r");
    char access_type;
    long long address; // address specifies a 64-bit hex value
    int hits = 0, misses = 0, evictions = 0;

    while(fscanf(tracefp, " %c %llx %*c %*d", &access_type, &address) == 2) {
        if (access_type == 'I') {
            /* Ignore instruction references since we are only interested in data references */
            continue;
        }

        int set_index;
        long long tag;
        cache_line *line_to_replace;
        
        set_index = (address >> b) & INDEX_MASK;
        tag = address >> (s+b);
        printf("%c, %llx, set = %d ", access_type, address, set_index);

        line_to_replace = cache_lookup(cache[set_index], assoc, tag);
        if (line_to_replace == NULL) {
            hits++;
            printf("hit\n");
        } else {
            misses++;
            printf("miss %d ",misses);
            if (line_to_replace -> valid) {
                evictions++;
                printf("eviction");
            } else {
                line_to_replace -> valid = 1;
            }
            printf("\n");

            line_to_replace -> tag = tag;
        }
        
        /* An 'M' or modify type of access reads a value and writes to the same location. So, irrespective of the result
         * of the read, the write is always a hit */
        if (access_type == 'M') {
            hits++;
        }
    }

    printSummary(hits, misses, evictions);
    return 0;
}
cache_line *cache_lookup(cache_line set[], int assoc, long long tag) {
/* cache_lookup searches all the lines in the set for a matching tag and returns the address where the incoming
 * block is to placed in case of a cache miss or a nullptr if it is a hit */

    /* Check if the required data is already present in the cache */
    for (int i = 0; i < assoc; i++) {
        if (set[i].valid && (set[i].tag == tag)) {
            update_lru_cntr(set, assoc, set[i].lru_cntr); 
            return NULL;
        }
    }

    /* If the data is not present, we need to find a slot for the incoming data */
    short max_lru_cntr = -1;
    cache_line *lru_line;
    for (int i = 0; i < assoc; i++) {
        if (set[i].valid == 0) {
            /* If an empty slot is found, no line needs to be evicted */
            update_lru_cntr(set, assoc, set[i].lru_cntr); 
            return &set[i];
        } else if (set[i].lru_cntr > max_lru_cntr) {
            /* If the set is full, evict the least recently used line */
            max_lru_cntr = set[i].lru_cntr;
            lru_line = &set[i];
        }
    }

    /* The lru line must necessarily have a lru_cntr value of assoc-1 */
    assert(max_lru_cntr == assoc-1);
    update_lru_cntr(set, assoc, assoc-1);
    return lru_line;
}

void update_lru_cntr(cache_line set[], int assoc, short lru_cntr_accessed) {
    /* All the lines in the set with lru_cntr values less than that of the accessed block must be incremented and the
     * accessed block must have a lru_cntr value of 0 */
    for (int i = 0; i < assoc; i++) {
        if (set[i].lru_cntr < lru_cntr_accessed) {
            set[i].lru_cntr++;
        } else if (set[i].lru_cntr == lru_cntr_accessed) {
            set[i].lru_cntr = 0;
        }
    }
}

void usage(char *argv[]) {
    printf("%s -s <num> -E <num> -b <num> -t <file>\n", argv[0]);
    printf("\nOptions:\n");
    printf("  -s <num>   Number of set index bits\n");
    printf("  -E <num>   Number of lines per set.\n");
    printf("  -b <num>   Number of block offset bits.\n");
    printf("  -t <file>  Trace file.\n");
    printf("\nExample : %s -s 4 -E 1 -b 4 -t traces/yi.trace\n", argv[0]);       
}
