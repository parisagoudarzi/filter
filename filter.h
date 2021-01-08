#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "hash_function.h"
struct qf_iterator {
	uint64_t qfi_index;
	uint64_t qfi_quotient;
	uint64_t qfi_visited;
};
struct hash_entry {
    uint32_t key;
    uint32_t status; /* AVAILABLE, OCCUPIED */
    uint32_t value;
};

struct hash_table {
    struct hash_entry **buckets;
    struct hash_entry *items;

    uint32_t bucket_num;
    uint32_t table_size;
    uint32_t num_items;
};
enum { AVAILABLE, OCCUPIED, };
/*
 * Initializes a quotient filter with capacity 2^q.
 * Increasing r improves the filter's accuracy but uses more space.
 * 
 * Returns false if q == 0, r == 0, q+r > 64, or on ENOMEM.
 */
bool qf_init(struct quotient_filter *qf, uint32_t q, uint32_t r);
/*
 * Inserts a hash into the QF.
 * Only the lowest q+r bits are actually inserted into the QF table.
 *
 * Returns false if the QF is full.
 */
bool qf_insert(struct quotient_filter *qf, uint64_t hash);

/*
 * Returns true if the QF may contain the hash.
 *
 * Returns false otherwise.
 * Returns true if the QF may contain the hash. Returns false otherwise.
 */
bool qf_may_contain(struct quotient_filter *qf, uint64_t hash);

/*
 * Removes a hash from the QF.
 *
 * Caution: If you plan on using this function, make sure that your hash
 * function emits no more than q+r bits. Consider the following scenario;
 *
 *	insert(qf, A:X)   # X is in the lowest q+r bits.
 *	insert(qf, B:X)   # This is a no-op, since X is already in the table.
 *	remove(qf, A:X)   # X is removed from the table.
 *
 * Now, may-contain(qf, B:X) == false, which is a ruinous false negative.
 *
 * Returns false if the hash uses more than q+r bits.
 */
bool qf_remove(struct quotient_filter *qf, uint64_t hash);
/*
 * Initializes qfout and copies over all elements from qf1 and qf2.
 * Caution: qfout holds twice as many entries as either qf1 or qf2.
 *
 * Returns false on ENOMEM.
 */
bool qf_merge(struct quotient_filter *qf1, struct quotient_filter *qf2,
	struct quotient_filter *qfout);

/*
 * Resets the QF table.
 *
 * This function does not deallocate any memory.
 * Resets the QF table. This function does not deallocate any memory.
 */
void qf_clear(struct quotient_filter *qf);

/*
 * Finds the size (in bytes) of a QF table.
 *
 * Caution: sizeof(struct quotient_filter) is not included.
 */
size_t qf_table_size(uint32_t q, uint32_t r);
/*
 * Deallocates the QF table.
 */
void qf_destroy(struct quotient_filter *qf);
/*
 * Initialize an iterator for the QF.
 */
void qfi_start(struct quotient_filter *qf, struct qf_iterator *i);
/*
 * Returns true if there are no elements left to visit.
 */
bool qfi_done(struct quotient_filter *qf, struct qf_iterator *i);
/*
 * Returns the next (q+r)-bit fingerprint in the QF.
 *
 * Caution: Do not call this routine if qfi_done() == true.
 */
uint64_t qfi_next(struct quotient_filter *qf, struct qf_iterator *i);
static inline int is_pow_of_2(uint32_t x)
{
    return !(x & (x-1));
}

static inline uint32_t upper_pow_of_2(uint32_t x)
{
    if (is_pow_of_2(x))
        return x;

    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}

/*
int cuckoo_hash_put(struct hash_table *table, uint32_t key, uint32_t value);
void cuckoo_hash_remove(struct hash_table *table, uint32_t key);
*/
void cuckoo_insert(uint32_t key, uint32_t value);
void cuckoo_remove(uint32_t key);
int cuckoo_lookup(uint32_t key);

double load_factor();
void show_hash_table();
int cuckoo_hashing_init(size_t size);

#endif

