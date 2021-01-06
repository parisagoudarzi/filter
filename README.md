# filter
cuckoo filter. quotient filter
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
struct qf_iterator {
	uint64_t qfi_index;
	uint64_t qfi_quotient;
	uint64_t qfi_visited;
};
#include "hash_function.h"

//#define CUCKOO_DEBUG
#define ASSOC_WAY       (4)  /* 4-way association */
#define MAX_KICK_OUT    (64) /* maximum number of cuckoo kicks before claiming failure */
#define HASH_FUNC_NUM   16
#define XIP_VXT_TABLE_SIZE  128
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

@@ -56,29 +54,27 @@
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

enum { AVAILABLE, OCCUPIED, };

/* The hash entries, including the key-value,
 * and the status of the buckets
 */
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

#include "cuckoo_hashing.h"
#include <stdlib.h>
#include <string.h>

#include "qf.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define LOW_MASK(n) ((1ULL << (n)) - 1ULL)

bool qf_init(struct quotient_filter *qf, uint32_t q, uint32_t r)
{
	if (q == 0 || r == 0 || q + r > 64) {
		return false;
	}

	qf->qf_qbits = q;
	qf->qf_rbits = r;
	qf->qf_elem_bits = qf->qf_rbits + 3;
	qf->qf_index_mask = LOW_MASK(q);
	qf->qf_rmask = LOW_MASK(r);
	qf->qf_elem_mask = LOW_MASK(qf->qf_elem_bits);
	qf->qf_entries = 0; 
	qf->qf_max_size = 1 << q;
	qf->qf_table = (uint64_t *) calloc(qf_table_size(q, r), 1);
	return qf->qf_table != NULL;
}

/* Return QF[idx] in the lower bits. */
static uint64_t get_elem(struct quotient_filter *qf, uint64_t idx)
{
	uint64_t elt = 0;
	size_t bitpos = qf->qf_elem_bits * idx;
	size_t tabpos = bitpos / 64;
	size_t slotpos = bitpos % 64;
	int spillbits = (slotpos + qf->qf_elem_bits) - 64;
	elt = (qf->qf_table[tabpos] >> slotpos) & qf->qf_elem_mask;
	if (spillbits > 0) {
		++tabpos;
		uint64_t x = qf->qf_table[tabpos] & LOW_MASK(spillbits);
		elt |= x << (qf->qf_elem_bits - spillbits);
	}
	return elt;
}

/* Store the lower bits of elt into QF[idx]. */
static void set_elem(struct quotient_filter *qf, uint64_t idx, uint64_t elt)
{
	size_t bitpos = qf->qf_elem_bits * idx;
	size_t tabpos = bitpos / 64;
	size_t slotpos = bitpos % 64;
	int spillbits = (slotpos + qf->qf_elem_bits) - 64;
	elt &= qf->qf_elem_mask;
	qf->qf_table[tabpos] &= ~(qf->qf_elem_mask << slotpos);
	qf->qf_table[tabpos] |= elt << slotpos;
	if (spillbits > 0) {
		++tabpos;
		qf->qf_table[tabpos] &= ~LOW_MASK(spillbits);
		qf->qf_table[tabpos] |= elt >> (qf->qf_elem_bits - spillbits);
	}
}

static inline uint64_t incr(struct quotient_filter *qf, uint64_t idx)
{
	return (idx + 1) & qf->qf_index_mask;
}

static inline uint64_t decr(struct quotient_filter *qf, uint64_t idx)
{
	return (idx - 1) & qf->qf_index_mask;
}

static inline int is_occupied(uint64_t elt)
{
	return elt & 1;
}

static inline uint64_t set_occupied(uint64_t elt)
{
	return elt | 1;
}

static inline uint64_t clr_occupied(uint64_t elt)
{
	return elt & ~1;
}

static inline int is_continuation(uint64_t elt)
{
	return elt & 2;
}

static inline uint64_t set_continuation(uint64_t elt)
{
	return elt | 2;
}

static inline uint64_t clr_continuation(uint64_t elt)
{
	return elt & ~2;
}

static inline int is_shifted(uint64_t elt)
{
	return elt & 4;
}

static inline uint64_t set_shifted(uint64_t elt)
{
	return elt | 4;
}

static inline uint64_t clr_shifted(uint64_t elt)
{
	return elt & ~4;
}

static inline uint64_t get_remainder(uint64_t elt)
{
	return elt >> 3;
}

static inline bool is_empty_element(uint64_t elt)
{
	return (elt & 7) == 0;
}

static inline bool is_cluster_start(uint64_t elt)
{
	return is_occupied(elt) && !is_continuation(elt) && !is_shifted(elt);
}

static inline bool is_run_start(uint64_t elt)
{
	return !is_continuation(elt) && (is_occupied(elt) || is_shifted(elt));
}

static inline uint64_t hash_to_quotient(struct quotient_filter *qf,
		uint64_t hash)
{
	return (hash >> qf->qf_rbits) & qf->qf_index_mask;
}

static inline uint64_t hash_to_remainder(struct quotient_filter *qf,
		uint64_t hash)
{
	return hash & qf->qf_rmask;
}

/* Find the start index of the run for fq (given that the run exists). */
static uint64_t find_run_index(struct quotient_filter *qf, uint64_t fq)
{
	/* Find the start of the cluster. */
	uint64_t b = fq;
	while (is_shifted(get_elem(qf, b))) {
		b = decr(qf, b);
	}

	/* Find the start of the run for fq. */
	uint64_t s = b;
	while (b != fq) {
		do {
			s = incr(qf, s);
		} while (is_continuation(get_elem(qf, s)));

		do {
			b = incr(qf, b);
		} while (!is_occupied(get_elem(qf, b)));
	}
	return s;
}

/* Insert elt into QF[s], shifting over elements as necessary. */
static void insert_into(struct quotient_filter *qf, uint64_t s, uint64_t elt)
{
	uint64_t prev;
	uint64_t curr = elt;
	bool empty;

	do {
		prev = get_elem(qf, s);
		empty = is_empty_element(prev);
		if (!empty) {
			/* Fix up `is_shifted' and `is_occupied'. */
			prev = set_shifted(prev);
			if (is_occupied(prev)) {
				curr = set_occupied(curr);
				prev = clr_occupied(prev);
			}
		}
		set_elem(qf, s, curr);
		curr = prev;
		s = incr(qf, s);
	} while (!empty);
}

bool qf_insert(struct quotient_filter *qf, uint64_t hash)
{
	if (qf->qf_entries >= qf->qf_max_size) {
		return false;
	}

	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);
	uint64_t entry = (fr << 3) & ~7;

	/* Special-case filling canonical slots to simplify insert_into(). */
	if (is_empty_element(T_fq)) {
		set_elem(qf, fq, set_occupied(entry));
		++qf->qf_entries;
		return true;
	}

	if (!is_occupied(T_fq)) {
		set_elem(qf, fq, set_occupied(T_fq));
	}

	uint64_t start = find_run_index(qf, fq);
	uint64_t s = start;

	if (is_occupied(T_fq)) {
		/* Move the cursor to the insert position in the fq run. */
		do {
			uint64_t rem = get_remainder(get_elem(qf, s));
			if (rem == fr) {
				return true;
			} else if (rem > fr) {
				break;
			}
			s = incr(qf, s);
		} while (is_continuation(get_elem(qf, s)));

		if (s == start) {
			/* The old start-of-run becomes a continuation. */
			uint64_t old_head = get_elem(qf, start);
			set_elem(qf, start, set_continuation(old_head));
		} else {
			/* The new element becomes a continuation. */
			entry = set_continuation(entry);
		}
	}

	/* Set the shifted bit if we can't use the canonical slot. */
	if (s != fq) {
		entry = set_shifted(entry);
	}

	insert_into(qf, s, entry);
	++qf->qf_entries;
	return true;
}

bool qf_may_contain(struct quotient_filter *qf, uint64_t hash)
{
	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);

	/* If this quotient has no run, give up. */
	if (!is_occupied(T_fq)) {
		return false;
	}

	/* Scan the sorted run for the target remainder. */
	uint64_t s = find_run_index(qf, fq);
	do {
		uint64_t rem = get_remainder(get_elem(qf, s));
		if (rem == fr) {
			return true;
		} else if (rem > fr) {
			return false;
		}
		s = incr(qf, s);
	} while (is_continuation(get_elem(qf, s)));
	return false;
}

/* Remove the entry in QF[s] and slide the rest of the cluster forward. */
static void delete_entry(struct quotient_filter *qf, uint64_t s, uint64_t quot)
{
	uint64_t next;
	uint64_t curr = get_elem(qf, s);
	uint64_t sp = incr(qf, s);
	uint64_t orig = s;

	/*
	 * FIXME(vsk): This loop looks ugly. Rewrite.
	 */
	while (true) {
		next = get_elem(qf, sp);
		bool curr_occupied = is_occupied(curr);

		if (is_empty_element(next) || is_cluster_start(next) || sp == orig) {
			set_elem(qf, s, 0);
			return;
		} else {
			/* Fix entries which slide into canonical slots. */
			uint64_t updated_next = next;
			if (is_run_start(next)) {
				do {
					quot = incr(qf, quot);
				} while (!is_occupied(get_elem(qf, quot)));

				if (curr_occupied && quot == s) {
					updated_next = clr_shifted(next);
				}
			}

			set_elem(qf, s, curr_occupied ?
					set_occupied(updated_next) :
					clr_occupied(updated_next));
			s = sp;
			sp = incr(qf, sp);
			curr = next;
		}
	}
}

bool qf_remove(struct quotient_filter *qf, uint64_t hash)
{
	uint64_t highbits = hash >> (qf->qf_qbits + qf->qf_rbits);
	if (highbits) {
		return false;
	}

	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);

	if (!is_occupied(T_fq) || !qf->qf_entries) {
		return true;
	}

	uint64_t start = find_run_index(qf, fq);
	uint64_t s = start;
	uint64_t rem;

	/* Find the offending table index (or give up). */
	do {
		rem = get_remainder(get_elem(qf, s));
		if (rem == fr) {
			break;
		} else if (rem > fr) {
			return true;
		}
		s = incr(qf, s);
	} while (is_continuation(get_elem(qf, s)));
	if (rem != fr) {
		return true;
	}

	uint64_t kill = (s == fq) ? T_fq : get_elem(qf, s);
	bool replace_run_start = is_run_start(kill);

	/* If we're deleting the last entry in a run, clear `is_occupied'. */
	if (is_run_start(kill)) {
		uint64_t next = get_elem(qf, incr(qf, s));
		if (!is_continuation(next)) {
			T_fq = clr_occupied(T_fq);
			set_elem(qf, fq, T_fq);
		}
	}

	delete_entry(qf, s, fq);

	if (replace_run_start) {
		uint64_t next = get_elem(qf, s);
		uint64_t updated_next = next;
		if (is_continuation(next)) {
			/* The new start-of-run is no longer a continuation. */
			updated_next = clr_continuation(next);
		}
		if (s == fq && is_run_start(updated_next)) {
			/* The new start-of-run is in the canonical slot. */
			updated_next = clr_shifted(updated_next);
		}
		if (updated_next != next) {
			set_elem(qf, s, updated_next);
		}
	}

	--qf->qf_entries;
	return true;
}

bool qf_merge(struct quotient_filter *qf1, struct quotient_filter *qf2,
		struct quotient_filter *qfout)
{
	uint32_t q = 1 + MAX(qf1->qf_qbits, qf2->qf_qbits);
	uint32_t r = MAX(qf1->qf_rbits, qf2->qf_rbits);
	if (!qf_init(qfout, q, r)) {
		return false;
	}

	struct qf_iterator qfi;
	qfi_start(qf1, &qfi);
	while (!qfi_done(qf1, &qfi)) {
		qf_insert(qfout, qfi_next(qf1, &qfi));
	}
	qfi_start(qf2, &qfi);
	while (!qfi_done(qf2, &qfi)) {
		qf_insert(qfout, qfi_next(qf2, &qfi));
	}
	return true;
}

void qf_clear(struct quotient_filter *qf)
{
	qf->qf_entries = 0;
	memset(qf->qf_table, 0, qf_table_size(qf->qf_qbits, qf->qf_rbits));
}

size_t qf_table_size(uint32_t q, uint32_t r)
{
	size_t bits = (1 << q) * (r + 3);
	size_t bytes = bits / 8;
	return (bits % 8) ? (bytes + 1) : bytes;
}

void qf_destroy(struct quotient_filter *qf)
{
	free(qf->qf_table);
}

void qfi_start(struct quotient_filter *qf, struct qf_iterator *i)
{
	/* Mark the iterator as done. */
	i->qfi_visited = qf->qf_entries;

	if (qf->qf_entries == 0) {
		return;
	}

	/* Find the start of a cluster. */
	uint64_t start;
	for (start = 0; start < qf->qf_max_size; ++start) {
		if (is_cluster_start(get_elem(qf, start))) {
			break;
		}
	}

	i->qfi_visited = 0;
	i->qfi_index = start;
}

bool qfi_done(struct quotient_filter *qf, struct qf_iterator *i)
{
	return qf->qf_entries == i->qfi_visited;
}

uint64_t qfi_next(struct quotient_filter *qf, struct qf_iterator *i)
{
	while (!qfi_done(qf, i)) {
		uint64_t elt = get_elem(qf, i->qfi_index);

		/* Keep track of the current run. */
		if (is_cluster_start(elt)) {
			i->qfi_quotient = i->qfi_index;
		} else {
			if (is_run_start(elt)) {
				uint64_t quot = i->qfi_quotient;
				do {
					quot = incr(qf, quot);
				} while (!is_occupied(get_elem(qf, quot)));
				i->qfi_quotient = quot;
			}
		}

		i->qfi_index = incr(qf, i->qfi_index);

		if (!is_empty_element(elt)) {
			uint64_t quot = i->qfi_quotient;
			uint64_t rem = get_remainder(elt);
			uint64_t hash = (quot << qf->qf_rbits) | rem;
			++i->qfi_visited;
			return hash;
		}
	}

	abort();
}

static struct hash_table cuckoo_hash_table;

static uint (* hash_func[HASH_FUNC_NUM])(const unsigned char * str, uint len) = 
{BOB1, JSHash, OCaml, OAAT, PJWHash, RSHash,  SDBM, Simple, SML, STL,
    APHash, BKDR, DEKHash, DJBHash, FNV32, Hsieh};

double load_factor() {
    return 1.0 * cuckoo_hash_table.num_items / (cuckoo_hash_table.table_size + 0.0);
}

void show_hash_table()
{
    int i, j;

    printf("Show all entries in the cuckoo hash table [key, value, status]:\n");
    for (i = 0; i < cuckoo_hash_table.bucket_num; i++) {
        printf("bucket[%02x]:", i);

        struct hash_entry *items = cuckoo_hash_table.buckets[i];

        for (j = 0; j < ASSOC_WAY; j++) {
            printf("\t[%04x, %02x, %x]", items[j].key, items[j].value, items[j].status);
        }

        printf("\n");
    }
}

static int cuckoo_hash_collide(struct hash_table *table, uint32_t *hash, uint32_t key, uint32_t value)
{
    int j, k, alt_cnt;
    uint32_t old_hash[2], alternate_hash;
    uint32_t old_val, old_key;
    struct hash_entry *items;

    /* Kick out the old bucket and move it to the alternative bucket. */
    //old_hash[0] = hash[0];
    //old_hash[1] = hash[1];

    items = table->buckets[hash[0]];

    /* Obtain the old (key, value), and calculate the corresponding hashes */
    old_key = items[0].key;
    old_val = items[0].value;
    
    //old_hash[0] = old_key & (table->bucket_num - 1);
    //old_hash[1] = (old_hash[0] + 11) & (table->bucket_num - 1);
    
    old_hash[0] = hash_func[0]((const unsigned char *)&old_key, sizeof(old_key)) % table->bucket_num;
    old_hash[1] = hash_func[1]((const unsigned char *)&old_key, sizeof(old_key)) % table->bucket_num;

    alternate_hash = (old_hash[0] == hash[0] ? old_hash[1] : old_hash[0]);

    /* insert the new (key, value) into the current slot */
    items[0].key = key;
    items[0].value = value;

    k = 0;
    alt_cnt = 0;

KICK_OUT:
    items = table->buckets[alternate_hash];

    for (j = 0; j < ASSOC_WAY; j++) {
        if (items[j].status == AVAILABLE) {
            items[j].status = OCCUPIED;
            items[j].key = old_key;
            items[j].value = old_val;
            break;
        }
    }

    if (j == ASSOC_WAY) {
        if (++alt_cnt > MAX_KICK_OUT) {
            if (k == ASSOC_WAY - 1) {
                /* Hash table is almost full and needs to be resized */
                return 1;
            } else {
                k++;
            }
        }

        uint32_t tmp_key = items[k].key;
        uint32_t tmp_val = items[k].value;
        items[k].key = old_key;
        items[k].value = old_val;

        old_key = tmp_key;
        old_val = tmp_val;
        //old_hash[0] = old_key & (table->bucket_num - 1);
        //old_hash[1] = (old_hash[0] + 11) & (table->bucket_num - 1);
        
        old_hash[0] = hash_func[0]((const unsigned char *)&old_key, sizeof(old_key)) % table->bucket_num;
        old_hash[1] = hash_func[1]((const unsigned char *)&old_key, sizeof(old_key)) % table->bucket_num;


        alternate_hash = (old_hash[0] == alternate_hash ? old_hash[1] : old_hash[0]);

        goto KICK_OUT;
    }

    return 0;
}

static int cuckoo_hash_get(struct hash_table *table, uint32_t key)
{
    int i, j;
    uint32_t hash[2];
    uint32_t val;

    struct hash_entry *items;

    //hash[0] = key & (table->bucket_num - 1);
    //hash[1] = (hash[0] + 11) & (table->bucket_num - 1);

    hash[0] = hash_func[0]((const unsigned char *)&key, sizeof(key)) % table->bucket_num;
    hash[1] = hash_func[1]((const unsigned char *)&key, sizeof(key)) % table->bucket_num;

#ifdef CUCKOO_DEBUG
    printf("get h0:%x h1:%x\n", hash[0], hash[1]);
#endif

    /* Filter the key. */
    items = table->buckets[hash[0]];

    for (i = 0; i < ASSOC_WAY; i++) {
        if (key == items[i].key) {
            if (items[i].status == OCCUPIED) {
                val = items[i].value;
                return val;
            }
        }
    }

    if (i == ASSOC_WAY) {
        items = table->buckets[hash[1]];
        for (j = 0; j < ASSOC_WAY; j++) {
            if (key == items[j].key) {
                if (items[j].status == OCCUPIED) {
                    val = items[j].value;
                    return val;
                }
            }
        }

        if (j == ASSOC_WAY) {
#ifdef CUCKOO_DEBUG
            printf("Key not exists!\n");
#endif
            return AVAILABLE;
        }
    }

    return -1;
}

static int cuckoo_hash_put(struct hash_table *table, uint32_t key, uint32_t value)
{
    int i, j;
    uint32_t hash[2];
    struct hash_entry *items;

    //hash[0] = key & (table->bucket_num - 1);
    //hash[1] = (hash[0] + 11) & (table->bucket_num - 1);

    hash[0] = hash_func[0]((const unsigned char *)&key, sizeof(key)) % table->bucket_num;
    hash[1] = hash_func[1]((const unsigned char *)&key, sizeof(key)) % table->bucket_num;

#ifdef CUCKOO_DEBUG
    printf("put: value:%x h0:%x h1:%x\n", value, hash[0], hash[1]);
#endif

    /* Insert new key into hash buckets. */
    items = table->buckets[hash[0]];
    for (i = 0; i < ASSOC_WAY; i++) {
        if (items[i].status == AVAILABLE) {
            items[i].status = OCCUPIED;
            items[i].key = key;
            items[i].value = value;
            break;
        }
    }

    if (i == ASSOC_WAY) {
        items = table->buckets[hash[1]];
        for (j = 0; j < ASSOC_WAY; j++) {
            if (items[j].status == AVAILABLE) {
                items[j].status = OCCUPIED;
                items[j].key = key;
                items[j].value = value;
                break;
            }
        }

        if (j == ASSOC_WAY) {
            if (cuckoo_hash_collide(table, hash, key, value)) {
#ifdef CUCKOO_DEBUG
                printf("Hash table collision!\n");
#endif
                return -1;
            }
        }
    }

#ifdef CUCKOO_DEBUG
    show_hash_table();
#endif

    return 0;
}

static void cuckoo_hash_remove(struct hash_table *table, uint32_t key)
{
    uint32_t i, j, hash[2];
    struct hash_entry *items;

    //hash[0] = key & (table->bucket_num - 1);
    //hash[1] = (hash[0] + 11) & (table->bucket_num - 1);

    hash[0] = hash_func[0]((const unsigned char *)&key, sizeof(key)) % table->bucket_num;
    hash[1] = hash_func[1]((const unsigned char *)&key, sizeof(key)) % table->bucket_num;

#ifdef CUCKOO_DEBUG
    printf("Remove: h0:%x h1:%x\n", hash[0], hash[1]);
#endif

    items = table->buckets[hash[0]];

    for (i = 0; i < ASSOC_WAY; i++) {
        if (key == items[i].key) {
            items[i].status = AVAILABLE;
            items[i].key = 0;
            items[i].value = 0;
            cuckoo_hash_table.num_items --;
            return;
        }
    }

    if (i == ASSOC_WAY) {
        items = table->buckets[hash[1]];
        for (j = 0; j < ASSOC_WAY; j++) {
            if (key == items[j].key) {
                items[j].status = AVAILABLE;
                items[j].key = 0;
                items[j].value = 0;
                cuckoo_hash_table.num_items --;
                return;
            }
        }

        if (j == ASSOC_WAY) {
#ifdef CUCKOO_DEBUG
            printf("Key not exists!\n");
#endif
        }
    }
}

static void cuckoo_rehash(struct hash_table *table)
{
    int i;
    struct hash_table old_table;

    /* Reallocate hash slots */
    old_table.items = table->items;
    old_table.table_size = table->table_size;

    table->table_size *= 2;
    table->items = calloc(table->table_size, sizeof(struct hash_entry));

    if (table->items == NULL) {
        table->items = old_table.items;
        return;
    }

    /* Reallocate hash buckets associated with slots */
    old_table.buckets = table->buckets;
    old_table.bucket_num = table->bucket_num;

    table->bucket_num *= 2;
    table->buckets = malloc(table->bucket_num * sizeof(struct hash_entry *));

    if (table->buckets == NULL) {
        free(table->items);
        table->items = old_table.items;
        table->buckets = old_table.buckets;
        return;
    }

    for (i = 0; i < table->bucket_num; i++) {
        table->buckets[i] = &table->items[i * ASSOC_WAY];
    }

    /* Rehash all hash slots */
    for (i = 0; i < old_table.bucket_num; i++) {
        /* Duplicated keys in hash table which can cause eternal
         * hashing collision! Be careful of that!
         */
        assert(!cuckoo_hash_put(table, old_table.buckets[i]->key, old_table.buckets[i]->value));
    }

    free(old_table.items);
    free(old_table.buckets);
}

void cuckoo_insert(uint32_t key, uint32_t value)
{
    if (value != 0) {
        /* Reject duplicated keys keeping from eternal collision */
        int status = cuckoo_hash_get(&cuckoo_hash_table, key);
        if (status == OCCUPIED) {
            return;
        } else {
            /* Insert into hash slots */
            if (cuckoo_hash_put(&cuckoo_hash_table, key, value) == -1) {
                //cuckoo_rehash(&cuckoo_hash_table);
                cuckoo_hash_put(&cuckoo_hash_table, key, value);
            }

            cuckoo_hash_table.num_items ++;
        }
    } else {
        /* Delete at the hash slot */
        cuckoo_hash_remove(&cuckoo_hash_table, key);
    }
}

void cuckoo_remove(uint32_t key) {
    cuckoo_hash_remove(&cuckoo_hash_table, key);
#ifdef CUCKOO_DEBUG
    show_hash_table();
#endif
}

int cuckoo_lookup(uint32_t key)
{
    return cuckoo_hash_get(&cuckoo_hash_table, key);
}

int cuckoo_hashing_init(size_t size)
{
    int i;

    /* Otain the hash table size, it must be the power of 2 */
    cuckoo_hash_table.table_size = upper_pow_of_2(size);

    /* Make rehashing happen */
    //hash_table.slot_num /= 4;
    cuckoo_hash_table.items = calloc(cuckoo_hash_table.table_size, sizeof(struct hash_entry));
    if (cuckoo_hash_table.items == NULL) {
        return -1;
    }

    /* Allocate hash buckets associated with slots */
    cuckoo_hash_table.bucket_num = cuckoo_hash_table.table_size / ASSOC_WAY;
    cuckoo_hash_table.buckets = malloc(cuckoo_hash_table.bucket_num * sizeof(struct hash_entry *));

    if (cuckoo_hash_table.buckets == NULL) {
        free(cuckoo_hash_table.items);
        return -1;
    }

    for (i = 0; i < cuckoo_hash_table.bucket_num; i++) {
        cuckoo_hash_table.buckets[i] = &cuckoo_hash_table.items[i * ASSOC_WAY];
    }

    return 0;
}




