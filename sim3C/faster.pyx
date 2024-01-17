#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, overflowcheck=False, initializedcheck=False, nonecheck=False

from libc.stdint cimport uint8_t, uint32_t, uint64_t, UINT32_MAX
from libc.stdlib cimport malloc, free
from libc.string cimport strlen
from cpython cimport array

#
# The following random number generator depends on the C-implementation of the PCG random number generator
#
# Website: http://www.pcg-random.org/
# Recent release: https://www.pcg-random.org/downloads/pcg-c-0.94.zip
# Github: https://github.com/imneme/pcg-c
#

cdef extern from "pcg_variants.h":

    struct pcg_state_setseq_64:
        uint64_t state
        uint64_t inc

    ctypedef pcg_state_setseq_64 pcg32_random_t

    void pcg_setseq_64_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq) nogil
    uint32_t pcg32_random_r(pcg32_random_t* rng) nogil
    uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound) nogil


cdef class PCGRandomState:
    # Python wrapper for holding the random state

    cdef pcg32_random_t rng
    cdef uint64_t state, inc

    def __init__(self, state=None, inc=None):
        if state is not None and inc is not None:
            self.state = state
            self.inc = inc
        self.seed(self.state, self.inc)

    cdef void seed(self, uint64_t state, uint64_t inc):
        pcg_setseq_64_srandom_r(&self.rng, state, inc)


cdef uint8_t DELETION = 45  # Int value for the deletion byte '-'

cdef class PCGRandom:
    """
    Main class for access to random number generation.
    """

    cdef PCGRandomState random_state

    def __init__(self, uint64_t state, uint64_t inc):
        self.random_state = PCGRandomState(state, inc)

    cpdef float uniform(self):
        """
        Generate a random float in the range [0, 1)
        :return: a random unsigned 32bit float
        """
        return random_uniform_ptr(&self.random_state.rng)

    cpdef uint32_t integer(self, uint32_t high):
        """
        Generate a random integer in the range [0, high)
        :param high: highest integer to return
        :return: a random unsigned 32bit int
        """
        return random_integer_ptr(&self.random_state.rng, high)

    cpdef uint32_t[::1] integers(self, uint32_t high, uint32_t size):
        """
        Generate an array of random integers in the range [0, high)
        :param high: highest integer to return
        :param size: the number of values to return as an array
        :return: a MemoryView on a Cython array
        """
        return _random_integer_array(&self.random_state.rng, high, size)

    cpdef bytes bytes(self, char* alphabet, uint32_t size):
        """
        Generate a random string of bytes from the given alphabet
        :param alphabet: the alphabet to use
        :param size: the number of bytes to generate
        :return: a bytes object
        """
        return random_bytes_ptr(&self.random_state.rng, alphabet, size)

    cpdef knockout_nucleotide(self, uint8_t nt):
        """
        Generate a random nucleotide that is not the given nucleotide
        :param nt: the excluded nucleotide
        :return: a single unsigned 8bit ordinal
        """
        return random_knockout_nucleotide_ptr(&self.random_state.rng, nt)

    cpdef nucleotide(self):
        """
        Generate a random nucleotide (A, C, G or T)
        :return: a single unsigned 8bit integer ordinal
        """
        return random_nucleotide_ptr(&self.random_state.rng)

    cpdef void parse_error(self,
                           uint8_t[::1] qual,
                           uint8_t[::1] seq,
                           double[::1] perr,
                           uint8_t n_symb):
        """
        When analyzed, sequences are potentially modified by the simulated quality scores.
        Beginning with the basic transcription from Art C/C++ code, this method has been reimplemented to use
        Numpy for speed improvements, but does not employ Numba as we must respect the existing random state.
    
        :param qual: quality scores, modified in place
        :param seq: bytearray of DNA sequence, to be modified in place
        :param perr: probability of an error
        :param n_symb: ambiguous symbol integer value
        """
        cdef uint32_t i
        cdef pcg32_random_t* rng = &self.random_state.rng

        for i in range(seq.shape[0]):
            if seq[i] == n_symb:
                qual[i] = 1
            elif self.uniform() < perr[qual[i]]:
                seq[i] = random_knockout_nucleotide_ptr(rng, seq[i])

    # cpdef int64_t[::1] get_quals(self,
    #                              int64_t[:, :, ::1] qual_dist_for_symb,
    #                              int64_t read_len,
    #                              int64_t max_dist_number):
    #     return get_qual_from_dist(&self.random_state.rng, qual_dist_for_symb, read_len, max_dist_number)


# cdef array.array INT64_ARRAY_TEMPLATE = array.array('l', [])
#
# cdef int64_t[::1] get_qual_from_dist(pcg32_random_t* random_state,
#                                      int64_t[:, :, ::1] qual_dist_for_symb,
#                                      int64_t read_len,
#                                      int64_t max_dist_number):
#     cdef uint32_t[::1] rv_list
#     cdef int64_t[::1] quals
#     cdef int64_t[:, ::1] qcdf
#     cdef uint32_t i, j
#
#     rvals = _random_integer_array(random_state, max_dist_number, size=read_len)
#     quals = array.clone(INT64_ARRAY_TEMPLATE, read_len, zero=<bint>True)
#     for i in range(read_len):
#         qcdf = qual_dist_for_symb[i]
#         j = np.searchsorted(qcdf[:, 0], 1+rvals[i])
#         quals[i] = qcdf[j, 1]
#
#     return quals


    # cpdef int32_t simulate_indels(self,
    #                               uint32_t[::1] indel_pos,
    #                               uint8_t[::1] indel_nuc,
    #                               double[::1] del_rate,
    #                               double[::1] ins_rate,
    #                               int32_t read_len,
    #                               int32_t temp_len):
    #
    #     cdef int32_t mut_len
    #
    #     # step 1
    #     mut_len = self.get_indel(indel_pos, indel_nuc, del_rate, ins_rate, read_len)
    #     # step 2: ensure that this read will fit within the extent of the template
    #     if read_len - mut_len > temp_len:
    #         mut_len = self.get_indel_2(indel_pos, indel_nuc, del_rate, ins_rate, read_len)
    #
    #     return mut_len
    #
    # cdef int32_t get_indel(self,
    #                        uint32_t[::1] indel_pos,
    #                        uint8_t[::1] indel_nuc,
    #                        double[::1] del_rate,
    #                        double[::1] ins_rate,
    #                        int32_t read_len):
    #     cdef uint32_t pos, n
    #     cdef int32_t i, j, ins_len, del_len
    #     cdef pcg32_random_t* rng = &self.random_state.rng
    #
    #     ins_len = 0
    #     del_len = 0
    #
    #     indel_pos[:] = 0
    #     indel_nuc[:] = 0
    #
    #     n = 0
    #     # deletion
    #     for i in range(del_rate.shape[0]-1, -1, -1):
    #         if del_rate[i] >= random_uniform_ptr(rng):
    #             del_len = i+1
    #             j = i
    #             while j >= 0:
    #                 # invalid deletion positions: 0 or read_len-1
    #                 pos = random_integer_ptr(rng, read_len)
    #                 if pos == 0:
    #                     continue
    #                 if not in_array(pos, indel_pos, n):
    #                     indel_pos[n] = pos
    #                     indel_nuc[n] = DELETION
    #                     n += 1
    #                     j -= 1
    #             break
    #
    #     # # insertion
    #     for i in range(ins_rate.shape[0]-1, -1, -1):
    #         # ensure that enough unchanged positions for mutation
    #         if read_len - del_len - ins_len < i+1:
    #             continue
    #         if ins_rate[i] >= random_uniform_ptr(rng):
    #             ins_len = i+1
    #             j = i
    #             while j >= 0:
    #                 pos = random_integer_ptr(rng, read_len)
    #                 if not in_array(pos, indel_pos, n):
    #                     indel_pos[n] = pos
    #                     indel_nuc[n] = random_nucleotide_ptr(rng)
    #                     n += 1
    #                     j -= 1
    #             break
    #
    #     return ins_len - del_len

    # cdef int32_t get_indel_2(self,
    #                         uint32_t[::1] indel_pos,
    #                         uint8_t[::1] indel_nuc,
    #                         double[::1] del_rate,
    #                         double[::1] ins_rate,
    #                         int32_t read_len):
    #     """
    #     Generate insertion and deletions
    #     :return: net change in length, i.e. insertion_length - deletion_length
    #     """
    #     cdef uint32_t pos, n
    #     cdef int32_t i, j, ins_len, del_len, size
    #     cdef pcg32_random_t* rng = &self.random_state.rng
    #
    #     ins_len = 0
    #     del_len = 0
    #
    #     indel_pos[:] = 0
    #     indel_nuc[:] = 0
    #
    #     n = 0
    #     # deletion
    #     for i in range(del_rate.shape[0]-1, -1, -1):
    #         if del_rate[i] >= random_uniform_ptr(rng):
    #             del_len = i+1
    #             j = i
    #             while j >= 0:
    #                 # invalid deletion positions: 0 or read_len-1
    #                 pos = random_integer_ptr(rng, read_len)
    #                 if pos == 0:
    #                     continue
    #                 if not in_array(pos, indel_pos, n):
    #                     indel_pos[n] = pos
    #                     indel_nuc[n] = -1
    #                     n += 1
    #                     j -= 1
    #             break
    #
    #     # insertion
    #     for i in range(ins_rate.shape[0]-1, -1, -1):
    #         # ensure that enough unchanged positions for mutation
    #         if read_len - del_len - ins_len < i+1:
    #             continue
    #         if ins_rate[i] >= random_uniform_ptr(rng):
    #             ins_len = i+1
    #             j = i
    #             while j >= 0:
    #                 pos = random_integer_ptr(rng, read_len)
    #                 if not in_array(pos, indel_pos, n):
    #                     indel_pos[n] = pos
    #                     indel_nuc[n] = random_nucleotide_ptr(rng)
    #                     n += 1
    #                     j -= 1
    #             break
    #
    #     return ins_len - del_len

    # cpdef void transfer_sequence(self,
    #                             uint32_t[::1] indel_pos,
    #                             uint8_t[::1] indel_nuc,
    #                             char[::1] seq_ref,
    #                             char[::1] seq_read,
    #                             uint32_t read_len):
    #     """
    #     From the reference (mother) sequence, generating the read's sequence along
    #     with the indels.
    #     """
    #     cdef uint32_t i, k, n, n_indels
    #     cdef char* seq_read_ptr = &seq_read[0]
    #     cdef char* seq_ref_ptr = &seq_ref[0]
    #
    #
    #     # count deletions and insertions
    #     n_indels = 0
    #     for i in range(read_len):
    #         if indel_nuc[i] > 0:
    #             n_indels += 1
    #         else:
    #             break
    #
    #     if n_indels == 0:
    #         strcpy(seq_read_ptr, seq_ref_ptr)
    #     else:
    #         n = 0
    #         k = 0
    #         i = 0
    #         while i < read_len:
    #             if not in_array(k, indel_pos, n_indels):
    #                 seq_read[n] = seq_ref[i]
    #                 n += 1
    #                 i += 1
    #                 k += 1
    #             elif indel_nuc[i] == DELETION:
    #                 # deletion
    #                 i += 1
    #                 k += 1
    #             else:
    #                 # insertion
    #                 seq_read[n] = indel_nuc[i]
    #                 n += 1
    #                 k += 1
    #
    #         while in_array(k, indel_pos, n_indels):
    #             seq_read[n] = indel_nuc[k]
    #             n += 1
    #             k += 1


cdef bint in_array(uint32_t val, uint32_t[::1] arr, uint32_t size) noexcept:
# cdef bint in_array(uint32_t val, uint32_t[::1] arr, uint32_t size):
    """
    Test if a value is contained in an array. 
    :param val: the value to find 
    :param arr: the array to search
    :param size: the length of the array
    :return: True the value was found, False otherwise
    """
    cdef uint32_t i
    for i in range(size):
        if arr[i] == val:
            return <bint>True
    return <bint>False


cdef float random_uniform_ptr(pcg32_random_t* random_state) noexcept:
# cdef float random_uniform_ptr(pcg32_random_t* random_state):
    """
    Generate a uniform random between [0, 1)
    :param random_state: a point to an initialised state
    :return: 32bit float
    """
    return <float>pcg32_random_r(random_state) / <float>UINT32_MAX


cdef uint32_t random_integer_ptr(pcg32_random_t* random_state, uint32_t high) noexcept:
# cdef uint32_t random_integer_ptr(pcg32_random_t* random_state, uint32_t high):
    """
    Generate a single random integer [0, high)
    :param random_state: a pointer to an initialised state
    :param high: highest value to return
    :return: the random intgger
    """
    return pcg32_boundedrand_r(random_state, high)


# template for array creation
cdef array.array UINT32_ARRAY_TEMPLATE = array.array('I', [])

cdef uint32_t[::1] _random_integer_array(pcg32_random_t* random_state, uint32_t high, uint32_t size) noexcept:
# cdef uint32_t[::1] _random_integer_array(pcg32_random_t* random_state, uint32_t high, uint32_t size):
    """
    Generate an array of random integers in the range [0, high)
    :param random_state: a pointer to an initialised state
    :param high: highest integer to return
    :param size: the number of values to return as an array
    :return: a MemoryView to a Cython array
    """
    cdef uint32_t i
    cdef uint32_t [::1] rvals

    rvals = array.clone(UINT32_ARRAY_TEMPLATE, size, zero=<bint>False)
    for i in range(size):
        rvals[i] = pcg32_boundedrand_r(random_state, high)
    return rvals


cdef bytes random_bytes_ptr(pcg32_random_t* random_state, char* alphabet, uint32_t size) noexcept:
# cdef bytes random_bytes_ptr(pcg32_random_t* random_state, char* alphabet, uint32_t size):
    """
    Generate a random string of bytes from the given alphabet
    :param random_state: a pointer to an initialised state 
    :param alphabet: the alphabet to use
    :param size: the number of bytes to generate
    :return: a bytes object
    """
    cdef uint32_t i
    cdef uint32_t high = strlen(alphabet)
    cdef bytes py_string

    cdef char* randoms = <char *> malloc((size+1) * sizeof(char))
    if not randoms:
        return None

    for i in range(size):
        randoms[i] = alphabet[pcg32_boundedrand_r(random_state, high)]
    randoms[size] = 0

    try:
        py_string = randoms[:size]
    finally:
        free(randoms)

    return py_string


cdef uint8_t random_knockout_nucleotide_ptr(pcg32_random_t* random_state, uint8_t nt) noexcept:
# cdef uint8_t random_knockout_nucleotide_ptr(pcg32_random_t* random_state, uint8_t nt):
    """
    Generate a random nucleotide that is not the given nucleotide
    :param random_state: a pointer to an initialised state
    :param nt: the nucleotide to exclude
    :return: a random nucleotide
    """
    cdef uint32_t[3] A_KO = [67, 71, 84]
    cdef uint32_t[3] C_KO = [65, 71, 84]
    cdef uint32_t[3] G_KO = [65, 67, 84]
    cdef uint32_t[3] T_KO = [65, 67, 71]
    cdef uint32_t ix

    ix = random_integer_ptr(random_state, 3)

    if nt == 65:
        rnt = A_KO[ix]
    elif nt == 67:
        rnt = C_KO[ix]
    elif nt == 71:
        rnt = G_KO[ix]
    elif nt == 84:
        rnt = T_KO[ix]
    else:
        rnt = 0
    return rnt


cdef uint8_t random_nucleotide_ptr(pcg32_random_t* random_state) noexcept:
# cdef uint8_t random_nucleotide_ptr(pcg32_random_t* random_state):
    """
    Generate a random nucleotide (A, C, G or T)
    :param random_state: a pointer to an initialised state
    :return: a random nucleotide
    """
    cdef uint32_t[4] ALL_NT = [65, 67, 71, 84]
    cdef uint32_t ix

    ix = random_integer_ptr(random_state, 4)
    return ALL_NT[ix]

#
# Miscellaneous functions
#

cpdef bytes qualities_to_bytes(uint8_t[::1] qualities):
    """
    Convert an array of quality scores to a bytes object
    :param qualities: array if ints
    :return: byte string representation
    """
    cdef uint32_t i
    cdef uint32_t n_qual = qualities.shape[0]
    cdef bytes py_out
    cdef char* phred_bytes = <char *> malloc((n_qual + 1) * sizeof(char))
    if not phred_bytes:
        return None
    for i in range(n_qual):
        phred_bytes[i] = qualities[i] + 33
    try:
        py_out = phred_bytes[:n_qual]
    finally:
        free(phred_bytes)
    return py_out
