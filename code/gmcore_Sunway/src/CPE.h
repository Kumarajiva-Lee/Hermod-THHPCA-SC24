#ifndef CPE_H_INCLUDED 
#define CPE_H_INCLUDED 1

#include "crts.hpp"

#define COL_MASK 0x07
#define ROW_MASK (COL_MASK << 3)
#define COL(myid) ((myid) & COL_MASK)
#define ROW(myid) (((myid) >> 3) & COL_MASK)
#define TID(rid, cid) (((rid) << 3) | (cid))

#define NCPES       64
#define NCPES_SHIFT 6
#define NCPES_MASK  0x3F
#define NCOLS       8
#define NCOLS_SHIFT 3
#define NCOLS_MASK  0x07
#define NROWS       8
#define NROWS_SHIFT 3
#define NROWS_MASK  0x07

#define SYNC_ROW  CRTS_ssync_row()  // athread_syn in 1.0
#define SYNC_COL  CRTS_ssync_col()
#define SYNC_ALL  CRTS_ssync_array()

#define CPE_ID CRTS_tid // _MYID in 1.0


/***********************  For 2.0 porting  *********************/
// declare DMA 
#define DMA_DECL_DESC(mode) \
    crts_rply_t dma_reply_ = 0; \
    long dma_counter_ = 0;

#define DMA_WAIT_ALL CRTS_dma_wait_value(&dma_reply_, dma_counter_)
#define DMA_BARRIER  CRTS_dma_barrier()

#ifdef __cplusplus // pointer cast
#define CRTS_DMA_IGET(ldm, mem, size, reply) CRTS_dma_iget(ldm, const_cast<void*>(static_cast<const void*>(mem)), size, reply)
#define CRTS_DMA_IPUT(mem, ldm, size, reply) CRTS_dma_iput(mem, const_cast<void*>(static_cast<const void*>(ldm)), size, reply)
#define CRTS_DMA_GET(ldm, mem, size)         CRTS_dma_get(ldm, const_cast<void*>(static_cast<const void*>(mem)), size)
#define CRTS_DMA_PUT(mem, ldm, size)         CRTS_dma_put(mem, const_cast<void*>(static_cast<const void*>(ldm)), size)
#else
#define CRTS_DMA_IGET(ldm, mem, size, reply) CRTS_dma_iget(ldm, (void*)((uintptr_t)(mem)), size, reply)
#define CRTS_DMA_IPUT(mem, ldm, size, reply) CRTS_dma_iput(mem, (void*)((uintptr_t)(ldm)), size, reply)
#define CRTS_DMA_GET(ldm, mem, size)         CRTS_dma_get(ldm, (void*)((uintptr_t)(mem)), size)
#define CRTS_DMA_PUT(mem, ldm, size)         CRTS_dma_put(mem, (void*)((uintptr_t)(ldm)), size)
#endif

#define DMA_IREAD(mem, ldm, size) {\
    CRTS_DMA_IGET(ldm, mem, size, &dma_reply_); \
    ++dma_counter_; \
}

#define DMA_READ(mem, ldm, size) { \
    DMA_IREAD(mem, ldm, size); \
    DMA_WAIT_ALL; \
}

#define DMA_SREAD(mem, ldm, size) { \
    DMA_IREAD(mem, ldm, size); \
    DMA_BARRIER; \
}

#define DMA_IWRITE(mem, ldm, size) {\
    CRTS_DMA_IPUT(mem, ldm, size, &dma_reply_); \
    ++dma_counter_; \
}

#define DMA_WRITE(mem, ldm, size) { \
    DMA_IWRITE(mem, ldm, size); \
    DMA_WAIT_ALL; \
}

#define DMA_SWRITE(mem, ldm, size) { \
    DMA_IWRITE(mem, ldm, size); \
    DMA_BARRIER; \
}

#define DMA_READ_NODESC(mem, ldm, size) CRTS_DMA_GET(ldm, mem, size)
#define DMA_WRITE_NODESC(mem, ldm, size) CRTS_DMA_PUT(mem, ldm, size)

/* separation of replies */
#define DMA_DECL_DESC_2(mode) \
    crts_rply_t dma_reply_2_ = 0; \
    long dma_counter_2_ = 0;

#define DMA_WAIT_2 CRTS_dma_wait_value(&dma_reply_2_, dma_counter_2_)

#define DMA_IREAD_2(mem, ldm, size) {\
    CRTS_DMA_IGET(ldm, mem, size, &dma_reply_2_); \
    ++dma_counter_2_; \
}

#define DMA_IWRITE_2(mem, ldm, size) {\
    CRTS_DMA_IPUT(mem, ldm, size, &dma_reply_2_); \
    ++dma_counter_2_; \
}

inline void mem_fence() {
    asm volatile("memb" ::: "memory");
}

#ifdef __cplusplus
template <typename T>
inline T *remote_ldm_addr(T *ptr, uintptr_t id) {
    uintptr_t x = reinterpret_cast<uintptr_t>(ptr), y = id << 20, z = 1ull << 45;
    return reinterpret_cast<T *>(x | y | z);
}
#endif

#endif
