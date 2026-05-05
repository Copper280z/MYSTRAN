#include "internal.h"
#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_host.h>
#elif defined(__linux__)
#include <sys/sysinfo.h>
#endif

#if defined(_MSC_VER)
__declspec(thread) static ldlt_error_trap *g_error_trap = NULL;
#else
static __thread ldlt_error_trap *g_error_trap = NULL;
#endif

void ldlt_push_error_trap(ldlt_error_trap *trap) {
    if (!trap) return;
    trap->status = LDLT_OK;
    trap->prev = g_error_trap;
    g_error_trap = trap;
}

void ldlt_pop_error_trap(ldlt_error_trap *trap) {
    if (trap && g_error_trap == trap)
        g_error_trap = trap->prev;
}

void ldlt_raise_status(ldlt_status status) {
    if (g_error_trap) {
        g_error_trap->status = status;
        longjmp(g_error_trap->env, 1);
    }
    fprintf(stderr, "ldlt: unrecoverable error: %s\n", ldlt_status_str(status));
    abort();
}

void *ldlt_xmalloc(size_t bytes) {
    if (bytes == 0) bytes = 1;
    void *p = malloc(bytes);
    if (!p) {
        fprintf(stderr, "ldlt: malloc(%zu) failed\n", bytes);
        ldlt_raise_status(LDLT_ERR_NOMEM);
    }
    return p;
}
void *ldlt_xcalloc(size_t n, size_t sz) {
    if (n == 0 || sz == 0) { n = 1; sz = 1; }
    void *p = calloc(n, sz);
    if (!p) {
        fprintf(stderr, "ldlt: calloc(%zu,%zu) failed\n", n, sz);
        ldlt_raise_status(LDLT_ERR_NOMEM);
    }
    return p;
}
void *ldlt_xrealloc(void *p, size_t bytes) {
    void *q = realloc(p, bytes ? bytes : 1);
    if (!q) {
        fprintf(stderr, "ldlt: realloc(%zu) failed\n", bytes);
        ldlt_raise_status(LDLT_ERR_NOMEM);
    }
    return q;
}

uint64_t ldlt_available_memory_bytes(void) {
#if defined(__APPLE__)
    mach_port_t host = mach_host_self();
    vm_size_t page_size = 0;
    vm_statistics64_data_t vmstat;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;

    if (host_page_size(host, &page_size) != KERN_SUCCESS)
        return 0;
    if (host_statistics64(host, HOST_VM_INFO64, (host_info64_t)&vmstat,
                          &count) != KERN_SUCCESS)
        return 0;

    uint64_t pages = (uint64_t)vmstat.free_count
                   + (uint64_t)vmstat.inactive_count
                   + (uint64_t)vmstat.speculative_count;
    return pages * (uint64_t)page_size;
#elif defined(__linux__)
    struct sysinfo si;
    if (sysinfo(&si) != 0)
        return 0;
    return (uint64_t)si.freeram * (uint64_t)si.mem_unit;
#else
    return 0;
#endif
}

void ldlt_options_default(ldlt_options *opt) {
    opt->ordering = LDLT_ORDER_AMD;
    opt->user_perm = NULL;
    opt->relax_supernodes = 16;
    opt->nthreads = 0;
    opt->pivot = LDLT_PIVOT_NONE;
    opt->pivot_threshold = 0.0;
    opt->require_pd = 0;
}

const char *ldlt_status_str(ldlt_status s) {
    switch (s) {
    case LDLT_OK: return "ok";
    case LDLT_ERR_INPUT: return "invalid input";
    case LDLT_ERR_NOMEM: return "out of memory";
    case LDLT_ERR_INDEFINITE: return "matrix indefinite";
    case LDLT_ERR_SINGULAR: return "matrix singular";
    case LDLT_ERR_UNSUPPORTED: return "feature not supported";
    }
    return "unknown";
}

int32_t ldlt_n(const ldlt_symbolic *S)      { return S ? S->n : 0; }
int32_t ldlt_nsuper(const ldlt_symbolic *S) { return S ? S->nsuper : 0; }
int64_t ldlt_lnz(const ldlt_symbolic *S)    { return S ? S->lnz : 0; }
double  ldlt_flops(const ldlt_symbolic *S)  { return S ? S->flops : 0.0; }
int32_t ldlt_nlevels(const ldlt_symbolic *S) { return S ? S->nlevels : 0; }
int32_t ldlt_max_level_width(const ldlt_symbolic *S) { return S ? S->max_level_width : 0; }

void ldlt_free_symbolic(ldlt_symbolic *S) {
    if (!S) return;
    free(S->perm); free(S->iperm);
    free(S->parent); free(S->postorder);
    free(S->col_to_super);
    if (S->super) {
        for (int32_t i = 0; i < S->nsuper; ++i) free(S->super[i].sn_perm);
        free(S->super);
    }
    free(S->rows_pool);
    free(S->super_parent);
    free(S->child_ptr); free(S->child_idx);
    free(S->edge_ptr); free(S->edge_pos);
    free(S->upd_ptr); free(S->upd_idx);
    free(S->panel_off); free(S->d_off);
    free(S);
}

void ldlt_free_numeric(ldlt_numeric *N) {
    if (!N) return;
    free(N->L); free(N->D);
    free(N->piv_size);
    free(N->piv_var_off);
    free(N->piv_vars);
    free(N->piv_d_off);
    free(N->piv_D);
    free(N->piv_l_row_off);
    free(N->piv_l_val_off);
    free(N->piv_l_rows);
    free(N->piv_l_vals);
    free(N);
}
