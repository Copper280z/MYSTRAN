/* Internal AMD stub: identity permutation. Selectable via -Dordering=internal,
   or used as last-resort fallback when SuiteSparse AMD/METIS aren't available. */

#include "internal.h"

ldlt_status ldlt_amd_stub(int32_t n, const int32_t *Ap, const int32_t *Ai, int32_t *perm) {
    (void)Ap; (void)Ai;
    for (int32_t i = 0; i < n; ++i) perm[i] = i;
    return LDLT_OK;
}
