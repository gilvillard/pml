#include <flint/ulong_extras.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/fft_small.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod.h>

#include "nmod_poly_mat_utils.h" // for rand
#include "nmod_poly_mat_multiply.h"



/***********************
*  bit reversed copy  *
***********************/

long RevInc(long a, long k)
{
    long j, m;

    j = k;
    m = 1L << (k-1);

    while (j && (m & a)) {
        a ^= m;
        m >>= 1;
        j--;
    }
    if (j) a ^= m;
    return a;
}

// indices initialized with length >= k
void brc_indices(mp_limb_t * indices, long k)
{
    const long n = (1L << k);
    for (long i = 0, j = 0; i < n; i++, j = RevInc(j, k))
        indices[i] = j;
}




// Discrete Fourier transform
//Input:
//    :coeffs: coefficients [mat_0,...,mat_{d-1}] of some polynomial matrix of degree < d
//    :lg: nonnegative integer (log, or length); below d == 2**lg
//    :ws: table of powers [1],[1,w**(d/2-1)],...,[1,w**2,w**4,...,w**(d-2)],[1,w,w**2,...,w**(d-1)] of principal d-th root of unity w
//    :winvs: precomputed inverses for the above powers
//Output:
//    evaluations [mat(1),mat(w),...,mat(w**(d-1))] (not in this order!!)
//Complexity:
//    3/2 d log_2(d)
//"""
void _dft(nmod_mat_struct * coeffs, ulong lg, mp_limb_t ** ws, mp_limb_t ** winvs, ulong m, ulong n, nmod_t mod)
{
    if (lg == 0)
        return;

    // two polynomials for recursive calls, of length < len == 2**(lg-1)
    ulong len = 1 << (lg-1);
    // gm = f rem x**len - 1, hence gm(w**(2*k)) == f(w**(2*k))
    // gp = (f rem x**len + 1)(w x), hence gp(w**(2*k)) == f(w * w**(2*k)) = f(w**(2*k+1))

    // in place, compute:
    // gm = [f[k] + f[k+len] for k in range(len)] in first half of mat
    // gp = [(f[k] - f[k+len]) * ws[k] for k in range(len)] in second half of mat
    mp_limb_t buf;
    for (ulong k = 0; k < len; k++)
    {
        for (ulong i = 0; i < m; i++)
        {
            for (ulong j = 0; j < n; j++)
            {
                buf = (coeffs+k)->rows[i][j];
                // nmod_mat_add(coeffs+k, buf, coeffs+(k+len)); -->
                (coeffs+k)->rows[i][j] = _nmod_add(buf, (coeffs+k+len)->rows[i][j], mod);
                // nmod_mat_sub(coeffs+(k+len), coeffs+(k+len), buf); -->
                (coeffs+k+len)->rows[i][j] = _nmod_sub(buf, (coeffs+k+len)->rows[i][j], mod);
                // nmod_mat_scalar_mul(coeffs+(k+len), coeffs+(k+len), ws[lg][k]); -->
                (coeffs+k+len)->rows[i][j] = n_mulmod_shoup(ws[lg][k], (coeffs+k+len)->rows[i][j], winvs[lg][k], mod.n);
            }
        }
    }

    _dft(coeffs,lg-1,ws,winvs,m,n,mod);
    _dft(coeffs+len,lg-1,ws,winvs,m,n,mod);
    // TODO reorganize here? or leave as such?
}


// matrix given as matrix polynomial, degree < len
// len is a power of two
// rt is a principal len-th root of unity
// computes dft w.r.t len
void dft(nmod_mat_poly_t mat, ulong len, mp_limb_t rt)
{
    // ceiling(log_2(len)) == log_2(len)
    const ulong order = n_clog(len, 2);

    // array of powers of rt
    mp_limb_t ** ws = flint_malloc((order+1) * sizeof(mp_limb_t *));
    mp_limb_t ** winvs = flint_malloc((order+1) * sizeof(mp_limb_t *));
        //mp_limb_t winv = n_mulmod_precomp_shoup(ws[lg][k], buf->mod.n);
    ws[order] = _nmod_vec_init(len);
    winvs[order] = _nmod_vec_init(len);
    ws[order][0] = 1UL;
    winvs[order][0] = n_mulmod_precomp_shoup(1UL, mat->mod.n);
    for (ulong k = 1; k < len; k++)
    {
        ws[order][k] = nmod_mul(ws[order][k-1], rt, mat->mod);
        winvs[order][k] = n_mulmod_precomp_shoup(ws[order][k], mat->mod.n);
    }
    for (slong i = order - 1; i >= 0; i--)
    {
        len = len/2;
        ws[i] = _nmod_vec_init(len);
        winvs[i] = _nmod_vec_init(len);
        for (ulong k = 0; k < len; k++)
        {
            ws[i][k] = ws[i+1][2*k];
            winvs[i][k] = n_mulmod_precomp_shoup(ws[i][k], mat->mod.n);
        }
    }

    // careful here len == 1 or 0, something like this

    _dft(mat->coeffs, order, ws, winvs, mat->r, mat->c, mat->mod);

    for (ulong i = 0; i <= order; i++)
        _nmod_vec_clear(ws[i]);
    for (ulong i = 0; i <= order; i++)
        _nmod_vec_clear(winvs[i]);
    flint_free(ws);
    flint_free(winvs);
}

int test_matrix_dft(ulong m, ulong n, ulong len)
{
    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, time(NULL), time(NULL));

    // prime field Z/nZ
    mp_limb_t modulus = 1108307720798209;  // 7 * 3^2 * 2**44 + 1
    nmod_t mod;
    nmod_init(&mod, modulus);
    // root of unity of max order (n-1); 2**44-th root of unity
    mp_limb_t gen = 1101616563820748;
    mp_limb_t rt = nmod_pow_ui(gen, 7*3*3, mod);

    mp_limb_t w = nmod_pow_ui(rt, (1UL<<(44-n_clog(len, 2))), mod);

    const ulong order = n_clog(len, 2);
    ulong length = len;
    mp_limb_t ** ws = flint_malloc((order+1) * sizeof(mp_limb_t *));
    ws[order] = _nmod_vec_init(length);
    ws[order][0] = 1UL;
    for (ulong k = 1; k < length; k++)
        ws[order][k] = nmod_mul(ws[order][k-1], w, mod);
    for (slong i = order - 1; i >= 0; i--)
    {
        length = length/2;
        ws[i] = _nmod_vec_init(length);
        for (ulong k = 0; k < length; k++)
            ws[i][k] = ws[i+1][2*k];
    }

    // root of unity in bit reverse order
    mp_limb_t * brc_ind = _nmod_vec_init(len);
    brc_indices(brc_ind, order);
    mp_limb_t * brws = flint_malloc(len * sizeof(mp_limb_t));
    for (ulong k = 0; k < len; k++)
        brws[k] = ws[order][brc_ind[k]];
    _nmod_vec_clear(brc_ind);

    // random matrix
    nmod_mat_poly_t matA;
    nmod_mat_poly_init(matA, m, n, mod.n);
    nmod_mat_poly_rand(matA, state, len);

    // dft, v1
    nmod_mat_poly_t matB;
    nmod_mat_poly_init(matB, m, n, mod.n);
    nmod_mat_poly_set(matB, matA);
    dft(matB, len, w);

    // naive Horner
    nmod_mat_poly_t matC;
    nmod_mat_poly_init(matC, m, n, mod.n);
    nmod_mat_poly_fit_length(matC, len);
    _nmod_mat_poly_set_length(matC, len);
    for (ulong k = 0; k < len; k++)
        nmod_mat_poly_evaluate_nmod(matC->coeffs+k, matA, brws[k]);

    // check
    for (ulong k = 0; k < len; k++)
        if (!nmod_mat_equal(matB->coeffs+k, matC->coeffs+k))
        {
            printf("\n!!!WRONG EVAL!!! %ld\n", k);
            return 0;
        }

    nmod_mat_poly_clear(matA);
    nmod_mat_poly_clear(matB);
    nmod_mat_poly_clear(matC);
    return 1;
}

/*--------------------------------------------------------------*/
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
void time_matrix_dft(ulong m, ulong n, ulong len)
{
    flint_rand_t state;
    flint_randinit(state);

    // prime field Z/nZ
    mp_limb_t modulus = 1108307720798209;  // 7 * 3^2 * 2**44 + 1
    nmod_t mod;
    nmod_init(&mod, modulus);
    // root of unity of max order (n-1); 2**44-th root of unity
    mp_limb_t gen = 1101616563820748;
    mp_limb_t rt = nmod_pow_ui(gen, 7*3*3, mod);

    double tt1 = 0.0;
    double tt_eval = 0.0;

    {
        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, m, n, modulus);
        nmod_poly_mat_rand(A, state, len);
        tt1 = 0.0;
        long nb_iter1 = 0;
        while (tt1 < 0.5)
        {
            nmod_mat_poly_t matA;
            nmod_mat_poly_init(matA, A->r, A->c, A->modulus);
            nmod_mat_poly_set_from_poly_mat(matA, A);
            clock_t t1 = clock();
            dft(matA, len, nmod_pow_ui(rt, (1UL<<(44-n_clog(len, 2))), mod));
            tt1 += (double)(clock() - t1) / CLOCKS_PER_SEC;
            nmod_mat_poly_clear(matA);
            nb_iter1++;
        }
        tt1 = tt1/nb_iter1;
        printf("%.1e\t", tt1);
    }

    {
        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, m, n, modulus);
        nmod_poly_mat_rand(A, state, len/2);
        nmod_poly_mat_t B;
        nmod_poly_mat_init(B, n, m, modulus);
        nmod_poly_mat_rand(B, state, len/2);
        nmod_poly_mat_t C;
        nmod_poly_mat_init(C, m, m, modulus);

        double tt2 = 0.0;
        long nb_iter2 = 0;
        while (tt2 < 0.5)
        {
            clock_t t2 = clock();
            tt_eval += nmod_poly_mat_mul_tft(C, A, B);
            tt2 += (double)(clock() - t2) / CLOCKS_PER_SEC;
            nb_iter2++;
        }
        tt2 = tt2 / nb_iter2;
        tt_eval = tt_eval / nb_iter2;
        printf("%.1e\t", tt_eval);
        printf("%.1e\t", tt2);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    printf("%.1e\n", tt_eval / tt1);

    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls tets                                              */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_set_num_threads(1);

    if (argc != 2)
    {
        printf("Usage: '%s i' with i in 1 (test only), 2 (bench only), 3 (both)\n", argv[0]);
        return 0;
    }

    if (atoi(argv[1]) == 1 || atoi(argv[1]) == 3)
    {
        printf("Testing (slow: has Horner eval):\n");
        for (ulong dim = 1; dim < 25; dim += 3)
        {
            printf("dim = %ld\nlen-s ok: ", dim);
            for (ulong len = 2; len < (1<<16)/(dim); len *= 2)
            {
                if ((len < (1L<<12)) & (! test_matrix_dft(dim, dim, len)))
                    return 0;
                printf("%ld, ", len);
                fflush(stdout);
            }
        }
        printf("\n");
        return 0;
    }
    else if (atoi(argv[1]) == 2 || atoi(argv[1]) == 3)
    {
        printf("Bench:\n");
        printf("dim\tlen\tnew\ttft\tfullmul\tratio tft/new\n");
        for (ulong dim = 1; dim < 5; dim += 1)
        {
            for (ulong len = 4; len < (1<<26)/(dim*dim); len *= 2)
            {
                printf("%ld\t%ld\t", dim, len);
                time_matrix_dft(dim, dim, len);
            }
        }
        for (ulong dim = 5; dim < 16; dim += 3)
        {
            for (ulong len = 4; len < (1<<26)/(dim*dim); len *= 2)
            {
                printf("%ld\t%ld\t", dim, len);
                time_matrix_dft(dim, dim, len);
            }
        }
        for (ulong dim = 16; dim < 256; dim *= 2)
        {
            for (ulong len = 4; len < (1<<26)/(dim*dim); len *= 2)
            {
                printf("%ld\t%ld\t", dim, len);
                time_matrix_dft(dim, dim, len);
            }
        }
        return 0;
    }
    else
    {
        printf("Usage: '%s i' with i in 1 (test only), 2 (bench only), 3 (both)\n", argv[0]);
        return 0;
    }
    return 0;
}
