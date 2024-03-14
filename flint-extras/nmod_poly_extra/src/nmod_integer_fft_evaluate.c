#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include "nmod_poly_fft.h"

/** fft dif evaluation, in place                                 
 * - len == 2**order, order >= 0
 * - input p has length at least len, seen as polynomial
 * - output p contains, in 0...len-1, the evaluations of trunc(p, len)
 *       at w**k, for 0 <= k < len, in bit reverse order
 * TODO : ASSUMES order suitable for w,                      
 *          and tables for 'order' already precomputed in F  
 *
 *  also assumes we can add a few terms before overflowing
 */
void _nmod_poly_dif_inplace_radix2_rec(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
    }
    else if (order == 2)
    {
        // for I == w**(2**(order-2)) square root of -1,
        //                            [1  1  1  1]
        // p <- [p[0] p[1] p[2] p[3]] [1 -1  I -I]
        //                            [1  1 -1 -1]
        //                            [1 -1 -I  I]
        // p[0] = (p[0] + p[2]) + (p[1] + p[3]);
        // p[1] = (p[0] + p[2]) - (p[1] + p[3]);
        // p[2] = (p[0] - p[2]) + I*(p[1] - p[3]);
        // p[3] = (p[0] - p[2]) - I*(p[1] - p[3]);
        mp_limb_t tmp02 = nmod_add(p[0], p[2], F->mod); // p[0] + p[2], [0,n-1]
        p[2] = nmod_sub(p[0], p[2], F->mod);    // p[0] - p[2], [0,n-1]
        mp_limb_t tmp13 = nmod_sub(p[1], p[3], F->mod); // p[1] - p[3], [0,n-1]
        p[1] = nmod_add(p[1], p[3], F->mod);  // p[1] + p[3], [0,2n-2]

        p[0] = nmod_add(tmp02, p[1], F->mod);  // p[0] now ok,  [0, n-1]
        p[1] = nmod_sub(tmp02, p[1], F->mod);  // p[1] now ok,  [0, n-1]

        // compute I*(p[1] - p[3]),  [0, n-1]
        NMOD_MUL_PRENORM(tmp13, tmp13 << F->mod.norm, F->tab_w[0][1], F->mod);
        p[3] = nmod_sub(p[2], tmp13, F->mod); // p[3] now ok, [0, n-1]
        p[2] = nmod_add(p[2], tmp13, F->mod); // p[2] now ok, [0, n-1]
    }
    else
    {
        const mp_ptr phi = p+(len/2);
        for (ulong k = 0; k < len/2; k++)
        {
            const ulong tmp = p[k];
            p[k] = nmod_add(p[k], phi[k], F->mod);
            phi[k] = nmod_sub(tmp, phi[k], F->mod);
            NMOD_MUL_PRENORM(phi[k], phi[k] << F->mod.norm, F->tab_w[order-2][k], F->mod);
        }
        _nmod_poly_dif_inplace_radix2_rec(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec(phi, len/2, order-1, F);
    }

}

void _nmod_poly_dif_inplace_radix2_rec_v2(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(p[0], p[1], F->mod);
    }
    else if (order == 2)
    {
        // for I == w**(2**(order-2)) square root of -1,
        //                            [1  1  1  1]
        // p <- [p[0] p[1] p[2] p[3]] [1 -1  I -I]
        //                            [1  1 -1 -1]
        //                            [1 -1 -I  I]
        // p[0] = (p[0] + p[2]) + (p[1] + p[3]);
        // p[1] = (p[0] + p[2]) - (p[1] + p[3]);
        // p[2] = (p[0] - p[2]) + I*(p[1] - p[3]);
        // p[3] = (p[0] - p[2]) - I*(p[1] - p[3]);
        mp_limb_t tmp02 = nmod_add(p[0], p[2], F->mod); // p[0] + p[2], [0,n-1]
        p[2] = nmod_sub(p[0], p[2], F->mod);    // p[0] - p[2], [0,n-1]
        mp_limb_t tmp13 = nmod_sub(p[1], p[3], F->mod); // p[1] - p[3], [0,n-1]
        p[1] = nmod_add(p[1], p[3], F->mod);  // p[1] + p[3], [0,2n-2]

        p[0] = nmod_add(tmp02, p[1], F->mod);  // p[0] now ok,  [0, n-1]
        p[1] = nmod_sub(tmp02, p[1], F->mod);  // p[1] now ok,  [0, n-1]

        // compute I*(p[1] - p[3]),  [0, n-1]
        NMOD_MUL_PRENORM(tmp13, tmp13 << F->mod.norm, F->tab_w[0][1], F->mod);
        p[3] = nmod_sub(p[2], tmp13, F->mod); // p[3] now ok, [0, n-1]
        p[2] = nmod_add(p[2], tmp13, F->mod); // p[2] now ok, [0, n-1]
    }
    else
    {
        // here order >= 3, len >= 8, try some unrolling
        const mp_ptr phi = p+(len/2);
        mp_limb_t tmp[4];
        for (ulong k = 0; k+3 < len/2; k+=4)
        {
            tmp[0] = p[k+0];
            tmp[1] = p[k+1];
            tmp[2] = p[k+2];
            tmp[3] = p[k+3];
            p[k+0] = nmod_add(tmp[0], phi[k+0], F->mod);
            p[k+1] = nmod_add(tmp[1], phi[k+1], F->mod);
            p[k+2] = nmod_add(tmp[2], phi[k+2], F->mod);
            p[k+3] = nmod_add(tmp[3], phi[k+3], F->mod);
            phi[k+0] = nmod_sub(tmp[0], phi[k+0], F->mod);
            phi[k+1] = nmod_sub(tmp[1], phi[k+1], F->mod);
            phi[k+2] = nmod_sub(tmp[2], phi[k+2], F->mod);
            phi[k+3] = nmod_sub(tmp[3], phi[k+3], F->mod);
            NMOD_MUL_PRENORM(phi[k+0], phi[k+0], F->tab_w[order-2][k+0] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+1], phi[k+1], F->tab_w[order-2][k+1] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+2], phi[k+2], F->tab_w[order-2][k+2] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+3], phi[k+3], F->tab_w[order-2][k+3] << F->mod.norm, F->mod);
        }
        _nmod_poly_dif_inplace_radix2_rec_v2(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_v2(phi, len/2, order-1, F);
    }
}


void _nmod_poly_dif_inplace_radix2_rec_v3(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(p[0], p[1], F->mod);
    }
    else if (order == 2)
    {
        // for I == w**(2**(order-2)) square root of -1,
        //                            [1  1  1  1]
        // p <- [p[0] p[1] p[2] p[3]] [1 -1  I -I]
        //                            [1  1 -1 -1]
        //                            [1 -1 -I  I]
        // p[0] = (p[0] + p[2]) + (p[1] + p[3]);
        // p[1] = (p[0] + p[2]) - (p[1] + p[3]);
        // p[2] = (p[0] - p[2]) + I*(p[1] - p[3]);
        // p[3] = (p[0] - p[2]) - I*(p[1] - p[3]);
        mp_limb_t tmp02 = nmod_add(p[0], p[2], F->mod); // p[0] + p[2], [0,n-1]
        p[2] = nmod_sub(p[0], p[2], F->mod);    // p[0] - p[2], [0,n-1]
        mp_limb_t tmp13 = nmod_sub(p[1], p[3], F->mod); // p[1] - p[3], [0,n-1]
        p[1] = nmod_add(p[1], p[3], F->mod);  // p[1] + p[3], [0,2n-2]

        p[0] = nmod_add(tmp02, p[1], F->mod);  // p[0] now ok,  [0, n-1]
        p[1] = nmod_sub(tmp02, p[1], F->mod);  // p[1] now ok,  [0, n-1]

        // compute I*(p[1] - p[3]),  [0, n-1]
        NMOD_MUL_PRENORM(tmp13, tmp13 << F->mod.norm, F->tab_w[0][1], F->mod);
        p[3] = nmod_sub(p[2], tmp13, F->mod); // p[3] now ok, [0, n-1]
        p[2] = nmod_add(p[2], tmp13, F->mod); // p[2] now ok, [0, n-1]
    }
    else if (order == 3)
    {
        // here order = 3, len = 8, try some unrolling
        const mp_limb_t tmp[4] = {p[0], p[1], p[2], p[3]};
        p[0] = nmod_add(tmp[0], p[4], F->mod);
        p[1] = nmod_add(tmp[1], p[5], F->mod);
        p[2] = nmod_add(tmp[2], p[6], F->mod);
        p[3] = nmod_add(tmp[3], p[7], F->mod);
        p[4] = nmod_sub(tmp[0], p[4], F->mod);
        p[5] = nmod_sub(tmp[1], p[5], F->mod);
        p[6] = nmod_sub(tmp[2], p[6], F->mod);
        p[7] = nmod_sub(tmp[3], p[7], F->mod);
        NMOD_MUL_PRENORM(p[4], p[4], F->tab_w[order-2][0] << F->mod.norm, F->mod);
        NMOD_MUL_PRENORM(p[5], p[5], F->tab_w[order-2][1] << F->mod.norm, F->mod);
        NMOD_MUL_PRENORM(p[6], p[6], F->tab_w[order-2][2] << F->mod.norm, F->mod);
        NMOD_MUL_PRENORM(p[7], p[7], F->tab_w[order-2][3] << F->mod.norm, F->mod);
        _nmod_poly_dif_inplace_radix2_rec_v2(p, 4, 2, F);
        _nmod_poly_dif_inplace_radix2_rec_v2(p+4, 4, 2, F);
    }
    else
    {
        // here order >= 4, len >= 16, try some unrolling
        const mp_ptr phi = p+(len/2);
        mp_limb_t tmp[8];
        for (ulong k = 0; k+7 < len/2; k+=8)
        {
            tmp[0] = p[k+0];
            tmp[1] = p[k+1];
            tmp[2] = p[k+2];
            tmp[3] = p[k+3];
            tmp[4] = p[k+4];
            tmp[5] = p[k+5];
            tmp[6] = p[k+6];
            tmp[7] = p[k+7];
            p[k+0] = nmod_add(tmp[0], phi[k+0], F->mod);
            p[k+1] = nmod_add(tmp[1], phi[k+1], F->mod);
            p[k+2] = nmod_add(tmp[2], phi[k+2], F->mod);
            p[k+3] = nmod_add(tmp[3], phi[k+3], F->mod);
            p[k+4] = nmod_add(tmp[4], phi[k+4], F->mod);
            p[k+5] = nmod_add(tmp[5], phi[k+5], F->mod);
            p[k+6] = nmod_add(tmp[6], phi[k+6], F->mod);
            p[k+7] = nmod_add(tmp[7], phi[k+7], F->mod);
            phi[k+0] = nmod_sub(tmp[0], phi[k+0], F->mod);
            phi[k+1] = nmod_sub(tmp[1], phi[k+1], F->mod);
            phi[k+2] = nmod_sub(tmp[2], phi[k+2], F->mod);
            phi[k+3] = nmod_sub(tmp[3], phi[k+3], F->mod);
            phi[k+4] = nmod_sub(tmp[4], phi[k+4], F->mod);
            phi[k+5] = nmod_sub(tmp[5], phi[k+5], F->mod);
            phi[k+6] = nmod_sub(tmp[6], phi[k+6], F->mod);
            phi[k+7] = nmod_sub(tmp[7], phi[k+7], F->mod);
            NMOD_MUL_PRENORM(phi[k+0], phi[k+0], F->tab_w[order-2][k+0] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+1], phi[k+1], F->tab_w[order-2][k+1] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+2], phi[k+2], F->tab_w[order-2][k+2] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+3], phi[k+3], F->tab_w[order-2][k+3] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+4], phi[k+4], F->tab_w[order-2][k+4] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+5], phi[k+5], F->tab_w[order-2][k+5] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+6], phi[k+6], F->tab_w[order-2][k+6] << F->mod.norm, F->mod);
            NMOD_MUL_PRENORM(phi[k+7], phi[k+7], F->tab_w[order-2][k+7] << F->mod.norm, F->mod);
        }
        _nmod_poly_dif_inplace_radix2_rec_v3(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_v3(phi, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_rec_v4(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(p[0], p[1], F->mod);
    }
    else if (order == 2)
    {
        // for I == w**(2**(order-2)) square root of -1,
        //                            [1  1  1  1]
        // p <- [p[0] p[1] p[2] p[3]] [1 -1  I -I]
        //                            [1  1 -1 -1]
        //                            [1 -1 -I  I]
        // p[0] = (p[0] + p[2]) + (p[1] + p[3]);
        // p[1] = (p[0] + p[2]) - (p[1] + p[3]);
        // p[2] = (p[0] - p[2]) + I*(p[1] - p[3]);
        // p[3] = (p[0] - p[2]) - I*(p[1] - p[3]);
        mp_limb_t tmp02 = nmod_add(p[0], p[2], F->mod); // p[0] + p[2], [0,n-1]
        p[2] = nmod_sub(p[0], p[2], F->mod);    // p[0] - p[2], [0,n-1]
        mp_limb_t tmp13 = nmod_sub(p[1], p[3], F->mod); // p[1] - p[3], [0,n-1]
        p[1] = nmod_add(p[1], p[3], F->mod);  // p[1] + p[3], [0,2n-2]

        p[0] = nmod_add(tmp02, p[1], F->mod);  // p[0] now ok,  [0, n-1]
        p[1] = nmod_sub(tmp02, p[1], F->mod);  // p[1] now ok,  [0, n-1]

        // compute I*(p[1] - p[3]),  [0, n-1]
        NMOD_MUL_PRENORM(tmp13, tmp13 << F->mod.norm, F->tab_w[0][1], F->mod);
        p[3] = nmod_sub(p[2], tmp13, F->mod); // p[3] now ok, [0, n-1]
        p[2] = nmod_add(p[2], tmp13, F->mod); // p[2] now ok, [0, n-1]
    }
    else
    {
        // here order >= 3, len >= 8, try some unrolling
        const mp_ptr phi = p+(len/2);
        mp_limb_t tmp[4];
        for (ulong k = 0; k+3 < len/2; k+=4)
        {
            tmp[0] = p[k+0];
            tmp[1] = p[k+1];
            tmp[2] = p[k+2];
            tmp[3] = p[k+3];
            p[k+0] = nmod_add(tmp[0], phi[k+0], F->mod);
            p[k+1] = nmod_add(tmp[1], phi[k+1], F->mod);
            p[k+2] = nmod_add(tmp[2], phi[k+2], F->mod);
            p[k+3] = nmod_add(tmp[3], phi[k+3], F->mod);
            phi[k+0] = nmod_sub(tmp[0], phi[k+0], F->mod);
            phi[k+1] = nmod_sub(tmp[1], phi[k+1], F->mod);
            phi[k+2] = nmod_sub(tmp[2], phi[k+2], F->mod);
            phi[k+3] = nmod_sub(tmp[3], phi[k+3], F->mod);
            phi[k+0] = n_mulmod_shoup(F->tab_w[order-2][k+0], phi[k+0], F->tab_w_pre[order-2][k+0], F->mod.n);
            phi[k+1] = n_mulmod_shoup(F->tab_w[order-2][k+1], phi[k+1], F->tab_w_pre[order-2][k+1], F->mod.n);
            phi[k+2] = n_mulmod_shoup(F->tab_w[order-2][k+2], phi[k+2], F->tab_w_pre[order-2][k+2], F->mod.n);
            phi[k+3] = n_mulmod_shoup(F->tab_w[order-2][k+3], phi[k+3], F->tab_w_pre[order-2][k+3], F->mod.n);
        }
        _nmod_poly_dif_inplace_radix2_rec_v4(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_v4(phi, len/2, order-1, F);
    }

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
