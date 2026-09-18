#![allow(unused)]

use crate::{
    params::*,
    poly::*,
    polyvec::*,
    packing::*,
    fips202::shake256,
};

/*
pub struct PublicParams<'a> {
    rho: &'a [u8],
    b_h: &'a [u8],
}

pub struct MasterKey<'a> {
    s: &'a [u8],
}

pub struct PartialPrivateKey<'a> {
    s:      &'a [u8],
    y:      &'a [u8],
    b_l:    &'a [u8],
}
*/

#[derive(Copy, Clone)]
pub struct Params {
    pub rho:   [u8; SEEDBYTES],
    pub b11_h: Polyveck,
}

#[derive(Copy, Clone)]
pub struct PartialPrivateKey {
    pub s11:   Polyvecl,
    pub y11:   Polyvecl,
    pub b11_l: Polyveck,
}

#[derive(Copy, Clone)]
pub struct MasterKey {
    pub s11: Polyvecl,
}

// 32-bytes ID size.
pub const ID_SIZE: usize = 32;

///
/// Input:
///     identity:   ID_A
///     rho:        security parameter
pub fn partial_private_key_generation(
    // params:     &mut [u8],
    // master_key: &mut [u8],
    // ppk:        &mut [u8],
    identity:   &[u8],
    rho:        &[u8],
) -> (Params, MasterKey, PartialPrivateKey) {
    assert_eq!(rho.len(), SEEDBYTES);
    assert_eq!(identity.len(), ID_SIZE);

    // Collision Resistant Hash function: shake256

    // Compute Seed: r and r_prime
    let mut seedbuf = [0u8; SEEDBYTES + CRHBYTES];
    // let mut r = [0u8; SEEDBYTES];

    let mut hash_input = [0u8; SEEDBYTES + ID_SIZE];
    hash_input[..SEEDBYTES].copy_from_slice(rho);
    hash_input[SEEDBYTES..].copy_from_slice(identity);

    // Write seed to seedbuf.
    shake256(&mut seedbuf, SEEDBYTES + CRHBYTES, &hash_input, SEEDBYTES + ID_SIZE);

    let mut r = [0u8; SEEDBYTES];
    let mut r_prime = [0u8; CRHBYTES];

    r.copy_from_slice(&seedbuf[..SEEDBYTES]);
    r_prime.copy_from_slice(&seedbuf[SEEDBYTES..]);

    // Public Matrix Generation: mat_a11, mat_a12

    let mut mat = [Polyvecl::default(); K];

    polyvec_matrix_expand(&mut mat, &r);

    let mut mat_a11 = mat;

    // Sample Secrets: s11, e11
    let mut s11 = Polyvecl::default();
    let mut e11 = Polyveck::default();

    // Sample short vectors s11 and e11.
    polyvecl_uniform_eta(&mut s11, &r_prime, 0);
    polyveck_uniform_eta(&mut e11, &r_prime, L_U16);

    let (mut b_h, mut b_l) = (
        Polyveck::default(),
        Polyveck::default(),
    );
    
    // Matrix-vector multiplication
    let mut s_hat = s11;
    polyvecl_ntt(&mut s_hat);
    polyvec_matrix_pointwise_montgomery(&mut b_h, &mat_a11, &s_hat);
    polyveck_reduce(&mut b_h);
    polyveck_invntt_tomont(&mut b_h);

    // Add error vector e11
    polyveck_add(&mut b_h, &e11);
    // Extract b_h
    polyveck_caddq(&mut b_h);
    polyveck_power2round(&mut b_h, &mut b_l);

    // Partial-private key generation.
    let mut y = Polyvecl::default();
    polyvecl_uniform_gamma1(&mut y, &r_prime, L_U16);

    /*
    let mut v = Polyveck::default();
    polyvecl_ntt(&mut y); // TODO: do not forget to invert y.
    polyvec_matrix_pointwise_montgomery(&mut v, &mat_a11, &y);
    polyveck_reduce(&mut v);
    polyveck_invntt_tomont(&mut v);
    */

    // Set/pack params, master_key and partial-private key.

    // Pack params.
    // pack_params(params, &rho, &b_h);
    let mut params = Params {
        rho: [0u8; SEEDBYTES],
        b11_h: b_h,
    };
    params.rho.copy_from_slice(&rho);
    
    // Pack master key.
    // pack_short_vecl(master_key, &s11);
    let msk = MasterKey {
        s11: s11,
    };

    // Pack partial private key
    // pack_ppk(ppk, &s11, &y, &b_l);
    let ppk = PartialPrivateKey {
        s11:   s11,
        y11:   y,
        b11_l: b_l
    };

    (params, msk, ppk)
}

// Inspired from pack_pk.
pub fn pack_params(params: &mut [u8], rho: &[u8], b_h: &Polyveck) {
    params[..SEEDBYTES].copy_from_slice(&rho[..SEEDBYTES]);
    for i in 0..K {
        polyt1_pack(&mut params[SEEDBYTES + i * POLYT1_PACKEDBYTES..], &b_h.vec[i]);
    }
}

// Inspired from unpack_pk.
pub fn unpack_params(rho: &mut [u8], b_h: &mut Polyveck, params: &[u8]) {
    rho[..SEEDBYTES].copy_from_slice(&params[..SEEDBYTES]);

    for i in 0..K {
        polyt1_unpack(&mut b_h.vec[i], &params[SEEDBYTES + i * POLYT1_PACKEDBYTES..]);
    }
}

pub fn pack_short_vecl(master_key: &mut [u8], s: &Polyvecl) {
    for i in 0..L {
        polyeta_pack(&mut master_key[i * POLYETA_PACKEDBYTES..], &s.vec[i]);
    }
}

// Inspired from pack_sk.
pub fn pack_ppk(ppk: &mut [u8], s: &Polyvecl, y: &Polyvecl, b_l: &Polyveck) {
    let mut idx = 0usize;

    for i in 0..L {
        polyeta_pack(&mut ppk[idx + i * POLYETA_PACKEDBYTES..], &s.vec[i]);
    }
    idx += L * POLYETA_PACKEDBYTES;

    for i in 0..L {
        polyz_pack(&mut ppk[idx + i * POLYZ_PACKEDBYTES..], &y.vec[i]);
    }
    idx += L * POLYZ_PACKEDBYTES;

    for i in 0..K {
        polyeta_pack(&mut ppk[idx + i * POLYETA_PACKEDBYTES..], &b_l.vec[i]);
    }
}

// Inspired from unpack_sk.
pub fn unpack_ppk(s: &mut Polyvecl, y: &mut Polyvecl, b_l: &mut Polyveck, ppk: &[u8]) {
    let mut idx = 0usize;

    for i in 0..L {
        polyeta_unpack(&mut s.vec[i], &ppk[idx + i * POLYETA_PACKEDBYTES..]);
    }
    idx += L * POLYETA_PACKEDBYTES;

    for i in 0..L {
        polyz_unpack(&mut y.vec[i], &ppk[idx + i * POLYZ_PACKEDBYTES..]);
    }
    idx += L * POLYZ_PACKEDBYTES;

    for i in 0..K {
        polyeta_unpack(&mut b_l.vec[i], &ppk[idx + i * POLYETA_PACKEDBYTES..]);
    }
}
