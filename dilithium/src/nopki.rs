#![allow(unused)]

use core::cell::Cell;

use crate::{
    params::*,
    poly::*,
    polyvec::*,
    packing::*,
    fips202::*,
    randombytes::*,
};

// Some slight differences between vector dimensions used by paper and code.
//
// paper -> code
// m     -> k
// l     -> n

// 32-bytes ID size.
pub const ID_SIZE: usize = CRHBYTES;

#[derive(Copy, Clone)]
pub struct Params {
    rho: [u8; SEEDBYTES],
}

pub struct Mpk {
    b11_h: Polyveck,
}

pub struct MasterKey {
    s11: Polyvecl,
}

pub struct PartialPrivateKey {
    e11: Polyveck,
    b11_l: Polyveck,
    user_secret: [u8; 32],
}

#[derive(Copy, Clone)]
pub struct UserKey {
    z1: Polyvecl,
    h1: Polyveck,
    c1_hashed: [u8; SEEDBYTES],
}

/// KGC Setup and Partial Private Key generation.
pub fn kgc_setup(identity: &[u8], rho: &[u8]) -> (Params, Mpk, MasterKey, 
    PartialPrivateKey, UserKey) {

    // Assert sizes on identity and rho.
    assert_eq!(rho.len(), SEEDBYTES);
    assert_eq!(identity.len(), ID_SIZE);

    // Compute r0, and r1.
    let mut seedbuf = [0u8; SEEDBYTES + CRHBYTES]; // 64 + 32.

    // This will be the input to hash function.
    // let mut hash_input = [0u8; SEEDBYTES + CRHBYTES];
    let mut hash_input = [0u8; SEEDBYTES];

    // Fill hash_input with rho and identity.
    hash_input[..SEEDBYTES].copy_from_slice(rho);
    // hash_input[SEEDBYTES..].copy_from_slice(identity);

    // Write hash to seedbuf.
    shake256(&mut seedbuf, SEEDBYTES + CRHBYTES, &hash_input, hash_input.len());

    let mut r0 = [0u8; SEEDBYTES]; // Used to sample matrix.
    let mut r1 = [0u8; CRHBYTES];  // Used during sampling secret vectors.

    r0.copy_from_slice(&seedbuf[..SEEDBYTES]);
    r1.copy_from_slice(&seedbuf[SEEDBYTES..]);

    // Public matrix generation.
    let mut mat_a11 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a11, &r0);

    // NTT transform the entire matrix, useful for later multiplication.
    let mut mat_a11_ntt = mat_a11;
    let rows = mat_a11_ntt.len();
    for i in 0..rows {
        polyvecl_ntt(&mut mat_a11_ntt[i]);
    }

    // Sample Secrets: s11, e11
    let mut s11 = Polyvecl::default();
    let mut e11 = Polyveck::default();
    
    // Sample short vectors s11 and e11.
    polyvecl_uniform_eta(&mut s11, &r1, 0);
    polyveck_uniform_eta(&mut e11, &r1, L_U16);

    // Create NTT form of s11 and e11.
    let mut s11_ntt = s11;
    polyvecl_ntt(&mut s11_ntt);
    
    let mut e11_ntt = e11;
    polyveck_ntt(&mut e11_ntt);

    // Compute b11.
    let mut b11 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut b11, &mat_a11_ntt, &s11_ntt);
    polyveck_invntt_tomont(&mut b11);
    polyveck_reduce(&mut b11);
    polyveck_add(&mut b11, &e11);

    // Power2round b11.
    let (mut b11_h, mut b11_l) = (Polyveck::default(), Polyveck::default());
    b11_h = b11;
    polyveck_caddq(&mut b11_h);
    polyveck_power2round(&mut b11_h, &mut b11_l);

    // NTT form of b11_l, will be used later to multiply.
    let mut b11_l_ntt = b11_l;
    polyveck_ntt(&mut b11_l_ntt);

    let mut iter_idx: usize = 0;

    // Pre-hash rho and identity to avoid repeated hashing of these two inside
    // the loop.
    let mut prehashed_rho_id = KeccakState::default();
    shake256_absorb(&mut prehashed_rho_id, &rho, rho.len());
    shake256_absorb(&mut prehashed_rho_id, &identity, identity.len());

    // Pre-hash r0 and identity to avoid repeated hashing inside the loop.
    let mut prehashed_r0_id = KeccakState::default();
    shake256_absorb(&mut prehashed_r0_id, &r0, r0.len());
    shake256_absorb(&mut prehashed_r0_id, &identity, identity.len());

    loop {
        // Size of user secret.
        const SECRET_SIZE: usize = 32; // in bytes, so 256-bits.

        // User secret.
        let mut user_secret = [0u8; SECRET_SIZE]; // 256-bits.


        // Sample a user secret based on the hash of iter_idx. Better randomness
        // could be used.
        shake256(&mut user_secret, SECRET_SIZE, &iter_idx.to_ne_bytes(),
            size_of::<usize>());

        // Complete hash of rho, identity and user_secret.
        let mut rho_id_user_secret_hash = [0u8; CRHBYTES];
        let mut hasher_state = prehashed_rho_id;
        shake256_absorb(&mut hasher_state, &user_secret, user_secret.len());
        shake256_finalize(&mut hasher_state);
        shake256_squeeze(&mut rho_id_user_secret_hash, CRHBYTES, &mut hasher_state);

        // Use the just finished hash to sample y11.
        let mut y11 = Polyvecl::default();
        polyvecl_uniform_gamma1(&mut y11, &rho_id_user_secret_hash, 0);

        // NTT form of y11.
        let mut y11_ntt = y11;
        polyvecl_ntt(&mut y11_ntt);

        // Compute v11.
        let mut v11 = Polyveck::default();
        polyvec_matrix_pointwise_montgomery(&mut v11, &mat_a11_ntt, &y11_ntt);
        polyveck_invntt_tomont(&mut v11);
        polyveck_reduce(&mut v11);

        // Decompose v11 into high and low.
        let (mut v11_h, mut v11_l) = (Polyveck::default(), Polyveck::default());
        v11_h = v11;
        polyveck_reduce(&mut v11_h);
        polyveck_caddq(&mut v11_h);
        polyveck_decompose(&mut v11_h, &mut v11_l);

        // Finish the hashing of r0, id, and v11_h.

        // Pack the v11_h bits.
        let mut v11_h_packed = [0u8; K * POLYW1_PACKEDBYTES];
        polyveck_pack_w1(v11_h_packed.as_mut_slice(), &v11_h);

        let mut hasher_state = prehashed_r0_id;
        shake256_absorb(&mut hasher_state, &v11_h_packed, v11_h_packed.len());
        shake256_finalize(&mut hasher_state);

        let mut c1_hashed = [0u8; SEEDBYTES];
        shake256_squeeze(&mut c1_hashed, SEEDBYTES, &mut hasher_state);

        // Sample in ball, using c1_hashed seed, immediately convert it to ntt
        // form.
        let mut c1_ntt = Poly::default();
        poly_challenge_nonced(&mut c1_ntt, &c1_hashed, 0);
        poly_ntt(&mut c1_ntt);

        // Compute z1.
        let mut z1 = Polyvecl::default();
        polyvecl_pointwise_poly_montgomery(&mut z1, &c1_ntt, &s11_ntt);
        polyvecl_invntt_tomont(&mut z1);
        polyvecl_reduce(&mut z1);
        polyvecl_add(&mut z1, &y11);

        // Compute r_poly, then decompose it into high and low.
        let mut r_poly = Polyveck::default();

        // Product of c1 and e11.
        let mut c1_e11 = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c1_e11, &c1_ntt, &e11_ntt);
        polyveck_invntt_tomont(&mut c1_e11);
        polyveck_reduce(&mut c1_e11);

        polyveck_sub(&mut r_poly, &c1_e11);
        polyveck_add(&mut r_poly, &v11);

        // Decompose r_poly into high and low.
        let (mut r_h, mut r_l) = (Polyveck::default(), Polyveck::default());
        r_h = r_poly;
        polyveck_reduce(&mut r_h);
        polyveck_caddq(&mut r_h);
        polyveck_decompose(&mut r_h, &mut r_l);

        // Check the following conditions before proceeding:
        // - z1 norm
        // - r_l norm
        // r_h == v11_h

        // Check z1 norm.
        if polyvecl_chknorm(&z1, (GAMMA1 - BETA) as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        // Check r_l norm.
        if polyveck_chknorm(&r_l, (GAMMA2 - BETA) as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        // Check r_h == v11_h
        let mut rh_v11h_diff = r_h;
        polyveck_sub(&mut rh_v11h_diff, &v11_h);
        polyveck_reduce(&mut rh_v11h_diff);

        // If they are same then their infinite norm must be less than 1.
        if polyveck_chknorm(&rh_v11h_diff, 1) > 0 {
            iter_idx += 1;
            continue;
        }

        // Compute c1_b11_l in order to later make hint.
        let mut c1_b11_l = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c1_b11_l, &c1_ntt, &b11_l_ntt);
        polyveck_invntt_tomont(&mut c1_b11_l);
        polyveck_reduce(&mut c1_b11_l);

        // Check norm of c1_b11_l.
        if polyveck_chknorm(&c1_b11_l, GAMMA2 as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        let mut neg_c1_b11_l = Polyveck::default();
        polyveck_sub(&mut neg_c1_b11_l, &c1_b11_l);

        let mut hint_from = v11;
        polyveck_sub(&mut hint_from, &c1_e11);
        polyveck_add(&mut hint_from, &c1_b11_l);
        polyveck_reduce(&mut hint_from);
        polyveck_caddq(&mut hint_from);

        // The hint.
        let mut h1 = Polyveck::default();

        // Make hint and get the number of 1's in it.
        let n = polyveck_make_hint_simple(&mut h1, &neg_c1_b11_l, &hint_from);

        if n > OMEGA as i32 {
            // println!("kgc: omega fail: n: {n}, omega: {OMEGA}");
            iter_idx += 1;
            continue;
        }

        // println!("kgc: loop ended after {} extra iterations!", iter_idx);

        use std::convert::TryInto;
        let params = Params {
            rho: rho.try_into().unwrap(),
        };

        let mpk = Mpk {
            b11_h: b11_h,
        };

        let msk = MasterKey {
            s11: s11,
        };

        let ppk = PartialPrivateKey {
            e11: e11,
            b11_l: b11_l,
            user_secret: user_secret,
        };
        
        let userkey = UserKey {
            z1: z1,
            h1: h1,
            c1_hashed: c1_hashed,
        };

        break (params, mpk, msk, ppk, userkey);
    }
}

pub struct SecretKey {
    b12_l: Polyveck,
    y11: Polyvecl,
    y12: Polyvecl,
    s12: Polyvecl,
    e11: Polyveck,
    e12: Polyveck,
}

#[derive(Copy, Clone)]
pub struct PublicKey {
    b12_h: Polyveck,
    z2: Polyvecl,
    h2: Polyveck,
    c2_hashed: [u8; SEEDBYTES],
    c3_hashed: [u8; SEEDBYTES],
}

/// User Key Generation.
pub fn user_keygen(
    params: Params,
    ppk: PartialPrivateKey,
    identity: &[u8],
    upk: UserKey,
    msk: MasterKey,
) -> (PublicKey, SecretKey) {

    // Unpack arguments.
    let Params { rho } = params;
    let PartialPrivateKey { e11, b11_l, user_secret } = ppk;
    let UserKey { z1, h1, c1_hashed } = upk;
    let MasterKey { s11 } = msk;

    // Some assertions about lengths.
    assert_eq!(rho.len(), SEEDBYTES);
    assert_eq!(identity.len(), ID_SIZE);

    // NTT forms for later use.
    let mut e11_ntt = e11;
    polyveck_ntt(&mut e11_ntt);

    let mut b11_l_ntt = b11_l;
    polyveck_ntt(&mut b11_l_ntt);

    // Compute r0, r1 <- CRH(rho)
    let mut seedbuf = [0u8; SEEDBYTES + CRHBYTES];
    shake256(&mut seedbuf, SEEDBYTES + CRHBYTES, &rho, rho.len());

    let mut r0 = [0u8; SEEDBYTES];
    let mut r1 = [0u8; CRHBYTES];

    r0.copy_from_slice(&seedbuf[..SEEDBYTES]);
    r1.copy_from_slice(&seedbuf[SEEDBYTES..]);

    // Public matrix generation.
    //
    // A11 <- Expand(r0)
    // A12 <- Expand(r1)
    let mut mat_a11 = [Polyvecl::default(); K];
    let mut mat_a12 = [Polyvecl::default(); K];

    // Use only first 32-bytes for expansion.
    polyvec_matrix_expand(&mut mat_a11, &r0[..SEEDBYTES]);
    polyvec_matrix_expand(&mut mat_a12, &r1[..SEEDBYTES]);

    // NTT transform matrices for later multiplication.
    let mut mat_a11_ntt = mat_a11;
    let mut mat_a12_ntt = mat_a12;

    for i in 0..K {
        polyvecl_ntt(&mut mat_a11_ntt[i]);
        polyvecl_ntt(&mut mat_a12_ntt[i]);
    }

    // Sample s12, e12
    let mut s12 = Polyvecl::default();
    let mut e12 = Polyveck::default();

    polyvecl_uniform_eta(&mut s12, &r1, 0);
    polyveck_uniform_eta(&mut e12, &r1, L_U16);

    // NTT forms of s12, e12.
    let mut s12_ntt = s12;
    polyvecl_ntt(&mut s12_ntt);

    let mut e12_ntt = e12;
    polyveck_ntt(&mut e12_ntt);

    // Compute b12 = A12*s12 + e12
    let mut b12 = Polyveck::default();

    polyvec_matrix_pointwise_montgomery(&mut b12, &mat_a12_ntt, &s12_ntt);
    polyveck_invntt_tomont(&mut b12);
    polyveck_reduce(&mut b12);
    polyveck_add(&mut b12, &e12);

    // Power2Round.
    let (mut b12_h, mut b12_l) = (Polyveck::default(), Polyveck::default());
    b12_h = b12;
    polyveck_caddq(&mut b12_h);
    polyveck_power2round(&mut b12_h, &mut b12_l);

    // NTT form of b12_l.
    let mut b12_l_ntt = b12_l;
    polyveck_ntt(&mut b12_l_ntt);

    // Sample y11 <- S^(n)_{gamma1-1}
    //
    // y11 <- CRH(rho || ID || user_secret)
    let mut hasher_state = KeccakState::default();

    shake256_absorb(&mut hasher_state, &rho, rho.len());
    shake256_absorb(&mut hasher_state, identity, identity.len());
    shake256_absorb(&mut hasher_state, &user_secret, user_secret.len());
    shake256_finalize(&mut hasher_state);

    let mut rho_id_user_secret_hash = [0u8; CRHBYTES];

    shake256_squeeze(&mut rho_id_user_secret_hash, CRHBYTES, &mut hasher_state);

    let mut y11 = Polyvecl::default();
    polyvecl_uniform_gamma1(&mut y11, &rho_id_user_secret_hash, 0);

    // Compute v11 = A11*y11
    let mut y11_ntt = y11;
    polyvecl_ntt(&mut y11_ntt);

    let mut v11 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut v11, &mat_a11_ntt, &y11_ntt);
    polyveck_invntt_tomont(&mut v11);
    polyveck_reduce(&mut v11);

    // v11_h = HighBits(v11, 2*gamma2)
    let (mut v11_h, mut v11_l) = (Polyveck::default(), Polyveck::default());
    v11_h = v11;
    polyveck_reduce(&mut v11_h);
    polyveck_caddq(&mut v11_h);
    polyveck_decompose(&mut v11_h, &mut v11_l);

    // c1 from the upk.
    let mut c1_ntt = Poly::default();
    poly_challenge_nonced(&mut c1_ntt, &c1_hashed, 0);
    poly_ntt(&mut c1_ntt);

    let mut iter_idx: usize = 0;

    loop {
        // Sample y12 <- S^(n)_{gamma1-1}
        let mut y12 = Polyvecl::default();
        polyvecl_uniform_gamma1(&mut y12, &r1, iter_idx as u16);

        // NTT form.
        let mut y12_ntt = y12;
        polyvecl_ntt(&mut y12_ntt);

        // v12 = A12*y12
        let mut v12 = Polyveck::default();
        polyvec_matrix_pointwise_montgomery(&mut v12, &mat_a12_ntt, &y12_ntt);
        polyveck_invntt_tomont(&mut v12);
        polyveck_reduce(&mut v12);

        // v12_h = HighBits(v12, 2*gamma2)
        let (mut v12_h, mut v12_l) = (Polyveck::default(), Polyveck::default());
        v12_h = v12;
        polyveck_reduce(&mut v12_h);
        polyveck_caddq(&mut v12_h);
        polyveck_decompose(&mut v12_h, &mut v12_l);

        // v1 = v11 + v12
        let mut v1 = v11;
        polyveck_add(&mut v1, &v12);
        polyveck_reduce(&mut v1);
        polyveck_caddq(&mut v1);

        // Decompose v1.
        let (mut v1_h, mut v1_l) = (Polyveck::default(), Polyveck::default());
        v1_h = v1;
        polyveck_decompose(&mut v1_h, &mut v1_l);

        // c2 = CRH(r1 || ID || v12_h)

        // Pack v12_h.
        let mut v12_h_packed = [0u8; K * POLYW1_PACKEDBYTES];

        polyveck_pack_w1(v12_h_packed.as_mut_slice(), &v12_h);

        // Hash r1, identity and v12_h_packed
        let mut hasher_state = KeccakState::default();

        shake256_absorb(&mut hasher_state, &r1, r1.len());
        shake256_absorb(&mut hasher_state, identity, identity.len());
        shake256_absorb(&mut hasher_state, &v12_h_packed, v12_h_packed.len());
        shake256_finalize(&mut hasher_state);

        let mut c2_hashed = [0u8; SEEDBYTES];
        shake256_squeeze(&mut c2_hashed, SEEDBYTES, &mut hasher_state);

        // Challenge polynomial c2.
        let mut c2_ntt = Poly::default();
        poly_challenge_nonced(&mut c2_ntt, &c2_hashed, 0);
        poly_ntt(&mut c2_ntt);

        // z2 = y12 + s12*c2
        let mut z2 = Polyvecl::default();
        polyvecl_pointwise_poly_montgomery(&mut z2, &c2_ntt, &s12_ntt);
        polyvecl_invntt_tomont(&mut z2);
        polyvecl_reduce(&mut z2);
        polyvecl_add(&mut z2, &y12);

        // c1*e11
        let mut c1_e11 = Polyveck::default();

        polyveck_pointwise_poly_montgomery(&mut c1_e11, &c1_ntt, &e11_ntt);
        polyveck_invntt_tomont(&mut c1_e11);
        polyveck_reduce(&mut c1_e11);

        // c2*e12
        let mut c2_e12 = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c2_e12, &c2_ntt, &e12_ntt);
        polyveck_invntt_tomont(&mut c2_e12);
        polyveck_reduce(&mut c2_e12);

        // r_poly = v1 - c2*e12 - c1*e11
        let mut r_poly = v1;
        polyveck_sub(&mut r_poly, &c2_e12);
        polyveck_sub(&mut r_poly, &c1_e11);

        // Decompose r1.
        let (mut r_h, mut r_l) = (Polyveck::default(), Polyveck::default());
        r_h = r_poly;
        polyveck_reduce(&mut r_h);
        polyveck_caddq(&mut r_h);
        polyveck_decompose(&mut r_h, &mut r_l);

        // Check:
        //
        // ||z2|| >= gamma1-beta
        // ||r_l|| >= gamma2-beta
        // r_h != v1_h
        if polyvecl_chknorm(&z2, (GAMMA1 - BETA) as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        if polyveck_chknorm(&r_l, (GAMMA2 - BETA) as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        let mut rh_v1h_diff = r_h;
        polyveck_sub(&mut rh_v1h_diff, &v1_h);
        polyveck_reduce(&mut rh_v1h_diff);

        if polyveck_chknorm(&rh_v1h_diff, 1) > 0 {
            iter_idx += 1;
            continue;
        }

        // c1*b11_l
        let mut c1_b11_l = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c1_b11_l, &c1_ntt, &b11_l_ntt);
        polyveck_invntt_tomont(&mut c1_b11_l);
        polyveck_reduce(&mut c1_b11_l);

        // c2*b12_l
        let mut c2_b12_l = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c2_b12_l, &c2_ntt, &b12_l_ntt);
        polyveck_invntt_tomont(&mut c2_b12_l);
        polyveck_reduce(&mut c2_b12_l);

        // c1*b11_l + c2*b12_l
        let mut c_b_l = c1_b11_l;
        polyveck_add(&mut c_b_l, &c2_b12_l);
        polyveck_reduce(&mut c_b_l);

        // Check: ||c1*b11_l + c2*b12_l|| >= gamma2
        if polyveck_chknorm(&c_b_l, GAMMA2 as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        // Compute: v1 - (c1*e11 + c2*e12) + (c1*b11_l + c2*b12_l)
        let mut hint_from = v1;

        let mut c_e = c1_e11;
        polyveck_add(&mut c_e, &c2_e12);
        polyveck_reduce(&mut c_e);

        polyveck_sub(&mut hint_from, &c_e);
        polyveck_add(&mut hint_from, &c_b_l);
        polyveck_reduce(&mut hint_from);
        polyveck_caddq(&mut hint_from);

        // h2 = MakeHint(-(c1*b11_l + c2*b12_l), hint_from)
        let mut neg_c_b_l = Polyveck::default();
        polyveck_sub(&mut neg_c_b_l, &c_b_l);

        let mut h2 = Polyveck::default();

        let n = polyveck_make_hint_simple(&mut h2, &neg_c_b_l, &hint_from);

        if n > OMEGA as i32 {
            // println!("userkeygen: omega fail: n: {n}, omega: {OMEGA}");
            iter_idx += 1;
            continue;
        }

        // c3 = CRH(c1 || c2 || v1_h)

        // Reduce v1_h before packing and hashing.
        polyveck_reduce(&mut v1_h);
        polyveck_caddq(&mut v1_h);
        let mut v1_h_packed = [0u8; K * POLYW1_PACKEDBYTES];
        polyveck_pack_w1(v1_h_packed.as_mut_slice(), &v1_h);

        let mut hasher_state = KeccakState::default();

        shake256_absorb(&mut hasher_state, &c1_hashed, c1_hashed.len());
        shake256_absorb(&mut hasher_state, &c2_hashed, c2_hashed.len());
        shake256_absorb(&mut hasher_state, &v1_h_packed, v1_h_packed.len());
        shake256_finalize(&mut hasher_state);

        let mut c3_hashed = [0u8; SEEDBYTES];
        shake256_squeeze(&mut c3_hashed, SEEDBYTES, &mut hasher_state);

        let sk = SecretKey {
            b12_l,
            y11,
            y12,
            s12,
            e11,
            e12,
        };

        let pk = PublicKey {
            b12_h,
            z2,
            h2,
            c2_hashed,
            c3_hashed,
        };

        // println!("user keygen: loop ended after {} iterations!", iter_idx);

        break (pk, sk);
    }
}

pub struct Signature {
    z: Polyvecl,
    h: Polyveck,
    c_i_hashed: [u8; SEEDBYTES],
}

/// Sign message m_i.
pub fn sign(params: Params, pk: PublicKey, sk: SecretKey, identity: &[u8], message: &[u8])
    -> Signature {

    // Unpack arguments.
    let Params { rho } = params;
    let PublicKey { b12_h, z2, h2, c2_hashed, c3_hashed } = pk;
    let SecretKey { b12_l, y11, y12, s12, e11, e12 } = sk;

    // Some assertions about lengths.
    assert_eq!(rho.len(), SEEDBYTES);
    assert_eq!(identity.len(), ID_SIZE);

    // NTT forms for later use.
    let mut b12_l_ntt = b12_l;
    polyveck_ntt(&mut b12_l_ntt);

    // Compute r0, r1 <- CRH(rho)
    let mut seedbuf = [0u8; SEEDBYTES + CRHBYTES];
    shake256(&mut seedbuf, SEEDBYTES + CRHBYTES, &rho, rho.len());

    let mut r0 = [0u8; SEEDBYTES];
    let mut r1 = [0u8; CRHBYTES];

    r0.copy_from_slice(&seedbuf[..SEEDBYTES]);
    r1.copy_from_slice(&seedbuf[SEEDBYTES..]);

    // Public matrix generation.
    //
    // A12 <- Expand(r1)
    let mut mat_a12 = [Polyvecl::default(); K];

    // Use only first 32-bytes for expansion.
    polyvec_matrix_expand(&mut mat_a12, &r1[..SEEDBYTES]);

    // NTT transform matrix for multiplication.
    let mut mat_a12_ntt = mat_a12;
    for i in 0..K {
        polyvecl_ntt(&mut mat_a12_ntt[i]);
    }

    // NTT forms of s12, e12.
    let mut s12_ntt = s12;
    polyvecl_ntt(&mut s12_ntt);

    let mut e12_ntt = e12;
    polyveck_ntt(&mut e12_ntt);

    let mut iter_idx: usize = 0;

    loop {

        // Sample y_i <- S^(n)_{gamma1-1}
        let mut y_i = Polyvecl::default();

        polyvecl_uniform_gamma1(&mut y_i, &r1, iter_idx as u16);

        // NTT form for matrix multiplication.
        let mut y_i_ntt = y_i;
        polyvecl_ntt(&mut y_i_ntt);

        // v_i = A12*y_i
        let mut v_i = Polyveck::default();

        polyvec_matrix_pointwise_montgomery(&mut v_i, &mat_a12_ntt, &y_i_ntt);
        polyveck_invntt_tomont(&mut v_i);
        polyveck_reduce(&mut v_i);

        // v_i^h = HighBits(v_i, 2*gamma2)
        let (mut v_i_h, mut v_i_l) = (Polyveck::default(), Polyveck::default());
        v_i_h = v_i;
        polyveck_reduce(&mut v_i_h);
        polyveck_caddq(&mut v_i_h);
        polyveck_decompose(&mut v_i_h, &mut v_i_l);

        // c_i := H(ID || v_i^h || m_i)
        //
        // Pack v_i^h before hashing.
        let mut v_i_h_packed = [0u8; K * POLYW1_PACKEDBYTES];

        polyveck_pack_w1(v_i_h_packed.as_mut_slice(), &v_i_h);

        let mut hasher_state = KeccakState::default();

        shake256_absorb(&mut hasher_state, identity, identity.len());
        shake256_absorb(&mut hasher_state, &v_i_h_packed, v_i_h_packed.len(),);
        shake256_absorb(&mut hasher_state, message, message.len());
        shake256_finalize(&mut hasher_state);

        let mut c_i_hashed = [0u8; SEEDBYTES];

        shake256_squeeze(&mut c_i_hashed, SEEDBYTES, &mut hasher_state);

        // c_i in B_60 := F(c_i)
        let mut c_i_ntt = Poly::default();

        poly_challenge_nonced(&mut c_i_ntt, &c_i_hashed, 0);
        poly_ntt(&mut c_i_ntt);

        // z_i = y_i + s12*c_i
        let mut z_i = Polyvecl::default();
        polyvecl_pointwise_poly_montgomery(&mut z_i, &c_i_ntt, &s12_ntt);
        polyvecl_invntt_tomont(&mut z_i);
        polyvecl_reduce(&mut z_i);
        polyvecl_add(&mut z_i, &y_i);

        // r_i = v_i - c_i*e12
        let mut c_i_e12 = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c_i_e12, &c_i_ntt, &e12_ntt);
        polyveck_invntt_tomont(&mut c_i_e12);
        polyveck_reduce(&mut c_i_e12);

        let mut r_i = v_i;
        polyveck_sub(&mut r_i, &c_i_e12);

        // (r_i^h, r_i^l) = Decompose_q(r_i, 2*gamma2)
        let (mut r_i_h, mut r_i_l) = (Polyveck::default(), Polyveck::default());
        r_i_h = r_i;
        polyveck_reduce(&mut r_i_h);
        polyveck_caddq(&mut r_i_h);
        polyveck_decompose(&mut r_i_h, &mut r_i_l);

        // Check:
        //
        // ||z_i|| >= gamma1 - beta
        // ||r_i^l|| >= gamma2 - beta
        // r_i^h != v_i^h

        if polyvecl_chknorm(&z_i, (GAMMA1 - BETA) as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        if polyveck_chknorm(&r_i_l, (GAMMA2 - BETA) as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        let mut ri_h_vi_h_diff = r_i_h;

        polyveck_sub(&mut ri_h_vi_h_diff, &v_i_h);

        polyveck_reduce(&mut ri_h_vi_h_diff);

        if polyveck_chknorm(&ri_h_vi_h_diff, 1) > 0 {
            iter_idx += 1;
            continue;
        }

        // h_i = MakeHint_q(-c_i*b12_l, v_i - c_i*e12 + c_i*b12_l, 2*gamma2)

        // c_i*b12_l
        let mut c_i_b12_l = Polyveck::default();

        polyveck_pointwise_poly_montgomery(&mut c_i_b12_l, &c_i_ntt, &b12_l_ntt);
        polyveck_invntt_tomont(&mut c_i_b12_l);
        polyveck_reduce(&mut c_i_b12_l);

        // Check ||c_i*b12_l|| < gamma2.
        if polyveck_chknorm(&c_i_b12_l, GAMMA2 as i32) > 0 {
            iter_idx += 1;
            continue;
        }

        // v_i - c_i*e12 + c_i*b12_l
        let mut hint_from = v_i;

        polyveck_sub(&mut hint_from, &c_i_e12);
        polyveck_add(&mut hint_from, &c_i_b12_l);
        polyveck_reduce(&mut hint_from);
        polyveck_caddq(&mut hint_from);

        // -c_i*b12_l
        let mut neg_c_i_b12_l = Polyveck::default();

        polyveck_sub(&mut neg_c_i_b12_l, &c_i_b12_l);

        // Make hint.
        let mut h_i = Polyveck::default();

        let n = polyveck_make_hint_simple(&mut h_i, &neg_c_i_b12_l, &hint_from);

        // Check number of hint bits.
        if n > OMEGA as i32 {
            // println!("sign: omega fail: n: {n}, omega: {OMEGA}");
            iter_idx += 1;
            continue;
        }

        // Signature found.
        // println!("sign: loop ended after {} iterations!", iter_idx);

        break Signature {
            z: z_i,
            h: h_i,
            c_i_hashed,
        };
    }
}

pub fn verify(
    params: Params,
    mpk: Mpk,
    pk: PublicKey,
    upk: UserKey,
    identity: &[u8],
    message: &[u8],
    signature: Signature,
) -> bool {

    // Unpack arguments.
    let Params { rho } = params;
    let Mpk { b11_h } = mpk;
    let PublicKey { b12_h, z2, h2, c2_hashed, c3_hashed } = pk;
    let UserKey { z1, h1, c1_hashed } = upk;
    let Signature { z: z_i, h: h_i, c_i_hashed } = signature;

    // Some assertions about lengths.
    assert_eq!(rho.len(), SEEDBYTES);
    assert_eq!(identity.len(), ID_SIZE);

    // Compute r0, r1 <- CRH(rho)
    let mut seedbuf = [0u8; SEEDBYTES + CRHBYTES];
    shake256(&mut seedbuf, SEEDBYTES + CRHBYTES, &rho, rho.len());

    let mut r0 = [0u8; SEEDBYTES];
    let mut r1 = [0u8; CRHBYTES];

    r0.copy_from_slice(&seedbuf[..SEEDBYTES]);
    r1.copy_from_slice(&seedbuf[SEEDBYTES..]);

    // Public matrix generation.
    //
    // A11 <- Expand(r0)
    // A12 <- Expand(r1)
    let mut mat_a11 = [Polyvecl::default(); K];
    let mut mat_a12 = [Polyvecl::default(); K];

    polyvec_matrix_expand(&mut mat_a11, &r0[..SEEDBYTES]);
    polyvec_matrix_expand(&mut mat_a12, &r1[..SEEDBYTES]);

    // NTT forms for matrix-vector multiplication.
    let mut mat_a11_ntt = mat_a11;
    let mut mat_a12_ntt = mat_a12;

    for i in 0..K {
        polyvecl_ntt(&mut mat_a11_ntt[i]);
        polyvecl_ntt(&mut mat_a12_ntt[i]);
    }

    // c1, c2, c_i in B_60 := F(c)

    // c1 = F(c1_hashed)
    let mut c1_ntt = Poly::default();
    poly_challenge_nonced(&mut c1_ntt, &c1_hashed, 0);
    poly_ntt(&mut c1_ntt);

    // c2 = F(c2_hashed)
    let mut c2_ntt = Poly::default();
    poly_challenge_nonced(&mut c2_ntt, &c2_hashed, 0);
    poly_ntt(&mut c2_ntt);

    // c_i = F(c_i_hashed)
    let mut c_i_ntt = Poly::default();
    poly_challenge_nonced(&mut c_i_ntt, &c_i_hashed, 0);
    poly_ntt(&mut c_i_ntt);

    // Compute:
    //
    // v1_prime = UseHint_q(h2, A11*z1 + A12*z2 - (c1*b11_h + c2*b12_h)*2^d, 2*gamma2)

    // A11 * z1
    let mut a11_z1 = Polyveck::default();

    let mut z1_ntt = z1;
    polyvecl_ntt(&mut z1_ntt);

    polyvec_matrix_pointwise_montgomery(&mut a11_z1, &mat_a11_ntt, &z1_ntt);
    polyveck_invntt_tomont(&mut a11_z1);
    polyveck_reduce(&mut a11_z1);

    // A12 * z2
    let mut a12_z2 = Polyveck::default();

    let mut z2_ntt = z2;
    polyvecl_ntt(&mut z2_ntt);

    polyvec_matrix_pointwise_montgomery(&mut a12_z2, &mat_a12_ntt, &z2_ntt);
    polyveck_invntt_tomont(&mut a12_z2);
    polyveck_reduce(&mut a12_z2);

    // A11*z1 + A12*z2
    let mut a_z = a11_z1;
    polyveck_add(&mut a_z, &a12_z2);
    polyveck_reduce(&mut a_z);

    // c1 * b11_h

    // b11_h is in normal representation here. Convert to NTT.
    let mut b11_h_ntt = b11_h;
    polyveck_ntt(&mut b11_h_ntt);

    let mut c1_b11_h = Polyveck::default();
    polyveck_pointwise_poly_montgomery(&mut c1_b11_h, &c1_ntt, &b11_h_ntt);
    polyveck_invntt_tomont(&mut c1_b11_h);
    polyveck_reduce(&mut c1_b11_h);

    // c2 * b12_h

    let mut b12_h_ntt = b12_h;
    polyveck_ntt(&mut b12_h_ntt);

    let mut c2_b12_h = Polyveck::default();
    polyveck_pointwise_poly_montgomery(&mut c2_b12_h, &c2_ntt, &b12_h_ntt);
    polyveck_invntt_tomont(&mut c2_b12_h);
    polyveck_reduce(&mut c2_b12_h);

    // c1*b11_h + c2*b12_h
    let mut c_b_h = c1_b11_h;
    polyveck_add(&mut c_b_h, &c2_b12_h);
    polyveck_reduce(&mut c_b_h);

    // Multiply by 2^d.
    let mut c_b_h_shifted = c_b_h;
    polyveck_shiftl(&mut c_b_h_shifted);
    polyveck_reduce(&mut c_b_h_shifted);

    // A11*z1 + A12*z2 - (c1*b11_h + c2*b12_h)*2^d
    let mut prehint_1 = a_z;
    polyveck_sub(&mut prehint_1, &c_b_h_shifted);
    polyveck_reduce(&mut prehint_1);
    polyveck_caddq(&mut prehint_1);

    // v1_prime = UseHint_q(h2, ..., 2*gamma2)
    let mut v1_prime = Polyveck::default();
    polyveck_use_hint_simple(&mut v1_prime, &h2, &prehint_1);

    // Check:
    //
    // ||z1||_inf < gamma1 - beta
    // ||z2||_inf < gamma1 - beta
    // c3 == H(c1 || c2 || v1')
    // #1's in h1, h2 <= omega

    if polyvecl_chknorm(&z1, (GAMMA1 - BETA) as i32) > 0 {
        // println!("z1 norm check failed");
        return false;
    }

    if polyvecl_chknorm(&z2, (GAMMA1 - BETA) as i32) > 0 {
        // println!("z2 norm check failed");
        return false;
    }

    // Check number of hint bits in h1.
    //
    // Assuming this to be true, as we just created the hints in above
    // functions.

    // c3 = H(c1 || c2 || v1_prime)

    // Reduce v1_prime before packing and hashing.
    polyveck_reduce(&mut v1_prime);
    polyveck_caddq(&mut v1_prime);

    let mut v1_prime_packed = [0u8; K * POLYW1_PACKEDBYTES];
    polyveck_pack_w1(v1_prime_packed.as_mut_slice(), &v1_prime);

    let mut hasher_state = KeccakState::default();

    shake256_absorb(&mut hasher_state, &c1_hashed, c1_hashed.len());
    shake256_absorb(&mut hasher_state, &c2_hashed, c2_hashed.len());
    shake256_absorb(&mut hasher_state, &v1_prime_packed, v1_prime_packed.len());
    shake256_finalize(&mut hasher_state);

    let mut computed_c3_hashed = [0u8; SEEDBYTES];

    shake256_squeeze(&mut computed_c3_hashed, SEEDBYTES, &mut hasher_state);

    // Compare c3 from public key against recomputed computed_c3_hashed.
    if c3_hashed != computed_c3_hashed {
        // println!("final hash(c3) did not matched");
        return false;
    }

    // v_i_prime = UseHint_q(h_i, A12*z_i - c_i*b12_h*2^d, 2*gamma2)

    // A12 * z_i
    let mut a12_zi = Polyveck::default();

    let mut z_i_ntt = z_i;
    polyvecl_ntt(&mut z_i_ntt);

    polyvec_matrix_pointwise_montgomery(&mut a12_zi, &mat_a12_ntt, &z_i_ntt);
    polyveck_invntt_tomont(&mut a12_zi);
    polyveck_reduce(&mut a12_zi);

    // c_i * b12_h
    let mut c_i_b12_h = Polyveck::default();

    polyveck_pointwise_poly_montgomery(&mut c_i_b12_h, &c_i_ntt, &b12_h_ntt);
    polyveck_invntt_tomont(&mut c_i_b12_h);
    polyveck_reduce(&mut c_i_b12_h);

    // Multiply c_i*b12_h by 2^d.
    let mut c_i_b12_h_shifted = c_i_b12_h;
    polyveck_shiftl(&mut c_i_b12_h_shifted);
    polyveck_reduce(&mut c_i_b12_h_shifted);

    // A12*z_i - c_i*b12_h*2^d
    let mut prehint_i = a12_zi;

    polyveck_sub( &mut prehint_i, &c_i_b12_h_shifted,);
    polyveck_reduce(&mut prehint_i);
    polyveck_caddq(&mut prehint_i);

    // v_i_prime = UseHint_q(h_i, ..., 2*gamma2)
    let mut v_i_prime = Polyveck::default();
    polyveck_use_hint_simple(&mut v_i_prime, &h_i, &prehint_i);

    // Check:
    //
    // ||z_i||_inf < gamma1 - beta
    // c_i == H(ID || v_i' || m_i)
    // #1's in h_i <= omega

    if polyvecl_chknorm(&z_i, (GAMMA1 - BETA) as i32) > 0 {
        // println!("z_i norm check failed");
        return false;
    }

    // #1's in h_i < omega
    // Assuming this to be true, as we just created the hints in above.

    // Pack v_i_prime for hashing.
    let mut v_i_prime_packed = [0u8; K * POLYW1_PACKEDBYTES];

    polyveck_pack_w1(v_i_prime_packed.as_mut_slice(), &v_i_prime);

    // c_i = H(ID || v_i_prime || m_i)
    let mut hasher_state = KeccakState::default();

    shake256_absorb(&mut hasher_state, identity, identity.len());
    shake256_absorb(&mut hasher_state, &v_i_prime_packed, v_i_prime_packed.len());
    shake256_absorb(&mut hasher_state, message, message.len());
    shake256_finalize(&mut hasher_state);

    let mut computed_c_i_hashed = [0u8; SEEDBYTES];
    shake256_squeeze(&mut computed_c_i_hashed, SEEDBYTES, &mut hasher_state);

    // Compare signature challenge against recomputed challenge.
    if computed_c_i_hashed != c_i_hashed {
        // println!("final hash(c_i) did not matched");
        return false;
    }

    // All verification conditions passed.
    true
}
