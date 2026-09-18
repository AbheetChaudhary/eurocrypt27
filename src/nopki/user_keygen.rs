#![allow(unused)]

pub const POLY_V_H_PACKEDBYTES: usize = POLYW1_PACKEDBYTES;

pub const SIG_RETRY_THRESHOLD: usize = 100;
pub const KEYGEN_RETRY_THRESHOLD: usize = 20;

const R_SIZE: usize = 2 * SEEDBYTES;
const R_PRIME_SIZE: usize = CRHBYTES;

#[derive(Debug, Default)]
struct FailReason {
    z_i1:      usize,
    z_i2:      usize,
    r_i_l:   usize,
    rih_vih: usize,
    cb:      usize,
    omega:   usize,
}

use crate::{
    params::*,
    poly::*,
    polyvec::*,
    packing::*,
    fips202::*,
    randombytes::*,

    nopki::kgc::*,
};

#[derive(Copy, Clone)]
pub struct PublicKey {
    pub b12_h: Polyveck,
    pub z1:    Polyvecl,
    pub z2:    Polyvecl,
    pub h:     Polyveck,
    pub c:     [u8; SEEDBYTES],
}

#[derive(Copy, Clone)]
pub struct SecretKey {
    pub b11_l:     Polyveck,
    pub b12_l:     Polyveck,
    pub y11:       Polyvecl,
    pub y12:       Polyvecl,
    pub s11:       Polyvecl,
    pub s12_prime: Polyvecl,
}

#[derive(Debug)]
pub enum KeygenError {
    ThresholdError,
}

pub fn user_generate_key(
    identity: &[u8],
    params:   Params,
    ppk:      PartialPrivateKey,
    seed:     Option<&[u8]>,
) -> Result<(PublicKey, SecretKey), KeygenError> {
    assert!(seed.is_none(), "No need to use this while benchmarking");

    /*
    // Buffer to contain seed.
    let mut init_seed = [0u8; SEEDBYTES];

    // If seed is Some then take that as the seed otherwise generate new random
    // seed.
    if let Some(seed) = seed {
        init_seed.copy_from_slice(seed);
    } else {
        randombytes(&mut init_seed, SEEDBYTES);
    }
    */

    // Unpack params.
    /*
    let mut rho = [0u8; SEEDBYTES];
    let mut b11_h = Polyveck::default();
    unpack_params(&mut rho, &mut b11_h, params);
    */

    let Params {mut rho, mut b11_h} = params;

    debug_assert_eq!(rho.len(), SEEDBYTES);
    debug_assert_eq!(identity.len(), ID_SIZE);

    // Unpack ppk into s11, y11, b11_l.
    /*
    let mut s11 = Polyvecl::default();
    let mut y11 = Polyvecl::default();
    let mut b11_l = Polyveck::default();
    unpack_ppk(&mut s11, &mut y11, &mut b11_l, ppk);
    */
    let PartialPrivateKey {mut s11, mut y11, mut b11_l} = ppk;

    // NTT representation of s11, y11, and b11_l
    let mut s11_hat = s11;
    polyvecl_ntt(&mut s11_hat);
    let s11_hat = s11_hat;

    let mut y11_hat = y11;
    polyvecl_ntt(&mut y11_hat);
    let y11_hat = y11_hat;

    let mut b11_l_hat = b11_l;
    polyveck_ntt(&mut b11_l_hat);
    let b11_l_hat = b11_l_hat;

    // Compute seed r and r_prime
    let mut seedbuf = [0u8; R_SIZE + R_PRIME_SIZE]; // r || r_prime

    let mut hash_input = [0u8; SEEDBYTES + ID_SIZE]; // rho || ID_A
    hash_input[..SEEDBYTES].copy_from_slice(&rho);
    hash_input[SEEDBYTES..].copy_from_slice(&identity);

    // Write seed to seedbuf.
    shake256(&mut seedbuf, R_SIZE + R_PRIME_SIZE, &hash_input, SEEDBYTES + ID_SIZE);

    let mut r = [0u8; R_SIZE];
    let mut r_prime = [0u8; R_PRIME_SIZE];

    r.copy_from_slice(&seedbuf[..R_SIZE]);
    r_prime.copy_from_slice(&seedbuf[R_SIZE..]);

    // Seed Expansion: mat_a11, mat_a12.

    // Expand seed for both matrices.
    // 
    // Use first SEEDBYTES bytes of r for first matrix.
    let mut mat_a11 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a11, &r[..SEEDBYTES]);

    // Use later SEEDBYTES bytes of r for second matrix.
    let mut mat_a12 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a12, &r[SEEDBYTES..]);

    // Secret Sampling.

    // Generate s12_prime.
    let mut s12_prime = Polyvecl::default();
    polyvecl_uniform_eta(&mut s12_prime, &r_prime, 0);

    // NTT form of s12_prime
    let mut s12_prime_hat = s12_prime;
    polyvecl_ntt(&mut s12_prime_hat);

    // Public Matrix Generation: b11 := mat_a12 x s12_prime + e12

    // Sample e12.
    let mut e12 = Polyveck::default();
    polyveck_uniform_eta(&mut e12, &r_prime, L_U16);

    let mut e12_hat = e12;
    polyveck_ntt(&mut e12_hat); // NTT form of e12

    let (mut b12_h, mut b12_l) = (
        Polyveck::default(),
        Polyveck::default(),
    );

    // Matrix-vector multiplication
    polyvec_matrix_pointwise_montgomery(&mut b12_h, &mat_a12, &s12_prime_hat);
    polyveck_reduce(&mut b12_h);
    polyveck_invntt_tomont(&mut b12_h);

    // Add error vector e12
    polyveck_add(&mut b12_h, &e12);

    // Extract b12_h
    polyveck_reduce(&mut b12_h);
    polyveck_caddq(&mut b12_h);
    polyveck_power2round(&mut b12_h, &mut b12_l);

    // NTT representation of b12_l for later use.
    let mut b12_l_hat = b12_l;
    polyveck_ntt(&mut b12_l_hat);

    // Compute v11_h = A11 x y11
    let (mut v11_h, mut v11_l) = (Polyveck::default(), Polyveck::default());
    polyvec_matrix_pointwise_montgomery(&mut v11_h, &mat_a11, &y11_hat);
    polyveck_reduce(&mut v11_h);
    polyveck_invntt_tomont(&mut v11_h);

    let v11 = v11_h;

    // Decompose v11 into (v11_h, v11_l)
    polyveck_reduce(&mut v11_h);
    polyveck_caddq(&mut v11_h);
    polyveck_decompose(&mut v11_h, &mut v11_l);

    // Compute: e11
    let mut e11 = b11_h;
    polyveck_shiftl(&mut e11);
    polyveck_add(&mut e11, &b11_l);

    // Result of the matrix multiplication of mat_a11 with s11.
    let mut v_a11_s11 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut v_a11_s11, &mat_a11, &s11_hat);
    polyveck_reduce(&mut v_a11_s11);
    polyveck_invntt_tomont(&mut v_a11_s11);

    // Update e11 by subtracting v_a11_s11 from it.
    polyveck_sub(&mut e11, &v_a11_s11);

    let mut e11_hat = e11;
    polyveck_ntt(&mut e11_hat);

    // z := nil
    // let mut z = [Polyvecl::default(); 2]; // (z1, z2)

    // Allocate intermediates.
    let mut y12 = Polyvecl::default();
    let (mut v12_h, mut v12_l) = (Polyveck::default(), Polyveck::default());
    let mut state = KeccakState::default();

    let mut nonce = 0u16;

    let mut retry = 0;

    loop {
        if retry == KEYGEN_RETRY_THRESHOLD {
            break Err(KeygenError::ThresholdError);
        }

        // Sample intermediate vector y12
        polyvecl_uniform_gamma1(&mut y12, &r_prime, nonce);
        nonce += 1;

        // Compute v12_h as matrix-vector multiplication of mat_a12 and y12
        let mut y12_hat = y12;
        polyvecl_ntt(&mut y12_hat);
        
        polyvec_matrix_pointwise_montgomery(&mut v12_h, &mat_a12, &y12_hat);
        polyveck_reduce(&mut v12_h);
        polyveck_invntt_tomont(&mut v12_h);

        let v12 = v12_h;

        // Create v1 = v11 + v12 before decomposing v12 into (v12_h, v12_l)
        let mut v1 = v11;
        polyveck_add(&mut v1, &v12);

        // Decompose v1 into v1_h and v1_l
        let (mut v1_h, mut v1_l) = (Polyveck::default(), Polyveck::default());
        v1_h = v1;
        polyveck_reduce(&mut v1_h);
        polyveck_caddq(&mut v1_h);
        polyveck_decompose(&mut v1_h, &mut v1_l);

        // Pack bits of v1_h for calculating c
        let mut v1_h_packed = [0u8; K * POLYW1_PACKEDBYTES];

        // Get coefficients in proper range before packing.
        polyveck_reduce(&mut v1_h);
        polyveck_caddq(&mut v1_h);

        polyveck_pack_w1(v1_h_packed.as_mut_slice(), &v1_h);

        let mut c = [0u8; SEEDBYTES];

        state.init();
        shake256_absorb(&mut state, &r, R_SIZE);
        shake256_absorb(&mut state, &v1_h_packed, K * POLYW1_PACKEDBYTES);
        shake256_finalize(&mut state);
        shake256_squeeze(&mut c, SEEDBYTES, &mut state);

        // Got c, now use it to sample from B_tau.

        // Sample c1 and c2 from B_tau and convert them to NTT form.
        let mut c1 = Poly::default();
        poly_challenge_nonced(&mut c1, &c, 0);
        poly_ntt(&mut c1);

        let mut c2 = Poly::default();
        poly_challenge_nonced(&mut c2, &c, 1);
        poly_ntt(&mut c2);

        // Create z1 and z2.
        let (mut z1, mut z2) = (Polyvecl::default(), Polyvecl::default());

        polyvecl_pointwise_poly_montgomery(&mut z1, &c1, &s11_hat);
        polyvecl_invntt_tomont(&mut z1);
        polyvecl_add(&mut z1, &y11);
        polyvecl_reduce(&mut z1);

        polyvecl_pointwise_poly_montgomery(&mut z2, &c2, &s12_prime_hat);
        polyvecl_invntt_tomont(&mut z2);
        polyvecl_add(&mut z2, &y12);
        polyvecl_reduce(&mut z2);

        // Check norm of z = (z1 z2)^T
        if polyvecl_chknorm(&z1, (GAMMA1 - BETA) as i32) > 0 {
            retry += 1;
            continue;
        }

        if polyvecl_chknorm(&z2, (GAMMA1 - BETA) as i32) > 0 {
            retry += 1;
            continue;
        }

        let mut rpoly_h = v1;
        let (mut c1_e11, mut c2_e12) = (Polyveck::default(), Polyveck::default());

        polyveck_pointwise_poly_montgomery(&mut c1_e11, &c1, &e11_hat);
        polyveck_reduce(&mut c1_e11);
        polyveck_invntt_tomont(&mut c1_e11);
        polyveck_sub(&mut rpoly_h, &c1_e11);

        polyveck_pointwise_poly_montgomery(&mut c2_e12, &c2, &e12_hat);
        polyveck_reduce(&mut c2_e12);
        polyveck_invntt_tomont(&mut c2_e12);
        polyveck_sub(&mut rpoly_h, &c2_e12);

        polyveck_reduce(&mut rpoly_h);
        polyveck_caddq(&mut rpoly_h);

        let rpoly = rpoly_h;

        // Decompose rpoly_h in to rpoly_h and rpoly_l.
        let mut rpoly_l = Polyveck::default();
        
        polyveck_decompose(&mut rpoly_h, &mut rpoly_l);

        // Check norm of rpoly_l.
        if polyveck_chknorm(&rpoly_l, (GAMMA2 - BETA) as i32) > 0 {
            retry += 1;
            continue;
        }

        // Check whether rpoly_h and v1 have same high bits.
        let mut rh_vh_diff = rpoly_h;
        polyveck_sub(&mut rh_vh_diff, &v1_h);

        // Polyveck cannot be compared with each other using == operator.
        /*
        if rpoly_h != v1_full_h {
            continue;
        }
        */

        // The infinite norm must be strictly less than 1, ie be 0.
        if polyveck_chknorm(&rh_vh_diff, 1 as i32) > 0 {
            retry += 1;
            continue;
        }

        let (mut c1_b11_l, mut c2_b12_l) = (Polyveck::default(), Polyveck::default());

        polyveck_pointwise_poly_montgomery(&mut c1_b11_l, &c1, &b11_l_hat);
        polyveck_reduce(&mut c1_b11_l);
        polyveck_invntt_tomont(&mut c1_b11_l);
        // polyveck_add(&mut cb, &c1_b11_l);

        polyveck_pointwise_poly_montgomery(&mut c2_b12_l, &c2, &b12_l_hat);
        polyveck_reduce(&mut c2_b12_l);
        polyveck_invntt_tomont(&mut c2_b12_l);
        // polyveck_add(&mut cb, &c2_b12_l);

        let mut cb = Polyveck::default(); // -(c1_b11_l + c2_b12_l)
        polyveck_sub(&mut cb, &c1_b11_l);
        polyveck_sub(&mut cb, &c2_b12_l);
        // polyveck_reduce(&mut cb);
        // polyveck_caddq(&mut cb);

        // Check if norm of cb is less than GAMMA2.
        if polyveck_chknorm(&cb, GAMMA2 as i32) > 0 {
            retry += 1;
            continue;
        }

        // Make hint.
        let mut h = Polyveck::default();

        let mut hint_from = rpoly;
        polyveck_sub(&mut hint_from, &cb);
        polyveck_reduce(&mut hint_from);
        polyveck_caddq(&mut hint_from);

        let n = polyveck_make_hint_simple(&mut h, &cb, &hint_from);
        if n > OMEGA as i32 {
            retry += 1;
            continue;
        }

        // No packing of pk, just pass plain b12_h and z = (z1 z2)^T
        let pk = PublicKey {
            b12_h, z1, z2, h, c,

        };

        // No packing of sk either, just pass the fields.
        let sk = SecretKey {
            b11_l, b12_l, y11, y12, s11, s12_prime
        };

        break Ok((pk, sk));
    }
}

#[derive(Debug)]
pub enum SigError {
    ThresholdError, // Too many attempts retrying within a loop.
}

#[derive(Copy, Clone)]
pub struct Signature {
    z_i1: Polyvecl,
    z_i2: Polyvecl,
    h_i:  Polyveck,
    c_i: [u8; SEEDBYTES], // TODO: is c_i a hash or a polynomial?
}

pub fn generate_signature(
    msg:        &[u8],
    identity:   &[u8],
    params:     Params, 
    public_key: PublicKey,
    secret_key: SecretKey,
) -> Result<Signature, SigError> {

    assert_eq!(identity.len(), ID_SIZE);

    // Unpack params, public key and secret key.
    let Params { rho, b11_h } = params;
    let PublicKey { b12_h, z1, z2, h, c: _ } = public_key;
    let SecretKey { b11_l, b12_l, y11, y12, s11, s12_prime } = secret_key;

    // NTT form of b11_l and b12_l for later use.
    let (mut b11_l_hat, mut b12_l_hat) = (b11_l, b12_l);
    polyveck_ntt(&mut b11_l_hat);
    polyveck_ntt(&mut b12_l_hat);

    // NTT form of y11 and y12 for later use.
    let mut y11_hat = y11;
    polyvecl_ntt(&mut y11_hat);

    let mut y12_hat = y12;
    polyvecl_ntt(&mut y12_hat);

    // NTT form of s11 and s12_prim for later use.
    let mut s11_hat = s11;
    polyvecl_ntt(&mut s11_hat);

    let mut s12_prime_hat = s12_prime;
    polyvecl_ntt(&mut s12_prime_hat);

    // Compute seed r and r_prime
    let mut seedbuf = [0u8; R_SIZE + R_PRIME_SIZE]; // r || r_prime

    let mut hash_input = [0u8; SEEDBYTES + ID_SIZE]; // rho || ID_A
    hash_input[..SEEDBYTES].copy_from_slice(&rho);
    hash_input[SEEDBYTES..].copy_from_slice(&identity);

    // Write seed to seedbuf.
    shake256(&mut seedbuf, R_SIZE + R_PRIME_SIZE, &hash_input, SEEDBYTES + ID_SIZE);

    let mut r = [0u8; R_SIZE];
    let mut r_prime = [0u8; R_PRIME_SIZE];

    r.copy_from_slice(&seedbuf[..R_SIZE]);
    r_prime.copy_from_slice(&seedbuf[R_SIZE..]);

    // Seed Expansion: mat_a11, mat_a12.

    // Expand seed for both matrices.
    // 
    // Use first SEEDBYTES bytes of r for first matrix.
    let mut mat_a11 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a11, &r[..SEEDBYTES]);

    // Use later SEEDBYTES bytes of r for second matrix.
    let mut mat_a12 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a12, &r[SEEDBYTES..]);

    // Compute v11 and v12, and decompose them into their high and low
    //       bits versions.
    let mut v11 = Polyveck::default();

    polyvec_matrix_pointwise_montgomery(&mut v11, &mat_a11, &y11_hat);
    polyveck_reduce(&mut v11);
    polyveck_invntt_tomont(&mut v11);

    // Do the same for v12.
    let mut v12 = Polyveck::default();

    polyvec_matrix_pointwise_montgomery(&mut v12, &mat_a12, &y12_hat);
    polyveck_reduce(&mut v12);
    polyveck_invntt_tomont(&mut v12);

    // Create v1 as sum of v11 and v12
    let mut v1 = v11;
    polyveck_add(&mut v1, &v12);
    polyveck_reduce(&mut v1);
    polyveck_caddq(&mut v1);

    // Decompose v1 into high and low parts.
    let mut v1_h = v1;
    let mut v1_l = Polyveck::default();
    polyveck_decompose(&mut v1_h, &mut v1_l);

    // Compute e11 and e12.
    let mut e11 = b11_h;
    polyveck_shiftl(&mut e11);
    polyveck_add(&mut e11, &b11_l);

    // Result of the matrix multiplication of mat_a11 with s11.
    let mut v_a11_s11 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut v_a11_s11, &mat_a11, &s11_hat);
    polyveck_reduce(&mut v_a11_s11);
    polyveck_invntt_tomont(&mut v_a11_s11);

    // Update e11 by subtracting v_a11_s11 from it.
    polyveck_sub(&mut e11, &v_a11_s11);

    // Convert e11 into NTT form.
    let mut e11_hat = e11;
    polyveck_ntt(&mut e11_hat);

    // Do the same for e12
    let mut e12 = b12_h;
    polyveck_shiftl(&mut e12);
    polyveck_add(&mut e12, &b12_l);

    // Result of the matrix multiplication of mat_a12 with s12.
    let mut v_a12_s12_prime = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut v_a12_s12_prime, &mat_a12, &s12_prime_hat);
    polyveck_reduce(&mut v_a12_s12_prime);
    polyveck_invntt_tomont(&mut v_a12_s12_prime);

    // Update e12 by subtracting v_a12_s12_prime from it.
    polyveck_sub(&mut e12, &v_a12_s12_prime);

    // Convert e12 into NTT form.
    let mut e12_hat = e12;
    polyveck_ntt(&mut e12_hat);

    let mut nonce = 0u16;
    let (mut y_i1, mut y_i2) = (Polyvecl::default(), Polyvecl::default());

    // Reset state, this will be used later to sample c_i
    let mut state = KeccakState::default();

    let mut retry = 0;

    let mut fail_reason = FailReason::default();

    let mut prehashed = KeccakState::default();
    shake256_absorb(&mut prehashed, &r, R_SIZE);
    shake256_absorb(&mut prehashed, msg, msg.len());

    loop {
        if retry == SIG_RETRY_THRESHOLD {
            break Err(SigError::ThresholdError);
        }

        // Sample intermediate vectors y_i1 and y_i2.
        polyvecl_uniform_gamma1(&mut y_i1, &r_prime, nonce);
        nonce += 1;

        // NTT transform y_i1 for future use.
        let mut y_i1_hat = y_i1;
        polyvecl_ntt(&mut y_i1_hat);

        polyvecl_uniform_gamma1(&mut y_i2, &r_prime, nonce);
        nonce += 1;

        // NTT transform y_i2 for future use.
        let mut y_i2_hat = y_i2;
        polyvecl_ntt(&mut y_i2_hat);

        // Compute v_i1 and v_i2
        let mut v_i1 = Polyveck::default();
        polyvec_matrix_pointwise_montgomery(&mut v_i1, &mat_a11, &y_i1_hat);
        polyveck_reduce(&mut v_i1);
        polyveck_invntt_tomont(&mut v_i1);

        let mut v_i2 = Polyveck::default();
        polyvec_matrix_pointwise_montgomery(&mut v_i2, &mat_a12, &y_i2_hat);
        polyveck_reduce(&mut v_i2);
        polyveck_invntt_tomont(&mut v_i2);

        // Compute v_i.
        let mut v_i = v_i1;
        polyveck_add(&mut v_i, &v_i2);

        // Decompose v_i into high and low parts.
        let (mut v_i_h, mut v_i_l) = (Polyveck::default(), Polyveck::default());
        v_i_h = v_i;
        // polyveck_reduce(&mut v_i_h);
        polyveck_caddq(&mut v_i_h);
        polyveck_decompose(&mut v_i_h, &mut v_i_l);

        // Compute c_i = CRH(r || v_... || msg)
        
        // Sum of high parts of v1 and v_i.
        let mut v_h_sum = v1_h;
        polyveck_add(&mut v_h_sum, &v_i_h);

        // polyveck_reduce(&mut v_h_sum);
        polyveck_caddq(&mut v_h_sum);

        let mut v_h_sum_packed = [0u8; K * POLYW1_PACKEDBYTES];
        polyveck_pack_w1(v_h_sum_packed.as_mut_slice(), &v_h_sum);

        let mut c_i = [0u8; SEEDBYTES];

        let mut hasher_state = prehashed;
        shake256_absorb(&mut hasher_state, &v_h_sum_packed, K * POLYW1_PACKEDBYTES);
        shake256_finalize(&mut hasher_state);
        shake256_squeeze(&mut c_i, SEEDBYTES, &mut hasher_state);

        // Sample in ball, using c_i seed.
        let mut c_i1 = Poly::default();
        // poly_challenge_nonced(&mut c_i1, &c_i, 0);
        poly_challenge_nonced(&mut c_i1, &c_i, 0);
        poly_ntt(&mut c_i1);

        let mut c_i2 = Poly::default();
        poly_challenge_nonced(&mut c_i2, &c_i, 1);
        poly_ntt(&mut c_i2);

        let (mut z_i1, mut z_i2) = (Polyvecl::default(), Polyvecl::default());

        polyvecl_pointwise_poly_montgomery(&mut z_i1, &c_i1, &s11_hat);
        polyvecl_invntt_tomont(&mut z_i1);
        polyvecl_add(&mut z_i1, &y_i1);
        polyvecl_reduce(&mut z_i1);

        // Check norm of z_i1.
        if polyvecl_chknorm(&z_i1, (GAMMA1 - BETA) as i32) > 0 {
            retry += 1;
            fail_reason.z_i1 += 1;
            continue;
        }

        polyvecl_pointwise_poly_montgomery(&mut z_i2, &c_i2, &s12_prime_hat);
        polyvecl_invntt_tomont(&mut z_i2);
        polyvecl_add(&mut z_i2, &y_i2);
        polyvecl_reduce(&mut z_i2);

        // Check norm of z_i2.
        if polyvecl_chknorm(&z_i2, (GAMMA1 - BETA) as i32) > 0 {
            retry += 1;
            fail_reason.z_i2 += 1;
            continue;
        }

        // Compute r_poly and decompose it to high and low parts.
        let mut r_i_poly = v_i;
        let (mut ci1_e11, mut ci2_e12) = (Polyveck::default(), Polyveck::default());

        polyveck_pointwise_poly_montgomery(&mut ci1_e11, &c_i1, &e11_hat);
        polyveck_invntt_tomont(&mut ci1_e11);
        // polyveck_reduce(&mut ci1_e11);
        polyveck_sub(&mut r_i_poly, &ci1_e11);

        polyveck_pointwise_poly_montgomery(&mut ci2_e12, &c_i2, &e12_hat);
        polyveck_invntt_tomont(&mut ci2_e12);
        // polyveck_reduce(&mut ci2_e12);
        polyveck_sub(&mut r_i_poly, &ci2_e12);

        polyveck_reduce(&mut r_i_poly);
        polyveck_caddq(&mut r_i_poly);

        let (mut r_i_poly_h, mut r_i_poly_l) = (Polyveck::default(), Polyveck::default());
        
        r_i_poly_h = r_i_poly;
        polyveck_decompose(&mut r_i_poly_h, &mut r_i_poly_l);

        // Check norm of r_i_poly_l
        if polyveck_chknorm(&r_i_poly_l, (GAMMA2 - BETA) as i32) > 0 {
            retry += 1;
            fail_reason.r_i_l += 1;
            continue;
        }

        // Check whether high parts of r_poly and v_i are same or not.
        let mut rih_vih_diff = r_i_poly_h;
        polyveck_sub(&mut rih_vih_diff, &v_i_h);
        polyveck_reduce(&mut rih_vih_diff);

        if polyveck_chknorm(&rih_vih_diff, 1 as i32) > 0 {
            retry += 1;
            fail_reason.rih_vih += 1;
            continue;
        }

        let mut c_i1_b11_l = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c_i1_b11_l, &c_i1, &b11_l_hat);
        polyveck_reduce(&mut c_i1_b11_l);
        polyveck_invntt_tomont(&mut c_i1_b11_l);

        let mut c_i2_b12_l = Polyveck::default();
        polyveck_pointwise_poly_montgomery(&mut c_i2_b12_l, &c_i2, &b12_l_hat);
        polyveck_reduce(&mut c_i2_b12_l);
        polyveck_invntt_tomont(&mut c_i2_b12_l);

        let mut cib = Polyveck::default();
        polyveck_sub(&mut cib, &c_i1_b11_l);
        polyveck_sub(&mut cib, &c_i2_b12_l);
        polyveck_reduce(&mut cib);
        // polyveck_caddq(&mut cib);

        // Check norm of cib
        if polyveck_chknorm(&cib, 2 * GAMMA2_I32) > 0 {
            retry += 1;
            fail_reason.cb += 1;
            continue;
        }

        // Make hint
        let mut hint_from = r_i_poly;
        polyveck_sub(&mut hint_from, &cib);
        polyveck_reduce(&mut hint_from);
        polyveck_caddq(&mut hint_from);

        let mut h_i = Polyveck::default();
        let n = polyveck_make_hint_simple(&mut h_i, &cib, &hint_from);

        if n > OMEGA as i32 {
            retry += 1;
            fail_reason.omega += 1;
            continue;
        }

        let signature = Signature {
            z_i1,
            z_i2,
            h_i,
            c_i,
        };

        break Ok(signature);
    }
}

pub fn verify_sign(
    msg:       &[u8],
    identity:  &[u8],
    signature: Signature,
    params:    Params,
    pk:        PublicKey,
) -> bool {
    assert_eq!(identity.len(), ID_SIZE);

    let Params { rho, b11_h } = params;
    let PublicKey { b12_h, z1, z2, h, c } = pk;
    let Signature { z_i1, z_i2, h_i, c_i } = signature;

    // NTT + shifted form of b11_h and b12_h for later use.
    let (mut b11_h_hat, mut b12_h_hat) = (b11_h, b12_h);
    polyveck_shiftl(&mut b11_h_hat);
    polyveck_ntt(&mut b11_h_hat);

    polyveck_shiftl(&mut b12_h_hat);
    polyveck_ntt(&mut b12_h_hat);

    // NTT versions of z's for later use.
    let mut z1_hat = z1;
    polyvecl_ntt(&mut z1_hat);

    let mut z2_hat = z2;
    polyvecl_ntt(&mut z2_hat);

    let mut z_i1_hat = z_i1;
    polyvecl_ntt(&mut z_i1_hat);

    let mut z_i2_hat = z_i2;
    polyvecl_ntt(&mut z_i2_hat);

    // Compute seed: r
    let mut r = [0u8; R_SIZE];

    let mut hash_input = [0u8; SEEDBYTES + ID_SIZE]; // rho || ID_A
    hash_input[..SEEDBYTES].copy_from_slice(&rho);
    hash_input[SEEDBYTES..].copy_from_slice(&identity);

    // Write seed to r,
    shake256(&mut r, R_SIZE, &hash_input, SEEDBYTES + ID_SIZE);

    // Seed Expansion: mat_a11, mat_a12.

    // Expand seed for both matrices.
    // 
    // Use first SEEDBYTES bytes of r for first matrix.
    let mut mat_a11 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a11, &r[..SEEDBYTES]);

    // Use later SEEDBYTES bytes of r for second matrix.
    let mut mat_a12 = [Polyvecl::default(); K];
    polyvec_matrix_expand(&mut mat_a12, &r[SEEDBYTES..]);

    // Sample c1 and c2 from B_tau using c and convert them to NTT form.
    let mut c1 = Poly::default();
    poly_challenge_nonced(&mut c1, &c, 0);
    poly_ntt(&mut c1);

    let mut c2 = Poly::default();
    poly_challenge_nonced(&mut c2, &c, 1);
    poly_ntt(&mut c2);

    // Do the same thing for c_i1 and c_i2 using c_i
    let mut c_i1 = Poly::default();
    poly_challenge_nonced(&mut c_i1, &c_i, 0);
    poly_ntt(&mut c_i1);

    let mut c_i2 = Poly::default();
    poly_challenge_nonced(&mut c_i2, &c_i, 1);
    poly_ntt(&mut c_i2);

    // Check whether z and z_i have correct norms.
    if polyvecl_chknorm(&z1, (GAMMA1 - BETA) as i32) > 0 {
        return false;
    }

    if polyvecl_chknorm(&z2, (GAMMA1 - BETA) as i32) > 0 {
        return false;
    }

    if polyvecl_chknorm(&z_i1, (GAMMA1 - BETA) as i32) > 0 {
        return false;
    }

    if polyvecl_chknorm(&z_i2, (GAMMA1 - BETA) as i32) > 0 {
        return false;
    }

    // Calculate cb and cib for v1 and v_i_prime.
    let mut c1_b11_h = Polyveck::default();
    polyveck_pointwise_poly_montgomery(&mut c1_b11_h, &c1, &b11_h_hat);
    polyveck_reduce(&mut c1_b11_h);
    polyveck_invntt_tomont(&mut c1_b11_h);

    let mut c2_b12_h = Polyveck::default();
    polyveck_pointwise_poly_montgomery(&mut c2_b12_h, &c2, &b12_h_hat);
    polyveck_reduce(&mut c2_b12_h);
    polyveck_invntt_tomont(&mut c2_b12_h);

    let mut cb = c1_b11_h;
    polyveck_add(&mut cb, &c2_b12_h);

    let mut c_i1_b11_h = Polyveck::default();
    polyveck_pointwise_poly_montgomery(&mut c_i1_b11_h, &c_i1, &b11_h_hat);
    polyveck_reduce(&mut c_i1_b11_h);
    polyveck_invntt_tomont(&mut c_i1_b11_h);

    let mut c_i2_b12_h = Polyveck::default();
    polyveck_pointwise_poly_montgomery(&mut c_i2_b12_h, &c_i2, &b12_h_hat);
    polyveck_reduce(&mut c_i2_b12_h);
    polyveck_invntt_tomont(&mut c_i2_b12_h);

    let mut cib = c_i1_b11_h;
    polyveck_add(&mut cib, &c_i2_b12_h);

    // Calculate product of A's and z's for v1 and v_i_prime.

    let mut a11_z1 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut a11_z1, &mat_a11, &z1_hat);
    polyveck_reduce(&mut a11_z1);
    polyveck_invntt_tomont(&mut a11_z1);

    let mut a12_z2 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut a12_z2, &mat_a12, &z2_hat);
    polyveck_reduce(&mut a12_z2);
    polyveck_invntt_tomont(&mut a12_z2);

    let mut a11_z_i1 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut a11_z_i1, &mat_a11, &z_i1_hat);
    polyveck_reduce(&mut a11_z_i1);
    polyveck_invntt_tomont(&mut a11_z_i1);

    let mut a12_z_i2 = Polyveck::default();
    polyvec_matrix_pointwise_montgomery(&mut a12_z_i2, &mat_a12, &z_i2_hat);
    polyveck_reduce(&mut a12_z_i2);
    polyveck_invntt_tomont(&mut a12_z_i2);

    // Calculate v1's pre-hint
    let mut prehint1 = a11_z1;
    polyveck_add(&mut prehint1, &a12_z2);
    // polyveck_reduce(&mut prehint1);
    polyveck_sub(&mut prehint1, &cb);
    polyveck_reduce(&mut prehint1);
    polyveck_caddq(&mut prehint1);

    // Calculate v_i_prime's pre-hint
    let mut prehint_i = a11_z_i1;
    polyveck_add(&mut prehint_i, &a12_z_i2);
    polyveck_sub(&mut prehint_i, &cib);
    polyveck_reduce(&mut prehint_i);
    polyveck_caddq(&mut prehint_i);

    // Use hint on v1
    let mut v1_prime = Polyveck::default();
    polyveck_use_hint_simple(&mut v1_prime, &h, &prehint1);

    // Use hint on v_i_prime
    let mut v_i_prime = Polyveck::default();
    polyveck_use_hint_simple(&mut v_i_prime, &h_i, &prehint_i);

    // Sum v1 and v_i_prime.
    let mut v_sum = v1_prime;
    polyveck_add(&mut v_sum, &v_i_prime);

    polyveck_reduce(&mut v_sum);
    polyveck_caddq(&mut v_sum);

    let mut v_sum_packed = [0u8; K * POLYW1_PACKEDBYTES];
    polyveck_pack_w1(v_sum_packed.as_mut_slice(), &v_sum);

    let mut final_hash = [0u8; SEEDBYTES];

    let mut state = KeccakState::default();
    shake256_absorb(&mut state, &r, R_SIZE);
    shake256_absorb(&mut state, msg, msg.len());
    shake256_absorb(&mut state, &v_sum_packed, K * POLYW1_PACKEDBYTES);
    shake256_finalize(&mut state);
    shake256_squeeze(&mut final_hash, SEEDBYTES, &mut state);

    if c_i != final_hash {
        return false;
    }

    return true
}
