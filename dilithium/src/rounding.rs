use crate::params::*;

/// For finite field element a, compute a0, a1 such that
/// a mod^+ Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
/// Assumes a to be standard representative.
///
/// Returns a1.
pub fn power2round(a: i32, a0: &mut i32) -> i32 {
  let a1 = (a + (1 << (D - 1)) - 1) >> D;
  *a0 = a - (a1 << D);
  return a1;
}

/// For finite field element a, compute high and low bits a0, a1 such
/// that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
/// if a1 = (Q-1)/ALPHA where we set a1 = 0 and
/// -ALPHA/2 <= a0 = a mod^+ Q - Q < 0. Assumes a to be standard
/// representative.
///
/// Returns a1.
pub fn decompose(a0: &mut i32, a: i32) -> i32 {
  let mut a1 = (a + 127) >> 7;
  if GAMMA2 == (Q - 1) / 32 {
    a1 = (a1 * 1025 + (1 << 21)) >> 22;
    a1 &= 15;
  } else if GAMMA2 == (Q - 1) / 88 {
    a1 = (a1 * 11275 + (1 << 23)) >> 24;
    a1 ^= ((43 - a1) >> 31) & a1;
  }
  *a0 = a - a1 * 2 * GAMMA2_I32;
  *a0 -= (((Q_I32 - 1) / 2 - *a0) >> 31) & Q_I32;
  a1
}

/// Compute hint bit indicating whether the low bits of the
/// input element overflow into the high bits.
///
/// Returns 1 if overflow.
pub fn make_hint(a0: i32, a1: i32) -> u8 {
  if a0 > GAMMA2_I32 || a0 < -GAMMA2_I32 || (a0 == -GAMMA2_I32 && a1 != 0) {
    return 1;
  }
  return 0;
}

/// Correct high bits according to hint.
///
/// Returns corrected high bits.
pub fn use_hint(a: i32, hint: u8) -> i32 {
  let mut a0 = 0i32;
  let a1 = decompose(&mut a0, a);
  if hint == 0 {
    return a1;
  }

  if GAMMA2 == (Q - 1) / 32 {
    if a0 > 0 {
      return (a1 + 1) & 15;
    } else {
      return (a1 - 1) & 15;
    }
  } else {
    if a0 > 0 {
      if a1 == 43 {
        return 0;
      } else {
        return a1 + 1;
      };
    } else {
      if a1 == 0 {
        return 43;
      } else {
        return a1 - 1;
      }
    }
  }
}

// Input: `a` must be positive.
//
// Returns r1
pub fn decompose_scaled(a0: &mut i32, a: i32) -> i32 {
    const ALPHA: i32 = 4 * GAMMA2_I32;

    let mut r0 = a % ALPHA;
    if r0 > ALPHA / 2 {
        r0 -= ALPHA;
    } else if r0 <= -ALPHA / 2 {
        r0 += ALPHA;
    }

    let r1: i32;
    if a - r0 == Q_I32 - 1 {
        r1 = 0;
        r0 -= 1;
    } else {
        r1 = (a - r0) / ALPHA;
    }

    *a0 = r0;

    r1
}

pub fn use_hint_scaled(a: i32, hint: u8) -> i32 {
    let mut a0 = 0i32;
    let a1 = decompose_scaled(&mut a0, a);

    if hint == 0 {
        return a1;
    }

    const M: i32 = if GAMMA2 == (Q - 1) / 32 {
        8
    } else {
        22
    };

    if a0 > 0 {
        return (((a1 + 1) % M) + M) % M;
    } else {
        return (((a1 - 1) % M) + M) % M;
    }
}

pub fn make_hint_scaled(a0: i32, a1: i32) -> u8 {
    /*
    if a0 > GAMMA2_I32 || a0 < -GAMMA2_I32 || (a0 == -GAMMA2_I32 && a1 != 0) {
        return 1;
    }
    return 0;
    */

    let z = a0;
    let r = a1;

    let r_z = (((r + z) % Q_I32) + Q_I32) % Q_I32;

    let mut r_z_l = 0;
    let r_z_h = decompose_scaled(&mut r_z_l, r_z);

    let mut r_l = 0;
    let r_h = decompose_scaled(&mut r_l, r);

    (r_h != r_z_h) as u8
}

pub fn decompose_simple(a0: &mut i32, a: i32) -> i32 {
    assert!(a >= 0);
    const ALPHA: i32 = 2 * GAMMA2_I32;

    let mut r0 = a % ALPHA;
    if r0 > ALPHA / 2 {
        r0 -= ALPHA;
    } else if r0 <= -ALPHA / 2 {
        r0 += ALPHA;
    }

    let r1: i32;
    if a - r0 == Q_I32 - 1 {
        r1 = 0;
        r0 -= 1;
    } else {
        r1 = (a - r0) / ALPHA;
    }

    *a0 = r0;

    r1
}

pub fn use_hint_simple(a: i32, hint: u8) -> i32 {
    assert!(hint == 0 || hint == 1);
    assert!(a >= 0);

    let mut a0 = 0i32;
    let a1 = decompose_simple(&mut a0, a);

    if hint == 0 {
        return a1;
    }

    const M: i32 = if GAMMA2 == (Q - 1) / 32 {
        16
    } else {
        44
    };

    if a0 > 0 {
        return (((a1 + 1) % M) + M) % M;
    } else {
        return (((a1 - 1) % M) + M) % M;
    }
}

pub fn make_hint_simple(a0: i32, a1: i32) -> u8 {
    assert!(a1 >= 0);
    /*
    if a0 > GAMMA2_I32 || a0 < -GAMMA2_I32 || (a0 == -GAMMA2_I32 && a1 != 0) {
        return 1;
    }
    return 0;
    */

    let z = a0;
    let r = a1;

    let r_z = (((r + z) % Q_I32) + Q_I32) % Q_I32;

    let mut r_z_l = 0;
    let r_z_h = decompose_simple(&mut r_z_l, r_z);

    let mut r_l = 0;
    let r_h = decompose_simple(&mut r_l, r);

    (r_h != r_z_h) as u8
}

