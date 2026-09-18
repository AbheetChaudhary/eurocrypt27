use rand::prelude::*;

pub fn randombytes(x: &mut [u8], len: usize) {
  rand::rng().fill_bytes(&mut x[..len])
}
