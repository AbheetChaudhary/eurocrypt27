#![allow(dead_code, unused)]

use pqc_dilithium::euro;
use std::time::{Instant, Duration};

/*
fn main() {
    let mut rho = [0u8; 32];
    rand::fill(&mut rho);

    let mut identity = [0u8; 64];
    rand::fill(&mut identity);

    let (params, mpk, msk, ppk, upk) = euro::kgc_setup(&identity, &rho);

    let (pk, sk) = euro::user_keygen(params, ppk, &identity, upk, msk);

    // let message: &[u8] = b"lorem ipsum dolor sit amet";

    let mut message = vec![0u8; 4096];
    rand::fill(&mut message);
    let sig = euro::sign(params, pk, sk, &identity, &message);

    let verify_result = euro::verify(params, mpk, pk, upk, &identity, &message, sig);

    if verify_result {
        println!("verify success!");
    } else {
        println!("verify failed!");
    }
}
*/

fn main() {
    let mode_name = if cfg!(feature = "mode2") {
        "mode2"
    } else if cfg!(feature = "mode3") {
        "mode3"
    } else if cfg!(feature = "mode5") {
        "mode5"
    } else {
        println!("Please specify dilithium mode: one of mode2, mode3 or mode5");
        println!("example: cargo run --release --features mode3");
        panic!("No mode selected!");
    };


    println!("------------------NOPKI - {}----------------", mode_name);
    print_timing_info(1000);
    println!("");
    println!("----------------Baseline Dilithium---------------");
    baseline_dilithium_timings(1000);
    println!("");
}

fn print_timing_info(count: usize) {
    let mut i = 0;

    // A message vector of 1KB size.
    let mut msg = vec![0u8; 1024];

    // Durations of differnet stages.
    let mut ppk_duration = Duration::ZERO;
    let mut keygen_duration = Duration::ZERO;
    let mut sig_duration = Duration::ZERO;
    let mut verify_duration = Duration::ZERO;

    // Number of times we successfully verified.
    let mut verify_success_count = 0;

    while i < count {
        let mut identity = [0u8; 64];
        let mut rho      = [0u8; 32];

        // Fill identity and rho with random bytes.
        rand::fill(&mut identity[..]);
        rand::fill(&mut rho[..]);

        // Fill the message with random bytes.
        rand::fill(&mut msg);

        let ppk_begin = Instant::now();
        let (params, mpk, msk, ppk, upk) = euro::kgc_setup(&identity, &rho);
        ppk_duration += ppk_begin.elapsed();

        let keygen_begin = Instant::now();
        let (pk, sk) = euro::user_keygen(params, ppk, &identity, upk, msk);
        keygen_duration += keygen_begin.elapsed();

        let sig_begin = Instant::now();
        let sig = euro::sign(params, pk, sk, &identity, &msg);
        sig_duration += sig_begin.elapsed();

        let verify_begin = Instant::now();
        let verify_result = euro::verify(params, mpk, pk, upk, &identity, &msg, sig);
        verify_duration += verify_begin.elapsed();

        if verify_result {
            verify_success_count += 1;
        }

        i += 1;
    }

    use std::convert::TryInto;

    println!("iterations: {}, verify success count: {}, message length: 1024bytes",
        count, verify_success_count);
    println!("average ppk time:       {:.3?}", ppk_duration / count.try_into().unwrap());
    println!("average keygen time:    {:.3?}", keygen_duration / count.try_into().unwrap());
    println!("average signature time: {:.3?}", sig_duration / count.try_into().unwrap());
    println!("average verify time:    {:.3?}", verify_duration / count.try_into().unwrap());
}

fn baseline_dilithium_timings(count: usize) {
    let mut keygen_duration = Duration::ZERO;
    let mut signature_duration = Duration::ZERO;
    let mut verify_duration = Duration::ZERO;

    let mut i = 0;

    let msg = b"lorem ipsum dolor sit amet";

    use pqc_dilithium as baseline;

    while i < count {
        let begin_keygen = Instant::now();
        let keypair = baseline::Keypair::generate();
        keygen_duration += begin_keygen.elapsed();

        let begin_signature = Instant::now();
        let signature = keypair.sign(msg);
        signature_duration += begin_signature.elapsed();

        let begin_verify = Instant::now();
        let verify_result = baseline::verify(&signature, msg, &keypair.public);
        verify_duration += begin_verify.elapsed();

        if verify_result.is_err() {
            println!("baseline dilithium failed once");
        }

        i += 1;
    }

    use std::convert::TryInto;

    println!("iterations: {count}");
    println!("average keygen time:    {:.3?}", keygen_duration / count.try_into().unwrap());
    println!("average signature time: {:.3?}", signature_duration / count.try_into().unwrap());
    println!("average verify time:    {:.3?}", verify_duration / count.try_into().unwrap());
}

