# EUROCRYPT

To implement our Certificateless Signature Scheme we have modified the
dilithium implementation present at [Argyle-Software/dilithium](https://github.com/Argyle-Software/dilithium).
The modified code is present in the `dilithium` directory.

We compare the running times of our proposed Certificateless Signature Scheme, as printed `NOPKI` in the console, with that of a version of dilithium that use PKI infrastructure, as printed `Baseline Dilithium` in the console.

The comparisons are happening for differnt modes of dilithium. The three modes of dilithium `mode2`, `mode3`, and `mode5` provide different levels of security, respectively 128-bit, 192-bit, and 256-bit of security.

## Running the benchmarks.
To run the benchmark for any level of dilithim just go to the `dilithium` directory.

```
cd dilithium
```

The run any of the following commands:

For mode2:
```
cargo run --release --features mode2
```

For mode3:
```
cargo run --release --features mode3
```

For mode5:
```
cargo run --release --features mode5
```

By default, just running

```shell
cargo run --release
```

...will print a helpful error message and ask to select a mode.

This might take a few seconds to compile and run.

NOTE: you will need `Rust` setup and its package manager `cargo` for this. To know how to
install Rust, follow [Install Rust](https://rust-lang.org/tools/install/).

For example, the above commands will print something like this.

```
------------------NOPKI - mode2----------------
iterations: 1000, message length: 1024bytes
average ppk time:       492.634µs
average keygen time:    727.057µs
average signature time: 499.767µs
average verify time:    338.064µs

----------------Baseline Dilithium---------------
iterations: 1000
average keygen time:    109.997µs
average signature time: 392.252µs
average verify time:    116.739µs

```

This compares our NO-PKI dilithium with one that uses PKI infrastructure.

This particular run is for dilithium2, i.e. it was ran with `--features mode2` CLI option.

`ppk time`: time taken to generate partial private key. \
`keygen time`: time taken to generate public and secret key \
`sig time`: time taken to sign a message \
`verify time`: time taken to verify the signature for a message

...and also print the number of times the whole process is done. By default 1000
iterations are done, each for a different message of length 1024 bytes.

