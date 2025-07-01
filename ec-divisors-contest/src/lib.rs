#![no_std]
#![allow(static_mut_refs)]

use ec_divisors::{
  DivisorCurve as DivisorCurveRef, Poly as PolyRef, ScalarDecomposition as ScalarDecompositionRef,
};
pub use ec_divisors_contest_src::{precompute, Precomp};
use ec_divisors_contest_src::{DivisorCurve, Poly, ScalarDecomposition};

use ciphersuite::{
  group::{ff::Field, Group},
  Ciphersuite, Ed25519,
};
use dalek_ff_group::{EdwardsPoint, FieldElement, Scalar};

use zeroize::Zeroizing;

use rand_chacha::ChaCha20Rng;
use rand_core::SeedableRng;

use std_shims::sync::OnceLock;

pub fn check_init_ref(point: &EdwardsPoint, scalar: &Scalar) {
  let scalar = ScalarDecompositionRef::new(*scalar).unwrap();
  let point = Zeroizing::new(*point * scalar.scalar());
  let (_, _) = <EdwardsPoint as DivisorCurveRef>::to_xy(*point).unwrap();
}

pub fn check_init_contest(point: &EdwardsPoint, scalar: &Scalar) {
  let scalar = ScalarDecomposition::new(*scalar).expect("failed scalar decomposition");
  let point = Zeroizing::new(*point * scalar.scalar());
  let (_, _) = <EdwardsPoint as DivisorCurve>::to_xy(*point).expect("zero scalar was decomposed");
}

pub fn run_bench_ref(
  point: &EdwardsPoint,
  scalar: &Scalar,
) -> (PolyRef<FieldElement>, ScalarDecompositionRef<Scalar>) {
  let scalar = ScalarDecompositionRef::new(*scalar).unwrap();
  (scalar.scalar_mul_divisor(*point), scalar)
}

pub fn run_bench_contest(
  point: &EdwardsPoint,
  scalar: &Scalar,
) -> (Poly<FieldElement>, ScalarDecomposition<Scalar>) {
  let scalar = ScalarDecomposition::new(*scalar).unwrap();
  (scalar.scalar_mul_divisor(*point), scalar)
}

// For error: no global memory allocator found but one is required; link to std or add `#[global_allocator]` to a static item that implements the GlobalAlloc trait
// dlmalloc is the allocator used on wasm32-unknown-unknown: https://doc.rust-lang.org/rustc/platform-support/wasm32-unknown-unknown.html
#[cfg(target_arch = "wasm32")]
#[global_allocator]
static ALLOCATOR: dlmalloc::GlobalDlmalloc = dlmalloc::GlobalDlmalloc;

#[cfg(target_arch = "wasm32")]
use core::{
  panic::PanicInfo,
  result::{
    Result,
    Result::{Err, Ok},
  },
  unimplemented,
};

#[cfg(target_arch = "wasm32")]
use getrandom::{register_custom_getrandom, Error};

// For error: `#[panic_handler]` function required, but not found
#[cfg(target_arch = "wasm32")]
#[panic_handler]
fn panic(_: &PanicInfo) -> ! {
  loop {}
}

// https://forum.dfinity.org/t/module-imports-function-wbindgen-describe-from-wbindgen-placeholder-that-is-not-exported-by-the-runtime/11545/8
#[cfg(target_arch = "wasm32")]
pub fn custom_getrandom(_buf: &mut [u8]) -> Result<(), Error> {
  unimplemented!("custom_getrandom not implemented");
}

#[cfg(target_arch = "wasm32")]
register_custom_getrandom!(custom_getrandom);

struct EcDivisorsParamsRef {
  point: EdwardsPoint,
  scalar: Scalar,
}

struct EcDivisorsParams {
  point: EdwardsPoint,
  scalar: Scalar,
}

impl EcDivisorsParamsRef {
  pub fn new(rng_seed: [u8; 32]) -> Self {
    let point = EdwardsPoint::generator();
    let mut rng = ChaCha20Rng::from_seed(rng_seed);
    let scalar = <Ed25519 as Ciphersuite>::F::random(&mut rng);
    check_init_ref(&point, &scalar);
    Self { point, scalar }
  }
}

impl EcDivisorsParams {
  pub fn new(rng_seed: [u8; 32]) -> Self {
    let point = EdwardsPoint::generator();
    let mut rng = ChaCha20Rng::from_seed(rng_seed);
    let scalar = <Ed25519 as Ciphersuite>::F::random(&mut rng);
    check_init_contest(&point, &scalar);
    Self { point, scalar }
  }
}

static mut EC_DIVISORS_PARAMS: OnceLock<EcDivisorsParams> = OnceLock::new();
static mut EC_DIVISORS_PARAMS_REF: OnceLock<EcDivisorsParamsRef> = OnceLock::new();

// Tests using different seed with for https://github.com/j-berman/wasm-cycles
#[no_mangle]
pub extern "C" fn case_scalar_mul_divisor_ref1() {
  unsafe {
    EC_DIVISORS_PARAMS_REF = OnceLock::new();
    let _ = EC_DIVISORS_PARAMS_REF.get_or_init(|| EcDivisorsParamsRef::new([0xff; 32]));
  };
}

#[no_mangle]
pub extern "C" fn case_scalar_mul_divisor_ref2() {
  unsafe {
    EC_DIVISORS_PARAMS_REF = OnceLock::new();
    let _ = EC_DIVISORS_PARAMS_REF.get_or_init(|| EcDivisorsParamsRef::new([0xde; 32]));
  };
}

#[no_mangle]
pub extern "C" fn case_scalar_mul_divisor_contest1() {
  unsafe {
    EC_DIVISORS_PARAMS = OnceLock::new();
    let _ = EC_DIVISORS_PARAMS.get_or_init(|| EcDivisorsParams::new([0xff; 32]));
  };
}

#[no_mangle]
pub extern "C" fn case_scalar_mul_divisor_contest2() {
  unsafe {
    EC_DIVISORS_PARAMS = OnceLock::new();
    let _ = EC_DIVISORS_PARAMS.get_or_init(|| EcDivisorsParams::new([0xde; 32]));
  };
}

#[no_mangle]
pub extern "C" fn test_scalar_mul_divisor_ref() {
  core::hint::black_box(unsafe {
    let params = EC_DIVISORS_PARAMS_REF.get().unwrap();
    run_bench_ref(&params.point, &params.scalar);
  });
}

#[no_mangle]
pub extern "C" fn test_scalar_mul_divisor_contest() {
  core::hint::black_box(unsafe {
    let params = EC_DIVISORS_PARAMS.get().unwrap();
    run_bench_contest(&params.point, &params.scalar);
  });
}
