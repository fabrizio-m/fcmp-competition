use ciphersuite::{
  group::{ff::Field, Group},
  Ciphersuite, Ed25519,
};
use dalek_ff_group::EdwardsPoint;
use ec_divisors_contest::{check_init_contest, check_init_ref, run_bench_contest, run_bench_ref};
use rand_core::OsRng;

use criterion::{criterion_group, criterion_main, Criterion};

// Benchmark the reference implementation and the implementation in
// ../ec-divisors-contest-src
fn bench_scalar_mul_divisors(c: &mut Criterion) {
  let mut group = c.benchmark_group("ec-divisors");
  group.measurement_time(core::time::Duration::from_secs(60));

  let point = EdwardsPoint::generator();
  let rand_scalar = <Ed25519 as Ciphersuite>::F::random(&mut OsRng);

  check_init_ref(&point, &rand_scalar);
  check_init_contest(&point, &rand_scalar);

  // Run the benchmark for the reference implementation
  group.bench_function("reference-impl", |b| {
    b.iter_with_large_drop(|| run_bench_ref(&point, &rand_scalar))
  });

  // Run the benchmark for the contest implementation
  group.bench_function("contest-impl", |b| {
    b.iter_with_large_drop(|| run_bench_contest(&point, &rand_scalar))
  });

  group.finish();
}

criterion_group!(benches, bench_scalar_mul_divisors);
criterion_main!(benches);
