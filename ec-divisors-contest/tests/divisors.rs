use ec_divisors_contest::{check_init_contest, check_init_ref, run_bench_contest, run_bench_ref};

use ciphersuite::{
    group::{ff::Field, Group},
    Ciphersuite, Ed25519,
};
use dalek_ff_group::EdwardsPoint;
use rand_core::OsRng;

// Test that the implementation in ../ec-divisors-contest-src produces the same
// result as the reference implementation
#[test]
fn divisors_contest_test() {
    static N_ITERS: usize = 20;

    for i in 0..N_ITERS {
        println!("Testing with random scalar {} / {}", i + 1, N_ITERS);

        let point = EdwardsPoint::generator();
        let scalar = <Ed25519 as Ciphersuite>::F::random(&mut OsRng);

        check_init_ref(&point, &scalar);
        check_init_contest(&point, &scalar);

        // Get divisors
        let ref_res = run_bench_ref(&point, &scalar).0;
        let res = run_bench_contest(&point, &scalar).0;

        assert_eq!(ref_res.y_coefficients, res.y_coefficients);
        assert_eq!(ref_res.yx_coefficients, res.yx_coefficients);
        assert_eq!(ref_res.x_coefficients, res.x_coefficients);
        assert_eq!(ref_res.zero_coefficient, res.zero_coefficient);
    }
}
