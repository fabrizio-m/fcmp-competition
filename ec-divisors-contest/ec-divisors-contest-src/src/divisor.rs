use core::ops::{Div, Mul};
use ff::PrimeField;
use std::rc::Rc;

/// Divisor of form f(x,y) = A(x) - yB(x), with A and B
/// represented as enough evaluations for their degree.
struct Divisor<F: PrimeField> {
    a: Evals<F>,
    b: Evals<F>,
    // to substitute y^2
    modulus: Rc<Evals<F>>,
}

/// Represented as coefficients as it can be efficiently evaluated
/// in O(n) additions.
struct SmallDivisor<F: PrimeField> {
    // (a,b) for ax + b
    a: (F, F),
    b: F,
}

impl<F: PrimeField> Divisor<F> {
    fn ab(&self, i: usize) -> (F, F) {
        let a = self.a.evals[i];
        let b = self.b.evals[i];
        (a, b)
    }

    fn ab_mut(&mut self, i: usize) -> (&mut F, &mut F) {
        let a = &mut self.a.evals[i];
        let b = &mut self.b.evals[i];
        (a, b)
    }

    /// New degree after mul.
    fn new_degree(&self, other: &Self) -> (usize, usize) {
        // f1 * f2 = A1A2 - y(A1B2 + A2B1) + (x^3 + ax + b) B1B2
        // A = A1A2 + (x^3 + ax + b) B1B2
        // B = A1B2 + A2B1
        // deg(A) = max(A1 + A2, 3 + B1 + B2)
        // deg(B) = max(A1 + B2, A2 + B1)
        let (a1, b1) = (self.a.degree, self.b.degree);
        let (a2, b2) = (other.a.degree, other.b.degree);
        let a = (a1 + a2).max(3 + b1 + b2);
        let b = (a1 + b2).max(a2 + b1);
        (a, b)
    }

    /// Remove 2 points by dividing by (x - x1) * (x - x2)
    fn remove_diff(self, x1: F, x2: F) -> Self {
        assert_eq!(self.a.len(), self.b.len());
        let mut denominator = Vec::with_capacity(self.a.len());
        for i in 0..self.a.len() {
            let x = F::from(i as u64);
            denominator.push((x - x1) * (x - x2));
        }
        let denominator = Evals {
            evals: denominator,
            degree: 2,
        };
        self / denominator
    }

    fn merge(divisors: [Self; 2], small: SmallDivisor<F>, denom: (F, F)) -> Self {
        let [d1, d2] = divisors;
        let numerator = d1 * &d2;
        //TODO: second operand can be reused for scratch space
        let numerator = numerator * small;
        let (x1, x2) = denom;
        numerator.remove_diff(x1, x2)
    }
}

impl<F: PrimeField> SmallDivisor<F> {
    fn new(a: (F, F), b: F) -> Self {
        Self { a, b }
    }
}

struct Evals<F: PrimeField> {
    evals: Vec<F>,
    degree: usize,
}

impl<F: PrimeField> Evals<F> {
    fn len(&self) -> usize {
        self.evals.len()
    }
}

/*
impl<F: PrimeField> Add<&Self> for Evals<F> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        debug_assert_eq!(self.evals.len(), rhs.evals.len());
        let Self { mut evals, degree } = self;
        let degree = degree.max(rhs.degree);
        for e in evals.iter_mut().zip(rhs.evals.iter()) {
            let (l, r) = e;
            *l += r;
        }
        Self { evals, degree }
    }
}

impl<F: PrimeField> Mul<&Self> for Evals<F> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        debug_assert_eq!(self.evals.len(), rhs.evals.len());
        let Self { mut evals, degree } = self;
        let degree = degree + rhs.degree;
        debug_assert!(evals.len() > degree);
        for e in evals.iter_mut().zip(rhs.evals.iter()) {
            let (l, r) = e;
            *l *= r;
        }
        Self { evals, degree }
    }
}

impl<F: PrimeField> Add<F> for Evals<F> {
    type Output = Self;

    fn add(mut self, rhs: F) -> Self::Output {
        for eval in self.evals.iter_mut() {
            *eval += rhs;
        }
        self
    }
}


impl<F: PrimeField> Mul<F> for Evals<F> {
    type Output = Self;

    fn mul(mut self, rhs: F) -> Self::Output {
        for eval in self.evals.iter_mut() {
            *eval *= rhs;
        }
        self
    }
}

*/
impl<F: PrimeField> Mul<&Self> for Divisor<F> {
    type Output = Self;

    fn mul(mut self, rhs: &Self) -> Self::Output {
        debug_assert_eq!(self.a.len(), rhs.a.len());
        debug_assert_eq!(self.b.len(), rhs.b.len());
        debug_assert_eq!(self.a.len(), self.b.len());
        let len = self.a.len();

        let new_degree = self.new_degree(&rhs);
        // f1 * f2 = A1A2 - y(A1B2 + A2B1) + y^2 B1B2
        // f1 * f2 = A1A2 - y(A1B2 + A2B1) + (x^3 + ax + b) B1B2
        // (A1+B1)(A2+B2)
        // A1A2 + A1B2 + B1A2 + B1B2
        for i in 0..len {
            let modulus = self.modulus.evals[i];
            let (a1, b1) = self.ab_mut(i);
            let (a2, b2) = rhs.ab(i);
            let a1a2 = *a1 * a2;
            let b1b2 = *b1 * b2;
            // (A1+B1)(A2+B2)
            let cross = (*a1 + *b1) * (a2 + b2);
            let b = cross - (a1a2 + b1b2);
            let a = a1a2 + b1b2 * modulus;
            *a1 = a;
            *b1 = b;
        }
        let (a, b) = new_degree;
        self.a.degree = a;
        self.b.degree = b;
        self
    }
}

impl<F: PrimeField> Mul<SmallDivisor<F>> for Divisor<F> {
    type Output = Self;

    fn mul(mut self, rhs: SmallDivisor<F>) -> Self::Output {
        debug_assert_eq!(self.a.len(), self.b.len());
        let len = self.a.len();

        // constant term for x = 0
        let mut a2 = rhs.a.1;
        let b2 = rhs.b;
        for i in 0..len {
            let modulus = self.modulus.evals[i];
            let (a1, b1) = self.ab_mut(i);
            let a1a2 = *a1 * a2;
            let b1b2 = *b1 * b2;
            // (A1+B1)(A2+B2)
            let cross = (*a1 + *b1) * (a2 + b2);
            let b = cross - (a1a2 + b1b2);
            let a = a1a2 + b1b2 * modulus;
            a2 += rhs.a.0;
            *a1 = a;
            *b1 = b;
        }
        let da = self.a.degree;
        let db = self.a.degree;
        self.a.degree = da.max(db + 3);
        self.b.degree = da.max(1 + db);
        self
    }
}

impl<F: PrimeField> Div<Evals<F>> for Divisor<F> {
    type Output = Self;

    fn div(mut self, mut rhs: Evals<F>) -> Self::Output {
        let len = rhs.len();
        debug_assert_eq!(self.a.len(), self.b.len());
        debug_assert_eq!(self.a.len(), len);
        // TODO: maybe reuse
        let mut scratch_space: Vec<F> = rhs.evals.clone();
        ff::BatchInverter::invert_with_external_scratch(&mut rhs.evals, &mut scratch_space);
        for i in 0..len {
            let denom = rhs.evals[i];
            self.a.evals[i] *= denom;
            self.b.evals[i] *= denom;
        }
        self.a.degree -= rhs.degree;
        self.b.degree -= rhs.degree;
        self
    }
}
