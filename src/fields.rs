use std::fmt;
use std::num::FromStrRadix;
use num::{BigUint, Zero, One};
use num::bigint::{ToBigUint};

// To do:  Follow Rust pattern of implementing .to_field() functions on a
//     handful of standard types.
// -------------------------------------------------------------------------
// Generic Field and FieldElement descriptions / implementations.
// -------------------------------------------------------------------------
//
// Prime field elements are constructed with a BigUint and a choice of
// field.
//
// All the field operations can be computed with the normal unary/binary
// operators.  Equality and exponentiation are constant time.
//
// Unary:  Negation.
// Binary:  Addition, Subtraction, Multiplication, Division.
//
// - Element `x` can be inverted by calling `x.invert()`.
// - Powers are computed by calling `x.pow(exp)`, where exp is a BigUint.
pub trait Field : Clone { // To do:  Implement binary fields.
    fn modulus(&self) -> BigUint;
    fn reduce(&self, z: BigUint) -> FieldElem<Self>;
    fn negate(&self, z: &FieldElem<Self>) -> FieldElem<Self>;
}

#[deriving(Clone)]
pub struct FieldElem<F: Field> {
    pub limbs: BigUint,
    pub field: F
}


impl<F: Field> PartialEq for FieldElem<F> {
    fn eq(&self, other: &FieldElem<F>) -> bool {
        if self.limbs.bits() != other.limbs.bits() { false }
        else { (self.limbs ^ other.limbs).bits() == 0 }
    }
}

// Notes:  Ordering is omitted because prime fields have no norm.

impl<F: Field> Add<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn add(&self, other: &FieldElem<F>) -> FieldElem<F> {
        self.field.reduce(self.limbs + other.limbs)
    }
}

impl<F: Field> Sub<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn sub(&self, other: &FieldElem<F>) -> FieldElem<F> {
        self.field.reduce(self.limbs + (-*other).limbs)
    }
}

impl<F: Field> Mul<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn mul(&self, other: &FieldElem<F>) -> FieldElem<F> {
        self.field.reduce(self.limbs * other.limbs)
    }
}

impl<F: Field> Div<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn div(&self, other: &FieldElem<F>) -> FieldElem<F> {
        *self * other.invert()
    }
}

impl<F: Field> Neg<FieldElem<F>> for FieldElem<F> {
    fn neg(&self) -> FieldElem<F> { self.field.negate(self) }
}

impl<F: Field> FieldElem<F> {
    pub fn pow(&self, exp: &BigUint) -> FieldElem<F> {
        // Montgomery Ladder.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();

        let m = exp.bits() + 1;
        let mut r0 = self.field.reduce(one.clone());
        let mut r1 = self.clone();

        for i in range(0u, m) {
            if ((one << (m - i - 1)) & *exp) == zer {
                r1 = r0 * r1;
                r0 = r0 * r0;
            } else {
                r0 = r0 * r1;
                r1 = r1 * r1;
            }
        }

        r0
    }

    pub fn invert(&self) -> FieldElem<F> {
        // Modulus must be odd and relatively prime to self.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();

        let p = self.field.modulus();
        let mut tmp: BigUint;

        let mut a: BigUint = One::one();
        let mut b: BigUint = Zero::zero();
        let mut x: BigUint = self.limbs.clone();
        let mut y: BigUint = self.field.modulus();

        while x != zer {
            if x & one == one {
                if x < y { // x < y; Swap everything
                    tmp = x; x = y; y = tmp;
                    tmp = a; a = b; b = tmp;
                }

                x = x - y;
                if a < b { a = a + p; }
                a = a - b;
            }

            x = x >> 1;
            if a & one == one { a = a + p }
            a = a >> 1;
        }

        self.field.reduce(b)
    }

    // To do:  Find a way to actually implement Zero and One instead of mimic it.
    pub fn zero(&self) -> FieldElem<F> { self.field.reduce(0i.to_biguint().unwrap()) }
    pub fn one(&self) -> FieldElem<F> { self.field.reduce(1i.to_biguint().unwrap()) }
    pub fn is_zero(&self) -> bool { self.limbs.is_zero() }
}

impl<F: Field> fmt::Show for FieldElem<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.limbs)
    }
}


// -------------------------------------------------------------------------
// Implementations of specific fields.
// -------------------------------------------------------------------------
#[deriving(Clone)]
pub struct P192;
impl Field for P192 {
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("6277101735386680763835789423207666416083908700390324961279", 10).unwrap()
    }

    fn reduce(&self, z: BigUint) -> FieldElem<P192> {
        FieldElem {
            limbs: z % self.modulus(),
            field: P192
        }
    }

    fn negate(&self, z: &FieldElem<P192>) -> FieldElem<P192> {
        self.reduce(self.modulus() - self.reduce(z.limbs.clone()).limbs)
    }
}

#[deriving(Clone)]
pub struct R192;
impl Field for R192 {
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("6277101735386680763835789423176059013767194773182842284081", 10).unwrap()
    }

    fn reduce(&self, z: BigUint) -> FieldElem<R192> {
        FieldElem { limbs: z % self.modulus(), field: R192 }
    }

    fn negate(&self, z: &FieldElem<R192>) -> FieldElem<R192> {
        self.reduce(self.modulus() - self.reduce(z.limbs.clone()).limbs)
    }
}

#[deriving(Clone)]
pub struct P521;
impl Field for P521 {
    fn modulus(&self) -> BigUint {
        let one: BigUint = One::one();

        (one << 521) - one
    }

    fn reduce(&self, z: BigUint) -> FieldElem<P521> {
        if z.bits() > 521 {
            let u = z >> 521;
            let v = z ^ (u << 521);

            self.reduce(u + v)
        } else {
            FieldElem {
                limbs: z,
                field: P521
            }
        }
    }

    fn negate(&self, z: &FieldElem<P521>) -> FieldElem<P521> {
        self.reduce(self.modulus() ^ self.reduce(z.limbs.clone()).limbs)
    }
}

#[deriving(Clone)]
pub struct R521;
impl Field for R521 {
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449", 10).unwrap()
    }

    fn reduce(&self, z: BigUint) -> FieldElem<R521> {
        FieldElem { limbs: z % self.modulus(), field: R521 }
    }

    fn negate(&self, z: &FieldElem<R521>) -> FieldElem<R521> {
        self.reduce(self.modulus() - self.reduce(z.limbs.clone()).limbs)
    }
}

// -------------------------------------------------------------------------
// Unit Tests
// -------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    extern crate test;
    use self::test::Bencher;

    use std::num::FromStrRadix;
    use num::bigint::{ToBigUint, BigUint};
    use super::{FieldElem, P192, P521};

    #[test]
    fn large_exponentiation() {
        let p = FieldElem {
            limbs: 3i.to_biguint().unwrap(),
            field: P192
        };

        let x = p.pow(&(300i.to_biguint().unwrap())).limbs;
        let y: Option<BigUint> = FromStrRadix::from_str_radix("301636b368b805e333093013fb588bb4f3bc9fd0578b2e6", 16);

        assert_eq!(x, y.unwrap())
    }

    #[bench]
    fn bench_exponentiation_p192(b: &mut Bencher) {
        let p = FieldElem {
            limbs: 3i.to_biguint().unwrap(),
            field: P192
        };
        let exp = 3000i.to_biguint().unwrap();

        b.iter(|| { p.pow(&exp) })
     }

     #[bench]
     fn bench_exponentiation_p521(b: &mut Bencher) {
         let p = FieldElem {
             limbs: 3i.to_biguint().unwrap(),
             field: P521
         };
         let exp = 3000i.to_biguint().unwrap();

         b.iter(|| { p.pow(&exp) })
      }
}
