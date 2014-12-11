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
}

#[deriving(Clone)]
pub struct FieldElem<F: Field> {
    pub limbs: BigUint,
    pub field: F
}


impl<F: Field> PartialEq for FieldElem<F> {
    fn eq(&self, other: &FieldElem<F>) -> bool {
        if self.limbs.bits() != other.limbs.bits() {
            false
        } else {
            let one: BigUint = One::one();
            let m = self.limbs.bits() - 1;
            let mut ok = true;

            for i in range(0u, m) {
                ok &= (self.limbs & (one << i)) ^ (other.limbs & (one << i)) != (one << i);
            }

            ok
        }
    }
}

// Notes:  Ordering is omitted because prime fields have no norm.

impl<F: Field> Add<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn add(&self, other: &FieldElem<F>) -> FieldElem<F> {
        FieldElem {
            limbs: (self.limbs + other.limbs) % self.field.modulus(),
            field: self.field.clone()
        }
    }
}

impl<F: Field> Sub<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn sub(&self, other: &FieldElem<F>) -> FieldElem<F> {
        FieldElem {
            limbs: (self.limbs + (-*other).limbs) % self.field.modulus(),
            field: self.field.clone()
        }
    }
}

impl<F: Field> Mul<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn mul(&self, other: &FieldElem<F>) -> FieldElem<F> {
        FieldElem {
            limbs: (self.limbs * other.limbs) % self.field.modulus(),
            field: self.field.clone()
        }
    }
}

impl<F: Field> Div<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn div(&self, other: &FieldElem<F>) -> FieldElem<F> {
        *self * other.invert()
    }
}

impl<F: Field> Neg<FieldElem<F>> for FieldElem<F> {
    fn neg(&self) -> FieldElem<F> {
        FieldElem {
            limbs: (self.field.modulus() - self.limbs) % self.field.modulus(),
            field: self.field.clone()
        }
    }
}

impl<F: Field> FieldElem<F> {
    pub fn pow(&self, exp: &BigUint) -> FieldElem<F> {
        // Montgomery Ladder.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();

        let m = exp.bits() + 1;
        let mut r0 = FieldElem {
            limbs: one.clone(),
            field: self.field.clone()
        };
        let mut r1 = FieldElem {
            limbs: self.limbs.clone(),
            field: self.field.clone()
        };

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

        FieldElem {
            limbs: b,
            field: self.field.clone()
        }
    }

    // To do:  Find a way to actually implement Zero instead of mimic it.
    pub fn zero(&self) -> FieldElem<F> {
        FieldElem {
            limbs: 0i.to_biguint().unwrap(),
            field: self.field.clone()
        }
    }

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
}

#[deriving(Clone)]
pub struct R192;
impl Field for R192 {
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("6277101735386680763835789423176059013767194773182842284081", 10).unwrap()
    }
}

// -------------------------------------------------------------------------
// Unit Tests
// -------------------------------------------------------------------------
#[cfg(test)]
mod test {
    use std::num::FromStrRadix;
    use num::bigint::{ToBigUint, BigUint};
    use super::{FieldElem, P192};

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
}
