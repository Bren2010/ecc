use std::fmt;
use std::default::Default;
use std::num::FromStrRadix;
use num::{BigUint, Zero, One};
use num::bigint::ToBigUint;

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
pub trait Field : Clone + Default { // To do:  Implement binary fields.
    fn factors(&self) -> (uint, Vec<int>);
    fn modulus(&self) -> BigUint;
    fn reduce(&self, z: BigUint) -> FieldElem<Self>;
    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<Self>;
}

#[deriving(Clone)]
pub struct FieldElem<F: Field> {
    pub limbs: BigUint
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
        let f: F = Default::default();
        f.reduce(self.limbs + other.limbs)
    }
}

impl<F: Field> Sub<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn sub(&self, other: &FieldElem<F>) -> FieldElem<F> {
        let f: F = Default::default();
        f.reduce(self.limbs + (-*other).limbs)
    }
}

impl<F: Field> Mul<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn mul(&self, other: &FieldElem<F>) -> FieldElem<F> {
        let f: F = Default::default();
        f.reduce(self.limbs * other.limbs)
    }
}

impl<F: Field> Div<FieldElem<F>, FieldElem<F>> for FieldElem<F> {
    fn div(&self, other: &FieldElem<F>) -> FieldElem<F> {
        *self * other.invert()
    }
}

impl<F: Field> Neg<FieldElem<F>> for FieldElem<F> {
    fn neg(&self) -> FieldElem<F> {
        let f: F = Default::default();

        f.reduce(f.modulus() - f.reduce(self.limbs.clone()).limbs)
    }
}

impl<F: Field> Zero for FieldElem<F> {
    fn zero() -> FieldElem<F> {
        let f: F = Default::default();
        f.reduce(0i.to_biguint().unwrap())
    }

    fn is_zero(&self) -> bool { self.limbs.is_zero() }
}

impl<F: Field> One for FieldElem<F> {
    fn one() -> FieldElem<F> {
        let f: F = Default::default();
        f.reduce(1i.to_biguint().unwrap())
    }
}

impl<F: Field> FieldElem<F> {
    pub fn serialize(&self) -> Vec<uint> {
        let f: F = Default::default();
        let mask = 255u.to_biguint().unwrap();
        let m = f.modulus().bits() / 8;
        let mut out: Vec<uint> = Vec::new();

        for i in range(0, m) {
            out.push(((self.limbs >> ((m - i - 1) * 8)) & mask).to_uint().unwrap())
        }

        out
    }

    pub fn pow(&self, cand_exp: &BigUint) -> FieldElem<F> {
        // Montgomery Ladder.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();
        let f: F = Default::default();

        let order: BigUint = f.modulus() - One::one();
        let mut exp: BigUint = *cand_exp + order;

        if exp.bits() == order.bits() { exp = exp + order; }

        let m = exp.bits() + 1;
        let mut r0 = f.reduce(one.clone());
        let mut r1 = (*self).clone();

        for i in range(0u, m) {
            if ((one << (m - i - 1)) & exp) == zer {
                r1 = r0 * r1;
                r0 = r0 * r0;
            } else {
                r0 = r0 * r1;
                r1 = r1 * r1;
            }
        }

        r0
    }

    pub fn sqrt(&self) -> (FieldElem<F>, FieldElem<F>) {
        let f: F = Default::default();

        let exp = (f.modulus() + One::one()) >> 2;
        let out = self.pow(&exp);

        (out.clone(), -out)
    }

    pub fn invert(&self) -> FieldElem<F> {
        // Modulus must be odd and relatively prime to self.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();
        let f: F = Default::default();

        let p = f.modulus();
        let mut tmp: BigUint;

        let mut a: BigUint = One::one();
        let mut b: BigUint = Zero::zero();
        let mut x: BigUint = self.limbs.clone();
        let mut y: BigUint = f.modulus();

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

        f.reduce(b)
    }
}

impl<F: Field> fmt::Show for FieldElem<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.limbs)
    }
}

fn reduce(modulus: (uint, Vec<int>), z: BigUint) -> BigUint {
    let (size, reducs) = modulus;

    let u = z >> size;
    let v = z ^ (u << size);

    let mut out = v;

    for part in reducs.iter() {
        if *part > 0 {
            out = out - (u << part.to_uint().unwrap())
        } else {
            out = out + (u << (*part * -1).to_uint().unwrap())
        }
    }

    out
}

fn unserialize(s: &Vec<uint>) -> BigUint {
    let mut out: BigUint = Zero::zero();

    for x in s.iter() {
        out = (out << 8) ^ (x.to_biguint().unwrap())
    }

    out
}

// -------------------------------------------------------------------------
// Implementations of specific fields.
// -------------------------------------------------------------------------
#[deriving(Clone, Default)]
pub struct P192;
impl Field for P192 {
    fn factors(&self) -> (uint, Vec<int>) { (192, vec![-64, 0]) }
    fn modulus(&self) -> BigUint {
        let one: BigUint = One::one();

        (one << 192) - (one << 64) - one
    }

    fn reduce(&self, z: BigUint) -> FieldElem<P192> {
        if z.bits() > 192 { self.reduce(reduce(self.factors(), z)) }
        else if z.bits() == 192 {
            let m = self.modulus();

            if z >= m { self.reduce(z - m) }
            else { FieldElem { limbs: z } }
        } else { FieldElem { limbs: z } }
    }

    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<P192> {
        self.reduce(unserialize(s))
    }
}

#[deriving(Clone, Default)]
pub struct R192;
impl Field for R192 {
    fn factors(&self) -> (uint, Vec<int>) { (192, vec![-64, 0]) }
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("6277101735386680763835789423176059013767194773182842284081", 10).unwrap()
    }

    fn reduce(&self, z: BigUint) -> FieldElem<R192> {
        FieldElem { limbs: z % self.modulus() }
    }

    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<R192> {
        self.reduce(unserialize(s))
    }
}

#[deriving(Clone, Default)]
pub struct P256;
impl Field for P256 {
    fn factors(&self) -> (uint, Vec<int>) { (256, vec![-224, 192, 96, 0]) }
    fn modulus(&self) -> BigUint {
        let one: BigUint = One::one();

        (one << 256) - (one << 224) + (one << 192) + (one << 96) - one
    }

    fn reduce(&self, z: BigUint) -> FieldElem<P256> {
        if z.bits() > 256 { self.reduce(reduce(self.factors(), z)) }
        else if z.bits() == 256 {
            let m = self.modulus();

            if z >= m { self.reduce(z - m) }
            else { FieldElem { limbs: z } }
        } else { FieldElem { limbs: z } }
    }

    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<P256> {
        self.reduce(unserialize(s))
    }
}

#[deriving(Clone, Default)]
pub struct R256;
impl Field for R256 {
    fn factors(&self) -> (uint, Vec<int>) { (256, vec![-224, 192, 96, 0]) }
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("115792089210356248762697446949407573529996955224135760342422259061068512044369", 10).unwrap()
    }

    fn reduce(&self, z: BigUint) -> FieldElem<R256> {
        FieldElem { limbs: z % self.modulus() }
    }

    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<R256> {
        self.reduce(unserialize(s))
    }
}

#[deriving(Clone, Default)]
pub struct P521;
impl Field for P521 {
    fn factors(&self) -> (uint, Vec<int>) { (521, vec![0]) }
    fn modulus(&self) -> BigUint {
        let one: BigUint = One::one();

        (one << 521) - one
    }

    fn reduce(&self, z: BigUint) -> FieldElem<P521> {
        if z.bits() > 521 { self.reduce(reduce(self.factors(), z)) }
        else if z.bits() == 521 {
            let m = self.modulus();

            if z >= m { self.reduce(z - m) }
            else { FieldElem { limbs: z } }
        } else { FieldElem { limbs: z } }
    }

    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<P521> {
        self.reduce(unserialize(s))
    }
}

#[deriving(Clone, Default)]
pub struct R521;
impl Field for R521 {
    fn factors(&self) -> (uint, Vec<int>) { (521, vec![0]) }
    fn modulus(&self) -> BigUint {
        FromStrRadix::from_str_radix("6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449", 10).unwrap()
    }

    fn reduce(&self, z: BigUint) -> FieldElem<R521> {
        FieldElem { limbs: z % self.modulus() }
    }

    fn unserialize(&self, s: &Vec<uint>) -> FieldElem<R521> {
        self.reduce(unserialize(s))
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
        let p: FieldElem<P192> = FieldElem { limbs: 3i.to_biguint().unwrap() };

        let x = p.pow(&(300i.to_biguint().unwrap())).limbs;
        let y: Option<BigUint> = FromStrRadix::from_str_radix("301636b368b805e333093013fb588bb4f3bc9fd0578b2e6", 16);

        assert_eq!(x, y.unwrap())
    }

    #[bench]
    fn bench_exponentiation_p192(b: &mut Bencher) {
        let p: FieldElem<P192> = FieldElem {
            limbs: FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap()
        };
        let exp = 3000i.to_biguint().unwrap();

        b.iter(|| { p.pow(&exp) })
    }

    #[bench]
    fn bench_exponentiation_p192_right_zeroes(b: &mut Bencher) {
        let p: FieldElem<P192> = FieldElem {
            limbs: FromStrRadix::from_str_radix("7e0000000000000000000000000000000000000000000000", 16).unwrap()
        };
        let exp = 3000i.to_biguint().unwrap();

        b.iter(|| { p.pow(&exp) })
    }

    #[bench]
    fn bench_exponentiation_p192_left_zeroes(b: &mut Bencher) {
        let p: FieldElem<P192> = FieldElem {
            limbs: 3i.to_biguint().unwrap()
        };
        let exp = 3000i.to_biguint().unwrap();

        b.iter(|| { p.pow(&exp) })
    }

    #[bench]
    fn bench_exponentiation_p521(b: &mut Bencher) {
        let p: FieldElem<P521> = FieldElem {
            limbs: FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap(),
        };
        let exp = 3000i.to_biguint().unwrap();

        b.iter(|| { p.pow(&exp) })
    }
}
