use std::fmt;
use std::default::Default;
use std::num::FromStrRadix;
use num::{BigUint, Zero, One};
use num::bigint::ToBigUint;

use fields::{Field, FieldElem, P192, R192, P256, R256, P521, R521};

#[allow(non_snake_case)]
// Weierstrass curve with large characteristic field.
pub trait Curve<F: Field, G: Field>: Clone + Default { // To do:  Implement binary curves.
    fn AN3(&self) -> bool; // If A = -3
    fn unserialize(&self, s: &Vec<uint>) -> AffinePoint<Self, F, G>;

    // Standard curve paraemters:
    fn A(&self) -> FieldElem<F>;
    fn B(&self) -> FieldElem<F>;
    fn G(&self) -> AffinePoint<Self, F, G>;
}

#[deriving(Clone)]
pub struct AffinePoint<C: Curve<F, G>, F: Field, G: Field> {
    pub x: FieldElem<F>,
    pub y: FieldElem<F>,
}

#[deriving(Clone)]
pub struct JacobianPoint<C: Curve<F, G>, F: Field, G: Field> {
    pub x: FieldElem<F>,
    pub y: FieldElem<F>,
    pub z: FieldElem<F>,
}

// Constant-time exponentiation/scalar-multiplication.
fn pow<T: Zero + Add<T, T> + Clone>(this: &T, cand_exp: &BigUint, order: &BigUint) -> T {
    // Montgomery Ladder.
    let zer: BigUint = Zero::zero();
    let one: BigUint = One::one();
    let mut exp: BigUint = *cand_exp + *order;

    if exp.bits() == order.bits() { exp = exp + *order; }

    let m = exp.bits() + 1;
    let mut r0: T = Zero::zero();
    let mut r1: T = (*this).clone();

    for i in range(0u, m) {
        if ((one << (m - i - 1)) & exp) == zer {
            r1 = r0 + r1;
            r0 = r0 + r0;
        } else {
            r0 = r0 + r1;
            r1 = r1 + r1;
        }
    }

    r0
}

impl<C: Curve<F, G>, F: Field, G: Field> PartialEq for AffinePoint<C, F, G> {
    fn eq(&self, other: &AffinePoint<C, F, G>) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> PartialEq for JacobianPoint<C, F, G> {
    fn eq(&self, other: &JacobianPoint<C, F, G>) -> bool {
        self.to_affine() == other.to_affine()
    }
}

// Notes:  Ordering is omitted because elliptic curve groups have no norm.

impl<C: Curve<F, G>, F: Field, G: Field> Add<AffinePoint<C, F, G>, AffinePoint<C, F, G>> for AffinePoint<C, F, G> {
    fn add(&self, other: &AffinePoint<C, F, G>) -> AffinePoint<C, F, G> {
        if (*self).is_zero() {
            (*other).clone()
        } else if (*other).is_zero() {
            (*self).clone()
        } else if self.x != other.x {
            let m = (other.y - self.y) / (other.x - self.x);

            let x3 = (m * m) - self.x - other.x;
            let y3 = m * (self.x - x3) - self.y;

            AffinePoint { x: x3, y: y3 }
        } else if self.y != other.y || self.y.limbs.bits() == 0 {
            Zero::zero()
        } else {
            let c: C = Default::default();
            let two: FieldElem<F> = FieldElem {
                limbs: 2i.to_biguint().unwrap(),
            };
            let three: FieldElem<F> = FieldElem {
                limbs: 3i.to_biguint().unwrap(),
            };

            let m = ((three * (self.x * self.x)) + c.A()) / (two * self.y);

            let x3 = (m * m) - (self.x * two);
            let y3 = m * (self.x - x3) - self.y;

            AffinePoint { x: x3, y: y3 }
        }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> Add<JacobianPoint<C, F, G>, JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn add(&self, other: &JacobianPoint<C, F, G>) -> JacobianPoint<C, F, G> {
        if (*self).is_zero() {
            (*other).clone()
        } else if (*other).is_zero() {
            (*self).clone()
        } else if self.x == other.x && self.y == other.y && self.z == other.z {
            self.double()
        } else {
            let z12 = self.z * self.z;
            let z13 = z12 * self.z;
            let z22 = other.z * other.z;
            let z23 = z22 * other.z;

            let u1 = self.x * z22;
            let u2 = other.x * z12;
            let s1 = self.y * z23;
            let s2 = other.y * z13;

            if u1 == u2 {
                if s1 != s2 { // P1 = +/- P2
                    Zero::zero()
                } else {
                    self.double()
                }
            } else {
                let h = u2 - u1;
                let h2 = h * h;
                let h3 = h2 * h;

                let r = s2 - s1;
                let u1h2 = u1 * h2;

                let x3 = (r * r) - h3 - u1h2 - u1h2;
                let y3 = r * (u1h2 - x3) - (s1 * h3);
                let z3 = h * self.z * other.z;

                JacobianPoint { x: x3, y: y3, z: z3 }
            }
        }
    }
}

fn add_jacobian_to_affine<C: Curve<F, G>, F: Field, G: Field>(this: &JacobianPoint<C, F, G>, other: &AffinePoint<C, F, G>) -> JacobianPoint<C, F, G> {
    if this.is_zero() {
        (*other).to_jacobian()
    } else if other.is_zero() {
        (*this).clone()
    } else {
        let z2 = this.z * this.z;
        let c = (other.x * z2) - this.x;

        if c.is_zero() {
            if (other.y * z2 * this.z) == this.y { // Same point
                this.double()
            } else { // Negatives
                Zero::zero()
            }
        } else {
            let d = (other.y * z2 * this.z) - this.y;
            let c2 = c * c;

            let x3 = (d * d) - (c2 * c) - (c2 * (this.x + this.x));
            let y3 = d * ((this.x * c2) - x3) - (this.y * c2 * c);
            let z3 = this.z * c;

            JacobianPoint { x: x3, y: y3, z: z3 }
        }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> Add<JacobianPoint<C, F, G>, JacobianPoint<C, F, G>> for AffinePoint<C, F, G> {
    fn add(&self, other: &JacobianPoint<C, F, G>) -> JacobianPoint<C, F, G> { add_jacobian_to_affine(other, self) }
}

impl<C: Curve<F, G>, F: Field, G: Field> Add<AffinePoint<C, F, G>, JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn add(&self, other: &AffinePoint<C, F, G>) -> JacobianPoint<C, F, G> { add_jacobian_to_affine(self, other) }
}

impl<C: Curve<F, G>, F: Field, G: Field> Sub<AffinePoint<C, F, G>, AffinePoint<C, F, G>> for AffinePoint<C, F, G> {
    fn sub(&self, other: &AffinePoint<C, F, G>) -> AffinePoint<C, F, G> { *self + (-*other) }
}

impl<C: Curve<F, G>, F: Field, G: Field> Sub<JacobianPoint<C, F, G>, JacobianPoint<C, F, G>> for AffinePoint<C, F, G> {
    fn sub(&self, other: &JacobianPoint<C, F, G>) -> JacobianPoint<C, F, G> { *self + (-*other) }
}

impl<C: Curve<F, G>, F: Field, G: Field> Sub<JacobianPoint<C, F, G>, JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn sub(&self, other: &JacobianPoint<C, F, G>) -> JacobianPoint<C, F, G> { *self + (-*other) }
}

impl<C: Curve<F, G>, F: Field, G: Field> Sub<AffinePoint<C, F, G>, JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn sub(&self, other: &AffinePoint<C, F, G>) -> JacobianPoint<C, F, G> { *self + (-*other) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<FieldElem<G>, AffinePoint<C, F, G>> for AffinePoint<C, F, G> {
    fn mul(&self, other: &FieldElem<G>) -> AffinePoint<C, F, G> { self.pow(&other.limbs) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<AffinePoint<C, F, G>, AffinePoint<C, F, G>> for FieldElem<G> {
    fn mul(&self, other: &AffinePoint<C, F, G>) -> AffinePoint<C, F, G> { other.pow(&self.limbs) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<FieldElem<G>, JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn mul(&self, other: &FieldElem<G>) -> JacobianPoint<C, F, G> { self.pow(&other.limbs) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<JacobianPoint<C, F, G>, JacobianPoint<C, F, G>> for FieldElem<G> {
    fn mul(&self, other: &JacobianPoint<C, F, G>) -> JacobianPoint<C, F, G> { other.pow(&self.limbs) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<BigUint, AffinePoint<C, F, G>> for AffinePoint<C, F, G> {
    fn mul(&self, other: &BigUint) -> AffinePoint<C, F, G> { self.pow(other) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<AffinePoint<C, F, G>, AffinePoint<C, F, G>> for BigUint {
    fn mul(&self, other: &AffinePoint<C, F, G>) -> AffinePoint<C, F, G> { other.pow(self) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<BigUint, JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn mul(&self, other: &BigUint) -> JacobianPoint<C, F, G> { self.pow(other) }
}

impl<C: Curve<F, G>, F:Field, G: Field> Mul<JacobianPoint<C, F, G>, JacobianPoint<C, F, G>> for BigUint {
    fn mul(&self, other: &JacobianPoint<C, F, G>) -> JacobianPoint<C, F, G> { other.pow(self) }
}

impl<C: Curve<F, G>, F: Field, G: Field> Neg<AffinePoint<C, F, G>> for AffinePoint<C, F, G> {
    fn neg(&self) -> AffinePoint<C, F, G> {
        AffinePoint { x: self.x.clone(), y: -(self.y.clone()) }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> Neg<JacobianPoint<C, F, G>> for JacobianPoint<C, F, G> {
    fn neg(&self) -> JacobianPoint<C, F, G> {
        JacobianPoint { x: self.x.clone(), y: -(self.y.clone()), z: self.z.clone() }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> Zero for AffinePoint<C, F, G> {
    fn zero() -> AffinePoint<C, F, G> {
        AffinePoint { x: Zero::zero(), y: Zero::zero() }
    }

    fn is_zero(&self) -> bool { self.x.is_zero() && self.y.is_zero() }
}

impl<C: Curve<F, G>, F: Field, G: Field> Zero for JacobianPoint<C, F, G> {
    fn zero() -> JacobianPoint<C, F, G> {
        JacobianPoint { x: One::one(), y: One::one(), z: Zero::zero() }
    }

    fn is_zero(&self) -> bool { self.z.is_zero() }
}

impl<C: Curve<F, G>, F: Field, G: Field> AffinePoint<C, F, G> {
    pub fn is_valid(&self) -> bool {
        if self.is_zero() {
            true
        } else {
            let c: C = Default::default();
            let y2 = self.y * self.y;
            let x = (self.x * self.x * self.x) + (c.A() * self.x) + c.B();

            y2 == x
        }
    }

    pub fn serialize(&self) -> Vec<uint> {
        let mut out: Vec<uint> = self.x.serialize();

        if self.y.limbs < (-self.y).limbs { out.push(0) }
        else { out.push(1) }

        out
    }

    fn pow(&self, exp: &BigUint) -> AffinePoint<C, F, G> {
        let g: G = Default::default();
        pow(self, exp, &g.modulus())
    }

    pub fn to_jacobian(&self) -> JacobianPoint<C, F, G> {
        let p = JacobianPoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: One::one()
        };

        if self.is_zero() {
            Zero::zero()
        } else {
            p
        }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> JacobianPoint<C, F, G> {
    pub fn is_valid(&self) -> bool {
        if self.is_zero() {
            true
        } else {
            let c: C = Default::default();
            let z4 = self.z * self.z * self.z * self.z;

            let y2 = self.y * self.y;
            let x = (self.x * self.x * self.x) + (c.A() * self.x * z4)
                + (c.B() * z4 * self.z * self.x);

            y2 == x
        }
    }

    fn pow(&self, exp: &BigUint) -> JacobianPoint<C, F, G> { // Replace with generic
        let g: G = Default::default();
        pow(self, exp, &g.modulus())
    }

    pub fn double(&self) -> JacobianPoint<C, F, G> {
        let c: C = Default::default();
        let f: F = Default::default();

        let y2 = self.y * self.y;
        let y4 = y2 * y2;

        let z2 = self.z * self.z;

        let xy2 = self.x * y2;
        let yz = self.y * self.z;

        let v = xy2 + xy2 + xy2 + xy2;
        let w: FieldElem<F>;

        if c.AN3() {
            let tmp = (self.x + z2) * (self.x - z2);
            w = tmp + tmp + tmp;
        } else {
            let x2 = self.x * self.x;
            let z4  = z2 * z2;
            w = x2 + x2 + x2 + (c.A() * z4);
        }

        let v2 = f.reduce(v.limbs << 1);
        let y84 = f.reduce(y4.limbs << 3);

        let xf = (w * w) - v2;
        let yf = (w * (v - xf)) - y84;
        let zf = yz + yz;

        JacobianPoint { x: xf, y: yf, z: zf }
    }

    pub fn to_affine(&self) -> AffinePoint<C, F, G> {
        if self.is_zero() {
            Zero::zero()
        } else {
            let zinv = self.z.invert();

            AffinePoint {
                x: self.x * zinv * zinv,
                y: self.y * zinv * zinv * zinv
            }
        }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> fmt::Show for AffinePoint<C, F, G> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        } else {
            return write!(f, "({}, {})", self.x.limbs, self.y.limbs)
        }
    }
}

impl<C: Curve<F, G>, F: Field, G: Field> fmt::Show for JacobianPoint<C, F, G> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        } else {
            return write!(f, "({} : {} : {})", self.x.limbs, self.y.limbs, self.z.limbs)
        }
    }
}

fn unserialize<C: Curve<F, G>, F: Field, G: Field>(c: &C, s: &Vec<uint>) -> (FieldElem<F>, FieldElem<F>) {
    let f: F = Default::default();

    let mut t = s.clone();
    let sign = t.pop().unwrap();

    let x = f.unserialize(&t);
    let (y1, y2) = ((x * x * x) + (c.A() * x) + c.B()).sqrt();

    if (sign == 0) ^ (y1.limbs < (y2).limbs) {
        (x, y2)
    } else {
        (x, y1)
    }
}

#[deriving(Clone, Default)]
pub struct C192<P192, R192>;
impl Curve<P192, R192> for C192<P192, R192> {
    fn AN3(&self) -> bool { true }
    fn A(&self) -> FieldElem<P192> {
        -(FieldElem { limbs: 3i.to_biguint().unwrap() })
    }

    fn B(&self) -> FieldElem<P192> {
        FieldElem { limbs: FromStrRadix::from_str_radix("64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1", 16).unwrap() }
    }

    fn G(&self) -> AffinePoint<C192<P192, R192>, P192, R192> {
        AffinePoint {
            x: FieldElem { // To do:  Implement FromStrRadix for FieldElem
                limbs: FromStrRadix::from_str_radix("188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012", 16).unwrap()
            },
            y: FieldElem {
                limbs: FromStrRadix::from_str_radix("07192b95ffc8da78631011ed6b24cdd573f977a11e794811", 16).unwrap()
            }
        }
    }

    fn unserialize(&self, s: &Vec<uint>) -> AffinePoint<C192<P192, R192>, P192, R192> {
        let (x, y) = unserialize(self, s);

        AffinePoint { x: x, y: y }
    }
}

#[deriving(Clone, Default)]
pub struct C256<P256, R256>;
impl Curve<P256, R256> for C256<P256, R256> {
    fn AN3(&self) -> bool { true }
    fn A(&self) -> FieldElem<P256> {
        -(FieldElem { limbs: 3i.to_biguint().unwrap() })
    }

    fn B(&self) -> FieldElem<P256> {
        FieldElem { limbs: FromStrRadix::from_str_radix("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b", 16).unwrap() }
    }

    fn G(&self) -> AffinePoint<C256<P256, R256>, P256, R256> {
        AffinePoint {
            x: FieldElem { // To do:  Implement FromStrRadix for FieldElem
                limbs: FromStrRadix::from_str_radix("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16).unwrap()
            },
            y: FieldElem {
                limbs: FromStrRadix::from_str_radix("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16).unwrap()
            }
        }
    }

    fn unserialize(&self, s: &Vec<uint>) -> AffinePoint<C256<P256, R256>, P256, R256> {
        let (x, y) = unserialize(self, s);

        AffinePoint { x: x, y: y }
    }}

#[deriving(Clone, Default)]
pub struct C521<P521, R521>;
impl Curve<P521, R521> for C521<P521, R521> {
    fn AN3(&self) -> bool { true }
    fn A(&self) -> FieldElem<P521> {
        -(FieldElem {
            limbs: 3i.to_biguint().unwrap(),
        })
    }

    fn B(&self) -> FieldElem<P521> {
        FieldElem {
            limbs: FromStrRadix::from_str_radix("051953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00", 16).unwrap()
        }
    }

    fn G(&self) -> AffinePoint<C521<P521, R521>, P521, R521> {
        AffinePoint {
            x: FieldElem { // To do:  Implement FromStrRadix for FieldElem
                limbs: FromStrRadix::from_str_radix("c6858e06b70404e9cd9e3ecb662395b4429c648139053fb521f828af606b4d3dbaa14b5e77efe75928fe1dc127a2ffa8de3348b3c1856a429bf97e7e31c2e5bd66", 16).unwrap()
            },
            y: FieldElem {
                limbs: FromStrRadix::from_str_radix("11839296a789a3bc0045c8a5fb42c7d1bd998f54449579b446817afbd17273e662c97ee72995ef42640c550b9013fad0761353c7086a272c24088be94769fd16650", 16).unwrap()
            }
        }
    }

    fn unserialize(&self, s: &Vec<uint>) -> AffinePoint<C521<P521, R521>, P521, R521> {
        let (x, y) = unserialize(self, s);

        AffinePoint { x: x, y: y }
    }}

// -------------------------------------------------------------------------
// Unit Tests
// -------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    extern crate test;
    use self::test::Bencher;

    use std::num::FromStrRadix;
    use num::{BigUint, One};
    use num::bigint::ToBigUint;

    use fields::{FieldElem, P192, R192, P256, R256, P521, R521};
    use super::{AffinePoint, Curve, C192, C256, C521};

    #[test]
    fn accept_valid_point() {
        let c: C192<P192, R192> = C192;
        assert_eq!(c.G().is_valid(), true)
    }

    #[test]
    fn reject_invalid_point() {
        let c: C192<P192, R192> = C192;
        let p = AffinePoint {
            x: c.G().x,
            y: FieldElem { limbs: c.G().y.limbs + One::one() }
        };

        assert_eq!(p.is_valid(), false)
    }

    #[test]
    fn base_point_field_is_r192() {
        let c: C192<P192, R192> = C192;
        let one: FieldElem<R192> = FieldElem { limbs: One::one() };

        let x: FieldElem<R192> = FieldElem { limbs: 3i.to_biguint().unwrap() };

        let y = x.invert();

        let a = x * c.G();
        let b = y * a;

        assert!(x != y);
        assert!((x * y) == one);
        assert!(c.G() != a);
        assert!(c.G() == b);
    }

    #[test]
    fn affine_point_multiplication() {
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
        let x: BigUint = FromStrRadix::from_str_radix("48162eae1116dbbd5b7a0d9494ff0c9b414a31ce3d8b273f", 16).unwrap();
        let y: BigUint = FromStrRadix::from_str_radix("4c221e09f96b3229c95af490487612c8e3bd81704724eeda", 16).unwrap();

        let c: C192<P192, R192> = C192;

        let a = sec * c.G();

        assert_eq!(a.x.limbs, x);
        assert_eq!(a.y.limbs, y);
    }

    #[test]
    fn jacobian_point_multiplication_c192() {
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
        let x: BigUint = FromStrRadix::from_str_radix("48162eae1116dbbd5b7a0d9494ff0c9b414a31ce3d8b273f", 16).unwrap();
        let y: BigUint = FromStrRadix::from_str_radix("4c221e09f96b3229c95af490487612c8e3bd81704724eeda", 16).unwrap();

        let c: C192<P192, R192> = C192;

        let a = (sec * c.G().to_jacobian()).to_affine();

        assert_eq!(a.x.limbs, x);
        assert_eq!(a.y.limbs, y);
    }

    #[test]
    fn mixed_point_addition_c192() {
        let c: C192<P192, R192> = C192;

        let a = c.G().to_jacobian().double() + c.G();
        let b = c.G() + c.G().to_jacobian().double();
        let c = c.G().to_jacobian().double() + c.G().to_jacobian();

        assert_eq!(a, c);
        assert_eq!(b, c);
    }

    #[test]
    fn jacobian_point_multiplication_c256() {
        let sec: BigUint = FromStrRadix::from_str_radix("26254a72f6a0ce35958ce62ff0cea754b84ac449b2b340383faef50606d03b01", 16).unwrap();
        let x: BigUint = FromStrRadix::from_str_radix("6ee12372b80bad6f5432d0e6a3f02199db2b1617414f7fd8fe90b6bcf8b7aa68", 16).unwrap();
        let y: BigUint = FromStrRadix::from_str_radix("5c71967784765fe888995bfd9a1fb76f329630018430e1d9aca7b59dc672cad8", 16).unwrap();

        let c: C256<P256, R256> = C256;

        let a = (sec * c.G().to_jacobian()).to_affine();

        assert_eq!(a.x.limbs, x);
        assert_eq!(a.y.limbs, y);
    }

    #[bench]
    fn bench_point_mult_c192(b: &mut Bencher) {
        let c: C192<P192, R192> = C192;
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
        let p = c.G().to_jacobian();

        b.iter(|| { sec * p })
    }

    #[bench]
    fn bench_point_mult_c192_right_zeroes(b: &mut Bencher) {
        let c: C192<P192, R192> = C192;
        let sec: BigUint = FromStrRadix::from_str_radix("7e0000000000000000000000000000000000000000000000", 16).unwrap();
        let p = c.G().to_jacobian();

        b.iter(|| { sec * p })
    }

    #[bench]
    fn bench_point_mult_c192_left_zeroes(b: &mut Bencher) {
        let c: C192<P192, R192> = C192;
        let sec: BigUint = 3i.to_biguint().unwrap();
        let p = c.G().to_jacobian();

        b.iter(|| { sec * p })
    }

    #[bench]
    fn bench_point_mult_c521(b: &mut Bencher) {
        let c: C521<P521, R521> = C521;
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
        let p = c.G().to_jacobian();

        b.iter(|| { sec * p })
    }
}
