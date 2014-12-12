use std::fmt;
use std::num::FromStrRadix;
use num::{BigUint, Zero, One};
use num::bigint::ToBigUint;

use fields::{Field, FieldElem, P192, P521};

#[allow(non_snake_case)]
pub trait Curve<F: Field>: Clone { // To do:  Implement binary curves.
    fn A(&self) -> FieldElem<F>;
    fn B(&self) -> FieldElem<F>;
    fn G(&self) -> AffinePoint<Self, F>;
}

#[deriving(Clone)]
pub struct AffinePoint<C: Curve<F>, F: Field> {
    pub x: FieldElem<F>,
    pub y: FieldElem<F>,
    pub curve: C,
}

#[deriving(Clone)]
pub struct JacobianPoint<C: Curve<F>, F: Field> {
    pub x: FieldElem<F>,
    pub y: FieldElem<F>,
    pub z: FieldElem<F>,
    pub curve: C,
}


impl<C: Curve<F>, F: Field> PartialEq for AffinePoint<C, F> {
    fn eq(&self, other: &AffinePoint<C, F>) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: Curve<F>, F: Field> PartialEq for JacobianPoint<C, F> {
    fn eq(&self, other: &JacobianPoint<C, F>) -> bool {
        self.to_affine() == other.to_affine()
    }
}

// Notes:  Ordering is omitted because elliptic curve groups have no norm.

impl<C: Curve<F>, F: Field> Add<AffinePoint<C, F>, AffinePoint<C, F>> for AffinePoint<C, F> {
    fn add(&self, other: &AffinePoint<C, F>) -> AffinePoint<C, F> {
        if (*self).is_zero() {
            (*other).clone()
        } else if (*other).is_zero() {
            (*self).clone()
        } else if self.x != other.x {
            let m = (other.y - self.y) / (other.x - self.x);

            let x3 = (m * m) - self.x - other.x;
            let y3 = m * (self.x - x3) - self.y;

            AffinePoint {
                x: x3,
                y: y3,
                curve: self.curve.clone()
            }
        } else if self.y != other.y || self.y.limbs.bits() == 0 {
            self.zero()
        } else {
            let two = FieldElem {
                limbs: 2i.to_biguint().unwrap(),
                field: self.x.field.clone()
            };
            let three = FieldElem {
                limbs: 3i.to_biguint().unwrap(),
                field: self.x.field.clone()
            };

            let m = ((three * (self.x * self.x)) + self.curve.A()) / (two * self.y);

            let x3 = (m * m) - (self.x * two);
            let y3 = m * (self.x - x3) - self.y;

            AffinePoint {
                x: x3,
                y: y3,
                curve: self.curve.clone()
            }
        }
    }
}

impl<C: Curve<F>, F: Field> Add<JacobianPoint<C, F>, JacobianPoint<C, F>> for JacobianPoint<C, F> {
    fn add(&self, other: &JacobianPoint<C, F>) -> JacobianPoint<C, F> {
        if (*self).is_zero() {
            (*other).clone()
        } else if (*other).is_zero() {
            (*self).clone()
        } else if self.x == other.x && self.y == other.y && self.z == other.z {
            return self.double();
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
                    return self.zero();
                } else {
                    return self.double();
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

                return JacobianPoint {
                    x: x3,
                    y: y3,
                    z: z3,
                    curve: self.curve.clone()
                }
            }
        }
    }
}

impl<C: Curve<F>, F: Field> Sub<AffinePoint<C, F>, AffinePoint<C, F>> for AffinePoint<C, F> {
    fn sub(&self, other: &AffinePoint<C, F>) -> AffinePoint<C, F> { *self + (-*other) }
}

impl<C: Curve<F>, F: Field> Sub<JacobianPoint<C, F>, JacobianPoint<C, F>> for JacobianPoint<C, F> {
    fn sub(&self, other: &JacobianPoint<C, F>) -> JacobianPoint<C, F> { *self + (-*other) }
}

impl<C: Curve<F>, F:Field> Mul<BigUint, AffinePoint<C, F>> for AffinePoint<C, F> {
    fn mul(&self, other: &BigUint) -> AffinePoint<C, F> { self.pow(other) }
}

impl<C: Curve<F>, F:Field> Mul<AffinePoint<C, F>, AffinePoint<C, F>> for BigUint {
    fn mul(&self, other: &AffinePoint<C, F>) -> AffinePoint<C, F> { other.pow(self) }
}

impl<C: Curve<F>, F:Field> Mul<BigUint, JacobianPoint<C, F>> for JacobianPoint<C, F> {
    fn mul(&self, other: &BigUint) -> JacobianPoint<C, F> { self.pow(other) }
}

impl<C: Curve<F>, F:Field> Mul<JacobianPoint<C, F>, JacobianPoint<C, F>> for BigUint {
    fn mul(&self, other: &JacobianPoint<C, F>) -> JacobianPoint<C, F> { other.pow(self) }
}

impl<C: Curve<F>, F: Field> Neg<AffinePoint<C, F>> for AffinePoint<C, F> {
    fn neg(&self) -> AffinePoint<C, F> {
        AffinePoint {
            x: self.x.clone(),
            y: -(self.y.clone()),
            curve: self.curve.clone()
        }
    }
}

impl<C: Curve<F>, F: Field> Neg<JacobianPoint<C, F>> for JacobianPoint<C, F> {
    fn neg(&self) -> JacobianPoint<C, F> {
        JacobianPoint {
            x: self.x.clone(),
            y: -(self.y.clone()),
            z: self.z.clone(),
            curve: self.curve.clone()
        }
    }
}

impl<C: Curve<F>, F: Field> AffinePoint<C, F> {
    pub fn is_valid(&self) -> bool {
        if self.is_zero() {
            true
        } else {
            let y2 = self.y * self.y;
            let x = (self.x * self.x * self.x)
                + (self.curve.A() * self.x)
                + self.curve.B();

            y2 == x
        }
    }

    fn pow(&self, exp: &BigUint) -> AffinePoint<C, F> {
        // Montgomery Ladder.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();

        let m = exp.bits() + 1;
        let mut r0 = self.zero();
        let mut r1 = self.clone();

        for i in range(0u, m) {
            if ((one << (m - i - 1)) & *exp) == zer {
                r1 = r0 + r1;
                r0 = r0 + r0;
            } else {
                r0 = r0 + r1;
                r1 = r1 + r1;
            }
        }

        r0
    }

    pub fn to_jacobian(&self) -> JacobianPoint<C, F> {
        let p = JacobianPoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.x.one(),
            curve: self.curve.clone()
        };

        if self.is_zero() {
            p.zero() // Jacobian 0 != (0 : 0 : 1)
        } else {
            p
        }
    }

    // To do:  Find a way to actually implement Zero instead of mimic it.
    pub fn zero(&self) -> AffinePoint<C, F> {
        AffinePoint {
            x: self.x.zero(),
            y: self.y.zero(),
            curve: self.curve.clone()
        }
    }

    pub fn is_zero(&self) -> bool { self.x.is_zero() && self.y.is_zero() }
}

impl<C: Curve<F>, F: Field> JacobianPoint<C, F> {
    pub fn is_valid(&self) -> bool {
        if self.is_zero() {
            true
        } else {
            let z4 = self.z * self.z * self.z * self.z;

            let y2 = self.y * self.y;
            let x = (self.x * self.x * self.x)
                + (self.curve.A() * self.x * z4)
                + (self.curve.B() * z4 * self.z * self.x);

            y2 == x
        }
    }

    fn pow(&self, exp: &BigUint) -> JacobianPoint<C, F> { // Replace with generic
        // Montgomery Ladder.
        let zer: BigUint = Zero::zero();
        let one: BigUint = One::one();

        let m = exp.bits() + 1;
        let mut r0 = self.zero();
        let mut r1 = self.clone();

        for i in range(0u, m) {
            if ((one << (m - i - 1)) & *exp) == zer {
                r1 = r0 + r1;
                r0 = r0 + r0;
            } else {
                r0 = r0 + r1;
                r1 = r1 + r1;
            }
        }

        r0
    }

    pub fn double(&self) -> JacobianPoint<C, F> {
        let x2 = self.x * self.x;

        let y2 = self.y * self.y;
        let y4 = y2 * y2;

        let z2 = self.z * self.z;
        let z4  = z2 * z2;

        let xy2 = self.x * y2;
        let yz = self.y * self.z;

        let v = xy2 + xy2 + xy2 + xy2;
        let w = x2 + x2 + x2 + (self.curve.A() * z4);

        let v2 = self.x.field.reduce(v.limbs << 1);
        let y84 = self.x.field.reduce(y4.limbs << 3);

        let xf = (w * w) - v2;
        let yf = (w * (v - xf)) - y84;
        let zf = yz + yz;

        JacobianPoint {
            x: xf,
            y: yf,
            z: zf,
            curve: self.curve.clone()
        }
    }

    pub fn to_affine(&self) -> AffinePoint<C, F> {
        if self.is_zero() {
            let p = AffinePoint {
                x: self.x.zero(),
                y: self.y.zero(),
                curve: self.curve.clone()
            };

            p.zero() // In the event that affine zero isn't (0, 0)
        } else {
            let zinv = self.z.invert();

            AffinePoint {
                x: self.x * zinv * zinv,
                y: self.y * zinv * zinv * zinv,
                curve: self.curve.clone()
            }
        }
    }

    // To do:  Find a way to actually implement Zero instead of mimic it.
    pub fn zero(&self) -> JacobianPoint<C, F> {
        JacobianPoint {
            x: self.x.one(),
            y: self.y.one(),
            z: self.z.zero(),
            curve: self.curve.clone()
        }
    }

    pub fn is_zero(&self) -> bool { self.z.is_zero() }
}

impl<C: Curve<F>, F: Field> fmt::Show for AffinePoint<C, F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        } else {
            return write!(f, "({}, {})", self.x.limbs, self.y.limbs)
        }
    }
}

impl<C: Curve<F>, F: Field> fmt::Show for JacobianPoint<C, F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        } else {
            return write!(f, "({} : {} : {})", self.x.limbs, self.y.limbs, self.z.limbs)
        }
    }
}


#[deriving(Clone)]
pub struct C192<P192>;
impl Curve<P192> for C192<P192> {
    fn A(&self) -> FieldElem<P192> {
        -(FieldElem {
            limbs: 3i.to_biguint().unwrap(),
            field: P192
        })
    }

    fn B(&self) -> FieldElem<P192> {
        FieldElem {
            limbs: FromStrRadix::from_str_radix("64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1", 16).unwrap(),
            field: P192
        }
    }

    fn G(&self) -> AffinePoint<C192<P192>, P192> {
        AffinePoint {
            x: FieldElem { // To do:  Implement FromStrRadix for FieldElem
                limbs: FromStrRadix::from_str_radix("188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012", 16).unwrap(),
                field: P192
            },
            y: FieldElem {
                limbs: FromStrRadix::from_str_radix("07192b95ffc8da78631011ed6b24cdd573f977a11e794811", 16).unwrap(),
                field: P192
            },
            curve: C192
        }
    }
}

impl C192<P192> {
    pub fn base_point(&self) -> AffinePoint<C192<P192>, P192> { self.G() }
}

#[deriving(Clone)]
pub struct C521<P521>;
impl Curve<P521> for C521<P521> {
    fn A(&self) -> FieldElem<P521> {
        -(FieldElem {
            limbs: 3i.to_biguint().unwrap(),
            field: P521
        })
    }

    fn B(&self) -> FieldElem<P521> {
        FieldElem {
            limbs: FromStrRadix::from_str_radix("051953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00", 16).unwrap(),
            field: P521
        }
    }

    fn G(&self) -> AffinePoint<C521<P521>, P521> {
        AffinePoint {
            x: FieldElem { // To do:  Implement FromStrRadix for FieldElem
                limbs: FromStrRadix::from_str_radix("c6858e06b70404e9cd9e3ecb662395b4429c648139053fb521f828af606b4d3dbaa14b5e77efe75928fe1dc127a2ffa8de3348b3c1856a429bf97e7e31c2e5bd66", 16).unwrap(),
                field: P521
            },
            y: FieldElem {
                limbs: FromStrRadix::from_str_radix("11839296a789a3bc0045c8a5fb42c7d1bd998f54449579b446817afbd17273e662c97ee72995ef42640c550b9013fad0761353c7086a272c24088be94769fd16650", 16).unwrap(),
                field: P521
            },
            curve: C521
        }
    }
}

impl C521<P521> {
    pub fn base_point(&self) -> AffinePoint<C521<P521>, P521> { self.G() }
}

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

    use fields::{FieldElem, P192, R192, P521};
    use super::{AffinePoint, C192, C521};

    #[test]
    fn accept_valid_point() {
        let c: C192<P192> = C192;
        assert_eq!(c.base_point().is_valid(), true)
    }

    #[test]
    fn reject_invalid_point() {
        let c: C192<P192> = C192;
        let p = AffinePoint {
            x: c.base_point().x,
            y: FieldElem {
                limbs: c.base_point().y.limbs + One::one(),
                field: c.base_point().y.field
            },
            curve: c
        };

        assert_eq!(p.is_valid(), false)
    }

    #[test]
    fn base_point_field_is_r192() {
        let c: C192<P192> = C192;
        let one = FieldElem {
            limbs: One::one(),
            field: R192
        };

        let x = FieldElem {
            limbs: 3i.to_biguint().unwrap(),
            field: R192
        };

        let y = x.invert();

        let a = x.limbs * c.base_point();
        let b = y.limbs * a;

        assert!(x != y);
        assert!((x * y) == one);
        assert!(c.base_point() != a);
        assert!(c.base_point() == b);
    }

    #[test]
    fn affine_point_multiplication() {
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
        let x: BigUint = FromStrRadix::from_str_radix("48162eae1116dbbd5b7a0d9494ff0c9b414a31ce3d8b273f", 16).unwrap();
        let y: BigUint = FromStrRadix::from_str_radix("4c221e09f96b3229c95af490487612c8e3bd81704724eeda", 16).unwrap();

        let c: C192<P192> = C192;

        let a = sec * c.base_point();

        assert_eq!(a.x.limbs, x);
        assert_eq!(a.y.limbs, y);
    }

    #[test]
    fn jacobian_point_multiplication() {
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
        let x: BigUint = FromStrRadix::from_str_radix("48162eae1116dbbd5b7a0d9494ff0c9b414a31ce3d8b273f", 16).unwrap();
        let y: BigUint = FromStrRadix::from_str_radix("4c221e09f96b3229c95af490487612c8e3bd81704724eeda", 16).unwrap();

        let c: C192<P192> = C192;

        let a = (sec * c.base_point().to_jacobian()).to_affine();

        assert_eq!(a.x.limbs, x);
        assert_eq!(a.y.limbs, y);
    }

    #[test]
    fn curve521() {
        let sec: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();

        let c: C521<P521> = C521;

        let a = (sec * c.base_point().to_jacobian()).to_affine();

        println!("{}", a)
    }

    #[bench]
    fn bench_point_mult_c192(b: &mut Bencher) {
        let c: C192<P192> = C192;
        let sec: BigUint = FromStrRadix::from_str_radix("712627dcc76d4df2099af7c8e5", 16).unwrap();
        let p = c.base_point().to_jacobian();

        b.iter(|| { sec * p })
    }

    #[bench]
    fn bench_point_mult_c521(b: &mut Bencher) {
        let c: C521<P521> = C521;
        let sec: BigUint = FromStrRadix::from_str_radix("712627dcc76d4df2099af7c8e5", 16).unwrap();
        let p = c.base_point().to_jacobian();

        b.iter(|| { sec * p })
    }
}
