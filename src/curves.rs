use std::fmt;
use std::num::FromStrRadix;
use num::{BigUint, Zero, One};
use num::bigint::ToBigUint;

use fields::{Field, FieldElem, P192};

#[allow(non_snake_case)]
pub trait Curve<F: Field>: Clone { // To do:  Implement binary curves.
    fn A(&self) -> FieldElem<F>;
    fn B(&self) -> FieldElem<F>;
    fn G(&self) -> Point<Self, F>;
}

#[deriving(Clone)]
pub struct Point<C: Curve<F>, F: Field> {
    pub x: FieldElem<F>,
    pub y: FieldElem<F>,
    pub curve: C,
}


impl<C: Curve<F>, F: Field> PartialEq for Point<C, F> {
    fn eq(&self, other: &Point<C, F>) -> bool {
        self.x == other.x && self.y == other.y
    }
}

// Notes:  Ordering is omitted because elliptic curve groups have no norm.

impl<C: Curve<F>, F: Field> Add<Point<C, F>, Point<C, F>> for Point<C, F> {
    fn add(&self, other: &Point<C, F>) -> Point<C, F> {
        if (*self).is_zero() {
            (*other).clone()
        } else if (*other).is_zero() {
            (*self).clone()
        } else if self.x != other.x {
            let m = (other.y - self.y) / (other.x - self.x);

            let x3 = (m * m) - self.x - other.x;
            let y3 = m * (self.x - x3) - self.y;

            Point {
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

            Point {
                x: x3,
                y: y3,
                curve: self.curve.clone()
            }
        }
    }
}

impl<C: Curve<F>, F: Field> Sub<Point<C, F>, Point<C, F>> for Point<C, F> {
    fn sub(&self, other: &Point<C, F>) -> Point<C, F> { *self + (-*other) }
}

impl<C: Curve<F>, F:Field> Mul<BigUint, Point<C, F>> for Point<C, F> {
    fn mul(&self, other: &BigUint) -> Point<C, F> { self.pow(other) }
}

impl<C: Curve<F>, F:Field> Mul<Point<C, F>, Point<C, F>> for BigUint {
    fn mul(&self, other: &Point<C, F>) -> Point<C, F> { other.pow(self) }
}

impl<C: Curve<F>, F: Field> Neg<Point<C, F>> for Point<C, F> {
    fn neg(&self) -> Point<C, F> {
        Point {
            x: self.x.clone(),
            y: -(self.y.clone()),
            curve: self.curve.clone()
        }
    }
}

impl<C: Curve<F>, F: Field> Point<C, F> {
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

    fn pow(&self, exp: &BigUint) -> Point<C, F> {
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

    // To do:  Find a way to actually implement Zero instead of mimic it.
    pub fn zero(&self) -> Point<C, F> {
        Point {
            x: self.x.zero(),
            y: self.y.zero(),
            curve: self.curve.clone()
        }
    }

    pub fn is_zero(&self) -> bool { self.x.is_zero() && self.y.is_zero() }
}

impl<C: Curve<F>, F: Field> fmt::Show for Point<C, F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        } else {
            return write!(f, "({}, {})", self.x.limbs, self.y.limbs)
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

    fn G(&self) -> Point<C192<P192>, P192> {
        Point {
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
    pub fn base_point(&self) -> Point<C192<P192>, P192> { self.G() }
}


// -------------------------------------------------------------------------
// Unit Tests
// -------------------------------------------------------------------------
#[cfg(test)]
mod test {
    use std::num::FromStrRadix;
    use num::{BigUint, One};
    use num::bigint::ToBigUint;

    use fields::{FieldElem, P192, R192};
    use super::{Point, C192};

    #[test]
    fn accept_valid_point() {
        let c: C192<P192> = C192;
        assert_eq!(c.base_point().is_valid(), true)
    }

    #[test]
    fn reject_invalid_point() {
        let c: C192<P192> = C192;
        let p = Point {
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
}
