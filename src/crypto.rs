use std::rand;
use std::default::Default;

use num::{BigUint, Zero};
use num::bigint::RandBigInt;

use fields::Field;
use curves::{Curve, AffinePoint};

pub struct DiffieHellman;

impl DiffieHellman {
    pub fn key_gen<C: Curve<F, G>, F: Field, G: Field>() -> (AffinePoint<C, F, G>, BigUint) {
        let mut rng = match rand::OsRng::new() {
            Ok(g) => g,
            Err(e) => panic!("Could not load the OS' RNG! Error: {}", e)
        };

        let c: C = Default::default();
        let f: F = Default::default();

        let x = rng.gen_biguint_below(&f.modulus());
        let X = x * c.G();

        (X, x)
    }

    pub fn shared<C: Curve<F, G>, F: Field, G: Field>(x: &BigUint, Y: &AffinePoint<C, F, G>) -> Option<AffinePoint<C, F, G>> {
        let r: G = Default::default();
        let n: BigUint = r.modulus();

        if Y.is_valid() && !Y.is_zero() && (n * *Y).is_zero() {
            Some((*x * *Y))
        } else {
            None
        }
    }
}


pub struct ElGamalEncryptor<C: Curve<F, G>, F: Field, G: Field> {
    pub pub_key: AffinePoint<C, F, G>,
}

impl<C: Curve<F, G>, F: Field, G: Field> ElGamalEncryptor<C, F, G> {
    pub fn encrypt(&self, msg: &AffinePoint<C, F, G>) -> (AffinePoint<C, F, G>, AffinePoint<C, F, G>) {
        let (K, k) = DiffieHellman::key_gen();
        let ct = *msg + (k * self.pub_key);

        (K, ct)
    }
}


pub struct ElGamalDecryptor<C: Curve<F, G>, F: Field, G: Field> {
    pub pub_key: AffinePoint<C, F, G>,
    pub priv_key: BigUint
}

impl<C: Curve<F, G>, F: Field, G: Field> ElGamalDecryptor<C, F, G> {
    pub fn new(sk: &BigUint) -> ElGamalDecryptor<C, F, G> {
        let c: C = Default::default();

        ElGamalDecryptor {
            pub_key: (*sk * c.G().to_jacobian()).to_affine(),
            priv_key: (*sk).clone()
        }
    }

    pub fn new_key_gen() -> ElGamalDecryptor<C, F, G> {
        let (A, a): (AffinePoint<C, F, G>, BigUint) = DiffieHellman::key_gen();

        ElGamalDecryptor {
            pub_key: A,
            priv_key: a
        }
    }

    pub fn decrypt(&self, nonce: &AffinePoint<C, F, G>, ct: &AffinePoint<C, F, G>) -> Option<AffinePoint<C, F, G>> {
        let s = DiffieHellman::shared(&self.priv_key, nonce);

        if s == None || !ct.is_valid() || ct.is_zero() {
            None
        } else {
            Some(*ct - (self.priv_key * *nonce))
        }
    }
}

// -------------------------------------------------------------------------
// Unit Tests
// -------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use num::Zero;
    use num::bigint::ToBigUint;

    use fields::{FieldElem, P192, R192};
    use curves::{Curve, AffinePoint, C192};
    use super::{DiffieHellman, ElGamalEncryptor, ElGamalDecryptor};

    #[test]
    fn valid_dh_exchange() {
        type Point = AffinePoint<C192<P192, R192>, P192, R192>;

        let (X, x): (Point, _) = DiffieHellman::key_gen();
        let (Y, y): (Point, _) = DiffieHellman::key_gen();

        let s1 = DiffieHellman::shared(&x, &Y).unwrap();
        let s2 = DiffieHellman::shared(&y, &X).unwrap();

        assert_eq!(s1, s2)
    }

    #[test]
    fn invalid_dh_exchange() {
        type Point = AffinePoint<C192<P192, R192>, P192, R192>;

        let (_, x): (Point, _) = DiffieHellman::key_gen();
        let Y: Point = AffinePoint {
            x: FieldElem { limbs: Zero::zero() },
            y: FieldElem { limbs: Zero::zero() }
        };

        let s1 = DiffieHellman::shared(&x, &Y);

        assert_eq!(s1, None)
    }

    #[test]
    fn elgamal_encryption() {
        type Enc = ElGamalEncryptor<C192<P192, R192>, P192, R192>;
        type Dec = ElGamalDecryptor<C192<P192, R192>, P192, R192>;

        let c: C192<P192, R192> = C192;
        let M = 7i.to_biguint().unwrap() * c.G();

        let dec: Dec = ElGamalDecryptor::new_key_gen();
        let enc: Enc = ElGamalEncryptor{ pub_key: dec.pub_key.clone() };

        let (nonce1, ct1) = enc.encrypt(&M);
        let (nonce2, ct2) = enc.encrypt(&M);

        assert!(nonce1 != nonce2);
        assert!(ct1 != ct2);

        let dec1 = dec.decrypt(&nonce1, &ct1);
        let dec2 = dec.decrypt(&nonce2, &ct2);

        assert_eq!(dec1.unwrap(), M);
        assert_eq!(dec2.unwrap(), M);
    }
}
