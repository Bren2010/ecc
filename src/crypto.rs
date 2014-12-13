use std::rand;

use num::BigUint;
use num::bigint::RandBigInt;

use fields::Field;
use curves::Curve;

pub struct DiffieHellman;

impl DiffieHellman {
    pub fn key_gen<C: Curve<F>, F: Field>(c: &C) -> (Vec<uint>, BigUint) {
        let mut rng = match rand::OsRng::new() {
            Ok(g) => g,
            Err(e) => panic!("Could not load the OS' RNG! Error: {}", e)
        };

        let x = rng.gen_biguint_below(&c.G().x.field.modulus());
        let X = x * c.G();

        (X.serialize(), x)
    }

    pub fn shared<C: Curve<F>, F: Field>(c: &C, x: &BigUint, Y: &Vec<uint>) -> Option<Vec<uint>> {
        let Y2 = c.unserialize(Y);

        if Y2.is_valid() && !Y2.is_zero() {
            Some((*x * Y2).x.serialize())
        } else {
            None
        }
    }
}


// -------------------------------------------------------------------------
// Unit Tests
// -------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use fields::P192;
    use curves::C192;
    use super::DiffieHellman;

    #[test]
    fn valid_dh_exchange() {
        let c: C192<P192> = C192;

        let (X, x) = DiffieHellman::key_gen(&c);
        let (Y, y) = DiffieHellman::key_gen(&c);

        let s1 = DiffieHellman::shared(&c, &x, &Y).unwrap();
        let s2 = DiffieHellman::shared(&c, &y, &X).unwrap();

        assert_eq!(s1, s2)
    }

    #[test]
    fn invalid_dh_exchange() {
        let c: C192<P192> = C192;

        let (_, x) = DiffieHellman::key_gen(&c);
        let Y = vec![0, 0, 1, 0];

        let s1 = DiffieHellman::shared(&c, &x, &Y);

        assert_eq!(s1, None)
    }
}
