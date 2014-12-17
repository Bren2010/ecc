Elliptic Curve Crypto
=====================

An elliptic curve arithmetic and cryptography library in Pure Rust&reg;.

Install by adding this to your `Cargo.toml`:
```
[dependencies.ecc]
version = "*"
```

The obligatory calling-my-code-terrible-and-warning-you-to-use-it-at-your-own-peril paragraph:  I've attempted to prevent any timing or invalid curve attacks, but this is the first thing I've ever written in Rust and I'm also probably the only person who's ever read the code.  There are tests that say the code is correct and benchmarks that say it's not horrendously slow, but I don't know if I believe either of them myself.  This is all a work-in-progress.

### To Do:
1.  ElGamal encryption
2.  ECDSA / ElGamal signatures


Basics
------

Currently the only high-level operation implemented is Diffie-Hellman key exchanges:
```rust
extern crate ecc;

use ecc::fields::{P256, R256}; // NIST's fields for P-256
use ecc::curves::C256; // NIST's curve P-256
use ecc::crypto::DiffieHellman;

fn main() {
  type Curve = C192<P192, R192>;
  type Point = AffinePoint<Curve, P192, R192>;

  let (X, x): (Point, _) = DiffieHellman::key_gen();
  let out = X.serialize();
  // out: Vec<uint>  -> Send to the other person.
  // x: BigUint      -> Kept secret.

  // X             -----> Other person
  // Y: Vec<uint>  <-----

  let c: Curve = C192;
  let in = c.unserialize(Y);

  let s: Option<Point> = DiffieHellman::shared(&x, &in);
  // Will return None if trickery occured.
  // Will return the shared secret Some(...), which should be serialized and
  // put through a KDF or something and then used in a cipher/MAC.
}
```

Fields
------

There are six fields implemented in total:  `P192`, `P256`, and `P521` are the base fields of their corresponding NIST curves.  `R192`, `R256`, and `R521`, however, are the fields that act on the base point of their corresponding curves.

Creating a field is easy, but not much use:
```rust
let f: P192 = P192;
```

Creating a field element is more interesting:
```rust
extern crate num;

use num::bigint::ToBigUint;
use ecc::fields::{FieldElem, P192};

...

let x: FieldElem<P192> = FieldElem { limbs: 3i.to_biguint().unwrap() }
// Or any other BigUint in the `limbs` field.
```

Field operations can be applied to field elements with the normal unary/binary operators (negation, addition, subtraction, multiplication, and division).  Calling `x.invert()` will return `x`'s inverse, and calling `x.pow(exp)`, where `exp` is a `BigUint`, will raise `x` to the power of `exp`.

Field elements don't support an order, but they can be tested for equality with the `==` operator.

A field element can be serialized to a byte vector by calling `x.serialize()`, or unserialized by calling `field.unserialize(...)`.

Exponentiation and testing for equality are constant-time, by default.

Curves
------

There are currently only three NIST curves implemented:  `C192`, `C256`, and `C521`.

Creating a curve is similar to creating a field:
```rust
let c: C192<P192, R192> = C192;
```

Calling `c.G()` will return the curve's base point in affine coordinates.  Calling `.to_jacobian()` will convert it to Jacobian coordinates, and alternately, calling `.to_affine()` on a Jacobian point will convert it back to affine.

Negation, addition, and subtraction are supported on points of the same coordinate system.  Points can also be multiplied by a `BigUint`.  To check if a point is valid, call `A.is_valid()`.  To check if a point is zero (the point at inifinity), call `A.is_zero()`.

A point can be serialized to a byte vector by calling `A.serialize()`, or unserialized by calling `curve.unserialize(...)`.

Point multiplication and testing for equality is constant-time, by default.

```rust
extern crate num;

use std::num::FromStrRadix;
use fields::{P192, R192};
use curves::{Curve, C192};

fn main() {
  let curve C192<P192, R192> = C192;
  let a: BigUint = FromStrRadix::from_str_radix("7e48c5ab7f43e4d9c17bd9712627dcc76d4df2099af7c8e5", 16).unwrap();
  let G = curve.G().to_jacobian();

  let A = a * G;
}
```
