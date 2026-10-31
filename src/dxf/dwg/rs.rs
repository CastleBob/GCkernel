//! Reed-Solomon (255, k) decoding over GF(2^8) for the R2007 file organisation (ODA specification 5.13).
//!
//! Written from the public mathematics of Reed-Solomon codes: syndromes, the Berlekamp-Massey algorithm for the
//! error locator, a Chien search for its roots and Forney's formula for the error values.
//!
//! The parameters: the specification gives the field polynomials as coefficient lists of x^0..x^7 (x^8 implied):
//! system pages and the file header (1, 0, 0, 1, 0, 1, 1, 0) = x^8 + x^6 + x^5 + x^3 + 1 (0x169) with 16 parity
//! bytes; data pages (1, 0, 1, 1, 1, 0, 0, 0) = x^8 + x^4 + x^3 + x^2 + 1 (0x11D) with 4 parity bytes. It does not
//! state the generator roots or how the code word bytes map to polynomial coefficients; those were found by checking
//! the parity of a real AC1021 file (tests/fixtures/dwg/fixture_2007.dwg) for zero syndromes: the field element
//! alpha = x, the code word's 255 bytes as stored (k data bytes, then the parity bytes) are the coefficients of
//! x^0..x^254, and the code's roots are alpha^1 .. alpha^(n-k). (Equivalently, after a cyclic shift: parity bytes
//! at x^0.. and data at x^(n-k).., the usual systematic form with generator roots alpha^1 ..; a cyclic shift of a
//! code word of this code is again a code word, so both views have the same syndromes up to a factor.) For the
//! file header every first root alpha^0 .. alpha^254 was tried with 0x169 and 0x11D: only alpha^1 verifies (or
//! alpha^239 with the bytes read in reverse order, the same code backwards). The data page choice verifies every
//! data code word of the sample, interleaved (encoding 4) or not (encoding 1); see r2007.rs's test.

/// Code word length.
pub(super) const N: usize = 255;

pub(super) struct Code {
    exp: [u8; 2 * N],
    log: [u8; 256],
    /// Parity bytes per code word (n - k); corrects up to half as many byte errors.
    parity: usize,
    /// The exponent of the first root, alpha^first.
    first: usize,
}

/// At most 16 parity bytes.
const MAX_PARITY: usize = 16;

impl Code {
    /// `poly` has the x^8 term (bit 8) set and must be primitive.
    pub(super) fn new(poly: u16, parity: usize, first: usize) -> Code {
        let mut exp = [0u8; 2 * N];
        let mut log = [0u8; 256];
        let mut x: u16 = 1;
        for (i, e) in exp.iter_mut().enumerate().take(N) {
            *e = x as u8;
            log[x as usize] = i as u8;
            x <<= 1;
            if x & 0x100 != 0 {
                x ^= poly;
            }
        }
        for i in N..2 * N {
            exp[i] = exp[i - N];
        }
        Code { exp, log, parity: parity.min(MAX_PARITY), first }
    }

    /// The system page and file header code: RS(255, 239).
    pub(super) fn system() -> Code {
        Code::new(0x169, 16, 1)
    }

    /// The data page code: RS(255, 251).
    pub(super) fn data() -> Code {
        Code::new(0x11D, 4, 1)
    }

    /// Data bytes per code word.
    pub(super) fn k(&self) -> usize {
        N - self.parity
    }

    fn mul(&self, a: u8, b: u8) -> u8 {
        if a == 0 || b == 0 {
            0
        } else {
            self.exp[self.log[a as usize] as usize + self.log[b as usize] as usize]
        }
    }

    fn div(&self, a: u8, b: u8) -> u8 {
        // b != 0 is the caller's guarantee.
        if a == 0 || b == 0 {
            0
        } else {
            self.exp[(self.log[a as usize] as usize + N - self.log[b as usize] as usize) % N]
        }
    }

    /// alpha^e for any exponent.
    fn pow(&self, e: usize) -> u8 {
        self.exp[e % N]
    }

    /// p(alpha^e), p given lowest coefficient first.
    fn eval(&self, p: &[u8], e: usize) -> u8 {
        // Horner from the top coefficient.
        let x = self.pow(e);
        p.iter().rev().fold(0u8, |acc, &c| self.mul(acc, x) ^ c)
    }

    /// The syndromes S_j = c(alpha^(first + j)), j = 0 .. parity - 1.
    pub(super) fn syndromes(&self, word: &[u8; N]) -> [u8; MAX_PARITY] {
        let mut s = [0u8; MAX_PARITY];
        for (j, sj) in s.iter_mut().enumerate().take(self.parity) {
            *sj = self.eval(word, self.first + j);
        }
        s
    }

    /// Corrects up to parity / 2 byte errors in place. Returns the number of bytes corrected, or None when the
    /// word has more errors than the code can correct (as far as that can be told: such a word is refused unless
    /// it happens to lie within reach of another code word).
    pub(super) fn correct(&self, word: &mut [u8; N]) -> Option<usize> {
        let t2 = self.parity;
        let s = self.syndromes(word);
        if s[..t2].iter().all(|&x| x == 0) {
            return Some(0);
        }
        // Berlekamp-Massey: the shortest linear recurrence (error locator lambda) generating the syndromes.
        let mut lambda = [0u8; MAX_PARITY + 1];
        let mut prev = [0u8; MAX_PARITY + 1];
        lambda[0] = 1;
        prev[0] = 1;
        let mut l = 0usize;
        let mut m = 1usize;
        let mut b = 1u8;
        for n in 0..t2 {
            let mut d = s[n];
            for i in 1..=l.min(n) {
                d ^= self.mul(lambda[i], s[n - i]);
            }
            if d == 0 {
                m += 1;
                continue;
            }
            let coef = self.div(d, b);
            let before = lambda;
            for i in 0..=MAX_PARITY {
                if i >= m && prev[i - m] != 0 {
                    lambda[i] ^= self.mul(coef, prev[i - m]);
                }
            }
            if 2 * l <= n {
                l = n + 1 - l;
                prev = before;
                b = d;
                m = 1;
            } else {
                m += 1;
            }
        }
        if l == 0 || 2 * l > t2 {
            return None;
        }
        // Chien search: position i is in error when lambda(alpha^-i) = 0.
        let mut places = [0usize; MAX_PARITY];
        let mut found = 0usize;
        for i in 0..N {
            if self.eval(&lambda[..=l], N - i) == 0 {
                if found == l {
                    return None;
                }
                places[found] = i;
                found += 1;
            }
        }
        if found != l {
            return None;
        }
        // Forney: e = X^(1 - first) * omega(X^-1) / lambda'(X^-1), omega = S * lambda mod x^(parity).
        let mut omega = [0u8; MAX_PARITY];
        for i in 0..t2 {
            for j in 0..=i.min(l) {
                omega[i] ^= self.mul(s[i - j], lambda[j]);
            }
        }
        let mut derivative = [0u8; MAX_PARITY];
        for k in (1..=l).step_by(2) {
            derivative[k - 1] = lambda[k];
        }
        for &i in &places[..found] {
            let inverse = N - i; // X^-1 = alpha^(-i)
            let denominator = self.eval(&derivative[..l], inverse);
            if denominator == 0 {
                return None;
            }
            let numerator = self.eval(&omega[..t2], inverse);
            // X^(1 - first) = alpha^(i * (1 - first)), taken mod 255.
            let shift = (i * ((1 + N - self.first % N) % N)) % N;
            let value = self.mul(self.pow(shift), self.div(numerator, denominator));
            word[i] ^= value;
        }
        // A corrected word is a code word; anything else was beyond repair.
        if self.syndromes(word)[..t2].iter().any(|&x| x != 0) {
            return None;
        }
        Some(found)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A code word of `code`: data at the low positions, then parity chosen so that every syndrome is zero (solved
    /// as a linear system over GF(2^8): the parity positions' contributions to each syndrome).
    fn encode(code: &Code, data: &[u8]) -> [u8; N] {
        let k = code.k();
        let p = code.parity;
        let mut word = [0u8; N];
        word[..k].copy_from_slice(&data[..k]);
        let s = code.syndromes(&word);
        // Matrix rows j: sum_i alpha^((first + j)(k + i)) * parity_i = s_j.
        let mut a = vec![vec![0u8; p + 1]; p];
        for j in 0..p {
            for i in 0..p {
                a[j][i] = code.pow((code.first + j) * (k + i));
            }
            a[j][p] = s[j];
        }
        for col in 0..p {
            let pivot = (col..p).find(|&r| a[r][col] != 0).unwrap();
            a.swap(col, pivot);
            let inv = code.div(1, a[col][col]);
            for c in 0..=p {
                a[col][c] = code.mul(a[col][c], inv);
            }
            for r in 0..p {
                if r != col && a[r][col] != 0 {
                    let f = a[r][col];
                    for c in 0..=p {
                        let v = code.mul(f, a[col][c]);
                        a[r][c] ^= v;
                    }
                }
            }
        }
        for i in 0..p {
            word[k + i] = a[i][p];
        }
        assert!(code.syndromes(&word)[..p].iter().all(|&x| x == 0));
        word
    }

    #[test]
    fn fields_are_primitive() {
        for code in [Code::system(), Code::data()] {
            let mut seen = [false; 256];
            for i in 0..N {
                seen[code.exp[i] as usize] = true;
            }
            assert_eq!(seen.iter().filter(|&&x| x).count(), 255);
        }
    }

    #[test]
    fn corrects_up_to_half_the_parity() {
        let mut seed = 0x2545F4914F6CDD1Du64;
        let mut next = move || { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed };
        for code in [Code::system(), Code::data()] {
            let data: Vec<u8> = (0..N).map(|_| next() as u8).collect();
            let good = encode(&code, &data);
            for errors in 1..=code.parity / 2 {
                for _ in 0..200 {
                    let mut w = good;
                    let mut places = Vec::new();
                    while places.len() < errors {
                        let at = (next() as usize) % N;
                        if !places.contains(&at) {
                            places.push(at);
                        }
                    }
                    for &at in &places {
                        w[at] ^= 1 + (next() % 255) as u8;
                    }
                    assert_eq!(code.correct(&mut w), Some(errors));
                    assert_eq!(w, good);
                }
            }
            // One more error than the code corrects: refused (or, rarely, taken to another code word, never to
            // the original one).
            let mut refused = 0;
            for _ in 0..200 {
                let mut w = good;
                let mut places = Vec::new();
                while places.len() < code.parity / 2 + 1 {
                    let at = (next() as usize) % N;
                    if !places.contains(&at) {
                        places.push(at);
                    }
                }
                for &at in &places {
                    w[at] ^= 1 + (next() % 255) as u8;
                }
                match code.correct(&mut w) {
                    None => refused += 1,
                    Some(_) => assert_ne!(w, good),
                }
            }
            if code.parity == 16 {
                assert_eq!(refused, 200);
            } else {
                assert!(refused > 50, "{refused}");
            }
        }
    }
}
