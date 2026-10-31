//! The standard security handler (ISO 32000-1 7.6: Algorithm 1 object keys, Algorithm 2 file key, Algorithms 4-6
//! user password check; ISO 32000-2 Algorithms 2.A and 2.B for AES-256, revisions 5 and 6), with our own MD5
//! (RFC 1321), SHA-256/384/512 (FIPS 180-4), RC4 and AES (FIPS 197). Only documents that open without a password
//! (an empty user password, or an empty owner password) are decrypted.

use super::lexer::{Dict, Obj};
use crate::KernelError;

// ---- MD5 (RFC 1321) ------------------------------------------------------------------------------------------------

const MD5_T: [u32; 64] = [
    0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501, 0x698098d8,
    0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821, 0xf61e2562, 0xc040b340,
    0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8, 0x21e1cde6, 0xc33707d6, 0xf4d50d87,
    0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a, 0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
    0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, 0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039,
    0xe6db99e5, 0x1fa27cf8, 0xc4ac5665, 0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92,
    0xffeff47d, 0x85845dd1, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb,
    0xeb86d391,
];

/// Per-step shift amounts of the four rounds (RFC 1321 3.4).
const MD5_S: [[u32; 4]; 4] = [[7, 12, 17, 22], [5, 9, 14, 20], [4, 11, 16, 23], [6, 10, 15, 21]];

pub fn md5(data: &[u8]) -> [u8; 16] {
    let mut h: [u32; 4] = [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476];
    let mut msg = data.to_vec();
    let bits = (data.len() as u64).wrapping_mul(8);
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bits.to_le_bytes());
    for block in msg.chunks_exact(64) {
        let mut x = [0u32; 16];
        for (j, w) in x.iter_mut().enumerate() {
            *w = u32::from_le_bytes([block[4 * j], block[4 * j + 1], block[4 * j + 2], block[4 * j + 3]]);
        }
        let (mut a, mut b, mut c, mut d) = (h[0], h[1], h[2], h[3]);
        for i in 0..64 {
            let round = i / 16;
            let (f, k) = match round {
                0 => ((b & c) | (!b & d), i),
                1 => ((b & d) | (c & !d), (1 + 5 * i) % 16),
                2 => (b ^ c ^ d, (5 + 3 * i) % 16),
                _ => (c ^ (b | !d), (7 * i) % 16),
            };
            let t = a.wrapping_add(f).wrapping_add(x[k]).wrapping_add(MD5_T[i]).rotate_left(MD5_S[round][i % 4]);
            let new_b = b.wrapping_add(t);
            a = d;
            d = c;
            c = b;
            b = new_b;
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
    }
    let mut out = [0u8; 16];
    for (i, w) in h.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&w.to_le_bytes());
    }
    out
}

// ---- SHA-2 (FIPS 180-4) --------------------------------------------------------------------------------------------

const SHA256_K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98,
    0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
    0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8,
    0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819,
    0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
    0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
    0xc67178f2,
];

pub fn sha256(data: &[u8]) -> [u8; 32] {
    let mut h: [u32; 8] = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    let mut msg = data.to_vec();
    let bits = (data.len() as u64).wrapping_mul(8);
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bits.to_be_bytes());
    for block in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for t in 0..16 {
            w[t] = u32::from_be_bytes([block[4 * t], block[4 * t + 1], block[4 * t + 2], block[4 * t + 3]]);
        }
        for t in 16..64 {
            let s0 = w[t - 15].rotate_right(7) ^ w[t - 15].rotate_right(18) ^ (w[t - 15] >> 3);
            let s1 = w[t - 2].rotate_right(17) ^ w[t - 2].rotate_right(19) ^ (w[t - 2] >> 10);
            w[t] = s1.wrapping_add(w[t - 7]).wrapping_add(s0).wrapping_add(w[t - 16]);
        }
        let mut v = h;
        for t in 0..64 {
            let (a, b, c, d, e, f, g, hh) = (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
            let big1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ (!e & g);
            let t1 = hh.wrapping_add(big1).wrapping_add(ch).wrapping_add(SHA256_K[t]).wrapping_add(w[t]);
            let big0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = big0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), a, b, c, d.wrapping_add(t1), e, f, g];
        }
        for i in 0..8 {
            h[i] = h[i].wrapping_add(v[i]);
        }
    }
    let mut out = [0u8; 32];
    for (i, w) in h.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&w.to_be_bytes());
    }
    out
}

const SHA512_K: [u64; 80] = [
    0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc, 0x3956c25bf348b538,
    0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118, 0xd807aa98a3030242, 0x12835b0145706fbe,
    0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2, 0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235,
    0xc19bf174cf692694, 0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
    0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5, 0x983e5152ee66dfab,
    0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4, 0xc6e00bf33da88fc2, 0xd5a79147930aa725,
    0x06ca6351e003826f, 0x142929670a0e6e70, 0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed,
    0x53380d139d95b3df, 0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
    0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30, 0xd192e819d6ef5218,
    0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8, 0x19a4c116b8d2d0c8, 0x1e376c085141ab53,
    0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8, 0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373,
    0x682e6ff3d6b2b8a3, 0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
    0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b, 0xca273eceea26619c,
    0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178, 0x06f067aa72176fba, 0x0a637dc5a2c898a6,
    0x113f9804bef90dae, 0x1b710b35131c471b, 0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc,
    0x431d67c49c100d4c, 0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817,
];

fn sha512_core(data: &[u8], init: [u64; 8], out_len: usize) -> Vec<u8> {
    let mut h = init;
    let mut msg = data.to_vec();
    let bits = (data.len() as u128).wrapping_mul(8);
    msg.push(0x80);
    while msg.len() % 128 != 112 {
        msg.push(0);
    }
    msg.extend_from_slice(&bits.to_be_bytes());
    for block in msg.chunks_exact(128) {
        let mut w = [0u64; 80];
        for t in 0..16 {
            let mut b8 = [0u8; 8];
            b8.copy_from_slice(&block[8 * t..8 * t + 8]);
            w[t] = u64::from_be_bytes(b8);
        }
        for t in 16..80 {
            let s0 = w[t - 15].rotate_right(1) ^ w[t - 15].rotate_right(8) ^ (w[t - 15] >> 7);
            let s1 = w[t - 2].rotate_right(19) ^ w[t - 2].rotate_right(61) ^ (w[t - 2] >> 6);
            w[t] = s1.wrapping_add(w[t - 7]).wrapping_add(s0).wrapping_add(w[t - 16]);
        }
        let mut v = h;
        for t in 0..80 {
            let (a, b, c, d, e, f, g, hh) = (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
            let big1 = e.rotate_right(14) ^ e.rotate_right(18) ^ e.rotate_right(41);
            let ch = (e & f) ^ (!e & g);
            let t1 = hh.wrapping_add(big1).wrapping_add(ch).wrapping_add(SHA512_K[t]).wrapping_add(w[t]);
            let big0 = a.rotate_right(28) ^ a.rotate_right(34) ^ a.rotate_right(39);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = big0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), a, b, c, d.wrapping_add(t1), e, f, g];
        }
        for i in 0..8 {
            h[i] = h[i].wrapping_add(v[i]);
        }
    }
    let mut out: Vec<u8> = h.iter().flat_map(|w| w.to_be_bytes()).collect();
    out.truncate(out_len);
    out
}

pub fn sha384(data: &[u8]) -> Vec<u8> {
    sha512_core(data, [
        0xcbbb9d5dc1059ed8, 0x629a292a367cd507, 0x9159015a3070dd17, 0x152fecd8f70e5939,
        0x67332667ffc00b31, 0x8eb44a8768581511, 0xdb0c2e0d64f98fa7, 0x47b5481dbefa4fa4,
    ], 48)
}

pub fn sha512(data: &[u8]) -> Vec<u8> {
    sha512_core(data, [
        0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
        0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179,
    ], 64)
}

// ---- RC4 -----------------------------------------------------------------------------------------------------------

/// RC4 is its own inverse: the key schedule permutes 0..255 by the key, then each output byte XORs the keystream.
pub fn rc4(key: &[u8], data: &[u8]) -> Vec<u8> {
    if key.is_empty() {
        return data.to_vec();
    }
    let mut s: [u8; 256] = std::array::from_fn(|i| i as u8);
    let mut j: u8 = 0;
    for i in 0..256 {
        j = j.wrapping_add(s[i]).wrapping_add(key[i % key.len()]);
        s.swap(i, j as usize);
    }
    let (mut i, mut j) = (0u8, 0u8);
    data.iter()
        .map(|&b| {
            i = i.wrapping_add(1);
            j = j.wrapping_add(s[i as usize]);
            s.swap(i as usize, j as usize);
            b ^ s[s[i as usize].wrapping_add(s[j as usize]) as usize]
        })
        .collect()
}

// ---- AES (FIPS 197) ------------------------------------------------------------------------------------------------

const SBOX: [u8; 256] = [
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
];

/// Multiplication in GF(2^8) modulo x^8 + x^4 + x^3 + x + 1 (FIPS 197 4.2).
fn gmul(mut a: u8, mut b: u8) -> u8 {
    let mut p = 0u8;
    while b != 0 {
        if b & 1 != 0 {
            p ^= a;
        }
        let high = a & 0x80;
        a <<= 1;
        if high != 0 {
            a ^= 0x1b;
        }
        b >>= 1;
    }
    p
}

pub struct Aes {
    /// Round keys, 16 bytes per round, Nr + 1 rounds.
    w: Vec<[u8; 16]>,
    inv: [u8; 256],
}

impl Aes {
    /// KEYEXPANSION for 16, 24 or 32 byte keys.
    pub fn new(key: &[u8]) -> Option<Aes> {
        let nk = key.len() / 4;
        if !matches!(key.len(), 16 | 24 | 32) {
            return None;
        }
        let nr = nk + 6;
        let total = 4 * (nr + 1);
        let mut words: Vec<[u8; 4]> = Vec::with_capacity(total);
        for i in 0..nk {
            words.push([key[4 * i], key[4 * i + 1], key[4 * i + 2], key[4 * i + 3]]);
        }
        let mut rcon = 1u8;
        for i in nk..total {
            let mut t = words[i - 1];
            if i % nk == 0 {
                t = [SBOX[t[1] as usize] ^ rcon, SBOX[t[2] as usize], SBOX[t[3] as usize], SBOX[t[0] as usize]];
                rcon = gmul(rcon, 2);
            } else if nk > 6 && i % nk == 4 {
                t = t.map(|b| SBOX[b as usize]);
            }
            let p = words[i - nk];
            words.push([p[0] ^ t[0], p[1] ^ t[1], p[2] ^ t[2], p[3] ^ t[3]]);
        }
        let w = words.chunks_exact(4).map(|c| {
            let mut k = [0u8; 16];
            for (j, word) in c.iter().enumerate() {
                k[4 * j..4 * j + 4].copy_from_slice(word);
            }
            k
        }).collect();
        let mut inv = [0u8; 256];
        for (i, s) in SBOX.iter().enumerate() {
            inv[*s as usize] = i as u8;
        }
        Some(Aes { w, inv })
    }

    /// The last round key's last word (for checking the key expansion against FIPS 197 Appendix A).
    #[cfg(test)]
    pub fn last_word(&self) -> [u8; 4] {
        let k = self.w[self.w.len() - 1];
        [k[12], k[13], k[14], k[15]]
    }

    /// CIPHER(): the state is the block column by column (byte r + 4c is row r, column c).
    pub fn encrypt_block(&self, block: &mut [u8; 16]) {
        let nr = self.w.len() - 1;
        xor16(block, &self.w[0]);
        for round in 1..=nr {
            for b in block.iter_mut() {
                *b = SBOX[*b as usize];
            }
            shift_rows(block);
            if round != nr {
                for c in 0..4 {
                    let s = [block[4 * c], block[4 * c + 1], block[4 * c + 2], block[4 * c + 3]];
                    block[4 * c] = gmul(s[0], 2) ^ gmul(s[1], 3) ^ s[2] ^ s[3];
                    block[4 * c + 1] = s[0] ^ gmul(s[1], 2) ^ gmul(s[2], 3) ^ s[3];
                    block[4 * c + 2] = s[0] ^ s[1] ^ gmul(s[2], 2) ^ gmul(s[3], 3);
                    block[4 * c + 3] = gmul(s[0], 3) ^ s[1] ^ s[2] ^ gmul(s[3], 2);
                }
            }
            xor16(block, &self.w[round]);
        }
    }

    /// INVCIPHER(): the rounds undone in reverse order.
    pub fn decrypt_block(&self, block: &mut [u8; 16]) {
        let nr = self.w.len() - 1;
        xor16(block, &self.w[nr]);
        for round in (0..nr).rev() {
            inv_shift_rows(block);
            for b in block.iter_mut() {
                *b = self.inv[*b as usize];
            }
            xor16(block, &self.w[round]);
            if round != 0 {
                for c in 0..4 {
                    let s = [block[4 * c], block[4 * c + 1], block[4 * c + 2], block[4 * c + 3]];
                    block[4 * c] = gmul(s[0], 14) ^ gmul(s[1], 11) ^ gmul(s[2], 13) ^ gmul(s[3], 9);
                    block[4 * c + 1] = gmul(s[0], 9) ^ gmul(s[1], 14) ^ gmul(s[2], 11) ^ gmul(s[3], 13);
                    block[4 * c + 2] = gmul(s[0], 13) ^ gmul(s[1], 9) ^ gmul(s[2], 14) ^ gmul(s[3], 11);
                    block[4 * c + 3] = gmul(s[0], 11) ^ gmul(s[1], 13) ^ gmul(s[2], 9) ^ gmul(s[3], 14);
                }
            }
        }
    }

    /// CBC encryption without padding (the input is a whole number of blocks).
    pub fn cbc_encrypt(&self, iv: &[u8; 16], data: &[u8]) -> Vec<u8> {
        let mut prev = *iv;
        let mut out = Vec::with_capacity(data.len());
        for chunk in data.chunks_exact(16) {
            let mut b = [0u8; 16];
            for i in 0..16 {
                b[i] = chunk[i] ^ prev[i];
            }
            self.encrypt_block(&mut b);
            out.extend_from_slice(&b);
            prev = b;
        }
        out
    }

    /// CBC decryption of whole blocks, no padding removed.
    pub fn cbc_decrypt(&self, iv: &[u8; 16], data: &[u8]) -> Vec<u8> {
        let mut prev = *iv;
        let mut out = Vec::with_capacity(data.len());
        for chunk in data.chunks_exact(16) {
            let mut b = [0u8; 16];
            b.copy_from_slice(chunk);
            let cipher = b;
            self.decrypt_block(&mut b);
            for i in 0..16 {
                b[i] ^= prev[i];
            }
            out.extend_from_slice(&b);
            prev = cipher;
        }
        out
    }
}

fn xor16(a: &mut [u8; 16], k: &[u8; 16]) {
    for i in 0..16 {
        a[i] ^= k[i];
    }
}

fn shift_rows(s: &mut [u8; 16]) {
    let t = *s;
    for r in 1..4 {
        for c in 0..4 {
            s[r + 4 * c] = t[r + 4 * ((c + r) % 4)];
        }
    }
}

fn inv_shift_rows(s: &mut [u8; 16]) {
    let t = *s;
    for r in 1..4 {
        for c in 0..4 {
            s[r + 4 * ((c + r) % 4)] = t[r + 4 * c];
        }
    }
}

/// AES-CBC as PDF stores it: a 16-byte IV, then the blocks, PKCS#5 padded (7.6.2). A damaged pad is kept as data.
fn aes_pdf_decrypt(key: &[u8], data: &[u8]) -> Vec<u8> {
    let Some(aes) = Aes::new(key) else { return data.to_vec() };
    if data.len() < 16 {
        return Vec::new();
    }
    let mut iv = [0u8; 16];
    iv.copy_from_slice(&data[..16]);
    let body = &data[16..];
    let whole = body.len() / 16 * 16;
    let mut out = aes.cbc_decrypt(&iv, &body[..whole]);
    if let Some(&pad) = out.last() {
        let pad = pad as usize;
        if (1..=16).contains(&pad) && pad <= out.len() && out[out.len() - pad..].iter().all(|&b| b as usize == pad) {
            out.truncate(out.len() - pad);
        }
    }
    out
}

// ---- The security handler ------------------------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Method {
    Identity,
    Rc4,
    AesV2,
    AesV3,
}

pub struct Crypt {
    key: Vec<u8>,
    /// Method for streams and strings without their own crypt filter.
    pub stream: Method,
    pub string: Method,
    /// Named crypt filters from /CF.
    filters: Vec<(Vec<u8>, Method)>,
    pub encrypt_metadata: bool,
    pub description: String,
}

const PAD: [u8; 32] = [
    0x28, 0xBF, 0x4E, 0x5E, 0x4E, 0x75, 0x8A, 0x41, 0x64, 0x00, 0x4E, 0x56, 0xFF, 0xFA, 0x01, 0x08, 0x2E, 0x2E, 0x00,
    0xB6, 0xD0, 0x68, 0x3E, 0x80, 0x2F, 0x0C, 0xA9, 0xFE, 0x64, 0x53, 0x69, 0x7A,
];

fn password_refusal() -> KernelError {
    KernelError::Unsupported("this PDF needs a password to open".into())
}

fn int(d: &Dict, key: &[u8]) -> Option<i64> {
    d.get(key).and_then(Obj::as_i64)
}

fn bytes<'a>(d: &'a Dict, key: &[u8]) -> &'a [u8] {
    d.get(key).and_then(Obj::as_str).unwrap_or(&[])
}

impl Crypt {
    /// Sets up decryption from the encryption dictionary (direct objects only; the caller resolves references) and
    /// the first element of the trailer's /ID.
    pub fn new(enc: &Dict, id0: &[u8]) -> Result<Crypt, KernelError> {
        let filter = enc.get(b"Filter").and_then(Obj::as_name).unwrap_or(b"");
        if filter != b"Standard" {
            return Err(KernelError::Unsupported(format!(
                "this PDF is encrypted with the '{}' security handler (only the standard password handler is read)",
                String::from_utf8_lossy(filter))));
        }
        let v = int(enc, b"V").unwrap_or(0);
        let r = int(enc, b"R").unwrap_or(2);
        let p = int(enc, b"P").unwrap_or(-1) as i32;
        let encrypt_metadata = !matches!(enc.get(b"EncryptMetadata"), Some(Obj::Bool(false)));
        let o = bytes(enc, b"O");
        let u = bytes(enc, b"U");

        // Crypt filters (V 4 and 5).
        let mut filters = Vec::new();
        let mut cf_length: Option<i64> = None;
        if let Some(Obj::Dict(cf)) = enc.get(b"CF") {
            for (name, value) in &cf.0 {
                let Some(fd) = value.as_dict() else { continue };
                let method = match fd.get(b"CFM").and_then(Obj::as_name) {
                    Some(b"V2") => Method::Rc4,
                    Some(b"AESV2") => Method::AesV2,
                    Some(b"AESV3") => Method::AesV3,
                    Some(b"None") | None => Method::Identity,
                    Some(other) => return Err(KernelError::Unsupported(format!(
                        "this PDF is encrypted with the unknown method '{}'", String::from_utf8_lossy(other)))),
                };
                if let Some(l) = int(fd, b"Length") {
                    cf_length = Some(l);
                }
                filters.push((name.clone(), method));
            }
        }
        let lookup = |key: &[u8]| -> Method {
            match enc.get(key).and_then(Obj::as_name) {
                None | Some(b"Identity") => Method::Identity,
                Some(n) => filters.iter().find(|(k, _)| k.as_slice() == n).map_or(Method::Identity, |(_, m)| *m),
            }
        };
        let (stream, string) = match v {
            1 | 2 | 3 => (Method::Rc4, Method::Rc4),
            4 | 5 => (lookup(b"StmF"), lookup(b"StrF")),
            _ => return Err(KernelError::Unsupported(format!("this PDF uses encryption version {v}, which is not read"))),
        };

        let key = if r >= 5 {
            aes256_key(r, enc, o, u)?
        } else {
            // Key length n in bytes: 5 for revision 2, else /Length bits / 8 (a crypt filter's Length is in bytes
            // when 16 or less, else bits).
            let mut n = if r == 2 { 5 } else {
                int(enc, b"Length").map(|l| l / 8).or(cf_length.map(|l| if l <= 16 { l } else { l / 8 })).unwrap_or(if v >= 4 { 16 } else { 5 })
            };
            if !(5..=16).contains(&n) {
                n = 16;
            }
            let n = n as usize;
            let compute = |password: &[u8]| -> Vec<u8> {
                let mut input = padded(password).to_vec();
                input.extend_from_slice(&o[..o.len().min(32)]);
                input.extend_from_slice(&(p as u32).to_le_bytes());
                input.extend_from_slice(id0);
                if r >= 4 && !encrypt_metadata {
                    input.extend_from_slice(&[0xFF; 4]);
                }
                let mut hash = md5(&input);
                if r >= 3 {
                    for _ in 0..50 {
                        hash = md5(&hash[..n]);
                    }
                }
                hash[..n].to_vec()
            };
            let user_ok = |key: &[u8]| -> bool {
                if r == 2 {
                    rc4(key, &PAD).as_slice() == &u[..u.len().min(32)]
                } else {
                    let mut input = PAD.to_vec();
                    input.extend_from_slice(id0);
                    let mut x = md5(&input).to_vec();
                    for i in 0..20u8 {
                        let k: Vec<u8> = key.iter().map(|b| b ^ i).collect();
                        x = rc4(&k, &x);
                    }
                    u.len() >= 16 && x[..16] == u[..16]
                }
            };
            let empty_user = compute(b"");
            if user_ok(&empty_user) {
                empty_user
            } else {
                // Algorithm 7 with an empty owner password: recover the user password from O.
                let mut hash = md5(&PAD);
                if r >= 3 {
                    for _ in 0..50 {
                        hash = md5(&hash);
                    }
                }
                let okey = &hash[..n];
                let mut user = o[..o.len().min(32)].to_vec();
                if r == 2 {
                    user = rc4(okey, &user);
                } else {
                    for i in (0..20u8).rev() {
                        let k: Vec<u8> = okey.iter().map(|b| b ^ i).collect();
                        user = rc4(&k, &user);
                    }
                }
                let key = compute(&user);
                if !user_ok(&key) {
                    return Err(password_refusal());
                }
                key
            }
        };
        let name = |m: Method| match m {
            Method::Identity => "none".to_string(),
            Method::Rc4 => format!("RC4 {}-bit", key.len() * 8),
            Method::AesV2 => "AES 128-bit".to_string(),
            Method::AesV3 => "AES 256-bit".to_string(),
        };
        let description = format!("the PDF is encrypted ({}, security revision {r}) but opens without a password; it was decrypted", name(stream));
        Ok(Crypt { key, stream, string, filters, encrypt_metadata, description })
    }

    pub fn filter(&self, name: &[u8]) -> Method {
        if name == b"Identity" {
            return Method::Identity;
        }
        self.filters.iter().find(|(k, _)| k.as_slice() == name).map_or(Method::Identity, |(_, m)| *m)
    }

    /// Algorithm 1: the object key, then RC4 or AES.
    pub fn decrypt(&self, method: Method, num: u32, gen: u16, data: &[u8]) -> Vec<u8> {
        match method {
            Method::Identity => data.to_vec(),
            Method::AesV3 => aes_pdf_decrypt(&self.key, data),
            Method::Rc4 | Method::AesV2 => {
                let mut input = self.key.clone();
                input.extend_from_slice(&num.to_le_bytes()[..3]);
                input.extend_from_slice(&gen.to_le_bytes());
                if method == Method::AesV2 {
                    input.extend_from_slice(b"sAlT");
                }
                let hash = md5(&input);
                let k = &hash[..(self.key.len() + 5).min(16)];
                if method == Method::Rc4 { rc4(k, data) } else { aes_pdf_decrypt(k, data) }
            }
        }
    }
}

fn padded(password: &[u8]) -> [u8; 32] {
    let mut out = [0u8; 32];
    let n = password.len().min(32);
    out[..n].copy_from_slice(&password[..n]);
    out[n..].copy_from_slice(&PAD[..32 - n]);
    out
}

/// ISO 32000-2 Algorithm 2.B (revision 6) or the plain SHA-256 of revision 5.
pub fn hash_r6(r: i64, password: &[u8], salt: &[u8], udata: &[u8]) -> [u8; 32] {
    let mut input = password.to_vec();
    input.extend_from_slice(salt);
    input.extend_from_slice(udata);
    let first = sha256(&input);
    if r < 6 {
        return first;
    }
    let mut k: Vec<u8> = first.to_vec();
    let mut round = 0u32;
    loop {
        let mut k1 = Vec::with_capacity(64 * (password.len() + k.len() + udata.len()));
        for _ in 0..64 {
            k1.extend_from_slice(password);
            k1.extend_from_slice(&k);
            k1.extend_from_slice(udata);
        }
        let Some(aes) = Aes::new(&k[..16]) else { break };
        let mut iv = [0u8; 16];
        iv.copy_from_slice(&k[16..32]);
        let e = aes.cbc_encrypt(&iv, &k1);
        let sum: u32 = e[..16].iter().map(|&b| b as u32).sum();
        k = match sum % 3 {
            0 => sha256(&e).to_vec(),
            1 => sha384(&e),
            _ => sha512(&e),
        };
        round += 1;
        let last = *e.last().unwrap_or(&0) as u32;
        if round >= 64 && last <= round - 32 {
            break;
        }
        if round > 10_000 {
            break;
        }
    }
    let mut out = [0u8; 32];
    out.copy_from_slice(&k[..32]);
    out
}

/// Algorithm 2.A: the file key from an empty user password (or an empty owner password).
fn aes256_key(r: i64, enc: &Dict, o: &[u8], u: &[u8]) -> Result<Vec<u8>, KernelError> {
    if o.len() < 48 || u.len() < 48 {
        return Err(KernelError::InvalidInput("the PDF's AES-256 encryption dictionary is damaged (short O or U)".into()));
    }
    let ue = bytes(enc, b"UE");
    let oe = bytes(enc, b"OE");
    let unwrap = |intermediate: [u8; 32], wrapped: &[u8]| -> Option<Vec<u8>> {
        if wrapped.len() < 32 {
            return None;
        }
        let aes = Aes::new(&intermediate)?;
        Some(aes.cbc_decrypt(&[0u8; 16], &wrapped[..32]))
    };
    if hash_r6(r, b"", &u[32..40], &[]) == u[..32] {
        return unwrap(hash_r6(r, b"", &u[40..48], &[]), ue)
            .ok_or_else(|| KernelError::InvalidInput("the PDF's AES-256 encryption dictionary has no UE key".into()));
    }
    let udata = &u[..48];
    if hash_r6(r, b"", &o[32..40], udata) == o[..32] {
        return unwrap(hash_r6(r, b"", &o[40..48], udata), oe)
            .ok_or_else(|| KernelError::InvalidInput("the PDF's AES-256 encryption dictionary has no OE key".into()));
    }
    Err(password_refusal())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(b: &[u8]) -> String {
        b.iter().map(|x| format!("{x:02x}")).collect()
    }

    #[test]
    fn md5_rfc_1321_test_suite() {
        assert_eq!(hex(&md5(b"")), "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(hex(&md5(b"abc")), "900150983cd24fb0d6963f7d28e17f72");
        assert_eq!(hex(&md5(b"12345678901234567890123456789012345678901234567890123456789012345678901234567890")),
            "57edf4a22be3c955ac49da2e2107b67a");
    }

    #[test]
    fn sha2_fips_examples() {
        assert_eq!(hex(&sha256(b"abc")), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        assert_eq!(hex(&sha384(b"abc")),
            "cb00753f45a35e8bb5a03d699ac65007272c32ab0eded1631a8b605a43ff5bed8086072ba1e7cc2358baeca134c825a7");
        assert_eq!(hex(&sha512(b"abc")), concat!("ddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a",
            "2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f"));
    }

    #[test]
    fn aes_fips_197_examples() {
        let key = [0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c];
        let aes = Aes::new(&key).unwrap();
        let mut b = [0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34];
        let plain = b;
        aes.encrypt_block(&mut b);
        assert_eq!(hex(&b), "3925841d02dc09fbdc118597196a0b32");
        aes.decrypt_block(&mut b);
        assert_eq!(b, plain);
        // A.3: the 256-bit key expansion ends with w59 = 706c631e.
        let k256: Vec<u8> = (0..32).map(|i| u8::from_str_radix(&"603deb1015ca71be2b73aef0857d77811f352c073b6108d72d9810a30914dff4"[2 * i..2 * i + 2], 16).unwrap()).collect();
        let aes = Aes::new(&k256).unwrap();
        assert_eq!(hex(&aes.last_word()), "706c631e");
        let mut b = plain;
        aes.encrypt_block(&mut b);
        aes.decrypt_block(&mut b);
        assert_eq!(b, plain);
    }

    #[test]
    fn rc4_round_trips() {
        let c = rc4(b"Key", b"Plaintext");
        assert_eq!(hex(&c), "bbf316e8d940af0ad3");
        assert_eq!(rc4(b"Key", &c), b"Plaintext");
    }
}
