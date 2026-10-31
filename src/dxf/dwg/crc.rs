//! The R2007 integrity values (ODA specification 5.12 and 5.4.1): 64-bit CRCs of the file header, the system pages
//! and the data pages, and the 32-bit checksum of each data page's decompressed data.
//!
//! The specification prints both 256-entry CRC tables. They are generated here from their polynomials (public
//! mathematics) and checked against the printed tables entry for entry (tests below, and in the spec text):
//! the normal CRC is ECMA-182, polynomial 0x42F0E1EBA9EA3693, most significant bit first; the "mirrored" table is
//! the least-significant-bit-first table of 0x95AC9329AC4BC9B5 (the bit reversal of 0xAD93D23594C935A9).
//!
//! What each value covers and its seed were found by checking the real AC1021 sample (fixture_2007.dwg), where all
//! of them match:
//! - file header (decompressed, 0x110 bytes, its CRC field at 0x108 taken as zero): normal CRC, seed
//!   UpdateSeed2(0, length) (spec 5.2.1.2);
//! - the file header's stored (compressed) data: normal CRC, seed UpdateSeed2(0, length) (5.2.1.3);
//! - page map and section map: mirrored CRC of the stored (compressed) bytes and of the decompressed data, seed
//!   UpdateSeed1(CrcSeed, length) with CrcSeed from the file header (5.3; the sample's CrcSeed is 0);
//! - data pages: "Page CRC" = mirrored CRC of the stored (compressed) bytes, seed UpdateSeed1(CrcSeed, length)
//!   (5.4); "Page Checksum" = GetCheckSum(0, decompressed page data) (5.4.1; the spec leaves the seed open, 0 fits).

use std::sync::OnceLock;

const NORMAL: u64 = 0x42F0_E1EB_A9EA_3693;
const MIRRORED: u64 = 0x95AC_9329_AC4B_C9B5;

fn tables() -> &'static ([u64; 256], [u64; 256]) {
    static TABLES: OnceLock<([u64; 256], [u64; 256])> = OnceLock::new();
    TABLES.get_or_init(|| {
        let mut normal = [0u64; 256];
        let mut mirrored = [0u64; 256];
        for i in 0..256u64 {
            let mut n = i << 56;
            let mut m = i;
            for _ in 0..8 {
                n = if n >> 63 != 0 { (n << 1) ^ NORMAL } else { n << 1 };
                m = if m & 1 != 0 { (m >> 1) ^ MIRRORED } else { m >> 1 };
            }
            normal[i as usize] = n;
            mirrored[i as usize] = m;
        }
        (normal, mirrored)
    })
}

/// The order the CRC and the checksum visit bytes in (spec 5.12's table, the same in 5.4.1): each 8-byte block as
/// bytes 6 7 4 5 2 3 0 1, then the 1-7 remaining bytes by the table for their count.
fn visit(data: &[u8], mut f: impl FnMut(u8)) {
    let mut blocks = data.chunks_exact(8);
    for b in &mut blocks {
        for i in [6, 7, 4, 5, 2, 3, 0, 1] {
            f(b[i]);
        }
    }
    let rest = blocks.remainder();
    let order: &[usize] = match rest.len() {
        0 => &[],
        1 => &[0],
        2 => &[0, 1],
        3 => &[0, 1, 2],
        4 => &[2, 3, 0, 1],
        5 => &[2, 3, 0, 1, 4],
        6 => &[2, 3, 0, 1, 4, 5],
        _ => &[2, 3, 0, 1, 4, 5, 6],
    };
    for &i in order {
        f(rest[i]);
    }
}

/// Normal 64-bit CRC (5.12.1), inverted at the end.
pub(super) fn normal(data: &[u8], init: u64) -> u64 {
    let table = &tables().0;
    let mut crc = init;
    visit(data, |b| crc = table[((u64::from(b) ^ (crc >> 56)) & 0xff) as usize] ^ (crc << 8));
    !crc
}

/// Mirrored 64-bit CRC (5.12.2), not inverted.
pub(super) fn mirrored(data: &[u8], init: u64) -> u64 {
    let table = &tables().1;
    let mut crc = init;
    visit(data, |b| crc = table[((crc ^ u64::from(b)) & 0xff) as usize] ^ (crc >> 8));
    crc
}

/// UpdateSeed1 (5.12) in 64-bit arithmetic, and the same with its first step in 32 bits: the specification's type
/// ("UInt3264") leaves that open and the sample cannot tell them apart (they differ only for lengths over 20116
/// bytes), so a value matching either is accepted.
pub(super) fn seed1(seed: u64, length: usize) -> [u64; 2] {
    let rest = |s: u64| !(s | s.wrapping_mul(0x343fd << 32).wrapping_add(0x269ec3 << 32));
    let s = seed.wrapping_add(length as u64).wrapping_mul(0x343fd).wrapping_add(0x269ec3);
    [rest(s), rest(s & 0xffff_ffff)]
}

/// UpdateSeed2 (5.12), 64-bit arithmetic (used for lengths under 0x400, where both readings agree).
pub(super) fn seed2(seed: u64, length: usize) -> u64 {
    let n = length as u64;
    let s = seed.wrapping_add(n).wrapping_mul(0x343fd).wrapping_add(0x269ec3);
    !s.wrapping_mul((1 << 32) + 0x343fd).wrapping_add(n + 0x269ec3)
}

/// True when the mirrored CRC of `data` with seed UpdateSeed1(`seed`, length) is `expected`.
pub(super) fn mirrored_matches(data: &[u8], seed: u64, expected: u64) -> bool {
    seed1(seed, data.len()).iter().any(|&init| mirrored(data, init) == expected)
}

/// The data page checksum (5.4.1): two running sums modulo 0xFFF1 started from the seed.
pub(super) fn checksum(seed: u64, data: &[u8]) -> u32 {
    let seed = seed.wrapping_add(data.len() as u64).wrapping_mul(0x343fd).wrapping_add(0x269ec3);
    let mut sum1 = seed & 0xffff;
    let mut sum2 = (seed >> 16) & 0xffff;
    for chunk in data.chunks(0x15b0) {
        visit(chunk, |b| {
            sum1 += u64::from(b);
            sum2 += sum1;
        });
        sum1 %= 0xfff1;
        sum2 %= 0xfff1;
    }
    ((sum2 << 16) | (sum1 & 0xffff)) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Entries of the tables as printed in the specification (5.12.1, 5.12.2).
    #[test]
    fn tables_match_the_specification() {
        let (n, m) = tables();
        for (i, v) in [(1, 0x42f0e1eba9ea3693u64), (2, 0x85e1c3d753d46d26), (3, 0xc711223cfa3e5bb5),
            (16, 0x663d78ff90e185ef), (255, 0x9afce626ce85b507)] {
            assert_eq!(n[i], v, "normal {i}");
        }
        for (i, v) in [(1, 0x7ad870c830358979u64), (2, 0xf5b0e190606b12f2), (3, 0x8f689158505e9b8b),
            (16, 0x7d08ff3b88be6f81), (64, 0xdf7adabd7a6e2d6f), (255, 0x29b7d047efec8728)] {
            assert_eq!(m[i], v, "mirrored {i}");
        }
    }

    #[test]
    fn byte_order_and_checksum() {
        // The visiting order on 11 bytes: one block of 8, then 3.
        let mut seen = Vec::new();
        visit(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], |b| seen.push(b));
        assert_eq!(seen, vec![6, 7, 4, 5, 2, 3, 0, 1, 8, 9, 10]);
        // Empty data: only the seeded sums.
        assert_eq!(checksum(0, &[]), 0x0026_9ec3);
    }
}
