//! Deflate decompression (RFC 1951) and the zlib wrapper (RFC 1950), written from the RFC texts: stored, fixed and
//! dynamic Huffman blocks, canonical codes decoded bit by bit, the Adler-32 check.

use crate::KernelError;

struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
    buf: u64,
    count: u32,
}

/// Input ran out in the middle of the data.
const SHORT: &str = "ended early";

impl Bits<'_> {
    /// n (<= 32) bits, least significant first (RFC 1951 3.1.1).
    fn take(&mut self, n: u32) -> Result<u32, &'static str> {
        while self.count < n {
            let Some(&b) = self.data.get(self.pos) else { return Err(SHORT) };
            self.pos += 1;
            self.buf |= (b as u64) << self.count;
            self.count += 8;
        }
        let v = (self.buf & ((1u64 << n) - 1)) as u32;
        self.buf >>= n;
        self.count -= n;
        Ok(v)
    }

    fn align(&mut self) {
        let drop = self.count % 8;
        self.buf >>= drop;
        self.count -= drop;
    }
}

/// A canonical Huffman code (RFC 1951 3.2.2): the number of codes of each length and the symbols in code order.
struct Huffman {
    counts: [u16; 16],
    symbols: Vec<u16>,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Result<Huffman, &'static str> {
        let mut counts = [0u16; 16];
        for &l in lengths {
            counts[l as usize] += 1;
        }
        counts[0] = 0;
        // No length may be over-subscribed.
        let mut left: i32 = 1;
        for &c in &counts[1..] {
            left = left * 2 - c as i32;
            if left < 0 {
                return Err("has an over-subscribed Huffman code");
            }
        }
        let mut offsets = [0u16; 16];
        for len in 1..15 {
            offsets[len + 1] = offsets[len] + counts[len];
        }
        let mut symbols = vec![0u16; lengths.len()];
        for (sym, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbols[offsets[l as usize] as usize] = sym as u16;
                offsets[l as usize] += 1;
            }
        }
        Ok(Huffman { counts, symbols })
    }

    /// Codes of one length are consecutive integers, and shorter codes come before longer ones; the code's bits
    /// arrive most significant first.
    fn decode(&self, bits: &mut Bits) -> Result<u16, &'static str> {
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for len in 1..16 {
            code |= bits.take(1)? as i32;
            let count = self.counts[len] as i32;
            if code - first < count {
                return self.symbols.get((index + code - first) as usize).copied().ok_or("has a bad Huffman code");
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        Err("has a bad Huffman code")
    }
}

const LENGTH_BASE: [u16; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LENGTH_EXTRA: [u8; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DIST_BASE: [u16; 30] = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
const DIST_EXTRA: [u8; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

/// How a decode ended.
pub struct Inflated {
    pub data: Vec<u8>,
    /// Bytes of input used (for finding the zlib checksum after the data).
    pub used: usize,
    /// Why the data stopped short of a proper end, if it did.
    pub damage: Option<&'static str>,
}

/// Raw Deflate data. Output past `limit` bytes is an error; damage after some output keeps the output and says why.
pub fn inflate(input: &[u8], limit: usize) -> Result<Inflated, KernelError> {
    let mut bits = Bits { data: input, pos: 0, buf: 0, count: 0 };
    let mut out: Vec<u8> = Vec::with_capacity(input.len().saturating_mul(4).min(limit).min(1 << 26));
    let result = blocks(&mut bits, &mut out, limit);
    let used = bits.pos - (bits.count / 8) as usize;
    match result {
        Ok(()) => Ok(Inflated { data: out, used, damage: None }),
        Err("LIMIT") => Err(KernelError::InvalidInput(format!(
            "a compressed stream in the PDF expands past {} MB; the file was refused", limit >> 20))),
        Err(why) => Ok(Inflated { data: out, used, damage: Some(why) }),
    }
}

fn blocks(bits: &mut Bits, out: &mut Vec<u8>, limit: usize) -> Result<(), &'static str> {
    loop {
        let last = bits.take(1)?;
        match bits.take(2)? {
            0 => {
                bits.align();
                let len = bits.take(16)?;
                let nlen = bits.take(16)?;
                if len != !nlen & 0xFFFF {
                    return Err("has a damaged stored block");
                }
                if out.len() + len as usize > limit {
                    return Err("LIMIT");
                }
                for _ in 0..len {
                    out.push(bits.take(8)? as u8);
                }
            }
            1 => {
                let mut lengths = [0u8; 320];
                lengths[..144].fill(8);
                lengths[144..256].fill(9);
                lengths[256..280].fill(7);
                lengths[280..288].fill(8);
                lengths[288..].fill(5);
                let lit = Huffman::new(&lengths[..288])?;
                let dist = Huffman::new(&lengths[288..])?;
                codes(bits, out, &lit, &dist, limit)?;
            }
            2 => {
                let hlit = bits.take(5)? as usize + 257;
                let hdist = bits.take(5)? as usize + 1;
                let hclen = bits.take(4)? as usize + 4;
                let mut cl = [0u8; 19];
                for &i in &ORDER[..hclen] {
                    cl[i] = bits.take(3)? as u8;
                }
                let clh = Huffman::new(&cl)?;
                let mut lengths = vec![0u8; hlit + hdist];
                let mut i = 0;
                while i < hlit + hdist {
                    let sym = clh.decode(bits)?;
                    let (value, repeat) = match sym {
                        0..=15 => (sym as u8, 1),
                        16 => {
                            if i == 0 {
                                return Err("repeats a code length before the first");
                            }
                            (lengths[i - 1], 3 + bits.take(2)? as usize)
                        }
                        17 => (0, 3 + bits.take(3)? as usize),
                        _ => (0, 11 + bits.take(7)? as usize),
                    };
                    if i + repeat > lengths.len() {
                        return Err("has too many code lengths");
                    }
                    lengths[i..i + repeat].fill(value);
                    i += repeat;
                }
                if lengths[256] == 0 {
                    return Err("has no end-of-block code");
                }
                let lit = Huffman::new(&lengths[..hlit])?;
                let dist = Huffman::new(&lengths[hlit..])?;
                codes(bits, out, &lit, &dist, limit)?;
            }
            _ => return Err("has an invalid block type"),
        }
        if last == 1 {
            return Ok(());
        }
    }
}

fn codes(bits: &mut Bits, out: &mut Vec<u8>, lit: &Huffman, dist: &Huffman, limit: usize) -> Result<(), &'static str> {
    loop {
        let sym = lit.decode(bits)? as usize;
        if sym < 256 {
            if out.len() >= limit {
                return Err("LIMIT");
            }
            out.push(sym as u8);
        } else if sym == 256 {
            return Ok(());
        } else {
            let i = sym - 257;
            if i >= 29 {
                return Err("has an invalid length code");
            }
            let len = LENGTH_BASE[i] as usize + bits.take(LENGTH_EXTRA[i] as u32)? as usize;
            let d = dist.decode(bits)? as usize;
            if d >= 30 {
                return Err("has an invalid distance code");
            }
            let distance = DIST_BASE[d] as usize + bits.take(DIST_EXTRA[d] as u32)? as usize;
            if distance > out.len() {
                return Err("refers back before its start");
            }
            if out.len() + len > limit {
                return Err("LIMIT");
            }
            let start = out.len() - distance;
            for k in 0..len {
                let b = out[start + k];
                out.push(b);
            }
        }
    }
}

fn adler32(data: &[u8]) -> u32 {
    let (mut a, mut b) = (1u32, 0u32);
    for chunk in data.chunks(5552) {
        for &x in chunk {
            a += x as u32;
            b += a;
        }
        a %= 65521;
        b %= 65521;
    }
    (b << 16) | a
}

/// zlib data (RFC 1950). A missing header is read as raw Deflate, a missing or wrong checksum is tolerated, and
/// damage after some output keeps that output; each case returns a note.
pub fn zlib(input: &[u8], limit: usize) -> Result<(Vec<u8>, Option<String>), KernelError> {
    let header_ok = input.len() >= 2 && input[0] & 0x0F == 8 && input[0] >> 4 <= 7
        && ((input[0] as u32) << 8 | input[1] as u32) % 31 == 0 && input[1] & 0x20 == 0;
    let body = if header_ok { &input[2..] } else { input };
    let r = inflate(body, limit)?;
    if let Some(why) = r.damage {
        if r.data.is_empty() {
            return Err(KernelError::InvalidInput(format!("a compressed stream in the PDF {why} and gave no data")));
        }
        let n = r.data.len();
        return Ok((r.data, Some(format!("a compressed stream in the PDF {why}; the {n} bytes before the damage were used"))));
    }
    if !header_ok {
        return Ok((r.data, Some("a compressed stream had no zlib header; it was read as raw deflate data".into())));
    }
    let tail = body.get(r.used..r.used + 4);
    let good = tail.is_some_and(|t| u32::from_be_bytes([t[0], t[1], t[2], t[3]]) == adler32(&r.data));
    let note = (!good).then(|| "a compressed stream's checksum was missing or wrong; its data was used as decoded".to_string());
    Ok((r.data, note))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stored_and_fixed_blocks() {
        // A stored block "abc" (BFINAL 1, BTYPE 00).
        let (d, note) = zlib(&[0x78, 0x01, 0x01, 0x03, 0x00, 0xFC, 0xFF, b'a', b'b', b'c', 0x02, 0x4D, 0x01, 0x27], 1 << 20).unwrap();
        assert_eq!(d, b"abc");
        assert!(note.is_none(), "{note:?}");
        // A fixed-Huffman block of "a": literal 'a' (0x61 -> 8-bit code 0x91), then end of block (7 zero bits).
        let r = inflate(&[0x4B, 0x04, 0x00], 1 << 20).unwrap();
        assert_eq!(r.data, b"a");
        assert!(r.damage.is_none());
    }

    #[test]
    fn raw_deflate_and_a_missing_checksum_are_read_with_a_note() {
        let (d, note) = zlib(&[0x4B, 0x04, 0x00], 1 << 20).unwrap();
        assert_eq!(d, b"a");
        assert!(note.unwrap().contains("no zlib header"));
        let (d, note) = zlib(&[0x78, 0x9C, 0x4B, 0x04, 0x00], 1 << 20).unwrap();
        assert_eq!(d, b"a");
        assert!(note.unwrap().contains("checksum"));
        // Damage after some output keeps that output.
        let (d, note) = zlib(&[0x78, 0x01, 0x00, 0x03, 0x00, 0xFC, 0xFF, b'a', b'b'], 1 << 20).unwrap();
        assert_eq!(d, b"ab");
        assert!(note.unwrap().contains("ended early"));
    }

    #[test]
    fn adler_of_wikipedia() {
        assert_eq!(adler32(b"Wikipedia"), 0x11E60398);
    }
}
