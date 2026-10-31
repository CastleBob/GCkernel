//! The DWG bit stream (ODA specification chapter 2, "Bit codes and data definitions"): values packed most
//! significant bit first, little-endian bytes, the compressed codes B, BB, BS, BL, BLL, BD, DD, MC, MS, H, T, TU.
//! Every read is checked against the stream's end; nothing here can panic on damaged data.

use crate::KernelError;

pub(super) type R<T> = Result<T, KernelError>;

pub(super) fn bad<T>(message: impl Into<String>) -> R<T> {
    Err(KernelError::InvalidInput(format!("DWG: {}", message.into())))
}

/// The file format generations the reader tells apart (each later one includes the earlier changes).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Ver {
    R13,
    R14,
    R2000,
    R2004,
    R2007,
    R2010,
    R2013,
    R2018,
}

/// A handle reference: |code (4 bits)|counter (4 bits)|value bytes|, with offset codes (6, 8, A, C) resolved.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Handle {
    pub code: u8,
    pub value: u64,
}

#[derive(Clone)]
pub(super) struct Bits<'a> {
    data: &'a [u8],
    /// Next bit to read, counted from the start of `data`.
    pub pos: usize,
    /// One past the last readable bit.
    pub end: usize,
    pub ver: Ver,
}

impl<'a> Bits<'a> {
    pub fn new(data: &'a [u8], ver: Ver) -> Self {
        Bits { data, pos: 0, end: data.len().saturating_mul(8), ver }
    }

    /// The same bytes read from bit `start` up to bit `end` (clamped to the data).
    pub fn window(&self, start: usize, end: usize) -> Bits<'a> {
        let limit = self.data.len().saturating_mul(8);
        let end = end.min(limit);
        Bits { data: self.data, pos: start.min(end), end, ver: self.ver }
    }

    pub fn left(&self) -> usize {
        self.end.saturating_sub(self.pos)
    }

    fn need(&self, n: usize) -> R<()> {
        if self.pos.checked_add(n).map_or(true, |e| e > self.end) {
            return bad(format!("data ends early (bit {} of {})", self.pos, self.end));
        }
        Ok(())
    }

    pub fn skip(&mut self, n: usize) -> R<()> {
        self.need(n)?;
        self.pos += n;
        Ok(())
    }

    pub fn b(&mut self) -> R<bool> {
        self.need(1)?;
        let byte = self.data[self.pos >> 3];
        let bit = (byte >> (7 - (self.pos & 7))) & 1;
        self.pos += 1;
        Ok(bit == 1)
    }

    /// `n` (at most 8) bits, most significant first.
    fn small(&mut self, n: usize) -> R<u8> {
        self.need(n)?;
        let mut v = 0u8;
        for _ in 0..n {
            let byte = self.data[self.pos >> 3];
            v = (v << 1) | ((byte >> (7 - (self.pos & 7))) & 1);
            self.pos += 1;
        }
        Ok(v)
    }

    pub fn bb(&mut self) -> R<u8> {
        self.small(2)
    }


    pub fn rc(&mut self) -> R<u8> {
        self.need(8)?;
        let at = self.pos >> 3;
        let shift = self.pos & 7;
        let v = if shift == 0 {
            self.data[at]
        } else {
            let hi = self.data[at] << shift;
            let lo = self.data.get(at + 1).copied().unwrap_or(0) >> (8 - shift);
            hi | lo
        };
        self.pos += 8;
        Ok(v)
    }

    fn bytes<const N: usize>(&mut self) -> R<[u8; N]> {
        self.need(8 * N)?;
        let mut out = [0u8; N];
        for b in out.iter_mut() {
            *b = self.rc()?;
        }
        Ok(out)
    }

    pub fn rs(&mut self) -> R<u16> {
        Ok(u16::from_le_bytes(self.bytes::<2>()?))
    }

    pub fn rl(&mut self) -> R<u32> {
        Ok(u32::from_le_bytes(self.bytes::<4>()?))
    }

    pub fn rd(&mut self) -> R<f64> {
        Ok(f64::from_le_bytes(self.bytes::<8>()?))
    }

    pub fn bs(&mut self) -> R<i16> {
        Ok(match self.bb()? {
            0 => self.rs()? as i16,
            1 => self.rc()? as i16,
            2 => 0,
            _ => 256,
        })
    }

    pub fn bl(&mut self) -> R<i32> {
        match self.bb()? {
            0 => Ok(self.rl()? as i32),
            1 => Ok(self.rc()? as i32),
            2 => Ok(0),
            _ => bad(format!("bit long with the unused code 11 at bit {}", self.pos)),
        }
    }

    /// A bit long used as a count of things that follow, each at least `min_bits` long: refused when the rest of
    /// the stream could not hold that many.
    pub fn count(&mut self, min_bits: usize, what: &str) -> R<usize> {
        let n = self.bl()?;
        self.check_count(n as i64, min_bits, what)
    }

    pub fn check_count(&self, n: i64, min_bits: usize, what: &str) -> R<usize> {
        if n < 0 || (n as u64).saturating_mul(min_bits.max(1) as u64) > self.left() as u64 {
            return bad(format!("{what} count {n} does not fit the data"));
        }
        Ok(n as usize)
    }

    pub fn bll(&mut self) -> R<u64> {
        // The byte count is a plain 3-bit number (files show "001" for one byte, which the 1-3 bit code would read as 0).
        let n = self.small(3)? as usize;
        let mut v = 0u64;
        for i in 0..n {
            v |= (self.rc()? as u64) << (8 * i);
        }
        Ok(v)
    }

    pub fn bd(&mut self) -> R<f64> {
        match self.bb()? {
            0 => self.rd(),
            1 => Ok(1.0),
            2 => Ok(0.0),
            _ => bad(format!("bit double with the unused code 11 at bit {}", self.pos)),
        }
    }

    /// Bit double with default.
    pub fn dd(&mut self, default: f64) -> R<f64> {
        let mut b = default.to_le_bytes();
        match self.bb()? {
            0 => {}
            1 => {
                for k in b.iter_mut().take(4) {
                    *k = self.rc()?;
                }
            }
            2 => {
                b[4] = self.rc()?;
                b[5] = self.rc()?;
                for k in b.iter_mut().take(4) {
                    *k = self.rc()?;
                }
            }
            _ => return self.rd(),
        }
        Ok(f64::from_le_bytes(b))
    }

    pub fn p3(&mut self) -> R<[f64; 3]> {
        Ok([self.bd()?, self.bd()?, self.bd()?])
    }

    /// Bit extrusion.
    pub fn be(&mut self) -> R<[f64; 3]> {
        if self.ver >= Ver::R2000 && self.b()? {
            return Ok([0.0, 0.0, 1.0]);
        }
        self.p3()
    }

    /// Bit thickness.
    pub fn bt(&mut self) -> R<f64> {
        if self.ver >= Ver::R2000 && self.b()? {
            return Ok(0.0);
        }
        self.bd()
    }

    /// Object type: a bit short until R2007, then a bit pair and one or two bytes.
    pub fn ot(&mut self) -> R<u16> {
        if self.ver < Ver::R2010 {
            return Ok(self.bs()? as u16);
        }
        Ok(match self.bb()? {
            0 => self.rc()? as u16,
            1 => self.rc()? as u16 + 0x1f0,
            _ => self.rs()?,
        })
    }

    /// A handle reference; offset codes are resolved against `reference` (the object's own handle).
    pub fn h(&mut self, reference: u64) -> R<Handle> {
        let first = self.rc()?;
        let code = first >> 4;
        let counter = (first & 0x0f) as usize;
        if counter > 8 {
            return bad(format!("a handle {counter} bytes long at bit {}", self.pos));
        }
        let mut v = 0u64;
        for _ in 0..counter {
            v = (v << 8) | self.rc()? as u64;
        }
        let value = match code {
            6 => reference.wrapping_add(1),
            8 => reference.wrapping_sub(1),
            0xA => reference.wrapping_add(v),
            0xC => reference.wrapping_sub(v),
            _ => v,
        };
        Ok(Handle { code, value })
    }

    /// Handle stream items whose value is not needed.
    pub fn skip_h(&mut self, n: usize) -> R<()> {
        for _ in 0..n {
            self.h(0)?;
        }
        Ok(())
    }

    /// Text: T before R2007 (bit short length, then code-page bytes), TU from R2007 (bit short length, then UTF-16).
    pub fn text(&mut self) -> R<String> {
        let n = self.bs()? as u16 as usize;
        if self.ver >= Ver::R2007 {
            self.need(16 * n)?;
            let units: Vec<u16> = (0..n).map(|_| self.rs()).collect::<R<_>>()?;
            let s = String::from_utf16_lossy(&units);
            return Ok(s.trim_end_matches('\0').to_string());
        }
        self.need(8 * n)?;
        let bytes: Vec<u8> = (0..n).map(|_| self.rc()).collect::<R<_>>()?;
        let bytes = match bytes.iter().position(|&b| b == 0) {
            Some(z) => &bytes[..z],
            None => &bytes[..],
        };
        Ok(crate::dxf::codes::decode(bytes))
    }

    /// CmColor as stored in tables and the header: a bit short index before R2004; from R2004 index, RGB, flags byte
    /// and optional names. Gives the AutoCAD colour index (negative indexes are kept: a layer that is off).
    pub fn cmc(&mut self, strings: &mut Option<Bits<'a>>) -> R<i32> {
        let index = self.bs()? as i32;
        if self.ver < Ver::R2004 {
            return Ok(index);
        }
        let rgb = self.bl()? as u32;
        let flags = self.rc()?;
        if flags & 1 != 0 {
            self.tv(strings)?;
        }
        if flags & 2 != 0 {
            self.tv(strings)?;
        }
        Ok(if index != 0 { index } else { method_index(rgb) })
    }

    /// TV: text in the object data before R2007, in the string stream from R2007.
    pub fn tv(&mut self, strings: &mut Option<Bits<'a>>) -> R<String> {
        if self.ver >= Ver::R2007 {
            return match strings.as_mut() {
                Some(s) => s.text(),
                None => Ok(String::new()),
            };
        }
        self.text()
    }
}

/// The colour index carried by a 32-bit colour value (its high byte is how the colour is given: 0xC0 by layer,
/// 0xC1 by block, 0xC2 true colour, 0xC3 colour index in the low byte).
pub(super) fn method_index(rgb: u32) -> i32 {
    match rgb >> 24 {
        0xC0 => 256,
        0xC1 => 0,
        0xC3 => (rgb & 0xff) as i32,
        0xC2 => primary(rgb & 0xff_ffff).unwrap_or(256),
        _ => 7,
    }
}

/// Modular char (object map, R2010 handle stream size): 7 bits per byte, low first, high bit set while more follow.
/// When `signed`, bit 0x40 of the last byte negates the value.
pub(super) fn mc(data: &[u8], at: &mut usize, signed: bool) -> R<i64> {
    let mut v: u64 = 0;
    let mut shift = 0u32;
    loop {
        let Some(&byte) = data.get(*at) else { return bad("modular char runs off the data") };
        *at += 1;
        if byte & 0x80 != 0 {
            if shift > 56 {
                return bad("modular char too long");
            }
            v |= ((byte & 0x7f) as u64) << shift;
            shift += 7;
            continue;
        }
        if signed {
            if shift > 57 {
                return bad("modular char too long");
            }
            v |= ((byte & 0x3f) as u64) << shift;
            let v = v as i64;
            return Ok(if byte & 0x40 != 0 { -v } else { v });
        }
        if shift > 57 {
            return bad("modular char too long");
        }
        v |= ((byte & 0x7f) as u64) << shift;
        return Ok(v as i64);
    }
}

/// Modular short (object sizes): 15 bits per little-endian short, high bit set while more follow.
pub(super) fn ms(data: &[u8], at: &mut usize) -> R<usize> {
    let mut v: u64 = 0;
    let mut shift = 0u32;
    loop {
        let (Some(&lo), Some(&hi)) = (data.get(*at), data.get(*at + 1)) else {
            return bad("object size runs off the data");
        };
        *at += 2;
        let word = u16::from_le_bytes([lo, hi]) as u64;
        v |= (word & 0x7fff) << shift;
        if word & 0x8000 == 0 {
            return Ok(v as usize);
        }
        shift += 15;
        if shift > 45 {
            return bad("object size too long");
        }
    }
}

pub(super) fn le32(data: &[u8], at: usize) -> R<u32> {
    match data.get(at..at.wrapping_add(4)) {
        Some(b) if at.checked_add(4).is_some() => Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
        _ => bad(format!("data ends early at byte {at}")),
    }
}

pub(super) fn le64(data: &[u8], at: usize) -> R<u64> {
    Ok(le32(data, at)? as u64 | (le32(data, at.wrapping_add(4))? as u64) << 32)
}


impl<'a> Bits<'a> {
    /// The start of every table record (layer, block header, ...): name, the 64 flag, xref index + 1 and the
    /// xref-dependent bit. Files from R2007 on store only the xref index (the two flags are not in the data).
    pub fn table_entry(&mut self, strings: &mut Option<Bits<'a>>) -> R<String> {
        let name = self.tv(strings)?;
        if self.ver < Ver::R2007 {
            self.b()?;
        }
        self.bs()?;
        if self.ver < Ver::R2007 {
            self.b()?;
        }
        Ok(name)
    }
}

/// The colour index of a true colour that is exactly one of the seven standard colours (1 red, 2 yellow, 3 green,
/// 4 cyan, 5 blue, 6 magenta, 7 white or black); other true colours have no exact index.
fn primary(rgb: u32) -> Option<i32> {
    Some(match rgb {
        0xFF0000 => 1,
        0xFFFF00 => 2,
        0x00FF00 => 3,
        0x00FFFF => 4,
        0x0000FF => 5,
        0xFF00FF => 6,
        0xFFFFFF | 0x000000 => 7,
        _ => return None,
    })
}
