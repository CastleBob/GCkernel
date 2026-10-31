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
        0xC2 => nearest_index(rgb & 0xff_ffff),
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

/// The RGB value of AutoCAD Colour Index `index` (1-255), from the palette's documented structure: 1-7 the named
/// colours (red, yellow, green, cyan, blue, magenta, white), 8 and 9 greys (128 and 192); 10-249 24 hues 15 degrees
/// apart (index / 10 - 1 gives the hue, red first), each in five brightnesses (100, 80, 60, 50, 30 %) at full and
/// half saturation (index % 10: even full, odd half; 0-1 brightest); 250-255 greys from 20 % to 100 % in steps of
/// 16 %. Channel values are the exact HSV values truncated to whole numbers (so 50 % is 127, 30 % is 76).
pub(crate) fn aci_rgb(index: u8) -> [u8; 3] {
    match index {
        0 => [0, 0, 0],
        1 => [255, 0, 0],
        2 => [255, 255, 0],
        3 => [0, 255, 0],
        4 => [0, 255, 255],
        5 => [0, 0, 255],
        6 => [255, 0, 255],
        7 => [255, 255, 255],
        8 => [128, 128, 128],
        9 => [192, 192, 192],
        250..=255 => {
            let v = (255 * (20 + 16 * (index as u32 - 250)) / 100) as u8;
            [v, v, v]
        }
        _ => {
            let hue = (index as u32 - 10) / 10; // 0..23, 15 degrees each
            let shade = index as u32 % 10;
            let brightness = [10, 8, 6, 5, 3][(shade / 2) as usize]; // tenths
            let half = shade % 2 == 1;
            // In 1/80ths of a unit so every value is a whole number before truncation.
            let max = 255 * brightness * 8;
            let min = if half { max / 2 } else { 0 };
            let step = (max - min) * (hue % 4) / 4;
            let (rise, fall) = (min + step, max - step);
            let [r, g, b] = match hue / 4 {
                0 => [max, rise, min],
                1 => [fall, max, min],
                2 => [min, max, rise],
                3 => [min, fall, max],
                4 => [rise, min, max],
                _ => [max, min, fall],
            };
            [(r / 80) as u8, (g / 80) as u8, (b / 80) as u8]
        }
    }
}

/// The AutoCAD Colour Index (1-255) nearest a true colour 0xRRGGBB by Euclidean distance in RGB; ties go to the
/// lower index (so the seven named colours keep their indexes 1-7, and black is 18).
pub(crate) fn nearest_index(rgb: u32) -> i32 {
    let want = [(rgb >> 16) as u8, (rgb >> 8) as u8, rgb as u8];
    let mut best = (u32::MAX, 7);
    for index in 1..=255u8 {
        let c = aci_rgb(index);
        let d: u32 = (0..3).map(|k| {
            let e = i32::from(c[k]) - i32::from(want[k]);
            (e * e) as u32
        }).sum();
        if d < best.0 {
            best = (d, index as i32);
        }
    }
    best.1
}

#[cfg(test)]
mod tests {
    use super::{aci_rgb, method_index, nearest_index};

    #[test]
    fn the_colour_index_palette() {
        let known: [(u8, [u8; 3]); 29] = [
            (1, [255, 0, 0]), (2, [255, 255, 0]), (3, [0, 255, 0]), (4, [0, 255, 255]), (5, [0, 0, 255]),
            (6, [255, 0, 255]), (7, [255, 255, 255]), (8, [128, 128, 128]), (9, [192, 192, 192]),
            (10, [255, 0, 0]), (11, [255, 127, 127]), (12, [204, 0, 0]), (13, [204, 102, 102]), (14, [153, 0, 0]),
            (15, [153, 76, 76]), (16, [127, 0, 0]), (17, [127, 63, 63]), (18, [76, 0, 0]), (19, [76, 38, 38]),
            (20, [255, 63, 0]), (21, [255, 159, 127]), (30, [255, 127, 0]), (40, [255, 191, 0]), (50, [255, 255, 0]),
            (96, [0, 127, 0]), (130, [0, 255, 255]), (170, [0, 0, 255]), (250, [51, 51, 51]), (251, [91, 91, 91]),
        ];
        for (index, rgb) in known {
            assert_eq!(aci_rgb(index), rgb, "ACI {index}");
        }
        let greys: Vec<[u8; 3]> = (250..=255).map(aci_rgb).collect();
        assert_eq!(greys, vec![[51; 3], [91; 3], [132; 3], [173; 3], [214; 3], [255; 3]]);
    }

    #[test]
    fn true_colours_take_the_nearest_index() {
        for index in 1..=255u8 {
            let [r, g, b] = aci_rgb(index);
            let rgb = (u32::from(r) << 16) | (u32::from(g) << 8) | u32::from(b);
            // Exact palette colours give their own index, or a lower one with the same colour (10 is 1, 255 is 7).
            let found = nearest_index(rgb);
            assert!(found <= index as i32 && aci_rgb(found as u8) == [r, g, b], "{index} -> {found}");
        }
        assert_eq!(nearest_index(0x000000), 18);
        assert_eq!(nearest_index(0x008000), 96);
        assert_eq!(nearest_index(0xFE0101), 1);
        assert_eq!(nearest_index(0x808080), 8);
        assert_eq!(method_index(0xC2_000000), 18);
        assert_eq!(method_index(0xC2_FFFFFF), 7);
        assert_eq!(method_index(0xC3_000005), 5);
        assert_eq!(method_index(0xC0_000000), 256);
    }
}
