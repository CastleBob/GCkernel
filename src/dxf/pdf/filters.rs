//! Stream filters (ISO 32000-1 7.4): FlateDecode (see inflate.rs), LZWDecode (7.4.4.2, with EarlyChange),
//! ASCIIHexDecode, ASCII85Decode, RunLengthDecode, and the TIFF 2 and PNG predictors of 7.4.4.4. Image-only filters
//! (DCT, JPX, CCITT, JBIG2) are refused: images are never decoded.

use super::inflate;
use super::lexer::{Dict, Obj};
use crate::KernelError;

fn bad(m: String) -> KernelError {
    KernelError::InvalidInput(m)
}

/// One filter with its (already resolved) parameters.
pub fn apply(name: &[u8], parms: Option<&Dict>, data: &[u8], limit: usize, notes: &mut Vec<String>) -> Result<Vec<u8>, KernelError> {
    let out = match name {
        b"FlateDecode" | b"Fl" => {
            let (d, note) = inflate::zlib(data, limit)?;
            if let Some(n) = note {
                notes.push(n);
            }
            d
        }
        b"LZWDecode" | b"LZW" => {
            let early = parms.and_then(|p| p.get(b"EarlyChange")).and_then(Obj::as_i64).unwrap_or(1);
            let (d, damaged) = lzw(data, early != 0, limit)?;
            if damaged {
                notes.push(format!("an LZW stream in the PDF was damaged; the {} bytes before the damage were used", d.len()));
            }
            d
        }
        b"ASCIIHexDecode" | b"AHx" => hex(data)?,
        b"ASCII85Decode" | b"A85" => ascii85(data)?,
        b"RunLengthDecode" | b"RL" => run_length(data, limit)?,
        b"Crypt" => data.to_vec(),
        b"DCTDecode" | b"DCT" | b"JPXDecode" | b"CCITTFaxDecode" | b"CCF" | b"JBIG2Decode" => {
            return Err(KernelError::Unsupported(format!(
                "the PDF uses the image filter {} on a drawing stream", String::from_utf8_lossy(name))))
        }
        other => return Err(KernelError::Unsupported(format!("the PDF uses the unknown filter '{}'", String::from_utf8_lossy(other)))),
    };
    if out.len() > limit {
        return Err(bad(format!("a stream in the PDF expands past {} MB; the file was refused", limit >> 20)));
    }
    match name {
        b"FlateDecode" | b"Fl" | b"LZWDecode" | b"LZW" => predict(parms, out),
        _ => Ok(out),
    }
}

// ---- LZW -----------------------------------------------------------------------------------------------------------

/// Returns the data and whether it stopped at damage instead of an EOD code or the end of input.
pub fn lzw(input: &[u8], early: bool, limit: usize) -> Result<(Vec<u8>, bool), KernelError> {
    let mut out: Vec<u8> = Vec::new();
    // Entry i (258 and up) is out[start .. start + len]: the previous string plus the next string's first byte,
    // which in the output always immediately follows it.
    let mut table: Vec<(usize, usize)> = Vec::with_capacity(4096);
    let (mut buf, mut count, mut pos) = (0u32, 0u32, 0usize);
    let mut prev: Option<(usize, usize)> = None;
    let early = early as usize;
    loop {
        let next_code = 258 + table.len();
        let width = if next_code + 1 + early >= 2048 { 12 } else if next_code + 1 + early >= 1024 { 11 }
            else if next_code + 1 + early >= 512 { 10 } else { 9 };
        while count < width {
            let Some(&b) = input.get(pos) else { return Ok((out, false)) };
            pos += 1;
            buf = (buf << 8) | b as u32;
            count += 8;
        }
        let code = ((buf >> (count - width)) & ((1 << width) - 1)) as usize;
        count -= width;
        buf &= (1u32 << count).wrapping_sub(1);
        if code == 256 {
            table.clear();
            prev = None;
            continue;
        }
        if code == 257 {
            return Ok((out, false));
        }
        let start = out.len();
        if code < 256 {
            out.push(code as u8);
        } else if code - 258 < table.len() {
            let (s, l) = table[code - 258];
            if out.len() + l > limit {
                return Err(bad(format!("an LZW stream in the PDF expands past {} MB", limit >> 20)));
            }
            for k in 0..l {
                let b = out[s + k];
                out.push(b);
            }
        } else if code - 258 == table.len() && prev.is_some() {
            // The code being defined right now: the previous string plus its own first byte.
            let (s, l) = prev.unwrap_or((0, 0));
            if out.len() + l + 1 > limit {
                return Err(bad(format!("an LZW stream in the PDF expands past {} MB", limit >> 20)));
            }
            for k in 0..l {
                let b = out[s + k];
                out.push(b);
            }
            let first = out[s];
            out.push(first);
        } else {
            return Ok((out, true));
        }
        if out.len() > limit {
            return Err(bad(format!("an LZW stream in the PDF expands past {} MB", limit >> 20)));
        }
        if let Some((s, l)) = prev {
            if table.len() < 4096 - 258 {
                table.push((s, l + 1));
            }
        }
        prev = Some((start, out.len() - start));
    }
}

// ---- ASCII filters and run length ----------------------------------------------------------------------------------

fn is_white(c: u8) -> bool {
    matches!(c, 0 | 9 | 10 | 12 | 13 | 32)
}

pub fn hex(data: &[u8]) -> Result<Vec<u8>, KernelError> {
    let mut out = Vec::with_capacity(data.len() / 2);
    let mut high: Option<u8> = None;
    for &c in data {
        if c == b'>' {
            break;
        }
        if is_white(c) {
            continue;
        }
        let v = match c {
            b'0'..=b'9' => c - b'0',
            b'a'..=b'f' => c - b'a' + 10,
            b'A'..=b'F' => c - b'A' + 10,
            _ => return Err(bad(format!("an ASCIIHex stream in the PDF has the character {:?}", c as char))),
        };
        match high.take() {
            Some(h) => out.push(h << 4 | v),
            None => high = Some(v),
        }
    }
    if let Some(h) = high {
        out.push(h << 4);
    }
    Ok(out)
}

pub fn ascii85(data: &[u8]) -> Result<Vec<u8>, KernelError> {
    let mut out = Vec::with_capacity(data.len() * 4 / 5 + 4);
    let mut group = [0u8; 5];
    let mut n = 0;
    let body = data.strip_prefix(b"<~").unwrap_or(data);
    for &c in body {
        match c {
            b'~' => break,
            b'z' if n == 0 => out.extend_from_slice(&[0, 0, 0, 0]),
            b'!'..=b'u' => {
                group[n] = c - b'!';
                n += 1;
                if n == 5 {
                    let v = group.iter().fold(0u64, |acc, &d| acc * 85 + d as u64);
                    if v > u32::MAX as u64 {
                        return Err(bad("an ASCII85 stream in the PDF has a group out of range".into()));
                    }
                    out.extend_from_slice(&(v as u32).to_be_bytes());
                    n = 0;
                }
            }
            c if is_white(c) => {}
            _ => return Err(bad(format!("an ASCII85 stream in the PDF has the character {:?}", c as char))),
        }
    }
    if n == 1 {
        return Err(bad("an ASCII85 stream in the PDF ends with a single character".into()));
    }
    if n > 1 {
        for g in group.iter_mut().skip(n) {
            *g = 84;
        }
        let v = group.iter().fold(0u64, |acc, &d| acc * 85 + d as u64).min(u32::MAX as u64);
        out.extend_from_slice(&(v as u32).to_be_bytes()[..n - 1]);
    }
    Ok(out)
}

pub fn run_length(data: &[u8], limit: usize) -> Result<Vec<u8>, KernelError> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let l = data[i] as usize;
        i += 1;
        if l == 128 {
            break;
        }
        if l < 128 {
            let end = (i + l + 1).min(data.len());
            out.extend_from_slice(&data[i..end]);
            i = end;
        } else {
            let Some(&b) = data.get(i) else { break };
            out.extend(std::iter::repeat_n(b, 257 - l));
            i += 1;
        }
        if out.len() > limit {
            return Err(bad(format!("a run-length stream in the PDF expands past {} MB", limit >> 20)));
        }
    }
    Ok(out)
}

// ---- Predictors ----------------------------------------------------------------------------------------------------

fn predict(parms: Option<&Dict>, data: Vec<u8>) -> Result<Vec<u8>, KernelError> {
    let Some(p) = parms else { return Ok(data) };
    let get = |k: &[u8], d: i64| p.get(k).and_then(Obj::as_i64).unwrap_or(d);
    let predictor = get(b"Predictor", 1);
    if predictor <= 1 {
        return Ok(data);
    }
    let colors = get(b"Colors", 1);
    let bpc = get(b"BitsPerComponent", 8);
    let columns = get(b"Columns", 1);
    if !(1..=64).contains(&colors) || ![1, 2, 4, 8, 16].contains(&bpc) || !(1..=1 << 24).contains(&columns) {
        return Err(bad(format!("a stream in the PDF has bad predictor parameters (Colors {colors}, BitsPerComponent {bpc}, Columns {columns})")));
    }
    let bits_per_pixel = (colors * bpc) as usize;
    let row = (bits_per_pixel * columns as usize).div_ceil(8);
    let bpp = bits_per_pixel.div_ceil(8).max(1);
    if predictor == 2 {
        return Ok(tiff(data, row, colors as usize, bpc as usize));
    }
    if predictor < 10 {
        return Err(KernelError::Unsupported(format!("the PDF uses the unknown predictor {predictor}")));
    }
    let mut out = Vec::with_capacity(data.len());
    let mut prior = vec![0u8; row];
    for chunk in data.chunks(row + 1) {
        let tag = chunk[0];
        let mut cur = chunk[1..].to_vec();
        cur.resize(row, 0);
        for i in 0..row {
            let a = if i >= bpp { cur[i - bpp] } else { 0 };
            let b = prior[i];
            let c = if i >= bpp { prior[i - bpp] } else { 0 };
            let add = match tag {
                0 => 0,
                1 => a,
                2 => b,
                3 => ((a as u16 + b as u16) / 2) as u8,
                4 => {
                    let p = a as i16 + b as i16 - c as i16;
                    let (pa, pb, pc) = ((p - a as i16).abs(), (p - b as i16).abs(), (p - c as i16).abs());
                    if pa <= pb && pa <= pc { a } else if pb <= pc { b } else { c }
                }
                t => return Err(bad(format!("a stream in the PDF has the unknown PNG row type {t}"))),
            };
            cur[i] = cur[i].wrapping_add(add);
        }
        out.extend_from_slice(&cur[..chunk.len() - 1]);
        prior = cur;
    }
    Ok(out)
}

/// TIFF Predictor 2: each component is the difference from the same component of the sample to its left.
fn tiff(mut data: Vec<u8>, row: usize, colors: usize, bpc: usize) -> Vec<u8> {
    if row == 0 {
        return data;
    }
    for line in data.chunks_mut(row) {
        match bpc {
            8 => {
                for i in colors..line.len() {
                    line[i] = line[i].wrapping_add(line[i - colors]);
                }
            }
            16 => {
                let n = line.len() / 2;
                for i in colors..n {
                    let v = u16::from_be_bytes([line[2 * i], line[2 * i + 1]])
                        .wrapping_add(u16::from_be_bytes([line[2 * (i - colors)], line[2 * (i - colors) + 1]]));
                    line[2 * i..2 * i + 2].copy_from_slice(&v.to_be_bytes());
                }
            }
            _ => {
                let per_byte = 8 / bpc;
                let n = line.len() * per_byte;
                let mask = (1u16 << bpc) - 1;
                let get = |l: &[u8], i: usize| ((l[i / per_byte] as u16) >> (8 - bpc * (i % per_byte + 1))) & mask;
                for i in colors..n {
                    let v = (get(line, i) + get(line, i - colors)) & mask;
                    let shift = 8 - bpc * (i % per_byte + 1);
                    let byte = &mut line[i / per_byte];
                    *byte = (*byte & !((mask as u8) << shift)) | ((v as u8) << shift);
                }
            }
        }
    }
    data
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lzw_example_from_the_specification() {
        let (d, damaged) = lzw(&[0x80, 0x0B, 0x60, 0x50, 0x22, 0x0C, 0x0C, 0x85, 0x01], true, 1 << 20).unwrap();
        assert_eq!(d, [45, 45, 45, 45, 45, 65, 45, 45, 45, 66]);
        assert!(!damaged);
    }

    #[test]
    fn ascii85_and_hex() {
        assert_eq!(ascii85(b"87cURD]j7BEbo80~>").unwrap(), b"Hello world!");
        assert_eq!(ascii85(b"z~>").unwrap(), [0, 0, 0, 0]);
        assert_eq!(hex(b"48 65 6C6c 7>").unwrap(), b"Hellp");
        assert_eq!(run_length(&[2, b'a', b'b', b'c', 254, b'x', 128], 100).unwrap(), b"abcxxx");
    }
}
