//! The R2004 file organisation, also used by R2010, R2013 and R2018 (ODA specification chapter 4, and chapters 6-8
//! for the later versions): the masked file header at 0x80, the system section page map and section map, masked
//! data page headers and the R2004 LZ77 variant; gives each named section as one contiguous byte buffer.

use std::collections::HashMap;

use super::bits::{bad, le32, le64, R};

/// Sections larger than this are refused rather than allocated (a damaged size field).
const MAX_SECTION: u64 = 512 << 20;

pub(super) struct Sections {
    pub by_name: HashMap<String, Vec<u8>>,
}

/// The 0x6c-byte mask for the file header data (spec 4.1: a linear congruential sequence seeded with 1).
fn header_mask() -> [u8; 0x6c] {
    let mut out = [0u8; 0x6c];
    let mut seed: u32 = 1;
    for b in out.iter_mut() {
        seed = seed.wrapping_mul(0x343fd).wrapping_add(0x269ec3);
        *b = (seed >> 16) as u8;
    }
    out
}

struct Description {
    name: String,
    size: u64,
    max_page: u64,
    compressed: u32,
    encrypted: u32,
    pages: Vec<(i32, u64)>, // (page number, start offset)
}

pub(super) fn sections(file: &[u8], wanted: &[&str]) -> R<Sections> {
    let raw = file.get(0x80..0x80 + 0x6c).ok_or_else(|| short())?;
    let mask = header_mask();
    let mut h = [0u8; 0x6c];
    for i in 0..0x6c {
        h[i] = raw[i] ^ mask[i];
    }
    if &h[0..11] != b"AcFssFcAJMB" {
        return bad("the R2004 file header is not where it should be (damaged file)");
    }
    let page_map_address = le64(&h, 0x54)?.checked_add(0x100).ok_or_else(|| short())?;
    let section_map_id = le32(&h, 0x5c)? as i32;

    // Section page map: the address of every page, numbered from 1.
    let page_map = system_page(file, page_map_address, 0x41630e3b)?;
    let mut addresses: HashMap<i32, u64> = HashMap::new();
    let mut address: u64 = 0x100;
    let mut at = 0usize;
    while at + 8 <= page_map.len() {
        let number = le32(&page_map, at)? as i32;
        let size = le32(&page_map, at + 4)? as u64;
        at += 8;
        if number < 0 {
            at += 16;
        } else {
            addresses.insert(number, address);
        }
        address = address.saturating_add(size);
    }
    let section_map_address = *addresses.get(&section_map_id)
        .ok_or_else(|| super::bits::bad::<()>("the section map page is missing").unwrap_err())?;
    let map = system_page(file, section_map_address, 0x4163003b)?;

    let count = le32(&map, 0)? as usize;
    let mut descriptions = Vec::new();
    let mut at = 0x14usize;
    for _ in 0..count {
        if at + 0x60 > map.len() {
            return bad("the section map ends early");
        }
        let size = le64(&map, at)?;
        let page_count = le32(&map, at + 8)? as usize;
        let max_page = le32(&map, at + 0x0c)? as u64;
        let compressed = le32(&map, at + 0x14)?;
        let encrypted = le32(&map, at + 0x1c)?;
        let name_bytes = &map[at + 0x20..at + 0x60];
        let name_end = name_bytes.iter().position(|&b| b == 0).unwrap_or(name_bytes.len());
        let name = String::from_utf8_lossy(&name_bytes[..name_end]).to_string();
        at += 0x60;
        if page_count.saturating_mul(16) > map.len().saturating_sub(at) {
            return bad(format!("section '{name}' lists more pages than the section map holds"));
        }
        let mut pages = Vec::with_capacity(page_count);
        for _ in 0..page_count {
            pages.push((le32(&map, at)? as i32, le64(&map, at + 8)?));
            at += 16;
        }
        descriptions.push(Description { name, size, max_page, compressed, encrypted, pages });
    }

    let mut by_name = HashMap::new();
    for d in descriptions {
        if !wanted.iter().any(|w| w.eq_ignore_ascii_case(&d.name)) || d.name.is_empty() {
            continue;
        }
        if d.encrypted == 1 {
            return Err(crate::KernelError::Unsupported(
                "this DWG is password protected (encrypted); save it without a password".into()));
        }
        if d.size > MAX_SECTION || d.size > (file.len() as u64).saturating_mul(64) + 0x10000 {
            return bad(format!("section '{}' claims {} bytes", d.name, d.size));
        }
        if d.max_page == 0 || d.max_page > 0x100000 {
            return bad(format!("section '{}' has page size {}", d.name, d.max_page));
        }
        let mut data = vec![0u8; d.size as usize];
        for &(number, start) in &d.pages {
            let Some(&address) = addresses.get(&number) else {
                return bad(format!("section '{}' uses page {number}, which is not in the page map", d.name));
            };
            if start >= d.size {
                continue;
            }
            let page = data_page(file, address, d.compressed == 2, d.max_page as usize)?;
            let take = (d.size - start).min(d.max_page).min(page.len() as u64) as usize;
            let start = start as usize;
            data[start..start + take].copy_from_slice(&page[..take]);
        }
        by_name.insert(d.name, data);
    }
    Ok(Sections { by_name })
}

fn short() -> crate::KernelError {
    crate::KernelError::InvalidInput("DWG: the file ends early".into())
}

/// A system page (page map or section map): 0x14-byte header, then compressed data.
fn system_page(file: &[u8], address: u64, kind: u32) -> R<Vec<u8>> {
    let at = usize::try_from(address).ok().filter(|&a| a < file.len()).ok_or_else(short)?;
    if le32(file, at)? != kind {
        return bad(format!("expected a system page at byte {at}"));
    }
    let size = le32(file, at + 4)? as usize;
    let compressed = le32(file, at + 8)? as usize;
    let method = le32(file, at + 12)?;
    if size > 64 << 20 {
        return bad("a system page is too large");
    }
    let body = file.get(at + 0x14..(at + 0x14).saturating_add(compressed)).ok_or_else(short)?;
    if method == 2 {
        let out = decompress(body, size)?;
        if out.len() != size {
            return bad(format!("a system page decompressed to {} bytes, not {size}", out.len()));
        }
        Ok(out)
    } else {
        Ok(body.to_vec())
    }
}

/// A data page: 32-byte header masked with 0x4164536b ^ its address, then the stored payload.
fn data_page(file: &[u8], address: u64, compressed: bool, max: usize) -> R<Vec<u8>> {
    let at = usize::try_from(address).ok().filter(|&a| a < file.len()).ok_or_else(short)?;
    let mask = 0x4164536bu32 ^ (address as u32);
    let mut header = [0u32; 8];
    for (i, word) in header.iter_mut().enumerate() {
        *word = le32(file, at + 4 * i)? ^ mask;
    }
    if header[0] != 0x4163043b {
        return bad(format!("expected a data page at byte {at}"));
    }
    let stored = header[2] as usize;
    let body = file.get(at + 32..(at + 32).saturating_add(stored)).ok_or_else(short)?;
    if compressed {
        decompress(body, max)
    } else {
        Ok(body.to_vec())
    }
}

/// R2004 LZ77 variant (spec 4.7). Output is limited to `limit` bytes.
pub(super) fn decompress(input: &[u8], limit: usize) -> R<Vec<u8>> {
    let mut out: Vec<u8> = Vec::with_capacity(limit.min(1 << 20));
    let mut i = 0usize;
    let byte = |i: &mut usize| -> R<u8> {
        let b = *input.get(*i).ok_or_else(|| crate::KernelError::InvalidInput(
            "DWG: compressed data ends early".into()))?;
        *i += 1;
        Ok(b)
    };
    // Literal length: 0x01-0x0F add 3; 0x00 starts a running total; a byte with its high nibble set is the next
    // opcode and means no literals.
    fn literal(input: &[u8], i: &mut usize, pending: &mut Option<u8>) -> R<usize> {
        let Some(&first) = input.get(*i) else { return bad("compressed data ends early") };
        *i += 1;
        if first & 0xf0 != 0 {
            *pending = Some(first);
            return Ok(0);
        }
        if first != 0 {
            return Ok(first as usize + 3);
        }
        let mut total = 0x0fusize;
        loop {
            let Some(&b) = input.get(*i) else { return bad("compressed data ends early") };
            *i += 1;
            if b == 0 {
                total += 0xff;
                if total > 1 << 28 {
                    return bad("literal run too long");
                }
            } else {
                return Ok(total + b as usize + 3);
            }
        }
    }
    fn long_count(input: &[u8], i: &mut usize) -> R<usize> {
        let Some(&first) = input.get(*i) else { return bad("compressed data ends early") };
        *i += 1;
        if first != 0 {
            return Ok(first as usize);
        }
        let mut total = 0xffusize;
        loop {
            let Some(&b) = input.get(*i) else { return bad("compressed data ends early") };
            *i += 1;
            if b == 0 {
                total += 0xff;
                if total > 1 << 28 {
                    return bad("copy run too long");
                }
            } else {
                return Ok(total + b as usize);
            }
        }
    }
    let mut pending: Option<u8> = None;
    let copy_literal = |n: usize, i: &mut usize, out: &mut Vec<u8>| -> R<()> {
        let src = input.get(*i..i.saturating_add(n)).ok_or_else(|| crate::KernelError::InvalidInput(
            "DWG: compressed literal runs past the data".into()))?;
        if out.len() + n > limit {
            return bad("compressed data expands past its page");
        }
        out.extend_from_slice(src);
        *i += n;
        Ok(())
    };
    let n = literal(input, &mut i, &mut pending)?;
    copy_literal(n, &mut i, &mut out)?;
    loop {
        let op = match pending.take() {
            Some(op) => op,
            None => {
                if i >= input.len() {
                    break;
                }
                byte(&mut i)?
            }
        };
        let (count, offset, mut lit);
        match op {
            0x11 => break,
            0x10 | 0x12..=0x1f | 0x20..=0x3f => {
                count = match op {
                    0x10 => long_count(input, &mut i)? + 9,
                    0x12..=0x1f => (op & 0x0f) as usize + 2,
                    0x20 => long_count(input, &mut i)? + 0x21,
                    _ => op as usize - 0x1e,
                };
                let first = byte(&mut i)?;
                let second = byte(&mut i)?;
                let base = (first >> 2) as usize | (second as usize) << 6;
                offset = if op < 0x20 { base + 0x3fff } else { base };
                lit = (first & 3) as usize;
            }
            0x40..=0xff => {
                count = ((op & 0xf0) >> 4) as usize - 1;
                let second = byte(&mut i)?;
                offset = (second as usize) << 2 | ((op & 0x0c) >> 2) as usize;
                lit = (op & 3) as usize;
            }
            _ => return bad(format!("bad compression opcode {op:#x}")),
        }
        if offset + 1 > out.len() {
            return bad("compressed data refers back before its start");
        }
        if out.len() + count > limit {
            return bad("compressed data expands past its page");
        }
        let from = out.len() - offset - 1;
        for k in 0..count {
            let b = out[from + k];
            out.push(b);
        }
        if lit == 0 {
            lit = literal(input, &mut i, &mut pending)?;
        }
        copy_literal(lit, &mut i, &mut out)?;
    }
    Ok(out)
}
