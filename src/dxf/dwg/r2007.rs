//! The R2007 file organisation (ODA specification chapter 5): a Reed-Solomon protected file header at 0x80, system
//! pages (page map, section map) and data pages that are Reed-Solomon coded (interleaved or not) and compressed with
//! the R21 LZ77 variant. The files read here are taken as undamaged: the parity bytes are set aside, not used to
//! correct errors.

use std::collections::HashMap;

use super::bits::{bad, R};
use super::r2004::Sections;

/// Sections larger than this are refused rather than allocated (a damaged size field).
const MAX_SECTION: u64 = 512 << 20;

fn le64(d: &[u8], at: usize) -> R<u64> {
    match d.get(at..at.saturating_add(8)) {
        Some(b) if b.len() == 8 => Ok(u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])),
        _ => bad(format!("R2007 data ends early at byte {at}")),
    }
}

/// The data bytes of `blocks` Reed-Solomon code words of `k` data bytes each (n = 255): interleaved code words put
/// byte i of word j at i * blocks + j; otherwise the data blocks come first, then all the parity.
fn rs_data(raw: &[u8], blocks: usize, k: usize, interleaved: bool) -> R<Vec<u8>> {
    if raw.len() < blocks.saturating_mul(255) {
        return bad("a Reed-Solomon coded page runs past the end of the file");
    }
    let mut out = Vec::with_capacity(blocks * k);
    if interleaved {
        for j in 0..blocks {
            for i in 0..k {
                out.push(raw[i * blocks + j]);
            }
        }
    } else {
        out.extend_from_slice(&raw[..blocks * k]);
    }
    Ok(out)
}

pub(super) fn sections(file: &[u8], wanted: &[&str]) -> R<Sections> {
    // File header: 3 interleaved (255, 239) code words at 0x80.
    let raw = file.get(0x80..0x80 + 3 * 255).ok_or_else(|| short())?;
    let decoded = rs_data(raw, 3, 239, true)?;
    let length = i64::from(i32::from_le_bytes([decoded[0x18], decoded[0x19], decoded[0x1a], decoded[0x1b]]));
    let header = if length < 0 {
        decoded.get(0x20..0x20 + (-length) as usize).ok_or_else(|| short())?.to_vec()
    } else {
        let body = decoded.get(0x20..0x20 + length as usize).ok_or_else(|| short())?;
        decompress(body, 0x110)?
    };
    let field = |at: usize| le64(&header, at);
    let pages_map_offset = field(0x38)?;
    let pages_map_compressed = field(0x50)?;
    let pages_map_size = field(0x58)?;
    let pages_map_factor = field(0x18)?;
    let sections_map_compressed = field(0xB0)?;
    let sections_map_id = field(0xC0)?;
    let sections_map_size = field(0xC8)?;
    let sections_map_factor = field(0xD8)?;

    let page_map = system_page(file, pages_map_offset.saturating_add(0x480), pages_map_compressed, pages_map_size,
        pages_map_factor)?;
    let mut pages: HashMap<i64, (u64, u64)> = HashMap::new(); // id -> (stream position, size)
    let mut offset = 0u64;
    let mut at = 0;
    while at + 16 <= page_map.len() {
        let size = le64(&page_map, at)?;
        let id = le64(&page_map, at + 8)? as i64;
        at += 16;
        pages.insert(id.saturating_abs(), (offset.saturating_add(0x480), size));
        offset = offset.saturating_add(size);
    }
    let Some(&(map_at, _)) = pages.get(&(sections_map_id as i64)) else {
        return bad("the R2007 section map page is missing");
    };
    let map = system_page(file, map_at, sections_map_compressed, sections_map_size, sections_map_factor)?;

    let mut by_name = HashMap::new();
    let mut at = 0usize;
    while at + 0x40 <= map.len() {
        let size = le64(&map, at)?;
        let encrypted = le64(&map, at + 0x10)?;
        // The name length is in bytes, its two terminating zero bytes included.
        let name_len = le64(&map, at + 0x20)? as usize;
        let encoding = le64(&map, at + 0x30)?;
        let page_count = le64(&map, at + 0x38)? as usize;
        at += 0x40;
        if name_len > 256 || page_count > map.len() / 56 + 1 {
            return bad("the R2007 section map is damaged");
        }
        let units: Vec<u16> = (0..name_len / 2)
            .map(|i| map.get(at + 2 * i..at + 2 * i + 2).map(|b| u16::from_le_bytes([b[0], b[1]])))
            .collect::<Option<_>>()
            .ok_or_else(|| short())?;
        let name = String::from_utf16_lossy(&units).trim_end_matches('\0').to_string();
        at += name_len;
        let mut list = Vec::with_capacity(page_count);
        for _ in 0..page_count {
            // offset, size, id, uncompressed size, compressed size, checksum, CRC
            list.push((le64(&map, at)?, le64(&map, at + 16)? as i64, le64(&map, at + 24)?, le64(&map, at + 32)?));
            at += 56;
        }
        if !wanted.iter().any(|w| w.eq_ignore_ascii_case(&name)) {
            continue;
        }
        if encrypted == 1 {
            return Err(crate::KernelError::Unsupported(
                "this DWG is password protected (encrypted); save it without a password".into()));
        }
        if size > MAX_SECTION || size > (file.len() as u64).saturating_mul(64) + 0x10000 {
            return bad(format!("section '{name}' claims {size} bytes"));
        }
        let mut data = vec![0u8; size as usize];
        for (start, id, full, compressed) in list {
            let Some(&(page_at, _)) = pages.get(&id) else {
                return bad(format!("section '{name}' uses page {id}, which is not in the page map"));
            };
            if start >= size || full > 0x100000 || compressed > 0x100000 {
                continue;
            }
            let aligned = (compressed as usize + 7) & !7;
            let blocks = aligned.div_ceil(251);
            let at = usize::try_from(page_at).map_err(|_| short())?;
            let raw = file.get(at..at.saturating_add(blocks * 255)).ok_or_else(|| short())?;
            let stored = rs_data(raw, blocks, 251, encoding == 4)?;
            let stored = &stored[..compressed as usize];
            let page = if compressed < full { decompress(stored, full as usize)? } else { stored.to_vec() };
            let take = ((size - start) as usize).min(page.len());
            let start = start as usize;
            data[start..start + take].copy_from_slice(&page[..take]);
        }
        by_name.insert(name, data);
    }
    Ok(Sections { by_name })
}

fn short() -> crate::KernelError {
    crate::KernelError::InvalidInput("DWG: the file ends early".into())
}

/// A system page: the (padded) data repeated `factor` times, (255, 239) Reed-Solomon coded and interleaved.
fn system_page(file: &[u8], address: u64, compressed: u64, size: u64, factor: u64) -> R<Vec<u8>> {
    if compressed > 0x400000 || size > 0x400000 || factor == 0 || factor > 0x10000 {
        return bad("an R2007 system page has impossible sizes");
    }
    let aligned = (compressed as usize + 7) & !7;
    let blocks = (aligned * factor as usize).div_ceil(239);
    let at = usize::try_from(address).map_err(|_| short())?;
    let raw = file.get(at..at.saturating_add(blocks * 255)).ok_or_else(|| short())?;
    let data = rs_data(raw, blocks, 239, true)?;
    let stored = data.get(..compressed as usize).ok_or_else(|| short())?;
    if compressed < size {
        decompress(stored, size as usize)
    } else {
        Ok(stored[..size as usize].to_vec())
    }
}

/// The order in which a literal run of 1-32 bytes lands in the output: (block size, source index) pairs, from the
/// specification's table (5.10.1).
const LITERAL_ORDER: [&[(usize, usize)]; 33] = [
    &[],
    &[(1, 0)],
    &[(1, 1), (1, 0)],
    &[(1, 2), (1, 1), (1, 0)],
    &[(1, 0), (1, 1), (1, 2), (1, 3)],
    &[(1, 4), (4, 0)],
    &[(1, 5), (4, 1), (1, 0)],
    &[(2, 5), (4, 1), (1, 0)],
    &[(4, 0), (4, 4)],
    &[(1, 8), (8, 0)],
    &[(1, 9), (8, 1), (1, 0)],
    &[(2, 9), (8, 1), (1, 0)],
    &[(4, 8), (8, 0)],
    &[(1, 12), (4, 8), (8, 0)],
    &[(1, 13), (4, 9), (8, 1), (1, 0)],
    &[(2, 13), (4, 9), (8, 1), (1, 0)],
    &[(8, 8), (8, 0)],
    &[(8, 9), (1, 8), (8, 0)],
    &[(1, 17), (16, 1), (1, 0)],
    &[(3, 16), (16, 0)],
    &[(4, 16), (16, 0)],
    &[(1, 20), (4, 16), (16, 0)],
    &[(2, 20), (4, 16), (16, 0)],
    &[(3, 20), (4, 16), (16, 0)],
    &[(8, 16), (16, 0)],
    &[(8, 17), (1, 16), (16, 0)],
    &[(1, 25), (8, 17), (1, 16), (16, 0)],
    &[(2, 25), (8, 17), (1, 16), (16, 0)],
    &[(4, 24), (8, 16), (16, 0)],
    &[(1, 28), (4, 24), (8, 16), (16, 0)],
    &[(2, 28), (4, 24), (8, 16), (16, 0)],
    &[(1, 30), (4, 26), (8, 18), (16, 2), (2, 0)],
    &[(16, 16), (16, 0)],
];

/// A literal piece of 1-32 bytes in output order: each block of the table is itself laid out by the table for its
/// size (so a 3-byte block is reversed and a 16-byte block has its halves swapped), down to single bytes.
fn shuffle(piece: &[u8], out: &mut Vec<u8>) {
    if piece.len() == 1 {
        out.push(piece[0]);
        return;
    }
    for &(size, from) in LITERAL_ORDER[piece.len()] {
        shuffle(&piece[from..from + size], out);
    }
}

fn literal(src: &[u8], i: &mut usize, n: usize, out: &mut Vec<u8>, limit: usize) -> R<()> {
    let chunk = src.get(*i..i.saturating_add(n)).ok_or_else(|| crate::KernelError::InvalidInput(
        "DWG: R2007 compressed literal runs past the data".into()))?;
    if out.len() + n > limit {
        return bad("R2007 compressed data expands past its page");
    }
    let mut rest = chunk;
    while !rest.is_empty() {
        let take = rest.len().min(32);
        let (piece, tail) = rest.split_at(take);
        shuffle(piece, out);
        rest = tail;
    }
    *i += n;
    Ok(())
}

fn byte(src: &[u8], i: &mut usize) -> R<usize> {
    let b = *src.get(*i).ok_or_else(|| crate::KernelError::InvalidInput("DWG: R2007 compressed data ends early".into()))?;
    *i += 1;
    Ok(b as usize)
}

/// One back-reference instruction (spec 5.10.2): gives (offset back, length); `op` becomes the opcode whose low
/// three bits count the literals that follow.
fn instruction(src: &[u8], i: &mut usize, op: &mut u8) -> R<(usize, usize)> {
    let o = *op as usize;
    let (offset, length);
    match o >> 4 {
        0 => {
            let mut len = (o & 0xf) + 0x13;
            let low = byte(src, i)?;
            let next = byte(src, i)?;
            *op = next as u8;
            len += (next >> 3) & 0x10;
            offset = ((next & 0x78) << 5) + 1 + low;
            length = len;
        }
        1 => {
            let len = (o & 0xf) + 3;
            let low = byte(src, i)?;
            let next = byte(src, i)?;
            *op = next as u8;
            offset = ((next & 0xf8) << 5) + 1 + low;
            length = len;
        }
        2 => {
            let mut off = byte(src, i)?;
            off |= (byte(src, i)? << 8) & 0xff00;
            let mut len = o & 7;
            if o & 8 == 0 {
                let next = byte(src, i)?;
                *op = next as u8;
                len += next & 0xf8;
            } else {
                off += 1;
                len += byte(src, i)? << 3;
                let next = byte(src, i)?;
                *op = next as u8;
                len += ((next & 0xf8) << 8) + 0x100;
            }
            offset = off;
            length = len;
        }
        _ => {
            let len = o >> 4;
            let low = o & 15;
            let next = byte(src, i)?;
            *op = next as u8;
            offset = ((next & 0xf8) << 1) + low + 1;
            length = len;
        }
    }
    Ok((offset, length))
}

/// R21 LZ77 variant (spec 5.10): literal runs (stored in a shuffled byte order) alternating with back references.
pub(super) fn decompress(src: &[u8], limit: usize) -> R<Vec<u8>> {
    let mut out = Vec::with_capacity(limit.min(1 << 20));
    let mut i = 0usize;
    let mut length = 0usize;
    let mut pending: Option<u8> = None;
    if src.first().is_some_and(|b| b >> 4 == 2) {
        // Past the opcode and two more bytes; the next byte gives the first literal count.
        i = 4;
        length = *src.get(3).ok_or_else(short)? as usize & 7;
    }
    let mut guard = 0usize;
    while i < src.len() || pending.is_some() {
        guard += 1;
        if guard > 1 << 26 {
            return bad("R2007 compressed data does not end");
        }
        if length == 0 {
            let op = match pending.take() {
                Some(op) => op as usize,
                None => byte(src, &mut i)?,
            };
            length = op + 8;
            if length == 0x17 {
                let n = byte(src, &mut i)?;
                length += n;
                if n == 0xff {
                    loop {
                        let n = byte(src, &mut i)? | byte(src, &mut i)? << 8;
                        length += n;
                        if n != 0xffff {
                            break;
                        }
                        if length > limit {
                            return bad("R2007 literal run too long");
                        }
                    }
                }
            }
        }
        literal(src, &mut i, length, &mut out, limit)?;
        if i >= src.len() {
            break;
        }
        // Back references until one carries a literal count or a literal opcode comes up.
        let mut op = byte(src, &mut i)? as u8;
        let (mut offset, mut count) = instruction(src, &mut i, &mut op)?;
        loop {
            if offset > out.len() || offset == 0 {
                return bad("R2007 compressed data refers back before its start");
            }
            if out.len() + count > limit {
                return bad("R2007 compressed data expands past its page");
            }
            let from = out.len() - offset;
            for k in 0..count {
                let b = out[from + k];
                out.push(b);
            }
            length = (op & 7) as usize;
            if length != 0 || i >= src.len() {
                break;
            }
            op = byte(src, &mut i)? as u8;
            if op >> 4 == 0 {
                pending = Some(op);
                break;
            }
            if op >> 4 == 15 {
                op &= 15;
            }
            (offset, count) = instruction(src, &mut i, &mut op)?;
        }
    }
    Ok(out)
}
