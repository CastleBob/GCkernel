//! DWG drawings in (Bob 2026-10-31: "make sure we get a DWG importer"). Written from the Open Design Alliance's
//! published "Open Design Specification for .dwg files" only; no third-party DWG source.
//!
//! The reader decodes the file's objects into the same entity records (DXF group codes) the DXF reader builds and
//! hands them to the DXF reader's assembler, so blocks, object coordinate systems, layer inheritance, paper space,
//! splines and polylines behave identically for both formats.

mod bits;
mod header;
mod objects;
mod r15;
mod r2004;
mod r2007;

use std::collections::HashMap;

use bits::{bad, Ver, R};
use super::Drawing;
use crate::KernelError;

/// True when the bytes start like a DWG file ("AC10" followed by two digits).
pub fn looks_like(bytes: &[u8]) -> bool {
    bytes.len() >= 6 && &bytes[0..4] == b"AC10" && bytes[4].is_ascii_digit() && bytes[5].is_ascii_digit()
}

/// Reads a DWG file's bytes.
pub fn read(bytes: &[u8]) -> Result<Drawing, KernelError> {
    read_parts(bytes, false)
}

/// Only the header and layer table (for asking about units before an import).
pub fn read_header(bytes: &[u8]) -> Result<Drawing, KernelError> {
    read_parts(bytes, true)
}

fn version(bytes: &[u8]) -> R<(String, Ver)> {
    if !looks_like(bytes) {
        return bad("this is not a DWG file");
    }
    let name = String::from_utf8_lossy(&bytes[0..6]).to_string();
    let ver = match name.as_str() {
        "AC1012" => Ver::R13,
        "AC1014" => Ver::R14,
        "AC1015" => Ver::R2000,
        "AC1018" => Ver::R2004,
        "AC1021" => Ver::R2007,
        "AC1024" => Ver::R2010,
        "AC1027" => Ver::R2013,
        "AC1032" => Ver::R2018,
        v if v < "AC1012" => return Err(KernelError::Unsupported(format!(
            "DWG {v} (R12 or earlier) files are not supported; save as DWG 2000 or later, or DXF"))),
        v => return Err(KernelError::Unsupported(format!(
            "DWG {v} is newer than this reader knows (2018 format, AC1032); save as DWG 2018 or DXF"))),
    };
    Ok((name, ver))
}

/// The parts of a DWG file the reader uses, whatever the file organisation.
struct Source {
    ver: Ver,
    /// Header variables, from the first data byte after the size field(s).
    header: Vec<u8>,
    /// Class definitions, from the first data byte after the size field(s).
    classes: Vec<u8>,
    /// Object handle -> byte offset in `objects`.
    map: HashMap<u64, usize>,
    /// Handles in object map order.
    order: Vec<u64>,
    objects: Vec<u8>,
    measurement: Option<i32>,
}

const HEADER_SENTINEL: [u8; 16] = [0xCF, 0x7B, 0x1F, 0x23, 0xFD, 0xDE, 0x38, 0xA9, 0x5F, 0x7C, 0x68, 0xB8, 0x4E, 0x6D, 0x33, 0x5F];
const CLASSES_SENTINEL: [u8; 16] = [0x8D, 0xA1, 0xC4, 0xB8, 0xC4, 0xA9, 0xF8, 0xC5, 0xC0, 0xDC, 0xF4, 0x5F, 0xE7, 0xCF, 0xB6, 0x8A];

fn load(bytes: &[u8], ver: Ver) -> R<Source> {
    let maintenance = bytes.get(0x0b).copied().unwrap_or(0);
    if ver <= Ver::R2000 {
        let loc = r15::locators(bytes)?;
        let header = after_sentinel(bytes, loc.header, &HEADER_SENTINEL, "header", 0)?;
        let classes = after_sentinel(bytes, loc.classes, &CLASSES_SENTINEL, "classes", 0)?;
        let end = loc.handles.saturating_add(loc.handles_size).min(bytes.len());
        let (map, order) = object_map(&bytes[loc.handles..end], bytes.len())?;
        let measurement = loc.template.and_then(|at| measurement(bytes.get(at..)?));
        return Ok(Source { ver, header, classes, map, order, objects: bytes.to_vec(), measurement });
    }
    let wanted = ["AcDb:Header", "AcDb:Classes", "AcDb:Handles", "AcDb:AcDbObjects", "AcDb:Template"];
    let mut s = if ver == Ver::R2007 { r2007::sections(bytes, &wanted)? } else { r2004::sections(bytes, &wanted)? };
    let mut take = |name: &str| -> R<Vec<u8>> {
        s.by_name.remove(name).ok_or_else(|| bad::<()>(format!("the {name} section is missing")).unwrap_err())
    };
    // R2010 and R2013 (maintenance version above 3) and R2018 add a long after the size.
    let extra = if ver >= Ver::R2018 || (ver >= Ver::R2010 && maintenance > 3) { 4 } else { 0 };
    let header_section = take("AcDb:Header")?;
    let header = after_sentinel(&header_section, 0, &HEADER_SENTINEL, "header", extra)?;
    let classes_section = take("AcDb:Classes")?;
    let classes = after_sentinel(&classes_section, 0, &CLASSES_SENTINEL, "classes", extra)?;
    let objects = take("AcDb:AcDbObjects")?;
    let handles = take("AcDb:Handles")?;
    let (map, order) = object_map(&handles, objects.len())?;
    let measurement = take("AcDb:Template").ok().and_then(|t| measurement(&t));
    Ok(Source { ver, header, classes, map, order, objects, measurement })
}

/// MEASUREMENT from template data (spec chapter 22): a description length, its bytes, then the value (0 or 1).
fn measurement(t: &[u8]) -> Option<i32> {
    let n = u16::from_le_bytes([*t.first()?, *t.get(1)?]) as usize;
    let m = u16::from_le_bytes([*t.get(2 + n)?, *t.get(3 + n)?]);
    (m <= 1).then_some(m as i32)
}

/// The data after a section's 16-byte sentinel and its size long (and `extra` bytes more).
fn after_sentinel(data: &[u8], at: usize, sentinel: &[u8; 16], what: &str, extra: usize) -> R<Vec<u8>> {
    if data.get(at..at + 16) != Some(&sentinel[..]) {
        return bad(format!("the {what} section does not start with its sentinel"));
    }
    let size = bits::le32(data, at + 16)? as usize;
    let start = at + 20 + extra;
    let end = start.saturating_add(size).min(data.len());
    Ok(data.get(start..end).unwrap_or(&[]).to_vec())
}

/// The object map (spec chapter 23): big-endian sized chunks of (handle offset, location offset) modular chars.
fn object_map(data: &[u8], limit: usize) -> R<(HashMap<u64, usize>, Vec<u64>)> {
    let mut map = HashMap::new();
    let mut order = Vec::new();
    let mut at = 0usize;
    loop {
        let (Some(&hi), Some(&lo)) = (data.get(at), data.get(at + 1)) else { break };
        let size = u16::from_be_bytes([hi, lo]) as usize;
        if size <= 2 {
            break;
        }
        let end = (at + size).min(data.len());
        let mut i = at + 2;
        // Each chunk starts again from handle 0 and offset 0.
        let (mut handle, mut location) = (0u64, 0i64);
        while i < end {
            let dh = bits::mc(data, &mut i, false)?;
            let dl = bits::mc(data, &mut i, true)?;
            handle = handle.wrapping_add(dh as u64);
            location = location.wrapping_add(dl);
            if location >= 0 && (location as usize) < limit && !map.contains_key(&handle) {
                map.insert(handle, location as usize);
                order.push(handle);
            }
        }
        at = end + 2; // the chunk's CRC
    }
    if map.is_empty() {
        return bad("the object map is empty");
    }
    Ok((map, order))
}

fn read_parts(bytes: &[u8], header_only: bool) -> R<Drawing> {
    let (name, ver) = version(bytes)?;
    let source = load(bytes, ver)?;
    let mut drawing = Drawing { version: name, ..Drawing::default() };
    drawing.insunits = header::insunits(&source);
    drawing.measurement = source.measurement;
    let mut reader = objects::Reader::new(&source)?;
    reader.layers(&mut drawing)?;
    if header_only {
        drawing.layer_index("0");
        return Ok(drawing);
    }
    let (blocks, model) = reader.entities(&mut drawing)?;
    super::read::assemble(drawing, &blocks, &model)
}
