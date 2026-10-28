//! Wire protocol v1 between GCmaker and gckernel-worker (agreed with Astra, battle plan 2026-10-27).
//! All little-endian.
//!
//! Request: magic "GCKQ", u32 version (1), u32 op, u32 units (1 = inch), f64 linear tolerance,
//!          u32 count, then count f64.
//! Reply:   magic "GCKR", u32 version (1), u32 status (0 ok, 1 unsupported, 2 invalid input, 3 internal),
//!          u32 count, then count f64 when ok, else count UTF-8 message bytes.
//!
//! Ops: 1 check (a body: validated and echoed back, followed by its V, E, F, L counts),
//!      2 extrude (a profile -> the prism body, followed by its exact volume),
//!      3 cut (two bodies back to back -> body less tool),
//!      4 measure (a body -> its volume, surface area, and 1 when both are closed form, 0 when by quadrature),
//!      5 rebuild prism (a body -> the kernel's own extrusion of it, when it is a straight prism; rebuild.rs).

use crate::neutral::{write_body, Reader};
use crate::{boolean, extrude, mass, rebuild, topology, KernelError};

pub const VERSION: u32 = 1;
pub const UNITS_INCH: u32 = 1;

pub const OP_CHECK: u32 = 1;
pub const OP_EXTRUDE: u32 = 2;
pub const OP_CUT: u32 = 3;
pub const OP_MEASURE: u32 = 4;
pub const OP_REBUILD_PRISM: u32 = 5;

pub struct Request {
    pub op: u32,
    pub units: u32,
    pub tolerance: f64,
    pub data: Vec<f64>,
}

fn u32_at(bytes: &[u8], at: usize) -> Result<u32, KernelError> {
    bytes.get(at..at + 4).map(|b| u32::from_le_bytes(b.try_into().unwrap()))
        .ok_or_else(|| KernelError::InvalidInput("request ended early".into()))
}

pub fn parse_request(bytes: &[u8]) -> Result<Request, KernelError> {
    if bytes.get(0..4) != Some(b"GCKQ") { return Err(KernelError::InvalidInput("not a GC Kernel request".into())); }
    let version = u32_at(bytes, 4)?;
    if version != VERSION { return Err(KernelError::Unsupported(format!("request version {version}"))); }
    let op = u32_at(bytes, 8)?;
    let units = u32_at(bytes, 12)?;
    let tolerance = bytes.get(16..24).map(|b| f64::from_le_bytes(b.try_into().unwrap()))
        .ok_or_else(|| KernelError::InvalidInput("request ended early".into()))?;
    let count = u32_at(bytes, 24)? as usize;
    let payload = &bytes[28..];
    if payload.len() != count * 8 {
        return Err(KernelError::InvalidInput(format!("request says {count} values but carries {} bytes", payload.len())));
    }
    let data = payload.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
    Ok(Request { op, units, tolerance, data })
}

pub fn encode_request(op: u32, tolerance: f64, data: &[f64]) -> Vec<u8> {
    let mut out = b"GCKQ".to_vec();
    for v in [VERSION, op, UNITS_INCH] { out.extend_from_slice(&v.to_le_bytes()); }
    out.extend_from_slice(&tolerance.to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    for v in data { out.extend_from_slice(&v.to_le_bytes()); }
    out
}

/// Runs one request.
pub fn handle(request: &Request) -> Result<Vec<f64>, KernelError> {
    if request.units != UNITS_INCH {
        return Err(KernelError::Unsupported(format!("units code {} (v1 is inches; convert at capture)", request.units)));
    }
    let tol = request.tolerance;
    if !tol.is_finite() || tol <= 0.0 { return Err(KernelError::InvalidInput("the tolerance must be positive".into())); }
    let mut reader = Reader::with_tolerance(&request.data, tol);
    match request.op {
        OP_CHECK => {
            let body = reader.body()?;
            reader.finish()?;
            let report = topology::check(&body, tol)?;
            let mut out = write_body(&body);
            out.extend_from_slice(&[report.vertices as f64, report.edges as f64, report.faces as f64, report.loops as f64]);
            Ok(out)
        }
        OP_EXTRUDE => {
            let profile = reader.profile()?;
            reader.finish()?;
            let body = extrude::extrude(&profile, tol)?;
            topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("extrude built an invalid solid: {e}")))?;
            let mut out = write_body(&body);
            out.push(extrude::volume(&profile));
            Ok(out)
        }
        OP_CUT => {
            let body = reader.body()?;
            let tool = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            topology::check(&tool, tol)?;
            let result = boolean::cut(&body, &tool, tol)?;
            // Complete removal is an explicit, successful empty result (Astra's acceptance cases).
            if !result.faces.is_empty() {
                topology::check(&result, tol).map_err(|e| KernelError::Internal(format!("cut built an invalid solid: {e}")))?;
            }
            Ok(write_body(&result))
        }
        OP_MEASURE => {
            let body = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            let m = mass::measure(&body, tol)?;
            Ok(vec![m.volume, m.area, if m.exact { 1.0 } else { 0.0 }])
        }
        OP_REBUILD_PRISM => {
            let body = reader.body()?;
            reader.finish()?;
            topology::check(&body, tol)?;
            let rebuilt = rebuild::rebuild_prism(&body, tol)?;
            topology::check(&rebuilt, tol).map_err(|e| KernelError::Internal(format!("rebuild built an invalid solid: {e}")))?;
            Ok(write_body(&rebuilt))
        }
        op => Err(KernelError::Unsupported(format!("op {op}"))),
    }
}

pub fn encode_reply(result: &Result<Vec<f64>, KernelError>) -> Vec<u8> {
    let mut out = b"GCKR".to_vec();
    out.extend_from_slice(&VERSION.to_le_bytes());
    match result {
        Ok(values) => {
            out.extend_from_slice(&0u32.to_le_bytes());
            out.extend_from_slice(&(values.len() as u32).to_le_bytes());
            for v in values { out.extend_from_slice(&v.to_le_bytes()); }
        }
        Err(e) => {
            let (status, message) = match e {
                KernelError::Unsupported(m) => (1u32, m),
                KernelError::InvalidInput(m) => (2u32, m),
                KernelError::Internal(m) => (3u32, m),
            };
            out.extend_from_slice(&status.to_le_bytes());
            out.extend_from_slice(&(message.len() as u32).to_le_bytes());
            out.extend_from_slice(message.as_bytes());
        }
    }
    out
}
