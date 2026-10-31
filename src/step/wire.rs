//! Op 16: [1, UTF8 path byte count, bytes] -> [1, body count,
//! (source id, UTF8 name count, bytes, length-to-inch factor, volume,
//!  neutral body count, neutral body)...]. No source or destination writes.
use super::{invalid, Result};
use crate::{neutral, KernelError};
pub fn read_op(data: &[f64], tolerance: f64) -> Result<Vec<f64>> {
    if data.len() < 2 || data[0] != 1.0 {
        return invalid("invalid STEP request version/envelope");
    }
    let n = data[1];
    if !n.is_finite() || n.fract() != 0.0 || n < 1.0 || n > 32768.0 || data.len() != n as usize + 2
    {
        return invalid("invalid STEP path length");
    }
    let mut bytes = Vec::with_capacity(n as usize);
    for &v in &data[2..] {
        if !v.is_finite() || v.fract() != 0.0 || v < 1.0 || v > 255.0 {
            return invalid("invalid STEP path byte");
        }
        bytes.push(v as u8);
    }
    let path = String::from_utf8(bytes)
        .map_err(|_| KernelError::InvalidInput("STEP path is not UTF-8".into()))?;
    let imported = super::read_file(std::path::Path::new(&path), tolerance)?;
    let mut out = vec![1.0, imported.bodies.len() as f64];
    for item in imported.bodies {
        out.push(item.entity_id as f64);
        out.push(item.name.len() as f64);
        out.extend(item.name.bytes().map(f64::from));
        out.extend_from_slice(&[item.source_length_to_inches, item.volume]);
        let body = neutral::write_body(&item.body);
        out.push(body.len() as f64);
        out.extend(body);
    }
    Ok(out)
}
