//! Op 16: [1, UTF8 path byte count, bytes] -> [1, body count,
//! (source id, UTF8 name count, bytes, length-to-inch factor, volume,
//!  neutral body count, neutral body)...]. No source or destination writes.
use super::{invalid, Result};
use crate::{neutral, KernelError};
/// Op 17: [1, count, (name byte count, UTF8 name bytes, body value count,
/// body)...] -> [1, byte count, ASCII STEP bytes]. No filesystem writes.
pub fn write_op(data: &[f64], tolerance: f64) -> Result<Vec<f64>> {
    let mut r = neutral::Reader::with_tolerance(data, tolerance);
    if r.next()? != 1. {
        return invalid("invalid STEP export envelope version");
    }
    let n = r.count()?;
    if n == 0 || n > 1024 {
        return invalid("invalid STEP export body count");
    }
    let mut items = Vec::with_capacity(n);
    for _ in 0..n {
        let count = r.count()?;
        if count > 32768 {
            return invalid("STEP export name exceeds 32768 bytes");
        }
        let bytes = (0..count)
            .map(|_| r.index(256, "name byte").map(|b| b as u8))
            .collect::<Result<Vec<_>>>()?;
        let name = String::from_utf8(bytes)
            .map_err(|_| KernelError::InvalidInput("STEP export name is not UTF-8".into()))?;
        let count = r.count()?;
        let values = (0..count).map(|_| r.next()).collect::<Result<Vec<_>>>()?;
        let mut body = neutral::Reader::with_tolerance(&values, tolerance);
        let parsed = body.body()?;
        body.finish()?;
        items.push((name, parsed));
    }
    r.finish()?;
    let pairs = items
        .iter()
        .map(|(name, body)| (name.as_str(), body))
        .collect::<Vec<_>>();
    let bytes = super::write(&pairs, tolerance)?.into_bytes();
    let mut out = vec![1., bytes.len() as f64];
    out.extend(bytes.into_iter().map(f64::from));
    Ok(out)
}
pub fn read_op(data: &[f64], tolerance: f64) -> Result<Vec<f64>> {
    if data.len() < 2 || (data[0] != 1.0 && data[0] != 2.0) {
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
    if data[0] == 2. {
        let instances = super::assembly::read_file(std::path::Path::new(&path), tolerance)?;
        let mut out = vec![2., instances.len() as f64];
        for instance in instances {
            out.push(instance.path.len() as f64);
            out.extend(instance.path.iter().map(|&id| id as f64));
            out.push(instance.representation as f64);
            let (matrix, t) = instance.placement.rows();
            out.extend(matrix.into_iter().flatten());
            out.extend([t.x, t.y, t.z]);
            let item = instance.item;
            let mut record = vec![item.entity_id as f64, item.name.len() as f64];
            record.extend(item.name.bytes().map(f64::from));
            record.extend([item.source_length_to_inches, item.volume]);
            let body = neutral::write_body(&item.body);
            record.push(body.len() as f64);
            record.extend(body);
            out.push(record.len() as f64);
            out.extend(record);
        }
        return Ok(out);
    }
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
