//! Op 22 request [1,pathByteCount,UTF8 bytes], reply is the same detached body
//! envelope as STEP schema 1. Source IDs are IGES directory sequence numbers.
pub fn read_op(data: &[f64], tolerance: f64) -> super::Result<Vec<f64>> {
    let mut r = crate::neutral::Reader::new(data);
    if r.next()? != 1. {
        return super::bad("unknown request version");
    }
    let n = r.count()?;
    if n == 0 || n > 32768 {
        return super::bad("invalid path length");
    }
    let bytes = (0..n)
        .map(|_| r.index(256, "IGES path byte").map(|v| v as u8))
        .collect::<super::Result<Vec<_>>>()?;
    r.finish()?;
    if bytes.contains(&0) {
        return super::bad("path contains a zero byte");
    }
    let path = String::from_utf8(bytes).map_err(|_| super::error("path is not UTF-8"))?;
    let bodies = super::read_file(std::path::Path::new(&path), tolerance)?;
    let mut out = vec![1., bodies.len() as f64];
    for b in bodies {
        out.extend([b.entity_id as f64, b.name.len() as f64]);
        out.extend(b.name.bytes().map(f64::from));
        out.extend([b.source_length_to_inches, b.volume]);
        let body = crate::neutral::write_body(&b.body);
        out.push(body.len() as f64);
        out.extend(body);
    }
    Ok(out)
}
