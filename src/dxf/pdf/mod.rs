//! PDF (and Adobe Illustrator files saved with PDF compatibility) to `Drawing`, written from ISO 32000-1:2008
//! (PDF 1.7) and, for AES-256 encryption, the revision 6 password algorithms of ISO 32000-2; our own Deflate, LZW,
//! MD5, SHA-2, RC4 and AES from their RFC / FIPS texts. No third-party PDF, compression or cryptography source.
//!
//! The first page (in page-tree order) is read: every stroked or filled path, through forms, in inches (PDF units are
//! 1/72 inch times the page's /UserUnit) with $INSUNITS 1, measured from the lower-left corner of the page's crop box
//! (else media box) after its /Rotate. Text, images, shadings and annotations are counted, not drawn. Optional
//! content groups become layers; content on a group that is off in the default configuration is left out.

mod content;
mod crypt;
mod filters;
mod inflate;
mod lexer;
mod shapes;
mod xref;

use std::collections::HashSet;
use std::rc::Rc;

use lexer::{Dict, Obj};
use xref::Doc;

use super::Drawing;
use crate::KernelError;

/// x -> [a c e; b d f] x (the PDF matrix [a b c d e f]).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Affine(pub [f64; 6]);

impl Affine {
    pub const IDENTITY: Affine = Affine([1.0, 0.0, 0.0, 1.0, 0.0, 0.0]);

    pub fn point(&self, p: [f64; 2]) -> [f64; 2] {
        let m = self.0;
        [m[0] * p[0] + m[2] * p[1] + m[4], m[1] * p[0] + m[3] * p[1] + m[5]]
    }

    /// self after other.
    pub fn then(&self, o: &Affine) -> Affine {
        let (a, b) = (self.0, o.0);
        Affine([
            a[0] * b[0] + a[2] * b[1],
            a[1] * b[0] + a[3] * b[1],
            a[0] * b[2] + a[2] * b[3],
            a[1] * b[2] + a[3] * b[3],
            a[0] * b[4] + a[2] * b[5] + a[4],
            a[1] * b[4] + a[3] * b[5] + a[5],
        ])
    }

    /// The largest stretch of a column (how far a unit step can move).
    pub fn scale(&self) -> f64 {
        let m = self.0;
        (m[0] * m[0] + m[1] * m[1]).max(m[2] * m[2] + m[3] * m[3]).sqrt()
    }
}

/// Whether the bytes are a PDF (also an Illustrator file with its PDF part): "%PDF-" in the first 1024 bytes.
pub fn looks_like(bytes: &[u8]) -> bool {
    bytes[..bytes.len().min(1024)].windows(5).any(|w| w == b"%PDF-")
}

/// A PostScript-only Illustrator or EPS file (no PDF part), including EPS files with a DOS binary header.
pub fn looks_like_postscript(bytes: &[u8]) -> bool {
    !looks_like(bytes) && (bytes.starts_with(b"%!PS") || bytes.starts_with(&[0xC5, 0xD0, 0xD3, 0xC6]))
}

pub fn postscript_refusal() -> KernelError {
    KernelError::Unsupported("this Illustrator file has no PDF part (saved without 'Create PDF Compatible File', or \
        Illustrator 8 or older); re-save it with PDF compatibility, or as PDF/SVG/DXF".into())
}

/// The first page's geometry.
pub fn read(bytes: &[u8]) -> Result<Drawing, KernelError> {
    read_parts(bytes, true)
}

/// Units, version and layers only (no curves).
pub fn read_header(bytes: &[u8]) -> Result<Drawing, KernelError> {
    read_parts(bytes, false)
}

fn bad(m: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(m.into())
}

/// The page attributes a page inherits from its ancestors (7.7.3.4).
#[derive(Clone, Default)]
struct Inherited {
    resources: Option<Rc<Obj>>,
    media: Option<Rc<Obj>>,
    crop: Option<Rc<Obj>>,
    rotate: Option<Rc<Obj>>,
}

fn read_parts(bytes: &[u8], curves: bool) -> Result<Drawing, KernelError> {
    if !looks_like(bytes) {
        return Err(if looks_like_postscript(bytes) { postscript_refusal() } else { bad("the file is not a PDF") });
    }
    let mut doc = Doc::open(bytes)?;
    let root = doc.trailer.get(b"Root").cloned().ok_or_else(|| bad("the PDF has no document catalog"))?;
    let catalog = doc.resolve(&root)?;
    let catalog = catalog.as_dict().cloned().ok_or_else(|| bad("the PDF's document catalog is not a dictionary"))?;
    let mut version = doc.version.clone();
    if let Some(v) = catalog.get(b"Version").and_then(Obj::as_name) {
        let v = String::from_utf8_lossy(v).into_owned();
        if v.as_str() > version.as_str() {
            version = v;
        }
    }
    let mut d = Drawing { version: format!("PDF {version}"), insunits: Some(1), ..Drawing::default() };
    d.layer_index("0");
    let oc = optional_content(&mut doc, &catalog, &mut d)?;

    let pages = doc.get(&catalog, b"Pages")?.ok_or_else(|| bad("the PDF has no pages"))?;
    let mut inherited = Inherited::default();
    let mut seen = HashSet::new();
    let mut visits = 0usize;
    let root_num = match catalog.get(b"Pages") { Some(Obj::Ref(n, _)) => Some(*n), _ => None };
    let page = first_page(&mut doc, pages.clone(), root_num, &mut inherited, &mut seen, &mut visits, 0)?
        .ok_or_else(|| bad("the PDF has no pages"))?;
    let count = match pages.as_dict().and_then(|p| p.get(b"Count")).and_then(Obj::as_i64) {
        Some(n) if n >= 1 => n as usize,
        _ => {
            let (mut seen, mut visits) = (HashSet::new(), 0);
            count_pages(&mut doc, pages.clone(), &mut seen, &mut visits, 0)
        }
    };
    if count > 1 {
        d.warn(format!("page 1 of {count} was read"));
    }
    let page = page.as_dict().cloned().unwrap_or_default();

    // The page's box, rotation and units.
    let rect = |o: Option<Rc<Obj>>, doc: &mut Doc| -> Option<[f64; 4]> {
        let o = o?;
        let a = o.as_array()?;
        let v: Vec<f64> = a.iter().filter_map(|x| doc.resolve(x).ok().and_then(|r| r.as_f64())).collect();
        (v.len() == 4).then(|| [v[0].min(v[2]), v[1].min(v[3]), v[0].max(v[2]), v[1].max(v[3])])
            .filter(|r| r[2] > r[0] && r[3] > r[1])
    };
    let own = |key: &[u8], doc: &mut Doc| doc.get(&page, key).ok().flatten();
    let media = own(b"MediaBox", &mut doc).or(inherited.media.clone());
    let crop = own(b"CropBox", &mut doc).or(inherited.crop.clone());
    let media = match rect(media, &mut doc) {
        Some(r) => r,
        None => {
            d.warn("the PDF page has no media box; a US Letter page was assumed".into());
            [0.0, 0.0, 612.0, 792.0]
        }
    };
    let bx = rect(crop, &mut doc).unwrap_or(media);
    let rotate = own(b"Rotate", &mut doc).or(inherited.rotate.clone()).and_then(|r| r.as_i64()).unwrap_or(0);
    let mut rotate = rotate.rem_euclid(360);
    if rotate % 90 != 0 {
        d.warn(format!("the PDF page's /Rotate {rotate} is not a multiple of 90; it was ignored"));
        rotate = 0;
    }
    let unit = own(b"UserUnit", &mut doc).and_then(|u| u.as_f64()).filter(|u| *u > 0.0).unwrap_or(1.0);
    let (w, h) = (bx[2] - bx[0], bx[3] - bx[1]);
    let turn = match rotate {
        90 => Affine([0.0, -1.0, 1.0, 0.0, 0.0, w]),
        180 => Affine([-1.0, 0.0, 0.0, -1.0, w, h]),
        270 => Affine([0.0, 1.0, -1.0, 0.0, h, 0.0]),
        _ => Affine::IDENTITY,
    };
    let s = unit / 72.0;
    let page_matrix = Affine([s, 0.0, 0.0, s, 0.0, 0.0]).then(&turn).then(&Affine([1.0, 0.0, 0.0, 1.0, -bx[0], -bx[1]]));

    if curves {
        let resources = match own(b"Resources", &mut doc).or(inherited.resources.clone()) {
            Some(r) => r.as_dict().cloned().unwrap_or_default(),
            None => Dict::default(),
        };
        let mut data: Vec<u8> = Vec::new();
        if let Some(c) = page.get(b"Contents") {
            let c = doc.resolve(c)?;
            let parts: Vec<Obj> = match &*c {
                Obj::Array(a) => a.clone(),
                other => vec![other.clone()],
            };
            for part in parts.iter().take(100_000) {
                let p = doc.resolve(part)?;
                if let Obj::Stream(s) = &*p {
                    data.extend_from_slice(&doc.decode(s)?);
                    data.push(b'\n');
                }
            }
        }
        let annotations = doc.get(&page, b"Annots")?.and_then(|a| a.as_array().map(<[Obj]>::len)).unwrap_or(0);
        let mut painter = content::Painter::new(&mut doc, &mut d, &oc);
        painter.run(&data, &resources, page_matrix, content::Mark { layer: 0, hidden_by: None }, 0)?;
        let hidden: Vec<usize> = painter.hidden_layers().to_vec();
        for l in hidden {
            let name = d.layers[l].name.clone();
            d.warn(format!("layer '{name}' is off in the PDF; its content was left out"));
        }
        for _ in 0..annotations.min(100_000) {
            d.skip("ANNOTATION");
        }
    }
    for n in std::mem::take(&mut doc.notes) {
        d.warn(n);
    }
    Ok(d)
}

/// Optional-content groups (8.11.4): each group's layer and its state in the default configuration /D.
fn optional_content(doc: &mut Doc, catalog: &Dict, d: &mut Drawing) -> Result<content::Oc, KernelError> {
    let mut oc = content::Oc { groups: Default::default() };
    let Some(props) = doc.get(catalog, b"OCProperties")? else { return Ok(oc) };
    let Some(props) = props.as_dict().cloned() else { return Ok(oc) };
    let refs = |o: Option<Rc<Obj>>| -> Vec<u32> {
        o.and_then(|a| a.as_array().map(|a| a.iter().filter_map(|x| if let Obj::Ref(n, _) = x { Some(*n) } else { None }).collect()))
            .unwrap_or_default()
    };
    let all = refs(doc.get(&props, b"OCGs")?);
    let config = doc.get(&props, b"D")?.and_then(|c| c.as_dict().cloned()).unwrap_or_default();
    let base_on = config.get(b"BaseState").and_then(Obj::as_name) != Some(b"OFF");
    let on = refs(doc.get(&config, b"ON")?);
    let off = refs(doc.get(&config, b"OFF")?);
    for n in all.into_iter().take(100_000) {
        let g = doc.object(n)?;
        let name = g.as_dict().and_then(|g| g.get(b"Name")).and_then(Obj::as_str).map(lexer::text).unwrap_or_default();
        let state = if off.contains(&n) { false } else if on.contains(&n) { true } else { base_on };
        let layer = d.layer_index(&name);
        if !state {
            d.layers[layer].off = true;
        }
        oc.groups.insert(n, (layer, state));
    }
    Ok(oc)
}

const MAX_TREE_DEPTH: usize = 64;
const MAX_TREE_NODES: usize = 1_000_000;

fn first_page(doc: &mut Doc, node: Rc<Obj>, num: Option<u32>, inherited: &mut Inherited, seen: &mut HashSet<u32>,
    visits: &mut usize, depth: usize) -> Result<Option<Rc<Obj>>, KernelError> {
    *visits += 1;
    if depth > MAX_TREE_DEPTH || *visits > MAX_TREE_NODES || num.is_some_and(|n| !seen.insert(n)) {
        return Ok(None);
    }
    let Some(dict) = node.as_dict() else { return Ok(None) };
    let kids = dict.get(b"Kids");
    let is_page = dict.get(b"Type").and_then(Obj::as_name) == Some(b"Page") || kids.is_none();
    if is_page {
        return Ok(Some(node.clone()));
    }
    let mut mine = inherited.clone();
    for (key, slot) in [(&b"Resources"[..], &mut mine.resources), (b"MediaBox", &mut mine.media),
        (b"CropBox", &mut mine.crop), (b"Rotate", &mut mine.rotate)] {
        if let Some(v) = doc.get(dict, key)? {
            *slot = Some(v);
        }
    }
    let kids = doc.resolve(kids.unwrap_or(&Obj::Null))?;
    for kid in kids.as_array().unwrap_or(&[]) {
        let n = if let Obj::Ref(r, _) = kid { Some(*r) } else { None };
        let k = doc.resolve(kid)?;
        let mut below = mine.clone();
        if let Some(p) = first_page(doc, k, n, &mut below, seen, visits, depth + 1)? {
            *inherited = below;
            return Ok(Some(p));
        }
    }
    Ok(None)
}

fn count_pages(doc: &mut Doc, node: Rc<Obj>, seen: &mut HashSet<*const Obj>, visits: &mut usize, depth: usize) -> usize {
    *visits += 1;
    if depth > MAX_TREE_DEPTH || *visits > MAX_TREE_NODES || !seen.insert(Rc::as_ptr(&node)) {
        return 0;
    }
    let Some(dict) = node.as_dict() else { return 0 };
    let Some(kids) = dict.get(b"Kids") else { return 1 };
    let Ok(kids) = doc.resolve(kids) else { return 0 };
    let mut n = 0;
    for kid in kids.as_array().unwrap_or(&[]) {
        if let Ok(k) = doc.resolve(kid) {
            n += count_pages(doc, k, seen, visits, depth + 1);
        }
    }
    n
}
