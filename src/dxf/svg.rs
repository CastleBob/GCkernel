//! SVG to `Drawing`, from the W3C SVG 1.1/2 specifications (document structure, coordinate systems and
//! transformations, path data grammar, basic shapes, and the elliptical arc implementation notes F.6.5/F.6.6). Our own
//! minimal XML reader; no third-party source.
//!
//! The drawing is turned right way up (SVG's y runs down). When the root's width and height carry an absolute unit
//! (mm, cm, in, pt, pc) coordinates come out in millimetres with $INSUNITS 4; otherwise they stay in SVG pixels and the
//! units are left unknown for the user to confirm (exporters disagree on what a pixel is).

use std::collections::HashMap;
use std::f64::consts::TAU;

use super::{fit, Curve, Drawing, Geometry, P2};
use crate::KernelError;

struct Node {
    name: String,
    attrs: Vec<(String, String)>,
    children: Vec<Node>,
}

impl Node {
    fn attr(&self, name: &str) -> Option<&str> {
        self.attrs.iter().find(|(k, _)| k == name).map(|(_, v)| v.as_str())
    }
}

fn bad(message: impl Into<String>) -> KernelError {
    KernelError::InvalidInput(message.into())
}

/// Whether the bytes look like an SVG (or other XML) document rather than a DXF.
pub fn looks_like(bytes: &[u8]) -> bool {
    let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    bytes.iter().find(|b| !b.is_ascii_whitespace()) == Some(&b'<')
}

// ---- XML -----------------------------------------------------------------------------------------------------------

fn parse_xml(text: &str) -> Result<Node, KernelError> {
    let b = text.as_bytes();
    let mut i = 0;
    let mut stack: Vec<Node> = vec![Node { name: String::new(), attrs: Vec::new(), children: Vec::new() }];
    while i < b.len() {
        if b[i] != b'<' {
            i += 1;
            continue;
        }
        let rest = &text[i..];
        if rest.starts_with("<!--") {
            i += rest.find("-->").map(|k| k + 3).ok_or_else(|| bad("SVG: unterminated comment"))?;
        } else if rest.starts_with("<![CDATA[") {
            i += rest.find("]]>").map(|k| k + 3).ok_or_else(|| bad("SVG: unterminated CDATA"))?;
        } else if rest.starts_with("<?") {
            i += rest.find("?>").map(|k| k + 2).ok_or_else(|| bad("SVG: unterminated processing instruction"))?;
        } else if rest.starts_with("<!") {
            // DOCTYPE, possibly with an internal subset in brackets.
            let mut depth = 0;
            let mut k = 2;
            while k < rest.len() {
                match rest.as_bytes()[k] {
                    b'[' => depth += 1,
                    b']' => depth -= 1,
                    b'>' if depth <= 0 => break,
                    _ => {}
                }
                k += 1;
            }
            i += k + 1;
        } else if rest.starts_with("</") {
            let end = rest.find('>').ok_or_else(|| bad("SVG: unterminated end tag"))?;
            let name = local(rest[2..end].trim());
            let node = stack.pop().ok_or_else(|| bad("SVG: unbalanced end tag"))?;
            if node.name != name || stack.is_empty() {
                return Err(bad(format!("SVG: </{name}> does not close <{}>", node.name)));
            }
            stack.last_mut().unwrap().children.push(node);
            i += end + 1;
        } else {
            let (node, closed, used) = start_tag(rest)?;
            i += used;
            if closed {
                stack.last_mut().unwrap().children.push(node);
            } else {
                stack.push(node);
            }
        }
    }
    while stack.len() > 1 {
        // Unclosed elements at the end of the file: keep what was read.
        let node = stack.pop().unwrap();
        stack.last_mut().unwrap().children.push(node);
    }
    let document = stack.pop().unwrap();
    document.children.into_iter().find(|n| n.name == "svg").ok_or_else(|| bad("this is not an SVG file (no <svg> element)"))
}

/// One start tag: the element, whether it closes itself, and the bytes it used.
fn start_tag(rest: &str) -> Result<(Node, bool, usize), KernelError> {
    let b = rest.as_bytes();
    let mut k = 1;
    while k < b.len() && !b[k].is_ascii_whitespace() && b[k] != b'>' && b[k] != b'/' {
        k += 1;
    }
    let name = local(&rest[1..k]);
    let mut attrs = Vec::new();
    loop {
        while k < b.len() && b[k].is_ascii_whitespace() {
            k += 1;
        }
        if k >= b.len() {
            return Err(bad(format!("SVG: unterminated <{name}>")));
        }
        if b[k] == b'>' {
            return Ok((Node { name, attrs, children: Vec::new() }, false, k + 1));
        }
        if b[k] == b'/' {
            let end = rest[k..].find('>').ok_or_else(|| bad(format!("SVG: unterminated <{name}>")))?;
            return Ok((Node { name, attrs, children: Vec::new() }, true, k + end + 1));
        }
        let start = k;
        while k < b.len() && b[k] != b'=' && !b[k].is_ascii_whitespace() && b[k] != b'>' && b[k] != b'/' {
            k += 1;
        }
        let key = rest[start..k].to_string();
        while k < b.len() && b[k].is_ascii_whitespace() {
            k += 1;
        }
        if k < b.len() && b[k] == b'=' {
            k += 1;
            while k < b.len() && b[k].is_ascii_whitespace() {
                k += 1;
            }
            let quote = *b.get(k).ok_or_else(|| bad(format!("SVG: unterminated <{name}>")))?;
            if quote != b'"' && quote != b'\'' {
                return Err(bad(format!("SVG: attribute {key} of <{name}> is not quoted")));
            }
            let end = rest[k + 1..].find(quote as char).ok_or_else(|| bad(format!("SVG: unterminated attribute {key}")))?;
            attrs.push((key, entities(&rest[k + 1..k + 1 + end])));
            k += end + 2;
        } else {
            attrs.push((key, String::new()));
        }
    }
}

fn local(name: &str) -> String {
    name.rsplit(':').next().unwrap_or(name).to_string()
}

fn entities(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(at) = rest.find('&') {
        out.push_str(&rest[..at]);
        rest = &rest[at..];
        let Some(end) = rest.find(';').filter(|e| *e <= 10) else {
            out.push('&');
            rest = &rest[1..];
            continue;
        };
        let name = &rest[1..end];
        let c = match name {
            "amp" => Some('&'),
            "lt" => Some('<'),
            "gt" => Some('>'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ if name.starts_with("#x") || name.starts_with("#X") => u32::from_str_radix(&name[2..], 16).ok().and_then(char::from_u32),
            _ if name.starts_with('#') => name[1..].parse().ok().and_then(char::from_u32),
            _ => None,
        };
        match c {
            Some(c) => out.push(c),
            None => out.push_str(&rest[..=end]),
        }
        rest = &rest[end + 1..];
    }
    out.push_str(rest);
    out
}

// ---- Numbers, lengths, transforms ----------------------------------------------------------------------------------

/// x -> [a c e; b d f] x.
#[derive(Clone, Copy, Debug)]
struct Affine([f64; 6]);

impl Affine {
    const IDENTITY: Affine = Affine([1.0, 0.0, 0.0, 1.0, 0.0, 0.0]);

    fn point(&self, p: P2) -> P2 {
        let m = self.0;
        [m[0] * p[0] + m[2] * p[1] + m[4], m[1] * p[0] + m[3] * p[1] + m[5]]
    }

    fn vector(&self, v: P2) -> P2 {
        let m = self.0;
        [m[0] * v[0] + m[2] * v[1], m[1] * v[0] + m[3] * v[1]]
    }

    /// self after other.
    fn then(&self, o: &Affine) -> Affine {
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

    fn translate(x: f64, y: f64) -> Affine {
        Affine([1.0, 0.0, 0.0, 1.0, x, y])
    }

    fn scale(x: f64, y: f64) -> Affine {
        Affine([x, 0.0, 0.0, y, 0.0, 0.0])
    }
}

/// Numbers in SVG's compact syntax: "1.5.5" is 1.5 then 0.5, "1-2" is 1 then -2, separators are spaces or commas.
struct Numbers<'a> {
    s: &'a [u8],
    at: usize,
}

impl<'a> Numbers<'a> {
    fn new(s: &'a str) -> Self {
        Numbers { s: s.as_bytes(), at: 0 }
    }

    fn skip_separators(&mut self) {
        while self.at < self.s.len() && (self.s[self.at].is_ascii_whitespace() || self.s[self.at] == b',') {
            self.at += 1;
        }
    }

    fn peek(&mut self) -> Option<u8> {
        self.skip_separators();
        self.s.get(self.at).copied()
    }

    fn number(&mut self) -> Option<f64> {
        self.skip_separators();
        let start = self.at;
        let s = self.s;
        let mut k = self.at;
        if k < s.len() && (s[k] == b'+' || s[k] == b'-') { k += 1; }
        let mut digits = false;
        while k < s.len() && s[k].is_ascii_digit() { k += 1; digits = true; }
        if k < s.len() && s[k] == b'.' {
            k += 1;
            while k < s.len() && s[k].is_ascii_digit() { k += 1; digits = true; }
        }
        if !digits {
            return None;
        }
        if k < s.len() && (s[k] == b'e' || s[k] == b'E') {
            let mut e = k + 1;
            if e < s.len() && (s[e] == b'+' || s[e] == b'-') { e += 1; }
            if e < s.len() && s[e].is_ascii_digit() {
                while e < s.len() && s[e].is_ascii_digit() { e += 1; }
                k = e;
            }
        }
        self.at = k;
        std::str::from_utf8(&s[start..k]).ok()?.parse::<f64>().ok().filter(|v| v.is_finite())
    }

    /// An arc flag: a single 0 or 1, which may run into the next number ("a1 1 0 00 1 1").
    fn flag(&mut self) -> Option<bool> {
        self.skip_separators();
        let c = *self.s.get(self.at)?;
        if c == b'0' || c == b'1' {
            self.at += 1;
            Some(c == b'1')
        } else {
            None
        }
    }
}

/// A length in user units (CSS pixels): 96 per inch.
fn length(s: Option<&str>) -> Option<f64> {
    let s = s?.trim();
    let mut n = Numbers::new(s);
    let v = n.number()?;
    let unit = s[n.at..].trim();
    let per = match unit {
        "" | "px" => 1.0,
        "in" => 96.0,
        "mm" => 96.0 / 25.4,
        "cm" => 96.0 / 2.54,
        "pt" => 96.0 / 72.0,
        "pc" => 16.0,
        _ => return None,
    };
    Some(v * per)
}

fn absolute(s: Option<&str>) -> bool {
    s.is_some_and(|s| ["mm", "cm", "in", "pt", "pc"].iter().any(|u| s.trim().ends_with(u)))
}

fn transform(s: &str) -> Result<Affine, KernelError> {
    let mut m = Affine::IDENTITY;
    let mut rest = s.trim();
    while !rest.is_empty() {
        let open = rest.find('(').ok_or_else(|| bad(format!("SVG: bad transform '{s}'")))?;
        let close = rest[open..].find(')').map(|k| k + open).ok_or_else(|| bad(format!("SVG: bad transform '{s}'")))?;
        let name = rest[..open].trim_matches(|c: char| c.is_whitespace() || c == ',');
        let mut n = Numbers::new(&rest[open + 1..close]);
        let mut v = Vec::new();
        while let Some(x) = n.number() {
            v.push(x);
        }
        let arg = |i: usize, default: f64| v.get(i).copied().unwrap_or(default);
        let t = match name {
            "matrix" if v.len() == 6 => Affine([v[0], v[1], v[2], v[3], v[4], v[5]]),
            "translate" if !v.is_empty() => Affine::translate(v[0], arg(1, 0.0)),
            "scale" if !v.is_empty() => Affine::scale(v[0], arg(1, v[0])),
            "rotate" if !v.is_empty() => {
                let (s, c) = v[0].to_radians().sin_cos();
                let r = Affine([c, s, -s, c, 0.0, 0.0]);
                let (cx, cy) = (arg(1, 0.0), arg(2, 0.0));
                Affine::translate(cx, cy).then(&r).then(&Affine::translate(-cx, -cy))
            }
            "skewX" if !v.is_empty() => Affine([1.0, 0.0, v[0].to_radians().tan(), 1.0, 0.0, 0.0]),
            "skewY" if !v.is_empty() => Affine([1.0, v[0].to_radians().tan(), 0.0, 1.0, 0.0, 0.0]),
            _ => return Err(bad(format!("SVG: bad transform '{s}'"))),
        };
        m = m.then(&t);
        rest = rest[close + 1..].trim_start_matches(|c: char| c.is_whitespace() || c == ',');
    }
    Ok(m)
}

fn style(node: &Node, key: &str) -> Option<String> {
    if let Some(style) = node.attr("style") {
        for part in style.split(';') {
            if let Some((k, v)) = part.split_once(':') {
                if k.trim() == key {
                    return Some(v.trim().to_string());
                }
            }
        }
    }
    node.attr(key).map(|v| v.trim().to_string())
}

// ---- Reading -------------------------------------------------------------------------------------------------------

pub fn read(bytes: &[u8]) -> Result<Drawing, KernelError> {
    let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    let text = String::from_utf8_lossy(bytes);
    let root = parse_xml(&text)?;

    let mut d = Drawing { version: "SVG".into(), ..Drawing::default() };
    d.layer_index("0");
    let view_box: Option<Vec<f64>> = root.attr("viewBox").map(|v| {
        let mut n = Numbers::new(v);
        std::iter::from_fn(|| n.number()).collect()
    }).filter(|v: &Vec<f64>| v.len() == 4 && v[2] > 0.0 && v[3] > 0.0);
    let width = length(root.attr("width")).filter(|w| *w > 0.0).or(view_box.as_ref().map(|v| v[2]));
    let height = length(root.attr("height")).filter(|h| *h > 0.0).or(view_box.as_ref().map(|v| v[3]));
    let mut viewport = Affine::IDENTITY;
    if let (Some(vb), Some(w), Some(h)) = (&view_box, width, height) {
        viewport = view_box_transform(vb, w, h, root.attr("preserveAspectRatio"));
    }
    let metric = absolute(root.attr("width")) && absolute(root.attr("height"));
    let to_mm = if metric { 25.4 / 96.0 } else { 1.0 };
    if metric {
        d.insunits = Some(4);
    } else {
        d.warn("the SVG gives no physical size; its coordinates are SVG pixels".into());
    }
    let page_height = height.unwrap_or(0.0);
    let flip = Affine([to_mm, 0.0, 0.0, -to_mm, 0.0, page_height * to_mm]);

    let mut ids: HashMap<String, &Node> = HashMap::new();
    collect_ids(&root, &mut ids);
    let mut w = Walk { ids: &ids, d: &mut d, uses: Vec::new(), group: 0, text: 0 };
    for child in &root.children {
        w.node(child, &flip.then(&viewport), 0)?;
        w.group += 1;
    }
    let text = w.text;
    for _ in 0..text {
        d.skip("TEXT");
    }
    Ok(d)
}

fn view_box_transform(vb: &[f64], w: f64, h: f64, aspect: Option<&str>) -> Affine {
    let (sx, sy) = (w / vb[2], h / vb[3]);
    let aspect = aspect.unwrap_or("xMidYMid meet").trim();
    if aspect.starts_with("none") {
        return Affine::scale(sx, sy).then(&Affine::translate(-vb[0], -vb[1]));
    }
    let s = if aspect.contains("slice") { sx.max(sy) } else { sx.min(sy) };
    let align = aspect.split_whitespace().next().unwrap_or("xMidYMid");
    let fx = if align.contains("xMin") { 0.0 } else if align.contains("xMax") { 1.0 } else { 0.5 };
    let fy = if align.contains("YMin") { 0.0 } else if align.contains("YMax") { 1.0 } else { 0.5 };
    let (tx, ty) = ((w - vb[2] * s) * fx, (h - vb[3] * s) * fy);
    Affine::translate(tx, ty).then(&Affine::scale(s, s)).then(&Affine::translate(-vb[0], -vb[1]))
}

fn collect_ids<'a>(node: &'a Node, ids: &mut HashMap<String, &'a Node>) {
    if let Some(id) = node.attr("id") {
        ids.entry(id.to_string()).or_insert(node);
    }
    for c in &node.children {
        collect_ids(c, ids);
    }
}

struct Walk<'a, 'b> {
    ids: &'b HashMap<String, &'a Node>,
    d: &'b mut Drawing,
    uses: Vec<String>,
    group: usize,
    text: usize,
}

impl Walk<'_, '_> {
    fn node(&mut self, node: &Node, parent: &Affine, layer: usize) -> Result<(), KernelError> {
        if style(node, "display").as_deref() == Some("none") || style(node, "visibility").as_deref() == Some("hidden") {
            return Ok(());
        }
        let mut m = *parent;
        if let Some(t) = node.attr("transform") {
            m = m.then(&transform(t)?);
        }
        let mut layer = layer;
        match node.name.as_str() {
            "g" | "a" | "switch" => {
                if node.attr("inkscape:groupmode") == Some("layer") {
                    if let Some(label) = node.attr("inkscape:label").or(node.attr("id")) {
                        layer = self.d.layer_index(label);
                    }
                }
                for c in &node.children {
                    self.node(c, &m, layer)?;
                }
            }
            "svg" => {
                // A nested viewport: x, y and its own viewBox.
                let (x, y) = (length(node.attr("x")).unwrap_or(0.0), length(node.attr("y")).unwrap_or(0.0));
                m = m.then(&Affine::translate(x, y));
                let vb: Option<Vec<f64>> = node.attr("viewBox").map(|v| {
                    let mut n = Numbers::new(v);
                    std::iter::from_fn(|| n.number()).collect()
                }).filter(|v: &Vec<f64>| v.len() == 4 && v[2] > 0.0 && v[3] > 0.0);
                if let (Some(vb), Some(w), Some(h)) = (vb, length(node.attr("width")), length(node.attr("height"))) {
                    m = m.then(&view_box_transform(&vb, w, h, node.attr("preserveAspectRatio")));
                }
                for c in &node.children {
                    self.node(c, &m, layer)?;
                }
            }
            "use" => {
                let href = node.attr("href").or(node.attr("xlink:href")).unwrap_or("");
                let Some(id) = href.strip_prefix('#') else { return Ok(()) };
                let Some(target) = self.ids.get(id).copied() else {
                    self.d.warn(format!("<use> refers to '#{id}', which is not in the file; it was left out"));
                    return Ok(());
                };
                if self.uses.iter().any(|u| u == id) || self.uses.len() > 64 {
                    self.d.warn(format!("'#{id}' uses itself; the repeated use was left out"));
                    return Ok(());
                }
                let (x, y) = (length(node.attr("x")).unwrap_or(0.0), length(node.attr("y")).unwrap_or(0.0));
                let placed = m.then(&Affine::translate(x, y));
                self.uses.push(id.to_string());
                if target.name == "symbol" {
                    let mut inner = placed;
                    if let Some(t) = target.attr("transform") {
                        inner = inner.then(&transform(t)?);
                    }
                    for c in &target.children {
                        self.node(c, &inner, layer)?;
                    }
                } else {
                    self.node(target, &placed, layer)?;
                }
                self.uses.pop();
            }
            "path" => {
                if let Some(data) = node.attr("d") {
                    self.path(data, &m, layer)?;
                }
            }
            "rect" => {
                let x = length(node.attr("x")).unwrap_or(0.0);
                let y = length(node.attr("y")).unwrap_or(0.0);
                let (w, h) = (length(node.attr("width")).unwrap_or(0.0), length(node.attr("height")).unwrap_or(0.0));
                if w <= 0.0 || h <= 0.0 {
                    return Ok(());
                }
                let (rx, ry) = match (length(node.attr("rx")), length(node.attr("ry"))) {
                    (Some(rx), Some(ry)) => (rx, ry),
                    (Some(r), None) | (None, Some(r)) => (r, r),
                    _ => (0.0, 0.0),
                };
                let (rx, ry) = (rx.clamp(0.0, w / 2.0), ry.clamp(0.0, h / 2.0));
                let data = if rx > 0.0 && ry > 0.0 {
                    format!(
                        "M{} {}H{}A{rx} {ry} 0 0 1 {} {}V{}A{rx} {ry} 0 0 1 {} {}H{}A{rx} {ry} 0 0 1 {} {}V{}A{rx} {ry} 0 0 1 {} {}Z",
                        x + rx, y, x + w - rx, x + w, y + ry, y + h - ry, x + w - rx, y + h, x + rx, x, y + h - ry, y + ry, x + rx, y)
                } else {
                    format!("M{x} {y}H{}V{}H{x}Z", x + w, y + h)
                };
                self.path(&data, &m, layer)?;
            }
            "circle" | "ellipse" => {
                let c = [length(node.attr("cx")).unwrap_or(0.0), length(node.attr("cy")).unwrap_or(0.0)];
                let (rx, ry) = if node.name == "circle" {
                    let r = length(node.attr("r")).unwrap_or(0.0);
                    (r, r)
                } else {
                    (length(node.attr("rx")).unwrap_or(0.0), length(node.attr("ry")).unwrap_or(0.0))
                };
                if rx > 0.0 && ry > 0.0 {
                    self.conic(&m, c, [rx, 0.0], [0.0, ry], 0.0, TAU, layer);
                }
            }
            "line" => {
                let a = [length(node.attr("x1")).unwrap_or(0.0), length(node.attr("y1")).unwrap_or(0.0)];
                let b = [length(node.attr("x2")).unwrap_or(0.0), length(node.attr("y2")).unwrap_or(0.0)];
                self.line(&m, a, b, layer);
            }
            "polyline" | "polygon" => {
                let mut n = Numbers::new(node.attr("points").unwrap_or(""));
                let mut pts = Vec::new();
                while let (Some(x), Some(y)) = (n.number(), n.number()) {
                    pts.push([x, y]);
                }
                for w in pts.windows(2) {
                    self.line(&m, w[0], w[1], layer);
                }
                if node.name == "polygon" && pts.len() > 2 {
                    self.line(&m, pts[pts.len() - 1], pts[0], layer);
                }
            }
            "text" => self.text += 1,
            "image" => self.d.skip("IMAGE"),
            // Definitions and paint servers are drawn only through <use>; metadata draws nothing.
            _ => {}
        }
        Ok(())
    }

    fn push(&mut self, geometry: Geometry, layer: usize) {
        self.d.curves.push(Curve { geometry, layer, color: 256, group: self.group });
    }

    fn line(&mut self, m: &Affine, a: P2, b: P2, layer: usize) {
        let (a, b) = (m.point(a), m.point(b));
        if fit::dist(a, b) > 1e-12 * (1.0 + fit::len(a).max(fit::len(b))) {
            self.push(Geometry::Line { a, b }, layer);
        }
    }

    /// c + u cos t + v sin t for t0 <= t <= t1, mapped through m.
    fn conic(&mut self, m: &Affine, c: P2, u: P2, v: P2, t0: f64, t1: f64, layer: usize) {
        let (c, u, v) = (m.point(c), m.vector(u), m.vector(v));
        let size = fit::len(u).max(fit::len(v));
        let area = fit::cross(u, v);
        if size <= 0.0 || area.abs() <= 1e-9 * fit::len(u) * fit::len(v) {
            return;
        }
        let full = t1 - t0 >= TAU - 1e-12;
        let circular = (fit::len(u) - fit::len(v)).abs() <= 1e-9 * size && fit::dot(u, v).abs() <= 1e-9 * size * size;
        let geometry = if circular && full {
            Geometry::Circle { c, r: fit::len(u) }
        } else if circular {
            let p0 = [u[0] * t0.cos() + v[0] * t0.sin(), u[1] * t0.cos() + v[1] * t0.sin()];
            Geometry::Arc { c, r: fit::len(u), start: p0[1].atan2(p0[0]), sweep: (t1 - t0) * area.signum() }
        } else {
            Geometry::Ellipse { c, u, v, t0, t1 }
        };
        self.push(geometry, layer);
    }

    fn bezier(&mut self, m: &Affine, points: &[P2], layer: usize) {
        let pts: Vec<P2> = points.iter().map(|p| m.point(*p)).collect();
        let degree = pts.len() - 1;
        if pts.windows(2).all(|w| fit::dist(w[0], w[1]) == 0.0) {
            return;
        }
        let mut knots = vec![0.0; degree + 1];
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        self.push(Geometry::Spline { degree, knots, points: pts, weights: Vec::new() }, layer);
    }

    /// SVG path data (W3C SVG "Paths": moveto, lineto, curveto, quadratic, smooth, elliptical arc, closepath).
    fn path(&mut self, data: &str, m: &Affine, layer: usize) -> Result<(), KernelError> {
        let mut n = Numbers::new(data);
        let (mut cur, mut start) = ([0.0, 0.0], [0.0, 0.0]);
        let mut last_cubic: Option<P2> = None;
        let mut last_quad: Option<P2> = None;
        let mut command = 0u8;
        let fail = || bad(format!("SVG: path data '{}' is not valid", data.chars().take(60).collect::<String>()));
        loop {
            match n.peek() {
                None => break,
                Some(c) if c.is_ascii_alphabetic() => {
                    command = c;
                    n.at += 1;
                }
                Some(_) if command == 0 => return Err(fail()),
                Some(_) => {} // the previous command repeats
            }
            let rel = command.is_ascii_lowercase();
            let at = |p: P2, cur: P2| if rel { [cur[0] + p[0], cur[1] + p[1]] } else { p };
            let pair = |n: &mut Numbers| -> Result<P2, KernelError> { Ok([n.number().ok_or_else(fail)?, n.number().ok_or_else(fail)?]) };
            let (mut cubic, mut quad) = (None, None);
            match command.to_ascii_uppercase() {
                b'M' => {
                    cur = at(pair(&mut n)?, cur);
                    start = cur;
                    command = if rel { b'l' } else { b'L' };
                }
                b'L' => {
                    let p = at(pair(&mut n)?, cur);
                    self.line(m, cur, p, layer);
                    cur = p;
                }
                b'H' => {
                    let x = n.number().ok_or_else(fail)?;
                    let p = [if rel { cur[0] + x } else { x }, cur[1]];
                    self.line(m, cur, p, layer);
                    cur = p;
                }
                b'V' => {
                    let y = n.number().ok_or_else(fail)?;
                    let p = [cur[0], if rel { cur[1] + y } else { y }];
                    self.line(m, cur, p, layer);
                    cur = p;
                }
                b'C' | b'S' => {
                    let c1 = if command.to_ascii_uppercase() == b'C' {
                        at(pair(&mut n)?, cur)
                    } else {
                        last_cubic.map_or(cur, |l| [2.0 * cur[0] - l[0], 2.0 * cur[1] - l[1]])
                    };
                    let c2 = at(pair(&mut n)?, cur);
                    let p = at(pair(&mut n)?, cur);
                    self.bezier(m, &[cur, c1, c2, p], layer);
                    cubic = Some(c2);
                    cur = p;
                }
                b'Q' | b'T' => {
                    let c1 = if command.to_ascii_uppercase() == b'Q' {
                        at(pair(&mut n)?, cur)
                    } else {
                        last_quad.map_or(cur, |l| [2.0 * cur[0] - l[0], 2.0 * cur[1] - l[1]])
                    };
                    let p = at(pair(&mut n)?, cur);
                    self.bezier(m, &[cur, c1, p], layer);
                    quad = Some(c1);
                    cur = p;
                }
                b'A' => {
                    let rx = n.number().ok_or_else(fail)?.abs();
                    let ry = n.number().ok_or_else(fail)?.abs();
                    let phi = n.number().ok_or_else(fail)?.to_radians();
                    let large = n.flag().ok_or_else(fail)?;
                    let sweep = n.flag().ok_or_else(fail)?;
                    let p = at(pair(&mut n)?, cur);
                    self.arc(m, cur, p, rx, ry, phi, large, sweep, layer);
                    cur = p;
                }
                b'Z' => {
                    self.line(m, cur, start, layer);
                    cur = start;
                    command = 0; // numbers may not follow a closepath
                }
                _ => return Err(fail()),
            }
            last_cubic = cubic;
            last_quad = quad;
        }
        Ok(())
    }

    /// The endpoint arc as centre form (SVG implementation notes F.6.5, radii scaled up when too small, F.6.6).
    #[allow(clippy::too_many_arguments)]
    fn arc(&mut self, m: &Affine, p1: P2, p2: P2, rx: f64, ry: f64, phi: f64, large: bool, sweep: bool, layer: usize) {
        if fit::dist(p1, p2) == 0.0 {
            return;
        }
        if rx == 0.0 || ry == 0.0 {
            self.line(m, p1, p2, layer);
            return;
        }
        let (s, c) = phi.sin_cos();
        let (dx, dy) = ((p1[0] - p2[0]) / 2.0, (p1[1] - p2[1]) / 2.0);
        let (x1, y1) = (c * dx + s * dy, -s * dx + c * dy);
        let lambda = (x1 * x1) / (rx * rx) + (y1 * y1) / (ry * ry);
        let (rx, ry) = if lambda > 1.0 { (rx * lambda.sqrt(), ry * lambda.sqrt()) } else { (rx, ry) };
        let num = rx * rx * ry * ry - rx * rx * y1 * y1 - ry * ry * x1 * x1;
        let den = rx * rx * y1 * y1 + ry * ry * x1 * x1;
        let k = (num / den).max(0.0).sqrt() * if large != sweep { 1.0 } else { -1.0 };
        let (cx1, cy1) = (k * rx * y1 / ry, -k * ry * x1 / rx);
        let centre = [c * cx1 - s * cy1 + (p1[0] + p2[0]) / 2.0, s * cx1 + c * cy1 + (p1[1] + p2[1]) / 2.0];
        let angle = |u: P2, v: P2| fit::cross(u, v).atan2(fit::dot(u, v));
        let a = [(x1 - cx1) / rx, (y1 - cy1) / ry];
        let b = [(-x1 - cx1) / rx, (-y1 - cy1) / ry];
        let theta1 = angle([1.0, 0.0], a);
        let mut delta = angle(a, b);
        if !sweep && delta > 0.0 { delta -= TAU; }
        if sweep && delta < 0.0 { delta += TAU; }
        let u = [rx * c, rx * s];
        let v = [-ry * s, ry * c];
        if delta > 0.0 {
            self.conic(m, centre, u, v, theta1, theta1 + delta, layer);
        } else {
            self.conic(m, centre, u, [-v[0], -v[1]], -theta1, -theta1 - delta, layer);
        }
    }
}
