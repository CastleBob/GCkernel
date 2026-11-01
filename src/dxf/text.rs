//! Text values to plain text (DXF Reference: TEXT group 1 special codes, "MTEXT inline codes", and the caret
//! notation for control characters).

/// A TEXT or ATTRIB value: %%d degree, %%p plus/minus, %%c diameter, %%% percent, %%nnn character nnn; %%u and %%o
/// (underline and overline toggles) are dropped; caret control characters (^J line feed, ^I tab, "^ " a caret).
pub fn text_plain(s: &str) -> String {
    percent_codes(&caret_codes(s))
}

/// An MTEXT value: paragraph and line breaks (\P, \N, \X) become '\n', \~ a space, \\ \{ \} the character itself,
/// stacked fractions (\S1/2; \S1#2; \S1^2;) "1/2"; formatting (font \f \F, height \H, width \W, obliquing \Q,
/// tracking \T, alignment \A, colour \C \c, paragraph \p) up to its ';' and the toggles \L \l \O \o \K \k are
/// dropped, as are grouping braces. Then the TEXT codes (%%d ...) as for `text_plain`.
pub fn mtext_plain(s: &str) -> String {
    let chars: Vec<char> = s.chars().collect();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        match c {
            '\\' if i + 1 < chars.len() => {
                let code = chars[i + 1];
                i += 2;
                match code {
                    'P' | 'N' | 'X' => out.push('\n'),
                    '~' => out.push(' '),
                    '\\' | '{' | '}' => out.push(code),
                    'L' | 'l' | 'O' | 'o' | 'K' | 'k' => {}
                    'S' => {
                        let end = until_semicolon(&chars, i);
                        let stacked: String = chars[i..end].iter().collect();
                        out.extend(stacked.chars().map(|c| if c == '^' || c == '#' { '/' } else { c }));
                        i = (end + 1).min(chars.len());
                    }
                    'f' | 'F' | 'H' | 'W' | 'Q' | 'T' | 'A' | 'C' | 'c' | 'p' => {
                        i = (until_semicolon(&chars, i) + 1).min(chars.len());
                    }
                    other => out.push(other),
                }
            }
            '{' | '}' => i += 1,
            _ => {
                out.push(c);
                i += 1;
            }
        }
    }
    text_plain(&out)
}

fn until_semicolon(chars: &[char], from: usize) -> usize {
    (from..chars.len()).find(|&k| chars[k] == ';').unwrap_or(chars.len())
}

fn caret_codes(s: &str) -> String {
    if !s.contains('^') {
        return s.to_string();
    }
    let chars: Vec<char> = s.chars().collect();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '^' && i + 1 < chars.len() && (chars[i + 1] == ' ' || ('@'..='_').contains(&chars[i + 1])) {
            match chars[i + 1] {
                ' ' => out.push('^'),
                'J' => out.push('\n'),
                'I' => out.push('\t'),
                _ => {}
            }
            i += 2;
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}

fn percent_codes(s: &str) -> String {
    if !s.contains("%%") {
        return s.to_string();
    }
    let chars: Vec<char> = s.chars().collect();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '%' && i + 2 < chars.len() && chars[i + 1] == '%' {
            let code = chars[i + 2];
            let digits = chars.get(i + 2..i + 5).filter(|d| d.iter().all(|c| c.is_ascii_digit()));
            match code.to_ascii_lowercase() {
                'd' => out.push('\u{00B0}'),
                'p' => out.push('\u{00B1}'),
                'c' => out.push('\u{00D8}'),
                '%' => out.push('%'),
                'u' | 'o' | 'k' => {}
                _ => {
                    if let Some(d) = digits {
                        let n: u32 = d.iter().collect::<String>().parse().unwrap_or(0);
                        if let Some(ch) = char::from_u32(n).filter(|ch| !ch.is_control()) {
                            out.push(ch);
                        }
                        i += 5;
                        continue;
                    }
                    out.push_str("%%");
                    i += 2;
                    continue;
                }
            }
            i += 3;
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}
