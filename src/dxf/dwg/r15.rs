//! The R13-R15 file organisation (ODA specification chapter 3): section-locator records at 0x15 give the header
//! variables, classes and object map at absolute file addresses; objects are read straight from the file.

use super::bits::{bad, le32, R};

pub(super) struct Locators {
    pub header: usize,
    pub classes: usize,
    pub handles: usize,
    pub handles_size: usize,
    /// Record 4, where the MEASUREMENT variable may be stored (template layout, spec chapter 22).
    pub template: Option<usize>,
}

pub(super) fn locators(file: &[u8]) -> R<Locators> {
    let count = le32(file, 0x15)? as usize;
    if !(3..=16).contains(&count) {
        return bad(format!("{count} section locators (the header is damaged)"));
    }
    let mut found = [None; 5];
    for i in 0..count {
        let at = 0x19 + 9 * i;
        let number = *file.get(at).ok_or_else(|| bad::<()>("the file ends early").unwrap_err())? as usize;
        let seeker = le32(file, at + 1)? as usize;
        let size = le32(file, at + 5)? as usize;
        if number < 5 {
            found[number] = Some((seeker, size));
        }
    }
    let get = |k: usize, what: &str| -> R<(usize, usize)> {
        let (seeker, size) = found[k].ok_or_else(|| bad::<()>(format!("the {what} locator is missing")).unwrap_err())?;
        if seeker >= file.len() {
            return bad(format!("the {what} lies past the end of the file"));
        }
        Ok((seeker, size))
    };
    let header = get(0, "header")?.0;
    let classes = get(1, "classes")?.0;
    let (handles, handles_size) = get(2, "object map")?;
    let template = found[4].map(|(seeker, _)| seeker);
    Ok(Locators { header, classes, handles, handles_size, template })
}
