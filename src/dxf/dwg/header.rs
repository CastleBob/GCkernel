//! Header variables (ODA specification chapter 9), walked as far as INSUNITS (R2000 and later; R13 and R14 have no
//! INSUNITS). From R2007 the header's text and handles live in their own streams, so only the data stream is walked
//! there (HANDSEED alone stays in the data stream).

use super::bits::{Bits, Ver, R};
use super::Source;

/// $INSUNITS when the header can be walked to it and gives a known unit code; None otherwise.
pub(super) fn insunits(src: &Source) -> Option<i32> {
    let v = walk(src).ok()?;
    (0..=24).contains(&v).then_some(v)
}

fn walk(src: &Source) -> R<i32> {
    let ver = src.ver;
    if ver < Ver::R2000 {
        return super::bits::bad("R13 and R14 headers have no INSUNITS");
    }
    let mut b = Bits::new(&src.header, ver);
    let mut none: Option<Bits> = None;
    // Text: inline before R2007, in the string stream after (nothing to skip here).
    macro_rules! tv { () => { if ver < Ver::R2007 { b.text()?; } } }
    // Handles: inline before R2007, in the handle stream after.
    macro_rules! h { () => { if ver < Ver::R2007 { b.h(0)?; } } }
    macro_rules! bits { ($n:expr) => { for _ in 0..$n { b.b()?; } } }
    macro_rules! bs { ($n:expr) => { for _ in 0..$n { b.bs()?; } } }
    macro_rules! bl { ($n:expr) => { for _ in 0..$n { b.bl()?; } } }
    macro_rules! bd { ($n:expr) => { for _ in 0..$n { b.bd()?; } } }
    macro_rules! p3 { ($n:expr) => { for _ in 0..$n { b.p3()?; } } }
    macro_rules! rd { ($n:expr) => { for _ in 0..$n { b.rd()?; } } }
    macro_rules! cmc { () => { b.cmc(&mut none)?; } }

    if ver >= Ver::R2007 {
        b.rl()?; // size in bits (the spec lists it for R2007; R2010-R2018 files carry it too)
    }
    if ver >= Ver::R2013 {
        b.bll()?; // REQUIREDVERSIONS
    }
    bd!(4);
    for _ in 0..4 {
        tv!();
    }
    bl!(2);
    if ver < Ver::R2004 {
        h!();
    }
    bits!(2); // DIMASO, DIMSHO
    bits!(7); // PLINEGEN .. LIMCHECK
    if ver >= Ver::R2004 {
        bits!(1);
    }
    bits!(4); // USRTIMER, SKPOLY, ANGDIR, SPLFRAME
    bits!(2); // MIRRTEXT, WORLDVIEW
    bits!(3); // TILEMODE, PLIMCHECK, VISRETAIN
    bits!(2); // DISPSILH, PELLIPSE
    bs!(1); // PROXYGRAPHICS
    bs!(5); // TREEDEPTH .. AUPREC
    bs!(1); // ATTMODE
    bs!(1); // PDMODE
    if ver >= Ver::R2004 {
        bl!(3);
    }
    bs!(5 + 14); // USERI1-5, SPLINESEGS .. TEXTQLTY
    bd!(21); // LTSCALE .. CELTSCALE
    if ver <= Ver::R2004 {
        tv!(); // MENUNAME
    }
    bl!(4);
    if ver >= Ver::R2004 {
        bl!(3);
    }
    bl!(4);
    cmc!(); // CECOLOR
    b.h(0)?; // HANDSEED, always in the data stream
    h!(); // CLAYER
    h!(); // TEXTSTYLE
    h!(); // CELTYPE
    if ver >= Ver::R2007 {
        h!(); // CMATERIAL
    }
    h!(); // DIMSTYLE
    h!(); // CMLSTYLE
    if ver >= Ver::R2000 {
        bd!(1); // PSVPSCALE
    }
    for _ in 0..2 {
        // Paper space, then model space: INSBASE, EXTMIN, EXTMAX, LIMMIN, LIMMAX, ELEVATION, UCSORG/X/Y, UCSNAME.
        p3!(3);
        rd!(4);
        bd!(1);
        p3!(3);
        h!();
        if ver >= Ver::R2000 {
            h!();
            bs!(1);
            h!();
            p3!(6);
        }
    }
    if ver >= Ver::R2000 {
        tv!(); // DIMPOST
        tv!(); // DIMAPOST
    }
    bd!(9); // DIMSCALE .. DIMTM
    if ver >= Ver::R2007 {
        bd!(2);
        bs!(1);
        cmc!();
    }
    bits!(6); // DIMTOL .. DIMSE2
    bs!(3); // DIMTAD, DIMZIN, DIMAZIN
    if ver >= Ver::R2007 {
        bs!(1);
    }
    bd!(8); // DIMTXT .. DIMGAP
    bd!(1); // DIMALTRND
    bits!(1);
    bs!(1);
    bits!(4);
    cmc!();
    cmc!();
    cmc!();
    bs!(11); // DIMADEC .. DIMJUST
    bits!(2);
    bs!(4);
    bits!(1);
    bs!(1); // DIMATFIT
    if ver >= Ver::R2007 {
        bits!(1);
    }
    if ver >= Ver::R2010 {
        bits!(1);
        bd!(1);
        tv!();
        bd!(1);
        tv!();
    }
    for _ in 0..5 {
        h!(); // DIMTXSTY, DIMLDRBLK, DIMBLK, DIMBLK1, DIMBLK2
    }
    if ver >= Ver::R2007 {
        for _ in 0..3 {
            h!();
        }
    }
    bs!(2); // DIMLWD, DIMLWE
    for _ in 0..9 {
        h!(); // table control objects
    }
    if ver <= Ver::R2000 {
        h!();
    }
    for _ in 0..3 {
        h!(); // dictionaries
    }
    bs!(2); // TSTACKALIGN, TSTACKSIZE
    tv!();
    tv!();
    for _ in 0..3 {
        h!();
    }
    if ver >= Ver::R2004 {
        h!();
        h!();
    }
    if ver >= Ver::R2007 {
        h!();
    }
    if ver >= Ver::R2013 {
        h!();
    }
    bl!(1); // flags
    let units = b.bs()? as i32;
    // CEPSNTYPE follows (0-3); anything else means the walk lost its place.
    if !(0..=3).contains(&b.bs()?) {
        return super::bits::bad("the header walk lost its place");
    }
    Ok(units)
}
