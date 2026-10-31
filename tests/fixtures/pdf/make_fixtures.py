"""Writes the PDF fixtures for tests/pdf.rs (run: python -I tests/fixtures/pdf/make_fixtures.py).

Standard library only. The PDF writer, the LZW and PNG-predictor encoders, RC4, AES (S-box computed from its
GF(2^8) definition in FIPS 197) and the standard security handler's algorithms (ISO 32000-1 7.6.3; ISO 32000-2
Algorithm 2.B for revision 6) are written here from the specifications, independently of the Rust reader.
Output is deterministic (fixed salts and IVs) so the fixtures are reproducible.
"""

import hashlib
import os
import struct
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))


# ---- Objects -------------------------------------------------------------------------------------------------------

class Name(str):
    pass


class Ref:
    def __init__(self, n):
        self.n = n


class Stream:
    def __init__(self, d, data, length=None, length_ref=None, encrypt=True):
        self.d, self.data, self.length, self.length_ref, self.encrypt = d, data, length, length_ref, encrypt


class Raw:
    """Bytes written as they are."""
    def __init__(self, b):
        self.b = b


def num(x):
    if isinstance(x, bool):
        return b"true" if x else b"false"
    if isinstance(x, int):
        return str(x).encode()
    s = ("%.6f" % x).rstrip("0").rstrip(".")
    return (s if s not in ("", "-0") else "0").encode()


def ser(o, crypt=None):
    """crypt: function bytes -> bytes for strings of this object (None: plain)."""
    if isinstance(o, Raw):
        return o.b
    if isinstance(o, Name):
        out = b"/"
        for ch in o.encode("latin-1"):
            out += bytes([ch]) if 0x21 <= ch <= 0x7E and ch not in b"#()<>[]{}/%" else b"#%02X" % ch
        return out
    if o is None:
        return b"null"
    if isinstance(o, (bool, int, float)):
        return num(o)
    if isinstance(o, bytes):
        b = crypt(o) if crypt else o
        return b"<" + b.hex().encode() + b">"
    if isinstance(o, Ref):
        return b"%d 0 R" % o.n
    if isinstance(o, list):
        return b"[" + b" ".join(ser(x, crypt) for x in o) + b"]"
    if isinstance(o, dict):
        return b"<<" + b"".join(ser(Name(k)) + b" " + ser(v, crypt) + b"\n" for k, v in o.items()) + b">>"
    raise TypeError(o)


# ---- Filters -------------------------------------------------------------------------------------------------------

def lzw_encode(data, early=1, clear_after=None):
    """ISO 32000-1 7.4.4.2: 9 to 12 bit codes, high bit first, clear 256, EOD 257."""
    bits = []

    def emit(code, width):
        bits.extend((code >> (width - 1 - i)) & 1 for i in range(width))

    def width_for(next_code):
        n = next_code + early
        return 12 if n >= 2048 else 11 if n >= 1024 else 10 if n >= 512 else 9

    table = {bytes([i]): i for i in range(256)}
    next_code = 258
    emit(256, 9)
    w = b""
    emitted = 0
    for byte in data:
        wc = w + bytes([byte])
        if wc in table:
            w = wc
            continue
        emit(table[w], width_for(next_code))
        emitted += 1
        if next_code + early < 4095:
            table[wc] = next_code
            next_code += 1
        w = bytes([byte])
        if clear_after and emitted == clear_after:
            # Flush the pending byte, then start over (away from a code-width boundary, so the decoder agrees).
            emit(table[w], width_for(next_code))
            emit(256, width_for(next_code))
            table = {bytes([i]): i for i in range(256)}
            next_code = 258
            w = b""
    if w:
        emit(table[w], width_for(next_code))
        if next_code + early < 4095:
            next_code += 1
    emit(257, width_for(next_code))
    while len(bits) % 8:
        bits.append(0)
    return bytes(int("".join(map(str, bits[i:i + 8])), 2) for i in range(0, len(bits), 8))


def ascii85(data):
    out = b""
    for i in range(0, len(data), 4):
        chunk = data[i:i + 4]
        n = len(chunk)
        v = int.from_bytes(chunk + b"\0" * (4 - n), "big")
        if n == 4 and v == 0:
            out += b"z"
            continue
        digits = []
        for _ in range(5):
            digits.append(v % 85)
            v //= 85
        out += bytes(33 + d for d in reversed(digits))[:n + 1]
    return out + b"~>"


def run_length(data):
    out = b""
    i = 0
    while i < len(data):
        j = i
        while j < len(data) and j - i < 128 and data[j] == data[i]:
            j += 1
        if j - i >= 3:
            out += bytes([257 - (j - i), data[i]])
            i = j
        else:
            j = i
            while j < len(data) and j - i < 128 and not (j + 2 < len(data) and data[j] == data[j + 1] == data[j + 2]):
                j += 1
            out += bytes([j - i - 1]) + data[i:j]
            i = j
    return out + b"\x80"


def png_predict(data, columns, types):
    """PNG row filters (RFC 2083) with 1-byte samples; `types` cycles the filter type per row."""
    data = data + b" " * (-len(data) % columns)
    out = b""
    prior = bytes(columns)
    for r in range(0, len(data), columns):
        row = data[r:r + columns]
        t = types[(r // columns) % len(types)]
        enc = bytearray()
        for i in range(columns):
            a = row[i - 1] if i >= 1 else 0
            b = prior[i]
            c = prior[i - 1] if i >= 1 else 0
            if t == 0:
                p = 0
            elif t == 1:
                p = a
            elif t == 2:
                p = b
            elif t == 3:
                p = (a + b) // 2
            else:
                pp = a + b - c
                pa, pb, pc = abs(pp - a), abs(pp - b), abs(pp - c)
                p = a if pa <= pb and pa <= pc else b if pb <= pc else c
            enc.append((row[i] - p) % 256)
        out += bytes([t]) + bytes(enc)
        prior = row
    return out


# ---- Ciphers -------------------------------------------------------------------------------------------------------

def rc4(key, data):
    s = list(range(256))
    j = 0
    for i in range(256):
        j = (j + s[i] + key[i % len(key)]) % 256
        s[i], s[j] = s[j], s[i]
    i = j = 0
    out = bytearray()
    for ch in data:
        i = (i + 1) % 256
        j = (j + s[i]) % 256
        s[i], s[j] = s[j], s[i]
        out.append(ch ^ s[(s[i] + s[j]) % 256])
    return bytes(out)


def gf_mul(a, b):
    p = 0
    while b:
        if b & 1:
            p ^= a
        a = ((a << 1) ^ 0x1B) & 0xFF if a & 0x80 else a << 1
        b >>= 1
    return p


def make_sbox():
    """FIPS 197 5.1.1: multiplicative inverse in GF(2^8), then the affine transformation with c = 0x63."""
    inv = [0] * 256
    for x in range(1, 256):
        for y in range(1, 256):
            if gf_mul(x, y) == 1:
                inv[x] = y
                break
    box = []
    for x in range(256):
        b = inv[x]
        r = 0
        for i in range(8):
            bit = ((b >> i) ^ (b >> ((i + 4) % 8)) ^ (b >> ((i + 5) % 8)) ^ (b >> ((i + 6) % 8)) ^ (b >> ((i + 7) % 8)) ^ (0x63 >> i)) & 1
            r |= bit << i
        box.append(r)
    return box


SBOX = make_sbox()
assert SBOX[0x00] == 0x63 and SBOX[0x53] == 0xED


class AES:
    def __init__(self, key):
        nk = len(key) // 4
        self.nr = nk + 6
        w = [list(key[4 * i:4 * i + 4]) for i in range(nk)]
        rcon = 1
        for i in range(nk, 4 * (self.nr + 1)):
            t = list(w[i - 1])
            if i % nk == 0:
                t = [SBOX[t[1]] ^ rcon, SBOX[t[2]], SBOX[t[3]], SBOX[t[0]]]
                rcon = gf_mul(rcon, 2)
            elif nk > 6 and i % nk == 4:
                t = [SBOX[x] for x in t]
            w.append([w[i - nk][j] ^ t[j] for j in range(4)])
        self.rk = [sum(w[4 * r:4 * r + 4], []) for r in range(self.nr + 1)]

    def block(self, b):
        s = [b[i] ^ self.rk[0][i] for i in range(16)]
        for rnd in range(1, self.nr + 1):
            s = [SBOX[x] for x in s]
            s = [s[r + 4 * ((c + r) % 4)] for c in range(4) for r in range(4)]
            if rnd != self.nr:
                t = []
                for c in range(4):
                    a = s[4 * c:4 * c + 4]
                    t += [gf_mul(a[0], 2) ^ gf_mul(a[1], 3) ^ a[2] ^ a[3],
                          a[0] ^ gf_mul(a[1], 2) ^ gf_mul(a[2], 3) ^ a[3],
                          a[0] ^ a[1] ^ gf_mul(a[2], 2) ^ gf_mul(a[3], 3),
                          gf_mul(a[0], 3) ^ a[1] ^ a[2] ^ gf_mul(a[3], 2)]
                s = t
            s = [s[i] ^ self.rk[rnd][i] for i in range(16)]
        return bytes(s)

    def cbc(self, iv, data):
        out = b""
        prev = iv
        for i in range(0, len(data), 16):
            prev = self.block(bytes(x ^ y for x, y in zip(data[i:i + 16], prev)))
            out += prev
        return out


assert AES(bytes.fromhex("2b7e151628aed2a6abf7158809cf4f3c")).block(
    bytes.fromhex("3243f6a8885a308d313198a2e0370734")).hex() == "3925841d02dc09fbdc118597196a0b32"


def aes_pdf(key, data, seed):
    """7.6.2: a 16-byte IV first, then CBC over PKCS#5-padded data."""
    iv = hashlib.md5(b"iv" + seed).digest()
    pad = 16 - len(data) % 16
    return iv + AES(key).cbc(iv, data + bytes([pad]) * pad)


PAD = bytes.fromhex("28BF4E5E4E758A4164004E56FFFA01082E2E00B6D0683E802F0CA9FE6453697A")


def padded(pw):
    return (pw + PAD)[:32]


class Standard:
    """Standard security handler, revisions 2 (RC4 40-bit), 3 and 4 (RC4 or AES-128) or 6 (AES-256)."""

    def __init__(self, revision, user=b"", owner=b"owner", aes=False, encrypt_metadata=True, id0=b""):
        self.r, self.aes, self.id0 = revision, aes, id0
        self.p = -3904
        self.meta = encrypt_metadata
        if revision == 6:
            self.key = hashlib.sha256(b"file key").digest()
            vs, ks = b"uvsalt01", b"ukslt001"
            self.U = hash6(user, vs, b"") + vs + ks
            self.UE = AES(hash6(user, ks, b"")).cbc(bytes(16), self.key)
            ovs, oks = b"ovsalt01", b"okslt001"
            self.O = hash6(owner, ovs, self.U) + ovs + oks
            self.OE = AES(hash6(owner, oks, self.U)).cbc(bytes(16), self.key)
            perms = struct.pack("<i", self.p) + b"\xff\xff\xff\xff" + (b"T" if encrypt_metadata else b"F") + b"adb" + b"0000"
            self.Perms = AES(self.key).block(perms)
            return
        n = 5 if revision == 2 else 16
        rounds = 50 if revision >= 3 else 0
        # Algorithm 3: O.
        h = hashlib.md5(padded(owner)).digest()
        for _ in range(rounds):
            h = hashlib.md5(h).digest()
        ok = h[:n]
        x = rc4(ok, padded(user))
        for i in range(1, 20 if revision >= 3 else 1):
            x = rc4(bytes(b ^ i for b in ok), x)
        self.O = x
        # Algorithm 2: the file key.
        m = padded(user) + self.O + struct.pack("<i", self.p) + id0
        if revision >= 4 and not encrypt_metadata:
            m += b"\xff\xff\xff\xff"
        h = hashlib.md5(m).digest()
        for _ in range(rounds):
            h = hashlib.md5(h[:n]).digest()
        self.key = h[:n]
        if revision == 2:
            # Algorithm 4: U.
            self.U = rc4(self.key, PAD)
            return
        # Algorithm 5: U.
        x = rc4(self.key, hashlib.md5(PAD + id0).digest())
        for i in range(1, 20):
            x = rc4(bytes(b ^ i for b in self.key), x)
        self.U = x + b"\x00" * 16

    def dictionary(self):
        if self.r == 6:
            return {"Filter": Name("Standard"), "V": 5, "R": 6, "Length": 256,
                    "CF": {"StdCF": {"CFM": Name("AESV3"), "AuthEvent": Name("DocOpen"), "Length": 32}},
                    "StmF": Name("StdCF"), "StrF": Name("StdCF"),
                    "O": self.O, "U": self.U, "OE": self.OE, "UE": self.UE, "Perms": self.Perms, "P": self.p}
        if self.r == 4:
            d = {"Filter": Name("Standard"), "V": 4, "R": 4, "Length": 128,
                 "CF": {"StdCF": {"CFM": Name("AESV2" if self.aes else "V2"), "AuthEvent": Name("DocOpen"), "Length": 16}},
                 "StmF": Name("StdCF"), "StrF": Name("StdCF"), "O": self.O, "U": self.U, "P": self.p}
            if not self.meta:
                d["EncryptMetadata"] = False
            return d
        if self.r == 2:
            return {"Filter": Name("Standard"), "V": 1, "R": 2, "O": self.O, "U": self.U, "P": self.p}
        return {"Filter": Name("Standard"), "V": 2, "R": 3, "Length": 128, "O": self.O, "U": self.U, "P": self.p}

    def crypt(self, n, data):
        """Algorithm 1."""
        if self.r == 6:
            return aes_pdf(self.key, data, b"%d" % n + data[:8])
        k = hashlib.md5(self.key + struct.pack("<I", n)[:3] + b"\x00\x00" + (b"sAlT" if self.aes else b"")).digest()[:min(len(self.key) + 5, 16)]
        return aes_pdf(k, data, b"%d" % n + data[:8]) if self.aes else rc4(k, data)


def hash6(pw, salt, udata):
    """ISO 32000-2 Algorithm 2.B."""
    k = hashlib.sha256(pw + salt + udata).digest()
    rnd = 0
    while True:
        k1 = (pw + k + udata) * 64
        e = AES(k[:16]).cbc(k[16:32], k1)
        k = [hashlib.sha256, hashlib.sha384, hashlib.sha512][sum(e[:16]) % 3](e).digest()
        rnd += 1
        if rnd >= 64 and e[-1] <= rnd - 32:
            break
    return k[:32]


# ---- Files ---------------------------------------------------------------------------------------------------------

def write_object(n, o, enc=None):
    sc = (lambda b: enc.crypt(n, b)) if enc else None
    if isinstance(o, Stream):
        data = o.data
        if enc and o.encrypt:
            data = enc.crypt(n, data)
        d = dict(o.d)
        if o.length_ref is not None:
            d["Length"] = Ref(o.length_ref)
        else:
            d["Length"] = o.length if o.length is not None else len(data)
        return b"%d 0 obj\n" % n + ser(d, sc) + b"\nstream\n" + data + b"\nendstream\nendobj\n"
    return b"%d 0 obj\n" % n + ser(o, sc) + b"\nendobj\n"


HEAD = b"%PDF-{v}\n%\xe2\xe3\xcf\xd3\n"


def classic(objs, root, version="1.4", trailer=None, enc=None, enc_num=None, startxref=None, head=None):
    out = head if head is not None else HEAD.replace(b"{v}", version.encode())
    offsets = {}
    for n in sorted(objs):
        offsets[n] = len(out)
        out += write_object(n, objs[n], None if n == enc_num else enc)
    size = max(objs) + 1
    xref = len(out)
    out += b"xref\n0 %d\n0000000000 65535 f \n" % size
    for n in range(1, size):
        out += (b"%010d 00000 n \n" % offsets[n]) if n in offsets else b"0000000000 00000 f \n"
    t = {"Size": size, "Root": Ref(root)}
    t.update(trailer or {})
    out += b"trailer\n" + ser(t) + b"\nstartxref\n%d\n%%%%EOF\n" % (startxref if startxref is not None else xref)
    return out, offsets, xref


def with_object_streams(objs, root, version="1.5", trailer=None, enc=None, enc_num=None, packed=None, hybrid=False):
    """Non-stream objects (or those in `packed`) in one object stream; a cross-reference stream (PNG Up predictor)
    or, when hybrid, a classic table plus /XRefStm."""
    packed = [n for n in sorted(objs) if not isinstance(objs[n], Stream) and n != enc_num] if packed is None else packed
    size = max(objs) + 3
    stm_num, xref_num = size - 2, size - 1
    header, body = b"", b""
    for n in packed:
        header += b"%d %d " % (n, len(body))
        body += ser(objs[n]) + b"\n"
    header += b"\n"
    stm = Stream({"Type": Name("ObjStm"), "N": len(packed), "First": len(header), "Filter": Name("FlateDecode")},
                 zlib.compress(header + body))
    out = HEAD.replace(b"{v}", version.encode())
    offsets = {}
    for n in sorted(objs):
        if n in packed:
            continue
        offsets[n] = len(out)
        out += write_object(n, objs[n], None if n == enc_num else enc)
    offsets[stm_num] = len(out)
    out += write_object(stm_num, stm, enc)
    xref_at = len(out)
    offsets[xref_num] = xref_at
    rows = b""
    for n in range(size):
        if n in packed:
            rows += struct.pack(">BIH", 2, stm_num, packed.index(n))
        elif n in offsets:
            rows += struct.pack(">BIH", 1, offsets[n], 0)
        else:
            rows += struct.pack(">BIH", 0, 0, 0)
    d = {"Type": Name("XRef"), "Size": size, "W": [1, 4, 2], "Filter": Name("FlateDecode"),
         "DecodeParms": {"Predictor": 12, "Columns": 7}}
    if not hybrid:
        d["Root"] = Ref(root)
        d.update(trailer or {})
    xs = Stream(d, zlib.compress(png_predict(rows, 7, [2])), encrypt=False)
    out += write_object(xref_num, xs)
    if not hybrid:
        return out + b"startxref\n%d\n%%%%EOF\n" % xref_at
    table = len(out)
    out += b"xref\n0 %d\n0000000000 65535 f \n" % size
    for n in range(1, size):
        out += (b"%010d 00000 n \n" % offsets[n]) if n in offsets else b"0000000000 00000 f \n"
    t = {"Size": size, "Root": Ref(root), "XRefStm": xref_at}
    t.update(trailer or {})
    return out + b"trailer\n" + ser(t) + b"\nstartxref\n%d\n%%%%EOF\n" % table


def save(name, data):
    with open(os.path.join(HERE, name), "wb") as f:
        f.write(data)


def page_doc(content, resources=None, page_extra=None, extra=None, media=(0, 0, 720, 720), filters=None):
    """1 catalog, 2 pages, 3 page, 4 contents, then `extra` objects from 5."""
    objs = {1: {"Type": Name("Catalog"), "Pages": Ref(2)},
            2: {"Type": Name("Pages"), "Kids": [Ref(3)], "Count": 1},
            3: {"Type": Name("Page"), "Parent": Ref(2), "MediaBox": list(media), "Contents": Ref(4),
                "Resources": resources or {}},
            4: Stream(filters or {}, content)}
    objs[3].update(page_extra or {})
    objs.update(extra or {})
    return objs


K = 0.5522847498307936


def fmt(*xs):
    return b" ".join(num(float(x)) for x in xs)


# ---- The fixtures --------------------------------------------------------------------------------------------------

def shapes():
    k36, k72 = 36 * 0.5523, 72 * K
    c = b"q\n1 0 0 1 0 0 cm\n"
    c += b"72 72 m 144 72 l S\n"
    c += b"72 144 m 108 180 144 108 180 144 c S\n"
    c += b"216 216 72 36 re f\n"
    # A circle of radius 36 at (360, 360): four quarters with the 0.5523 handle.
    c += fmt(396, 360) + b" m " + fmt(396, 360 + k36, 360 + k36, 396, 360, 396) + b" c "
    c += fmt(360 - k36, 396, 324, 360 + k36, 324, 360) + b" c " + fmt(324, 360 - k36, 360 - k36, 324, 360, 324) + b" c "
    c += fmt(360 + k36, 324, 396, 360 - k36, 396, 360) + b" c h S\n"
    # A quarter arc of radius 72 about (72, 576), from (144, 576) to (72, 648).
    c += fmt(144, 576) + b" m " + fmt(144, 576 + k72, 72 + k72, 648, 72, 648) + b" c S\n"
    # An ellipse 72 by 36 about (540, 180).
    a, b = 72 * K, 36 * K
    c += fmt(612, 180) + b" m " + fmt(612, 180 + b, 540 + a, 216, 540, 216) + b" c " + fmt(540 - a, 216, 468, 180 + b, 468, 180)
    c += b" c " + fmt(468, 180 - b, 540 - a, 144, 540, 144) + b" c " + fmt(540 + a, 144, 612, 180 - b, 612, 180) + b" c S\n"
    c += b"Q\n"
    # Nested transforms: (0, 324) -> (36, 324) -> doubled.
    c += b"q 2 0 0 2 0 0 cm q 1 0 0 1 36 0 cm 0 324 m 36 324 l S Q Q\n"
    # A clip that is not painted.
    c += b"0 0 m 720 720 l W n\n"
    c += b"BT /F1 12 Tf 72 700 Td (Hello) Tj ET\n"
    # Inline images: unfiltered data that contains " EI ", and an ASCIIHex one found by searching.
    c += b"BI /W 4 /H 2 /BPC 8 /CS /G ID \n EI Q\x00\xff" + b"\nEI\n"
    c += b"BI /W 2 /H 1 /CS /RGB /F /AHx ID 00FF00FF0000> EI\n"
    c += b"/Fm1 Do /Im1 Do /Sh1 sh\n"
    form = Stream({"Type": Name("XObject"), "Subtype": Name("Form"), "BBox": [0, 0, 720, 720], "Matrix": [1, 0, 0, 1, 0, 72]},
                  b"288 72 m 360 72 l S")
    image = Stream({"Type": Name("XObject"), "Subtype": Name("Image"), "Width": 1, "Height": 1, "ColorSpace": Name("DeviceGray"),
                    "BitsPerComponent": 8}, b"\x80")
    res = {"XObject": {"Fm1": Ref(5), "Im1": Ref(6)}, "Shading": {"Sh1": Ref(7)}}
    shading = {"ShadingType": 2, "ColorSpace": Name("DeviceGray"), "Coords": [0, 0, 1, 0],
               "Function": {"FunctionType": 2, "Domain": [0, 1], "N": 1}}
    objs = page_doc(c, res, extra={5: form, 6: image, 7: shading})
    save("shapes.pdf", classic(objs, 1)[0])


def rotated():
    objs = page_doc(b"72 72 m 144 72 l S", media=(0, 0, 612, 792),
                    page_extra={"CropBox": [36, 36, 576, 756], "Rotate": 90, "UserUnit": 2})
    save("rotated.pdf", classic(objs, 1, version="1.6")[0])


def pages():
    form = Stream({"Type": Name("XObject"), "Subtype": Name("Form"), "BBox": [0, 0, 100, 100]}, b"0 0 m 72 72 l S")
    objs = {1: {"Type": Name("Catalog"), "Pages": Ref(2)},
            2: {"Type": Name("Pages"), "Kids": [Ref(3), Ref(6)], "Count": 3, "MediaBox": [0, 0, 144, 144], "Rotate": 0},
            3: {"Type": Name("Pages"), "Parent": Ref(2), "Kids": [Ref(4), Ref(5)], "Count": 2,
                "Resources": {"XObject": {"F0": Ref(9)}}},
            4: {"Type": Name("Page"), "Parent": Ref(3), "Contents": [Ref(7), Ref(10)]},
            5: {"Type": Name("Page"), "Parent": Ref(3), "Contents": Ref(8)},
            6: {"Type": Name("Page"), "Parent": Ref(2), "Contents": Ref(8)},
            7: Stream({}, b"/F0"),
            8: Stream({}, b"0 0 m 1 1 l S"),
            9: form,
            10: Stream({}, b"Do")}
    save("pages.pdf", classic(objs, 1)[0])


def layers_objs(name_cut=b"Cut"):
    c = b"0 0 m 72 0 l S\n"
    c += b"/OC /oc1 BDC 0 72 m 72 72 l S EMC\n"
    c += b"/OC /oc2 BDC 0 144 m 72 144 l S EMC\n"
    c += b"/OC /oc3 BDC 0 216 m 72 216 l S EMC\n"
    c += b"/OC /oc4 BDC 0 252 m 72 252 l S EMC\n"
    c += b"/Fm Do\n"
    c += b"/Span << /ActualText (x) >> BDC 0 288 m 72 288 l S EMC\n"
    form = Stream({"Type": Name("XObject"), "Subtype": Name("Form"), "BBox": [0, 0, 720, 720], "OC": Ref(7)},
                  b"0 400 m 72 400 l S")
    res = {"Properties": {"oc1": Ref(5), "oc2": Ref(6), "oc3": Ref(7), "oc4": Ref(9)}, "XObject": {"Fm": Ref(8)}}
    extra = {5: {"Type": Name("OCG"), "Name": name_cut},
             6: {"Type": Name("OCG"), "Name": b"\xfe\xff" + "Ätch".encode("utf-16-be")},
             7: {"Type": Name("OCG"), "Name": b"Hidden"},
             8: form,
             9: {"Type": Name("OCMD"), "OCGs": [Ref(5), Ref(7)], "P": Name("AllOn")}}
    objs = page_doc(c, res, extra=extra)
    objs[1]["OCProperties"] = {"OCGs": [Ref(5), Ref(6), Ref(7)], "D": {"Order": [Ref(5), Ref(6), Ref(7)], "OFF": [Ref(7)]}}
    return objs


def layers():
    save("layers.pdf", classic(layers_objs(), 1, version="1.5")[0])


def xref_stream():
    content = b"".join(b"%d 0 m %d 72 l S\n" % (10 * i, 10 * i) for i in range(40))
    objs = page_doc(b"")
    objs[4] = Stream({"Filter": Name("FlateDecode"), "DecodeParms": {"Predictor": 15, "Columns": 16}},
                     zlib.compress(png_predict(content, 16, [0, 1, 2, 3, 4])))
    save("xrefstream.pdf", with_object_streams(objs, 1))


def hybrid():
    objs = page_doc(zlib.compress(b"0 0 m 72 0 l S"), filters={"Filter": Name("FlateDecode")})
    save("hybrid.pdf", with_object_streams(objs, 1, packed=[2, 3], hybrid=True))


def updated():
    objs = page_doc(b"0 0 m 72 0 l S")
    base, offsets, xref = classic(objs, 1)
    # An incremental update replaces the contents with a different line.
    add = write_object(4, Stream({}, b"0 0 m 144 0 l S"))
    at = len(base)
    out = base + add
    table = len(out)
    out += b"xref\n0 1\n0000000000 65535 f \n4 1\n%010d 00000 n \n" % at
    out += b"trailer\n" + ser({"Size": 5, "Root": Ref(1), "Prev": xref}) + b"\nstartxref\n%d\n%%%%EOF\n" % table
    save("updated.pdf", out)


def damaged():
    objs = page_doc(b"0 0 m 72 0 l S", extra={})
    objs[3]["Contents"] = [Ref(4), Ref(5), Ref(6)]
    objs[4] = Stream({}, b"0 0 m 72 0 l S", length=500)          # /Length past the data
    objs[5] = Stream({}, b"0 72 m 72 72 l S", length_ref=7)      # indirect /Length, correct
    objs[6] = Stream({}, b"0 144 m 72 144 l S", length=3)        # /Length short of the data
    objs[7] = len(b"0 72 m 72 72 l S")
    data, _, _ = classic(objs, 1, startxref=99999)
    save("damaged.pdf", data)


def filters():
    lines = b"".join(b"%d %d m %d %d l S\n" % (i * 7 % 500, i * 13 % 300, i * 11 % 400, i * 17 % 600) for i in range(300))
    objs = page_doc(b"")
    objs[3]["Contents"] = [Ref(4), Ref(5), Ref(6), Ref(7)]
    objs[4] = Stream({"Filter": Name("LZWDecode")}, lzw_encode(lines, early=1, clear_after=700))
    objs[5] = Stream({"Filter": Name("LZWDecode"), "DecodeParms": {"EarlyChange": 0}}, lzw_encode(b"0 700 m 72 700 l S\n" * 40, early=0))
    objs[6] = Stream({"Filter": [Name("ASCII85Decode"), Name("FlateDecode")]}, ascii85(zlib.compress(b"0 710 m 72 710 l S")))
    objs[7] = Stream({"Filter": [Name("ASCIIHexDecode"), Name("RunLengthDecode")]}, run_length(b"0 720 m    72 720 l S").hex().encode() + b">")
    save("filters.pdf", classic(objs, 1)[0])
    with open(os.path.join(HERE, "filters_lines.txt"), "wb") as f:
        f.write(lines)


ID0 = hashlib.md5(b"gc kernel pdf fixture").digest()


def encrypted(name, enc, objstm=False, user=None):
    objs = layers_objs()
    # A stream left in the clear by the Identity crypt filter (V4 and up).
    if enc.r >= 4:
        objs[3]["Contents"] = [Ref(4), Ref(10)]
        objs[10] = Stream({"Filter": Name("Crypt"), "DecodeParms": {"Name": Name("Identity")}}, b"0 500 m 72 500 l S", encrypt=False)
        objs[11] = Stream({"Type": Name("Metadata"), "Subtype": Name("XML")}, b"<x:xmpmeta/>", encrypt=enc.meta)
        objs[1]["Metadata"] = Ref(11)
    enc_num = max(objs) + 1
    objs[enc_num] = enc.dictionary()
    trailer = {"Encrypt": Ref(enc_num), "ID": [ID0, ID0]}
    if objstm:
        data = with_object_streams(objs, 1, version="1.7", trailer=trailer, enc=enc, enc_num=enc_num)
    else:
        data = classic(objs, 1, version="1.6", trailer=trailer, enc=enc, enc_num=enc_num)[0]
    save(name, data)


def illustrator():
    objs = page_doc(b"0 0 m 72 0 l 72 72 l h S")
    objs[1]["PieceInfo"] = {"Illustrator": {"LastModified": b"D:20261031", "Private": Ref(5)}}
    objs[5] = {"AIMetaData": Ref(6), "ContainerVersion": 11, "CreatorVersion": 28, "NumBlock": 1, "RoundtripVersion": 24}
    objs[6] = Stream({"Length1": 99}, bytes(range(256)) + b"%%EndProlog endstream junk")
    head = b"%PDF-1.6\n%\xe2\xe3\xcf\xd3\n%AI illustrator private data follows\n"
    save("illustrator.ai", classic(objs, 1, head=head)[0])


if __name__ == "__main__":
    shapes()
    rotated()
    pages()
    layers()
    xref_stream()
    hybrid()
    updated()
    damaged()
    filters()
    encrypted("rc4_40_r2.pdf", Standard(2, id0=ID0))
    encrypted("rc4_r3.pdf", Standard(3, id0=ID0))
    encrypted("aes128_r4.pdf", Standard(4, aes=True, encrypt_metadata=False, id0=ID0))
    encrypted("aes256_r6.pdf", Standard(6, id0=ID0), objstm=True)
    encrypted("password.pdf", Standard(3, user=b"secret", id0=ID0))
    illustrator()
    print("fixtures written to", HERE)
