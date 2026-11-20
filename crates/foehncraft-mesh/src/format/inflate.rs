//! A small zlib/deflate decoder (RFC 1950 and 1951) for compressed VTK data.
//!
//! Written from the RFCs, so that reading the files ParaView saves needs no
//! compression library. Only decompression is provided.

const LEN_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LEN_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];

struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
    bit: u32,
    nbits: u32,
}

impl Bits<'_> {
    fn need(&mut self, n: u32) -> Result<(), String> {
        while self.nbits < n {
            let b = *self
                .data
                .get(self.pos)
                .ok_or("the compressed data ends too early")?;
            self.pos += 1;
            self.bit |= (b as u32) << self.nbits;
            self.nbits += 8;
        }
        Ok(())
    }

    fn take(&mut self, n: u32) -> Result<u32, String> {
        if n == 0 {
            return Ok(0);
        }
        self.need(n)?;
        let v = self.bit & ((1u32 << n) - 1);
        self.bit >>= n;
        self.nbits -= n;
        Ok(v)
    }
}

/// A canonical Huffman code: the number of codes of each length and the
/// symbols in code order.
struct Huffman {
    count: [u16; 16],
    symbol: Vec<u16>,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Result<Self, String> {
        let mut count = [0u16; 16];
        for &l in lengths {
            count[l as usize] += 1;
        }
        count[0] = 0;
        let mut offs = [0u16; 16];
        for i in 1..16 {
            offs[i] = offs[i - 1] + count[i - 1];
        }
        let mut symbol = vec![0u16; lengths.len()];
        for (s, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbol[offs[l as usize] as usize] = s as u16;
                offs[l as usize] += 1;
            }
        }
        Ok(Self { count, symbol })
    }

    fn decode(&self, b: &mut Bits) -> Result<u16, String> {
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for len in 1..16 {
            code |= b.take(1)? as i32;
            let c = self.count[len] as i32;
            if code - c < first {
                return Ok(self.symbol[(index + (code - first)) as usize]);
            }
            index += c;
            first += c;
            first <<= 1;
            code <<= 1;
        }
        Err("a bad Huffman code in the compressed data".into())
    }
}

fn block(b: &mut Bits, out: &mut Vec<u8>, lit: &Huffman, dist: &Huffman) -> Result<(), String> {
    loop {
        let s = lit.decode(b)? as usize;
        match s {
            0..=255 => out.push(s as u8),
            256 => return Ok(()),
            257..=285 => {
                let i = s - 257;
                let len = LEN_BASE[i] as usize + b.take(LEN_EXTRA[i] as u32)? as usize;
                let d = dist.decode(b)? as usize;
                if d >= 30 {
                    return Err("a bad distance in the compressed data".into());
                }
                let dist = DIST_BASE[d] as usize + b.take(DIST_EXTRA[d] as u32)? as usize;
                if dist > out.len() {
                    return Err("a distance reaches before the start of the data".into());
                }
                for _ in 0..len {
                    out.push(out[out.len() - dist]);
                }
            }
            _ => return Err("a bad symbol in the compressed data".into()),
        }
    }
}

/// Decompresses raw deflate data.
pub fn inflate(data: &[u8]) -> Result<Vec<u8>, String> {
    let mut b = Bits {
        data,
        pos: 0,
        bit: 0,
        nbits: 0,
    };
    let mut out = Vec::new();
    loop {
        let last = b.take(1)?;
        match b.take(2)? {
            0 => {
                b.bit = 0;
                b.nbits = 0;
                let p = b.pos;
                let h = data
                    .get(p..p + 4)
                    .ok_or("the compressed data ends too early")?;
                let len = u16::from_le_bytes([h[0], h[1]]) as usize;
                let body = data
                    .get(p + 4..p + 4 + len)
                    .ok_or("the compressed data ends too early")?;
                out.extend_from_slice(body);
                b.pos = p + 4 + len;
            }
            1 => {
                let mut l = [8u8; 288];
                l[144..256].fill(9);
                l[256..280].fill(7);
                block(
                    &mut b,
                    &mut out,
                    &Huffman::new(&l)?,
                    &Huffman::new(&[5u8; 30])?,
                )?;
            }
            2 => {
                let nlen = b.take(5)? as usize + 257;
                let ndist = b.take(5)? as usize + 1;
                let ncode = b.take(4)? as usize + 4;
                const ORDER: [usize; 19] = [
                    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
                ];
                let mut cl = [0u8; 19];
                for &o in ORDER.iter().take(ncode) {
                    cl[o] = b.take(3)? as u8;
                }
                let clh = Huffman::new(&cl)?;
                let mut lengths = vec![0u8; nlen + ndist];
                let mut i = 0;
                while i < nlen + ndist {
                    let s = clh.decode(&mut b)?;
                    if s < 16 {
                        lengths[i] = s as u8;
                        i += 1;
                    } else {
                        let (v, rep) = match s {
                            16 => {
                                if i == 0 {
                                    return Err("a bad code length repeat".into());
                                }
                                (lengths[i - 1], 3 + b.take(2)? as usize)
                            }
                            17 => (0, 3 + b.take(3)? as usize),
                            _ => (0, 11 + b.take(7)? as usize),
                        };
                        if i + rep > nlen + ndist {
                            return Err("a bad code length repeat".into());
                        }
                        lengths[i..i + rep].fill(v);
                        i += rep;
                    }
                }
                block(
                    &mut b,
                    &mut out,
                    &Huffman::new(&lengths[..nlen])?,
                    &Huffman::new(&lengths[nlen..])?,
                )?;
            }
            _ => return Err("a bad block type in the compressed data".into()),
        }
        if last == 1 {
            return Ok(out);
        }
    }
}

/// Decompresses zlib data (a two-byte header, deflate, a checksum that is
/// not checked).
pub fn zlib(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.len() < 2 || data[0] & 0x0f != 8 {
        return Err("not zlib data".into());
    }
    inflate(&data[2..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stored_and_fixed_blocks() {
        // stored: final block, 5 bytes
        let s = [0x01, 5, 0, 0xfa, 0xff, b'h', b'e', b'l', b'l', b'o'];
        assert_eq!(inflate(&s).unwrap(), b"hello");
        // fixed Huffman, from zlib: "aaaaaaaaaa"
        let f = [
            0x78, 0x9c, 0x4b, 0x4c, 0x84, 0x01, 0x00, 0x14, 0xe1, 0x03, 0xcb,
        ];
        assert_eq!(zlib(&f).unwrap(), b"aaaaaaaaaa");
        assert!(inflate(&[0x07]).is_err());
    }
}
