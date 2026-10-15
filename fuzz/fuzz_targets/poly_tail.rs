#![no_main]
//! Poly(A/T) tail detection: the first bytes set the options, the rest is
//! the clip.

use libfuzzer_sys::fuzz_target;
use mei_rs::reader::signals::{ClipOptions, poly_tail};

fuzz_target!(|data: &[u8]| {
    let Some((&first, clip)) = data.split_first() else {
        return;
    };
    let options = ClipOptions {
        polya_min_length: u32::from(first & 0x1f),
        polya_min_fraction: f64::from(first >> 5) / 7.0,
        ..ClipOptions::default()
    };
    for junction_at_end in [false, true] {
        if let Some(tail) = poly_tail(clip, junction_at_end, &options) {
            assert!(tail.length as usize <= clip.len());
            assert!(tail.length >= options.polya_min_length);
        }
    }
});
