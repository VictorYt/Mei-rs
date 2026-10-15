#![no_main]
//! Local alignment of two sequences cut from the input.

use libfuzzer_sys::fuzz_target;
use mei_rs::typer::align;

fuzz_target!(|data: &[u8]| {
    let data = &data[..data.len().min(600)];
    let (query, target) = data.split_at(data.len() / 2);
    if let Some(a) = align::local(query, target) {
        assert!(a.query_end <= query.len() && a.target_end <= target.len());
        assert!(a.matches <= a.columns);
    }
});
