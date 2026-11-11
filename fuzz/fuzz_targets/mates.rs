#![no_main]
//! The mate logic on candidates built from the input (see `common`):
//! sampling, typing by the mates, tail mates, then the filters and the score.

mod common;

use libfuzzer_sys::fuzz_target;
use mei_rs::genotype::score;
use mei_rs::typer::mates::MAX_CHECKED_MATES;

fuzz_target!(|data: &[u8]| {
    let set = common::consensus();
    let candidates = common::candidates(&mut common::Input(data), set);
    for c in &candidates {
        assert!(c.mate_checks.len() <= MAX_CHECKED_MATES);
        assert!(c.mate_checks.iter().all(|m| m.signal < c.signals.len()));
        assert!(c.mei_mates <= c.checked_mates);
        assert!(c.tail_mates <= c.checked_mates);
        assert!(c.mate_identity.is_none_or(|i| (0.0..=1.0).contains(&i)));
        assert_eq!(c.mate_identity.is_some(), c.mei_mates > 0);
        assert!(c.qual.is_some_and(|q| q <= score::max_quality()));
        if let Some(t) = &c.typing {
            assert!(t.element_start <= t.element_end);
            assert!((0.0..=1.0).contains(&t.identity));
        }
        if let Some(g) = &c.genotype {
            assert!(g.log_likelihoods.iter().all(|ll| ll.is_finite()));
        }
    }
});
