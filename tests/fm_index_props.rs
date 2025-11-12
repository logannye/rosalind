use proptest::prelude::*;
use rosalind::genomics::{BaseCode, BlockedFMIndex, FmSymbol};

fn fm_symbols() -> [FmSymbol; 6] {
    [
        FmSymbol::Sentinel,
        FmSymbol::Base(BaseCode::A),
        FmSymbol::Base(BaseCode::C),
        FmSymbol::Base(BaseCode::G),
        FmSymbol::Base(BaseCode::T),
        FmSymbol::Base(BaseCode::N),
    ]
}

proptest! {
    #[test]
    fn rank_totals_are_consistent(
        reference in proptest::collection::vec(prop_oneof![
            Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')
        ], 1..64),
        block_size in 1usize..32,
    ) {
        let index = BlockedFMIndex::build(&reference, block_size).expect("index build succeeds");

        for symbol in fm_symbols() {
            let total = index.total(symbol);
            let rank_at_end = index.rank(symbol, index.len());
            prop_assert_eq!(total, rank_at_end, "rank at end should equal total");

            let mut previous = 0;
            for pos in 0..=index.len() {
                let rank = index.rank(symbol, pos);
                prop_assert!(rank >= previous, "rank must be monotonic");
                previous = rank;
            }
        }

        let base_total: u32 = [BaseCode::A, BaseCode::C, BaseCode::G, BaseCode::T, BaseCode::N]
            .into_iter()
            .map(|code| index.total(FmSymbol::Base(code)))
            .sum();
        let sentinel_total = index.total(FmSymbol::Sentinel);
        prop_assert_eq!(base_total + sentinel_total, index.len() as u32, "totals should sum to BWT length");
    }
}
