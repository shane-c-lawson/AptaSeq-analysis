# AptaSeq Analysis

This pipeline treats each unique 40mer down to the base as a separate aptamer. An alternative approach I tried was to treat sequences with 90%+ identity, or <= 4 mismatches, as the same aptamer in order to increase the tallied number of hits for that group and better deconvolute which varieties perform the best across the entire pool. Unfortunately, with hundreds of thousands of unique aptamers present across the reads, calculating the Hamming distance between every single possible combination in order to create those groupings leads to a ridiculously high programming complexity that I estimate could take weeks to process. It could be done, and theoritically in a situation where we aim to identify ideal aptamers for potential therapeutics, it might be ideal. If you want to look at how I approached the analysis, it is left in the script though bypassed and commented out.

Because here we are more interested in the performance of specific aptamers (namely the controls) in each condition, similar sequences can be treated separately. Interestingly, looking at the results, the top hits do not actually appear to have much similarity. However, that does not factor in any aptamers that might share high sequence identity and thus increase their score, which could aid with potential improper basecalling. (I've neglected that possibility in general because it can't really be predicted without an alignment)

Judging from the results, the best combination is unsurprisingly 50k cells with 100 copies of aptamer per cell, the worst being any 10k cells samples or any 5 copies/cell. The diversity of aptamers statistic is insignificant. Though 100k/50 had some of the lowest total reads (see indexing-qc), it still produced the 2 positive controls as top 2 hits. 50k/100 had only marginally more total reads, but the most control hits. Incidentally, those samples were treated with the same level amplified aptamer pool. I would personally choose the 50k/100, though 100k might be as good.
