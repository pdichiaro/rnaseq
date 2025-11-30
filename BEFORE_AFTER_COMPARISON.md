# Before & After Comparison: Channel Refactoring

## The Problem We Solved

### BEFORE: Mix Then Filter (Fragile) ❌

```
┌──────────────────────────────────────────────────────┐
│  NORMALIZATION PROCESSES                              │
│                                                       │
│  ┌─────────────────┐     ┌─────────────────┐        │
│  │ Invariant Genes │     │   All Genes     │        │
│  │                 │     │                 │        │
│  │  scaling_factors│     │  scaling_factors│        │
│  └────────┬────────┘     └────────┬────────┘        │
│           │                       │                  │
│           └───────────┬───────────┘                  │
│                       ↓                              │
│          ┌────────────────────────┐                  │
│          │  SINGLE MIXED CHANNEL  │                  │
│          │                        │                  │
│          │ ch_scaling_factors_    │                  │
│          │     individual         │                  │
│          │                        │                  │
│          │ Contains ALL files     │                  │
│          │ from both methods!     │                  │
│          └───────────┬────────────┘                  │
└──────────────────────┼───────────────────────────────┘
                       │
                       │  ⚠️  PROBLEM: Mixed files!
                       │
┌──────────────────────┼───────────────────────────────┐
│  DEEPTOOLS BIGWIG GENERATION                         │
│                      │                               │
│                      ↓                               │
│         ┌────────────────────────┐                   │
│         │   .filter { file ->    │                   │
│         │  Check directory names │                   │
│         │  parent_dir.contains   │                   │
│         │    ('invariant')       │                   │
│         └─────┬────────┬─────────┘                   │
│               │        │                             │
│     ┌─────────┘        └─────────┐                   │
│     ↓                            ↓                   │
│  INVARIANT                    ALL GENES              │
│  BigWigs                      BigWigs                │
│                                                      │
│  ⚠️  ISSUES:                                         │
│  • Fragile: Depends on directory naming             │
│  • Error-prone: Could mix files                     │
│  • Hard to debug: Complex filtering logic           │
│  • Performance: Unnecessary filtering overhead      │
└──────────────────────────────────────────────────────┘
```

### AFTER: Separate at Source (Robust) ✅

```
┌──────────────────────────────────────────────────────┐
│  NORMALIZATION PROCESSES                              │
│                                                       │
│  ┌─────────────────┐     ┌─────────────────┐        │
│  │ Invariant Genes │     │   All Genes     │        │
│  │                 │     │                 │        │
│  │  scaling_factors│     │  scaling_factors│        │
│  └────────┬────────┘     └────────┬────────┘        │
│           │                       │                  │
│           │    ✓ Separate at      │                  │
│           │      source!          │                  │
│           │                       │                  │
│      ┌────┼────┐             ┌────┼────┐            │
│      ↓    │    ↓             ↓    │    ↓            │
│   ┌────┐ │ ┌────────┐   ┌────┐ │ ┌────────┐       │
│   │Mix │ │ │Invariant│  │Mix │ │ │AllGenes│       │
│   │Ch. │ │ │ Chan.  │   │Ch. │ │ │ Chan.  │       │
│   └────┘ │ └────────┘   └────┘ │ └────────┘       │
│          │                      │                   │
│   For    │ For DeepTools        │ For DeepTools    │
│   other  │ (Invariant)          │ (All Genes)      │
│   uses   │                      │                   │
└──────────┼──────────────────────┼───────────────────┘
           │                      │
           │ ✓ Clean separation   │
           │                      │
┌──────────┼──────────────────────┼───────────────────┐
│  DEEPTOOLS BIGWIG GENERATION    │                   │
│           │                     │                   │
│           ↓                     ↓                   │
│    ┌─────────────┐       ┌─────────────┐           │
│    │  INVARIANT  │       │  ALL GENES  │           │
│    │   Channel   │       │   Channel   │           │
│    │             │       │             │           │
│    │  .flatten() │       │  .flatten() │           │
│    │  .map { }   │       │  .map { }   │           │
│    │             │       │             │           │
│    │ ✓ Direct    │       │ ✓ Direct    │           │
│    │   use!      │       │   use!      │           │
│    └──────┬──────┘       └──────┬──────┘           │
│           ↓                     ↓                   │
│     INVARIANT                ALL GENES              │
│     BigWigs                  BigWigs                │
│                                                     │
│  ✅ BENEFITS:                                       │
│  • Reliable: No directory dependencies             │
│  • Safe: No cross-contamination possible           │
│  • Clear: Obvious data provenance                  │
│  • Fast: No filtering overhead                     │
└─────────────────────────────────────────────────────┘
```

## Code Comparison

### Channel Creation

#### BEFORE:
```groovy
// Only one channel for everything
ch_scaling_factors_individual = Channel.empty()

// Later, try to separate by filtering...
ch_scaling_per_sample_invariant = ch_scaling_factors_individual
    .flatten()
    .filter { file ->
        // ⚠️ Fragile! Depends on directory structure
        def parent_dir = file.getParent()?.getName() ?: ""
        def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
        parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
    }
    .map { file -> ... }
```

#### AFTER:
```groovy
// Three channels for clear separation
ch_scaling_factors_individual = Channel.empty()           // Mixed (compatibility)
ch_scaling_factors_individual_invariant = Channel.empty() // Dedicated
ch_scaling_factors_individual_all_genes = Channel.empty() // Dedicated

// Later, use directly with no filtering needed
ch_scaling_per_sample_invariant = ch_scaling_factors_individual_invariant
    .flatten()
    // ✅ Direct use - clean and simple!
    .map { file -> ... }
```

### Population

#### BEFORE:
```groovy
// Everything goes into one channel
if (normalization_methods.contains('invariant_genes')) {
    ch_scaling_factors_individual = ch_scaling_factors_individual.mix(
        NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.scaling_factors_individual
    )
}

if (normalization_methods.contains('all_genes')) {
    ch_scaling_factors_individual = ch_scaling_factors_individual.mix(
        NORMALIZE_DESEQ2_QC_ALL_GENES.out.scaling_factors_individual
    )
}

// ⚠️ Now they're all mixed together!
```

#### AFTER:
```groovy
// Separate at the source
if (normalization_methods.contains('invariant_genes')) {
    // Add to both mixed and dedicated channels
    ch_scaling_factors_individual = ch_scaling_factors_individual.mix(
        NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.scaling_factors_individual
    )
    ch_scaling_factors_individual_invariant = ch_scaling_factors_individual_invariant.mix(
        NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.scaling_factors_individual
    )
}

if (normalization_methods.contains('all_genes')) {
    // Add to both mixed and dedicated channels
    ch_scaling_factors_individual = ch_scaling_factors_individual.mix(
        NORMALIZE_DESEQ2_QC_ALL_GENES.out.scaling_factors_individual
    )
    ch_scaling_factors_individual_all_genes = ch_scaling_factors_individual_all_genes.mix(
        NORMALIZE_DESEQ2_QC_ALL_GENES.out.scaling_factors_individual
    )
}

// ✅ Clean separation maintained throughout
```

## Failure Scenarios

### BEFORE: What Could Go Wrong ❌

1. **Directory Structure Changes**
   ```groovy
   // If someone renames 'invariant' to 'invariant_genes_method'
   parent_dir.contains('invariant') // ⚠️ Still matches! False positive!
   ```

2. **Nested Directory Confusion**
   ```groovy
   // If paths become:
   // /results/all_genes_method/sample1_invariant/file.txt
   parent_dir.contains('invariant') // ⚠️ False positive!
   ```

3. **Multiple Filtering Points**
   ```groovy
   // Filtering logic duplicated and could diverge
   .filter { parent_dir.contains('invariant') || grandparent_dir.contains('invariant') }
   // vs
   .filter { !parent_dir.contains('invariant') && !grandparent_dir.contains('invariant') }
   // ⚠️ Logic errors possible!
   ```

### AFTER: Failure-Proof ✅

1. **No Directory Dependencies**
   ```groovy
   // Files go to correct channel at creation time
   // No chance of misrouting later
   ch_scaling_factors_individual_invariant.mix(...)
   ```

2. **Clear Data Flow**
   ```groovy
   // You can trace exactly where each file comes from
   Source Process → Dedicated Channel → Consumption
   ```

3. **Type Safety**
   ```groovy
   // Each channel has a clear, single purpose
   // No ambiguity about what files it contains
   ```

## Performance Comparison

### BEFORE: Inefficient
```
1. Mix all files together           (N files → 1 channel)
2. Flatten channel                  (1 channel → N items)
3. Filter by directory              (N items → scan each)
4. Extract sample info              (Filtered items → tuples)
5. Combine with BAM files           (Tuples → joined)

Total operations: Mix + Flatten + Filter + Map + Combine
```

### AFTER: Efficient
```
1. Separate at source               (N files → correct channel)
2. Flatten channel                  (1 channel → N items)
3. Extract sample info              (N items → tuples, no filtering!)
4. Combine with BAM files           (Tuples → joined)

Total operations: Mix + Flatten + Map + Combine
Savings: One Filter operation eliminated per run
```

## Maintainability Comparison

### BEFORE: Hard to Maintain ❌

```groovy
// To add a new normalization method:
// 1. Add files to mixed channel ✓
// 2. Update filtering logic in MULTIPLE places ⚠️
// 3. Ensure filtering doesn't conflict with existing ⚠️
// 4. Test all combinations ⚠️

.filter { file ->
    def parent_dir = file.getParent()?.getName() ?: ""
    // Need to add more conditions here...
    parent_dir.contains('invariant') || 
    parent_dir.contains('new_method') ||  // ⚠️ Getting complex!
    grandparent_dir.contains('invariant') ||
    grandparent_dir.contains('new_method')
}
```

### AFTER: Easy to Extend ✅

```groovy
// To add a new normalization method:
// 1. Create new dedicated channel
ch_scaling_factors_individual_new_method = Channel.empty()

// 2. Populate at source (one line per quantifier)
ch_scaling_factors_individual_new_method = 
    ch_scaling_factors_individual_new_method.mix(...)

// 3. Use directly (one section)
ch_scaling_per_sample_new_method = 
    ch_scaling_factors_individual_new_method.flatten().map { ... }

// ✅ Clean, predictable, no conflicts possible!
```

## Testing Verification

### Test Coverage

| Test Case | Before | After |
|-----------|--------|-------|
| Single method (invariant only) | ⚠️ Filtering needed | ✅ Direct use |
| Single method (all_genes only) | ⚠️ Filtering needed | ✅ Direct use |
| Both methods simultaneously | ⚠️ Complex filtering | ✅ Clean separation |
| Multiple quantifiers (STAR+RSEM) | ⚠️ Path detection fragile | ✅ Path detection + separation |
| Directory rename | ❌ Breaks | ✅ Still works |
| New normalization method | ⚠️ Update filters | ✅ Add channel |

### Verification Commands

```bash
# Before: Check if filtering logic is correct
grep -A 10 "filter { file ->" workflows/rnaseq/main.nf
# (Complex, multiple occurrences)

# After: Verify clean separation
./verify_channel_refactoring.sh
# ✅ ALL CHECKS PASSED!
```

## Conclusion

The refactoring transforms:
- **Fragile directory-based filtering** → **Robust source-level separation**
- **Complex filtering logic** → **Simple direct channel use**
- **Hard to maintain** → **Easy to extend**
- **Risk of cross-contamination** → **Guaranteed separation**

**Result**: More reliable, maintainable, and performant pipeline! 🎉
