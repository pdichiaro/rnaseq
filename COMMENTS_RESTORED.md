# Comments Restored to main.nf

## Summary

Important comments have been restored to the `workflows/rnaseq/main.nf` file to improve code documentation and maintainability. This document summarizes all the key comments added.

---

## 📍 Sections with Restored Comments

### 1. Module Imports Section (Lines ~1-50)

#### **DESeq2 Normalization Modules**
```groovy
//
// MODULE: Local modules for DESeq2 normalization with multiple methods
//
include { NORMALIZE_DESEQ2_QC_INVARIANT_GENES } from '../../modules/local/normalize_deseq2_qc_invariant_genes'
include { NORMALIZE_DESEQ2_QC_ALL_GENES      } from '../../modules/local/normalize_deseq2_qc_all_genes'
```

**Purpose**: Clarifies that these modules handle DESeq2 normalization with support for multiple normalization methods (all_genes and invariant_genes).

---

#### **DeepTools BigWig Normalization Module**
```groovy
// 
// MODULE: DeepTools BigWig normalization with DESeq2 scaling factors
//
include { DEEPTOOLS_BIGWIG_NORM              } from '../../modules/local/deeptools_bw_norm'
include { DEEPTOOLS_BIGWIG_NORM as DEEPTOOLS_BIGWIG_NORM_INVARIANT } from '../../modules/local/deeptools_bw_norm'
include { DEEPTOOLS_BIGWIG_NORM as DEEPTOOLS_BIGWIG_NORM_ALL_GENES } from '../../modules/local/deeptools_bw_norm'
```

**Purpose**: Documents the new feature that generates normalized BigWig files using DESeq2 scaling factors.

---

#### **Local Subworkflows**
```groovy
//
// SUBWORKFLOW: Local subworkflows
//
include { ALIGN_STAR                            } from '../../subworkflows/local/align_star'
```

**Purpose**: Separates and labels local subworkflow imports.

---

#### **nf-core Modules**
```groovy
//
// MODULE: Installed modules from nf-core/modules
//
include { DUPRADAR                   } from '../../modules/nf-core/dupradar'
```

**Purpose**: Distinguishes modules installed from nf-core/modules repository.

---

#### **nf-core Subworkflows**
```groovy
//
// SUBWORKFLOW: Installed subworkflows from nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
```

**Purpose**: Distinguishes subworkflows from nf-core repository.

---

### 2. RSEM Quantification Section (Lines ~340-360)

```groovy
    //
    // SUBWORKFLOW: Quantify gene and transcript abundance with RSEM
    //
    if (params.aligner == 'star' && !params.skip_quantification_method && quantification_methods.size() > 0) {
        
        if (quantification_methods.contains('rsem')) {
            //
            // QUANTIFICATION: STAR alignment + RSEM quantification
            //
            QUANTIFY_RSEM (
```

**Purpose**: Documents the RSEM quantification workflow for STAR-aligned reads.

---

### 3. DESeq2 Normalization Section (Lines ~360-400)

#### **Normalization Method Parsing**
```groovy
            if (!params.skip_qc & !params.skip_deseq2_qc) {
                //
                // DESeq2 normalization: Parse normalization methods (all_genes, invariant_genes, or both)
                //
                def normalization_methods = params.normalization_method instanceof List ? 
                    params.normalization_method : params.normalization_method.split(',').collect{it.trim()}
```

**Purpose**: Explains how the pipeline parses and handles multiple normalization methods.

---

#### **Invariant Genes Normalization**
```groovy
                //
                // MODULE: Invariant genes normalization (stable genes only)
                //
                if (normalization_methods.contains('invariant_genes')) {
                    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_ALIGNMENT (
```

**Purpose**: Clarifies that this branch handles normalization using only stable, non-differentially expressed genes.

---

#### **All Genes Normalization**
```groovy
                //
                // MODULE: All genes normalization (default DESeq2 method)
                //
                if (normalization_methods.contains('all_genes') || (!normalization_methods.contains('invariant_genes') && !normalization_methods.contains('all_genes'))) {
                    NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT (
```

**Purpose**: Documents that this is the default DESeq2 normalization method using all genes.

---

### 4. DeepTools BigWig Normalization Section (Lines ~1120-1200) - **MOST CRITICAL**

#### **Section Header**
```groovy
    //
    // MODULE: Generate normalized BigWig files using DESeq2 scaling factors
    // Uses mmrnaseq strategy: read scaling factors once, combine with BAM files, apply normalization
    //
    if (!params.skip_deeptools_norm && !params.skip_qc && !params.skip_deseq2_qc) {
        def normalization_methods = params.normalization_method instanceof List ? 
            params.normalization_method : params.normalization_method.split(',').collect{it.trim()}
        
        // Prepare BAM channel with BAI index files
        ch_bam_for_deeptools = ch_genome_bam
            .join(ch_genome_bam_index, by: [0])
```

**Purpose**: 
- Documents the new feature for generating normalized BigWig files
- Highlights the use of mmrnaseq strategy
- Explains the overall workflow

---

#### **Invariant Genes Channel Operations**
```groovy
        if (normalization_methods.contains('invariant_genes')) {
            //
            // CHANNEL OPERATION: Extract scaling factors from invariant genes normalization
            // Read individual *_scaling_factor.txt files, filter by 'invariant' directory
            // Extract sample name and scaling value as tuple: [sample_name, scaling_value]
            //
            ch_scaling_per_sample_invariant = ch_scaling_factors_individual
                .flatten()
                .filter { file ->
                    def parent_dir = file.getParent()?.getName() ?: ""
                    def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
                    parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
                }
                .map { file ->
                    def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
                    def scaling_value = file.text.trim()
                    [sample_name, scaling_value]
                }
            
            //
            // CHANNEL OPERATION: Combine BAM files with scaling factors (mmrnaseq strategy)
            // Use .combine() to create cartesian product, then filter by matching sample IDs
            // Result: [meta, bam, bai, scaling_value]
            //
            ch_combined_input_invariant = ch_bam_for_deeptools
                .combine(ch_scaling_per_sample_invariant)
                .map { meta, bam, bai, sample_id, scaling -> 
                    meta.id == sample_id ? [meta, bam, bai, scaling] : null
                }
                .filter { it != null }
```

**Purpose**: 
- Explains the two-step channel operation
- Documents file filtering by directory structure
- Describes the mmrnaseq combine strategy
- Shows the expected channel structure at each step

---

#### **All Genes Channel Operations**
```groovy
        if (normalization_methods.contains('all_genes') || (!normalization_methods.contains('invariant_genes') && !normalization_methods.contains('all_genes'))) {
            //
            // CHANNEL OPERATION: Extract scaling factors from all genes normalization
            // Read individual *_scaling_factor.txt files, exclude 'invariant' directory
            // Extract sample name and scaling value as tuple: [sample_name, scaling_value]
            //
            ch_scaling_per_sample_all_genes = ch_scaling_factors_individual
                .flatten()
                .filter { file ->
                    def parent_dir = file.getParent()?.getName() ?: ""
                    def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
                    !(parent_dir.contains('invariant') || grandparent_dir.contains('invariant'))
                }
                .map { file ->
                    def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
                    def scaling_value = file.text.trim()
                    [sample_name, scaling_value]
                }
            
            //
            // CHANNEL OPERATION: Combine BAM files with scaling factors (mmrnaseq strategy)
            // Use .combine() to create cartesian product, then filter by matching sample IDs
            // Result: [meta, bam, bai, scaling_value]
            //
            ch_combined_input_all_genes = ch_bam_for_deeptools
                .combine(ch_scaling_per_sample_all_genes)
                .map { meta, bam, bai, sample_id, scaling -> 
                    meta.id == sample_id ? [meta, bam, bai, scaling] : null
                }
                .filter { it != null }
```

**Purpose**: Same as invariant genes, but documents the filtering logic that **excludes** invariant directory files.

---

### 5. Software Versions Collection Section (Lines ~1200+)

```groovy
    //
    // MODULE: Collate software versions from all tools used in the pipeline
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }
```

**Purpose**: Documents the collection and formatting of software versions for MultiQC report.

---

### 6. MultiQC Section (Lines ~1210+)

```groovy
    //
    // MODULE: MultiQC - Aggregate all QC reports into a single HTML report
    //
    ch_multiqc_report = Channel.empty()

    if (!params.skip_multiqc) {
```

**Purpose**: Documents the MultiQC aggregation of all QC metrics.

---

## 🎯 Key Benefits of Restored Comments

### 1. **Code Navigation**
- Clear section headers make it easy to find specific functionality
- Comments act as landmarks in the 1,300+ line file

### 2. **Understanding New Features**
- DeepTools normalization is clearly documented as a new feature
- mmrnaseq strategy is explicitly mentioned
- Channel operations are explained step-by-step

### 3. **Maintainability**
- Future developers can understand the logic without deep analysis
- Channel transformations are documented with input/output formats
- Filtering logic is explained (e.g., why we check parent_dir and grandparent_dir)

### 4. **Training and Onboarding**
- New team members can understand the workflow structure quickly
- Comments explain **why** certain operations are done, not just **what**

### 5. **Debugging Support**
- Expected channel structures are documented
- Filtering logic is explained, making it easier to debug issues
- Strategy choices (mmrnaseq) are documented

---

## 📝 Comment Style Guidelines Used

### **Section Headers**
```groovy
//
// MODULE: Description of what this section does
//
```

### **Subsection Headers**
```groovy
//
// CHANNEL OPERATION: Specific operation description
// Additional details on multiple lines
//
```

### **Inline Comments**
```groovy
// Brief explanation of next line
ch_bam_for_deeptools = ch_genome_bam.join(ch_genome_bam_index, by: [0])
```

---

## 🔍 Most Critical Comments (Top 5)

### 1. **DeepTools BigWig Normalization Section Header**
- **Why critical**: Documents entirely new feature
- **Location**: Line ~1120
- **Impact**: Explains high-level strategy (mmrnaseq)

### 2. **Invariant Genes Channel Operations**
- **Why critical**: Complex channel logic with filtering
- **Location**: Line ~1135
- **Impact**: Documents file filtering and channel structure

### 3. **All Genes Channel Operations**
- **Why critical**: Similar complexity, opposite filtering logic
- **Location**: Line ~1165
- **Impact**: Clarifies the difference from invariant genes

### 4. **DESeq2 Normalization Method Parsing**
- **Why critical**: Explains multi-method support
- **Location**: Line ~362
- **Impact**: Documents how parameters are parsed

### 5. **Module Import Sections**
- **Why critical**: Organizes 30+ imports
- **Location**: Lines 1-100
- **Impact**: Makes navigation much easier

---

## 📊 Statistics

- **Total comments added**: ~15 major comment blocks
- **Lines of comments added**: ~50 lines
- **Sections documented**: 6 major sections
- **Critical sections**: 3 (DeepTools, DESeq2, Channel ops)

---

## ✅ Verification

To verify comments are in place:

```bash
# Check for DeepTools section comments
grep -A 5 "MODULE: Generate normalized BigWig" workflows/rnaseq/main.nf

# Check for channel operation comments
grep -A 3 "CHANNEL OPERATION:" workflows/rnaseq/main.nf

# Check for all section headers
grep "MODULE:\|SUBWORKFLOW:" workflows/rnaseq/main.nf
```

---

## 🎉 Result

The `main.nf` file now has comprehensive documentation that:
- ✅ Explains the purpose of each major section
- ✅ Documents complex channel operations
- ✅ Highlights new features (DeepTools normalization)
- ✅ References implementation strategies (mmrnaseq)
- ✅ Shows expected data structures
- ✅ Improves code maintainability

**The code is now production-ready with professional documentation!**
