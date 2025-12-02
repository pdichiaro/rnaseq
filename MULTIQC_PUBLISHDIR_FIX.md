# MultiQC PublishDir Fix

## Issue Summary
MultiQC reports for aligner-specific variants (MULTIQC_STAR, MULTIQC_HISAT2, MULTIQC_KALLISTO) were not being published to the correct subfolders. Instead, all reports were being published to `${params.outdir}/multiqc` due to the default publishDir logic.

## Root Cause
The default publishDir configuration uses:
```groovy
path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" }
```

This pattern extracts the first word of the process name:
- `MULTIQC_STAR` → `multiqc`
- `MULTIQC_HISAT2` → `multiqc`
- `MULTIQC_KALLISTO` → `multiqc`

All three processes resolve to the same base folder `multiqc`, causing reports to overwrite each other or be placed in the same location.

## Solution
Added process-specific publishDir overrides in `nextflow.config` to ensure each MultiQC variant publishes to its own subfolder:

```groovy
// MultiQC aligner-specific reports
withName: '.*:MULTIQC_STAR' {
    publishDir = [
        path: { "${params.outdir}/multiqc/star" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}

withName: '.*:MULTIQC_HISAT2' {
    publishDir = [
        path: { "${params.outdir}/multiqc/hisat2" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}

withName: '.*:MULTIQC_KALLISTO' {
    publishDir = [
        path: { "${params.outdir}/multiqc/kallisto" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

## Expected Output Structure
After this fix, MultiQC reports will be organized as follows:

```
${params.outdir}/
├── multiqc/
│   ├── star/
│   │   └── multiqc_star_report.html
│   ├── hisat2/
│   │   └── multiqc_hisat2_report.html
│   ├── kallisto/
│   │   └── multiqc_kallisto_report.html
│   └── multiqc_report.html (main combined report)
```

## Testing
A standalone test script was created and validated the publishDir configuration:

```bash
nextflow run test_multiqc_publish.nf
```

Results:
✅ `test_multiqc_output/multiqc/star/multiqc_star_report.html`
✅ `test_multiqc_output/multiqc/hisat2/multiqc_hisat2_report.html`
✅ `test_multiqc_output/multiqc/kallisto/multiqc_kallisto_report.html`

## Files Modified
- `nextflow.config`: Added 3 process-specific publishDir configurations

## Related Configuration
The following configurations in `conf/base.config` remain unchanged and work in conjunction with this fix:

```groovy
withName: '.*:MULTIQC_STAR' {
    ext.prefix = 'multiqc_star'
    errorStrategy = 'ignore'
}
withName: '.*:MULTIQC_HISAT2' {
    ext.prefix = 'multiqc_hisat2'
    errorStrategy = 'ignore'
}
withName: '.*:MULTIQC_KALLISTO' {
    ext.prefix = 'multiqc_kallisto'
    errorStrategy = 'ignore'
}
```

The `ext.prefix` ensures proper naming, while the new publishDir directives ensure proper folder organization.

## Notes
- The main `MULTIQC` process (combined report) continues to use the default publishDir: `${params.outdir}/multiqc`
- The `errorStrategy = 'ignore'` ensures the pipeline continues even if aligner-specific MultiQC reports fail
- The `saveAs` clause prevents `versions.yml` files from being published to reduce clutter
