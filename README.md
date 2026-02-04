# `lrfq`

CLI tool to generate summary statistics of long-read FASTQ files.

## Usage

```bash
# for uncompressed FASTQ:
lrfq < bc1015.fq
# for gzipped FASTQ:
gunzip -c bc1015.fq.gz | lrfq
```

## Output

```json
{
  "n": 2762,
  "n50": 1535,
  "median_length": 3693.0,
  "mean_length": 3944.025343953657,
  "longest_read": 24704,
  "shortest_read": 1088,
  "total_bases": 10893398,
  "gc_prop": 0.4792368735632353,
  "mean_read_qual": 90.03900824780611,
  "median_read_qual": 92.99046709870389
}
```

### Keys

* `n`: number of reads
* `n50`: read N50
* `median_length`: median read length
* `mean_length`: mean read length
* `longest_read`: longest read length
* `shortest_read`: shortest read length
* `total_bases`: sum of all read lengths
* `gc_prop`: GC content proportion
* `mean_read_qual`: mean PHRED read quality
* `median_read_qual`: median PHRED read quality
