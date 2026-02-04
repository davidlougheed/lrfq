use noodles_fastq as fastq;
use serde::Serialize;
// use std::fs::File;
use std::io::{BufReader, Error, Stdin, stdin};

struct FastqStatsAgg {
    read_lengths: Vec<usize>,
    read_quals: Vec<f64>,
    total_bases: usize,
    total_gc: usize,
}

impl FastqStatsAgg {
    fn update(&mut self, seq: &[u8], quals: &[u8]) {
        let seq_len = seq.len();

        // TODO: optimize by not storing every read length (use vec histogram for hifi?) / keeping running averages

        self.read_lengths.push(seq_len);
        self.read_quals.push(quals.iter().map(|&q| (q - 33) as f64).sum::<f64>() / seq_len as f64);
        self.total_bases += seq_len;
        self.total_gc += seq.iter().filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c').count();

        // TODO
    }
}

#[derive(Serialize)]
struct FastqStatsReport {
    n: usize,
    // read lengths
    n50: usize,
    median_length: f64,
    mean_length: f64,
    longest_read: usize,
    shortest_read: usize,
    // bases
    total_bases: usize,
    gc_prop: f64,
    // qualities
    mean_read_qual: f64,
    median_read_qual: f64,
}

fn n50(x: &Vec<usize>, total_bases: usize) -> usize {
    let half_total_bases = total_bases / 2;
    // x: already sorted smallest-to-largest; start at the largest end
    let mut sum: usize = 0;
    for i in(0..x.len()).rev() {
        sum += x[i];
        if sum >= half_total_bases {
            return x[i];
        }
    }
    0  // unreachable
}

fn median_usize(x: &Vec<usize>) -> f64 {
    if x.len() % 2 == 1 {
        x[x.len() / 2] as f64
    } else {
        ((x[x.len() / 2] + x[(x.len() / 2) + 1]) as f64) / 2.0f64
    }
}

fn median_f64(x: &Vec<f64>) -> f64 {
    if x.len() % 2 == 1 {
        x[x.len() / 2]
    } else {
        (x[x.len() / 2] + x[(x.len() / 2) + 1]) / 2.0f64
    }
}

fn fastq_report() -> Result<FastqStatsReport, Error> {
    let mut stats_agg = FastqStatsAgg {
        read_lengths: Vec::new(),
        read_quals: Vec::new(),
        total_bases: 0,
        total_gc: 0,
    };

    let mut reader: fastq::io::Reader<BufReader<Stdin>> = fastq::io::Reader::new(BufReader::new(stdin()));
    let mut record = fastq::Record::default();

    loop {
        let bytes_read: usize = reader.read_record(&mut record)?;

        if bytes_read == 0 {
            // EOF
            break;
        }

        let seq = record.sequence();
        let quals = record.quality_scores();

        stats_agg.update(seq, quals);
    }

    stats_agg.read_lengths.sort_unstable();

    // ------------------------------------------------------------------------

    let n_reads = stats_agg.read_lengths.len();
    let total_bases = stats_agg.total_bases;

    // TODO: histogram

    Ok(FastqStatsReport {
        n: n_reads,
        n50: n50(&stats_agg.read_lengths, total_bases),
        median_length: median_usize(&stats_agg.read_lengths),
        mean_length: total_bases as f64 / n_reads as f64,
        longest_read: stats_agg.read_lengths[n_reads - 1],
        shortest_read: stats_agg.read_lengths[0],
        total_bases,
        gc_prop: stats_agg.total_gc as f64 / total_bases as f64,
        mean_read_qual: stats_agg.read_quals.iter().sum::<f64>() / n_reads as f64,
        median_read_qual: median_f64(&stats_agg.read_quals)
    })
}

fn main() -> Result<(), Error> {
    let report = fastq_report()?;

    println!("{}", serde_json::to_string_pretty(&report)?);

    Ok(())
}
