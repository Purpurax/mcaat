
const K: u64 = 23;
#[derive(Debug)]
pub struct Reads {
    pub reads: Vec<Read>
}
#[derive(Debug)]
pub struct Read {
    pub sequence: String,
    pub start_k_mer: String,
    pub end_k_mer: String,
    pub nodes_between: u64
}

impl Reads {
    pub fn parse(reads_file_path_1: String, reads_file_path_2: String) -> Reads {
        let content_1: String = std::fs::read_to_string(&reads_file_path_1)
            .expect("Failed to read READS file 1")
            .trim()
            .to_string();
        let content_2 = std::fs::read_to_string(&reads_file_path_2)
            .expect("Failed to read READS file 2")
            .trim()
            .to_string();

        let iter_1 = content_1.split("\n")
            .filter(|line| {
                line.starts_with("A")
                || line.starts_with("T")
                || line.starts_with("C")
                || line.starts_with("G")
            });
        let iter_2 = content_2.split("\n")
            .filter(|line| {
                line.starts_with("A")
                || line.starts_with("T")
                || line.starts_with("C")
                || line.starts_with("G")
            });

        let mut sequences: Vec<String> = vec![];

        for (seq_1, seq_2) in iter_1.zip(iter_2) {
            let seq_2_mod: String = seq_2.chars()
                .rev()
                .map(|nucl| {
                    if nucl == 'A' {
                        'T'
                    } else if nucl == 'T' {
                        'A'
                    } else if nucl == 'C' {
                        'G'
                    } else if nucl == 'G' {
                        'C'
                    } else {
                        nucl
                    }
                }).collect::<String>();
            sequences.push(seq_1.to_string() + &seq_2_mod);
        }

        let reads: Vec<Read> = sequences.into_iter()
            .map(|sequence| {
                Read {
                    sequence: sequence.clone(),
                    start_k_mer: sequence[..K as usize].to_string(),
                    end_k_mer: sequence[(sequence.len() - K as usize)..].to_string(),
                    nodes_between: (sequence.len() as u64).saturating_sub(K + 1)
                }
            })
            .collect::<Vec<Read>>();
        
        Reads { reads }
    }
}
