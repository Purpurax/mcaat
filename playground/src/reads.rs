use std::collections::HashSet;

use crate::graph::Graph;


const K: usize = 23;
#[derive(Clone, Debug)]
pub struct Reads {
    pub reads: Vec<Read>
}
#[derive(Clone, Debug)]
pub struct Read {
    pub sequence: String,
    pub start_k_mer: String,
    pub end_k_mer: String,
    pub nodes_between: usize
}

impl Reads {
    pub fn new() -> Reads {
        Reads {
            reads: vec![]
        }
    }

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
                line.chars()
                    .all(|char| char == 'A'
                        || char == 'T'
                        || char == 'C'
                        || char == 'G')
            });
        let iter_2 = content_2.split("\n")
            .filter(|line| {
                line.chars()
                    .all(|char| char == 'A'
                        || char == 'T'
                        || char == 'C'
                        || char == 'G')
            });

        let sequences: Vec<String> = iter_1.map(|s| s.to_string())
            .chain(
                iter_2.map(|seq| {
                    seq.chars()
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
                        }).collect::<String>()
                })
            ).collect::<Vec<String>>();

        let reads: Vec<Read> = sequences.into_iter()
            .map(|sequence| {
                Read {
                    sequence: sequence.clone(),
                    start_k_mer: sequence[..K as usize].to_string(),
                    end_k_mer: sequence[(sequence.len() - K as usize)..].to_string(),
                    nodes_between: (sequence.len()).saturating_sub(K + 1)
                }
            })
            .collect::<Vec<Read>>();
        
        Reads { reads }
    }

    pub fn export_as_desired_input(&self, file_path: Option<String>) {
        let header: String = "read_start_k_mer,read_end_k_mer,nodes_between\n".to_string();
        
        let content: String = self.reads.clone().into_iter()
            .map(|read| {
                format!("{},{},{}", read.start_k_mer, read.end_k_mer, read.nodes_between)
            }).fold(String::new(), |a, b| a + "\n" + &b)
            .trim()
            .to_string();

        let _ = std::fs::write(
            file_path.unwrap_or(
                "/home/master/Documents/UNI/Informatik/Semester-4/Bachelor/mcaat/data/desired_inputs/reads.csv".to_string()
            ), header + &content
        );
    }

    pub fn get_relevant(&self, graph: &Graph) -> Reads {
        let unique_nodes_with_seqs: HashSet<String> = graph.get_unique_node_seq();

        let mut relevant_reads = Reads::new();
        for read in self.reads.clone().into_iter() {
            if unique_nodes_with_seqs.contains(&read.start_k_mer)
            && unique_nodes_with_seqs.contains(&read.end_k_mer) {
                relevant_reads.reads.push(read);
            }
        }

        relevant_reads
    }
}
