use std::collections::HashSet;

use crate::graph::Graph;
use fastq::{Parser, Record};


pub const K: usize = 23;
#[derive(Clone, Debug)]
pub struct Reads {
    pub reads: Vec<Read>
}
#[derive(Clone, Debug)]
pub struct Read {
    pub sequence: String,
    pub start_k_mer: String,
    pub middle_k_mer: String,
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
        let mut sequences = Vec::new();
        
        let parser = Parser::new(std::fs::File::open(&reads_file_path_1)
            .expect("Failed to open READS file 1"));
        parser.each(|record| {
            sequences.push(std::str::from_utf8(record.seq()).unwrap().to_string());
            true
        }).expect("Failed to parse READS file 1");
        
        let parser = Parser::new(std::fs::File::open(&reads_file_path_2)
            .expect("Failed to open READS file 2"));
        parser.each(|record| {
            sequences.push(std::str::from_utf8(record.seq()).unwrap().to_string());
            true
        }).expect("Failed to parse READS file 2");

        let reads: Vec<Read> = sequences.into_iter().map(|sequence| {
            Read {
                sequence: sequence.clone(),
                start_k_mer: sequence[..K].to_string(),
                middle_k_mer: sequence[(sequence.len() - K) / 2..(sequence.len() + K) / 2].to_string(),
                end_k_mer: sequence[(sequence.len() - K)..].to_string(),
                nodes_between: sequence.len().saturating_sub(K + 1)
            }
        }).collect();
        
        Reads { reads }
    }

    pub fn export_as_desired_input(&self, file_path: String) {
        let header: String = "read_start_k_mer,read_middle_k_mer,read_end_k_mer,nodes_between\n".to_string();
        
        let content: String = self.reads.clone().into_iter()
            .map(|read| {
                format!("{},{},{},{}", read.start_k_mer, read.middle_k_mer, read.end_k_mer, read.nodes_between)
            }).fold(String::new(), |a, b| a + "\n" + &b)
            .trim()
            .to_string();

        let _ = std::fs::write(file_path, header + &content);
    }

    pub fn get_relevant(&self, graph: &Graph) -> Reads {
        let unique_nodes_with_seqs: HashSet<String> = graph.get_unique_node_seq();

        let mut relevant_reads = Reads::new();
        for read in self.reads.clone().into_iter() {
            if unique_nodes_with_seqs.contains(&read.start_k_mer)
            || unique_nodes_with_seqs.contains(&read.middle_k_mer)
            || unique_nodes_with_seqs.contains(&read.end_k_mer) {
                relevant_reads.reads.push(read);
            }
        }

        relevant_reads
    }
}
