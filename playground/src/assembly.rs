use std::collections::HashMap;

use crate::{graph::Graph, reads::Reads};
use crate::reads::K;

fn get_unique_path(paths: Vec<Vec<u64>>) -> Vec<Option<u64>> {
    paths.first()
        .unwrap()
        .iter()
        .enumerate()
        .map(|(i, node)| {
            let all_paths_share_same_node: bool = paths.iter()
                .all(|node_to_check| {
                    node == node_to_check.get(i).unwrap()
                });
            if all_paths_share_same_node {
                Some(*node)
            } else {
                None
            }
        })
        .collect::<Vec<Option<u64>>>()
}

fn extend_sequence(locked_sequence: &mut String, sequence_to_add: &String) -> bool {
    if locked_sequence.len() < K || sequence_to_add.len() < K {
        return false
    }
    
    let mut index: usize = locked_sequence.len().saturating_sub(K) + 1;
    let mut left_matches: Vec<usize> = vec![];
    let mut center_matches: Vec<usize> = vec![];
    let mut right_matches: Vec<usize> = vec![];

    // right
    while index > 0 {
        index -= 1;

        let early_termination: bool = locked_sequence.len() - index
        > sequence_to_add.len();
        if early_termination {
            index = 0;
            continue
        }
        
        let is_matching: bool = locked_sequence.chars()
            .skip(index)
            .enumerate()
            .all(|(i, character)| {
                let add_seq_char = sequence_to_add.chars().nth(i).unwrap();
                add_seq_char == 'N'
                || character == 'N'
                || character == add_seq_char
            });
        if is_matching {
            right_matches.push(index)
        }
    }

    // middle
    index = locked_sequence.len();

    while index > 0 {
        index -= 1;
        
        let is_matching: bool = locked_sequence.chars()
            .skip(index)
            .enumerate()
            .all(|(i, character)| {
                if i < sequence_to_add.len() {
                    let locked_seq_char = sequence_to_add.chars().nth(i).unwrap();
                    locked_seq_char == 'N'
                    || character == 'N'
                    || character == locked_seq_char
                } else {
                    true
                }
            });
        if is_matching {
            center_matches.push(index)
        }
    }
    
    // left
    index = sequence_to_add.len().saturating_sub(K) + 1;

    while index > 0 {
        index -= 1;
        
        let is_matching: bool = sequence_to_add.chars()
            .skip(index)
            .enumerate()
            .all(|(i, character)| {
                let locked_seq_char = locked_sequence.chars().nth(i).unwrap();
                locked_seq_char == 'N'
                || character == 'N'
                || character == locked_seq_char
            });
        if is_matching {
            left_matches.push(index)
        }
    }

    // apply matches
    if right_matches.len() > 0 {
        let match_index = *right_matches.last().unwrap();
        let inside_chunk: String = sequence_to_add.chars().take(locked_sequence.len() - match_index).collect();
        let outside_chunk: String = sequence_to_add.chars().skip(locked_sequence.len() - match_index).collect();
        
        *locked_sequence = locked_sequence.chars()
            .enumerate()
            .map(|(i, char)| {
                if i >= match_index && char == 'N' {
                    inside_chunk.chars().nth(i - match_index).unwrap()
                } else {
                    char
                }
            })
            .collect();
        
        locked_sequence.push_str(&outside_chunk);
    }
    for match_index in center_matches {
        *locked_sequence = locked_sequence.chars()
            .enumerate()
            .map(|(i, char)| {
                if i >= match_index
                && i - match_index < sequence_to_add.len()
                && char == 'N' {
                    sequence_to_add.chars().nth(i - match_index).unwrap()
                } else {
                    char
                }
            }).collect();
    }
    if left_matches.len() > 0 {
        let match_index = *left_matches.last().unwrap();
        let outside_chunk: String = sequence_to_add.chars().take(match_index).collect();
        let inside_chunk: String = sequence_to_add.chars().skip(match_index).collect();
        
        *locked_sequence = locked_sequence.chars()
            .enumerate()
            .map(|(i, char)| {
                if i < inside_chunk.len() && char == 'N' {
                    inside_chunk.chars().nth(i).unwrap()
                } else {
                    char
                }
            })
            .collect();
        locked_sequence.insert_str(0, &outside_chunk);
    }


    left_matches.len() > 0 || right_matches.len() > 0
}

pub fn greedy_assembly(graph: Graph, reads: Reads, debug: bool) -> String {
    let sequence_to_node_id_map: HashMap<String, u64> = graph.nodes.iter()
        .map(|(node_id, sequence)| (sequence.clone(), *node_id))
        .collect();
    
    let mut sorted_reads = reads.reads.clone()
        .into_iter()
        .map(|read| {
            let start_node_id: u64 = *sequence_to_node_id_map.get(&read.start_k_mer).unwrap();
            let end_node_id: u64 = *sequence_to_node_id_map.get(&read.end_k_mer).unwrap();
            let paths: Vec<Vec<u64>> = graph.find_path(start_node_id, end_node_id, read.nodes_between + 1);
            
            let unique_paths: Vec<Option<u64>> = get_unique_path(paths);
            let paths_seq = unique_paths.into_iter()
                .map(|node_id_option| {
                    if let Some(node_id) = node_id_option {
                        if let Some(seq) = graph.nodes.get(&node_id) {
                            Some(seq.clone())
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect::<Vec<Option<String>>>();

            let path_sequence: String = paths_seq.first().unwrap().as_ref().unwrap().clone()
            + &paths_seq.into_iter()
                .skip(1)
                .map(|seq_optional| {
                    if let Some(seq) = seq_optional {
                        seq.chars().last().unwrap()
                    } else {
                        'N'
                    }
                })
                .collect::<String>();

            path_sequence
        }).collect::<Vec<String>>();
    
    sorted_reads.sort_by(|seq_a, seq_b| {
        seq_b.chars().filter(|&c| c == 'N').count()
        .cmp(
            &seq_a.chars().filter(|&c| c == 'N').count()
        )
    });

    // println!("{:?}", sorted_reads);
    let mut final_sequence: String = sorted_reads.pop().unwrap();

    let mut reads_cycling_counter: usize = 0;

    while sorted_reads.len() > 0 && reads_cycling_counter < sorted_reads.len() {
        let read_seq: &String = sorted_reads.get(reads_cycling_counter).unwrap();
        
        let successful: bool = extend_sequence(&mut final_sequence, read_seq);
        if successful {
            sorted_reads.remove(reads_cycling_counter);
            reads_cycling_counter = 0;
        } else {
            reads_cycling_counter += 1;
        }
    }

    if debug && sorted_reads.len() > 0 {
        println!("Only {} out of {} reads were considered",
            reads.reads.len().saturating_sub(sorted_reads.len()),
            reads.reads.len()
        );
    }
    
    final_sequence
}
