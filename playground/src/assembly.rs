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

fn sequences_match_at(
    seq_a: &String,
    seq_b: &String,
    index_in_a: usize
) -> bool {
    seq_a.chars()
        .skip(index_in_a)
        .enumerate()
        .all(|(i, character)| {
            if i < seq_b.len() {
                let other_character = seq_b.chars().nth(i).unwrap();
                
                other_character == 'N'
                || character == 'N'
                || character == other_character
            } else {
                true
            }
        })
}

fn overlay_onto_sequence_at(
    seq_a: &mut String,
    seq_b: &String,
    index_in_a: usize
) {
    *seq_a = seq_a.chars()
        .enumerate()
        .map(|(i, char)| {
            if i >= index_in_a
            && i - index_in_a < seq_b.len()
            && char == 'N' {
                seq_b.chars().nth(i - index_in_a).unwrap()
            } else {
                char
            }
        }).collect();
}

fn extend_sequence_center(locked_sequence: &mut String, sequence_to_add: &String) -> bool {
    let matching_indices = (0..(locked_sequence.len() - 1)).into_iter()
        .filter(|index| {
            sequences_match_at(&locked_sequence, sequence_to_add, *index)
        })
        .collect::<Vec<usize>>();
    
    for index in matching_indices.iter() {
        overlay_onto_sequence_at(locked_sequence, sequence_to_add, *index);
    }
    
    !matching_indices.is_empty()
}

fn extend_sequence_left(locked_sequence: &mut String, sequence_to_add: &String) -> bool {
    let matching_indices = (0..(sequence_to_add.len().saturating_sub(K) + 1))
        .into_iter()
        .filter(|index| {
            sequences_match_at(sequence_to_add, &locked_sequence, *index)
        })
        .collect::<Vec<usize>>();
    
    if matching_indices.len() > 0 {
        let match_index = *matching_indices.last().unwrap();
        let outside_chunk: String = sequence_to_add.chars().take(match_index).collect();
        let inside_chunk: String = sequence_to_add.chars().skip(match_index).collect();
        
        overlay_onto_sequence_at(locked_sequence, &inside_chunk, 0);
        locked_sequence.insert_str(0, &outside_chunk);

        true
    } else {
        false
    }
}

fn extend_sequence_right(locked_sequence: &mut String, sequence_to_add: &String) -> bool {
    let matching_indices = (0..(locked_sequence.len().saturating_sub(K) + 1))
        .into_iter()
        .filter(|index| {
            locked_sequence.len() - index <= sequence_to_add.len()
        })
        .filter(|index| {
            sequences_match_at(&locked_sequence, sequence_to_add, *index)
        })
        .collect::<Vec<usize>>();

    if matching_indices.len() > 0 {
        let match_index = *matching_indices.last().unwrap();
        let inside_chunk: String = sequence_to_add.chars().take(locked_sequence.len() - match_index).collect();
        let outside_chunk: String = sequence_to_add.chars().skip(locked_sequence.len() - match_index).collect();
        
        overlay_onto_sequence_at(locked_sequence, &inside_chunk, match_index);
        locked_sequence.push_str(&outside_chunk);

        true
    } else {
        false
    }
}

fn extend_sequence(locked_sequence: &mut String, sequence_to_add: &String) -> bool {
    if locked_sequence.len() < K || sequence_to_add.len() < K {
        return false
    }

    extend_sequence_right(locked_sequence, sequence_to_add)
    || extend_sequence_center(locked_sequence, sequence_to_add)
    || extend_sequence_left(locked_sequence, sequence_to_add)
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
