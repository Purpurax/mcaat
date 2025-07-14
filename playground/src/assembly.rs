use std::collections::{HashMap, HashSet};

use crate::{graph::Graph, reads::Reads};

fn get_node_to_cycle_map(raw_cycles: Vec<Vec<u64>>) -> HashMap<u64, Vec<usize>> {
    let mut node_to_cycle: HashMap<u64, Vec<usize>> = HashMap::new();
    
    for (node_id, cycle_index) in raw_cycles.iter()
        .enumerate()
        .flat_map(|(i, cycle)| {
            cycle.into_iter()
                .map(move |node_id| {
                    (*node_id, i)
                })
    }) {
        if node_to_cycle.contains_key(&node_id) {
            node_to_cycle.get_mut(&node_id).unwrap().push(cycle_index);
        } else {
            node_to_cycle.insert(node_id, vec![cycle_index]);
        }
    }

    node_to_cycle
}

fn map_reads_to_node_id_pairs(reads: Reads, graph: &Graph) -> Vec<(u64, u64, usize)> {
    let sequence_to_node_id_map: HashMap<String, u64> = graph.nodes.iter()
        .map(|(node_id, sequence)| (sequence.clone(), *node_id))
        .collect();

    reads.reads.into_iter()
        .map(|read| {
            let start_node_id: u64 = *sequence_to_node_id_map.get(&read.start_k_mer).unwrap();
            let end_node_id: u64 = *sequence_to_node_id_map.get(&read.end_k_mer).unwrap();
            
            (start_node_id, end_node_id, read.nodes_between)
        })
        .collect::<Vec<(u64, u64, usize)>>()
}

fn get_jumps(graph: &Graph, reads: Reads, raw_cycles: Vec<Vec<u64>>) -> Vec<(usize, usize)> {
    let node_to_cycle: HashMap<u64, Vec<usize>> = get_node_to_cycle_map(raw_cycles);

    let mut jumps: Vec<(usize, usize)> = vec![];

    for (start_node_id, end_node_id, node_in_between) in map_reads_to_node_id_pairs(reads, &graph) {
        let paths = graph.find_path(start_node_id, end_node_id, node_in_between + 1);

        let mut unique_cycles_on_path = paths.first()
            .unwrap()
            .iter()
            .enumerate()
            .filter_map(|(i, _)| {
                let nodes_cycles = paths.iter()
                    .map(|path| {
                        path.get(i).unwrap()
                    })
                    .filter(|node| {
                        node_to_cycle.contains_key(node)
                    })
                    .flat_map(|node| {
                        node_to_cycle.get(&node).unwrap()
                    })
                    .map(|x| *x)
                    .collect::<HashSet<usize>>()
                    .into_iter()
                    .collect::<Vec<usize>>();
                
                if nodes_cycles.len() == 1 {
                    Some(*nodes_cycles.first().unwrap())
                } else {
                    None
                }
            })
            .collect::<Vec<usize>>();
        
        if unique_cycles_on_path.len() < 2 {
            println!("The read does not contain a good jump");
        }

        while unique_cycles_on_path.len() >= 2 {
            let first_cycle_index: usize = unique_cycles_on_path.remove(0);
            jumps.push((
                first_cycle_index,
                *unique_cycles_on_path.first().unwrap()
            ))
        }
    }

    jumps.into_iter()
        .filter(|(start_cycle, end_cycle)| {
            start_cycle != end_cycle
        }).collect::<Vec<(usize, usize)>>()
}

fn left_hand_reduce_jumps(jumps: Vec<(usize, usize)>) -> HashMap<usize, Vec<usize>> {
    let mut reduced_jumps: HashMap<usize, Vec<usize>> = HashMap::new();

    for (start_cycle, end_cycle) in jumps.into_iter() {
        if reduced_jumps.contains_key(&start_cycle) {
            let already_contained_end_cycles = reduced_jumps.get_mut(&start_cycle).unwrap();

            if !already_contained_end_cycles.contains(&end_cycle) {
                already_contained_end_cycles.push(end_cycle);
            }
        } else {
            reduced_jumps.insert(start_cycle, vec![end_cycle]);
        }
    }

    reduced_jumps
}

fn get_possible_start_cycles(jumps: &HashMap<usize, Vec<usize>>) -> usize {
    let all_start_cycles = jumps.keys().into_iter().map(|x| *x).collect::<Vec<usize>>();

    let mut sorted_jumps = jumps.into_iter().map(|(k, v)| (*k, v.clone())).collect::<Vec<(usize, Vec<usize>)>>();
    
    sorted_jumps.sort_by(|(_, landing_cycle_indices_a), (_, landing_cycle_indices_b)| {
            landing_cycle_indices_a.len().cmp(&landing_cycle_indices_b.len())
        });

    if sorted_jumps.len() == 0 {
        panic!("There are not enough reads, as no start_cycle was found");
    }
    if sorted_jumps.first().unwrap().1.len() > 1 {
        panic!("The start_cycle is not clear as every node has an incoming edge");
    }
    if sorted_jumps.len() >= 1
    && sorted_jumps.get(1).unwrap().1.len() == 1 {
        println!("WARNING: multiple start cycles are possible, and one of them has been chosen");
    }

    sorted_jumps.first().unwrap().0
}

fn reconstruct_cycle_order_from_start(start_cycle: usize, jumps: HashMap<usize, Vec<usize>>) -> Vec<usize> {
    let mut next_cycles: Vec<usize> = jumps.get(&start_cycle).unwrap().clone();

    let mut cycle_order: Vec<usize> = vec![start_cycle];

    while !next_cycles.is_empty() {
        if next_cycles.len() == 1 {
            let next_cycle = next_cycles.pop().unwrap();
            cycle_order.push(next_cycle);

            next_cycles = jumps.get(&next_cycle).unwrap().clone();
        } else {
            panic!("WARNING: The next cycle is ambiguous");
        }
    }

    cycle_order
}

pub fn assembly(graph: Graph, reads: Reads, raw_cycles: Vec<Vec<u64>>, debug: bool) -> String {
    let full_jumps = get_jumps(&graph, reads, raw_cycles);
    let jumps = left_hand_reduce_jumps(full_jumps);
    println!("jumps: {:?}", jumps);
    
    let start_cycle: usize = get_possible_start_cycles(&jumps);
    
    let cycle_order: Vec<usize> = reconstruct_cycle_order_from_start(start_cycle, jumps);
    println!("{:?}", cycle_order);

    "".to_string()
}
