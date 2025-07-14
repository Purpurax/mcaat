use std::{collections::{HashMap, HashSet}, hash::Hash};

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

fn get_jumps(graph: &Graph, reads: Reads, raw_cycles: Vec<Vec<u64>>) -> Vec<(Vec<usize>, Vec<usize>)> {
    let node_to_cycle: HashMap<u64, Vec<usize>> = get_node_to_cycle_map(raw_cycles);

    let mut jumps: Vec<(Vec<usize>, Vec<usize>)> = vec![];

    for (start_node_id, end_node_id, node_in_between) in map_reads_to_node_id_pairs(reads, &graph) {
        let paths = graph.find_path(start_node_id, end_node_id, node_in_between + 1);

        let mut cycle_indices_order = paths.first()
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, _)| {
                paths.iter()
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
                    .collect::<Vec<usize>>()
            })
            .collect::<Vec<Vec<usize>>>();

        while cycle_indices_order.len() >= 2 {
            let first_cycle_index: Vec<usize> = cycle_indices_order.remove(0);

            for landing_cycle_index in cycle_indices_order.iter() {
                jumps.push((
                    first_cycle_index.clone(),
                    landing_cycle_index.clone()
                ))
            }
        }
    }

    jumps = jumps.into_iter()
        .filter(|(start_cycles, end_cycles)| {
            !start_cycles.into_iter()
                .all(|start_cycle| {
                    end_cycles.contains(start_cycle)
                })
        })
        .collect::<HashSet<(Vec<usize>, Vec<usize>)>>()
        .into_iter()
        .collect::<Vec<(Vec<usize>, Vec<usize>)>>();

    jumps.sort_by(|(start_a, end_a), (start_b, end_b)| {
        let a_score = start_a.len() + end_a.len();
        let b_score = start_b.len() + end_b.len();

        a_score.cmp(&b_score)
    });

    jumps
}

// fn left_hand_reduce_jumps(jumps: Vec<(usize, usize)>) -> HashMap<usize, Vec<usize>> {
//     let mut reduced_jumps: HashMap<usize, Vec<usize>> = HashMap::new();

//     for (start_cycle, end_cycle) in jumps.into_iter() {
//         if reduced_jumps.contains_key(&start_cycle) {
//             let already_contained_end_cycles = reduced_jumps.get_mut(&start_cycle).unwrap();

//             if !already_contained_end_cycles.contains(&end_cycle) {
//                 already_contained_end_cycles.push(end_cycle);
//             }
//         } else {
//             reduced_jumps.insert(start_cycle, vec![end_cycle]);
//         }
//     }

//     reduced_jumps
// }

fn get_possible_start_cycles(jumps: &Vec<(Vec<usize>, Vec<usize>)>) -> Vec<usize> {
    let possible_cycles = jumps.into_iter()
        .flat_map(|(start_cycles, _)| {
            start_cycles.clone()
        })
        .collect::<HashSet<usize>>()
        .into_iter()
        .collect::<Vec<usize>>();
    
    possible_cycles
}

fn reconstruct_cycle_order_from_start(start_cycle: usize, jumps: &Vec<(Vec<usize>, Vec<usize>)>) -> Vec<usize> {
    let mut cycle_order: Vec<usize> = vec![];
    let mut next_cycle_opt: Option<usize> = Some(start_cycle);

    while let Some(next_cycle) = next_cycle_opt {
        cycle_order.push(next_cycle);

        // remove cycle from jumps
        let filtered_jumps = jumps.into_iter()
            .map(|(start_cycles, end_cycles)| {
                (
                    start_cycles.into_iter()
                        .filter(|start_cycle| {
                            !cycle_order.contains(*start_cycle)
                            || **start_cycle == next_cycle
                        })
                        .map(|x| *x)
                        .collect::<Vec<usize>>(),
                    end_cycles.into_iter()
                        .filter(|end_cycle| {
                            !cycle_order.contains(*end_cycle)
                            || **end_cycle == next_cycle
                        })
                        .map(|x| *x)
                        .collect::<Vec<usize>>()
                )
            })
            .filter(|(start, end)| {
                start.len() != 0 && end.len() != 0
            })
            .collect::<Vec<(Vec<usize>, Vec<usize>)>>();
        // println!("f_j: {:?}", filtered_jumps.clone()
        //     .into_iter()
        //     .filter(|(start_cycles, end_cycles)| {
        //         start_cycles.len() == 1 || end_cycles.len() == 1
        //     })
        //     .collect::<Vec<(Vec<usize>, Vec<usize>)>>()
        // );

        let mut cycle_to_cycles_map: HashMap<usize, Vec<usize>> = HashMap::new();
        
        for (start_cycle, end_cycle) in  filtered_jumps.into_iter()
            // .filter(|(start_cycles, end_cycles)| {
            //     start_cycles.len() == 1 && end_cycles.len() == 1
            // })
            // .map(|(start_cycles, end_cycles)| {
            //     (start_cycles.first().unwrap().clone(), end_cycles.first().unwrap().clone())
            // }) {
            .filter(|(start_cycles, end_cycles)| {
                start_cycles.len() == 1 && end_cycles.len() == 1
            })
            .map(|(start_cycles, end_cycles)| {
                (start_cycles.first().unwrap().clone(), end_cycles.first().unwrap().clone())
            }) {

            if cycle_to_cycles_map.contains_key(&start_cycle) {
                cycle_to_cycles_map.get_mut(&start_cycle).unwrap().push(end_cycle);
            } else {
                cycle_to_cycles_map.insert(start_cycle, vec![end_cycle]);
            }
        }

        let cycle_to_cycle_map = cycle_to_cycles_map.clone().into_iter()
            .filter(|(_, values)| {
                values.len() == 1
            })
            .map(|(key, values)| {
                (key, *values.first().unwrap())
            })
            .collect::<HashMap<usize, usize>>();
        // println!("{:?}", cycle_to_cycles_map);

        next_cycle_opt = cycle_to_cycle_map.get(&next_cycle).copied();
    }

    cycle_order
}

pub fn assembly(graph: Graph, reads: Reads, raw_cycles: Vec<Vec<u64>>, debug: bool) -> String {
    let jumps = get_jumps(&graph, reads, raw_cycles);
    // let jumps = left_hand_reduce_jumps(full_jumps);
    // println!("jumps: {:?}", jumps);
    
    let start_cycles: Vec<usize> = get_possible_start_cycles(&jumps);
    
    for start_cycle in start_cycles {
        let cycle_order: Vec<usize> = reconstruct_cycle_order_from_start(start_cycle, &jumps);
        println!("{} has cycle order: {:?}", start_cycle, cycle_order);
    }

    "".to_string()
}
