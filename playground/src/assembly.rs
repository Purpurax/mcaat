use std::{collections::{HashMap, HashSet}, fs::File, io};
use std::io::Write;

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

    jumps.into_iter()
        .filter(|(start_cycles, end_cycles)| {
            !start_cycles.into_iter()
                .all(|start_cycle| {
                    end_cycles.contains(start_cycle)
                })
        })
        .collect::<Vec<(Vec<usize>, Vec<usize>)>>()
}

fn get_nodes(jumps: &Vec<(Vec<usize>, Vec<usize>)>) -> Vec<usize> {
    let mut all_cycles = HashSet::new();

    jumps.clone().into_iter()
        .for_each(|(start_cycles, end_cycles)| {
            for start_cycle in start_cycles {
                all_cycles.insert(start_cycle);
            }
            for end_cycle in end_cycles {
                all_cycles.insert(end_cycle);
            }
        });

    all_cycles.into_iter()
        .collect::<Vec<usize>>()
}

// start, end, weight
fn get_edges(jumps: &Vec<(Vec<usize>, Vec<usize>)>) -> HashMap<(usize, usize), f64> {
    let all_cycles: Vec<usize> = get_nodes(jumps);

    let mut edges = all_cycles.clone().into_iter()
        .flat_map(|node| {
            all_cycles.clone().into_iter()
                .map(move |other_node| {
                    ((node, other_node), 0.0)
                })
        }).collect::<HashMap<(usize, usize), f64>>();

    for (start_cycles, end_cycles) in jumps.clone() {
        for start_cycle in start_cycles.clone() {
            for end_cycle in end_cycles.clone() {
                let confidence = (start_cycles.len() * end_cycles.len()) as f64;
                *edges.get_mut(&(start_cycle, end_cycle)).unwrap() += confidence;
            }
        }
    }

    edges
}

pub fn export_to_dot(nodes: &Vec<usize>, edges: &HashMap<(usize, usize), f64>, file_path: &str) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    writeln!(file, "digraph G {{")?;

    for source in nodes.clone() {
        for destination in nodes.clone() {
            let weight = edges.get(&(source, destination)).unwrap();
            writeln!(file, "  \"{}\" -> \"{}\" [label=\"{}\"];", source, destination, weight)?;
        }
    }

    writeln!(file, "}}")?;
    Ok(())
}

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

fn reconstruct_cycle_order_using_constraint_graph(
    nodes: &Vec<usize>,
    edges: &HashMap<(usize, usize), f64>,
    start_node: usize
) -> (Vec<usize>, Vec<f64>) {
    let mut order = vec![];
    let mut overall_confidence: Vec<f64> = vec![];
    let mut current_node_opt: Option<usize> = Some(start_node);

    while let Some(current_node) = current_node_opt {
        order.push(current_node);
        // greedy by confidence
        let remaining_nodes = nodes.clone()
            .into_iter()
            .filter(|node| {
                !order.contains(node)
            });
        
        let mut best_node: Option<usize> = None;
        let mut best_confidence: f64 = 0.0;

        let normalize_confidence: f64 = nodes.clone()
            .into_iter()
            .map(|node| {
                *edges.get(&(current_node, node)).unwrap()
            })
            .sum::<f64>();

        for destination in remaining_nodes {
            let edge_weight = *edges.get(&(current_node, destination)).unwrap();
            if edge_weight > best_confidence {
                best_confidence = edge_weight;
                best_node = Some(destination);
            }
        }

        if best_node.is_some() {
            overall_confidence.push(best_confidence / normalize_confidence);
        }

        current_node_opt = best_node;
    }

    (order, overall_confidence)
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
    let nodes = get_nodes(&jumps);
    let edges = get_edges(&jumps);

    let _ = export_to_dot(&nodes, &edges, "./results/assemply-graph.dot");
    
    let start_cycles: Vec<usize> = get_possible_start_cycles(&jumps);
    let orders_with_confidence = start_cycles.into_iter().map(|start_cycle| {
        // let cycle_order: Vec<usize> = reconstruct_cycle_order_from_start(start_cycle, &jumps);
        let (cycle_order, partial_confidences): (Vec<usize>, Vec<f64>) =
            reconstruct_cycle_order_using_constraint_graph(&nodes, &edges, start_cycle);
        let confidence = partial_confidences.into_iter().sum::<f64>();
        println!("{} has cycle order {:?} with a confidence of {:?}", start_cycle, cycle_order, confidence);

        (cycle_order, confidence)
    });

    let best_order = orders_with_confidence.max_by(|(_, conf_a), (_, conf_b)| {
        conf_a.partial_cmp(conf_b).unwrap()
    })
    .map(|(order, _)| order)
    .unwrap();

    println!("Best order is {:?}", best_order);

    "".to_string()
}
