use std::{collections::{HashMap, HashSet}, fs::File, io};
use std::io::Write;

use crate::{graph::Graph, reads::Reads};
use rand::seq::SliceRandom;
use itertools::Itertools;

// Set cover problem:
//   All node_ids have to be covered with as little cycles possible
fn get_node_to_cycle_map(raw_cycles: &Vec<Vec<u64>>) -> HashMap<u64, usize> {
    let sets = raw_cycles.clone()
        .into_iter()
        .map(|cycle| cycle.into_iter().collect::<HashSet<u64>>())
        .collect::<Vec<HashSet<u64>>>();
    let total_set = raw_cycles.clone()
        .into_iter()
        .flatten()
        .collect::<HashSet<u64>>();

    // let cycle_indices = set_cover_solver(sets, total_set);
    let cycle_indices = set_cover_brute_force(sets.clone(), total_set.clone());

    raw_cycles.clone()
        .into_iter()
        .enumerate()
        .filter(|(cycle_index, _cycle)| cycle_indices.contains(cycle_index))
        .flat_map(|(cycle_index, cycle)| {
            cycle.into_iter()
                .filter({
                    let value = sets.clone();
                    let cycle_indices_copy = cycle_indices.clone();
                    move |node_id| {
                        value.clone().into_iter()
                            .enumerate()
                            .filter(|(set_index, _set)| cycle_indices_copy.contains(set_index))
                            .filter(|(set_index, _set)| *set_index != cycle_index)
                            .all(|(_, set)| !set.contains(&node_id))
                    }
                })
                .map(move |node_id| {
                    (node_id, cycle_index)
                })
        })
        .collect::<HashMap<u64, usize>>()

    // raw_cycles.clone()
    //     .into_iter()
    //     .enumerate()
    //     .filter(|(cycle_index, _cycle)| cycle_indices.contains(cycle_index))
    //     .flat_map(|(cycle_index, cycle)| {
    //         cycle.into_iter()
    //             .map(move |node_id| {
    //                 (node_id, cycle_index)
    //             })
    //     })
    //     .collect::<HashMap<u64, usize>>()
}

fn set_cover_solver(
    sets: Vec<HashSet<u64>>,
    universe: HashSet<u64>
) -> Vec<usize> {
    let uniquely_covering_sets = sets.clone().into_iter()
        .enumerate()
        .filter(|(set_index, set)| {
            let mut unique_elements: HashSet<u64> = set.clone();

            sets.clone().into_iter()
                .enumerate()
                .filter(|(other_set_index, _)| set_index != other_set_index)
                .for_each(|(_, other_set)| {
                    unique_elements = unique_elements.difference(&other_set).map(|x| *x).collect();
                });

            set.len() > 0
        });
    
    let mut result = uniquely_covering_sets.map(|(index, _)| index).collect::<Vec<usize>>();
    let current_coverage = sets.clone()
        .into_iter()
        .enumerate()
        .filter(|(set_index, _)| result.contains(set_index))
        .map(|(_, set)| set)
        .flatten()
        .collect::<HashSet<u64>>();
    let mut missing_elements = universe.difference(&current_coverage)
        .map(|x| *x)
        .collect::<Vec<u64>>();

    while missing_elements.len() > 0
    && result.len() < sets.len() {
        let possible_indices = (0..sets.len() - 1).filter(|i| !result.contains(i));

        // TODO

    }

    result
}

fn set_cover_brute_force(
    sets: Vec<HashSet<u64>>,
    universe: HashSet<u64>
) -> Vec<usize> {
    let n = sets.len();
    let universe_vec: Vec<u64> = universe.iter().cloned().collect();

    for size in 1..=n {
        let indices: Vec<usize> = (0..n).collect();
        let combos = indices.iter().cloned().combinations(size);
        for combo in combos {
            let mut covered: HashSet<u64> = HashSet::new();
            for &i in &combo {
                covered.extend(&sets[i]);
            }
            if universe_vec.iter().all(|u| covered.contains(u)) {
                return combo.clone()
            }
        }
    }

    return vec![]
}

fn get_cycles_with_unique_nodes(raw_cycles: &Vec<Vec<u64>>) -> HashMap<u64, usize> {
    raw_cycles.iter() 
        .enumerate()
        .map(|(i, raw_cycle)| {
            let cycle_specific_nodes = raw_cycle.iter()
                .map(|node_id| *node_id)
                .filter(|node_id| {
                    raw_cycles.iter()
                        .enumerate()
                        .filter(|(j, _)| i != *j)
                        .all(|(_, other_cycle)| {
                            other_cycle.iter()
                                .all(|other_node_id| {
                                    node_id != other_node_id
                                })
                        })
                })
                .collect::<Vec<u64>>();
            (i, cycle_specific_nodes)
        })
        .filter(|(_, cycle_specific_nodes)| cycle_specific_nodes.len() > 0)
        .flat_map(|(i, cycle_specific_nodes)| {
            cycle_specific_nodes.into_iter()
                .map(move |node_id| {
                    (node_id, i)
                })
        })
        .collect::<HashMap<u64, usize>>()
}

fn map_reads_to_node_id_pairs(reads: Reads, graph: &Graph) -> Vec<(u64, u64, u64, usize)> {
    let sequence_to_node_id_map: HashMap<String, u64> = graph.nodes.iter()
        .map(|(node_id, sequence)| (sequence.clone(), *node_id))
        .collect();

    reads.reads.into_iter()
        .map(|read| {
            let start_node_id: u64 = *sequence_to_node_id_map.get(&read.start_k_mer).unwrap();
            let middle_node_id: u64 = *sequence_to_node_id_map.get(&read.middle_k_mer).unwrap();
            let end_node_id: u64 = *sequence_to_node_id_map.get(&read.end_k_mer).unwrap();
            
            (start_node_id, middle_node_id, end_node_id, read.nodes_between)
        })
        .collect::<Vec<(u64, u64, u64, usize)>>()
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

fn get_reads_path(graph: &Graph, reads: &Reads) -> Vec<Vec<u64>> {
    const K: usize = 23;

    reads.reads.clone().into_iter()
        .map(|read| {
            let sequence = read.sequence;

            let mut node_seq = vec![];
            let mut counter = 0;
            while counter + K < sequence.chars().count() {
                node_seq.push(
                    sequence.chars()
                        .skip(counter)
                        .take(K)
                        .collect::<String>()
                );

                counter += 1;
            }

            let node_id_to_seq: HashMap<String, u64> = graph.nodes
                .iter()
                .map(|(node_id, seq)| (seq.clone(), *node_id))
                .collect();

            node_seq.into_iter()
                .map(|node_s| {
                    node_id_to_seq.get(&node_s)
                })
                .filter_map(|res| res)
                .map(|x| *x)
                .collect::<Vec<u64>>()
        })
        .collect::<Vec<Vec<u64>>>()
}

fn generate_constraints(graph: &Graph, reads: Reads, node_to_cycle: HashMap<u64, usize>) -> Vec<(usize, usize)> {
    let mut constraints: Vec<(usize, usize)> = vec![];

    for path in get_reads_path(graph, &reads).into_iter() {
        let cycles_on_path = path.into_iter()
            .filter(|node_id| {
                node_to_cycle.contains_key(node_id)
            })
            .map(|node_id| {
                *node_to_cycle.get(&node_id).unwrap()
            })
            .collect::<Vec<usize>>();
        
        if cycles_on_path.len() == 0 {
            println!("The read goes through no cycle specific node");
        } else if cycles_on_path.len() == 1 {
            println!("The read is visiting only one cycle specific node: {}", cycles_on_path.first().unwrap());
        }
        
        cycles_on_path.clone()
            .into_iter()
            .enumerate()
            .for_each(|(i, cycle_index)| {
                cycles_on_path.iter()
                    .enumerate()
                    .filter(|(j, _)| i < *j)
                    .filter(|(_, other_cycle_index)| cycle_index != **other_cycle_index)
                    .for_each(|(_, other_cycle_index)| {
                        constraints.push((cycle_index, *other_cycle_index));
                    })
            });
    }

    constraints
}

fn extend_constraints_by_transitive(constraints: &mut Vec<(usize, usize)>) {
    
}

fn sort_topologically(constraints: &Vec<(usize, usize)>) -> Vec<Vec<usize>> {
    let mut edges: HashMap<(usize, usize), usize> = HashMap::new();

    constraints.into_iter()
        .for_each(|(source, destination)| {
            let weight_opt = edges.get_mut(&(*source, *destination));
            if let Some(weight) = weight_opt {
                *weight += 1;
            } else {
                edges.insert((*source, *destination), 0);
            }
        });
    let nodes = edges.clone().into_iter()
        .flat_map(|((source, destination), _)| {
            vec![source, destination]
        })
        .collect::<HashSet<usize>>()
        .into_iter()
        .collect::<Vec<usize>>();

    println!("Nodes: {:?}", nodes);
    println!("Edges: {:?}", edges);
    resolve_cycles(&mut edges);

    let possible_start_nodes = nodes.clone()
        .into_iter()
        .filter(|node| {
            nodes.iter()
                .filter(|source_node| {
                    edges.contains_key(&(**source_node, *node))
                })
                .count() == 0
        })
        .collect::<Vec<usize>>();

    apply_topological_sort(possible_start_nodes, &nodes, edges)
}

fn resolve_cycles(edges_with_weights: &mut HashMap<(usize, usize), usize>) {
    let mut cycles = cycle_finder(&edges_with_weights);

    while cycles.len() > 0 {
        let mut edges_sorted = edges_with_weights.clone().into_iter()
            .collect::<Vec<((usize, usize), usize)>>();
        edges_sorted.sort_by(|((_, _), weight_a), ((_, _), weight_b)| {
            weight_a.cmp(weight_b)
        });

        let removed_edge = edges_sorted.remove(0);
        edges_with_weights.remove(&(removed_edge.0.0, removed_edge.0.1));
        cycles = cycle_finder(&edges_with_weights);
    }
}

fn cycle_finder(edges_with_weights: &HashMap<(usize, usize), usize>) -> Vec<Vec<usize>> {
    let mut all_cycles = Vec::new();
    let nodes: Vec<usize> = edges_with_weights.clone().into_iter()
        .flat_map(|((source, destination), _)| {
            vec![source, destination].into_iter()
        })
        .collect::<HashSet<usize>>()
        .into_iter()
        .collect::<Vec<usize>>();
    let mut edges: HashMap<usize, HashSet<usize>> = HashMap::new();
    for ((source, destination), _) in edges_with_weights.clone() {
        let destinations_opt = edges.get_mut(&source);
        if let Some(destinations) = destinations_opt {
            destinations.insert(destination);
        } else {
            let mut hashy = HashSet::new();
            hashy.insert(destination);
            edges.insert(source, hashy);
        }
    }

    for &start_node in &nodes {
        let mut stack = Vec::new();
        let mut visited = HashSet::new();
        dfs_cycles(
            &edges,
            start_node,
            start_node,
            &mut stack,
            &mut visited,
            &mut all_cycles,
        );
    }

    all_cycles
}

fn dfs_cycles(
    edges: &HashMap<usize, HashSet<usize>>,
    start_node: usize,
    current_node: usize,
    stack: &mut Vec<usize>,
    visited: &mut HashSet<usize>,
    all_cycles: &mut Vec<Vec<usize>>,
) {
    stack.push(current_node);
    visited.insert(current_node);

    if let Some(neighbors) = edges.get(&current_node) {
        for &neighbor in neighbors {
            if neighbor == start_node && stack.len() > 1 {
                // Found a cycle
                all_cycles.push(stack.clone());
            } else if !visited.contains(&neighbor) {
                dfs_cycles(edges, start_node, neighbor, stack, visited, all_cycles);
            }
        }
    }

    stack.pop();
    visited.remove(&current_node);
}

fn apply_topological_sort(
    possible_start_nodes: Vec<usize>,
    nodes: &Vec<usize>,
    edges: HashMap<(usize, usize), usize>
) -> Vec<Vec<usize>> {
    if possible_start_nodes.len() == 0 {
        return vec![vec![]]
    }

    possible_start_nodes.clone().into_iter()
        .flat_map(|start_node| {
            let new_possible_start_nodes = possible_start_nodes.clone()
                .into_iter()
                .filter(|node| *node != start_node)
                .chain(
                    edges.clone().into_iter()
                        .filter(|((source, _), _)| {
                            *source == start_node
                        })
                        .map(|((_, destination), _)| destination)
                )
                .collect::<Vec<usize>>();
            let new_edges = edges.clone()
                .into_iter()
                .filter(|((source, destination), _)| {
                    *source != start_node
                    && start_node != *destination
                })
                .collect::<HashMap<(usize, usize), usize>>();
            
            apply_topological_sort(new_possible_start_nodes, nodes, new_edges)
                .into_iter()
                .map(|possible_order| {
                    let mut new_order = vec![start_node];
                    new_order.extend(possible_order);
                    new_order
                })
                .collect::<Vec<Vec<usize>>>()
        })
        .collect::<Vec<Vec<usize>>>()


    // if edges.keys().len() == 0 {
    //     let missing_nodes = nodes.clone()
    //         .into_iter()
    //         .filter(|node| {
    //             !resulting_order.contains(&node)
    //         })
    //         .collect::<Vec<usize>>();

    //     for missing_node in missing_nodes {
    //         let possible_spots = nodes.len() - resulting_order.len();

    //         *resulting_confidence *= 1.0 / possible_spots as f64;
    //         resulting_order.push(missing_node);
    //     }

    //     return
    // }
    
    // let possible_start_nodes = edges.clone()
    //     .into_iter()
    //     .flat_map(|((source, destination), _)| {
    //         vec![source, destination]
    //     })
    //     .filter(|node| {
    //         nodes.iter()
    //             .filter(|source_node| {
    //                 edges.contains_key(&(**source_node, *node))
    //             })
    //             .count() == 0
    //     })
    //     .collect::<Vec<usize>>();

    // let node_to_relevancy_map = possible_start_nodes
    //     .clone()
    //     .into_iter()
    //     .map(|node| {
    //         let mut summed_node_relevancy: usize = 0;
    //         nodes.iter()
    //             .for_each(|other_node| {
    //                 if let Some(weight) = edges.get(&(node, *other_node)) {
    //                     summed_node_relevancy += *weight;
    //                 }
    //                 if let Some(weight) = edges.get(&(*other_node, node)) {
    //                     summed_node_relevancy += *weight;
    //                 }
    //             });
    //         (node, summed_node_relevancy)
    //     })
    //     .collect::<HashMap<usize, usize>>();
    // let next_node = possible_start_nodes
    //     .clone()
    //     .into_iter()
    //     .max_by(|node_a, node_b| {
    //         let relevancy_a = node_to_relevancy_map.get(node_a).unwrap();
    //         let relevancy_b = node_to_relevancy_map.get(node_b).unwrap();
    //         relevancy_a.cmp(relevancy_b)
    //     })
    //     .unwrap();
    
    // let total_relevancy = possible_start_nodes
    //     .into_iter()
    //     .map(|node| {
    //         node_to_relevancy_map.get(&node).unwrap()
    //     })
    //     .sum::<usize>();
    // let next_node_relevancy = *node_to_relevancy_map.get(&next_node).unwrap();

    // resulting_order.push(next_node);
    // *resulting_confidence *= next_node_relevancy as f64 / total_relevancy as f64;

    // *edges = edges.clone()
    //     .into_iter()
    //     .filter(|((source, destination), _)| {
    //         *source != next_node
    //         && *next_node != *destination
    //     })
    //     .collect::<HashMap<(usize, usize), usize>>();

    // apply_topological_sort(nodes, edges, resulting_order, resulting_confidence);
}

fn cycle_order_to_node_order(cycle_order: &Vec<usize>, raw_cycles: &Vec<Vec<u64>>) -> Vec<u64> {
    cycle_order.clone()
        .into_iter()
        .flat_map(|cycle_index| {
            raw_cycles.get(cycle_index).unwrap()
        })
        .map(|x| *x)
        .collect::<Vec<u64>>()
}

fn get_sequence_from_node_order(graph: &Graph, node_order: &Vec<u64>) -> String {
    let mut sequence: String = graph.nodes.get(node_order.first().unwrap()).unwrap().to_string();
    
    for node_id in node_order.into_iter().skip(1) {
        let seq_of_node = graph.nodes.get(node_id).unwrap();
        sequence.push(seq_of_node.chars().last().unwrap());
    }

    sequence
}

// Almost good enough
// Remember having all cycles in the path is another constraint
// so check the filters, so that they value unique nodes (not in the primitive constraints) as extra
// BEFORE doing anything, implement a function that performs a constraint check on a possible path
// using this path, you can test the solution 5, 7, 2, 1, 0, 4, 8, 6, 3 (which is valid under the fully reduced constraints)
// IF it is valid and the reductions all are valid, the paths can be considered correct and additional heuristical or other approaches are required to get the perfect path
pub fn assembly(graph: Graph, reads: Reads, raw_cycles: Vec<Vec<u64>>, debug: bool) -> String {
    let node_to_cycle: HashMap<u64, usize> = get_node_to_cycle_map(&raw_cycles);
    println!("Kept cycles: {:?} with a total of {} cycles",
        node_to_cycle.clone().into_iter()
            .map(|(_, index)| index)
            .collect::<HashSet<usize>>(),
        raw_cycles.len());
    let constraints = generate_constraints(&graph, reads, node_to_cycle);
    let all_cycle_orders = sort_topologically(&constraints);
    let overall_confidence: f64 = 1.0 / all_cycle_orders.len() as f64;

    let sequences = all_cycle_orders.clone().into_iter()
        .map(|cycle_order| {
            let node_order = cycle_order_to_node_order(&cycle_order, &raw_cycles);
            let sequence = get_sequence_from_node_order(&graph, &node_order);
            sequence
        })
        .collect::<Vec<String>>();
    let exact_sequence = sequences.first().unwrap()
        .chars()
        .enumerate()
        .map(|(index, character)| {
            if sequences.iter()
                .map(|sequence| sequence.chars().nth(index).unwrap())
                .all(|other_character| other_character == character)
            {
                character
            } else {
                'N'
            }
        })
        .collect::<String>();

    println!("TopoSort with confidence {:.2}% and order of {:?}", overall_confidence * 100.0, all_cycle_orders.first().unwrap());
    
    exact_sequence
}
