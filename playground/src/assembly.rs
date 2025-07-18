use std::{cmp::Ordering, collections::{HashMap, HashSet}, fs::File, io};
use std::io::Write;

use crate::{graph::Graph, reads::Reads};
use rand::seq::SliceRandom;

fn get_node_to_cycle_map(raw_cycles: &Vec<Vec<u64>>) -> HashMap<u64, Vec<usize>> {
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

fn generate_constraints(graph: &Graph, reads: Reads, raw_cycles: &Vec<Vec<u64>>) -> HashSet<(Vec<usize>, Vec<usize>)> {
    let node_to_cycle: HashMap<u64, Vec<usize>> = get_node_to_cycle_map(raw_cycles);

    println!("nodes_to_cycles_map: {:?}", node_to_cycle);

    let mut constraints: HashSet<(Vec<usize>, Vec<usize>)> = HashSet::new();

    for (start_node_id, middle_node_id, end_node_id, node_in_between)
    in map_reads_to_node_id_pairs(reads, graph) {
        let paths = graph.find_path(start_node_id, middle_node_id, end_node_id, node_in_between);
        println!("{}", paths.len());
        let cycles_on_path = paths.first()
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
    
        // get every combination and store the constraint in hashset
        cycles_on_path.clone()
            .into_iter()
            .filter(|cycles| cycles.len() > 0)
            .enumerate()
            .for_each(|(i, cycles)| {
                cycles_on_path.clone()
                    .into_iter()
                    .filter(|other_cycles| other_cycles.len() > 0)
                    .enumerate()
                    .filter(|(j, _)| i < *j)
                    .for_each(|(_, other_cycles)| {
                        let mut left_side = cycles.clone();
                        let mut right_side = other_cycles.clone();
                        left_side.sort();
                        right_side.sort();

                        constraints.insert((left_side, right_side));
                    });
            });
    }
    
    constraints
}

fn filter_satisfied_constraints(constraints: &mut HashSet<(Vec<usize>, Vec<usize>)>) {
    *constraints = constraints.clone().into_iter()
        .filter(|(from_cycles, to_cycles)| {
            let from = from_cycles.clone().into_iter().collect::<HashSet<usize>>();
            let to = to_cycles.clone().into_iter().collect::<HashSet<usize>>();

            from.intersection(&to).count() == 0
        })
        .collect::<HashSet<(Vec<usize>, Vec<usize>)>>();
}

fn filter_loose_constraints(constraints: &mut HashSet<(Vec<usize>, Vec<usize>)>) {
    let indexable_constraints = constraints.clone().into_iter().collect::<Vec<(Vec<usize>, Vec<usize>)>>();

    let mut left_cycles_occurances: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    indexable_constraints.clone().into_iter()
        .enumerate()
        .for_each(|(i, (left_cycles, _right_cycles))| {
            if left_cycles_occurances.contains_key(&left_cycles) {
                left_cycles_occurances.get_mut(&left_cycles).unwrap().push(i);
            } else {
                left_cycles_occurances.insert(left_cycles, vec![i]);
            }
        });

    let mut constraints_to_delete: HashSet<usize> = HashSet::new();
    
    left_cycles_occurances.clone()
        .into_iter()
        .for_each(|(_left_cycles, indices_of_occurance)| {
            // go through every index and check whether one is a loose constraint
            indices_of_occurance.clone().into_iter()
                .for_each(|index| {
                    let index_right_side = indexable_constraints.get(index).unwrap().1.clone();

                    let is_loose_constraint = indices_of_occurance.clone()
                        .into_iter()
                        .filter(|other_index| index != *other_index)
                        .map(|other_index| {
                            indexable_constraints.get(other_index).unwrap().1.clone()
                        }) // other_indexes right side
                        .filter(|other_right_side| other_right_side.len() < index_right_side.len())
                        .any(|other_right_side| {
                            let index_right_side_hash = index_right_side.clone().into_iter().collect::<HashSet<usize>>();
                            let other_right_side_hash = other_right_side.clone().into_iter().collect::<HashSet<usize>>();

                            index_right_side_hash.intersection(&other_right_side_hash).count()
                            == other_right_side.len()
                            // other right side has as many elements as the intersection
                            //  ergo other is an index for a stricter constraint
                        });
                    if is_loose_constraint {
                        constraints_to_delete.insert(index);
                    }
                })
        });

    constraints_to_delete.into_iter()
        .for_each(|index| {
            let constraint = indexable_constraints.get(index).unwrap();
            constraints.remove(constraint);
        });
}

fn display_constraints(con: &HashSet<(Vec<usize>, Vec<usize>)>) {
    let mut constraints = con.clone().into_iter().collect::<Vec<(Vec<usize>, Vec<usize>)>>();
    constraints.sort_by(|(left_a, right_a), (left_b, right_b)| {
        let ordering = left_a.len().cmp(&left_b.len());
        if ordering == Ordering::Equal {
            right_a.len().cmp(&right_b.len())
        } else {
            ordering
        }
    });

    println!("Displaying constraints: {:?}", constraints);
}

fn filter_transitively_satisfied_constraints(constraints: &mut HashSet<(Vec<usize>, Vec<usize>)>) {
    let mut constraints_to_delete: HashSet<usize> = HashSet::new();
    let mut delete_at_least_one: bool = true;
        
    while delete_at_least_one {
        delete_at_least_one = false;

        let mut primitive_constraints: HashMap<usize, usize> = HashMap::new();

        constraints.clone().into_iter()
            .for_each(|(left, right)| {
                if left.len() == 1 && right.len() == 1 {
                    primitive_constraints.insert(left.first().unwrap().clone(), right.first().unwrap().clone());
                }
            });
        let indexable_constraints = constraints.clone().into_iter().collect::<Vec<(Vec<usize>, Vec<usize>)>>();

        indexable_constraints.clone().into_iter()
            .enumerate()
            .filter(|(_, (left, right))| left.len() > 1 || right.len() > 1)
            .for_each(|(i, (left, right))| {
                let transitive_duplicate = left.into_iter()
                    .any(|left_element| {
                        let mut visited: HashSet<usize> = HashSet::new();
                        let mut current_index: usize = left_element;

                        while !visited.contains(&current_index) {
                            visited.insert(current_index);

                            if let Some(next_index) = primitive_constraints.get(&current_index) {
                                current_index = *next_index
                            }
                        }

                        let right_hashed = right.clone().into_iter().collect::<HashSet<usize>>();

                        visited.intersection(&right_hashed).count()
                        == right_hashed.len()
                    });
                
                if transitive_duplicate {
                    constraints_to_delete.insert(i);
                    delete_at_least_one = true;
                }
            });
        
        constraints_to_delete.drain()
            .for_each(|index| {
                let constraint = indexable_constraints.get(index).unwrap();
                constraints.remove(constraint);
            });
    }
}

fn get_amount_of_violated_contraints(
    path: Vec<usize>,
    constraints: &HashSet<(Vec<usize>, Vec<usize>)>
) -> usize {
    constraints.clone().into_iter()
        .filter(|(left, right)| {
            !left.into_iter()
                .filter(|left_index| path.contains(&left_index))
                .any(|left_index| {
                    right.clone().into_iter()
                        .filter(|right_index| {
                            path.contains(right_index)
                        })
                        .any(|right_index| {
                            let path_index_of_left_cycle = path.iter().position(|ele| ele == left_index);
                            let path_index_of_right_cycle = path.iter().position(|ele| *ele == right_index);

                            path_index_of_left_cycle <= path_index_of_right_cycle
                        })
                })
        })
        .count()
}

// Almost good enough
// Remember having all cycles in the path is another constraint
// so check the filters, so that they value unique nodes (not in the primitive constraints) as extra
// BEFORE doing anything, implement a function that performs a constraint check on a possible path
// using this path, you can test the solution 5, 7, 2, 1, 0, 4, 8, 6, 3 (which is valid under the fully reduced constraints)
// IF it is valid and the reductions all are valid, the paths can be considered correct and additional heuristical or other approaches are required to get the perfect path
pub fn assembly(graph: Graph, reads: Reads, raw_cycles: Vec<Vec<u64>>, debug: bool) -> String {
    let mut constraints = generate_constraints(&graph, reads, &raw_cycles);

    filter_satisfied_constraints(&mut constraints);
    filter_loose_constraints(&mut constraints);
    filter_transitively_satisfied_constraints(&mut constraints);

    let mut rng = rand::rng();
    let mut cycle_indices: Vec<usize> = (0..raw_cycles.len()).collect();
    let mut satisfying_permutations = 0;
    let mut first_working_permutation: Option<Vec<usize>> = None;
    let sample_size = 1000;
    
    for _ in 0..sample_size {
        cycle_indices.shuffle(&mut rng);
        
        if get_amount_of_violated_contraints(cycle_indices.clone(), &constraints) == 0 {
            satisfying_permutations += 1;
            if first_working_permutation.is_none() {
                first_working_permutation = Some(cycle_indices.clone());
            }
        }
    }
    
    if let Some(first_permutation) = first_working_permutation {
        println!("First working permutation found: {:?}", first_permutation);
        println!("Satisfying permutations: {}/{} samples", satisfying_permutations, sample_size);
    } else {
        println!("No satisfying permutation found in {} samples", sample_size);
    }

    display_constraints(&constraints);


    "".to_string()
}
