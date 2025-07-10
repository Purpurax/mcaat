use std::collections::HashSet;
use std::{collections::HashMap, fs::File, io};
use std::io::Write;

use crate::cycle::Cycle;
use crate::reads::Reads;
use crate::tree::Tree;


pub struct Graph {
    pub nodes: HashMap<u64, String>,
    pub edges: HashMap<u64, Vec<u64>>,
    pub crispr_nodes: HashSet<u64>
}

impl Graph {
    pub fn new() -> Graph {
        Graph {
            nodes: HashMap::new(),
            edges: HashMap::new(),
            crispr_nodes: HashSet::new()
        }
    }

    /// Assumes the nodes have ids from 0 to n
    pub fn parse(edges_file_path: String, nodes_file_path: String) -> Graph {
        let edges_file_content = std::fs::read_to_string(&edges_file_path)
            .expect("Failed to read edges file");
        let nodes_file_content = std::fs::read_to_string(&nodes_file_path)
            .expect("Failed to read nodes file");

        let nodes: HashMap<u64, String> = nodes_file_content.split("\n")
            .filter_map(|line| {
                let line_splitted = line.split(": ").collect::<Vec<&str>>();

                if let Ok(key) = line_splitted.get(0).unwrap().parse::<u64>() {
                    if let Some(value) = line_splitted.get(1) {
                        return Some((key, value.to_string()))
                    }
                }
                None
            }).collect::<HashMap<u64, String>>();
        let edges: HashMap<u64, Vec<u64>> = edges_file_content.split("\n")
            .map(|line| {
                let splitted: Vec<&str> = line.split(":").collect();
                let key: u64 = splitted.get(0)
                    .unwrap_or(&"")
                    .parse()
                    .unwrap_or(u64::MAX);
                let mut values_str: &str = splitted.get(1)
                    .unwrap_or(&"[]");
                if values_str.len() >= 2 {
                    values_str = {
                        let mut chars = values_str.chars();
                        chars.next();
                        chars.next_back();
                        chars.as_str()
                    };
                } else {
                    values_str = "";
                }

                let values: Vec<u64> = values_str.split(",")
                    .map(|to_node_id_str| {
                        to_node_id_str.parse().unwrap_or(u64::MAX)
                    }).collect::<Vec<u64>>();

                (key, values)
            }).collect::<HashMap<u64, Vec<u64>>>();

        Graph {
            nodes,
            edges,
            crispr_nodes: HashSet::new()
        }
    }

    pub fn get_unique_node_seq(&self) -> HashSet<String> {
        let unique_ids = self.nodes.keys()
            .map(|key| *key)
            .collect::<HashSet<u64>>()
            .into_iter()
            .collect::<Vec<u64>>();

        unique_ids.into_iter()
            .map(|id| {
                self.nodes.get(&id).unwrap().clone()
            })
            .collect::<HashSet<String>>()
    }

    pub fn keep_crispr_cycles(&self, cycles: Vec<Cycle>) -> Graph {
        let unique_cycle_nodes: HashSet<u64> = cycles.into_iter()
            .flat_map(|cycle| cycle.into_iter())
            .collect::<HashSet<u64>>();

        let new_edges = self.edges.clone()
                .into_iter()
                .filter(|edge|
                    unique_cycle_nodes.contains(&edge.0)
                    || edge.1.iter().any(|end_node| unique_cycle_nodes.contains(end_node))
                ).collect::<HashMap<u64, Vec<u64>>>();

        let all_crispr_nodes = new_edges.iter()
                .map(|edge| {
                    [*edge.0].into_iter().chain(
                        edge.1.iter().map(|val| *val))
                }).flatten()
                .collect::<HashSet<u64>>();

        let new_nodes = self.nodes.clone()
                .into_iter()
                .filter(|(node_id, _node_seq)| {
                    all_crispr_nodes.contains(node_id)
                })
                .collect::<HashMap<u64, String>>();

        Graph {
            nodes: new_nodes,
            edges: new_edges,
            crispr_nodes: all_crispr_nodes
        }
    }

    pub fn get_strongly_connected_components(&self) -> Vec<Graph> {
        let mut visited: HashSet<u64> = HashSet::new();
        let mut groups: Vec<HashSet<u64>> = vec![];

        for &start_node in self.edges.keys() {
            if visited.contains(&start_node) {
                continue;
            }

            let mut group: HashSet<u64> = HashSet::new();
            let mut walking_nodes: Vec<u64> = vec![start_node];

            while let Some(walking_node) = walking_nodes.pop() {
                if visited.contains(&walking_node) {
                    continue;
                }
                
                visited.insert(walking_node);
                group.insert(walking_node);

                if let Some(landing_nodes) = self.edges.get(&walking_node) {
                    for &landing_node in landing_nodes {
                        if !visited.contains(&landing_node) {
                            walking_nodes.push(landing_node);
                        }
                    }
                }
            }

            if !group.is_empty() {
                groups.push(group);
            }
        }

        groups.into_iter().map(|group: HashSet<u64>| {
            group.into_iter().collect::<Vec<u64>>()
        })
        .filter(|group_nodes: &Vec<u64>| group_nodes.len() > 1)
        .map(|group_nodes: Vec<u64>| {
            let mut graph: Graph = Graph::new();
            for node in group_nodes {
                if let Some(node_name) = self.nodes.get(&node) {
                    graph.nodes.insert(node, node_name.clone());
                }
                if let Some(node_edges) = self.edges.get(&node) {
                    graph.edges.insert(node, node_edges.clone());
                }
            }
            graph
        }).collect::<Vec<Graph>>()
    }

    pub fn extend_by_m_steps(&mut self, m: usize, full_graph: &Graph) {
        let mut nodes_to_check: HashSet<u64> = self.edges.keys()
            .copied().collect();

        let mut incoming_into_relevant: HashSet<u64> = HashSet::new();
        let mut outgoing_from_relevant: HashSet<u64> = HashSet::new();
        
        let mut counter = 0;
        while counter < m && !nodes_to_check.is_empty() {
            // Check specified node for extensions
            for node in nodes_to_check.drain() {
                if let Some(to_nodes) = full_graph.edges.get(&node) {
                    for to_node in to_nodes {
                        if !self.edges.contains_key(to_node) {
                            outgoing_from_relevant.insert(*to_node);
                        }
                    }
                }
                
                // BRUTE FORCE
                for (from_node, to_nodes) in &full_graph.edges {
                    if to_nodes.contains(&node)
                    && !self.edges.contains_key(from_node) {
                        incoming_into_relevant.insert(*from_node);
                    }
                }
            }

            // Add the extensions as the new nodes_to_check for further extensions
            nodes_to_check.extend(incoming_into_relevant.drain());
            nodes_to_check.extend(outgoing_from_relevant.drain());

            // Add as new edges to relevant graph
            for node in nodes_to_check.iter() {
                if let Some(edges) = full_graph.edges.get(node) {
                    self.edges.insert(*node, edges.clone());
                }
            }

            counter += 1;
        }
    }

    fn find_path(&self, start: u64, end: u64, node_count: usize) -> Vec<Vec<u64>> {
        if node_count == 0 && start == end {
            return vec![vec![end]]
        } else if node_count == 0 {
            return vec![]
        }
        
        // Safe handling of missing nodes
        let Some(neighbors) = self.edges.get(&start) else {
            return vec![];
        };
        
        neighbors.iter()
            .flat_map(|&next_node| {
                let mut sub_paths = self.find_path(next_node, end, node_count - 1);
                for sub_path in &mut sub_paths {
                    sub_path.insert(0, start);
                }
                sub_paths
            })
            .collect()
    }

    pub fn reconstruct_crispr_sequence(&self, reads: Reads) -> String {
        let sequence_to_node_id_map: HashMap<String, u64> = self.nodes.iter()
            .map(|(node_id, sequence)| (sequence.clone(), *node_id))
            .collect();

        let mut order_tree: Tree = Tree::new();

        for read in reads.reads {
            let start_node_id: u64 = *sequence_to_node_id_map.get(&read.start_k_mer).unwrap();
            let end_node_id: u64 = *sequence_to_node_id_map.get(&read.end_k_mer).unwrap();
            let paths: Vec<Vec<u64>> = self.find_path(start_node_id, end_node_id, read.nodes_between);
        
            if paths.len() == 0 {
                continue
            } else if paths.len() == 1 && paths.first().unwrap().len() > 0 {
                let path = paths.first().unwrap();

                order_tree.add_ordered_path(path);
            } else {
                // TODO can still be used for ordering, by seeking similarities in the paths
                continue
            }
        }

        "".to_string()
    }

    pub fn export_to_dot(&self, path: &str) -> io::Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "digraph G {{")?;

        // Example: Assuming self.nodes is a collection of node IDs
        // for node_id in self.iter_nodes() {
        //     writeln!(file, "  \"{}\";", node_id)?;
        // }

        // Example: Assuming self.adj is an adjacency list like HashMap<u32, Vec<u32>>
        for (source, destinations) in &self.edges {
            for destination in destinations {
                writeln!(file, "  \"{}\" -> \"{}\";", source, destination)?;
            }
        }

        writeln!(file, "}}")?;
        Ok(())
    }
}


