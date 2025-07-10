use std::collections::HashSet;
use std::iter;
use std::{collections::HashMap, fs::File, io};
use std::io::Write;

use crate::cycle::Cycle;


pub struct Graph {
    pub nodes: Vec<String>,
    pub edges: HashMap<u64, Vec<u64>>
}

impl Graph {
    /// Assumes the nodes have ids from 0 to n
    pub fn parse(edges_file_path: String, nodes_file_path: String) -> Graph {
        let edges_file_content = std::fs::read_to_string(&edges_file_path)
            .expect("Failed to read edges file");
        let nodes_file_content = std::fs::read_to_string(&nodes_file_path)
            .expect("Failed to read nodes file");

        let nodes: Vec<String> = nodes_file_content.split("\n")
            .map(|line| {
                line.split(" ")
                    .collect::<Vec<&str>>()
                    .get(1)
                    .unwrap_or(&"NO SEQUENCE")
                    .to_string()
            }).collect::<Vec<String>>();
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

        Graph { nodes, edges }
    }

    pub fn keep_relevant(&self, cycles: Vec<Cycle>) -> Graph {
        let unique_nodes: HashSet<u64> = cycles.into_iter()
            .flat_map(|cycle| cycle.into_iter())
            .collect::<HashSet<u64>>();

        let new_edges = self.edges.clone().into_iter()
            .filter(|edge|
                unique_nodes.contains(&edge.0)
                || edge.1.iter().any(|end_node| unique_nodes.contains(end_node))
            ).collect::<HashMap<u64, Vec<u64>>>();

        let new_nodes = self.nodes.clone();

        Graph { nodes: new_nodes, edges: new_edges }
    }

    pub fn extend_by_m_steps(&mut self, m: usize, full_graph: Graph) {
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


