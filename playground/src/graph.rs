use std::{collections::HashMap, fs::File, io};
use std::io::Write;


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


