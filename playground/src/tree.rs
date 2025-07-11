use std::{collections::{HashMap, HashSet}, fs::File, io};

use std::io::Write;

pub struct Tree {
    pub roots: Vec<usize>,
    pub nodes: Vec<(u64, Vec<usize>)>,
    pub node_location_map: HashMap<u64, Vec<usize>>
}

enum AddPathStage {
    Unclear,
    NewNodes(usize, usize), // root_index, location
    ContinueAtNode(Vec<usize>), // locations
    TrenchNewPath(Vec<usize>) // locations
}

impl Tree {
    pub fn new() -> Tree {
        Tree {
            roots: vec![],
            nodes: vec![],
            node_location_map: HashMap::new()
        }
    }

    pub fn add_ordered_path(&mut self, path: &Vec<u64>) {
        let mut add_path_stage = AddPathStage::Unclear;

        for node in path {
            add_path_stage = match add_path_stage {
                AddPathStage::Unclear =>
                    self.add_path_node_in_stage_unclear(*node),
                AddPathStage::NewNodes(root_index, location) =>
                    self.add_path_node_in_stage_new_nodes(*node, root_index, location),
                AddPathStage::ContinueAtNode(locations) =>
                    self.add_path_node_in_stage_continue_at_node(*node, locations),
                AddPathStage::TrenchNewPath(locations) =>
                    self.add_path_node_in_stage_trench_new_path(*node, locations),
            };
            // println!("\n---\nROOTS\n>{:?}\nNODES\n>{:?}\nMAP\n>{:?}", self.roots, self.nodes, self.node_location_map);
        }
    }

    pub fn export_to_dot(&self, path: &str) -> io::Result<()> {
        let mut file = File::create(path)?;
        
        writeln!(file, "digraph Tree {{")?;

        // Add node styling - roots with different appearance
        for (index, (node_id, _)) in self.nodes.iter().enumerate() {
            if self.roots.contains(&index) {
            writeln!(file, "  \"{}\" [label=\"{}\", shape=box, style=filled, fillcolor=lightblue];", index, node_id)?;
            } else {
            writeln!(file, "  \"{}\" [label=\"{}\"];", index, node_id)?;
            }
        }
        
        // Add edges
        for (index, (_, children)) in self.nodes.iter().enumerate() {
            for &child_index in children {
            writeln!(file, "  \"{}\" -> \"{}\";", index, child_index)?;
            }
        }

        writeln!(file, "}}")?;

        Ok(())
    }

    pub fn get_longest_path(&self) -> Vec<u64> {
        self.roots.iter()
            .map(|node_index| {
                self.get_longest_path_rec(*node_index)
            })
            .max_by(|path_a, path_b| path_a.len().cmp(&path_b.len()))
            .unwrap_or(vec![])
    }

    fn get_longest_path_rec(&self, start_index: usize) -> Vec<u64> {
        let mut path = vec![self.nodes.get(start_index).unwrap().0];
        path.extend(
            self.nodes.get(start_index)
                .unwrap().1
                .iter()
                .map(|next_index| {
                    self.get_longest_path_rec(*next_index)
                })
                .max_by(|path_a, path_b| path_a.len().cmp(&path_b.len()))
                .unwrap_or(vec![])
                .into_iter()
        );
        path
    }

    fn add_location_map_entry(&mut self, node_id: u64, location: usize) {
        let mut old_links = (*self.node_location_map
            .get(&node_id)
            .unwrap_or(&vec![])).clone();
        old_links.push(location);

        self.node_location_map.insert(
            node_id,
            old_links
        );
    }

    fn add_root(&mut self, node_id: u64) -> usize {
        self.nodes.push((node_id, vec![]));
        
        let new_location: usize = self.nodes.len() - 1;
        self.add_location_map_entry(node_id, new_location);
        self.roots.push(new_location);
        new_location
    }

    fn add_normal_node(&mut self, parent_location: usize, node_id: u64) -> usize {
        self.nodes.push((node_id, vec![]));
        
        let new_location: usize = self.nodes.len() - 1;
        self.add_location_map_entry(node_id, new_location);

        let parent = self.nodes.get_mut(parent_location).unwrap();
        parent.1.push(new_location);

        new_location
    }

    fn add_node(&mut self, parent_location: Option<usize>, node_id: u64) -> usize {
        if let Some(parent_loc) = parent_location {
            self.add_normal_node(parent_loc, node_id)
        } else {
            self.add_root(node_id)
        }
    }

    fn attach_root(&mut self, roots_index: usize, parent_location: usize) -> usize {
        let index_in_nodes_of_root: usize = self.roots.remove(roots_index);

        let parent = self.nodes.get_mut(parent_location).unwrap();
        parent.1.push(index_in_nodes_of_root);
        index_in_nodes_of_root
    }

    fn add_path_node_in_stage_unclear(&mut self, node: u64) -> AddPathStage {
        if let Some(location) = self.node_location_map.get(&node) {
            AddPathStage::ContinueAtNode(location.clone())
        } else {
            let location = self.add_node(None, node);
            let root_index: usize = self.roots.len() - 1;
            AddPathStage::NewNodes(root_index, location)
        }
    }

    fn add_path_node_in_stage_new_nodes(&mut self, node: u64, parent_root_index: usize, location: usize) -> AddPathStage {
        let matching_root: Option<usize> = self.roots.iter()
            .enumerate()
            .filter(|(_, root_location)| {
                let root_id: u64 = (*self.nodes.get(**root_location).unwrap()).0;
                root_id == node
            })
            .filter(|(in_roots_index, _)| *in_roots_index != parent_root_index)
            .map(|(in_roots_index, _)| in_roots_index)
            .next();
        if let Some(roots_index) = matching_root {
            let new_location = self.attach_root(roots_index, location);
            
            AddPathStage::ContinueAtNode(vec![new_location])
        } else {
            let new_location = self.add_node(Some(location), node);

            AddPathStage::NewNodes(parent_root_index, new_location)
        }
    }

    fn add_path_node_in_stage_continue_at_node(
        &mut self,
        node: u64,
        previous_locations: Vec<usize>
    ) -> AddPathStage {
        let mut total_future_locations: Vec<usize> = vec![];
        
        if let Some(node_locations) = self.node_location_map.get(&node) {
            // Node has some locations
            // check if any future_location starting from one of the previous_location
            //    contains a node_locations
            //  i.e. a path to that node exists and can be walked
            // if NONE of the future_locations match, a new path has to be created (see else)
            for previous_location in previous_locations {
                let future_locations: &Vec<usize> = &self.nodes.get(previous_location).unwrap().1;
                
                let matching_future_locations: Vec<usize> =
                    intersect_vectors(node_locations, future_locations);

                if matching_future_locations.len() > 0 {
                    total_future_locations.extend(matching_future_locations)
                }
            }

            if total_future_locations.len() > 0 {
                AddPathStage::ContinueAtNode(total_future_locations)
            } else {
                AddPathStage::TrenchNewPath(node_locations.clone())
            }
        } else {
            // Node does not exist and has to be added
            // after that the rest of the nodes will be newly added as well
            let new_locations = previous_locations.into_iter()
                .map(|previous_location| {
                    self.add_node(Some(previous_location), node)
                })
                .collect::<Vec<usize>>();
            AddPathStage::TrenchNewPath(new_locations)
        }
    }

    fn add_path_node_in_stage_trench_new_path(
        &mut self,
        node: u64,
        locations: Vec<usize>
    ) -> AddPathStage {
        let new_locations = locations.into_iter()
            .map(|location| {
                self.add_node(Some(location), node)
            })
            .collect::<Vec<usize>>();

        AddPathStage::TrenchNewPath(new_locations)
    }
}

pub fn intersect_vectors<T>(vec_1: &Vec<T>, vec_2: &Vec<T>) -> Vec<T> 
where 
    T: Eq + std::hash::Hash + Clone,
{
    let set1: HashSet<T> = vec_1.clone().into_iter().collect();
    let set2: HashSet<T> = vec_2.clone().into_iter().collect();
    
    set1.intersection(&set2)
        .cloned()
        .collect()
}
