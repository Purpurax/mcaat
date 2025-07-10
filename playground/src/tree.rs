use std::collections::{HashMap, HashSet};


pub struct Tree {
    pub roots: Vec<u64>,
    pub nodes: Vec<(u64, Vec<usize>)>,
    pub node_location_map: HashMap<u64, Vec<usize>>
}

enum AddPathStage {
    Unclear,
    NewNodes,
    ContinueAtNode(Vec<usize>),
    TrenchNewPath(Vec<usize>)
}

impl Tree {
    pub fn new() -> Tree {
        Tree {
            roots: vec![],
            nodes: vec![],
            node_location_map: HashMap::new()
        }
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
        self.roots.push(node_id);
        self.nodes.push((node_id, vec![]));

        let new_location: usize = self.nodes.len() - 1;
        self.add_location_map_entry(node_id, new_location);
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

    pub fn add_node(&mut self, parent_location: Option<usize>, node_id: u64) -> usize {
        if let Some(parent_loc) = parent_location {
            self.add_normal_node(parent_loc, node_id)
        } else {
            self.add_root(node_id)
        }
    }

    pub fn add_ordered_path(&mut self, path: &Vec<u64>) {
        let mut add_path_stage = AddPathStage::Unclear;

        for node in path {
            add_path_stage = match add_path_stage {
                AddPathStage::Unclear =>
                    self.add_path_node_in_stage_unclear(*node),
                AddPathStage::NewNodes =>
                    self.add_path_node_in_stage_new_nodes(*node),
                AddPathStage::ContinueAtNode(locations) =>
                    self.add_path_node_in_stage_continue_at_node(*node, locations),
                AddPathStage::TrenchNewPath(locations) =>
                    self.add_path_node_in_stage_trench_new_path(*node, locations),
            }
        }
    }

    fn add_path_node_in_stage_unclear(&mut self, node: u64) -> AddPathStage {
        if let Some(location) = self.node_location_map.get(&node) {
            AddPathStage::ContinueAtNode(location.clone())
        } else {
            self.add_node(None, node);
            AddPathStage::NewNodes
        }
    }

    fn add_path_node_in_stage_new_nodes(&mut self, node: u64) -> AddPathStage {
        
        0
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
