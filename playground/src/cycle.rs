use std::collections::HashSet;

pub type Cycle = Vec<u64>;

pub fn parse(file_path: String) -> Vec<Cycle> {
    let content: String = std::fs::read_to_string(&file_path)
        .expect("Failed to read file for cycles");

    content.split("\n")
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .map(|line| {
            line.split(" ")
                .map(|id_str| id_str.parse::<u64>().unwrap_or(u64::MAX))
                .collect::<Cycle>()
        }).collect::<Vec<Cycle>>()
}

pub fn export_as_desired_input(cycles: Vec<Cycle>, file_path: String) {
    let mut all_nodes: Vec<u64> = cycles.into_iter()
        .flat_map(|cycle| cycle.into_iter())
        .collect::<HashSet<u64>>()
        .into_iter()
        .collect();
    all_nodes.sort();

    let content: String = all_nodes.into_iter()
        .map(|distinct_cycle| {
            distinct_cycle.to_string()
        }).fold(String::new(), |a, b| a + " " + &b)
        .trim()
        .to_string();

    let _ = std::fs::write(file_path, content);
}
