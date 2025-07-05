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
