use cycle::Cycle;
use graph::Graph;
use reads::Reads;

use clap::Parser;

pub mod graph;
pub mod reads;
pub mod cycle;

/// Program which basically filters out only relevant data for CRISPR array ordering
#[derive(Parser, Debug)]
struct Args {
    /// Graph structure file path
    #[arg(long, default_value = "/home/master/Documents/UNI/Informatik/Semester-4/Bachelor/mcaat/data/graph/graph_structure.txt")]
    graph_structure: String,

    /// Graph nodes file path
    #[arg(long, default_value = "/home/master/Documents/UNI/Informatik/Semester-4/Bachelor/mcaat/data/graph/nodes.txt")]
    graph_nodes: String,

    /// READS file path
    #[arg(long, default_value = "/home/master/Documents/UNI/Informatik/Semester-4/Bachelor/mcaat/data/reads_input/basic_read_R1.fastq")]
    reads: String,

    /// Second READS file path (optional)
    #[arg(long, default_value = "/home/master/Documents/UNI/Informatik/Semester-4/Bachelor/mcaat/data/reads_input/basic_read_R2.fastq")]
    reads_2: Option<String>,

    /// File path of the cycles
    #[arg(long, default_value = "/home/master/Documents/UNI/Informatik/Semester-4/Bachelor/mcaat/data/cycles/cycles.txt")]
    cycles: String,

    /// Output dot file for visualization of the graph
    #[arg(long, action = clap::ArgAction::SetTrue)]
    output: bool
}

fn main() {
    let args = Args::parse();

    let graph = Graph::parse(args.graph_structure, args.graph_nodes);

    let reads = Reads::parse(args.reads, args.reads_2.unwrap());

    let cycles: Vec<Cycle> = cycle::parse(args.cycles);

    println!("{:?}", reads);
    println!("\n\n\n\n{:?}", cycles);

    if args.output {
        graph.export_to_dot("output.dot").expect("Failed outputing graph into dot file");
    }
}
