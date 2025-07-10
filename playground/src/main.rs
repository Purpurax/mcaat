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
    let m: usize = reads.reads.first().map_or(0, |read| read.nodes_between as usize + 1);

    let cycles: Vec<Cycle> = cycle::parse(args.cycles);

    // println!("{:?}", reads);
    // println!("\n\n\n\n{:?}", cycles);

    // reads.export_as_desired_input(None);
    // cycle::export_as_desired_input(cycles, None);

    if args.output {
        let result_folder: String = "./results/".to_string();

        graph.export_to_dot(&(result_folder.clone() + "output-full.dot")).expect("Failed outputing full graph into dot file");
        
        let mut relevant_graph: Graph = graph.keep_relevant(cycles.clone());
        relevant_graph.export_to_dot(&(result_folder.clone() + "output-only-crispr.dot")).expect("Failed outputing only crispr graph into dot file");
        
        relevant_graph.extend_by_m_steps(m, graph);
        relevant_graph.export_to_dot(&(result_folder + "output-relevant.dot")).expect("Failed outputing relevant graph into dot file");
    }
}
