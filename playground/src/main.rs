use cycle::Cycle;
use graph::Graph;
use reads::Reads;

use clap::Parser;

pub mod cycle;
pub mod graph;
pub mod reads;
pub mod tree;

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

    let cycles: Vec<Cycle> = cycle::parse(args.cycles);

    let all_reads = Reads::parse(args.reads, args.reads_2.unwrap());
    let m: usize = all_reads.reads.first().map_or(0, |read| read.nodes_between as usize + 1);

    let result_folder: String = "./results/".to_string();

    let full_graph: Graph = Graph::parse(args.graph_structure, args.graph_nodes);
    let relevant_graph: Graph = full_graph.keep_crispr_cycles(cycles.clone());
    let mut crispr_subgraphs: Vec<Graph> = relevant_graph.get_strongly_connected_components();
    crispr_subgraphs.iter_mut().for_each(|sub_graph: &mut Graph| {
        sub_graph.extend_by_m_steps(m, &full_graph)
    });

    let problem_cases = crispr_subgraphs.into_iter()
        .map(|subgraph| {
            let relevant_reads: Reads = all_reads.get_relevant(&subgraph);
            (subgraph, relevant_reads)
        });
    
    for (i, (subgraph, reads)) in problem_cases.enumerate() {
        if args.output {
            let folder = format!("{}{}{}", result_folder.clone(), "crispr_", i);
            let crispr_sequence: String = subgraph.reconstruct_crispr_sequence(reads, Some(folder));
            println!("{}: {}", i, crispr_sequence);
        } else {
            let crispr_sequence: String = subgraph.reconstruct_crispr_sequence(reads, None);
            println!("{}: {}", i, crispr_sequence);
        }
    }

    if args.output {
        all_reads.export_as_desired_input(result_folder.clone() + "/desired-reads-output.csv");
        cycle::export_as_desired_input(cycles, result_folder.clone() + "/desired-cycles-output.csv");

        full_graph.export_to_dot(&(result_folder.clone() + "full-graph.dot")).expect("Failed outputing full graph into dot file");
        relevant_graph.export_to_dot(&(result_folder.clone() + "all-crispr-arrays.dot")).expect("Failed outputing only crispr graph into dot file");
        
        for (i, mut subgraph) in relevant_graph.get_strongly_connected_components().into_iter().enumerate() {
            subgraph.export_to_dot(&(result_folder.clone() + format!("one-crispr-array-{}.dot", i).as_str())).expect("Failed outputing relevant graph into dot file");
    
            subgraph.extend_by_m_steps(m, &full_graph);
            subgraph.export_to_dot(&(result_folder.clone() + format!("one-crispr-array-extended-{}.dot", i).as_str())).expect("Failed outputing relevant graph into dot file");
        }
    }
}
