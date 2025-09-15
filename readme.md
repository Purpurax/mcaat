## 😸 MCAAT - Metagenomic CRISPR Array Analysis Tool

- CRISPR-Cas is a bacterial immune system also famous for its use in genome editing. The diversity of known systems could be significantly increased by metagenomic data. 
- Here we present the Metagenomic CRISPR Array Analysis Tool MCAAT, a highly sensitive algorithm for finding CRISPR Arrays in un-assembled metagenomic data. 
- It takes advantage of the properties of CRISPR arrays that form multicycles in de Bruijn graphs. 
- MCAAT's assembly-free graph-based strategy outperforms assembly-based workflows and other assembly-free methods on synthetic and real metagenomes. 
---

### 🥳 NEWS
- Docker container available under: https://hub.docker.com/r/feeka94/mcaat
- Version 0.3 makes use of following optimization techniques:
  - Better data structures for preprocessing, `phmap::flat_hash_set`
  - Added compiler intrinsics to guide the hardware in the right direction
  - Reserving the capacity to prevent rehashing
In depth technical details: [educational resource](./src/z_educational_guide.md) and [optimization developer notes](./src/z_optimization_dev_notes.md). As a result of the above optimizations we achieved __17-25__ times speedup in __1billion__ node graph(from 3 <span style="color:red">days</span> to 3 <span style="color:green">hours</span>). Considering the complexity of the graphs, this is a huge improvement.


### Installation using docker
#### Docker Build

```bash
docker build -t mcaat .
```

---

#### Run the Tool Using Docker

Mount your working directory to access input/output files:

```bash
docker run --rm -v $(pwd):/data mcaat \
  --input_files /data/reads_R1.fastq /data/reads_R2.fastq \
  --output-folder /data/results
```

---

#### Final Image Size

The final image is based on `debian:bookworm-slim` and includes only:

- The `mcaat` binary
- Runtime libraries: `libomp5`, `zlib1g`

This keeps the image small and portable.

---

#### Clean Up

To remove the image:

```bash
docker rmi mcaat
```

### Compiling the project

#### 🔧 Build the Project
To allow ./install.sh make changes, we execute following command:
```bash
chmod +x ./install.sh
```
You can build the project and the working version will be saved in the build folder.
```bash
./install.sh
```
It is also possible to install the library by simply putting the --install flag.
```bash
./install.sh --install
```
To clean up you can use --clean flag.

---


### Usage

```bash
./mcaat --input-files <file1> [file2] [--ram <amount>] [--threads <num>] [--output-folder <path>] [--help]
```

---
### 🧾 Command-Line Arguments

#### ✅ Required

| Argument                  | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| `--input_files <file1> [file2]` | One or two input FASTA/FASTQ files. If one file is provided, it is treated as single-end data. If two files are provided, they are treated as paired-end reads. |

#### ⚙️ Optional

| Argument                  | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| `--ram <amount>`          | Maximum RAM to use. Units: `B`, `K`, `M`, `G`. <br>**Default:** 95% of system RAM <br>**Example:** `--ram 4G` |
| `--threads <num>`         | Number of threads to use. <br>**Default:** total CPU cores minus 2          |
| `--output-folder <path>`  | Output directory for results. <br>If not provided, a timestamped folder will be created automatically. If provided, the folder is used exactly as given. |
| `--help`, `-h`            | Show usage information and exit                                            |

---

### 📁 Output Structure

The tool creates the following directory structure inside the specified output folder:

```
<output-folder>/
├── CRISPR_Arrays.txt         # Raw CRISPR array output
```

---

### 🧪 Example Usage

| Scenario                     | Command                                                                 |
|-----------------------------|-------------------------------------------------------------------------|
| Paired-end input with custom output | `./mcaat --input_files reads_R1.fastq reads_R2.fastq --ram 8G --threads 12 --output-folder results/my_run` |
| Single-end input with default output | `./mcaat --input_files reads.fastq` <br>Creates a folder like `mcaat_run_2025-07-07_15-30-00/` |


---

#### Notes

- Input files must exist and be accessible.
- If RAM is set below 1 GB or above system capacity, the program will exit with an error.
- If only one input file is provided, the tool assumes single-end data.

---

#### Requirements

- C++17 compiler
- [RapidFuzz](https://github.com/maxbachmann/rapidfuzz-cpp) (for fuzzy string matching)
- Filesystem support (`<filesystem>`)

---

## Support

If you encounter issues or have questions, feel free to open an issue or write us an email: fikrat.talibli@ibmg.uni-stuttgart.de. If you are using this software please cite this paper: https://academic.oup.com/microlife/article/doi/10.1093/femsml/uqaf016/8205558.
