use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead};
use std::path::Path;
use clap::Parser;
use clap::Subcommand;
use std::error::Error;
use std::process::{Command};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {

   #[clap(subcommand)]
   subcommand: SubCommands,

}

#[derive(Debug, Subcommand)]
enum SubCommands {

   /// After CheckV, append information from CheckV database
   /// Optionally, search PHAGE database to populate the abundance file.
   Populate {
       /// Path to abundance file (row, contig name; col, sample name)
       #[clap(short, long, value_parser)]
       abundance: String,

       /// Path to completeness.tsv from CheckV
       #[clap(short, long, value_parser)]
       completeness: String,

       /// Path to CheckV database directory, only for v1.2
       #[clap(short, long, value_parser)]
       db: String,

       /// Path to TSV file from PHAGE (Millard lab)
       #[clap(short, long, value_parser)]
       phage: Option<String>,
    },
    /// Search similar reference from PHAGE database using mash from queried contigs
    /// Mash executable needs to be in PATH.
    /// Need to provide mash database build from PHAGE database.
    Host {
        #[clap(short, long, value_parser)]
        mash_db: String,
        #[clap(short, long, value_parser)]
        phage_genome_db: Option<String>,
        #[clap(short, long, value_parser)]
        contigs: String,
        #[clap(short, long, value_parser)]
        all_info: bool,
        abundance: Option<String>,
    }
}



#[tokio::main]
async fn obtain() -> Result<String, Box<dyn Error>> {
    // println!("Obtaining through web ...");
    let res = reqwest::get("http://inphared.s3.climb.ac.uk/1Aug2022_data.tsv").await?;
    let body = res.text().await?;
    Ok(body)
}

fn main() {
    let args = Args::parse();

    match args.subcommand {
        SubCommands::Populate {abundance, completeness, db, phage} => {
            let abundance_file = abundance;
            let completeness_file = completeness;
            let db_path = db;

            let mil_contents;
            if phage.is_none() {
                mil_contents = obtain().unwrap();
            } else {
                let mil_path = phage.unwrap();
                mil_contents = fs::read_to_string(mil_path)
                    .expect("Cannot read the LR or NC data file");
            }

            // Read contig ids from abundance file
            let mut contig_ids = Vec::<String>::new();
            if let Ok(lines) = read_lines(abundance_file.clone()) {
                for line in lines {
                    if let Ok(vi) = line {
                        let split = vi.split("\t")
                                    .map(|s| s.to_string())
                                    .collect::<Vec<String>>();
                        contig_ids.push(split[0].clone());
                    }
                }
            }
            // Read completeness.tsv
            let mut contig_to_checkvdb = Vec::<&str>::new();
            let completeness_contents = fs::read_to_string(completeness_file)
                .expect("Cannot read the completeness file");
            for contig_id in contig_ids {
                for line in completeness_contents.lines() {
                    if line.contains(&contig_id) {
                        contig_to_checkvdb.push(line);
                    }
                }
            }
            
            let mut ref_ids = Vec::<String>::new();
            for line in contig_to_checkvdb {
                let ctc_line = line.split("\t")
                                   .map(|s| s.to_string())
                                   .collect::<Vec<String>>();
                ref_ids.push(ctc_line[8].clone());
            }

            let start1 = "DTR";
            let start2 = "LR";
            let start3 = "NC";
            let start4 = "GCA";

            let circular_path = db_path.to_owned() + "/genome_db/checkv_circular.tsv";
            let circular_contents = fs::read_to_string(circular_path)
                .expect("Cannot read the circular file");

            let gb_path = db_path.to_owned() + "/genome_db/checkv_genbank.tsv";
            let gb_contents = fs::read_to_string(gb_path)
                .expect("Cannot read the genbank file");

            let mut habitat = Vec::<String>::new();
            let mut lineage = Vec::<String>::new();
            let mut source = Vec::<String>::new();
            let mut ncbi_name = Vec::<String>::new();
            let mut ncbi_id = Vec::<String>::new();
            let mut vog_clade = Vec::<String>::new();
            let mut phage_desc = Vec::<String>::new();
            let mut phage_class = Vec::<String>::new();
            let mut phage_host = Vec::<String>::new();
            let mut ids = Vec::<String>::new();

            for ref_id in ref_ids {
                ids.push(ref_id.clone());
                if ref_id.starts_with(start1) {
                    for line in circular_contents.lines() {
                        if line.contains(&ref_id) {
                            let split = line.split("\t")
                                            .map(|s| s.to_string())
                                            .collect::<Vec<String>>();
                            habitat.push(split[6].clone());
                            lineage.push(split[5].clone());
                            source.push(split[1].clone());

                            phage_desc.push("".to_string());
                            phage_class.push("".to_string());
                            phage_host.push("".to_string());
                            ncbi_name.push("".to_string());
                            ncbi_id.push("".to_string());
                            vog_clade.push("".to_string());

                        }
                    }
                } else if ref_id.starts_with(start4) {
                    for line in gb_contents.lines() {
                        if line.contains(&ref_id) {
                            let split = line.split("\t")
                                            .map(|s| s.to_string())
                                            .collect::<Vec<String>>();

                            ncbi_name.push(split[3].clone());
                            ncbi_id.push(split[2].clone());
                            vog_clade.push(split[6].clone());
                            lineage.push(split[7].clone());

                            phage_desc.push("".to_string());
                            phage_class.push("".to_string());
                            phage_host.push("".to_string());
                            habitat.push("".to_string());
                            source.push("".to_string());
                        }
                    }
                } else if ref_id.starts_with(start2) || ref_id.starts_with(start3) {
                    for line in mil_contents.lines() {
                        if line.contains(&ref_id) {
                            let split = line.split("\t")
                                            .map(|s| s.to_string())
                                            .collect::<Vec<String>>();

                            phage_desc.push(split[1].clone());
                            phage_class.push(split[2].clone());
                            phage_host.push(split[14].clone());
                            lineage.push(split[19].clone()); // Parse Class

                            ncbi_name.push("".to_string());
                            ncbi_id.push("".to_string());
                            vog_clade.push("".to_string());
                            habitat.push("".to_string());
                            source.push("".to_string());
                        }
                    }
                } else {
                    ncbi_name.push("".to_string());
                    ncbi_id.push("".to_string());
                    vog_clade.push("".to_string());
                    habitat.push("".to_string());
                    source.push("".to_string());
                    phage_desc.push("".to_string());
                    phage_class.push("".to_string());
                    phage_host.push("".to_string());
                    lineage.push("".to_string());
                }
                if lineage.len() != ids.len() {
                    ncbi_name.push("".to_string());
                    ncbi_id.push("".to_string());
                    vog_clade.push("".to_string());
                    habitat.push("".to_string());
                    source.push("".to_string());
                    phage_desc.push("".to_string());
                    phage_class.push("".to_string());
                    phage_host.push("".to_string());
                    lineage.push("".to_string());
                }
            }

            ids.insert(0, "ref_id".to_string());
            ncbi_name.insert(0, "ncbi_name".to_string());
            ncbi_id.insert(0, "ncbi_id".to_string());
            vog_clade.insert(0, "vog_clade".to_string());
            habitat.insert(0, "habitat".to_string());
            source.insert(0, "source".to_string());
            phage_desc.insert(0, "phage_desc".to_string());
            phage_class.insert(0, "phage_class".to_string());
            phage_host.insert(0, "phage_host".to_string());
            lineage.insert(0, "lineage".to_string());

            // Output
            if let Ok(lines) = read_lines(abundance_file) {
                for (val, line) in lines.enumerate() {
                    let new_line = line.unwrap() + "\t" + 
                    &ids[val] + "\t" +
                    &ncbi_name[val] + "\t" + 
                    &ncbi_id[val] + "\t" + 
                    &vog_clade[val] + "\t" + 
                    &habitat[val] + "\t" + 
                    &source[val] + "\t" + 
                    &phage_desc[val] + "\t" + 
                    &phage_class[val] + "\t" + 
                    &phage_host[val] + "\t" + 
                    &lineage[val];
                    println!("{}",new_line);
                }
            }
        },
        SubCommands::Host {mash_db, contigs, phage_genome_db, abundance, all_info} => {

            let mil_contents;
            if phage_genome_db.is_none() {
                mil_contents = obtain().unwrap();
            } else {
                let mil_path = phage_genome_db.unwrap();
                mil_contents = fs::read_to_string(mil_path)
                    .expect("Cannot read phage genome db file");
            }


            let output = Command::new("mash")
            .arg("dist")
            .arg("-i")
            .arg(mash_db)
            .arg(contigs)
            .arg("-v")
            .arg("1")
            .output()
            .expect("Failed to execute mash");

            let stdout = String::from_utf8(output.stdout).unwrap();
            let mut res: HashMap<String, Vec<String>> = HashMap::new();

            for line in stdout.lines() {
                let split = line.split("\t").map(|s| s.to_string()).collect::<Vec<String>>();
                // contig_id: dist, db_id, pvalue
                if !(res.contains_key(&split[1])) {
                    res.entry(split[1].clone())
                       .or_default().push(split[2].clone());
                    res.entry(split[1].clone())
                       .or_default().push(split[0].clone());
                    res.entry(split[1].clone())
                       .or_default().push(split[3].clone());
                } else {
                    let cur_d = split[2].clone().parse::<f32>().unwrap();
                    let past_d = res.get(&split[1]).unwrap()[0].parse::<f32>().unwrap();
                    if cur_d < past_d {
                        res.remove(&split[1]);
                        res.entry(split[1].clone()).or_default().push(split[2].clone());
                        res.entry(split[1].clone()).or_default().push(split[0].clone());
                        res.entry(split[1].clone()).or_default().push(split[3].clone());
                    }
                }
            }


            if all_info {
                let first_line_all = "contig_id\tphage_id\tmash_distance\tmash_p\t".to_owned()+mil_contents.lines().nth(0).unwrap();
                println!("{}", first_line_all);
            } else {
                println!("contig_id\tphage_id\tmash_distance\tmash_p\thost");
            }

            for (k, v) in &res {
                let mut host = "".to_string();
                for line in mil_contents.lines() {
                        if line.contains(&v[1]) {
                            let split = line.split("\t")
                                            .map(|s| s.to_string())
                                            .collect::<Vec<String>>();
                            if all_info {
                                host = split[0].clone();
                                for n in 1..split.len() {
                                    host = host + "\t" + &split[n];
                                }
                            } else {
                                host = split[14].clone();
                            }
                        }
                }

                let new_line = k.to_owned() + "\t" + &v[1] + "\t" + &v[0] + "\t" + &v[2] + "\t" + &host;
                println!("{}", new_line);
            }
        }
    }
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<fs::File>>>
where P: AsRef<Path>, {
    let file = fs::File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}


