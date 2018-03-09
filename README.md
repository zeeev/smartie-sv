# smartie-sv

Have contigs and a reference genome? Smartie-sv will align query contigs against a reference genome and call structural variants.  

### Dependencies

1. Anaconda
2. Snakemake
3. Bedtools

### Running

Clone the repo:
``` 
git clone --recursive https://github.com/zeeev/smartie-sv.git
```
Export HDF5LIB if it is NOT in a global path. If the build fails check this first.

Run make:
```
cd smartie-sv && make
```
Edit the config in the pipeline folder

Index the reference genome:

```
bin/sawriter target.fasta
```
Run the pipeline (Sun Grid Engine): 
```
snakemake -j 10 --cluster-config cluster.config.json --cluster "qsub -l {cluster.h_rt} -l {cluster.mfree} -pe {cluster.pe} -q {cluster.q}" -s Snakefile --verbose -p
```

Run the pipeline (SLURM):

```
snakemake -p -w 25 -j 3 --verbose -s Snakefile --cluster-config cluster.config.slurm.json --cluster "sbatch -p {cluster.partition} -n {cluster.n}  -t {cluster.time} -c {cluster.c} -o out -e error" -w 30
```

Run the pipeline on a local machine:

```
snakemake -s Snakefile -w 50  -p -k -j 20
```
 
