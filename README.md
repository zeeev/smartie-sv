# smartie-sv

Have contigs and a reference genome? Smartie-sv will align query contigs against a reference genome and call structural variants.  

### Dependancies

1. Ananconda
2. Snakemake
3. Bedtools

### Running

1. Edit the config.json.
2. Edit the config.sh (you need to export the path to hdf5lib, see the example in config.sh)
2. If you don't have "modules" make sure your path contains bedtools. 
3. Index target genomes with Blasr:

```
blasr/alignment/bin/sawriter target.fasta
```


If you're on an SGE cluster you can use the following:

```
snakemake -j 10 --cluster-config cluster.config.json --cluster "qsub -l {cluster.h_rt} -l {cluster.mfree} -pe {cluster.pe} -q {cluster.q}" -s Snakefile --verbose -p
```
 
If you're on a SLUM cluster you can use the following: 
 
 Otherwise:

 On a single node:

```
snakemake -s Snakefile -w 50  -p -k -j 20
```
 
 Or modidfy the cluster settings in the config.json for your cluster.
