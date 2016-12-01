# smartie-sv

Have contigs and a reference genome?  Smartie-sv will align query contigs against a reference genome and call strucutral variants. 

### Dependancies

1. Ananconda
2. Snakemake
3. Bedtools

### Running

1. Edit the config.jason.
2. If you don't have "modules" make sure your path contains bedtools. 
3. Index target genomes with Blasr:

```
blasr/alignment/bin/sawriter target.fasta
```


If you're on an SGE cluster you can use the following:

```
 snakemake  -c "qsub {params.sge_opts}" -s Snakefile -w 50  -p -k -j 20
```
 
 Otherwise:

 On a single node:

```
snakemake -s Snakefile -w 50  -p -k -j 20
```
 
 Or modidfy the cluster settings in the config.json for your cluster.
