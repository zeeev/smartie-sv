# smartie-sv

Have contigs and a reference genome?  Smartie-sv will align target contigs against a reference genome and call strucutral variants. 

### Dependancies

1. Ananconda
2. Snakemake
3. Bedtools

### Running

1. Edit the config.jason.
2. If you don't have "modules" make sure your path contains bedtools. 

If you're on an SGE cluster you can use the following:

```
 snakemake  -c "qsub {params.sge_opts}" -s Snakefile -w 50  -p -k -j 20
```
 
 Otherwise:

```
snakemake -s Snakefile -w 50  -p -k -j 20
```
 
