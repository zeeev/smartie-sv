shell.prefix("source config.sh; set -eo pipefail ; ")

configfile: "config.json"

def _get_target_files(wildcards):
    return config["targets"][wildcards.target]

def _get_query_files(wildcards):
        return config["queries"][wildcards.query]

rule dummy:
     input: expand("coverage/{target}-{query}-cov.bed", query=config["queries"].keys(), target=config["targets"].keys()), expand("counts_per_window/{target}-{query}-counts.bed", query=config["queries"].keys(), target=config["targets"].keys())


rule counts:
     input: WINS="data/hg38.windows.bed", SVS="dist_anno_svs/{target}-{query}-svs.bed"
     output: "counts_per_window/{target}-{query}-counts.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell:  """

     bedtools intersect -a {input.WINS} -b {input.SVS} -wao | bedtools groupby -c 1 -o count > $TMPDIR/counts.bed
     bedtools coverage -a $TMPDIR/counts.bed -b {input.SVS} > {output}

     """

rule badRegions:
     input: FA=_get_target_files
     output: "data/{target}_low_confidence_regions.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell: """

     cat {input.FA}.fai | perl -lane 'if($F[1] > 1000000){{$e = $F[1] - 1000000 ; print "$F[0]\t0\t1000000\ttelomere\n$F[0]\t$e\t$F[1]\ttelomere"}}' > $TMPDIR/low_confidence_regions.bed

     curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz" | gunzip -c | grep cen | bedtools merge | perl -lane '$s = $F[1] - 2000000; $e = $F[2] + 2000000; print "$F[0]\t$s\t$e\tcentromere"' >> $TMPDIR/low_confidence_regions.bed

     sort -k1,1 -k2,2n $TMPDIR/low_confidence_regions.bed > {output}

"""

rule filt:
     input: BED="raw_svs/{target}-{query}-svs.bed", TARGET=_get_query_files, COV="coverage/{target}-{query}-cov.bed", CEN_TEL="data/{target}_low_confidence_regions.bed"
     output: "dist_anno_svs/{target}-{query}-svs.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell: """

     perl helper_scripts/annotate_size_dist.pl -f {input.TARGET}.fai -s {input.BED} | perl -lane 'print if $F[-3] > 200000 && $F[-1] > 0.05 && $F[-1] < 0.95' > $TMPDIR/filtA.bed

     perl -lane 'print if $F[-1] < 0.5 || $F[3] > 50' {input.COV} | sort -k1,1 -k2,2n > $TMPDIR/low_cov_regions.bed

     bedtools subtract -A -a $TMPDIR/filtA.bed -b $TMPDIR/low_cov_regions.bed >  $TMPDIR/filtB.bed

     bedtools subtract -A -a $TMPDIR/filtB.bed -b {input.CEN_TEL} > {output}


     """

rule coverage:
     input: WIN="data/hg38.windows.bed", BEDS="beds/{target}-{query}-aligned.bed", OVER="contig_overlap/{target}-{query}-overlap.bed"
     output: "coverage/{target}-{query}-cov.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell: """
     bedtools intersect -a {input.WIN} -b {input.BEDS} -wao | bedtools groupby -c 1 -o count > $TMPDIR/counts.bed
     bedtools coverage -a $TMPDIR/counts.bed -b {input.BEDS} > {output}
     """

rule makeWin:
     input : GENOME=config["bed_genome"]
     output: "data/hg38.windows.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell: """

     bedtools makewindows -g {input.GENOME} -w 1000000 -s 250000 > {output}

     """

rule contigDepth:
     input:  "beds/{target}-{query}-aligned.bed"
     output: "contig_overlap/{target}-{query}-overlap.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell:   """

     bedtools merge -d -100 -c 1 -o count -i {input} | sort -k1,1 -k2,2n > {output}

     """

rule bed:
     input: SAM="mappings/{target}-{query}-aligned.sam", STB=config["samTo_bed_bin"],  SVS="raw_svs/{target}-{query}-svs.bed"
     output: "beds/{target}-{query}-aligned.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell:  """

     {input.STB} {input.SAM} | sort -k1,1 -k2,2n > {output}

     """

rule call:
     message: "Calling SVs."
     input: SAM="mappings/{target}-{query}-aligned.sam", PG=config["print_gaps_bin"], TARGET=_get_target_files
     output: "raw_svs/{target}-{query}-svs.bed"
     params:  sge_opts=config["cluster_settings"]["lite"]
     shell:  """
          {input.PG} --minq 20 --condense 20 --outFile {output} {input.TARGET} {input.SAM}
          """

rule runBlasr:
     message: "Aligning query to target"
     input:   BL="helper_scripts/blasr", TARGET=_get_target_files, QUERY=_get_query_files
     output:  "mappings/{target}-{query}-aligned.sam", "unmappings/{target}-{query}-unaligned.sam"
     params:  sge_opts=config["cluster_settings"]["heavy"], THR=config["threads"]
     shell:   """
              {input.BL} -alignContigs -sam -minMapQV 30 -nproc {params.THR} -minPctIdentity 50 -unaligned {output[1]} {input.QUERY} {input.TARGET} -out {output[0]}
     """

rule installBlasr:
     message: "Installing blasr"
     output:  "helper_scripts/blasr"
     shell:   """
            git clone https://github.com/mchaisso/blasr.git
            cd blasr
            git checkout 2956caed6de86518c48af78c9e0cd80aa307e412
            make
            cp alignment/bin/blasr ../helper_scripts/
     """