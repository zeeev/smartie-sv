#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <map>
#include "faidx.h"
#include "split.hpp"

std::map<char, std::string > keys;

struct event{
    char         type;
    uint32_t      len;
    uint32_t t_offset;
    uint32_t q_offset;
    int          tlen;
    int          qlen;
};

bool parseCigar(std::string & data,
                std::vector<event *> & cigar,
                long int * match,
                long int * bases,
                long int * tConsumed,
                long int * qConsumed ){
    uint64_t start = 0,   nchars  = 0;
    uint64_t tOffset = 0, qOffset = 0;
    for(uint64_t i = 0; i < data.size(); i++){
        if(data[i] > 60){
            event * tmp = new event;
            tmp->t_offset = tOffset;
            tmp->q_offset = qOffset;
            tmp->len  = atol(data.substr(start, nchars).c_str());
            tmp->type = data[i];
            switch(data[i]){
            case 'M':
                {
                    *qConsumed += tmp->len;
                    *tConsumed += tmp->len;
                    *bases     += tmp->len;
                    *match     += tmp->len;
                    tOffset    += tmp->len;
                    qOffset    += tmp->len;
                    tmp->tlen  = tmp->len;
                    tmp->qlen  = tmp->len;
                    break;
                }
            case 'X':
                {
                    *qConsumed += tmp->len;
                    *tConsumed += tmp->len;
                    *bases     += tmp->len;
                    tOffset    += tmp->len;
                    qOffset    += tmp->len;
                    tmp->tlen  = tmp->len;
                    tmp->qlen  = tmp->len;
                    break;
                }
            case '=':
                {
                    *qConsumed += tmp->len;
                    *tConsumed += tmp->len;
                    *bases     += tmp->len;
                    *match     += tmp->len;
                    tOffset    += tmp->len;
                    qOffset    += tmp->len;
                    tmp->tlen  = tmp->len;
                    tmp->qlen  = tmp->len;
                    break;
                }
            case 'I':
                {
                    *qConsumed += tmp->len;
                    *bases     += tmp->len;
                    qOffset    += tmp->len;
                    tmp->tlen  = 0       ;
                    tmp->qlen  = tmp->len;
                    break;
                }
            case 'D':
                {
                    *tConsumed += tmp->len;
                    *bases     += tmp->len;
                    tOffset    += tmp->len;
                    tmp->tlen  = tmp->len;
                    tmp->qlen  = 0;
                    break;
                }
            case 'S':
                {
                    *qConsumed += tmp->len;
                    qOffset    += tmp->len;
                    tmp->tlen  = 0       ;
                    tmp->qlen  = tmp->len;
                    break;
                }
            case 'H':
                {
                    break;
                }
            case 'N':
                {
                    *tConsumed += tmp->len;
                    *bases     += tmp->len;
                    tOffset    += tmp->len;
                    tmp->tlen  = tmp->len;
                    tmp->qlen  = 0       ;
                    break;
                }
            case 'P':
                {
                    break;
                }
            default:
                {
                    std::cerr << "FATAL: bad cigar" << data[i] << std::endl;
                    exit(1);
                }
            }

            start += nchars +1 ;
            nchars = 0;
            cigar.push_back(tmp);
        }
        else{
            nchars++;
        }
    }
    return true;
}

bool freeCigar(std::vector<event *> & data){
    for(uint64_t i = 0; i < data.size(); i++){
        delete data[i];
    }
}

bool printGaps(std::vector<event *> & cigars,
               uint64_t  tStart,
               std::string & qSeq,
               std::string & qname,
               std::string & rname,
               int alignmentFlag,
               faidx_t  * FA,
               long int * match,
               long int * bases,
               long int * tConsumed,
               long int * qConsumed,
               std::ofstream * bed,
               std::ofstream * snps,
               std::ofstream * svs,
               std::ofstream * indels){

    uint64_t queryStart  = 0;
    uint64_t queryOffset = 0;
    long int queryLen    = *qConsumed;

    if(cigars.front()->type == 'S') *qConsumed -= cigars.front()->len;
    if(cigars.back()->type == 'S') *qConsumed -= cigars.back()->len;

    if(cigars.front()->type == 'H') queryOffset += cigars.front()->len;
    if(cigars.front()->type == 'H' || cigars.front()->type == 'S') queryStart += cigars.front()->len;

    char strand = '+';
    if(alignmentFlag & 16) strand = '-';
    if(cigars.front()->type == 'H') queryLen += cigars.front()->len;
    if(cigars.back()->type  == 'H') queryLen += cigars.back()->len;

    *bed << rname  << "\t"
         << tStart - 1 << "\t"
         << tStart + *tConsumed - 1 << "\t"
         << qname << "\t"
         << queryStart - 1 << "\t"
         << queryStart + *qConsumed - 1<< "\t"
         << strand <<  "\t"
         << queryLen << "\t"
         << *match << "\t"
         << *bases << "\t"
         << double(*match)/double(*bases)
         << std::endl;

    for(uint32_t i = 0; i < cigars.size(); i++){
        if(cigars[i]->type == '=' || cigars[i]->type == 'M' || cigars[i]->type == 'H' || cigars[i]->type == 'S') continue;

        std::string dna =".";
        if(cigars[i]->type == 'D' || cigars[i]->type == 'X'){

            int len = 0;

            char * tmp = faidx_fetch_seq(FA,
                                         rname.c_str(),
                                         (cigars[i]->t_offset + tStart - 1),
                                         (cigars[i]->t_offset + tStart + cigars[i]->len - 2),
                                         &len);
            dna = (std::string)tmp;

        }
        if(cigars[i]->type == 'I'){
            dna = qSeq.substr(cigars[i]->q_offset, cigars[i]->len);
        }

        std::stringstream ss;

        if(keys.find(cigars[i]->type) == keys.end()){
            std::cerr << "FATAL: unknown cigar op for gap printing" << cigars[i]->type << std::endl;
            exit(1);
        }
        ss << rname << "\t"
           << cigars[i]->t_offset + tStart - 1 << "\t"
           << cigars[i]->t_offset + tStart + cigars[i]->tlen -1 <<  "\t"
           << keys[cigars[i]->type] << "\t"
           << cigars[i]->len  << "\t"
           << strand << "\t"
           << qname << "\t"
           << cigars[i]->q_offset - 1 << "\t"
           << cigars[i]->q_offset + cigars[i]->qlen << "\t"
           << queryLen << "\t"
           << *match << "\t" << *bases << "\t" << double(*match)/double(*bases) << "\t"
           << dna;

        if(cigars[i]->type == 'X'){
            dna = qSeq.substr(cigars[i]->q_offset, cigars[i]->len);
            ss << "\t" << dna << std::endl;
        }
        else{
            ss << std::endl;
        }





        if(cigars[i]->type == 'X') *snps << ss.str();
        if((cigars[i]->type != 'X') & (cigars[i]->len > 49) ) *svs    << ss.str();
        if((cigars[i]->type != 'X') & (cigars[i]->len < 50) ) *indels << ss.str();
    }
}

int main(int argc, char *argv[]){

    if(argc != 3){
        std::cerr << "Usage: cat <sam> | printgaps <fasta> <prefix> " << std::endl << std::endl;
        std::cerr << "Required: <fasta>  - The target fasta the query was aligned against, matching input." <<std::endl;
        std::cerr << "          <prefix> - A prefix for output files.                                     " <<std::endl;
        std::cerr << std::endl;
        std::cerr << "Details:  This program expects a sam file from stdin." << std::endl;
        std::cerr << "          The input can be piped from samtools view if you have a bam." << std::endl << std::endl;
        std::cerr << "Output:   The program outputs four files, each with a header line:" << std::endl << std::endl;
        std::cerr << "          1. A bed file describing the alignment blocks." << std::endl;
        std::cerr << "          2. A SNV file describing mismatches from the CIGAR." << std::endl;
        std::cerr << "          3. An indel file for events greater than 2bp and less than 50bp." << std::endl;
        std::cerr << "          4. A Structural variant file. " << std::endl;
        exit(1);
    }


    keys['X'] = "snv";
    keys['D'] = "deletion";
    keys['I'] = "insertion";

    std::ofstream bed      ;
    std::ofstream snp      ;
    std::ofstream indel    ;
    std::ofstream svs      ;

    std::string prefix = argv[2];

    bed.open( (prefix   + ".aln.coords.bed").c_str());
    snp.open( (prefix   + ".snp.bed").c_str());
    svs.open( (prefix   + ".svs.bed").c_str());
    indel.open( (prefix + ".indel.bed").c_str());

    bed << "#"
        << "target_name\t"
        << "target_start_0\t"
        << "target_end_0\t"
        << "query_name\t"
        << "query_start_0\t"
        << "query_end_0\t"
        << "query_strand\t"
        << "query_length\t"
        << "matching_bases_count\t"
        << "consumed_bases_count\t"
        << "percent_identity\n";

    std::stringstream header;
    header << "#"
           << "target_name\t"
           << "target_start_0\t"
           << "target_end_0\t"
           << "sv_type\t"
           << "sv_len\t"
           << "query_strand\t"
           << "query_name\t"
           << "query_start_0\t"
           << "query_end_0\t"
           << "query_length\t"
           << "matching_bases_count\t"
           << "consumed_bases_count\t"
           << "percent_identity\t"
           << "sequence";

    indel << header.str() << std::endl;
    snp   << header.str() << "\talt_sequence" << std::endl;
    svs   << header.str() << std::endl;


    std::string line;

    uint64_t nlines = 0;

    const char * fasta = argv[1];
    std::cerr << "INFO: fasta: " << fasta << std::endl;

    faidx_t * FA = fai_load(argv[1]);

    while(std::getline(std::cin, line)){
        if(line[0] == '@'){
            continue;
        }
        std::vector<std::string> lineDat = split(line, "\t");
        if(atol(lineDat[4].c_str()) < 20 || lineDat[5] == "*" ){
            continue;
        }

        long int match = 0, bases = 0, tConsumed = 0, qConsumed = 0;

        std::vector<event *> cigarData;
        parseCigar(lineDat[5], cigarData, &match, &bases, &tConsumed, &qConsumed);

        nlines++;
        if((nlines % 100) == 0 ) std::cerr << "parsed n: " << nlines << " alignments." << std::endl;

        printGaps(cigarData, atol( lineDat[3].c_str() ),
                  lineDat[9],
                  lineDat[0],
                  lineDat[2],
                  atoi(lineDat[1].c_str()),
                  FA,
                  &match,
                  &bases,
                  &tConsumed,
                  &qConsumed,
                  &bed,
                  &snp,
                  &svs,
                  &indel);

        freeCigar(cigarData);

    }


    return 0;
}
