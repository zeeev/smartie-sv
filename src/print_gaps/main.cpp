#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <map>
#include "faidx.h"
#include "split.hpp"

std::map<char, std::string > keys;

struct event{
    char         type;
    uint32_t      len;
    uint32_t t_offset;
    uint32_t q_offset;
};

bool parseCigar(std::string & data, std::vector<event *> & cigar, long int * match, long int * bases){
    uint64_t start = 0,   nchars  = 0;
    uint64_t tOffset = 0, qOffset = 0;
    for(uint64_t i = 0; i < data.size(); i++){
        if(data[i] > 60){
            event * tmp = new event;
            tmp->len  = atol(data.substr(start, nchars).c_str());
            tmp->type = data[i];
            switch(data[i]){
            case 'M':
                {
                    *bases  += tmp->len;
                    *match  += tmp->len;
                    tOffset += tmp->len;
                    qOffset += tmp->len;
                    break;
                }
            case 'X':
                {
                    *bases  += tmp->len;
                    tOffset += tmp->len;
                    qOffset += tmp->len;
                    break;
                }
            case '=':
                {
                    *bases  += tmp->len;
                    *match  += tmp->len;
                    tOffset += tmp->len;
                    qOffset += tmp->len;
                    break;
                }
            case 'I':
                {
                    *bases  += tmp->len;
                    qOffset += tmp->len;
                    break;
                }
            case 'D':
                {
                    *bases  += tmp->len;
                    tOffset += tmp->len;
                    break;
                }
            case 'S':
                {
                    qOffset += tmp->len;
                    break;
                }
            case 'H':
                {
                    break;
                }
            case 'N':
                {
                    *bases  += tmp->len;
                    tOffset += tmp->len;
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
            tmp->t_offset = tOffset;
            tmp->q_offset = qOffset;

            start += nchars +1 ;// skip the cigar op too
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
               std::string & rname,
               faidx_t * FA,
               long int * match,
               long int * bases ){

    uint64_t queryStart = 0;

    if(cigars.front()->type == 'S'){
        queryStart += cigars.front()->len;
    }

    for(uint32_t i = 0; i < cigars.size(); i++){
        std::string dna =".";
        if(cigars[i]->type == 'D' || cigars[i]->type == 'X'){

            int len = 0;

            char * tmp = faidx_fetch_seq(FA, rname.c_str(),
                                         cigars[i]->t_offset + tStart - 1,
                                        cigars[i]->t_offset + tStart - 1 + cigars[i]->len -1,
                                         &len);
            dna = (std::string)tmp;

        }
        if(cigars[i]->type == 'M' || cigars[i]->type == '=' ){
            continue;
        }
        if(cigars[i]->type == 'I'){
            dna = qSeq.substr(queryStart + cigars[i]->q_offset - 1, cigars[i]->len);
        }


        std::cout << rname << "\t"
                  << cigars[i]->t_offset + tStart << "\t"
                  << cigars[i]->t_offset + tStart + cigars[i]->len <<  "\t"
                  << keys[cigars[i]->type] << "\t"
                  << cigars[i]->len  << "\t"
                  << *match << "\t" << *bases << "\t" << double(*match)/double(*bases) << "\t"
                  << dna
                  << std::endl;
            }
}


int main(int argc, char *argv[]){

    keys['X'] = "snv";
    keys['D'] = "deletion";
    keys['I'] = "insertion";

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

        long int match = 0, bases = 0;

        std::vector<event *> cigarData;
        parseCigar(lineDat[5], cigarData, &match, &bases);

        nlines++;
        if((nlines % 100) == 0 ) std::cerr << "parsed n: " << nlines << " alignments." << std::endl;

        printGaps(cigarData, atol( lineDat[3].c_str() ), lineDat[9], lineDat[2], FA, &match, &bases);
        freeCigar(cigarData);

    }


    return 0;
}
