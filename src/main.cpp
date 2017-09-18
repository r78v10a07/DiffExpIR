/* 
 * File:   main.cpp
 * Author: veraalva
 *
 * Created on July 10, 2017, 11:58 AM
 */
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <ctime>
#include <cstdint>
#include <random>
#include <chrono>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <map>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <set>

#include "api/BamReader.h"

#include "Global.h"
#include "Exceptions.h"
#include "TimeUtils.h"
#include "bstring.h"
#include "TextParser.h"
#include "Sequence.h"
#include "GenomeFactory.h"
#include "ReadFactory.h"
#include "DiffExpIR.h"
#include "Stats.h"

using namespace std;
using namespace ngs;
using namespace sequence;
using namespace genome;
using namespace BamTools;

Global *Global::s_instance = 0;

void print_usage(char *program_name, int exit_code) {
    cerr << "\n********************************************************************************\n";
    cerr << "\nUsage: " << program_name;
    cerr << "\n\n" << program_name << " options:\n\n";
    cerr << "-h    Display this usage information.\n";
    cerr << "-g    GTF file\n";
    cerr << "-o    Output file name\n";
    cerr << "-d    Directory with the TPM output files\n";
    cerr << "-k    Gene key to use from GTF file. Default: gene_id\n";
    cerr << "-t    Transcript key to use from GTF file. Default: transcript_id\n";
    cerr << "-c    Smaller size allowed for an intron created for genes. Default: 16. We recommend to use the reads length\n";
    cerr << "-p    Prefix for grouping samples. (sample_1,sample_2)\n";
    cerr << "-s    Stat method: ttest (default), wilcox\n";
    cerr << "-f    Minimum fold change to filter out (default value: 2.0), wilcox\n";
    cerr << "-v    Minimum P-Value to filter out (default value: 1.0E-6)\n";
    cerr << "-r    Minimum fold change between intron and neighboring exons (default value: -1.0)\n";
    cerr << "-fdr    FRD Correction on the P-Values\n";
    cerr << "\n********************************************************************************\n";
    cerr << "\n                        Roberto Vera Alvarez, PhD\n";
    cerr << "            Emails: veraalva@ncbi.nlm.nih.gov, r78v10a07@gmail.com\n\n";
    cerr << "********************************************************************************\n";
    exit(exit_code);
}

int main(int argc, char *argv[]) {
    TimeUtils uTime;
    string gtfFileName;
    string tpmDirName;
    string geneNameKey = "gene_id";
    string transcriptNameKey = "transcript_id";
    string method = "ttest";
    string output_name;
    double fc_cutoff = 2.0;
    double pvalue_cutoff = 1.0E-6;
    double r_cutoff = -1.0;
    bool useFDR = false;
    int intronCutOff = 16;
    vector<string> samples;
    set<string>features = {"exon"};
    unordered_map<string, string> featuresToCreate = {
        {"exon", "intron"}
    };
    ReadFactory readFactory;
    DiffExpIR diffExpIR;

    if (argc == 1) {
        print_usage(argv[0], 0);
    }

    for (int i = 1; i < argc; i++) {
        string option(argv[i]);
        if (option.compare(0, 1, "-") == 0 && option.compare(1, 1, "-") != 0 && option.size() == 2) {
            if (option.compare(1, 1, "h") == 0) {
                print_usage(argv[0], 0);
            } else if (option.compare(1, 1, "f") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option f require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    fc_cutoff = atof(argument.c_str());
                } else {
                    cerr << "Option f require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 3, "fdr") == 0) {
                useFDR = true;                
            } else if (option.compare(1, 1, "v") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option v require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    pvalue_cutoff = atof(argument.c_str());
                } else {
                    cerr << "Option v require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "r") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    r_cutoff = atof(argument.c_str());
                } else {
                    cerr << "Option r require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "g") == 0) {
                i++;
                if (i < argc) {
                    gtfFileName = argv[i];
                    if (gtfFileName.compare(0, 1, "-") == 0) {
                        cerr << "Option g require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option g require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "o") == 0) {
                i++;
                if (i < argc) {
                    output_name = argv[i];
                    if (output_name.compare(0, 1, "-") == 0) {
                        cerr << "Option o require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option o require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "s") == 0) {
                i++;
                if (i < argc) {
                    method = argv[i];
                    if (method.compare(0, 1, "-") == 0) {
                        cerr << "Option s require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option s require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "d") == 0) {
                i++;
                if (i < argc) {
                    tpmDirName = argv[i];
                    if (tpmDirName.compare(0, 1, "-") == 0) {
                        cerr << "Option d require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option d require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "k") == 0) {
                i++;
                if (i < argc) {
                    geneNameKey = argv[i];
                    if (geneNameKey.compare(0, 1, "-") == 0) {
                        cerr << "Option k require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option k require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "t") == 0) {
                i++;
                if (i < argc) {
                    transcriptNameKey = argv[i];
                    if (transcriptNameKey.compare(0, 1, "-") == 0) {
                        cerr << "Option t require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option t require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "p") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option p require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    BString::split(argument, ",", samples);
                    if (samples.size() != 2) {
                        cerr << "Please, use two samples prefix coma separated" << endl;
                        print_usage(argv[0], -1);
                    }
                } else {
                    cerr << "Option p require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else if (option.compare(1, 1, "c") == 0) {
                i++;
                if (i < argc) {
                    string argument(argv[i]);
                    if (argument.compare(0, 1, "-") == 0) {
                        cerr << "Option c require an argument" << endl;
                        print_usage(argv[0], -1);
                    }
                    intronCutOff = atoi(argv[i]);
                } else {
                    cerr << "Option t require an argument" << endl;
                    print_usage(argv[0], -1);
                }
            } else {
                cerr << "Unsupported option: " << option << endl;
                print_usage(argv[0], -1);
            }
        } else {
            cerr << "Unsupported option: " << option << endl;
            print_usage(argv[0], -1);
        }
    }

    if (gtfFileName.empty()) {
        cerr << "\nGTF is required. See -g option" << endl;
        print_usage(argv[0], -1);
    }

    if (tpmDirName.empty()) {
        cerr << "\nDirectory with the BTPM output files is required. See -d or -b options" << endl;
        print_usage(argv[0], -1);
    }

    if (samples.empty()) {
        cerr << "\nSamples prefixes are required. See -p option" << endl;
        print_usage(argv[0], -1);
    }

    if (method != "ttest" && method != "wilcox") {
        cerr << "\nStat method should be ttest or wilcox. See -s option" << endl;
        print_usage(argv[0], -1);
    }

    if (output_name.empty()) {
        cerr << "\nThe output prefix are required. See -o option" << endl;
        print_usage(argv[0], -1);
    }

    uTime.setTime();
    cerr << "Reading GTF file ... " << endl;
    readFactory.getGenomeFactory().setIntronCutOff(intronCutOff);
    readFactory.getGenomeFactory().processGTFFile(gtfFileName, geneNameKey, transcriptNameKey, features, featuresToCreate);
    cerr << "Done in " << uTime.getElapseTimeSec() << " seconds" << endl;


    uTime.setTime();
    cerr << "Reading TMP output data ... " << endl;
    readFactory.loadTPMCalculatorGenesOutput(tpmDirName);
    cerr << "Done in " << uTime.getElapseTimeSec() << " seconds" << endl;

    uTime.setTime();
    cerr << "Processing TMP data ... " << endl;
    diffExpIR.calculateDiffExpIR(readFactory, samples, method, useFDR);
    cerr << "Done in " << uTime.getElapseTimeSec() << " seconds" << endl;

    uTime.setTime();
    cerr << "Printing output data ... " << endl;
    diffExpIR.printDiffExpIR(output_name, fc_cutoff, pvalue_cutoff, r_cutoff);
    cerr << "Done in " << uTime.getElapseTimeSec() << " seconds" << endl;

    cerr << "Total time: " << uTime.getTotalTimeSec() << " seconds" << endl;
    return 0;
}