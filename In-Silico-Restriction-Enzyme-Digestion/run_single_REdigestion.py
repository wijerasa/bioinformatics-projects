#Saranga Wijeratne, MCIC, OARDC, The Ohio State University
#Copyright (c) 2015 Saranga Wijeratne.
#Permission is hereby granted, free of charge, to any person obtaining a copy of this
#software and associated documentation files (the "Software"), to deal in the Software
#without restriction, including without limitation the rights to use, copy, modify,
#merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to the following
#conditions:
#
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
import itertools
import difflib as df
import os,sys,subprocess
import shutil,glob,re
import ConfigParser,glob
import pprint
import argparse


def making_directory(path, dirname):
    path_to_dir=os.path.join(path,dirname)
    try:
        print "Removing exsiting {0} Directory".format(dirname)
        shutil.rmtree(path_to_dir)
        print "Creating New {0} Directory".format(dirname)
        os.mkdir(path_to_dir)
    except:
        print "Creating new {0} Directory".format(dirname)
        os.makedirs(path_to_dir)


def enzyme_database():
    enzyme_DB={}
    read_enzymes_file=open("Enzymes_Tassel_mod.cfg", "r").read().splitlines()
    headers=read_enzymes_file[0]
    for lines in read_enzymes_file[1:]:
        enzyme=lines.split()[0]
        enzyme_site=lines.split()[2]
        enzyme_DB[enzyme]=enzyme_site

    return enzyme_DB





def singleDenzyme(genome, ouputfolder, enzyme, enzyme_site):
    proc = subprocess.Popen(['perl', 'singleDigestion.pl', genome, enzyme, enzyme_site,ouputfolder],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=False)
    (out, err) = proc.communicate()
    if not err:
        print "Single Digestion Restriction Enzyme analysis complete for enzyme {0}.....".format(enzyme)
    else:
        print " Erro {0}".format(err)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This scritp will calculate number of restriction fragment of digestons (Single Enzyme Digestion) with and wthout size selection process Author Saranga \n Contact wijeratne.3 AT osue.edu')
    parser.add_argument('-i','--input', help='Path to genome file in fasta format', required=True, metavar="FILE")
    parser.add_argument('-o', "--output", help="Path to Output folder", required=True, metavar="FILE")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0 ')
    args = parser.parse_args()
    if args.input:
        genome_fasta = args.input
    if args.output:
        analysis_folder = args.output
    single_dEnzyme_db=enzyme_database()
    fragment_information=[]
    for enzyme,enzyme_site in single_dEnzyme_db.items():
        making_directory(analysis_folder, enzyme)
        output_folder=os.path.join(analysis_folder,enzyme)
        singleDenzyme(genome_fasta, output_folder, enzyme, enzyme_site)
        summary_file=os.path.join(output_folder,"digestion_summary_{0}_by_{1}".format(genome_fasta,enzyme))
        open_summary_files=open(summary_file,"r").read()
        open_summary_files_data=filter(None,open_summary_files.split('\n',2)[2].split('\n'))
        fragment_data_information_header_index=open_summary_files_data.index('all_fragments_number\tall_fragments_coverage\tfragments_in_range_number\tfragments_in_range_coverage')
        fragment_data_information=open_summary_files_data[fragment_data_information_header_index+1]
        fragment_histro_data_header_index=open_summary_files_data.index('Length_range\tfragments_number\tfragments_ratio')
        fragment_histro_data_header=open_summary_files_data[fragment_histro_data_header_index]
        fragment_histro_data=open_summary_files_data[fragment_histro_data_header_index+1:]
        fragment_information.append([enzyme,fragment_data_information])
        



    with open("Summary_Single_Digestion.txt", "w") as outfile:
        outfile.write('Enzyme\tall_fragments_number\tall_fragments_coverage\tfragments_in_(300-600)\tfragments_in_(300-600)_coverage\n')
        for enzyme_name,data in fragment_information:
            outfile.write("{0}\t{1}\n".format(enzyme_name,data))
    #pprint.pprint(single_dEnzyme_db)