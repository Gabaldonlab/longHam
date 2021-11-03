#!/usr/bin/env python


"""
  longHAm v1.0 - automated pipeline for a hybrid long read assembly
  Copyright (C) 2021 - Marina Marcet-Houben, Toni Gabaldon
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import subprocess as sp

cluster_path = "/gpfs/projects/bsc40/current/"
canu_path = cluster_path+"mmarcet/anaconda_old/envs/assembly/bin/canu"
trimmomatic_path = cluster_path+"mmarcet/anaconda_old/envs/assembly/bin/trimmomatic"
platanus_path = cluster_path+"mmarcet/nanopore/programs/platanus"
sparse_path = cluster_path+"mmarcet/anaconda_old/envs/assembly/bin/SparseAssembler"
DBG2OLC_path = cluster_path+"mmarcet/anaconda_old/envs/assembly/bin/DBG2OLC"
masurca_path = cluster_path+"mmarcet/anaconda/envs/masurca/bin/masurca"
bwa_path = cluster_path+"mmarcet/anaconda_old/envs/assembly/bin/bwa"
samtools_path = cluster_path+"mmarcet/anaconda_old/envs/assembly/bin/samtools"
pilon_path = cluster_path+"mmarcet/nanopore/programs/pilon-1.22.jar"
adapters_file = cluster_path+"mmarcet/nanopore/GABALDON02/scripts/adapters_trimmonmatic/TruSeq3-PE.fa"
masurca_master_config=cluster_path+"mmarcet/nanopore/GABALDON02/scripts/master_config_file.masurca.txt"
ragout_path = cluster_path+"mmarcet/nanopore/programs/ragout/ragout-2.0-linux-x86_64/ragout.py"
wtdbg2_path = cluster_path+"mmarcet/nanopore/programs/wtdbg2/"

########################################################################
# Standard modules
########################################################################

#Checks if a folder exists and if not creates it		
def create_folder(name):
    if not os.path.exists(name):
        cmd = "mkdir "+name
        try:
            run_command(cmd,False)
        except:
            print("Unable to create directory ",name)

#This function will execute a command in bash
def run_command(cmd,ommit):
    if ommit:
        try: process = sp.Popen(cmd,shell=True)
        except: pass
        process.communicate("Y\n")
        if process.wait() != 0: print("Error ocurred, but you chose to ommit it")
    else:
        try: process = sp.Popen(cmd,shell=True)
        except OSError: sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0: sys.exit("ERROR: Execution cmd failed")


########################################################################
# Illumina modules
########################################################################

def trim_illumina_reads(reads1,reads2):
    #cmd = "java -jar "+trimmomatic_path+" PE "+reads1+" "+reads2+" reads.paired.1.fastq reads.unpaired.1.fastq reads.paired.2.fastq reads.unpaired.2.fastq ILLUMINACLIP:"+adapters_file+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    cmd = trimmomatic_path+" PE "+reads1+" "+reads2+" reads.paired.1.fastq reads.unpaired.1.fastq reads.paired.2.fastq reads.unpaired.2.fastq ILLUMINACLIP:"+adapters_file+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    print(cmd)
    run_command(cmd,False)
    
def platanus_assembly():
    create_folder("platanus")
    os.chdir("platanus")
    cmd = platanus_path+" assemble -f ../reads.paired* -t "+str(args.threads)
    print(cmd)
    run_command(cmd,False)
    cmd = "cp out_contig.fa ../../assemblies/platanus_assembly.fasta"
    run_command(cmd,False)
    os.chdir("../")

def sparse_assembly():
    create_folder("sparse")
    os.chdir("sparse")
    cmd = sparse_path+" GS "+args.genome_size+" LD 0 k 51 g 15 NodeCovTh 1 EdgeCovTh 0 f ../reads.paired.1.fastq f ../reads.paired.2.fastq"
    print(cmd)
    run_command(cmd,False)
    cmd = "cp Contigs.txt ../../assemblies/sparse_assembly.fasta"
    run_command(cmd,False)
    os.chdir("../")

def canu_whole(nanoporeFile,genomeSize,tag,outAssem):
    if tag == "nanopore":
        cmd = canu_path+" -d '.' -p canu -useGrid=False genomeSize="+genomeSize+" maxMemory=32G maxThreads=25 -nanopore-raw "+nanoporeFile
    elif tag == "pacbio":
        cmd = canu_path+" -correct -d '.' -p canu -useGrid=False genomeSize="+genomeSize+" maxMemory=32G maxThreads=25 -pacbio-raw "+nanoporeFile
    print(cmd)
    run_command(cmd,False)
    cmd = "gunzip canu.correctedReads.fasta.gz"
    run_command(cmd,False)
    cmd = "cp canu.contigs.fasta "+outAssem
    run_command(cmd,False)

def assemble_DBG2OLC(nanoporeFile,illumina_assembly_file,outputfolder,outputName):
    os.chdir(outputfolder)
    cmd = DBG2OLC_path+" k 17 AdaptiveTh 0.01 KmerCovTh 10 MinOverlap 100 RemoveChimera 2 ContigTh 2 Contigs "+illumina_assembly_file+" f "+nanoporeFile
    print(cmd)
    run_command(cmd,False)
    cmd = "cp backbone_raw.fasta ../../assemblies/DBG2OLC_"+outputName+".fasta"
    print(cmd)
    run_command(cmd,False)

def adapt_config_file_masurca(reads_ill1,reads_ill2,readsLong,threads,read_type,coverage):
    cmd = "cp "+masurca_master_config+" config.txt"
    run_command(cmd,False)
    outfile = open("config2.txt","w")
    for line in open("config.txt"):
        line = line.strip()
        if "ILLUMINA_READS_FILE" in line:
            tag = reads_ill1+" "+reads_ill2
            line = line.replace("ILLUMINA_READS_FILE",tag)
        if "LONGREADS" in line:
            print("ENTRO",tag)
            if read_type == "nanopore":
                line = "NANOPORE="+readsLong
            elif read_type == "pacbio":
                line = "PACBIO="+readsLong
        if "NUM_THREADS" in line:
            line = "NUM_THREADS = "+str(threads)
        if "SELECTED_COVERAGE" in line:
            if not coverage:
                line = "LHE_COVERAGE = 30"
            else:
                line = "LHE_COVERAGE = "+str(coverage)
        print(line,file=outfile)
    outfile.close()
    cmd = "mv config2.txt config.txt"
    run_command(cmd,False)

def assemble_masurca():
    cmd = masurca_path+" config.txt"
    run_command(cmd,False)
    cmd = "./assemble.sh"
    run_command(cmd,False)
    cmd = "cp CA.mr.*/final.genome.scf.fasta ../assemblies/masurca_assembly.fasta"
    run_command(cmd,False)

def align_illumina_reads(assembly,alignment,threads,illReads1,illReads2):
    cmd = bwa_path+" index "+assembly
    run_command(cmd,False)
    if not os.path.exists(alignment+".sam"):
        cmd = bwa_path+" mem -t "+threads+" "+assembly+" "+illReads1+" "+illReads2+" >"+alignment+".sam"
        print(cmd)
        run_command(cmd,False)
    if not os.path.exists(alignment+".sorted.bam"):
        cmd = samtools_path+" view -@ "+threads+" -Sb "+alignment+".sam | "+samtools_path+" sort -@ "+threads+" -o "+alignment+".sorted.bam -"
        print(cmd)
        run_command(cmd,False)
        cmd = samtools_path+" index "+alignment+".sorted.bam"
        print(cmd)
        run_command(cmd,False)
        cmd = "rm "+alignment+".sam"
        run_command(cmd,True)

def run_pilon(alignmentFile,genome,folderName,threads):
    cmd = "java -Xmx96G -jar "+pilon_path+" --genome "+genome+" --bam "+alignmentFile+".sorted.bam --outdir "+folderName+" --vcf --changes --threads " +threads
    run_command(cmd,False)

def rename_contigs_pilon(folderName):
    genome = folderName+"/pilon_corrected.fasta"
    outfile = open(genome,"w")
    for line in open(folderName+"/pilon.fasta"):
        line = line.strip()
        if ">" in line:
            line = line.replace("_pilon","")
        print(line,file=outfile)
    outfile.close()
    return genome

def correct_assembly(genome,outDir,repetitions,outfileName,threads,illReads1,illReads2,tag):
    create_folder(outDir)
    original_genome = genome
    for a in range(0,repetitions):
        n = a + 1
        folderName = outDir+"/run_"+str(n)
        create_folder(folderName)
        alignmentFile = folderName+"/alignment"
        print("ALIGN ILLUMINA READS")
        align_illumina_reads(genome,alignmentFile,threads,illReads1,illReads2)
        print("CORRECT WITH PILON")
        run_pilon(alignmentFile,genome,folderName,threads)
        genome = rename_contigs_pilon(folderName)
    rename_contigs(outfileName,folderName+"/pilon_corrected.fasta",tag)
    cmd = "rm "+original_genome+".*"
    run_command(cmd,False)
    
def rename_contigs(outfileName,infileName,tag):
    outfile = open(outfileName,"w")
    seqs = load_sequences(infileName,"")
    contigs = list(seqs.keys())
    contigs = sorted(contigs,key=lambda x:len(seqs[x]),reverse=True)
    for num,contig in enumerate(contigs):
        name = tag+"_"+str(num+1)
        print_sequence(name,seqs[contig],outfile)
    outfile.close()
    

def create_receipt(main_assembly,reference1,reference2,reference_genome,ragout_folder):
    outfile = open(ragout_folder+"/receipt.txt","w")
    if reference_genome:
        print(".references = ref_genome,ref1,ref2",file=outfile)
    else:
        print(".references = ref1,ref2",file=outfile)
    print(".target = target\n",file=outfile)
    
    if reference_genome:
        print("ref_genome.fasta = "+reference_genome,file=outfile)
    print("ref1.fasta = "+reference1,file=outfile)
    print("ref2.fasta = "+reference2,file=outfile)
    print("target.fasta = "+main_assembly,file=outfile)
    outfile.close()
    
def run_ragout(ragout_folder):
    cmd = "python2 "+ragout_path+" -o "+ragout_folder+"/scaffolds -t "+str(args.threads)+" "+ragout_folder+"/receipt.txt"
    run_command(cmd,False)

def calculate_coverage(nanoporeFile,genome_size):
    size = 0
    seqs = load_sequences(nanoporeFile,"")
    for code in seqs:
        size += len(seqs[code])
    coverage = float(size) / genome_size
    return seqs,coverage

def get_nanopore_subset(seqs,reads_names,genome_size,coverage):
    seqs2 = {}
    counting_size = 0
    counting_coverage = 0
    for code in reads_names:
        seqs2[code] = seqs[code]
        counting_size += len(seqs2[code])
        counting_coverage = int(float(counting_size) / genome_size)
        if counting_coverage >= coverage:
            break
    print(len(seqs2),"sequences selected with a coverage of ",counting_coverage)
    return seqs2

def load_sequences(contigFile,delimiter):
    seqs = {}
    name = ""
    s = []
    for line in open(contigFile):
        line = line.strip()
        if ">" in line:
            if name != "":
                seqs[name] = "".join(s)
            if delimiter == "":
                name = line.replace(">","")
            else:
                name = line.replace(">","").split(delimiter)[0]
            s = []
        else:
            s.append(line.upper())
    if len(s) != 0:
        seqs[name] = "".join(s)
    return seqs

def print_sequence(code,sequence,outfile):
    print(">"+code,file=outfile)
    i = 0
    if sequence[-1] == "*":
        sequence = sequence[:-1]
    while i < len(sequence):
        print(sequence[i:i+60],file=outfile)
        i += 60

def assemble_WTDBG2(corrected_nanopore,outputFolder,assemblyWTDBG2):
    os.chdir(outputFolder)
    cmd = wtdbg2_path+"/wtdbg2 -t "+str(args.threads)+" -i "+corrected_nanopore+" -fo prefix"
    run_command(cmd,False)
    cmd = wtdbg2_path+"/wtpoa-cns -t "+str(args.threads)+" -i prefix.ctg.lay -fo assembly.fa"
    run_command(cmd,False)
    cmd = "mv assembly.fa "+assemblyWTDBG2
    run_command(cmd,False)

def get_genome_size(genome_size):
    if "M" in genome_size or "m" in genome_size:
        genome_size = float(genome_size.replace("M","").replace("m",""))*1000000.0
    elif "G" in genome_size or "g" in genome_size:
        genome_size = float(genome_size.replace("G","").replace("g",""))*1000000000.0
    else:
        genome_size = float(genome_size)
    return genome_size

def calculate_genome_stats(contigsFile, toprint=True):
    seqs = load_sequences(contigsFile," ")
    lengths = []
    for code in seqs:
        lengths.append(len(seqs[code]))
    #GENOME SIZE
    size = 0
    for l in lengths:
        size += l
    #Num contigs
    NContigs = len(seqs)
    NContigs1kb = 0
    NContigsValid = 0
    for code in seqs:
        if len(seqs[code]) > 100000:
            NContigs1kb += 1
        if len(seqs[code]) > 1000:
            NContigsValid += 1
    #N50
    median = size/2
    lengths.sort()
    lengths.reverse()
    i = 0
    N50 = False
    L50 = False
    for num,c in enumerate(lengths):
        i += c
        if i > median and not N50:
            N50 = c
            L50 = num + 1
    #GC content
    GC = 0
    total = 0
    for code in seqs:
        seq = seqs[code].upper()
        GC += seq.count("C") + seq.count("G")
        total += seq.count("C") + seq.count("G")+seq.count("A") + seq.count("T")
    percGC = float(GC)/float(total)

    #Number of N
    n = 0
    for code in seqs:
        n += seqs[code].count("N")
    if toprint:
        print("GENOME SIZE: ",size)
        print("CONTIGS: ",NContigs)
        print("CONTIGS > 1000: ",NContigsValid)
        print("CONTIGS > 100kb: ",NContigs1kb)
        print("N50: ",N50)
        print("L50: ",L50)
        print("GC content: ",percGC)
        print("Number of N: ",n)
    return size,NContigs1kb,NContigs,N50,L50

def build_illumina_assemblies():
    #Create illumina assemblies
    illOutput = args.outputFolder+"/illumina"
    create_folder(illOutput)
    os.chdir(illOutput)
    if not os.path.exists(illOutput+"/reads.paired.1.fastq"):
        trim_illumina_reads(args.reads_ill1,args.reads_ill2)

    reads_ill1 = illOutput+"/reads.paired.1.fastq"
    reads_ill2 = illOutput+"/reads.paired.2.fastq"
    
    platanusFile = assemblyOutput+"/platanus_assembly.fasta"
    if not os.path.exists(platanusFile):
        platanus_assembly()
    sparseFile = assemblyOutput+"/sparse_assembly.fasta"
    if not os.path.exists(sparseFile):
        sparse_assembly()
    return reads_ill1,reads_ill2,platanusFile,sparseFile

def build_canu_assembly(primary_assemblies,genome_size):
    #Create CANU assembly    
    canu_output = args.outputFolder+"/canu"
    create_folder(canu_output)
    os.chdir(canu_output)
    assemblyCANU = assemblyOutput+"/canu_assembly.fasta"
    print(assemblyCANU)
    if not os.path.exists(assemblyCANU):
        canu_whole(args.reads_nanopore,str(genome_size),args.reads_type,assemblyCANU)

    primary_assemblies["canu"] = assemblyCANU

    corrected_nanopore = canu_output+"/canu.correctedReads.fasta"
    
    return primary_assemblies,corrected_nanopore,canu_output

def build_wtdbg2_assembly(corrected_nanopore,primary_assemblies):

    #Build long-reads only assembly
    wtdbg2_output = args.outputFolder+"/wtdbg2"
    create_folder(wtdbg2_output)
    try:
        assemblyWTDBG2 = assemblyOutput+"/wtdbg2_assembly.fasta"
        if not os.path.exists(assemblyWTDBG2):
            assemble_WTDBG2(corrected_nanopore,wtdbg2_output,assemblyWTDBG2)
        primary_assemblies["wtdbg2"] = assemblyWTDBG2
    except:
        assemblyWTDBG2 = None
    return primary_assemblies

def subsample_long_reads(long_reads_corrected,genome_size,canu_output):
    #Cut corrected reads and only take 30X coverage
    seqs,coverage = calculate_coverage(long_reads_corrected,genome_size)

    if int(coverage) < args.nanopore_coverage:
        pass
    else:
        reads_names = seqs.keys()
        reads_names = sorted(reads_names,key=lambda x: len(seqs[x]))
        reads_names.reverse()
        seqs2 = get_nanopore_subset(seqs,reads_names,genome_size,args.nanopore_coverage)
        corrected_nanopore = canu_output+"/canu.correctedReads."+str(args.nanopore_coverage)+"X.fasta"
        outfile = open(corrected_nanopore,"w")
        for code in seqs2:
            print_sequence(code,seqs2[code],outfile)
        outfile.close()
    return corrected_nanopore


def build_scaffolding_assemblies(corrected_nanopore,platanusFile,sparseFile,primary_assemblies):
    #Create nanopore + illumina assemblies
    DBG2OLC_output = args.outputFolder+"/DBG2OLC"
    create_folder(DBG2OLC_output)
    DBG2OLC_output_platanus = DBG2OLC_output+"/platanus/"
    outputName = "platanus"
    create_folder(DBG2OLC_output_platanus)
    try:
        assemblyPL = assemblyOutput+"/DBG2OLC_platanus.fasta"
        if not os.path.exists(assemblyPL):
            assemble_DBG2OLC(corrected_nanopore,platanusFile,DBG2OLC_output_platanus,outputName)
        primary_assemblies["DBG2OLC_platanus"] = assemblyPL
    except:
        assemblyPL = None
        print("Platanus + DBG2OLC assembly failed")
        
    DBG2OLC_output_sparse = DBG2OLC_output+"/sparse/"
    outputName = "sparse"
    create_folder(DBG2OLC_output_sparse)
    try:
        assemblySP = assemblyOutput+"/DBG2OLC_sparse.fasta"
        if not os.path.exists(assemblySP):
            assemble_DBG2OLC(corrected_nanopore,sparseFile,DBG2OLC_output_sparse,outputName)
        primary_assemblies["DBG2OLC_sparse"] = assemblySP
    except:
        assemblySP = None
        print("Sparse + DBG2OLC assembly failed")
    return primary_assemblies

def build_masurca_assembly(primary_assemblies):
    masurca_output = args.outputFolder+"/masurca/"
    create_folder(masurca_output)
    os.chdir(masurca_output)
    assemblyM = assemblyOutput+"/masurca_assembly.fasta"
    if not os.path.exists(assemblyM):
        adapt_config_file_masurca(args.reads_ill1,args.reads_ill2,args.reads_nanopore,args.threads,args.reads_type,args.nanopore_coverage)
        assemble_masurca()
    primary_assemblies["masurca"] = assemblyM
    return primary_assemblies

def correct_primary_assemblies(primary_assemblies,reads_ill1,reads_ill2):
    #Correct genome assemblies
    os.chdir(args.outputFolder)
    correction_output = args.outputFolder+"/correction"
    create_folder(correction_output)

    corrected_assemblies = {}
    for assem in primary_assemblies:
        correction_output1 = correction_output+"/"+assem+"/"
        corrected_assembly = args.outputFolder+"/assemblies/"+assem+"_assembly.corrected.fasta"
        corrected_assemblies[assem] = corrected_assembly
        if not os.path.exists(corrected_assembly):
            correct_assembly(primary_assemblies[assem],correction_output1,1,corrected_assembly,str(args.threads),reads_ill1,reads_ill2,assem)
    return corrected_assemblies

def assess_assemblies(corrected_assemblies):
    #Assess assemblies
    info = {}
    for assem in corrected_assemblies:
        size,NContigs1kb,Ncontigs,N50,L50 = calculate_genome_stats(corrected_assemblies[assem],toprint=False)
        diff = abs(genome_size - size)
        info[assem] = [diff,size,NContigs1kb,Ncontigs,N50,L50]

    assems = list(info.keys())
    assems = sorted(assems,key=lambda x:info[x][0])
    valid_assemblies = []
    for assem in assems:
        p = info[assem][1] / genome_size * 100.0
        if p < 70.0 or p > 130.0:
            print("Assembly ",assem,"has more than 30% variation compared to the reference genome size provided. This assembly will not be used")
        else:
            valid_assemblies.append(assem)

    return valid_assemblies,info

def ragout_analysis(valid_assemblies,info,corrected_assemblies):
    correction_output = args.outputFolder+"/correction"
    if len(valid_assemblies) < 3:
        exit("There were not enough valid assemblies to use Ragout. The pipeline will finish here")
    else:
        valid_assemblies = sorted(valid_assemblies,key=lambda x:info[x][2])
        best1 = valid_assemblies[0]
        refAssem1 = corrected_assemblies[best1]
        best2 = valid_assemblies[1]
        refAssem2 = corrected_assemblies[best2]
        others = valid_assemblies[2:]
        ragout_folder = args.outputFolder+"/ragout/"
        create_folder(ragout_folder)
        os.chdir(ragout_folder)
        ragout_assemblies = []
        for assem in others:
            try:
                ragout1_folder = ragout_folder + "/target_"+assem
                create_folder(ragout1_folder)
                if not os.path.exists(ragout1_folder+"/scaffolds/target_scaffolds.fasta"):
                    main_assembly = corrected_assemblies[assem]
                    create_receipt(main_assembly,refAssem1,refAssem2,args.reference_genome,ragout1_folder)
                    run_ragout(ragout1_folder)
                ragout_assembly = assemblyOutput+"/ragout_"+str(assem)+".fasta"
                rename_contigs(ragout_assembly,ragout1_folder+"/scaffolds/target_scaffolds.fasta","Ragout_"+assem)
                correction_output1 = correction_output+"/ragout_"+assem+"/"
                corrected_assembly = assemblyOutput+"/ragout_"+str(assem)+".corrected.fasta"
                correct_assembly(ragout_assembly,correction_output1,3,corrected_assembly,str(args.threads),reads_ill1,reads_ill2,"Ragout_"+assem)
            except:
                print("Ragout with "+assem+" failed to produce an assembly")

parser = argparse.ArgumentParser(description="Needed to calculate different statistics from the clusters")
parser.add_argument("-s1",dest="reads_ill1",action="store",required=True,help="Illumina reads 1")
parser.add_argument("-s2",dest="reads_ill2",action="store",required=True,help="Illumina reads 2")
parser.add_argument("-n",dest="reads_nanopore",action="store",required=True,help="Nanopore reads. Adapter trimming needs to be performed beforhand")
parser.add_argument("-o",dest="outputFolder",action="store",required=True,help="Folder where the results will be printed")
parser.add_argument("-t",dest="threads",action="store",type=int,default=8,help="Number of threads")
parser.add_argument("-g",dest="genome_size",action="store",required=True,help="Estimated genome size")
parser.add_argument("-r",dest="reference_genome",action="store",default=None,help="Reference genome path if it's used for scaffolding")
parser.add_argument("-c",dest="nanopore_coverage",action="store",type=int,default=30,help="Limits the nanopore coverage, set to 0 if you don't want to use this filter")
parser.add_argument("--reads_type",dest="reads_type",action="store",default="nanopore",help="Which kind of long reads we have")
args = parser.parse_args()

#Create output folder
create_folder(args.outputFolder)

#Create assemblies folder
assemblyOutput = args.outputFolder+"/assemblies/"
create_folder(assemblyOutput)

genome_size = get_genome_size(args.genome_size)

primary_assemblies = {}

reads_ill1,reads_ill2,platanusFile,sparseFile = build_illumina_assemblies()

primary_assemblies,corrected_nanopore,canu_output = build_canu_assembly(primary_assemblies,genome_size)

if not os.path.exists(corrected_nanopore):
    exit("For some reason corrected nanopore reads do not exists, please, check out that CANU ran correctly and that the corrected reads are uncompressed")

primary_assemblies = build_wtdbg2_assembly(corrected_nanopore,primary_assemblies)

corrected_nanopore = subsample_long_reads(corrected_nanopore,genome_size,canu_output)

primary_assemblies = build_scaffolding_assemblies(corrected_nanopore,platanusFile,sparseFile,primary_assemblies)

primary_assemblies = build_masurca_assembly(primary_assemblies)

corrected_assemblies = correct_primary_assemblies(primary_assemblies,reads_ill1,reads_ill2)

valid_assemblies,info = assess_assemblies(corrected_assemblies)
 
ragout_analysis(valid_assemblies,info,corrected_assemblies)


