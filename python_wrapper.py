# imports python packages
import os
import subprocess as sp
from Bio import SeqIO


# all variable names, in a list of list. [name, assc1, assc1, ftp address, SRA asscession]
var_names = [
        ['HM27', 'GCF_000387825.2', 'ASM38782v2', 
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz', 'SRR1278956'],
        ['HM46', 'GCF_000387845.2', 'ASM38784v2', 
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz', 'SRR1278960'], 
        ['HM65', 'GCF_000387785.2', 'ASM38778v2', 
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz', 'SRR1283106'], 
        ['HM69', 'GCF_000387865.2', 'ASM38786v2', 
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/latest_assembly_versions/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz', 'SRR1278963']]


# to fetch fasta file from ftp address
def fetch_files():
    out_file.write('Fetching fasta files.\n \n')
    
    for name in var_names:
        # make directory for downloading fasta file
        sp.call(['mkdir', name[0]])
        sp.call(['mkdir', name[0] + '/fasta'])
        
        # fetches ftp address with wget
        sp.call(['wget', '-O', name[0] + '/fasta/' + name[0] + '_fasta.fna.gz', name[3]])
        
        # unzip .gz file into fna file
        sp.call(['gunzip', name[0] + "/fasta/" + name[0] + '_fasta.fna.gz'])
        out_file.write(str(name[0]) + '.fna file downloaded.\n')
    
    out_file.write('\n \n \n')


# counts contings in each file
def num_contigs():
    out_file.write('Counting contigs in each assembly.\n \n')
    
    for name in var_names:
        # opens fasta file with SeqIO into a list (records)
        with open(name[0] + '/fasta/' + name[0] + '_fasta.fna', 'rU') as f:
            records = list(SeqIO.parse(f, 'fasta'))
            
        # prints the length of records
        out_file.write('There are {} contigs in the assembly {}.\n'.format(len(records), name[0]))
    
    out_file.write('\n \n \n')
    
    
def assembly_length():
    out_file.write('Getting length of each assembly.\n \n')
    
    for name in var_names:
        length = 0
        # opens fasta file with SeqIO into list (records)
        with open(name[0] + '/fasta/' + name[0] + '_fasta.fna', 'rU') as f:
            records = list(SeqIO.parse(f, 'fasta'))
        
        # loops records and adds length of record to length if > 1000
        for each in records:
            if len(each.seq) > 1000:
                length += len(each.seq)
        out_file.write('There are {} bp in the assembly {}.\n'.format(length, name[0]))
            
    out_file.write('\n \n \n')
    
    
def run_prokka():
    # runns prokka to annotate assembly
    out_file.write('Running Prokka.\n \n')
    
    for name in var_names:
        sp.call(['prokka', '--outdir', name[0] + '/annotation', '--prefix', name[0], '--genus', 'Escherichia', name[0] + '/fasta/' + name[0] + '_fasta.fna'])
        
        out_file.write('Running: prokka --outdir' + name[0] + '/annotation --prefix ' + name[0] + ' --genus Escherichia ' + name[0] + '/fasta/' + name[0] + '_fasta.fna\n \n')
        
    out_file.write('\n \n \n')
    
    
def diff_prokka():
    # calculates the difference from the reference
    for name in var_names:
        CDS = 0
        tRNA = 0
        # opens annotation file
        with open(name[0] + '/annotation/' + name[0] + '.txt', 'r') as f:
            out_file.write('Annotation of ' + name[0] + '.\n')
            # finds CDS and tRNA results
            for text in f:
                if 'CDS' in text:
                    CDS = int(text.split()[1])
                if 'tRNA' in text:
                    tRNA = int(text.split()[1])
                out_file.write(text)
            out_file.write('\n')
        
        # prints based on differences
        if CDS > 4140:
            if tRNA > 89:
                out_file.write('Prokka found {} additional CDS and {} more tRNA than the RefSeq in assembly {}'.format(CDS-4140, tRNA-89, name[0]))
            elif tRNA < 89:
                out_file.write('Prokka found {} additional CDS and {} less tRNA than the RefSeq in assembly {}'.format(CDS-4140, 89-tRNA, name[0]))
            else:
                out_file.write('Prokka found {} additional CDS and the same tRNA than the RefSeq in assembly {}'.format(CDS-4140, name[0]))
        elif CDS < 4140:
            if tRNA > 89:
                out_file.write('Prokka found {} fewer CDS and {} more tRNA than the RefSeq in assembly {}'.format(4140-CDS, tRNA-89, name[0]))
            elif tRNA < 89:
                out_file.write('Prokka found {} fewer CDS and {} less tRNA than the RefSeq in assembly {}'.format(4140-CDS, 89-tRNA, name[0]))
            else:
                out_file.write('Prokka found {} fewer CDS and the same tRNA than the RefSeq in assembly {}'.format(4140-CDS, name[0]))
        else:
            if tRNA > 89:
                out_file.write('Prokka found the same CDS and {} more tRNA than the RefSeq in assembly {}'.format(tRNA-89, name[0]))
            elif tRNA < 89:
                out_file.write('Prokka found the same CDS and {} less tRNA than the RefSeq in assembly {}'.format(89-tRNA, name[0]))
            else:
                out_file.write('Prokka found the same CDS and the same tRNA than the RefSeq in assembly {}'.format(name[0]))
        
        out_file.write('\n \n')

    out_file.write('\n \n \n')


def get_map():
    # gets map of assembly using bowtie2 and and tophat2
    for name in var_names:
        
        # fetches SRA files
        out_file.write('Fetching ' + name[0] + ' SRA fastq files.\n \n')
        sp.call(['mkdir', name[0] + '/transcriptomes/'])
        sp.call(['prefetch', name[4]])
        sp.call(['fastq-dump', '-I', '--split-files', '-O', name[0] + '/transcriptomes/', name[4]])
        
        # renames all files into propper location with _index naming scheme
        sp.call(['mv', name[0] + '/transcriptomes/' + name[4] + '_1.fastq', name[0] + '/transcriptomes/' + name[0] + '_1.fastq'])
        sp.call(['mv', name[0] + '/transcriptomes/' + name[4] + '_2.fastq', name[0] + '/transcriptomes/' + name[0] + '_2.fastq'])
        sp.call(['mv', name[0] + '/annotation/' + name[0] + '.gff', name[0] + '/transcriptomes/' + name[0] + '_index.gff'])
        sp.call(['mv', name[0] + '/fasta/' + name[0] + '_fasta.fna', name[0] + '/transcriptomes/' + name[0] + '_index.fa'])
        
        out_file.write('Mapping expression of ' + name[0] + '\n')
        # runs bowtie2 and tophat2
        sp.call(['bowtie2-build', '--threads', '20', '-f', name[0] + '/transcriptomes/' + name[0] + "_index.fa", name[0] + '/transcriptomes/' + name[0] + '_index'])
        sp.call(['tophat2', '-p', '20', '-o', name[0] + '/transcriptomes/' + name[0] + '_out', name[0] + '/transcriptomes/' + name[0] + '_index', name[0] + '/transcriptomes/' + name[0] + '_1.fastq', name[0] + '/transcriptomes/' + name[0] + '_2.fastq'])
        
        # runs samtools to sort
        out_file.write('running Cufflinks.\n')
        sp.call(['samtools', 'sort', name[0] + '/transcriptomes/' + name[0] + 'out/accepted_hits.bam', name[0] + '/transcriptomes/' + name[0] + 'out/accepted_hits.sorted.bam'])
        # runs cufflinks to get expression
        sp.call(['cufflinks', '-p', '4', '-G', name[0] + '/transcriptomes/' + name[0] + '_index.gff', '-o', name[0] + '/transcriptomes/cufflinks', name[0] + '/transcriptomes/' + name[0] + '_out/' + 'accepted_hits.sorted.bam'])
        out_file.write('Sorted expression to ' + name[0] + '/transcriptomes/cufflinks' + name[0] + 'genes.fpkm_tracking')


# main method
if __name__ == '__main__':
    # opens log file
    out_file = open("pythonwrapper.log", 'w')
    fetch_files()
    num_contigs()
    assembly_length()
    run_prokka()
    diff_prokka()
    get_map()
    # closes log file
    out_file.close()
