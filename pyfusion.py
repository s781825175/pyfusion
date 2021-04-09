#coding=utf-8
import argparse
import os
import time
import subprocess
import re

parser = argparse.ArgumentParser()
parser.add_argument("bam", type=str, help="minimum number of breakpoint-spanning reads required for output")
parser.add_argument("exon", type=str, help="minimum number of breakpoint-spanning reads required for output")
parser.add_argument("twobit", type=str, help="minimum number of breakpoint-spanning reads required for output")
parser.add_argument("important-genelist", type=str, help="the most important gene that you don't want to filter it")
parser.add_argument("--targets", type=str, default='0', help="minimum number of breakpoint-spanning reads required for output")
parser.add_argument("-o", nargs='?', action="store" , type=str, required=True, help="output directory (tumor.bam directory)")
parser.add_argument("-r", nargs='?', const=5, default=5 , type=int, help="minimum number of breakpoint-spanning reads required for output")
parser.add_argument("-m", nargs='?', const=2, default=2 , type=int, help="minimum number of discordant reads required for each candidate fusion")
parser.add_argument("-x", nargs='?', const=5, default=5 , type=int, help="maximum number of breakpoints to examine for each pair of genomic regions")
parser.add_argument("-s", nargs='?', const=1, default=1 , type=int, help="minimum number of reads with the same breakpoint")
parser.add_argument("-f", nargs='?', const=0.9, default=0.9 , type=float, help="minimum fraction of read bases required for alignment to fusion template")
parser.add_argument("-S", nargs='?', const=0.95, default=0.95 , type=float, help="minimum similarity required for read to match fusion template")
parser.add_argument("-k", nargs='?', const=12, default=12 , type=int, help="k-mer size for fragment comparison")
parser.add_argument("-c", nargs='?', const=16, default=16 , type=int, help="minimum size of soft-clipped region to consider")
parser.add_argument("-b", nargs='?', const=500, default=500 , type=int, help="number of bases flanking breakpoint for fusion template")
parser.add_argument("-p", nargs='?', const=4, default=4 , type=int, help="number of threads for blastn search")
parser.add_argument("-a", nargs='?', const=50, default=50 , type=int, help="number of bases flanking breakpoint to provide in output")
parser.add_argument("-e", nargs='?', const=True, default=False , help="disable grouping of input coordinates by column 4 of exons.bed")
parser.add_argument("-v", nargs='?', const=True, default=False , help="disable verbose output")
parser.add_argument("-t", nargs='?', const=True, default=False , help="disable running time output")
parser.add_argument("-F", nargs='?', const=True, default=False , help="force remake of BLAST database for a particular input")
args = parser.parse_args()

#===========================================================================================================
#retrieve depth of region surrounding fusion(s)
def important_gene(important_genelist,gene1,gene2):
    with open(important_genelist,'r',encoding='gbk') as f:
        data = f.readlines()
        if gene1 + '-' + gene2 + '\n' in data or gene2 + '-' + gene1 + '\n' in data:
            return True
        else:
            return False
    

def getNormDepth(bestfusions, buffer, targets, OUTPUTDIR, bam):
    coors = {} #store surrounding depth for each translocation
    for bp_depth in bestfusions:
        for bps in bestfusions[bp_depth]:
            data = bestfusions[bp_depth][bps][0] + bestfusions[bp_depth][bps][1]
            tokens = bps.split(' ')
            bp1 = tokens[0]
            bp2 = tokens[1]
            chr1 = bp1.split(':')[0]
            bp1 = int(bp1.split(':')[1])
            chr2 = bp2.split(':')[0]
            bp2 = int(bp2.split(':')[1])
            var = data.split('\t')
            orient1 = var[7]
            cut = orient1.split(' ')
            if '-' in cut[0]:
                bp1 -= buffer
            if chr1 not in coors:
                coors[chr1] = {}
            coors[chr1][bp1] = tokens[0] 
            if '-' in  cut[1]:
                bp2 -= buffer
            if chr2 not in coors:
                coors[chr2] = {}
            coors[chr2][bp2] = tokens[1]

            insert = "-l " + targets
            if targets == '0':
                insert = ''
            (tmp_error,tmp_data) = subprocess.getstatusoutput("samtools mpileup -ABr {S1}:{S2}-{S3} -Q20 {S4} -d 10000000 {S5}".format(S1=chr1,S2=bp1,S3=bp1+buffer,S4=insert,S5=bam))
            d1 = []
            for line in tmp_data.split('\n'):
                if 'mpileup' in line or 'older than' in line:
                    continue
                depth = line.split('\t')[3]
                d1.append(depth)
            if len(d1) == 0:
                med_depth1 = 0
            else:
                med_depth1 = sorted(d1)[int(len(d1)/2)]
            (tmp_error,tmp_data) = subprocess.getstatusoutput("samtools mpileup -ABr {S1}:{S2}-{S3} -Q20 {S4} -d 10000000 {S5}".format(S1=chr2,S2=bp2,S3=bp2+buffer,S4=insert,S5=bam))
            d2 = []
            for line in tmp_data.split('\n'):
                if 'mpileup' in line or 'older than' in line:
                    continue
                depth = line.split('\t')[3]
                d2.append(depth)
            if len(d2) == 0:
                med_depth2 = 0
            else:
                med_depth2 = sorted(d2)[int(len(d2)/2)]

            final_depth = max(int(med_depth1), int(med_depth2))
            coors[bps] = final_depth

    return coors




#=========================================================================================================== 
def makeBLASTdb(now,storeids4blast,VERBOSE,blast,OUTPUTDIR,bam,blastitor,fusetargets):
    currtime = time.time() - now
    
    #print all unmapped reads to blast database fasta file 
    (blastreads_error,blastreads) = subprocess.getstatusoutput("samtools view -F 2 -L {S1} {S2}".format(S1=fusetargets,S2=bam))
    
    count = 0
    progress = ['|','/','-','\\']
    progressItor = 0 
    tokens_1 = [73,133,89,121,165,181,101,117,153,185,69,137,77,141]
    for line in blastreads.split('\n'):
        tokens = line.split("\t")
        if tokens[0] in storeids4blast:
            count += 1
            if count % 1000 == 0 and VERBOSE == 1:
                progressItor += 1
                if progressItor > 3:
                    progressItor = 0
            Ns = tokens[9].replace('N','')
            if float(len(Ns)/len(tokens[9])) < 0.75:
                continue
            blastitor += 1
            blast.write(">{S1}_{S2}_IP\n{S3}\n".format(S1=tokens[0],S2=blastitor,S3=tokens[9]))
        elif tokens[1] in tokens_1:
            count += 1
            if count % 1000 == 0 and VERBOSE == 1:
                if progressItor > 3:
                    progressItor = 0
            #so makeblastdb does not complain about >40% Ns
            Ns = tokens[9].replace('N','')
            if float(len(Ns)/len(tokens[9])) < 0.75:
                continue
            blastitor += 1
            blast.write(">{S1}_{S2}_UM\n{S3}\n".format(S1=tokens[0],S2=blastitor,S3=tokens[9]))
    print("\n")
    blast.close()

    bdbname = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.blastreads.fa')
    currtime = time.time() - now
    print(currtime)
    print(" Creating blast database...")

    os.system("makeblastdb -in {S1} -dbtype 'nucl' >/dev/null".format(S1=bdbname))

    if VERBOSE == 1:
        print("\n           - Done\n")

#===========================================================================================================       
#calculate fusion coordinates and retrieve genomic sequences surrounding breakpoint using user-definable pad length (default 200)

def getFusionSeq(clipOrder, clipOrder2, bp_, bp2_, bp, bp2, offset, buffer, offset1,twobit,OUTPUTDIR,bam):
    g1start = 0
    g1end = 0
    g2start = 0 
    g2end = 0 
    if clipOrder == clipOrder2:
        if clipOrder == "NC":
            bp_ += offset1
            g1start = bp_ - buffer
            g1end = bp_
            bp2_ += offset
            g2start = bp2_ - buffer
            g2end = bp2_
        else:
            bp_ += offset1
            g1start = bp_ - 1
            g1end = bp_ + buffer - 1
            bp2_ += offset
            g2start = bp2_ - 1
            g2end = bp2_ + buffer - 1
    else:
        if clipOrder == "NC":
            bp_ += offset1
            g1start = bp_ - buffer
            g1end = bp_
            bp2_ += offset
            g2end = bp2_ + buffer - 1
            g2start = bp2_ - 1
        else:
            bp_ += offset1
            g1start = bp_ - 1
            g1end = bp_ + buffer - 1
            bp2_ += offset
            g2end = bp2_
            g2start = bp2_ - buffer
    index1 = bp.split(':')[0]
    index1 = index1 + ':' + str(g1start) + '-' + str(g1end)
    index2 = bp2.split(':')[0]
    index2 = index2 + ':' + str(g2start) + '-' + str(g2end)
    if g1start < 0 or g2start < 0:
        return ">{S1},{S2}\nNNNNNNN".format(S1=index1,S2=index2)
    #extract genomic sequence for putative translocation and write to temporary file
    output = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.tmp.fa')
    os.system("twoBitToFa -noMask {S1}:{S2},{S3} {S4}".format(S1=twobit,S2=index1,S3=index2,S4=output))
    FILE = open(output,'r')
    seqCount = 0 
    seq1 = ''
    seq2 = ''
    data = FILE.readlines()
    for line in data:
        if '>' in line:
            seqCount += 1
            continue
        if '\n' in line:
            line = line[:-1]
        if seqCount == 1:
            seq1 += line
        else:
            seq2 += line
    FILE.close()
    os.system('rm -rf ' + output)

    #merge both sequences to create putative fusion gene sequence
    fusionseq = seq1
    if clipOrder == clipOrder2: #need to flip fragment 2 into reverse complemen t(seq2)
        seq2 = doRevComp(seq2)
    if clipOrder == "NC":
        fusionseq += seq2  #order=fragment 1, fragment 2
    else:
        fusionseq = seq2 + fusionseq #order=fragment 2, fragment 1
    
    return ">{S1},{S2}\n{S3}".format(S1=index1,S2=index2,S3=fusionseq)
#===========================================================================================================       
#BLAST all reads that map to genes involved in fusion region (properly paired, improperly paired, and not paired)
def doBLAST(fusionSeqFa,gene,gene2,buffer,min_len,bam,OUTPUTDIR,BLASTTHREADS,minsimilarity):
    #write fusion gene to temp file for blast search
    query = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.blastquery.fa')
    bquery = open(query,'w')
    bquery.write(fusionSeqFa)
    bquery.close()

    bdbname = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.blastreads.fa')
    (blast_error,blastout) = subprocess.getstatusoutput("blastn -task 'megablast' -query {S1} -db {S2} -outfmt 6 -num_threads {S3} -max_target_seqs 9999999".format(S1=query,S2=bdbname,S3=BLASTTHREADS))
    if blastout == '':
        return ''
    bp_depth = 0 #breakpoint depth
    bpd_SC = 0  #count soft-clipped properly paired read support
    bpd_UM = 0 #count unpaired read support
    bpd_IP = 0 #count soft-clipped improperly paired read support

    proper_pair_depth = 0  #extra support due to properly paired reads that flank breakpoint
    properpairs = {}  #store proper pairs

    storeused = {} #store read ids to ensure each fragment only counted once
    
    percentsim = 100 * minsimilarity

    for line in blastout.split('\n'):
        tokens = line.split('\t')
        id_ = tokens[1]
        readtype = id_.split('_')[-1]
        id_ = id_.split('_')[0] #store read type code (SC=soft-clipped, UM=unmapped, IP=improperly paired)
        percent = tokens[2]
        align_len = tokens[3]
        left = tokens[6]
        right = tokens[7]

        if float(percent) >= percentsim and float(align_len) >= min_len:
            if id_ in storeused:
                properpairs[id_][left] = right
                continue
            if int(min(left,right)) < buffer - 15 and int(max(left,right)) > buffer + 15:
                bp_depth += 1
                if readtype == 'SC':
                    bpd_SC += 1
                if readtype == 'UM':
                    bpd_UM += 1
                if readtype == 'IP':
                    bpd_IP += 1
                storeused[id_] = id_ #store id
            if id_ not in properpairs:
                properpairs[id_] = {}
            properpairs[id_][left] = right #group proper pairs   
    #=====================================  
    #count proper pairs that flank breakpoint  
    for id_ in properpairs:
        size = properpairs[id_]
        if len(size) != 2:
            continue
        coors = []
        for left in properpairs[id_]:
            coors.append(int(left))
            coors.append(int(properpairs[id_][left]))
        if min(coors) < (buffer - 85) and max(coors) > (buffer + 85):
            proper_pair_depth += 1
    #=====================================   
      
    #print break-point depth
    data = "{S1}\t{S2}\t{S3}\t{S4}\t{S5}\t".format(S1=bp_depth,S2=bpd_SC,S3=bpd_UM,S4=bpd_IP,S5=proper_pair_depth)
    #print proper pair depth
    return data   

            


        



#=========================================================================================================== 
#parse cigar string    
def parsecigar(cigar):
    cigarItor = 0
    var = cigar[cigarItor:cigarItor+1]
    num = ""
    #extract number and cigar code  
    while re.match(r'^[+-]?\d+$',var) and cigarItor < len(cigar):
        cigarItor += 1
        num = num + var
        var = cigar[cigarItor:cigarItor+1]
    
    cigar = cigar[cigarItor+1:]

    data = []
    data.append(num)
    data.append(var)
    data.append(cigar)

    return data


#============================================================================================================
#given start coordinate of a read, find and return the closest gene----------
def findGene(pos,exonstart,coors2gene,chr,getothercoor):
    start_coor = BinSearch(pos,exonstart)
    if int(start_coor) - start_coor != 0:
        p1 = int(start_coor - 0.5)
        p2 = int(start_coor + 0.5)
        if p2 == len(exonstart):
            pos2 = 0
        else:
            pos2 = int(exonstart[p2])
        pos1 = int(exonstart[p1])
        if abs(pos1 - pos) < abs(pos2 - pos):
            start_coor = pos1
        else:
            start_coor = pos2
    else:
        start_coor = exonstart[start_coor]
    start_coor_end = getothercoor[chr + ' ' + str(start_coor)]
    gene_start = start_coor
    gene_end = start_coor_end
    if (chr + ' ' + str(start_coor) + ' ' + start_coor_end) in coors2gene:
        gene = coors2gene[chr + ' ' + str(start_coor) + ' ' + start_coor_end]
    else:
        gene = 'no gene'
    return [gene_start,gene_end,gene]

def BinSearch(target,exonstart):
    posmin = 0
    posmax = len(exonstart)-1
    if cmpFunc(0,exonstart,target) > 0:
        return -0.5
    if cmpFunc(posmax,exonstart,target) < 0:
        return posmax + 0.5
    while 1:
        mid = int((posmin + posmax) / 2)
        result = cmpFunc(mid,exonstart,target)
        if result < 0:
            if mid == posmin and posmax != posmin:
                posmin = posmax
                continue
            if mid == posmin:
                return mid + 0.5
            posmin = mid
        elif result > 0:
            if mid == posmax and posmax != posmin:
                posmax = posmin
                continue
            if mid == posmax:
                return mid - 0.5
            posmax = mid
        else:
            return mid



def cmpFunc(index,arrayRef,target):
    item = int(arrayRef[index])
    if item > target:
        return 1
    elif item < target:
        return -1
    else:
        return 0

#store consensus sequence for soft-clipped reads for a given breakpoint
def getSeqArray(consensus, gene, bp, num1, adjustlen, num2, clipped, screads, clipped2):
    i = num1
    leng = adjustlen
    j = num2
    leng2 = clipped

    seqarray = []
    
    if bp in screads[gene]:
        if gene not in consensus:
            consensus[gene] = {}
        if bp not in consensus[gene]:
            consensus[gene][bp] = []
        seqarray = consensus[gene][bp]
    while i < leng:
        if i < 0:
            continue
        ch = clipped2[j:j+1]
        if i + 1 > len(seqarray):
            for hhh in range(i+1-len(seqarray)):
                seqarray.append([0,0,0,0])
        if ch == "A":
            seqarray[i][0] += 1
        if ch == "T":
            seqarray[i][1] += 1
        if ch == "C":
            seqarray[i][2] += 1
        if ch == "G":
            seqarray[i][3] += 1
        j+=1
        if j > leng2:
            break
        i += 1
    return seqarray

#===========================================================================================================       
#compare clipped read to non-clipped read using k-mer hashtable
def compKmers(kmers1,kmers1rc,read2,k,threshold,read1,doOffset,clip1,clip2,clipsize):
    sumk = 0 #count forward hits
    sumkrc = 0 #count reverse complement hits
    read1len = len(read1)
    offset = 0
    offsetrc = 0
    firstmatch = 0
    firstmatchrc = 0

    for itor in range(len(read2) - k):
        km = read2[itor:itor+k]
        if km in kmers1:
            sumk += 1
            if firstmatch == 0:
                if clip1 == "CN" and clip2 == "NC":
                    offset = itor - kmers1[km]
                elif clip1 == "NC" and clip2 == "CN":
                    offset = read1len - kmers1[km] - (len(read2) - itor) 
                firstmatch = 1
        if km in kmers1rc:
            sumkrc += 1
            if firstmatch == 0:
                if clip1 == "NC" and clip2 == "NC":
                    offset = itor - kmers1rc[km]
                elif clip1 == "CN" and clip2 == "CN":
                    offset = read1len - kmers1rc[km] - (len(read2) - itor) 
                firstmatchrc = 1
    
    if sumk >= max(threshold,clipsize-k):   #threshold = min(25,0.5*(min(len(read2),len(read1))- k))
        tmp = [offset,'F']
        return tmp
    elif sumkrc >= max(threshold,clipsize-k):
        tmp = [offsetrc,'RC']
        return tmp
    else:
        return [9999999]

#perform bp correction by checking subsequence against reference
def bp_correction(read2,same,offsetbp2,bp2_,chr2,OUTPUTDIR,bam,twobit):
    offsetbp1 = 0
    index = ""
    r2seq = ""
    r2 = read2.split(' ')
    if same[0] > 0:
        #grab from read 2 downstream of cutpoint
        r2seq = r2[1][:same[0]]
        index = chr2 + ":" + str(bp2_ - 1) + "-" + str(bp2_ + same[0] - 1)
        return [offsetbp1,offsetbp2]

    else:
        #grab from read 2 upstream of cutpoint
        r2seq = r2[0][len(r2[0])+same[0]:len(r2[0])]
        index = chr2 + ":" + str(bp2_ + same[0])
        return [offsetbp1,offsetbp2]
    
    if len(r2seq) < 3:
        return [offsetbp1,offsetbp2]
    output = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.tmp.fa')
    os.system("twoBitToFa -noMask {S1}:{S2} {S3}".format(S1=twobit,S2=index,S3=output))

    #open temporary fasta file containing genomic sequences
    FILE = open(output,'r')
    seq1 = ''
    data1 = FILE.readlines()
    for line in data1:
        line = line[:-1]
        if '>' in line:
            continue
        seq1 = seq1 + line
    FILE.close()
    #remove temporary file
    os.system("rm -rf {S1}".format(S1=output))

    if r2seq == seq1:
        #adjust breakpoint 1
        if same[1] == 'F':
            offsetbp1 = offsetbp2 * -1
        else:
            offsetbp1 = offsetbp2
        offsetbp2 = 0 
    data = [offsetbp1,offsetbp2]
    return data 
        




#===========================================================================================================       
#convert read into kmers; store as hashtable
def getKmers(read1,k):
    kmers = {}
    for itor in range(len(read1) - k):
        km = read1[itor:itor+k]
        if km not in kmers:
            kmers[km] = itor
    return kmers


#===========================================================================================================
#return reverse complement of read
def doRevComp(read1):
    rc = ""
    read2 = read1[::-1]
    for c in read2:
        if c == "A":
            c = "T"
        elif c == "T":
            c = "A"
        elif c == "C":
            c = "G"
        elif c == "G":
            c = "C"
        rc += c
    return rc

#get consensus sequence for soft-clipped reads for a given breakpoint
def getConsensus(consensus,gene,bp_):
    seqarray = consensus[gene][bp_]
    sq = ''
    for j in range(len(seqarray)):
        mx_val = 0 
        mx_itor = 0
        mx_char = 'A'
        for r in range(4):
            if seqarray[j][r] > mx_val:
                mx_val = seqarray[j][r]
                mx_itor = r
        if mx_val == 0:
            continue
        if mx_itor == 1:
            mx_char = "T"
        if mx_itor == 2:
            mx_char = "C"
        if mx_itor == 3:
            mx_char = "G"
        sq = sq + mx_char
    return sq
        

def main():
    opts = vars(args)

    bam = opts['bam'] #bam file
    exon = opts['exon'] #exon coordinate bed file [chr <tab> start <tab> end <tab> genename]
    twobit = opts['twobit'] #2bit genome file (e.g., hg19.2bit)
    important_genelist = opts['important-genelist'] #important-genelist
    targets = 0 #restrict analysis to targeted regions (.bed) [use 0 if non-targeted]
    for opt in opts:
        if opt == 'targets':
            targets = opts['targets']
            print("Restrict results to regions specified in: $targets\n")
        if opt == 'o':
            if os.path.isdir(opts['o']):
                print('Output directory: ' + opts['o'])
            else:   print(opts['o'] + "is not a directory. Abort!\n")
        if opt == 'r':
            if opts['r'] < 1:
                print("Minimum breakpoint-spanning reads needed must be >0. Abort!\n")
            else:   print("Minimum breakpoint-spanning reads needed for output set to: " + str(opts['r']))
        if opt == 'm':
            if opts['m'] < 1:
                print("Minimum discordant pairs for candidate fusion must be >0. Abort!\n")
            else:   print("Minimum discordant pairs for candidate fusion set to: " + str(opts['m']))
        if opt == 'x':
            if opts['x'] < 1:
                print("Maximum breakpoints for each gene pair must be >0. Abort!\n")
            else:   print("Maximum breakpoints for each gene pair set to: " + str(opts['x']))
        if opt == 's':
            if opts['s'] < 1:
                print("Minimum number of reads with same breakpoint must be >0. Abort!\n")
            else:   print("Minimum number of reads with same breakpoint set to: " + str(opts['s']))
        if opt == 'f':
            if opts['f'] < 1 and opts['f'] > 0:
                print("Minimum fraction of matching bases for fusion template must be a fraction (0-1). Abort!\n")
            else:   print("Minimum fraction of matching bases for fusion template set to: " + str(opts['f']))
        if opt == 'S':
            if opts['S'] < 1 and opts['S'] > 0:
                print("Minimum similarity between reads and fusion template must be a fraction (0-1). Abort!\n")
            else:   print("Minimum similarity between reads and fusion template set to: " + str(opts['S']))
        if opt == 'k':
            if opts['k'] < 1:
                print("k-mer size for fragment comparison must be >0. Abort!\n")
            else:   print("k-mer size for fragment comparison set to: " + str(opts['k']))
        if opt == 'c':
            if opts['c'] < 1:
                print("Minimum size of soft-clipped read must be >0. Abort!\n")
            else:   print("Minimum size of soft-clipped read set to: " + str(opts['c']))
        if opt == 'b':
            if opts['b'] < 1:
                print("Number of bases flanking breakpoint of fusion template must be >0. Abort!\n")
            else:   print("Number of bases flanking breakpoint of fusion template set to: " + str(opts['b']))
        if opt == 'p':
            if opts['p'] < 1:
                print("Number of threads for blastn search must be >0. Abort!\n")
            else:   print("Number of threads for blastn search set to: " + str(opts['p']))
        if opt == 'a':
            if opts['a'] < 1:
                print("Number of bases flanking breakpoint in output must be >0. Abort!\n")
            else:   print("Number of bases flanking breakpoint in output set to: " + str(opts['a']))
        if opt == 'e':  print("Disable grouping of input coordinates (e.g., allow fusions between exons)\n")
        if opt == 'v':  print("Suppress verbose output\n")
        if opt == 't':  print("Suppress running time output\n")
        if opt == 'C':  print("Disable addition of \'chr\' prefix to chromosome names\n")
        if opt == 'F':  print("Force remake of blast database\n")

    #PARAMETERS============================================================================================

    MINSPANNINGREADS = args.r #minimum number of breakpoint-spanning reads needed to output a fusion
    MINIMPROPERREADS = args.m #minimum number of discordant reads needed to consider putative fusion
    MAXBPS2EXAMINE = args.x #maximum number of putative breakpoints to consider for each unique gene pair
    MINBPSUPPORT = args.s #minimum number of reads spanning putative breakpoint to proceed
    MINBLASTFRAC = args.f #fraction of read length to determine minimum alignment length for matching sequences
    VERBOSE = 1
    TIME = 1
    USE_ALL_COORS = 0
    if args.v:
        VERBOSE = 0 #if 1, verbose output, including running time; o.w. only results will be printed
    if args.t: 
        TIME = 0    #if 1, print running time 
    if args.e:
        USE_ALL_COORS = 1 #if 0, cluster input coordinates (second argument) by 4th column of bed file (e.g., gene symbol); o.w. treat all coordinates separately (i.e., each line in bed file will be distinct; e.g., to find fusions between exons, rather than between genes)
    makeblastdb = 1 #will set to 0 if a $bam.nhr file is detected
    FORCEREMAKE = 0 #if 1, force remake of blast database if one already exists
    if args.F:
        FORCEREMAKE = 1
    disablechrprefix = 0 #if 1, do not add "chr" prefix to chromosome names
    k = args.k
    buffer = args.b
    pad = args.a
    clipsize = args.c
    minsimilarity = args.S
    BLASTTHREADS = args.p
    READLEN = 101
    READLENBIT = 0
    if args.o:
        OUTPUTDIR = args.o
    if k > clipsize:
        print("k-mer size {S1} cannot be greater than minimum size of soft-clipped read {S1}. Abort!\n".format(S1=k,S2=clipsize))
    #==============================================================================================================
    #function


    #write parameters used and arguments to disk
    params = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.parameters.txt')
    with open(params,'w') as output0:
        v_ = VERBOSE
        if v_ == 1:
            v_ = 0
        t_ = TIME
        if t_ == 1:
            t_ = 0
        Ta_ = targets
        if Ta_ == '0':
            Ta_ = 'null'
        output0_info = ["Parameter (symbolic argument)\tValue", "\nInput BAM:\t" + bam, "\nGenomic coordinates:\t" + exon, "\nReference genome:\t" + twobit, "\nTargeted regions:\t" + Ta_,
        "\nOutput directory (o):\t" + OUTPUTDIR, "\nMinimum breakpoint-spanning reads for each fusion (r):\t" + str(MINSPANNINGREADS), "\nMinimum discordant reads for each candidate fusion (m):\t"+ str(MINIMPROPERREADS),
        "\nMaximum breakpoints to examine for each read pair (x):\t" + str(MAXBPS2EXAMINE), "\nMinimum reads with same breakpoint required (s):\t" + str(MINBPSUPPORT), "\nMinimum fraction of read bases for fusion template alignment (f):\t" + str(MINBLASTFRAC),
        "\nMinimum similarity for read to match fusion template (S):\t" + str(minsimilarity), "\nK-mer size for fragment comparison (k):\t" + str(k), "\nMinimum size of soft-clipped read region (c):\t" + str(clipsize),
        "\nNumber of bases flanking breakpoint for fusion template (b):\t" + str(buffer), "\nNumber of threads for blastn search (p):\t" + str(BLASTTHREADS), "\nNumber of bases flanking breakpoint to provide in output (a):\t" + str(pad),
        "\nDisable grouping of input coordinates by column 4 of input bed (e):\t" + str(USE_ALL_COORS), "\nDisable verbose output (v):\t" + str(v_), "\nDisable running time output (t):\t" + str(t_), "\nForce remake of blast database (F):\t" + str(FORCEREMAKE)]
        output0.writelines(output0_info)

    #if blast database already exists, don't make it again...

    #write reads in fasta format to create blast database of improperly paired, soft-clipped, and unmapped reads
    blastdb = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.blastreads.fa')
    if os.path.isfile(blastdb) and FORCEREMAKE == 0:
        makeblastdb = 0
        print("Old blast database detected. Use -F argument to force remake.\n")
    #==============================================================================================================
    exon_start = [] #store exon start coordinates
    coors2gene = {} #store gene symbol for every unique exon [format: chr_start_end]
    chrom = {} #store index for each chromosome
    getothercoor = {} #get end coordinate if input is start and vice versa [format: gene_start or gene_end]
    blastitor = 0 #iterate so every subject has unique identifier
    now = time.time() #running time

    #open(NULL, ">", File::Spec->devnull); #suppress standard error for mpileup

    #read in exon coors: store coordinates in sorted arrays for each chromosome
    with open(exon,'r') as FILE:
        data = FILE.readlines()
        for line in data:
            tokens = line.split('\t')
            chr = tokens[0]
            start = tokens[1]
            end = tokens[2]
            gene = tokens[3][:-1]
            if chr not in chrom:
                chrom[chr] = len(chrom)
            chrindex = chrom[chr] #retrieve chromosome index
            if chrindex >= len(exon_start):
                exon_start.append([])
            exon_start[chrindex].append(int(start)) #add start coordinates
            if (chr + ' ' + start + ' ' + end) not in coors2gene:
                insert = gene
                if USE_ALL_COORS == 1:
                    insert = gene + ":" + start
                coors2gene[chr + ' ' + start + ' ' + end] = insert
            getothercoor[chr + ' ' + start] = end
            getothercoor[chr + ' ' + end] = start
        if VERBOSE == 1:
            print("Done loading genomic coordinates\n")
        #sort coordinates for each chromosome

    for itor in exon_start:
        itor.sort()

    if VERBOSE == 1:
        print("Done sorting genomic coordinates\n")

    impropname = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.discordantpair.details.txt')
    depthname =  OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.discordantpair.depth.txt')
    output = open(impropname,'w')
    output2 = open(depthname,'w')
    insert = "-L" + targets
    if targets == '0':
        insert = ''
    (imp_error,imp_pairs) = subprocess.getstatusoutput("samtools view -F 2 {S1} {S2} | awk '$2 ~ /81|161|97|145|65|129|113|177/'".format(S1=insert,S2=bam))
    output.write("GENE1\tGENE2\tFUSION\tCHR1\tCHR2\tGENE1_START\tGENE2_START\tDIST\tREAD1_POS\tREAD2_POS\n")
    output2.write("FUSION\tDEPTH\n")
    if makeblastdb == 1:
        blast = open(blastdb,'w')
    #collect all ids of mapped improper pairs in targeted region 
    insert = "-L " + targets
    if targets == '0':
        insert = ''


    


    count = 0

    prevgene = {} #store genes already seen [format: chr_start]
    prevread = {} #store reads already seen
    fusepairs = {} #store unique read pairs per fusion
    genes2exons = {} #store exon coordinates and chromosome for each candidate fusion gene
    storeids4blast = {} #store read ids of properly paired reads that passed filtration criteria; used to build blast database
    progress = ['|','/','-','\\']
    progressItor = 0 

    if VERBOSE == 1:
        currtime = time.time() - now
        print("\n   Time                 Current step\n")
        print("---------- --------------------------------------\n")
        print(currtime)
        print(" Analyzing improperly paired reads...\n")

    for line in imp_pairs.split('\n'):
        tokens = line.split('\t')
        count += 1
        if count % 1000 == 0 and VERBOSE == 1:
            progressItor += 1
            if progressItor > 3:
                progressItor = 0
        Q = 0
        quality = tokens[10]
        for i in quality:
            Q += ord(i) - 33
        if Q/len(quality) < 30:
            continue
        id = tokens[0]
        #skip read ids that have already been processed
        if id in prevread:
            continue
        else:
            prevread[id] = id
        
        #read1
        chr = tokens[2]
        pos = int(tokens[3])
        #read2
        chrb = tokens[6]
        if chrb == '=':
            chrb = chr
        posb = int(tokens[7])
        gene = ''
        gene_start = ''
        gene_end = ''
        geneb = ''
        gene_startb = ''
        gene_endb = ''
        if (chr + ' ' + str(pos)) not in prevgene:
            if chr in chrom:
                exonstart = exon_start[chrom[chr]]
            else:
                continue
            [gene_start,gene_end,gene] = findGene(pos,exonstart,coors2gene,chr,getothercoor)
            if gene == 'no gene':
                continue
            prevgene[chr + ' ' + str(pos)] = gene + '|' + str(gene_start) + '|' + str(gene_end)
        else:
            tmp = prevgene[chr + ' ' + str(pos)]
            var = tmp.split('|')
            gene = var[0]
            gene_start = var[1]
            gene_end = var[2]
#----------------find closest start coordinate to read2 (and corresponding end)------------
        if (chrb + ' ' + str(posb)) not in prevgene:
            if chrb in chrom:
                exonstart = exon_start[chrom[chrb]]
            else:
                continue
            [gene_startb,gene_endb,geneb] = findGene(posb,exonstart,coors2gene,chrb,getothercoor)
            if geneb == 'no gene':
                continue
            prevgene[chrb + ' ' + str(posb)] = geneb + '|' + str(gene_startb) + '|' + str(gene_endb)
        else:
            tmp = prevgene[chrb + ' ' + str(posb)]
            var = tmp.split('|')
            geneb = var[0]
            gene_startb = var[1]
            gene_endb = var[2]       
#-----------------if candidate translocation---------------------------------
        if gene != geneb:  #discordant gene pair...
            if gene == '' or geneb == '':
                continue
            dist = 0
            if chr == chrb:
                dist = abs(int(gene_start) - int(gene_startb))
            
            fusion = gene + '-' + geneb
            output.write(gene+'\t'+geneb+'\t'+fusion+'\t'+chr+'\t'+chrb+'\t'+str(gene_start)+'\t'+str(gene_startb)+'\t'+str(dist)+'\t'+str(pos)+'\t'+str(posb)+'\n')
            
            #store chromosome and start position of each exon corresponding to candidate fusion gene
            minstart = max(1,min(int(gene_start),int(pos)-300))
            minstartb = max(1,min(int(gene_startb),int(posb)-300))
            maxend = max(int(gene_end),int(pos)+300)
            maxendb = max(int(gene_endb),int(posb)+300)
            if gene not in genes2exons:
                genes2exons[gene] = {}
            if geneb not in genes2exons:
                genes2exons[geneb] = {}
            if chr + ' ' + str(minstart) not in genes2exons[gene]:
                genes2exons[gene][chr + ' ' + str(minstart)] = maxend
            elif genes2exons[gene][chr + ' ' + str(minstart)] < maxend:
                genes2exons[gene][chr + ' ' + str(minstart)] = maxend
            if chrb + ' ' + str(minstartb) not in genes2exons[geneb]:
                genes2exons[geneb][chrb + ' ' + str(minstartb)] = maxendb
            elif genes2exons[geneb][chrb + ' ' + str(minstartb)] < maxendb:
                genes2exons[geneb][chrb + ' ' + str(minstartb)] = maxendb
            #eliminate redundant fragments (same start for both reads) for depth statistics 
            if gene + ' ' + geneb not in fusepairs:
                fusepairs[gene + ' ' + geneb] = {}
            if str(pos) + ' ' + str(posb) not in fusepairs[gene + ' ' + geneb]:
                fusepairs[gene + ' ' + geneb][str(pos) + ' ' + str(posb)] = str(pos) + ' ' + str(posb)
                storeids4blast[id] = id
    output.close()
    if VERBOSE == 1:
        print('\n')

#-------------print out depth statistics for mapped improper pairs-----------
    genes = {}
    pairs = {}
    depth = {}
    depth_sort = []
    fusetargets = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.fusiontargets.bed')

    fusiontargets = open(fusetargets,'w')

    for g in fusepairs:
        d = fusepairs[g]
        depth_sort.append([g,len(d)])
        depth[g] = d
    count = 0


    depth_sort = sorted(depth_sort,key = lambda x:x[1],reverse=True)
    for g in depth_sort:
        count += 1
        if g[1] < MINIMPROPERREADS:
            continue
        gs = g[0].split(' ')
        if gs[0] not in genes:
            genes[gs[0]] = gs[0]
            f_chr = ''
            f_min = 9999999999
            f_max = 0
            for key in genes2exons[gs[0]]:
                f_vars = key.split(' ')
                f_chr = f_vars[0]
                start = int(f_vars[1])
                if start < f_min:
                    f_min = start
                end = genes2exons[gs[0]][key]
                if end > f_max:
                    f_max = end
            fusiontargets.write(f_chr+'\t'+str(f_min)+'\t'+str(f_max)+'\t'+g[0]+'\n')
        if gs[1] not in genes:
            genes[gs[1]] = gs[1]
            f_chr = ''
            f_min = 9999999999
            f_max = 0
            for key in genes2exons[gs[1]]:
                f_vars = key.split(' ')
                f_chr = f_vars[0]
                start = int(f_vars[1])
                if start < f_min:
                    f_min = start
                end = genes2exons[gs[1]][key]
                if end > f_max:
                    f_max = end
            fusiontargets.write(f_chr+'\t'+str(f_min)+'\t'+str(f_max)+'\t'+g[0]+'\n')
        output2.write(g[0]+'\t'+str(len(depth[g[0]]))+'\n')
        if gs[0] not in pairs:
            pairs[gs[0]] = {}
        if gs[1] not in pairs:
            pairs[gs[1]] = {}
        pairs[gs[0]][gs[1]] = depth[g[0]]
        pairs[gs[1]][gs[0]] = depth[g[0]]
    fusiontargets.close()
    output2.close()
#===============================================================================================
#------------infer break-point by soft-clipping, and increment depth of fusion pairs accordingly
    count = 0
#open my $sclipped, '-|', "samtools view -F 12 -L $fusetargets $bam | awk '\$2 ~ /99|147|83|163|67|131|115|179/' | awk '\$6 ~ /S/'" or die;
    (sclipped_error,sclipped) = subprocess.getstatusoutput("samtools view -F 12 -L {S1} {S2} | awk '$2 ~ /99|147|83|163|67|131|115|179|81|161|97|145|65|129|113|177/' | awk '$6 ~ /S/'".format(S1=fusetargets,S2=bam))
    breakpoints = {}
    screads = {}
    consensus = {}
    if VERBOSE == 1:
        currtime = time.time() - now
        print(currtime)
        print(" Analyzing soft-clipped reads...\n")
    for line in sclipped.split('\n'):
        tokens = line.split('\t')
        cigar = tokens[5]
        if 'D' in cigar or 'I' in cigar:
            continue
        
        count += 1
        if count % 100 == 0 and VERBOSE == 1:
            progressItor += 1
            if progressItor > 3:
                progressItor = 0
        
        tmp = cigar.replace('S','')
        
        if len(tmp) != len(cigar) - 1:
            continue
        
        parsed = parsecigar(cigar)
        parsed2 = parsecigar(parsed[2])
        matched = 0 #cigar string 'M'
        skipped = 0 #cigar string 'S'

        #if forward orientation, properly paired, and with match followed by clip, go to next..
        if tokens[1] == '99' or tokens[1] == '163':
            if parsed[1] == "M":
                continue
        #if reverse orientation, properly paired, and with clip followed by match, go to next..
        if tokens[1] == '147' or tokens[1] == '83':
            if parsed[1] == "S":
                continue
        #only keep clipped regions at least clipsize (default 16 bases; only 1 in 4.3B by random chance)
        if parsed[1] == 'S':
            skipped = parsed[0]
            matched = int(parsed2[0])
            if int(parsed[0]) < clipsize:
                continue
        elif parsed2[1] == 'S':
            skipped = parsed2[0]
            matched = int(parsed[0])
            if int(parsed2[0]) < clipsize:
                continue
        s_chr = tokens[2]
        exonstart = exon_start[chrom[s_chr]]
        pos = int(tokens[3])
        gene = ''
        gene_start = ''
        gene_end = ''
        [gene_start,gene_end,gene] = findGene(pos,exonstart,coors2gene,s_chr,getothercoor)
        if gene == 'no gene':
            continue
        if gene not in genes:
            continue
        if makeblastdb == 1:
            #so makeblastdb does not complain about >40% Ns
            Ns = tokens[9].replace('N','')
            if float(len(Ns)/len(tokens[9])) >= 0.75:
                blast.write(">{S1}_{S2}_SC\n{S3}\n".format(S1=tokens[0],S2=blastitor,S3=tokens[9]))
                blastitor += 1
        #print to blast database
        bp = pos
        if parsed[1] == "M":
            bp += matched - 1
            if gene not in screads:
                screads[gene] = {}
            if bp not in screads[gene]:
                screads[gene][bp] = []
            if len(screads[gene][bp]) > 0:
                prevLen = screads[gene][bp][0]
            else:
                prevLen = 0
            clipped = tokens[9][matched:]
            notclipped = tokens[9][:matched]
            middle = len(tokens[9])/2

            if abs(len(clipped)) -  middle <= abs(prevLen - middle):
                temp = [len(clipped),notclipped + ' ' + clipped, "NC", tokens[0]]
                screads[gene][bp] = temp
                adjustlen = 500
                #store consensus sequence for soft-clipped portion of a given breakpoint
                seqarray = getSeqArray(consensus, gene, bp, 0, adjustlen - clipsize, 0, len(clipped), screads, clipped)
                consensus[gene][bp] = seqarray
        else:
            if gene not in screads:
                screads[gene] = {}
            if bp not in screads[gene]:
                screads[gene][bp] = []
            if len(screads[gene][bp]) > 0:
                prevLen = screads[gene][bp][0]
            else:
                prevLen = 0
            clipped = tokens[9][:len(tokens[9]) - matched]
            notclipped = tokens[9][(len(tokens[9]) - matched):]
            middle = len(tokens[9])/2
            if abs(len(clipped)) -  middle <= abs(prevLen - middle):
                temp = [len(clipped),clipped + ' ' + notclipped, "CN", tokens[0]]
                screads[gene][bp] = temp
                adjustlen = 500
                #store consensus sequence for soft-clipped portion of a given breakpoint
                seqarray = getSeqArray(consensus, gene, bp, adjustlen - clipsize - len(clipped), adjustlen - clipsize, 0, len(clipped), screads, clipped)
                consensus[gene][bp] = seqarray
        if gene not in breakpoints:
            breakpoints[gene] = {}
        if s_chr + ':' + str(bp) not in breakpoints[gene]:
            breakpoints[gene][s_chr + ':' + str(bp)] = 0
        breakpoints[gene][s_chr + ':' + str(bp)] += 1
    
    if VERBOSE == 1:
        print("\n")
#------------------------------------------------------------------------------
#finish writing reads to disk, then make blast database if not already made
    if makeblastdb == 1:
        makeBLASTdb(now,storeids4blast,VERBOSE,blast,OUTPUTDIR,bam,blastitor,fusetargets)
#------------------------------------------------------------------------------
#sort soft-clipped reads by breakpoint, identify "best" breakpoints and analyze
    usedpair = {}
    bestfusions = {}

    fusions = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.fusions.txt')

    output3 = open(fusions, 'w')

    fusionseqs = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.fusionseqs.fa')

    output4 = open(fusionseqs, 'w')

    fusions_bed = OUTPUTDIR + '/' + bam[::-1].split('/',1)[0][::-1].replace('bam','factera.fusions.bed')

    output5 = open(fusions_bed, 'w')

    if VERBOSE == 1:
        currtime = time.time() - now
        print(currtime)
        print(" Validating candidate fusions...\n")
#==================================================================================================

    progress_itor = 0
    progressItor = 0
    prev_prog = -1
    fusioncons = ''
    for gene in breakpoints:
        rankedlist = {}
        #Rank all breakpoints for partner genes wrt read depth
        for gene2 in pairs[gene]:
            if gene2 not in breakpoints:
                continue
            for bp2 in sorted(breakpoints[gene2], key=breakpoints[gene2].__getitem__, reverse=True):
                if len(breakpoints[gene2]) == 0:
                    continue
                depth = breakpoints[gene2][bp2]
                if depth not in rankedlist:
                    rankedlist[depth] = {}
                if gene2 not in rankedlist[depth]:
                    rankedlist[depth][gene2] = {}
                rankedlist[depth][gene2][bp2] = bp2
        
        itor = 0
        for bp in sorted(breakpoints[gene], key=breakpoints[gene].__getitem__, reverse=True):
            itor += 1
            switch = 0
            if itor > MAXBPS2EXAMINE:
                break #only examine the x most abundant putative breakpoints for gene 1
            count = breakpoints[gene][bp]
            if count < MINBPSUPPORT:
                continue #skip if <x reads support
            bp_ = bp
            bp_ = bp.split(':')[1]
            chr_ = bp.split(':')[0]
            if gene not in screads:
                continue
            bp_ = int(bp_)
            if bp_ not in screads[gene]:
                continue
            c_id = screads[gene][bp_][3]
            clipOrder = screads[gene][bp_][2]
            read1 = screads[gene][bp_][1]
            read1b = read1
            readLen = len(read1) - 1
            if READLENBIT == 0:
                READLEN = readLen
                READLENBIT = 1
            read1cut = screads[gene][bp_][0]
            if clipOrder == "NC" :
                read1 = read1[:len(read1)-read1cut-1]
                read1btemp = read1b[len(read1b)-read1cut:]
                read1b = getConsensus(consensus,gene,bp_)
                read1b = read1b[:len(read1btemp)]
            else:
                read1 = read1[read1cut + 1:]
                read1btemp = read1b[:read1cut]
                read1b = getConsensus(consensus,gene,bp_)
                read1b = read1b[:1 + len(read1b) - len(read1btemp)]
            kmers1 = getKmers(read1,k)
            kmers1rc = getKmers(doRevComp(read1),k)
            kmers1b = getKmers(read1b,k)
            kmers1rcb = getKmers(doRevComp(read1b),k)

            itor2 = 0
            for depth in sorted(rankedlist, reverse=True):
                for gene2 in rankedlist[depth]:
                    for bp2 in rankedlist[depth][gene2]:
                        itor2 += 1
                        if important_gene(important_genelist,gene,gene2) == False:
                            if itor2 > MAXBPS2EXAMINE:
                                break #examine x most abundant putative breakpoints paired with gene 1
                            count2 = breakpoints[gene2][bp2]
                            if count2 < MINBPSUPPORT:
                                continue #skip if <x reads support
                        
                        progress_itor += 1
                        if progress_itor % 100 == 0 and VERBOSE == 1:
                            progressItor += 1
                            if progressItor > 3:
                                progressItor = 0

                        bp2_ = bp2.split(':')[1]
                        chr2 = bp2.split(':')[0]
                        if gene2 not in screads:
                            continue
                        bp2_ = int(bp2_)
                        if bp2_ not in screads[gene2]:
                            continue
                        id2 = screads[gene2][bp2_][3]
                        clipOrder2 = screads[gene2][bp2_][2]
                        read2 = screads[gene2][bp2_][1]
                        read2b = read2
                        read2cut = screads[gene2][bp2_][0]
                        if clipOrder2 == "NC": #excise clipped read using order bit
                            read2temp = read2[len(read2) - read2cut:]
                            read2b = read2b[:len(read2b) - read2cut - 1]
                            read2 = getConsensus(consensus,gene2,bp2_)
                            read2 = read2[:len(read2temp)]
                        else:
                            read2temp = read2[:read2cut]
                            read2b = read2b[read2cut+1:]
                            read2 = getConsensus(consensus,gene2,bp2_)
                            read2 = read2[1+len(read2)-len(read2temp):]
                        
                        cmpthresh = 0.5
                        threshold = min(25,cmpthresh*(min(len(read2),len(read1))- k))
                        thresholdb = min(25,cmpthresh*(min(len(read2b),len(read1b))- k))
                        doOffset = 1 #if 1, allow offset adjustment

                        same = compKmers(kmers1,kmers1rc,read2,k,threshold,read1,doOffset,clipOrder,clipOrder2,clipsize)
                        clipOrderb = "NC"
                        clipOrder2b = "NC"
                        if clipOrder == "NC":
                            clipOrderb = "CN"
                        if clipOrder2 == "NC":
                            clipOrder2b = "CN"
                        sameb = compKmers(kmers1b,kmers1rcb,read2b,k,thresholdb,read1b,doOffset,clipOrderb,clipOrder2b,clipsize)
                        target_gene = ['ALK','RET','ROS1']
                        if same[0] != 9999999 and sameb[0] != 9999999:
                            if '{S1}_{S2}_{S3}_{S4}_{S5}_{S6}'.format(S1=gene,S2=gene2,S3=bp,S4=bp2,S5=clipOrder,S6=clipOrder2) in usedpair:
                                continue
                            if '{S1}_{S2}_{S3}_{S4}_{S5}_{S6}'.format(S1=gene2,S2=gene,S3=bp2,S4=bp,S5=clipOrder2,S6=clipOrder) in usedpair:
                                continue
                            if clipOrder == clipOrder2 and same[1] == 'F': #'F' = forward read
                                continue
                            if clipOrder != clipOrder2 and same[1] == "RC": #'RC' = reverse complement read 
                                continue
                            usedpair['{S1}_{S2}_{S3}_{S4}_{S5}_{S6}'.format(S1=gene2,S2=gene,S3=bp2,S4=bp,S5=clipOrder2,S6=clipOrder)] = 1 #store unique gene pair

                            #assign strand orientation to fusion
                            orientation1 = "2-"
                            if clipOrder == "NC":
                                orientation1 = "1+"
                            orientation2 = "2-"
                            if clipOrder == "CN":
                                orientation1 = "1+"
                            if clipOrder != clipOrder2 and clipOrder == "CN":
                                orientation1 = "2+"
                                orientation2 = "1+"
                            if clipOrder != clipOrder2 and clipOrder == "NC":
                                orientation2 = "2+"
                                orientation1 = "1+"

                            #do bp correction by comparison to reference=================
                            
                            offsetbp1 = 0
                            offsetbp2 = same[0]

                            if same[0] != 0:
                                tmp = bp_correction(screads[gene2][bp2_][1],same,offsetbp2,bp2_,chr2,OUTPUTDIR,bam,twobit)
                                offsetbp1 = tmp[0]
                                offsetbp2 = tmp[1]
                            #=============================================================

                            #Find and store insert (non-templated or microhomologous) sequence (if any)
                            if offsetbp2 != 0:
                                if offsetbp2 > 0 and clipOrder2 == "NC":
                                    fusioncons = read2[:offsetbp2]
                                elif offsetbp2 < 0 and clipOrder == "NC":
                                    fusioncons = doRevComp(read1b[:(offsetbp2 * -1)])
                            
                            #return and print adjusted fragment 2 breakpoint
                            fusionSeqFa = getFusionSeq(clipOrder, clipOrder2, bp_, bp2_, bp, bp2, offsetbp2, buffer, offsetbp1,twobit,OUTPUTDIR,bam)
                            vals = fusionSeqFa.split("\n")
                            fusionSeq = vals[1]
                            part1 = fusionSeq[:buffer]
                            part2 = fusionSeq[buffer:buffer+buffer]
                            if len(fusioncons) > 0:
                                nt_len = len(fusioncons)
                                if '1' in orientation1 :
                                    fusioncons = doRevComp(fusioncons)
                                    if part2[:len(fusioncons)] == fusioncons:
                                        fusioncons = ""
                                    else:
                                        part2 = fusioncons + part2[len(fusioncons):]
                                        if '+' in orientation2:
                                            bp2_ += nt_len
                                        else:
                                            bp2_ -= nt_len
                                else:
                                    if part1[len(part1)-len(fusioncons):len(part1)] == fusioncons:
                                        fusioncons = ''
                                    else:
                                        part1 = part1[:len(part1)-len(fusioncons)] + fusioncons + part1[:len(part1)] 
                                        if '+' in orientation1:
                                            bp2_ -= nt_len
                                        else:
                                            bp2_ += nt_len
                            
                            #print "$fusioncons\t$bp_\t$bp2_\n"; 
                            if fusioncons == "":
                                fusioncons = '-'
                                fusioncons_bp = 0
                            else:
                                fusioncons = "[{S1}]".format(S1=fusioncons)
                            nonref_insert = ''
                            fusionSeqFa = vals[0] + '\n' + part1 + "" + part2
                            #=============================================================

                            data = gene + '\t' + gene2 + '\t' + chr_ + ':' + str(bp_ + offsetbp1) + '\t' + chr2 + ':' + str(bp2_ + offsetbp2)
                            data += '\t{S1}\t{S2}\t{S3}\t{S4} {S5}\t{S6}\t{S7}\t'.format(S1=count,S2=count2,S3=same[0],S4=orientation1,S5=orientation2,S6=clipOrder,S7=clipOrder2)
                            #do BLAST search
                            datatmp = doBLAST(fusionSeqFa,gene,gene2,buffer,(MINBLASTFRAC * READLEN),bam,OUTPUTDIR,BLASTTHREADS,minsimilarity)
                            data += datatmp

                            #print close-up of breakpoint
                            vals = fusionSeqFa.split('\n')
                            fusionSeq = vals[1]
                            part1 = fusionSeq[buffer-pad:buffer]
                            part2 = fusionSeq[buffer:buffer+pad]

                            #add non-templated sequence into estimated fusion sequence
                            data = data + part1 + " <> " + part2 + "{S1}".format(S1=fusioncons)
                            fusioncons = ''
                            #print "$data\n";
                            if datatmp == '':
                                continue
                            vars_ = datatmp.split('\t')
                            bp_depth = vars_[0]
                            if float(bp_depth) < MINSPANNINGREADS: #don't save putative translocations with too few reads spanning breakpoint
                                continue
                            tmp = [data , fusionSeqFa]
                            code = chr_ + ":" + str(bp_ + offsetbp1) + ' ' + chr2 + ":" + str(bp2_ + offsetbp2)
                            if bp_depth not in bestfusions:
                                bestfusions[bp_depth] = {}
                            if code not in bestfusions[bp_depth]:
                                bestfusions[bp_depth][code] = tmp
                            else:
                                tmp2 = bestfusions[bp_depth][code]

    #if(($progress_itor == $total_comp || $progress_itor == 0) && $VERBOSE == 1) {print "\n";}
    #print header
    if len(bestfusions) > 0:
        output3.write("Est_Type\tRegion1\tRegion2\tBreak1\tBreak2\tBreak_support1\tBreak_support2\tBreak_offset\tOrientation\tOrder1\tOrder2\tBreak_depth\tProper_pair_support\tUnmapped_support\tImproper_pair_support\tPaired_end_depth\tTotal_depth\tFusion_seq\tNon-templated_seq\n")
        if VERBOSE == 1:
            print("\n-------------------------------------------------\n\nFACTERA results:\nEst_Type\tRegion1\tRegion2\tBreak1\tBreak2\tBreak_support1\tBreak_support2\tBreak_offset\tOrientation\tOrder1\tOrder2\tBreak_depth\tProper_pair_support\tUnmapped_support\tImproper_pair_support\tPaired_end_depth\tTotal_depth\tFusion_seq\tNon-templated_seq\n")
    #===========================================================================================================       
    #print final fusion predictions, ordered by decreasing breakpoint depth (remove redundancies with 20bp window)
    bestfusions_ = {} #store best fusion events and use to remove redundancies (+/- 20bp padding)

    normdepth = getNormDepth(bestfusions, buffer, targets, OUTPUTDIR, bam) #total depth of genomic regions flanking candidate fusion(s) 

    transitor = 0     #count translocations/fusions
    for bp_depth in bestfusions:
        for bps in bestfusions[bp_depth]:
            orig_bps = bps 
            tmp = bestfusions[bp_depth][bps]
            data = tmp[0]
            tokens = bps.split(' ')                   
            bp1 = int(tokens[0].split(':')[1])
            bp2 = int(tokens[1].split(':')[1])
            chr1 = tokens[0].split(':')[0]
            chr2 = tokens[1].split(':')[0]
            var = data.split('\t')
            orient1 = var[7]

            #Estimate rearrangement type
            first = orient1[:1]
            second = orient1[3:4]
            first_polarity = orient1[1:2]
            sec_polarity = orient1[4:5]
            SR_type = ' - '
            if chr1 != chr2:
                SR_type = 'TRA'
            elif first_polarity != sec_polarity:
                SR_type = 'INV'
            elif (bp1 < bp2 and first == 1) or (bp1 > bp2 and first == 2):
                SR_type = "DEL"
            #elsif(abs(1+$bp1-$bp2)>1000000){$SR_type = "TRA";}
        
            #Switch order for instances where region 2 is first (5' of fusion)
            if first == 2:
                var[0],var[1] = var[1],var[0]
                var[2],var[3] = var[3],var[2]
                var[4],var[5] = var[5],var[4]
                var[8],var[9] = var[9],var[8]
                var[7] = var[7].replace('1','0')
                orient1 = var[7].replace('2','0')
                bp1,bp2 = bp2,bp1
                chr1,chr2 = chr2,chr1
                bps = chr1 + ":" + str(bp1) + ' ' + chr2 + ':' +str(bp2)
            for bp_depth_ in bestfusions:
                for bps_ in bestfusions[bp_depth]:
                    tokens2 = bps_.split(' ')
                    bp1_ = int(tokens2[0].split(':')[1])
                    bp2_ = int(tokens2[1].split(':')[1])
                    chr1_ = tokens2[0].split(':')[0]
                    chr2_ = tokens2[1].split(':')[0]
                    tmp_ = bestfusions[bp_depth][bps_]
                    data_ = tmp_[0]
                    var_ = data_.split('\t')
                    orient2 = var_[7]
                    first = orient2[:1]
                    if first == 2:
                        var_[7] = var_[7].replace('1','0')
                        orient2 = var_[7].replace('2','0')
                        bp1_ , bp2_ = bp2_ , bp1_
                        chr1_ , chr2_ = chr2_ , chr1_
                        bps_ = chr1_ + ":" + bp1_ + ' ' + chr2_ + ':' + bp2_
                    if bps == bps_:
                        continue
                    if chr1 != chr1_ or chr2 != chr2_:
                        continue
                    if orient1 != orient2:
                        continue
                    
                    if bp1_ <= bp1 + 20 and bp1_ >= bp1 - 20 and bp2_ <= bp2 + 20 and bp2_ >= bp2 - 20:
                        bestfusions_[bps_] = bps_
            
            if bps not in bestfusions_:
                bestfusions_[bps] = bps
                seq = var[-1]
                tokens = seq.split(' ')
                seq = '\t'.join(tokens)
                var[-1] = str(normdepth[orig_bps]) #retrieve genomic depth
                data = '\t'.join(var)
                data = SR_type + '\t' + data + '\t' + seq + '\n'
                output3.write(data)
                if VERBOSE == 1:
                    print(data+'\n')
                fusionSeq = tmp[1] + '\n'
                output4.write(fusionSeq)
                t1 = var[2].split(":")
                t2 = var[3].split(":")
                wr5 = '\t'.join([t1[0],t1[1],t1[1],'']) + '_'.join([var[0],var[1],var[2],var[3]]) + '\n'
                wr5_ = '\t'.join([t2[0],t2[1],t2[1],'']) + '_'.join([var[0],var[1],var[2],var[3]]) + '\n'
                output5.write(wr5)
                output5.write(wr5_)
                transitor += 1
    
    if VERBOSE == 1:
        if transitor > 0:
            ins = 'fusion'
            if transitor > 1:
                ins = 'fusions'
            print(str(transitor) + ' ' + ins + 'found!\n')
        else:
            print('\nNo fusions found\n')
    output3.close()
    output4.close()
    output5.close()

    now = time.time() - now
    if TIME == 1:
        print(now)





if __name__ == "__main__":
    main()





