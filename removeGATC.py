from Bio.Seq import Seq
from sys import argv
from sys import exit
from Bio import SeqIO
from copy import deepcopy
from time import strftime

if len(argv) == 1:
    exit('input file is required')
elif len(argv) == 2:
    input_file = argv[1]
elif len(argv) == 3:
    input_file = argv[1]
    output_file = argv[2]
elif len(argv) > 3:
    input_file = argv[1]
    output_file = argv[2]
    print("Warning: more than one argument given, only using first two")


record_list = [record for record in SeqIO.parse(input_file, 'fasta')]
print("Processing %i sequences \n\n" % (len(record_list)))


new_records = deepcopy(record_list)
for tfn, record in enumerate(record_list):
    print("Sequence %i of %i: %s" % (tfn+1, len(record_list), record.name))
    # load the CDS sequence of the TF. Should start with ATG (or at least be in reading frame 1)
    tfseq = record.seq

    # convert sequence to list of characters
    newtfseq = [x for x in tfseq]

    # get the starting position and coding frame for each GATC site
    pos = [x+1 for x, y in enumerate(tfseq) if tfseq[x:(x+4)] == 'GATC']
    frame = [round(x/3 % 1, 1) for x in pos]
    print("  is %i bp long and has %i GATC sites " % (len(tfseq), len(pos)))

    # create empty lists for codons around GATC sites
    oldcodons = []
    newcodons = []

    # go through sequence and replace GATC sites depending on reading frame (tried to use common codons)
    for x, y in enumerate(pos):
        # if GATC is in frame 3, the Ile codon ATC is replaced with ATT
        if frame[x] == 0:
            oldcodons.append(tfseq[y-3:(y+3)])
            newcodons.append(tfseq[y-3:y+2]+'T')
            newtfseq[y-3:(y+3)] = newtfseq[y-3:y+2] + ['T']
        # if the GATC is in frame 1, there are three possibilities
        # not sure why I didn't just change GATC to GACC, maybe had to do with potential to create new GATC sites?
        elif frame[x] == 0.3:
            oldcodons.append(tfseq[y-1:y+5])
            # if nucleotide following GATC is a T, the second codon must be Leu, replace CTN with TTA
            if tfseq[y+3] == 'T':
                newcodons.append('GATTTA')
                newtfseq[y-1:y+5] = [z for z in 'GATTTA']
            # if nucleotide following GATC is a G, the second codon must be Arg, replace CGN with AGA
            elif tfseq[y+3] == 'G':
                newcodons.append('GATAGA')
                newtfseq[y-1:y+5] = [z for z in 'GATAGA']
            # if neither of the above, replace GATC with GAC
            else:
                newcodons.append('GAC'+tfseq[y+2:y+5])
                newtfseq[y-1:y+5] = ['G'] + ['A'] + ['C'] + newtfseq[y+2:y+5]
        # if GATC is in frame 2, the Ser codon TCC is replaced by AGT
        elif frame[x] == 0.7:
            oldcodons.append(tfseq[y-2:y+4])
            newcodons.append(tfseq[y-2:y+1]+'AGT')
            newtfseq[y-2:y+4] = newtfseq[y-2:y+1] + ['A'] + ['G'] + ['T']

    if len(pos) > 0:
        print('  GATC sites removed')

    # get the starting postion for each of the BsaI and BsmBI sites
    posBsaIF = [x+1 for x, y in enumerate(tfseq) if tfseq[x:(x+6)] == 'GGTCTC']
    posBsaIR = [x+1 for x, y in enumerate(tfseq) if tfseq[x:(x+6)] == 'GAGACC']
    frameBsaIF = [round(x/3 % 1, 1) for x in posBsaIF]
    frameBsaIR = [round(x/3 % 1, 1) for x in posBsaIR]
    posBsmBIF = [x+1 for x, y in enumerate(tfseq) if tfseq[x:(x+6)] == 'CGTCTC']
    posBsmBIR = [x+1 for x, y in enumerate(tfseq) if tfseq[x:(x+6)] == 'GAGACG']
    frameBsmBIF = [round(x/3 % 1, 1) for x in posBsmBIF]
    frameBsmBIR = [round(x/3 % 1, 1) for x in posBsmBIR]
    print("  it has %i BsaI and %i BsmBI sites " % (len(posBsaIF + posBsaIR), len(posBsmBIF + posBsmBIR)))
    if len(posBsaIF + posBsaIR) == 0:
        print("  Clone into DamID vector with BsaI")
    elif len(posBsmBIF + posBsmBIR) == 0:
        print("  Clone into DamID vector with BsmBI")
    else:
        if len(posBsaIF + posBsaIR) <= len(posBsmBIF + posBsmBIR):
            # go through sequence and replace Fwd BsaI sites (GGTCTC) depending on reading frame
            for x, y in enumerate(posBsaIF):
                # if BsaI is in frame 1, the Gly codon GGT is replaced with GGA
                if frame[x] == 0.3:
                    newtfseq[y + 1] = 'A'
                # if the BsaI site is in frame 2, replace Val codon GTC with GTG
                elif frame[x] == 0.7:
                    newtfseq[y + 2] = 'G'
                # if the BsaI site is in frame 3, replace Ser codon TCT with TCA
                elif frame[x] == 0:
                    newtfseq[y + 3] = 'A'
            # replace Rev BsaI sites (GAGACC) depending on frame
            for x, y in enumerate(posBsaIR):
                # if BsaI is in frame 1, the Glu codon GAG is replaced with GAA
                if frame[x] == 0.3:
                    newtfseq[y + 1] = 'A'
                # if the BsaI site is in frame 2, replace Arg codon AGA with AGG
                elif frame[x] == 0.7:
                    newtfseq[y + 2] = 'G'
                # if the BsaI site is in frame 3, can't replace Asp because it makes GATC
                # therefore need to change codon before Asp which is NGA. N can't be T because that is STOP
                # if N is C or A, it is Arg and can change to AGG
                # if N is G, it is Gly and change to GGT
                elif frame[x] == 0:
                    if newtfseq[y - 2] in ['C', 'A']:
                        newtfseq[y-2:y+1] = ['A', 'G', 'G']
                    else:
                        newtfseq[y] = 'T'
            print("  Removed %i BsaI sites\n  Clone into DamID vector with BsaI" % len(posBsaIF + posBsaIR))
            if [x + 1 for x, y in enumerate(newtfseq) if newtfseq[x:(x+4)] == ['G', 'A', 'T', 'C']] == [] and \
               [x + 1 for x, y in enumerate(newtfseq) if newtfseq[x:(x + 6)] in [['G', 'G', 'T', 'C', 'T', 'C'],
                                                                                 ['G', 'A', 'G', 'A', 'C', 'C']]] == []:
                print("  Confirmed no GATC or BsaI sites remaining")
            else:
                print("  WARNING: GATC or BsaI sites remain")
        else:
            # go through sequence and replace Fwd BsmBI sites (CGTCTC) depending on reading frame
            for x, y in enumerate(posBsmBIF):
                # if BsmBI is in frame 1, the Leu codon CTC is replaced with CTT
                if frame[x] == 0.3:
                    newtfseq[y + 4] = 'T'
                # if the BsmBI site is in frame 2, replace Val codon GTC with GTG
                elif frame[x] == 0.7:
                    newtfseq[y + 2] = 'G'
                # if the BsmBI site is in frame 3, replace Ser codon TCT with TCA
                elif frame[x] == 0:
                    newtfseq[y + 3] = 'A'
            # replace Rev BsmBI sites (GAGACG) depending on frame
            for x, y in enumerate(posBsmBIR):
                # if BsaI is in frame 1, the Thr codon ACG is replaced with ACA
                if frame[x] == 0.3:
                    newtfseq[y + 4] = 'A'
                # if the BsmBI site is in frame 2, replace Arg codon AGA with AGG
                elif frame[x] == 0.7:
                    newtfseq[y + 2] = 'G'
                # if the BsmBI site is in frame 3, replace Asp GAC with GAT
                elif frame[x] == 0:
                    newtfseq[y + 3] = 'T'
            print("  Removed %i BsmBI sites\n  Clone into DamID vector with BsmBI" % len(posBsaIF + posBsaIR))
            if [x + 1 for x, y in enumerate(newtfseq) if newtfseq[x:(x + 4)] == ['G', 'A', 'T', 'C']] == [] and \
               [x + 1 for x, y in enumerate(newtfseq) if newtfseq[x:(x + 6)] in [['C', 'G', 'T', 'C', 'T', 'C'],
                                                                                 ['G', 'A', 'G', 'A', 'C', 'G']]] == []:
                print("  Confirmed no GATC or BsmBI sites remaining")
            else:
                print("  WARNING: GATC or BsmBI sites remain")

    # convert the list of characters into a string.
    newtfseqstr = Seq(''.join(newtfseq))

    # confirm new translated sequence is same as old sequence
    if tfseq.translate() == newtfseqstr.translate():
        print("  new sequence translation matches old: Confirmed\n")
    else:
        print("  somebody (Matt) done fucked up")

    new_records[tfn].seq = newtfseqstr
# write the output file
if 'output_file' in locals():
    SeqIO.write(new_records, output_file, 'fasta')
else:
    timestr = strftime("%Y%m%d_%H%M%S")
    SeqIO.write(new_records, timestr + "_RemoveGATC_out.fa", 'fasta')
    
    
