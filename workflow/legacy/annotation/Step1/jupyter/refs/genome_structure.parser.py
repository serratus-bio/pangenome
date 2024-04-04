from Bio import SearchIO
import sys, os

filename='genome_structure/SRR9967737.transeq.domtbl' # incomplete rdrp
filename="genome_structure/SRR11939968.transeq.domtbl" # complete genome

"""
another example of annotations where everything is there but VADR either marks it as PASS or FAIL.
pass:
SRR11826966.darth/SRR11826966/SRR11826966.vadr.pass.tbl
fail:
SRR11953795.darth/SRR11953795/SRR11953795.vadr.ftr
"""

if len(sys.argv) > 1:
    filename = sys.argv[1]

def inv_transeq_coords(ctg_name, ctg_len, start, end):
    frame = int(ctg_name[-1])
    if frame <= 3:
        return (start*3, end*3, "+") # rough approximation of coordinates (+/- 3bp precise?)
    else:
        return (ctg_len*3-end*3, ctg_len*3-start*3, "-") # rough approximation here too

res = []
for qresult in SearchIO.parse(filename,'hmmsearch3-domtab'):
    for hit in qresult.hits:
        hmm_len = qresult.seq_len
        hmm_id = qresult.id
        contig_len = hit.seq_len
        #print("newhit",hit.query_id,qresult.seq_len,qresult.id)
        for item in hit.hsps:
            hit_name      = item.query_id
            contig_name   = item.hit_id
            transeq_start = item.hit_start
            transeq_end   = item.hit_end
            contig_start, contig_end, contig_strand = inv_transeq_coords(contig_name, contig_len, transeq_start, transeq_end)
            #print(hit_name,contig_name,contig_start, contig_end)
            complete = item.hit_span >= 0.75 * hmm_len 
            # FIXME: unsure that will properly take care of incomplete genes. Verify that result agrees with Robert's suggestion: " You have to check for an incomplete RdRP alignment, you can't see that without counting gaps in the PFAM a2m. I use 25% gaps as the cutoff.""
            res += [(contig_start, contig_end, contig_strand, hit_name, contig_name, complete)]

# have ".nido" annotations
pfam_hmm = ".nido-pfam" if "nido" in filename else "" # "" is shorthand for "sars-cov-2"

accession = '.'.join(os.path.basename(filename).split('.')[:2])
print("(\"%s%s\", (\\" % (accession, pfam_hmm))
for e in sorted(res):
    print(str(e)+",")
print(")),")

