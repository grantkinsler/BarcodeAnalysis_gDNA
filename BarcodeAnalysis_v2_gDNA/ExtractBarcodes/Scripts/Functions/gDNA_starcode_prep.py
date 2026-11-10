# This function takes in Fastq files taken from gDNA barcodes and reverse complements the sequences + trims them
# so that starcode can run properly

from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
from rapidfuzz import fuzz
from time import sleep

def gDNA_starcode_prep(FastqfoldergDNA,sc_in,GSAMP,bclen,filt_haveStart_gDNA,strtseq,startseqMatch,sc_mm,tag_seqs=None,tagMatch=75,include_tag_in_clustering=False,minQual=None,exclude_tag_seqs=None,reverse_complement=True,assume_tag_present=False):
    # tag_seqs: dict of {tag_name: forward-orientation tag sequence immediately
    # following the barcode in the construct}. Reads are reverse-complemented
    # before strtseq is located, so in that search orientation the tag actually
    # sits BETWEEN strtseq and the barcode.
    #
    # include_tag_in_clustering=False (default): skip the tag_len-nt window before
    # slicing out the true bclen-nt barcode, so the tag never reaches starcode /
    # the Hamming-distance clustering (tag is still called and reported separately).
    #
    # include_tag_in_clustering=True: don't skip anything -- always carry the
    # tag_len-nt window immediately before the barcode through into what's fed to
    # starcode, so clustering sees (tag_or_whatever_is_there + barcode) as one unit.
    # Applied uniformly regardless of whether that window fuzzy-matches a known tag,
    # so every read contributes a fixed-length window and reads genuinely lacking a
    # tag (e.g. ladder/spike-in) just carry whatever sequence sits there instead.
    #
    # minQual (mean Phred quality, or None to disable): applied separately to the tag
    # window and the barcode window. Low tag-window quality doesn't change bc_offset
    # (a tag that fuzzy-matched is still structurally skipped so the barcode stays in
    # register) -- it only downgrades the tag call itself to "None", since we're no
    # longer confident WHICH tag it is. Low barcode-window quality drops the read
    # entirely (same as failing the anchor match), since an unreliable barcode call
    # would otherwise get fed straight into starcode clustering.
    #
    # exclude_tag_seqs: dict of {tag_name: forward-orientation sequence}, same shape
    # as tag_seqs, for tags that are known NOT to belong to this experiment (e.g.
    # Tag3-6 confirmed 2026-07-31 as contamination from an unrelated multi-tag library
    # sharing sequencing capacity -- see eGRK1059_lod_and_hamming_clustering memory).
    # The tag window is matched against tag_seqs and exclude_tag_seqs together, taking
    # the single best-scoring candidate across both; if that best match is an excluded
    # tag, the read is dropped entirely (never reaches starcode), same as failing the
    # anchor match -- not merely relabeled "None" the way a genuinely tag-less
    # (e.g. ladder) read would be.
    #
    # reverse_complement: bool (default True), or {sample_name: bool} for a per-sample
    # override, with True as the default for any sample not listed. True -- the
    # historical/only behavior before this parameter existed, and confirmed to also
    # be correct for an R2 read from a library with read orientation deliberately
    # swapped relative to the normal R1 setup -- reverse-complements the raw read
    # before searching for strtseq. Both R1 and R2 reads get the same treatment here;
    # the fastq file searched for a sample is whichever of *_R1*.fastq / *_R2*.fastq
    # exists in that sample's folder (see barcode_read in run_barcode_pipeline.py),
    # but which one it physically is doesn't change this default. Set False only if a
    # specific construct/read genuinely reads directly into strtseq -> tag -> barcode
    # without needing a flip.
    #
    # assume_tag_present=False (default): bc_offset (whether the barcode slice skips
    # the tag_len-nt tag window) is only set when a read's tag window fuzzy-matches a
    # declared tag at >=tagMatch. A read whose tag window has enough sequencing error
    # to miss that threshold falls back to bc_offset=0 -- correct for a construct where
    # some reads genuinely have NO tag at that position (e.g. a tag-less ladder/spike-in
    # mixed in with tagged samples), but WRONG for a single-tag library where the tag is
    # always physically present: those misses silently shift the reported "barcode" to
    # [[the noisy tag][first (bclen-tag_len) nt of the real barcode]] instead of the
    # real barcode, which starcode then can't merge back with the correctly-offset
    # majority (identical real barcode, totally different noisy prefix each time) --
    # this fabricates one spurious extra "barcode" per distinct tag-read-error pattern
    # and can inflate the apparent distinct-barcode count by an order of magnitude.
    # assume_tag_present=True: skip the tag_len-nt window unconditionally whenever
    # tag_len>0, regardless of fuzzy-match score -- use this when every read in this
    # sample is known to carry a tag at that fixed position (e.g. tag_seqs has exactly
    # one entry and there's no genuinely tag-less read population). The fuzzy-match
    # result still determines the reported per-read tag LABEL (for QC visibility into
    # how often the tag region itself is noisy), it just no longer gates bc_offset.
    if tag_seqs is None:
        tag_seqs = {}
    if exclude_tag_seqs is None:
        exclude_tag_seqs = {}
    tag_names = list(tag_seqs.keys())
    tag_seqs_rc = {name: str(Seq(seq).reverse_complement()) for name, seq in tag_seqs.items()}
    exclude_tag_names = list(exclude_tag_seqs.keys())
    exclude_tag_seqs_rc = {name: str(Seq(seq).reverse_complement()) for name, seq in exclude_tag_seqs.items()}
    tag_len = len(next(iter(tag_seqs.values()))) if tag_seqs else 0
    sc_bclen = bclen + tag_len if include_tag_in_clustering else bclen
   #get all gDNA barcode-read fastq file paths. Matches *_R1* or *_R2* -- a sample's
   #barcode-bearing read is normally R1 for this construct, but may be staged as R2
   #for a sample whose read orientation was swapped (see reverse_complement above).
   #Each sample is expected to have exactly one such file in its own folder; the
   #variable names below still say "R1" for historical reasons but cover both.
    all_R1_gDNA_unfilt = glob(FastqfoldergDNA + "/**/*_R1*.fastq", recursive = True)
    all_R1_gDNA_unfilt += glob(FastqfoldergDNA + "/**/*_R2*.fastq", recursive = True)
    all_R1_gDNA_unfilt.sort()
    
    # Remove any Read1 fastq files that you dont care about
    all_R1_gDNA_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_R1_gDNA_unfilt):
                all_R1_gDNA_temp.append([s for s in all_R1_gDNA_unfilt if str(smp) in s])

    all_R1_gDNA = [item for sublist in all_R1_gDNA_temp for item in sublist]

    #define samples Read1
    samples_R1_gDNA = []
    for paths in all_R1_gDNA:

        samples_R1_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_R1_gDNA = list(set(samples_R1_gDNA))
    samples_R1_gDNA.sort()
    
    
    #loop through all the Read 1 fastq files from gDNA (which contain barcode sequences) and write the textfile for starcode inputs
    readsKeep_gDNA = []
    counter = 0
    file_endpoints_gDNA = []
    for sample in samples_R1_gDNA:
        file_endpoints_gDNA.append([])
        file_endpoints_gDNA[counter].insert(0,0)

        s_fastq = []
        for path in all_R1_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)

        # resolve this sample's reverse_complement setting (bool, or {sample: bool}
        # with True as the default for any sample not listed)
        if isinstance(reverse_complement, dict):
            sample_rc = reverse_complement.get(sample, True)
        else:
            sample_rc = reverse_complement

        # get JUST the sequences in all the fastqs contains all sequences for a sample
        seqs = []
        quals = []  # per-base Phred quality, same reversal as seqs (no complementing -- quality has no complement)
        for fsmp in s_fastq:
            for record in SeqIO.parse(fsmp, "fastq"):
                if sample_rc:
                    seqs.append(str(record.seq.reverse_complement()))
                    quals.append(record.letter_annotations["phred_quality"][::-1])
                else:
                    seqs.append(str(record.seq))
                    quals.append(record.letter_annotations["phred_quality"])
            file_endpoints_gDNA[counter].append(len(seqs))

        #determine which seqeunces are actual barcodes by determining which start with the constant sequence before the barcode

        # variables: the constant sequence before the barcode, percent match to call it a barcode (here it is 70)

        # since we stagger priming sight, here we determine what the correct stagger is for a given fastq file.
        # Reads whose anchor is at this file's most common ("mode") stagger position dominate, but some
        # reads (e.g. a different insert, like a spike-in/ladder construct with no tag and a different
        # length) have a genuinely different, consistent stagger of their own. Rather than forcing every
        # read to the single global mode (which silently drops any read staggered differently), find each
        # read's own exact anchor position first and only fall back to the file's mode when a read has no
        # exact match (e.g. because of a substitution error within the anchor itself).
        strt_ind = []
        strt_ind_per_read = []
        for lines in seqs:
            c_ind = len(lines.split(strtseq)[0])
            if c_ind < len(lines):
                strt_ind.append(c_ind)
                strt_ind_per_read.append(c_ind)
            else:
                strt_ind_per_read.append(None)

        bc_strt = statistics.mode(strt_ind)  # fallback for reads with no exact anchor match, and for the (inconsequential) quality-line slice below

        # modify sequences so that those that have a start sequence (where there can be erros) get replaces with a pefect start sequence, and get rid of stagger so that all sequence start with the perfect start sequences
        modseq1 = [] # For writing just the good sequences
        modseq2 = [] # For writing good and bad sequnces
        tag_calls = [] # Tag call for each read written to sc_input, same order/length as modseq1/sc_input
        readsKeep_gDNA.append([])
        for read_idx, i in enumerate(seqs):
            this_bc_strt = strt_ind_per_read[read_idx]
            if this_bc_strt is None:
                this_bc_strt = bc_strt
            strt = i[this_bc_strt:(len(strtseq)+this_bc_strt)]
            pctmatch = (fuzz.ratio(strtseq,strt))

            if pctmatch >= startseqMatch:
                trim = i[this_bc_strt + len(strtseq) :]
                trim_q = quals[read_idx][this_bc_strt + len(strtseq) :]

                # Check for a tag at the fixed offset immediately after strtseq
                # (i.e. before the barcode, in this reverse-complemented search
                # orientation). Fuzzy-matched the same way strtseq is, so a
                # sequencing error in the tag doesn't just get treated as "no tag".
                called_tag = "None"
                bc_offset = 0
                dropped_as_excluded_tag = False
                if tag_len > 0:
                    candidate = trim[0:tag_len]
                    best_name, best_score, best_is_excluded = None, -1, False
                    for name in tag_names:
                        score = fuzz.ratio(tag_seqs_rc[name], candidate)
                        if score > best_score:
                            best_name, best_score, best_is_excluded = name, score, False
                    for name in exclude_tag_names:
                        score = fuzz.ratio(exclude_tag_seqs_rc[name], candidate)
                        if score > best_score:
                            best_name, best_score, best_is_excluded = name, score, True
                    tag_matched = best_score >= tagMatch
                    if tag_matched and best_is_excluded:
                        dropped_as_excluded_tag = True
                    else:
                        # assume_tag_present: skip the tag window unconditionally (the
                        # tag is known to always be physically present at this fixed
                        # offset for this sample), independent of fuzzy-match confidence.
                        # Otherwise (default): only skip when the tag confidently matched
                        # -- structurally a tag is present (sequence matched well enough),
                        # so the barcode window starts after it regardless of quality;
                        # quality only gates whether we trust WHICH tag it is.
                        if assume_tag_present or tag_matched:
                            bc_offset = tag_len
                        if tag_matched:
                            candidate_q = trim_q[0:tag_len]
                            if minQual is None or statistics.mean(candidate_q) >= minQual:
                                called_tag = best_name

                if dropped_as_excluded_tag:
                    readsKeep_gDNA[counter].append(False)
                    modseq2.append(i)
                    continue

                # tag calling above is unaffected either way; this only controls whether
                # the detected/skippable window actually gets dropped from what's clustered
                effective_offset = 0 if include_tag_in_clustering else bc_offset
                bc_trim = trim[effective_offset:]

                # Drop the read entirely if the true bclen-nt barcode window (always
                # starting at bc_offset, regardless of include_tag_in_clustering) is
                # too low quality to trust feeding into starcode clustering.
                bc_window_q = trim_q[bc_offset:bc_offset + bclen]
                if minQual is not None and (len(bc_window_q) < bclen or statistics.mean(bc_window_q) < minQual):
                    readsKeep_gDNA[counter].append(False)
                    modseq2.append(i)
                    continue

                # NB: strtseq is deliberately NOT included here. It's identical on every read
                # (this literal constant is prepended regardless of the read's actual anchor
                # bases), so it can never affect which reads cluster together in starcode --
                # prepending an identical string to two sequences never changes their edit
                # distance. Including it anyway just makes every starcode comparison and every
                # input/output file bigger for no benefit, so we drop it before writing sc_input.
                modseq1.append(bc_trim)
                modseq2.append(bc_trim)
                tag_calls.append(called_tag)
                readsKeep_gDNA[counter].append(True)
            else :
                readsKeep_gDNA[counter].append(False)
                modseq2.append(i)

        #trim barcodes
        # variable: how long do you want your barcode (sc_bclen = bclen, or bclen+tag_len
        # when include_tag_in_clustering so the tag window rides along into starcode)
        sc_input = []
        for i in modseq1:
            sc_input.append(i[0:sc_bclen])

        #write files with these edited barcodes ( these are used as the input into starcode)
        f = open(sc_in + "/" "sc_input" + "_" + sample +".txt","w")
        f.write('\n'.join(sc_input))
        f.close()

        # write the tag call for each read, in the same line order as sc_input,
        # so it can be rejoined after starcode clustering (which never sees the tag)
        f = open(sc_in + "/" "tag_calls" + "_" + sample + ".txt","w")
        f.write('\n'.join(tag_calls))
        f.close()
        sleep(20)
    
    
        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_gDNA + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_gDNA[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_gDNA[counter][lane_counter+1]))
                for j in range(file_endpoints_gDNA[counter][lane_counter],file_endpoints_gDNA[counter][lane_counter+1]):
                    if readsKeep_gDNA[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(modseq2[j][0:sc_bclen]+'\n') # pulls the trimmed sequence (strtseq-free)
                        f.write(lines[(4*j)+2])
                        # NB: this quality slice doesn't shift for a per-read tag offset; harmless because
                        # Step3.py overwrites this read's sequence+quality with the starcode canonical
                        # barcode before anything downstream reads it. Length MUST match the sequence
                        # line above (sc_bclen) or gDNA_after_starcode's fastq parser errors out on a
                        # sequence/quality length mismatch. Start is shifted past strtseq's own quality
                        # scores since the sequence line above no longer includes strtseq either.
                        #
                        # KNOWN GAP: this [::-1] unconditionally assumes sample_rc=True (matches the
                        # sequence being written on the strtseq-trimmed REVERSE-COMPLEMENT, consistent
                        # with every call site below). It is NOT fixed up for a reverse_complement=False
                        # sample -- this CellRanger-output path (Step3.py) was not exercised for the
                        # orientation-swap feature. run_barcode_pipeline.py never reads this output
                        # (it only reads sc_input_<sample>.txt / tag_calls_<sample>.txt, which ARE
                        # correct for reverse_complement=False), so this only matters if you use
                        # Step3.py's CellRanger feature-barcoding output directly with a swapped sample.
                        f.write(lines[(4*j)+3][::-1][bc_strt+len(strtseq)+1:(sc_bclen+len(strtseq)+bc_strt)+1]+'\n') # Also reverses the phred score sequence
        
            lane_counter = lane_counter + 1

    
        counter = counter+1
        
    
    
    #get all gDNA Index1 fastq file paths
    all_I1_gDNA_unfilt = glob(FastqfoldergDNA + "/**/*_I1*.fastq", recursive = True)
    all_I1_gDNA_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I1_gDNA_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I1_gDNA_unfilt):
                all_I1_gDNA_temp.append([s for s in all_I1_gDNA_unfilt if str(smp) in s])

    all_I1_gDNA = [item for sublist in all_I1_gDNA_temp for item in sublist]

    #define samples Read1
    samples_I1_gDNA = []
    for paths in all_I1_gDNA:

        samples_I1_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_I1_gDNA = list(set(samples_I1_gDNA))
    samples_I1_gDNA.sort()
    
    # Loop through the corresponding index1 files and remove the bad reads
    counter = 0
    for sample in samples_I1_gDNA:
        s_fastq = []
        for path in all_I1_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_gDNA + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_gDNA[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_gDNA[counter][lane_counter+1]))
                for j in range(file_endpoints_gDNA[counter][lane_counter],file_endpoints_gDNA[counter][lane_counter+1]):
                    if readsKeep_gDNA[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3])
        
            lane_counter = lane_counter + 1

    
        counter = counter+1
        
    #get all gDNA Index2 fastq file paths
    all_I2_gDNA_unfilt = glob(FastqfoldergDNA + "/**/*_I2*.fastq", recursive = True)
    all_I2_gDNA_unfilt.sort()

    # Remove any Read1 fastq files that you dont care about
    all_I2_gDNA_temp = []
    for grp in GSAMP:
        for smp in grp:
            if str(smp) in str(all_I2_gDNA_unfilt):
                all_I2_gDNA_temp.append([s for s in all_I2_gDNA_unfilt if str(smp) in s])

    all_I2_gDNA = [item for sublist in all_I2_gDNA_temp for item in sublist]

    #define samples Read1
    samples_I2_gDNA = []
    for paths in all_I2_gDNA:

        samples_I2_gDNA.append(paths.split("/")[-1].split("_L0")[0])
    samples_I2_gDNA = list(set(samples_I2_gDNA))
    samples_I2_gDNA.sort()
    
    # Loop through the corresponding index1 files and remove the bad reads
    counter = 0
    for sample in samples_I2_gDNA:
        s_fastq = []
        for path in all_I2_gDNA:
            if "/"+ sample in path:
                s_fastq.append(path)

        # get all of the fastq lines from the fastqs
        lines = []
        for fsmp in s_fastq:
            with open(fsmp) as f:
                for line in f:
                    lines.append(line)
                
        # Write the edited sequence and trimmed file  to a fastq folder
        lane_counter = 0
        for fsamp in s_fastq:

            #write files with these edited barcodes ( these are used as the input into starcode)
            with open(filt_haveStart_gDNA + "/" + fsamp.split('/')[-1],"w") as f:
                # print('Beginning count: ' + str(file_endpoints_gDNA[counter][lane_counter]))
                # print('End count: '+ str(file_endpoints_gDNA[counter][lane_counter+1]))
                for j in range(file_endpoints_gDNA[counter][lane_counter],file_endpoints_gDNA[counter][lane_counter+1]):
                    if readsKeep_gDNA[counter][j] == 1:
                        f.write(lines[(4*j)])
                        f.write(lines[(4*j)+1]) # pulls the trimmed sequence
                        f.write(lines[(4*j)+2])
                        f.write(lines[(4*j)+3])
        
            lane_counter = lane_counter + 1

    
        counter = counter+1

