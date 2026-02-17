from collections import defaultdict

from collections import defaultdict

def build_full_length_alignments(sam_path, ref_fasta_path):
    """
    Efficiently builds fully gapped reference-anchored alignments.
    Insertions relative to reference create new columns; deletions are '-'.
    Returns a tuple of (reference_string, dict of read_name -> aligned_read_string)
    """
    # Load reference
    with open(ref_fasta_path) as f:
        ref = "".join(line.strip() for line in f if not line.startswith(">")).upper()
    ref_len = len(ref)

    reads_data = []
    read_names = []
    insertions_at_pos = defaultdict(list)

    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            read_name = fields[0]
            flag = int(fields[1])
            if flag & 4:  # unmapped
                continue

            pos = int(fields[3]) - 1
            cigar = fields[5]
            seq = fields[9].upper()
            if cigar == "*":
                continue

            ref_i = pos
            read_i = 0
            num = ""
            read_tuples = []  # list of (ref_pos, base_or_insertion)

            for char in cigar:
                if char.isdigit():
                    num += char
                    continue
                length = int(num)
                num = ""

                if char in "M=X":
                    for _ in range(length):
                        read_tuples.append((ref_i, seq[read_i]))
                        ref_i += 1
                        read_i += 1
                elif char == "D":
                    for _ in range(length):
                        read_tuples.append((ref_i, "-"))
                        ref_i += 1
                elif char == "I":
                    ins_seq = seq[read_i:read_i+length]
                    # associate insertion with the *previous* reference position (i.e. the insertion
                    # appears after the last consumed reference base). If ref_i==0 (insertion at very
                    # start), keep it at ref_i==0 so we still have a slot.
                    ins_pos = ref_i - 1 if ref_i > 0 else 0
                    insertions_at_pos[ins_pos].append(ins_seq)
                    # Mark it as insertion explicitly (we'll handle placement later)
                    read_tuples.append((ins_pos, ("INS", ins_seq)))
                    read_i += length
                elif char == "S":
                    read_i += length
                elif char == "H":
                    pass

            reads_data.append(read_tuples)
            read_names.append(read_name)

    # Determine final alignment length
    max_inserts = {}
    for pos, ins_lists in insertions_at_pos.items():
        max_len = max(len(ins_seq) for ins_seq in ins_lists)
        max_inserts[pos] = max_len

    alignment_length = ref_len + sum(max_inserts.values())

    # Map reference positions to final alignment positions
    ref_to_aln = {}
    aln_pos = 0
    for i in range(ref_len):
        ref_to_aln[i] = aln_pos
        aln_pos += 1
        if i in max_inserts:
            aln_pos += max_inserts[i]

    # Build reference row
    ref_row = ["-"] * alignment_length
    for i, base in enumerate(ref):
        ref_row[ref_to_aln[i]] = base
    ref_str = "".join(ref_row)

    # Build read rows
    read_dict = {}
    for read_name, read_tuples in zip(read_names, reads_data):
        row = ["-"] * alignment_length
        for ref_pos, base in read_tuples:
            # deletions and matches: single-character bases or '-'
            if isinstance(base, str) and len(base) == 1:
                aln_idx = ref_to_aln[ref_pos]
                row[aln_idx] = base
            else:
                # insertion entry stored as ("INS", ins_seq)
                tag, ins_seq = base
                # If insertion associated with ref_pos beyond last base (shouldn't normally happen),
                # put it after the last reference base:
                if ref_pos not in ref_to_aln:
                    # place at end of alignment (safest fallback)
                    aln_idx = alignment_length - 1
                else:
                    # insertion columns follow the reference base column
                    aln_idx = ref_to_aln[ref_pos] + 1

                # write insertion sequence into the insertion columns reserved for this ref_pos
                for j, b in enumerate(ins_seq):
                    # safety: guard against overrunning alignment (should not happen if max_inserts computed)
                    if aln_idx + j < len(row):
                        row[aln_idx + j] = b
        read_dict[read_name] = "".join(row)

    return ref_str, read_dict


def add_pipes(piped_reference, ref_aln, read_alns_dict):
    """
    Add pipe characters from piped_reference to reference and reads.
    read_alns_dict: dict of read_name -> aligned read string
    Returns: (piped_ref_aln, dict of read_name -> piped aligned string)
    """
    read_names = list(read_alns_dict.keys())
    read_alns = [read_alns_dict[name] for name in read_names]
    read_count = len(read_alns)

    piped_ref_aln = []
    alns_transposed = []

    ref_iterator = iter(piped_reference)

    for ref_ali_base, read_ali_bases in zip(ref_aln, zip(*read_alns)):
        if ref_ali_base != '-':
            ref_base = next(ref_iterator)
            if ref_base == '|':
                piped_ref_aln.append('|')
                alns_transposed.append(['|'] * read_count)
                ref_base = next(ref_iterator)
            piped_ref_aln.append(ref_ali_base)
            alns_transposed.append(read_ali_bases)
        else:
            piped_ref_aln.append(ref_ali_base)
            alns_transposed.append(read_ali_bases)

    piped_ref_aln = ''.join(piped_ref_aln)
    piped_read_alns = {
        name: ''.join(bases)
        for name, bases in zip(read_names, zip(*alns_transposed))
    }

    return piped_ref_aln, piped_read_alns
