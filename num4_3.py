def get_reverse_primer(seq, seq_f, n, amplicon_min, amplicon_max, gc_min=50, gc_max=60, tm_min=50, tm_max=60, max_self_compliment=4):
    def complement(seq):
        dct = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        return "".join([dct[k] for k in seq])

    def get_subseq_adress(seq, subseq):
        for i in range(0, len(seq)):
            if seq[i:i+len(subseq)] == subseq:
                return [i, i+len(subseq)-1]
        return False

    def gc_content(seq):
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq)

    def moving_window(seq, n):
        for i in range(0, len(seq)):
            if i + n <= len(seq):
                yield seq[i:i+n]
            else:
                break

    def calculate_tm(string):
        string = string.upper()
        w = string.count('A') 
        x = string.count('T')
        y = string.count('G') 
        z = string.count('C')
        Tm = 64.9 + 41 * (y + z - 16.4) / (w + x + y + z)
        return round(Tm)

    def is_self_compliment(string, max_n):
        string = string.upper()
        dct = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        string_complement = "".join([dct[k] for k in string[-max_n:]])
        string_complement = string_complement[::-1]
        return string_complement in string
    result = []

    seq_f_adress = get_subseq_adress(seq, seq_f)
    seq_f_adress_comp = [len(seq) - seq_f_adress[1] - 1, len(seq) - seq_f_adress[0] - 1]
    margin = amplicon_min - 2 * n
    seq_comp_back = complement(seq)[::-1]
    mbl = amplicon_max - amplicon_min
    seq_1 = seq_comp_back[0:seq_f_adress_comp[0] - margin-1]
    seq_2 = seq_comp_back[seq_f_adress_comp[1] + margin+1:]
    if len(seq_1) > mbl:
        seq_1 = seq_1[len(seq_1) - mbl:]
    
    if len(seq_2) > mbl:
        seq_2 = seq_2[:mbl]
    
    for subseq in moving_window(seq_1, n):
        if gc_content(subseq) >= gc_min / 100 and gc_content(subseq) <= gc_max / 100 and calculate_tm(subseq) >= tm_min and calculate_tm(subseq) <= tm_max and not is_self_compliment(subseq, max_self_compliment):
            result.append(subseq)
            return subseq
    if not result:
        for subseq in moving_window(seq_2, n):
            if gc_content(subseq) >= gc_min / 100 and gc_content(subseq) <= gc_max / 100 and calculate_tm(subseq) >= tm_min and calculate_tm(subseq) <= tm_max and not is_self_compliment(subseq, max_self_compliment):
                result.append(subseq)
            
            return subseq
    
    if not result:
        return False

