def get_primer(seq, n, gc_min=50, gc_max=60, tm_min=50, tm_max=60, max_self_compliment=4):
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

    for subseq in moving_window(seq, n):
        if gc_content(subseq) >= gc_min / 100 and gc_content(subseq) <= gc_max / 100 and calculate_tm(subseq) >= tm_min and calculate_tm(subseq) <= tm_max and not is_self_compliment(subseq, max_self_compliment):
            result.append(subseq)
            return subseq
    if not result:
        return False