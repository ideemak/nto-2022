def sliding_window(seq, n):
    def gc_content(seq):
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq)

    def moving_window(seq, n):
        for i in range(0, len(seq)):
            if i + n <= len(seq):
                yield seq[i:i+n]
            else:
                break

    result = []
    for subseq in moving_window(seq, n):
        result.append(gc_content(subseq))
    return result