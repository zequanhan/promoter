def encode_sequences(sequences, max_len=99):
    # 核苷酸编码，加上一个用于填充的0编码
    nucleotide_to_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    encoded_seqs = []
    for seq in sequences:
        # 编码序列
        encoded_seq = [nucleotide_to_int.get(nuc, 0) for nuc in seq]
        # 如果序列长度小于max_len，则用0填充
        padded_seq = encoded_seq + [0] * (max_len - len(encoded_seq))
        encoded_seqs.append(padded_seq[:max_len])  # 确保序列长度不超过max_len
    return np.array(encoded_seqs)
