def split_fasta_by_length(fasta_file_path, num_splits):
    # 读取FASTA文件并计算每个序列的长度和总长度
    sequences = []
    total_length = 0
    with open(fasta_file_path, 'r') as fasta:
        seq = ''
        for line in fasta:
            if line.startswith('>'):
                if seq:
                    sequences.append((header, seq))
                    total_length += len(seq)
                    seq = ''
                header = line.strip()
            else:
                seq += line.strip()
        # 添加最后一个序列
        if seq:
            sequences.append((header, seq))
            total_length += len(seq)

    # 计算每个子文件应含有的大致总序列长度
    length_per_split = total_length // num_splits

    # 分割序列并写入新的子文件
    split_idx = 1
    current_length = 0
    split_sequences = []
    split_files = []

    for header, seq in sequences:
        if current_length < length_per_split:
            split_sequences.append((header, seq))
            current_length += len(seq)
        else:
            # 当前子文件已满，写入子文件
            with open(f'G:\kfwang\miniTest\database\split_{split_idx}.fasta', 'w') as out:
                for split_header, split_seq in split_sequences:
                    out.write(f'{split_header}\n{split_seq}\n')
            split_files.append(f'split_{split_idx}.fasta')
            
            # 重置计数并开始下一个子文件
            split_idx += 1
            split_sequences = [(header, seq)]  # 开始新子文件的序列列表包含当前序列
            current_length = len(seq)  # 重置长度计数为当前序列长度

    # 处理最后一个子文件（包含所有剩余的序列）
    if split_sequences:
        with open(f'G:\kfwang\miniTest\database\split_{split_idx}.fasta', 'w') as out:
            for split_header, split_seq in split_sequences:
                out.write(f'{split_header}\n{split_seq}\n')
        split_files.append(f'split_{split_idx}.fasta')

    return split_files


fasta_file_path = r'G:\kfwang\miniTest\database\target_group_pace2.fasta'
num_splits = 3
split_files = split_fasta_by_length(fasta_file_path, num_splits)
print(split_files)