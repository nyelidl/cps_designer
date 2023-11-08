#!/usr/bin/python3

from tqdm import tqdm



def generate_combinations():
    amino_acids = [chr(i) for i in range(65, 85)]  # 使用A到T代表20种氨基酸
    all_combinations = set()
    

    # 生成所有可能的五氨基酸组合
    for comb in tqdm(itertools.product(amino_acids, repeat=5)):
        # 将组合转换为字符串
        comb_str = ''.join(comb)
        if 'CC' in comb_str:
            if_new = False
            break

        
        # 检查所有旋转
        is_new = True
        for i in range(5):
            if comb_str[i:] + comb_str[:i] in all_combinations:
                is_new = False
                break
        
        # 如果考虑翻转，检查镜像对称的组合
        reversed_str = comb_str[::-1]
        for i in range(5):
            if reversed_str[i:] + reversed_str[:i] in all_combinations:
                is_new = False
                break
        
        if is_new:
            all_combinations.add(comb_str)

    return all_combinations

import itertools
combinations = generate_combinations()
print(len(combinations))
#print(combinations)

