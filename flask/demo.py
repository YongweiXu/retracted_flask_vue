from flask import jsonify
from collections import defaultdict
import itertools

def cooccurrence_matrix_with_max_combination_json():
    global df

    # 获取 MH 数据
    mh_data = df['MH'].dropna().tolist()

    # 初始化共现矩阵
    matrix = defaultdict(int)

    # 遍历每个子列表
    for terms in mh_data:
        unique_terms = sorted(set(terms))
        # 更新共现次数
        for term1, term2 in itertools.combinations(unique_terms, 2):
            matrix[(term1, term2)] += 1

    # 找到共现次数前 20 的组合
    top_20_combinations = sorted(matrix.items(), key=lambda x: x[1], reverse=True)[:20]

    # 创建二维矩阵
    matrix_data = [[""] * (len(top_20_combinations) + 1) for _ in range(len(top_20_combinations) + 1)]

    # 添加行和列标签
    matrix_data[0][0] = "Terms"
    for i, (terms, _) in enumerate(top_20_combinations):
        matrix_data[0][i+1] = f"{terms[0]} 和 {terms[1]}"
        matrix_data[i+1][0] = f"{terms[0]} 和 {terms[1]}"

    # 填充矩阵数据
    for i, (terms, value) in enumerate(top_20_combinations):
        term1, term2 = terms
        matrix_data[i+1][i+1] = str(value)

    return matrix_data
x = cooccurrence_matrix_with_max_combination_json()
print(x)