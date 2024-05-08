def most_frequent_combination(lists_of_strings):
    def co_occurrence_matrix(lists_of_strings):
        # 构建词汇表
        vocabulary = set()
        for string_list in lists_of_strings:
            vocabulary.update(string_list)
        vocabulary = sorted(list(vocabulary))

        # 初始化共现矩阵
        matrix = [[0] * len(vocabulary) for _ in range(len(vocabulary))]
        # 填充共现矩阵
        for string_list in lists_of_strings:
            for i, word1 in enumerate(vocabulary):
                if word1 in string_list:
                    for j, word2 in enumerate(vocabulary):
                        if word2 in string_list and word1 != word2:
                            matrix[i][j] += 1

        return matrix, vocabulary

    matrix, vocabulary = co_occurrence_matrix(lists_of_strings)

    max_count = 0
    max_combinations = []
    for i in range(len(vocabulary)):
        for j in range(i + 1, len(vocabulary)):
            if matrix[i][j] > max_count:
                max_count = matrix[i][j]
                max_combinations = [(vocabulary[i], vocabulary[j])]
            elif matrix[i][j] == max_count:
                max_combinations.append((vocabulary[i], vocabulary[j]))

    return max_count, max_combinations