import itertools
from collections import Counter

from flask import Flask, request, jsonify, send_file
from flask_pymongo import PyMongo
import pandas as pd
import re
from scipy.sparse import dok_matrix
from collections import defaultdict
from wordcloud import WordCloud

app = Flask(__name__, static_url_path='/', static_folder='./../flask-dist', template_folder='./../flask-dist')
app.config["MONGO_URI"] = "mongodb://localhost:27017/pubmed"
mongo = PyMongo(app)

df = None
##连接数据库并获取数据
def get_data_from_mongodb():
    global df  # 声明 df 为全局变量
    try:
        # 从名为 'retracted_papers' 的集合中读取数据
        data = mongo.db.pubmed_collection.find()
        # 转换为 DataFrame
        df = pd.DataFrame(list(data))
        # 删除 '_id' 字段
        df.drop('_id', axis=1, inplace=True)
        print(df)
        return df
    except Exception as e:
        print("Error while getting data from MongoDB:", e)
        return None

##存储数据到全局变量
def load_data():
    # 加载数据并保存到全局变量 df 中
    global df
    df = get_data_from_mongodb()
    if df is not None:
        return 'Data loaded successfully'
    else:
        return 'Failed to load data'

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response

##数据库记录计数
@app.route('/counts')
def data_counts():
    global df

    #退稿数据计数
    original_count = mongo.db.pubmed_collection.count_documents({})
    #通过数据计数
    pass_count = mongo.db.pubmed_review_collection.count_documents({})
    # 构建返回的 JSON 对象
    counts = {
        're_count': original_count,
        'pass_count' : pass_count
    }

    return jsonify(counts)



##国家计数
@app.route('/country')
def data_country():
    global df

    # 获取数据，暂存并处理仅保留 PL 不为空的数据
    df_country = df[df['PL'].notna()]

    # 对 PL 中的内容进行处理：先将所有 'England' 替换为 'United Kingdom'，然后去除括号及括号内的内容
    df_country['PL'] = df_country['PL'].apply(
        lambda x: re.sub(r'\(.*?\)', '', x).strip().replace('England', 'United Kingdom'))

    # 对 PL 中的内容计数
    country_counts = df_country['PL'].value_counts().to_dict()

    # 将数据格式化为你需要的格式
    data = []
    for country, count in country_counts.items():
        data.append({'name': country, 'value': count})

    return jsonify(data)

##共现矩阵
@app.route('/Coc')
def cooccurrence_matrix_with_max_combination_json():
    global df

    # 获取 MH 数据
    mh_data = df['MH'].dropna().tolist()

    # 初始化共现矩阵字典
    matrix_dict = {}

    # 遍历每个子列表
    for terms in mh_data:
        unique_terms = set(terms)
        for term1 in unique_terms:
            for term2 in unique_terms:
                if term1 != term2:
                    # 更新共现次数
                    if term1 < term2:
                        if (term1, term2) not in matrix_dict:
                            matrix_dict[(term1, term2)] = 0
                        matrix_dict[(term1, term2)] += 1
                    else:
                        if (term2, term1) not in matrix_dict:
                            matrix_dict[(term2, term1)] = 0
                        matrix_dict[(term2, term1)] += 1

    # 找到共现次数前400的组合
    top_400_combinations = sorted(matrix_dict.items(), key=lambda x: x[1], reverse=True)[:42]

    # 初始化输出格式的数据
    output_data = {'term1': [], 'term2': [], 'matrix': []}

    # 添加所有术语到输出数据中
    all_terms = set()
    for (term1, term2), _ in top_400_combinations:
        all_terms.add(term1)
        all_terms.add(term2)

    output_data['term1'] = sorted(all_terms)
    output_data['term2'] = sorted(all_terms)

    # 构建共现次数矩阵
    for term1 in output_data['term1']:
        matrix_row = []
        for term2 in output_data['term2']:
            if term1 < term2:
                matrix_row.append(matrix_dict.get((term1, term2), 0))
            else:
                matrix_row.append(matrix_dict.get((term2, term1), 0))
        output_data['matrix'].append(matrix_row)

    return jsonify(output_data)

##年份计数
@app.route('/PY')
def count_by_year():
    global df
    # 检查数据是否已加载
    if df is None:
        return 'Data not loaded yet', 400
    # 对 'PY' 字段进行计数统计
    df_cleaned = df.dropna(subset=['DCOM'])
    year_counts = df_cleaned['DCOM'].astype(str).str[:4].value_counts().to_dict()
    # 构建包含年份和计数的字典列表
    result = [{'年份': year, '计数': count} for year, count in year_counts.items()]
    return jsonify(result)


##期刊年份退稿关系
@app.route('/PU')
#读取数据并构建出版社年份退稿关系
def generate_chart_data():
    global df
    df = df.dropna(subset=['DCOM'])
    df = df.dropna(subset=['TA'])
    # 计算每个出版社每年的退稿次数
    publisher_year_counts = defaultdict(lambda: defaultdict(int))
    for index, row in df.iterrows():
        publishers = re.sub(r'\([^)]*\)', '', row['TA']).split(' & ')  # 删除括号内的内容
        year = str(row['DCOM'])[:4]  # 只取年份的前四位数字
        for publisher in publishers:
            publisher_year_counts[publisher][year] += 1

    # 计算总退稿次数和阈值
    total_count = sum(sum(counts.values()) for counts in publisher_year_counts.values())
    threshold = total_count * 0.02

    # 构建数据
    years_sorted = sorted(set(str(row['DCOM'])[:4] for _, row in df.iterrows()))  # 年份递增排序
    chart_data = ['product,' + ','.join(years_sorted)]
    other_counts = defaultdict(int)
    for publisher, year_counts in publisher_year_counts.items():
        if sum(year_counts.values()) >= threshold:
            data = [str(year_counts.get(year, 0)) for year in years_sorted]
            chart_data.append(f'{publisher},{",".join(data)}')
        else:
            for year, count in year_counts.items():
                other_counts[year] += count

    # 添加 "other" 行
    other_data = [str(other_counts[year]) for year in years_sorted]
    chart_data.append(f'other,{",".join(other_data)}')

    return jsonify(data=chart_data)

#关键词年份关系
@app.route('/DE')
def get_de_dcom_relationship():
    global df
    df = df.dropna(subset=['DCOM'])
    df = df.dropna(subset=['MH'])
    # 创建一个空的字典，用于存储关键词和每年的计数
    de_dcom_relationship = {}

    # 遍历 DataFrame 中的每一行
    for index, row in df.iterrows():
        # 获取关键词列表
        keywords = row['MH']
        # 更新关键词和每年的计数
        for keyword in keywords:
            de_dcom_relationship.setdefault(keyword, Counter())
            year = str(row['DCOM'])[:4]  # 只取年份的前四位数字
            de_dcom_relationship[keyword][year] += 1

    # 取出计数值总数最高的10个关键词
    top_10_keywords = dict(
        sorted(de_dcom_relationship.items(), key=lambda item: sum(item[1].values()), reverse=True)[:15])

    # 返回每个关键词每一年的变化情况
    result = {}
    for keyword, year_counter in top_10_keywords.items():
        result[keyword] = dict(year_counter)

    # 返回关系变化数据
    return jsonify(result)

#标题和被引用次数
@app.route('/TIZ9')
def get_echarts_data():
    global df
    df = df.dropna(subset=['VI'])
    df = df.dropna(subset=['MH'])
    # 将Z9列作为纵坐标轴的数据（dataAxi）
    dataAxis = df['VI'].tolist()
    # 将TI列作为横坐标轴的数据（data）
    data = df['TI'].tolist()

    # 创建一个字典，将标题作为键，引用次数作为值
    chart_data = {}
    for title, reference_count in zip(dataAxis, data):
        chart_data[title] = reference_count

    # 按照引用次数降序排序
    sorted_chart_data = dict(sorted(chart_data.items(), key=lambda item: item[1], reverse=True))

    # 提取前20个标题和对应的引用次数
    top_20_chart_data = dict(list(sorted_chart_data.items())[:20])

    # 返回ECharts数据
    return jsonify(top_20_chart_data)


##关键词词云图
@app.route('/wordcloud')
def generate_wordcloud():
    global df
    df = df.dropna(subset=['MH'])
    # 合并 "MH" 列中的列表内容为一个字符串
    mh_text = ''
    for index, row in df.iterrows():
        # 将列表内容拼接成一个字符串
        mh_text += ', '.join(row['MH']) + ' '

    # 使用 WordCloud 生成词云图像，并设置背景为透明
    wordcloud = WordCloud(width=800, height=400, background_color='rgba(0, 0, 0, 0)').generate(mh_text)

    # 保存词云图像到临时文件
    wordcloud_file = 'wordcloud.png'
    wordcloud.to_file(wordcloud_file)

    # 返回词云图像的 URL
    return send_file(wordcloud_file, mimetype='image/png')

##被退稿最多的作者
@app.route('/AU')
def top_rejected_authors():
    global df
    df = df.dropna(subset=['FAU'])
    # 将FAU列展开为单独的行
    authors = df.explode('FAU')

    # 计算每位作者的退稿次数
    rejected_authors = authors['FAU'].value_counts()

    # 按退稿次数降序排序
    sorted_rejected_authors = rejected_authors.sort_values(ascending=False)

    # 返回退稿次数最高的前10位作者及其退稿次数
    top_rejected_authors = sorted_rejected_authors.head(15).to_dict()

    return jsonify(top_rejected_authors)




# @app.route('/query')
# def query_data():
#     global df
#     # 获取查询参数
#     query_params = request.args.to_dict()
#     # 如果查询参数为空，则返回原始数据
#     if not query_params:
#         df = get_data_from_mongodb()
#         return jsonify(df.to_dict(orient='records'))
#     # 如果只有一个查询字段，则使用该字段进行查询
#     elif len(query_params) == 1:
#         field, value = list(query_params.items())[0]
#         df = get_data_from_mongodb()
#         query_result = df[df[field] == value]
#         return jsonify(query_result.to_dict(orient='records'))
#     # 如果有多个查询字段，则使用联合查询的方式
#     else:
#         df = get_data_from_mongodb()
#         query = " & ".join([f"{field} == '{value}'" for field, value in query_params.items()])
#         query_result = df.query(query)
#         return jsonify(query_result.to_dict(orient='records'))


load_data()
if __name__ == '__main__':
    load_data()
    app.run(debug=True)
