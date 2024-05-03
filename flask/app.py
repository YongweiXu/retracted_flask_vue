from collections import Counter

from flask import Flask, request, jsonify, send_file
from flask_pymongo import PyMongo
import pandas as pd
from collections import defaultdict
from wordcloud import WordCloud

app = Flask(__name__, static_url_path='/', static_folder='./../flask-dist', template_folder='./../flask-dist')
app.config["MONGO_URI"] = "mongodb://localhost:27017/pubmed"
mongo = PyMongo(app)

df = None

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

@app.route('/counts')
def data_counts():
    global df

    # 获取清洗前的数据量
    original_count = mongo.db.pubmed_collection.count_documents({})

    # 获取清洗后的数据量
    cleaned_count = len(df)

    # 构建返回的 JSON 对象
    counts = {
        'original_count': original_count,
        'cleaned_count': cleaned_count
    }

    return jsonify(counts)

import re

@app.route('/country')
def data_country():
    global df

    # 获取数据，暂存并处理仅保留 PL 不为空的数据
    df_country = df[df['PL'].notna()]
    # 对 PL 中的内容计数
    country_counts = df_country['PL'].value_counts().to_dict()

    # 将数据格式化为你需要的格式
    data = []
    for country, count in country_counts.items():
        # 使用正则表达式替换括号及括号内的内容为空
        cleaned_country = re.sub(r'\(.*?\)', '', country).strip()
        data.append({'name': cleaned_country, 'value': count})

    return jsonify(data)




# @app.route('/wordcloud')
# def generate_wordcloud():
#     global df
#
#     # 合并 "DE" 字段和 "WC" 字段中的关键词
#     de_wc_text = ''
#     for index, row in df.iterrows():
#         de_wc_text += row['DE'].replace(' ', ',') + ' ' + row['WC'].replace(' ', ',') + ' '
#
#     # 使用 WordCloud 生成词云图像，并设置背景为透明
#     wordcloud = WordCloud(width=800, height=400, background_color='rgba(0, 0, 0, 0)').generate(de_wc_text)
#
#     # 保存词云图像到临时文件
#     wordcloud_file = 'wordcloud.png'
#     wordcloud.to_file(wordcloud_file)
#
#     # 返回词云图像的 URL
#     return send_file(wordcloud_file, mimetype='image/png')








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


# @app.route('/PY')
# def count_by_year():
#     global df
#     # 检查数据是否已加载
#     if df is None:
#         return 'Data not loaded yet', 400
#     # 对 'PY' 字段进行计数统计
#     year_counts = df['PY'].value_counts().to_dict()
#     # 构建包含年份和计数的字典列表
#     result = [{'年份': year, '计数': count} for year, count in year_counts.items()]
#     return jsonify(result)


# @app.route('/PU')
# 读取数据并构建出版社年份退稿关系

# def generate_chart_data():
#     global df
#
#     # 计算每个出版社每年的退稿次数
#     publisher_year_counts = defaultdict(lambda: defaultdict(int))
#     for index, row in df.iterrows():
#         publishers = row['PU'].split(' & ')
#         year = row['PY']
#         for publisher in publishers:
#             publisher_year_counts[publisher][year] += 1
#
#     # 计算总退稿次数和阈值
#     total_count = sum(sum(counts.values()) for counts in publisher_year_counts.values())
#     threshold = total_count * 0.02
#
#     # 构建数据
#     years_sorted = sorted(df['PY'].unique())  # 年份递增排序
#     chart_data = ['product,' + ','.join(str(year) for year in years_sorted)]
#     other_counts = defaultdict(int)
#     for publisher, year_counts in publisher_year_counts.items():
#         if sum(year_counts.values()) >= threshold:
#             data = [str(year_counts.get(year, 0)) for year in years_sorted]
#             chart_data.append(f'{publisher},{",".join(data)}')
#         else:
#             for year, count in year_counts.items():
#                 other_counts[year] += count
#
#     # 添加 "other" 行
#     other_data = [str(other_counts[year]) for year in years_sorted]
#     chart_data.append(f'other,{",".join(other_data)}')
#
#     return jsonify(data=chart_data)


#关键词年份关系
# @app.route('/DE')
# def get_de_py_relationship():
#     global df
#     # 创建一个空的字典，用于存储关键词和每年的计数
#     de_py_relationship = {}
#
#     # 遍历 DataFrame 中的每一行
#     for index, row in df.iterrows():
#         # 将 DE 列拆分成单独的关键词
#         keywords = row['DE'].split('; ')
#         # 更新关键词和每年的计数
#         for keyword in keywords:
#             de_py_relationship.setdefault(keyword, Counter())
#             year = str(row['PY'])
#             de_py_relationship[keyword][year] += 1
#
#     # 取出计数值总数最高的10个关键词
#     top_10_keywords = dict(
#         sorted(de_py_relationship.items(), key=lambda item: sum(item[1].values()), reverse=True)[:10])
#
#     # 返回每个关键词每一年的变化情况
#     result = {}
#     for keyword, year_counter in top_10_keywords.items():
#         result[keyword] = dict(year_counter)
#
#     # 返回关系变化数据
#     return jsonify(result)


# #标题和被引用次数
# @app.route('/TIZ9')
# def get_echarts_data():
#     global df
#
#     # 将Z9列作为纵坐标轴的数据（dataAxi）
#     dataAxis = df['Z9'].tolist()
#     # 将TI列作为横坐标轴的数据（data）
#     data = df['TI'].tolist()
#
#     # 创建一个字典，将标题作为键，引用次数作为值
#     chart_data = {}
#     for title, reference_count in zip(dataAxis, data):
#         chart_data[title] = reference_count
#
#     # 按照引用次数降序排序
#     sorted_chart_data = dict(sorted(chart_data.items(), key=lambda item: item[1], reverse=True))
#
#     # 提取前20个标题和对应的引用次数
#     top_20_chart_data = dict(list(sorted_chart_data.items())[:20])
#
#     # 返回ECharts数据
#     return jsonify(top_20_chart_data)





# @app.route('/WC')
# def wc_counts():
#     # 创建一个字典来统计每个学科的出现次数
#     category_counts = defaultdict(int)
#     total_count = 0
#
#     # 统计所有学科的出现次数
#     for categories in df['WC']:
#         unique_categories = set()
#         for category in categories.split('; '):
#             for sub_category in category.split(', '):
#                 unique_categories.add(sub_category.strip())
#         total_count += len(unique_categories)
#         # 统计每个学科的出现次数
#         for category in unique_categories:
#             category_counts[category] += 1
#
#     # 筛选频次高于百分之三的学科
#     filtered_category_counts = {category: count for category, count in category_counts.items() if
#                                 count / total_count > 0.03}
#
#     # 计算 "other" 的频次
#     other_count = total_count - sum(filtered_category_counts.values())
#     filtered_category_counts['other'] = other_count
#
#     # 将结果格式化为所需的格式
#     result = [{'name': category, 'count': count} for category, count in filtered_category_counts.items()]
#
#     return jsonify(result)


# @app.route('/PI')
# def count_by_PI():
#     global df
#     # 检查数据是否已加载
#     if df is None:
#         return 'Data not loaded yet', 400
#     # 对 'PI' 字段进行计数统计
#     PI_counts = Counter(df['PI'])
#     # 取出前十个最常出现的PI值及其计数值，并将PI值转换为首字母大写
#     top_10_PI = {pi.capitalize(): count for pi, count in PI_counts.most_common(10)}
#     return jsonify(top_10_PI)








load_data()
if __name__ == '__main__':
    load_data()
    app.run(debug=True)
