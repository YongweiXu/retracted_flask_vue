from collections import Counter
from flask import Flask, request, jsonify, send_file
from flask_pymongo import PyMongo
import pandas as pd
import re
from Bio import Entrez, Medline
import time
from http.client import IncompleteRead
from collections import defaultdict
from wordcloud import WordCloud
from predict_text import predict_text




app = Flask(__name__, static_url_path='/', static_folder='./../flask-dist', template_folder='./../flask-dist')
app.config["MONGO_URI"] = "mongodb://localhost:27017/pubmed"
mongo = PyMongo(app)


API_KEY = "c61d19c529b4ca05f8dc340f265d9f7d4408"
EMAIL = "xyw0206070@outlook.com"

#连接数据库并获取数据
def get_data_from_mongodb():
    try:

        data = mongo.db.pubmed_collection.find()
        df = pd.DataFrame(list(data))
        df.drop('_id', axis=1, inplace=True)
        return df
    except Exception as e:
        return "Error while getting data from MongoDB:", e


@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response


##数据库记录计数
@app.route('/counts', methods=['GET'])
def data_counts():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500


    original_count = mongo.db.pubmed_collection.count_documents({})
    pass_count = mongo.db.pubmed_review_collection.count_documents({})
    counts = {
        're_count': original_count,
        'pass_count': pass_count
    }

    return jsonify(counts)


#国家计数
@app.route('/country', methods=['GET'])
def data_country():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500

    df_country = df[df['PL'].notna()]

    #England替换为United Kingdom
    df_country['PL'] = df_country['PL'].apply(
        lambda x: re.sub(r'\(.*?\)', '', x).strip().replace('England', 'United Kingdom'))

    country_counts = df_country['PL'].value_counts().to_dict()

    data = []
    for country, count in country_counts.items():
        data.append({'name': country, 'value': count})

    return jsonify(data)


##共现矩阵
@app.route('/Coc', methods=['GET'])
def cooccurrence_matrix_with_max_combination_json():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500

    mh_data = df['MH'].dropna().tolist()

    # 初始化共现矩阵字典
    matrix_dict = {}

    for terms in mh_data:
        unique_terms = set(terms)
        for term1 in unique_terms:
            for term2 in unique_terms:
                if term1 != term2:
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

    output_data = {'term1': [], 'term2': [], 'matrix': []}

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


#年份计数
@app.route('/PY', methods=['GET'])
def count_by_year():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500
    df_cleaned = df.dropna(subset=['DCOM'])
    year_counts = df_cleaned['DCOM'].astype(str).str[:4].value_counts().to_dict()
    result = [{'年份': year, '计数': count} for year, count in year_counts.items()]
    return jsonify(result)


#期刊年份退稿关系
@app.route('/PU', methods=['GET'])
def generate_chart_data():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500
    df = df.dropna(subset=['DCOM'])
    df = df.dropna(subset=['TA'])
    publisher_year_counts = defaultdict(lambda: defaultdict(int))
    for index, row in df.iterrows():
        publishers = re.sub(r'\([^)]*\)', '', row['TA']).split(' & ')  # 删除括号内的内容
        year = str(row['DCOM'])[:4]  # 只取年份的前四位数字
        for publisher in publishers:
            publisher_year_counts[publisher][year] += 1

    # 计算总退稿次数和阈值
    total_count = sum(sum(counts.values()) for counts in publisher_year_counts.values())
    threshold = total_count * 0.02

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

    # 添加other
    other_data = [str(other_counts[year]) for year in years_sorted]
    chart_data.append(f'other,{",".join(other_data)}')

    return jsonify(data=chart_data)


#关键词年份关系
@app.route('/DE', methods=['GET'])
def get_de_dcom_relationship():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500
    df = df.dropna(subset=['DCOM'])
    df = df.dropna(subset=['MH'])

    de_dcom_relationship = {}

    for index, row in df.iterrows():
        keywords = row['MH']
        for keyword in keywords:
            de_dcom_relationship.setdefault(keyword, Counter())
            year = str(row['DCOM'])[:4]
            de_dcom_relationship[keyword][year] += 1

    top_10_keywords = dict(
        sorted(de_dcom_relationship.items(), key=lambda item: sum(item[1].values()), reverse=True)[:15])

    result = {}
    for keyword, year_counter in top_10_keywords.items():
        years_data = {}
        for year in range(1990, 2025):
            year_str = str(year)
            years_data[year_str] = year_counter.get(year_str, 0)  # 如果年份不存在，频次用0补充
        result[keyword] = years_data

    return jsonify(result)


#标题和被引用次数
@app.route('/TIZ9', methods=['GET'])
def get_echarts_data():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500
    df = df.dropna(subset=['VI'])
    df = df.dropna(subset=['MH'])
    dataAxis = df['VI'].tolist()
    data = df['TI'].tolist()

    chart_data = {}
    for title, reference_count in zip(dataAxis, data):
        chart_data[title] = reference_count

    sorted_chart_data = dict(sorted(chart_data.items(), key=lambda item: item[1], reverse=True))

    top_20_chart_data = dict(list(sorted_chart_data.items())[:20])

    return jsonify(top_20_chart_data)


#关键词词云图
@app.route('/wordcloud', methods=['GET'])
def generate_wordcloud():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500
    df = df.dropna(subset=['MH'])
    mh_text = ''
    for index, row in df.iterrows():
        mh_text += ', '.join(row['MH']) + ' '

    # 使用 WordCloud 生成词云图像，并设置背景为透明
    wordcloud = WordCloud(width=800, height=400, background_color='rgba(0, 0, 0, 0)').generate(mh_text)
    # 保存词云图像到临时文件
    wordcloud_file = 'wordcloud.png'
    wordcloud.to_file(wordcloud_file)
    # 返回词云图像的 URL
    return send_file(wordcloud_file, mimetype='image/png')


#被退稿最多的作者
@app.route('/AU', methods=['GET'])
def top_rejected_authors():
    df = get_data_from_mongodb()
    if df is None:
        return 'Failed to load data', 500
    df = df.dropna(subset=['FAU'])
    authors = df.explode('FAU')

    rejected_authors = authors['FAU'].value_counts()
    sorted_rejected_authors = rejected_authors.sort_values(ascending=False)
    top_rejected_authors = sorted_rejected_authors.head(15).to_dict()

    return jsonify(top_rejected_authors)


#数据轮播图
@app.route('/data', methods=['GET'])
def get_data():
    global sampled_data

    # 将 sampled_data 重置为 None，以便刷新
    sampled_data = None
    # 如果抽样数据未被设置或数据量发生变化，重新进行随机抽样
    if sampled_data is None or len(sampled_data) != 15:
        df = get_data_from_mongodb()
        if df is None:
            return '无法加载数据', 500

        df['DCOM'] = pd.to_datetime(df['DCOM'], format='%Y%m%d', errors='coerce')
        df['DCOM'] = df['DCOM'].dt.strftime('%Y-%m-%d')
        df['AU'] = df['AU'].apply(lambda x: ', '.join(x) if isinstance(x, list) else x)
        # 随机抽样
        if len(df) >= 15:
            sampled_data = df[['PMID', 'STAT', 'DCOM', 'IS', 'TI', 'AU']].dropna().sample(n=15,
                                                                                          replace=False).values.tolist()
        else:
            sampled_data = df[['PMID', 'STAT', 'DCOM', 'IS', 'TI', 'AU']].dropna().values.tolist()  # 如果数据不足15条，返回所有数据

    return jsonify(sampled_data)



@app.route('/predict', methods=['POST'])
def prediction():
    data = request.get_json()
    text_to_predict = data['text']
    prediction_result = predict_text(text_to_predict)
    return jsonify(prediction_result)



@app.route('/get_pubmed_records', methods=['POST'])
def get_pubmed_records():
    # 从前端获取关键词
    search_term = request.json['term']
    Entrez.api_key = API_KEY
    Entrez.email = EMAIL
    # 最多获取记录数
    max_records = 10000
    retries = 0
    while retries < 3:
        try:
            # 在 PubMed 中搜索
            handle = Entrez.esearch(db='pubmed', term=search_term, usehistory='y')
            search_results = Entrez.read(handle)
            webenv = search_results['WebEnv']
            query_key = search_results['QueryKey']
            total_records = min(int(search_results['Count']), max_records)

            # 获取所有记录
            batch_size = min(10000, max_records)
            records = []
            for start in range(0, total_records, batch_size):
                end = min(total_records, start + batch_size)
                handle = Entrez.efetch(db='pubmed', webenv=webenv, query_key=query_key, retstart=start,
                                       retmax=batch_size, retmode='text', rettype='medline')
                batch_records = list(Medline.parse(handle))

                # 提取所需字段内容
                for record in batch_records:
                    record_info = {}
                    for field in ['PMID', 'STAT', 'DCOM', 'IS', 'TI', 'AU', 'AB']:
                        record_info[field] = record.get(field, '')
                    records.append(record_info)

            return jsonify({'records': records})

        except IncompleteRead:
            retries += 1
            time.sleep(2)
            continue

        except Exception as e:
            retries += 1
            time.sleep(2)
            print(f"错误：{e}")
            continue

    return jsonify({'error': 'Failed to retrieve PubMed records'})


if __name__ == '__main__':
    app.run(debug=True)
