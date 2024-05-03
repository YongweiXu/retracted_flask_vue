import time
from datetime import datetime

from Bio import Entrez, Medline
from pymongo import MongoClient
from urllib3.exceptions import IncompleteRead

# 设置你的 API 密钥和电子邮件
API_KEY = "c61d19c529b4ca05f8dc340f265d9f7d4408"
EMAIL = "xyw0206070@outlook.com"

Entrez.api_key = API_KEY
Entrez.email = EMAIL

# 定义开始年份和结束年份
start_year = 1951
current_year = datetime.now().year

# 定义查询之间的睡眠时间和重试延迟
retry_delay = 600
max_retries = 3
timeout = 30  # 设置超时时间为30秒

# 初始化变量
records = []
pmids = set()

# 连接到 MongoDB
MONGO_URI = "mongodb://localhost:27017/"
DATABASE_NAME = "pubmed"
client = MongoClient(MONGO_URI)
db = client[DATABASE_NAME]
collection = db.pubmed_review_collection

# 主循环
query_start_year = start_year
while query_start_year <= current_year:
    retries = 0
    while retries < max_retries:
        try:
            # 构建搜索条件
            query_end_year = query_start_year + 1
            term = f'Review[pt] AND {query_start_year}:{query_end_year}[dp]'

            # 在 PubMed 中搜索
            handle = Entrez.esearch(db='pubmed', term=term, usehistory='y', timeout=timeout)
            search_results = Entrez.read(handle)
            webenv = search_results['WebEnv']
            query_key = search_results['QueryKey']
            total_records = int(search_results['Count'])

            # 获取所有记录
            batch_size = 2000  # 每年最多获取2000条记录
            for start in range(0, total_records, batch_size):
                end = min(total_records, start + batch_size)
                handle = Entrez.efetch(db='pubmed', webenv=webenv, query_key=query_key, retstart=start, retmax=batch_size, retmode='text', rettype='medline', timeout=timeout)
                batch_records = list(Medline.parse(handle))

                # 存储记录
                for record in batch_records:
                    pmid = record.get('PMID')
                    if pmid not in pmids:
                        records.append(record)
                        pmids.add(pmid)

                # 打印进度
                print(f"已发现 {len(records)} 条有效记录，{query_start_year}-{query_end_year} 年。此次发现 {len(batch_records)} 条记录。")

                # 清除本次迭代的记录列表，以免重复添加
                batch_records.clear()

                if len(records) >= 2000:  # 如果超过2000条记录，跳出循环
                    break

            break

        except IncompleteRead:
            retries += 1
            time.sleep(retry_delay)
            continue

        except Exception as e:
            retries += 1
            time.sleep(retry_delay)
            print(f"错误：{e}")
            continue

    if retries == max_retries:
        print(f"尝试 {max_retries} 次仍然失败，放弃查询 {query_start_year}-{query_end_year} 年。")

    # 移动到下一年
    query_start_year = query_end_year + 1

print("查询结束。")

# 将数据写入 MongoDB
try:
    if records:
        collection.insert_many(records)
        print(f"成功插入 {len(records)} 条记录到数据库。")
    else:
        print("没有记录需要插入。")
except Exception as e:
    print(f"插入数据库时出错：{e}")
