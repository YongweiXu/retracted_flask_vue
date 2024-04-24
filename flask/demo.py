import requests

# 定义 Flask 应用的地址
base_url = 'http://127.0.0.1:5000'


# 访问所有路由的函数
def access_all_routes():
    # 查询所有数据
    response = requests.get(f'{base_url}/query')
    print('查询所有数据:', response.json())

    # 按年份统计文章数量
    response = requests.get(f'{base_url}/PY')
    print('按年份统计文章数量:', response.json())

    # 获取出版社年份退稿关系数据
    response = requests.get(f'{base_url}/PU')
    print('获取出版社年份退稿关系数据:', response.json())

    # 获取关键词年份关系数据
    response = requests.get(f'{base_url}/DE')
    print('获取关键词年份关系数据:', response.json())

    # 获取标题和被引用次数数据
    response = requests.get(f'{base_url}/TIZ9')
    print('获取标题和被引用次数数据:', response.json())

    # 获取每个学科的文章数量统计
    response = requests.get(f'{base_url}/WC')
    print('获取每个学科的文章数量统计:', response.json())


# 调用函数访问所有路由
access_all_routes()
