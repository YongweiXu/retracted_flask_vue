from pymongo import MongoClient

# 连接到 MongoDB
client = MongoClient('localhost', 27017)
db = client['pubmed']
collection = db['pubmed_collection']

# 获取集合中的所有文档的字段名
all_documents = collection.find({})
field_names = set()
for doc in all_documents:
    field_names.update(doc.keys())

print("All field names in the collection:", field_names)
