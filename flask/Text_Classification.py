import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Embedding, Conv1D, GlobalMaxPooling1D, Dense, Dropout
from pymongo import MongoClient
import tensorflow as tf
import os

# 连接到MongoDB
client = MongoClient('mongodb://localhost:27017/')

# 选择数据库
db = client['pubmed']

# 选择集合
collection1 = db['pubmed_collection']
collection2 = db['pubmed_review_collection']

# 获取集合中的数据并转换为DataFrame
data1 = list(collection1.find({}, {'_id': 0, 'AB': 1}))
data2 = list(collection2.find({}, {'_id': 0, 'AB': 1}))

# 转换为DataFrame
df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)

# 为DataFrame添加label列
df1['label'] = 'rejection'
df2['label'] = 'pass'

# 重命名AB列
df1.rename(columns={'AB': 'text'}, inplace=True)
df2.rename(columns={'AB': 'text'}, inplace=True)

# 去除空值
df1.dropna(inplace=True)
df2.dropna(inplace=True)

print("训练数据量:", len(df1) + len(df2))

# 打印df1和df2的行数
print("df1行数:", len(df1))
print("df2行数:", len(df2))

# 数据预处理
texts = df1['text'].tolist() + df2['text'].tolist()
labels = df1['label'].tolist() + df2['label'].tolist()

# 标签编码
label_map = {'rejection': 0, 'pass': 1}
labels = [label_map[label] for label in labels]

# 将文本转换为序列
tokenizer = Tokenizer()
tokenizer.fit_on_texts(texts)
sequences = tokenizer.texts_to_sequences(texts)

# 序列填充
max_length = 500  # 设定序列最大长度
sequences_padded = pad_sequences(sequences, maxlen=max_length, padding='post')

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(sequences_padded, labels, test_size=0.2, random_state=42)
print("sequences_padded shape:", sequences_padded.shape)
print("X_train shape:", X_train.shape)
print("X_test shape:", X_test.shape)

model_path = 'text_classification_model'

if os.path.exists(model_path):
    model = load_model(model_path)
    print("已加载现有模型。")
else:
    # 创建新模型
    model = Sequential([
        Embedding(input_dim=len(tokenizer.word_index) + 1, output_dim=100, input_length=max_length),
        Conv1D(128, 5, activation='relu'),
        GlobalMaxPooling1D(),
        Dropout(0.5),  # 添加dropout层，丢弃率为0.5
        Dense(64, activation='relu'),
        Dropout(0.5),  # 添加dropout层，丢弃率为0.5
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    print("已创建新模型。")

# 训练模型
model.fit(np.array(X_train), np.array(y_train), epochs=10, batch_size=8, validation_data=(np.array(X_test), np.array(y_test)))

# 评估模型
loss, accuracy = model.evaluate(np.array(X_test), np.array(y_test))
print(f'Test Accuracy: {accuracy:.4f}')

# 保存模型
tf.saved_model.save(model, model_path)
print("模型已保存。")
