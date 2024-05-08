import torch
from torch import nn
from gensim.models import Word2Vec
import numpy as np
import pandas as pd
from pymongo import MongoClient
from sklearn.model_selection import train_test_split
from nltk.tokenize import word_tokenize
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
from CharCNN import CharCNN
import pickle

#检查是否有GPU可用
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def generate_embedding_matrix(word_to_index):
    # 加载预训练的 Word2Vec 模型
    word2vec_model = Word2Vec.load("model/word2vec.model")

    embedding_dim = word2vec_model.vector_size
    vocab_size = len(word_to_index)
    embedding_matrix = np.zeros((vocab_size, embedding_dim))

    for word, index in word_to_index.items():
        if word in word2vec_model.wv:
            embedding_matrix[index] = word2vec_model.wv[word]

    return embedding_matrix


#连接MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client['pubmed']

#从集合中获取数据并转换为DataFrame
data1 = list(db['pubmed_collection'].find({}, {'AB': 1}))
data2 = list(db['pubmed_review_collection'].find({}, {'AB': 1}))

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df1.drop('_id', axis=1, inplace=True)
df2.drop('_id', axis=1, inplace=True)
df1.dropna(subset=['AB'], inplace=True)
df2.dropna(subset=['AB'], inplace=True)
df1.rename(columns={'AB': 'text'}, inplace=True)
df2.rename(columns={'AB': 'text'}, inplace=True)
df1['label'] = 'rejection'
df2['label'] = 'pass'

df = pd.concat([df1, df2], ignore_index=True)

#对文本进行分词
texts = [word_tokenize(text.lower()) for text in df['text']]

#构建词汇表和词汇索引映射
word_set = set(word for text in texts for word in text)
word_to_index = {word: idx for idx, word in enumerate(word_set)}
vocab_size = len(word_to_index)

#构建词嵌入矩阵
embeddings = generate_embedding_matrix(word_to_index)

torch.save(embeddings, 'model/embeddings.pt')

# 数据预处理
sequences = [[word_to_index[word] for word in text] for text in texts]
max_seq_len = max(len(seq) for seq in sequences)
with open('model/max_seq_len.pkl', 'wb') as f:
    pickle.dump(max_seq_len, f)

X = np.array([seq + [0] * (max_seq_len - len(seq)) for seq in sequences])
y = np.array(df['label'].map({'rejection': 0, 'pass': 1}))

#划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

#调整批量大小
batch_size = 4

#创建数据加载器
train_dataset = TensorDataset(torch.LongTensor(X_train), torch.LongTensor(y_train))
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

#模型配置
class Config:
    num_channels = 256
    seq_len = max_seq_len
    linear_size = 256  # 与卷积层输出大小匹配
    output_size = 2
    dropout_keep = 0.5

#创建模型实例并移到GPU
model = CharCNN(Config(), embeddings)
model.to(device)

#损失函数
criterion = nn.CrossEntropyLoss()

#选择优化器
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

#学习率调度器
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.1)

#训练模型
num_epochs = 20
for epoch in range(num_epochs):
    model.train()
    loop = tqdm(train_loader, leave=True)
    for inputs, labels in loop:
        optimizer.zero_grad()
        inputs, labels = inputs.to(device), labels.to(device)  # 将数据移动到设备上
        inputs = inputs.long()  # 将输入张量的数据类型转换为 Long 类型
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        loop.set_description(f'Epoch [{epoch + 1}/{num_epochs}]')
        loop.set_postfix(loss=loss.item())
    #调用学习率调度器
    scheduler.step()

torch.save(model.state_dict(), 'model/char_cnn_model.pth')
