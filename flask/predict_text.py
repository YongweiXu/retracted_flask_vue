import pickle
from flask import Flask, request, jsonify
import torch
from nltk import word_tokenize
from CharCNN import CharCNN  # 导入模型定义

app = Flask(__name__)

# 加载预训练的词嵌入矩阵
embeddings = torch.load('model/embeddings.pt')


# 加载模型配置
class Config:
    num_channels = 256
    seq_len = 1000
    linear_size = 256  # 与卷积层输出大小匹配
    output_size = 2
    dropout_keep = 0.5


# 创建模型实例
model = CharCNN(Config(), embeddings)

# 加载模型参数
model.load_state_dict(torch.load('model/char_cnn_model.pth'))

# 设置模型为评估模式
model.eval()

# 加载保存的最大序列长度
with open('model/max_seq_len.pkl', 'rb') as f:
    max_seq_len = pickle.load(f)


# 定义文本预处理函数
def preprocess_text(text, word_to_index):
    tokens = word_tokenize(text.lower())
    seq = [word_to_index.get(word, 0) for word in tokens]  # 使用0作为未知词的索引
    seq = seq[:max_seq_len]  # 截断或填充序列以匹配模型期望的长度
    seq += [0] * (max_seq_len - len(seq))  # 补0对齐
    return torch.LongTensor(seq).unsqueeze(0)  # 不需要将张量移动到 GPU


# 定义预测函数
def predict_text(text):
    # 加载保存的词汇表
    with open('model/word_to_index.pkl', 'rb') as f:
        word_to_index = pickle.load(f)

    with torch.no_grad():
        model.eval()  # 确保处于评估模式
        input_tensor = preprocess_text(text, word_to_index)
        output = model(input_tensor)
        probabilities = torch.softmax(output, dim=1)  # 对模型输出进行 softmax 处理

        # 获取类别标签列表
        class_labels = ['rejected', 'pass']

        # 获取预测的类别和概率
        _, predicted_class = torch.max(probabilities, dim=1)
        predicted_class_label = class_labels[predicted_class.item()]
        predicted_probability = probabilities.squeeze().tolist()

        # 构建结果字典
        result = {
            'probabilities': {
                class_labels[0]: predicted_probability[0],
                class_labels[1]: predicted_probability[1]
            }
        }

        # 返回预测结果
        return result