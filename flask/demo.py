from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.preprocessing.text import Tokenizer

# 加载模型
model = load_model('text_classification_model')

# 加载文本预处理器
tokenizer = Tokenizer()
# 这里定义并传递texts变量，该变量应该包含用于训练模型的文本数据
# texts = [
#     "This is the first text for training the model.",
#     "This is the second text for training the model.",
# ]
# tokenizer.fit_on_texts(texts)
max_length = 500  # 与训练时一致

# 定义函数来对文本进行分类
def classify_text(text):
    # 将文本转换为序列
    sequence = tokenizer.texts_to_sequences([text])
    # 序列填充
    sequence_padded = pad_sequences(sequence, maxlen=max_length, padding='post')
    # 进行预测
    prediction = model.predict(sequence_padded)[0][0]
    # 返回概率值
    positive_probability = prediction
    negative_probability = 1 - prediction
    return {'positive_probability': positive_probability, 'negative_probability': negative_probability}

# 使用示例
text_to_classify = "This is a sample text to classify."
result = classify_text(text_to_classify)
print(result)
