import torch
from torch import nn


class CharCNN(nn.Module):
    def __init__(self, config,  embeddings):
        super(CharCNN, self).__init__()
        self.config = config
        embed_size = embeddings.shape[1]  # 使用提供的嵌入维度

        # 将嵌入转换为PyTorch张量
        embeddings = torch.tensor(embeddings)

        # 嵌入层
        self.embeddings = nn.Embedding.from_pretrained(embeddings, freeze=True)

        # 卷积层
        self.conv1 = nn.Sequential(
            nn.Conv1d(in_channels=embed_size, out_channels=self.config.num_channels, kernel_size=7),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=3)
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(in_channels=self.config.num_channels, out_channels=self.config.num_channels, kernel_size=7),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=3)
        )

        # 线性层
        # 调整输入大小以确保与卷积层输出的大小匹配
        self.flatten = nn.Flatten()
        self.linear1 = nn.Sequential(
            nn.Linear(48384, self.config.linear_size),
            nn.ReLU(),
            nn.Dropout(self.config.dropout_keep)
        )
        self.linear2 = nn.Linear(self.config.linear_size, self.config.output_size)

    def forward(self, x):
        embedded_sent = self.embeddings(x).permute(0, 2, 1)
        conv_out1 = self.conv1(embedded_sent.float())
        conv_out2 = self.conv2(conv_out1)
        conv_out = self.flatten(conv_out2)
        linear_output = self.linear1(conv_out)
        linear_output = self.linear2(linear_output)
        return linear_output
