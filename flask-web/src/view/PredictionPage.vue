<template>
  <div class="head clearfix">
      <h1 class="pulll_left">生物医学期刊撤稿大数据可视化看板</h1>
      <div class="menu menu2 pulll_left">
        <ul>
          <li><router-link to="/">首页</router-link></li>
          <li><router-link to="/homepage">数据看板</router-link></li>
          <li><router-link to="/homepage2">影响力看板</router-link></li>
          <li><router-link to="/predictionpage">摘要预测</router-link></li>
          <li><router-link to="/query">数据查询</router-link></li>
        </ul>
      </div>
      <div class="time">{{ currentTime }}</div>
  </div>
  <div class="container">
    <span style="margin-bottom: 10px;">本预测依赖本次收录数据训练的CharCNN模型（即字符级卷积神经网络）进行处理，结果仅供参考</span>
    <textarea v-model="textInput" placeholder="请输入要预测的文献的摘要" rows="15" cols="100"></textarea>
    <button @click="predict" style="width: 80px;height:60px;border-radius: 10px;">预测</button>
    <pie-chart v-if="showChartDialog" :pass-probability="passProbability" :rejected-probability="rejectedProbability" />
  </div>
</template>

<script>
import axios from 'axios';
import PieChart from './PieChart.vue';

export default {
  data() {
    return {
      textInput: '',
      predictionResult: null,
      passProbability: 0,
      rejectedProbability: 0,
      showChartDialog: false
    };
  },
  components: {
    PieChart
  },
  methods: {
    async predict() {
      // 发送POST请求，提交文本框内容
      try {
        const response = await axios.post('http://localhost:5000/predict', {text: this.textInput});
        // 清空文本框
        this.textInput = '';
        // 获取预测结果
        this.predictionResult = response.data;
        // 更新通过概率和拒绝概率
        this.passProbability = this.predictionResult.probabilities.pass;
        this.rejectedProbability = this.predictionResult.probabilities.rejected;

        // 直接显示图表
        this.showChartDialog = true;
      } catch (error) {
        console.error('预测请求出错：', error);
        // 处理错误
      }
    }
  }
};
</script>
<style>
.container {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100vh; /* 使内容垂直居中 */
}

textarea {
  border-radius: 10px;
  margin-bottom: 10px;
}
</style>