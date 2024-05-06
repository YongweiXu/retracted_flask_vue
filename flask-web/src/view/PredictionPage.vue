<template>
  <div>
    <textarea v-model="textInput" placeholder="请输入文本" rows="10" cols="50"></textarea>
    <button @click="predict">预测</button>
    <button @click="showChart">显示图表</button>
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
        const response = await axios.post('http://localhost:5000/predict', { term: this.textInput });
        // 清空文本框
        this.textInput = '';
        // 获取预测结果
        this.predictionResult = response.data;
        // 更新通过概率和拒绝概率
        this.passProbability = this.predictionResult.probabilities.pass;
        this.rejectedProbability = this.predictionResult.probabilities.rejected;
      } catch (error) {
        console.error('预测请求出错：', error);
        // 处理错误
      }
    },
    showChart() {
      this.showChartDialog = true;
    }
  }
};
</script>
