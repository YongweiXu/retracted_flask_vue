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
    <div class="t_container">
    </div>

    <section class="section">
      <div class="row">
        <div class="column">
          <div class="t_box">
            <header class="t_title">
              <span>撤稿论文的被引用次数</span>
            </header>
            <div class="stats" id="myChart1" style="height: 500px;width: 100%"></div>
          </div>
        </div>

        <div class="column">
          <div class="t_box">
            <header class="t_title">
              <span>撤稿文献关键词及学科词云</span>
            </header>
            <div class="image-container">
              <img :src="'http://localhost:5000/wordcloud'" alt="" class="image">
            </div>
          </div>
        </div>

        <div class="column">
          <div class="t_box">
            <header class="t_title">
              <span>综合年度被撤稿关键词TOP15</span>
            </header>
            <div class="stats" id="myChart6" style="height: 500px;width: 100%"></div>
          </div>
        </div>
      </div>

    </section>
  </div>
</template>

<script>
import * as echarts from 'echarts';
import { Chart1, Chart6 } from './js/Chart';

export default {
  data() {
    return {};
  },
  mounted() {
    const myChart1 = echarts.init(document.getElementById('myChart1'));
    const myChart6 = echarts.init(document.getElementById('myChart6'));
    Chart1(myChart1);
    Chart6(myChart6);
  },
  methods: {
    // 返回操作
    goBack() {
      // 这里可以添加你的返回逻辑，比如返回上一页或者执行其他操作
      console.log('返回按钮被点击');
    }
  }
};
</script>

<style scoped>
* {
  margin: 0;
  padding: 0;
  box-sizing: border-box;
}

.container {
  background-size: cover;
  line-height: 1rem;
  display: flex;
  flex-direction: column;
  height: 100vh;
}

.row {
  display: flex;
  text-align: center;
  justify-content: center;
  margin-top: auto;
}

.column {
  flex: 1;
  display: flex;
  justify-content: center;
  align-items: center;
  text-align: center;
}


.section {
  flex: 1;
  display: flex;
  flex-direction: column;
}

.t_container {
  width: 100%;
}

.t_h_bg {
  width: 100%;
  height: 80px;
  line-height: 80px;
  background: url('static/images/t_header.png') no-repeat;
  background-size: 100% 100%;
  text-align: center;
}

.t_h_bg span {
  font-size: 40px;
  color: #CCE5FF;
  display: inline-block;
}

.t_box {
  margin: 0.125rem;
  width: 100%;
  height: 550px;
  background: url('static/images/t_bg.png') no-repeat;
  background-size: 100% 100%;
  position: relative;
}

.stats {
  display: flex;
  justify-content: center;
  align-items: center;
}


.t_title {
  text-align: center;
  font-size: 0.3rem;
  color: #CCE5FF;
  font-weight: bold;
  height: 0.4rem;
  line-height: 0.4rem;
  margin-top: 0.4rem;
}

.image-container {
  width: 100%;
  height: 100%;
  display: flex;
  justify-content: center;
  align-items: center;
}

.image {
  max-width: 96%; /* 让图片宽度充满容器 */
  max-height: 96%; /* 让图片高度充满容器 */
}

.button-row {
  display: flex;
  justify-content: center;
  align-items: center;
  height: 40px;
  margin-top: 10px;
}

.back-button {
  text-align: center;
  padding: 10px 20px;
  font-size: 16px;
  background-color: #007bff;
  color: #fff;
  border: none;
  border-radius: 5px;
  cursor: pointer;
}

</style>
