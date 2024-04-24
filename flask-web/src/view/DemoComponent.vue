<template>
  <div>
    <header>
      <h1>撤稿数据可视化大屏</h1>
      <div class="showtime">{{ currentTime }}</div>
    </header>
    <section>
      <div class="mainbox">
        <div class="column">
          <div class="panel bar">
            <h2>撤稿论文的被引用次数</h2>
            <div id="myChart1" style="width: 110%; height: 110%; text-align: center"></div>
            <div class="panel-footer"></div>
          </div>
          <div class="panel line">
            <h2>各出版社每年撤稿占比</h2>
            <div id="myChart2" style="width: 100%; height: 100%"></div>
            <div class="panel-footer"></div>
          </div>
          <div class="panel pie">
            <h2>撤稿出版社地区前十</h2>
            <div id="myChart3" style="width: 100%; height: 100%;"></div>
            <div class="panel-footer"></div>
          </div>
        </div>

        <div class="column">
          <div class="no">
            <div class="no-hd">
              <ul>
                <li>{{ frontendDemand }}</li>
                <li>{{ marketSupply }}</li>
              </ul>
            </div>
            <div class="no-bd">
              <ul>
                <li style="list-style-type:none;">总数据数</li>
                <li style="list-style-type:none;">有效数据数</li>
              </ul>
            </div>
          </div>
          <div class="map">
            <div class="chart"></div>

            <div class="map2"></div>
            <div>
              <img :src="'http://localhost:5000/wordcloud'" alt="" style="width: 500px; height: 300px; text-align: center">
            </div>
            <div class="map3"></div>
          </div>
        </div>
        <div class="column">
          <div class="panel bar2">
            <h2>年份撤稿关系</h2>
            <div id="myChart4" style="width: 110%; height: 110%; text-align: center"></div>
            <div class="panel-footer"></div>
          </div>
          <div class="panel line2">
            <h2>撤稿学科前十占比</h2>
            <div id="myChart5" style="width:110%; height: 110%; text-align: center"></div>
            <div class="panel-footer"></div>
          </div>
          <div class="panel pie2">
            <h2>综合年度被撤稿关键词TOP10</h2>
            <div id="myChart6" style="width:110%; height: 110%; text-align: center"></div>
            <div class="panel-footer"></div>
          </div>
        </div>
      </div>
    </section>
  </div>
</template>

<script>
import * as echarts from 'echarts';
import axios from 'axios';
import { Chart1 } from './static/js/Chart.js';
import { Chart2 } from "./static/js/Chart.js"; // 导入 drawBarChart 函数
import { Chart3 } from "./static/js/Chart.js"; // 导入 drawLineChart 函数
import { Chart4 } from "./static/js/Chart.js"; // 导入 drawPieChart 函数
import { Chart5 } from "./static/js/Chart.js"; // 导入 drawLineChart 函数
import { Chart6 } from "./static/js/Chart.js"; // 导入 drawPieChart 函数
export default {
  data() {
    return {
      currentTime: '',
      frontendDemand: '',
      marketSupply: ''
    };
  },
  mounted() {
    axios.get('http://localhost:5000/counts')
    .then(response => {
      // 更新原始数据量和清洗后数据量
      this.frontendDemand = response.data.original_count;
      this.marketSupply = response.data.cleaned_count;
    })
    .catch(error => {
      console.error('Error fetching data:', error);
    });
    this.updateTime();
    setInterval(this.updateTime, 1000);
        // 初始化图表
    const myChart1 = echarts.init(document.getElementById('myChart1'));
    const myChart2 = echarts.init(document.getElementById('myChart2'));
    const myChart3 = echarts.init(document.getElementById('myChart3'));
    const myChart4 = echarts.init(document.getElementById('myChart4'));
    const myChart5 = echarts.init(document.getElementById('myChart5'));
    const myChart6 = echarts.init(document.getElementById('myChart6'));

    Chart1(myChart1);
    Chart2(myChart2);
    Chart3(myChart3);
    Chart4(myChart4);
    Chart5(myChart5);
    Chart6(myChart6);
  },
  methods: {
    updateTime() {
      const now = new Date();
      const year = now.getFullYear();
      const month = now.getMonth() + 1;
      const day = now.getDate();
      const hour = now.getHours();
      const minute = now.getMinutes();
      const second = now.getSeconds();
      this.currentTime = `${year}年${month}月${day}日-${hour}时${minute}分${second}秒`;
    }
  }
};
</script>
<style>
#app {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
  margin-top: 60px;
}


</style>