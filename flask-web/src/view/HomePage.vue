<template>
  <div>
    <div class="head clearfix">
      <h1 class="pulll_left">生物医学期刊撤稿大数据可视化看板</h1>
      <div class="menu menu2 pulll_left">
        <ul>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">首页</a></li>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">导航标题样式</a></li>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">导航标题</a></li>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">导航标题</a></li>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">导航标题</a></li>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">导航标题</a></li>
        </ul>
      </div>
      <div class="time">{{ currentTime }}</div>
    </div>
    <div class="mainbox">
      <ul class="clearfix nav1">
        <li style="width: 22%">
          <div class="box">
            <div class="tit">数据信息一览</div>
            <div class="boxnav" style="height: 330px;">
              <div class="yqlist">
                <ul class="clearfix">
                  <li>
                    <div class="yq" id="yq">{{Original + Pass}}</div>
                    <span>总数据量</span>
                  </li>
                  <li>
                    <div class="yq">{{Pass}}</div>
                    <span>未被退稿数据量数据量</span>
                  </li>
                  <li>
                    <div class="yq">{{Original}}</div>
                    <span>被退稿数据量</span>
                  </li>
                  <li>
                    <div class="yq">57241</div>
                    <span>模型训练数据量</span>
                  </li>
                </ul>
              </div>
            </div>
          </div>
          <div class="box">
            <div class="tit">退稿作者top15</div>
            <div  class="boxnav" style="height: 430px">
              <div id="myRanking" style="height: 430px "></div>
            </div>
          </div>
        </li>
        <li style="width: 56%">
          <div class="box">
            <div class="boxnav mapc" style="height: 550px; position: relative">
              <div class="map" id="map">
                <iframe ref="iframe"
                        :src="htmlSrc"
                        width="100%"
                        height="100%"
                        frameborder="0">
                </iframe>
              </div>
            </div>
          </div>
          <div class="box">
            <div class="tit">撤稿信息轮播图</div>
            <div class="boxnav" style="height: 250px;" id="echart3">
              <iframe ref="iframe"
                      :src="htmlSrc2"
                      width="100%"
                      height="100%"
                      frameborder="0"
                      style="height: 100%;">
              </iframe>
            </div>
          </div>
        </li>
        <li style="width: 22%">
          <div class="box">
            <div class="tit">关键词共现矩阵</div>
            <div  class="boxnav" style="height: 200px">
              <div id="myMatrix" style="height: 140% ;width: 100%;" ></div>
            </div>
          </div>
          <div class="box">
            <div class="tit">年份撤稿关系</div>
            <div  class="boxnav" style="height: 250px">
              <div  id="myChart4" style="height: 100% ;width: 100%;"></div>
            </div>
          </div>
          <div class="box">
            <div class="tit">期刊年份退稿关系</div>
            <div  class="boxnav" style="height: 250px">
              <div  id="myChart2" style="height: 100% ;width: 100%;"></div>
            </div>
          </div>
        </li>
      </ul>
    </div>
  </div>
</template>

<script>
import * as echarts from 'echarts';
import { Ranking } from './js/newChart';
import { matrix } from './js/newChart';
import { Chart4 } from './js/Chart';
import { Chart2 } from './js/Chart';
import axios from "axios";
export default {
  data() {
    return {
      currentTime: '',
      htmlSrc: 'static/js/test2.html',
      htmlSrc2: 'static/js/test.html',
      salesTotal: 2634, // Example data, replace with actual data
      Original: '',
      Pass: '',
    };
  },
  mounted() {
    axios.get("http://localhost:5000/counts")
        .then(response =>{
          this.Original = response.data.re_count;
          this.Pass = response.data.pass_count;
        })
    this.updateTime();
    setInterval(this.updateTime, 1000);
    window.addEventListener('message', this.getiframeMsg);
    const myRanking = echarts.init(document.getElementById('myRanking'))
    const myChart4 = echarts.init(document.getElementById('myChart4'))
    const myChart2 = echarts.init(document.getElementById('myChart2'))
    // 在组件挂载后调用 matrix 函数来初始化 echarts 实例
    this.initMatrixChart();
    Ranking(myRanking);
    Chart4(myChart4);
    Chart2(myChart2);
  },
  methods: {
    updateTime() {
      const now = new Date();
      const y = now.getFullYear();
      const mt = now.getMonth() + 1;
      const day = now.getDate();
      const h = now.getHours();
      const m = now.getMinutes();
      const s = now.getSeconds();
      this.currentTime = `${y}/${mt}/${day} ${h}:${m}:${s}`;
    },
    // 初始化 echarts 实例的方法
    initMatrixChart() {
      const myMatrix = echarts.init(document.getElementById('myMatrix'));

      // 假设 matrix 函数接受一个参数来初始化 echarts 实例
      matrix(myMatrix);
    },
    // vue获取iframe传递过来的信息
    getiframeMsg(event) {
      const res = event.data;
      console.log(event);
      if (res.cmd == 'myIframe') {
        console.log(res);
      }
    },
    // vue向iframe传递信息
    vueSendMsg() {
      const iframeWindow = this.$refs.iframe.contentWindow;
      iframeWindow.postMessage({
        cmd: 'myVue',
        params: {
          info: 'Vue向iframe传递的消息',
        }
      }, '*');
    },
    // 触发iframe中的方法
    iframeMethods() {
      this.$refs.iframe.contentWindow.triggerByVue('通过Vue触发iframe中的方法');
    }
  }
};
</script>

<style scoped>
/* 这里可以添加一些局部样式 */
</style>
