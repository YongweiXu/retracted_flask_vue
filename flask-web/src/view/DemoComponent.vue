<template>
  <div>
    <div class="head clearfix">
      <h1 class="pulll_left">生物医学期刊撤稿大数据可视化看板</h1>
      <div class="menu menu2 pulll_left">
        <ul>
          <li><a href="https://gitee.com/iGaoWei/big-data-view">导航标题</a></li>
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
            <div class="tit">模块标题</div>
            <div class="boxnav" style="height: 330px;">
              <div class="yqlist">
                <ul class="clearfix">
                  <li>
                    <div class="yq" id="yq">2634</div>
                    <span>原始数据量</span>
                  </li>
                  <li>
                    <div class="yq">567</div>
                    <span>有效数据量</span>
                  </li>
                  <li>
                    <div class="yq">56345</div>
                    <span>重新发表</span>
                  </li>
                  <li>
                    <div class="yq">721</div>
                    <span>数据展示(4)</span>
                  </li>
                </ul>
              </div>
            </div>
          </div>
          <div class="box">
            <div class="tit">模块标题</div>
            <div class="boxnav" style="height: 430px;">
              <div class="" style="height: 406px" id="echart2"></div>
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
            <div class="tit">模块标题</div>
            <div class="boxnav" style="height: 250px;" id="echart3"></div>
          </div>
        </li>
        <li style="width: 22%">
          <div class="box">
            <div class="tit">模块标题</div>
            <div class="boxnav" id="echart4" style="height: 200px;"></div>
          </div>
          <div class="box">
            <div class="tit">模块标题</div>
            <div class="boxnav" style="height: 250px;" id="echart5"></div>
          </div>
          <div class="box">
            <div class="tit">模块标题</div>
            <div class="boxnav" style="height: 250px;" id="echart6"></div>
          </div>
        </li>
      </ul>
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      currentTime: '',
      htmlSrc: 'static/js/test2.html',
      salesTotal: 2634 // Example data, replace with actual data
    };
  },
  mounted() {
    this.updateTime();
    setInterval(this.updateTime, 1000);
    window.addEventListener('message', this.getiframeMsg)
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

    // vue获取iframe传递过来的信息
    getiframeMsg(event) {
      const res = event.data;
      console.log(event)
      if (res.cmd == 'myIframe') {
        console.log(res)
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
      }, '*')
    },
    // 触发iframe中的方法
    iframeMethods() {
      this.$refs.iframe.contentWindow.triggerByVue('通过Vue触发iframe中的方法');
    }
  }
};
</script>

<style scoped>

</style>
