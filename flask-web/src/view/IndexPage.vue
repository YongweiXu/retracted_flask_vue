<template>
  <div @click="redirectToHomePage" class="container">
    <img src="./static/images/1.gif" alt="" class="bg-image">
    <div class="content">
      <h1>生物医学撤稿大数据看板</h1>
      <p>
        <span>这是一个基于PubMed数据库的项目</span><br>
        <span>旨在为研究人员和学者提供一个集中展示</span><br>
        <span>生物医学领域撤稿文献信息的平台</span><br>
        <span>该项目的目标是为生物医学研究者提供一个</span><br>
        <span>高效、便捷的信息管理和分析工具</span><br>
        <span>帮助他们更好地了解领域动态、指导研究方向</span><br/>
      </p>
      <br><br>
      <h2>点击首页任意位置进入</h2>
      <h2>Biomedical Retraction Big Data Visualization Dashboard Design</h2>
    </div>
  </div>

<router-view></router-view>
</template>

<script>
export default {
  methods: {
    redirectToHomePage() {
      if (this.$router) {
        this.$router.push('/homepage');
      } else {
        console.error('Router is not available');
      }
}

  }
};
</script>

<style scoped>
.container {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100%;
}
.bg-image {
  width: 100%;
  height: 100%;
  object-fit: cover;
  position: absolute;
  top: 0;
  left: 0;
  z-index: -1;
}
.content {
  text-align: right;
  color: #ADD8E6;
  z-index: 1;
  position: relative;
  margin-left: auto;
  margin-right: 150px;
  transform: scale(1);
  text-shadow: -1px -1px 0 #3b677d,
                1px -1px 0 #3b677d,
                -1px 1px 0 #3b677d,
                1px 1px 0 #3b677d;
}
h1{
  font-size: 4em;
  margin-top: 100px;
  margin-bottom: 110px;
}
p {
  font-size: 1.6em;
}
</style>


