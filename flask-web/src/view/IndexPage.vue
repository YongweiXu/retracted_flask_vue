<template>
  <div @click="redirectToHomePage" class="container">
    <img src="./static/images/bg2.png" alt="" class="bg-image">
    <div class="content">
      <h1>生物医学撤稿大数据看板</h1>
      <p>
        <span>"生物医学撤稿大数据看板"是一个基于PubMed数据库的项目，旨在为研究人员和学者提供一个集中</span><br/>
        <span>展示生物医学领域撤稿文献信息的平台。该项目的目标是为生物医学研究者提供一个高效、</span><br/>
        <span>便捷的信息管理和分析工具，帮助他们更好地了解领域动态、指导研究方向，</span><br/>
        <span>并促进学术交流和合作。通过该平台，用户可以从首页直接进</span><br/>
        <span>入看板页面，即可对生物医学撤稿文献信息进行了解</span>
      </p>
    </div>
  </div>
</template>

<script>
export default {
  methods: {
    redirectToHomePage() {
      // 使用编程式导航跳转到HomePage页面
      this.$router.push({ name: 'HomePage' });
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
  height: 100vh;
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
  text-align: center;
  color: white;
  z-index: 1;
  position: relative;
}

h1 {
  font-size: 4em;
  margin-top: 200px;
  margin-bottom: 100px;
}

p {
  font-size: 1.6em;
}
</style>
