<template>
  <div @click="redirectToHomePage">
    <img src="././static/images/bg2.png" alt="" class="bg-image">
    <div>
      <h1>生物医学撤稿大数据看板</h1>
    </div>
    <div>
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
body {
  margin: 0;
  padding: 0;
  height: 100vh;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  cursor: pointer;
  position: relative;
}

.bg-image {
  position: fixed;
  width: 100%;
  height: 100vh;
  opacity: 0.84;
}

h1 {
  color: white;
  font-family: Arial, sans-serif;
  font-size: 4em;
  text-align: center;
  margin-top: 200px;
  margin-bottom: 100px;
  position: relative;
  z-index: 1;
}

div {
  position: relative;
  left: 50px;
  bottom: 200px;
  text-shadow: 1px 1px 2px;
  opacity: 1;
  z-index: 1;
}

p {
  color: white;
  font-family: main_Arial, sans-serif;
  font-size: 1.6em;
  text-align: center;
  margin-top: auto;
  margin-bottom: auto;
}
</style>
