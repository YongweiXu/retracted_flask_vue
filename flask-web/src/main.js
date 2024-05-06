import App from './App.vue';
import '../public/static/css/index.css';
import '../public/static/js/flexible.js';
import './util/rem';
import router from "./router/index";
import { createApp } from 'vue';

// 设置默认页面标题
document.title = '退稿数据可视化大屏';

createApp({
  router, // 使用路由
  render: () => App // 直接返回 App 组件
}).mount('#app');
