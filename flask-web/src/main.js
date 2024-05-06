import { createApp } from 'vue';
import App from './App.vue';
import '../public/static/css/index.css';
import '../public/static/js/flexible.js';
import './util/rem';
import router from "./router/index";

// 设置默认页面标题
document.title = '退稿数据可视化大屏';

const app = createApp(App); // 创建 Vue 应用程序实例
app.use(router); // 使用路由器
app.mount('#app'); // 挂载应用程序到 #app 元素上
