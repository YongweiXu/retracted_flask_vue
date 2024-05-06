import { createApp } from 'vue';
import App from './App.vue';
import router from "./router";

import '../public/static/css/index.css';
import '../public/static/js/flexible.js';
import './util/rem';

// 设置默认页面标题
document.title = '退稿数据可视化大屏';

const app = createApp(App);

app.use(router);

app.mount('#app');
