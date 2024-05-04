import { createApp } from 'vue';
import App from '../src/view/DemoComponent.vue';

import '../public/static/css/index.css';
import '../public/static/js/flexible.js';
import './util/rem';

const app = createApp(App);

// 设置默认页面标题
document.title = '退稿数据可视化大屏';

app.mount('#app');
