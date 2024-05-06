import { createRouter, createWebHistory } from 'vue-router';
import Index from "@/view/IndexPage.vue";
import HomePage from "@/view/HomePage.vue";
import Prediction from "@/view/PredictionPage.vue";

const router = createRouter({
  history: createWebHistory(), // 使用 createWebHistory 创建历史模式
  routes: [
    {
      path: '/',
      name: 'Index',
      component: Index
    },
    {
      path: '/homepage',
      name: 'HomePage',
      component: HomePage
    },
    {
      path: '/prediction',
      name: 'Prediction',
      component: Prediction
    }
  ]
});

export default router;

