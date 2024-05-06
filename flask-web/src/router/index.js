import {createRouter,createWebHistory} from "vue-router";

const router = createRouter({
    routes:[
        {
            path:"/",
            name:"index",
            component:()=>import("../view/IndexPage.vue"),
        },
        {
            path:"/homepage",
            name:"homepage",
            component:()=>import("../view/HomePage.vue"),
        },
        {
          path:"/predictionpage",
          name:"predictionpage",
          component:()=>import("../view/PredictionPage.vue"),
        },
        {
            path:"/query",
            name:"query",
            component:()=>import("../view/QueryPage.vue"),
        },
        {
            path:"/homepage2",
            name:"homepage2",
            component:()=>import("../view/HomePage2.vue"),
        }
    ],
    history:createWebHistory()
})
export default router;