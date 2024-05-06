<template>
  <div id="app">
    <router-view></router-view>
  </div>
</template>

<script>
import { defineComponent } from 'vue'
import router from './router'

export default defineComponent({
  name: 'App',
  router
})

</script>

