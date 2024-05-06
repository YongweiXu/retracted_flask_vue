<template>
  <div ref="pieChart" style="width: 400px; height: 400px;"></div>
</template>

<script>
import * as echarts from 'echarts';

export default {
  props: ['passProbability', 'rejectedProbability'],
  mounted() {
    this.drawPieChart();
  },
  methods: {
    drawPieChart() {
      const pieChart = echarts.init(this.$refs.pieChart);

      const option = {
        title: {
          text: '预测结果',
          subtext: '通过概率 vs 拒绝概率',
          left: 'center',
          subtextStyle:{
            color: '#efe5e3' // 标题颜色为黑色
          },
          textStyle: {
            color: '#ee6b3f' // 标题颜色为黑色
    }
        },
        tooltip: {
          trigger: 'item',
          formatter: '{a} <br/>{b} : {c} ({d}%)'
        },
        legend: {
          orient: 'vertical',
          left: 'left',
          data: ['通过概率', '拒绝概率'],
          textStyle: {
            color: '#efe5e3' // 标题颜色为黑色
    }
        },
        series: [
          {
            name: '预测结果',
            type: 'pie',
            radius: '55%',
            center: ['50%', '60%'],
            textStyle: {
            color: '#efe5e3' // 标题颜色为黑色
            },
            data: [
              { value: this.passProbability, name: '通过概率', itemStyle: { color: '#67C23A' } }, // 通过概率的颜色为绿色
              { value: this.rejectedProbability, name: '拒绝概率', itemStyle: { color: '#F56C6C' } } // 拒绝概率的颜色为红色
            ],
            emphasis: {
              itemStyle: {
                shadowBlur: 10,
                shadowOffsetX: 0,
                shadowColor: 'rgba(0, 0, 0, 0.5)'
              }
            }
          }
        ]
      };

      pieChart.setOption(option);
    }
  }
};
</script>

<style scoped>
/* 可以添加一些样式来调整图表的外观 */
</style>
