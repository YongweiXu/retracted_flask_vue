// static/js/Chart.js
import * as echarts from 'echarts';
import axios from 'axios';

export function Chart1(myChart1) {
  fetch('http://localhost:5000/TIZ9')
    .then(response => response.json())
    .then(data => {
      const titles = Object.values(data);
      const referenceCounts = Object.keys(data);

      const top20Titles = titles.slice(0, 20);
      const top20ReferenceCounts = referenceCounts.slice(0, 20);

      let option = {
        color: ['#3398DB'],
        xAxis: {
          show: false,
          type: 'category',
          data: top20Titles,
          axisTick: {
            alignWithLabel: true
          }
        },
        yAxis: {
          type: 'value',
          axisLine: {
            lineStyle: {
              color: '#fff'
            }
          },
          axisLabel: {
            color: '#fff'
          }
        },
        series: [{
          data: top20ReferenceCounts.map(count => parseInt(count)),
          type: 'bar',
          barWidth: 'auto',
          label: {
            show: true,
            position: 'top',
            color: '#fff',
            emphasis: {
              textBorderColor: 'transparent'
            }
          }
        }],
        tooltip: {
          position: function (pos, params, el, elRect, size) {
            var obj = { top: 10 };
            obj[['left', 'right'][+(pos[0] < size.viewSize[0] / 2)]] = 10;
            return obj;
          },
          formatter: function(params) {
            const title = params.name;
            const referenceCount = params.value;
            return `${title}: ${referenceCount}`;
          }
        },
        dataZoom: [{
          type: 'slider',
          show: false,
          start: 0,
          end: 100
        }, {
          type: 'inside',
          start: 0,
          end: 100
        }]
      };

      myChart1.setOption(option);

      const zoomSize = 6;
      myChart1.on('click', function (params) {
        const dataAxis = top20Titles;
        myChart1.dispatchAction({
          type: 'dataZoom',
          startValue: dataAxis[Math.max(params.dataIndex - zoomSize / 2, 0)],
          endValue: dataAxis[Math.min(params.dataIndex + zoomSize / 2, top20Titles.length - 1)]
        });
      });
    })
    .catch(error => {
      console.error('Error fetching data:', error);
    });
}




export function Chart2(myChart2) {
    // 从后端API获取数据
    fetch('http://localhost:5000/PU')
        .then(response => response.json())
        .then(data => {
            const sourceData = data.data.slice(1).map(row => row.split(','));
            const productNames = sourceData.map(row => row[0]);
            const years = data.data[0].split(',').slice(1); // 不倒置年份

            const seriesData = sourceData.map((row, index) => ({
                type: 'line',
                smooth: true,
                seriesLayoutBy: 'row',
                emphasis: { focus: 'series' },
                showSymbol: true,
                name: productNames[index],
                data: row.slice(1).map(value => parseFloat(value)) // 不倒置数据
            }));

            const option = {
                legend: false, // 隐藏图例
                tooltip: {
                    trigger: 'axis',
                    formatter: function (params) {
                        const data = params[0];
                        return data.seriesName + '<br/>' + data.name + ' : ' + data.value;
                    }
                },
                dataset: {
                    source: [
                        ['product', ...years], // 年份在第一行
                        ...sourceData // 不需要倒置数据
                    ]
                },
                xAxis: {
                    type: 'category',
                    data: years, // 设置x轴的数据为年份
                    axisPointer: {
                        type: 'shadow' // 设置悬停时的指示器类型为阴影
                    },
                    axisLabel: {
                        color: '#fff' // 设置 x 轴文字颜色为白色
                    }
                },
                yAxis: {
                    gridIndex: 0,
                    axisLabel: {
                        show: false // 隐藏 y 轴文字
                    }
                },
                grid: { top: '55%' },
                dataZoom: [ // 配置dataZoom组件
                    {
                        type: 'inside', // 内置的数据区域缩放组件
                        xAxisIndex: 0, // 表示作用在第一个x轴上
                        start: 0, // 缩放起始位置
                        end: 100 // 缩放结束位置
                    }
                ],
                series: seriesData.concat({
                    type: 'pie',
                    id: 'pie',
                    radius: '30%',
                    center: ['50%', '25%'],
                    emphasis: { focus: 'self' },
                    label: {
                        normal: {
                            show: false // 正常状态下不显示标签
                        },
                        emphasis: {
                            show: true, // 悬停时显示标签
                            formatter: '{b}: {@2012} ({d}%)'
                        }
                    },
                    encode: { itemName: 'product', value: '2012', tooltip: '2012' }
                })
            };

            myChart2.on('updateAxisPointer', function (event) {
                const xAxisInfo = event.axesInfo[0];
                if (xAxisInfo) {
                    const dimension = xAxisInfo.value + 1;
                    myChart2.setOption({
                        series: {
                            id: 'pie',
                            label: { formatter: '{b}: {@' + dimension + '} ({d}%)' },
                            encode: { value: dimension, tooltip: dimension }
                        }
                    });
                }
            });

            myChart2.setOption(option);
        })
        .catch(error => {
            console.error('获取数据时出错：', error);
        });
}



export function Chart3(myChart3) {
  axios.get('http://localhost:5000/PI')
    .then(response => {
      const data = response.data;
      renderChart(myChart3, data); // 将数据传递给渲染函数
    })
    .catch(error => {
      console.error('获取数据时出错:', error);
    });

  function renderChart(myChart, data) {
    const option = {
      polar: {
        radius: [30, '80%']
      },
      angleAxis: {
        startAngle: 75
      },
      radiusAxis: {
        type: 'category',
        data: Object.keys(data) // 使用数据的键作为角度轴的数据
      },
      tooltip: {
        show: true, // 显示提示框
        trigger: 'axis', // 设置触发类型为 axis，鼠标悬停时触发
        formatter: function(params) {
          const value = params[0].value; // 获取数据值
          const result = value * 8; // 计算乘以 8 的结果
          return params[0].name + ': ' + result; // 返回格式化后的字符串
        }
      },
      series: [
        {
          type: 'bar',
          data: Object.entries(data).map(([city, value]) => ({
            value: value / 8, // 将高度降低为原来的1/10
            label: {
              show: false, // 不显示标签
              position: 'top', // 设置标签位置为顶部
              formatter: function() {
                return `{b}: ${city}: ${value}`;
              }
            }
          })),
          coordinateSystem: 'polar'
        }
      ]
    };

    myChart.setOption(option);
  }
}

// 在 Vue 组件中使用 Chart3 函数
export default {
  mounted() {
    const myChart3 = echarts.init(document.getElementById('myChart3'));
    Chart3(myChart3);
  }
};






export function Chart4(myChart4) {
  function fetchDataFromBackend() {
    return fetch('http://localhost:5000/PY')
      .then(response => {
        if (!response.ok) {
          throw new Error('Network response was not ok');
        }
        return response.json();
      })
      .catch(error => {
        console.error('There was a problem with your fetch operation:', error);
      });
  }

  const option = {
    tooltip: {
      trigger: 'axis',
      axisPointer: {
        type: 'cross'
      }
    },
    toolbox: {
      show: true,
      feature: {
        saveAsImage: {}
      }
    },
    xAxis: {
      type: 'category',
      boundaryGap: false,
      data: [], // 根据后端数据动态设置
      axisLabel: {
        color: '#fff' // 设置 x 轴文字颜色为白色
      }
    },
    yAxis: {
      type: 'value',
      name: 'Counts', // y轴名称
      axisPointer: {
        snap: true
      },
      axisLabel: {
        show: false // 隐藏 y 轴文字
      }
    },
    dataZoom: [{ // 添加dataZoom组件配置
      type: 'inside', // 内置缩放
      start: 0,
      end: 100
    }],
    series: [
      {
        name: 'Rejection Counts',
        type: 'line',
        smooth: true,
        data: [], // 根据后端数据动态设置
        animation: false // 禁用折线动画效果
      }
    ]
  };

  fetchDataFromBackend().then(data => {
    console.log('Received data from backend:', data); // 添加这行来打印数据，以便检查数据是否正确
    if (data && data.length > 0) {
      const sortedData = data.sort((a, b) => a['年份'] - b['年份']);
      const xAxisData = sortedData.map(item => item['年份']);
      const seriesData = sortedData.map(item => item['计数']);
      option.xAxis.data = xAxisData;
      option.series[0].data = seriesData.map((count, index) => ({
        value: count,
        itemStyle: index > 0 && count > seriesData[index - 1] ? { color: 'green' } : { color: 'red' }
      }));

      console.log('Prepared option:', option); // 添加这行来打印准备好的选项，以便检查选项是否正确
      myChart4.setOption(option);
    } else {
      console.error('Invalid data received from backend:', data);
    }
  });
}









export function Chart5(myChart5, myChart6) {
    // 发起 HTTP 请求从后端获取数据
    fetch('http://localhost:5000/WC')
        .then(response => response.json())
        .then(data => {
            // 处理数据，确保 count 属性值是数字类型
            data.forEach(item => {
                item.count = parseInt(item.count); // 将 count 属性值转换为整数类型
            });

            // 筛选出 "other" 数据项并将其值缩小
            const otherItem = data.find(item => item.name === "other");
            if (otherItem) {
                otherItem.originalCount = otherItem.count; // 保留 "other" 的原始值
                otherItem.count = Math.round(otherItem.count * 0.08); // 将 "other" 数据项的值缩小
            }

            // 对数据进行排序
            data.sort((a, b) => b.count - a.count);

            // 定义数据项
            var seriesData = data.map(item => ({
                name: item.name,
                value: item.count
            }));

            // 设置图表配置
            var option = {
                toolbox: {
                    show: true,
                    feature: {
                        mark: { show: true },
                        dataView: { show: true, readOnly: false },
                        restore: { show: true },
                        saveAsImage: { show: true }
                    }
                },
                tooltip: {
                    trigger: 'item',
                    formatter: function(params) {
                        if (params.name === "other") {
                            return params.seriesName + "<br/>" + params.name + " : " + otherItem.originalCount + " (原始值)";
                        } else {
                            return params.seriesName + "<br/>" + params.name + " : " + params.value + " (" + params.percent + "%)";
                        }
                    }
                },
                series: [{
                    name: 'Nightingale Chart',
                    type: 'pie',
                    radius: '55%',
                    center: ['50%', '50%'],
                    roseType: 'area',
                    label: {
                        show: true
                    },
                    emphasis: {
                        label: {
                            show: true
                        }
                    },
                    data: seriesData
                }]
            };

            // 设置图表配置并渲染图表
            myChart5.setOption(option);

            // 创建另一个图表容器，并渲染 "other" 数据项
            var option2 = {
                toolbox: {
                    show: true,
                    feature: {
                        mark: { show: true },
                        dataView: { show: true, readOnly: false },
                        restore: { show: true },
                        saveAsImage: { show: true }
                    }
                },
                tooltip: {
                    trigger: 'item',
                    formatter: '{a} <br/>{b} : {c} ({d}%)'
                },
                series: [{
                    type: 'pie',
                    radius: '55%',
                    center: ['50%', '50%'],
                    roseType: 'area',
                    label: {
                        show: true
                    },
                    emphasis: {
                        label: {
                            show: true
                        }
                    },
                    data: [{
                        name: otherItem.name,
                        value: otherItem.count
                    }]
                }]
            };

            // 渲染第二个图表
            myChart6.setOption(option2);

        })
        .catch(error => {
            console.error('Error fetching data:', error);
        });
}


export function Chart6(myChart) {
  // 发送请求获取关键词数据
  fetch('http://localhost:5000/DE')
    .then(response => {
      if (!response.ok) {
        throw new Error('网络响应异常');
      }
      return response.json();
    })
    .then(data => {
      const keywords = Object.keys(data);

      const seriesList = [];
      const years = [];

      for (let year = 1990; year <= 2025; year++) {
        years.push(year.toString());
      }

      years.forEach(year => {
        const dataItem = { Year: year };
        keywords.forEach(keyword => {
          dataItem[keyword] = data[keyword][year] || 0;
        });
        seriesList.push(dataItem);
      });

      const option = {
        animationDuration: 10000,
        dataset: [{
          source: seriesList
        }],
        tooltip: {
          order: 'valueDesc',
          trigger: 'axis'
        },
        xAxis: {
          type: 'category',
          data: years,
          nameLocation: 'middle',
          boundaryGap: false,
          axisLine: {
            lineStyle: {
              color: '#fff'
            }
          },
          axisLabel: {
            color: '#fff'
          }
        },
        yAxis: {
          name: '数量',
          axisLine: {
            lineStyle: {
              color: '#fff'
            }
          },
          axisLabel: {
            color: '#fff'
          }
        },
        grid: {
          right: 140
        },
        dataZoom: [{
          type: 'inside',
          start: 0,
          end: 100
        }],
        series: keywords.map(keyword => ({
          type: 'line',
          showSymbol: false,
          name: keyword,
          endLabel: {
            show: false,
            formatter: function (params) {
              return params.value[3] + ': ' + params.value[0];
            }
          },
          labelLayout: {
            moveOverlap: 'shiftY'
          },
          emphasis: {
            focus: 'series'
          },
          encode: {
            x: 'Year',
            y: keyword,
            label: {
              show: false,
              position: 'top',
              formatter: function(params) {
                return params.value[1];
              }
            },
            tooltip: [keyword]
          }
        }))
      };

      myChart.setOption(option);
    })
    .catch(error => {
      console.error('发生了问题:', error);
    });
}
