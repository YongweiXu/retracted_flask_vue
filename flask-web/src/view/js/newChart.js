export function matrix(myChart) {
    // 从指定的 URL 获取数据
    fetch('http://127.0.0.1:5000/Coc')
    .then(response => response.json()) // 解析 JSON 响应
    .then(data => {
        // 从响应中提取数据
        const matrix = data.matrix;
        const term1 = data.term1;
        const term2 = data.term2;
        const matrix1 = [];
        for (let i = 0; i < term1.length; i++) {
            for (let j = 0; j < term2.length; j++) {
                // 注意调换了 j 和 i 的位置，以便正确显示在热力图中
                matrix1.push([term2[j], term1[i], matrix[i][j]]);
            }
        }

        // 使用提取的数据渲染图表
        const option = {
            tooltip: {
                position: 'top',
                formatter: function(params) {
                    return  params.data[1] + '<br>' + // 修改这里以适应调换后的数据格式
                            params.data[0] + '<br>' + // 修改这里以适应调换后的数据格式
                            params.data[2];
                }
            },
            xAxis: {
                type: 'category', // 设置 x 轴为分类轴
                data: term2, // 使用 term2 作为 x 轴数据
                show: false, // 关闭 x 轴显示
            },
            yAxis: {
                type: 'category', // 设置 y 轴为分类轴
                data: term1, // 使用 term1 作为 y 轴数据
                show: false, // 关闭 y 轴显示
            },
            visualMap: {
                min: 0,
                max: 5000, // 加大颜色取值范围
                calculable: false,
                orient: 'horizontal',
                left: 'center',
                bottom: '12%'
            },
            series: [{
                name: '共现矩阵',
                type: 'heatmap',
                data: matrix1,
                label: {
                    show: false,
                },
                emphasis: {
                    itemStyle: {
                        shadowBlur: 20,
                        shadowColor: 'rgba(0, 0, 0, 0.5)'
                    }
                }
            }]
        };

        // 设置图表配置
        myChart.setOption(option);

        // 在上一级容器中垂直居中
        const parentContainer = myChart.getDom().parentNode;
        parentContainer.style.display = 'flex';
        parentContainer.style.justifyContent = 'center';
        parentContainer.style.alignItems = 'center';
    })
    .catch(error => {
        console.error('Error fetching data:', error);
    });
}

export function Ranking(myChart) {
    fetch('http://localhost:5000/AU')
        .then(response => response.json())
        .then(data => {
            const titlename = Object.keys(data);
            const dataValues = Object.values(data);

            const option = {
                grid: {
                    left: '0',
                    top: '0',
                    right: '0',
                    bottom: '0%',
                    containLabel: true
                },
                xAxis: {
                    show: false
                },
                yAxis: [{
                    show: true,
                    data: titlename,
                    inverse: true,
                    axisLine: { show: false },
                    splitLine: { show: false },
                    axisTick: { show: false },
                    axisLabel: {
                        textStyle: { color: '#fff' },
                    },
                }, {
                    show: false,
                    inverse: true,
                    data: dataValues,
                    axisLabel: { textStyle: { color: '#fff' } },
                    axisLine: { show: false },
                    splitLine: { show: false },
                    axisTick: { show: false },
                }],
                series: [{
                    name: '条',
                    type: 'bar',
                    yAxisIndex: 0,
                    data: dataValues,
                    barWidth: 15,
                    itemStyle: {
                        normal: {
                            barBorderRadius: 50,
                            color: '#1089E7',
                        }
                    },
                    label: {
                        normal: {
                            show: true,
                            position: 'right',
                            formatter: '{c}',
                            textStyle: { color: 'rgba(255,255,255,.5)' }
                        }
                    },
                }]
            };

            myChart.setOption(option);
            window.addEventListener("resize", function () {
                myChart.resize();
            });
        })
        .catch(error => {
            console.error('获取数据时出错：', error);
        });
}
