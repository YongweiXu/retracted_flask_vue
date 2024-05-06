<template>
  <div class="head clearfix">
      <h1 class="pulll_left">生物医学期刊撤稿大数据可视化看板</h1>
      <div class="menu menu2 pulll_left">
        <ul>
          <li><router-link to="/">首页</router-link></li>
          <li><router-link to="/homepage">数据看板</router-link></li>
          <li><router-link to="/homepage2">影响力看板</router-link></li>
          <li><router-link to="/predictionpage">摘要预测</router-link></li>
          <li><router-link to="/query">数据查询</router-link></li>
        </ul>
      </div>
      <div class="time">{{ currentTime }}</div>
  </div>
  <div class="container">
    <form @submit.prevent="handleSubmit" class="search-form">
      <label for="searchTerm">搜索词:</label>
      <input type="text" id="searchTerm" v-model="searchTerm" placeholder="请输入搜索词">
      <button type="submit">搜索</button>
    </form>

    <div class="table-container">
      <table v-if="pagedRecords.length > 0 || loading" class="data-table">
        <thead>
          <tr>
            <th>PMID</th>
            <th>ISSN</th>
            <th>收录状态</th>
            <th>时间</th>
            <th>标题</th>
            <th>摘要</th>
            <th>作者</th>
          </tr>
        </thead>
        <tbody>
          <tr v-for="(record, index) in pagedRecords" :key="index" class="table-row">
            <td>{{ record.PMID }}</td>
            <td>{{ record.ISSN }}</td>
            <td>{{ record.STAT }}</td>
            <td>{{ record.DOCM }}</td>
            <td>{{ record.TI }}</td>
            <td>{{ record.AB }}</td>
            <td>
              <span v-for="(author, authorIndex) in record.AU" :key="authorIndex">
                {{ author }}
                <template v-if="authorIndex !== record.AU.length - 1">, </template>
              </span>
            </td>
          </tr>
          <tr v-if="loading && pagedRecords.length === 0">
            <td colspan="7">加载中...</td>
          </tr>
          <tr v-if="!loading && pagedRecords.length === 0">
            <td colspan="7">暂无数据</td>
          </tr>
        </tbody>
      </table>
    </div>

    <div class="pagination">
      <button :disabled="currentPage === 1" @click="previousPage">上一页</button>
      <span>{{ currentPage }}</span>
      <button :disabled="currentPage === totalPages" @click="nextPage">下一页</button>
    </div>
  </div>
</template>

<script>
import axios from 'axios';

export default {
  data() {
    return {
      searchTerm: '',
      records: [],
      loading: false,
      currentPage: 1,
      pageSize: 10 // 每页显示的记录数
    };
  },
  computed: {
    totalPages() {
      return Math.ceil(this.records.length / this.pageSize);
    },
    pagedRecords() {
      const startIndex = (this.currentPage - 1) * this.pageSize;
      const endIndex = startIndex + this.pageSize;
      return this.records.slice(startIndex, endIndex);
    }
  },
  methods: {
    handleSubmit() {
      this.loading = true;
      axios.post('http://localhost:5000/get_pubmed_records', {term: this.searchTerm})
          .then(response => {
            this.loading = false;
            this.records = response.data.records;
            this.currentPage = 1; // 重置到第一页
          })
          .catch(error => {
            this.loading = false;
            console.error('Error fetching PubMed records:', error);
          });
    },
    previousPage() {
      if (this.currentPage > 1) {
        this.currentPage--;
      }
    },
    nextPage() {
      if (this.currentPage < this.totalPages) {
        this.currentPage++;
      }
    }
  }
};
</script>

<style>
.container {
  text-align: center;
}

.search-form {
  margin-bottom: 20px;
}

.table-container {
  max-width: 100%;
  overflow-x: auto;
}

.data-table {
  width: 80%;
  border-collapse: collapse;
  margin: auto;
}

.data-table th,
.data-table td {
  border: 1px solid #ddd;
  padding: 8px;
}

.data-table th {
  background-color: #5da6eb;
  color: #333; /* 避免与背景色冲突 */
}

.table-row:hover {
  background-color: #63cce2;
}

.pagination {
  margin-top: 20px;
}
</style>
