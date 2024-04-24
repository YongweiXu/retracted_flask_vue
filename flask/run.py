from flask import Flask, render_template  #引入包中要使用的类

app = Flask(__name__, static_url_path='/', static_folder='./../flask-dist', template_folder='./../flask-dist')


#生成app对象  static_folder 设置资源位置  就是 js文件夹 css文件夹的目录   template_folder html文件位置


@app.route('/')
def index():
    return render_template('index.html')


# template_folder 设置的路径下的 index.html

if __name__ == '__main__':
    app.run(debug=True)
