# HCA_script
### How to use the git command
https://blog.csdn.net/qq_37623429/article/details/80649266

1.通过下面2条命令，设置用户名和邮箱地址（相当于注册一个账号来代表这台电脑）
```
$ git config --global user.name "haojiang@HCA"
$ git config --global user.email "haojiang9999@outlook.com"
```
2.创建本地仓库（repository）
```
$ pwd       //显示当前目录
$ cd e: /gitTest    //我打算将git的仓库建立在e盘的gitTest文件夹下。
$ git init        //将当前路径初始化为git仓库
```
3、将文件添加到仓库
```
$ vim hello.txt ## create a test file
$ git add hello.txt //将hello.txt文件添加到仓库
$ gti commit -m "add hello.txt"  //把文件提交到仓库，-m 之后的"add hello.txt"是对这次操作自己添加的描述，可以随便写.
$ git status   //查看仓库状态
```
4、设置ssh（仅首次使用需要配置）
```
## （注意其中C是大写的）然后一路回车到底
$ ssh-keygen -t rsa -C "youremail@example.com"
```
在提示的路径中的.ssh文件夹中有2个文件，分别为id_rsa(私钥)和id_rsa.pub（公钥）

用记事本打开id_rsa.pub,复制其中的内容，打开并登陆github。

在个人设置（Settings）中找到SSH设置，点击添加sshkey（new ssh key）

5、在github上新建一个仓库（repository）

6、将本地仓库同步到空远程仓库（repository）

https://stackoverflow.com/questions/10904339/github-fatal-remote-origin-already-exists
```
## 1.If the repository is empty
$ git remote add origin git@github.com:mozhilei/test1.git
$ git push -u origin master
## 2.Error: “fatal: remote origin already exists”
# you should just update the existing remote:
$ git remote set-url origin git@github.com:ppreyer/first_app.git
## 3.获取远程库与本地同步合并（如果远程库不为空必须做这一步，否则后面的提交会失败）
$ git pull --rebase origin master
$ git push -u origin master
```

7、状态查询命令
```
git status
```
