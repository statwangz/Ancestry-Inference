# R语言高性能数据处理包：data.table

（更新中）

王志玮

## 介绍

R中的`data.table`包是一个`data.frame`的高级版本，它是一个高性能数据处理包，快速高效，适合处理大型数据集，例如基因组等。但是两者的用法又有一些区别，如果直接套用可能会出错，所以这里简单介绍一下`data.table`的基本语法。

## 下载安装

```r
install.packages("data.table")
library(data.table)
```

## 创建

```r
name <- c("wang", "zhao", "peng", "zhou")
chinese <- c(85, 79, 92, 90)
math <- c(96, 98, 81, 84)
english <- c(78, 87, 89, 99)
dt <- data.table(name, chinese, math, english)
```

## 读写文件

使用`fread()`和`fwrite()`函数可以处理`.txt` `.csv` `.dat` `.raw`等格式的文件，下面是基本用法：
```r
fread("filename")    # 读
fwrite(dt, "dt.csv") # 写
```

一些参数：
```r
fread("filename", stringsAsFactors = T)) # 将字符型转换为因子型，默认为F
fread("filename", data.table = F)        # 生成data.frame，默认为T
```

## 转换

使用`setDF()`函数将数据框从`data.table`转换为`data.frame`：
```r
df <- setDF(dt)
```

使用`setDT()`函数将数据框从`data.frame`转换为`data.table`：
```r
dt <- setDT(df)
```

## 选取

选取数据时注意与`data.frame`的区别。

### 选取点

```r
dt[1, 2]        # 返回"data.table"，注意与"data.frame"的区别
dt[[1, 2]]      # 返回一个值
dt[1:3, 3]      # 返回"data.table"
dt[1:3, "math"] # 返回"data.table"
dt[1:3, math]   # 返回一个向量
```

### 选取行

```r
dt[1]     # []里可以只写一个数，不需要加逗号
dt[1, ]   # 也可以加逗号
dt[1:3]   # 多行
dt[1:3, ] # 多行
#以上都是返回"data.table"
```

### 选取列

根据坐标选取：
```r
dt[ , 2]   # 返回"data.table"，注意与"data.frame"的区别
dt[[2]]    # [[]]里只写一个数选取列，返回一个向量
dt[ , 1:3] # 多列
```

根据列名选取：
```r
dt[["math"]]             # 返回一个向量
dt$math                  # 返回一个向量

dt[ , "math"]            # 列名用""括起来，返回"data.table"
dt[ , c("name", "math")] # 多列，返回"data.table"
dt[ , name:math]         # 多列，返回"data.table"

dt[ , math]              # 直接写列名表达式（expression，没有""），返回一个向量
dt[ , c(name, math)]     # 将多列合并成一个向量输出，数据降维

# 如果想使用表达式且数据不降维，即返回"data.table"，可使用list()
dt[ , list(math)]
dt[ , .(math)]           # .()是list()的简写
dt[ , .(name, math)]     # 多列
```
如果将列名存储在一个变量中，并且想要通过这个变量选取指定列，因为用`data.table[]`选取时可以接受一个不加`""`的对象，所以会将这个变量名当作列名，产生混淆，这时就需要使用`with`参数：
```r
subject <- c("chinese", "math")
dt[ , subject]           # 报错
dt[ , subject, with = F] # 选取出chinese和math列，with默认为T
```
这里`with`参数的作用就类似于`with()`函数。

### 筛选

```r
dt[chinese>80 & math<85]
dt[name=="wang" | name=="zhao"]

# 使用on参数选取某一列是某一个值的行
dt[c("wang", "zhao"), on = "name"]

# 按照多列条件筛选
dt[.("wang", 96), on = .(name, math)]

# 找不到这样的行则生成满足列条件的新行，其他列的值为NA，但不改变原数据集（dt）
dt[.("wang", c(89, 96)), on = .(name, math)]

# 找不到则生成新行并填充
dt[.("wang", c(89, 96)), on = .(name, math), roll = -Inf]

# 找不到也不生成新行
dt[.("wang", c(89, 96)), on = .(name, math), nomatch = 0]
```

## 删除

删除行:
```r
# 删除第2、3行
dt[-c(2, 3)]
dt[!c(2, 3)]
```

删除列:
```r
# 删除第2、4列
dt[ , -c(2, 4)]
dt[ , !c(2, 4)]
dt[ , -c("chinese", "english")]
dt[ , !c("chinese", "english")]
```

## 排序

```r
# 按照"math"排序，默认升序
dt[order(math)]
dt[order(math), ]
dt[order(math, decreasing = T)] # 降序
```

## []内操作

`data.table`可以通过直接在`[]`中输入指令来进行计算，`[]`里可以有三个参数，第一个表示对哪些行进行操作，第二个表示对哪些列进行什么操作，第三个表示筛选条件。下面是一些简单例子：
```r
# 对所有行的"math"值求平均值
dt[ , mean(math)]
# 对前三行的"math"值求平均值
dt[1:3, mean(math)]

# 同时进行多种操作，结果长度不同时自动循环补齐，返回data.table
dt[ , .(sort(chinese), sum(math))]
# 返回向量
dt[ , c(sort(chinese), sum(math))]

# 自定义结果列名
dt[ , .(ch=mean(chinese), ma=sum(math))]
dt[ , c(ch=mean(chinese), ma=sum(math))] # 注意结果列名与上面的区别

# 筛选
dt[c("wang", "zhao"), math+2, on = "name"]
```

分组操作：
```r
sex <- c("M", "F", "F", "M")
class <- c(1, 2, 2, 3)
dt <- data.table(name, class, sex, chinese, math, english)

dt[ , mean(math), by=class]             # 按照"class"值进行分组计算
dt[ , mean(math), keyby=sex]            # 分组计算并按照"class"值排序
dt[ , mean(math), by=chinese<90]        # 按照"chinese"值条件分组
dt[ , mean(math), by=.(class, sex)]     # 按照多个条件分组
dt[ , mean(math), by=c("class", "sex")]
```

## key

`data.table`的行名默认为`"1" "2" "3"`，且无法像`data.frame`一样通过`rownames()`函数有效地修改:
```r
rownames(dt) # "1" "2" "3" "4"
rownames(dt) <- letters[1:4]
rownames(dt) # "a" "b" "c" "d"
```
从输出结果来看修改成功了，但实际上：
```r
dt      # 输出数据框时行名仍是"1" "2" "3" "4"
dt["a"] # 报错
dt[1]   # 成功
```
即无法通过“修改后”的行名选取数据。

但是`data.table`有更强大的工具`key`，我们可以通过设置`key`来指定任意一列作为“行名”:
```r
haskey(dt)            # 返回TF，检查是否有key
setkey(dt, name)      # 指定name作为key
haskey(dt)
key(dt)               # 查看key

dt["zhao"]            # 选取"zhao"这一行
dt[c("wang", "zhao")] # 选取多行
dt["zhao", math]
dt["zhao", "math"]

# 当指定的行名为数字格式时
setkey(dt, class)     # 指定class作为key
dt[3]                 # 选取第3行
dt[.(3)]              # 选取"class"为3的行
dt[J(3)]              # 选取"class"为3的行
dt[.(c(1, 3))]        # 选取多行
dt[.(J(1, 3))]        # 选取多行

# 当指定的key有重复时，根据key选取会将同名的行都选取出来
dt[.(2)]
```

当用`as.data.table()`函数将`data.frame`转换成`data.table`时，默认抛弃行名，但是也可以通过设置参数`keep.rownames`保留行名成为新的一列：
```r
df <- data.frame(class, sex, chinese, math, english, row.names = name)
dt1 <- as.data.table(df)
dt2 <- as.data.table(df, keep.rownames = T)         # 保留原数据框行名作为一列，列名为"rn"
dt3 <- as.data.table(df, keep.rownames = "name") # 自定义列名
```