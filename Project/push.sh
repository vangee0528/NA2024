#!/bin/bash

# 将所有新增的内容加入暂存区
git add .

# 提交更改，备注默认为 "update"
git commit -m "update"

# 推送到远程仓库
git push