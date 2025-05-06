#!/bin/bash

echo "==================================="
echo "GCH1 多腫瘤基因網路視覺化生成器"
echo "==================================="
echo ""

echo "正在生成網路視覺化..."
python multi_tumor_network.py

echo ""
echo "網路視覺化已成功生成！"
echo ""
echo "是否在瀏覽器中打開索引頁面？(Y/N)"
read -p "> " choice

if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        open "index.html"
    else
        # Linux
        xdg-open "index.html" &>/dev/null || firefox "index.html" &>/dev/null || google-chrome "index.html" &>/dev/null
    fi
    echo "已打開索引頁面。"
else
    echo "您可以稍後手動打開 index.html 查看結果。"
fi

echo ""
echo "按Enter鍵退出..."
read 