import akshare as ak
import pandas as pd
# 获取 A 股股票列表
print("Start fetching stock list...")
stocks = ak.stock_info_a_code_name()
# 剔除 ST 和 PT 股票
print("Filtering ST and PT stocks...")
filtered_stocks = stocks[~stocks['name'].str.contains('ST|PT')]
# 日期范围
start_date = "20190531"
end_date = "20200104"
# 存储数据
all_data = []

print("Start fetching historical data...")
print("Total stocks:", len(filtered_stocks))
for code in filtered_stocks['code']:
    I = 1
    try:
        # 获取单个股票的行情数据
        df = ak.stock_zh_a_hist(symbol=code, period="daily", start_date=start_date, end_date=end_date, adjust="qfq")
        df['code'] = code  # 添加股票代码
        all_data.append(df)
        print(I,"Fetched data for", code)
        I += 1
        if I > 500:
            break
    except Exception as e:
        print(I,f"Error fetching data for {code}: {e}")
    
    

# 合并所有股票的数据
print("Merging data...")
result = pd.concat(all_data, ignore_index=True)

# 保存数据
print("Saving data...")
result.to_csv("a_shares_historical_data.csv", index=False)