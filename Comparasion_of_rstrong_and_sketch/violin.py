import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 定义参数范围
k_values = [20, 30, 40, 50, 60]
r_values = ["0.001", "0.01", "0.05", "0.2"]

# 读取数据函数
def read_data(file_name):
    try:
        # 每行两个值，分别为 r_hat 和 r_prime
        data = pd.read_csv(file_name, header=None, names=["r_hat", "r_prime"])
        return data
    except Exception as e:
        print(f"Error reading file {file_name}: {e}")
        return None

# 合并所有数据到一个 DataFrame
data_list = []
for r in r_values:
    for k in k_values:
        file_name = f"results/r{r}_k{k}.output"  # 文件名格式
        if os.path.exists(file_name):
            df = read_data(file_name)
            if df is not None:
                df["r"] = float(r)  # 添加 r 值（转换为浮点数）
                df["k"] = k  # 添加 k 值
                data_list.append(df)
        else:
            print(f"Warning: File {file_name} not found.")

# 合并数据
if len(data_list) > 0:
    all_data = pd.concat(data_list, ignore_index=True)
else:
    raise ValueError("No data available. Check your file paths or inputs.")

# 创建画布
fig, axes = plt.subplots(
    nrows=len(r_values),
    ncols=len(k_values) * 2,  # 每个 r 和 k 有两个子图，分别绘制 r_hat 和 r_prime
    figsize=(20, 15),
    sharey="row"
)

# 遍历 r 和 k 的组合
for i, r in enumerate(r_values):
    for j, k in enumerate(k_values):
        # 筛选当前 r 和 k 的数据
        subset = all_data[(all_data["r"] == float(r)) & (all_data["k"] == k)]

        # 选择两个子图：一个绘制 r_hat，一个绘制 r_prime
        ax_r_hat = axes[i, j * 2]  # 左边的子图绘制 r_hat
        ax_r_prime = axes[i, j * 2 + 1]  # 右边的子图绘制 r_prime

        # 绘制 r_hat 的小提琴图
        sns.violinplot(
            data=subset,
            x=["r_hat"] * len(subset),  # 固定 x 为 r_hat
            y="r_hat",
            ax=ax_r_hat,
            inner="box",
            density_norm="width",
            color="blue"
        )
        ax_r_hat.axhline(float(r), color="red", linestyle="--", linewidth=1)  # 水平线
        ax_r_hat.set_title(f"r_hat (r={r}, k={k})", fontsize=10)

        # 绘制 r_prime 的小提琴图
        sns.violinplot(
            data=subset,
            x=["r_prime"] * len(subset),  # 固定 x 为 r_prime
            y="r_prime",
            ax=ax_r_prime,
            inner="box",
            density_norm="width",
            color="orange"
        )
        ax_r_prime.axhline(float(r), color="red", linestyle="--", linewidth=1)  # 水平线
        ax_r_prime.set_title(f"r_prime (r={r}, k={k})", fontsize=10)

# 调整整体布局
fig.tight_layout(pad=3.0)

plt.savefig("violin_r_hat_vs_r_prime.png")

# import os
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt

# # 定义参数范围
# k_values = [20, 30, 40, 50, 60]
# r_values = ["0.001", "0.01", "0.05", "0.2"]

# # 读取数据函数
# def read_data(file_name):
#     try:
#         # 每行两个值，分别为 r_hat 和 r_prime
#         data = pd.read_csv(file_name, header=None, names=["r_hat", "r_prime"])
#         return data
#     except Exception as e:
#         print(f"Error reading file {file_name}: {e}")
#         return None

# # 合并所有数据到一个 DataFrame
# data_list = []
# for r in r_values:
#     for k in k_values:
#         file_name = f"results/r{r}_k{k}.output"  # 文件名格式
#         if os.path.exists(file_name):
#             df = read_data(file_name)
#             if df is not None:
#                 df["r"] = float(r)  # 添加 r 值（转换为浮点数）
#                 df["k"] = k  # 添加 k 值
#                 data_list.append(df)
#         else:
#             print(f"Warning: File {file_name} not found.")

# # 合并数据
# if len(data_list) > 0:
#     all_data = pd.concat(data_list, ignore_index=True)
# else:
#     raise ValueError("No data available. Check your file paths or inputs.")

# # 创建画布
# fig, axes = plt.subplots(
#     nrows=len(r_values), ncols=len(k_values),
#     figsize=(20, 10),
#     sharex=False, sharey=False  # 不共享 y 轴，允许自动调整
# )

# # 遍历 r 和 k 的组合
# for i, r in enumerate(r_values):
#     for j, k in enumerate(k_values):
#         # 筛选当前 r 和 k 的数据
#         subset = all_data[(all_data["r"] == float(r)) & (all_data["k"] == k)]

#         # 获取当前子图
#         ax = axes[i, j]

#         # 绘制小提琴图
#         if not subset.empty:
#             sns.violinplot(
#                 data=subset,
#                 y="r_hat",  # 只绘制 r_hat
#                 ax=ax,
#                 inner="box",
#                 density_norm="width",  # 替代 scale 参数
#                 color="blue"
#             )
#             # 自动调整 y 轴范围
#             y_min, y_max = subset["r_hat"].min(), subset["r_hat"].max()
#             padding = (y_max - y_min) * 0.3  # 留 10% 空白
#             ax.set_ylim(y_min - padding, y_max + padding)

#             # 在子图中添加水平线表示 r 的值
#             ax.axhline(float(r), color="red", linestyle="--", linewidth=1, label=f"r = {r}")
#             ax.legend(loc="upper right", fontsize=8)

#         # 设置标题
#         ax.set_title(f"r={r}, k={k}", fontsize=10)

# # 调整整体布局
# fig.tight_layout(pad=3.0)
# plt.savefig("violin_r_hat.png")