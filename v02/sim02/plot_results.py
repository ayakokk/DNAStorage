"""
論文Figure 6, 7, 8のプロット生成スクリプト (Python版)
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import glob

def load_data(filename):
    """データファイルを読み込む"""
    try:
        data = np.loadtxt(filename, skiprows=2)  # ヘッダー行をスキップ
        return data
    except Exception as e:
        print(f"Warning: Could not load {filename}: {e}")
        return None

def plot_channel_results(channel_name, figure_num, ylabel="Mutual Information (bits/channel use)", 
                        column=1, ylim=None):
    """指定されたチャネルの結果をプロット"""
    
    # データファイルを探す
    pattern = f"results/fig{figure_num}_*.dat"
    files = glob.glob(pattern)
    
    if not files:
        print(f"No data files found for {pattern}")
        return
    
    plt.figure(figsize=(12, 8))
    
    # 重み設定と色の対応（修正版）
    weight_configs = {
        'uniform': {'color': '#1f77b4', 'marker': 'o', 'label': 'Uniform [1,1,1,1]'},
        'light': {'color': '#ff7f0e', 'marker': 's', 'label': 'Light [6,5,4,3]'},
        'medium': {'color': '#2ca02c', 'marker': '^', 'label': 'Medium [5,4,3,2]'},
        'strong': {'color': '#d62728', 'marker': 'D', 'label': 'Strong [4,3,2,1]'},
        'extreme': {'color': '#9467bd', 'marker': 'v', 'label': 'Extreme [10,6,3,1]'},
        'light_reverce': {'color': '#8c564b', 'marker': 'P', 'label': 'Light_reverse [3,4,5,6]'},
        'medium_reverce': {'color': '#e377c2', 'marker': 'X', 'label': 'Medium_reverse [2,3,4,5]'},
        'strong_reverce': {'color': '#7f7f7f', 'marker': 'h', 'label': 'Strong_reverse [1,2,3,4]'},
        'extreme_reverce': {'color': '#bcbd22', 'marker': 'H', 'label': 'Extreme_reverse [1,3,6,10]'}
    }
    
    plotted_count = 0
    
    for file in sorted(files):
        # ファイル名から重み設定を抽出
        basename = os.path.basename(file)
        print(f"Processing file: {basename}")
        
        # ファイル名解析の改善
        parts = basename.split('_')
        if len(parts) >= 2:
            weight_type = '_'.join(parts[1:]).split('.')[0]  # 複数のアンダースコアに対応
        else:
            print(f"  Skipping {basename}: unexpected filename format")
            continue
        
        print(f"  Weight type extracted: '{weight_type}'")
        
        if weight_type not in weight_configs:
            print(f"  Skipping {weight_type}: not in weight_configs")
            print(f"  Available configs: {list(weight_configs.keys())}")
            continue
            
        data = load_data(file)
        if data is None:
            print(f"  Failed to load data from {file}")
            continue
        
        if data.shape[0] == 0:
            print(f"  No data points in {file}")
            continue
            
        config = weight_configs[weight_type]
        
        print(f"  Plotting {len(data)} data points for {weight_type}")
        
        # プロット
        plt.plot(data[:, 0], data[:, column], 
                marker=config['marker'], 
                color=config['color'],
                linewidth=2, 
                markersize=8,
                label=config['label'])
        plotted_count += 1
    
    print(f"Total plotted series: {plotted_count}")
    
    plt.xlabel('Error Rate', fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(f'Figure {figure_num}: {channel_name}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10, loc='best')
    
    if ylim:
        plt.ylim(ylim)
    
    # 保存
    output_file = f'results/figure{figure_num}_{channel_name.lower().replace(" ", "_")}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated: {output_file}")


def main():
    """メイン関数"""
    print("Generating PNG plots using Python/matplotlib...")
    
    # results ディレクトリの確認
    if not os.path.exists('results'):
        print("Error: results directory not found")
        return
    
    print("\n=== Generating Plots ===")
    
    # 各図をプロット
    
    # Figure 6: IDS Channel - Mutual Information
    plot_channel_results("IDS Channel", 6, 
                        ylabel="Mutual Information (bits/channel use)",
                        column=1)
    
    # Figure 7: ID Channel - Mutual Information  
    plot_channel_results("ID Channel", 7,
                        ylabel="Mutual Information (bits/channel use)",
                        column=1)
    
    # Figure 8: Deletion Channel - Mutual Information
    plot_channel_results("Deletion Channel", 8,
                        ylabel="Mutual Information (bits/channel use)", 
                        column=1)
    
    print("\nAll Python plots generated successfully!")
    print("Check results/ directory for PNG files.")

if __name__ == "__main__":
    main()