#!/usr/bin/env python3
"""
CSVファイルから正規化されたIxy値をプロットしてPDFに保存するスクリプト
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 日本語フォントの設定（必要に応じて）
# plt.rcParams['font.family'] = 'DejaVu Sans'

def plot_normalized_ixy():
    # CSVファイルを読み込む
    csv_file = 'result/ex01b_normalized_ixy_table.csv'
    df = pd.read_csv(csv_file)
    
    # データの準備
    pi_values = df['Pi'].values
    
    # プロットの設定
    plt.figure(figsize=(10, 6))
    
    # マーカーとラインスタイルの定義
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', 'h']
    linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-']
    colors = plt.cm.tab10(np.linspace(0, 1, 9))
    
    # 各条件をプロット
    conditions = ['Normal', 'Extreme', 'Extreme_reverse', 'Light', 'Light_reverse', 
                  'Medium', 'Medium_reverse', 'Strong', 'Strong_reverse']
    
    for i, condition in enumerate(conditions):
        if condition in df.columns:
            # N/Aを除外してプロット
            values = df[condition].replace('N/A', np.nan)
            
            # 数値に変換（科学的記数法に対応）
            numeric_values = pd.to_numeric(values, errors='coerce')
            
            # プロット
            plt.plot(pi_values, numeric_values, 
                    marker=markers[i], 
                    linestyle=linestyles[i],
                    color=colors[i],
                    label=condition,
                    markersize=8,
                    linewidth=2)
    
    # グラフの装飾
    plt.xlabel('Pi=Pd', fontsize=14)
    plt.ylabel('Mutual information (per channel bit)', fontsize=14)
    plt.title('Normalized Ixy values for different conditions', fontsize=16)
    
    # グリッドを追加
    plt.grid(True, alpha=0.3, linestyle='--')
    
    # 凡例を追加（グラフの外側に配置）
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # x軸の範囲を設定
    plt.xlim(0.001, 0.010)
    
    # y軸の範囲を自動調整（データに基づいて）
    all_values = []
    for condition in conditions:
        if condition in df.columns:
            values = pd.to_numeric(df[condition].replace('N/A', np.nan), errors='coerce')
            all_values.extend(values.dropna().tolist())
    
    if all_values:
        y_min = min(all_values) * 0.95
        y_max = max(all_values) * 1.05
        plt.ylim(y_min, y_max)
    
    # レイアウトを調整
    plt.tight_layout()
    
    # PDFとして保存
    output_pdf = 'result/ex01b_normalized_ixy_plot.pdf'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_pdf}")
    
    # PNGとしても保存（プレビュー用）
    output_png = 'result/ex01b_normalized_ixy_plot.png'
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    print(f"Plot also saved as PNG: {output_png}")
    
    # グラフを表示（オプション）
    # plt.show()

def plot_comparison_subplots():
    """
    正逆ペアごとに比較するサブプロット版
    """
    csv_file = 'result/ex01b_normalized_ixy_table.csv'
    df = pd.read_csv(csv_file)
    
    pi_values = df['Pi'].values
    
    # 2x2のサブプロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    pairs = [
        ('Extreme', 'Extreme_reverse'),
        ('Light', 'Light_reverse'),
        ('Medium', 'Medium_reverse'),
        ('Strong', 'Strong_reverse')
    ]
    
    for idx, (forward, reverse) in enumerate(pairs):
        ax = axes[idx]
        
        # Forward条件
        if forward in df.columns:
            values = pd.to_numeric(df[forward].replace('N/A', np.nan), errors='coerce')
            ax.plot(pi_values, values, 'o-', label=forward, markersize=8, linewidth=2)
        
        # Reverse条件
        if reverse in df.columns:
            values = pd.to_numeric(df[reverse].replace('N/A', np.nan), errors='coerce')
            ax.plot(pi_values, values, 's--', label=reverse, markersize=8, linewidth=2)
        
        ax.set_xlabel('Pi=Pd=Ps')
        ax.set_ylabel('Normalized Ixy')
        ax.set_title(f'{forward.split("_")[0]} comparison')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlim(0.003, 0.010)
    
    plt.suptitle('Normalized Ixy: Forward vs Reverse conditions', fontsize=16)
    plt.tight_layout()
    
    # PDFとして保存
    output_pdf = 'result/ex01b_comparison_subplots.pdf'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved to: {output_pdf}")

if __name__ == "__main__":
    # メインプロット
    plot_normalized_ixy()
    
    # 比較用サブプロット
    plot_comparison_subplots()