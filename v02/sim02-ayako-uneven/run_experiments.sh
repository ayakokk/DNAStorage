#!/bin/bash

# 論文 Figure 6, 7, 8 の実験を実行するスクリプト

echo "=== Running experiments for Figure 6, 7, 8 reproduction ==="

# パラメータ設定
ICB_DIR="./ICB/ex01"
N=100
SEED=1

# 結果ディレクトリ作成
mkdir -p results

echo "Running comprehensive experiment..."
echo "  ICB Directory: ${ICB_DIR}"
echo "  Block length N: ${N}"
echo "  Random seed: ${SEED}"
echo ""

# 実行
./sim_all ${ICB_DIR} ${N} ${SEED} > results/comprehensive_results.txt 2> results/experiment_log.txt

if [ $? -eq 0 ]; then
    echo "✓ Comprehensive experiment completed successfully"
else
    echo "✗ Experiment failed"
    exit 1
fi

echo ""
echo "=== Extracting results for plotting ==="

# 結果ファイルから各チャネル・重み設定の組み合わせを抽出
python3 << 'EOF'
import re
import os

# 結果ファイルを読み込み
with open('results/comprehensive_results.txt', 'r') as f:
    content = f.read()

# チャネルタイプと重み設定を抽出
channels = ['IDS', 'ID', 'DEL']
weight_configs = ['uniform', 'light', 'medium', 'strong', 'extreme', 'light_reverce', 'medium_reverce', 'strong_reverce', 'extreme_reverce']

for channel in channels:
    for weight in weight_configs:
        # 該当する章を探す
        pattern = f"# ========== {channel} Channel Experiment ==========.*?## Weight Config: {weight}"
        match = re.search(pattern, content, re.DOTALL)
        
        if match:
            # データ行を抽出（数字で始まる行）
            start_pos = match.end()
            
            # 次の## Weight Config:または次の# ==========まで探す
            next_pattern = r"(## Weight Config:|# ==========)"
            next_match = re.search(next_pattern, content[start_pos:])
            
            if next_match:
                data_section = content[start_pos:start_pos + next_match.start()]
            else:
                data_section = content[start_pos:]
            
            # 数値行を抽出
            data_lines = []
            for line in data_section.split('\n'):
                if re.match(r'^\d+\.\d+', line.strip()):
                    data_lines.append(line.strip())
            
            # ファイルに保存
            filename = f'results/fig{channels.index(channel)+6}_{weight}.dat'
            with open(filename, 'w') as f:
                f.write(f"# {channel} Channel - {weight} weighting\n")
                f.write("# Error_Rate\tMutual_Info\tEntropy_X\tEntropy_XY\tSER\n")
                for line in data_lines:
                    f.write(line + '\n')
            
            print(f"Created: {filename} ({len(data_lines)} data points)")

print("Data extraction completed!")
EOF

echo ""
echo "=== Results Summary ==="
echo "Files created in results/ directory:"
echo "  - comprehensive_results.txt: Full experiment output"
echo "  - experiment_log.txt: Progress log and error messages"
echo "  - fig6_*.dat: IDS channel data (Figure 6) for each weight config"
echo "  - fig7_*.dat: ID channel data (Figure 7) for each weight config"
echo "  - fig8_*.dat: Deletion channel data (Figure 8) for each weight config"

echo ""
echo "=== Available data files ==="
ls -la results/*.dat 2>/dev/null || echo "No .dat files found"

echo ""
echo "=== Sample plotting commands (gnuplot) ==="
cat << 'EOF'
# Figure 6 (IDS Channel)
set xlabel "Error Rate"
set ylabel "Mutual Information (bits)"
set title "Figure 6: IDS Channel"
plot 'results/fig6_uniform.dat' u 1:2 w lp title 'Uniform', \
     'results/fig6_light.dat' u 1:2 w lp title 'Light', \
     'results/fig6_medium.dat' u 1:2 w lp title 'Medium', \
     'results/fig6_strong.dat' u 1:2 w lp title 'Strong', \
     'results/fig6_extreme.dat' u 1:2 w lp title 'Extreme' \
     'results/fig6_light_reverce.dat' u 1:2 w lp title 'Light_reverce', \
     'results/fig6_medium_reverce.dat' u 1:2 w lp title 'Medium_reverce', \
     'results/fig6_strong_reverce.dat' u 1:2 w lp title 'Strong_reverce', \
     'results/fig6_extreme_reverce.dat' u 1:2 w lp title 'Extreme_reverce'

# Figure 7 (ID Channel)
set title "Figure 7: ID Channel"
plot 'results/fig7_uniform.dat' u 1:2 w lp title 'Uniform', \
     'results/fig7_light.dat' u 1:2 w lp title 'Light', \
     'results/fig7_medium.dat' u 1:2 w lp title 'Medium', \
     'results/fig7_strong.dat' u 1:2 w lp title 'Strong', \
     'results/fig7_extreme.dat' u 1:2 w lp title 'Extreme' \
     'results/fig7_light_reverce.dat' u 1:2 w lp title 'Light_reverce', \
     'results/fig7_medium_reverce.dat' u 1:2 w lp title 'Medium_reverce', \
     'results/fig7_strong_reverce.dat' u 1:2 w lp title 'Strong_reverce', \
     'results/fig7_extreme_reverce.dat' u 1:2 w lp title 'Extreme_reverce' \

# Figure 8 (Deletion Channel)
set title "Figure 8: Deletion Channel"
plot 'results/fig8_uniform.dat' u 1:2 w lp title 'Uniform', \
     'results/fig8_light.dat' u 1:2 w lp title 'Light', \
     'results/fig8_medium.dat' u 1:2 w lp title 'Medium', \
     'results/fig8_strong.dat' u 1:2 w lp title 'Strong', \
     'results/fig8_extreme.dat' u 1:2 w lp title 'Extreme' \
     'results/fig8_light_reverce.dat' u 1:2 w lp title 'Light_reverce', \
     'results/fig8_medium_reverce.dat' u 1:2 w lp title 'Medium_reverce', \
     'results/fig8_strong_reverce.dat' u 1:2 w lp title 'Strong_reverce', \
     'results/fig8_extreme_reverce.dat' u 1:2 w lp title 'Extreme_reverce'
EOF

echo ""
echo "Experiments completed successfully!"