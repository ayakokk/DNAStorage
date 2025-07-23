#!/bin/bash

# 全ての条件を定義
conditions=(
    "Normal"
    "Extreme"
    "Extreme_reverse"
    "Light"
    "Light_reverse"
    "Medium"
    "Medium_reverse"
    "Strong"
    "Strong_reverse"
)

# 結果を保存するディレクトリを作成
mkdir -p result

# メインの結果ファイル
OUTPUT_FILE="result/ex01b_normalized_ixy_table.txt"
CSV_FILE="result/ex01b_normalized_ixy_table.csv"

# 一時ファイルでデータを収集
temp_data=$(mktemp)

echo "=== Extracting normalized Ixy values from all conditions ==="

# 各条件からデータを抽出
for condition in "${conditions[@]}"; do
    echo "Processing $condition..."
    
    OUTPUT_DIR="log/${condition}"
    
    for i in 002 003 004 005 006 007 008 009; do
        LOGFILE="${OUTPUT_DIR}/ex01b_${i}_${condition}.log"
        
        if [ -f "$LOGFILE" ]; then
            # 最終行（100000で始まる行）から最後のフィールド（Ixy）を抽出
            FINAL_IXY=$(awk '/^100000 / {print $NF}' "$LOGFILE")
            
            if [ -n "$FINAL_IXY" ]; then
                PI_VALUE=$(echo "scale=3; $i / 1000" | bc -l)
                
                # 正規化された値を計算
                NORMALIZED_IXY=$(echo "$FINAL_IXY" | awk '{nu=6; printf "%.6E", $1/nu}')
                
                # 一時ファイルに保存
                echo "$PI_VALUE $condition $NORMALIZED_IXY" >> "$temp_data"
            else
                # データが見つからない場合はN/Aを記録
                PI_VALUE=$(echo "scale=3; $i / 1000" | bc -l)
                echo "$PI_VALUE $condition N/A" >> "$temp_data"
            fi
        else
            # ファイルが見つからない場合はN/Aを記録
            PI_VALUE=$(echo "scale=3; $i / 1000" | bc -l)
            echo "$PI_VALUE $condition N/A" >> "$temp_data"
        fi
    done
done

# テーブル形式のテキストファイルを作成
echo "# Normalized Ixy values (Ixy/nu) for all conditions" > $OUTPUT_FILE
echo "# nu=6, log(2)=0.693147" >> $OUTPUT_FILE
echo "#" >> $OUTPUT_FILE

# ヘッダー行を作成（固定幅で整形）
printf "%-8s" "Pi" >> $OUTPUT_FILE
for condition in "${conditions[@]}"; do
    printf "%-18s" "$condition" >> $OUTPUT_FILE
done
echo "" >> $OUTPUT_FILE

# 区切り線
printf "%-8s" "--------" >> $OUTPUT_FILE
for condition in "${conditions[@]}"; do
    printf "%-18s" "------------------" >> $OUTPUT_FILE
done
echo "" >> $OUTPUT_FILE

# データ行を作成
for i in 002 003 004 005 006 007 008 009; do
    PI_VALUE=$(echo "scale=3; $i / 1000" | bc -l)
    printf "%-8s" "$PI_VALUE" >> $OUTPUT_FILE
    
    for condition in "${conditions[@]}"; do
        # 該当するデータを検索
        VALUE=$(awk -v pi="$PI_VALUE" -v cond="$condition" '$1==pi && $2==cond {print $3}' "$temp_data")
        if [ -n "$VALUE" ]; then
            printf "%-18s" "$VALUE" >> $OUTPUT_FILE
        else
            printf "%-18s" "N/A" >> $OUTPUT_FILE
        fi
    done
    echo "" >> $OUTPUT_FILE
done

# CSV形式のファイルも作成
echo "Pi,Normal,Extreme,Extreme_reverse,Light,Light_reverse,Medium,Medium_reverse,Strong,Strong_reverse" > $CSV_FILE

# CSVデータ行
for i in 002 003 004 005 006 007 008 009; do
    PI_VALUE=$(echo "scale=3; $i / 1000" | bc -l)
    echo -n "$PI_VALUE" >> $CSV_FILE
    
    for condition in "${conditions[@]}"; do
        VALUE=$(awk -v pi="$PI_VALUE" -v cond="$condition" '$1==pi && $2==cond {print $3}' "$temp_data")
        if [ -n "$VALUE" ]; then
            echo -n ",$VALUE" >> $CSV_FILE
        else
            echo -n ",N/A" >> $CSV_FILE
        fi
    done
    echo "" >> $CSV_FILE
done

# LaTeX形式のテーブルも作成（オプション）
LATEX_FILE="result/ex01b_normalized_ixy_table.tex"
echo "% LaTeX table of normalized Ixy values" > $LATEX_FILE
echo "\\begin{table}[htbp]" >> $LATEX_FILE
echo "\\centering" >> $LATEX_FILE
echo "\\caption{Normalized Ixy values (Ixy/nu) for all conditions}" >> $LATEX_FILE
echo "\\begin{tabular}{|c|cccccccc|}" >> $LATEX_FILE
echo "\\hline" >> $LATEX_FILE
echo -n "Pi " >> $LATEX_FILE
for condition in "${conditions[@]}"; do
    # アンダースコアをエスケープ
    escaped_condition=$(echo "$condition" | sed 's/_/\\_/g')
    echo -n "& $escaped_condition " >> $LATEX_FILE
done
echo "\\\\" >> $LATEX_FILE
echo "\\hline" >> $LATEX_FILE

# LaTeXデータ行
for i in 002 003 004 005 006 007 008 009; do
    PI_VALUE=$(echo "scale=3; $i / 1000" | bc -l)
    echo -n "$PI_VALUE " >> $LATEX_FILE
    
    for condition in "${conditions[@]}"; do
        VALUE=$(awk -v pi="$PI_VALUE" -v cond="$condition" '$1==pi && $2==cond {print $3}' "$temp_data")
        if [ -n "$VALUE" ] && [ "$VALUE" != "N/A" ]; then
            # 科学的記数法を維持
            echo -n "& $VALUE " >> $LATEX_FILE
        else
            echo -n "& N/A " >> $LATEX_FILE
        fi
    done
    echo "\\\\" >> $LATEX_FILE
done

echo "\\hline" >> $LATEX_FILE
echo "\\end{tabular}" >> $LATEX_FILE
echo "\\end{table}" >> $LATEX_FILE

# 一時ファイルを削除
rm "$temp_data"

# 結果を表示
echo ""
echo "=== Normalized Ixy Table ==="
cat $OUTPUT_FILE

echo ""
echo "=== Files created ==="
echo "Text format: $OUTPUT_FILE"
echo "CSV format: $CSV_FILE"
echo "LaTeX format: $LATEX_FILE"

# 簡易統計情報
echo ""
echo "=== Summary Statistics ==="
echo "Average normalized Ixy by condition:"
for condition in "${conditions[@]}"; do
    # CSVファイルから該当列を抽出して平均を計算
    col_index=$(($(echo "${conditions[@]}" | tr ' ' '\n' | grep -n "^$condition$" | cut -d: -f1) + 1))
    avg=$(tail -n +2 $CSV_FILE | cut -d, -f$col_index | grep -v "N/A" | \
          awk '{if($1!="N/A" && $1!="") {sum+=$1; count++}} END {if(count>0) printf "%.6E", sum/count; else print "N/A"}')
    printf "  %-20s: %s\n" "$condition" "$avg"
done

echo ""
echo "=== Processing complete ==="