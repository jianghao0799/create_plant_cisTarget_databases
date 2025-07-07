#!/usr/bin/env python3
import os
from pathlib import Path

def generate_cb_files():
    # 创建输出目录
    output_dir = Path("motif_dir")
    output_dir.mkdir(exist_ok=True)
    
    # 读取motifs_id.txt文件
    with open("motifs_id.txt", "r") as f:
        motif_lines = f.read().strip().split('\n')
    
    processed_count = 0
    missing_count = 0
    
    for line in motif_lines:
        line = line.strip()
        if not line:
            continue
            
        # 提取motif_id (例如: M00775_3.00)
        parts = line.split('_')
        if len(parts) >= 3:
            motif_id = f"{parts[0]}_{parts[1]}"  # M00775_3.00
        else:
            print(f"Warning: Invalid format in line: {line}")
            continue
        
        # 检查对应的PWM文件是否存在
        pwm_file = Path(f"pwms_all_motifs/{motif_id}.txt")
        
        if pwm_file.exists():
            # 创建输出文件
            output_file = output_dir / f"{line}.cb"
            
            try:
                with open(pwm_file, "r") as infile:
                    lines = infile.readlines()
                
                # 写入输出文件
                with open(output_file, "w") as outfile:
                    # 写入头部
                    outfile.write(f">{line}\n")
                    
                    # 处理PWM数据：跳过第一行，删除第一列
                    for pwm_line in lines[1:]:  # 跳过第一行表头
                        if pwm_line.strip():
                            # 分割并去掉第一列（位置信息）
                            cols = pwm_line.strip().split('\t')
                            if len(cols) > 1:
                                data_cols = cols[1:]  # 去掉第一列
                                outfile.write('\t'.join(data_cols) + '\n')
                
                print(f"Generated: {output_file}")
                processed_count += 1
                
            except Exception as e:
                print(f"Error processing {pwm_file}: {e}")
                
        else:
            print(f"Warning: PWM file not found for {motif_id}")
            missing_count += 1
    
    print(f"\nSummary:")
    print(f"Successfully processed: {processed_count} files")
    print(f"Missing PWM files: {missing_count}")
    print(f"All .cb files generated in {output_dir}/")

if __name__ == "__main__":
    generate_cb_files()
