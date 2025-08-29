#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BER/FER Performance Comparison Plot
Legacy 2D Lattice vs Dec2 3D Lattice Decoder
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# Set Japanese font
plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['font.size'] = 12

# Performance data
decoders = ['Legacy\n(2D Lattice)', 'Dec2\n(3D Lattice)']
ber_values = [7.13, 5.64]
accuracy_values = [92.87, 94.36]

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# BER Comparison
bars1 = ax1.bar(decoders, ber_values, color=['#ff6b6b', '#4ecdc4'], alpha=0.8, width=0.6)
ax1.set_ylabel('Bit Error Rate (%)', fontsize=12, fontweight='bold')
ax1.set_title('BER Performance Comparison', fontsize=14, fontweight='bold')
ax1.set_ylim(0, 8)
ax1.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, value in zip(bars1, ber_values):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
             f'{value:.2f}%', ha='center', va='bottom', fontweight='bold')

# Accuracy Comparison
bars2 = ax2.bar(decoders, accuracy_values, color=['#ff6b6b', '#4ecdc4'], alpha=0.8, width=0.6)
ax2.set_ylabel('Accuracy (%)', fontsize=12, fontweight='bold')
ax2.set_title('Accuracy Performance Comparison', fontsize=14, fontweight='bold')
ax2.set_ylim(90, 96)
ax2.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, value in zip(bars2, accuracy_values):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.05,
             f'{value:.2f}%', ha='center', va='bottom', fontweight='bold')

# Add improvement annotations
# BER improvement
improvement_ber = (ber_values[0] - ber_values[1]) / ber_values[0] * 100
ax1.annotate(f'{improvement_ber:.1f}% improvement', 
            xy=(1, ber_values[1]), xytext=(0.5, ber_values[0] - 1),
            arrowprops=dict(arrowstyle='->', color='green', lw=2),
            fontsize=11, fontweight='bold', color='green', ha='center')

# Accuracy improvement
improvement_acc = accuracy_values[1] - accuracy_values[0]
ax2.annotate(f'+{improvement_acc:.2f}% points', 
            xy=(1, accuracy_values[1]), xytext=(0.5, accuracy_values[1] + 0.5),
            arrowprops=dict(arrowstyle='->', color='green', lw=2),
            fontsize=11, fontweight='bold', color='green', ha='center')

# Add technical details as text box
technical_details = """
Technical Details:
• Block length N = 50
• Codebook: ICB/ex01-new-4gen (4-ary extended)
• Channel parameters: Pi=Pd=Ps=0.004
• Test blocks: 100,000
• Legacy: 2D lattice (drift states only)
• Dec2: 3D lattice (drift + error state memory)
"""

fig.text(0.02, 0.02, technical_details, fontsize=10, 
         bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8),
         verticalalignment='bottom')

# Overall title
fig.suptitle('Dec2 3D Lattice Decoder Performance Evaluation\nError State Memory P(e_{i+1}|e_i) Implementation', 
             fontsize=16, fontweight='bold', y=0.95)

plt.tight_layout()
plt.subplots_adjust(top=0.85, bottom=0.25)

# Save the plot
plt.savefig('decoder_performance_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print("=== Performance Comparison Summary ===")
print(f"Legacy Decoder (2D): BER = {ber_values[0]:.2f}%, Accuracy = {accuracy_values[0]:.2f}%")
print(f"Dec2 Decoder (3D):   BER = {ber_values[1]:.2f}%, Accuracy = {accuracy_values[1]:.2f}%")
print(f"BER Improvement: {improvement_ber:.1f}%")
print(f"Accuracy Improvement: +{improvement_acc:.2f} percentage points")
print("\n✅ Dec2 decoder shows significant performance improvement!")
print("📊 Plot saved as: decoder_performance_comparison.png")