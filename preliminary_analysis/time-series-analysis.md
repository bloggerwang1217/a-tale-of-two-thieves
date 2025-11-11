# Time Series Analysis

1. Drum-to-Drum Variance 是什麼？
2. AR(1) 自相關是考慮前一個和自己嗎？
3. Durbin-Watson test 是什麼？ 

我用你們的 Tablet 資料來解釋這三個概念：

## 1. Drum-to-Drum Variance（鼓間變異）

情境回顧：
- 你們有 30 個 drums（桶），隨機抽 10 桶
- 每桶抽 3 個 tablets 測量

Drum-to-Drum Variance 就是「桶和桶之間的差異」

分解變異來源：
總變異 = 桶間變異 + 桶內變異

桶間變異 (Between-Drum Variance) = 0.4394
→ 不同桶的平均值差異
→ 例如：Drum 1 平均 35.5, Drum 5 平均 36.2

桶內變異 (Within-Drum Variance) = 1.3729
→ 同一桶內 3 個 tablets 的差異
→ 例如：Drum 1 的三個值是 35.77, 39.44, 36.43

實務意義：
- 桶間變異小 (0.4394) = 每桶的品質很一致 ✓
- CV = 1.85% = 變異係數很低，品質控制良好 ✓
- 對比 Location Variance (0.9976) =
壓錠後的一致性比混合階段好很多！

### CV (Coefficient of Variation) 計算方式：

```r
# 1. 計算 between-drum 標準差
sd_drum <- sqrt(drum_variance)
sd_drum <- sqrt(0.4394) = 0.6629 mg/100mg

# 2. 計算總平均數
mean_tablet <- mean(tablet_data$ASSAY) = 35.82 mg/100mg

# 3. 計算變異係數
cv_drum <- (sd_drum / mean_tablet) × 100
cv_drum <- (0.6629 / 35.82) × 100 = 1.85%
```

**公式：**
$$\text{CV} = \frac{\text{SD}_{\text{between-drum}}}{\text{Mean}_{\text{total}}} \times 100\%$$

**為什麼用 between-drum SD？**
- 這個 CV 專門用來衡量「drum-to-drum」的相對變異性
- 分子：between-drum variance 的標準差（0.6629）
- 分母：總平均數（35.82）
- 目的：單獨評估桶與桶之間的一致性

**如果要計算 total CV（包含 within-drum + between-drum）：**
```r
total_sd <- sqrt(drum_variance + residual_variance)
          = sqrt(0.4394 + 1.3729)
          = sqrt(1.8123)
          = 1.346 mg/100mg
total_cv <- (1.346 / 35.82) × 100 = 3.76%
```

---
## 2. AR(1) 自相關（Autoregressive Model of Order 1）

AR(1) = 只考慮「前一個」觀測值對「現在」的影響

### 數學表達：

$y_t = \mu + \rho × (y_{t-1} - \mu) + \epsilon_t$

$y_t$ = 第 t 個觀測值
$y_{t-1}$ = 前一個觀測值
$\rho$ = 自相關係數（你們的結果 $\rho$ = 0.251）

### 實際例子（你們的 Tablet 資料）：

時間序列：

Tablet 1: 35.77   ← 第 1 個
Tablet 2: 39.44   ← 第 2 個（前一個是 35.77）
Tablet 3: 36.43   ← 第 3 個（前一個是 39.44）
...
Tablet 30: 35.60  ← 第 30 個（前一個是 36.63）

AR(1) 的意思：
- Tablet 2 的值會不會受到 Tablet 1 影響？
- Tablet 3 的值會不會受到 Tablet 2 影響？
- 只看前一個，不看更早的（例如不看 Tablet 1 對 Tablet 3
的影響）

你們的結果：$\rho$ = 0.251
- 0 = 完全沒相關（獨立）
- 1 = 完全相關（前一個完全決定下一個）
- 0.251 = 弱正相關，表示前一個 tablet
值稍微會影響下一個，但影響很小

為什麼要檢查？
如果 tablets 是按順序從同一批壓出來的：
- 可能機器狀態會連續（例如溫度慢慢升高）
- 前後 tablets 可能比較接近
- 但你們的結果顯示：沒有顯著相關 ✓

---
## 3. Durbin-Watson Test

目的：檢測殘差（residuals）是否有自相關

DW 統計量範圍：

DW = 0    → 完全正自相關（前一個大，下一個也大）
DW = 2    → 沒有自相關（獨立） ← 理想狀態
DW = 4    → 完全負自相關（前一個大，下一個就小）

你們的結果：DW = 1.565, p = 0.077

判讀：
- 1.565 稍微 < 2 → 輕微正自相關的跡象
- p = 0.077 > 0.05 → 但不顯著
- 結論：雖然有一點點正相關傾向，但統計上可以當作獨立

---
AR(1) vs Durbin-Watson：差在哪？

| 比較項目 | Durbin-Watson Test | AR(1) Model  |
|------|--------------------|--------------|
| 性質   | 診斷工具               | 建模方法         |
| 時機   | 模型跑完後檢查            | 建模時直接納入|
| 功能   | 告訴你「有沒有」自相關        |直接「處理」自相關    |
| 結果   | DW 統計量 + p-value   | 估計 ρ + 調整標準誤 |

### 實務應用邏輯：

步驟 1：先用 DW test 檢查是否有自相關
    ↓
    DW ≈ 2 或 p > 0.05 → 沒問題，用一般模型
    DW 偏離 2 且 p < 0.05 → 有自相關，需要處理
    ↓
步驟 2：如果有自相關，改用 AR(1) model
    → 直接在模型中估計 ρ
    → 調整標準誤和信賴區間
    → 得到更準確的推論

### 你們的情況：

Durbin-Watson Test:
- DW = 1.565, p = 0.077
- 結論：無顯著自相關 ✓

AR(1) Model:
- $\rho$ = 0.251, LRT p = 0.235
- LRT 顯示：加入 AR(1) 結構沒有顯著改善模型 ✓

兩個結果一致！
→ 你們的 tablet 資料可以當作獨立觀測值
→ 不需要複雜的時間序列模型
→ 標準的統計方法就夠用了

---
## 總結：你們的 Tablet 資料

1. Drum-to-Drum Variance = 0.4394
    - 每桶品質很一致，CV 只有 1.85%
2. AR(1) ρ = 0.251 (p = 0.235)
    - 前後 tablets 沒有顯著相關
3. Durbin-Watson = 1.565 (p = 0.077)
    - 殘差沒有顯著自相關

結論：壓錠和裝桶的品質控制做得很好！✓