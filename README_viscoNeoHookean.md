# ViscoNeoHookeanElastic 材料模型使用说明

## 概述
`viscoNeoHookeanElastic` 是一个简化的粘弹性-超弹性本构模型，结合了 Neo-Hookean 超弹性行为和简化的粘弹性松弛。

## 参数配置

### 在 `constant/mechanicalProperties` 文件中添加以下参数：

```cpp
// 材料模型类型
type viscoNeoHookeanElastic;

// 超弹性参数 (Neo-Hookean 模型)
mu0             1e6;        // 主剪切模量 [Pa]
mu1             1e5;        // 次剪切模量 [Pa] 
mu2             1e4;        // 第三剪切模量 [Pa]

// 粘弹性参数
relaxationTime  1.0;        // 松弛时间 [s]

// 传统线性弹性参数 (兼容性)
E               1e7;        // 杨氏模量 [Pa]
nu              0.3;        // 泊松比 [-]

// 密度
rho
{
    type        uniform;
    value       1000;       // 密度 [kg/m³]
}

// 几何类型
planeStress     false;      // false = 平面应变, true = 平面应力
```

## 参数说明

### 超弹性参数
- **`mu0`**: 主剪切模量，控制材料的主要刚度
- **`mu1`**: 次剪切模量，影响非线性行为
- **`mu2`**: 第三剪切模量，提供额外的非线性响应

### 粘弹性参数
- **`relaxationTime`**: 松弛时间，控制粘弹性响应的速度
  - 较小值 (0.1-1.0 s): 快速松弛，接近弹性行为
  - 较大值 (10-100 s): 慢速松弛，明显的粘弹性行为

### 传统参数 (兼容性)
- **`E`**: 杨氏模量，用于数值稳定性
- **`nu`**: 泊松比，用于数值稳定性
- **`rho`**: 材料密度

## 材料行为

### 超弹性部分 (Neo-Hookean)
```
τ_hyper = μ₀(B - I) + (μ₁ + μ₂)(I₁ - 3)I
```
其中：
- `B` = 左 Cauchy-Green 张量
- `I₁` = 第一不变量
- `I` = 单位张量

### 粘弹性部分
```
s_new = α·s_old + (1-α)·τ_hyper
```
其中：
- `α = exp(-dt/relaxationTime)`
- `dt` = 时间步长

### 总应力
```
σ = (τ_hyper + s) / J
```
其中：
- `J` = 变形梯度行列式

## 使用建议

1. **参数选择**：
   - 从线性弹性参数开始 (`E`, `nu`)
   - 根据材料特性调整超弹性参数
   - 通过松弛时间控制粘弹性程度

2. **数值稳定性**：
   - 确保 `mu0` 足够大以提供数值稳定性
   - 松弛时间不应过小 (建议 > 0.1s)

3. **边界条件**：
   - 使用 `fixedValue` 固定位移
   - 使用 `tractionDisplacement` 混合边界条件

## 示例配置

### 软材料 (如橡胶)
```cpp
type viscoNeoHookeanElastic;
mu0             1e5;        // 较小的剪切模量
mu1             1e4;
mu2             1e3;
relaxationTime  0.5;         // 快速松弛
```

### 硬材料 (如金属)
```cpp
type viscoNeoHookeanElastic;
mu0             1e9;        // 较大的剪切模量
mu1             1e8;
mu2             1e7;
relaxationTime  10.0;        // 慢速松弛
```
