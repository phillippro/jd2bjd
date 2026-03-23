# JD ↔ BJD 转换工具

将儒略日 (JD) 转换为质心儒略日 (BJD-TDB)，考虑地球轨道运动和恒星运动。

## 功能

- ✅ JD → BJD-TDB 转换
- ✅ 地球质心位置计算
- ✅ 光行时 (LTT) 修正
- ✅ 多观测站支持
- ✅ 视差和自行修正（可选）
- ✅ 批量转换 API

## 安装

```bash
pip install -r requirements.txt
```

## 运行

```bash
python app.py
```

然后浏览器访问 http://localhost:5000

## API 使用

```bash
# 单次转换
curl -X POST http://localhost:5000/convert \
  -H "Content-Type: application/json" \
  -d '{
    "jd": 2460000.5,
    "ra": "83.633",
    "dec": "-5.391",
    "observatory": "Shanghai",
    "parallax": 76.29
  }'

# 批量转换
curl -X POST http://localhost:5000/batch \
  -H "Content-Type: application/json" \
  -d '{
    "jd_list": [2460000.5, 2460001.5, 2460002.5],
    "ra": "83.633",
    "dec": "-5.391"
  }'
```

## 计算原理

1. **地球轨道修正**: 使用 astropy 计算地球在太阳系质心坐标系中的位置
2. **光行时修正**: 计算从观测者到质心的光行时间
3. **坐标变换**: 从地心坐标系转换到质心坐标系

```
BJD_TDB = JD_UTC + ΔT_UTC→TDB + Δt_LTT
```

其中:
- ΔT_UTC→TDB: UTC → TDB 时间转换 (~分钟级)
- Δt_LTT: 光行时修正 (~±16分钟，取决于地球位置)