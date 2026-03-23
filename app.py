"""
JD → BJD 转换服务
使用 Astropy 进行高精度天体测量计算

考虑:
1. 地球轨道运动 (质心修正)
2. 恒星视差和自行 (可选)
3. 光行时效应
"""

from flask import Flask, render_template_string, request, jsonify
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, get_body
from astropy import units as u
import math

app = Flask(__name__)

# 默认观测站点 (上海天文台)
DEFAULT_OBSERVATORY = {
    'name': 'Shanghai',
    'lon': 121.54,  # 经度
    'lat': 31.22,   # 纬度
    'height': 7     # 海拔 (米)
}

# 著名观测站
OBSERVATORIES = {
    'Shanghai': {'lon': 121.54, 'lat': 31.22, 'height': 7},
    'Mauna Kea': {'lon': -155.48, 'lat': 19.82, 'height': 4207},
    'Paranal': {'lon': -70.40, 'lat': -24.63, 'height': 2635},
    'La Silla': {'lon': -70.73, 'lat': -29.26, 'height': 2400},
    'McDonald': {'lon': -104.02, 'lat': 30.67, 'height': 2070},
    'Siding Spring': {'lon': 149.04, 'lat': -31.27, 'height': 1164},
    'Cerro Tololo': {'lon': -70.81, 'lat': -30.17, 'height': 2215},
}


def jd_to_bjd(jd, ra, dec, obs_name='Shanghai', parallax=0, pmra=0, pmdec=0, rv_star=0):
    """
    将 JD 转换为 BJD
    
    参数:
    - jd: 儒略日 (可以是单个值或列表)
    - ra: 赤经 (度或 'HH:MM:SS' 格式)
    - dec: 赤纬 (度或 'DD:MM:SS' 格式)
    - obs_name: 观测站名称
    - parallax: 视差 (mas)
    - pmra: 赤经自行 (mas/yr)
    - pmdec: 赤纬自行 (mas/yr)
    - rv_star: 恒星视向速度 (km/s)
    
    返回:
    - BJD_TDB (质心儒略日，质心坐标系)
    """
    
    # 获取观测站位置
    obs = OBSERVATORIES.get(obs_name, DEFAULT_OBSERVATORY)
    location = EarthLocation(lon=obs['lon']*u.deg, 
                            lat=obs['lat']*u.deg, 
                            height=obs['height']*u.m)
    
    # 创建目标坐标
    if isinstance(ra, str):
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
    else:
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    
    # 添加视差和自行
    if parallax > 0:
        coord.distance = (1000/parallax) * u.pc
    if pmra != 0 or pmdec != 0:
        coord.pm_ra_cosdec = pmra * u.mas/u.yr
        coord.pm_dec = pmdec * u.mas/u.yr
    
    # 转换为时间对象
    if isinstance(jd, (list, tuple)):
        t = Time(jd, format='jd', scale='utc', location=location)
    else:
        t = Time(jd, format='jd', scale='utc', location=location)
    
    # 转换为 TDB (质心动力学时) - 这已经包含了地球轨道运动的修正
    t_tdb = t.tdb
    
    # 计算光行时 (Light Travel Time Correction)
    # 这是由于地球-质心距离变化引起的
    ltt = t_tdb.light_travel_time(coord)
    
    # BJD_TDB = TDB + 光行时修正
    bjd_tdb = t_tdb + ltt
    
    return {
        'bjd_tdb': bjd_tdb.value,
        'tdb': t_tdb.value,
        'ltt_correction_days': ltt.to(u.day).value,
        'ltt_correction_sec': ltt.to(u.second).value,
    }


def get_earth_position(jd):
    """获取地球在质心坐标系中的位置"""
    t = Time(jd, format='jd', scale='utc')
    earth = get_body('earth', t)
    return {
        'x_au': earth.x.to(u.AU).value,
        'y_au': earth.y.to(u.AU).value,
        'z_au': earth.z.to(u.AU).value,
    }


@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)


@app.route('/convert', methods=['POST'])
def convert():
    data = request.json
    
    try:
        # 解析输入
        jd = float(data.get('jd', 0))
        ra = data.get('ra', '0')
        dec = data.get('dec', '0')
        obs_name = data.get('observatory', 'Shanghai')
        parallax = float(data.get('parallax', 0))
        pmra = float(data.get('pmra', 0))
        pmdec = float(data.get('pmdec', 0))
        
        if jd == 0:
            return jsonify({'error': '请输入有效的 JD'})
        
        # 转换
        result = jd_to_bjd(jd, ra, dec, obs_name, parallax, pmra, pmdec)
        
        # 获取地球位置（用于教学展示）
        earth_pos = get_earth_position(jd)
        
        return jsonify({
            'success': True,
            'input': {
                'jd': jd,
                'ra': ra,
                'dec': dec,
                'observatory': obs_name,
                'parallax': parallax,
                'pmra': pmra,
                'pmdec': pmdec,
            },
            'result': {
                'bjd_tdb': result['bjd_tdb'],
                'tdb': result['tdb'],
                'ltt_correction_days': result['ltt_correction_days'],
                'ltt_correction_sec': result['ltt_correction_sec'],
            },
            'earth_position_au': earth_pos,
        })
        
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/batch', methods=['POST'])
def batch_convert():
    """批量转换多个 JD"""
    data = request.json
    
    try:
        jd_list = data.get('jd_list', [])
        ra = data.get('ra', '0')
        dec = data.get('dec', '0')
        obs_name = data.get('observatory', 'Shanghai')
        
        results = []
        for jd in jd_list:
            result = jd_to_bjd(float(jd), ra, dec, obs_name)
            results.append({
                'jd': jd,
                'bjd_tdb': result['bjd_tdb'],
                'ltt_correction_sec': result['ltt_correction_sec'],
            })
        
        return jsonify({'success': True, 'results': results})
        
    except Exception as e:
        return jsonify({'error': str(e)})


HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>JD → BJD 转换器</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            min-height: 100vh;
            padding: 20px;
            color: #eee;
        }
        .container {
            max-width: 800px;
            margin: 0 auto;
            background: rgba(255,255,255,0.05);
            border-radius: 16px;
            padding: 30px;
            backdrop-filter: blur(10px);
        }
        h1 {
            text-align: center;
            margin-bottom: 10px;
            background: linear-gradient(90deg, #00d4ff, #7b2ff7);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .subtitle {
            text-align: center;
            color: #888;
            margin-bottom: 30px;
            font-size: 14px;
        }
        .form-group {
            margin-bottom: 20px;
        }
        label {
            display: block;
            margin-bottom: 8px;
            color: #aaa;
            font-size: 14px;
        }
        input, select {
            width: 100%;
            padding: 12px;
            background: rgba(255,255,255,0.1);
            border: 1px solid rgba(255,255,255,0.1);
            border-radius: 8px;
            color: #fff;
            font-size: 16px;
        }
        input:focus, select:focus {
            outline: none;
            border-color: #00d4ff;
        }
        .row {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }
        .btn {
            width: 100%;
            padding: 14px;
            background: linear-gradient(90deg, #00d4ff, #7b2ff7);
            border: none;
            border-radius: 8px;
            color: white;
            font-size: 16px;
            font-weight: bold;
            cursor: pointer;
            transition: transform 0.2s;
        }
        .btn:hover {
            transform: scale(1.02);
        }
        .btn:active {
            transform: scale(0.98);
        }
        .result {
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,212,255,0.1);
            border-radius: 12px;
            border: 1px solid rgba(0,212,255,0.3);
            display: none;
        }
        .result.show { display: block; }
        .result h3 { color: #00d4ff; margin-bottom: 15px; }
        .result-item {
            display: flex;
            justify-content: space-between;
            padding: 10px 0;
            border-bottom: 1px solid rgba(255,255,255,0.1);
        }
        .result-item:last-child { border-bottom: none; }
        .result-label { color: #aaa; }
        .result-value { 
            color: #00d4ff; 
            font-family: 'Monaco', monospace;
            font-size: 18px;
        }
        .info {
            margin-top: 30px;
            padding: 20px;
            background: rgba(255,255,255,0.03);
            border-radius: 12px;
            font-size: 13px;
            color: #888;
            line-height: 1.6;
        }
        .info h4 { color: #aaa; margin-bottom: 10px; }
        .loading {
            text-align: center;
            padding: 20px;
            display: none;
        }
        .loading.show { display: block; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🔭 JD → BJD 转换器</h1>
        <p class="subtitle">儒略日 (JD) → 质心儒略日 (BJD-TDB)</p>
        
        <div class="row">
            <div class="form-group">
                <label>JD (儒略日)</label>
                <input type="number" id="jd" placeholder="例如: 2460000.5" step="0.0001">
            </div>
            <div class="form-group">
                <label>观测站</label>
                <select id="observatory">
                    <option value="Shanghai">上海天文台</option>
                    <option value="Mauna Kea">Mauna Kea (夏威夷)</option>
                    <option value="Paranal">Paranal (智利)</option>
                    <option value="La Silla">La Silla (智利)</option>
                    <option value="McDonald">McDonald (德州)</option>
                    <option value="Siding Spring">Siding Spring (澳洲)</option>
                    <option value="Cerro Tololo">Cerro Tololo (智利)</option>
                </select>
            </div>
        </div>
        
        <div class="row">
            <div class="form-group">
                <label>赤经 RA (度 或 HH:MM:SS)</label>
                <input type="text" id="ra" placeholder="例如: 83.633 或 05:34:32" value="83.633">
            </div>
            <div class="form-group">
                <label>赤纬 Dec (度 或 DD:MM:SS)</label>
                <input type="text" id="dec" placeholder="例如: -5.391 或 -05:23:28" value="-5.391">
            </div>
        </div>
        
        <div class="row">
            <div class="form-group">
                <label>视差 parallax (mas)</label>
                <input type="number" id="parallax" placeholder="例如: 76.29" step="0.01" value="0">
            </div>
            <div class="form-group">
                <label>恒星视向速度 (km/s)</label>
                <input type="number" id="rv" placeholder="例如: 5.2" step="0.1" value="0">
            </div>
        </div>
        
        <div class="row">
            <div class="form-group">
                <label>赤经自行 pmRA (mas/yr)</label>
                <input type="number" id="pmra" placeholder="例如: 3.5" step="0.01" value="0">
            </div>
            <div class="form-group">
                <label>赤纬自行 pmDec (mas/yr)</label>
                <input type="number" id="pmdec" placeholder="例如: -5.8" step="0.01" value="0">
            </div>
        </div>
        
        <button class="btn" onclick="convert()">转换</button>
        
        <div class="loading" id="loading">计算中...</div>
        
        <div class="result" id="result">
            <h3>🎯 转换结果</h3>
            <div class="result-item">
                <span class="result-label">输入 JD (UTC)</span>
                <span class="result-value" id="res-jd"></span>
            </div>
            <div class="result-item">
                <span class="result-label">BJD-TDB (质心)</span>
                <span class="result-value" id="res-bjd"></span>
            </div>
            <div class="result-item">
                <span class="result-label">TDB</span>
                <span class="result-value" id="res-tdb"></span>
            </div>
            <div class="result-item">
                <span class="result-label">光行时修正 (天)</span>
                <span class="result-value" id="res-ltt-d"></span>
            </div>
            <div class="result-item">
                <span class="result-label">光行时修正 (秒)</span>
                <span class="result-value" id="res-ltt-s"></span>
            </div>
            <div class="result-item">
                <span class="result-label">地球位置 X (AU)</span>
                <span class="result-value" id="res-ex"></span>
            </div>
            <div class="result-item">
                <span class="result-label">地球位置 Y (AU)</span>
                <span class="result-value" id="res-ey"></span>
            </div>
            <div class="result-item">
                <span class="result-label">地球位置 Z (AU)</span>
                <span class="result-value" id="res-ez"></span>
            </div>
        </div>
        
        <div class="info">
            <h4>📖 说明</h4>
            <p>• <b>JD</b>: 儒略日，从公元前4713年1月1日正午开始计算的日数</p>
            <p>• <b>BJD-TDB</b>: 质心儒略日，以太阳系质心为参考点的时间</p>
            <p>• <b>光行时修正</b>: 地球到质心距离变化导致的时间差，最大约 ±16分钟</p>
            <p>• <b>视差</b>: 输入视差可以获得更精确的恒星距离修正</p>
            <p>• <b>自行</b>: 恒星的固有运动，对长时间跨度的观测很重要</p>
        </div>
    </div>
    
    <script>
        async function convert() {
            const jd = document.getElementById('jd').value;
            if (!jd) {
                alert('请输入 JD');
                return;
            }
            
            document.getElementById('loading').classList.add('show');
            document.getElementById('result').classList.remove('show');
            
            try {
                const response = await fetch('/convert', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({
                        jd: parseFloat(jd),
                        ra: document.getElementById('ra').value,
                        dec: document.getElementById('dec').value,
                        observatory: document.getElementById('observatory').value,
                        parallax: parseFloat(document.getElementById('parallax').value) || 0,
                        pmra: parseFloat(document.getElementById('pmra').value) || 0,
                        pmdec: parseFloat(document.getElementById('pmdec').value) || 0,
                    })
                });
                
                const data = await response.json();
                
                if (data.error) {
                    alert(data.error);
                    return;
                }
                
                document.getElementById('res-jd').textContent = data.input.jd.toFixed(6);
                document.getElementById('res-bjd').textContent = data.result.bjd_tdb.toFixed(6);
                document.getElementById('res-tdb').textContent = data.result.tdb.toFixed(6);
                document.getElementById('res-ltt-d').textContent = data.result.ltt_correction_days.toFixed(6);
                document.getElementById('res-ltt-s').textContent = data.result.ltt_correction_sec.toFixed(3);
                document.getElementById('res-ex').textContent = data.earth_position_au.x_au.toFixed(4);
                document.getElementById('res-ey').textContent = data.earth_position_au.y_au.toFixed(4);
                document.getElementById('res-ez').textContent = data.earth_position_au.z_au.toFixed(4);
                
                document.getElementById('result').classList.add('show');
                
            } catch (e) {
                alert('错误: ' + e);
            } finally {
                document.getElementById('loading').classList.remove('show');
            }
        }
        
        // 默认使用当前时间
        document.getElementById('jd').value = Math.floor(Date.now()/86400000 + 2440587.5);
    </script>
</body>
</html>
'''


if __name__ == '__main__':
    print("启动 JD→BJD 转换服务...")
    print("请访问: http://localhost:5000")
    app.run(host='0.0.0.0', port=5000, debug=True)