"""Raw (station-adjusted) triangle closures, Colorado-Utah, Schott 1900 pp. 555-566.
No sea-level reduction, no figure adjustment. Fit closure = K * plane area."""
import math
# station: {target: (deg, min, sea-level reduction ")}  -- pp. 555-559
D = {
 'Tushar': {'Mount Nebo':(155,33,.155),'Wasatch':(182,45,.228),'Mount Ellen':(238,41,-.102)},
 'Mount Nebo': {'Patmos Head':(99,26,-.096),'Wasatch':(155,13,-.137),'Tushar':(194,36,.155)},
 'Patmos Head': {'Tavaputs':(98,42,.012),'Mount Waas':(149,29,-.243),'Mount Ellen':(207,9,.120),'Wasatch':(257,55,.161),'Mount Nebo':(297,5,-.125)},
 'Wasatch': {'Mount Nebo':(56,15,-.149),'Patmos Head':(141,19,.147),'Mount Ellen':(228,22,-.186),'Tushar':(302,49,.244)},
 'Mount Waas': {'Mount Ellen':(57,49,.138),'Patmos Head':(124,44,-.199),'Tavaputs':(175,33,.060),'Treasury Mountain':(239,14,.150),'Uncompahgre':(273,50,-.172)},
 'Mount Ellen': {'Tushar':(121,30,-.117),'Wasatch':(171,6,-.184),'Patmos Head':(213,17,.105),'Mount Waas':(268,43,.156),'Uncompahgre':(287,44,0.0)},
 'Treasury Mountain': {'Mount Elbert':(137,13,.115),'Mount Ouray':(189,27,-.278),'Uncompahgre':(255,51,.161),'Mount Waas':(313,40,.126),'Tavaputs':(349,2,-.114)},
 'Tavaputs': {'Treasury Mountain':(87,15,-.164),'Uncompahgre':(118,24,-.280),'Mount Waas':(168,13,.083),'Patmos Head':(246,38,.012)},
 'Mount Elbert': {'Bison':(176,0,.082),'Pikes Peak':(199,22,-.131),'Mount Ouray':(261,34,-.132),'Uncompahgre':(313,14,.278),'Treasury Mountain':(354,19,.104)},
 'Uncompahgre': {'Mount Ellen':(17,57,-.017),'Mount Waas':(34,57,-.158),'Tavaputs':(66,53,-.177),'Treasury Mountain':(122,33,.153),'Mount Elbert':(142,52,.286),'Mount Ouray':(175,40,.186)},
 'Mount Ouray': {'Uncompahgre':(73,31,.184),'Treasury Mountain':(134,1,-.273),'Mount Elbert':(169,2,-.138),'Bison':(217,35,.238),'Pikes Peak':(248,16,.219),'Plateau':(273,44,-.003)},
 'Pikes Peak': {'Plateau':(0,24,-.108),'Mount Ouray':(107,11,.210),'Mount Elbert':(145,46,-.143),'Bison':(179,36,-.250),'Divide':(281,54,.128),'Big Springs':(319,1,-.035)},
 'Bison': {'Pikes Peak':(8,5,-.281),'Mount Ouray':(84,58,.263),'Mount Elbert':(130,53,.089),'Divide':(331,53,-.060)},
 'Plateau': {'Pikes Peak':(0,0,-.286),'Big Springs':(73,43,.118),'Mount Ouray':(312,14,-.019)},
 'Divide': {'Big Springs':(33,19,-.114),'Pikes Peak':(126,59,.240),'Bison':(168,29,-.107)},
 'Big Springs': {'Divide':(33,35,-.137),'Plateau':(279,28,.101),'Pikes Peak':(344,22,-.083)},
}
# No, (A,B,C), printed angle degrees at A,B,C, observed closure ", spherical excess ", sides (m)  -- pp. 564-566
T = [
 (1,('Wasatch','Tushar','Mount Nebo'),(113,27,39),22.325,21.643,(164256.13,81792.99,113592.39)),
 (2,('Mount Ellen','Tushar','Wasatch'),(49,55,74),34.111,34.332,(113592.39,123556.31,143690.32)),
 (3,('Patmos Head','Wasatch','Mount Nebo'),(39,85,55),22.417,22.153,(81792.99,129043.98,107094.02)),
 (4,('Patmos Head','Mount Ellen','Wasatch'),(50,42,87),32.472,33.549,(123556.31,107094.02,159283.42)),
 (5,('Mount Waas','Mount Ellen','Patmos Head'),(66,55,57),48.581,48.709,(159283.42,142552.33,146290.19)),
 (6,('Tavaputs','Mount Waas','Patmos Head'),(78,50,50),30.352,31.622,(142552.33,112788.37,112740.81)),
 (7,('Treasury Mountain','Mount Waas','Tavaputs'),(35,63,80),49.531,49.345,(112740.81,174569.80,192356.71)),
 (8,('Uncompahgre','Mount Waas','Tavaputs'),(31,98,49),45.494,46.152,(112740.81,211050.40,162926.74)),
 (9,('Uncompahgre','Mount Waas','Treasury Mountain'),(87,34,57),44.951,45.184,(192356.71,109320.60,162926.74)),
 (10,('Treasury Mountain','Uncompahgre','Tavaputs'),(93,55,31),48.988,48.378,(211050.40,174569.80,109320.60)),
 (11,('Uncompahgre','Mount Ellen','Mount Waas'),(17,19,143),36.064,35.592,(146290.19,162926.74,294111.77)),
 (12,('Mount Ouray','Uncompahgre','Treasury Mountain'),(60,53,66),26.159,25.554,(109320.60,100470.89,115110.46)),
 (13,('Mount Elbert','Uncompahgre','Treasury Mountain'),(41,20,118),15.420,14.064,(109320.60,57730.33,146049.96)),
 (14,('Mount Elbert','Mount Ouray','Treasury Mountain'),(92,35,52),13.527,11.639,(100470.89,57730.33,79504.24)),
 (15,('Mount Ouray','Uncompahgre','Mount Elbert'),(95,32,51),24.266,23.128,(146049.96,79504.24,115110.46)),
 (16,('Bison','Mount Ouray','Mount Elbert'),(45,48,85),15.692,16.695,(79504.24,82965.96,110377.41)),
 (17,('Pikes Peak','Mount Ouray','Mount Elbert'),(38,79,62),20.933,22.362,(79504.24,125253.81,112774.75)),
 (18,('Pikes Peak','Mount Ouray','Bison'),(72,30,76),16.471,16.130,(110377.41,59101.18,112774.75)),
 (19,('Bison','Pikes Peak','Mount Elbert'),(122,33,23),11.230,10.463,(125253.81,82965.96,59101.18)),
 (20,('Plateau','Mount Ouray','Pikes Peak'),(47,25,106),18.705,17.954,(112774.75,65495.32,145858.19)),
 (21,('Divide','Pikes Peak','Bison'),(41,102,36),6.770,7.721,(59101.18,87141.54,52673.50)),
 (22,('Big Springs','Plateau','Pikes Peak'),(64,73,41),9.941,7.630,(65495.32,69423.13,47805.01)),
 (23,('Big Springs','Pikes Peak','Divide'),(49,37,93),5.415,5.602,(52673.50,41981.60,69423.13)),
]
RHO = 180 * 3600 / math.pi
rows = []
for no, st, degs, obs, se, s in T:
    dsl = 0.0
    for i, A in enumerate(st):
        B, C = st[(i + 1) % 3], st[(i + 2) % 3]
        (db, mb, sb), (dc, mc, sc) = D[A][B], D[A][C]
        a = ((dc + mc / 60) - (db + mb / 60)) % 360
        sl = sc - sb
        if a > 180: a, sl = 360 - a, -sl
        assert abs(a - degs[i]) < 1.01, (no, A, a, degs[i])   # orientation check against printed angle
        dsl += sl
    raw = obs - dsl
    a, b, c = s; p = (a + b + c) / 2
    area = math.sqrt(p * (p - a) * (p - b) * (p - c))          # plane (Heron) area, m^2
    rows.append((no, st, raw, obs, se, area))

def fit(y):  # least squares through origin: y = k * area
    k = sum(r[5] * v for r, v in zip(rows, y)) / sum(r[5] ** 2 for r in rows)
    res = [v - k * r[5] for r, v in zip(rows, y)]
    sig = math.sqrt(sum(e * e for e in res) / (len(rows) - 1) / sum(r[5] ** 2 for r in rows))
    return k, sig, res

def fit2(y):  # with intercept: y = c + k * area
    n = len(rows); X = [r[5] for r in rows]; mx, my = sum(X) / n, sum(y) / n
    sxx = sum((x - mx) ** 2 for x in X); k = sum((x - mx) * (v - my) for x, v in zip(X, y)) / sxx
    c = my - k * mx; res = [v - c - k * x for x, v in zip(X, y)]
    s2 = sum(e * e for e in res) / (n - 2)
    return c, math.sqrt(s2 * (1 / n + mx * mx / sxx)), k, math.sqrt(s2 / sxx)

if __name__ == '__main__':
    print('No | Triangle | Area km2 | Raw closure | After sea level | Schott SE | Raw - SE')
    for no, st, raw, obs, se, area in rows:
        print(f"{no} | {' - '.join(st)} | {area/1e6:,.0f} | {raw:.3f} | {obs:.3f} | {se:.3f} | {raw-se:+.3f}")
    for lab, y in [('raw', [r[2] for r in rows]), ('sea level', [r[3] for r in rows]), ('Schott SE', [r[4] for r in rows])]:
        k, sk, res = fit(y)
        K = k / RHO; R = 1 / math.sqrt(K) / 1000; sR = R * sk / k / 2
        c, sc, k2, sk2 = fit2(y)
        rms = math.sqrt(sum(e * e for e in res) / len(res))
        print(f'{lab}: k = {k*1e10:.4f} "/1e4 km2 +- {sk*1e10:.4f}  R = {R:,.0f} +- {sR:,.0f} km  rms {rms:.3f}"  '
              f'| with intercept: c = {c:+.3f} +- {sc:.3f}", slope {k2*1e10:.4f} +- {sk2*1e10:.4f}')
    d = [r[2] - r[4] for r in rows]
    m = sum(d) / len(d); sd = math.sqrt(sum((x - m) ** 2 for x in d) / (len(d) - 1))
    print(f'raw - SE: mean {m:+.3f}, sd {sd:.3f}, min {min(d):+.3f}, max {max(d):+.3f}')
    print('min raw closure', min(r[2] for r in rows))
