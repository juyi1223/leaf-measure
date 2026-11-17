/* 香榧叶片测量 —— 核心算法（无外部依赖，浏览器内离线运行）
 * 流程：找四角定位方块 → 透视拉正 → 叶片/纸张分割 → 每片叶子沿中脉测量
 * 坐标约定：
 *   原图坐标：像素 i 覆盖 [i, i+1)，中心在 i+0.5
 *   底板坐标：毫米（mm），与打印底板 1:1（名义尺寸，打印缩放用 sx/sy 校正）
 *   拉正图坐标：像素中心在整数处，mm = region.x0 + (px + 0.5) / S
 */
(function (root) {
  'use strict';

  const BOARD = {
    markers: [[20, 35], [190, 35], [190, 215], [20, 215]], // 左上 右上 右下 左下（中心，mm）
    markerArea: [128, 144, 144, 144],                     // 左上带 4×4 白心
    region: { x0: 27, y0: 44, x1: 183, y1: 172 },         // 拉正、分析的区域
    grid: { x0: 30, y0: 45, cw: 50, rh: 12, cols: 3, rows: 10 },
    num: { dx0: 0.5, dx1: 5.4, dy0: 4.0, dy1: 7.9 },      // 每格编号数字所在的小框（相对格子左上角）
  };

  /* ---------- 线性代数 ---------- */
  function solveLinear(A, b) {
    const n = b.length;
    const M = A.map((row, i) => row.slice().concat([b[i]]));
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      if (Math.abs(M[p][c]) < 1e-12) throw new Error('矩阵奇异');
      const t = M[c]; M[c] = M[p]; M[p] = t;
      for (let r = c + 1; r < n; r++) {
        const f = M[r][c] / M[c][c];
        if (f) for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    const x = new Array(n);
    for (let r = n - 1; r >= 0; r--) {
      let s = M[r][n];
      for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k];
      x[r] = s / M[r][r];
    }
    return x;
  }

  // 求 src(底板 mm) → dst(图像 px) 的单应矩阵
  function homography(src, dst) {
    const A = [], b = [];
    for (let k = 0; k < 4; k++) {
      const [x, y] = src[k], [X, Y] = dst[k];
      A.push([x, y, 1, 0, 0, 0, -x * X, -y * X]); b.push(X);
      A.push([0, 0, 0, x, y, 1, -x * Y, -y * Y]); b.push(Y);
    }
    const h = solveLinear(A, b);
    return [h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], 1];
  }
  function applyH(H, x, y) {
    const w = H[6] * x + H[7] * y + H[8];
    return [(H[0] * x + H[1] * y + H[2]) / w, (H[3] * x + H[4] * y + H[5]) / w];
  }

  // 多项式最小二乘 v = p(u)
  function fitPoly(us, vs, deg) {
    const n = us.length;
    let um = 0; for (let i = 0; i < n; i++) um += us[i]; um /= n;
    let sc = 1e-9; for (let i = 0; i < n; i++) sc = Math.max(sc, Math.abs(us[i] - um));
    const m = deg + 1;
    const A = [], b = new Array(m).fill(0);
    for (let r = 0; r < m; r++) A.push(new Array(m).fill(0));
    const pw = new Array(2 * m);
    for (let i = 0; i < n; i++) {
      const t = (us[i] - um) / sc;
      pw[0] = 1; for (let k = 1; k < 2 * m; k++) pw[k] = pw[k - 1] * t;
      for (let r = 0; r < m; r++) { b[r] += pw[r] * vs[i]; for (let c = 0; c < m; c++) A[r][c] += pw[r + c]; }
    }
    for (let r = 0; r < m; r++) A[r][r] += 1e-9;
    const co = solveLinear(A, b);
    return {
      f(u) { const t = (u - um) / sc; let s = 0, p = 1; for (let k = 0; k < m; k++) { s += co[k] * p; p *= t; } return s; },
      d(u) { const t = (u - um) / sc; let s = 0, p = 1; for (let k = 1; k < m; k++) { s += k * co[k] * p; p *= t; } return s / sc; },
    };
  }

  function quantile(arr, q) {
    if (!arr.length) return 0;
    const a = Float32Array.from(arr); a.sort();
    return a[Math.min(a.length - 1, Math.max(0, Math.floor(q * (a.length - 1))))];
  }

  /* ---------- 图像采样 ---------- */
  // 原图双线性采样（输入原图连续坐标），结果写入 out[0..2]
  function sampleRGB(img, X, Y, out) {
    const W = img.width, Hh = img.height, d = img.data;
    let x = X - 0.5, y = Y - 0.5;
    if (x < 0) x = 0; else if (x > W - 1.001) x = W - 1.001;
    if (y < 0) y = 0; else if (y > Hh - 1.001) y = Hh - 1.001;
    const i0 = x | 0, j0 = y | 0, fx = x - i0, fy = y - j0;
    const p00 = (j0 * W + i0) * 4, p10 = p00 + 4, p01 = p00 + W * 4, p11 = p01 + 4;
    const w00 = (1 - fx) * (1 - fy), w10 = fx * (1 - fy), w01 = (1 - fx) * fy, w11 = fx * fy;
    out[0] = d[p00] * w00 + d[p10] * w10 + d[p01] * w01 + d[p11] * w11;
    out[1] = d[p00 + 1] * w00 + d[p10 + 1] * w10 + d[p01 + 1] * w01 + d[p11 + 1] * w11;
    out[2] = d[p00 + 2] * w00 + d[p10 + 2] * w10 + d[p01 + 2] * w01 + d[p11 + 2] * w11;
  }

  /* ---------- 连通域 ---------- */
  function labelComponents(mask, w, h, conn8, labels) {
    labels = labels || new Int32Array(w * h);
    const stack = new Int32Array(w * h);
    const comps = [null];
    let id = 0;
    const N = w * h;
    for (let p0 = 0; p0 < N; p0++) {
      if (!mask[p0] || labels[p0]) continue;
      id++;
      let sp = 0; stack[sp++] = p0; labels[p0] = id;
      let area = 0, sx = 0, sy = 0, minx = w, maxx = -1, miny = h, maxy = -1;
      while (sp) {
        const p = stack[--sp];
        const y = (p / w) | 0, x = p - y * w;
        area++; sx += x; sy += y;
        if (x < minx) minx = x; if (x > maxx) maxx = x; if (y < miny) miny = y; if (y > maxy) maxy = y;
        for (let dy = -1; dy <= 1; dy++) {
          const yy = y + dy; if (yy < 0 || yy >= h) continue;
          for (let dx = -1; dx <= 1; dx++) {
            if (!dx && !dy) continue;
            if (!conn8 && dx && dy) continue;
            const xx = x + dx; if (xx < 0 || xx >= w) continue;
            const q = yy * w + xx;
            if (mask[q] && !labels[q]) { labels[q] = id; stack[sp++] = q; }
          }
        }
      }
      comps.push({ id, area, cx: sx / area, cy: sy / area, minx, maxx, miny, maxy });
    }
    return { labels, comps };
  }

  /* ---------- 1. 找四角定位方块 ---------- */
  function detectMarkers(img) {
    const W = img.width, H = img.height, d = img.data;
    const k = Math.max(1, Math.ceil(Math.max(W, H) / 1000));
    const sw = Math.floor(W / k), sh = Math.floor(H / k);
    const L = new Float32Array(sw * sh), G = new Float32Array(sw * sh);
    const nn = k * k;
    for (let j = 0; j < sh; j++) for (let i = 0; i < sw; i++) {
      let r = 0, g = 0, b = 0;
      for (let yy = 0; yy < k; yy++) {
        let p = ((j * k + yy) * W + i * k) * 4;
        for (let xx = 0; xx < k; xx++, p += 4) { r += d[p]; g += d[p + 1]; b += d[p + 2]; }
      }
      r /= nn; g /= nn; b /= nn;
      L[j * sw + i] = (r + g + b) / 3; G[j * sw + i] = 2 * g - r - b;
    }
    // 局部"纸张亮度"（分块 90 分位），用来适应光照不均
    const nb = 10, bw = Math.ceil(sw / nb), bh = Math.ceil(sh / nb);
    const nbx = Math.ceil(sw / bw), nby = Math.ceil(sh / bh);
    const P = new Float32Array(nbx * nby);
    for (let by = 0; by < nby; by++) for (let bx = 0; bx < nbx; bx++) {
      const vals = [];
      for (let j = by * bh; j < Math.min(sh, (by + 1) * bh); j += 2)
        for (let i = bx * bw; i < Math.min(sw, (bx + 1) * bw); i += 2) vals.push(L[j * sw + i]);
      P[by * nbx + bx] = Math.max(40, quantile(vals, 0.9));
    }
    const paperAt = (i, j) => {
      let fx = (i + 0.5) / bw - 0.5, fy = (j + 0.5) / bh - 0.5;
      fx = Math.max(0, Math.min(nbx - 1, fx)); fy = Math.max(0, Math.min(nby - 1, fy));
      const x0 = Math.min(nbx - 2, fx | 0), y0 = Math.min(nby - 2, fy | 0);
      if (nbx < 2 || nby < 2) return P[0];
      const ax = fx - x0, ay = fy - y0;
      return P[y0 * nbx + x0] * (1 - ax) * (1 - ay) + P[y0 * nbx + x0 + 1] * ax * (1 - ay) +
        P[(y0 + 1) * nbx + x0] * (1 - ax) * ay + P[(y0 + 1) * nbx + x0 + 1] * ax * ay;
    };
    const mask = new Uint8Array(sw * sh);
    for (let j = 0; j < sh; j++) for (let i = 0; i < sw; i++) {
      const p = j * sw + i, pp = paperAt(i, j);
      if (L[p] < 0.5 * pp && G[p] < 0.15 * pp) mask[p] = 1;
    }
    const { comps } = labelComponents(mask, sw, sh, false);
    const cands = [];
    for (let c = 1; c < comps.length; c++) {
      const o = comps[c];
      const bwid = o.maxx - o.minx + 1, bhei = o.maxy - o.miny + 1;
      if (o.area < 25 || bwid < 5 || bhei < 5) continue;
      if (o.minx === 0 || o.miny === 0 || o.maxx === sw - 1 || o.maxy === sh - 1) continue;
      if (Math.max(bwid, bhei) / Math.min(bwid, bhei) > 1.8) continue;
      if (o.area / (bwid * bhei) < 0.38) continue;
      // 中心是否为白（左上角方块有白心）
      const ci = Math.round(o.cx), cj = Math.round(o.cy);
      let s = 0, n = 0;
      for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) {
        const ii = ci + di, jj = cj + dj;
        if (ii >= 0 && jj >= 0 && ii < sw && jj < sh) { s += L[jj * sw + ii]; n++; }
      }
      const hole = s / n > 0.55 * paperAt(ci, cj);
      cands.push({ x: (o.cx + 0.5) * k, y: (o.cy + 0.5) * k, area: o.area * k * k, hole });
    }
    cands.sort((a, b) => b.area - a.area);
    const top = cands.slice(0, 12);
    let best = null;
    const n = top.length;
    for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++) for (let e = c + 1; e < n; e++) {
      const set = [top[a], top[b], top[c], top[e]];
      const holes = set.filter(s => s.hole);
      if (holes.length !== 1) continue;
      const areas = set.map(s => s.area);
      if (Math.max(...areas) / Math.min(...areas) > 4) continue;
      const mx = set.reduce((s, p) => s + p.x, 0) / 4, my = set.reduce((s, p) => s + p.y, 0) / 4;
      const ordered = set.slice().sort((p, q) => Math.atan2(p.y - my, p.x - mx) - Math.atan2(q.y - my, q.x - mx));
      const ti = ordered.indexOf(holes[0]);
      const q = [0, 1, 2, 3].map(i => ordered[(ti + i) % 4]);
      const dist = (p, r) => Math.hypot(p.x - r.x, p.y - r.y);
      const t = dist(q[0], q[1]), r = dist(q[1], q[2]), bo = dist(q[2], q[3]), l = dist(q[3], q[0]);
      if (Math.min(t, r, bo, l) < 1) continue;
      if (t / bo > 2 || bo / t > 2 || r / l > 2 || l / r > 2) continue;
      const asp = (r + l) / (t + bo);
      if (asp < 0.5 || asp > 2.2) continue;
      const size = areas.reduce((s, v) => s + Math.sqrt(v), 0) / 4;
      const rel = size / ((t + bo) / 2);
      if (rel < 0.035 || rel > 0.14) continue;
      let sign = 0, convex = true;
      for (let i = 0; i < 4; i++) {
        const p0 = q[i], p1 = q[(i + 1) % 4], p2 = q[(i + 2) % 4];
        const cr = (p1.x - p0.x) * (p2.y - p1.y) - (p1.y - p0.y) * (p2.x - p1.x);
        const sg = Math.sign(cr);
        if (!sign) sign = sg; else if (sg !== sign) convex = false;
      }
      if (!convex) continue;
      const score = areas.reduce((s, v) => s + v, 0);
      if (!best || score > best.score) best = { score, pts: q.map(p => [p.x, p.y]) };
    }
    if (!best) {
      const holeCnt = cands.filter(c => c.hole).length;
      if (cands.length >= 4 && holeCnt === 0)
        throw new Error('找不到左上角带白心的定位方块。请确认底板四个角的黑方块都完整拍进画面，且没有反光。');
      throw new Error('没有找到底板四角的 4 个黑色定位方块。请让 4 个方块都完整进入画面，不要被手指、阴影或反光挡住。');
    }
    return best.pts;
  }

  /* ---------- 2. 精确定位（在拉正坐标中求方块质心，迭代） ---------- */
  function refineMarkers(img, pts) {
    let H = homography(BOARD.markers, pts);
    const half = 9, sp = 10, n = half * 2 * sp;
    const rgb = new Float64Array(3);
    let info = [];
    for (let it = 0; it < 3; it++) {
      const np = [];
      info = [];
      for (let m = 0; m < 4; m++) {
        const [cx, cy] = BOARD.markers[m];
        const Ls = new Float32Array(n * n);
        for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
          const [X, Y] = applyH(H, cx - half + (i + 0.5) / sp, cy - half + (j + 0.5) / sp);
          sampleRGB(img, X, Y, rgb);
          Ls[j * n + i] = (rgb[0] + rgb[1] + rgb[2]) / 3;
        }
        const Pp = quantile(Ls, 0.85), Dk = quantile(Ls, 0.03);
        const span = Math.max(1, Pp - Dk);
        let sw = 0, sx = 0, sy = 0;
        for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
          let w = (Pp - Ls[j * n + i]) / span;
          w = w < 0.15 ? 0 : w > 0.85 ? 1 : (w - 0.15) / 0.7;
          if (w) { sw += w; sx += w * (cx - half + (i + 0.5) / sp); sy += w * (cy - half + (j + 0.5) / sp); }
        }
        if (sw < 1) throw new Error('定位方块对比度太低，请在光线更均匀的地方重拍。');
        const c = [sx / sw, sy / sw];
        np.push(applyH(H, c[0], c[1]));
        info.push({ area: sw / (sp * sp), contrast: span, offset: Math.hypot(c[0] - cx, c[1] - cy) });
      }
      pts = np;
      H = homography(BOARD.markers, pts);
    }
    return { H, pts, info };
  }

  /* ---------- 3. 拉正 + 计算"叶片分数" ---------- */
  function rectify(img, H, S, Sd) {
    const R = BOARD.region, rgb = new Float64Array(3);
    const W = Math.round((R.x1 - R.x0) * S), Hh = Math.round((R.y1 - R.y0) * S);
    // 纸张亮度图：8mm 分块的 90 分位，再做 3×3 中值
    const B = 8, nbx = Math.ceil((R.x1 - R.x0) / B), nby = Math.ceil((R.y1 - R.y0) / B);
    const pm0 = new Float32Array(nbx * nby);
    for (let by = 0; by < nby; by++) for (let bx = 0; bx < nbx; bx++) {
      const vals = [];
      for (let j = 0; j < 16; j++) for (let i = 0; i < 16; i++) {
        const x = R.x0 + bx * B + (i + 0.5) * B / 16, y = R.y0 + by * B + (j + 0.5) * B / 16;
        if (x > R.x1 || y > R.y1) continue;
        const [X, Y] = applyH(H, x, y);
        sampleRGB(img, X, Y, rgb); vals.push((rgb[0] + rgb[1] + rgb[2]) / 3);
      }
      pm0[by * nbx + bx] = Math.max(30, quantile(vals, 0.9));
    }
    const pm = new Float32Array(nbx * nby);
    for (let by = 0; by < nby; by++) for (let bx = 0; bx < nbx; bx++) {
      const v = [];
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const xx = bx + dx, yy = by + dy;
        if (xx >= 0 && yy >= 0 && xx < nbx && yy < nby) v.push(pm0[yy * nbx + xx]);
      }
      v.sort((a, b) => a - b); pm[by * nbx + bx] = v[v.length >> 1];
    }
    const paperAt = (x, y) => {
      let fx = (x - R.x0) / B - 0.5, fy = (y - R.y0) / B - 0.5;
      fx = Math.max(0, Math.min(nbx - 1, fx)); fy = Math.max(0, Math.min(nby - 1, fy));
      const x0 = Math.min(nbx - 2, fx | 0), y0 = Math.min(nby - 2, fy | 0);
      const ax = fx - x0, ay = fy - y0;
      return pm[y0 * nbx + x0] * (1 - ax) * (1 - ay) + pm[y0 * nbx + x0 + 1] * ax * (1 - ay) +
        pm[(y0 + 1) * nbx + x0] * (1 - ax) * ay + pm[(y0 + 1) * nbx + x0 + 1] * ax * ay;
    };
    const g = BOARD.grid, nm = BOARD.num;
    const sc = new Uint8Array(W * Hh);
    for (let j = 0; j < Hh; j++) {
      const y = R.y0 + (j + 0.5) / S;
      const gy = y - g.y0, ry = gy - Math.floor(gy / g.rh) * g.rh;
      const rowNum = gy >= 0 && gy < g.rh * g.rows && ry >= nm.dy0 && ry < nm.dy1;
      for (let i = 0; i < W; i++) {
        const x = R.x0 + (i + 0.5) / S;
        const [X, Y] = applyH(H, x, y);
        sampleRGB(img, X, Y, rgb);
        const P = paperAt(x, y);
        const L = (rgb[0] + rgb[1] + rgb[2]) / 3, exg = (2 * rgb[1] - rgb[0] - rgb[2]) / P;
        let s = 1 - L / P + 0.5 * Math.max(0, exg);
        if (rowNum) {
          const gx = x - g.x0;
          if (gx >= 0 && gx < g.cw * g.cols) {
            const rx = gx - Math.floor(gx / g.cw) * g.cw;
            if (rx >= nm.dx0 && rx < nm.dx1 && exg < 0.08) s = 0; // 印刷的灰色编号，不算叶片
          }
        }
        const v = s * 200;
        sc[j * W + i] = v < 0 ? 0 : v > 255 ? 255 : v;
      }
    }
    // 显示用的彩色拉正图
    const Wd = Math.round((R.x1 - R.x0) * Sd), Hd = Math.round((R.y1 - R.y0) * Sd);
    const disp = new Uint8ClampedArray(Wd * Hd * 4);
    for (let j = 0; j < Hd; j++) for (let i = 0; i < Wd; i++) {
      const [X, Y] = applyH(H, R.x0 + (i + 0.5) / Sd, R.y0 + (j + 0.5) / Sd);
      sampleRGB(img, X, Y, rgb);
      const p = (j * Wd + i) * 4;
      disp[p] = rgb[0]; disp[p + 1] = rgb[1]; disp[p + 2] = rgb[2]; disp[p + 3] = 255;
    }
    return { sc, W, Hh, S, disp: { data: disp, width: Wd, height: Hd, S: Sd } };
  }

  /* ---------- 4. 测量单片叶子 ---------- */
  function measureLeaf(ctx, comp) {
    const { sc, labels, W, Hh, S, paperLev } = ctx;
    const id = comp.id;
    const ox = comp.minx, oy = comp.miny, bw = comp.maxx - ox + 1, bh = comp.maxy - oy + 1;
    // 填洞（叶面反光造成的亮斑）
    const loc = new Uint8Array(bw * bh);
    for (let y = 0; y < bh; y++) for (let x = 0; x < bw; x++) if (labels[(y + oy) * W + x + ox] === id) loc[y * bw + x] = 1;
    const st = new Int32Array(bw * bh); let sp = 0;
    const pushIf = (x, y) => { const p = y * bw + x; if (!loc[p]) { loc[p] = 2; st[sp++] = p; } };
    for (let x = 0; x < bw; x++) { pushIf(x, 0); pushIf(x, bh - 1); }
    for (let y = 0; y < bh; y++) { pushIf(0, y); pushIf(bw - 1, y); }
    while (sp) {
      const p = st[--sp], y = (p / bw) | 0, x = p - y * bw;
      if (x > 0) pushIf(x - 1, y); if (x < bw - 1) pushIf(x + 1, y);
      if (y > 0) pushIf(x, y - 1); if (y < bh - 1) pushIf(x, y + 1);
    }
    let n = 0, mx = 0, my = 0;
    const vals = [];
    for (let y = 0; y < bh; y++) for (let x = 0; x < bw; x++) {
      const p = y * bw + x, q = (y + oy) * W + x + ox;
      if (loc[p] === 0 && labels[q] === 0) { labels[q] = id; loc[p] = 1; }
      if (loc[p] === 1 && labels[q] === id) { n++; mx += x + ox; my += y + oy; vals.push(sc[q]); }
    }
    mx /= n; my /= n;
    const leafLev = quantile(vals, 0.6);
    const Te = 0.5 * (leafLev + paperLev);
    // 主方向
    let cxx = 0, cyy = 0, cxy = 0;
    for (let y = 0; y < bh; y++) for (let x = 0; x < bw; x++) {
      if (labels[(y + oy) * W + x + ox] !== id) continue;
      const dx = x + ox - mx, dy = y + oy - my;
      cxx += dx * dx; cyy += dy * dy; cxy += dx * dy;
    }
    const th = 0.5 * Math.atan2(2 * cxy, cxx - cyy);
    let e1 = [Math.cos(th), Math.sin(th)];
    if (e1[0] < 0) e1 = [-e1[0], -e1[1]];
    const e2 = [-e1[1], e1[0]];
    // 按 u 分箱求中线初值
    let umin = Infinity, umax = -Infinity;
    const us = [], vs = [];
    for (let y = 0; y < bh; y++) for (let x = 0; x < bw; x++) {
      if (labels[(y + oy) * W + x + ox] !== id) continue;
      const dx = x + ox - mx, dy = y + oy - my;
      const u = dx * e1[0] + dy * e1[1], v = dx * e2[0] + dy * e2[1];
      us.push(u); vs.push(v);
      if (u < umin) umin = u; if (u > umax) umax = u;
    }
    const nbin = Math.floor(umax - umin) + 1;
    const vmin = new Float32Array(nbin).fill(Infinity), vmax = new Float32Array(nbin).fill(-Infinity);
    for (let i = 0; i < us.length; i++) {
      const b = Math.min(nbin - 1, Math.floor(us[i] - umin));
      if (vs[i] < vmin[b]) vmin[b] = vs[i]; if (vs[i] > vmax[b]) vmax[b] = vs[i];
    }
    const fu = [], fv = [];
    for (let b = Math.floor(nbin * 0.08); b < Math.ceil(nbin * 0.92); b++)
      if (vmin[b] <= vmax[b]) { fu.push(umin + b + 0.5); fv.push(0.5 * (vmin[b] + vmax[b])); }
    let poly = fitPoly(fu.length >= 4 ? fu : us, fu.length >= 4 ? fv : vs, fu.length >= 8 ? 2 : 1);

    const inside = (x, y) => {
      const xi = Math.round(x), yi = Math.round(y);
      return xi >= 0 && yi >= 0 && xi < W && yi < Hh && labels[yi * W + xi] === id;
    };
    const samp = (x, y) => {
      if (x < 0) x = 0; else if (x > W - 1.001) x = W - 1.001;
      if (y < 0) y = 0; else if (y > Hh - 1.001) y = Hh - 1.001;
      const i0 = x | 0, j0 = y | 0, fx = x - i0, fy = y - j0, p = j0 * W + i0;
      return sc[p] * (1 - fx) * (1 - fy) + sc[p + 1] * fx * (1 - fy) + sc[p + W] * (1 - fx) * fy + sc[p + W + 1] * fx * fy;
    };
    const maxR = 4 * S;
    const edge = (px, py, nx, ny, dir) => {
      let r = 0;
      while (r < maxR) {
        const rr = r + 0.5;
        if (!inside(px + dir * nx * rr, py + dir * ny * rr)) break;
        r = rr;
      }
      if (r >= maxR) return r;
      let prevR = r + 3, prevV = samp(px + dir * nx * prevR, py + dir * ny * prevR);
      for (let rr = r + 2.9; rr >= Math.max(0, r - 3); rr -= 0.1) {
        const v = samp(px + dir * nx * rr, py + dir * ny * rr);
        if (v >= Te && prevV < Te) return rr + (v - Te) / (v - prevV) * (prevR - rr);
        prevV = v; prevR = rr;
      }
      return r + 0.25;
    };

    const uB = umin - 0.5;
    let uT = umax + 0.5;
    let N, du;
    const stations = (poly) => {
      N = Math.max(10, Math.ceil((uT - uB) / (0.05 * S)) + 1);
      du = (uT - uB) / (N - 1);
      const P = new Float64Array(2 * N), E = new Float64Array(4 * N), w = new Float64Array(N), nrm = new Float64Array(2 * N);
      for (let k = 0; k < N; k++) {
        const u = uB + k * du, v = poly.f(u), dv = poly.d(u);
        const px = mx + u * e1[0] + v * e2[0], py = my + u * e1[1] + v * e2[1];
        let tx = e1[0] + dv * e2[0], ty = e1[1] + dv * e2[1];
        const tl = Math.hypot(tx, ty); tx /= tl; ty /= tl;
        const nx = -ty, ny = tx;
        P[2 * k] = px; P[2 * k + 1] = py; nrm[2 * k] = nx; nrm[2 * k + 1] = ny;
        if (inside(px, py)) {
          const rp = edge(px, py, nx, ny, 1), rn = edge(px, py, nx, ny, -1);
          E[4 * k] = px + nx * rp; E[4 * k + 1] = py + ny * rp;
          E[4 * k + 2] = px - nx * rn; E[4 * k + 3] = py - ny * rn;
          w[k] = rp + rn;
        } else {
          E[4 * k] = E[4 * k + 2] = px; E[4 * k + 1] = E[4 * k + 3] = py;
        }
      }
      return { P, E, w, nrm };
    };
    let st1 = stations(poly);
    // 用左右边缘中点重新拟合中线
    let wmax0 = 0; for (let k = 0; k < N; k++) wmax0 = Math.max(wmax0, st1.w[k]);
    const ru = [], rv = [];
    for (let k = 0; k < N; k++) {
      if (st1.w[k] < 0.4 * wmax0) continue;
      const cx = 0.5 * (st1.E[4 * k] + st1.E[4 * k + 2]) - mx, cy = 0.5 * (st1.E[4 * k + 1] + st1.E[4 * k + 3]) - my;
      ru.push(cx * e1[0] + cy * e1[1]); rv.push(cx * e2[0] + cy * e2[1]);
    }
    if (ru.length >= 12) poly = fitPoly(ru, rv, (uT - uB) / S > 10 && ru.length >= 30 ? 3 : 2);
    // 叶尖：沿中线向外找，用较低的阈值（细尖在模糊下对比度会变低）
    {
      const levTip = paperLev + 0.3 * (leafLev - paperLev);
      let found = null, prevU = null, prevV = null;
      for (let u = uT - 1.0 * S; u <= uT + 1.0 * S; u += 0.1) {
        const v = poly.f(u), dv = poly.d(u);
        const px = mx + u * e1[0] + v * e2[0], py = my + u * e1[1] + v * e2[1];
        let tx = e1[0] + dv * e2[0], ty = e1[1] + dv * e2[1];
        const tl = Math.hypot(tx, ty), nx = -ty / tl, ny = tx / tl;
        let best = 0;
        for (let r = -1.5; r <= 1.5; r += 0.5) best = Math.max(best, samp(px + nx * r, py + ny * r));
        if (prevV !== null && prevV >= levTip && best < levTip) { found = prevU + (prevV - levTip) / (prevV - best) * (u - prevU); break; }
        prevU = u; prevV = best;
      }
      if (found !== null && found > uT - 0.6 * S) uT = found;
    }
    st1 = stations(poly);
    const { P, E, w } = st1;

    // 阴影/翘起检查、模糊度
    const off = 0.3 * S;
    let sPos = 0, sNeg = 0, cntS = 0;
    let wm = 0; for (let k = 0; k < N; k++) wm = Math.max(wm, w[k]);
    const blurs = [];
    const lev25 = paperLev + 0.25 * (leafLev - paperLev), lev75 = paperLev + 0.75 * (leafLev - paperLev);
    for (let k = 0; k < N; k++) {
      if (w[k] < 0.6 * wm) continue;
      const nx = st1.nrm[2 * k], ny = st1.nrm[2 * k + 1];
      const ax = E[4 * k] + nx * off, ay = E[4 * k + 1] + ny * off;
      const bx = E[4 * k + 2] - nx * off, by = E[4 * k + 3] - ny * off;
      const la = labels[Math.round(ay) * W + Math.round(ax)], lb = labels[Math.round(by) * W + Math.round(bx)];
      if (la === 0 && lb === 0) { sPos += samp(ax, ay); sNeg += samp(bx, by); cntS++; }
      if (k % 7 === 0) {
        // 边缘 25%→75% 过渡宽度
        let r25 = null, r75 = null, prev = samp(E[4 * k] + nx * 3, E[4 * k + 1] + ny * 3);
        for (let r = 2.9; r >= -3; r -= 0.1) {
          const v = samp(E[4 * k] + nx * r, E[4 * k + 1] + ny * r);
          if (r25 === null && v >= lev25 && prev < lev25) r25 = r;
          if (r75 === null && v >= lev75 && prev < lev75) { r75 = r; break; }
          prev = v;
        }
        if (r25 !== null && r75 !== null) blurs.push((r25 - r75) / S);
      }
    }
    const shadow = cntS ? (Math.max(sPos, sNeg) / cntS - paperLev) / 200 : 0;

    const toMM = (px, py) => [BOARD.region.x0 + (px + 0.5) / S, BOARD.region.y0 + (py + 0.5) / S];
    const ptsNom = new Float64Array(2 * N), edgeNom = new Float64Array(4 * N);
    for (let k = 0; k < N; k++) {
      const a = toMM(P[2 * k], P[2 * k + 1]); ptsNom[2 * k] = a[0]; ptsNom[2 * k + 1] = a[1];
      const b = toMM(E[4 * k], E[4 * k + 1]), c = toMM(E[4 * k + 2], E[4 * k + 3]);
      edgeNom[4 * k] = b[0]; edgeNom[4 * k + 1] = b[1]; edgeNom[4 * k + 2] = c[0]; edgeNom[4 * k + 3] = c[1];
    }
    const bb = [toMM(comp.minx - 0.5, comp.miny - 0.5), toMM(comp.maxx + 0.5, comp.maxy + 0.5)];
    return {
      N, ptsNom, edgeNom, bboxMM: [bb[0][0], bb[0][1], bb[1][0], bb[1][1]],
      shadow, blur: blurs.length ? quantile(blurs, 0.5) : null,
      reversed: false, cutMode: 'auto', cutIdx: 0, deleted: false, flags: [],
    };
  }

  /* ---------- 5. 由测量数据算长度/宽度（可反复调用：改切点、翻转、改校准） ---------- */
  function movAvg(a, half) {
    const n = a.length, out = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      let s = 0, c = 0;
      for (let j = Math.max(0, k - half); j <= Math.min(n - 1, k + half); j++) { s += a[j]; c++; }
      out[k] = s / c;
    }
    return out;
  }
  function reverseLeaf(leaf) {
    const N = leaf.N, p = leaf.ptsNom, e = leaf.edgeNom;
    const np = new Float64Array(2 * N), ne = new Float64Array(4 * N);
    for (let k = 0; k < N; k++) {
      const j = N - 1 - k;
      np[2 * k] = p[2 * j]; np[2 * k + 1] = p[2 * j + 1];
      ne[4 * k] = e[4 * j + 2]; ne[4 * k + 1] = e[4 * j + 3]; ne[4 * k + 2] = e[4 * j]; ne[4 * k + 3] = e[4 * j + 1];
    }
    leaf.ptsNom = np; leaf.edgeNom = ne; leaf.reversed = !leaf.reversed;
    leaf.cutMode = 'auto';
  }
  function compute(leaf, opt) {
    const sx = opt.sx || 1, sy = opt.sy || 1, frac = opt.frac == null ? 0.3 : opt.frac;
    const N = leaf.N, p = leaf.ptsNom, e = leaf.edgeNom;
    const w = new Float64Array(N), arc = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      w[k] = Math.hypot((e[4 * k] - e[4 * k + 2]) * sx, (e[4 * k + 1] - e[4 * k + 3]) * sy);
      if (k) arc[k] = arc[k - 1] + Math.hypot((p[2 * k] - p[2 * k - 2]) * sx, (p[2 * k + 1] - p[2 * k - 1]) * sy);
    }
    const wS = movAvg(w, 3), wC = movAvg(w, 2);
    let wmax = 0; for (let k = 0; k < N; k++) wmax = Math.max(wmax, wS[k]);
    // 自动叶柄切点：从叶柄端起，宽度第一次达到最大宽度 frac 处
    let auto = 0;
    for (let k = 0; frac > 0 && k < N; k++) {
      if (wC[k] >= frac * wmax) {
        auto = k === 0 ? 0 : k - 1 + (frac * wmax - wC[k - 1]) / Math.max(1e-9, wC[k] - wC[k - 1]);
        break;
      }
    }
    leaf.cutAuto = auto;
    if (leaf.cutMode === 'auto') leaf.cutIdx = auto;
    const c = Math.max(0, Math.min(N - 2, leaf.cutIdx));
    const k0 = Math.floor(c), f = c - k0;
    const cx = p[2 * k0] + f * (p[2 * k0 + 2] - p[2 * k0]), cy = p[2 * k0 + 1] + f * (p[2 * k0 + 3] - p[2 * k0 + 1]);
    const arcC = arc[k0] + f * (arc[k0 + 1] - arc[k0]);
    let maxIdx = Math.ceil(c), bw = 0;
    for (let k = Math.ceil(c); k < N; k++) if (wS[k] > bw) { bw = wS[k]; maxIdx = k; }
    // 最大宽度取 ±0.15mm 内的平均（降噪）
    leaf.curved = arc[N - 1] - arcC;
    leaf.straight = Math.hypot((p[2 * N - 2] - cx) * sx, (p[2 * N - 1] - cy) * sy);
    leaf.width = bw;
    leaf.maxIdx = maxIdx;
    leaf.wmaxAll = wmax;
    // 方向检查：叶柄端（起点）到 70% 宽度的距离应短于叶尖端
    let kb = 0; while (kb < N && wS[kb] < 0.7 * wmax) kb++;
    let kt = N - 1; while (kt > 0 && wS[kt] < 0.7 * wmax) kt--;
    leaf.dirSuspect = frac > 0 && (arc[kb] > 1.2 * (arc[N - 1] - arc[kt]) && arc[kb] > 1.5);
    return leaf;
  }

  /* ---------- 主流程 ---------- */
  async function process(img, opt, onStep) {
    opt = opt || {};
    const step = async (msg) => { if (onStep) { onStep(msg); await new Promise(r => setTimeout(r, 0)); } };
    const t0 = Date.now();
    const warnings = [];
    await step('寻找定位方块…');
    const rough = detectMarkers(img);
    const ref = refineMarkers(img, rough);
    const H = ref.H;
    ref.info.forEach((m, i) => {
      const exp = BOARD.markerArea[i];
      if (m.area < exp * 0.8 || m.area > exp * 1.25)
        warnings.push(`第 ${i + 1} 个定位方块看起来不完整（被挡住或有反光），尺寸可能不准，建议重拍。`);
    });
    // 每毫米像素数
    const R = BOARD.region;
    const cx = (R.x0 + R.x1) / 2, cy = (R.y0 + R.y1) / 2;
    const a = applyH(H, cx, cy), b = applyH(H, cx + 1, cy), c2 = applyH(H, cx, cy + 1);
    const native = Math.sqrt(Math.hypot(b[0] - a[0], b[1] - a[1]) * Math.hypot(c2[0] - a[0], c2[1] - a[1]));
    const S = Math.max(10, Math.min(opt.maxS || 30, native));
    if (native < 12) warnings.push(`照片分辨率偏低（每毫米约 ${native.toFixed(1)} 像素），建议手机靠近一些或打开高像素模式。`);

    await step('拉正照片…');
    const rect = rectify(img, H, S, opt.dispS || 6);
    await step('识别叶片…');
    const { sc, W, Hh } = rect;
    // 阈值：纸张水平与叶片水平的中点
    const hist = new Uint32Array(256);
    for (let i = 0; i < sc.length; i++) hist[sc[i]]++;
    const histQ = (lo, q) => { let tot = 0; for (let v = lo; v < 256; v++) tot += hist[v]; let acc = 0; for (let v = lo; v < 256; v++) { acc += hist[v]; if (acc >= q * tot) return v; } return 255; };
    const paperLev = histQ(0, 0.5);
    const leafLev = histQ(70, 0.5);
    const T = Math.max(40, 0.5 * (paperLev + leafLev));
    const mask = new Uint8Array(W * Hh);
    for (let i = 0; i < sc.length; i++) mask[i] = sc[i] > T ? 1 : 0;
    const { labels, comps } = labelComponents(mask, W, Hh, true);
    const g = BOARD.grid;
    const cells = new Map();
    const minArea = 3 * S * S;
    const areas = [];
    for (let k = 1; k < comps.length; k++) {
      const o = comps[k];
      if (o.area < minArea) continue;
      const xm = R.x0 + (o.cx + 0.5) / S, ym = R.y0 + (o.cy + 0.5) / S;
      const col = Math.floor((xm - g.x0) / g.cw), row = Math.floor((ym - g.y0) / g.rh);
      if (col < 0 || col >= g.cols || row < 0 || row >= g.rows) continue;
      const no = col * g.rows + row + 1;
      if (!cells.has(no)) cells.set(no, []);
      cells.get(no).push(o);
      areas.push(o.area);
    }
    const medArea = quantile(areas, 0.5);
    await step('测量中…');
    const ctx = { sc, labels, W, Hh, S, paperLev };
    const leaves = [];
    for (const [no, list] of cells) {
      list.sort((p, q) => q.area - p.area);
      const o = list[0];
      const leaf = measureLeaf(ctx, o);
      leaf.no = no;
      const col = Math.floor((no - 1) / g.rows), row = (no - 1) % g.rows;
      const flags = [];
      if (list.length > 1) flags.push('格内有多个物体');
      const touches = o.minx === 0 || o.miny === 0 || o.maxx === W - 1 || o.maxy === Hh - 1;
      const [bx0, by0, bx1, by1] = leaf.bboxMM;
      const cellX0 = g.x0 + col * g.cw, cellY0 = g.y0 + row * g.rh;
      if (touches) flags.push('超出拍摄区域');
      else if (o.area > 2.2 * medArea || by0 < cellY0 - 3 || by1 > cellY0 + g.rh + 3 || bx1 > cellX0 + g.cw + 4) flags.push('可能和别的叶子粘连');
      if (leaf.shadow > 0.12) flags.push('边缘有阴影/可能翘起');
      leaf.flags = flags;
      compute(leaf, opt);
      if (leaf.dirSuspect) leaf.flags.push('叶柄方向可疑');
      leaves.push(leaf);
    }
    leaves.sort((p, q) => p.no - q.no);
    // 长度异常
    const lens = leaves.map(l => l.curved);
    const mL = quantile(lens, 0.5);
    leaves.forEach(l => { if (l.curved < 0.6 * mL || l.curved > 1.5 * mL) l.flags.push('长度和同组差别大'); });
    const blurs = leaves.map(l => l.blur).filter(v => v != null);
    const blur = blurs.length ? quantile(blurs, 0.5) : null;
    if (blur != null && blur > 0.13) warnings.push(`照片偏模糊（边缘过渡约 ${blur.toFixed(2)} mm），建议对好焦、拿稳后重拍。`);
    return { H, S, native, leaves, warnings, blur, disp: rect.disp, ms: Date.now() - t0, markers: ref.pts };
  }

  // 从原图按底板坐标渲染一块区域（用于放大查看单片叶子）
  function renderRegion(img, H, x0, y0, x1, y1, S) {
    const W = Math.max(1, Math.round((x1 - x0) * S)), Hh = Math.max(1, Math.round((y1 - y0) * S));
    const out = new Uint8ClampedArray(W * Hh * 4), rgb = new Float64Array(3);
    for (let j = 0; j < Hh; j++) for (let i = 0; i < W; i++) {
      const [X, Y] = applyH(H, x0 + (i + 0.5) / S, y0 + (j + 0.5) / S);
      sampleRGB(img, X, Y, rgb);
      const p = (j * W + i) * 4;
      out[p] = rgb[0]; out[p + 1] = rgb[1]; out[p + 2] = rgb[2]; out[p + 3] = 255;
    }
    return { data: out, width: W, height: Hh };
  }

  root.LeafMeasure = { BOARD, process, compute, reverseLeaf, renderRegion, homography, applyH };
})(typeof window !== 'undefined' ? window : globalThis);
