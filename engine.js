/* Liquid Layers Mod - physics engine (no DOM).
 * Grid/particle (MPM-flavoured) liquid with several density layers, plus rigid
 * bodies that float on it via sampled buoyancy + drag, with reaction on the fluid.
 */
(function (root) {
'use strict';

const REST = 2.0;

// bottom -> top. m = mass per particle, heavier layers sink.
const LAYERS = [
  { m: 2.00, color: '#6f3bd6' },
  { m: 1.65, color: '#2f6fe4' },
  { m: 1.35, color: '#1fb8c9' },
  { m: 1.12, color: '#3ec46d' },
  { m: 1.00, color: '#f2c230' },
];

class Liquid {
  constructor(GX, GY, layers) {
    this.GX = GX; this.GY = GY; this.layers = layers;
    const maxP = GX * GY * 2;
    this.maxP = maxP; this.n = 0;
    this.x = new Float32Array(maxP); this.y = new Float32Array(maxP);
    this.u = new Float32Array(maxP); this.v = new Float32Array(maxP);
    this.mat = new Uint8Array(maxP);
    this.cx = new Int32Array(maxP); this.cy = new Int32Array(maxP);
    this.px = new Float32Array(maxP * 3); this.py = new Float32Array(maxP * 3);
    this.gx = new Float32Array(maxP * 3); this.gy = new Float32Array(maxP * 3);
    const N = GX * GY;
    this.nm = new Float32Array(N); this.nd = new Float32Array(N);
    this.nu = new Float32Array(N); this.nv = new Float32Array(N);
    this.ax = new Float32Array(N); this.ay = new Float32Array(N);
    this.flag = new Uint8Array(N); this.act = new Int32Array(N); this.na = 0;
    this.layerM = Float32Array.from(layers.map(l => l.m));
    this.k = 2.6;        // stiffness
    this.visc = 0.05;    // viscosity
    this.flip = 0.9;     // FLIP/PIC blend
    this.maxSpeed = 0.9; // cells per time unit
    this.fill(0.42);
  }

  // Fill the tank bottom with equal-thickness layers (heaviest at the bottom).
  fill(fraction) {
    const GX = this.GX, GY = this.GY, L = this.layers.length;
    const x0 = 3.6, x1 = GX - 4.6, yBot = GY - 4.6, yTop = yBot - (GY - 8) * fraction;
    const s = Math.SQRT1_2;
    this.n = 0;
    let row = 0;
    for (let y = yBot; y > yTop; y -= s, row++) {
      const frac = (yBot - y) / (yBot - yTop);
      const li = Math.min(L - 1, (frac * L) | 0);
      const off = (row & 1) ? s / 2 : 0;
      for (let x = x0 + off; x < x1; x += s) {
        const p = this.n++;
        this.x[p] = x + (Math.random() - 0.5) * 0.04;
        this.y[p] = y + (Math.random() - 0.5) * 0.04;
        this.u[p] = 0; this.v[p] = 0; this.mat[p] = li;
      }
    }
    for (let a = 0; a < this.na; a++) {
      const i = this.act[a];
      this.nm[i] = this.nd[i] = this.nu[i] = this.nv[i] = this.ax[i] = this.ay[i] = 0; this.flag[i] = 0;
    }
    this.na = 0;
  }

  // Passes 1-3: nodes get mass/velocity, forces, then the node velocity update.
  prepare(grx, gry, dt) {
    const GX = this.GX, GY = this.GY, n = this.n;
    const X = this.x, Y = this.y, U = this.u, V = this.v, MAT = this.mat, CX = this.cx, CY = this.cy;
    const PX = this.px, PY = this.py, GXa = this.gx, GYa = this.gy;
    const nm = this.nm, nd = this.nd, nu = this.nu, nv = this.nv, ax = this.ax, ay = this.ay;
    const flag = this.flag, act = this.act, LM = this.layerM;

    for (let a = 0; a < this.na; a++) {
      const i = act[a];
      nm[i] = nd[i] = nu[i] = nv[i] = ax[i] = ay[i] = 0; flag[i] = 0;
    }
    let na = 0;

    // pass 1: weights + node mass/density/momentum
    for (let p = 0; p < n; p++) {
      const x = X[p], y = Y[p];
      const cx = (x - 0.5) | 0, cy = (y - 0.5) | 0;
      CX[p] = cx; CY[p] = cy;
      const o = p * 3;
      let t = cx - x;
      PX[o] = 0.5 * t * t + 1.5 * t + 1.125; GXa[o] = t + 1.5;
      t += 1; PX[o + 1] = 0.75 - t * t; GXa[o + 1] = -2 * t;
      t += 1; PX[o + 2] = 0.5 * t * t - 1.5 * t + 1.125; GXa[o + 2] = t - 1.5;
      t = cy - y;
      PY[o] = 0.5 * t * t + 1.5 * t + 1.125; GYa[o] = t + 1.5;
      t += 1; PY[o + 1] = 0.75 - t * t; GYa[o + 1] = -2 * t;
      t += 1; PY[o + 2] = 0.5 * t * t - 1.5 * t + 1.125; GYa[o + 2] = t - 1.5;
      const pm = LM[MAT[p]], pu = U[p], pv = V[p];
      for (let i = 0; i < 3; i++) {
        const base = (cx + i) * GY + cy;
        for (let j = 0; j < 3; j++) {
          const idx = base + j;
          if (!flag[idx]) { flag[idx] = 1; act[na++] = idx; }
          const phi = PX[o + i] * PY[o + j];
          nm[idx] += phi * pm; nd[idx] += phi;
          nu[idx] += phi * pm * pu; nv[idx] += phi * pm * pv;
        }
      }
    }
    this.na = na;
    for (let a = 0; a < na; a++) {
      const i = act[a], m = nm[i];
      if (m > 1e-6) { const inv = 1 / m; nu[i] *= inv; nv[i] *= inv; } else { nu[i] = 0; nv[i] = 0; }
    }

    // pass 2: pressure + viscous stress -> node forces
    const k = this.k, mu = this.visc;
    for (let p = 0; p < n; p++) {
      const o = p * 3, cx = CX[p], cy = CY[p];
      let dens = 0, dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
      for (let i = 0; i < 3; i++) {
        const base = (cx + i) * GY + cy;
        for (let j = 0; j < 3; j++) {
          const idx = base + j;
          const phi = PX[o + i] * PY[o + j];
          const gxv = GXa[o + i] * PY[o + j], gyv = PX[o + i] * GYa[o + j];
          dens += phi * nd[idx];
          const uu = nu[idx], vv = nv[idx];
          dudx -= uu * gxv; dudy -= uu * gyv; dvdx -= vv * gxv; dvdy -= vv * gyv;
        }
      }
      let pr = k * (dens - REST) / REST;
      if (pr > 2.5) pr = 2.5; else if (pr < 0) pr *= 0.25;
      const sxx = -pr + 2 * mu * dudx, syy = -pr + 2 * mu * dvdy, sxy = mu * (dudy + dvdx);
      const x = X[p], y = Y[p];
      let fx = 0, fy = 0;
      if (x < 4) fx = 0.6 * (4 - x); else if (x > GX - 5) fx = 0.6 * (GX - 5 - x);
      if (y < 4) fy = 0.6 * (4 - y); else if (y > GY - 5) fy = 0.6 * (GY - 5 - y);
      for (let i = 0; i < 3; i++) {
        const base = (cx + i) * GY + cy;
        for (let j = 0; j < 3; j++) {
          const idx = base + j;
          const phi = PX[o + i] * PY[o + j];
          const gxv = GXa[o + i] * PY[o + j], gyv = PX[o + i] * GYa[o + j];
          ax[idx] += sxx * gxv + sxy * gyv + fx * phi;
          ay[idx] += sxy * gxv + syy * gyv + fy * phi;
        }
      }
    }

    // pass 3: node velocity update (ax/ay become the velocity change)
    for (let a = 0; a < na; a++) {
      const i = act[a], m = nm[i];
      let du = grx * dt, dv = gry * dt;
      if (m > 1e-6) { const inv = dt / m; du += ax[i] * inv; dv += ay[i] * inv; }
      ax[i] = du; ay[i] = dv; nu[i] += du; nv[i] += dv;
    }
  }

  // Pass 4: gather velocities back to particles and move them.
  advect(dt) {
    const GX = this.GX, GY = this.GY, n = this.n;
    const X = this.x, Y = this.y, U = this.u, V = this.v, CX = this.cx, CY = this.cy;
    const PX = this.px, PY = this.py;
    const nu = this.nu, nv = this.nv, ax = this.ax, ay = this.ay;
    const f = this.flip, ms = this.maxSpeed;
    const xl = 3, xh = GX - 4, yl = 3, yh = GY - 4;
    for (let p = 0; p < n; p++) {
      const o = p * 3, cx = CX[p], cy = CY[p];
      let gu = 0, gv = 0, du = 0, dv = 0;
      for (let i = 0; i < 3; i++) {
        const base = (cx + i) * GY + cy;
        for (let j = 0; j < 3; j++) {
          const idx = base + j, phi = PX[o + i] * PY[o + j];
          gu += phi * nu[idx]; gv += phi * nv[idx]; du += phi * ax[idx]; dv += phi * ay[idx];
        }
      }
      let pu = f * (U[p] + du) + (1 - f) * gu, pv = f * (V[p] + dv) + (1 - f) * gv;
      let s = pu * pu + pv * pv;
      if (s > ms * ms) { s = ms / Math.sqrt(s); pu *= s; pv *= s; }
      U[p] = pu; V[p] = pv;
      let mx = gu, my = gv;
      s = mx * mx + my * my;
      if (s > ms * ms) { s = ms / Math.sqrt(s); mx *= s; my *= s; }
      let x = X[p] + mx * dt, y = Y[p] + my * dt;
      if (x < xl) { x = xl; if (U[p] < 0) U[p] = 0; } else if (x > xh) { x = xh; if (U[p] > 0) U[p] = 0; }
      if (y < yl) { y = yl; if (V[p] < 0) V[p] = 0; } else if (y > yh) { y = yh; if (V[p] > 0) V[p] = 0; }
      X[p] = x; Y[p] = y;
    }
  }

  // Read the node fields at a point (call between prepare and advect).
  sample(x, y, out) {
    const GX = this.GX, GY = this.GY;
    if (x < 2.5) x = 2.5; else if (x > GX - 3.5) x = GX - 3.5;
    if (y < 2.5) y = 2.5; else if (y > GY - 3.5) y = GY - 3.5;
    const cx = (x - 0.5) | 0, cy = (y - 0.5) | 0;
    const wx = this._w0, wy = this._w1;
    let t = cx - x;
    wx[0] = 0.5 * t * t + 1.5 * t + 1.125; t += 1; wx[1] = 0.75 - t * t; t += 1; wx[2] = 0.5 * t * t - 1.5 * t + 1.125;
    t = cy - y;
    wy[0] = 0.5 * t * t + 1.5 * t + 1.125; t += 1; wy[1] = 0.75 - t * t; t += 1; wy[2] = 0.5 * t * t - 1.5 * t + 1.125;
    let rho = 0, dens = 0, su = 0, sv = 0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const idx = (cx + i) * GY + cy + j, phi = wx[i] * wy[j], d = this.nd[idx];
        rho += phi * this.nm[idx]; dens += phi * d;
        su += phi * d * this.nu[idx]; sv += phi * d * this.nv[idx];
      }
    }
    out.rho = rho; out.dens = dens;
    out.u = dens > 1e-6 ? su / dens : 0; out.v = dens > 1e-6 ? sv / dens : 0;
  }

  // Apply an equal-and-opposite force (per unit time) on the grid at a point.
  push(x, y, fx, fy, dt) {
    const GX = this.GX, GY = this.GY;
    if (x < 2.5) x = 2.5; else if (x > GX - 3.5) x = GX - 3.5;
    if (y < 2.5) y = 2.5; else if (y > GY - 3.5) y = GY - 3.5;
    const cx = (x - 0.5) | 0, cy = (y - 0.5) | 0;
    const wx = this._w0, wy = this._w1;
    let t = cx - x;
    wx[0] = 0.5 * t * t + 1.5 * t + 1.125; t += 1; wx[1] = 0.75 - t * t; t += 1; wx[2] = 0.5 * t * t - 1.5 * t + 1.125;
    t = cy - y;
    wy[0] = 0.5 * t * t + 1.5 * t + 1.125; t += 1; wy[1] = 0.75 - t * t; t += 1; wy[2] = 0.5 * t * t - 1.5 * t + 1.125;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const idx = (cx + i) * GY + cy + j, m = this.nm[idx];
        if (m < 1e-3) continue;
        const phi = wx[i] * wy[j], inv = phi * dt / m;
        let du = fx * inv, dv = fy * inv;
        if (du > 0.05) du = 0.05; else if (du < -0.05) du = -0.05;
        if (dv > 0.05) dv = 0.05; else if (dv < -0.05) dv = -0.05;
        this.nu[idx] += du; this.nv[idx] += dv; this.ax[idx] += du; this.ay[idx] += dv;
      }
    }
  }

  // Pointer stirring: drag nearby particles along.
  stir(x, y, dx, dy, r) {
    const r2 = r * r, X = this.x, Y = this.y, U = this.u, V = this.v, ms = this.maxSpeed;
    for (let p = 0; p < this.n; p++) {
      const ex = X[p] - x, ey = Y[p] - y, d2 = ex * ex + ey * ey;
      if (d2 > r2) continue;
      const w = 1 - Math.sqrt(d2) / r, ww = w * w * 0.12;
      let u = U[p] + dx * ww, v = V[p] + dy * ww;
      const s = u * u + v * v;
      if (s > ms * ms) { const q = ms / Math.sqrt(s); u *= q; v *= q; }
      U[p] = u; V[p] = v;
    }
  }
}
Liquid.prototype._w0 = new Float32Array(3);
Liquid.prototype._w1 = new Float32Array(3);

/* ---------------- rigid bodies ---------------- */

function polyArea(pts) {
  let a = 0;
  for (let i = 0, n = pts.length; i < n; i++) {
    const p = pts[i], q = pts[(i + 1) % n];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return a / 2;
}
function polyCentroid(pts) {
  let a = 0, cx = 0, cy = 0;
  for (let i = 0, n = pts.length; i < n; i++) {
    const p = pts[i], q = pts[(i + 1) % n], c = p[0] * q[1] - q[0] * p[1];
    a += c; cx += (p[0] + q[0]) * c; cy += (p[1] + q[1]) * c;
  }
  a /= 2;
  return [cx / (6 * a), cy / (6 * a)];
}
function inPoly(pts, x, y) {
  let ins = false;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
    const xi = pts[i][0], yi = pts[i][1], xj = pts[j][0], yj = pts[j][1];
    if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) ins = !ins;
  }
  return ins;
}
function circlePoly(r, n) {
  const pts = [];
  for (let i = 0; i < n; i++) { const a = i / n * Math.PI * 2; pts.push([Math.cos(a) * r, Math.sin(a) * r]); }
  return pts;
}

// dens = mass per area (fluid is ~2.0 at the top layer .. ~4.0 at the bottom).
// com = centre-of-mass offset from the centroid (positive y = lower, makes it self-right).
const SPECS = {
  ball:   { poly: circlePoly(3.2, 18), dens: 0.80, com: [0, 0.0] },
  duck:   { poly: circlePoly(3.0, 18), dens: 0.60, com: [0, 1.1] },
  crate:  { poly: [[-3, -3], [3, -3], [3, 3], [-3, 3]], dens: 1.5, com: [0, 0] },
  boat:   { poly: [[-7, -2.5], [7, -2.5], [4.5, 2.5], [-4.5, 2.5]], mass: 38, com: [0, 1.3] },
  anchor: { poly: [[0, -3.2], [1.8, -1], [3.2, 2.6], [0, 1.4], [-3.2, 2.6], [-1.8, -1]], dens: 6.0, com: [0, 0.4] },
};

class Body {
  constructor(type, x, y) {
    const spec = SPECS[type];
    this.type = type;
    const art = spec.poly;
    const area = Math.abs(polyArea(art));
    const cen = polyCentroid(art);
    this.comX = cen[0] + spec.com[0]; this.comY = cen[1] + spec.com[1]; // COM in art space
    this.hull = art.map(p => [p[0] - this.comX, p[1] - this.comY]);
    this.area = area;
    this.mass = spec.mass || spec.dens * area;
    // samples for fluid coupling
    let minx = 1e9, maxx = -1e9, miny = 1e9, maxy = -1e9, rad = 0;
    for (const p of this.hull) {
      minx = Math.min(minx, p[0]); maxx = Math.max(maxx, p[0]);
      miny = Math.min(miny, p[1]); maxy = Math.max(maxy, p[1]);
      rad = Math.max(rad, Math.hypot(p[0], p[1]));
    }
    this.radius = rad;
    const step = 0.8, sm = [];
    for (let sy = miny + step / 2; sy < maxy; sy += step)
      for (let sx = minx + step / 2; sx < maxx; sx += step)
        if (inPoly(this.hull, sx, sy)) sm.push(sx, sy);
    this.samples = Float32Array.from(sm);
    this.ns = sm.length / 2;
    this.sv = area / this.ns;
    let I = 0;
    for (let i = 0; i < this.ns; i++) I += (this.mass / this.ns) * (sm[2 * i] ** 2 + sm[2 * i + 1] ** 2);
    this.inertia = I;
    this.cd = Math.min(0.1, 0.5 * this.mass / (4.4 * area)); // drag, kept stable for light bodies
    this.x = x; this.y = y; this.a = 0; this.vx = 0; this.vy = 0; this.w = 0;
    this.Fx = 0; this.Fy = 0; this.T = 0;
    this.grab = null;
  }

  toWorld(lx, ly, out) {
    const c = Math.cos(this.a), s = Math.sin(this.a);
    out[0] = this.x + c * lx - s * ly; out[1] = this.y + s * lx + c * ly;
  }
  contains(x, y) {
    const dx = x - this.x, dy = y - this.y;
    if (dx * dx + dy * dy > (this.radius + 0.8) ** 2) return false;
    const c = Math.cos(-this.a), s = Math.sin(-this.a);
    return inPoly(this.hull, c * dx - s * dy, s * dx + c * dy);
  }

  // Sample the fluid, accumulate forces, and push back on the fluid.
  couple(liq, grx, gry, dt, tmp, wp) {
    let Fx = 0, Fy = 0, T = 0;
    const c = Math.cos(this.a), s = Math.sin(this.a), sv = this.sv, cd = this.cd;
    for (let i = 0; i < this.ns; i++) {
      const lx = this.samples[2 * i], ly = this.samples[2 * i + 1];
      const rx = c * lx - s * ly, ry = s * lx + c * ly;
      liq.sample(this.x + rx, this.y + ry, tmp);
      const rho = tmp.rho;
      if (rho < 1e-4) continue;
      const svx = this.vx - this.w * ry, svy = this.vy + this.w * rx;
      const fx = -rho * sv * grx + cd * rho * sv * (tmp.u - svx);
      const fy = -rho * sv * gry + cd * rho * sv * (tmp.v - svy);
      Fx += fx; Fy += fy; T += rx * fy - ry * fx;
      liq.push(this.x + rx, this.y + ry, -fx, -fy, dt);
    }
    this.Fx = Fx; this.Fy = Fy; this.T = T;
  }

  integrate(grx, gry, dt, B) {
    const M = this.mass;
    let ax = this.Fx / M + grx, ay = this.Fy / M + gry, aw = this.T / this.inertia;
    if (this.grab) {
      const g = this.grab, c = Math.cos(this.a), s = Math.sin(this.a);
      const rx = c * g.lx - s * g.ly, ry = s * g.lx + c * g.ly;
      const px = this.x + rx, py = this.y + ry;
      const pvx = this.vx - this.w * ry, pvy = this.vy + this.w * rx;
      const fx = 0.04 * (g.tx - px) - 0.25 * pvx, fy = 0.04 * (g.ty - py) - 0.25 * pvy;
      ax += fx; ay += fy; aw += (rx * fy - ry * fx) * M / this.inertia;
    }
    this.vx += ax * dt; this.vy += ay * dt; this.w += aw * dt;
    this.w *= Math.pow(0.992, dt);
    let sp = Math.hypot(this.vx, this.vy);
    if (sp > 1.5) { this.vx *= 1.5 / sp; this.vy *= 1.5 / sp; }
    if (this.w > 0.3) this.w = 0.3; else if (this.w < -0.3) this.w = -0.3;
    this.x += this.vx * dt; this.y += this.vy * dt; this.a += this.w * dt;
    this.collideWalls(B);
  }

  collideWalls(B) {
    const M = this.mass, I = this.inertia;
    for (let wall = 0; wall < 4; wall++) {
      for (let k = 0; k < this.hull.length; k++) {
        const c = Math.cos(this.a), s = Math.sin(this.a);
        const lx = this.hull[k][0], ly = this.hull[k][1];
        const rx = c * lx - s * ly, ry = s * lx + c * ly;
        const wx = this.x + rx, wy = this.y + ry;
        let nx = 0, ny = 0, pen = 0;
        if (wall === 0) { pen = wy - B.ymax; ny = -1; }
        else if (wall === 1) { pen = B.ymin - wy; ny = 1; }
        else if (wall === 2) { pen = B.xmin - wx; nx = 1; }
        else { pen = wx - B.xmax; nx = -1; }
        if (pen <= 0) continue;
        this.x += nx * pen; this.y += ny * pen;
        const vpx = this.vx - this.w * ry, vpy = this.vy + this.w * rx;
        const vn = vpx * nx + vpy * ny;
        if (vn >= 0) continue;
        const rn = rx * ny - ry * nx, invM = 1 / M + rn * rn / I;
        const j = -(1.2) * vn / invM;
        this.vx += j * nx / M; this.vy += j * ny / M; this.w += j * rn / I;
        const tx = -ny, ty = nx, vt = vpx * tx + vpy * ty, rt = rx * ty - ry * tx;
        let jt = -vt / (1 / M + rt * rt / I);
        const lim = 0.4 * j;
        if (jt > lim) jt = lim; else if (jt < -lim) jt = -lim;
        this.vx += jt * tx / M; this.vy += jt * ty / M; this.w += jt * rt / I;
      }
    }
  }
}

function collideBodies(bodies) {
  for (let i = 0; i < bodies.length; i++) {
    for (let j = i + 1; j < bodies.length; j++) {
      const a = bodies[i], b = bodies[j];
      let dx = b.x - a.x, dy = b.y - a.y;
      const d = Math.hypot(dx, dy), rr = (a.radius + b.radius) * 0.9;
      if (d >= rr || d < 1e-6) continue;
      dx /= d; dy /= d;
      const pen = rr - d, tm = a.mass + b.mass;
      a.x -= dx * pen * b.mass / tm; a.y -= dy * pen * b.mass / tm;
      b.x += dx * pen * a.mass / tm; b.y += dy * pen * a.mass / tm;
      const vn = (b.vx - a.vx) * dx + (b.vy - a.vy) * dy;
      if (vn < 0) {
        const jn = -1.2 * vn / (1 / a.mass + 1 / b.mass);
        a.vx -= jn * dx / a.mass; a.vy -= jn * dy / a.mass;
        b.vx += jn * dx / b.mass; b.vy += jn * dy / b.mass;
      }
    }
  }
}

class World {
  constructor(GX, GY) {
    this.GX = GX; this.GY = GY;
    this.liq = new Liquid(GX, GY, LAYERS);
    this.bodies = [];
    this.g = 0.02;      // gravity strength
    this.tilt = 0;      // radians
    this.dt = 0.5;
    this.bounds = { xmin: 3.2, xmax: GX - 4.2, ymin: 3.2, ymax: GY - 4.2 };
    this.stirReq = null;
    this._tmp = { rho: 0, dens: 0, u: 0, v: 0 };
  }
  add(type, x, y) {
    const b = new Body(type, x, y);
    this.bodies.push(b);
    return b;
  }
  step() {
    const grx = this.g * Math.sin(this.tilt), gry = this.g * Math.cos(this.tilt), dt = this.dt;
    if (this.stirReq) {
      const r = this.stirReq;
      this.liq.stir(r.x, r.y, r.dx, r.dy, r.r);
      r.dx = 0; r.dy = 0;
    }
    this.liq.prepare(grx, gry, dt);
    for (const b of this.bodies) b.couple(this.liq, grx, gry, dt, this._tmp);
    this.liq.advect(dt);
    for (const b of this.bodies) b.integrate(grx, gry, dt, this.bounds);
    collideBodies(this.bodies);
  }
}

const api = { Liquid, Body, World, LAYERS, SPECS, REST };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
else root.LL = api;
})(typeof window !== 'undefined' ? window : globalThis);
