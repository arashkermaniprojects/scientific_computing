/* BCS 311 — Modules 4 & 5 | practice-gen45.js
   Endless random practice for CLO-4 (interpolation) and CLO-5 (integration, ODEs).
   Plugs into the widget in practice-gen.js (load AFTER it).
   Keys: interp1 · lagrange · newton   (Module 4)
         trap · simp · gaussq · euler · rk4   (Module 5)
   Every displayed intermediate value is the one carried forward. */
(function () {
  "use strict";
  const GEN = window.BCS311_GEN; if (!GEN) return;
  const rnd = (a, b) => a + Math.floor(Math.random() * (b - a + 1));
  const pick = a => a[rnd(0, a.length - 1)];
  const R = (x, d) => { const p = Math.pow(10, d); return Math.round(x * p) / p; };
  const F = (x, d) => { const s = R(x, d).toFixed(d); return /^-0\.0*$/.test(s) ? s.slice(1) : s; };
  const P = v => (String(v)[0] === "-" ? "(" + v + ")" : String(v));
  const sg = (x, d) => (x < 0 ? "-" : "+") + F(Math.abs(x), d);

  function distinctX(n, lo, hi) { const s = new Set(); while (s.size < n) s.add(rnd(lo, hi)); return [...s].sort((a, b) => a - b); }
  const ptsTeX = (xs, ys) => xs.map((x, i) => "(" + x + ",\\," + ys[i] + ")").join(",\\ ");

  /* ---------------- Module 4 ---------------- */
  function genInterp1() {
    const x0 = rnd(0, 4), dx = rnd(1, 3), x1 = x0 + dx, y0 = rnd(-3, 9), y1 = y0 + rnd(-6, 8) + (Math.random() < .5 ? 1 : 0);
    const xq = x0 + pick([0.25, 0.5, 0.75]) * dx;
    const sl = (y1 - y0) / (x1 - x0), v = y0 + sl * (xq - x0);
    return {
      given: `Given the two points \\((${x0},${y0})\\) and \\((${x1},${y1})\\), use linear interpolation to estimate \\(y\\) at \\(x=${xq}\\).`,
      steps: [
        { title: "Formula", html: `\\(y=y_0+\\dfrac{y_1-y_0}{x_1-x_0}(x-x_0)\\) with \\(x_0=${x0},\\,y_0=${y0},\\,x_1=${x1},\\,y_1=${y1}\\).` },
        { title: "Slope", html: `\\(\\dfrac{y_1-y_0}{x_1-x_0}=\\dfrac{${y1}-${P(y0)}}{${x1}-${x0}}=${F(sl, 4)}\\)` },
        { title: "Final answer (a number)", html: `\\(y(${xq})=${y0}+${P(F(sl, 4))}(${xq}-${x0})=\\boxed{${F(v, 4)}}\\)` }
      ]
    };
  }
  genInterp1.title = "linear interpolation";

  function genLagrange() {
    const xs = distinctX(3, 0, 5), ys = [rnd(-3, 8), rnd(-3, 8), rnd(-3, 8)];
    let xq; do { xq = rnd(0, 5); } while (xs.includes(xq));
    const L = [], steps = [];
    for (let i = 0; i < 3; i++) {
      const o = [0, 1, 2].filter(j => j !== i);
      const num = (xq - xs[o[0]]) * (xq - xs[o[1]]), den = (xs[i] - xs[o[0]]) * (xs[i] - xs[o[1]]);
      L.push(num / den);
      steps.push({ title: `Basis \\(L_{${i}}(${xq})\\)`, html:
        `\\(L_{${i}}(x)=\\dfrac{(x-${P(xs[o[0]])})(x-${P(xs[o[1]])})}{(${xs[i]}-${P(xs[o[0]])})(${xs[i]}-${P(xs[o[1]])})}\\Rightarrow L_{${i}}(${xq})=\\dfrac{(${xq - xs[o[0]]})(${xq - xs[o[1]]})}{(${xs[i] - xs[o[0]]})(${xs[i] - xs[o[1]]})}=\\dfrac{${num}}{${den}}=${F(L[i], 4)}\\)` });
    }
    const Lr = L.map(v => R(v, 4)), v = ys.reduce((s, y, i) => s + y * Lr[i], 0);
    steps.push({ title: "Final answer (a number)", html:
      `\\(P(${xq})=\\sum y_iL_i=${ys.map((y, i) => `${i ? "+" : ""}${P(y)}(${F(Lr[i], 4)})`).join("")}=\\boxed{${F(v, 4)}}\\)<br>(Check: \\(L_0+L_1+L_2=${F(Lr[0] + Lr[1] + Lr[2], 4)}\\approx 1\\).)` });
    return { given: `Data: \\(${ptsTeX(xs, ys)}\\). Use the Lagrange form to estimate \\(P(${xq})\\). (Each \\(L_i\\) uses only the \\(x\\)-values.)`, steps };
  }
  genLagrange.title = "Lagrange polynomial";

  function genNewton() {
    const xs = distinctX(3, 0, 5), ys = [rnd(-3, 8), rnd(-3, 8), rnd(-3, 8)];
    let xq; do { xq = rnd(0, 6); } while (xs.includes(xq));
    const d01 = R((ys[1] - ys[0]) / (xs[1] - xs[0]), 4), d12 = R((ys[2] - ys[1]) / (xs[2] - xs[1]), 4);
    const d012 = R((d12 - d01) / (xs[2] - xs[0]), 4);
    const v = ys[0] + d01 * (xq - xs[0]) + d012 * (xq - xs[0]) * (xq - xs[1]);
    return {
      given: `Data: \\(${ptsTeX(xs, ys)}\\). Build the Newton divided-difference polynomial and estimate \\(P(${xq})\\).`,
      steps: [
        { title: "First divided differences", html:
          `\\(f[x_0,x_1]=\\dfrac{${ys[1]}-${P(ys[0])}}{${xs[1]}-${xs[0]}}=${F(d01, 4)}\\)<br>\\(f[x_1,x_2]=\\dfrac{${ys[2]}-${P(ys[1])}}{${xs[2]}-${xs[1]}}=${F(d12, 4)}\\)` },
        { title: "Second divided difference", html:
          `\\(f[x_0,x_1,x_2]=\\dfrac{${F(d12, 4)}-${P(F(d01, 4))}}{${xs[2]}-${xs[0]}}=${F(d012, 4)}\\)` },
        { title: "The polynomial (an equation)", html:
          `\\(P(x)=${ys[0]}${sg(d01, 4)}(x-${P(xs[0])})${sg(d012, 4)}(x-${P(xs[0])})(x-${P(xs[1])})\\)` },
        { title: "Final answer (a number)", html:
          `\\(P(${xq})=${ys[0]}${sg(d01, 4)}(${xq - xs[0]})${sg(d012, 4)}(${xq - xs[0]})(${xq - xs[1]})=\\boxed{${F(v, 4)}}\\)` }
      ]
    };
  }
  genNewton.title = "Newton divided differences";

  /* ---------------- Module 5 ---------------- */
  function quadFam() {
    const c2 = rnd(1, 4), c1 = pick([0, 0, 1, 2]), c0 = rnd(1, 5), a = rnd(0, 2), b = a + 4;
    const f = x => c2 * x * x + c1 * x + c0;
    const tex = `${c2 === 1 ? "" : c2}x^2${c1 ? "+" + (c1 === 1 ? "" : c1) + "x" : ""}+${c0}`;
    const exact = c2 * (b ** 3 - a ** 3) / 3 + c1 * (b * b - a * a) / 2 + c0 * (b - a);
    return { f, tex, a, b, exact };
  }

  function genTrap() {
    const q = quadFam(), n = 4, h = (q.b - q.a) / n, xs = [0, 1, 2, 3, 4].map(i => q.a + i * h), ys = xs.map(q.f);
    const s = ys[0] + 2 * (ys[1] + ys[2] + ys[3]) + ys[4], T = h / 2 * s;
    return {
      given: `Estimate \\(\\displaystyle\\int_{${q.a}}^{${q.b}}(${q.tex})\\,dx\\) with the trapezoidal rule, \\(n=4\\). Compare with the exact value.`,
      header: ["i", "\\(x_i\\)", "\\(f(x_i)\\)"],
      steps: xs.map((x, i) => ({ title: `Node \\(x_${i}=${x}\\)`, html: `\\(f(${x})=${ys[i]}\\)`, row: [i, x, ys[i]] })).concat([
        { title: "Final answer (a number)", html: `\\(h=${h}\\).<br>\\(T=\\dfrac h2\\,[f_0+2(f_1+f_2+f_3)+f_4]=\\dfrac{${h}}{2}[${ys[0]}+2(${ys[1]}+${ys[2]}+${ys[3]})+${ys[4]}]=\\dfrac{${h}}{2}(${s})=\\boxed{${F(T, 4)}}\\)<br>Exact \\(=${F(q.exact, 4)}\\); error \\(=${F(T - q.exact, 4)}\\) (the trapezoid overestimates because \\(f\\) curves upward).` }])
    };
  }
  genTrap.title = "trapezoidal rule";

  function genSimp() {
    const q = quadFam(), n = 4, h = (q.b - q.a) / n, xs = [0, 1, 2, 3, 4].map(i => q.a + i * h), ys = xs.map(q.f);
    const s = ys[0] + 4 * (ys[1] + ys[3]) + 2 * ys[2] + ys[4], S = h / 3 * s;
    return {
      given: `Estimate \\(\\displaystyle\\int_{${q.a}}^{${q.b}}(${q.tex})\\,dx\\) with Simpson's 1/3 rule, \\(n=4\\). Compare with the exact value.`,
      header: ["i", "\\(x_i\\)", "\\(f(x_i)\\)", "weight"],
      steps: xs.map((x, i) => ({ title: `Node \\(x_${i}=${x}\\)`, html: `\\(f(${x})=${ys[i]}\\), weight ${[1, 4, 2, 4, 1][i]}`, row: [i, x, ys[i], [1, 4, 2, 4, 1][i]] })).concat([
        { title: "Final answer (a number)", html: `\\(h=${h}\\).<br>\\(S=\\dfrac h3[f_0+4f_1+2f_2+4f_3+f_4]=\\dfrac{${h}}{3}[${ys[0]}+4(${ys[1]})+2(${ys[2]})+4(${ys[3]})+${ys[4]}]=\\dfrac{${h}}{3}(${s})=\\boxed{${F(S, 4)}}\\)<br>Exact \\(=${F(q.exact, 4)}\\). Simpson is exact for polynomials up to degree 3, so the two agree.` }])
    };
  }
  genSimp.title = "Simpson's 1/3 rule";

  function genGauss() {
    const q = quadFam(), m = (q.b - q.a) / 2, c = (q.b + q.a) / 2, t = 1 / Math.sqrt(3);
    const x1 = R(c - m * t, 4), x2 = R(c + m * t, 4), f1 = R(q.f(x1), 4), f2 = R(q.f(x2), 4), G = R(m * (f1 + f2), 4);
    return {
      given: `Use 2-point Gaussian quadrature to estimate \\(\\displaystyle\\int_{${q.a}}^{${q.b}}(${q.tex})\\,dx\\). (Nodes \\(t=\\pm1/\\sqrt3\\), weights 1 on \\([-1,1]\\).)`,
      steps: [
        { title: "Map [a,b] to [−1,1]", html: `\\(x=\\dfrac{b-a}{2}t+\\dfrac{b+a}{2}=${m}t+${c}\\), and \\(dx=${m}\\,dt\\).` },
        { title: "Nodes in x", html: `\\(x_1=${m}(-0.5774)+${c}=${F(x1, 4)}\\), \\(x_2=${m}(0.5774)+${c}=${F(x2, 4)}\\)` },
        { title: "Function values", html: `\\(f(x_1)=${F(f1, 4)}\\), \\(f(x_2)=${F(f2, 4)}\\)` },
        { title: "Final answer (a number)", html: `\\(G=\\dfrac{b-a}{2}[w_1f(x_1)+w_2f(x_2)]=${m}(${F(f1, 4)}+${F(f2, 4)})=\\boxed{${F(G, 4)}}\\)<br>Exact \\(=${F(q.exact, 4)}\\): two nodes are exact for polynomials up to degree 3 (tiny differences are rounding).` }
      ]
    };
  }
  genGauss.title = "Gaussian quadrature";

  const ODEF = [
    { tex: "x+y", f: (x, y) => x + y, fs: (x, y) => `${x}+${P(y)}` },
    { tex: "y-x", f: (x, y) => y - x, fs: (x, y) => `${P(y)}-${P(x)}` },
    { tex: "2y", f: (x, y) => 2 * y, fs: (x, y) => `2(${y})` },
    { tex: "x-y", f: (x, y) => x - y, fs: (x, y) => `${x}-${P(y)}` }
  ];

  function genEuler() {
    const o = pick(ODEF), h = pick([0.1, 0.2, 0.25, 0.5]), y0 = rnd(1, 3), n = 3;
    let x = 0, y = y0; const steps = [];
    for (let i = 0; i < n; i++) {
      const fv = R(o.f(x, y), 4), yn = R(y + h * fv, 4);
      steps.push({ title: `Step ${i + 1}: \\(x_${i}=${F(x, 2)}\\to x_${i + 1}=${F(x + h, 2)}\\)`, html:
        `\\(f(x_${i},y_${i})=${o.fs(F(x, 2), F(y, 4))}=${F(fv, 4)}\\)<br>\\(y_${i + 1}=y_${i}+h f=${F(y, 4)}+${h}(${F(fv, 4)})=${F(yn, 4)}\\)`,
        row: [i, F(x, 2), F(y, 4), F(fv, 4), F(yn, 4)] });
      x = R(x + h, 4); y = yn;
    }
    steps[n - 1].title = "Final answer (a number): " + steps[n - 1].title;
    return { given: `\\(\\dfrac{dy}{dx}=${o.tex},\\ y(0)=${y0}\\). Use Euler's method with \\(h=${h}\\) to estimate \\(y(${F(n * h, 2)})\\) in ${n} steps.`,
      header: ["n", "\\(x_n\\)", "\\(y_n\\)", "\\(f(x_n,y_n)\\)", "\\(y_{n+1}\\)"], steps };
  }
  genEuler.title = "Euler's method";

  function genRK4() {
    const o = pick(ODEF.slice(0, 2).concat(ODEF.slice(3))), h = pick([0.2, 0.5]), y0 = rnd(1, 3), x0 = 0;
    const k1 = R(o.f(x0, y0), 4), k2 = R(o.f(x0 + h / 2, y0 + h / 2 * k1), 4), k3 = R(o.f(x0 + h / 2, y0 + h / 2 * k2), 4), k4 = R(o.f(x0 + h, y0 + h * k3), 4);
    const y1 = R(y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 4);
    return {
      given: `\\(\\dfrac{dy}{dx}=${o.tex},\\ y(0)=${y0}\\). Take <b>one</b> step of classical RK4 with \\(h=${h}\\) to estimate \\(y(${h})\\).`,
      steps: [
        { title: "Slope at the start", html: `\\(k_1=f(0,${y0})=${F(k1, 4)}\\)` },
        { title: "Slope at the midpoint (using k₁)", html: `\\(k_2=f(${h / 2},\\,${y0}+${h / 2}(${F(k1, 4)}))=f(${h / 2},${F(y0 + h / 2 * k1, 4)})=${F(k2, 4)}\\)` },
        { title: "Midpoint again (using k₂)", html: `\\(k_3=f(${h / 2},\\,${y0}+${h / 2}(${F(k2, 4)}))=f(${h / 2},${F(y0 + h / 2 * k2, 4)})=${F(k3, 4)}\\)` },
        { title: "Slope at the end (using k₃)", html: `\\(k_4=f(${h},\\,${y0}+${h}(${F(k3, 4)}))=f(${h},${F(y0 + h * k3, 4)})=${F(k4, 4)}\\)` },
        { title: "Final answer (a number)", html: `\\(y_1=y_0+\\dfrac h6(k_1+2k_2+2k_3+k_4)=${y0}+\\dfrac{${h}}{6}(${F(k1, 4)}+2(${F(k2, 4)})+2(${F(k3, 4)})+${F(k4, 4)})=\\boxed{${F(y1, 4)}}\\)` }
      ]
    };
  }
  genRK4.title = "Runge–Kutta 4";

  Object.assign(GEN, { interp1: genInterp1, lagrange: genLagrange, newton4: genNewton, trap: genTrap, simp: genSimp, gaussq: genGauss, euler: genEuler, rk4: genRK4 });
})();
