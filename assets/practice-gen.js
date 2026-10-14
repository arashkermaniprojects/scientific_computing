/* ===========================================================
   BCS 311 — Module 2 | practice-gen.js
   Random, step-by-step practice problems.
   Usage:  <div class="gen-practice" data-gen="newton"></div>
   Keys:   fixed · bisect · newton · secant · gauss · lu · chol
   Every displayed intermediate value is the one carried forward,
   so a student rounding to the stated decimals reproduces it exactly.
   =========================================================== */
(function () {
  "use strict";

  /* ---------- small helpers ---------- */
  const rnd  = (a, b) => a + Math.floor(Math.random() * (b - a + 1));
  const pick = arr => arr[rnd(0, arr.length - 1)];
  const R = (x, d) => { const p = Math.pow(10, d); return Math.round(x * p) / p; };
  const F = (x, d) => R(x, d).toFixed(d);                               // "1.4142", never "-0.0000"
  const S = (x, d) => { const s = F(x, d); return s[0] === "-" ? s : "+" + s; }; // signed
  const P = v => (String(v)[0] === "-" ? "(" + v + ")" : String(v));   // parenthesise a negative factor
  const EPS = t => (t >= 0.01 ? String(t) : "10^{" + Math.round(Math.log10(t)) + "}");

  /* TeX builders */
  const bm  = M => "\\begin{bmatrix}" + M.map(r => r.join("&")).join("\\\\") + "\\end{bmatrix}";
  const col = v => "\\begin{bmatrix}" + v.join("\\\\") + "\\end{bmatrix}";
  const aug = (M, b) => "\\left[\\begin{array}{ccc|c}" + M.map((r, i) => r.concat([b[i]]).join("&")).join("\\\\") + "\\end{array}\\right]";
  const mul = (A, B) => A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
  const mv  = (A, x) => A.map(r => r.reduce((s, v, k) => s + v * x[k], 0));
  const tr  = A => A[0].map((_, j) => A.map(r => r[j]));
  function lin(cs, vs) {                      // [2,3,-1] → "2x+3y-z"
    let s = "";
    cs.forEach((c, i) => {
      if (c === 0) return;
      const t = (Math.abs(c) === 1 ? "" : Math.abs(c)) + vs[i];
      s += s === "" ? (c < 0 ? "-" : "") + t : (c < 0 ? "-" : "+") + t;
    });
    return s || "0";
  }
  function subst(cs, vals) {                  // "2(1)+3(-2)-(3)" with numbers in place of unknowns
    let s = "";
    cs.forEach((c, i) => {
      const t = (Math.abs(c) === 1 ? "" : Math.abs(c)) + "(" + vals[i] + ")";
      s += s === "" ? (c < 0 ? "-" : "") + t : (c < 0 ? "-" : "+") + t;
    });
    return s;
  }
  const TRUE = Fm => Fm.truth ? `<br>(True value \\(${Fm.truth}=${F(Fm.tv, 5)}\\).)` : "";
  const CONT = (E, tol, n) => E < tol
    ? `\\(E_{${n}}=${F(E, 4)}\\lt\\varepsilon=${EPS(tol)}\\) ⇒ <b>stop</b> ✓`
    : `\\(E_{${n}}=${F(E, 4)}\\ge\\varepsilon=${EPS(tol)}\\) ⇒ continue`;

  /* ---------- random f(x) families for bisection / Newton / secant ---------- */
  function rootFamily() {
    const k = rnd(0, 4);
    if (k === 0) {
      let c; do { c = rnd(2, 40); } while (Number.isInteger(Math.sqrt(c)));
      const a = Math.floor(Math.sqrt(c));
      return { fTex: `x^{2}-${c}`, dTex: "2x", f: x => x * x - c, d: x => 2 * x, a, b: a + 1,
               fSub: x => `(${x})^{2}-${c}`, dSub: x => `2(${x})`, truth: `\\sqrt{${c}}`, tv: Math.sqrt(c) };
    }
    if (k === 1) {
      let c; do { c = rnd(2, 80); } while (Math.round(Math.cbrt(c)) ** 3 === c);
      const a = Math.floor(Math.cbrt(c));
      return { fTex: `x^{3}-${c}`, dTex: "3x^{2}", f: x => x ** 3 - c, d: x => 3 * x * x, a, b: a + 1,
               fSub: x => `(${x})^{3}-${c}`, dSub: x => `3(${x})^{2}`, truth: `\\sqrt[3]{${c}}`, tv: Math.cbrt(c) };
    }
    if (k === 2) {
      const c = rnd(3, 40); const f = x => x ** 3 + x - c;
      let a = 0; while (f(a + 1) <= 0) a++;
      return { fTex: `x^{3}+x-${c}`, dTex: "3x^{2}+1", f, d: x => 3 * x * x + 1, a, b: a + 1,
               fSub: x => `(${x})^{3}+(${x})-${c}`, dSub: x => `3(${x})^{2}+1` };
    }
    if (k === 3) {
      const c = rnd(2, 20); const a = Math.floor(Math.log(c));
      return { fTex: `e^{x}-${c}`, dTex: "e^{x}", f: x => Math.exp(x) - c, d: x => Math.exp(x), a, b: a + 1,
               fSub: x => `e^{${x}}-${c}`, dSub: x => `e^{${x}}`, truth: `\\ln ${c}`, tv: Math.log(c) };
    }
    const m = pick([1, 2]); const f = x => m * Math.cos(x) - x;
    let a = 0; while (f(a + 1) >= 0) a++;
    const ms = m === 1 ? "" : m;
    return { fTex: `${ms}\\cos x-x`, dTex: `-${ms}\\sin x-1`, f, d: x => -m * Math.sin(x) - 1, a, b: a + 1, rad: true,
             fSub: x => `${ms}\\cos(${x})-(${x})`, dSub: x => `-${ms}\\sin(${x})-1` };
  }

  /* ========== 2.1 FIXED-POINT ========== */
  function genFixed() {
    const fam = rnd(0, 2); const d = 4, tol = pick([0.001, 0.0005]), nmax = 12;
    let a, b, g, dg, fTex, gTex, dgTex, gSub, dgSub, f;
    if (fam === 0) {                                  // x^3 - a x - b = 0  →  x = cbrt(a x + b)
      a = rnd(1, 3); b = rnd(2, 9); const as = a === 1 ? "" : a;
      f = x => x ** 3 - a * x - b; g = x => Math.cbrt(a * x + b); dg = x => a / (3 * Math.cbrt(a * x + b) ** 2);
      fTex = `x^{3}-${as}x-${b}`; gTex = `\\sqrt[3]{${as}x+${b}}`; dgTex = `\\dfrac{${a}}{3\\,(${as}x+${b})^{2/3}}`;
      gSub = x => `\\sqrt[3]{${a}(${x})+${b}}=\\sqrt[3]{${F(a * (+x) + b, d)}}`;
      dgSub = x => `\\dfrac{${a}}{3\\,(${F(a * (+x) + b, d)})^{2/3}}`;
    } else if (fam === 1) {                           // x^2 - a x - b = 0  →  x = sqrt(a x + b)
      a = rnd(1, 4); b = rnd(1, 9); const as = a === 1 ? "" : a;
      f = x => x * x - a * x - b; g = x => Math.sqrt(a * x + b); dg = x => a / (2 * Math.sqrt(a * x + b));
      fTex = `x^{2}-${as}x-${b}`; gTex = `\\sqrt{${as}x+${b}}`; dgTex = `\\dfrac{${a}}{2\\sqrt{${as}x+${b}}}`;
      gSub = x => `\\sqrt{${a}(${x})+${b}}=\\sqrt{${F(a * (+x) + b, d)}}`;
      dgSub = x => `\\dfrac{${a}}{2\\sqrt{${F(a * (+x) + b, d)}}}`;
    } else {                                          // x^2 + a x - b = 0  →  x = b/(x + a)
      a = rnd(3, 5); b = rnd(1, 6);
      f = x => x * x + a * x - b; g = x => b / (x + a); dg = x => -b / ((x + a) ** 2);
      fTex = `x^{2}+${a}x-${b}`; gTex = `\\dfrac{${b}}{x+${a}}`; dgTex = `-\\dfrac{${b}}{(x+${a})^{2}}`;
      gSub = x => `\\dfrac{${b}}{${x}+${a}}=\\dfrac{${b}}{${F(+x + a, d)}}`;
      dgSub = x => `-\\dfrac{${b}}{(${F(+x + a, d)})^{2}}`;
    }
    let xs = 1; for (let i = 0; i < 300; i++) xs = g(xs);             // the fixed point
    let x0 = pick([Math.floor(xs), Math.ceil(xs), Math.floor(xs) + 0.5]);
    if (Math.abs(g(x0) - x0) < 0.02) x0 = Math.floor(xs) + 1;
    const x0s = Number.isInteger(x0) ? String(x0) : F(x0, 1);

    const given = `Solve \\(${fTex}=0\\) by <b>fixed-point iteration</b> with the rearrangement \\(x=g(x)=${gTex}\\).<br>
      initial value x₀ = ${x0s},   tolerance ε = ${tol.toFixed(4).replace(/0+$/,'')},   max ${nmax} iterations,   work to ${d} decimals (round \\(x_{n+1}\\) and carry the rounded value).<br>
      Stop when \\(E_n=|x_{n+1}-x_n|\\lt\\varepsilon\\).`;
    const steps = [];
    const g0 = dg(x0);
    steps.push({ title: "Set up and check the convergence condition", html:
      `\\(${fTex}=0 \\;\\Leftrightarrow\\; x=${gTex}\\), so \\(g(x)=${gTex}\\) and \\(g'(x)=${dgTex}\\).<br>
       At the start: \\(|g'(x_0)|=\\left|${dgSub(x0s)}\\right|=${F(Math.abs(g0), d)}\\lt 1\\) ⇒ the iteration is expected to converge near \\(x_0\\). ✓` });
    let x = x0, n = 0, x1 = x0, E = 0, ok = false;
    for (; n < nmax; n++) {
      x1 = R(g(x), d); E = R(Math.abs(x1 - x), d);
      const xsn = Number.isInteger(x) && n === 0 ? x0s : F(x, d);
      steps.push({ title: `Iteration n = ${n}`, html:
        `\\(x_{${n + 1}}=g(x_{${n}})=${gSub(xsn)}=${F(x1, d)}\\)<br>
         \\(E_{${n}}=|x_{${n + 1}}-x_{${n}}|=|${F(x1, d)}-${xsn}|=${F(E, d)}\\)<br>${CONT(E, tol, n)}`,
        row: [n, xsn, F(x1, d), F(E, d), E < tol ? "✓ stop" : "no"], rowClass: E < tol ? "converged" : "" });
      if (E < tol) { ok = true; break; }
      x = x1;
    }
    steps.push({ title: ok ? "Answer" : "Maximum iterations reached", html:
      (ok ? `Stopped at \\(n=${n}\\) because \\(E_{${n}}\\lt\\varepsilon\\).` : `\\(E_n\\) never dropped below \\(\\varepsilon\\) within ${nmax} iterations — report the latest value.`) +
      ` <b>Root \\(x\\approx${F(x1, d)}\\).</b><br>Check in the original equation: \\(f(${F(x1, d)})=${F(f(x1), d)}\\approx0\\). ✓` });
    return { given, header: ["n", "\\(x_n\\)", "\\(x_{n+1}=g(x_n)\\)", "\\(E_n=|x_{n+1}-x_n|\\)", "\\(E_n\\lt\\varepsilon?\\)"], steps };
  }
  genFixed.title = "fixed-point iteration";

  /* ========== 2.2 BISECTION ========== */
  function genBisect() {
    const Fm = rootFamily(); const d = 4, tol = pick([0.05, 0.02, 0.01]), nmax = 12;
    let a = Fm.a, b = Fm.b;
    const given = `Solve \\(f(x)=${Fm.fTex}=0\\) by <b>bisection</b>.<br>
      initial interval [a, b] = [${a}, ${b}],   c = (a+b)/2,   tolerance ε = ${tol},   max ${nmax} iterations,   work to ${d} decimals${Fm.rad ? " (x in radians)" : ""}.<br>
      Stop when \\(E_n=\\tfrac{b-a}{2}\\lt\\varepsilon\\). Keep the half in which \\(f\\) changes sign.`;
    const steps = [];
    let fa = R(Fm.f(a), d), fb = R(Fm.f(b), d);
    steps.push({ title: "Check the starting bracket", html:
      `\\(f(a)=f(${a})=${Fm.fSub(a)}=${S(fa, d)}\\),  \\(f(b)=f(${b})=${Fm.fSub(b)}=${S(fb, d)}\\).<br>
       \\(f(a)\\,f(b)\\lt 0\\) ⇒ opposite signs ⇒ a root lies in \\([${a},${b}]\\). ✓` });
    let n = 1, c = a, fc = 0, E = 1, ok = false;
    for (; n <= nmax; n++) {
      c = R((a + b) / 2, d); fc = R(Fm.f(c), d); E = R((b - a) / 2, d);
      const left = fa * fc < 0;                       // sign change in [a, c]?
      const na = left ? a : c, nb = left ? c : b;
      steps.push({ title: `Iteration n = ${n}`, html:
        `\\(c=\\tfrac{a+b}{2}=\\tfrac{${F(a, d)}+${F(b, d)}}{2}=${F(c, d)}\\)<br>
         \\(f(c)=${Fm.fSub(F(c, d))}=${S(fc, d)}\\)<br>
         \\(f(a)\\,f(c)=(${S(fa, d)})(${S(fc, d)})\\) is ${left ? "<b>negative</b> ⇒ the root is in \\([a,c]\\) ⇒ set \\(b\\leftarrow c\\)" : "<b>positive</b> ⇒ the root is in \\([c,b]\\) ⇒ set \\(a\\leftarrow c\\)"}<br>
         New interval \\([${F(na, d)},\\,${F(nb, d)}]\\);  error bound \\(E_{${n}}=\\tfrac{b-a}{2}=\\tfrac{${F(b, d)}-${F(a, d)}}{2}=${F(E, d)}\\).<br>${CONT(E, tol, n)}`,
        row: [n, F(a, d), F(b, d), F(c, d), S(fc, d), left ? "[a, c]" : "[c, b]", F(E, d), E < tol ? "✓ stop" : "no"], rowClass: E < tol ? "converged" : "" });
      if (left) { b = c; fb = fc; } else { a = c; fa = fc; }
      if (E < tol) { ok = true; break; }
    }
    steps.push({ title: ok ? "Answer" : "Maximum iterations reached", html:
      (ok ? `Stopped at \\(n=${n}\\) because \\(E_{${n}}=${F(E, d)}\\lt\\varepsilon\\).` : `Stopped after ${nmax} iterations.`) +
      ` <b>Root \\(x\\approx c=${F(c, d)}\\)</b>, guaranteed to lie in \\([${F(a, d)},\\,${F(b, d)}]\\) (error ≤ ${F(E, d)}).` + TRUE(Fm) });
    return { given, header: ["n", "a", "b", "\\(c=\\tfrac{a+b}{2}\\)", "f(c)", "keep", "\\(E_n=\\tfrac{b-a}{2}\\)", "\\(E_n\\lt\\varepsilon?\\)"], steps };
  }
  genBisect.title = "bisection";

  /* ========== 2.3 NEWTON ========== */
  function genNewton() {
    const Fm = rootFamily(); const d = 4, tol = pick([0.001, 0.0005]), nmax = 8;
    const x0 = pick([Fm.a, Fm.b]);
    const given = `Solve \\(f(x)=${Fm.fTex}=0\\) by <b>Newton's method</b>, \\(x_{n+1}=x_n-\\dfrac{f(x_n)}{f'(x_n)}\\).<br>
      initial value x₀ = ${x0},   tolerance ε = ${tol.toFixed(4).replace(/0+$/, '')},   max ${nmax} iterations,   work to ${d} decimals (round \\(f\\), \\(f'\\) and \\(x_{n+1}\\); carry the rounded values)${Fm.rad ? ",   x in radians" : ""}.<br>
      Stop when \\(E_n=|x_{n+1}-x_n|\\lt\\varepsilon\\).`;
    const steps = [];
    steps.push({ title: "Differentiate and write the iteration", html:
      `\\(f(x)=${Fm.fTex}\\) ⇒ \\(f'(x)=${Fm.dTex}\\).<br>
       Newton step: \\(x_{n+1}=x_n-\\dfrac{${Fm.fTex}}{${Fm.dTex}}\\) evaluated at \\(x=x_n\\).` });
    let x = x0, n = 0, x1 = x0, E = 0, ok = false, fx = 0, dfx = 0;
    for (; n < nmax; n++) {
      fx = R(Fm.f(x), d); dfx = R(Fm.d(x), d);
      if (dfx === 0) { steps.push({ title: "Breakdown", html: `\\(f'(x_{${n}})=0\\): the tangent is horizontal and Newton cannot continue. Choose a different \\(x_0\\).` }); break; }
      x1 = R(x - fx / dfx, d); E = R(Math.abs(x1 - x), d);
      const xs = n === 0 ? String(x0) : F(x, d);
      steps.push({ title: `Iteration n = ${n}`, html:
        `\\(f(x_{${n}})=${Fm.fSub(xs)}=${F(fx, d)}\\)<br>
         \\(f'(x_{${n}})=${Fm.dSub(xs)}=${F(dfx, d)}\\)<br>
         \\(x_{${n + 1}}=x_{${n}}-\\dfrac{f(x_{${n}})}{f'(x_{${n}})}=${xs}-\\dfrac{${F(fx, d)}}{${F(dfx, d)}}=${F(x1, d)}\\)<br>
         \\(E_{${n}}=|x_{${n + 1}}-x_{${n}}|=|${F(x1, d)}-${xs}|=${F(E, d)}\\)<br>${CONT(E, tol, n)}`,
        row: [n, xs, F(fx, d), F(dfx, d), F(x1, d), F(E, d), E < tol ? "✓ stop" : "no"], rowClass: E < tol ? "converged" : "" });
      if (E < tol) { ok = true; break; }
      x = x1;
    }
    steps.push({ title: ok ? "Answer" : "Maximum iterations reached", html:
      (ok ? `Stopped at \\(n=${n}\\) because \\(E_{${n}}\\lt\\varepsilon\\).` : `Stopped after ${nmax} iterations.`) +
      ` <b>Root \\(x\\approx${F(x1, d)}\\).</b><br>Check: \\(f(${F(x1, d)})=${F(Fm.f(x1), d)}\\approx0\\). ✓` +
      ` Notice how the error shrinks faster and faster (quadratic convergence).` + TRUE(Fm) });
    return { given, header: ["n", "\\(x_n\\)", "\\(f(x_n)\\)", "\\(f'(x_n)\\)", "\\(x_{n+1}\\)", "\\(E_n\\)", "\\(E_n\\lt\\varepsilon?\\)"], steps };
  }
  genNewton.title = "Newton's method";

  /* ========== 2.4 SECANT ========== */
  function genSecant() {
    const Fm = rootFamily(); const d = 4, tol = 0.001, nmax = 8;
    let x0 = Fm.a, x1 = Fm.b;
    const given = `Solve \\(f(x)=${Fm.fTex}=0\\) by the <b>secant method</b>, \\(x_{k+1}=x_k-f(x_k)\\,\\dfrac{x_k-x_{k-1}}{f(x_k)-f(x_{k-1})}\\).<br>
      initial values x₀ = ${x0}, x₁ = ${x1},   tolerance ε = 0.001,   max ${nmax} iterations,   work to ${d} decimals (round \\(f\\) and \\(x_{k+1}\\); carry the rounded values)${Fm.rad ? ",   x in radians" : ""}.<br>
      Stop when \\(E_{k+1}=|x_{k+1}-x_k|\\lt\\varepsilon\\).`;
    const steps = [];
    let f0 = R(Fm.f(x0), d), f1 = R(Fm.f(x1), d);
    steps.push({ title: "Evaluate f at the two starting points", html:
      `\\(f(x_0)=${Fm.fSub(x0)}=${F(f0, d)}\\),  \\(f(x_1)=${Fm.fSub(x1)}=${F(f1, d)}\\).<br>Each step uses the <b>two most recent</b> points.` });
    let k = 1, x2 = x1, f2 = f1, E = 0, ok = false;
    for (; k <= nmax; k++) {
      if (f1 - f0 === 0) { steps.push({ title: "Breakdown", html: `\\(f(x_k)=f(x_{k-1})\\): the secant is horizontal — the method cannot continue.` }); break; }
      x2 = R(x1 - f1 * (x1 - x0) / (f1 - f0), d); f2 = R(Fm.f(x2), d); E = R(Math.abs(x2 - x1), d);
      const s0 = k === 1 ? String(x0) : F(x0, d), s1 = k === 1 ? String(x1) : F(x1, d);
      steps.push({ title: `Iteration k = ${k}  (from \\(x_{${k - 1}},x_{${k}}\\) to \\(x_{${k + 1}}\\))`, html:
        `\\(x_{${k + 1}}=x_{${k}}-f(x_{${k}})\\,\\dfrac{x_{${k}}-x_{${k - 1}}}{f(x_{${k}})-f(x_{${k - 1}})}\\)<br>
         \\(\\phantom{x_{${k + 1}}}=${s1}-${P(F(f1, d))}\\cdot\\dfrac{${s1}-${P(s0)}}{${F(f1, d)}-${P(F(f0, d))}}=${F(x2, d)}\\)<br>
         \\(f(x_{${k + 1}})=${Fm.fSub(F(x2, d))}=${F(f2, d)}\\)<br>
         \\(E_{${k + 1}}=|x_{${k + 1}}-x_{${k}}|=|${F(x2, d)}-${s1}|=${F(E, d)}\\)<br>${CONT(E, tol, k + 1)}`,
        row: [k, s0, s1, F(f0, d), F(f1, d), F(x2, d), F(f2, d), F(E, d), E < tol ? "✓ stop" : "no"], rowClass: E < tol ? "converged" : "" });
      if (E < tol) { ok = true; break; }
      x0 = x1; f0 = f1; x1 = x2; f1 = f2;
    }
    steps.push({ title: ok ? "Answer" : "Maximum iterations reached", html:
      (ok ? `Stopped at \\(k=${k}\\) because \\(E_{${k + 1}}\\lt\\varepsilon\\).` : `Stopped after ${nmax} iterations.`) +
      ` <b>Root \\(x\\approx${F(x2, d)}\\).</b><br>Check: \\(f(${F(x2, d)})=${F(f2, d)}\\approx0\\). ✓ No derivative was needed.` + TRUE(Fm) });
    return { given, header: ["k", "\\(x_{k-1}\\)", "\\(x_k\\)", "\\(f(x_{k-1})\\)", "\\(f(x_k)\\)", "\\(x_{k+1}\\)", "\\(f(x_{k+1})\\)", "\\(E_{k+1}\\)", "\\(E\\lt\\varepsilon?\\)"], steps };
  }
  genSecant.title = "secant method";

  /* ---------- integer 3×3 systems built as A = L·U (so no fractions appear) ---------- */
  function luSystem() {
    const nz = () => pick([-3, -2, -1, 1, 2, 3]);
    for (;;) {
      const L = [[1, 0, 0], [nz(), 1, 0], [nz(), nz(), 1]];
      const U = [[pick([1, 2, 3, -1, -2]), rnd(-4, 4), rnd(-4, 4)], [0, pick([1, 2, 3, -1, -2, -3]), rnd(-4, 4)], [0, 0, pick([1, 2, 3, 4, -1, -2, -3])]];
      const A = mul(L, U);
      if (A.flat().filter(v => v === 0).length > 2) continue;
      const x = [rnd(-4, 4), rnd(-4, 4), rnd(-4, 4)]; if (x.every(v => v === 0)) x[1] = 2;
      return { L, U, A, x, b: mv(A, x) };
    }
  }
  /* one elimination row-operation, written out entry by entry */
  function rowOp(rowTarget, m, rowPivot, label, pivotLabel) {
    const parts = rowTarget.map((v, j) => `${v}-(${m})(${rowPivot[j]})`);
    const res = rowTarget.map((v, j) => v - m * rowPivot[j]);
    const n = rowTarget.length;
    const wrap = arr => `\\big[\\,${arr.slice(0, n - 1).join(",\\ ")}\\ \\big|\\ ${arr[n - 1]}\\,\\big]`;
    return { tex: `${label}\\leftarrow ${label}-(${m})${pivotLabel}:\\)<br>\\(${wrap(parts)}\\)<br>\\(=${wrap(res)}`, res };
  }

  /* ========== 2.6 GAUSSIAN ELIMINATION ========== */
  function genGauss() {
    const { A, b, x } = luSystem(); const V = ["x", "y", "z"];
    const given = `Solve by <b>Gaussian elimination</b> (forward elimination, then back substitution; no row swaps are needed):<br>
      \\(${lin(A[0], V)}=${b[0]}\\),   \\(${lin(A[1], V)}=${b[1]}\\),   \\(${lin(A[2], V)}=${b[2]}\\).<br>All arithmetic is exact (integers).`;
    const steps = [];
    const M = A.map((r, i) => r.concat([b[i]]));
    steps.push({ title: "Write the augmented matrix [A | b]", html: `\\[${aug(A, b)}\\]` });
    const m21 = M[1][0] / M[0][0];
    const o1 = rowOp(M[1], m21, M[0], "R_2", "R_1"); M[1] = o1.res;
    steps.push({ title: "Column 1 — eliminate row 2", html:
      `Multiplier \\(m_{21}=\\dfrac{a_{21}}{a_{11}}=\\dfrac{${P(A[1][0])}}{${P(A[0][0])}}=${m21}\\).<br>\\(${o1.tex}\\)<br>\\[${aug(M.map(r => r.slice(0, 3)), M.map(r => r[3]))}\\]` });
    const m31 = M[2][0] / M[0][0];
    const o2 = rowOp(M[2], m31, M[0], "R_3", "R_1"); M[2] = o2.res;
    steps.push({ title: "Column 1 — eliminate row 3", html:
      `Multiplier \\(m_{31}=\\dfrac{a_{31}}{a_{11}}=\\dfrac{${P(A[2][0])}}{${P(A[0][0])}}=${m31}\\).<br>\\(${o2.tex}\\)<br>\\[${aug(M.map(r => r.slice(0, 3)), M.map(r => r[3]))}\\]` });
    const m32 = M[2][1] / M[1][1];
    const o3 = rowOp(M[2], m32, M[1], "R_3", "R_2"); M[2] = o3.res;
    steps.push({ title: "Column 2 — eliminate row 3 (now upper triangular)", html:
      `Multiplier \\(m_{32}=\\dfrac{a'_{32}}{a'_{22}}=\\dfrac{${P(o2.res[1])}}{${P(M[1][1])}}=${m32}\\).<br>\\(${o3.tex}\\)<br>\\[${aug(M.map(r => r.slice(0, 3)), M.map(r => r[3]))}\\]Forward elimination is complete.` });
    const U = M.map(r => r.slice(0, 3)), c = M.map(r => r[3]);
    const z = c[2] / U[2][2];
    steps.push({ title: "Back substitution — row 3 gives z", html: `\\(${lin([0, 0, U[2][2]], V)}=${c[2]}\\;\\Rightarrow\\;z=\\dfrac{${c[2]}}{${U[2][2]}}=${z}\\)` });
    const y = (c[1] - U[1][2] * z) / U[1][1];
    steps.push({ title: "Back substitution — row 2 gives y", html:
      `\\(${lin([0, U[1][1], U[1][2]], V)}=${c[1]}\\;\\Rightarrow\\;y=\\dfrac{${c[1]}-(${U[1][2]})(${z})}{${U[1][1]}}=\\dfrac{${c[1] - U[1][2] * z}}{${U[1][1]}}=${y}\\)` });
    const xx = (c[0] - U[0][1] * y - U[0][2] * z) / U[0][0];
    steps.push({ title: "Back substitution — row 1 gives x", html:
      `\\(${lin(U[0], V)}=${c[0]}\\;\\Rightarrow\\;x=\\dfrac{${c[0]}-(${U[0][1]})(${y})-(${U[0][2]})(${z})}{${U[0][0]}}=\\dfrac{${c[0] - U[0][1] * y - U[0][2] * z}}{${U[0][0]}}=${xx}\\)` });
    steps.push({ title: "Answer and check", html:
      `<b>\\((x,\\,y,\\,z)=(${xx},\\,${y},\\,${z})\\)</b><br>Check in the original equations: ` +
      A.map((r, i) => `\\(${subst(r, [xx, y, z])}=${b[i]}\\)`).join(",  ") + ` ✓` });
    return { given, steps };
  }
  genGauss.title = "Gaussian elimination";

  /* ========== 2.7 LU FACTORIZATION ========== */
  function genLU() {
    const { A, b, x, L: Lt, U: Ut } = luSystem();
    const b2 = mv(A, [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
    const given = `Factor \\(A=LU\\) (Doolittle: ones on the diagonal of \\(L\\), no pivoting), then solve \\(A\\mathbf{x}=\\mathbf{b}\\) by two triangular solves, and finally re-use \\(L,U\\) for a second right-hand side \\(\\mathbf{b}'\\).<br>
      \\(A=${bm(A)},\\quad \\mathbf{b}=${col(b)},\\quad \\mathbf{b}'=${col(b2)}\\).<br>All arithmetic is exact (integers).`;
    const steps = [];
    const M = A.map(r => r.slice());
    const l21 = M[1][0] / M[0][0], l31 = M[2][0] / M[0][0];
    steps.push({ title: "Column 1 — the multipliers", html:
      `Pivot \\(u_{11}=a_{11}=${A[0][0]}\\).  \\(l_{21}=\\dfrac{a_{21}}{u_{11}}=\\dfrac{${P(A[1][0])}}{${P(A[0][0])}}=${l21}\\),  \\(l_{31}=\\dfrac{a_{31}}{u_{11}}=\\dfrac{${P(A[2][0])}}{${P(A[0][0])}}=${l31}\\).<br>Unlike plain elimination, we <b>keep</b> these — they become column 1 of \\(L\\).` });
    const wrap3 = arr => `\\big[\\,${arr.join(",\\ ")}\\,\\big]`;
    const op = (row, m, piv) => ({ parts: row.map((v, j) => `${v}-(${m})(${piv[j]})`), res: row.map((v, j) => v - m * piv[j]) });
    const o1 = op(M[1], l21, M[0]); M[1] = o1.res; const o2 = op(M[2], l31, M[0]); M[2] = o2.res;
    steps.push({ title: "Column 1 — row operations", html:
      `\\(R_2\\leftarrow R_2-(${l21})R_1:\\)<br>\\(${wrap3(o1.parts)}=${wrap3(o1.res)}\\)<br>
       \\(R_3\\leftarrow R_3-(${l31})R_1:\\)<br>\\(${wrap3(o2.parts)}=${wrap3(o2.res)}\\)<br>\\[${bm(M)}\\]` });
    const l32 = M[2][1] / M[1][1];
    const o3 = op(M[2], l32, M[1]); M[2] = o3.res;
    steps.push({ title: "Column 2 — multiplier and row operation", html:
      `Pivot \\(u_{22}=${M[1][1]}\\).  \\(l_{32}=\\dfrac{${P(o2.res[1])}}{${P(M[1][1])}}=${l32}\\).<br>
       \\(R_3\\leftarrow R_3-(${l32})R_2:\\)<br>\\(${wrap3(o3.parts)}=${wrap3(o3.res)}\\)<br>\\[U=${bm(M)}\\]` });
    const L = [[1, 0, 0], [l21, 1, 0], [l31, l32, 1]], U = M;
    steps.push({ title: "Assemble L and U, and verify LU = A", html:
      `\\[L=${bm(L)},\\qquad U=${bm(U)}\\]\\[LU=${bm(L)}${bm(U)}=${bm(mul(L, U))}=A\\ ✓\\]` });
    const y1 = b[0], y2 = b[1] - l21 * y1, y3 = b[2] - l31 * y1 - l32 * y2;
    steps.push({ title: "Forward solve L y = b — row 1", html: `\\(y_1=b_1=${y1}\\)` });
    steps.push({ title: "Forward solve L y = b — row 2", html: `\\(${P(l21)}y_1+y_2=${b[1]}\\;\\Rightarrow\\;y_2=${b[1]}-(${l21})(${y1})=${y2}\\)` });
    steps.push({ title: "Forward solve L y = b — row 3", html: `\\(${P(l31)}y_1+${P(l32)}y_2+y_3=${b[2]}\\;\\Rightarrow\\;y_3=${b[2]}-(${l31})(${y1})-(${l32})(${y2})=${y3}\\)<br>So \\(\\mathbf{y}=${col([y1, y2, y3])}\\).` });
    const x3 = y3 / U[2][2], x2 = (y2 - U[1][2] * x3) / U[1][1], x1 = (y1 - U[0][1] * x2 - U[0][2] * x3) / U[0][0];
    steps.push({ title: "Backward solve U x = y — row 3", html: `\\(${lin([0, 0, U[2][2]], ["x_1", "x_2", "x_3"])}=${y3}\\;\\Rightarrow\\;x_3=\\dfrac{${y3}}{${U[2][2]}}=${x3}\\)` });
    steps.push({ title: "Backward solve U x = y — row 2", html: `\\(${lin([0, U[1][1], U[1][2]], ["x_1", "x_2", "x_3"])}=${y2}\\;\\Rightarrow\\;x_2=\\dfrac{${y2}-(${U[1][2]})(${x3})}{${U[1][1]}}=${x2}\\)` });
    steps.push({ title: "Backward solve U x = y — row 1", html: `\\(${lin(U[0], ["x_1", "x_2", "x_3"])}=${y1}\\;\\Rightarrow\\;x_1=\\dfrac{${y1}-(${U[0][1]})(${x2})-(${U[0][2]})(${x3})}{${U[0][0]}}=${x1}\\)` });
    steps.push({ title: "Answer and check", html:
      `<b>\\(\\mathbf{x}=${col([x1, x2, x3])}\\)</b>.  Check \\(A\\mathbf{x}\\): ` + A.map((r, i) => `\\(${subst(r, [x1, x2, x3])}=${b[i]}\\)`).join(",  ") + ` ✓` });
    const w1 = b2[0], w2 = b2[1] - l21 * w1, w3 = b2[2] - l31 * w1 - l32 * w2;
    const v3 = w3 / U[2][2], v2 = (w2 - U[1][2] * v3) / U[1][1], v1 = (w1 - U[0][1] * v2 - U[0][2] * v3) / U[0][0];
    steps.push({ title: "New right-hand side b′ — forward solve L y′ = b′ (no new elimination!)", html:
      `\\(y'_1=${w1}\\),  \\(y'_2=${b2[1]}-(${l21})(${w1})=${w2}\\),  \\(y'_3=${b2[2]}-(${l31})(${w1})-(${l32})(${w2})=${w3}\\)` });
    steps.push({ title: "New right-hand side b′ — backward solve U x′ = y′", html:
      `\\(x'_3=\\dfrac{${w3}}{${U[2][2]}}=${v3}\\),  \\(x'_2=\\dfrac{${w2}-(${U[1][2]})(${v3})}{${U[1][1]}}=${v2}\\),  \\(x'_1=\\dfrac{${w1}-(${U[0][1]})(${v2})-(${U[0][2]})(${v3})}{${U[0][0]}}=${v1}\\)<br>
       <b>\\(\\mathbf{x}'=${col([v1, v2, v3])}\\)</b> — two cheap triangular solves (≈ \\(2n^2\\) operations) instead of a fresh elimination (≈ \\(\\tfrac23n^3\\)).` });
    return { given, steps };
  }
  genLU.title = "LU factorization";

  /* ========== 2.8 CHOLESKY ========== */
  function genChol() {
    let L;
    for (;;) {
      L = [[pick([1, 2, 2, 3, 3, 4]), 0, 0], [rnd(-3, 3), pick([1, 2, 2, 3]), 0], [rnd(-3, 3), rnd(-3, 3), pick([1, 2, 2, 3])]];
      if ([L[1][0], L[2][0], L[2][1]].filter(v => v !== 0).length >= 2) break;
    }
    const A = mul(L, tr(L));
    const x = [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]; if (x.every(v => v === 0)) x[0] = 1;
    const b = mv(A, x);
    const given = `Solve \\(A\\mathbf{x}=\\mathbf{b}\\) by <b>Cholesky factorization</b> \\(A=LL^{\\mathsf T}\\), then the two triangular solves \\(L\\mathbf{y}=\\mathbf{b}\\), \\(L^{\\mathsf T}\\mathbf{x}=\\mathbf{y}\\).<br>
      \\(A=${bm(A)},\\qquad \\mathbf{b}=${col(b)}\\).<br>All square roots are exact here (integers).`;
    const steps = [];
    steps.push({ title: "Check that A is symmetric", html:
      `\\(a_{12}=a_{21}=${A[0][1]}\\),  \\(a_{13}=a_{31}=${A[0][2]}\\),  \\(a_{23}=a_{32}=${A[1][2]}\\) ⇒ \\(A=A^{\\mathsf T}\\). ✓  (Positive-definiteness is confirmed if every square root below is of a positive number.)` });
    const l11 = L[0][0], l21 = L[1][0], l31 = L[2][0], l22 = L[1][1], l32 = L[2][1], l33 = L[2][2];
    steps.push({ title: "Column 1 — the diagonal entry", html: `\\(l_{11}=\\sqrt{a_{11}}=\\sqrt{${A[0][0]}}=${l11}\\)  (positive ✓)` });
    steps.push({ title: "Column 1 — below the diagonal", html:
      `\\(l_{21}=\\dfrac{a_{21}}{l_{11}}=\\dfrac{${P(A[1][0])}}{${l11}}=${l21}\\),  \\(l_{31}=\\dfrac{a_{31}}{l_{11}}=\\dfrac{${P(A[2][0])}}{${l11}}=${l31}\\)` });
    steps.push({ title: "Column 2 — the diagonal entry", html:
      `\\(l_{22}=\\sqrt{a_{22}-l_{21}^{2}}=\\sqrt{${A[1][1]}-(${l21})^{2}}=\\sqrt{${A[1][1] - l21 * l21}}=${l22}\\)  (positive ✓)` });
    steps.push({ title: "Column 2 — below the diagonal", html:
      `\\(l_{32}=\\dfrac{a_{32}-l_{31}l_{21}}{l_{22}}=\\dfrac{${A[2][1]}-(${l31})(${l21})}{${l22}}=\\dfrac{${A[2][1] - l31 * l21}}{${l22}}=${l32}\\)` });
    steps.push({ title: "Column 3 — the last diagonal entry", html:
      `\\(l_{33}=\\sqrt{a_{33}-l_{31}^{2}-l_{32}^{2}}=\\sqrt{${A[2][2]}-(${l31})^{2}-(${l32})^{2}}=\\sqrt{${A[2][2] - l31 * l31 - l32 * l32}}=${l33}\\)  (positive ✓ ⇒ \\(A\\) is positive definite)` });
    steps.push({ title: "Assemble L and verify L Lᵀ = A", html:
      `\\[L=${bm(L)},\\qquad L^{\\mathsf T}=${bm(tr(L))}\\]\\[LL^{\\mathsf T}=${bm(mul(L, tr(L)))}=A\\ ✓\\]` });
    const y1 = b[0] / l11, y2 = (b[1] - l21 * y1) / l22, y3 = (b[2] - l31 * y1 - l32 * y2) / l33;
    steps.push({ title: "Forward solve L y = b", html:
      `\\(y_1=\\dfrac{b_1}{l_{11}}=\\dfrac{${b[0]}}{${l11}}=${y1}\\)<br>
       \\(y_2=\\dfrac{b_2-l_{21}y_1}{l_{22}}=\\dfrac{${b[1]}-(${l21})(${y1})}{${l22}}=${y2}\\)<br>
       \\(y_3=\\dfrac{b_3-l_{31}y_1-l_{32}y_2}{l_{33}}=\\dfrac{${b[2]}-(${l31})(${y1})-(${l32})(${y2})}{${l33}}=${y3}\\)` });
    const x3 = y3 / l33, x2 = (y2 - l32 * x3) / l22, x1 = (y1 - l21 * x2 - l31 * x3) / l11;
    steps.push({ title: "Backward solve Lᵀ x = y", html:
      `Rows of \\(L^{\\mathsf T}\\): \\((l_{11},l_{21},l_{31})\\), \\((0,l_{22},l_{32})\\), \\((0,0,l_{33})\\).<br>
       \\(x_3=\\dfrac{y_3}{l_{33}}=\\dfrac{${y3}}{${l33}}=${x3}\\)<br>
       \\(x_2=\\dfrac{y_2-l_{32}x_3}{l_{22}}=\\dfrac{${y2}-(${l32})(${x3})}{${l22}}=${x2}\\)<br>
       \\(x_1=\\dfrac{y_1-l_{21}x_2-l_{31}x_3}{l_{11}}=\\dfrac{${y1}-(${l21})(${x2})-(${l31})(${x3})}{${l11}}=${x1}\\)` });
    steps.push({ title: "Answer and check", html:
      `<b>\\(\\mathbf{x}=${col([x1, x2, x3])}\\)</b>.  Check \\(A\\mathbf{x}\\): ` + A.map((r, i) => `\\(${subst(r, [x1, x2, x3])}=${b[i]}\\)`).join(",  ") + ` ✓` });
    return { given, steps };
  }
  genChol.title = "Cholesky factorization";

  const GEN = { fixed: genFixed, bisect: genBisect, newton: genNewton, secant: genSecant, gauss: genGauss, lu: genLU, chol: genChol };
  window.BCS311_GEN = GEN;   // exposed for testing

  /* ---------- MathJax (loaded async) ---------- */
  function typeset(el) {
    const go = () => MathJax.typesetPromise([el]).catch(() => {});
    if (window.MathJax && MathJax.typesetPromise) { go(); return; }
    let tries = 0;
    const t = setInterval(() => {
      if (window.MathJax && MathJax.typesetPromise) { clearInterval(t); go(); }
      else if (++tries > 80) clearInterval(t);
    }, 250);
  }

  /* ---------- the widget ---------- */
  function mount(host) {
    const gen = GEN[host.dataset.gen];
    if (!gen) return;
    host.classList.add("genp");
    host.innerHTML =
      '<div class="genp-head"><span class="ptag">🎲 Random practice · ' + gen.title + '</span>' +
      '<button class="btn accent genp-new">🎲 New numbers</button></div>' +
      '<p class="genp-hint">Work the problem on paper. After each step, press <b>Check next step</b>: only that step is revealed, so you can compare before going on.</p>' +
      '<div class="given genp-given"></div>' +
      '<div class="tablewrap genp-table" hidden><table class="data"><thead></thead><tbody></tbody></table></div>' +
      '<div class="genp-steps"></div>' +
      '<div class="btnrow"><button class="btn genp-next"></button><button class="btn ghost genp-all">Reveal all steps</button><button class="btn ghost genp-hide">Hide steps</button></div>';
    const given = host.querySelector(".genp-given"), stepsEl = host.querySelector(".genp-steps"),
          tableWrap = host.querySelector(".genp-table"), thead = tableWrap.querySelector("thead"), tbody = tableWrap.querySelector("tbody"),
          btnNext = host.querySelector(".genp-next"), btnAll = host.querySelector(".genp-all"), btnHide = host.querySelector(".genp-hide"), btnNew = host.querySelector(".genp-new");
    let prob, shown = 0, headTypeset = false;

    function updateBtns() {
      const N = prob.steps.length, done = shown >= N;
      btnNext.textContent = done ? "All " + N + " steps shown ✓" : "Check next step (" + (shown + 1) + " of " + N + ") ▸";
      btnNext.disabled = done; btnAll.disabled = done; btnHide.disabled = shown === 0;
    }
    function fresh() {
      prob = gen(); shown = 0; headTypeset = false;
      given.innerHTML = prob.given; stepsEl.innerHTML = ""; tbody.innerHTML = ""; tableWrap.hidden = true;
      thead.innerHTML = prob.header ? "<tr>" + prob.header.map(h => "<th>" + h + "</th>").join("") + "</tr>" : "";
      typeset(given); updateBtns();
    }
    function reveal(scroll) {
      if (shown >= prob.steps.length) return;
      const st = prob.steps[shown], last = shown === prob.steps.length - 1;
      const card = document.createElement("div");
      card.className = "genp-step" + (last ? " final" : "");
      card.innerHTML = '<div class="t"><span class="k">Step ' + (shown + 1) + '</span><span>' + st.title + '</span></div><div class="b">' + st.html + '</div>';
      stepsEl.appendChild(card);
      if (st.row) {
        if (tableWrap.hidden) { tableWrap.hidden = false; }
        if (!headTypeset) { headTypeset = true; typeset(thead); }
        const trEl = document.createElement("tr");
        if (st.rowClass) trEl.className = st.rowClass;
        trEl.innerHTML = st.row.map(c => "<td>" + c + "</td>").join("");
        tbody.appendChild(trEl);
      }
      shown++; typeset(card); updateBtns();
      if (scroll) card.scrollIntoView({ block: "nearest", behavior: "smooth" });
    }
    btnNext.addEventListener("click", () => reveal(true));
    btnAll.addEventListener("click", () => { while (shown < prob.steps.length) reveal(false); });
    btnHide.addEventListener("click", () => { shown = 0; stepsEl.innerHTML = ""; tbody.innerHTML = ""; tableWrap.hidden = true; updateBtns(); });
    btnNew.addEventListener("click", fresh);
    fresh();
  }

  document.addEventListener("DOMContentLoaded", () => {
    document.querySelectorAll(".gen-practice[data-gen]").forEach(mount);
  });
})();
