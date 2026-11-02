/* ===========================================================
   BCS 311 — Modules 2 & 3 | practice-gen.js
   Random, step-by-step practice problems.
   Usage:  <div class="gen-practice" data-gen="newton"></div>
   Keys:   fixed · bisect · newton · secant · gauss · lu · chol   (Module 2)
           linreg · rsq · linz · quad · polyfit                (Module 3)
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
      let c, a;                                       // skip c = a³+a (10, 30): the root would be the integer endpoint itself
      do { c = rnd(3, 40); a = 0; while ((a + 1) ** 3 + (a + 1) <= c) a++; } while (a ** 3 + a === c);
      const f = x => x ** 3 + x - c;
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
      const left = fa * fc <= 0;                      // sign change in [a, c]?  (fc = 0 exactly: keep c as the right end)
      const na = left ? a : c, nb = left ? c : b;
      const verdict = fc === 0 ? "<b>zero</b> ⇒ \\(c\\) is a root ⇒ keep \\([a,c]\\), set \\(b\\leftarrow c\\)"
                    : left ? "<b>negative</b> ⇒ the root is in \\([a,c]\\) ⇒ set \\(b\\leftarrow c\\)"
                    : "<b>positive</b> ⇒ the root is in \\([c,b]\\) ⇒ set \\(a\\leftarrow c\\)";
      steps.push({ title: `Iteration n = ${n}`, html:
        `\\(c=\\tfrac{a+b}{2}=\\tfrac{${F(a, d)}+${F(b, d)}}{2}=${F(c, d)}\\)<br>
         \\(f(c)=${Fm.fSub(F(c, d))}=${S(fc, d)}\\)<br>
         \\(f(a)\\,f(c)=(${S(fa, d)})(${S(fc, d)})\\) is ${verdict}<br>
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

  /* =====================================================================
     MODULE 3 — least-squares curve fitting
     Keys: linreg · rsq · linz · quad · polyfit
     ===================================================================== */
  const tbl = (head, rows, foot) =>
    "\\[\\begin{array}{" + "r".repeat(head.length) + "}" + head.join("&") + "\\\\\\hline " +
    rows.map(r => r.join("&")).join("\\\\") + (foot ? "\\\\\\hline " + foot.join("&") : "") + "\\end{array}\\]";
  const dataGiven = (xs, ys, xl, yl) =>
    "\\[\\begin{array}{c|" + "r".repeat(xs.length) + "}" + xl + "&" + xs.join("&") + "\\\\" + yl + "&" + ys.join("&") + "\\end{array}\\]";
  const sum = a => a.reduce((s, v) => s + v, 0);
  const lineTex = (a0, a1, d) => `${F(a0, d)}${a1 < 0 ? "-" : "+"}${F(Math.abs(a1), d)}\\,x`;

  /* ========== 3.2 LINEAR REGRESSION ========== */
  function genLinReg() {
    const d = 6;
    const T = v => String(parseFloat(F(v, d)));   // trimmed: 8, 1.714286
    let n, xs, ys, a1, num, den;
    for (;;) {
      n = rnd(5, 7); const x0 = pick([0, 1, 1, 2]), dx = pick([1, 1, 2]);
      xs = Array.from({ length: n }, (_, i) => x0 + dx * i);
      const m = pick([0.5, 1, 1.5, 2, 3, -1, -1.5]), c = m > 0 ? rnd(0, 5) : rnd(12, 18);
      ys = xs.map(x => Math.round(c + m * x + rnd(-2, 2)));
      num = n * sum(xs.map((x, i) => x * ys[i])) - sum(xs) * sum(ys);
      den = n * sum(xs.map(x => x * x)) - sum(xs) ** 2;
      a1 = R(num / den, d);
      if (ys.every(y => y >= 0) && new Set(ys).size > 2 && Math.abs(a1) >= 0.3) break;
    }
    const Sx = sum(xs), Sy = sum(ys), Sxy = sum(xs.map((x, i) => x * ys[i])), Sxx = sum(xs.map(x => x * x));
    const given = `Fit a straight line \\(y=a_0+a_1x\\) to the data by <b>linear least squares</b> (normal-equation formulas), working to ${d} decimals (round \\(a_1\\), \\(\\bar x\\), \\(\\bar y\\), \\(a_0\\) and carry the rounded values). Then predict \\(y\\) at \\(x=${xs[n - 1] + (xs[1] - xs[0])}\\).
      ${dataGiven(xs, ys, "x", "y")}`;
    const steps = [];
    steps.push({ title: "Tabulate the four sums", html:
      tbl(["i", "x_i", "y_i", "x_iy_i", "x_i^2"], xs.map((x, i) => [i + 1, x, ys[i], x * ys[i], x * x]), ["", Sx, Sy, Sxy, Sxx]) +
      `\\(n=${n},\\ \\sum x_i=${Sx},\\ \\sum y_i=${Sy},\\ \\sum x_iy_i=${Sxy},\\ \\sum x_i^2=${Sxx}\\).` });
    steps.push({ title: "Slope a₁", html:
      `\\(a_1=\\dfrac{n\\sum x_iy_i-\\sum x_i\\sum y_i}{n\\sum x_i^2-\\left(\\sum x_i\\right)^2}=\\dfrac{${n}(${Sxy})-(${Sx})(${Sy})}{${n}(${Sxx})-(${Sx})^2}=\\dfrac{${num}}{${den}}=${F(a1, d)}\\)` });
    const xb = R(Sx / n, d), yb = R(Sy / n, d), a0 = R(yb - a1 * xb, d);
    steps.push({ title: "Means and intercept a₀", html:
      `\\(\\bar x=\\dfrac{${Sx}}{${n}}=${T(xb)}\\),  \\(\\bar y=\\dfrac{${Sy}}{${n}}=${T(yb)}\\)<br>
       \\(a_0=\\bar y-a_1\\bar x=${T(yb)}-(${F(a1, d)})(${T(xb)})=${F(a0, d)}\\)` });
    const xq = xs[n - 1] + (xs[1] - xs[0]), yq = R(a0 + a1 * xq, 4);
    steps.push({ title: "Answer and prediction", html:
      `<b>\\(y=${lineTex(a0, a1, d)}\\)</b>.<br>
       Sanity check: the line passes through the centroid, \\(a_0+a_1\\bar x=${F(a0, d)}+(${F(a1, d)})(${T(xb)})=${T(R(a0 + a1 * xb, d))}\\approx\\bar y=${T(yb)}\\). ✓<br>
       Prediction: \\(y(${xq})=${F(a0, d)}+(${F(a1, d)})(${xq})=${F(yq, 4)}\\).` });
    return { given, steps };
  }
  genLinReg.title = "linear regression";

  /* ========== 3.3 HOW GOOD IS THE FIT? ========== */
  function genRsq() {
    const d = 4;
    let n, xs, ys, a0, a1, yh, e, e2, Sr, St, yb, r2, s;
    for (;;) {
      n = pick([4, 5]); const x0 = pick([0, 1]);
      xs = Array.from({ length: n }, (_, i) => x0 + i);
      const m = pick([1, 1.5, 2, 3, -1, -2]), c = m > 0 ? rnd(0, 4) : rnd(10, 14);
      ys = xs.map(x => Math.round(c + m * x + rnd(-2, 2)));
      const Sx = sum(xs), Sy = sum(ys), Sxy = sum(xs.map((x, i) => x * ys[i])), Sxx = sum(xs.map(x => x * x));
      a1 = R((n * Sxy - Sx * Sy) / (n * Sxx - Sx * Sx), 2);        // exact: den is 20 or 50
      a0 = R((Sy - a1 * Sx) / n, d);
      yh = xs.map(x => R(a0 + a1 * x, d));
      e = ys.map((y, i) => R(y - yh[i], d)); e2 = e.map(v => R(v * v, d)); Sr = R(sum(e2), d);
      yb = R(Sy / n, 4); const dv = ys.map(y => R(y - yb, d)), dv2 = dv.map(v => R(v * v, d));
      St = R(sum(dv2), d); r2 = R((St - Sr) / St, d);
      if (St > 0 && Sr > 0 && r2 >= 0.6 && r2 < 0.999 && ys.every(y => y >= 0)) { var DV = dv, DV2 = dv2; break; }
    }
    s = R(Math.sqrt(R(Sr / (n - 2), d)), d);
    const given = `For the data below, the least-squares line is \\(\\hat y=${lineTex(a0, a1, 4)}\\). Find the sum of squared residuals \\(S_r\\), the spread about the mean \\(S_t\\), the coefficient of determination \\(r^2\\) and the standard error \\(s_{y/x}\\). Work to ${d} decimals and carry the rounded values.
      ${dataGiven(xs, ys, "x", "y")}`;
    const steps = [];
    steps.push({ title: "Predictions and residuals about the line", html:
      tbl(["x_i", "y_i", "\\hat y_i=" + F(a0, d) + (a1 < 0 ? "-" : "+") + F(Math.abs(a1), d) + "x_i", "e_i=y_i-\\hat y_i", "e_i^2"],
          xs.map((x, i) => [x, ys[i], F(yh[i], d), F(e[i], d), F(e2[i], d)]), ["", "", "", "", F(Sr, d)]) +
      `<b>\\(S_r=\\sum e_i^2=${F(Sr, d)}\\)</b>` });
    steps.push({ title: "Spread about the mean — S_t", html:
      `\\(\\bar y=\\dfrac{${sum(ys)}}{${n}}=${F(yb, d)}\\)<br>` +
      tbl(["y_i", "y_i-\\bar y", "(y_i-\\bar y)^2"], ys.map((y, i) => [y, F(DV[i], d), F(DV2[i], d)]), ["", "", F(St, d)]) +
      `<b>\\(S_t=\\sum (y_i-\\bar y)^2=${F(St, d)}\\)</b>` });
    steps.push({ title: "Coefficient of determination r²", html:
      `\\(r^2=\\dfrac{S_t-S_r}{S_t}=\\dfrac{${F(St, d)}-${F(Sr, d)}}{${F(St, d)}}=${F(r2, d)}\\)<br>
       The line explains ${F(r2 * 100, 2)}% of the variation in \\(y\\); correlation \\(r=\\sqrt{r^2}=${F(Math.sqrt(r2), d)}\\).` });
    steps.push({ title: "Standard error of the estimate", html:
      `\\(s_{y/x}=\\sqrt{\\dfrac{S_r}{n-2}}=\\sqrt{\\dfrac{${F(Sr, d)}}{${n - 2}}}=\\sqrt{${F(R(Sr / (n - 2), d), d)}}=${F(s, d)}\\)<br>
       <b>Answer:</b> \\(S_r=${F(Sr, d)}\\), \\(S_t=${F(St, d)}\\), \\(r^2=${F(r2, d)}\\), \\(s_{y/x}=${F(s, d)}\\). A typical point lies about ${F(s, 2)} units from the line.` });
    return { given, steps };
  }
  genRsq.title = "r² and standard error";

  /* ========== 3.4 LINEARIZING (power / exponential) ========== */
  function genLinz() {
    const d = 4, pow = Math.random() < 0.5;
    let xs, ys, X, Y, XY, XX, a, b;
    const yfmt = v => v >= 100 ? Math.round(v) : v >= 10 ? R(v, 1) : R(v, 2);
    for (;;) {
      if (pow) {
        const A = pick([0.5, 1, 2, 3, 5]), B = pick([0.5, 1, 1.5, 2, 2.5]); xs = pick([[1, 2, 3, 4, 5], [2, 4, 6, 8, 10], [10, 20, 30, 40, 50]]);
        ys = xs.map(x => yfmt(A * x ** B * (1 + rnd(-10, 10) / 100)));
      } else {
        const A = pick([2, 3, 5, 8]), B = pick([0.2, 0.3, 0.4, -0.2, -0.3]); xs = pick([[0, 1, 2, 3, 4], [1, 2, 3, 4, 5]]);
        ys = xs.map(x => yfmt(A * Math.exp(B * x) * (1 + rnd(-10, 10) / 100)));
      }
      if (!ys.every(y => y >= 0.5)) continue;
      /* reject draws where carrying 4-decimal values drifts from the exact least-squares fit (ill-conditioned rounding) */
      const n0 = xs.length, lg0 = pow ? Math.log10 : Math.log, X0 = pow ? xs.map(Math.log10) : xs, Y0 = ys.map(lg0);
      const mx = sum(X0) / n0, my = sum(Y0) / n0;
      const bx = sum(X0.map((v, i) => (v - mx) * (Y0[i] - my))) / sum(X0.map(v => (v - mx) ** 2)), ax = pow ? Math.pow(10, my - bx * mx) : Math.exp(my - bx * mx);
      const Xr = pow ? xs.map(x => R(Math.log10(x), d)) : xs, Yr = ys.map(y => R(lg0(y), d));
      const sX = R(sum(Xr), d), sY = R(sum(Yr), d), sXY = R(sum(Xr.map((v, i) => R(v * Yr[i], d))), d), sXX = R(sum(Xr.map(v => R(v * v, d))), d);
      const bb = R(R(n0 * sXY - sX * sY, d) / R(n0 * sXX - sX * sX, d), d);
      const aa = R(pow ? Math.pow(10, R(R(sY / n0, d) - bb * R(sX / n0, d), d)) : Math.exp(R(R(sY / n0, d) - bb * R(sX / n0, d), d)), d);
      if (Math.abs(bb - bx) < 5e-4 && Math.abs(aa - ax) < 1e-3 * Math.max(1, ax)) break;
    }
    const n = xs.length, lg = pow ? Math.log10 : Math.log;
    X = pow ? xs.map(x => R(Math.log10(x), d)) : xs.slice();
    Y = ys.map(y => R(lg(y), d));
    XY = X.map((v, i) => R(v * Y[i], d)); XX = X.map(v => R(v * v, d));
    const SX = R(sum(X), d), SY = R(sum(Y), d), SXY = R(sum(XY), d), SXX = R(sum(XX), d);
    const num = R(n * SXY - SX * SY, d), den = R(n * SXX - SX * SX, d);
    b = R(num / den, d);
    const Xb = R(SX / n, d), Yb = R(SY / n, d), A0 = R(Yb - b * Xb, d);
    a = R(pow ? Math.pow(10, A0) : Math.exp(A0), d);
    const L = pow ? "\\log_{10}" : "\\ln", Xn = pow ? "X" : "x", base = pow ? "10" : "e";
    const model = pow ? "y=a\\,x^{b}" : "y=a\\,e^{bx}";
    const given = `Fit the model \\(${model}\\) to the data by <b>linearizing</b> it and using straight-line least squares. Work to ${d} decimals (round every log, product and sum to ${d} decimals and carry the rounded values). Then predict \\(y\\) at \\(x=${pow ? xs[n - 1] + (xs[1] - xs[0]) : xs[n - 1] + 1}\\).
      ${dataGiven(xs, ys, "x", "y")}`;
    const steps = [];
    steps.push({ title: "Linearize: take logs", html: pow
      ? `\\(y=ax^{b}\\Rightarrow \\log_{10}y=\\log_{10}a+b\\,\\log_{10}x\\). Put \\(X=\\log_{10}x\\), \\(Y=\\log_{10}y\\): a straight line \\(Y=A+bX\\) with slope \\(b\\) and intercept \\(A=\\log_{10}a\\).`
      : `\\(y=ae^{bx}\\Rightarrow \\ln y=\\ln a+b\\,x\\). Put \\(Y=\\ln y\\): a straight line \\(Y=A+bx\\) in the original \\(x\\), with slope \\(b\\) and intercept \\(A=\\ln a\\).` });
    steps.push({ title: "Transform the data and tabulate", html: pow
      ? tbl(["x_i", "y_i", "X_i=\\log x_i", "Y_i=\\log y_i", "X_iY_i", "X_i^2"], xs.map((x, i) => [x, ys[i], F(X[i], d), F(Y[i], d), F(XY[i], d), F(XX[i], d)]), ["", "", F(SX, d), F(SY, d), F(SXY, d), F(SXX, d)])
      : tbl(["x_i", "y_i", "Y_i=\\ln y_i", "x_iY_i", "x_i^2"], xs.map((x, i) => [x, ys[i], F(Y[i], d), F(XY[i], d), F(XX[i], d)]), [F(SX, 0), "", F(SY, d), F(SXY, d), F(SXX, d)]) });
    steps.push({ title: "Slope b", html:
      `\\(b=\\dfrac{n\\sum ${Xn}Y-\\sum ${Xn}\\sum Y}{n\\sum ${Xn}^2-\\left(\\sum ${Xn}\\right)^2}=\\dfrac{${n}(${F(SXY, d)})-(${F(SX, d)})(${F(SY, d)})}{${n}(${F(SXX, d)})-(${F(SX, d)})^2}=\\dfrac{${F(num, d)}}{${F(den, d)}}=${F(b, d)}\\)` });
    steps.push({ title: "Intercept A, then undo the log to get a", html:
      `\\(\\bar ${Xn}=${F(Xb, d)}\\),  \\(\\bar Y=${F(Yb, d)}\\) ⇒ \\(A=\\bar Y-b\\bar ${Xn}=${F(Yb, d)}-(${F(b, d)})(${F(Xb, d)})=${F(A0, d)}\\)<br>
       \\(a=${base}^{A}=${base}^{${F(A0, d)}}=${F(a, d)}\\)` });
    const xq = pow ? xs[n - 1] + (xs[1] - xs[0]) : xs[n - 1] + 1;
    const yq = pow ? a * xq ** b : a * Math.exp(b * xq);
    steps.push({ title: "Answer and prediction", html:
      `<b>\\(y=${F(a, d)}\\,${pow ? `x^{${F(b, d)}}` : `e^{${F(b, d)}x}`}\\)</b><br>
       Prediction: \\(y(${xq})=${F(a, d)}${pow ? `(${xq})^{${F(b, d)}}` : `\\,e^{${F(b, d)}\\cdot ${xq}}`}\\approx ${F(yq, 2)}\\).<br>
       Note the fit is least-squares on \\(\\${pow ? "log" : "ln"}\\,y\\), not on \\(y\\) itself, so it weights small and large \\(y\\) differently from a direct fit.` });
    return { given, steps };
  }
  genLinz.title = "linearizing a non-linear model";

  /* ========== 3.5 POLYNOMIAL (QUADRATIC) REGRESSION ========== */
  function genQuad() {
    const d = 4, dc = 6, n = 5;
    let xs, ys;
    for (;;) {
      const s = pick([-2, 0, 1, 2]);
      xs = Array.from({ length: n }, (_, i) => s + i);
      const c2 = pick([1, 1, 2, -1, -2, 0.5]), c1 = rnd(-3, 3), c0 = rnd(0, 6);
      ys = xs.map(x => Math.round(c0 + c1 * x + c2 * x * x + rnd(-1, 1)));
      const Sx = sum(xs), Sxx = sum(xs.map(x => x * x));
      if (new Set(ys).size >= 4 && Sx % n === 0 && Sxx % n === 0) break;
    }
    const S = k => sum(xs.map(x => x ** k)), SY = sum(ys), SXY = sum(xs.map((x, i) => x * ys[i])), SXXY = sum(xs.map((x, i) => x * x * ys[i]));
    const M = [[n, S(1), S(2)], [S(1), S(2), S(3)], [S(2), S(3), S(4)]], B = [SY, SXY, SXXY];
    const given = `Fit a parabola \\(y=a_0+a_1x+a_2x^2\\) to the data by <b>polynomial least squares</b>: build the normal equations, solve them by Gaussian elimination, and give \\(a_0,a_1,a_2\\) to ${d} decimals (carry ${dc} decimals through the back-substitution, then round the answer). Then find the sum of squared residuals \\(S_r\\) and \\(r^2\\).
      ${dataGiven(xs, ys, "x", "y")}`;
    const steps = [];
    steps.push({ title: "Tabulate the sums", html:
      tbl(["x_i", "y_i", "x_i^2", "x_i^3", "x_i^4", "x_iy_i", "x_i^2y_i"], xs.map((x, i) => [x, ys[i], x ** 2, x ** 3, x ** 4, x * ys[i], x * x * ys[i]]),
          [S(1), SY, S(2), S(3), S(4), SXY, SXXY]) });
    steps.push({ title: "Write the normal equations", html:
      `\\[${bm([[n, S(1), S(2)], [S(1), S(2), S(3)], [S(2), S(3), S(4)]])}\\begin{bmatrix}a_0\\\\a_1\\\\a_2\\end{bmatrix}=${col(B)}\\]
       (Linear in the unknowns \\(a_0,a_1,a_2\\) — a 3×3 system exactly like Module 2.)` });
    const m21 = M[1][0] / M[0][0], m31 = M[2][0] / M[0][0];
    const r2 = [0, M[1][1] - m21 * M[0][1], M[1][2] - m21 * M[0][2]], b2 = B[1] - m21 * B[0];
    const r3 = [0, M[2][1] - m31 * M[0][1], M[2][2] - m31 * M[0][2]], b3 = B[2] - m31 * B[0];
    steps.push({ title: "Eliminate column 1", html:
      `\\(m_{21}=\\dfrac{${M[1][0]}}{${n}}=${m21}\\),  \\(m_{31}=\\dfrac{${M[2][0]}}{${n}}=${m31}\\)<br>
       \\(R_2\\leftarrow R_2-(${m21})R_1:\\ [0,\\ ${r2[1]},\\ ${r2[2]}\\ |\\ ${b2}]\\)<br>
       \\(R_3\\leftarrow R_3-(${m31})R_1:\\ [0,\\ ${r3[1]},\\ ${r3[2]}\\ |\\ ${b3}]\\)` });
    const m32 = r3[1] / r2[1], u33 = r3[2] - m32 * r2[2], c3 = b3 - m32 * b2;
    steps.push({ title: "Eliminate column 2", html:
      `\\(m_{32}=\\dfrac{${r3[1]}}{${r2[1]}}=${m32}\\)<br>
       \\(R_3\\leftarrow R_3-(${m32})R_2:\\ [0,\\ 0,\\ ${u33}\\ |\\ ${c3}]\\)` });
    const a2c = R(c3 / u33, dc), a1c = R((b2 - r2[2] * a2c) / r2[1], dc), a0c = R((B[0] - M[0][1] * a1c - M[0][2] * a2c) / n, dc);
    const a2 = R(a2c, d), a1 = R(a1c, d), a0 = R(a0c, d);
    steps.push({ title: "Back-substitution", html:
      `\\(a_2=\\dfrac{${c3}}{${u33}}=${F(a2c, dc)}\\)<br>
       \\(a_1=\\dfrac{${b2}-(${r2[2]})(${F(a2c, dc)})}{${r2[1]}}=${F(a1c, dc)}\\)<br>
       \\(a_0=\\dfrac{${B[0]}-(${M[0][1]})(${F(a1c, dc)})-(${M[0][2]})(${F(a2c, dc)})}{${n}}=${F(a0c, dc)}\\)<br>
       Rounded to ${d} decimals: \\(a_0=${F(a0, d)},\\ a_1=${F(a1, d)},\\ a_2=${F(a2, d)}\\).` });
    const yh = xs.map(x => R(a0c + a1c * x + a2c * x * x, d)), e = ys.map((y, i) => R(y - yh[i], d)), e2 = e.map(v => R(v * v, d));
    const Sr = R(sum(e2), d), yb = R(SY / n, d), dv2 = ys.map(y => R((y - yb) ** 2, d)), St = R(sum(dv2), d), r2v = R((St - Sr) / St, d);
    const poly = `${F(a0, d)}${a1 < 0 ? "-" : "+"}${F(Math.abs(a1), d)}x${a2 < 0 ? "-" : "+"}${F(Math.abs(a2), d)}x^{2}`;
    steps.push({ title: "Answer, residuals and r²", html:
      `<b>\\(y=${poly}\\)</b><br>` +
      tbl(["x_i", "y_i", "\\hat y_i", "e_i=y_i-\\hat y_i", "e_i^2"], xs.map((x, i) => [x, ys[i], F(yh[i], d), F(e[i], d), F(e2[i], d)]), ["", "", "", "", F(Sr, d)]) +
      `\\(\\bar y=${F(yb, d)}\\), \\(S_t=\\sum(y_i-\\bar y)^2=${F(St, d)}\\), so \\(r^2=\\dfrac{S_t-S_r}{S_t}=\\dfrac{${F(St, d)}-${F(Sr, d)}}{${F(St, d)}}=${F(r2v, d)}\\).` });
    return { given, steps };
  }
  genQuad.title = "polynomial regression";

  /* ========== 3.6 MATLAB polyfit / polyval ========== */
  function genPolyfit() {
    const m = pick([1, 2]);
    let xs, c;
    const x0 = pick([0, 1]);
    xs = Array.from({ length: 5 }, (_, i) => x0 + i);
    if (m === 1) c = [pick([1, 2, 3, -1, -2, 4]), rnd(-4, 6)];                      // [a1 a0]
    else c = [pick([1, 2, -1, -2, 3]), rnd(-3, 3), rnd(-4, 6)];                     // [a2 a1 a0]
    const ev = x => c.reduce((s, v) => s * x + v, 0);
    const ys = xs.map(ev);
    const xq = x0 + 6;
    const vec = a => "[" + a.join(" ") + "]";
    const pname = m === 1 ? "line" : "parabola";
    const given = `A student runs this MATLAB code. The data lie <em>exactly</em> on a ${pname}, so there is no scatter.
      <pre><code>x = ${vec(xs)};
y = ${vec(ys)};
p  = polyfit(x, y, ${m});
yq = polyval(p, ${xq});</code></pre>
      What is the vector <code>p</code> (in MATLAB's order), what does <code>yq</code> equal, and what would the sum of squared residuals \\(S_r\\) be?`;
    const dif = a => a.slice(1).map((v, i) => v - a[i]);
    const d1 = dif(ys), d2 = dif(d1);
    const steps = [];
    steps.push({ title: "What kind of curve is it? (finite differences)", html: m === 1
      ? `First differences of \\(y\\): \\(${d1.join(",\\ ")}\\) — constant ⇒ the data are a straight line, slope \\(a_1=${d1[0]}\\) per unit \\(x\\).`
      : `First differences: \\(${d1.join(",\\ ")}\\) (not constant). Second differences: \\(${d2.join(",\\ ")}\\) — constant ⇒ a parabola with \\(2a_2=${d2[0]}\\), so \\(a_2=${d2[0] / 2}\\).` });
    if (m === 1) {
      steps.push({ title: "Intercept from one point", html:
        `At \\(x=${xs[0]}\\): \\(y=${ys[0]}=a_0+${P(c[0])}(${xs[0]})\\) ⇒ \\(a_0=${c[1]}\\).<br>So \\(y=${c[1]}${c[0] < 0 ? "-" : "+"}${Math.abs(c[0])}x\\).` });
    } else {
      const a2 = c[0], r = x => ev(x) - a2 * x * x;
      steps.push({ title: "Peel off the x² term, then read the line", html:
        `Subtract \\(${P(a2)}x^2\\): at \\(x=${xs[0]}\\) get \\(${ys[0]}-${P(a2)}(${xs[0] ** 2})=${r(xs[0])}\\); at \\(x=${xs[1]}\\) get \\(${ys[1]}-${P(a2)}(${xs[1] ** 2})=${r(xs[1])}\\).<br>
         Those two values fix the line \\(a_0+a_1x\\): slope \\(a_1=${r(xs[1]) - r(xs[0])}\\), then \\(a_0=${c[2]}\\).<br>So \\(y=${c[2]}${c[1] < 0 ? "-" : "+"}${Math.abs(c[1])}x${c[0] < 0 ? "-" : "+"}${Math.abs(c[0])}x^2\\).` });
    }
    steps.push({ title: "MATLAB's order: highest power first", html:
      `<code>polyfit</code> returns <b>highest power first</b>: <code>p = ${m === 1 ? "[a1 a0]" : "[a2 a1 a0]"} = ${vec(c)}</code>.<br>
       (MATLAB may print it as <code>${c.map(v => v.toFixed(4)).join("  ")}</code> because floating-point rounding leaves tiny errors around exact integers.)` });
    const yq = ev(xq);
    steps.push({ title: "polyval and S<sub>r</sub>", html:
      `<code>polyval(p, ${xq})</code> evaluates \\(${m === 1 ? `${c[0]}(${xq})+${P(c[1])}` : `${c[0]}(${xq})^2+${P(c[1])}(${xq})+${P(c[2])}`}=${yq}\\).<br>
       <b>Answer:</b> <code>p = ${vec(c)}</code>, <code>yq = ${yq}</code>, and \\(S_r=0\\) (up to floating-point round-off) because every point lies exactly on the fitted ${pname}, so \\(r^2=1\\).` });
    return { given, steps };
  }
  genPolyfit.title = "MATLAB polyfit / polyval";

  const GEN = { fixed: genFixed, bisect: genBisect, newton: genNewton, secant: genSecant, gauss: genGauss, lu: genLU, chol: genChol,
                 linreg: genLinReg, rsq: genRsq, linz: genLinz, quad: genQuad, polyfit: genPolyfit };
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
