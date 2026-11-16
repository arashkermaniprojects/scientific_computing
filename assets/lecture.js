/* BCS 311 — lecture.js : turns a module page into a self-explanatory narrated lecture.
   Reads assets/lecture_m<N>.json  ->  goals box, "Narrate this chapter", layout switch,
   per-section "big picture", "common mistake" and "check your understanding". */
(function () {
  "use strict";
  const m = (location.pathname.match(/module(\d)/) || [])[1]; if (!m) return;
  const $ = (s, r) => (r || document).querySelector(s), $$ = (s, r) => [...(r || document).querySelectorAll(s)];
  const esc = s => s;  // JSON content is authored by us; may contain inline HTML / TeX
  const typeset = el => { let t = 0; const iv = setInterval(() => { if (window.MathJax && MathJax.typesetPromise) { clearInterval(iv); MathJax.typesetPromise([el]).catch(() => {}); } else if (++t > 80) clearInterval(iv); }, 250); };

  fetch("assets/lecture_m" + m + ".json?v=" + Date.now()).then(r => r.json()).then(build).catch(() => {});

  function build(D) {
    const hero = $(".hero"); if (!hero) return;
    const secs = $$("main.content > section.section").filter(s => s.id && s.id !== "quiz" || s.id === "quiz");
    const narr = $$(".voicebar audio");

    /* ---------- chapter intro ---------- */
    const box = document.createElement("section");
    box.className = "lec-intro";
    box.innerHTML =
      (D.oneIdea ? '<div class="lec-one"><div class="lec-k">THE ONE IDEA OF THIS CHAPTER</div><p>' + D.oneIdea + '</p></div>' : '') +
      '<div class="lec-goals"><div class="lec-k">AFTER THIS CHAPTER YOU CAN</div><ul>' + (D.goals || []).map(g => "<li>" + g + "</li>").join("") + '</ul></div>' +
      '<div class="lec-bar"><button class="btn lec-narr" ' + (narr.length ? "" : "hidden") + '>▶ Narrate this chapter</button>' +
      '<span class="lec-meta">' + narr.length + ' narrated parts' + (D.minutes ? ' · about ' + D.minutes + ' min' : '') + '</span>' +
      '<span class="lec-layout">Layout <button data-l="text">Text</button><button data-l="even" class="on">Even</button><button data-l="diagram">Diagram</button></span></div>' +
      '<div class="lec-now" hidden></div>';
    hero.insertAdjacentElement("afterend", box);
    typeset(box);

    /* ---------- layout switch ---------- */
    const setLayout = l => {
      document.body.classList.remove("lay-text", "lay-diagram"); if (l !== "even") document.body.classList.add("lay-" + l);
      $$(".lec-layout button", box).forEach(b => b.classList.toggle("on", b.dataset.l === l));
      try { localStorage.setItem("bcs311-layout", l); } catch (e) {}
    };
    $$(".lec-layout button", box).forEach(b => b.addEventListener("click", () => setLayout(b.dataset.l)));
    try { const l = localStorage.getItem("bcs311-layout"); if (l) setLayout(l); } catch (e) {}

    /* ---------- per-section material ---------- */
    let qn = 0;
    secs.forEach(sec => {
      const S = (D.sections || {})[sec.id]; if (!S) return;
      const panel = $(".panel", sec) || sec, h2 = $("h2", panel);
      if (S.idea && h2) {
        const d = document.createElement("div"); d.className = "lec-idea";
        d.innerHTML = '<span class="lec-tag">Big picture</span> ' + S.idea;
        const after = $(".voicebar", panel) || h2; after.insertAdjacentElement("afterend", d); typeset(d);
      }
      const tail = document.createElement("div"); tail.className = "lec-tail";
      if (S.mistake) tail.innerHTML += '<div class="lec-mistake"><div class="lec-k">COMMON MISTAKE</div>' + S.mistake + '</div>';
      if (S.check) {
        const c = S.check, id = "q" + (++qn);
        tail.innerHTML += '<div class="lec-check"><div class="lec-k">CHECK YOUR UNDERSTANDING</div><p class="lec-q">' + c.q + '</p>' +
          c.opts.map((o, i) => '<button class="lec-opt" data-i="' + i + '">' + o + '</button>').join("") +
          '<div class="lec-why" hidden></div></div>';
        tail.addEventListener("click", e => {
          const b = e.target.closest(".lec-opt"); if (!b) return;
          const ok = +b.dataset.i === c.ans, chk = b.closest(".lec-check");
          $$(".lec-opt", chk).forEach(x => x.classList.remove("right", "wrong"));
          b.classList.add(ok ? "right" : "wrong");
          const w = $(".lec-why", chk); w.hidden = false;
          w.innerHTML = (ok ? "✅ Correct. " : "❌ Not quite. ") + (c.why || ""); typeset(w);
        });
      }
      if (tail.innerHTML) { panel.appendChild(tail); typeset(tail); }
    });

    /* ---------- follow-along: keep the page focused on what is being said ---------- */
    const TOP = 90;                                   // sticky top bar clearance
    let userScrollAt = 0;
    ["wheel", "touchmove", "keydown"].forEach(ev => window.addEventListener(ev, () => { userScrollAt = Date.now(); }, { passive: true }));
    const weight = el => {
      if (el.matches("img,figure,.figout,.mini-svg") || el.querySelector("img,canvas,.mini-svg")) return 160;
      if (el.matches(".demo") || el.querySelector(".demo")) return 120;
      if (el.matches(".codewrap,pre,.tablewrap,table,.eqn")) return 70;
      return Math.max(40, (el.textContent || "").trim().length);
    };
    function blocksFor(a) {
      const bar = a.closest(".voicebar"), panel = bar.closest(".panel") || bar.parentElement, out = [];
      let top = bar; while (top.parentElement && top.parentElement !== panel) top = top.parentElement;
      for (let e = top.nextElementSibling; e; e = e.nextElementSibling) {
        if (e.matches(".voicebar,.lec-tail,.practice,.gen-practice") || e.querySelector(".voicebar")) break;
        if (e.hidden || getComputedStyle(e).display === "none") continue;
        out.push(e);
      }
      return out.length ? out : [bar];
    }
    const clearFocus = () => $$(".lec-focus").forEach(x => x.classList.remove("lec-focus"));
    narr.forEach(a => {
      let blocks = [], cum = [], cur = -1;
      a.addEventListener("play", () => {
        blocks = blocksFor(a); let t = 0; cum = blocks.map(b => (t += weight(b))); cum = cum.map(v => v / t); cur = -1;
      });
      a.addEventListener("timeupdate", () => {
        if (!a.duration || !blocks.length) return;
        const p = a.currentTime / a.duration;
        let i = cum.findIndex(v => p <= v); if (i < 0) i = blocks.length - 1;
        if (i === cur) return; cur = i;
        clearFocus(); blocks[i].classList.add("lec-focus");
        if (Date.now() - userScrollAt < 4000) return;            // the reader is scrolling: do not fight them
        const r = blocks[i].getBoundingClientRect();
        if (r.top < TOP || r.top > innerHeight * .55 || r.bottom > innerHeight) scrollTo({ top: scrollY + r.top - TOP - 20, behavior: "smooth" });
      });
      a.addEventListener("pause", clearFocus); a.addEventListener("ended", clearFocus);
    });

    /* ---------- narrate the whole chapter ---------- */
    const btn = $(".lec-narr", box), now = $(".lec-now", box);
    let idx = -1, running = false;
    const stop = () => { running = false; narr.forEach(a => { a.pause(); }); btn.textContent = "▶ Narrate this chapter"; now.hidden = true; $$(".lec-playing").forEach(x => x.classList.remove("lec-playing")); };
    function play(i) {
      if (i >= narr.length) { stop(); return; }
      idx = i; const a = narr[i], sec = a.closest("section");
      narr.forEach(x => { if (x !== a) x.pause(); });
      $$(".lec-playing").forEach(x => x.classList.remove("lec-playing")); sec.classList.add("lec-playing");
      userScrollAt = 0; const bar = a.closest(".voicebar"); (bar || sec).scrollIntoView({ behavior: "smooth", block: "start" });
      now.hidden = false; now.textContent = "Playing part " + (i + 1) + " of " + narr.length + " · " + ((sec.querySelector("h2") || {}).textContent || "");
      a.currentTime = 0; a.play().catch(() => {});
    }
    narr.forEach((a, i) => a.addEventListener("ended", () => { if (running) play(i + 1); }));
    btn.addEventListener("click", () => { if (running) { stop(); } else { running = true; btn.textContent = "⏸ Stop narration"; play(0); } });
  }
})();
