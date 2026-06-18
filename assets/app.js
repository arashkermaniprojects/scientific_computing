/* ===========================================================
   BCS 311 — Scientific Computing | shared app.js
   Progress tracking, quiz engine, code copy, navigation
   =========================================================== */
(function () {
  "use strict";

  /* ---------- Module registry ---------- */
  const MODULES = [
    { id: "m1", file: "module1.html", clo: "CLO-1", title: "Computing with MATLAB" },
    { id: "m2", file: "module2.html", clo: "CLO-2", title: "Linear & Non-linear Systems" },
    { id: "m3", file: "module3.html", clo: "CLO-3", title: "Least-Squares Curve Fitting" },
    { id: "m4", file: "module4.html", clo: "CLO-4", title: "Interpolation & Approximation" },
    { id: "m5", file: "module5.html", clo: "CLO-5", title: "Integration & ODEs" }
  ];
  const STORE_KEY = "bcs311_progress_v1";

  /* ---------- Progress storage (localStorage, graceful fallback) ---------- */
  function loadProgress() {
    try { return JSON.parse(localStorage.getItem(STORE_KEY)) || {}; }
    catch (e) { return {}; }
  }
  function saveProgress(p) {
    try { localStorage.setItem(STORE_KEY, JSON.stringify(p)); } catch (e) {}
  }
  function isDone(id) { return !!loadProgress()[id]; }
  function setDone(id, v) {
    const p = loadProgress();
    if (v) p[id] = Date.now(); else delete p[id];
    saveProgress(p);
  }
  function completedCount() {
    const p = loadProgress();
    return MODULES.filter(m => p[m.id]).length;
  }

  window.BCS311 = { MODULES, loadProgress, isDone, setDone, completedCount };

  /* ---------- Code copy buttons ---------- */
  function wireCopyButtons() {
    document.querySelectorAll(".copybtn").forEach(btn => {
      btn.addEventListener("click", () => {
        const wrap = btn.closest(".codewrap");
        const src = wrap.querySelector("code, textarea");
        const text = src.value !== undefined ? src.value : src.innerText;
        navigator.clipboard.writeText(text).then(() => {
          const old = btn.textContent;
          btn.textContent = "Copied!";
          btn.classList.add("copied");
          setTimeout(() => { btn.textContent = old; btn.classList.remove("copied"); }, 1400);
        });
      });
    });
  }

  /* ---------- Quiz engine ----------
     Usage: renderQuiz("quizContainerId", "moduleId", [ {q, options:[...], answer:Index, explain} ]) */
  window.renderQuiz = function (containerId, moduleId, questions) {
    const root = document.getElementById(containerId);
    if (!root) return;
    let answered = 0, correct = 0;
    const total = questions.length;

    questions.forEach((item, qi) => {
      const div = document.createElement("div");
      div.className = "qitem";
      const qEl = document.createElement("div");
      qEl.className = "qq";
      qEl.innerHTML = '<span class="n">Q' + (qi + 1) + '.</span>' + item.q;
      div.appendChild(qEl);

      const explain = document.createElement("div");
      explain.className = "explain";
      explain.innerHTML = "<b>Why:</b> " + item.explain;

      item.options.forEach((opt, oi) => {
        const label = document.createElement("label");
        label.className = "opt";
        label.innerHTML = '<input type="radio" name="' + containerId + "_" + qi + '">' + opt;
        label.querySelector("input").addEventListener("change", function () {
          if (div.dataset.locked) return;
          div.dataset.locked = "1";
          answered++;
          div.querySelectorAll(".opt").forEach((o, k) => {
            o.classList.add("disabled");
            if (k === item.answer) o.classList.add("correct");
          });
          if (oi === item.answer) { label.classList.add("correct"); correct++; }
          else { label.classList.add("wrong"); }
          explain.classList.add("show");
          updateScore();
        });
        div.appendChild(label);
      });
      div.appendChild(explain);
      root.appendChild(div);
    });

    const actions = document.createElement("div");
    actions.className = "quiz-actions";
    const scoreEl = document.createElement("span");
    scoreEl.className = "score";
    scoreEl.textContent = "Answered 0 / " + total;
    const reset = document.createElement("button");
    reset.className = "btn ghost";
    reset.textContent = "Reset quiz";
    reset.addEventListener("click", () => {
      root.innerHTML = "";
      renderQuiz(containerId, moduleId, questions);
    });
    actions.appendChild(scoreEl);
    actions.appendChild(reset);
    root.appendChild(actions);

    function updateScore() {
      scoreEl.textContent = "Score: " + correct + " / " + total +
        (answered === total ? " — complete!" : "  (answered " + answered + "/" + total + ")");
      scoreEl.className = "score " + (correct / total >= 0.8 ? "good" : correct / total >= 0.5 ? "mid" : "low");
    }
  };

  /* ---------- Mark-complete bar (per module page) ---------- */
  window.wireCompleteBar = function (moduleId) {
    const bar = document.getElementById("completeBar");
    if (!bar) return;
    const btn = bar.querySelector("button");
    const txt = bar.querySelector(".txt");
    function refresh() {
      const done = isDone(moduleId);
      bar.classList.toggle("done", done);
      txt.textContent = done ? "✓ Module marked complete — great work!" : "Finished this module? Mark it complete to track your progress.";
      btn.textContent = done ? "Mark as not complete" : "Mark module complete";
      btn.className = done ? "btn ghost" : "btn";
    }
    btn.addEventListener("click", () => { setDone(moduleId, !isDone(moduleId)); refresh(); buildSidebarProgress(); });
    refresh();
  };

  /* ---------- Sidebar progress widget ---------- */
  function buildSidebarProgress() {
    const host = document.getElementById("progressWidget");
    if (!host) return;
    const done = completedCount(), total = MODULES.length;
    const pct = Math.round((done / total) * 100);
    host.innerHTML =
      '<h4>Your progress</h4>' +
      '<div class="bar"><span style="width:' + pct + '%"></span></div>' +
      '<div class="pct">' + done + ' of ' + total + ' modules complete (' + pct + '%)</div>';
  }

  /* ---------- Active-link highlighting in TOC via scrollspy ---------- */
  function wireScrollSpy() {
    const links = Array.from(document.querySelectorAll(".toc a[href^='#']"));
    if (!links.length) return;
    const map = {};
    links.forEach(a => { const t = document.getElementById(a.getAttribute("href").slice(1)); if (t) map[a.getAttribute("href").slice(1)] = a; });
    const obs = new IntersectionObserver(entries => {
      entries.forEach(e => {
        if (e.isIntersecting) {
          links.forEach(l => l.classList.remove("active"));
          if (map[e.target.id]) map[e.target.id].classList.add("active");
        }
      });
    }, { rootMargin: "-80px 0px -70% 0px", threshold: 0 });
    Object.keys(map).forEach(id => { const el = document.getElementById(id); if (el) obs.observe(el); });
  }

  /* ---------- Mobile menus ---------- */
  function wireMobileMenus() {
    const navToggle = document.querySelector(".menu-toggle");
    if (navToggle) navToggle.addEventListener("click", () => {
      document.querySelector(".topbar nav").classList.toggle("open");
    });
  }

  /* ---------- Home page module-card completion ticks ---------- */
  function markHomeCards() {
    document.querySelectorAll(".modcard[data-module]").forEach(card => {
      if (isDone(card.dataset.module)) card.classList.add("done");
    });
  }

  document.addEventListener("DOMContentLoaded", () => {
    wireCopyButtons();
    buildSidebarProgress();
    wireScrollSpy();
    wireMobileMenus();
    markHomeCards();
  });
})();
