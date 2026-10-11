/* Symbolic (number-free) step-by-step animations for Module 2:
   Gaussian elimination, LU factorization, Cholesky factorization.
   Usage: <div class="demo sd-demo" data-sd="gauss|lu|chol"></div>
   Each algorithm is simulated once into a list of snapshots (one per step);
   the renderer just diffs the DOM against the chosen snapshot, so Back/Scrub are free. */
(function(){
'use strict';

/* ---------- symbol helpers ---------- */
// S('a',2,1,0) -> a_21 with stage superscript (0 -> none).  idx parts are concatenated.
function S(name,i,j,s){
  const idx=(j==null||j==='')?String(i):''+i+j;
  const sub=`<sub>${idx}</sub>`;
  if(s==null||s===0) return `<i>${name}</i>${sub}`;
  return `<i>${name}</i><span class="ss"><sup>(${s})</sup>${sub}</span>`;
}
const grid=(r,c,f)=>Array.from({length:r},(_,i)=>Array.from({length:c},(_,j)=>f(i,j)));
const range=(a,b)=>{const o=[];for(let i=a;i<=b;i++)o.push(i);return o;};

/* ---------- generic snapshot builder ---------- */
function Builder(meta){
  const steps=[]; let hl={}, lab={};
  return {
    meta, steps,
    clear(){hl={};lab={};},
    set(id,r,c,k){hl[id+':'+r+','+c]=k;},
    setLab(id,r,k){lab[id+':'+r]=k;},
    get(id,r,c){return hl[id+':'+r+','+c]||'';},
    getLab(id,r){return lab[id+':'+r]||'';},
    emit(phase,cap,mats){steps.push({phase,cap,mats});}
  };
}
const cap=(t,w,g,i)=>({t,w,g,i:i||''});

/* =========================================================
   GAUSSIAN ELIMINATION   (0-based internals, 1-based labels)
   ========================================================= */
function buildGauss(n){
  const aug=n+1;
  const B=Builder({
    mats:[{id:'M',title:'Augmented matrix [ <i>A</i> | <b>b</b> ]',rows:n,cols:aug,bar:n-1,rowLabels:true},
          {id:'X',title:'Unknowns <b>x</b>',rows:n,cols:1}],
    layout:[['M','X']],
    phases:['Forward elimination','Back substitution'],
    legend:[['k-pivot','pivot'],['k-src','entries used'],['k-tgt','entry / row being changed'],['k-done','new zero · solved']]
  });
  const W=grid(n,aug,()=>({s:0,z:false}));
  const xd=Array(n).fill(false);
  const lbl=(i,j)=>{const c=W[i][j]; if(c.z) return '0'; return j===n?S('b',i+1,null,c.s):S('a',i+1,j+1,c.s);};
  const emit=(ph,c)=>B.emit(ph,c,{
    M:{cells:grid(n,aug,(i,j)=>({t:lbl(i,j),c:B.get('M',i,j)||(W[i][j].z?'k-zero':'')})),
       lab:range(0,n-1).map(i=>B.getLab('M',i))},
    X:{cells:grid(n,1,(i)=>({t:S('x',i+1),c:B.get('X',i,0)||(xd[i]?'k-done':'k-pend')}))}
  });

  B.clear();
  emit(0,cap('Start',
    'Every entry is a symbol: <i>a<sub>ij</sub></i> from <i>A</i>, <i>b<sub>i</sub></i> from <b>b</b>. The goal is to turn everything below the diagonal into zeros, then read off the unknowns from the bottom up.',
    'forward elimination → upper-triangular → back substitution',
    `size <i>n</i> = ${n}  ·  press <b>Step ▸</b> or <b>Auto ▶</b>`));

  for(let k=1;k<=n-1;k++){
    // --- pivot selection
    B.clear();
    B.set('M',k-1,k-1,'k-pivot'); B.setLab('M',k-1,'k-pivot');
    for(let j=k;j<aug;j++) B.set('M',k-1,j,'k-src');
    for(let i=k+1;i<=n;i++){ B.set('M',i-1,k-1,'k-soft'); B.setLab('M',i-1,'k-soft'); }
    emit(0,cap(`Stage <i>k</i> = ${k}: take the pivot`,
      `The pivot ${S('a',k,k,k-1)} sits on the diagonal. Every entry below it in column ${k} has to be cleared, one row at a time, using row ${k}.`,
      `for <i>i</i> = <i>k</i>+1 … <i>n</i>:  ${S('m','i','k')} = ${S('a','i','k','k−1')} / ${S('a','k','k','k−1')},   <i>R<sub>i</sub></i> ← <i>R<sub>i</sub></i> − ${S('m','i','k')}·<i>R<sub>k</sub></i>`,
      `<i>k</i> = ${k}  →  rows ${k+1}${k+1<n?'…'+n:''} still to clear`));

    for(let i=k+1;i<=n;i++){
      // --- multiplier
      B.clear();
      B.set('M',k-1,k-1,'k-pivot'); B.set('M',i-1,k-1,'k-tgt');
      B.setLab('M',k-1,'k-pivot'); B.setLab('M',i-1,'k-tgt');
      emit(0,cap(`Multiplier for row ${i}`,
        `Ratio of the entry to remove over the pivot — exactly the factor that will cancel it.`,
        `${S('m','i','k')} = ${S('a','i','k','k−1')} / ${S('a','k','k','k−1')}`,
        `${S('m',i,k)} = ${S('a',i,k,k-1)} / ${S('a',k,k,k-1)}`));
      // --- row update
      for(let j=k-1;j<aug;j++) W[i-1][j].s=k;
      W[i-1][k-1].z=true;
      B.clear();
      B.set('M',k-1,k-1,'k-pivot'); B.setLab('M',k-1,'k-pivot'); B.setLab('M',i-1,'k-tgt');
      for(let j=k;j<aug;j++){ B.set('M',k-1,j,'k-src'); B.set('M',i-1,j,'k-tgt'); }
      B.set('M',i-1,k-1,'k-done');
      const j1=k+1;
      emit(0,cap(`Row ${i} ← row ${i} − ${S('m',i,k)} · row ${k}`,
        `Subtract the multiple of the pivot row. The entry under the pivot becomes <b>0</b>; the rest of row ${i} (and its right-hand side) change to stage (${k}).`,
        `${S('a','i','j','k')} = ${S('a','i','j','k−1')} − ${S('m','i','k')}·${S('a','k','j','k−1')}   (<i>j</i> = <i>k</i>+1 … <i>n</i>, and the same for <i>b</i>)`,
        j1<=n?`e.g. ${S('a',i,j1,k)} = ${S('a',i,j1,k-1)} − ${S('m',i,k)}·${S('a',k,j1,k-1)}   …`:`${S('b',i,null,k)} = ${S('b',i,null,k-1)} − ${S('m',i,k)}·${S('b',k,null,k-1)}`));
    }
  }
  // --- triangular
  B.clear();
  for(let i=0;i<n;i++){ B.set('M',i,i,'k-pivot'); }
  emit(0,cap('Upper-triangular',
    `Everything below the diagonal is <b>0</b>. The system now reads row by row from the bottom: the last row has one unknown, the one above has two, and so on.`,
    `<i>U x</i> = <b>c</b>,   with ${S('a','i','j','i−1')} on and above the diagonal`,
    `forward elimination cost grows like <i>n</i>³ — back substitution only like <i>n</i>²`));

  // --- back substitution
  for(let i=n;i>=1;i--){
    B.clear();
    B.set('M',i-1,i-1,'k-pivot'); B.setLab('M',i-1,'k-tgt');
    for(let j=i+1;j<=n;j++){ B.set('M',i-1,j-1,'k-src'); B.set('X',j-1,0,'k-src'); }
    B.set('M',i-1,n,'k-src');
    B.set('X',i-1,0,'k-tgt'); xd[i-1]=true;
    const terms=range(i+1,n).map(j=>`${S('a',i,j,i-1)}·${S('x',j)}`);
    emit(1,cap(i===n?`Solve the last row for ${S('x',n)}`:`Back-substitute for ${S('x',i)}`,
      i===n?`The last row has a single unknown, so divide the right-hand side by the pivot.`
           :`Row ${i} has one new unknown, ${S('x',i)}; everything to its right is already known (blue), so move it across and divide by the pivot.`,
      `${S('x','i')} = ( ${S('b','i',null,'i−1')} − Σ<sub>j=i+1…n</sub> ${S('a','i','j','i−1')}·${S('x','j')} ) / ${S('a','i','i','i−1')}`,
      `${S('x',i)} = ${terms.length?`( ${S('b',i,null,i-1)} − ${terms.join(' − ')} )`:S('b',i,null,i-1)} / ${S('a',i,i,i-1)}`));
  }
  B.clear();
  emit(1,cap('Solved',
    `All unknowns are found, bottom to top. Notice the pattern: <b>every step needed only quantities already computed</b> — that is why the triangular shape matters.`,
    `${S('x',1)}, …, ${S('x','n')} determined`,
    `Next: LU factorization keeps the multipliers ${S('m','i','k')} so this work can be reused.`));
  return B;
}

/* =========================================================
   LU FACTORIZATION
   ========================================================= */
function buildLU(n){
  const B=Builder({
    mats:[{id:'W',title:'<i>A</i> → <i>U</i>  (working matrix)',rows:n,cols:n,rowLabels:true},
          {id:'L',title:'<i>L</i>  (multipliers collected)',rows:n,cols:n},
          {id:'b',title:'<b>b</b>',rows:n,cols:1},
          {id:'y',title:'<b>y</b>  ( <i>L y</i> = <b>b</b> )',rows:n,cols:1},
          {id:'x',title:'<b>x</b>  ( <i>U x</i> = <b>y</b> )',rows:n,cols:1}],
    layout:[['W','L'],['b','y','x']],
    phases:['Factor  A = L U','Forward solve  L y = b','Backward solve  U x = y'],
    legend:[['k-pivot','pivot'],['k-src','entries used'],['k-tgt','entry / row being changed'],['k-done','finished']]
  });
  const W=grid(n,n,()=>({s:0,z:false,f:false}));
  const Ld=grid(n,n,()=>false);
  const yd=Array(n).fill(false), xd=Array(n).fill(false);
  let live=0; // 0 = factoring (vectors greyed), 1 = solving
  const wl=(i,j)=>{const c=W[i][j]; if(c.z) return '0'; return c.f?S('u',i+1,j+1):S('a',i+1,j+1,c.s);};
  const wc=(i,j)=>B.get('W',i,j)||(W[i][j].z?'k-zero':(W[i][j].f?'k-done':''));
  const ll=(i,j)=>i===j?'1':(i<j?'0':(Ld[i][j]?S('l',i+1,j+1):'·'));
  const lc=(i,j)=>B.get('L',i,j)||(i===j||i<j?'k-dim':(Ld[i][j]?'k-done':'k-pend'));
  const vec=(id,sym,done)=>grid(n,1,(i)=>({t:S(sym,i+1),
     c:B.get(id,i,0)||(!live&&id!=='b'?'k-pend k-dim':(id==='b'?(live?'':'k-dim'):(done[i]?'k-done':'k-pend')))}));
  const emit=(ph,c)=>B.emit(ph,c,{
    W:{cells:grid(n,n,(i,j)=>({t:wl(i,j),c:wc(i,j)})),lab:range(0,n-1).map(i=>B.getLab('W',i))},
    L:{cells:grid(n,n,(i,j)=>({t:ll(i,j),c:lc(i,j)}))},
    b:{cells:vec('b','b',null)},
    y:{cells:vec('y','y',yd)},
    x:{cells:vec('x','x',xd)}
  });
  const finRow=k=>{for(let j=k-1;j<n;j++) W[k-1][j].f=true;};

  B.clear();
  emit(0,cap('Start',
    `Same elimination as before, but now we <b>keep</b> each multiplier: it goes into <i>L</i>. The working matrix will end up as <i>U</i>. The right-hand side <b>b</b> is not touched at all.`,
    `<i>A</i> = <i>L U</i>,   <i>L</i> unit lower-triangular,  <i>U</i> upper-triangular`,
    `size <i>n</i> = ${n}  ·  press <b>Step ▸</b> or <b>Auto ▶</b>`));

  for(let k=1;k<=n-1;k++){
    finRow(k);
    B.clear();
    B.set('W',k-1,k-1,'k-pivot'); B.setLab('W',k-1,'k-pivot');
    for(let i=k+1;i<=n;i++){ B.set('W',i-1,k-1,'k-soft'); B.setLab('W',i-1,'k-soft'); }
    emit(0,cap(`Stage <i>k</i> = ${k}: row ${k} is final`,
      `Row ${k} will not change any more, so it becomes row ${k} of <i>U</i> (entries renamed ${S('u',k,'j')}). Its pivot ${S('u',k,k)} is used to clear column ${k} below.`,
      `${S('u','k','j')} = ${S('a','k','j','k−1')}   (<i>j</i> ≥ <i>k</i>)`,
      `<i>k</i> = ${k}`));
    for(let i=k+1;i<=n;i++){
      Ld[i-1][k-1]=true;
      B.clear();
      B.set('W',k-1,k-1,'k-pivot'); B.set('W',i-1,k-1,'k-src');
      B.setLab('W',k-1,'k-pivot'); B.setLab('W',i-1,'k-tgt');
      B.set('L',i-1,k-1,'k-tgt');
      emit(0,cap(`Multiplier ${S('l',i,k)} → store it in <i>L</i>`,
        `The ratio that cancels the entry is computed and written straight into position (${i},${k}) of <i>L</i> — this is the bookkeeping that makes the factorization reusable.`,
        `${S('l','i','k')} = ${S('a','i','k','k−1')} / ${S('u','k','k')}`,
        `${S('l',i,k)} = ${S('a',i,k,k-1)} / ${S('u',k,k)}`));
      for(let j=k-1;j<n;j++) W[i-1][j].s=k;
      W[i-1][k-1].z=true;
      B.clear();
      B.set('W',k-1,k-1,'k-pivot'); B.setLab('W',k-1,'k-pivot'); B.setLab('W',i-1,'k-tgt');
      for(let j=k;j<n;j++){ B.set('W',k-1,j,'k-src'); B.set('W',i-1,j,'k-tgt'); }
      B.set('W',i-1,k-1,'k-done'); B.set('L',i-1,k-1,'k-src');
      const j1=k+1;
      emit(0,cap(`Row ${i} ← row ${i} − ${S('l',i,k)} · row ${k}`,
        `Apply the stored multiplier to the rest of row ${i}. The entry under the pivot vanishes (it is now remembered by ${S('l',i,k)}).`,
        `${S('a','i','j','k')} = ${S('a','i','j','k−1')} − ${S('l','i','k')}·${S('u','k','j')}   (<i>j</i> = <i>k</i>+1 … <i>n</i>)`,
        j1<=n?`e.g. ${S('a',i,j1,k)} = ${S('a',i,j1,k-1)} − ${S('l',i,k)}·${S('u',k,j1)}   …`:''));
    }
  }
  finRow(n);
  B.clear();
  for(let i=0;i<n;i++){B.set('W',i,i,'k-pivot');}
  emit(0,cap('Factorization complete: A = L U',
    `The working matrix is now <i>U</i> and the multipliers fill <i>L</i> below a diagonal of ones. The expensive part is <b>done once</b>.`,
    `<i>A</i> = <i>L U</i>`,
    `cost ≈ <i>n</i>³ operations — never repeated for new right-hand sides`));

  // ---- forward solve
  live=1;
  for(let i=1;i<=n;i++){
    B.clear();
    B.set('L',i-1,i-1,'k-pivot'); B.set('b',i-1,0,'k-src'); B.set('y',i-1,0,'k-tgt');
    for(let j=1;j<i;j++){ B.set('L',i-1,j-1,'k-src'); B.set('y',j-1,0,'k-src'); }
    yd[i-1]=true;
    const terms=range(1,i-1).map(j=>`${S('l',i,j)}·${S('y',j)}`);
    emit(1,cap(i===1?`Forward solve: ${S('y',1)}`:`Forward solve: ${S('y',i)}`,
      i===1?`<i>L</i> has ones on the diagonal, so the first unknown is just the first entry of <b>b</b>.`
           :`Row ${i} of <i>L</i> multiplies the already-known ${S('y',1)}…${S('y',i-1)} (blue). Move those terms across; the diagonal entry is 1, so no division is needed.`,
      `${S('y','i')} = ${S('b','i')} − Σ<sub>j&lt;i</sub> ${S('l','i','j')}·${S('y','j')}`,
      `${S('y',i)} = ${S('b',i)}${terms.length?' − '+terms.join(' − '):''}`));
  }
  // ---- backward solve
  for(let i=n;i>=1;i--){
    B.clear();
    B.set('W',i-1,i-1,'k-pivot'); B.setLab('W',i-1,'k-tgt'); B.set('y',i-1,0,'k-src'); B.set('x',i-1,0,'k-tgt');
    for(let j=i+1;j<=n;j++){ B.set('W',i-1,j-1,'k-src'); B.set('x',j-1,0,'k-src'); }
    xd[i-1]=true;
    const terms=range(i+1,n).map(j=>`${S('u',i,j)}·${S('x',j)}`);
    emit(2,cap(`Backward solve: ${S('x',i)}`,
      i===n?`The last row of <i>U</i> has one unknown: divide ${S('y',n)} by ${S('u',n,n)}.`
           :`Subtract the known terms to the right of the diagonal (blue) from ${S('y',i)}, then divide by the pivot ${S('u',i,i)}.`,
      `${S('x','i')} = ( ${S('y','i')} − Σ<sub>j&gt;i</sub> ${S('u','i','j')}·${S('x','j')} ) / ${S('u','i','i')}`,
      `${S('x',i)} = ${terms.length?`( ${S('y',i)} − ${terms.join(' − ')} )`:S('y',i)} / ${S('u',i,i)}`));
  }
  B.clear();
  emit(2,cap('Solved — and <i>L</i>, <i>U</i> are reusable',
    `For a new right-hand side <b>b′</b>, skip the whole factorization: reuse <i>L</i> and <i>U</i> and repeat only the two triangular solves.`,
    `factor once: ~ <i>n</i>³   ·   each new <b>b</b>: ~ <i>n</i>²`,
    ''));
  return B;
}

/* =========================================================
   CHOLESKY FACTORIZATION
   ========================================================= */
function buildChol(n){
  const B=Builder({
    mats:[{id:'A',title:'<i>A</i>  (symmetric, positive definite)',rows:n,cols:n},
          {id:'L',title:'<i>L</i>',rows:n,cols:n},
          {id:'T',title:'<i>L</i><sup>T</sup>',rows:n,cols:n},
          {id:'b',title:'<b>b</b>',rows:n,cols:1},
          {id:'y',title:'<b>y</b>  ( <i>L y</i> = <b>b</b> )',rows:n,cols:1},
          {id:'x',title:'<b>x</b>  ( <i>L</i><sup>T</sup><b>x</b> = <b>y</b> )',rows:n,cols:1}],
    layout:[['A','=','L','·','T'],['b','y','x']],
    phases:['Factor  A = L Lᵀ','Forward solve  L y = b','Backward solve  Lᵀ x = y'],
    legend:[['k-pivot','divisor'],['k-src','entries used'],['k-tgt','entry being computed'],['k-done','finished']]
  });
  const Ld=grid(n,n,()=>false);
  const yd=Array(n).fill(false), xd=Array(n).fill(false);
  let live=0;
  const ac=(i,j)=>B.get('A',i,j)||(i<j?'k-dim':'');
  const lc=(i,j)=>B.get('L',i,j)||(i<j?'k-dim':(Ld[i][j]?'k-done':'k-pend'));
  const ll=(i,j)=>i<j?'0':(Ld[i][j]?S('l',i+1,j+1):'·');
  const tc=(i,j)=>B.get('T',i,j)||(i>j?'k-dim':(Ld[j][i]?'k-done':'k-pend'));
  const tl=(i,j)=>i>j?'0':(Ld[j][i]?S('l',j+1,i+1):'·');
  const vec=(id,sym,done)=>grid(n,1,(i)=>({t:S(sym,i+1),
     c:B.get(id,i,0)||(id==='b'?(live?'':'k-dim'):(!live?'k-pend k-dim':(done[i]?'k-done':'k-pend')))}));
  const emit=(ph,c)=>B.emit(ph,c,{
    A:{cells:grid(n,n,(i,j)=>({t:S('a',i+1,j+1),c:ac(i,j)}))},
    L:{cells:grid(n,n,(i,j)=>({t:ll(i,j),c:lc(i,j)}))},
    T:{cells:grid(n,n,(i,j)=>({t:tl(i,j),c:tc(i,j)}))},
    b:{cells:vec('b','b',null)}, y:{cells:vec('y','y',yd)}, x:{cells:vec('x','x',xd)}
  });

  B.clear();
  emit(0,cap('Start',
    `Because <i>A</i> is symmetric, only its lower triangle (and diagonal) is ever read — the faded upper part is a mirror image. We build <i>L</i> one column at a time, and <i>L</i><sup>T</sup> is just its mirror.`,
    `<i>A</i> = <i>L L</i><sup>T</sup>,   <i>L</i> lower-triangular with positive diagonal`,
    `size <i>n</i> = ${n}  ·  press <b>Step ▸</b> or <b>Auto ▶</b>`));

  for(let j=1;j<=n;j++){
    // diagonal
    B.clear();
    B.set('A',j-1,j-1,'k-src'); B.set('L',j-1,j-1,'k-tgt'); B.set('T',j-1,j-1,'k-tgt');
    for(let k=1;k<j;k++){ B.set('L',j-1,k-1,'k-src'); B.set('T',k-1,j-1,'k-src'); }
    Ld[j-1][j-1]=true;
    const sq=range(1,j-1).map(k=>`${S('l',j,k)}²`);
    emit(0,cap(`Column ${j}: diagonal ${S('l',j,j)}`,
      j===1?`The first diagonal entry is simply the square root of ${S('a',1,1)}.`
           :`Subtract the squares of the entries already placed to the left in row ${j} (blue) from ${S('a',j,j)}, then take the square root. The radicand <b>must be positive</b> — that is the positive-definite test.`,
      `${S('l','j','j')} = √( ${S('a','j','j')} − Σ<sub>k&lt;j</sub> ${S('l','j','k')}² )`,
      `${S('l',j,j)} = √( ${S('a',j,j)}${sq.length?' − '+sq.join(' − '):''} )`));
    // off-diagonal in this column
    for(let i=j+1;i<=n;i++){
      B.clear();
      B.set('A',i-1,j-1,'k-src'); B.set('L',j-1,j-1,'k-pivot'); B.set('T',j-1,j-1,'k-pivot');
      B.set('L',i-1,j-1,'k-tgt'); B.set('T',j-1,i-1,'k-tgt');
      for(let k=1;k<j;k++){ B.set('L',i-1,k-1,'k-src'); B.set('L',j-1,k-1,'k-src'); }
      Ld[i-1][j-1]=true;
      const prods=range(1,j-1).map(k=>`${S('l',i,k)}·${S('l',j,k)}`);
      emit(0,cap(`Column ${j}: entry ${S('l',i,j)}`,
        j===1?`First column: divide ${S('a',i,1)} by the diagonal ${S('l',1,1)} (teal).`
             :`Rows ${i} and ${j} of <i>L</i> overlap to the left of column ${j} (blue pairs). Subtract their dot product from ${S('a',i,j)}, then divide by the diagonal ${S('l',j,j)} (teal).`,
        `${S('l','i','j')} = ( ${S('a','i','j')} − Σ<sub>k&lt;j</sub> ${S('l','i','k')}·${S('l','j','k')} ) / ${S('l','j','j')}`,
        `${S('l',i,j)} = ( ${S('a',i,j)}${prods.length?' − '+prods.join(' − '):''} ) / ${S('l',j,j)}`));
    }
  }
  B.clear();
  emit(0,cap('Factorization complete: A = L Lᵀ',
    `Only one triangle was computed and stored, and the second factor is free: <i>L</i><sup>T</sup> is the mirror image. That is why Cholesky costs about <b>half</b> of LU.`,
    `<i>A</i> = <i>L L</i><sup>T</sup>`,
    `cost ≈ <i>n</i>³/3 operations (LU: 2<i>n</i>³/3)`));

  live=1;
  for(let i=1;i<=n;i++){
    B.clear();
    B.set('L',i-1,i-1,'k-pivot'); B.set('b',i-1,0,'k-src'); B.set('y',i-1,0,'k-tgt');
    for(let j=1;j<i;j++){ B.set('L',i-1,j-1,'k-src'); B.set('y',j-1,0,'k-src'); }
    yd[i-1]=true;
    const terms=range(1,i-1).map(j=>`${S('l',i,j)}·${S('y',j)}`);
    emit(1,cap(`Forward solve: ${S('y',i)}`,
      i===1?`Only the diagonal of row 1 is present: divide ${S('b',1)} by ${S('l',1,1)}.`
           :`Use the known ${S('y',1)}…${S('y',i-1)} (blue), then divide by the diagonal ${S('l',i,i)} (teal).`,
      `${S('y','i')} = ( ${S('b','i')} − Σ<sub>j&lt;i</sub> ${S('l','i','j')}·${S('y','j')} ) / ${S('l','i','i')}`,
      `${S('y',i)} = ${terms.length?`( ${S('b',i)} − ${terms.join(' − ')} )`:S('b',i)} / ${S('l',i,i)}`));
  }
  for(let i=n;i>=1;i--){
    B.clear();
    B.set('T',i-1,i-1,'k-pivot'); B.set('y',i-1,0,'k-src'); B.set('x',i-1,0,'k-tgt');
    for(let j=i+1;j<=n;j++){ B.set('T',i-1,j-1,'k-src'); B.set('x',j-1,0,'k-src'); }
    xd[i-1]=true;
    const terms=range(i+1,n).map(j=>`${S('l',j,i)}·${S('x',j)}`);
    emit(2,cap(`Backward solve: ${S('x',i)}`,
      i===n?`Last row of <i>L</i><sup>T</sup>: divide ${S('y',n)} by ${S('l',n,n)}.`
           :`Row ${i} of <i>L</i><sup>T</sup> is column ${i} of <i>L</i>. Subtract the known terms (blue) from ${S('y',i)} and divide by ${S('l',i,i)}.`,
      `${S('x','i')} = ( ${S('y','i')} − Σ<sub>j&gt;i</sub> ${S('l','j','i')}·${S('x','j')} ) / ${S('l','i','i')}`,
      `${S('x',i)} = ${terms.length?`( ${S('y',i)} − ${terms.join(' − ')} )`:S('y',i)} / ${S('l',i,i)}`));
  }
  B.clear();
  emit(2,cap('Solved',
    `Two triangular solves with the <i>same</i> stored factor <i>L</i> (once as is, once transposed) finish the job.`,
    `<i>L y</i> = <b>b</b>,   <i>L</i><sup>T</sup><b>x</b> = <b>y</b>`,
    `If a radicand had been ≤ 0, <i>A</i> would not be positive definite — fall back to LU.`));
  return B;
}

const BUILDERS={gauss:buildGauss,lu:buildLU,chol:buildChol};

/* =========================================================
   RENDERER
   ========================================================= */
function mount(root){
  const kind=root.dataset.sd, build=BUILDERS[kind]; if(!build) return;
  let n=+(root.dataset.n||4), data, k=0, timer=null, speed=1100, first=true;

  root.innerHTML=`
    <div class="demo-head"><span class="pill">Live</span> <span class="sd-title"></span></div>
    <div class="sd-phases"></div>
    <div class="sd-stage"></div>
    <div class="sd-legend"></div>
    <div class="sd-prog"><input type="range" class="sd-scrub" min="0" value="0" aria-label="Progress through the algorithm"><span class="sd-count"></span></div>
    <div class="sd-cap" aria-live="polite"><div class="sd-ct"></div><div class="sd-cw"></div><div class="sd-cg"></div><div class="sd-ci"></div></div>
    <div class="btnrow sd-ctl">
      <button class="btn ghost sd-back">◂ Back</button>
      <button class="btn sd-step">Step ▸</button>
      <button class="btn accent sd-auto">Auto ▶</button>
      <button class="btn ghost sd-reset">Reset</button>
      <label class="sd-sel">size <i>n</i> <select class="sd-n"><option>3</option><option>4</option><option>5</option></select></label>
      <label class="sd-sel">speed <select class="sd-sp"><option value="1700">slow</option><option value="1100" selected>normal</option><option value="600">fast</option></select></label>
    </div>`;
  const $=s=>root.querySelector(s);
  const stage=$('.sd-stage'), scrub=$('.sd-scrub'), nSel=$('.sd-n'), autoBtn=$('.sd-auto');
  $('.sd-title').textContent={gauss:'Gaussian elimination — symbolic walkthrough',lu:'LU factorization — symbolic walkthrough',chol:'Cholesky factorization — symbolic walkthrough'}[kind];
  nSel.value=String(n);

  function setup(){
    stop(); data=build(n); k=0; first=true;
    stage.innerHTML=''; stage.style.setProperty('--sd',(n>=5?(kind==='chol'?36:40):46)+'px');
    const byId={}; data.meta.mats.forEach(m=>byId[m.id]=m);
    data.meta.layout.forEach(row=>{
      const r=document.createElement('div'); r.className='sd-row';
      row.forEach(id=>{
        if(!byId[id]){ const g=document.createElement('div'); g.className='sd-glue'; g.textContent=id; r.appendChild(g); return; }
        const m=byId[id];
        const box=document.createElement('div'); box.className='sd-mat';
        box.innerHTML=`<div class="sd-mt">${m.title}</div>`;
        const g=document.createElement('div'); g.className='sd-grid';
        g.style.gridTemplateColumns=(m.rowLabels?'auto ':'')+`repeat(${m.cols},var(--sd))`;
        m.labEls=[]; m.cellEls=[];
        for(let i=0;i<m.rows;i++){
          if(m.rowLabels){ const l=document.createElement('div'); l.className='sd-rl'; l.innerHTML=`<i>R</i><sub>${i+1}</sub>`; g.appendChild(l); m.labEls.push(l); }
          m.cellEls.push([]);
          for(let j=0;j<m.cols;j++){
            const c=document.createElement('div'); c.className='sd-c'; c._t=null; c._c=null;
            if(m.bar!=null&&j===m.bar+1) c.dataset.bar='1';
            g.appendChild(c); m.cellEls[i].push(c);
          }
        }
        box.appendChild(g); r.appendChild(box);
      });
      stage.appendChild(r);
    });
    // phases + legend
    $('.sd-phases').innerHTML=data.meta.phases.map((p,i)=>`<span class="sd-ph" data-i="${i}"><b>${i+1}</b> ${p}</span>`).join('<span class="sd-arrow">›</span>');
    $('.sd-legend').innerHTML=data.meta.legend.map(l=>`<span><span class="sd-sw ${l[0]}"></span> ${l[1]}</span>`).join('');
    scrub.max=data.steps.length-1; scrub.value=0;
    render(); first=false;
  }

  function render(){
    const st=data.steps[k];
    data.meta.mats.forEach(m=>{
      const snap=st.mats[m.id];
      for(let i=0;i<m.rows;i++){
        for(let j=0;j<m.cols;j++){
          const el=m.cellEls[i][j], cell=snap.cells[i][j];
          const cls='sd-c '+(cell.c||'');
          if(el._c!==cls){ el.className=cls; el._c=cls; }
          if(el._t!==cell.t){
            el.innerHTML=cell.t; el._t=cell.t;
            if(!first){ el.classList.remove('pop'); void el.offsetWidth; el.classList.add('pop'); }
          }
        }
        if(m.labEls.length){ m.labEls[i].className='sd-rl '+((snap.lab&&snap.lab[i])||''); }
      }
    });
    root.querySelectorAll('.sd-ph').forEach(p=>{
      const i=+p.dataset.i; p.classList.toggle('active',i===st.phase); p.classList.toggle('done',i<st.phase||(k===data.steps.length-1));
    });
    $('.sd-ct').innerHTML=st.cap.t; $('.sd-cw').innerHTML=st.cap.w;
    $('.sd-cg').innerHTML=st.cap.g?`<span class="sd-tag">general</span> ${st.cap.g}`:'';
    $('.sd-ci').innerHTML=st.cap.i?`<span class="sd-tag now">this step</span> ${st.cap.i}`:'';
    const last=data.steps.length-1;
    scrub.value=k; scrub.style.setProperty('--p',(last?100*k/last:0)+'%');
    $('.sd-count').textContent=`step ${k} / ${last}`;
    $('.sd-back').disabled=k===0; $('.sd-step').disabled=k===last;
  }

  function go(i){ k=Math.max(0,Math.min(data.steps.length-1,i)); render(); }
  function stop(){ if(timer){clearInterval(timer);timer=null;} autoBtn.textContent='Auto ▶'; }
  function play(){
    if(k>=data.steps.length-1) go(0);
    autoBtn.textContent='Pause ⏸';
    timer=setInterval(()=>{ if(k>=data.steps.length-1){stop();return;} go(k+1); },speed);
  }
  $('.sd-step').onclick=()=>{ stop(); go(k+1); };
  $('.sd-back').onclick=()=>{ stop(); go(k-1); };
  $('.sd-reset').onclick=()=>{ stop(); go(0); };
  autoBtn.onclick=()=>{ timer?stop():play(); };
  scrub.oninput=()=>{ stop(); go(+scrub.value); };
  nSel.onchange=()=>{ n=+nSel.value; setup(); };
  $('.sd-sp').onchange=e=>{ speed=+e.target.value; if(timer){stop();play();} };
  setup();
}

function init(){ document.querySelectorAll('[data-sd]').forEach(mount); }
if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',init); else init();
})();
