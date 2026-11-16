/* Module 5, section 5.4: interactive Gauss-Legendre demo and adaptive-Simpson demo (plain SVG, no libraries). */
(function(){
'use strict';
var NS='http://www.w3.org/2000/svg', W=560, H=300, L=46, R=16, T=14, B=34;
var C={teal:'#0e7c86',org:'#f0913a',ink:'#1d2b3a',grid:'#dbe3ee',soft:'#7b8aa0',red:'#d8553b',blue:'#2f5fa7'};
function el(t,a,p){var e=document.createElementNS(NS,t);for(var k in a)e.setAttribute(k,a[k]);if(p)p.appendChild(e);return e;}
function $(id){return document.getElementById(id);}
function fmt(x,d){return (Math.abs(x)<1e-3&&x!==0)?x.toExponential(2):x.toFixed(d===undefined?6:d);}

/* Gauss-Legendre nodes/weights on [-1,1] by Newton's method on P_n */
function gauss(n){
  var x=[],w=[];
  for(var i=0;i<n;i++){
    var z=Math.cos(Math.PI*(i+0.75)/(n+0.5)),pp=1;
    for(var it=0;it<100;it++){
      var p1=1,p2=0;
      for(var j=1;j<=n;j++){var p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j;}
      pp=n*(z*p1-p2)/(z*z-1);
      var dz=p1/pp;z-=dz;if(Math.abs(dz)<1e-15)break;
    }
    x.push(z);w.push(2/((1-z*z)*pp*pp));
  }
  var idx=x.map(function(_,i){return i;}).sort(function(a,b){return x[a]-x[b];});
  return {t:idx.map(function(i){return x[i];}),w:idx.map(function(i){return w[i];})};
}
function frame(svg,x0,x1,y0,y1,f,xt){
  var X=function(x){return L+(x-x0)/(x1-x0)*(W-L-R);},Y=function(y){return H-B-(y-y0)/(y1-y0)*(H-B-T);};
  el('rect',{width:W,height:H,rx:10,fill:'#fff'},svg);
  el('line',{x1:L,y1:Y(0),x2:W-R,y2:Y(0),stroke:C.ink,'stroke-width':1.2},svg);
  el('line',{x1:L,y1:T,x2:L,y2:H-B,stroke:C.ink,'stroke-width':1.2},svg);
  (xt||[]).forEach(function(v){el('line',{x1:X(v),y1:H-B,x2:X(v),y2:H-B+5,stroke:C.ink},svg);
    el('text',{x:X(v),y:H-B+19,'font-size':11.5,fill:C.ink,'text-anchor':'middle'},svg).textContent=+v.toFixed(2);});
  for(var k=0;k<=4;k++){var yv=y0+(y1-y0)*k/4;el('text',{x:L-6,y:Y(yv)+4,'font-size':11,fill:C.ink,'text-anchor':'end'},svg).textContent=+yv.toFixed(2);
    el('line',{x1:L,y1:Y(yv),x2:W-R,y2:Y(yv),stroke:C.grid,'stroke-width':.7},svg);}
  var d='',N=400;
  for(var i=0;i<=N;i++){var x=x0+(x1-x0)*i/N;d+=(i?'L':'M')+X(x).toFixed(1)+' '+Y(f(x)).toFixed(1);}
  return {X:X,Y:Y,curve:function(){el('path',{d:d,fill:'none',stroke:C.ink,'stroke-width':2.2},svg);}};
}
function range(f,a,b){var lo=0,hi=-1e9;for(var i=0;i<=400;i++){var v=f(a+(b-a)*i/400);if(v<lo)lo=v;if(v>hi)hi=v;}
  var pad=(hi-lo)*.08||1;return [Math.min(0,lo-pad*0.2),hi+pad];}
function refSimpson(f,a,b){var n=200000,h=(b-a)/n,s=f(a)+f(b);for(var i=1;i<n;i++)s+=f(a+i*h)*(i%2?4:2);return s*h/3;}

/* ---------------- Gauss demo ---------------- */
function initGauss(root){
  var FN=[
    {name:'eˣ on [0, 2]',f:function(x){return Math.exp(x);},a:0,b:2,exact:Math.exp(2)-1},
    {name:'sin x on [0, π]',f:Math.sin,a:0,b:Math.PI,exact:2},
    {name:'1/(1+x²) on [0, 4]',f:function(x){return 1/(1+x*x);},a:0,b:4,exact:Math.atan(4)},
    {name:'x⁵ on [0, 2]  (a polynomial of degree 5)',f:function(x){return Math.pow(x,5);},a:0,b:2,exact:64/6}
  ];
  var svg=root.querySelector('svg'),sel=root.querySelector('.g-fn'),sl=root.querySelector('.g-n'),nv=root.querySelector('.g-nval'),out=root.querySelector('.readout');
  sel.innerHTML=FN.map(function(o,i){return '<option value="'+i+'">'+o.name+'</option>';}).join('');
  function draw(){
    var o=FN[+sel.value],n=+sl.value;nv.textContent=n;
    var g=gauss(n),a=o.a,b=o.b,xs=g.t.map(function(t){return (b-a)/2*t+(a+b)/2;});
    var G=0;xs.forEach(function(x,i){G+=g.w[i]*o.f(x);});G*=(b-a)/2;
    var h=(b-a)/n,M=0;for(var i=0;i<n;i++)M+=o.f(a+(i+.5)*h);M*=h;
    svg.innerHTML='';var r=range(o.f,a,b);
    var fr=frame(svg,a,b,r[0],r[1],o.f,[a,(a+b)/2,b]);
    for(var i=0;i<n;i++){var xm=a+(i+.5)*h;el('line',{x1:fr.X(xm),y1:fr.Y(0),x2:fr.X(xm),y2:fr.Y(0)-10,stroke:C.soft,'stroke-width':2.5},svg);}
    fr.curve();
    xs.forEach(function(x,i){
      el('line',{x1:fr.X(x),y1:fr.Y(0),x2:fr.X(x),y2:fr.Y(o.f(x)),stroke:C.org,'stroke-width':1.3,'stroke-dasharray':'4 3'},svg);
      el('circle',{cx:fr.X(x),cy:fr.Y(o.f(x)),r:6,fill:C.org,stroke:'#fff','stroke-width':1.5},svg);
      el('rect',{x:fr.X(x)-4,y:fr.Y(0)-g.w[i]*14,width:8,height:g.w[i]*14,fill:C.teal,opacity:.85},svg);
    });
    el('text',{x:W-R-4,y:T+12,'font-size':11.5,fill:C.org,'text-anchor':'end','font-weight':700},svg).textContent='● Gauss nodes ('+n+')';
    el('text',{x:W-R-4,y:T+27,'font-size':11.5,fill:C.teal,'text-anchor':'end','font-weight':700},svg).textContent='▮ height = weight';
    el('text',{x:W-R-4,y:T+42,'font-size':11.5,fill:C.soft,'text-anchor':'end','font-weight':700},svg).textContent='| midpoint-rule sample points ('+n+')';
    var rows=g.t.map(function(t,i){return '<tr><td>'+(i+1)+'</td><td>'+t.toFixed(6)+'</td><td>'+g.w[i].toFixed(6)+'</td><td>'+xs[i].toFixed(5)+'</td></tr>';}).join('');
    out.innerHTML='<b>Exact</b> = '+o.exact.toFixed(8)+'<br><span style="color:'+C.org+'"><b>Gauss, '+n+' evaluations:</b></span> '+G.toFixed(8)+' &nbsp; error = '+fmt(Math.abs(G-o.exact))+
      '<br><span style="color:'+C.soft+'"><b>Midpoint rule, '+n+' evaluations:</b></span> '+M.toFixed(8)+' &nbsp; error = '+fmt(Math.abs(M-o.exact))+
      '<table class="data" style="margin:.6rem 0 0;font-size:.82rem;text-align:center"><thead><tr><th>i</th><th>node t<sub>i</sub> on [−1,1]</th><th>weight w<sub>i</sub></th><th>node x<sub>i</sub> on [a,b]</th></tr></thead><tbody>'+rows+'</tbody></table>';
  }
  sel.addEventListener('change',draw);sl.addEventListener('input',draw);draw();
}

/* ---------------- adaptive Simpson demo ---------------- */
function initAdapt(root){
  var FN=[
    {name:'e^(−x/3)·sin(0.9x²) on [0, 4]  (wiggles speed up)',f:function(x){return Math.exp(-x/3)*Math.sin(0.9*x*x);},a:0,b:4},
    {name:'√x on [0, 1]  (steep at the left end)',f:Math.sqrt,a:0,b:1},
    {name:'1/(1+100x²) on [−2, 2]  (narrow spike)',f:function(x){return 1/(1+100*x*x);},a:-2,b:2},
    {name:'sin x on [0, π]  (smooth everywhere)',f:Math.sin,a:0,b:Math.PI}
  ];
  var svg=root.querySelector('svg'),sel=root.querySelector('.a-fn'),sl=root.querySelector('.a-tol'),tv=root.querySelector('.a-tolval'),out=root.querySelector('.readout');
  sel.innerHTML=FN.map(function(o,i){return '<option value="'+i+'">'+o.name+'</option>';}).join('');
  var refs={};
  function adapt(f,a,b,tol){
    var evals=0,strips=[];
    function F(x){evals++;return f(x);}
    function rec(a,b,fa,fm,fb,tol,depth){
      var m=(a+b)/2,lm=(a+m)/2,rm=(m+b)/2,flm=F(lm),frm=F(rm);
      var S1=(b-a)/6*(fa+4*fm+fb),S2=(m-a)/6*(fa+4*flm+fm)+(b-m)/6*(fm+4*frm+fb);
      if(Math.abs(S2-S1)<=15*tol||depth>=40){strips.push({a:a,b:b,depth:depth});return S2+(S2-S1)/15;}
      return rec(a,m,fa,flm,fm,tol/2,depth+1)+rec(m,b,fm,frm,fb,tol/2,depth+1);
    }
    var fa=F(a),fb=F(b),fm=F((a+b)/2);
    var I=rec(a,b,fa,fm,fb,tol,0);
    strips.sort(function(p,q){return p.a-q.a;});
    return {I:I,evals:evals,strips:strips};
  }
  function uniformNeeded(f,a,b,exact,target){
    var n=2;
    for(;n<=1048576;n*=2){
      var h=(b-a)/n,s=f(a)+f(b);for(var i=1;i<n;i++)s+=f(a+i*h)*(i%2?4:2);
      if(Math.abs(s*h/3-exact)<=target)return n+1;
    }
    return null;
  }
  function draw(){
    var k=+sel.value,o=FN[k],e=+sl.value;tv.textContent='10^'+e;
    var tol=Math.pow(10,e);
    if(!(k in refs))refs[k]=refSimpson(o.f,o.a,o.b);
    var res=adapt(o.f,o.a,o.b,tol),exact=refs[k],err=Math.abs(res.I-exact);
    svg.innerHTML='';var r=range(o.f,o.a,o.b),lo=Math.min(r[0],-0.0),hi=r[1];
    var fr=frame(svg,o.a,o.b,lo,hi,o.f,[o.a,(o.a+o.b)/2,o.b]);
    var maxd=0;res.strips.forEach(function(s){if(s.depth>maxd)maxd=s.depth;});
    res.strips.forEach(function(s){
      var t=maxd?s.depth/maxd:0,col=t>0.66?C.org:(t>0.33?'#c8b25a':C.teal);
      var pts='M'+fr.X(s.a)+' '+fr.Y(0),N=Math.max(2,Math.min(30,Math.round((fr.X(s.b)-fr.X(s.a))/3)));
      for(var i=0;i<=N;i++){var x=s.a+(s.b-s.a)*i/N;pts+='L'+fr.X(x).toFixed(1)+' '+fr.Y(o.f(x)).toFixed(1);}
      pts+='L'+fr.X(s.b)+' '+fr.Y(0)+'Z';
      el('path',{d:pts,fill:col,opacity:.35,stroke:'#fff','stroke-width':.8},svg);
      el('line',{x1:fr.X(s.a),y1:H-B,x2:fr.X(s.a),y2:H-B+9,stroke:col,'stroke-width':2},svg);
    });
    fr.curve();
    el('text',{x:W-R-4,y:T+12,'font-size':11.5,fill:C.teal,'text-anchor':'end','font-weight':700},svg).textContent='■ wide strip (smooth part)';
    el('text',{x:W-R-4,y:T+27,'font-size':11.5,fill:C.org,'text-anchor':'end','font-weight':700},svg).textContent='■ narrow strip (wiggly part)';
    var uni=uniformNeeded(o.f,o.a,o.b,exact,Math.max(err,1e-13));
    var widths=res.strips.map(function(s){return s.b-s.a;}),wmin=Math.min.apply(null,widths),wmax=Math.max.apply(null,widths);
    out.innerHTML='<b>Tolerance</b> 10<sup>'+e+'</sup> &nbsp;→&nbsp; <b>'+res.strips.length+' strips</b>, <b>'+res.evals+' function evaluations</b>, error = '+fmt(err)+
      '<br>Strip widths range from '+wmin.toPrecision(3)+' (narrowest) to '+wmax.toPrecision(3)+' (widest).'+
      '<br>Uniform Simpson (equal strips) needs about <b>'+(uni===null?'more than 10⁶':uni)+'</b> evaluations for the same accuracy'+
      (uni?' (<b>'+(uni/res.evals).toFixed(1)+'×</b> more).':'.');
  }
  sel.addEventListener('change',draw);sl.addEventListener('input',draw);draw();
}
function boot(){var g=$('gaussDemo'),a=$('adaptDemo');if(g)initGauss(g);if(a)initAdapt(a);}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',boot);else boot();
})();
