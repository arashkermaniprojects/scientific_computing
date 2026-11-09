/* Module 5 SVG animations: adaptive quadrature (5.4) and slope field / ODE (5.5).
   Each figure is a pure function of time t in [0,1]; play / pause / scrub just move t. */
(function(){
'use strict';
var NS='http://www.w3.org/2000/svg';
var C={teal:'#0e7c86',org:'#f0913a',ink:'#1d2b3a',grid:'#c9d4e3',soft:'#7b8aa0',red:'#d8553b'};
var W=560,H0=300,L=48,R=20,T=16;
function el(tag,attrs,parent){var e=document.createElementNS(NS,tag);for(var k in attrs)e.setAttribute(k,attrs[k]);if(parent)parent.appendChild(e);return e;}
function halo(a){a['stroke']='#fff';a['stroke-width']=3;a['paint-order']='stroke';return a;}
function axes(svg,x0,x1,y0,y1,xt,yt,hh,bb){
  var H=hh||H0, B=(bb===undefined?40:bb);
  var X=function(x){return L+(x-x0)/(x1-x0)*(W-L-R);}, Y=function(y){return H-B-(y-y0)/(y1-y0)*(H-B-T);};
  el('rect',{width:W,height:H,rx:10,fill:'#fff'},svg);
  el('line',{x1:L,y1:H-B,x2:W-R,y2:H-B,stroke:C.ink,'stroke-width':1.4},svg);
  el('line',{x1:L,y1:T,x2:L,y2:H-B,stroke:C.ink,'stroke-width':1.4},svg);
  xt.forEach(function(t){el('line',{x1:X(t[0]),y1:H-B,x2:X(t[0]),y2:H-B+5,stroke:C.ink},svg);
    el('text',{x:X(t[0]),y:H-B+19,'font-size':11.5,fill:C.ink,'text-anchor':'middle'},svg).textContent=t[1];});
  yt.forEach(function(t){el('line',{x1:L-5,y1:Y(t[0]),x2:L,y2:Y(t[0]),stroke:C.ink},svg);
    el('text',{x:L-8,y:Y(t[0])+4,'font-size':11.5,fill:C.ink,'text-anchor':'end'},svg).textContent=t[1];});
  return {X:X,Y:Y,pw:W-L-R,ph:H-B-T};
}
function player(fig,draw,dur){
  var btn=fig.querySelector('button.anim-play'), sl=fig.querySelector('input.anim-scrub');
  var t=0,raf=0,last=0,playing=false;
  function set(v){t=Math.max(0,Math.min(1,v));sl.value=Math.round(t*1000);draw(t);}
  function label(){btn.textContent=playing?'❚❚ Pause':(t>=1?'⟳ Replay':'▶ Play');}
  function tick(ts){if(!playing)return;var dt=(ts-last)/1000;last=ts;set(t+dt/dur);
    if(t>=1){playing=false;label();}else raf=requestAnimationFrame(tick);}
  function play(){if(t>=1)set(0);playing=true;label();last=performance.now();raf=requestAnimationFrame(tick);}
  function stop(){playing=false;cancelAnimationFrame(raf);label();}
  btn.addEventListener('click',function(){playing?stop():play();});
  sl.addEventListener('input',function(){stop();set(sl.value/1000);label();});
  var reduce=window.matchMedia&&window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  set(reduce?1:0);label();
  if(!reduce&&'IntersectionObserver' in window){
    var io=new IntersectionObserver(function(es){es.forEach(function(e){if(e.isIntersecting){io.disconnect();if(t===0)play();}});},{threshold:.5});
    io.observe(fig);
  }
  return {seek:function(v){stop();set(v);label();}};
}

/* ---------------- 5.4 adaptive quadrature ---------------- */
function initQuad(fig){
  var f=function(x){return Math.exp(-x/3)*Math.sin(0.9*x*x);};
  var a=0,b=4;
  function simp(p,q){var m=(p+q)/2;return (q-p)/6*(f(p)+4*f(m)+f(q));}
  function ref(){var n=20000,h=(b-a)/n,s=f(a)+f(b);for(var i=1;i<n;i++)s+=f(a+i*h)*(i%2?4:2);return s*h/3;}
  var exact=ref(), tol=1e-5;
  var frames=[], ivs=[{a:a,b:b,born:false}];
  for(var lvl=0;lvl<12;lvl++){
    var est=0,open=0;
    ivs.forEach(function(v){var m=(v.a+v.b)/2,S1=simp(v.a,v.b),S2=simp(v.a,m)+simp(m,v.b);
      var allowed=tol*(v.b-v.a)/(b-a);
      v.ratio=(Math.abs(S2-S1)/15)/allowed; v.est=S2; v.split=v.ratio>1; est+=S2; if(v.split)open++;});
    frames.push({ivs:ivs.map(function(v){return {a:v.a,b:v.b,born:v.born,split:v.split,ratio:v.ratio};}),est:est,open:open});
    if(!open)break;
    var next=[];
    ivs.forEach(function(v){if(v.split){var m=(v.a+v.b)/2;next.push({a:v.a,b:m,born:true},{a:m,b:v.b,born:true});}
      else next.push({a:v.a,b:v.b,born:false});});
    ivs=next;
  }
  var svg=fig.querySelector('svg'), rd=fig.querySelector('.anim-read');
  var HH=390, BB=140, PB=HH-BB;                         /* function panel bottom edge */
  svg.setAttribute('viewBox','0 0 '+W+' '+HH);
  var ax=axes(svg,0,4.2,-1.2,1.9,[[0,'a'],[4,'b']],[],HH,BB);
  var fills=el('g',{},svg);
  el('line',{x1:ax.X(0),y1:ax.Y(0),x2:ax.X(4.2),y2:ax.Y(0),stroke:C.grid,'stroke-width':1},svg);
  var fine=el('g',{},svg), coarse=el('g',{},svg), divs=el('g',{},svg);
  var d='';for(var i=0;i<=240;i++){var x=4*i/240;d+=(i?'L':'M')+ax.X(x).toFixed(1)+','+ax.Y(f(x)).toFixed(1);}
  el('path',{d:d,fill:'none',stroke:C.ink,'stroke-width':1.6,'stroke-linejoin':'round',opacity:.55},svg);
  /* legend (function panel) */
  var lg=el('g',{},svg), ly=T+6;
  el('line',{x1:L+10,y1:ly+4,x2:L+34,y2:ly+4,stroke:C.blue||'#2563a8','stroke-width':2.4,'stroke-dasharray':'5 3'},lg);
  el('text',halo({x:L+40,y:ly+8,'font-size':11.5,fill:C.ink,'font-weight':600}),lg).textContent='ONE parabola over the whole strip';
  el('line',{x1:L+10,y1:ly+22,x2:L+34,y2:ly+22,stroke:C.ink,'stroke-width':2.4},lg);
  el('text',halo({x:L+40,y:ly+26,'font-size':11.5,fill:C.ink,'font-weight':600}),lg).textContent='TWO parabolas, one per half (colour = verdict)';
  el('text',halo({x:L+10,y:ly+44,'font-size':11.5,fill:C.ink,'font-weight':700}),lg).textContent='Do the two agree? yes \u2192 settled (teal)   no \u2192 still changing (orange)';
  /* bar panel */
  var bt=PB+46, bb=HH-14, bh=bb-bt;                    /* bars from baseline bb up to bt */
  el('text',{x:L,y:PB+40,'font-size':11.5,fill:C.ink,'font-weight':700},svg).textContent='Change when each strip is halved';
  el('line',{x1:L,y1:bb,x2:W-R,y2:bb,stroke:C.ink,'stroke-width':1.2},svg);
  var thr=bb-bh*0.5;
  var bars=el('g',{},svg);
  el('line',{x1:L,y1:thr,x2:W-R,y2:thr,stroke:C.red,'stroke-width':1.4,'stroke-dasharray':'6 4'},svg);
  el('line',{x1:W-R-118,y1:PB+36,x2:W-R-98,y2:PB+36,stroke:C.red,'stroke-width':1.6,'stroke-dasharray':'5 3'},svg);
  el('text',{x:W-R-94,y:PB+40,'font-size':11,fill:C.red,'font-weight':700},svg).textContent='allowed change';
  var minw=b-a;frames[frames.length-1].ivs.forEach(function(v){minw=Math.min(minw,v.b-v.a);});
  var uniform=Math.round((b-a)/minw);
  function para(x0,x1,x2,steps){var X=[x0,x1,x2],Y=X.map(f),d='';
    for(var j=0;j<=steps;j++){var x=x0+(x2-x0)*j/steps,y=0;
      for(var i=0;i<3;i++){var l=1;for(var k=0;k<3;k++)if(k!==i)l*=(x-X[k])/(X[i]-X[k]);y+=Y[i]*l;}
      d+=(j?'L':'M')+ax.X(x).toFixed(1)+','+ax.Y(y).toFixed(1);}
    return d;}
  function draw(t){
    var n=frames.length, s=t*(n-1), k=Math.min(n-1,Math.floor(s+1e-9)), fr=Math.min(1,s-k);
    if(t>=1){k=n-1;fr=1;}
    var F=frames[k]; fills.textContent='';fine.textContent='';coarse.textContent='';divs.textContent='';bars.textContent='';
    F.ivs.forEach(function(v){
      var col=v.split?C.org:C.teal, flash=v.born?.3*(1-Math.min(1,fr*2.5)):0;
      var p='M'+ax.X(v.a).toFixed(1)+','+ax.Y(0).toFixed(1), m=Math.max(4,Math.ceil((v.b-v.a)*40));
      for(var j=0;j<=m;j++){var x=v.a+(v.b-v.a)*j/m;p+='L'+ax.X(x).toFixed(1)+','+ax.Y(f(x)).toFixed(1);}
      p+='L'+ax.X(v.b).toFixed(1)+','+ax.Y(0).toFixed(1)+'Z';
      el('path',{d:p,fill:col,'fill-opacity':(v.split?.22:.14)+flash,stroke:'none'},fills);
      var mid=(v.a+v.b)/2, q1=(v.a+mid)/2, q3=(mid+v.b)/2;
      el('path',{d:para(v.a,mid,v.b,24),fill:'none',stroke:'#2563a8','stroke-width':2.2,'stroke-dasharray':'5 3','stroke-linecap':'round'},coarse);
      el('path',{d:para(v.a,q1,mid,12)+para(mid,q3,v.b,12).replace('M','L'),fill:'none',stroke:col,'stroke-width':2.2},fine);
      el('line',{x1:ax.X(v.a),y1:ax.Y(0),x2:ax.X(v.a),y2:ax.Y(f(v.a)),stroke:col,'stroke-width':1,opacity:v.born?Math.min(1,fr*3):.8},divs);
      /* bar: log scale, limit (ratio 1) sits at half height */
      var h=Math.max(0.03,Math.min(1,(Math.log(Math.max(v.ratio,1e-9))/Math.LN10+2)/4));
      var x0=ax.X(v.a)+1, w=Math.max(1.5,ax.X(v.b)-ax.X(v.a)-2);
      el('rect',{x:x0,y:bb-bh*h,width:w,height:bh*h,fill:col,'fill-opacity':.85,rx:1.5},bars);
    });
    var last=F.ivs[F.ivs.length-1];
    el('line',{x1:ax.X(last.b),y1:ax.Y(0),x2:ax.X(last.b),y2:ax.Y(f(last.b)),stroke:C.teal,'stroke-width':1,opacity:.8},divs);
    var msg='Stage '+k+' of '+(n-1)+' \u00B7 '+F.ivs.length+' strip'+(F.ivs.length>1?'s':'')+' \u00B7 estimate '+F.est.toFixed(6)+' (reference '+exact.toFixed(6)+')';
    msg+= F.open? ' \u00B7 '+F.open+' strip'+(F.open>1?'s':'')+' above the limit (orange bars) \u2192 split next' : ' \u00B7 every bar is below the limit \u2014 done. A uniform grid this fine would need '+uniform+' strips.';
    rd.textContent=msg;
  }
  var p=player(fig,draw,Math.max(6,frames.length*1.3));
  p.frames=frames;p.exact=exact;
  return p;
}

/* ---------------- 5.5 slope field, family, initial condition ---------------- */
function initODE(fig){
  var X0=0,X1=2.2,Y0=0,Y1=5.6;
  var svg=fig.querySelector('svg'), rd=fig.querySelector('.anim-read');
  var ax=axes(svg,X0,X1,Y0,Y1,[[0,'0'],[1,'1'],[2,'2']],[[1,'1'],[2,'2'],[4,'4']]);
  var cid='odeclip'+Math.floor(Math.random()*1e6);
  var cp=el('clipPath',{id:cid},svg);el('rect',{x:L,y:T,width:ax.pw,height:ax.ph},cp);
  var kx=(ax.ph/(Y1-Y0))/(ax.pw/(X1-X0));          // pixel slope = -slope*kx
  function seg(cx,cy,s,half){var dx=1,dy=-s*kx,n=Math.hypot(dx,dy);return [cx-dx/n*half,cy-dy/n*half,cx+dx/n*half,cy+dy/n*half];}
  var segs=[];
  var field=el('g',{'clip-path':'url(#'+cid+')'},svg);
  for(var i=0;i<7;i++)for(var j=0;j<8;j++){
    var gx=0.15+i*0.3, gy=0.4+0.7*j, q=seg(ax.X(gx),ax.Y(gy),gy,9);
    segs.push({o:i*0.5+j*0.5, e:el('line',{x1:q[0],y1:q[1],x2:q[2],y2:q[3],stroke:C.soft,'stroke-width':1.6,'stroke-linecap':'round',opacity:0},field)});
  }
  var cs=[0.35,0.6,1.0,1.7], curves=[], g2=el('g',{'clip-path':'url(#'+cid+')'},svg);
  cs.forEach(function(c){
    var xe=Math.min(X1,Math.log(Y1*1.02/c)),d='';
    for(var i=0;i<=120;i++){var x=xe*i/120;d+=(i?'L':'M')+ax.X(x).toFixed(1)+','+ax.Y(c*Math.exp(x)).toFixed(1);}
    curves.push({c:c,e:el('path',{d:d,fill:'none',stroke:c===1?C.teal:C.grid,'stroke-width':c===1?3.4:2.2,pathLength:1,'stroke-dasharray':'1 1','stroke-dashoffset':1,'stroke-linecap':'round'},g2)});
  });
  var lab=[[0.35,2.15,0.35],[0.6,2.0,0.6],[1,1.55,1],[1.7,1.0,1.7]], labs=lab.map(function(l){
    var x=l[1],y=l[2]*Math.exp(x);
    return el('text',halo({x:ax.X(x)-6,y:ax.Y(y)-6,'font-size':11.5,fill:C.soft,'font-weight':700,'text-anchor':'end',opacity:0}),svg);}); 
  labs.forEach(function(t,i){t.textContent='c = '+lab[i][0];});
  var dot=el('circle',{cx:ax.X(0),cy:ax.Y(1),r:6,fill:C.org,stroke:'#fff','stroke-width':1.5,opacity:0},svg);
  var ring=el('circle',{cx:ax.X(0),cy:ax.Y(1),r:6,fill:'none',stroke:C.org,'stroke-width':2,opacity:0},svg);
  var tan=el('line',{x1:0,y1:0,x2:0,y2:0,stroke:C.org,'stroke-width':3.4,'stroke-linecap':'round',opacity:0},svg);
  var mk=el('circle',{r:5.5,fill:C.org,stroke:'#fff','stroke-width':1.5,opacity:0},svg);
  var XM=1.7;
  function ease(u){return u<0?0:u>1?1:u;}
  function draw(t){
    var p1=ease(t/.18), p2=ease((t-.18)/.28), p3=ease((t-.46)/.12), p4=ease((t-.58)/.42);
    segs.forEach(function(s){var u=ease(p1*3.2-s.o/10);s.e.setAttribute('opacity',(0.9*u).toFixed(2));});
    curves.forEach(function(cv,i){
      var u=ease(p2*cs.length-i), teal=cv.c===1;
      cv.e.setAttribute('stroke-dashoffset',(1-u).toFixed(3));
      cv.e.setAttribute('opacity',teal?1:(1-0.55*p3).toFixed(2));
      cv.e.setAttribute('stroke',teal&&p3===0?C.grid:(teal?C.teal:C.grid));
      cv.e.setAttribute('stroke-width',teal?(2.2+1.2*p3).toFixed(2):2.2);
    });
    labs.forEach(function(l,i){l.setAttribute('opacity',(ease(p2*cs.length-i)*(1-0.4*p3)).toFixed(2));});
    dot.setAttribute('opacity',p3>0?1:0);
    ring.setAttribute('opacity',p3>0&&p3<1?(1-p3).toFixed(2):0);ring.setAttribute('r',(6+16*p3).toFixed(1));
    var msg;
    if(t<=0)msg='① The ODE gives a slope at every point (x, y).';
    else if(p4<=0&&p2<1&&p3<=0)msg=(p2<=0?'① The ODE dy/dx = y gives a slope at every point (x, y): the short grey segments.':'② Following the slopes gives a whole family of solutions, y = c eˣ (one curve per constant c).');
    else if(p4<=0)msg='③ The initial condition y(0) = 1 picks one curve from the family (c = 1).';
    else{
      var x=XM*p4,y=Math.exp(x);
      var q=seg(ax.X(x),ax.Y(y),y,22);
      tan.setAttribute('x1',q[0]);tan.setAttribute('y1',q[1]);tan.setAttribute('x2',q[2]);tan.setAttribute('y2',q[3]);tan.setAttribute('opacity',1);
      mk.setAttribute('cx',ax.X(x));mk.setAttribute('cy',ax.Y(y));mk.setAttribute('opacity',1);
      msg='④ Walk along the chosen curve: at x = '+x.toFixed(2)+', y = '+y.toFixed(2)+', the slope the ODE demands is dy/dx = y = '+y.toFixed(2)+' — the orange tangent always lies along the grey field.';
    }
    if(p4<=0){tan.setAttribute('opacity',0);mk.setAttribute('opacity',0);}
    rd.textContent=msg;
  }
  return player(fig,draw,16);
}

window.BCS_ANIM={};
function boot(){
  var q=document.getElementById('animQuad'), o=document.getElementById('animODE');
  if(q)window.BCS_ANIM.quad=initQuad(q);
  if(o)window.BCS_ANIM.ode=initODE(o);
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',boot);else boot();
})();
