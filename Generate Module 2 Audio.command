#!/bin/bash
# Self-contained: writes its own synth script to /tmp, runs it there (no reads of the
# iCloud/Documents-protected course folder), then copies the finished MP3s into the site.
SRC="/Users/arashkermanikolankeh/Documents/Work/Teaching/CUD/FA_2026/Scientific Computing/BCS311_Course_Website"
LOG="$SRC/run_tts.log"
{
  echo "=== $(date) ==="
  PY="$HOME/Documents/AudioBooks/WhyTeamsDontWork/.venv/bin/python"
  echo "python: $PY"
  rm -rf /tmp/m2tts && mkdir -p /tmp/m2tts/out
  cat > /tmp/m2tts/synth_inline.py <<'PYEOF'
import os, sys, time
from google.cloud import texttospeech as tts
VOICE='Charon'; MODEL='gemini-2.5-flash-tts'; LANG='en-US'
PROMPT=('Read aloud as a professional English audiobook narrator for a university '
        'numerical-methods lecture: calm, clear, unhurried, natural intonation. '
        'Pronounce mathematical phrases naturally. Read the text exactly as written, '
        'without adding or skipping anything.')
OUT='/tmp/m2tts/out'; os.makedirs(OUT, exist_ok=True)
SECTIONS=[('intro', "This module is about solving equations. Many problems reduce to finding a root: a value of x for which f of x equals zero. Graphically, that is where the curve y equals f of x crosses the x axis. For simple linear or quadratic equations we can solve by hand, but most non-linear equations, like cosine x equals x, or x cubed minus two x minus one equals zero, have no closed-form solution. So we approximate the root with an iterative method: start from a guess, improve it repeatedly, and stop when the change falls below a small tolerance. The four methods that follow, fixed-point iteration, bisection, Newton's method, and the secant method, trade off robustness, meaning will it converge, against speed, meaning how fast."), ('fixed', 'Fixed-point iteration rewrites f of x equals zero into the form x equals g of x. A solution is called a fixed point: a value the map returns to itself, so g of x-star equals x-star. We then iterate. x sub k plus one equals g of x sub k, and stop when the change is below the tolerance epsilon. Two pictures are worth holding at once. In the fixed-point view, we plot y equals g of x together with the line y equals x; the fixed point is where they cross, and the iteration traces a cobweb. In the root view, we plot y equals f of x, which is g minus x; the root is where that curve meets zero. The fixed point and the root are the same number. One step of the cobweb is two operations: go up to the curve to evaluate g, then across to the line to feed the result back as the next input. Crucially, the method converges only when the absolute value of g-prime at the root is less than one.'), ('fixed_theory', 'Why does it converge, and from where? Write the error at step k as e sub k. A one-term Taylor expansion shows that the next error is approximately g-prime at the root times the current error. So each step multiplies the error by g-prime. If the absolute value of g-prime is less than one, the error shrinks and the fixed point attracts; if it is greater than one, the point repels, and no starting guess will reach it. The interval, or contraction, condition tells us which starting points are safe: if g maps an interval into itself, and the slope stays below one throughout, the iteration converges for every start in that interval. The set of all starting points that reach a given root is called its basin of attraction. And because a rearrangement makes some roots attracting and others repelling, the choice of g decides which root you can reach, and from where.'), ('bisect', 'Bisection is the reliable bracket method. It needs two guesses, a and b, that bracket the root, meaning f of a and f of b have opposite signs; then, by the intermediate value theorem, a continuous function must cross zero between them. Compute the midpoint c equals a plus b over two, check the sign of f of c, and keep whichever half still changes sign. Repeat. The bracket halves every step, so the error after n cuts is at most b minus a over two to the n, about one binary digit per step. Bisection always converges for a sign-changing continuous function, but only linearly, so it is slow. It is the dependable fallback when faster methods misbehave.'), ('newton', "Newton's method uses the slope to leap toward the root. From a guess x sub k, draw the tangent line to the curve; where that tangent crosses the x axis is the next guess. The formula is: x sub k plus one equals x sub k minus f of x sub k, divided by f-prime of x sub k. It needs the derivative and a reasonable starting guess, but when those are available it converges quadratically, meaning the number of correct digits roughly doubles each step. The warning: if the derivative is near zero, the step explodes, and a poor start can diverge or jump to a different root. That is why robust solvers combine Newton with a bracket."), ('secant', "The secant method is Newton's method with the derivative replaced by a finite-difference estimate from the two most recent points, so no derivative is needed. Draw the straight line through the last two points on the curve, and take its x-intercept as the next guess. It needs two starting values, and converges superlinearly, with order about one point six one eight, slower than Newton but faster than bisection. It is the practical compromise when the derivative is awkward or expensive to compute."), ('compare', "How do you choose? Bisection is guaranteed but linear; fixed-point is linear and only conditional; Newton is quadratic but needs the derivative and a good start; the secant is superlinear and derivative-free. A common strategy is to use bisection to safely locate a bracket, then switch to Newton or the secant for fast final convergence, which is exactly what hybrid methods, like MATLAB's f-zero, do internally. For polynomials, MATLAB's roots function returns all the roots at once from the coefficient vector."), ('gauss', "Now we turn from one equation to systems: A x equals b. Gaussian elimination is the workhorse direct method. Using elementary row operations, swapping rows, scaling a row, and adding a multiple of one row to another, it transforms the augmented matrix into upper-triangular form. First, forward elimination introduces zeros below each pivot, column by column. Then back-substitution solves the last equation for the last unknown and works upward. One caution: if a pivot is zero or tiny, the multiplier blows up; partial pivoting, swapping in the row with the largest pivot, keeps elimination stable. MATLAB's backslash operator does this automatically."), ('lu', "L-U decomposition factors a square matrix into a lower-triangular L, with ones on its diagonal, times an upper-triangular U. This is Gaussian elimination's bookkeeping made reusable. Once you have L and U, solving A x equals b is two cheap triangular solves: forward-solve L y equals b, then back-solve U x equals y. The advantage is that the expensive elimination is done once; to solve the same matrix against many different right-hand sides, you reuse L and U and only repeat the two fast solves."), ('chol', 'When the matrix is symmetric and positive definite, meaning all its eigenvalues are positive, it has a special, cheaper factorization: A equals L times L transpose. Cholesky is about twice as fast as L-U and very stable, so it is the method of choice for the symmetric positive-definite systems that arise in least squares and many engineering problems. Solving uses the same two-step idea: L y equals b, then L transpose x equals y. If the Cholesky routine errors, the matrix is not positive definite, and you fall back to L-U.')]
client=tts.TextToSpeechClient()
def synth(text,path):
    for a in range(6):
        try:
            r=client.synthesize_speech(input=tts.SynthesisInput(text=text,prompt=PROMPT),
                voice=tts.VoiceSelectionParams(language_code=LANG,name=VOICE,model_name=MODEL),
                audio_config=tts.AudioConfig(audio_encoding=tts.AudioEncoding.MP3,speaking_rate=1.0),timeout=180)
            open(path,'wb').write(r.audio_content); return True
        except Exception as e:
            print('  retry',a+1,str(e)[:120],flush=True); time.sleep(min(30,2**a))
    return False
ok=0
for i,(sid,text) in enumerate(SECTIONS,1):
    p=os.path.join(OUT,'m2_%s.mp3'%sid)
    if os.path.exists(p) and os.path.getsize(p)>0:
        print('[%d/%d] skip %s'%(i,len(SECTIONS),sid),flush=True); ok+=1; continue
    print('[%d/%d] %s (%d chars)'%(i,len(SECTIONS),sid,len(text)),flush=True)
    if synth(text,p): ok+=1; print('   ok ->',os.path.getsize(p),'bytes',flush=True)
print('DONE %d/%d'%(ok,len(SECTIONS)),flush=True)
PYEOF
  echo "--- running Gemini-TTS (English, Charon) ---"
  "$PY" /tmp/m2tts/synth_inline.py
  echo "--- copying MP3s into the site ---"
  mkdir -p "$SRC/assets/audio"
  cp -v /tmp/m2tts/out/*.mp3 "$SRC/assets/audio/" 2>&1
  echo "--- site audio dir now ---"
  ls -la "$SRC/assets/audio/" 2>&1
  echo "=== done ==="
} > "$LOG" 2>&1
echo "Finished. See run_tts.log in the course folder."
