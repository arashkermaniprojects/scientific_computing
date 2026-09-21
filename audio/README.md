# Module 2 voice tracks

Per-section English narration for `module2.html`, generated with Google Cloud
Text-to-Speech (Gemini `gemini-2.5-flash-tts`, voice **Charon**) — the same API
used for the *Why Teams Don't Work* audiobook.

## Files
- `narration.json` — the narration text for each section (edit here to change wording).
- `synth.py` — generates one MP3 per section into `../assets/audio/m2_<id>.mp3`.

## Generate the audio (run on your Mac, where your Google credentials live)

```bash
cd "BCS311_Course_Website"

python3 -m venv .venv && source .venv/bin/activate
pip install google-cloud-texttospeech

# authenticate ONE of these ways:
gcloud auth application-default login
#   – or –
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json

python3 audio/synth.py                 # all 10 sections
python3 audio/synth.py --only fixed    # regenerate just one section
```

The script is **resumable**: it skips any `m2_<id>.mp3` that already exists, so
to redo a section, delete its file (or edit the text) and run again.

## How it shows up on the page
`module2.html` has a small "🔊 Listen to this section" player under each main
section heading, pointing at `assets/audio/m2_<id>.mp3`. Until the files exist
those players hide themselves automatically, so the page looks clean before the
audio is generated and starts working the moment the MP3s appear.

## Section → file map
| id | file | section |
|----|------|---------|
| intro | m2_intro.mp3 | Root finding: the problem |
| fixed | m2_fixed.mp3 | 2.1 Fixed-point iteration |
| fixed_theory | m2_fixed_theory.mp3 | 2.1.1 Why it converges — and from where |
| bisect | m2_bisect.mp3 | 2.2 Bisection |
| newton | m2_newton.mp3 | 2.3 Newton's method |
| secant | m2_secant.mp3 | 2.4 The secant method |
| compare | m2_compare.mp3 | 2.5 Choosing a method · polynomial roots |
| gauss | m2_gauss.mp3 | 2.6 Gaussian elimination |
| lu | m2_lu.mp3 | 2.7 LU factorization |
| chol | m2_chol.mp3 | 2.8 Cholesky factorization |

After generating, commit `assets/audio/*.mp3` and push — the site will play them.
