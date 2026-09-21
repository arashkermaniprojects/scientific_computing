"""
Synthesize the Module 2 section narrations into per-section MP3 files with
Google Cloud Text-to-Speech (the same Gemini-TTS setup used for the audiobook).

Reads   audio/narration.json   ->   writes   assets/audio/m2_<id>.mp3
Resumable: sections whose MP3 already exists are skipped.

Run it on YOUR machine (where your Google credentials live), from the
BCS311_Course_Website folder:

    # one-time setup
    python3 -m venv .venv && source .venv/bin/activate
    pip install google-cloud-texttospeech
    #   auth: either
    #     gcloud auth application-default login
    #   or
    #     export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json

    python3 audio/synth.py            # generate all sections
    python3 audio/synth.py --only fixed   # just one section (for testing)
"""
import json
import os
import sys
import time

from google.cloud import texttospeech as tts

# ---- voice / model (same family as the audiobook, now in English) ----
VOICE = "Charon"
MODEL = "gemini-2.5-flash-tts"
LANG  = "en-US"
PROMPT = ("Read aloud as a professional English audiobook narrator for a university "
          "numerical-methods lecture: calm, clear, unhurried, natural intonation. "
          "Pronounce mathematical phrases naturally. Read the text exactly as written, "
          "without adding or skipping anything.")

HERE       = os.path.dirname(os.path.abspath(__file__))
ROOT       = os.path.dirname(HERE)                       # BCS311_Course_Website
NARRATION  = os.path.join(HERE, "narration.json")
OUT_DIR    = os.path.join(ROOT, "assets", "audio")

only = None
if "--only" in sys.argv:
    only = sys.argv[sys.argv.index("--only") + 1]

client = tts.TextToSpeechClient()


def synth(text, path):
    for attempt in range(60):
        try:
            r = client.synthesize_speech(
                input=tts.SynthesisInput(text=text, prompt=PROMPT),
                voice=tts.VoiceSelectionParams(language_code=LANG, name=VOICE, model_name=MODEL),
                audio_config=tts.AudioConfig(audio_encoding=tts.AudioEncoding.MP3, speaking_rate=1.0),
                timeout=180,
            )
            tmp = path + ".part"
            open(tmp, "wb").write(r.audio_content)
            os.replace(tmp, path)
            return
        except Exception as e:                          # quota / transient errors
            wait = min(60, 2 ** attempt)
            print(f"  retry {attempt+1} for {os.path.basename(path)} in {wait}s: {str(e)[:90]}", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed: {path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    sections = json.load(open(NARRATION, encoding="utf-8"))
    todo = [s for s in sections if only is None or s["id"] == only]
    print(f"{len(sections)} sections; generating {len(todo)}", flush=True)
    t0 = time.time()
    for i, s in enumerate(todo, 1):
        path = os.path.join(OUT_DIR, f"m2_{s['id']}.mp3")
        if os.path.exists(path):
            print(f"  [{i}/{len(todo)}] skip {os.path.basename(path)} (exists)", flush=True)
            continue
        print(f"  [{i}/{len(todo)}] {s['id']}: {s['title']} ({len(s['text'])} chars)", flush=True)
        synth(s["text"], path)
    print(f"done in {time.time()-t0:.0f}s -> {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
