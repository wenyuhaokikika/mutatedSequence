# Narration voice options

The first cut used Microsoft Edge TTS with `zh-CN-XiaoyiNeural`. Xiaoyi is
free and lively, but its prosody is deliberately cartoon-like: short clauses
land with even, almost metronomic pauses. For a seven-minute explainer aimed at
non-specialists, that reads as “robot reading a script”.

The committed cut now uses **`zh-CN-XiaoxiaoNeural` at -2%**, which is a warmer
adult female voice with much smoother sentence prosody. A local **CosyVoice2**
backend is also implemented for a less synthetic timbre; it is optional because
it needs a free CUDA GPU and external model weights.

## Current cut: Edge TTS Xiaoxiao

- **Why:** immediate, reliable, no GPU; Xiaoxiao is markedly more natural than
  Xiaoyi for long Chinese paragraphs, without the cartoon cadence.
- **How:** `scripts/generate_tts_xiaoxiao.py` writes
  `public/narration/*.mp3` and `src/timing.json`.
- **Run:** `python -m pip install -r requirements-tts.txt` and
  `python scripts/generate_tts_xiaoxiao.py`.
- **Cut length:** 7:14 (433.6 s).

## Optional: CosyVoice2-0.5B

- **Why:** natural sentence prosody, 24 kHz output, zero-shot voice from a
  short reference clip, Apache-2.0. This is the best-sounding option tested.
- **How:** `scripts/generate_tts_cosyvoice.py` uses the female reference clip
  shipped with CosyVoice in zero-shot mode. It synthesizes each scene, measures
  the MP3 with `ffprobe`, and rewrites `src/timing.json`.
- **Cost:** model loading takes about two minutes; the first chunk pays a
  one-time CUDA warm-up. It needs roughly 10 GB of VRAM in fp16.
- **Status:** the pipeline works and was used to generate test audio, but the
  full run was not committed. On the author’s machine ComfyUI/Wan2.2 was
  occupying the GPU, and long Chinese paragraphs occasionally stalled. If the
  GPU is free, run it and re-render; the voice is noticeably less synthetic.
- **Run:** see the path variables at the top of `generate_tts_cosyvoice.py`.

## Alternatives considered

| Option | Strength | Why it was not used for this cut |
|---|---|---|
| **Edge TTS / Xiaoyi** | Free, instant, no GPU | The original voice; its clipped, even rhythm was the problem |
| **Edge TTS / Xiaoxiao** | Free, warm adult voice, reliable | Still a cloud Edge voice, not a local neural model |
| **CosyVoice2** | Natural local voice, zero-shot, 24 kHz | Needs GPU; long-run stability depends on free VRAM |
| **IndexTTS2** | Excellent emotion and duration control; strong Chinese naturalness | 1.5B model, ~8 GB VRAM, reference audio required; heavier setup |
| **Spark-TTS** | 0.5B, lightweight, 1–2 GB VRAM | 16 kHz output and less expressive long-form prosody |
| **Kokoro-82M** | Very fast, CPU-friendly | Chinese support is more limited; less expressive for a long narrative |
| **MiniMax Speech / ElevenLabs / Azure styles** | Often the most polished commercial voices; style control | Requires an API key, per-character cost, and sending text to a third party |
| **F5-TTS / GPT-SoVITS / Fish Speech** | Strong zero-shot cloning | Needs a clean reference recording and more tuning; not necessary for one narrator |

## If the voice still sounds uneven

1. Change `RATE` in `generate_tts_xiaoxiao.py` from `-2%` to `+0%`; the text is
   already written with short clauses, so the default rate may flow better.
2. For CosyVoice2, use a calm, warm 5–10 s female reference clip. The model
   clones timbre and speaking style from the prompt; an animated prompt gives
   an animated narrator.
3. Keep the text but shorten the first chunk. CosyVoice2 splits at about 80
   tokens; a long first sentence pays most of the warm-up and can sound
   different from the rest.
4. Apply a very mild `atempo=1.05–1.10` in post if a cut exceeds eight minutes;
   this changes speed without changing pitch.

## Reproducibility

`public/narration/*.mp3` and `src/timing.json` are committed, so the video can
be rendered without downloading a TTS model. Regenerating the voice is optional.
