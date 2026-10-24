#!/usr/bin/env python3
"""Generate the explainer narration with CosyVoice2 (local, GPU).

This is the higher-quality alternative to the Edge TTS backend. It uses the
Apache-2.0 CosyVoice2-0.5B model in zero-shot mode with the female reference
clip shipped with the CosyVoice repository. Zero-shot mode was more stable than
the instruct mode on long Chinese paragraphs: the instruct prompt occasionally
caused a chunk to loop, while the same text completed in seconds in zero-shot.

Prerequisites (outside this repository):
  * Python 3.10 environment with torch + torchaudio (CUDA build)
  * CosyVoice source at COSYVOICE_ROOT (default D:/workplace/_tts/CosyVoice-main)
  * Model weights at COSYVOICE_MODEL (default D:/workplace/_tts/models/CosyVoice2-0.5B)

Run:
  python scripts/generate_tts_cosyvoice.py
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Windows does not always add torch/lib to the DLL search path when a conda
# environment is invoked by absolute path. Do it before importing torch.
if os.name == "nt":
    _torch_lib = Path(sys.prefix) / "Lib" / "site-packages" / "torch" / "lib"
    if _torch_lib.exists():
        os.add_dll_directory(str(_torch_lib))
        os.environ["PATH"] = str(_torch_lib) + os.pathsep + os.environ.get("PATH", "")

import torch
import torchaudio

from narration import FPS, SCENES, TAIL_SECONDS

ROOT = Path(__file__).resolve().parent.parent
AUDIO_DIR = ROOT / "public" / "narration"
TIMING_PATH = ROOT / "src" / "timing.json"

COSY_ROOT = Path(os.environ.get("COSYVOICE_ROOT", r"D:\workplace\_tts\CosyVoice-main"))
MODEL_DIR = Path(os.environ.get("COSYVOICE_MODEL", r"D:\workplace\_tts\models\CosyVoice2-0.5B"))
PROMPT_WAV = COSY_ROOT / "asset" / "zero_shot_prompt.wav"
PROMPT_TEXT = "希望你以后能够做的比我还好呦。"
SPEED = 0.95
VOICE = "CosyVoice2-0.5B / zero_shot_prompt"


def audio_duration(path: Path) -> float:
    if shutil.which("ffprobe") is None:
        raise RuntimeError("ffprobe was not found; install ffmpeg")
    out = subprocess.check_output(
        [
            "ffprobe",
            "-v",
            "error",
            "-show_entries",
            "format=duration",
            "-of",
            "default=nw=1:nk=1",
            str(path),
        ],
        text=True,
    ).strip()
    return float(out)


def main() -> None:
    if not COSY_ROOT.exists():
        raise SystemExit(f"CosyVoice source not found: {COSY_ROOT}")
    if not MODEL_DIR.exists():
        raise SystemExit(f"CosyVoice model not found: {MODEL_DIR}")
    if not PROMPT_WAV.exists():
        raise SystemExit(f"Reference audio not found: {PROMPT_WAV}")

    sys.path.insert(0, str(COSY_ROOT / "third_party" / "Matcha-TTS"))
    sys.path.insert(0, str(COSY_ROOT))
    from cosyvoice.cli.cosyvoice import AutoModel

    print(f"loading {MODEL_DIR}", flush=True)
    model = AutoModel(model_dir=str(MODEL_DIR), fp16=torch.cuda.is_available())
    AUDIO_DIR.mkdir(parents=True, exist_ok=True)

    scenes = []
    total_frames = 0
    for scene in SCENES:
        wav_path = AUDIO_DIR / f"{scene['id']}.wav"
        mp3_path = AUDIO_DIR / f"{scene['id']}.mp3"
        pieces = []
        print(f"  {scene['id']}: {len(scene['narration'])} chars", flush=True)
        for chunk in model.inference_zero_shot(
            scene["narration"],
            PROMPT_TEXT,
            str(PROMPT_WAV),
            stream=False,
            speed=SPEED,
        ):
            speech = chunk["tts_speech"].detach().cpu()
            pieces.append(speech)
            print(f"    chunk {speech.shape[-1] / model.sample_rate:.2f}s", flush=True)
        if not pieces:
            raise RuntimeError(f"no audio produced for {scene['id']}")
        speech = torch.cat(pieces, dim=1)
        torchaudio.save(str(wav_path), speech, model.sample_rate)
        subprocess.run(
            [
                "ffmpeg",
                "-v",
                "error",
                "-y",
                "-i",
                str(wav_path),
                "-codec:a",
                "libmp3lame",
                "-b:a",
                "192k",
                str(mp3_path),
            ],
            check=True,
        )
        wav_path.unlink(missing_ok=True)
        seconds = audio_duration(mp3_path)
        frames = math.ceil((seconds + TAIL_SECONDS) * FPS)
        scenes.append(
            {
                **scene,
                "audio": f"narration/{scene['id']}.mp3",
                "seconds": round(seconds, 3),
                "durationInFrames": frames,
                "from": total_frames,
            }
        )
        total_frames += frames
        print(f"  {scene['id']}: {seconds:.2f}s -> {frames} frames", flush=True)

    timing = {
        "fps": FPS,
        "width": 1920,
        "height": 1080,
        "voice": VOICE,
        "rate": f"speed={SPEED}",
        "pitch": "natural",
        "backend": "CosyVoice2-0.5B",
        "reference": "asset/zero_shot_prompt.wav",
        "mode": "zero_shot",
        "totalFrames": total_frames,
        "totalSeconds": round(total_frames / FPS, 2),
        "scenes": scenes,
    }
    TIMING_PATH.parent.mkdir(parents=True, exist_ok=True)
    TIMING_PATH.write_text(
        json.dumps(timing, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"total: {total_frames} frames / {total_frames / FPS:.1f}s", flush=True)
    print(f"wrote {TIMING_PATH}", flush=True)


if __name__ == "__main__":
    main()
