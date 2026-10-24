#!/usr/bin/env python3
"""Fallback narration using Edge TTS Xiaoxiao (fast, CPU-only).

This is the committed narration backend: it needs no GPU and overwrites
public/narration/*.mp3 and src/timing.json.
"""

from __future__ import annotations

import asyncio
import json
import math
import shutil
import subprocess
from pathlib import Path

import edge_tts

from narration import FPS, SCENES, TAIL_SECONDS

ROOT = Path(__file__).resolve().parent.parent
AUDIO_DIR = ROOT / "public" / "narration"
TIMING_PATH = ROOT / "src" / "timing.json"

VOICE = "zh-CN-XiaoxiaoNeural"
RATE = "-2%"
PITCH = "+0Hz"


def audio_duration(path: Path) -> float:
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


async def synthesize(scene: dict) -> None:
    out = AUDIO_DIR / f"{scene['id']}.mp3"
    communicate = edge_tts.Communicate(
        scene["narration"], VOICE, rate=RATE, pitch=PITCH
    )
    await communicate.save(str(out))


async def main() -> None:
    AUDIO_DIR.mkdir(parents=True, exist_ok=True)
    for scene in SCENES:
        await synthesize(scene)
        print(f"  {scene['id']}: {audio_duration(AUDIO_DIR / (scene['id'] + '.mp3')):.2f}s")

    scenes = []
    total_frames = 0
    for scene in SCENES:
        seconds = audio_duration(AUDIO_DIR / f"{scene['id']}.mp3")
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

    timing = {
        "fps": FPS,
        "width": 1920,
        "height": 1080,
        "voice": VOICE,
        "rate": RATE,
        "pitch": PITCH,
        "backend": "edge-tts",
        "totalFrames": total_frames,
        "totalSeconds": round(total_frames / FPS, 2),
        "scenes": scenes,
    }
    TIMING_PATH.write_text(
        json.dumps(timing, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"total: {total_frames} frames / {total_frames / FPS:.1f}s")


if __name__ == "__main__":
    asyncio.run(main())
