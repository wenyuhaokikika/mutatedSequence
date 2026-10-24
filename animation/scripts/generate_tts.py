#!/usr/bin/env python3
"""Generate the narration and timing for the Remotion explainer.

Voice: zh-CN-XiaoyiNeural (female, lively) at -6% rate and +2 Hz pitch.
The result is intentionally a little playful and gentle, not a newsreader.
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

VOICE = "zh-CN-XiaoyiNeural"
RATE = "-6%"
PITCH = "+2Hz"
def audio_duration(path: Path) -> float:
    if shutil.which("ffprobe") is None:
        raise RuntimeError("ffprobe was not found; install ffmpeg or edit the timings manually")
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
        scene["narration"],
        VOICE,
        rate=RATE,
        pitch=PITCH,
    )
    await communicate.save(str(out))
    print(f"  {scene['id']}: {out.name}")


async def main() -> None:
    AUDIO_DIR.mkdir(parents=True, exist_ok=True)
    print(f"voice={VOICE} rate={RATE} pitch={PITCH}")
    for scene in SCENES:
        await synthesize(scene)

    scenes = []
    total_frames = 0
    for scene in SCENES:
        path = AUDIO_DIR / f"{scene['id']}.mp3"
        seconds = audio_duration(path)
        # Keep a short breathing tail after each line so the edit does not feel
        # rushed, and so the next scene can fade in before the voice starts.
        frames = math.ceil((seconds + TAIL_SECONDS) * FPS)
        start = total_frames
        scenes.append(
            {
                **scene,
                "audio": f"narration/{scene['id']}.mp3",
                "seconds": round(seconds, 3),
                "durationInFrames": frames,
                "from": start,
            }
        )
        total_frames += frames
        print(f"  {scene['id']}: {seconds:.2f}s -> {frames} frames")

    timing = {
        "fps": FPS,
        "width": 1920,
        "height": 1080,
        "voice": VOICE,
        "rate": RATE,
        "pitch": PITCH,
        "totalFrames": total_frames,
        "totalSeconds": round(total_frames / FPS, 2),
        "scenes": scenes,
    }
    TIMING_PATH.parent.mkdir(parents=True, exist_ok=True)
    TIMING_PATH.write_text(
        json.dumps(timing, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"total: {total_frames} frames / {total_frames / FPS:.1f}s")
    print(f"wrote {TIMING_PATH}")


if __name__ == "__main__":
    asyncio.run(main())
