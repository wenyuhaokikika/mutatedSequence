#!/usr/bin/env python3
"""Generate the narration and timing for the Remotion explainer.

Voice: zh-CN-XiaoyiNeural (female, lively) at -4% rate and +2 Hz pitch.
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

ROOT = Path(__file__).resolve().parent.parent
AUDIO_DIR = ROOT / "public" / "narration"
TIMING_PATH = ROOT / "src" / "timing.json"

VOICE = "zh-CN-XiaoyiNeural"
RATE = "-4%"
PITCH = "+2Hz"
FPS = 30
TAIL_SECONDS = 1.0

SCENES = [
    {
        "id": "hook",
        "title": "一个字母，一整条蛋白",
        "narration": "如果把 VCF 里的一个小字母，变成一整条蛋白，会发生什么？",
    },
    {
        "id": "problem",
        "title": "只按外显子翻译，会迷路",
        "narration": "以前的流程，逐个外显子翻译。可一遇到负链、跨外显子，或者结构变异，它就容易迷路。",
    },
    {
        "id": "inputs",
        "title": "三样东西，就够了",
        "narration": "现在，只要三样东西：参考基因组、转录本注释，还有你的 VCF。",
    },
    {
        "id": "ref",
        "title": "先对暗号",
        "narration": "第一步，先对暗号：VCF 的 REF，和参考基因组对得上吗？对不上，就温柔提醒；用了 strict，它才会立刻喊停。",
    },
    {
        "id": "allele",
        "title": "选一个等位基因",
        "narration": "第二步，选等位基因。零比一、一比一、一比二，还有 phased 的基因型；一个共识只选一个 ALT，别让所有等位基因挤在一起。",
    },
    {
        "id": "blocks",
        "title": "把变异画成积木",
        "narration": "第三步，把变异画成基因组上的积木。参考块、替换块、插入块；就算 deletion 横跨内含子，也能准确切回外显子。",
    },
    {
        "id": "protein",
        "title": "从外显子到蛋白",
        "narration": "第四步，沿着注释切出 cDNA、CDS，再翻译成蛋白。五撇端的 N 补齐、终止密码子、三撇端延伸，都会分开处理。",
    },
    {
        "id": "outputs",
        "title": "每一步都能审计",
        "narration": "第五步，写出可审计的 TSV：原始 REF 和 ALT、选中的等位基因、突变前后的序列，一条都不少。",
    },
    {
        "id": "limits",
        "title": "它也会诚实说不",
        "narration": "当然，它也有做不到的事：真正的 phasing、剪接跳跃、TRA、BND、CNV。遇到这些，它会诚实地说：我处理不了。",
    },
    {
        "id": "outro",
        "title": "mutatedSequence",
        "narration": "把 VCF 变成完整、可解释的蛋白。mutatedSequence，开源，Docker 一键运行。",
    },
]


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
