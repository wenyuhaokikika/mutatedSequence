# mutatedSequence · Remotion explainer

A **7:14**, 1920×1080 / 30 fps explainer for the VCF → transcript → protein
pipeline. The narration is Chinese, performed by the warmer, more natural
`zh-CN-XiaoxiaoNeural` voice.

- **Storyboard / direction:** [STORYBOARD.md](STORYBOARD.md)
- **Composition:** `src/index.ts` → `Explainer`
- **Narration + timing:** `src/timing.json` (generated from the actual MP3s)
- **Rendered output:** `out/explainer.mp4`
- **Online cut:** https://wenyuhaokikika.github.io/2026/10/24/mutatedsequence-vcf-to-transcript-protein/

## The cut

| # | Scene | Time |
|---:|---|---|
| 01 | Why — WGS data + AI structures, and the missing bridge | 0:00–0:54 |
| 02 | Map — genome → transcript → protein, in plain language | 0:54–1:29 |
| 03 | Three inputs and the outputs | 1:29–2:04 |
| 04 | Check the address (REF) | 2:04–2:31 |
| 05 | Two copies, choose one version | 2:31–3:04 |
| 06 | Edit genomic blocks first | 3:04–3:31 |
| 07 | cDNA → CDS → protein | 3:31–3:59 |
| 08 | Stop-codon loss and readthrough | 3:59–4:41 |
| 09 | Auditable TSV / FASTA | 4:41–5:07 |
| 10 | Docker in one command | 5:07–5:26 |
| 11 | Where it fits among VEP/Haplosaurus/VCF2Prot/ORFan | 5:26–6:17 |
| 12 | Honest limits | 6:17–6:56 |
| 13 | Outro | 6:56–7:14 |

## Voice

`scripts/generate_tts_xiaoxiao.py` uses:

```text
voice: zh-CN-XiaoxiaoNeural
rate:  -2%
pitch: +0Hz
tail:  0.9s per scene
```

Xiaoxiao is an adult female voice with much smoother sentence prosody than
Xiaoyi, which sounded clipped and metronomic in the first cut. The script calls
`ffprobe` to measure each MP3 and writes `src/timing.json`, so the edit is cut to
the voice rather than to guessed durations. The narration is split into 2–5
short caption chunks per scene.

Regenerate the current audio (requires `ffmpeg` + `edge-tts`):

```bash
python -m pip install -r requirements-tts.txt
python scripts/generate_tts_xiaoxiao.py
```

A local CosyVoice2 backend is also available and produces a more natural,
less synthetic timbre. It needs a CUDA GPU and the external CosyVoice checkout
and weights; see [TTS_OPTIONS.md](TTS_OPTIONS.md) and
`scripts/generate_tts_cosyvoice.py`. The older Xiaoyi backend is retained as
`scripts/generate_tts.py`.

The generated MP3s live in `public/narration/` and are committed, so a normal
render does not need network access or TTS credentials.

## Preview

```bash
npm install
npm start
```

The Remotion Studio opens the composition. If you only need stills:

```bash
npx remotion still src/index.ts Explainer out/still.png --frame=60
```

## Render (recommended: Docker)

```bash
./render.sh
# or
./render.sh out/my-cut.mp4
```

`render.sh` builds `Dockerfile.render` (Node 20 + Chromium Headless Shell
libraries + Noto CJK) and runs Remotion in the container. The project is
bind-mounted; no source is baked into the render image.

On Git Bash the script sets `MSYS_NO_PATHCONV=1` and uses `pwd -W` for the
volume path. On Linux, outputs are owned by the container's root user; add
`--user "$(id -u):$(id -g)"` to the `docker run` line if you want host-owned
files.

Direct local render (if you already have a Chrome/Chromium executable):

```bash
npx remotion render src/index.ts Explainer out/explainer.mp4
```

## Project layout

```text
animation/
├── src/
│   ├── index.ts              # registerRoot
│   ├── Root.tsx              # Composition metadata
│   ├── Explainer.tsx         # Sequences + audio
│   ├── timing.json           # generated voice timings
│   ├── theme.ts              # palette / type / card styles
│   ├── types.ts
│   ├── components/
│   │   ├── ui.tsx            # SceneShell, backdrop, pills, animation helpers
│   │   └── visuals.tsx       # VCF, chromosome, protein, sequence visuals
│   └── scenes/               # 10 scene components
├── scripts/generate_tts.py
├── public/narration/*.mp3
├── Dockerfile.render
├── render.sh
└── STORYBOARD.md
```

## Editing the story

1. Change the `narration` in `scripts/generate_tts.py`.
2. Run `python scripts/generate_tts.py` — this rewrites the MP3s and
   `src/timing.json`.
3. Adjust the matching scene component in `src/scenes/`.
4. Preview, then render.

Keep the narration lines short. The caption box is designed for **one or two
lines at 1920px**; very long lines will crowd the lower safe area. Technical
terms (`VCF`, `REF`, `ALT`, `allele`, `phasing`, `TRA/BND/CNV`) are intentionally
kept in English and should match the on-screen copy.

## Notes

- The video uses a cream paper, coral/teal/lavender palette and vector/CSS
  only—no stock footage, no external images.
- No BGM in v1; the voice and the small motion beats are the rhythm section.
- `out/` and `node_modules/` are gitignored. The narration MP3s are kept so the
  project is renderable without regenerating TTS.
