# mutatedSequence · Remotion explainer

A **1:52**, 1920×1080 / 30 fps explainer for the VCF → transcript → protein
pipeline. The narration is Chinese, performed by the playful-but-gentle
`zh-CN-XiaoyiNeural` voice.

- **Storyboard / direction:** [STORYBOARD.md](STORYBOARD.md)
- **Composition:** `src/index.ts` → `Explainer`
- **Narration + timing:** `src/timing.json` (generated from the actual MP3s)
- **Rendered output:** `out/explainer.mp4`

## The cut

| # | Scene | Time |
|---:|---|---|
| 01 | Hook — one letter, one protein | 0:00–0:07 |
| 02 | Old per-exon pipeline gets lost | 0:07–0:18 |
| 03 | Three inputs | 0:18–0:26 |
| 04 | REF check + `--strict` | 0:26–0:39 |
| 05 | Choose one allele | 0:39–0:52 |
| 06 | Genomic blocks | 0:52–1:05 |
| 07 | cDNA → CDS → protein | 1:05–1:18 |
| 08 | Auditable TSV / FASTA | 1:18–1:30 |
| 09 | Honest limits | 1:30–1:43 |
| 10 | Outro | 1:43–1:52 |

## Voice

`scripts/generate_tts.py` uses:

```text
voice: zh-CN-XiaoyiNeural
rate:  -4%
pitch: +2Hz
tail:  1.0s per scene
```

Xiaoyi is the lively/cartoon female voice in Edge TTS; the slightly slower rate
and small pitch lift keep it playful without sounding rushed. The script calls
`ffprobe` to measure each MP3 and writes `src/timing.json`, so the edit is cut to
the voice rather than to guessed durations.

Regenerate the audio (requires `ffmpeg` + `edge-tts`):

```bash
python -m pip install -r requirements-tts.txt
python scripts/generate_tts.py
```

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
