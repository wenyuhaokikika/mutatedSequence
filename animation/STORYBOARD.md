# 分镜与节奏｜mutatedSequence 流程动画

## 一句话

把“VCF 里一个字母，如何变成一条完整、可审计的蛋白”讲成一支 **1 分 52 秒** 的温柔科普短片；不堆术语，先把流程讲顺，再把边界讲清楚。

## 片子规格

| 项目 | 设定 |
|---|---|
| 片名 | 把 VCF 变成一条蛋白 |
| 时长 | 112.2 秒（3,366 帧） |
| 画幅 | 1920 × 1080，16:9 |
| 帧率 | 30 fps |
| 语言 | 普通话，技术名词保留英文 |
| 旁白 | `zh-CN-XiaoyiNeural`，女声，lively / 俏皮 |
| 语速 / 音高 | `rate -4%`，`pitch +2Hz`；温柔一点，但不拖 |
| 字幕 | 每场完整旁白，底部奶油色胶囊 |
| 配乐 | v1 不加 BGM，保留人声呼吸感；如果要加，建议 -22 LUFS 以下 |

## 视觉方向

- **底色**：奶油纸感 `#FFF8F0`，像一本摊开的实验记录本。
- **主色**：珊瑚 `#FF6F61`（变异/重点）、青绿 `#2A9D8F`（reference/通过）、薰衣草 `#CDB4DB`（边界/未建模）、金色 `#F4C95D`（CDS/stop）。
- **形状**：大圆角卡片、柔和投影、细网格、漂浮的 A/C/G/T、少量星星。
- **字体**：Noto Sans CJK SC；代码/序列用 JetBrains Mono / Consolas。
- **运动**：所有元素 8–14 帧 fade + 20px 上浮；spring 阻尼 14–18，不要弹得过火。每 1.5–2 秒给一个视觉变化点。
- **镜头感**：不追求 3D；用“信息卡片 + 轻漂浮”的扁平动画，像在纸上推演。
- **安全区**：左右 100px；标题下方 142px 开始内容；字幕底部 42px；进度条 7px。

## 叙事节奏

```text
0–26s    钩子 + 旧流程痛点 + 输入      快 → 慢
26–52s   REF 校验 + 等位基因           稳，两个“关卡”
52–78s   genomic blocks + 中心法则     全片信息密度最高
78–103s  输出审计 + 诚实边界           放慢，留白
103–112s 收束 + CTA                    暖、短、干净
```

节奏原则：

1. **先给结果，再给机制**：开场 3 秒就出现“小字母 → 一整条蛋白”。
2. **每场只讲一个动作**：REF、allele、blocks、protein、TSV 各自独立。
3. **关键转折留半拍**：`--strict`、unresolved、TRA/BND/CNV 出现时，动画放慢 6–8 帧。
4. **不把限制藏到最后**：在“看起来很顺”之后立刻讲“它也会说不”，可信度更高。
5. **结尾不喊口号**：用 GitHub / Docker 两个胶囊落地。

## 分镜表

| # | 时间 | 场景 | 画面 | 旁白 | 动效 / 节奏 |
|---:|---|---|---|---|---|
| 01 | 00:00–00:07.4 | Hook | 奶油底；左边大字“把小字母变成一整条蛋白”；右边 VCF 卡 + 蛋白珠链 | 如果把 VCF 里的一个小字母，变成一整条蛋白，会发生什么？ | 标题 2 帧入场；VCF 卡 20f；珊瑚箭头 38f；蛋白珠 34f 逐个弹入；金色下划线 18–42f 拉开 |
| 02 | 00:07.4–00:17.6 | Problem | 三张“迷路卡”：负链 / 跨外显子 / 结构变异；虚线路径；右上“REF 对不上？”气泡 | 以前的流程，逐个外显子翻译。可一遇到负链、跨外显子，或者结构变异，它就容易迷路。 | 卡片 12/22/32f 错峰；虚线 20–60f 从左到右；气泡 40f 浮入 |
| 03 | 00:17.6–00:26.1 | Inputs | genome.fa / annotation.gtf / sample.vcf 三张文件卡；三支箭头汇入 mutatedSequence 引擎卡 | 现在，只要三样东西：参考基因组、转录本注释，还有你的 VCF。 | 文件卡 10/18/26f；三支箭头 28f 同时向下；引擎卡 24f 缩放进入；optional 标签 36f |
| 04 | 00:26.1–00:38.6 | REF check | chr1 染色体条 + 变异红点；VCF 卡 REF C / genome A；右侧 green match + red mismatch + `--strict` 闸门 | 第一步，先对暗号：VCF 的 REF，和参考基因组对得上吗？对不上，就温柔提醒；用了 strict，它才会立刻喊停。 | 染色体 12f；VCF 卡 16f；绿/红两行 24/32f 错峰；`--strict` 46f 落下，exit 2 胶囊脉冲 |
| 05 | 00:38.6–00:52.1 | Allele | 四枚 GT 胶囊（0/1、1/1、1/2、1\|0）；选中 0/1；右侧 ALT1 selected / ALT2 not applied；底部“所有 ALT 一起上？不行” | 第二步，选等位基因。零比一、一比一、一比二，还有 phased 的基因型；一个共识只选一个 ALT，别让所有等位基因挤在一起。 | GT 12/18/24/30f 从左滑入；choose 箭头 30f；ALT 卡 34/42f；红色叉号 52f 轻弹 |
| 06 | 00:52.1–01:04.9 | Blocks | 基因组带：exon 1/2/3；红色虚线 bracket 跨 intron；INS 三角；下方 reference / ALT / insertion 三种 block 芯片；提取卡“exon1 + exon2 → mutated cDNA” | 第三步，把变异画成基因组上的积木。参考块、替换块、插入块；就算 deletion 横跨内含子，也能准确切回外显子。 | 基因组线 8–38f 拉开；exon 12/16/20f；bracket 28f 下压；block 芯片 40f；提取卡 58f；箭头 64f |
| 07 | 01:04.9–01:18.1 | Protein | 四栏：exons → cDNA → CDS → protein；序列条；N pad / stop 标签；annotated + 3′ extended 两条蛋白珠链 | 第四步，沿着注释切出 cDNA、CDS，再翻译成蛋白。五撇端的 N 补齐、终止密码子、三撇端延伸，都会分开处理。 | 四栏 10/18/26/34f 错峰；序列条逐行；蛋白珠 38/44f 两轮；N pad/stop 标签 30f |
| 08 | 01:18.1–01:29.6 | Outputs | mutation_report.tsv 表格：variant / REF-ALT / allele / protein；下方 before → after 序列条；右侧 FASTA / normalized.vcf / run_summary；AUDITABLE 印章 | 第五步，写出可审计的 TSV：原始 REF 和 ALT、选中的等位基因、突变前后的序列，一条都不少。 | 表格 10f；行 14–26f 错峰；before/after 38/44f；右侧卡 32f；印章 52f 旋转盖下 |
| 09 | 01:29.6–01:42.6 | Limits | 三张 “not modeled” 卡：phasing / splicing / TRA-BND-CNV；下方 supported DEL/INS/DUP/INV 芯片；右上诚实气泡 | 当然，它也有做不到的事：真正的 phasing、剪接跳跃、TRA、BND、CNV。遇到这些，它会诚实地说：我处理不了。 | 卡片 12/20/28f；气泡 40f；supported 芯片 52f；语气放慢，画面不抖 |
| 10 | 01:42.6–01:52.2 | Outro | 中心 DNA 图标 + “把 VCF 变成完整、可解释的蛋白”；mutatedSequence 字标；GitHub / Docker Hub 胶囊；open source / 10 tests / chr21 validated | 把 VCF 变成完整、可解释的蛋白。mutatedSequence，开源，Docker 一键运行。 | 图标 2f；标题 12f；链接 30f；金色线 10–44f；星星轻闪；最后 8f 淡出 |

## 旁白与 TTS

- 生成脚本：`scripts/generate_tts.py`
- 音色：`zh-CN-XiaoyiNeural`
- 参数：`--rate=-4% --pitch=+2Hz`
- 每场尾部留 **1.0 秒** 呼吸；`timing.json` 按实际 mp3 时长计算画面帧数。
- 字幕直接使用旁白文本；不做逐字高亮（v1 以完整语义为主）。
- 所有技术词保留英文：`VCF`、`REF`、`ALT`、`allele`、`consensus`、`phasing`、`splicing`、`TRA/BND/CNV`。

## 制作与渲染

```bash
cd animation
python scripts/generate_tts.py      # 重新生成旁白 + timing.json
npm install
npm start                           # Remotion Studio 预览
./render.sh                         # Docker 渲染 1080p MP4
```

渲染环境包含 Chromium Headless Shell 所需库和 Noto CJK 字体；视频输出在 `out/explainer.mp4`。

## 刻意不做的

- 不做 3D DNA 旋转：信息密度已经够高，3D 会抢旁白。
- 不做逐字卡拉 OK 字幕：中文断词容易碎，改用整句字幕。
- 不把全部限制塞进画面：只保留 phasing / splicing / TRA-BND-CNV 三类，细节留给 README。
- 不用悲伤/警告红铺满画面：这是“温柔而诚实”，不是故障演示。
