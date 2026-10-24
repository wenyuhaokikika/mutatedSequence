import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, Sparkle, fadeUp, useEnter} from '../components/ui';
import {ArrowLine} from '../components/visuals';
import type {SceneProps} from './types';

const ReadBars: React.FC<{delay: number}> = ({delay}) => {
  const p = useEnter(delay);
  const rows = [
    'ACGTTGCAAGCTTAGCGATCG',
    'TTGACCGTAGCTAGCTAACGT',
    'GGCTAACGTTAGCCGATCGTA',
    'CATGCGTAGCTAGCATCGATG',
    'AACGTTGCAAGCTTAGCGATC',
  ];
  return (
    <div style={{...card, padding: '16px 18px', width: 290, ...fadeUp(p, 18)}}>
      <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: 12}}>
        <span style={{fontFamily: MONO, fontSize: 12, fontWeight: 800}}>WGS reads</span>
        <span style={{fontSize: 11, color: COLORS.muted}}>sequencing</span>
      </div>
      {rows.map((r, i) => (
        <div
          key={r}
          style={{
            display: 'flex',
            gap: 2,
            marginBottom: 5,
            opacity: 0.4 + i * 0.1,
          }}
        >
          {r.split('').map((b, j) => (
            <span
              key={j}
              style={{
                width: 8,
                height: 15,
                borderRadius: 2,
                background:
                  b === 'A'
                    ? COLORS.coral
                    : b === 'C'
                      ? COLORS.teal
                      : b === 'G'
                        ? COLORS.gold
                        : COLORS.lavender,
                opacity: 0.55,
              }}
            />
          ))}
        </div>
      ))}
    </div>
  );
};

const VcfMini: React.FC<{delay: number}> = ({delay}) => {
  const p = useEnter(delay);
  const rows = [
    ['1', '12345', 'C', 'T'],
    ['1', '23456', 'A', 'G'],
    ['2', '88102', 'G', 'A'],
    ['2', '99120', 'T', 'C'],
  ];
  return (
    <div style={{...card, padding: '16px 18px', width: 300, ...fadeUp(p, 18)}}>
      <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: 10}}>
        <span style={{fontFamily: MONO, fontSize: 12, fontWeight: 800}}>sample.vcf</span>
        <span style={{fontSize: 11, color: COLORS.muted}}>mutations</span>
      </div>
      {rows.map((r, i) => (
        <div
          key={i}
          style={{
            display: 'grid',
            gridTemplateColumns: '22px 60px 22px 22px',
            gap: 7,
            padding: '5px 0',
            borderTop: i ? `1px dashed ${COLORS.line}` : 'none',
            fontFamily: MONO,
            fontSize: 11,
            color: i === 0 ? COLORS.coral : COLORS.ink,
            fontWeight: i === 0 ? 800 : 600,
          }}
        >
          <span>{r[0]}</span>
          <span>{r[1]}</span>
          <span>{r[2]}</span>
          <span style={{color: COLORS.teal}}>{r[3]}</span>
        </div>
      ))}
    </div>
  );
};

const Structure: React.FC<{delay: number}> = ({delay}) => {
  const p = useEnter(delay, 14);
  const frame = useCurrentFrame();
  const glow = 0.5 + Math.sin(frame / 14) * 0.08;
  return (
    <div
      style={{
        ...card,
        width: 350,
        height: 300,
        padding: '18px 20px',
        background:
          'linear-gradient(145deg, rgba(242,234,248,0.95), rgba(221,243,238,0.7))',
        ...fadeUp(p, 20),
        position: 'relative',
        overflow: 'hidden',
      }}
    >
      <div style={{display: 'flex', justifyContent: 'space-between'}}>
        <span style={{fontFamily: MONO, fontSize: 12, fontWeight: 800}}>protein structure</span>
        <span style={{fontSize: 11, color: COLORS.muted}}>AI</span>
      </div>
      <svg viewBox="0 0 310 230" style={{width: '100%', marginTop: 8}}>
        <defs>
          <filter id="glowWhy">
            <feGaussianBlur stdDeviation="7" result="blur" />
            <feMerge>
              <feMergeNode in="blur" />
              <feMergeNode in="SourceGraphic" />
            </feMerge>
          </filter>
        </defs>
        <path
          d="M45 150 C80 40, 120 45, 140 120 S195 210, 225 110 S270 55, 285 135"
          stroke={COLORS.teal}
          strokeWidth={10}
          fill="none"
          strokeLinecap="round"
          opacity={glow}
          filter="url(#glowWhy)"
        />
        <path
          d="M60 75 C100 160, 150 190, 180 100 S235 20, 270 75"
          stroke={COLORS.lavender}
          strokeWidth={7}
          fill="none"
          strokeLinecap="round"
          opacity={0.85}
        />
        <path
          d="M90 185 C130 110, 175 105, 205 175 S250 205, 285 190"
          stroke={COLORS.gold}
          strokeWidth={6}
          fill="none"
          strokeLinecap="round"
          opacity={0.8}
        />
        {[
          [45, 150],
          [140, 120],
          [225, 110],
          [285, 135],
          [100, 175],
          [205, 175],
        ].map(([x, y], i) => (
          <circle key={i} cx={x} cy={y} r={6} fill={COLORS.paper} stroke={COLORS.coral} strokeWidth={3} />
        ))}
      </svg>
      <div style={{display: 'flex', gap: 6, flexWrap: 'wrap', marginTop: 2}}>
        {['AlphaFold3', 'RoseTTAFold', 'ESM'].map((x) => (
          <span
            key={x}
            style={{
              fontFamily: MONO,
              fontSize: 10,
              fontWeight: 700,
              color: COLORS.teal,
              background: 'rgba(255,255,255,0.7)',
              border: `1px solid ${COLORS.teal}33`,
              borderRadius: 6,
              padding: '3px 6px',
            }}
          >
            {x}
          </span>
        ))}
      </div>
    </div>
  );
};

export const WhyScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const title = useEnter(2, 18);
  const reads = useEnter(10);
  const vcf = useEnter(20);
  const gap = useEnter(34);
  const structure = useEnter(48);
  const phase = Math.min(1, frame / Math.max(1, scene.durationInFrames));
  const gapPulse = gap * (0.8 + Math.sin(frame / 12) * 0.2);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.coral}>
      <AbsoluteFill>
        <div style={{display: 'flex', gap: 44, height: '100%', alignItems: 'center'}}>
          <div style={{width: 680, ...fadeUp(title, 22)}}>
            <Pill color={COLORS.coral} bg={COLORS.coralSoft} style={{marginBottom: 22}}>
              why this matters
            </Pill>
            <h2
              style={{
                margin: 0,
                fontFamily: FONT,
                fontSize: 66,
                lineHeight: 1.12,
                letterSpacing: -3,
                fontWeight: 900,
                color: COLORS.ink,
              }}
            >
              我们有海量突变，
              <br />
              也有了蛋白 AI，
              <br />
              <span style={{color: COLORS.coral}}>但中间缺一座桥</span>
            </h2>
            <p style={{fontSize: 20, lineHeight: 1.7, color: COLORS.muted, margin: '24px 0 0'}}>
              测序告诉我们「哪里变了」；
              <br />
              结构模型想要的是「一条干净的蛋白序列」。
            </p>
            <div style={{display: 'flex', gap: 10, marginTop: 22}}>
              <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.teal}}>WGS → VCF</span>
              <span style={{color: COLORS.muted}}>·</span>
              <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.coral}}>VCF → protein</span>
              <span style={{color: COLORS.muted}}>·</span>
              <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.green}}>protein → structure</span>
            </div>
          </div>

          <div style={{flex: 1, position: 'relative', height: 520}}>
            <div style={{position: 'absolute', left: 0, top: 0, ...fadeUp(reads, 18)}}>
              <ReadBars delay={12} />
            </div>
            <div style={{position: 'absolute', left: 0, top: 265, ...fadeUp(vcf, 18)}}>
              <VcfMini delay={22} />
            </div>
            <div style={{position: 'absolute', left: 250, top: 205, opacity: vcf}}>
              <ArrowLine width={90} delay={30} color={COLORS.coral} />
            </div>
            <div
              style={{
                position: 'absolute',
                left: 342,
                top: 198,
                width: 120,
                height: 120,
                borderRadius: '50%',
                border: `2px dashed ${COLORS.coral}`,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontFamily: MONO,
                fontSize: 50,
                fontWeight: 900,
                color: COLORS.coral,
                background: 'rgba(255,253,252,0.5)',
                opacity: gapPulse,
                transform: `scale(${0.9 + gap * 0.1})`,
              }}
            >
              ?
            </div>
            <div
              style={{
                position: 'absolute',
                left: 295,
                top: 325,
                fontFamily: MONO,
                fontSize: 12,
                color: COLORS.coral,
                fontWeight: 800,
                opacity: gap,
                whiteSpace: 'nowrap',
              }}
            >
              missing bridge
            </div>
            <div style={{position: 'absolute', right: 0, top: 115, ...fadeUp(structure, 20)}}>
              <Structure delay={50} />
            </div>
            <div
              style={{
                position: 'absolute',
                left: 455,
                top: 110,
                opacity: interpolate(phase, [0.35, 0.55], [0, 1], {
                  extrapolateLeft: 'clamp',
                  extrapolateRight: 'clamp',
                }),
              }}
            >
              <Sparkle size={18} color={COLORS.gold} />
            </div>
            <div
              style={{
                position: 'absolute',
                right: 30,
                top: 50,
                opacity: 0.65,
                fontFamily: MONO,
                fontSize: 11,
                color: COLORS.muted,
                textAlign: 'right',
              }}
            >
              sequence → structure
              <br />
              function / properties
            </div>
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
