import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import type {SceneProps} from './types';

const rows = [
  {
    tool: 'VEP / SnpEff',
    tag: '像查字典',
    body: '告诉你变异落在哪个基因、可能什么后果；不保证给你完整 cDNA、CDS 和蛋白。',
    color: COLORS.lavender,
  },
  {
    tool: 'bcftools csq',
    tag: '很快的注释',
    body: '能一起考虑同一段 DNA 上的多个变异，但输出仍是注释，不是完整序列。',
    color: COLORS.teal,
  },
  {
    tool: 'Haplosaurus / VCF2Prot',
    tag: '两套蛋白',
    body: '能在知道哪份拷贝是哪份时拼出蛋白；通常要先做这一步，主要做 CDS / 蛋白。',
    color: COLORS.gold,
  },
  {
    tool: 'ORFan',
    tag: '转录本 + 短肽',
    body: '直接出突变转录本和短肽；比较新，主要做编码区。',
    color: COLORS.peach,
  },
];

const CompareRow: React.FC<{
  tool: string;
  tag: string;
  body: string;
  color: string;
  delay: number;
}> = ({tool, tag, body, color, delay}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        ...card,
        padding: '14px 18px',
        display: 'grid',
        gridTemplateColumns: '250px 160px 1fr',
        gap: 18,
        alignItems: 'center',
        background: 'rgba(255,253,252,0.76)',
        opacity: p,
        transform: `translateY(${(1 - p) * 12}px)`,
      }}
    >
      <span style={{fontFamily: MONO, fontSize: 14, fontWeight: 800}}>{tool}</span>
      <span
        style={{
          justifySelf: 'start',
          fontFamily: MONO,
          fontSize: 11,
          fontWeight: 800,
          color,
          background: `${color}18`,
          borderRadius: 6,
          padding: '4px 7px',
        }}
      >
        {tag}
      </span>
      <span style={{fontSize: 14, color: COLORS.muted, lineHeight: 1.5}}>{body}</span>
    </div>
  );
};

export const CompareScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const note = useEnter(100);
  const highlight = useEnter(52);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.gold}>
      <AbsoluteFill>
        <div style={fadeUp(title, 18)}>
          <Pill color={COLORS.gold} bg={COLORS.goldSoft} style={{marginBottom: 14}}>
            where it fits
          </Pill>
          <h2
            style={{
              margin: 0,
              fontFamily: FONT,
              fontSize: 54,
              lineHeight: 1.15,
              letterSpacing: -2.2,
              fontWeight: 900,
            }}
          >
            和市面上的工具比，
            <span style={{color: COLORS.teal}}>它站在哪？</span>
          </h2>
        </div>

        <div style={{display: 'flex', flexDirection: 'column', gap: 10, marginTop: 28}}>
          {rows.map((row, i) => (
            <CompareRow
              key={row.tool}
              tool={row.tool}
              tag={row.tag}
              body={row.body}
              color={row.color}
              delay={16 + i * 8}
            />
          ))}

          <div
            style={{
              ...card,
              padding: '16px 18px',
              display: 'grid',
              gridTemplateColumns: '250px 160px 1fr',
              gap: 18,
              alignItems: 'center',
              background:
                'linear-gradient(110deg, rgba(255,225,219,0.9), rgba(221,243,238,0.8))',
              border: `2px solid ${COLORS.coral}`,
              opacity: highlight,
              transform: `translateY(${(1 - highlight) * 12}px)`,
            }}
          >
            <span style={{fontFamily: MONO, fontSize: 15, fontWeight: 900, color: COLORS.coral}}>
              mutatedSequence
            </span>
            <span
              style={{
                justifySelf: 'start',
                fontFamily: MONO,
                fontSize: 11,
                fontWeight: 900,
                color: COLORS.paper,
                background: COLORS.coral,
                borderRadius: 6,
                padding: '4px 7px',
              }}
            >
              完整序列 + 审计
            </span>
            <span style={{fontSize: 14, color: COLORS.ink, lineHeight: 1.5, fontWeight: 600}}>
              cDNA / CDS / 蛋白 + 前后对照表；负链、跨外显子、不完整转录本、常见结构变异。
            </span>
          </div>
        </div>

        <div
          style={{
            marginTop: 14,
            fontFamily: MONO,
            fontSize: 12,
            color: COLORS.muted,
            textAlign: 'center',
            ...fadeUp(note, 10),
          }}
        >
          不是全面更好，是解决另一段路：从 VCF 到可复查的完整序列。
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
