import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import type {SceneProps} from './types';

const ProblemCard: React.FC<{
  index: number;
  title: string;
  code: string;
  body: string;
  color: string;
  delay: number;
}> = ({index, title, code, body, color, delay}) => {
  const p = useEnter(delay, 14);
  return (
    <div
      style={{
        ...card,
        width: 470,
        minHeight: 310,
        padding: '28px 28px 26px',
        background: 'rgba(255,253,252,0.9)',
        ...fadeUp(p, 26),
        position: 'relative',
        overflow: 'hidden',
      }}
    >
      <div
        style={{
          position: 'absolute',
          right: -32,
          top: -34,
          width: 120,
          height: 120,
          borderRadius: '50%',
          background: `${color}1c`,
        }}
      />
      <div style={{display: 'flex', alignItems: 'center', justifyContent: 'space-between'}}>
        <span
          style={{
            fontFamily: MONO,
            fontSize: 13,
            fontWeight: 800,
            color,
            background: `${color}16`,
            borderRadius: 8,
            padding: '5px 9px',
          }}
        >
          {code}
        </span>
        <span style={{fontFamily: MONO, fontSize: 13, color: COLORS.muted}}>0{index}</span>
      </div>
      <div style={{fontSize: 30, fontWeight: 900, marginTop: 22, letterSpacing: -0.8}}>{title}</div>
      <div style={{fontSize: 17, lineHeight: 1.65, color: COLORS.muted, marginTop: 12}}>{body}</div>
      <div
        style={{
          marginTop: 22,
          height: 8,
          borderRadius: 8,
          background: `linear-gradient(90deg, ${color}, transparent)`,
          opacity: 0.55,
        }}
      />
    </div>
  );
};

export const ProblemScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const title = useEnter(2, 18);
  const path = interpolate(frame, [20, 60], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.lavender}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 20)}>
            <Pill color={COLORS.coral} bg={COLORS.coralSoft} style={{marginBottom: 16}}>
              旧流程 · per-exon bcftools
            </Pill>
            <h2
              style={{
                margin: 0,
                fontFamily: FONT,
                fontSize: 58,
                lineHeight: 1.15,
                letterSpacing: -2.4,
                fontWeight: 900,
              }}
            >
              只按外显子翻译，
              <br />
              <span style={{color: COLORS.coral}}>三个地方会迷路</span>
            </h2>
          </div>
          <div
            style={{
              width: 360,
              padding: '16px 20px',
              borderRadius: 18,
              background: 'rgba(255,253,252,0.75)',
              border: `1px dashed ${COLORS.coral}66`,
              color: COLORS.muted,
              fontSize: 14,
              lineHeight: 1.6,
              transform: `translateY(${(1 - path) * 12}px)`,
              opacity: path,
            }}
          >
            <span style={{fontFamily: MONO, color: COLORS.coral, fontWeight: 800}}>REF</span>{' '}
            对不上？
            <br />
            不是标签画错了，是坐标和链方向没对齐。
          </div>
        </div>

        <div
          style={{
            display: 'flex',
            gap: 34,
            justifyContent: 'center',
            marginTop: 54,
            position: 'relative',
          }}
        >
          <div
            style={{
              position: 'absolute',
              left: 60,
              right: 60,
              top: 154,
              height: 2,
              background: `repeating-linear-gradient(90deg, ${COLORS.lavender} 0 8px, transparent 8px 16px)`,
              opacity: 0.55,
              transform: `scaleX(${path})`,
              transformOrigin: 'left',
            }}
          />
          <ProblemCard
            index={1}
            code="strand"
            title="负链"
            color={COLORS.coral}
            body="cDNA 是反向互补的，VCF 的 REF/ALT 却是 forward strand；直接套，碱基当然对不上。"
            delay={12}
          />
          <ProblemCard
            index={2}
            code="junction"
            title="跨外显子"
            color={COLORS.teal}
            body="VCF 的 REF 可能跨过 intron；per-exon FASTA 里没有这段，deletion 就无处安放。"
            delay={22}
          />
          <ProblemCard
            index={3}
            code="SV"
            title="结构变异"
            color={COLORS.lavender}
            body="<INS>、<DUP>、<INV> 没有显式序列；bcftools consensus 也不支持这些符号型 ALT。"
            delay={32}
          />
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
