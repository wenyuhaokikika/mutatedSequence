import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {ArrowLine} from '../components/visuals';
import type {SceneProps} from './types';

const AnalogyCard: React.FC<{
  label: string;
  title: string;
  note: string;
  color: string;
  delay: number;
  children: React.ReactNode;
}> = ({label, title, note, color, delay, children}) => {
  const p = useEnter(delay, 15);
  return (
    <div
      style={{
        ...card,
        width: 470,
        height: 290,
        padding: '22px 24px',
        background: 'rgba(255,253,252,0.9)',
        ...fadeUp(p, 22),
        position: 'relative',
      }}
    >
      <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
        <span style={{fontFamily: MONO, fontSize: 11, fontWeight: 800, color}}>{label}</span>
        <span style={{fontSize: 12, color: COLORS.muted}}>{note}</span>
      </div>
      <div style={{marginTop: 16, height: 120}}>{children}</div>
      <div style={{fontSize: 22, fontWeight: 900, marginTop: 8, letterSpacing: -0.4}}>{title}</div>
    </div>
  );
};

export const MapScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const title = useEnter(2, 18);
  const arrows = interpolate(frame, [24, 52], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.teal}>
      <AbsoluteFill>
        <div style={fadeUp(title, 18)}>
          <Pill color={COLORS.teal} bg={COLORS.tealSoft} style={{marginBottom: 14}}>
            the gap, in plain language
          </Pill>
          <h2
            style={{
              margin: 0,
              fontFamily: FONT,
              fontSize: 56,
              lineHeight: 1.15,
              letterSpacing: -2.2,
              fontWeight: 900,
            }}
          >
            VCF 只是地址，
            <span style={{color: COLORS.coral}}>不是蛋白本身</span>
          </h2>
        </div>

        <div style={{display: 'flex', alignItems: 'center', gap: 18, marginTop: 40}}>
          <AnalogyCard label="GENOME" title="一本很厚的菜谱" note="参考基因组" color={COLORS.coral} delay={10}>
            <div style={{display: 'flex', gap: 8, alignItems: 'flex-end', height: 110}}>
              <div
                style={{
                  width: 120,
                  height: 100,
                  borderRadius: 12,
                  background: COLORS.coralSoft,
                  border: `2px solid ${COLORS.coral}`,
                  position: 'relative',
                }}
              >
                <div style={{position: 'absolute', left: 14, right: 14, top: 18, height: 5, borderRadius: 4, background: COLORS.coral, opacity: 0.6}} />
                <div style={{position: 'absolute', left: 14, right: 30, top: 32, height: 5, borderRadius: 4, background: COLORS.coral, opacity: 0.35}} />
                <div style={{position: 'absolute', left: 14, right: 20, top: 46, height: 5, borderRadius: 4, background: COLORS.coral, opacity: 0.35}} />
                <div
                  style={{
                    position: 'absolute',
                    left: 74,
                    top: 56,
                    width: 18,
                    height: 18,
                    borderRadius: '50%',
                    background: COLORS.coral,
                    border: '3px solid white',
                    boxShadow: '0 0 0 5px rgba(255,111,97,0.14)',
                  }}
                />
              </div>
              <div style={{fontFamily: MONO, fontSize: 11, color: COLORS.muted, lineHeight: 1.5}}>
                A C G T
                <br />
                第 12345 位
              </div>
            </div>
          </AnalogyCard>

          <div style={{opacity: arrows}}>
            <ArrowLine width={60} delay={24} color={COLORS.teal} />
          </div>

          <AnalogyCard label="TRANSCRIPT" title="抄出来的一页" note="剪掉 intron" color={COLORS.teal} delay={20}>
            <div style={{display: 'flex', gap: 6, alignItems: 'center', height: 100}}>
              <div style={{width: 54, height: 30, borderRadius: 7, background: COLORS.coral}} />
              <div style={{width: 58, height: 2, background: 'rgba(45,42,50,0.2)'}} />
              <div style={{width: 74, height: 30, borderRadius: 7, background: COLORS.teal}} />
              <div style={{width: 58, height: 2, background: 'rgba(45,42,50,0.2)'}} />
              <div style={{width: 62, height: 30, borderRadius: 7, background: COLORS.lavender}} />
            </div>
          </AnalogyCard>

          <div style={{opacity: arrows}}>
            <ArrowLine width={60} delay={30} color={COLORS.gold} />
          </div>

          <AnalogyCard label="PROTEIN" title="照着做出来的菜" note="氨基酸链" color={COLORS.gold} delay={30}>
            <div style={{display: 'flex', alignItems: 'center', gap: 4, height: 100}}>
              {['M', 'A', 'L', 'W', 'R', 'L', 'P'].map((aa, i) => (
                <div
                  key={i}
                  style={{
                    width: 34,
                    height: 34,
                    borderRadius: '50%',
                    background: i === 3 ? COLORS.coral : [COLORS.teal, COLORS.lavender, COLORS.gold, COLORS.peach][i % 4],
                    color: i === 3 ? COLORS.paper : COLORS.ink,
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    fontFamily: MONO,
                    fontWeight: 800,
                    fontSize: 13,
                    border: `1px solid ${COLORS.line}`,
                  }}
                >
                  {aa}
                </div>
              ))}
            </div>
          </AnalogyCard>
        </div>

        <div
          style={{
            display: 'flex',
            gap: 10,
            marginTop: 26,
            justifyContent: 'center',
            opacity: arrows,
          }}
        >
          {['正反方向', '剪接 / intron', 'indel / stop', '不同起始位点'].map((x, i) => (
            <span
              key={x}
              style={{
                fontFamily: MONO,
                fontSize: 12,
                color: i === 2 ? COLORS.coral : COLORS.muted,
                background: 'rgba(255,253,252,0.72)',
                border: `1px dashed ${i === 2 ? COLORS.coral : COLORS.line}`,
                borderRadius: 8,
                padding: '7px 10px',
              }}
            >
              {x}
            </span>
          ))}
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
