import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {Base, ChromosomeRibbon, VcfCard} from '../components/visuals';
import type {SceneProps} from './types';

const CheckRow: React.FC<{
  vcf: string;
  genome: string;
  ok: boolean;
  delay: number;
  note: string;
}> = ({vcf, genome, ok, delay, note}) => {
  const p = useEnter(delay);
  const color = ok ? COLORS.green : COLORS.red;
  return (
    <div
      style={{
        ...card,
        padding: '20px 22px',
        display: 'flex',
        alignItems: 'center',
        gap: 18,
        background: ok ? 'rgba(221,243,238,0.8)' : 'rgba(255,225,219,0.65)',
        border: `1px solid ${color}33`,
        ...fadeUp(p, 16),
      }}
    >
      <div
        style={{
          width: 42,
          height: 42,
          borderRadius: 14,
          background: ok ? COLORS.tealSoft : COLORS.coralSoft,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          color,
          fontWeight: 900,
          fontSize: 22,
        }}
      >
        {ok ? '✓' : '!'}
      </div>
      <div style={{flex: 1}}>
        <div style={{display: 'flex', alignItems: 'center', gap: 10}}>
          <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted}}>VCF 说</span>
          <Base value={vcf} bg="rgba(255,255,255,0.8)" color={COLORS.ink} size={28} />
          <span style={{fontSize: 14, color: COLORS.muted}}>vs</span>
          <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted}}>参考</span>
          <Base value={genome} bg="rgba(255,255,255,0.8)" color={color} size={28} />
        </div>
        <div style={{fontSize: 14, color: COLORS.muted, marginTop: 7, lineHeight: 1.4}}>{note}</div>
      </div>
      <span
        style={{
          fontFamily: MONO,
          fontSize: 12,
          fontWeight: 800,
          color,
          background: 'rgba(255,255,255,0.65)',
          borderRadius: 7,
          padding: '5px 8px',
        }}
      >
        {ok ? '对上了' : '对不上'}
      </span>
    </div>
  );
};

export const RefCheckScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const left = useEnter(10, 15);
  const right = useEnter(18);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.teal}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.teal} bg={COLORS.tealSoft} style={{marginBottom: 14}}>
              step 01 · 地址检查
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
              先对暗号：
              <span style={{color: COLORS.teal}}>原来的碱基对得上吗？</span>
            </h2>
          </div>
          <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.6, textAlign: 'right'}}>
            先确认版本、坐标、染色体名
            <br />
            对不上，就不要硬翻译
          </div>
        </div>

        <div style={{display: 'flex', gap: 46, marginTop: 38, alignItems: 'flex-start'}}>
          <div style={{width: 710, ...fadeUp(left, 20)}}>
            <ChromosomeRibbon width={680} highlight={0.42} delay={12} markerLabel="变异位置" />
            <div style={{display: 'flex', gap: 16, marginTop: 22, alignItems: 'center'}}>
              <VcfCard refBase="C" altBase="T" gt="0/1" delay={16} mismatch />
              <div style={{fontSize: 30, color: COLORS.coral, fontWeight: 800}}>≠</div>
              <div
                style={{
                  ...card,
                  padding: '18px 20px',
                  width: 220,
                  background: 'rgba(255,253,252,0.9)',
                }}
              >
                <div style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted}}>参考基因组这一位</div>
                <div style={{display: 'flex', alignItems: 'center', gap: 8, marginTop: 10}}>
                  <Base value="A" bg="rgba(255,255,255,0.8)" color={COLORS.teal} />
                  <span style={{fontSize: 13, color: COLORS.muted}}>实际是</span>
                </div>
              </div>
            </div>
            <div style={{display: 'flex', gap: 8, marginTop: 18, flexWrap: 'wrap'}}>
              {['同一版本', '染色体名', '大小写'].map((x) => (
                <span
                  key={x}
                  style={{
                    fontFamily: MONO,
                    fontSize: 11,
                    color: COLORS.muted,
                    background: 'rgba(255,253,252,0.7)',
                    border: `1px solid ${COLORS.line}`,
                    borderRadius: 6,
                    padding: '4px 7px',
                  }}
                >
                  {x}
                </span>
              ))}
            </div>
          </div>

          <div style={{flex: 1, ...fadeUp(right, 18)}}>
            <div style={{fontFamily: MONO, fontSize: 13, color: COLORS.muted, marginBottom: 12}}>
              两种结果
            </div>
            <div style={{display: 'flex', flexDirection: 'column', gap: 14}}>
              <CheckRow
                vcf="A"
                genome="A"
                ok
                delay={24}
                note="对上了，继续往下走。"
              />
              <CheckRow
                vcf="C"
                genome="A"
                ok={false}
                delay={32}
                note="对不上：先别翻译，检查地图。"
              />
            </div>
            <div
              style={{
                marginTop: 18,
                padding: '14px 18px',
                borderRadius: 15,
                background: 'rgba(255,253,252,0.78)',
                border: `1px dashed ${COLORS.red}33`,
                fontSize: 13,
                color: COLORS.muted,
                lineHeight: 1.55,
                ...fadeUp(right, 10),
              }}
            >
              对不上，不一定是工具坏了：先看 assembly、染色体名和坐标，再决定要不要继续。
            </div>
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
