import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {ArrowLine} from '../components/visuals';
import type {SceneProps} from './types';

const GenotypeChip: React.FC<{
  gt: string;
  active?: boolean;
  delay: number;
  label: string;
}> = ({gt, active, delay, label}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        padding: '13px 16px',
        borderRadius: 14,
        background: active ? COLORS.tealSoft : 'rgba(255,253,252,0.82)',
        border: `1px solid ${active ? COLORS.teal : COLORS.line}`,
        transform: `translateX(${(1 - p) * -18}px)`,
        opacity: p,
      }}
    >
      <div style={{display: 'flex', alignItems: 'center', gap: 12}}>
        <span style={{fontFamily: MONO, fontSize: 18, fontWeight: 900, color: active ? COLORS.teal : COLORS.ink}}>
          {gt}
        </span>
        <span style={{fontSize: 13, color: COLORS.muted}}>{label}</span>
      </div>
      <span
        style={{
          width: 9,
          height: 9,
          borderRadius: 9,
          background: active ? COLORS.teal : 'rgba(45,42,50,0.16)',
        }}
      />
    </div>
  );
};

const AlleleCard: React.FC<{
  title: string;
  detail: string;
  color: string;
  active?: boolean;
  delay: number;
}> = ({title, detail, color, active, delay}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        ...card,
        width: 340,
        padding: '22px 24px',
        background: active ? 'rgba(255,253,252,0.98)' : 'rgba(255,253,252,0.65)',
        border: `1.5px solid ${active ? color : COLORS.line}`,
        boxShadow: active ? COLORS.shadow : COLORS.shadowSoft,
        transform: `translateY(${(1 - p) * 20}px) scale(${0.97 + p * 0.03})`,
        opacity: p,
      }}
    >
      <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
        <span style={{fontFamily: MONO, fontSize: 13, color, fontWeight: 800}}>{title}</span>
        {active ? (
          <span
            style={{
              fontFamily: MONO,
              fontSize: 11,
              fontWeight: 800,
              color: COLORS.teal,
              background: COLORS.tealSoft,
              padding: '4px 8px',
              borderRadius: 6,
            }}
          >
            selected
          </span>
        ) : null}
      </div>
      <div style={{fontSize: 17, color: COLORS.muted, lineHeight: 1.55, marginTop: 12}}>{detail}</div>
    </div>
  );
};

export const AlleleScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const left = useEnter(8);
  const arrow = useEnter(30);
  const right = useEnter(24);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.coral}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'center', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.coral} bg={COLORS.coralSoft} style={{marginBottom: 14}}>
              step 02 · 选版本
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
              两份拷贝，
              <br />
              <span style={{color: COLORS.coral}}>只选一个版本</span>
            </h2>
          </div>
          <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.6, textAlign: 'right'}}>
            0/1 · 1/1 · 1/2 · phased
            <br />
            不要把两个等位基因挤进一条序列
          </div>
        </div>

        <div style={{display: 'flex', alignItems: 'center', gap: 34, marginTop: 42}}>
          <div style={{width: 540, ...fadeUp(left, 18)}}>
            <div style={{fontFamily: MONO, fontSize: 13, color: COLORS.muted, marginBottom: 12}}>
              基因型
            </div>
            <div style={{display: 'flex', flexDirection: 'column', gap: 10}}>
              <GenotypeChip gt="0/1" label="一份参考，一份突变" active delay={12} />
              <GenotypeChip gt="1/1" label="两份都突变" delay={18} />
              <GenotypeChip gt="1/2" label="两份突变不一样" delay={24} />
              <GenotypeChip gt="1|0" label="phased：知道哪份是哪份" delay={30} />
            </div>
            <div
              style={{
                marginTop: 16,
                fontFamily: MONO,
                fontSize: 12,
                color: COLORS.muted,
                lineHeight: 1.6,
                opacity: left,
              }}
            >
              两条拷贝不是魔法：
              <br />
              没 phasing 时，只是按基因型顺序取。
            </div>
          </div>

          <div style={{display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 10, opacity: arrow}}>
            <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.coral}}>choose</span>
            <ArrowLine width={100} delay={30} />
          </div>

          <div style={{flex: 1, ...fadeUp(right, 16)}}>
            <div style={{display: 'flex', gap: 18}}>
              <AlleleCard
                title="版本 1 · 0/1"
                detail="C → A，选中的突变版本；原始 REF/ALT 与基因型一起保留。"
                color={COLORS.coral}
                active
                delay={34}
              />
              <AlleleCard
                title="版本 2 · 1/2"
                detail="另一个版本不会被偷偷应用；需要时单独跑 hap2。"
                color={COLORS.lavender}
                delay={42}
              />
            </div>
            <div
              style={{
                marginTop: 20,
                padding: '14px 18px',
                borderRadius: 15,
                background: 'rgba(255,253,252,0.75)',
                border: `1px dashed ${COLORS.coral}55`,
                display: 'flex',
                alignItems: 'center',
                gap: 12,
                ...fadeUp(right, 10),
              }}
            >
              <span style={{fontSize: 17, color: COLORS.coral, fontWeight: 800}}>×</span>
              <span style={{fontSize: 14, color: COLORS.muted, lineHeight: 1.5}}>
                两个版本混在一起？不行。一条序列只能有一个版本。
              </span>
            </div>
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
