import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import type {SceneProps} from './types';

const LimitCard: React.FC<{
  k: string;
  title: string;
  body: string;
  delay: number;
}> = ({k, title, body, delay}) => {
  const p = useEnter(delay, 15);
  return (
    <div
      style={{
        ...card,
        width: 485,
        minHeight: 270,
        padding: '24px 26px',
        background: 'rgba(255,253,252,0.84)',
        border: `1px solid ${COLORS.lavender}66`,
        ...fadeUp(p, 20),
        position: 'relative',
      }}
    >
      <div style={{display: 'flex', alignItems: 'center', justifyContent: 'space-between'}}>
        <span
          style={{
            fontFamily: MONO,
            fontSize: 12,
            fontWeight: 800,
            color: COLORS.lavender,
            background: COLORS.lilac,
            borderRadius: 7,
            padding: '4px 8px',
          }}
        >
          {k}
        </span>
        <span
          style={{
            fontFamily: MONO,
            fontSize: 10,
            fontWeight: 800,
            color: COLORS.muted,
            border: `1px solid ${COLORS.line}`,
            borderRadius: 6,
            padding: '3px 6px',
          }}
        >
          not modeled
        </span>
      </div>
      <div style={{fontSize: 26, fontWeight: 900, marginTop: 18, letterSpacing: -0.5}}>{title}</div>
      <div style={{fontSize: 15, lineHeight: 1.6, color: COLORS.muted, marginTop: 10}}>{body}</div>
    </div>
  );
};

export const LimitsScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const bubble = useEnter(40);
  const supported = useEnter(52);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.lavender}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.lavender} bg={COLORS.lilac} style={{marginBottom: 14}}>
              honest limits
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
              它也会
              <span style={{color: COLORS.coral}}>诚实地说</span>
              <br />
              “这个我处理不了”
            </h2>
          </div>
          <div
            style={{
              maxWidth: 380,
              padding: '14px 18px',
              borderRadius: 18,
              background: 'rgba(255,253,252,0.8)',
              border: `1px dashed ${COLORS.lavender}`,
              color: COLORS.muted,
              fontSize: 14,
              lineHeight: 1.6,
              ...fadeUp(bubble, 12),
            }}
          >
            不求“看起来什么都行”。
            <br />
            不确定的地方，明确标出来。
          </div>
        </div>

        <div style={{display: 'flex', gap: 24, marginTop: 42, justifyContent: 'center'}}>
          <LimitCard
            k="phasing"
            title="两份菜谱混在一起"
            body="没有 phasing，两个副本的改动会变成一个 mosaic；看起来像，但不等于真实的两条单倍型。"
            delay={12}
          />
          <LimitCard
            k="splicing"
            title="跳页 / 剪接"
            body="剪接位点让页面被跳过或重排；工具仍按原目录翻，不会自动生成新版本。"
            delay={20}
          />
          <LimitCard
            k="SV"
            title="搬家 / 复制 / 倒装"
            body="TRA、BND、CNV 把段落搬到别的书里，或整章复制；一条线性转录本装不下。"
            delay={28}
          />
        </div>

        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            gap: 16,
            marginTop: 28,
            ...fadeUp(supported, 14),
          }}
        >
          <span style={{fontSize: 13, color: COLORS.muted}}>supported:</span>
          {['<DEL>', '<INS>', '<DUP>', '<INV>'].map((x) => (
            <span
              key={x}
              style={{
                fontFamily: MONO,
                fontSize: 13,
                fontWeight: 800,
                color: COLORS.green,
                background: COLORS.tealSoft,
                borderRadius: 8,
                padding: '6px 10px',
              }}
            >
              {x}
            </span>
          ))}
          <span style={{fontSize: 13, color: COLORS.muted}}>with documented caveats</span>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
