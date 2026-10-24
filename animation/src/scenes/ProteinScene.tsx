import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {ProteinChain, SequenceStrip} from '../components/visuals';
import type {SceneProps} from './types';

const StageCard: React.FC<{
  step: string;
  title: string;
  color: string;
  delay: number;
  children: React.ReactNode;
}> = ({step, title, color, delay, children}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        ...card,
        width: 396,
        height: 430,
        padding: '22px 22px 20px',
        background: 'rgba(255,253,252,0.88)',
        border: `1px solid ${color}33`,
        ...fadeUp(p, 22),
        position: 'relative',
        overflow: 'hidden',
      }}
    >
      <div style={{display: 'flex', alignItems: 'center', justifyContent: 'space-between'}}>
        <span
          style={{
            fontFamily: MONO,
            fontSize: 12,
            fontWeight: 800,
            color,
            background: `${color}16`,
            borderRadius: 7,
            padding: '4px 8px',
          }}
        >
          {step}
        </span>
        <span style={{fontSize: 14, fontWeight: 800, color: COLORS.ink}}>{title}</span>
      </div>
      <div style={{marginTop: 18}}>{children}</div>
    </div>
  );
};

const ExonMini: React.FC<{label: string; color: string; width: number}> = ({label, color, width}) => (
  <div
    style={{
      width,
      height: 34,
      borderRadius: 8,
      background: color,
      color: COLORS.paper,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      fontFamily: MONO,
      fontSize: 10,
      fontWeight: 800,
    }}
  >
    {label}
  </div>
);

export const ProteinScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.lavender}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.lavender} bg={COLORS.lilac} style={{marginBottom: 14}}>
              step 04 · 转录本 → 蛋白
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
              外显子拼起来，
              <span style={{color: COLORS.teal}}>蛋白就出来了</span>
            </h2>
          </div>
          <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.6, textAlign: 'right'}}>
            正常蛋白停在终止密码子；
            <br />
            3′ 延伸负责「停不下来」的情况。
          </div>
        </div>

        <div style={{display: 'flex', gap: 20, marginTop: 30, alignItems: 'stretch'}}>
          <StageCard step="外显子" title="剪接外显子" color={COLORS.coral} delay={10}>
            <div style={{display: 'flex', flexDirection: 'column', gap: 10, marginTop: 32}}>
              <div style={{display: 'flex', gap: 6, alignItems: 'center'}}>
                <ExonMini label="E1" color={COLORS.coral} width={84} />
                <div style={{height: 2, flex: 1, background: 'rgba(45,42,50,0.18)'}} />
                <ExonMini label="E2" color={COLORS.teal} width={92} />
                <div style={{height: 2, flex: 1, background: 'rgba(45,42,50,0.18)'}} />
                <ExonMini label="E3" color={COLORS.lavender} width={94} />
              </div>
              <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 14}}>
                按注释的外显子顺序切取；负链整体反向互补。
              </div>
              <div style={{display: 'flex', gap: 6, marginTop: 12}}>
                {['5′ UTR', 'CDS', '3′ UTR'].map((x) => (
                  <span
                    key={x}
                    style={{
                      fontFamily: MONO,
                      fontSize: 10,
                      color: COLORS.muted,
                      border: `1px solid ${COLORS.line}`,
                      borderRadius: 5,
                      padding: '3px 5px',
                    }}
                  >
                    {x}
                  </span>
                ))}
              </div>
            </div>
          </StageCard>

          <StageCard step="转录本" title="完整抄页" color={COLORS.teal} delay={18}>
            <SequenceStrip
              label="转录本 · 5′ → 3′"
              sequence="ATGGCACGTACGTAGCTAGCTAGC"
              color={COLORS.teal}
              highlight={5}
              delay={22}
            />
            <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 16}}>
              完整转录本：两端非编码区、外显子、剪接点都在。
            </div>
            <div
              style={{
                marginTop: 18,
                padding: '10px 12px',
                borderRadius: 10,
                background: COLORS.tealSoft,
                fontFamily: MONO,
                fontSize: 11,
                color: COLORS.teal,
                fontWeight: 700,
              }}
            >
              转录本 = 外显子拼接
            </div>
          </StageCard>

          <StageCard step="编码区" title="编码区" color={COLORS.gold} delay={26}>
            <SequenceStrip
              label="编码区"
              sequence="NNATGGCACGTACGTAGCTAGCTAA"
              color={COLORS.gold}
              highlight={2}
              delay={30}
            />
            <div style={{display: 'flex', gap: 6, marginTop: 14, flexWrap: 'wrap'}}>
              <span style={{fontSize: 11, color: COLORS.coral, fontFamily: MONO, fontWeight: 800}}>
                5′ 缺 N
              </span>
              <span style={{fontSize: 11, color: COLORS.muted}}>·</span>
              <span style={{fontSize: 11, color: COLORS.green, fontFamily: MONO, fontWeight: 800}}>
                stop TAA
              </span>
            </div>
            <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 14}}>
              有些注释不把终止密码子写进编码区，要单独补上；5′ 缺失时用 N 补齐。
            </div>
          </StageCard>

          <StageCard step="蛋白" title="翻译" color={COLORS.lavender} delay={34}>
            <div style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted, marginBottom: 8}}>
              正常翻译
            </div>
            <ProteinChain letters="MALWMRLLP" highlight={4} size={27} gap={4} delay={38} />
            <div style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted, marginTop: 18, marginBottom: 8}}>
              3′ 延伸
            </div>
            <ProteinChain letters="MALWMRLLP" highlight={4} size={27} gap={4} delay={44} />
            <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 14}}>
              终止密码子丢失或移码时，继续翻译到下一个终止密码子。
            </div>
          </StageCard>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
