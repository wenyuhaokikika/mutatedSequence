import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {ArrowLine} from '../components/visuals';
import type {SceneProps} from './types';

const ExonBox: React.FC<{
  left: string;
  width: string;
  label: string;
  color: string;
  delay: number;
}> = ({left, width, label, color, delay}) => {
  const p = useEnter(delay, 14);
  return (
    <div
      style={{
        position: 'absolute',
        left,
        width,
        top: 0,
        height: 58,
        borderRadius: 13,
        background: color,
        border: `1px solid ${COLORS.line}`,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        fontFamily: MONO,
        fontSize: 12,
        fontWeight: 800,
        color: COLORS.paper,
        transform: `scaleX(${p})`,
        transformOrigin: 'left',
        opacity: p,
      }}
    >
      {label}
    </div>
  );
};

export const BlocksScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const title = useEnter(2, 18);
  const line = interpolate(frame, [8, 38], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });
  const bracket = useEnter(28);
  const extraction = useEnter(58);
  const blockRail = useEnter(40);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.teal}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.teal} bg={COLORS.tealSoft} style={{marginBottom: 14}}>
              step 03 · genomic blocks
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
              先在基因组上改，
              <br />
              <span style={{color: COLORS.teal}}>再切回外显子</span>
            </h2>
          </div>
          <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.6, textAlign: 'right'}}>
            reference block / ALT block / insertion block
            <br />
            deletion 跨 intron 也不会丢
          </div>
        </div>

        <div
          style={{
            marginTop: 34,
            padding: '26px 34px 18px',
            ...card,
            background: 'rgba(255,253,252,0.86)',
            position: 'relative',
          }}
        >
          <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: 16}}>
            <span style={{fontFamily: MONO, fontSize: 13, color: COLORS.muted}}>reference genome</span>
            <span style={{fontSize: 12, color: COLORS.muted}}>5′ → 3′ · plus strand</span>
          </div>
          <div style={{position: 'relative', height: 98, margin: '0 10px'}}>
            <div
              style={{
                position: 'absolute',
                left: 0,
                right: 0,
                top: 29,
                height: 9,
                borderRadius: 9,
                background: 'rgba(45,42,50,0.12)',
                transform: `scaleX(${line})`,
                transformOrigin: 'left',
              }}
            />
            <div
              style={{
                position: 'absolute',
                left: 0,
                right: 0,
                top: 0,
                height: 58,
                opacity: line,
              }}
            >
              <ExonBox left="7%" width="14%" label="exon 1" color={COLORS.coral} delay={12} />
              <ExonBox left="33%" width="15%" label="exon 2" color={COLORS.teal} delay={16} />
              <ExonBox left="65%" width="17%" label="exon 3" color={COLORS.lavender} delay={20} />
              <div
                style={{
                  position: 'absolute',
                  left: '21%',
                  width: '12%',
                  top: 28,
                  borderTop: `2px solid ${COLORS.muted}`,
                  opacity: 0.35,
                }}
              />
              <div
                style={{
                  position: 'absolute',
                  left: '48%',
                  width: '17%',
                  top: 28,
                  borderTop: `2px solid ${COLORS.muted}`,
                  opacity: 0.35,
                }}
              />
              <div
                style={{
                  position: 'absolute',
                  left: '82%',
                  width: '12%',
                  top: 28,
                  borderTop: `2px solid ${COLORS.muted}`,
                  opacity: 0.35,
                }}
              />
            </div>

            <div
              style={{
                position: 'absolute',
                left: '24.5%',
                top: 5,
                width: '15%',
                height: 28,
                borderLeft: `2px dashed ${COLORS.coral}`,
                borderRight: `2px dashed ${COLORS.coral}`,
                borderTop: `2px dashed ${COLORS.coral}`,
                borderTopLeftRadius: 8,
                borderTopRightRadius: 8,
                opacity: bracket,
                transform: `translateY(${(1 - bracket) * -8}px)`,
              }}
            />
            <div
              style={{
                position: 'absolute',
                left: '27.5%',
                top: -24,
                fontSize: 12,
                fontFamily: MONO,
                color: COLORS.coral,
                fontWeight: 800,
                opacity: bracket,
              }}
            >
              DEL spans intron
            </div>
            <div
              style={{
                position: 'absolute',
                left: '77%',
                top: 17,
                width: 0,
                height: 0,
                borderLeft: '8px solid transparent',
                borderRight: '8px solid transparent',
                borderBottom: `12px solid ${COLORS.gold}`,
                opacity: bracket,
              }}
            />
            <div
              style={{
                position: 'absolute',
                left: '74.5%',
                top: -10,
                fontSize: 12,
                fontFamily: MONO,
                color: COLORS.green,
                fontWeight: 800,
                opacity: bracket,
              }}
            >
              INS
            </div>
            <div
              style={{
                position: 'absolute',
                left: '36%',
                top: 23,
                width: 12,
                height: 12,
                borderRadius: 12,
                background: COLORS.coral,
                border: '3px solid white',
                boxShadow: '0 0 0 5px rgba(255,111,97,0.14)',
                opacity: bracket,
              }}
            />
          </div>
        </div>

        <div
          style={{
            display: 'flex',
            gap: 16,
            marginTop: 20,
            alignItems: 'center',
            ...fadeUp(blockRail, 12),
          }}
        >
          <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted, width: 70}}>blocks →</span>
          {[
            {label: 'reference', color: COLORS.teal, detail: '从 genome 取'},
            {label: 'ALT', color: COLORS.coral, detail: '替换/删除'},
            {label: 'insertion', color: COLORS.gold, detail: '插在中间'},
          ].map((b) => (
            <div
              key={b.label}
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: 10,
                padding: '10px 14px',
                borderRadius: 12,
                background: 'rgba(255,253,252,0.75)',
                border: `1px solid ${COLORS.line}`,
              }}
            >
              <span
                style={{
                  width: 18,
                  height: 18,
                  borderRadius: 5,
                  background: b.color,
                  display: 'inline-block',
                }}
              />
              <div>
                <div style={{fontFamily: MONO, fontSize: 12, fontWeight: 800}}>{b.label}</div>
                <div style={{fontSize: 12, color: COLORS.muted}}>{b.detail}</div>
              </div>
            </div>
          ))}
        </div>

        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: 20,
            marginTop: 20,
            ...fadeUp(extraction, 16),
          }}
        >
          <div
            style={{
              ...card,
              padding: '14px 18px',
              display: 'flex',
              alignItems: 'center',
              gap: 12,
              background: 'rgba(221,243,238,0.7)',
              border: `1px solid ${COLORS.teal}33`,
            }}
          >
            <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.teal, fontWeight: 800}}>
              exon 1 + exon 2
            </span>
            <span style={{fontSize: 14, color: COLORS.muted}}>intersect blocks</span>
            <span style={{fontSize: 18, color: COLORS.teal}}>→</span>
            <span style={{fontSize: 14, fontWeight: 800, color: COLORS.ink}}>mutated cDNA</span>
          </div>
          <ArrowLine width={80} delay={64} color={COLORS.coral} />
          <span style={{fontSize: 14, color: COLORS.muted}}>
            哪怕 deletion 横跨 intron，也只切走 exonic 部分。
          </span>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
