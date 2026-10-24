import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {ArrowLine} from '../components/visuals';
import type {SceneProps} from './types';

const Codon: React.FC<{
  text: string;
  color?: string;
  bg?: string;
  strike?: boolean;
  delay: number;
}> = ({text, color = COLORS.ink, bg = 'rgba(255,253,252,0.6)', strike, delay}) => {
  const p = useEnter(delay);
  return (
    <span
      style={{
        fontFamily: MONO,
        fontSize: 20,
        fontWeight: 800,
        color,
        background: bg,
        border: `1px solid ${color}44`,
        borderRadius: 8,
        padding: '9px 11px',
        textDecoration: strike ? 'line-through' : 'none',
        opacity: p,
        transform: `scale(${0.9 + p * 0.1})`,
      }}
    >
      {text}
    </span>
  );
};

export const StopScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const title = useEnter(2, 18);
  const rail = useEnter(12);
  const broken = useEnter(34);
  const travel = interpolate(frame, [50, 150], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });
  const outcome = useEnter(76);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.coral}>
      <AbsoluteFill>
        <div style={fadeUp(title, 18)}>
          <Pill color={COLORS.coral} bg={COLORS.coralSoft} style={{marginBottom: 14}}>
            the stop-loss question
          </Pill>
          <h2
            style={{
              margin: 0,
              fontFamily: FONT,
              fontSize: 52,
              lineHeight: 1.15,
              letterSpacing: -2,
              fontWeight: 900,
            }}
          >
            终止密码子失灵，
            <span style={{color: COLORS.coral}}>会发生什么？</span>
          </h2>
        </div>

        <div style={{display: 'flex', gap: 34, marginTop: 34, alignItems: 'flex-start'}}>
          <div
            style={{
              width: 1000,
              ...card,
              padding: '22px 24px 20px',
              background: 'rgba(255,253,252,0.9)',
              ...fadeUp(rail, 18),
            }}
          >
            <div style={{display: 'flex', justifyContent: 'space-between', marginBottom: 18}}>
              <span style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted}}>CDS · reading frame</span>
              <span style={{fontSize: 12, color: COLORS.muted}}>三个碱基一个密码子</span>
            </div>
            <div style={{display: 'flex', alignItems: 'center', gap: 8, flexWrap: 'wrap'}}>
              {['ATG', 'GCT', 'TGG', 'CGT', 'TAT', 'GCT', 'TAA'].map((codon, i) => (
                <Codon
                  key={codon + i}
                  text={codon}
                  color={i === 6 ? COLORS.coral : COLORS.ink}
                  bg={i === 6 ? COLORS.coralSoft : 'rgba(255,253,252,0.6)'}
                  delay={16 + i * 3}
                />
              ))}
              <span style={{fontFamily: MONO, fontSize: 16, color: COLORS.muted}}>...</span>
            </div>

            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: 12,
                marginTop: 20,
                ...fadeUp(broken, 14),
              }}
            >
              <span style={{fontFamily: MONO, fontSize: 13, color: COLORS.muted}}>mutation</span>
              <span style={{fontFamily: MONO, fontSize: 18, fontWeight: 800, color: COLORS.coral}}>
                TAA → CAA
              </span>
              <span
                style={{
                  fontFamily: MONO,
                  fontSize: 12,
                  fontWeight: 800,
                  color: COLORS.paper,
                  background: COLORS.coral,
                  borderRadius: 6,
                  padding: '4px 7px',
                  transform: `rotate(${Math.sin(frame / 7) * 2}deg)`,
                }}
              >
                STOP 坏了
              </span>
            </div>

            <div
              style={{
                marginTop: 22,
                position: 'relative',
                height: 92,
                borderRadius: 16,
                background: 'rgba(221,243,238,0.55)',
                border: `1px solid ${COLORS.teal}22`,
                overflow: 'hidden',
              }}
            >
              <div
                style={{
                  position: 'absolute',
                  left: 18,
                  right: 18,
                  top: 44,
                  height: 8,
                  borderRadius: 8,
                  background: 'rgba(45,42,50,0.12)',
                }}
              />
              <div
                style={{
                  position: 'absolute',
                  left: 18,
                  top: 26,
                  fontSize: 11,
                  fontFamily: MONO,
                  color: COLORS.red,
                  fontWeight: 800,
                }}
              >
                original stop
              </div>
              <div
                style={{
                  position: 'absolute',
                  left: 18,
                  right: 18,
                  top: 63,
                  height: 8,
                  borderRadius: 8,
                  background: 'rgba(42,157,143,0.3)',
                }}
              />
              <div
                style={{
                  position: 'absolute',
                  left: 18,
                  bottom: 4,
                  fontSize: 11,
                  fontFamily: MONO,
                  color: COLORS.teal,
                  fontWeight: 800,
                }}
              >
                3′ UTR · keep translating
              </div>
              <div
                style={{
                  position: 'absolute',
                  left: `${15 + travel * 68}%`,
                  top: 33,
                  width: 36,
                  height: 36,
                  borderRadius: '50%',
                  background: COLORS.gold,
                  border: '5px solid white',
                  boxShadow: '0 8px 18px rgba(244,201,93,0.4)',
                  transform: 'translateX(-50%)',
                }}
              />
            </div>
            <div
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                marginTop: 10,
                fontSize: 11,
                fontFamily: MONO,
                color: COLORS.muted,
              }}
            >
              <span>annotated stop lost</span>
              <span>next in-frame stop?</span>
            </div>
          </div>

          <div style={{flex: 1, ...fadeUp(outcome, 16)}}>
            <div style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted, marginBottom: 12}}>
              two honest outcomes
            </div>
            <div style={{display: 'flex', flexDirection: 'column', gap: 14}}>
              <div
                style={{
                  ...card,
                  padding: '18px 20px',
                  background: 'rgba(221,243,238,0.75)',
                  border: `1px solid ${COLORS.teal}33`,
                }}
              >
                <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
                  <span style={{fontFamily: MONO, fontSize: 12, fontWeight: 800, color: COLORS.green}}>
                    next stop found
                  </span>
                  <span style={{fontFamily: MONO, fontSize: 10, color: COLORS.green}}>
                    protein_ext_stop=true
                  </span>
                </div>
                <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.55, marginTop: 8}}>
                  继续翻译，到下一个 in-frame stop 停下。蛋白比参考更长。
                </div>
              </div>
              <div
                style={{
                  ...card,
                  padding: '18px 20px',
                  background: 'rgba(255,225,219,0.65)',
                  border: `1px solid ${COLORS.red}33`,
                }}
              >
                <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
                  <span style={{fontFamily: MONO, fontSize: 12, fontWeight: 800, color: COLORS.red}}>
                    no stop found
                  </span>
                  <span style={{fontFamily: MONO, fontSize: 10, color: COLORS.red}}>
                    stop_lost=true
                  </span>
                </div>
                <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.55, marginTop: 8}}>
                  一直翻译到可用序列终点，明确标记「没有正常终止」。
                </div>
              </div>
            </div>
            <div
              style={{
                marginTop: 16,
                fontSize: 12,
                lineHeight: 1.6,
                color: COLORS.muted,
                fontFamily: MONO,
              }}
            >
              默认：剪接 3′ UTR，最多 5 kb
              <br />
              不会无限翻译整条染色体
            </div>
            <div style={{marginTop: 16, display: 'flex', justifyContent: 'center'}}>
              <ArrowLine width={70} delay={84} color={COLORS.coral} />
            </div>
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
