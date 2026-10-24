import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {SceneShell, Sparkle, fadeUp, useEnter} from '../components/ui';
import {ArrowLine, ProteinChain, VcfCard} from '../components/visuals';
import type {SceneProps} from './types';

export const HookScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const title = useEnter(2, 18);
  const sub = useEnter(12);
  const left = useEnter(20);
  const right = useEnter(28);
  const underline = interpolate(frame, [18, 42], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });

  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.coral}>
      <AbsoluteFill>
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            height: '100%',
            gap: 50,
          }}
        >
          <div style={{width: 770, ...fadeUp(title, 26)}}>
            <div
              style={{
                display: 'inline-flex',
                alignItems: 'center',
                gap: 10,
                marginBottom: 24,
                ...fadeUp(sub, 12),
              }}
            >
              <span
                style={{
                  width: 10,
                  height: 10,
                  borderRadius: 10,
                  background: COLORS.coral,
                }}
              />
              <span
                style={{
                  fontFamily: MONO,
                  fontSize: 15,
                  fontWeight: 700,
                  letterSpacing: 2,
                  color: COLORS.muted,
                }}
              >
                VCF → TRANSCRIPT → PROTEIN
              </span>
            </div>
            <h1
              style={{
                margin: 0,
                fontFamily: FONT,
                fontSize: 88,
                lineHeight: 1.1,
                letterSpacing: -4,
                fontWeight: 900,
                color: COLORS.ink,
              }}
            >
              把一个
              <span style={{color: COLORS.coral}}>小字母</span>
              <br />
              变成
              <span
                style={{
                  position: 'relative',
                  color: COLORS.teal,
                  whiteSpace: 'nowrap',
                }}
              >
                一整条蛋白
                <span
                  style={{
                    position: 'absolute',
                    left: 4,
                    right: 4,
                    bottom: 3,
                    height: 10,
                    borderRadius: 8,
                    background: COLORS.gold,
                    opacity: 0.5,
                    transform: `scaleX(${underline})`,
                    transformOrigin: 'left',
                  }}
                />
              </span>
            </h1>
            <p
              style={{
                margin: '30px 0 0',
                maxWidth: 620,
                fontSize: 25,
                lineHeight: 1.65,
                color: COLORS.muted,
                ...fadeUp(sub, 14),
              }}
            >
              不是只给一个 <span style={{fontFamily: MONO, color: COLORS.ink}}>p.</span> 注释，
              而是把突变后的 cDNA、CDS 和蛋白，一条一条写清楚。
            </p>
            <div style={{display: 'flex', gap: 10, marginTop: 30, ...fadeUp(sub, 10)}}>
              {['genome-first', 'allele-aware', 'auditable'].map((tag, i) => (
                <span
                  key={tag}
                  style={{
                    fontFamily: MONO,
                    fontSize: 13,
                    color: i === 0 ? COLORS.coral : i === 1 ? COLORS.teal : COLORS.green,
                    background:
                      i === 0
                        ? COLORS.coralSoft
                        : i === 1
                          ? COLORS.tealSoft
                          : COLORS.goldSoft,
                    padding: '7px 11px',
                    borderRadius: 8,
                    fontWeight: 700,
                  }}
                >
                  {tag}
                </span>
              ))}
            </div>
          </div>

          <div
            style={{
              width: 760,
              height: 540,
              position: 'relative',
              ...fadeUp(right, 22),
            }}
          >
            <div
              style={{
                position: 'absolute',
                inset: 20,
                borderRadius: 34,
                background:
                  'linear-gradient(145deg, rgba(255,255,255,0.75), rgba(255,240,235,0.45))',
                border: `1px solid ${COLORS.line}`,
                boxShadow: COLORS.shadow,
              }}
            />
            <div style={{position: 'absolute', left: 44, top: 46, ...fadeUp(left, 18)}}>
              <VcfCard refBase="A" altBase="G" gt="0/1" delay={24} highlight />
            </div>
            <div
              style={{
                position: 'absolute',
                left: 360,
                top: 224,
                transform: 'rotate(-8deg)',
              }}
            >
              <ArrowLine width={130} delay={38} color={COLORS.coral} />
            </div>
            <div
              style={{
                position: 'absolute',
                left: 480,
                top: 110,
                width: 230,
                ...fadeUp(right, 18),
              }}
            >
              <div style={{fontSize: 13, fontFamily: MONO, color: COLORS.muted, marginBottom: 10}}>
                mutated protein
              </div>
              <ProteinChain letters="MALWMRL" highlight={4} size={36} gap={5} delay={34} />
              <div style={{fontSize: 13, fontFamily: MONO, color: COLORS.muted, marginTop: 12}}>
                p.Leu5Arg
              </div>
            </div>
            <div
              style={{
                position: 'absolute',
                left: 480,
                top: 280,
                display: 'flex',
                gap: 8,
                alignItems: 'center',
                ...fadeUp(right, 12),
              }}
            >
              <span style={{fontSize: 14, color: COLORS.muted}}>one letter</span>
              <span style={{fontSize: 22, color: COLORS.coral}}>→</span>
              <span style={{fontSize: 14, color: COLORS.muted}}>whole protein</span>
            </div>
            <Sparkle size={22} style={{left: 58, top: 300}} />
            <Sparkle size={14} style={{right: 48, top: 70}} color={COLORS.coral} />
            <Sparkle size={18} style={{right: 90, bottom: 54}} color={COLORS.teal} />
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
