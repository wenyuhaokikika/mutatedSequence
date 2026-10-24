import React from 'react';
import {AbsoluteFill, interpolate, useCurrentFrame} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {DnaIcon, SceneShell, Sparkle, fadeUp, useEnter} from '../components/ui';
import type {SceneProps} from './types';

export const OutroScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const frame = useCurrentFrame();
  const logo = useEnter(2, 14);
  const title = useEnter(12, 18);
  const links = useEnter(30);
  const glow = 0.55 + Math.sin(frame / 13) * 0.08;
  const line = interpolate(frame, [10, 44], [0, 1], {
    extrapolateLeft: 'clamp',
    extrapolateRight: 'clamp',
  });
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.coral}>
      <AbsoluteFill style={{alignItems: 'center', justifyContent: 'center'}}>
        <div style={{display: 'flex', flexDirection: 'column', alignItems: 'center', textAlign: 'center'}}>
          <div
            style={{
              width: 112,
              height: 112,
              borderRadius: 34,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              background: 'rgba(255,253,252,0.8)',
              border: `1px solid ${COLORS.line}`,
              boxShadow: `0 20px 60px rgba(255,111,97,${glow * 0.3})`,
              transform: `scale(${logo}) rotate(${(1 - logo) * -8}deg)`,
            }}
          >
            <DnaIcon size={68} color={COLORS.coral} />
          </div>
          <h1
            style={{
              margin: '24px 0 0',
              fontFamily: FONT,
              fontSize: 74,
              letterSpacing: -3.5,
              fontWeight: 900,
              color: COLORS.ink,
              ...fadeUp(title, 20),
            }}
          >
            把 VCF 变成
            <span style={{color: COLORS.teal}}>完整、可解释的蛋白</span>
          </h1>
          <div
            style={{
              fontFamily: MONO,
              fontSize: 22,
              fontWeight: 800,
              letterSpacing: 1.5,
              color: COLORS.coral,
              marginTop: 12,
              opacity: title,
            }}
          >
            mutatedSequence
          </div>
          <div
            style={{
              width: 280,
              height: 6,
              borderRadius: 6,
              background: `linear-gradient(90deg, ${COLORS.coral}, ${COLORS.gold})`,
              marginTop: 20,
              transform: `scaleX(${line})`,
              opacity: 0.7,
            }}
          />
          <div
            style={{
              display: 'flex',
              gap: 12,
              marginTop: 30,
              ...fadeUp(links, 16),
            }}
          >
            <span
              style={{
                ...card,
                padding: '12px 18px',
                borderRadius: 12,
                fontFamily: MONO,
                fontSize: 14,
                fontWeight: 800,
                color: COLORS.ink,
                background: 'rgba(255,253,252,0.86)',
              }}
            >
              github.com/wenyuhaokikika/mutatedSequence
            </span>
            <span
              style={{
                ...card,
                padding: '12px 18px',
                borderRadius: 12,
                fontFamily: MONO,
                fontSize: 14,
                fontWeight: 800,
                color: COLORS.teal,
                background: COLORS.tealSoft,
                border: `1px solid ${COLORS.teal}33`,
              }}
            >
              Docker Hub · wenyuhaobuemoji/mutatedsequence
            </span>
          </div>
          <div
            style={{
              display: 'flex',
              gap: 18,
              marginTop: 26,
              fontSize: 13,
              color: COLORS.muted,
              ...fadeUp(links, 10),
            }}
          >
            <span>open source</span>
            <span>·</span>
            <span>10 tests</span>
            <span>·</span>
            <span>chr21 validated</span>
            <span>·</span>
            <span>no runtime Python deps</span>
          </div>
        </div>
        <Sparkle size={26} style={{left: 260, top: 160}} color={COLORS.coral} />
        <Sparkle size={16} style={{right: 320, top: 210}} color={COLORS.teal} />
        <Sparkle size={20} style={{left: 380, bottom: 150}} color={COLORS.gold} />
        <Sparkle size={14} style={{right: 380, bottom: 190}} color={COLORS.lavender} />
      </AbsoluteFill>
    </SceneShell>
  );
};
