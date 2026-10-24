import React from 'react';
import {
  AbsoluteFill,
  Easing,
  interpolate,
  spring,
  useCurrentFrame,
  useVideoConfig,
} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import type {SceneTiming} from '../types';

export const useEnter = (delay = 0, damping = 16) => {
  const frame = useCurrentFrame();
  const {fps} = useVideoConfig();
  return spring({
    frame: frame - delay,
    fps,
    config: {damping, mass: 0.85, stiffness: 120},
  });
};

export const fadeUp = (p: number, distance = 26): React.CSSProperties => ({
  opacity: p,
  transform: `translateY(${(1 - p) * distance}px)`,
});

export const clamp = {
  extrapolateLeft: 'clamp' as const,
  extrapolateRight: 'clamp' as const,
};

export const ProgressDots: React.FC<{active: number; total: number}> = ({
  active,
  total,
}) => (
  <div style={{display: 'flex', gap: 8, alignItems: 'center'}}>
    {Array.from({length: total}).map((_, i) => (
      <div
        key={i}
        style={{
          width: i === active ? 26 : 8,
          height: 8,
          borderRadius: 8,
          background:
            i === active
              ? COLORS.coral
              : i < active
                ? COLORS.peach
                : 'rgba(45,42,50,0.16)',
          transition: 'none',
        }}
      />
    ))}
  </div>
);

export const DnaIcon: React.FC<{size?: number; color?: string}> = ({
  size = 34,
  color = COLORS.coral,
}) => (
  <svg width={size} height={size} viewBox="0 0 64 64" fill="none">
    <path
      d="M18 8c16 10 16 38 0 48M46 8c-16 10-16 38 0 48"
      stroke={color}
      strokeWidth={4}
      strokeLinecap="round"
    />
    <path
      d="M21 18h22M19 28h26M19 38h26M21 48h22"
      stroke={COLORS.teal}
      strokeWidth={3}
      strokeLinecap="round"
      opacity={0.85}
    />
  </svg>
);

export const Wordmark: React.FC<{compact?: boolean}> = ({compact = false}) => (
  <div style={{display: 'flex', alignItems: 'center', gap: 10}}>
    <DnaIcon size={compact ? 24 : 28} />
    <span
      style={{
        fontFamily: MONO,
        fontSize: compact ? 16 : 18,
        fontWeight: 700,
        letterSpacing: -0.4,
        color: COLORS.ink,
        opacity: 0.82,
      }}
    >
      mutatedSequence
    </span>
  </div>
);

const FloatingBases: React.FC<{accent: string}> = ({accent}) => {
  const frame = useCurrentFrame();
  const {width, height} = useVideoConfig();
  const bases = [
    {b: 'A', x: 0.07, y: 0.2, size: 30, delay: 0, speed: 47},
    {b: 'C', x: 0.91, y: 0.26, size: 22, delay: 20, speed: 53},
    {b: 'G', x: 0.13, y: 0.76, size: 24, delay: 8, speed: 61},
    {b: 'T', x: 0.87, y: 0.72, size: 34, delay: 15, speed: 49},
    {b: 'A', x: 0.28, y: 0.1, size: 18, delay: 26, speed: 57},
    {b: 'G', x: 0.72, y: 0.88, size: 20, delay: 4, speed: 44},
    {b: 'T', x: 0.48, y: 0.06, size: 16, delay: 12, speed: 63},
    {b: 'C', x: 0.58, y: 0.94, size: 18, delay: 22, speed: 51},
  ];
  return (
    <AbsoluteFill style={{pointerEvents: 'none'}}>
      {bases.map((item, i) => {
        const y = Math.sin((frame + item.delay) / item.speed) * 16;
        const x = Math.cos((frame + item.delay) / (item.speed * 1.7)) * 9;
        return (
          <span
            key={i}
            style={{
              position: 'absolute',
              left: item.x * width,
              top: item.y * height,
              transform: `translate(${x}px, ${y}px) rotate(${y / 5}deg)`,
              fontFamily: MONO,
              fontSize: item.size,
              fontWeight: 700,
              color: i % 2 ? accent : COLORS.teal,
              opacity: 0.2,
            }}
          >
            {item.b}
          </span>
        );
      })}
    </AbsoluteFill>
  );
};

const Blob: React.FC<{
  left: number;
  top: number;
  size: number;
  color: string;
  speed: number;
  depth: number;
}> = ({left, top, size, color, speed, depth}) => {
  const frame = useCurrentFrame();
  const {width, height} = useVideoConfig();
  const drift = Math.sin(frame / speed) * 18;
  return (
    <div
      style={{
        position: 'absolute',
        left: left * width - size / 2,
        top: top * height - size / 2 + drift,
        width: size,
        height: size,
        borderRadius: '50%',
        background: color,
        filter: `blur(${22 + depth * 16}px)`,
        opacity: 0.42,
      }}
    />
  );
};

const Backdrop: React.FC<{accent: string}> = ({accent}) => {
  const frame = useCurrentFrame();
  return (
    <AbsoluteFill
      style={{
        background:
          'linear-gradient(135deg, #FFF9F2 0%, #FFF5EB 46%, #F7F0FF 100%)',
        overflow: 'hidden',
      }}
    >
      <div
        style={{
          position: 'absolute',
          inset: 0,
          backgroundImage:
            'linear-gradient(rgba(45,42,50,0.035) 1px, transparent 1px), linear-gradient(90deg, rgba(45,42,50,0.035) 1px, transparent 1px)',
          backgroundSize: '72px 72px',
          maskImage:
            'radial-gradient(circle at 50% 42%, black 0%, transparent 75%)',
        }}
      />
      <Blob left={0.82} top={0.18} size={420} color={accent} speed={90} depth={1} />
      <Blob left={0.13} top={0.82} size={360} color={COLORS.lavender} speed={110} depth={0} />
      <Blob left={0.52} top={0.06} size={280} color={COLORS.gold} speed={80} depth={0} />
      <div
        style={{
          position: 'absolute',
          inset: 24,
          border: '1px solid rgba(45,42,50,0.08)',
          borderRadius: 34,
        }}
      />
      <FloatingBases accent={accent} />
      <div
        style={{
          position: 'absolute',
          left: 0,
          right: 0,
          bottom: 0,
          height: 140,
          background: 'linear-gradient(transparent, rgba(255,248,240,0.65))',
        }}
      />
      <div
        style={{
          position: 'absolute',
          right: 84,
          top: 110,
          opacity: 0.12,
          transform: `rotate(${8 + Math.sin(frame / 120) * 2}deg)`,
        }}
      >
        <DnaIcon size={120} color={accent} />
      </div>
    </AbsoluteFill>
  );
};

export const SceneShell: React.FC<{
  scene: SceneTiming;
  index: number;
  total: number;
  accent?: string;
  children: React.ReactNode;
}> = ({scene, index, total, accent = COLORS.coral, children}) => {
  const frame = useCurrentFrame();
  const {fps} = useVideoConfig();
  const duration = scene.durationInFrames;
  const fadeIn = interpolate(frame, [0, 9], [0, 1], {
    ...clamp,
    easing: Easing.out(Easing.cubic),
  });
  const fadeOut = interpolate(frame, [duration - 10, duration], [1, 0], clamp);
  const opacity = Math.min(fadeIn, fadeOut);
  const caption = spring({
    frame: frame - 7,
    fps,
    config: {damping: 20, mass: 0.7, stiffness: 110},
  });
  const progress = ((index + Math.min(frame / duration, 1)) / total) * 100;
  return (
    <AbsoluteFill style={{opacity, color: COLORS.ink, fontFamily: FONT}}>
      <Backdrop accent={accent} />
      <div
        style={{
          position: 'absolute',
          top: 42,
          left: 64,
          right: 64,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
        }}
      >
        <div style={{display: 'flex', alignItems: 'center', gap: 18}}>
          <div
            style={{
              display: 'flex',
              alignItems: 'baseline',
              gap: 5,
              padding: '9px 14px',
              borderRadius: 14,
              background: 'rgba(255,253,252,0.78)',
              border: `1px solid ${COLORS.line}`,
              boxShadow: COLORS.shadowSoft,
            }}
          >
            <span style={{fontFamily: MONO, fontSize: 17, fontWeight: 800, color: accent}}>
              {String(index + 1).padStart(2, '0')}
            </span>
            <span style={{fontSize: 12, color: COLORS.muted}}>/ {String(total).padStart(2, '0')}</span>
          </div>
          <span style={{fontSize: 21, fontWeight: 700, color: COLORS.muted, letterSpacing: 0.2}}>
            {scene.title}
          </span>
        </div>
        <div style={{display: 'flex', alignItems: 'center', gap: 22}}>
          <ProgressDots active={index} total={total} />
          <Wordmark compact />
        </div>
      </div>
      <div
        style={{
          position: 'absolute',
          left: 100,
          right: 100,
          top: 142,
          bottom: 204,
        }}
      >
        {children}
      </div>
      <div
        style={{
          position: 'absolute',
          left: 0,
          right: 0,
          bottom: 42,
          display: 'flex',
          justifyContent: 'center',
          padding: '0 120px',
          opacity: caption,
          transform: `translateY(${(1 - caption) * 18}px)`,
        }}
      >
        <div
          style={{
            maxWidth: 1420,
            padding: '16px 26px',
            ...card,
            borderRadius: 20,
            background: 'rgba(255,253,252,0.82)',
            backdropFilter: 'blur(12px)',
            display: 'flex',
            alignItems: 'flex-start',
            gap: 14,
          }}
        >
          <div
            style={{
              width: 7,
              height: 7,
              borderRadius: 7,
              background: accent,
              marginTop: 12,
              flexShrink: 0,
            }}
          />
          <div
            style={{
              fontSize: 24,
              lineHeight: 1.55,
              color: COLORS.ink,
              letterSpacing: 0.15,
              textWrap: 'pretty',
            }}
          >
            {scene.narration}
          </div>
        </div>
      </div>
      <div
        style={{
          position: 'absolute',
          left: 0,
          right: 0,
          bottom: 0,
          height: 7,
          background: 'rgba(45,42,50,0.08)',
        }}
      >
        <div
          style={{
            width: `${progress}%`,
            height: '100%',
            background: `linear-gradient(90deg, ${COLORS.coral}, ${COLORS.gold})`,
          }}
        />
      </div>
    </AbsoluteFill>
  );
};

export const Pill: React.FC<{
  children: React.ReactNode;
  color?: string;
  bg?: string;
  style?: React.CSSProperties;
}> = ({children, color = COLORS.ink, bg = 'rgba(255,253,252,0.9)', style}) => (
  <span
    style={{
      display: 'inline-flex',
      alignItems: 'center',
      gap: 7,
      padding: '7px 13px',
      borderRadius: 999,
      background: bg,
      color,
      border: `1px solid ${COLORS.line}`,
      fontSize: 15,
      fontWeight: 700,
      whiteSpace: 'nowrap',
      ...style,
    }}
  >
    {children}
  </span>
);

export const Sparkle: React.FC<{
  size?: number;
  color?: string;
  opacity?: number;
  style?: React.CSSProperties;
}> = ({size = 18, color = COLORS.gold, opacity = 0.8, style}) => (
  <svg
    width={size}
    height={size}
    viewBox="0 0 24 24"
    fill="none"
    style={{position: 'absolute', ...style}}
  >
    <path
      d="M12 2l2.1 7.9L22 12l-7.9 2.1L12 22l-2.1-7.9L2 12l7.9-2.1L12 2z"
      fill={color}
      opacity={opacity}
    />
  </svg>
);
