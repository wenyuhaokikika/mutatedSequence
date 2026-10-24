import React from 'react';
import {interpolate, spring, useCurrentFrame, useVideoConfig} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {clamp, fadeUp, useEnter} from './ui';

const AA_COLORS = [COLORS.teal, COLORS.lavender, COLORS.gold, COLORS.peach];

export const ProteinChain: React.FC<{
  letters?: string;
  highlight?: number;
  delay?: number;
  size?: number;
  gap?: number;
}> = ({
  letters = 'MALWMRLLPLLALLALWGPDPAAA',
  highlight = 7,
  delay = 10,
  size = 44,
  gap = 8,
}) => {
  const frame = useCurrentFrame();
  const {fps} = useVideoConfig();
  return (
    <div style={{display: 'flex', alignItems: 'center', gap}}>
      {letters.split('').map((aa, i) => {
        const p = spring({
          frame: frame - delay - i * 2,
          fps,
          config: {damping: 15, mass: 0.6, stiffness: 170},
        });
        const isMutation = i === highlight;
        return (
          <React.Fragment key={i}>
            <div
              style={{
                width: size,
                height: size,
                borderRadius: '50%',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontFamily: MONO,
                fontWeight: 800,
                fontSize: size * 0.38,
                color: isMutation ? COLORS.paper : COLORS.ink,
                background: isMutation ? COLORS.coral : AA_COLORS[i % AA_COLORS.length],
                border: isMutation ? `3px solid ${COLORS.coral}` : `1px solid ${COLORS.line}`,
                boxShadow: isMutation
                  ? '0 8px 22px rgba(255,111,97,0.32)'
                  : '0 5px 12px rgba(45,42,50,0.06)',
                transform: `scale(${p})`,
                opacity: p,
              }}
            >
              {aa}
            </div>
            {i < letters.length - 1 ? (
              <div
                style={{
                  width: 9,
                  height: 3,
                  borderRadius: 3,
                  background: 'rgba(45,42,50,0.2)',
                  transform: `scaleX(${p})`,
                  transformOrigin: 'left',
                }}
              />
            ) : null}
          </React.Fragment>
        );
      })}
    </div>
  );
};

export const ChromosomeRibbon: React.FC<{
  width?: number;
  highlight?: number;
  delay?: number;
  label?: string;
  markerLabel?: string;
}> = ({width = 660, highlight = 0.42, delay = 0, label = 'chr1', markerLabel = 'REF'}) => {
  const frame = useCurrentFrame();
  const p = useEnter(delay);
  const pulse = 1 + Math.sin(frame / 11) * 0.08;
  const bands = [
    {left: '8%', w: '9%', c: COLORS.peach, o: 0.55},
    {left: '22%', w: '5%', c: COLORS.teal, o: 0.35},
    {left: '31%', w: '13%', c: COLORS.lavender, o: 0.55},
    {left: '51%', w: '6%', c: COLORS.gold, o: 0.5},
    {left: '62%', w: '15%', c: COLORS.peach, o: 0.42},
    {left: '83%', w: '7%', c: COLORS.teal, o: 0.45},
  ];
  return (
    <div style={{width, opacity: p, transform: `translateY(${(1 - p) * 16}px)`}}>
      <div
        style={{
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
          marginBottom: 12,
        }}
      >
        <span style={{fontFamily: MONO, fontSize: 16, color: COLORS.muted}}>{label}</span>
        <span style={{fontSize: 13, color: COLORS.muted}}>reference genome</span>
      </div>
      <div
        style={{
          position: 'relative',
          height: 66,
          borderRadius: 33,
          overflow: 'hidden',
          background:
            'linear-gradient(90deg, rgba(255,183,165,0.55), rgba(183,228,199,0.5), rgba(205,180,219,0.55))',
          border: `1px solid ${COLORS.line}`,
          boxShadow: COLORS.shadowSoft,
        }}
      >
        {bands.map((band, i) => (
          <div
            key={i}
            style={{
              position: 'absolute',
              left: band.left,
              top: i % 2 ? 12 : 0,
              width: band.w,
              height: i % 2 ? 42 : 22,
              background: band.c,
              opacity: band.o,
              borderRadius: 8,
            }}
          />
        ))}
        <div
          style={{
            position: 'absolute',
            left: `${highlight * 100}%`,
            top: '50%',
            width: 18,
            height: 18,
            borderRadius: '50%',
            background: COLORS.coral,
            border: '4px solid white',
            boxShadow: '0 0 0 7px rgba(255,111,97,0.15)',
            transform: `translate(-50%, -50%) scale(${pulse})`,
          }}
        />
      </div>
      <div style={{display: 'flex', justifyContent: 'flex-end', marginTop: 9}}>
        <span
          style={{
            fontSize: 13,
            color: COLORS.coral,
            fontWeight: 700,
            letterSpacing: 0.5,
          }}
        >
          {markerLabel}
        </span>
      </div>
    </div>
  );
};

export const FileCard: React.FC<{
  kind: 'dna' | 'gtf' | 'vcf' | 'protein';
  title: string;
  subtitle: string;
  color: string;
  delay?: number;
  compact?: boolean;
}> = ({kind, title, subtitle, color, delay = 0, compact = false}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        ...card,
        width: compact ? 230 : 300,
        padding: compact ? '18px 20px' : '24px 26px',
        background: 'rgba(255,253,252,0.9)',
        ...fadeUp(p, 24),
      }}
    >
      <div
        style={{
          width: compact ? 38 : 46,
          height: compact ? 38 : 46,
          borderRadius: 14,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          background: `${color}22`,
          marginBottom: compact ? 12 : 16,
        }}
      >
        {kind === 'dna' ? (
          <svg width={compact ? 22 : 26} height={compact ? 22 : 26} viewBox="0 0 24 24" fill="none">
            <path d="M7 3c6 4 4 14-2 18M17 3c-6 4-4 14 2 18M6 8h12M6 16h12" stroke={color} strokeWidth={2} strokeLinecap="round" />
          </svg>
        ) : null}
        {kind === 'gtf' ? (
          <svg width={compact ? 22 : 26} height={compact ? 22 : 26} viewBox="0 0 24 24" fill="none">
            <rect x="4" y="3" width="16" height="18" rx="2" stroke={color} strokeWidth={2} />
            <path d="M8 8h8M8 12h8M8 16h5" stroke={color} strokeWidth={2} strokeLinecap="round" />
          </svg>
        ) : null}
        {kind === 'vcf' ? (
          <svg width={compact ? 22 : 26} height={compact ? 22 : 26} viewBox="0 0 24 24" fill="none">
            <rect x="3" y="4" width="18" height="16" rx="2" stroke={color} strokeWidth={2} />
            <path d="M3 9h18M8 4v16M13 4v16" stroke={color} strokeWidth={1.5} opacity={0.7} />
            <path d="M10 6l1.2 1.2L13 5.4" stroke={color} strokeWidth={2} strokeLinecap="round" />
          </svg>
        ) : null}
        {kind === 'protein' ? (
          <svg width={compact ? 22 : 26} height={compact ? 22 : 26} viewBox="0 0 24 24" fill="none">
            <circle cx="7" cy="9" r="3" stroke={color} strokeWidth={2} />
            <circle cx="16" cy="7" r="2.5" stroke={color} strokeWidth={2} />
            <circle cx="15" cy="16" r="3.5" stroke={color} strokeWidth={2} />
            <path d="M9.5 10.5l3 3M13.8 8.6l-1.2 4.5" stroke={color} strokeWidth={1.8} strokeLinecap="round" />
          </svg>
        ) : null}
      </div>
      <div style={{fontFamily: MONO, fontSize: compact ? 14 : 16, fontWeight: 800, color: COLORS.ink}}>
        {title}
      </div>
      <div style={{fontSize: compact ? 13 : 14, lineHeight: 1.45, color: COLORS.muted, marginTop: 6}}>
        {subtitle}
      </div>
    </div>
  );
};

export const Base: React.FC<{
  value: string;
  color?: string;
  bg?: string;
  size?: number;
  strike?: boolean;
}> = ({value, color = COLORS.ink, bg = COLORS.coralSoft, size = 32, strike = false}) => (
  <span
    style={{
      display: 'inline-flex',
      alignItems: 'center',
      justifyContent: 'center',
      width: size,
      height: size,
      borderRadius: 9,
      background: bg,
      color,
      fontFamily: MONO,
      fontWeight: 800,
      fontSize: size * 0.46,
      textDecoration: strike ? 'line-through' : 'none',
      border: `1px solid ${COLORS.line}`,
    }}
  >
    {value}
  </span>
);

export const VcfCard: React.FC<{
  refBase: string;
  altBase: string;
  gt: string;
  delay?: number;
  highlight?: boolean;
  mismatch?: boolean;
}> = ({refBase, altBase, gt, delay = 0, highlight = false, mismatch = false}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        ...card,
        padding: '18px 20px',
        width: 290,
        background: 'rgba(255,253,252,0.94)',
        border: highlight
          ? `1.5px solid ${mismatch ? COLORS.red : COLORS.coral}`
          : `1px solid ${COLORS.line}`,
        ...fadeUp(p, 18),
      }}
    >
      <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 14}}>
        <span style={{fontFamily: MONO, fontSize: 13, color: COLORS.muted}}>sample.vcf</span>
        <span
          style={{
            fontFamily: MONO,
            fontSize: 12,
            fontWeight: 800,
            padding: '4px 8px',
            borderRadius: 7,
            background: mismatch ? 'rgba(196,69,54,0.1)' : COLORS.tealSoft,
            color: mismatch ? COLORS.red : COLORS.green,
          }}
        >
          {gt}
        </span>
      </div>
      <div style={{display: 'flex', alignItems: 'center', gap: 10}}>
        <span style={{fontSize: 12, color: COLORS.muted, width: 34}}>REF</span>
        <Base value={refBase} bg={mismatch ? COLORS.coralSoft : COLORS.tealSoft} color={mismatch ? COLORS.red : COLORS.green} />
        <span style={{fontSize: 16, color: COLORS.muted}}>→</span>
        <span style={{fontSize: 12, color: COLORS.muted, width: 34}}>ALT</span>
        <Base value={altBase} bg={COLORS.goldSoft} color={COLORS.ink} />
      </div>
      <div
        style={{
          marginTop: 14,
          paddingTop: 12,
          borderTop: `1px dashed ${COLORS.line}`,
          fontSize: 13,
          color: COLORS.muted,
          lineHeight: 1.45,
        }}
      >
        {mismatch ? 'REF 对不上，先别急着翻译' : '位置、REF、ALT、GT 一起看'}
      </div>
    </div>
  );
};

export const SequenceStrip: React.FC<{
  sequence: string;
  highlight?: number;
  color?: string;
  label?: string;
  delay?: number;
  width?: number;
}> = ({sequence, highlight = -1, color = COLORS.teal, label, delay = 0, width}) => {
  const p = useEnter(delay);
  const chars = sequence.split('');
  return (
    <div style={{width, opacity: p}}>
      {label ? (
        <div style={{fontSize: 12, color: COLORS.muted, marginBottom: 7, fontFamily: MONO}}>{label}</div>
      ) : null}
      <div style={{display: 'flex', flexWrap: 'wrap', gap: 4}}>
        {chars.map((ch, i) => (
          <span
            key={i}
            style={{
              width: 22,
              height: 28,
              borderRadius: 5,
              display: 'inline-flex',
              alignItems: 'center',
              justifyContent: 'center',
              fontFamily: MONO,
              fontSize: 13,
              fontWeight: 700,
              color: i === highlight ? COLORS.paper : COLORS.ink,
              background:
                i === highlight
                  ? COLORS.coral
                  : i % 8 < 4
                    ? `${color}1c`
                    : 'rgba(255,253,252,0.7)',
              border: `1px solid ${i === highlight ? COLORS.coral : COLORS.line}`,
              transform: `translateY(${(1 - p) * 8}px)`,
            }}
          >
            {ch}
          </span>
        ))}
      </div>
    </div>
  );
};

export const ArrowLine: React.FC<{
  width?: number;
  color?: string;
  dashed?: boolean;
  delay?: number;
  vertical?: boolean;
}> = ({width = 90, color = COLORS.coral, dashed = false, delay = 0, vertical = false}) => {
  const frame = useCurrentFrame();
  const p = interpolate(frame - delay, [0, 16], [0, 1], clamp);
  return (
    <svg
      width={vertical ? 24 : width}
      height={vertical ? width : 24}
      viewBox={vertical ? `0 0 24 ${width}` : `0 0 ${width} 24`}
      fill="none"
      style={{opacity: p}}
    >
      <path
        d={vertical ? `M12 0v${width - 10}` : `M0 12h${width - 10}`}
        stroke={color}
        strokeWidth={3}
        strokeLinecap="round"
        strokeDasharray={dashed ? '7 8' : undefined}
      />
      <path
        d={vertical ? `M5 ${width - 14}l7 9 7-9` : `M${width - 14} 5l9 7-9 7`}
        stroke={color}
        strokeWidth={3}
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
};

export const SampleTag: React.FC<{text: string; color?: string}> = ({text, color = COLORS.teal}) => (
  <span
    style={{
      display: 'inline-flex',
      alignItems: 'center',
      gap: 6,
      fontFamily: MONO,
      fontSize: 13,
      fontWeight: 700,
      color,
      background: `${color}14`,
      border: `1px solid ${color}33`,
      borderRadius: 7,
      padding: '4px 8px',
    }}
  >
    {text}
  </span>
);

export const MiniTable: React.FC<{
  rows: Array<{left: string; right: string; accent?: boolean}>;
  width?: number;
  delay?: number;
}> = ({rows, width = 420, delay = 0}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        width,
        ...card,
        padding: '16px 18px',
        ...fadeUp(p, 18),
      }}
    >
      {rows.map((row, i) => (
        <div
          key={i}
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            gap: 20,
            padding: '9px 0',
            borderBottom: i < rows.length - 1 ? `1px dashed ${COLORS.line}` : 'none',
            fontFamily: MONO,
            fontSize: 13,
            color: row.accent ? COLORS.coral : COLORS.ink,
            fontWeight: row.accent ? 800 : 600,
          }}
        >
          <span style={{color: COLORS.muted, fontWeight: 600}}>{row.left}</span>
          <span style={{textAlign: 'right'}}>{row.right}</span>
        </div>
      ))}
    </div>
  );
};

export const bigTitle: React.CSSProperties = {
  fontFamily: FONT,
  fontWeight: 900,
  fontSize: 72,
  lineHeight: 1.15,
  letterSpacing: -2.2,
  color: COLORS.ink,
  margin: 0,
};
