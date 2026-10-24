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
              step 04 · transcript → protein
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
            annotated protein 停在注释 stop；
            <br />
            3′ extension 负责 readthrough。
          </div>
        </div>

        <div style={{display: 'flex', gap: 20, marginTop: 30, alignItems: 'stretch'}}>
          <StageCard step="exons" title="cut & splice" color={COLORS.coral} delay={10}>
            <div style={{display: 'flex', flexDirection: 'column', gap: 10, marginTop: 32}}>
              <div style={{display: 'flex', gap: 6, alignItems: 'center'}}>
                <ExonMini label="E1" color={COLORS.coral} width={84} />
                <div style={{height: 2, flex: 1, background: 'rgba(45,42,50,0.18)'}} />
                <ExonMini label="E2" color={COLORS.teal} width={92} />
                <div style={{height: 2, flex: 1, background: 'rgba(45,42,50,0.18)'}} />
                <ExonMini label="E3" color={COLORS.lavender} width={94} />
              </div>
              <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 14}}>
                按 GTF 的 exon 顺序切取；负链整体反向互补。
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

          <StageCard step="cDNA" title="all exons" color={COLORS.teal} delay={18}>
            <SequenceStrip
              label="transcript · 5′ → 3′"
              sequence="ATGGCACGTACGTAGCTAGCTAGC"
              color={COLORS.teal}
              highlight={5}
              delay={22}
            />
            <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 16}}>
              完整 transcript：UTR、 exon、 exon junction 都在。
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
              cDNA = exon concatenation
            </div>
          </StageCard>

          <StageCard step="CDS" title="frame matters" color={COLORS.gold} delay={26}>
            <SequenceStrip
              label="coding sequence"
              sequence="NNATGGCACGTACGTAGCTAGCTAA"
              color={COLORS.gold}
              highlight={2}
              delay={30}
            />
            <div style={{display: 'flex', gap: 6, marginTop: 14, flexWrap: 'wrap'}}>
              <span style={{fontSize: 11, color: COLORS.coral, fontFamily: MONO, fontWeight: 800}}>
                5′ N pad
              </span>
              <span style={{fontSize: 11, color: COLORS.muted}}>·</span>
              <span style={{fontSize: 11, color: COLORS.green, fontFamily: MONO, fontWeight: 800}}>
                stop TAA
              </span>
            </div>
            <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 14}}>
              GTF CDS 不含 stop，就补上 stop_codon；partial 5′ 保留 N pad。
            </div>
          </StageCard>

          <StageCard step="protein" title="translate" color={COLORS.lavender} delay={34}>
            <div style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted, marginBottom: 8}}>
              annotated
            </div>
            <ProteinChain letters="MALWMRLLP" highlight={4} size={27} gap={4} delay={38} />
            <div style={{fontFamily: MONO, fontSize: 12, color: COLORS.muted, marginTop: 18, marginBottom: 8}}>
              3′ extended
            </div>
            <ProteinChain letters="MALWMRLLP" highlight={4} size={27} gap={4} delay={44} />
            <div style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.55, marginTop: 14}}>
              stop-lost / frameshift 时，继续翻译到下一个 in-frame stop。
            </div>
          </StageCard>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
