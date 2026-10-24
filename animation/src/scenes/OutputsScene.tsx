import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import {SampleTag, SequenceStrip} from '../components/visuals';
import type {SceneProps} from './types';

const rows = [
  {variant: 'A_missense', refalt: 'C → A', allele: '0/1', effect: 'p.Leu11Met'},
  {variant: 'A_nonsense', refalt: 'A → T', allele: '1/1', effect: 'p.Lys30*'},
  {variant: 'B_missense', refalt: 'G → T', allele: '1/1', effect: 'p.Leu11Met'},
  {variant: 'G_multi', refalt: 'C → A', allele: '1/2', effect: 'p.Gln95Lys'},
];

export const OutputsScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const table = useEnter(10);
  const side = useEnter(32);
  const stamp = useEnter(52);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.gold}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.green} bg={COLORS.goldSoft} style={{marginBottom: 14}}>
              step 05 · audit trail
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
              每个变异，
              <span style={{color: COLORS.coral}}>都能对上号</span>
            </h2>
          </div>
          <div style={{fontSize: 15, color: COLORS.muted, lineHeight: 1.6, textAlign: 'right'}}>
            REF/ALT、allele、consequence、
            <br />
            before/after sequence，一条 TSV 里。
          </div>
        </div>

        <div style={{display: 'flex', gap: 30, marginTop: 34, alignItems: 'flex-start'}}>
          <div
            style={{
              width: 1030,
              ...card,
              padding: '20px 24px',
              background: 'rgba(255,253,252,0.92)',
              ...fadeUp(table, 20),
            }}
          >
            <div
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                borderBottom: `1px solid ${COLORS.line}`,
                paddingBottom: 12,
                marginBottom: 4,
              }}
            >
              <span style={{fontFamily: MONO, fontSize: 14, fontWeight: 800}}>
                mutation_report.tsv
              </span>
              <SampleTag text="one row per variant × transcript" color={COLORS.green} />
            </div>
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1.1fr 1fr 0.6fr 1.2fr',
                gap: 10,
                padding: '10px 4px',
                fontFamily: MONO,
                fontSize: 11,
                color: COLORS.muted,
                textTransform: 'uppercase',
                letterSpacing: 0.5,
              }}
            >
              <span>variant</span>
              <span>REF/ALT</span>
              <span>allele</span>
              <span>protein</span>
            </div>
            {rows.map((row, i) => (
              <div
                key={row.variant}
                style={{
                  display: 'grid',
                  gridTemplateColumns: '1.1fr 1fr 0.6fr 1.2fr',
                  gap: 10,
                  padding: '11px 4px',
                  borderTop: `1px dashed ${COLORS.line}`,
                  fontFamily: MONO,
                  fontSize: 14,
                  color: COLORS.ink,
                  background:
                    i === 0
                      ? 'rgba(255,183,165,0.12)'
                      : i === 1
                        ? 'rgba(244,201,93,0.1)'
                        : 'transparent',
                }}
              >
                <span style={{fontWeight: 700}}>{row.variant}</span>
                <span style={{color: COLORS.coral, fontWeight: 800}}>{row.refalt}</span>
                <span>{row.allele}</span>
                <span style={{color: COLORS.green, fontWeight: 800}}>{row.effect}</span>
              </div>
            ))}
            <div
              style={{
                marginTop: 16,
                padding: '12px 14px',
                borderRadius: 12,
                background: COLORS.lilac,
                fontFamily: MONO,
                fontSize: 12,
                color: COLORS.muted,
                lineHeight: 1.55,
              }}
            >
              transcript_ref → transcript_alt · cds_ref → cds_alt · protein_ref → protein_alt
            </div>
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: 12,
                marginTop: 16,
              }}
            >
              <SequenceStrip
                label="before"
                sequence="ATGGCACGTAGC"
                highlight={5}
                width={330}
                color={COLORS.teal}
                delay={38}
              />
              <span style={{fontSize: 20, color: COLORS.coral, fontWeight: 800}}>→</span>
              <SequenceStrip
                label="after"
                sequence="ATGGCCCGTAGC"
                highlight={5}
                width={330}
                color={COLORS.coral}
                delay={44}
              />
            </div>
          </div>

          <div style={{flex: 1, ...fadeUp(side, 16)}}>
            <div style={{display: 'flex', flexDirection: 'column', gap: 12}}>
              {[
                {name: 'transcripts.alt.fa', detail: 'cDNA + CDS + protein', color: COLORS.teal},
                {name: 'normalized.vcf', detail: 'one ALT, explicit REF/ALT', color: COLORS.coral},
                {name: 'run_summary.json', detail: 'warnings + counts', color: COLORS.gold},
              ].map((f) => (
                <div
                  key={f.name}
                  style={{
                    ...card,
                    padding: '14px 16px',
                    display: 'flex',
                    alignItems: 'center',
                    gap: 12,
                    background: 'rgba(255,253,252,0.8)',
                  }}
                >
                  <span style={{width: 9, height: 9, borderRadius: 9, background: f.color}} />
                  <div>
                    <div style={{fontFamily: MONO, fontSize: 12, fontWeight: 800}}>{f.name}</div>
                    <div style={{fontSize: 12, color: COLORS.muted, marginTop: 3}}>{f.detail}</div>
                  </div>
                </div>
              ))}
            </div>
            <div
              style={{
                marginTop: 24,
                display: 'flex',
                justifyContent: 'center',
                ...fadeUp(stamp, 12),
              }}
            >
              <div
                style={{
                  border: `2px solid ${COLORS.coral}`,
                  color: COLORS.coral,
                  borderRadius: 14,
                  padding: '12px 20px',
                  fontFamily: MONO,
                  fontSize: 16,
                  fontWeight: 900,
                  letterSpacing: 2,
                  transform: `rotate(${stamp * -7 - 2}deg)`,
                  background: 'rgba(255,253,252,0.55)',
                }}
              >
                AUDITABLE
              </div>
            </div>
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
