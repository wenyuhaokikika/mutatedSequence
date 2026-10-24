import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {DnaIcon, Pill, SceneShell, Sparkle, fadeUp, useEnter} from '../components/ui';
import {ArrowLine, FileCard} from '../components/visuals';
import type {SceneProps} from './types';

export const InputsScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const engine = useEnter(24, 15);
  const optional = useEnter(36);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.teal}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'center', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 20)}>
            <Pill color={COLORS.teal} bg={COLORS.tealSoft} style={{marginBottom: 15}}>
              输入 · 只要三样
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
              参考基因组、注释、
              <br />
              <span style={{color: COLORS.teal}}>还有你的 VCF</span>
            </h2>
          </div>
          <div
            style={{
              fontSize: 16,
              lineHeight: 1.7,
              color: COLORS.muted,
              maxWidth: 430,
              textAlign: 'right',
            }}
          >
            可选再给 cDNA / CDS / pep FASTA，
            <br />
            用来验证重建结果，不强制。
          </div>
        </div>

        <div style={{display: 'flex', gap: 40, justifyContent: 'center', marginTop: 42}}>
          <FileCard
            kind="dna"
            title="genome.fa"
            subtitle="参考基因组 · 快速随机读取"
            color={COLORS.coral}
            delay={10}
          />
          <FileCard
            kind="gtf"
            title="annotation.gtf"
            subtitle="外显子 / 编码区 / 终止密码子"
            color={COLORS.teal}
            delay={18}
          />
          <FileCard
            kind="vcf"
            title="sample.vcf"
            subtitle="单碱基 / 插入缺失 / 基因型 / 结构变异"
            color={COLORS.gold}
            delay={26}
          />
        </div>

        <div
          style={{
            display: 'flex',
            gap: 260,
            justifyContent: 'center',
            marginTop: 16,
            opacity: engine,
          }}
        >
          <ArrowLine width={44} vertical delay={28} color={COLORS.coral} />
          <ArrowLine width={44} vertical delay={28} color={COLORS.teal} />
          <ArrowLine width={44} vertical delay={28} color={COLORS.gold} />
        </div>

        <div
          style={{
            display: 'flex',
            justifyContent: 'center',
            marginTop: 4,
            ...fadeUp(engine, 16),
          }}
        >
          <div
            style={{
              ...card,
              width: 900,
              padding: '24px 30px',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'space-between',
              background:
                'linear-gradient(120deg, rgba(255,253,252,0.95), rgba(221,243,238,0.8))',
              border: `1px solid ${COLORS.teal}33`,
            }}
          >
            <div style={{display: 'flex', alignItems: 'center', gap: 18}}>
              <div
                style={{
                  width: 54,
                  height: 54,
                  borderRadius: 18,
                  background: COLORS.tealSoft,
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                }}
              >
                <DnaIcon size={34} color={COLORS.teal} />
              </div>
              <div>
                <div style={{fontFamily: MONO, fontSize: 20, fontWeight: 800}}>
                  mutatedSequence
                </div>
                <div style={{fontSize: 15, color: COLORS.muted, marginTop: 5}}>
                  先在基因组上改，再切转录本和蛋白
                </div>
              </div>
            </div>
            <div style={{display: 'flex', gap: 8, flexWrap: 'wrap', justifyContent: 'flex-end', maxWidth: 330}}>
              {['cDNA', '编码区', '正常蛋白', '3′ 延伸'].map((x) => (
                <span
                  key={x}
                  style={{
                    fontFamily: MONO,
                    fontSize: 12,
                    color: COLORS.teal,
                    background: 'rgba(255,255,255,0.75)',
                    border: `1px solid ${COLORS.teal}33`,
                    borderRadius: 7,
                    padding: '5px 8px',
                  }}
                >
                  {x}
                </span>
              ))}
            </div>
          </div>
        </div>

        <div
          style={{
            display: 'flex',
            justifyContent: 'center',
            marginTop: 14,
            gap: 8,
            opacity: optional,
          }}
        >
          <span style={{fontSize: 13, color: COLORS.muted}}>optional validation:</span>
          {['转录本 FASTA', '编码区 FASTA', '蛋白 FASTA'].map((x) => (
            <span
              key={x}
              style={{
                fontFamily: MONO,
                fontSize: 12,
                color: COLORS.muted,
                background: 'rgba(255,253,252,0.7)',
                border: `1px solid ${COLORS.line}`,
                borderRadius: 6,
                padding: '3px 7px',
              }}
            >
              {x}
            </span>
          ))}
        </div>
        <Sparkle size={16} style={{left: 40, bottom: 10}} color={COLORS.teal} />
        <Sparkle size={12} style={{right: 80, top: 18}} color={COLORS.gold} />
      </AbsoluteFill>
    </SceneShell>
  );
};
